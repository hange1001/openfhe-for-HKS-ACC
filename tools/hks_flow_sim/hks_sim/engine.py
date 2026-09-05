"""编排：配置 -> trace -> 容量 pre-pass -> 调度 -> 指标。

（方案 §8 的文件清单里没有这一层；把编排从 report.py 里分出来，是为了让
report.py 只负责输出格式，测试可以直接调用 run_all 而不经过 CLI。）
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

from . import trace_dc, trace_mp, trace_oc
from .config import (SimConfig, hardware_config_hash, strategy_config_hash,
                     workload_config_hash)
from .cost_model import CostModel
from .memory import CapacityError, MemoryStats, apply_capacity
from .op import COMPUTE_KINDS, TRANSFER_KINDS, Op, OpKind
from .resources import BCONV, TRANSFORM, ResourceSet, build_resources
from .scheduler import ScheduleResult, schedule
from .trace_common import TraceContext
from .workload import Workload, build_workload

STRATEGY_BUILDERS: Dict[str, Callable[[TraceContext], List[Op]]] = {
    "dc": trace_dc.build,
    "mp": trace_mp.build,
    "oc": trace_oc.build,
}
#: 默认对比集。`oc` = `oc-w1` = 原始真 OC 基线。
STRATEGIES = ("dc", "mp", "oc")

_OC_PREFIX = "oc-w"


def parse_strategy(name: str) -> tuple[str, int]:
    """把 `oc-w3` 解析成 (builder_key, tile_width)。

    `oc` 等价于 `oc-w1`。tile width 是**调度**参数，绝不进 hardware hash——
    OC-w1..w5 必须能互相比较，所以它们的硬件配置必须判定为相同。
    """
    if name.startswith(_OC_PREFIX):
        suffix = name[len(_OC_PREFIX):]
        if not suffix.isdigit() or int(suffix) < 1:
            raise KeyError(f"非法策略名 {name!r}；应形如 oc-w1 .. oc-w5")
        return "oc", int(suffix)
    if name not in STRATEGY_BUILDERS:
        raise KeyError(f"未知策略 {name!r}；可选 {STRATEGIES} 或 oc-w1..oc-w<bconv_cols>")
    return name, 1


def oc_tile_strategies(bconv_cols: int) -> tuple[str, ...]:
    """给定阵列列数，列出全部合法的 OC tile 宽度策略名。"""
    return tuple(f"{_OC_PREFIX}{w}" for w in range(1, bconv_cols + 1))


class FairnessError(RuntimeError):
    """三种策略没有跑在同一硬件配置上——对比无效。"""


@dataclass
class RunResult:
    strategy: str
    cfg: SimConfig
    wl: Workload
    res: ResourceSet
    ops: List[Op]
    sched: Optional[ScheduleResult]
    mem: MemoryStats
    hardware_hash: str
    workload_hash: str
    invocations: int
    tile_width: int = 1
    strategy_hash: str = ""
    infeasible_reason: str = ""
    metrics: Dict[str, object] = field(default_factory=dict)

    @property
    def feasible(self) -> bool:
        return not self.infeasible_reason


def run_strategy(cfg: SimConfig, strategy: str) -> RunResult:
    builder_key, tile_width = parse_strategy(strategy)
    cfg.validate()
    wl = build_workload(cfg.workload)
    cm = CostModel(wl=cfg.workload, hw=cfg.hardware)
    res = build_resources(cfg.hardware)

    ctx = TraceContext(strategy=strategy, wl=wl, cm=cm, res=res,
                       boundary=cfg.boundary, invocation=cfg.invocation,
                       tile_width=tile_width)
    ops = STRATEGY_BUILDERS[builder_key](ctx)

    if cfg.include_init:
        _prepend_init(ops, cm.init_cycles())

    hw_hash = hardware_config_hash(cfg)
    wl_hash = workload_config_hash(cfg.workload)
    st_hash = strategy_config_hash(strategy, tile_width)

    try:
        ops, mem = apply_capacity(ops, cfg.hardware.sram_capacity_bytes,
                                  cm.spill_cycles, strategy)
    except CapacityError as exc:
        from .memory import analyze_lifetimes
        return RunResult(strategy=strategy, cfg=cfg, wl=wl, res=res, ops=ops,
                         sched=None, mem=analyze_lifetimes(ops),
                         hardware_hash=hw_hash, workload_hash=wl_hash,
                         invocations=ctx.invocations, tile_width=tile_width,
                         strategy_hash=st_hash,
                         infeasible_reason=str(exc))

    sched = schedule(ops, res, tie_break=cfg.tie_break,
                     dma_compute_overlap=cfg.hardware.dma_compute_overlap,
                     allow_engine_overlap=cfg.hardware.allow_engine_overlap)

    out = RunResult(strategy=strategy, cfg=cfg, wl=wl, res=res, ops=ops,
                    sched=sched, mem=mem, hardware_hash=hw_hash,
                    workload_hash=wl_hash, invocations=ctx.invocations,
                    tile_width=tile_width, strategy_hash=st_hash)
    out.metrics = _metrics(out)
    return out


def _prepend_init(ops: List[Op], init_cycles: int) -> None:
    """把 OP_INIT 作为所有事件的前置串进去（上下文级固定段）。"""
    from .resources import CONTROL
    init = Op(seq=-1, kind=OpKind.INIT, resource=CONTROL, cycles=init_cycles,
              note="OP_INIT (context-level fixed segment)")
    for op in ops:
        op.seq += 1
        op.deps = tuple(d + 1 for d in op.deps)
    init.seq = 0
    for op in ops:
        if not op.deps:
            op.deps = (0,)
    ops.insert(0, init)


def _metrics(r: RunResult) -> Dict[str, object]:
    s, hw = r.sched, r.cfg.hardware
    assert s is not None
    h2d = sum(op.bytes_moved for op in r.ops if op.kind is OpKind.H2D)
    d2h = sum(op.bytes_moved for op in r.ops if op.kind is OpKind.D2H)
    km_res = r.res.keymult_resource

    m: Dict[str, object] = {
        "strategy": r.strategy,
        "hardware_config_hash": r.hardware_hash,
        "workload_config_hash": r.workload_hash,
        "strategy_config_hash": r.strategy_hash,
        "oc_output_tile_width": r.tile_width,
        "workload": r.wl.describe(),
        "boundary": r.cfg.boundary,
        "invocation_granularity": r.cfg.invocation,
        "invocations": r.invocations,
        "evidence": hw.cost.evidence_level(),
        "total_cycles": s.makespan,
        "latency_us": round(hw.cycles_to_us(s.makespan), 6),
        "compute_cycles": s.compute_cycles,
        "memory_stall_cycles": s.transfer_cycles,
        "control_cycles": s.control_cycles,
        "dependency_stall_cycles": s.idle_cycles,
        "resource_stall_cycles": s.resource_stall_cycles,
        "h2d_bytes": h2d,
        "d2h_bytes": d2h,
        "ddr_read_bytes": r.mem.spill_read_bytes,
        "ddr_write_bytes": r.mem.spill_write_bytes,
        "spill_read_bytes": r.mem.spill_read_bytes,
        "spill_write_bytes": r.mem.spill_write_bytes,
        "peak_live_bytes": r.mem.peak_live_bytes,
        "peak_resident_bytes": r.mem.peak_resident_bytes,
        "transform_utilization": round(s.utilization(TRANSFORM, r.res), 6),
        "bconv_utilization": round(s.utilization(BCONV, r.res), 6),
        "keymult_utilization": round(s.utilization(km_res, r.res), 6),
    }
    for kind, cnt in sorted(s.counts_by_kind.items()):
        m[f"n_{kind}"] = cnt
    for kind, cyc in sorted(s.cycles_by_kind.items()):
        m[f"cycles_{kind}"] = cyc
    return m


def run_all(cfg: SimConfig,
            strategies: Sequence[str] = STRATEGIES) -> List[RunResult]:
    results = [run_strategy(cfg, s) for s in strategies]
    check_fairness(results)
    return results


def check_fairness(results: Sequence[RunResult]) -> None:
    """方案 §3：每份结果输出 hardware_config_hash，证明用的是同一配置。"""
    if not results:
        return
    hashes = {r.strategy: r.hardware_hash for r in results}
    if len(set(hashes.values())) != 1:
        raise FairnessError(f"hardware_config_hash 不一致，对比无效：{hashes}")
    wl_hashes = {r.strategy: r.workload_hash for r in results}
    if len(set(wl_hashes.values())) != 1:
        raise FairnessError(f"workload_config_hash 不一致，对比无效：{wl_hashes}")


def arithmetic_identity(results: Sequence[RunResult]) -> Dict[str, Dict[str, int]]:
    """返回三种策略的算术事件计数，供 trace count test 断言。

    DC/MP/OC 的 INTT / SCALE / NTT / KeyMult / Accumulate 计数必须逐项相等；
    BConv 的**调用次数**允许不同（OC 是 single-tower，调用更多），但它产出的
    输出塔总数必须相等。把这两件事分开，才不会把「调用次数多」误读成「算得多」。
    """
    out: Dict[str, Dict[str, int]] = {}
    for r in results:
        counts = {k.value: 0 for k in (OpKind.INTT, OpKind.SCALE, OpKind.NTT,
                                       OpKind.BCONV, OpKind.KEYMULT,
                                       OpKind.ACCUMULATE)}
        produced = 0
        for op in r.ops:
            if op.kind.value in counts:
                counts[op.kind.value] += 1
            if op.kind is OpKind.BCONV:
                produced += op.work.get("beta", 0)
        counts["bconv_output_towers"] = produced
        out[r.strategy] = counts
    return out
