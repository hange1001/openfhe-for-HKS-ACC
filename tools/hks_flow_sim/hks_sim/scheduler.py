"""资源约束下的事件驱动调度器（M2）。

不用 CPU wall-clock 模拟 FPGA：按事件 DAG + 资源容量 + 存储条件算 makespan。

**tie-break 是公平性的一部分。** 方案 §5.2 只规定「依赖完成且资源可用即可启动」，
但当多个 ready 事件争同一个引擎时选谁并没有规定，而这个自由度会系统性偏袒某种
dataflow（LIFO 偏袒 DC 的深度优先，FIFO 偏袒 MP 的广度优先）。默认
tie_break="trace_order"：严格按 trace 发射序，因为发射序本身就是 dataflow 的
定义，这是唯一不引入额外自由度的选择。该字段计入 hardware_config_hash。

makespan 的四分解（互斥且求和等于 makespan）：
    compute  -- 至少一个计算事件在跑
    transfer -- 没有计算在跑，但有搬运在跑
    control  -- 只有控制事件在跑
    idle     -- 什么都没跑（纯依赖等待）
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .op import (BOOKKEEPING_KINDS, COMPUTE_KINDS, CONTROL_KINDS,
                 TRANSFER_KINDS, Op, OpKind)
from .resources import ResourceSet


class ScheduleError(RuntimeError):
    """DAG 有环，或存在永远无法满足的依赖。"""


@dataclass
class ScheduleResult:
    ops: List[Op]
    makespan: int = 0
    compute_cycles: int = 0
    transfer_cycles: int = 0
    control_cycles: int = 0
    idle_cycles: int = 0
    resource_busy: Dict[str, int] = field(default_factory=dict)
    resource_stall_cycles: int = 0
    cycles_by_kind: Dict[str, int] = field(default_factory=dict)
    counts_by_kind: Dict[str, int] = field(default_factory=dict)

    def utilization(self, name: str, res: ResourceSet) -> float:
        if self.makespan <= 0:
            return 0.0
        return self.resource_busy.get(name, 0) / (self.makespan * res.cap(name))

    def as_dict(self) -> Dict[str, float]:
        return {
            "total_cycles": self.makespan,
            "compute_cycles": self.compute_cycles,
            "memory_stall_cycles": self.transfer_cycles,
            "control_cycles": self.control_cycles,
            "dependency_stall_cycles": self.idle_cycles,
            "resource_stall_cycles": self.resource_stall_cycles,
        }


def _class_of(op: Op) -> Optional[str]:
    if op.cycles <= 0 or op.kind in BOOKKEEPING_KINDS:
        return None
    if op.kind in COMPUTE_KINDS:
        return "compute"
    if op.kind in TRANSFER_KINDS:
        return "transfer"
    if op.kind in CONTROL_KINDS:
        return "control"
    return None


def schedule(
    ops: List[Op],
    res: ResourceSet,
    *,
    tie_break: str = "trace_order",
    dma_compute_overlap: bool = False,
    allow_engine_overlap: bool = False,
) -> ScheduleResult:
    n = len(ops)
    remaining = [len(op.deps) for op in ops]
    dependents: List[List[int]] = [[] for _ in range(n)]
    for op in ops:
        for d in op.deps:
            dependents[d].append(op.seq)

    ready: List[int] = [i for i in range(n) if remaining[i] == 0]
    ready_at: Dict[int, int] = {i: 0 for i in ready}
    free_cap = dict(res.capacity)
    running: List[Tuple[int, int]] = []      # (end_cycle, seq)
    active_compute = 0
    active_transfer = 0
    done = 0
    now = 0
    resource_stall = 0

    def sort_key(i: int):
        if tie_break == "fifo":
            return (ready_at.get(i, 0), i)
        if tie_break == "lifo":
            return (-ready_at.get(i, 0), -i)
        return (i,)                          # trace_order

    def can_start(op: Op) -> bool:
        if free_cap.get(op.resource, 0) <= 0:
            return False
        # opcode-RPC dispatcher（ADR-008）：顶层顺序调用各计算引擎，
        # 不同引擎之间也不并发，除非显式假设已做 DATAFLOW 拆分。
        if (not allow_engine_overlap
                and op.kind in COMPUTE_KINDS and active_compute > 0):
            return False
        if dma_compute_overlap:
            return True
        # 不允许 DMA 与计算重叠：两类互斥
        if op.kind in COMPUTE_KINDS and active_transfer > 0:
            return False
        if op.kind in TRANSFER_KINDS and active_compute > 0:
            return False
        return True

    def retire(seq: int) -> None:
        nonlocal done
        done += 1
        for nxt in dependents[seq]:
            remaining[nxt] -= 1
            if remaining[nxt] == 0:
                ready.append(nxt)
                ready_at[nxt] = now

    while done < n:
        progressed = True
        while progressed:
            progressed = False
            # 零成本记账事件立刻完成，不占资源、不推进时间
            for i in sorted([i for i in ready if ops[i].cycles <= 0], key=sort_key):
                op = ops[i]
                op.dep_ready_cycle = ready_at.get(i, now)
                op.start_cycle = op.end_cycle = now
                ready.remove(i)
                retire(i)
                progressed = True

            for i in sorted(list(ready), key=sort_key):
                op = ops[i]
                if not can_start(op):
                    continue
                op.dep_ready_cycle = ready_at.get(i, now)
                op.start_cycle = now
                op.end_cycle = now + op.cycles
                op.stall_reason = ("resource" if op.start_cycle > op.dep_ready_cycle
                                   else "")
                resource_stall += op.start_cycle - op.dep_ready_cycle
                free_cap[op.resource] -= 1
                if op.kind in COMPUTE_KINDS:
                    active_compute += 1
                elif op.kind in TRANSFER_KINDS:
                    active_transfer += 1
                heapq.heappush(running, (op.end_cycle, i))
                ready.remove(i)
                progressed = True

        if done >= n:
            break
        if not running:
            stuck = [i for i in range(n) if remaining[i] > 0 and i not in ready]
            raise ScheduleError(
                f"调度死锁：已完成 {done}/{n}，ready={len(ready)}，"
                f"仍有 {len(stuck)} 个事件依赖未满足（DAG 可能有环）"
            )

        now = running[0][0]
        while running and running[0][0] == now:
            _, seq = heapq.heappop(running)
            op = ops[seq]
            free_cap[op.resource] += 1
            if op.kind in COMPUTE_KINDS:
                active_compute -= 1
            elif op.kind in TRANSFER_KINDS:
                active_transfer -= 1
            retire(seq)

    return _summarize(ops, res, now, resource_stall)


def _summarize(ops: List[Op], res: ResourceSet, makespan: int,
               resource_stall: int) -> ScheduleResult:
    out = ScheduleResult(ops=ops, makespan=makespan,
                         resource_stall_cycles=resource_stall)

    for op in ops:
        out.counts_by_kind[op.kind.value] = out.counts_by_kind.get(op.kind.value, 0) + 1
        if op.cycles > 0:
            out.cycles_by_kind[op.kind.value] = (
                out.cycles_by_kind.get(op.kind.value, 0) + op.cycles)
            out.resource_busy[op.resource] = (
                out.resource_busy.get(op.resource, 0) + op.cycles)

    # 区间并集扫描：compute > transfer > control，同一拍只归一类
    events: List[Tuple[int, int, str]] = []
    for op in ops:
        cls = _class_of(op)
        if cls is not None and op.end_cycle > op.start_cycle:
            events.append((op.start_cycle, op.end_cycle, cls))
    if events:
        bounds = sorted({b for s, e, _ in events for b in (s, e)})
        for lo, hi in zip(bounds, bounds[1:]):
            span = hi - lo
            classes = {c for s, e, c in events if s <= lo and e >= hi}
            if "compute" in classes:
                out.compute_cycles += span
            elif "transfer" in classes:
                out.transfer_cycles += span
            elif "control" in classes:
                out.control_cycles += span
    out.idle_cycles = (makespan - out.compute_cycles
                       - out.transfer_cycles - out.control_cycles)
    return out
