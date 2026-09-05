"""三种 trace 共用的发射器。

DC / MP / OC 的区别只允许出现在**发射顺序和 buffer 生命周期**上，
不允许出现在算子成本或资源绑定上——后两者由 cost_model 和 resources 统一决定。
任何一种 dataflow 想「少做一次运算」都必须体现为 trace 里少一个 Op，
而不是同一个 Op 变便宜。

buffer 命名约定：
    in.d{d}.t{t}     digit d 的 native 塔（EVAL 域输入，同时是旁路输出的来源）
    coeff.d{d}.t{t}  INTT+SCALE 之后的系数域塔
    ext.d{d}.t{t}    partsCtExt 的第 t 塔（EVAL 域），KeyMult 的输入
    evka/evkb.d{d}.t{t}  evaluation key 的 A/B 塔
    acc0/acc1.t{t}   cTilda0 / cTilda1 的累加器
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence

from .cost_model import CostModel
from .op import Buffer, BufferKind, Op, OpKind, TraceBuilder
from .resources import BCONV, CONTROL, DMA_D2H, DMA_H2D, TRANSFORM, ResourceSet
from .workload import Workload

#: 仿真边界。modup_only 复现 P4 已实测的 OP_HKS_DIGIT；full_hks 是方案 §2 的范围。
BOUNDARY_MODUP = "modup_only"
BOUNDARY_FULL = "full_hks"
BOUNDARIES = (BOUNDARY_MODUP, BOUNDARY_FULL)

#: kernel 调用粒度。方案 §12 点名的「OC 每个 target tower 一次调用」风险，
#: 在这里是一个可配置、可证伪的假设，而不是写死的结论。
INVOKE_PER_HKS = "per_hks"
INVOKE_PER_DIGIT = "per_digit"
INVOKE_PER_OUTPUT_TOWER = "per_output_tower"
INVOKE_GRANULARITIES = (INVOKE_PER_HKS, INVOKE_PER_DIGIT, INVOKE_PER_OUTPUT_TOWER)


@dataclass
class TraceContext:
    """一次 trace 生成的全部上下文。"""

    strategy: str
    wl: Workload
    cm: CostModel
    res: ResourceSet
    boundary: str = BOUNDARY_MODUP
    invocation: str = INVOKE_PER_DIGIT
    #: OC 的 output tile 宽度（一次处理多少个输出塔）。1 = 原始真 OC。
    #: 这是**调度**参数不是硬件参数：同一套 3x5 阵列，只改调度粒度。
    tile_width: int = 1
    tb: TraceBuilder = field(init=False)
    #: buffer id -> 最后写它的 op seq，用来串真实数据依赖
    producer: Dict[str, int] = field(default_factory=dict)
    invocations: int = 0

    def __post_init__(self) -> None:
        if self.boundary not in BOUNDARIES:
            raise ValueError(f"未知 boundary: {self.boundary}")
        if self.invocation not in INVOKE_GRANULARITIES:
            raise ValueError(f"未知 invocation granularity: {self.invocation}")
        if not 1 <= self.tile_width <= self.res.bconv_cols:
            raise ValueError(
                f"oc_output_tile_width={self.tile_width} 超出 [1, bconv_cols="
                f"{self.res.bconv_cols}]：一次调用最多利用 bconv_cols 个输出列"
            )
        self.tb = TraceBuilder(self.strategy)

    # ---------- 基础工具 ----------
    @property
    def tower_bytes(self) -> int:
        return self.wl.bytes_per_tower

    def dep_of(self, *bids: str) -> List[int]:
        """把 buffer id 翻译成它们各自 producer 的 seq。"""
        out: List[int] = []
        for bid in bids:
            seq = self.producer.get(bid)
            if seq is not None and seq not in out:
                out.append(seq)
        return out

    def _emit(self, kind: OpKind, resource: str, cycles: int,
              deps: Sequence[int] = (), **kw) -> int:
        seq = self.tb.emit(kind, resource, deps, **kw)
        op = self.tb.ops[seq]
        op.cycles = cycles
        for bid in op.writes:
            self.producer[bid] = seq
        for buf in op.allocs:
            self.producer.setdefault(buf.bid, seq)
        # 释放过的 buffer 立刻退出 producer 表，后续 release/dep_of 才不会重复引用
        for bid in op.frees:
            self.producer.pop(bid, None)
        return seq

    def buf(self, bid: str, kind: BufferKind,
            tower: Optional[int] = None, digit: Optional[int] = None) -> Buffer:
        return Buffer(bid=bid, kind=kind, nbytes=self.tower_bytes,
                      tower=tower, digit=digit)

    # ---------- 事件发射 ----------
    def invoke(self, deps: Sequence[int] = (), note: str = "") -> int:
        """一次 kernel 调用的固定控制开销（calibrated 656 @P4）。

        **调用边界是串行的。** XRT 路径完全阻塞（bo.sync -> run.wait -> bo.sync，
        无 async、无双缓冲，见 task.yaml q8），所以下一次 kernel 调用不能在上一次
        返回之前开始。这里把 invoke 实现成对**此前全部事件**的 barrier。

        不加这道 barrier 的话，调度器会让下一个 digit 的 H2D 与上一个 digit 的
        D2H 重叠，P4 两个 digit 就会算成 89533 而不是实测的 91242。
        """
        self.invocations += 1
        barrier = list(deps) + [op.seq for op in self.tb.ops]
        return self._emit(OpKind.INVOKE, CONTROL, self.cm.invocation_cycles(),
                          sorted(set(barrier)),
                          note=note or f"invocation#{self.invocations}")

    def h2d_tower(self, bid: str, kind: BufferKind, deps: Sequence[int] = (),
                  tower: Optional[int] = None, digit: Optional[int] = None,
                  note: str = "") -> int:
        return self._emit(
            OpKind.H2D, DMA_H2D, self.cm.transfer_cycles(1), deps,
            digit=digit, tower=tower, writes=[bid],
            allocs=[self.buf(bid, kind, tower, digit)],
            bytes_moved=self.tower_bytes, work={"towers": 1}, note=note,
        )

    def d2h_tower(self, bid: str, deps: Sequence[int] = (),
                  tower: Optional[int] = None, digit: Optional[int] = None,
                  free: bool = True, note: str = "") -> int:
        return self._emit(
            OpKind.D2H, DMA_D2H, self.cm.transfer_cycles(1),
            list(deps) + self.dep_of(bid),
            digit=digit, tower=tower, reads=[bid],
            frees=[bid] if free else [],
            bytes_moved=self.tower_bytes, work={"towers": 1}, note=note,
        )

    def intt(self, digit: int, tower: int, deps: Sequence[int] = ()) -> int:
        src = f"in.d{digit}.t{tower}"
        dst = f"coeff.d{digit}.t{tower}"
        return self._emit(
            OpKind.INTT, TRANSFORM, self.cm.transform_cycles(),
            list(deps) + self.dep_of(src),
            digit=digit, tower=tower, reads=[src], writes=[dst],
            allocs=[self.buf(dst, BufferKind.COEFF_TOWER, tower, digit)],
            work={"N": self.wl.N},
        )

    def scale(self, digit: int, tower: int, deps: Sequence[int] = ()) -> int:
        """QHatInv 预乘。P4 之后原位读写，复用变换引擎的 4 路模乘。"""
        bid = f"coeff.d{digit}.t{tower}"
        return self._emit(
            OpKind.SCALE, TRANSFORM, self.cm.scale_cycles(),
            list(deps) + self.dep_of(bid),
            digit=digit, tower=tower, reads=[bid], writes=[bid],
            work={"N": self.wl.N},
        )

    def bconv(self, digit: int, out_towers: Sequence[int],
              deps: Sequence[int] = ()) -> int:
        """一次 BConv 调用，产出 out_towers 全部列。

        DC/MP 每 digit 调用一次（out_towers = 全部 complement 塔）；
        OC 每个 target tower 调用一次（out_towers 只有一个元素）。
        成本按 N*ceil(alpha/rows)*ceil(beta/cols)，因此在固定 3x5 阵列上
        产出 1 列与产出 5 列同价——OC 的 single-tower BConv 不省时间。
        """
        dg = self.wl.digits[digit]
        srcs = [f"coeff.d{digit}.t{t}" for t in dg.native_towers]
        dsts = [f"bconv.d{digit}.t{t}" for t in out_towers]
        cycles = self.cm.bconv_cycles(dg.alpha, len(out_towers))
        return self._emit(
            OpKind.BCONV, BCONV, cycles,
            list(deps) + self.dep_of(*srcs),
            digit=digit, tower=out_towers[0] if len(out_towers) == 1 else None,
            reads=srcs, writes=dsts,
            allocs=[self.buf(b, BufferKind.BCONV_OUT, t, digit)
                    for b, t in zip(dsts, out_towers)],
            work={"N": self.wl.N, "alpha": dg.alpha, "beta": len(out_towers)},
        )

    def ntt(self, digit: int, tower: int, deps: Sequence[int] = ()) -> int:
        src = f"bconv.d{digit}.t{tower}"
        dst = f"ext.d{digit}.t{tower}"
        return self._emit(
            OpKind.NTT, TRANSFORM, self.cm.transform_cycles(),
            list(deps) + self.dep_of(src),
            digit=digit, tower=tower, reads=[src], writes=[dst],
            allocs=[self.buf(dst, BufferKind.PARTS_CT_EXT, tower, digit)],
            frees=[src],
            work={"N": self.wl.N},
        )

    def bypass(self, digit: int, tower: int) -> str:
        """native 塔直接复用原 EVAL 系数，不做 BConv/NTT。

        返回该塔在 partsCtExt 中的 buffer id。不产生任何 Op——
        「不做运算」在模型里就是 trace 里少一个事件。
        """
        return f"in.d{digit}.t{tower}"

    def ext_buffer_id(self, digit: int, tower: int) -> str:
        if self.wl.is_native(digit, tower):
            return self.bypass(digit, tower)
        return f"ext.d{digit}.t{tower}"

    def keymult(self, digit: int, tower: int, deps: Sequence[int] = ()) -> int:
        """c_j 与 evk 的两路模乘。projected：硬件未实现。"""
        src = self.ext_buffer_id(digit, tower)
        ka, kb = f"evka.d{digit}.t{tower}", f"evkb.d{digit}.t{tower}"
        dst = f"km.d{digit}.t{tower}"
        return self._emit(
            OpKind.KEYMULT, self.res.keymult_resource, self.cm.keymult_cycles(),
            list(deps) + self.dep_of(src, ka, kb),
            digit=digit, tower=tower, reads=[src, ka, kb], writes=[dst],
            allocs=[self.buf(dst, BufferKind.PARTS_CT_EXT, tower, digit)],
            work={"N": self.wl.N, "passes": self.cm.c.keymult_mul_passes},
        )

    def accumulate(self, digit: int, tower: int, deps: Sequence[int] = (),
                   free_km: bool = True) -> int:
        """把 KeyMult 结果累加进 cTilda0/cTilda1。projected。"""
        km = f"km.d{digit}.t{tower}"
        a0, a1 = f"acc0.t{tower}", f"acc1.t{tower}"
        return self._emit(
            OpKind.ACCUMULATE, self.res.keymult_resource, self.cm.accumulate_cycles(),
            list(deps) + self.dep_of(km, a0, a1),
            digit=digit, tower=tower, reads=[km, a0, a1], writes=[a0, a1],
            frees=[km] if free_km else [],
            work={"N": self.wl.N},
        )

    def alloc_accumulators(self, tower: int, deps: Sequence[int] = ()) -> int:
        """给一个输出塔开 cTilda0/cTilda1 两个累加器（零初始化，不走 AXI）。"""
        a0, a1 = f"acc0.t{tower}", f"acc1.t{tower}"
        return self._emit(
            OpKind.RELEASE, CONTROL, 0, deps,
            tower=tower, writes=[a0, a1],
            allocs=[self.buf(a0, BufferKind.ACCUMULATOR, tower),
                    self.buf(a1, BufferKind.ACCUMULATOR, tower)],
            note=f"zero-init accumulators t{tower}",
        )

    def release(self, bids: Sequence[str], deps: Sequence[int] = ()) -> Optional[int]:
        """释放若干 buffer。零成本事件，只为让 lifetime 精确到事件级。

        已经不在场的 buffer（例如被 NTT 顺手 free 掉的 bconv 中间塔）自动跳过。
        """
        alive = [b for b in dict.fromkeys(bids) if b in self.producer]
        if not alive:
            return None
        return self._emit(OpKind.RELEASE, CONTROL, 0,
                          list(deps) + self.dep_of(*alive),
                          frees=alive, note="release")

    def load_evk(self, digit: int, tower: int, deps: Sequence[int] = ()) -> List[int]:
        out = []
        for tag in ("evka", "evkb"):
            bid = f"{tag}.d{digit}.t{tower}"
            out.append(self.h2d_tower(bid, BufferKind.EVK_TOWER, deps,
                                      tower=tower, digit=digit, note=tag))
        return out

    def free_evk(self, digit: int, tower: int) -> List[str]:
        return [f"evka.d{digit}.t{tower}", f"evkb.d{digit}.t{tower}"]

    # ---------- 结果 ----------
    def finish(self) -> List[Op]:
        return self.tb.ops


def load_digit_inputs(ctx: TraceContext, digit: int,
                      deps: Sequence[int] = ()) -> List[int]:
    """H2D 装入 digit 的 native 塔。它同时是 INTT 输入和旁路输出的来源。"""
    dg = ctx.wl.digits[digit]
    return [ctx.h2d_tower(f"in.d{digit}.t{t}", BufferKind.DIGIT_INPUT, deps,
                          tower=t, digit=digit, note="digit input")
            for t in dg.native_towers]


def modup_front_half(ctx: TraceContext, digit: int,
                     deps: Sequence[int] = ()) -> List[int]:
    """INTT -> SCALE，处理 digit 的全部 native 塔。返回各塔完成的 seq。"""
    dg = ctx.wl.digits[digit]
    out = []
    for t in dg.native_towers:
        s = ctx.intt(digit, t, deps)
        out.append(ctx.scale(digit, t, [s]))
    return out


def free_digit_coeffs(ctx: TraceContext, digit: int) -> List[str]:
    dg = ctx.wl.digits[digit]
    return [f"coeff.d{digit}.t{t}" for t in dg.native_towers]
