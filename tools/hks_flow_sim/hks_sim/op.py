"""硬件事件与缓冲区的数据结构。

trace 生成器只负责产生 Op 的 DAG 和 buffer 生命周期；周期数由 cost_model 填，
makespan 由 scheduler 算。三者分开，是为了让 M1（不依赖周期参数的算子/存储报告）
能独立成立。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Sequence, Tuple


class OpKind(str, Enum):
    INIT = "INIT"
    INVOKE = "Invoke"     # 一次 kernel 调用的固定控制开销（descriptor/参数/dispatch）
    RELEASE = "Release"   # 零成本记账事件，只承载 buffer 的 alloc/free，不占时间
    H2D = "H2D"
    D2H = "D2H"
    INTT = "INTT"
    SCALE = "SCALE"
    BCONV = "BConv"
    NTT = "NTT"
    KEYMULT = "KeyMult"
    ACCUMULATE = "Accumulate"
    SPILL_STORE = "SpillStore"
    SPILL_LOAD = "SpillLoad"


#: 占用计算引擎的事件。makespan 分解时算作 compute。
COMPUTE_KINDS = frozenset({
    OpKind.INTT, OpKind.SCALE, OpKind.BCONV,
    OpKind.NTT, OpKind.KEYMULT, OpKind.ACCUMULATE,
})

#: 搬运事件。makespan 分解时算作 transfer。
TRANSFER_KINDS = frozenset({
    OpKind.H2D, OpKind.D2H, OpKind.SPILL_STORE, OpKind.SPILL_LOAD,
})

#: 控制事件。既不是算术也不是搬运，单独一档，避免混进 compute 抬高利用率。
CONTROL_KINDS = frozenset({OpKind.INIT, OpKind.INVOKE})

#: 零成本记账事件。不占用任何资源时间，只用来精确标记 buffer 生命周期。
BOOKKEEPING_KINDS = frozenset({OpKind.RELEASE})


class BufferKind(str, Enum):
    """§7 要求跟踪的存储对象类别。"""

    DIGIT_INPUT = "digit_input"           # digit native 塔（EVAL 域输入 / 旁路来源）
    COEFF_TOWER = "coeff_tower"           # INTT/SCALE 之后的系数域塔
    BCONV_OUT = "bconv_out"               # BConv 产出、待 NTT 的塔
    PARTS_CT_EXT = "parts_ct_ext"         # ModUp 输出（EVAL 域），KeyMult 的输入
    EVK_TOWER = "evk_tower"               # evaluation key A/B
    ACCUMULATOR = "accumulator"           # cTilda0/cTilda1
    TWIDDLE = "twiddle"                   # twiddle / 常数 / metadata
    DMA_BUFFER = "dma_buffer"             # H2D/D2H 暂存


@dataclass(frozen=True)
class Buffer:
    """一块片上缓冲。生命周期由 Op 的 allocs/frees 驱动。"""

    bid: str
    kind: BufferKind
    nbytes: int
    tower: Optional[int] = None
    digit: Optional[int] = None
    note: str = ""


@dataclass
class Op:
    """一个硬件事件。

    seq 同时是 id 和 trace 发射序。tie_break="trace_order" 时调度器用它打破平局，
    因此 seq 必须严格按 trace 语义顺序分配，不能事后重排。
    """

    seq: int
    kind: OpKind
    resource: str
    deps: Tuple[int, ...] = ()
    cycles: int = 0                              # 由 cost_model 填
    digit: Optional[int] = None
    tower: Optional[int] = None
    reads: Tuple[str, ...] = ()
    writes: Tuple[str, ...] = ()
    allocs: Tuple[Buffer, ...] = ()
    frees: Tuple[str, ...] = ()
    bytes_moved: int = 0                         # 仅 TRANSFER_KINDS 有意义
    work: Dict[str, int] = field(default_factory=dict)   # N/alpha/beta/tower 数等
    note: str = ""

    # -- 调度结果，由 scheduler 回填 --
    dep_ready_cycle: int = 0
    start_cycle: int = 0
    end_cycle: int = 0
    stall_reason: str = ""

    @property
    def is_compute(self) -> bool:
        return self.kind in COMPUTE_KINDS

    @property
    def is_transfer(self) -> bool:
        return self.kind in TRANSFER_KINDS

    def label(self) -> str:
        bits = [self.kind.value]
        if self.digit is not None:
            bits.append(f"d{self.digit}")
        if self.tower is not None:
            bits.append(f"t{self.tower}")
        return ".".join(bits)


class TraceBuilder:
    """按发射序累积 Op，保证 seq 单调、依赖只指向更早的事件。"""

    def __init__(self, strategy: str) -> None:
        self.strategy = strategy
        self.ops: List[Op] = []

    def emit(
        self,
        kind: OpKind,
        resource: str,
        deps: Sequence[int] = (),
        *,
        digit: Optional[int] = None,
        tower: Optional[int] = None,
        reads: Sequence[str] = (),
        writes: Sequence[str] = (),
        allocs: Sequence[Buffer] = (),
        frees: Sequence[str] = (),
        bytes_moved: int = 0,
        work: Optional[Dict[str, int]] = None,
        note: str = "",
    ) -> int:
        seq = len(self.ops)
        for d in deps:
            if not 0 <= d < seq:
                raise ValueError(f"op {seq} 依赖了非法或未来的事件 {d}")
        self.ops.append(Op(
            seq=seq, kind=kind, resource=resource, deps=tuple(deps),
            digit=digit, tower=tower,
            reads=tuple(reads), writes=tuple(writes),
            allocs=tuple(allocs), frees=tuple(frees),
            bytes_moved=bytes_moved, work=dict(work or {}), note=note,
        ))
        return seq

    def counts_by_kind(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for op in self.ops:
            out[op.kind.value] = out.get(op.kind.value, 0) + 1
        return out
