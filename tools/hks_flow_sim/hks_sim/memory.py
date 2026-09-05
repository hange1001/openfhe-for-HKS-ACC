"""片上存储模型：lifetime、peak 与 spill。

方案 §3 要求：live data 超过片上容量时必须插入 SpillStore/SpillLoad，
不能为某种策略隐式增加存储。这里的做法是在调度之前做一次容量 pre-pass，
把 spill 事件真实插进 trace，因此它们会出现在 Gantt 和 stall 统计里，
而不是变成一个事后加上去的修正项。

peak 分两个口径，不要混用：
  peak_live_bytes      -- 调度**要求**的片上驻留量，与容量无关（M1 报告用）
  peak_resident_bytes  -- spill 之后实际驻留量，受 capacity 约束（M4 报告用）
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from .op import Buffer, Op, OpKind, TraceBuilder
from .resources import DMA_SPILL


class CapacityError(RuntimeError):
    """容量不足且无法通过 spill 化解 —— 该配置点标记为 infeasible。"""


@dataclass
class MemoryStats:
    peak_live_bytes: int = 0
    peak_resident_bytes: int = 0
    spill_read_bytes: int = 0
    spill_write_bytes: int = 0
    spill_store_ops: int = 0
    spill_load_ops: int = 0
    peak_live_by_kind: Dict[str, int] = field(default_factory=dict)
    #: peak 出现时刻的 buffer 明细，便于回答「峰值是谁撑起来的」
    peak_live_breakdown: Dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, int]:
        return {
            "peak_live_bytes": self.peak_live_bytes,
            "peak_resident_bytes": self.peak_resident_bytes,
            "spill_read_bytes": self.spill_read_bytes,
            "spill_write_bytes": self.spill_write_bytes,
            "spill_store_ops": self.spill_store_ops,
            "spill_load_ops": self.spill_load_ops,
        }


def analyze_lifetimes(ops: List[Op]) -> MemoryStats:
    """静态 pass：只算 lifetime 与 peak，不做容量约束，不插入任何事件。

    M1 的完成标准「不依赖周期参数即可生成 operation-count / memory 报告」
    就靠这个函数。
    """
    stats = MemoryStats()
    live: Dict[str, Buffer] = {}
    live_bytes = 0
    for op in ops:
        for buf in op.allocs:
            if buf.bid in live:
                raise ValueError(f"buffer {buf.bid} 被重复分配（op {op.seq}）")
            live[buf.bid] = buf
            live_bytes += buf.nbytes
            if live_bytes > stats.peak_live_bytes:
                stats.peak_live_bytes = live_bytes
                by_kind: Dict[str, int] = {}
                for b in live.values():
                    by_kind[b.kind.value] = by_kind.get(b.kind.value, 0) + b.nbytes
                stats.peak_live_breakdown = by_kind
        for bid in op.frees:
            buf = live.pop(bid, None)
            if buf is None:
                raise ValueError(f"释放了未分配的 buffer {bid}（op {op.seq}）")
            live_bytes -= buf.nbytes
    if live:
        # 末尾仍驻留的是输出/常驻对象，不算泄漏，但要计入 peak 口径说明。
        pass
    stats.peak_resident_bytes = stats.peak_live_bytes
    stats.peak_live_by_kind = dict(stats.peak_live_breakdown)
    return stats


def apply_capacity(
    ops: List[Op],
    capacity_bytes: int,
    spill_cycles_fn,
    strategy: str,
) -> Tuple[List[Op], MemoryStats]:
    """容量 pre-pass：按 LRU 逐出，把 SpillStore/SpillLoad 插进 trace。

    capacity_bytes <= 0 表示不限容量，直接返回原 trace 与静态 lifetime 统计。
    """
    static = analyze_lifetimes(ops)
    if capacity_bytes <= 0:
        return ops, static

    stats = MemoryStats(peak_live_bytes=static.peak_live_bytes,
                        peak_live_by_kind=dict(static.peak_live_by_kind),
                        peak_live_breakdown=dict(static.peak_live_breakdown))

    out = TraceBuilder(strategy)
    remap: Dict[int, int] = {}          # 旧 seq -> 新 seq
    live: Dict[str, Buffer] = {}
    resident: Set[str] = set()
    spilled: Set[str] = set()
    last_touch: Dict[str, int] = {}
    #: 每个 buffer 最后一次写它的新 seq，spill/load 靠它串依赖
    producer: Dict[str, int] = {}
    resident_bytes = 0
    clock = 0

    def evict_until(need: int, pinned: Set[str], deps: List[int]) -> None:
        nonlocal resident_bytes
        while resident_bytes + need > capacity_bytes:
            victims = [b for b in resident if b not in pinned]
            if not victims:
                raise CapacityError(
                    f"[{strategy}] 容量 {capacity_bytes} B 装不下当前工作集："
                    f"需要 {need} B，已驻留 {resident_bytes} B，且无可逐出对象"
                )
            victim = min(victims, key=lambda b: last_touch.get(b, 0))
            buf = live[victim]
            dep = [producer[victim]] if victim in producer else []
            seq = out.emit(
                OpKind.SPILL_STORE, DMA_SPILL, dep,
                tower=buf.tower, digit=buf.digit, reads=[victim],
                bytes_moved=buf.nbytes, note=f"evict {victim}",
            )
            out.ops[seq].cycles = spill_cycles_fn(buf.nbytes)
            producer[victim] = seq
            deps.append(seq)
            resident.discard(victim)
            spilled.add(victim)
            resident_bytes -= buf.nbytes
            stats.spill_write_bytes += buf.nbytes
            stats.spill_store_ops += 1

    for op in ops:
        clock += 1
        extra_deps: List[int] = []
        touched = set(op.reads) | set(op.writes)
        pinned = {b for b in touched if b in live}

        # 1) 被读/写但已被逐出的 buffer 要先取回
        for bid in sorted(touched & spilled):
            buf = live[bid]
            evict_until(buf.nbytes, pinned, extra_deps)
            dep = [producer[bid]] if bid in producer else []
            seq = out.emit(
                OpKind.SPILL_LOAD, DMA_SPILL, dep,
                tower=buf.tower, digit=buf.digit, writes=[bid],
                bytes_moved=buf.nbytes, note=f"reload {bid}",
            )
            out.ops[seq].cycles = spill_cycles_fn(buf.nbytes)
            producer[bid] = seq
            extra_deps.append(seq)
            spilled.discard(bid)
            resident.add(bid)
            resident_bytes += buf.nbytes
            stats.spill_read_bytes += buf.nbytes
            stats.spill_load_ops += 1

        # 2) 为新分配腾地方
        for buf in op.allocs:
            evict_until(buf.nbytes, pinned, extra_deps)
            live[buf.bid] = buf
            resident.add(buf.bid)
            resident_bytes += buf.nbytes
            stats.peak_resident_bytes = max(stats.peak_resident_bytes, resident_bytes)

        deps = [remap[d] for d in op.deps] + extra_deps
        seq = out.emit(
            op.kind, op.resource, deps,
            digit=op.digit, tower=op.tower,
            reads=op.reads, writes=op.writes,
            allocs=op.allocs, frees=op.frees,
            bytes_moved=op.bytes_moved, work=op.work, note=op.note,
        )
        out.ops[seq].cycles = op.cycles
        remap[op.seq] = seq

        for bid in touched | {b.bid for b in op.allocs}:
            last_touch[bid] = clock
        for bid in op.writes:
            producer[bid] = seq

        # 3) 释放
        for bid in op.frees:
            buf = live.pop(bid, None)
            if buf is None:
                raise ValueError(f"释放了未分配的 buffer {bid}（op {op.seq}）")
            if bid in resident:
                resident.discard(bid)
                resident_bytes -= buf.nbytes
            spilled.discard(bid)
            producer.pop(bid, None)

    return out.ops, stats
