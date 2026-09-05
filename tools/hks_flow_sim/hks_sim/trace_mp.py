"""Maximum-Parallel trace（方案 §6.2）。

    all digits INTT/SCALE
    barrier
    all digits BConv
    barrier
    all digits NTT
    barrier
    all digits KeyMult
    Accumulate

MP 与 DC 的算术量必须一致；差异来自 barrier、并行机会和 live working set。
单 transform engine + 单 BConv array 下，MP 不得通过虚假并行缩短时间——
这由 scheduler 的资源容量保证，trace 这边只负责如实写出 barrier 依赖。

barrier 的实现：下一阶段的每个事件依赖上一阶段的**全部**事件。这会把
所有中间塔的生命周期拉长到跨阶段，正是 MP 峰值存储偏高的来源。
"""

from __future__ import annotations

from typing import List

from .op import Op
from .trace_common import (
    BOUNDARY_FULL, INVOKE_PER_HKS,
    TraceContext, free_digit_coeffs, load_digit_inputs, modup_front_half,
)


def build(ctx: TraceContext) -> List[Op]:
    wl = ctx.wl
    full = ctx.boundary == BOUNDARY_FULL

    head: List[int] = []
    if ctx.invocation == INVOKE_PER_HKS:
        head = [ctx.invoke(note="mp: single fused invocation")]

    if full:
        for t in range(wl.total_towers):
            ctx.alloc_accumulators(t, head)

    # ---- 阶段 1：全部 digit 的 H2D + INTT + SCALE ----
    phase1: List[int] = []
    for dg in wl.digits:
        d = dg.index
        gate = list(head)
        if ctx.invocation != INVOKE_PER_HKS:
            gate = [ctx.invoke(head, note=f"mp: invocation for digit {d}")]
        loads = load_digit_inputs(ctx, d, gate)
        phase1.extend(modup_front_half(ctx, d, loads))

    # ---- barrier ----
    # ---- 阶段 2：全部 digit 的 BConv ----
    phase2: List[int] = []
    for dg in wl.digits:
        phase2.append(ctx.bconv(dg.index, dg.complement_towers, phase1))
    for dg in wl.digits:
        ctx.release(free_digit_coeffs(ctx, dg.index), phase2)

    # ---- barrier ----
    # ---- 阶段 3：全部 digit 的 NTT ----
    phase3: List[int] = []
    for dg in wl.digits:
        for t in dg.complement_towers:
            phase3.append(ctx.ntt(dg.index, t, phase2))

    if not full:
        for dg in wl.digits:
            d = dg.index
            for t in range(wl.total_towers):
                ctx.d2h_tower(ctx.ext_buffer_id(d, t), phase3, tower=t, digit=d,
                              note="partsCtExt writeback")
        return ctx.finish()

    # ---- barrier ----
    # ---- 阶段 4：全部 KeyMult，随后统一 Accumulate ----
    kms: List[int] = []
    for dg in wl.digits:
        d = dg.index
        for t in range(wl.total_towers):
            evk = ctx.load_evk(d, t, head)
            kms.append(ctx.keymult(d, t, phase3 + evk))

    accs: List[int] = []
    for dg in wl.digits:
        d = dg.index
        for t in range(wl.total_towers):
            accs.append(ctx.accumulate(d, t, kms))

    stale: List[str] = []
    for dg in wl.digits:
        d = dg.index
        stale.extend(f"in.d{d}.t{t}" for t in dg.native_towers)
        for t in range(wl.total_towers):
            stale.extend(ctx.free_evk(d, t))
            if not wl.is_native(d, t):
                stale.append(f"ext.d{d}.t{t}")
    ctx.release(stale, accs)

    for t in range(wl.total_towers):
        ctx.d2h_tower(f"acc0.t{t}", tower=t, note="cTilda0")
        ctx.d2h_tower(f"acc1.t{t}", tower=t, note="cTilda1")

    return ctx.finish()
