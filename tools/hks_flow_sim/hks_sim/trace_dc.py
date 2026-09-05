"""Digit-Centric trace（方案 §6.1）。

    for each digit d:
        INTT/SCALE native towers
        BConv all complement towers      <- 一次调用产出全部 beta_d 列
        NTT required towers
        KeyMult all output towers
        Accumulate into cTilda0/cTilda1
        release digit temporaries

DC 可以释放当前 digit 的 ModUp 临时数据，但 output accumulator 必须跨 digit 保活
（方案 §6.1）。这就是 DC 的存储特征：临时数据浅，累加器全程常驻。
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
        head = [ctx.invoke(note="dc: single fused invocation")]

    # 全局累加器：DC 逐 digit 累加，必须在所有 digit 之前开、全部结束后才释放
    if full:
        for t in range(wl.total_towers):
            ctx.alloc_accumulators(t, head)

    for dg in wl.digits:
        d = dg.index
        gate = list(head)
        if ctx.invocation != INVOKE_PER_HKS:
            gate = [ctx.invoke(head, note=f"dc: invocation for digit {d}")]

        loads = load_digit_inputs(ctx, d, gate)
        scaled = modup_front_half(ctx, d, loads)

        # 一次 BConv 产出该 digit 的全部 complement 塔
        bc = ctx.bconv(d, dg.complement_towers, scaled)
        ctx.release(free_digit_coeffs(ctx, d), [bc])

        ntts = [ctx.ntt(d, t, [bc]) for t in dg.complement_towers]

        if not full:
            # modup_only：partsCtExt 直接写回，复现 P4 的 OP_HKS_DIGIT 事务
            for t in range(wl.total_towers):
                ctx.d2h_tower(ctx.ext_buffer_id(d, t), ntts, tower=t, digit=d,
                              note="partsCtExt writeback")
            continue

        for t in range(wl.total_towers):
            evk = ctx.load_evk(d, t, gate)
            km = ctx.keymult(d, t, ntts + evk)
            acc = ctx.accumulate(d, t, [km])
            # 该 digit 对该塔的贡献已并入累加器，临时数据立即释放
            stale = list(ctx.free_evk(d, t))
            if not wl.is_native(d, t):
                stale.append(f"ext.d{d}.t{t}")
            ctx.release(stale, [acc])

        ctx.release([f"in.d{d}.t{t}" for t in dg.native_towers])

    if full:
        for t in range(wl.total_towers):
            ctx.d2h_tower(f"acc0.t{t}", tower=t, note="cTilda0")
            ctx.d2h_tower(f"acc1.t{t}", tower=t, note="cTilda1")

    return ctx.finish()
