"""Output-Centric trace（方案 §6.3）。

    retain/reload digit inputs under the common SRAM limit

    for each target output tower:
        allocate one cTilda0/cTilda1 tower accumulator
        for each required digit:
            bypass native Q tower when valid
            otherwise BConv only the target tower
            NTT if required
            KeyMult immediately
            Accumulate immediately
        write completed output tower
        release accumulator

两点必须看清楚，否则会把 OC 的收益记错：

1. **input residency 是 OC 的前提。** 输出塔循环在最外层，意味着所有 digit 的
   系数域塔（INTT+SCALE 之后）必须全程驻留。装不下就会被 memory 模型逐出并
   产生 SpillLoad——方案 §6.3 说的「每个 target tower 的重新读取成本必须计入」
   在这里是自动发生的，不需要额外建模。

2. **single-tower BConv 在固定阵列上不省时间。** cost_model 的
   ceil(beta/cols) 决定了：3x5 阵列产出 1 列和产出 5 列同价。因此 OC 的
   BConv 调用次数从 D 次涨到 L*(D-1)+K*D 次，每次仍是满价。这正是
   ADR-010 / OC_strategy_gap_analysis 说的「sizeP 倍冗余」的周期形式。
   要让 OC 真正省下 BConv 时间，必须缩窄阵列或开列分块（allow_tiling）。
"""

from __future__ import annotations

from typing import List

from .op import Op
from .trace_common import (
    BOUNDARY_FULL, INVOKE_PER_DIGIT, INVOKE_PER_HKS, INVOKE_PER_OUTPUT_TOWER,
    TraceContext, free_digit_coeffs, load_digit_inputs, modup_front_half,
)


def build(ctx: TraceContext) -> List[Op]:
    wl = ctx.wl
    full = ctx.boundary == BOUNDARY_FULL

    head: List[int] = []
    if ctx.invocation == INVOKE_PER_HKS:
        head = [ctx.invoke(note="oc: single fused invocation")]

    # ---- 前置：所有 digit 的输入装载与 INTT/SCALE，结果全程驻留 ----
    ready: List[int] = []
    for dg in wl.digits:
        d = dg.index
        gate = list(head)
        if ctx.invocation == INVOKE_PER_DIGIT:
            gate = [ctx.invoke(head, note=f"oc: invocation for digit {d} input")]
        loads = load_digit_inputs(ctx, d, gate)
        ready.extend(modup_front_half(ctx, d, loads))

    # ---- 主循环：以输出塔为外层 ----
    for t in range(wl.total_towers):
        gate = list(head)
        if ctx.invocation == INVOKE_PER_OUTPUT_TOWER:
            gate = [ctx.invoke(head, note=f"oc: invocation for output tower {t}")]

        if full:
            ctx.alloc_accumulators(t, gate)

        done_this_tower: List[int] = []
        for dg in wl.digits:
            d = dg.index
            if wl.is_native(d, t):
                # 旁路：直接复用原 EVAL 系数，不产生 BConv/NTT 事件
                src_ready: List[int] = list(gate)
            else:
                bc = ctx.bconv(d, [t], ready + gate)      # 只算这一个目标塔
                nt = ctx.ntt(d, t, [bc])
                src_ready = [nt]

            if not full:
                ctx.d2h_tower(ctx.ext_buffer_id(d, t), src_ready,
                              tower=t, digit=d, note="partsCtExt writeback")
                continue

            evk = ctx.load_evk(d, t, gate)
            km = ctx.keymult(d, t, src_ready + evk)
            acc = ctx.accumulate(d, t, [km])
            done_this_tower.append(acc)
            stale = list(ctx.free_evk(d, t))
            if not wl.is_native(d, t):
                stale.append(f"ext.d{d}.t{t}")
            ctx.release(stale, [acc])

        if full:
            # 该输出塔已集齐所有 digit 的贡献，立刻写出并释放累加器
            w0 = ctx.d2h_tower(f"acc0.t{t}", done_this_tower, tower=t, note="cTilda0")
            w1 = ctx.d2h_tower(f"acc1.t{t}", done_this_tower, tower=t, note="cTilda1")
            del w0, w1

    # ---- 收尾：释放常驻的 digit 输入与系数域塔 ----
    stale: List[str] = []
    for dg in wl.digits:
        d = dg.index
        stale.extend(free_digit_coeffs(ctx, d))
        stale.extend(f"in.d{d}.t{tt}" for tt in dg.native_towers)
    ctx.release(stale)

    return ctx.finish()
