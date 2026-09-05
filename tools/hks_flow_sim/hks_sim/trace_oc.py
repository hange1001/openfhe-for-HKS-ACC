"""Output-Centric trace，带 output tiling（方案 §6.3 + M4）。

    retain/reload digit inputs under the common SRAM limit

    for each output tile (w 个连续输出塔):
        allocate 2*w 个 cTilda0/cTilda1 accumulator
        for each digit:
            native Q towers 继续 bypass
            把该 tile 中需要 BConv 的 target towers 合并成**一次**调用
            NTT each
            KeyMult / Accumulate 仍逐 tower 执行
        write completed output towers
        release accumulators

## 为什么是 tiling 而不是把阵列变窄

固定 3x5 阵列的 BConv 成本是 `N * ceil(alpha/rows) * ceil(beta/cols) + 49`。
只要 `beta <= cols`，产出 1 列和产出 5 列同价 4145 拍。所以：

  w=1（原始真 OC）：每个 (目标塔, digit) 一次调用 -> L*(D-1)+K*D 次，每次满价
  w=cols          ：每个 (tile, digit) 一次调用   -> 退化到与 DC 相同的调用次数

同一套硬件，只改调度粒度，就能把 OC 的 BConv 冗余从「sizeP 倍」压回去。
代价是一次要为 w 个输出塔各留一对 accumulator，即 2*w 个塔的常驻存储。
这就是 M4 要画的 Pareto：tile width ↑ -> BConv 调用 ↓、latency ↓、accumulator 存储 ↑。

## 两点仍然成立

1. **input residency 是 OC 的前提。** 输出塔循环在最外层，所有 digit 的系数域塔
   必须全程驻留；装不下就会被 memory 模型逐出并产生 SpillLoad，
   「每个 target tower 的重新读取成本必须计入」是自动发生的。
2. **native 塔一律 bypass**，不产生 BConv/NTT 事件，与 w 无关。
"""

from __future__ import annotations

from typing import List, Sequence

from .op import Op
from .trace_common import (
    BOUNDARY_FULL, INVOKE_PER_DIGIT, INVOKE_PER_HKS, INVOKE_PER_OUTPUT_TOWER,
    TraceContext, free_digit_coeffs, load_digit_inputs, modup_front_half,
)


def output_tiles(total_towers: int, width: int) -> List[List[int]]:
    """把 0..total_towers-1 切成宽度 <= width 的连续 tile。末块可能更短。"""
    return [list(range(t, min(t + width, total_towers)))
            for t in range(0, total_towers, width)]


def build(ctx: TraceContext) -> List[Op]:
    wl = ctx.wl
    full = ctx.boundary == BOUNDARY_FULL
    w = ctx.tile_width

    head: List[int] = []
    if ctx.invocation == INVOKE_PER_HKS:
        head = [ctx.invoke(note=f"oc-w{w}: single fused invocation")]

    # ---- 前置：所有 digit 的输入装载与 INTT/SCALE，结果全程驻留 ----
    ready: List[int] = []
    for dg in wl.digits:
        d = dg.index
        gate = list(head)
        if ctx.invocation == INVOKE_PER_DIGIT:
            gate = [ctx.invoke(head, note=f"oc-w{w}: invocation for digit {d} input")]
        loads = load_digit_inputs(ctx, d, gate)
        ready.extend(modup_front_half(ctx, d, loads))

    # ---- 主循环：以 output tile 为外层 ----
    for tile in output_tiles(wl.total_towers, w):
        gate = list(head)
        if ctx.invocation == INVOKE_PER_OUTPUT_TOWER:
            gate = [ctx.invoke(head, note=f"oc-w{w}: invocation for tile {tile}")]

        if full:
            for t in tile:
                ctx.alloc_accumulators(t, gate)

        tile_done: List[int] = []
        for dg in wl.digits:
            d = dg.index
            # 该 tile 中需要 BConv 的目标塔；native 的走旁路
            targets = [t for t in tile if not wl.is_native(d, t)]
            ready_by_tower = {t: list(gate) for t in tile if wl.is_native(d, t)}
            if targets:
                # 合并成一次 BConv 调用，成本按实际 len(targets)
                bc = ctx.bconv(d, targets, ready + gate)
                for t in targets:
                    ready_by_tower[t] = [ctx.ntt(d, t, [bc])]

            for t in tile:
                src_ready = ready_by_tower[t]
                if not full:
                    ctx.d2h_tower(ctx.ext_buffer_id(d, t), src_ready,
                                  tower=t, digit=d, note="partsCtExt writeback")
                    continue
                evk = ctx.load_evk(d, t, gate)
                km = ctx.keymult(d, t, src_ready + evk)
                acc = ctx.accumulate(d, t, [km])
                tile_done.append(acc)
                stale = list(ctx.free_evk(d, t))
                if not wl.is_native(d, t):
                    stale.append(f"ext.d{d}.t{t}")
                ctx.release(stale, [acc])

        if full:
            # 该 tile 的输出塔已集齐全部 digit 的贡献，写出并释放 accumulator
            for t in tile:
                ctx.d2h_tower(f"acc0.t{t}", tile_done, tower=t, note="cTilda0")
                ctx.d2h_tower(f"acc1.t{t}", tile_done, tower=t, note="cTilda1")

    # ---- 收尾：释放常驻的 digit 输入与系数域塔 ----
    stale: List[str] = []
    for dg in wl.digits:
        d = dg.index
        stale.extend(free_digit_coeffs(ctx, d))
        stale.extend(f"in.d{d}.t{tt}" for tt in dg.native_towers)
    ctx.release(stale)

    return ctx.finish()
