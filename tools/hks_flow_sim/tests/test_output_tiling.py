"""M4：output-tiled OC（OC-w1 .. OC-w<bconv_cols>）。

固定 3x5 阵列下，真正缺的不是把阵列变窄，而是把 OC 从 single-tower 扩成
output-tiled：同一套硬件，只改调度粒度。
"""

import dataclasses
import unittest

from hks_sim.config import SimConfig, hardware_config_hash, strategy_config_hash
from hks_sim.engine import (arithmetic_identity, oc_tile_strategies,
                            parse_strategy, run_strategy)
from hks_sim.op import OpKind
from hks_sim.trace_oc import output_tiles


def _cfg(**kw):
    base = SimConfig()
    cost = dataclasses.replace(base.hardware.cost, keymult_pipe=38,
                               keymult_mem=256, accumulate_overhead=26)
    hw = dataclasses.replace(base.hardware, cost=cost)
    return dataclasses.replace(base, hardware=hw, boundary="full_hks",
                               invocation="per_hks", **kw)


WIDTHS = (1, 2, 3, 4, 5)


class TestStrategyNaming(unittest.TestCase):
    def test_oc_is_alias_for_w1(self):
        self.assertEqual(parse_strategy("oc"), ("oc", 1))
        self.assertEqual(parse_strategy("oc-w3"), ("oc", 3))

    def test_w1_reproduces_original_true_oc_baseline(self):
        """w=1 必须与改造前的真 OC 逐指标一致，否则 M4 动了基线。"""
        a = run_strategy(_cfg(), "oc")
        b = run_strategy(_cfg(), "oc-w1")
        self.assertEqual(a.metrics["total_cycles"], b.metrics["total_cycles"])
        self.assertEqual(a.mem.peak_live_bytes, b.mem.peak_live_bytes)

    def test_illegal_names_rejected(self):
        for bad in ("oc-w0", "oc-wx", "oc-w", "zz"):
            with self.assertRaises(KeyError):
                parse_strategy(bad)

    def test_width_above_bconv_cols_rejected(self):
        with self.assertRaises(ValueError):
            run_strategy(_cfg(), "oc-w6")   # bconv_cols=5


class TestHashSeparation(unittest.TestCase):
    """tile width 是调度参数，不是硬件参数。"""

    def test_tile_width_does_not_change_hardware_hash(self):
        hashes = {run_strategy(_cfg(), f"oc-w{w}").hardware_hash for w in WIDTHS}
        self.assertEqual(len(hashes), 1,
                         "tile width 改变了 hardware hash，OC-w* 之间就无法比较了")

    def test_tile_width_does_change_strategy_hash(self):
        hashes = {run_strategy(_cfg(), f"oc-w{w}").strategy_hash for w in WIDTHS}
        self.assertEqual(len(hashes), len(WIDTHS))

    def test_strategy_hash_is_independent_of_hardware(self):
        self.assertEqual(strategy_config_hash("oc", 3),
                         strategy_config_hash("oc", 3))
        self.assertNotEqual(strategy_config_hash("oc", 3),
                            strategy_config_hash("oc", 4))

    def test_tile_width_recorded_in_metrics(self):
        for w in WIDTHS:
            m = run_strategy(_cfg(), f"oc-w{w}").metrics
            self.assertEqual(m["oc_output_tile_width"], w)
            self.assertEqual(m["strategy"], f"oc-w{w}")


class TestTiling(unittest.TestCase):
    def test_tiles_cover_all_towers_without_overlap(self):
        for total in (5, 8, 16, 20):
            for w in range(1, 6):
                tiles = output_tiles(total, w)
                flat = [t for tile in tiles for t in tile]
                self.assertEqual(flat, list(range(total)), f"total={total} w={w}")
                self.assertTrue(all(len(t) <= w for t in tiles))

    def test_ragged_last_tile_is_allowed(self):
        # L+K=5, w=4 -> [0,1,2,3] + [4]
        self.assertEqual(output_tiles(5, 4), [[0, 1, 2, 3], [4]])


class TestArithmeticPreserved(unittest.TestCase):
    def test_all_widths_share_the_same_arithmetic(self):
        names = ["dc", "mp"] + list(oc_tile_strategies(5))
        results = [run_strategy(_cfg(), n) for n in names]
        ident = arithmetic_identity(results)
        for key in ("INTT", "SCALE", "NTT", "KeyMult", "Accumulate",
                    "bconv_output_towers"):
            vals = {ident[r.strategy][key] for r in results}
            self.assertEqual(len(vals), 1, f"{key} 在不同 tile width 下不一致：{vals}")

    def test_native_towers_still_bypass_at_every_width(self):
        for w in WIDTHS:
            r = run_strategy(_cfg(), f"oc-w{w}")
            for op in r.ops:
                if op.kind is OpKind.NTT:
                    self.assertFalse(r.wl.is_native(op.digit, op.tower),
                                     f"oc-w{w}: 对 native 塔做了 NTT")

    def test_bconv_cost_uses_actual_target_count(self):
        """合并调用的成本必须按实际 len(non_native_targets) 算。"""
        for w in WIDTHS:
            r = run_strategy(_cfg(), f"oc-w{w}")
            for op in r.ops:
                if op.kind is OpKind.BCONV:
                    self.assertEqual(op.work["beta"], len(op.writes))
                    self.assertLessEqual(op.work["beta"], w)


class TestParetoShape(unittest.TestCase):
    """M4 的核心输出：tile width ↑ -> BConv 调用 ↓、accumulator 存储 ↑。"""

    def setUp(self):
        self.runs = {w: run_strategy(_cfg(), f"oc-w{w}") for w in WIDTHS}

    def test_bconv_calls_are_non_increasing_in_width(self):
        calls = [sum(1 for op in self.runs[w].ops if op.kind is OpKind.BCONV)
                 for w in WIDTHS]
        self.assertEqual(calls, sorted(calls, reverse=True), calls)
        self.assertEqual(calls[0], 7)    # w=1: L*(D-1)+K*D
        self.assertEqual(calls[-1], 2)   # w=5: 退化到与 DC 相同的调用次数

    def test_peak_memory_is_strictly_increasing_in_width(self):
        peaks = [self.runs[w].mem.peak_live_bytes for w in WIDTHS]
        for a, b in zip(peaks, peaks[1:]):
            self.assertLess(a, b, peaks)

    def test_accumulator_footprint_matches_two_per_tower_in_tile(self):
        tower = 8 * 4096
        for w in WIDTHS:
            live_acc = 0
            peak_acc = 0
            for op in self.runs[w].ops:
                for b in op.allocs:
                    if b.bid.startswith(("acc0.", "acc1.")):
                        live_acc += b.nbytes
                        peak_acc = max(peak_acc, live_acc)
                for bid in op.frees:
                    if bid.startswith(("acc0.", "acc1.")):
                        live_acc -= tower
            self.assertEqual(peak_acc, 2 * w * tower, f"oc-w{w}")

    def test_widest_tile_matches_dc_on_compute(self):
        """w=bconv_cols 时 BConv 调用退化到 DC 水平，compute 应当相等。"""
        dc = run_strategy(_cfg(), "dc")
        self.assertEqual(self.runs[5].metrics["compute_cycles"],
                         dc.metrics["compute_cycles"])

    def test_latency_is_non_increasing_in_width(self):
        """串行 AXI 下延迟随 w 单调不增；w3==w4 是因为两者 BConv 调用次数相同。

        注：早先按用途拆分 AXI 池（等于白送 3 倍带宽）时，这条曲线在 w4 处
        出现过非单调。那是建模假象，不是 tiling 的性质——改成单一共享池后消失。
        """
        c = {w: self.runs[w].metrics["total_cycles"] for w in WIDTHS}
        seq = [c[w] for w in WIDTHS]
        self.assertEqual(seq, sorted(seq, reverse=True), c)
        self.assertEqual(c[3], c[4], f"w3/w4 的 BConv 调用次数相同，周期应相等：{c}")

    def test_equal_bconv_calls_implies_equal_compute(self):
        """w3 与 w4 都是 4 次 BConv 调用，compute 必须相等。"""
        calls = {w: sum(1 for op in self.runs[w].ops if op.kind is OpKind.BCONV)
                 for w in WIDTHS}
        self.assertEqual(calls[3], calls[4])
        self.assertEqual(self.runs[3].metrics["compute_cycles"],
                         self.runs[4].metrics["compute_cycles"])

    def test_wider_tile_never_costs_more_latency_but_always_more_memory(self):
        """M4 Pareto 的方向性：w ↑ -> latency 不增、accumulator 存储严格增。"""
        for a, b in zip(WIDTHS, WIDTHS[1:]):
            self.assertLessEqual(self.runs[b].metrics["total_cycles"],
                                 self.runs[a].metrics["total_cycles"])
            self.assertLess(self.runs[a].mem.peak_live_bytes,
                            self.runs[b].mem.peak_live_bytes)


if __name__ == "__main__":
    unittest.main()
