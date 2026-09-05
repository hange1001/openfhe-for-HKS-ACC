"""Capacity test（方案 §12）：容量不足时必须 spill 或判 infeasible。

以及 §3 的硬性要求：不能为某种策略隐式增加存储。
"""

import dataclasses
import unittest

from hks_sim.config import SimConfig
from hks_sim.engine import run_all, run_strategy
from hks_sim.memory import analyze_lifetimes
from hks_sim.op import OpKind


def _cfg(**hw_kw):
    base = SimConfig()
    return dataclasses.replace(base, hardware=dataclasses.replace(base.hardware,
                                                                  **hw_kw))


class TestLifetimeAccounting(unittest.TestCase):
    def test_peak_is_multiple_of_tower_size(self):
        r = run_strategy(SimConfig(), "dc")
        tower = r.wl.bytes_per_tower
        self.assertEqual(tower, 8 * 4096)
        self.assertEqual(r.mem.peak_live_bytes % tower, 0)

    def test_every_allocation_is_eventually_freed_or_is_output(self):
        """allocs/frees 必须配平；未释放的只允许是最终输出。"""
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(_cfg(), strategy)
            live = {}
            for op in r.ops:
                for b in op.allocs:
                    live[b.bid] = b
                for bid in op.frees:
                    live.pop(bid, None)
            for bid in live:
                self.assertTrue(bid.startswith(("acc0.", "acc1.", "ext.", "in.")),
                                f"{strategy}: 未释放的中间 buffer {bid}")

    def test_double_free_is_rejected(self):
        r = run_strategy(SimConfig(), "dc")
        ops = list(r.ops)
        ops.append(dataclasses.replace(ops[-1], seq=len(ops), allocs=(),
                                       frees=("in.d0.t0",)))
        with self.assertRaises(ValueError):
            analyze_lifetimes(ops)


class TestMpCostsMorePeakThanDc(unittest.TestCase):
    def test_barriers_inflate_working_set(self):
        """MP 的 barrier 把生命周期拉长，峰值必然不低于 DC。

        方案 §12 的待验证判断之一：MP 可能不减少 compute cycles，反而增加
        peak memory。
        """
        results = {r.strategy: r for r in run_all(SimConfig())}
        self.assertGreater(results["mp"].mem.peak_live_bytes,
                           results["dc"].mem.peak_live_bytes)


class TestCapacityPressure(unittest.TestCase):
    def test_tight_capacity_produces_spill_events(self):
        tower = 8 * 4096
        r = run_strategy(_cfg(sram_capacity_bytes=4 * tower), "dc")
        if r.feasible:
            kinds = {op.kind for op in r.ops}
            self.assertTrue(
                {OpKind.SPILL_STORE, OpKind.SPILL_LOAD} & kinds,
                "容量收紧后既没 spill 也没判 infeasible")
            self.assertGreater(r.mem.spill_write_bytes, 0)

    def test_impossible_capacity_is_reported_not_silently_ignored(self):
        """连一个塔都放不下时必须显式 infeasible，不能悄悄跑出个数来。"""
        r = run_strategy(_cfg(sram_capacity_bytes=1024), "dc")
        self.assertFalse(r.feasible)
        self.assertIn("容量", r.infeasible_reason)

    def test_spill_never_reduces_total_cycles(self):
        """Monotonicity（方案 §12）：减小 SRAM 不应让它变快。"""
        tower = 8 * 4096
        big = run_strategy(_cfg(sram_capacity_bytes=0), "dc")
        small = run_strategy(_cfg(sram_capacity_bytes=5 * tower), "dc")
        if small.feasible:
            self.assertGreaterEqual(small.metrics["total_cycles"],
                                    big.metrics["total_cycles"])

    def test_unbounded_capacity_inserts_no_spill(self):
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(_cfg(sram_capacity_bytes=0), strategy)
            self.assertEqual(r.mem.spill_read_bytes, 0)
            self.assertEqual(r.mem.spill_write_bytes, 0)


if __name__ == "__main__":
    unittest.main()
