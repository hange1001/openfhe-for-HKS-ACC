"""锁死一句容易说错的话，并锁死它的**理由**。

正确说法：DC 与 MP 的算子数量和理论 arithmetic workload 相同；在当前建模下
两者的周期指标也相同；**但 MP 需要约 2.8 倍 peak memory**，因为 barrier 把
中间塔的生命周期拉长到跨阶段。

不能说成「两者等价」——peak 差 2.8 倍是实打实的代价。
也不能说成「两者各项性能指标必然相同」——那只在当前的资源模型下成立，
一旦 AXI 通道数或 barrier 语义变化就可能分开（见 test_delta_appears_with_more_axi_channels）。

历史注记：早先按用途把 AXI 拆成 dma_h2d/dma_d2h/dma_spill 三个各自容量 1 的池，
等于白送 3 倍并发带宽，当时 DC 能靠增量回写藏掉 8424 拍而 MP 不能，于是两者
差 8424。那个差是**建模假象**：P4 实测 46671 只有在 AXI 全串行时才对得上。
改为单一共享池后差值归零。
"""

import dataclasses
import unittest
from collections import defaultdict

from hks_sim.config import SimConfig
from hks_sim.engine import arithmetic_identity, run_strategy
from hks_sim.op import COMPUTE_KINDS, TRANSFER_KINDS

MP_PEAK_RATIO_MIN = 2.5     # MP / DC 的 peak 倍数下限


def _full_cfg(**hw):
    base = SimConfig()
    cost = dataclasses.replace(base.hardware.cost, keymult_pipe=38,
                               keymult_mem=256, accumulate_overhead=26)
    return dataclasses.replace(
        base,
        hardware=dataclasses.replace(base.hardware, cost=cost, **hw),
        boundary="full_hks", invocation="per_hks")


def _per_kind(r):
    cyc, cnt = defaultdict(int), defaultdict(int)
    for op in r.ops:
        if op.cycles > 0:
            cyc[op.kind.value] += op.cycles
            cnt[op.kind.value] += 1
    return dict(cyc), dict(cnt)


class TestArithmeticIsIdentical(unittest.TestCase):
    def setUp(self):
        cfg = _full_cfg()
        self.dc = run_strategy(cfg, "dc")
        self.mp = run_strategy(cfg, "mp")

    def test_operator_counts_are_identical(self):
        ident = arithmetic_identity([self.dc, self.mp])
        self.assertEqual(ident["dc"], ident["mp"])

    def test_per_kind_cycle_sums_are_identical(self):
        self.assertEqual(_per_kind(self.dc), _per_kind(self.mp))

    def test_transfer_bytes_are_identical(self):
        for key in ("h2d_bytes", "d2h_bytes"):
            self.assertEqual(self.dc.metrics[key], self.mp.metrics[key])


class TestPeakMemoryIsWhereTheyDiffer(unittest.TestCase):
    """这才是 MP 的真实代价，也是「不能说两者等价」的理由。"""

    def setUp(self):
        cfg = _full_cfg()
        self.dc = run_strategy(cfg, "dc")
        self.mp = run_strategy(cfg, "mp")

    def test_mp_peak_is_much_larger(self):
        ratio = self.mp.mem.peak_live_bytes / self.dc.mem.peak_live_bytes
        self.assertGreaterEqual(ratio, MP_PEAK_RATIO_MIN,
                                f"MP/DC peak 倍数只有 {ratio:.2f}")

    def test_peak_gap_survives_in_modup_only_too(self):
        cfg = SimConfig()
        self.assertGreater(run_strategy(cfg, "mp").mem.peak_live_bytes,
                           run_strategy(cfg, "dc").mem.peak_live_bytes)

    def test_mp_needs_more_capacity_to_avoid_spill(self):
        """把容量卡在 DC 的 peak 上：DC 不 spill，MP 必须 spill。"""
        dc_peak = self.dc.mem.peak_live_bytes
        cfg = _full_cfg(sram_capacity_bytes=dc_peak)
        dc = run_strategy(cfg, "dc")
        mp = run_strategy(cfg, "mp")
        self.assertEqual(dc.mem.spill_write_bytes, 0)
        self.assertGreater(mp.mem.spill_write_bytes, 0)
        self.assertGreater(mp.metrics["total_cycles"], dc.metrics["total_cycles"])


class TestCycleParityIsConditional(unittest.TestCase):
    """周期相等**不是**必然的，它依赖当前的 AXI 资源模型。"""

    def test_cycles_are_equal_under_serial_axi(self):
        cfg = _full_cfg()
        self.assertEqual(run_strategy(cfg, "dc").metrics["total_cycles"],
                         run_strategy(cfg, "mp").metrics["total_cycles"])

    def test_delta_appears_with_more_axi_channels(self):
        """放宽到多通道，DC 的增量回写就能藏掉传输，两者立刻分开。

        这条断言存在的意义：提醒「DC=MP」只在串行 AXI 下成立，
        换了资源假设必须重新报数，不能把它当成 dataflow 的性质。
        """
        cfg = _full_cfg(axi_channels=3)
        dc = run_strategy(cfg, "dc").metrics["total_cycles"]
        mp = run_strategy(cfg, "mp").metrics["total_cycles"]
        self.assertLess(dc, mp)


class TestSchedulingInvariants(unittest.TestCase):
    def setUp(self):
        cfg = _full_cfg()
        self.dc = run_strategy(cfg, "dc")
        self.mp = run_strategy(cfg, "mp")

    def test_no_transfer_ever_overlaps_compute(self):
        for r in (self.dc, self.mp):
            xfer = [(op.start_cycle, op.end_cycle) for op in r.ops
                    if op.kind in TRANSFER_KINDS and op.cycles > 0]
            comp = [(op.start_cycle, op.end_cycle) for op in r.ops
                    if op.kind in COMPUTE_KINDS and op.cycles > 0]
            for xs, xe in xfer:
                for cs, ce in comp:
                    self.assertFalse(xs < ce and cs < xe,
                                     f"{r.strategy}: 传输与计算重叠")

    def test_no_two_transfers_overlap_on_single_channel(self):
        """axi_channels=1：任意两笔搬运不得同时在飞。"""
        for r in (self.dc, self.mp):
            spans = sorted((op.start_cycle, op.end_cycle) for op in r.ops
                           if op.kind in TRANSFER_KINDS and op.cycles > 0)
            for (s1, e1), (s2, e2) in zip(spans, spans[1:]):
                self.assertLessEqual(e1, s2, f"{r.strategy}: 两笔搬运重叠")


class TestModupOnlyStillMatchesMeasurement(unittest.TestCase):
    def test_calibration_survives_the_shared_axi_pool(self):
        cfg = SimConfig()
        for s in ("dc", "mp"):
            self.assertEqual(run_strategy(cfg, s).metrics["total_cycles"], 91242)


if __name__ == "__main__":
    unittest.main()
