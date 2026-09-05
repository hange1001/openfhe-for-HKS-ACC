"""M3 标定测试：P4 零残差 + P3 hold-out。

这是整个仿真器最重要的一组断言。方案 M3 的完成标准是「P1-P3 路径预测与可信
RTL 数据的目标偏差不超过 5%」；这里实际做到 P4 零残差、P3 hold-out 0.009%。
"""

import unittest

from hks_sim.config import SimConfig
from hks_sim.cost_model import CostModel
from hks_sim.engine import run_strategy

# --- measured，全部来自 docs/reports/hls/hks_mem_p4_20260905/ ---
P4_TRANSFORM = 6481      # transform-work-csynth.rpt: WORK_STAGE_LOOP 6*1080 + 1
P4_SCALE = 1047          # 同上，latency min（SCALE 模式）
P4_BCONV = 4145          # top-csynth.rpt: Compute_BConv_Systolic
P4_DIGIT_A2 = 46671      # perf-transactions.rpt transaction 1
P4_DIGIT_A1 = 44571      # perf-transactions.rpt transaction 2
P4_TWO_DIGITS = 91242
P4_INIT = 295063         # perf-transactions.rpt transaction 0

# --- hold-out：P3，docs/reports/hls/hks_mem_p3_20260904/perf-transactions.rpt ---
P3_DIGIT_A2 = 52779
P3_DIGIT_A1 = 47643
P3_TWO_DIGITS = 100422


class TestP4Calibration(unittest.TestCase):
    def setUp(self):
        cfg = SimConfig()
        self.cm = CostModel(wl=cfg.workload, hw=cfg.hardware)

    def test_primitive_cycles_match_hls_reports(self):
        self.assertEqual(self.cm.transform_cycles(), P4_TRANSFORM)
        self.assertEqual(self.cm.scale_cycles(), P4_SCALE)
        self.assertEqual(self.cm.bconv_cycles(2, 3), P4_BCONV)
        self.assertEqual(self.cm.bconv_cycles(1, 4), P4_BCONV)
        self.assertEqual(self.cm.init_cycles(), P4_INIT)

    def test_axi_beats_follow_bus_width(self):
        # N=4096 个 64-bit 系数走 256-bit 总线 = 1024 拍
        self.assertEqual(self.cm.axi_beats_per_tower(), 1024)
        self.assertEqual(self.cm.axi_tower_cycles(), 1024 + 29)

    def _digit_total(self, alpha, beta):
        """按模型算一个 digit 的事务周期。5 = alpha + beta 次变换，恒定。"""
        return (5 * self.cm.transform_cycles()
                + alpha * self.cm.scale_cycles()
                + self.cm.bconv_cycles(alpha, beta)
                + (alpha + 5) * self.cm.axi_tower_cycles()   # load alpha + store L+K
                + self.cm.invocation_cycles())

    def test_zero_residual_on_both_p4_transactions(self):
        self.assertEqual(self._digit_total(2, 3), P4_DIGIT_A2)
        self.assertEqual(self._digit_total(1, 4), P4_DIGIT_A1)

    def test_single_tower_bconv_is_not_cheaper_on_fixed_array(self):
        """固定 3x5 阵列上，产出 1 列与产出 5 列同价。

        这不是模型偷懒，而是 ADR-010 / OC_strategy_gap_analysis 描述的
        「sizeP 倍冗余」的周期形式。OC 想省下 BConv 时间必须缩窄阵列或做列分块。
        """
        self.assertEqual(self.cm.bconv_cycles(2, 1), self.cm.bconv_cycles(2, 5))


class TestP3HoldOut(unittest.TestCase):
    """把 P4 标定出的 axi/invocation 参数原样搬到结构不同的 P3。

    P3 的预乘是独立硬件、1 系数/拍（P4 才改成复用 4 路模乘），所以 P3 的
    per-tower 预乘成本是唯一的未知量。用两个 P3 事务各解一次，两次结果
    必须互相吻合——这才叫 hold-out，而不是再拟合一次。
    """

    def setUp(self):
        cfg = SimConfig()
        self.cm = CostModel(wl=cfg.workload, hw=cfg.hardware)

    def _solve_p3_prescale(self, alpha, beta, measured):
        fixed = (5 * self.cm.transform_cycles()
                 + self.cm.bconv_cycles(alpha, beta)
                 + (alpha + 5) * self.cm.axi_tower_cycles()
                 + self.cm.invocation_cycles())
        return (measured - fixed) / alpha

    def test_two_transactions_agree_on_prescale_cost(self):
        s2 = self._solve_p3_prescale(2, 3, P3_DIGIT_A2)
        s1 = self._solve_p3_prescale(1, 4, P3_DIGIT_A1)
        self.assertAlmostEqual(s2, 4101.0, places=6)
        self.assertAlmostEqual(s1, 4119.0, places=6)
        # 两个独立方程对同一未知量的解，相差必须远小于 M3 的 5% 门槛
        self.assertLess(abs(s2 - s1) / max(s2, s1), 0.005)

    def test_holdout_prediction_error_under_one_percent(self):
        s = (self._solve_p3_prescale(2, 3, P3_DIGIT_A2)
             + self._solve_p3_prescale(1, 4, P3_DIGIT_A1)) / 2.0

        def predict(alpha, beta):
            return (5 * self.cm.transform_cycles() + alpha * s
                    + self.cm.bconv_cycles(alpha, beta)
                    + (alpha + 5) * self.cm.axi_tower_cycles()
                    + self.cm.invocation_cycles())

        total = predict(2, 3) + predict(1, 4)
        err = abs(total - P3_TWO_DIGITS) / P3_TWO_DIGITS
        self.assertLess(err, 0.001, f"P3 hold-out 误差 {err:.5%} 超标")

    def test_p3_prescale_is_one_coefficient_per_cycle(self):
        """解出来的 P3 预乘应当接近 N 拍（1 系数/拍），而 P4 是 N/4 拍。

        这是对模型结构的独立印证：P4 报告说预乘从独立硬件改成复用 4 路 lane，
        比值就应该是 lane 数 4。
        """
        s = (self._solve_p3_prescale(2, 3, P3_DIGIT_A2)
             + self._solve_p3_prescale(1, 4, P3_DIGIT_A1)) / 2.0
        self.assertAlmostEqual(s / self.cm.scale_cycles(), 4.0, delta=0.1)


class TestTraceReproducesMeasurement(unittest.TestCase):
    """端到端：DC 的 modup_only trace 调度后必须精确等于实测两 digit 周期。"""

    def test_dc_modup_only_equals_measured(self):
        r = run_strategy(SimConfig(), "dc")
        self.assertEqual(r.metrics["total_cycles"], P4_TWO_DIGITS)
        self.assertAlmostEqual(r.metrics["latency_us"], 638.694, places=3)

    def test_no_idle_cycles_in_blocking_model(self):
        """XRT 完全阻塞 + 引擎顺序调用下，不应出现纯依赖空转。"""
        r = run_strategy(SimConfig(), "dc")
        self.assertEqual(r.metrics["dependency_stall_cycles"], 0)


if __name__ == "__main__":
    unittest.main()
