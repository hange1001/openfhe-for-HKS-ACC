"""M0 完成标准：正确推导 (alpha_d, beta_d, D)，拒绝非法参数。"""

import unittest

from hks_sim.config import ConfigError, WorkloadConfig
from hks_sim.workload import build_workload


class TestDigitLayout(unittest.TestCase):
    def test_p4_shape_matches_documented_reality(self):
        """L=3,K=2,D=2 必须给出 alpha=2/beta=3 与 alpha=1/beta=4。

        这两个形状在 docs/reports/hls/hks_mem_p4_20260905/ 里有独立事务周期，
        推错了后面全错。
        """
        wl = build_workload(WorkloadConfig())
        self.assertEqual(wl.D, 2)
        self.assertEqual(wl.alpha, 2)
        self.assertEqual([(d.alpha, d.beta) for d in wl.digits], [(2, 3), (1, 4)])

    def test_incomplete_last_digit_formula(self):
        """alpha_d = min(alpha, L - d*alpha)，beta_d = L + K - alpha_d。"""
        for L, K, D in [(3, 2, 2), (12, 4, 4), (7, 3, 3), (5, 2, 2), (24, 8, 8)]:
            wl = build_workload(WorkloadConfig(q_towers=L, p_towers=K, num_part_q=D))
            alpha = wl.alpha
            for d, dg in enumerate(wl.digits):
                self.assertEqual(dg.alpha, min(alpha, L - d * alpha),
                                 f"L={L} K={K} D={D} d={d}")
                self.assertEqual(dg.beta, L + K - dg.alpha)

    def test_digits_partition_q_towers_exactly(self):
        wl = build_workload(WorkloadConfig(q_towers=7, p_towers=3, num_part_q=3))
        covered = sorted(t for dg in wl.digits for t in dg.native_towers)
        self.assertEqual(covered, list(range(7)))

    def test_native_and_complement_are_disjoint_and_complete(self):
        wl = build_workload(WorkloadConfig(q_towers=12, p_towers=4, num_part_q=4))
        for dg in wl.digits:
            self.assertEqual(set(dg.native_towers) & set(dg.complement_towers), set())
            self.assertEqual(sorted(dg.native_towers + dg.complement_towers),
                             list(range(wl.total_towers)))


class TestOperationCounts(unittest.TestCase):
    def test_intt_count_equals_q_towers(self):
        wl = build_workload(WorkloadConfig(q_towers=12, p_towers=4, num_part_q=4))
        self.assertEqual(wl.operation_counts()["intt"], 12)

    def test_bconv_closed_form_matches_expansion(self):
        """方案 §6.3 的 N_BConv_OC = L*(D-1) + K*D 必须等于 sum_d beta_d。

        这条恒等式说明 DC 与 OC 的 BConv **算术量**相同，差别只在调用粒度。
        """
        for L, K, D in [(3, 2, 2), (12, 4, 4), (24, 8, 8), (6, 3, 3)]:
            wl = build_workload(WorkloadConfig(q_towers=L, p_towers=K, num_part_q=D))
            self.assertEqual(wl.operation_counts()["bconv"], wl.oc_bconv_closed_form())
            self.assertEqual(wl.oc_bconv_closed_form(), L * (D - 1) + K * D)

    def test_p4_counts_match_documented_attribution(self):
        """项目自己的周期归因写的是「3 次 INTT + 7 次 NTT」。"""
        counts = build_workload(WorkloadConfig()).operation_counts()
        self.assertEqual(counts["intt"], 3)
        self.assertEqual(counts["ntt"], 7)
        self.assertEqual(counts["bconv"], 7)


class TestRejectsIllegalConfigs(unittest.TestCase):
    def test_non_power_of_two_ring_dimension(self):
        with self.assertRaises(ConfigError):
            build_workload(WorkloadConfig(ring_dimension=5000))

    def test_num_part_q_exceeding_q_towers(self):
        with self.assertRaises(ConfigError):
            build_workload(WorkloadConfig(q_towers=3, num_part_q=4))

    def test_alpha_inconsistent_with_num_part_q(self):
        # L=12, alpha=3 -> 4 个 digit；却声明 num_part_q=3，自相矛盾
        with self.assertRaises(ConfigError):
            build_workload(WorkloadConfig(q_towers=12, num_part_q=3, alpha=3))

    def test_zero_p_towers(self):
        with self.assertRaises(ConfigError):
            build_workload(WorkloadConfig(p_towers=0))

    def test_limb_bits_out_of_range(self):
        with self.assertRaises(ConfigError):
            build_workload(WorkloadConfig(limb_bits=128))


if __name__ == "__main__":
    unittest.main()
