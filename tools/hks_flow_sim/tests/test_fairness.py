"""Fairness test（方案 §12）：三种策略必须跑在同一硬件配置上。

以及一条方案没写、但会静默毁掉结论的：调度器的 tie-break 与引擎并发假设
必须计入 hash，否则换一个假设就能换一个赢家而 hash 毫无变化。
"""

import dataclasses
import unittest

from hks_sim.config import (SimConfig, hardware_config_hash,
                            workload_config_hash)
from hks_sim.engine import FairnessError, check_fairness, run_all, run_strategy


class TestHashEquality(unittest.TestCase):
    def test_all_strategies_share_hardware_hash(self):
        results = run_all(SimConfig())
        self.assertEqual(len({r.hardware_hash for r in results}), 1)
        self.assertEqual(len({r.workload_hash for r in results}), 1)

    def test_check_fairness_rejects_mismatch(self):
        cfg_a = SimConfig()
        cfg_b = dataclasses.replace(
            cfg_a, hardware=dataclasses.replace(cfg_a.hardware, transform_lanes=8))
        results = [run_strategy(cfg_a, "dc"), run_strategy(cfg_b, "oc")]
        with self.assertRaises(FairnessError):
            check_fairness(results)


class TestHashSensitivity(unittest.TestCase):
    """凡是能改变赢家的自由度，都必须改变 hash。"""

    def setUp(self):
        self.base = SimConfig()

    def _hash(self, **kw):
        return hardware_config_hash(dataclasses.replace(self.base, **kw))

    def _hw_hash(self, **kw):
        hw = dataclasses.replace(self.base.hardware, **kw)
        return hardware_config_hash(dataclasses.replace(self.base, hardware=hw))

    def test_tie_break_changes_hash(self):
        self.assertNotEqual(self._hash(), self._hash(tie_break="fifo"))

    def test_boundary_and_invocation_change_hash(self):
        self.assertNotEqual(self._hash(), self._hash(boundary="full_hks"))
        self.assertNotEqual(self._hash(), self._hash(invocation="per_hks"))

    def test_engine_overlap_changes_hash(self):
        """OC 的胜负就悬在这一条上，它绝不能悄悄改变而 hash 不动。"""
        self.assertNotEqual(self._hash(), self._hw_hash(allow_engine_overlap=True))

    def test_dma_overlap_changes_hash(self):
        self.assertNotEqual(self._hash(), self._hw_hash(dma_compute_overlap=True))

    def test_cost_params_change_hash(self):
        cost = dataclasses.replace(self.base.hardware.cost, keymult_pipe=19)
        self.assertNotEqual(self._hash(), self._hw_hash(cost=cost))

    def test_hash_is_stable_across_runs(self):
        self.assertEqual(self._hash(), self._hash())

    def test_workload_hash_tracks_shape(self):
        wl = dataclasses.replace(self.base.workload, ring_dimension=8192)
        self.assertNotEqual(workload_config_hash(self.base.workload),
                            workload_config_hash(wl))


class TestNoStrategySpecificHardware(unittest.TestCase):
    def test_resource_capacities_identical_across_strategies(self):
        results = run_all(SimConfig())
        caps = [r.res.capacity for r in results]
        for c in caps[1:]:
            self.assertEqual(c, caps[0])

    def test_op_costs_identical_across_strategies(self):
        """同一类算子在三种策略里必须同价——省时间只能靠少发事件。"""
        results = run_all(SimConfig())
        by_kind = {}
        for r in results:
            for op in r.ops:
                if op.cycles <= 0 or op.kind.value == "BConv":
                    continue          # BConv 依 beta 变价，另有断言
                by_kind.setdefault(op.kind.value, {}).setdefault(r.strategy, set()
                                                                ).add(op.cycles)
        for kind, per_strategy in by_kind.items():
            sets = list(per_strategy.values())
            for s in sets[1:]:
                self.assertEqual(s, sets[0], f"{kind} 在不同策略下成本不同")


if __name__ == "__main__":
    unittest.main()
