"""Trace count test（方案 §12）：检查算子数量和 bypass。

比方案要求的更强一条：三种 dataflow 的 INTT/SCALE/NTT/KeyMult/Accumulate
计数必须**逐项相等**，BConv 的输出塔总数也必须相等。允许不同的只有 BConv
的调用次数（OC 是 single-tower，调用更多但每次产出更少）。

理由：sum_d beta_d = D*(L+K) - L = L*(D-1) + K*D，DC 与 OC 恒等。所以任何
一方算子数不同，就是 trace 写错了，而不是 dataflow 的性质。这与推导v1 §2.3
「每算子一次 offload 下三策略 AI 恒等」是同一件事。
"""

import dataclasses
import unittest

from hks_sim.config import SimConfig, WorkloadConfig
from hks_sim.engine import arithmetic_identity, run_all, run_strategy
from hks_sim.op import OpKind

SHAPES = [(3, 2, 2), (12, 4, 4), (7, 3, 3), (6, 3, 3)]


def _cfg(L=3, K=2, D=2, **kw):
    wl = WorkloadConfig(q_towers=L, p_towers=K, num_part_q=D)
    return dataclasses.replace(SimConfig(workload=wl), **kw)


class TestArithmeticIdentity(unittest.TestCase):
    def test_shared_counts_are_identical_modup_only(self):
        for L, K, D in SHAPES:
            ident = arithmetic_identity(run_all(_cfg(L, K, D)))
            for key in ("INTT", "SCALE", "NTT", "bconv_output_towers"):
                values = {s: c[key] for s, c in ident.items()}
                self.assertEqual(len(set(values.values())), 1,
                                 f"L={L},K={K},D={D} 的 {key} 不一致：{values}")

    def test_shared_counts_are_identical_full_hks(self):
        for L, K, D in SHAPES:
            cfg = _cfg(L, K, D, boundary="full_hks", invocation="per_hks")
            ident = arithmetic_identity(run_all(cfg))
            for key in ("INTT", "SCALE", "NTT", "KeyMult", "Accumulate",
                        "bconv_output_towers"):
                values = {s: c[key] for s, c in ident.items()}
                self.assertEqual(len(set(values.values())), 1,
                                 f"L={L},K={K},D={D} 的 {key} 不一致：{values}")

    def test_counts_match_closed_form(self):
        ident = arithmetic_identity(run_all(_cfg()))
        for s, c in ident.items():
            self.assertEqual(c["INTT"], 3, s)
            self.assertEqual(c["NTT"], 7, s)
            self.assertEqual(c["bconv_output_towers"], 7, s)

    def test_bconv_call_counts_differ_as_designed(self):
        """DC/MP 每 digit 一次；OC 每个 (目标塔, digit) 一次。"""
        ident = arithmetic_identity(run_all(_cfg()))
        self.assertEqual(ident["dc"]["BConv"], 2)
        self.assertEqual(ident["mp"]["BConv"], 2)
        self.assertEqual(ident["oc"]["BConv"], 7)   # = L*(D-1) + K*D


class TestBypass(unittest.TestCase):
    def test_native_towers_never_get_bconv_or_ntt(self):
        """旁路的定义：native 塔不产生 BConv/NTT 事件。"""
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(_cfg(), strategy)
            for op in r.ops:
                if op.kind in (OpKind.NTT,) and op.tower is not None:
                    self.assertFalse(
                        r.wl.is_native(op.digit, op.tower),
                        f"{strategy}: 对 native 塔 d{op.digit}.t{op.tower} 做了 NTT")

    def test_bypass_count_equals_q_towers(self):
        """旁路塔总数 = sum_d alpha_d = L。"""
        r = run_strategy(_cfg(), "oc")
        n_ntt = sum(1 for op in r.ops if op.kind is OpKind.NTT)
        self.assertEqual(n_ntt, r.wl.operation_counts()["ntt"])
        self.assertEqual(r.wl.D * r.wl.total_towers - n_ntt, r.wl.L)


class TestDependencies(unittest.TestCase):
    def test_keymult_never_starts_before_its_tower_is_ready(self):
        """Dependency test（方案 §12）：KeyMult 不得早于对应 tower 的 NTT。"""
        cfg = _cfg(boundary="full_hks", invocation="per_hks")
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(cfg, strategy)
            ntt_end = {(op.digit, op.tower): op.end_cycle
                       for op in r.ops if op.kind is OpKind.NTT}
            for op in r.ops:
                if op.kind is not OpKind.KEYMULT:
                    continue
                key = (op.digit, op.tower)
                if key in ntt_end:
                    self.assertGreaterEqual(
                        op.start_cycle, ntt_end[key],
                        f"{strategy}: KeyMult d{op.digit}.t{op.tower} 早于其 NTT")

    def test_accumulate_follows_its_keymult(self):
        cfg = _cfg(boundary="full_hks", invocation="per_hks")
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(cfg, strategy)
            km_end = {(op.digit, op.tower): op.end_cycle
                      for op in r.ops if op.kind is OpKind.KEYMULT}
            for op in r.ops:
                if op.kind is OpKind.ACCUMULATE:
                    self.assertGreaterEqual(op.start_cycle,
                                            km_end[(op.digit, op.tower)])


if __name__ == "__main__":
    unittest.main()
