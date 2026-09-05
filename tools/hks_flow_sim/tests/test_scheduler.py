"""M2 完成标准：单 engine 配置下，MP 不会因虚假并行获得额外吞吐。

外加 makespan 四分解的自洽性，和 tie-break 的公平性证据。
"""

import dataclasses
import unittest

from hks_sim.config import SimConfig
from hks_sim.engine import run_all, run_strategy
from hks_sim.op import BOOKKEEPING_KINDS, COMPUTE_KINDS


def _cfg(**kw):
    return dataclasses.replace(SimConfig(), **kw)


def _hw(**kw):
    base = SimConfig()
    return dataclasses.replace(base, hardware=dataclasses.replace(base.hardware, **kw))


class TestNoFalseParallelism(unittest.TestCase):
    def test_mp_does_not_beat_dc_on_single_engine(self):
        results = {r.strategy: r for r in run_all(SimConfig())}
        self.assertGreaterEqual(results["mp"].metrics["total_cycles"],
                                results["dc"].metrics["total_cycles"])

    def test_compute_ops_never_overlap_when_dispatcher_is_serial(self):
        """allow_engine_overlap=False 时，任意两个计算事件不得在时间上重叠。"""
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(SimConfig(), strategy)
            spans = sorted((op.start_cycle, op.end_cycle) for op in r.ops
                           if op.kind in COMPUTE_KINDS and op.cycles > 0)
            for (s1, e1), (s2, e2) in zip(spans, spans[1:]):
                self.assertLessEqual(e1, s2,
                                     f"{strategy}: 计算事件重叠 {(s1,e1)} {(s2,e2)}")

    def test_engine_overlap_flag_actually_changes_something(self):
        """这条开关必须真的有效，否则默认值就成了摆设。"""
        serial = run_strategy(SimConfig(), "oc").metrics["total_cycles"]
        overlap = run_strategy(_hw(allow_engine_overlap=True),
                               "oc").metrics["total_cycles"]
        self.assertLess(overlap, serial)

    def test_more_transform_engines_never_slower(self):
        base = run_strategy(SimConfig(), "mp").metrics["total_cycles"]
        wide = run_strategy(_hw(transform_engines=2, allow_engine_overlap=True),
                            "mp").metrics["total_cycles"]
        self.assertLessEqual(wide, base)


class TestMakespanDecomposition(unittest.TestCase):
    def test_four_buckets_sum_to_makespan(self):
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(SimConfig(), strategy)
            m = r.metrics
            total = (m["compute_cycles"] + m["memory_stall_cycles"]
                     + m["control_cycles"] + m["dependency_stall_cycles"])
            self.assertEqual(total, m["total_cycles"], strategy)

    def test_no_bucket_is_negative(self):
        for strategy in ("dc", "mp", "oc"):
            m = run_strategy(SimConfig(), strategy).metrics
            for key in ("compute_cycles", "memory_stall_cycles",
                        "control_cycles", "dependency_stall_cycles"):
                self.assertGreaterEqual(m[key], 0, f"{strategy}.{key}")

    def test_bookkeeping_events_take_no_time(self):
        r = run_strategy(SimConfig(), "dc")
        for op in r.ops:
            if op.kind in BOOKKEEPING_KINDS:
                self.assertEqual(op.start_cycle, op.end_cycle)


class TestTieBreak(unittest.TestCase):
    def test_default_is_trace_order(self):
        self.assertEqual(SimConfig().tie_break, "trace_order")

    def test_all_tie_breaks_run_and_preserve_arithmetic(self):
        """换 tie-break 不能改变做了多少运算，只能改变什么时候做。"""
        from hks_sim.engine import arithmetic_identity
        ref = None
        for tb in ("trace_order", "fifo", "lifo"):
            ident = arithmetic_identity(run_all(_cfg(tie_break=tb)))
            if ref is None:
                ref = ident
            else:
                self.assertEqual(ident, ref, f"tie_break={tb} 改变了算子计数")


class TestDependencyOrder(unittest.TestCase):
    def test_every_op_starts_after_all_its_deps_end(self):
        for strategy in ("dc", "mp", "oc"):
            r = run_strategy(_cfg(boundary="full_hks", invocation="per_hks"),
                             strategy)
            end = {op.seq: op.end_cycle for op in r.ops}
            for op in r.ops:
                for d in op.deps:
                    self.assertGreaterEqual(op.start_cycle, end[d],
                                            f"{strategy}: op{op.seq} 早于依赖 {d}")

    def test_invocations_are_serializing(self):
        """XRT 阻塞模型：下一次调用不得在上一次的任何事件结束前开始。"""
        r = run_strategy(SimConfig(), "dc")
        from hks_sim.op import OpKind
        invokes = [op for op in r.ops
                   if op.kind is OpKind.INVOKE and op.cycles > 0]
        self.assertGreaterEqual(len(invokes), 2)
        for prev, cur in zip(invokes, invokes[1:]):
            earlier_end = max(op.end_cycle for op in r.ops
                              if op.seq < cur.seq and op.cycles > 0)
            self.assertGreaterEqual(cur.start_cycle, earlier_end)


if __name__ == "__main__":
    unittest.main()
