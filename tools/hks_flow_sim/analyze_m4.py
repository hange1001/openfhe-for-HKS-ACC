#!/usr/bin/env python3
"""M4 分析：output-tile Pareto + 容量边界曲线 + 调度优先级敏感性。

主结论基线固定为 allow_engine_overlap=False / tie_break=trace_order；
engine overlap 与 FIFO/LIFO 只作为敏感性结果输出。

    python3 tools/hks_flow_sim/analyze_m4.py --output results/m4
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import os
import sys
from typing import Dict, List, Sequence

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hks_sim.config import SimConfig, WorkloadConfig
from hks_sim.engine import (RunResult, arithmetic_identity, oc_tile_strategies,
                            run_strategy)
from hks_sim.op import OpKind

KB = 1024


def band(cfg: SimConfig, pipe: int = 38, mem: int = 256, acc: int = 26) -> SimConfig:
    """套上 KeyMult 档位（默认 nominal）。"""
    cost = dataclasses.replace(cfg.hardware.cost, keymult_pipe=pipe,
                               keymult_mem=mem, accumulate_overhead=acc)
    return dataclasses.replace(
        cfg, hardware=dataclasses.replace(cfg.hardware, cost=cost))


def p4_config() -> SimConfig:
    return band(dataclasses.replace(SimConfig(), boundary="full_hks",
                                    invocation="per_hks"))


def industrial_config() -> SimConfig:
    wl = WorkloadConfig(ring_dimension=16384, q_towers=12, p_towers=4,
                        num_part_q=4, label="industrial")
    base = dataclasses.replace(SimConfig(workload=wl), boundary="full_hks",
                               invocation="per_hks")
    return band(base)


def all_strategies(cfg: SimConfig) -> List[str]:
    return ["dc", "mp"] + list(oc_tile_strategies(cfg.hardware.bconv_cols))


def with_capacity(cfg: SimConfig, nbytes: int) -> SimConfig:
    return dataclasses.replace(
        cfg, hardware=dataclasses.replace(cfg.hardware,
                                          sram_capacity_bytes=nbytes))


def bconv_calls(r: RunResult) -> int:
    return sum(1 for op in r.ops if op.kind is OpKind.BCONV)


# ---------------------------------------------------------------- Pareto

def pareto_front(points: Dict[str, tuple]) -> List[str]:
    """返回未被支配的策略名。支配 = latency 与 peak 都不差且至少一项更好。"""
    front = []
    for a, (la, pa) in points.items():
        dominated = any(lb <= la and pb <= pa and (lb < la or pb < pa)
                        for b, (lb, pb) in points.items() if b != a)
        if not dominated:
            front.append(a)
    return front


def run_pareto(cfg: SimConfig, label: str) -> List[dict]:
    rows = []
    results = [run_strategy(cfg, s) for s in all_strategies(cfg)]
    assert len({r.hardware_hash for r in results}) == 1, "hardware hash 不一致"
    ident = arithmetic_identity(results)
    pts = {r.strategy: (r.metrics["total_cycles"], r.mem.peak_live_bytes)
           for r in results}
    front = set(pareto_front(pts))

    print(f"\n=== {label}：output-tile Pareto ===")
    print(f"workload: {results[0].wl.describe()}")
    print(f"hardware_config_hash={results[0].hardware_hash}（全体一致）  "
          f"tower={results[0].wl.bytes_per_tower // KB} KB")
    print(f"{'策略':<8}{'cycles':>10}{'us':>10}{'compute':>10}{'xfer':>9}"
          f"{'BConv调用':>10}{'peak_KB':>9}{'acc塔':>7}  Pareto")
    print("-" * 82)
    for r in results:
        m = r.metrics
        acc_towers = 2 * r.tile_width if r.strategy.startswith("oc") else \
            2 * r.wl.total_towers
        mark = "front" if r.strategy in front else "dominated"
        print(f"{r.strategy:<8}{m['total_cycles']:>10}{m['latency_us']:>10.2f}"
              f"{m['compute_cycles']:>10}{m['memory_stall_cycles']:>9}"
              f"{bconv_calls(r):>10}{m['peak_live_bytes'] // KB:>9}"
              f"{acc_towers:>7}  {mark}")
        rows.append({
            "case": label, "strategy": r.strategy,
            "oc_output_tile_width": r.tile_width,
            "hardware_config_hash": r.hardware_hash,
            "strategy_config_hash": r.strategy_hash,
            "total_cycles": m["total_cycles"], "latency_us": m["latency_us"],
            "compute_cycles": m["compute_cycles"],
            "memory_stall_cycles": m["memory_stall_cycles"],
            "bconv_calls": bconv_calls(r),
            "peak_live_bytes": m["peak_live_bytes"],
            "accumulator_towers": acc_towers,
            "pareto": mark,
        })

    keys = ("INTT", "SCALE", "NTT", "KeyMult", "Accumulate", "bconv_output_towers")
    bad = [k for k in keys if len({ident[r.strategy][k] for r in results}) != 1]
    print(f"算术恒等式：{'全部一致' if not bad else '不一致 ' + str(bad)}")
    print(f"Pareto front：{sorted(front)}")
    return rows


# ---------------------------------------------------------- capacity curve

def run_capacity_curve(cfg: SimConfig, label: str) -> List[dict]:
    """以一个 tower 为步长扫描容量，覆盖各策略 peak 的转折点。"""
    strategies = all_strategies(cfg)
    unbounded = {s: run_strategy(cfg, s) for s in strategies}
    peaks = sorted({r.mem.peak_live_bytes for r in unbounded.values()})
    tower = cfg.workload.bytes_per_tower

    lo = max(tower * 2, peaks[0] - 4 * tower)
    hi = peaks[-1] + 4 * tower
    caps = list(range(lo, hi + tower, tower))

    print(f"\n=== {label}：容量边界扫描 ===")
    print(f"tower={tower // KB} KB，步长 1 tower，范围 "
          f"{lo // KB}..{hi // KB} KB（{len(caps)} 点）")
    print(f"无容量约束时各策略 peak：" +
          ", ".join(f"{s}={unbounded[s].mem.peak_live_bytes // KB}KB"
                    for s in strategies))

    rows = []
    transitions: Dict[str, List[str]] = {}
    for cap in caps:
        c = with_capacity(cfg, cap)
        for s in strategies:
            r = run_strategy(c, s)
            feasible = r.feasible
            rows.append({
                "case": label, "strategy": s, "oc_output_tile_width": r.tile_width,
                "sram_capacity_bytes": cap, "sram_capacity_KB": cap // KB,
                "feasible": feasible,
                "total_cycles": r.metrics["total_cycles"] if feasible else "",
                "spill_read_bytes": r.mem.spill_read_bytes if feasible else "",
                "spill_write_bytes": r.mem.spill_write_bytes if feasible else "",
                "peak_resident_bytes": r.mem.peak_resident_bytes if feasible else "",
            })
            state = ("infeasible" if not feasible
                     else "spill" if r.mem.spill_write_bytes else "clean")
            transitions.setdefault(s, []).append(state)

    print(f"\n{'策略':<8}{'无 spill 的最小容量':>22}{'可行的最小容量':>18}"
          f"{'该容量下 cycles':>18}")
    print("-" * 68)
    for s in strategies:
        states = transitions[s]
        clean = next((caps[i] for i, st in enumerate(states) if st == "clean"), None)
        alive = next((caps[i] for i, st in enumerate(states) if st != "infeasible"),
                     None)
        cyc = ""
        if alive is not None:
            rr = run_strategy(with_capacity(cfg, alive), s)
            cyc = rr.metrics["total_cycles"] if rr.feasible else ""
        print(f"{s:<8}{(str(clean // KB) + ' KB') if clean else 'n/a':>22}"
              f"{(str(alive // KB) + ' KB') if alive else 'n/a':>18}{str(cyc):>18}")
    return rows


# ------------------------------------------------------- tie-break sensitivity

def run_tie_break_sensitivity(cfg: SimConfig, label: str) -> List[dict]:
    strategies = all_strategies(cfg)
    rankings: Dict[str, List[str]] = {}
    rows = []
    print(f"\n=== {label}：调度优先级敏感性 ===")
    for tb in ("trace_order", "fifo", "lifo"):
        c = dataclasses.replace(cfg, tie_break=tb)
        res = [(s, run_strategy(c, s).metrics["total_cycles"]) for s in strategies]
        order = [s for s, _ in sorted(res, key=lambda kv: kv[1])]
        rankings[tb] = order
        for s, cyc in res:
            rows.append({"case": label, "tie_break": tb, "strategy": s,
                         "total_cycles": cyc})
        print(f"  {tb:<12} 排序: {' < '.join(order)}")

    base = rankings["trace_order"]
    flipped = [tb for tb, o in rankings.items() if o != base]
    if flipped:
        print(f"  ⚠ 调度优先级敏感：{flipped} 下排序与基线不同，"
              f"不能只挑一个赢家汇报")
    else:
        print("  排序在三种 tie-break 下一致，结论对调度优先级不敏感")
    return rows


# ------------------------------------------------------ engine-overlap sensitivity

def run_engine_overlap_sensitivity(cfg: SimConfig, label: str) -> List[dict]:
    strategies = all_strategies(cfg)
    rows = []
    print(f"\n=== {label}：allow_engine_overlap 敏感性（仅敏感性，非主结论） ===")
    for flag in (False, True):
        c = dataclasses.replace(
            cfg, hardware=dataclasses.replace(cfg.hardware,
                                              allow_engine_overlap=flag))
        res = [(s, run_strategy(c, s).metrics["total_cycles"]) for s in strategies]
        order = [s for s, _ in sorted(res, key=lambda kv: kv[1])]
        print(f"  overlap={str(flag):<6} 排序: {' < '.join(order)}")
        for s, cyc in res:
            rows.append({"case": label, "allow_engine_overlap": flag,
                         "strategy": s, "total_cycles": cyc})
    return rows


def _dump(rows: Sequence[dict], path: str) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cols = list(rows[0])
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", default="results/m4")
    ap.add_argument("--skip-industrial", action="store_true")
    args = ap.parse_args(argv)

    cases = [("P4 (N=4096,L=3,K=2,D=2)", p4_config())]
    if not args.skip_industrial:
        cases.append(("industrial (N=16384,L=12,K=4,D=4)", industrial_config()))

    pareto, capacity, tie, overlap = [], [], [], []
    for label, cfg in cases:
        pareto += run_pareto(cfg, label)
        capacity += run_capacity_curve(cfg, label)
        tie += run_tie_break_sensitivity(cfg, label)
        overlap += run_engine_overlap_sensitivity(cfg, label)

    _dump(pareto, os.path.join(args.output, "pareto.csv"))
    _dump(capacity, os.path.join(args.output, "capacity_curve.csv"))
    _dump(tie, os.path.join(args.output, "tie_break_sensitivity.csv"))
    _dump(overlap, os.path.join(args.output, "engine_overlap_sensitivity.csv"))
    print(f"\n输出目录：{args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
