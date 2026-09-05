"""输出：CSV / JSON / 事件 trace / 文本摘要（方案 §9）。

CSV+JSON 是必选输出，只用标准库。matplotlib/pandas 仅在画图时可选引入。
"""

from __future__ import annotations

import csv
import json
import os
from typing import Dict, Iterable, List, Sequence

from .engine import RunResult, arithmetic_identity
from .op import BOOKKEEPING_KINDS
from .resources import BCONV, TRANSFORM

#: §9 要求的指标列，顺序固定，便于跨次运行 diff
METRIC_COLUMNS = [
    "strategy", "hardware_config_hash", "workload_config_hash", "workload",
    "boundary", "invocation_granularity", "invocations", "evidence",
    "total_cycles", "latency_us",
    "compute_cycles", "memory_stall_cycles", "control_cycles",
    "dependency_stall_cycles", "resource_stall_cycles",
    "h2d_bytes", "d2h_bytes", "ddr_read_bytes", "ddr_write_bytes",
    "spill_read_bytes", "spill_write_bytes",
    "peak_live_bytes", "peak_resident_bytes",
    "transform_utilization", "bconv_utilization", "keymult_utilization",
]

TRACE_COLUMNS = [
    "event_id", "strategy", "kind", "digit", "tower", "resource",
    "start_cycle", "end_cycle", "bytes_read", "bytes_written", "stall_reason",
]


def _rows(results: Sequence[RunResult]) -> List[Dict[str, object]]:
    rows = []
    for r in results:
        if not r.feasible:
            rows.append({"strategy": r.strategy, "total_cycles": "infeasible",
                         "hardware_config_hash": r.hardware_hash,
                         "workload_config_hash": r.workload_hash,
                         "evidence": r.infeasible_reason})
            continue
        rows.append(dict(r.metrics))
    return rows


def write_results(results: Sequence[RunResult], outdir: str,
                  emit_trace: bool = False) -> Dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    written: Dict[str, str] = {}

    rows = _rows(results)
    extra = sorted({k for row in rows for k in row} - set(METRIC_COLUMNS))
    columns = METRIC_COLUMNS + extra

    csv_path = os.path.join(outdir, "results.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    written["results.csv"] = csv_path

    json_path = os.path.join(outdir, "results.json")
    payload = {
        "fairness": {
            "hardware_config_hash": results[0].hardware_hash if results else None,
            "workload_config_hash": results[0].workload_hash if results else None,
            "all_equal": len({r.hardware_hash for r in results}) <= 1,
        },
        "arithmetic_identity": arithmetic_identity(results),
        "runs": rows,
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, default=str)
    written["results.json"] = json_path

    if emit_trace:
        for r in results:
            if not r.feasible:
                continue
            path = os.path.join(outdir, f"trace_{r.strategy}.csv")
            with open(path, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=TRACE_COLUMNS)
                w.writeheader()
                for op in r.ops:
                    if op.kind in BOOKKEEPING_KINDS:
                        continue
                    w.writerow({
                        "event_id": op.seq, "strategy": r.strategy,
                        "kind": op.kind.value, "digit": op.digit,
                        "tower": op.tower, "resource": op.resource,
                        "start_cycle": op.start_cycle, "end_cycle": op.end_cycle,
                        "bytes_read": op.bytes_moved if op.reads else 0,
                        "bytes_written": op.bytes_moved if op.writes else 0,
                        "stall_reason": op.stall_reason,
                    })
            written[f"trace_{r.strategy}.csv"] = path
    return written


def format_summary(results: Sequence[RunResult]) -> str:
    if not results:
        return "(无结果)"
    lines: List[str] = []
    head = results[0]
    lines.append(f"workload : {head.wl.describe()}")
    lines.append(f"boundary : {head.cfg.boundary}   "
                 f"invocation: {head.cfg.invocation}   "
                 f"tie_break: {head.cfg.tie_break}")
    hw = head.cfg.hardware
    lines.append(f"hardware : clk={hw.clock_period_ns}ns lanes={hw.transform_lanes} "
                 f"bconv={hw.bconv_rows}x{hw.bconv_cols} "
                 f"dma_overlap={hw.dma_compute_overlap} "
                 f"engine_overlap={hw.allow_engine_overlap} "
                 f"sram={hw.sram_capacity_bytes or 'unbounded'}")
    ok = len({r.hardware_hash for r in results}) == 1
    lines.append(f"fairness : hardware_config_hash={head.hardware_hash} "
                 f"{'OK（三策略一致）' if ok else '不一致 —— 对比无效'}")
    lines.append(f"evidence : {hw.cost.evidence_level()}")
    lines.append("")

    hdr = (f"{'策略':<6}{'cycles':>10}{'latency_us':>12}{'compute':>10}"
           f"{'xfer':>9}{'ctrl':>7}{'idle':>8}{'peak_KB':>10}{'BConv调用':>10}")
    lines.append(hdr)
    lines.append("-" * len(hdr))
    ident = arithmetic_identity(results)
    for r in results:
        if not r.feasible:
            lines.append(f"{r.strategy:<6}{'INFEASIBLE':>10}  {r.infeasible_reason}")
            continue
        m = r.metrics
        lines.append(
            f"{r.strategy:<6}{m['total_cycles']:>10}{m['latency_us']:>12.3f}"
            f"{m['compute_cycles']:>10}{m['memory_stall_cycles']:>9}"
            f"{m['control_cycles']:>7}{m['dependency_stall_cycles']:>8}"
            f"{m['peak_live_bytes'] // 1024:>10}"
            f"{ident[r.strategy]['BConv']:>10}")

    feasible = [r for r in results if r.feasible]
    if len(feasible) > 1:
        best = min(feasible, key=lambda r: r.metrics["total_cycles"])
        lines.append("")
        lines.append(f"最快：{best.strategy}（{best.metrics['total_cycles']} cycles）")
        for r in feasible:
            if r is best:
                continue
            ratio = r.metrics["total_cycles"] / best.metrics["total_cycles"]
            lines.append(f"  {r.strategy} 为其 {ratio:.4f}x")
        lines.append("")
        lines.append("算术恒等式（三策略必须相等，BConv 调用次数除外）：")
        for s, c in ident.items():
            lines.append(f"  {s}: INTT={c['INTT']} SCALE={c['SCALE']} NTT={c['NTT']} "
                         f"KeyMult={c['KeyMult']} Accum={c['Accumulate']} | "
                         f"BConv调用={c['BConv']} 产出塔={c['bconv_output_towers']}")
    return "\n".join(lines)


def write_summary(results: Sequence[RunResult], outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "summary.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(format_summary(results) + "\n")
    return path
