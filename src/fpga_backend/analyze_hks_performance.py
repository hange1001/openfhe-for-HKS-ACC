#!/usr/bin/env python3
"""Compare identical OpenFHE ModUp fixtures in two HLS solutions and CPU JSON.

Reports RTL cycles and nominal-clock estimates, never measured FPGA speedup.
"""
import argparse
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def read_solution(path):
    report = ET.parse(path / "syn/report/Top_csynth.xml").getroot()
    cosim = (path / "sim/report/Top_cosim.rpt").read_text()
    if not re.search(r"\|\s*Verilog\s*\|\s*Pass\s*\|", cosim):
        raise ValueError("Verilog co-simulation did not pass: " + str(path))
    rows = re.findall(r"transaction\s+(\d+):\s+(\d+)\s+(\d+|x)",
                      (path / "sim/report/verilog/result.transaction.rpt").read_text())
    if [int(r[0]) for r in rows] != [0, 1, 2]:
        raise ValueError("Expected fixture transactions INIT, digit[2,0], digit[1,2]")
    cycles = [int(r[1]) for r in rows]
    period = float(report.findtext("UserAssignments/TargetClockPeriod"))
    estimated = float(report.findtext("PerformanceEstimates/SummaryOfTimingAnalysis/EstimatedClockPeriod"))
    uncertainty = float(report.findtext("UserAssignments/ClockUncertainty"))
    rtl = (path / "syn/verilog/Top.v").read_text()
    widths = {port: int(width) for port, width in
              re.findall(r"C_M_AXI_(GMEM\d)_DATA_WIDTH\s*=\s*(\d+)", rtl)}
    return {
        "resources": {v.tag: int(v.text) for v in report.find("AreaEstimates/Resources")},
        "part": report.findtext("UserAssignments/Part"),
        "clock_ns": period, "uncertainty_ns": uncertainty,
        "estimated_period_ns": estimated,
        "budget_slack_ns": round(period - uncertainty - estimated, 3),
        "axi_data_width_bits": widths,
        "transactions": [
            {"opcode": 0, "label": "INIT", "cycles": cycles[0]},
            {"opcode": 8, "alpha": 2, "start": 0, "cycles": cycles[1]},
            {"opcode": 8, "alpha": 1, "start": 2, "cycles": cycles[2]},
        ],
        "two_digit_cycles": sum(cycles[1:]),
        "two_digit_nominal_us": sum(cycles[1:]) * period / 1000,
        "init_nominal_us": cycles[0] * period / 1000,
        "cosim": "PASS",
    }


def analyze(args):
    old = read_solution(args.checkpoint)
    shared = read_solution(args.shared)
    keys = ["clock_ns", "uncertainty_ns", "part"]
    if not args.allow_width_change:
        keys.append("axi_data_width_bits")
    for key in keys:
        if old[key] != shared[key]:
            raise ValueError("Unmatched experiment configuration: " + key)
    cpu = [json.loads(path.read_text()) for path in args.cpu]
    for result in cpu:
        if (result["ring_dim"], result["Q"], result["P"], result["digits"]) != (4096, 3, 2, [2, 1]):
            raise ValueError("CPU workload shape differs from RTL fixture")
    return {
        "scope": "two ModUp digits; cached context; RTL includes testbench AXI, no PCIe/driver/host packing",
        "hardware_performance_measured": False,
        "clock_assumption_met_by_timing_signoff": False,
        "interface_width_change_allowed": args.allow_width_change,
        "checkpoint": old, "shared": shared,
        "resource_delta_shared_minus_checkpoint": {
            key: shared["resources"][key] - value for key, value in old["resources"].items()},
        "shared_latency_increase_percent":
            100 * (shared["two_digit_cycles"] / old["two_digit_cycles"] - 1),
        "cpu_comparison": [
            {"omp_max_threads": r["omp_max_threads"], "cpu_median_us": r["median_us"],
             "cpu_p95_us": r["p95_us"],
             "idealized_cpu_over_checkpoint_rtl_ratio": r["median_us"] / old["two_digit_nominal_us"],
             "idealized_cpu_over_shared_rtl_ratio": r["median_us"] / shared["two_digit_nominal_us"]}
            for r in cpu],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--shared", type=Path, required=True)
    parser.add_argument("--cpu", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--allow-width-change", action="store_true",
                        help="Explicitly compare a bandwidth optimization; still require same clock/device/workload")
    args = parser.parse_args()
    result = json.dumps(analyze(args), indent=2) + "\n"
    print(result, end="")
    if args.output:
        args.output.write_text(result)
