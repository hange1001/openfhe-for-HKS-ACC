#!/usr/bin/env python3
"""Audit generated Vitis 2023.2 RTL hierarchy, not C++ function names.

Usage: python3 check_shared_transform.py Solution/<project>/solution1
This is a structural synthesis check, not equivalence or timing signoff.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET


def audit(solution, lanes, axi_width=None, no_auto=False):
    rtl = solution / "syn" / "verilog"
    modules = {}
    for path in rtl.glob("*.v"):
        source = path.read_text()
        declaration = re.search(r"^\s*module\s+(\w+)\b", source, re.M)
        if declaration:
            modules[declaration[1]] = source
    if "Top" not in modules:
        raise ValueError("No generated Top module in " + str(rtl))
    # Vitis emits the relevant function instances without parameter overrides.
    # Filter by declared module names to exclude ports, tasks and expressions.
    edges = {name: [m[1] for m in re.finditer(r"^\s*(\w+)\s+(\w+)\s*\(", src, re.M)
                    if m[1] in modules]
             for name, src in modules.items()}

    def hierarchy(name, stack=()):
        if name in stack:
            raise ValueError("Recursive RTL hierarchy: " + name)
        counts = Counter({name: 1})
        for child in edges[name]:
            counts.update(hierarchy(child, stack + (name,)))
        return counts

    counts = hierarchy("Top")
    auto_instances = {name: n for name, n in counts.items()
                      if re.search(r"Compute_Auto|Load_Auto_Meta|^Top_Auto(?:_|$)", name)}
    if no_auto and auto_instances:
        raise ValueError("Retired automorphism hardware remains reachable: " + str(auto_instances))
    core = "CG_Transform_Banks" if counts["Top_CG_Transform_Banks"] else "CG_Transform_Kernel"
    chain = ["Top_Execute_Transform_Operation", "Top_Execute_Transform", "Top_" + core]
    for name in chain:
        if counts[name] != 1:
            raise ValueError(f"Expected one {name}, found {counts[name]}")
    old = {name: n for name, n in counts.items()
           if re.search(r"CG_NTT_Kernel|Execute_INTT|Execute_NTT", name)}
    if old:
        raise ValueError("Direction-specialized engines still reachable: " + str(old))
    core_counts = hierarchy(chain[-1])
    multipliers = sum(n for name, n in core_counts.items()
                      if re.fullmatch(r"Top_MultMod(?:_\d+)?", name))
    multiplier_location = "inside transform"
    if multipliers == 0:
        # Vitis may hoist a pure pipelined function to Top and share it with
        # standalone OP_MULT. Count distinct external lane interfaces AND the
        # matching physical pool; do not silently accept a missing multiplier.
        ports = re.findall(r"^output\s+\[63:0\]\s+(grp_(MultMod(?:_\d+)?)_fu_\d+)_p_din1;",
                           modules[chain[-1]], re.M)
        pools = Counter(kind for _, kind in ports)
        if len({name for name, _ in ports}) != lanes or len(pools) != 1:
            raise ValueError("Unrecognized hoisted multiplier interfaces: " + str(ports))
        kind = next(iter(pools))
        physical = edges["Top"].count("Top_" + kind)
        if physical != lanes or counts["Top_" + kind] != lanes:
            raise ValueError("Hoisted physical multiplier pool does not match lane interfaces")
        multipliers = physical
        multiplier_location = "Top: hoisted physical pool, external transform lane interfaces checked"
    if multipliers != lanes:
        raise ValueError(f"Expected {lanes} shared lane multipliers, found {multipliers}")
    transforms = {name: n for name, n in counts.items()
                  if re.fullmatch(r"Top_CG_Transform_(?:Banks|Kernel)(?:_\d+)?", name)}
    if transforms != {chain[-1]: 1}:
        raise ValueError("Additional transform engines/adapters reachable: " + str(transforms))
    total_multipliers = sum(n for name, n in counts.items()
                            if re.fullmatch(r"Top_MultMod(?:_\d+)?", name))
    widths = {port: int(width) for port, width in
              re.findall(r"C_M_AXI_(GMEM\d)_DATA_WIDTH\s*=\s*(\d+)", modules["Top"])}
    if axi_width is not None:
        if widths != {"GMEM0": axi_width, "GMEM1": axi_width, "GMEM2": axi_width}:
            raise ValueError("Unexpected actual AXI widths: " + str(widths))
        copies = {name: n for name, n in counts.items() if "Copy_HKS_Tower" in name}
        if copies:
            raise ValueError("Redundant tower copy engines remain: " + str(copies))

    def report(name):
        root = ET.parse(solution / "syn" / "report" / (name + "_csynth.xml")).getroot()
        resources = {v.tag: int(v.text) for v in root.find("AreaEstimates/Resources")}
        perf = root.find("PerformanceEstimates")
        return root, resources, perf.findtext("SummaryOfOverallLatency/Best-caseLatency")

    top, resources, latency = report("Top")
    _, wrapper_resources, wrapper_latency = report("Execute_Transform")
    _, core_resources, core_latency = report(core)
    butterfly_ii = {}
    for path in (solution / "syn" / "report").glob(core + "_Pipeline_BUTTERFLY_LOOP*_csynth.xml"):
        loops = ET.parse(path).getroot().find("PerformanceEstimates/SummaryOfLoopLatency")
        if loops is not None:
            for loop in loops:
                butterfly_ii[path.stem + ":" + loop.tag] = int(loop.findtext("PipelineII"))
    if len(butterfly_ii) != 2 or any(ii != 1 for ii in butterfly_ii.values()):
        raise ValueError("Expected both butterfly parity loops at II=1: " + str(butterfly_ii))
    target = float(top.findtext("UserAssignments/TargetClockPeriod"))
    uncertainty = float(top.findtext("UserAssignments/ClockUncertainty"))
    estimated = float(top.findtext("PerformanceEstimates/SummaryOfTimingAnalysis/EstimatedClockPeriod"))
    return {
        "status": "PASS", "scope": "Vitis generated RTL structural audit only",
        "chain_instances": {name: counts[name] for name in chain},
        "auto_instances": auto_instances,
        "core_shared_MultMod_instances": multipliers,
        "multiplier_location": multiplier_location,
        "total_MultMod_instances": total_multipliers,
        "axi_data_width_bits": widths,
        "butterfly_loop_ii": butterfly_ii,
        "top_resources": resources, "transform_including_adapters_resources": wrapper_resources,
        "core_resources": core_resources, "core_latency_cycles": core_latency,
        "transform_including_adapters_latency_cycles": wrapper_latency,
        "top_latency_cycles": latency,
        "target_period_ns": target, "uncertainty_ns": uncertainty,
        "estimated_period_ns": estimated,
        "estimated_budget_slack_ns": round(target - uncertainty - estimated, 3),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solution", type=Path)
    parser.add_argument("--lanes", type=int, default=4)
    parser.add_argument("--axi-width", type=int, help="Require all actual AXI data widths and no tower-copy modules")
    parser.add_argument("--total-multipliers", type=int, help="Require unchanged/reduced whole-Top compute capacity")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-auto", action="store_true", help="Reject reachable retired automorphism hardware")
    args = parser.parse_args()
    try:
        result = audit(args.solution, args.lanes, args.axi_width, args.no_auto)
        if args.total_multipliers is not None and result["total_MultMod_instances"] != args.total_multipliers:
            raise ValueError("Unexpected whole-Top multiplier count")
    except (ValueError, OSError, ET.ParseError) as exc:
        parser.exit(1, "FAIL: " + str(exc) + "\n")
    payload = json.dumps(result, indent=2) + "\n"
    print(payload, end="")
    if args.output:
        args.output.write_text(payload)
