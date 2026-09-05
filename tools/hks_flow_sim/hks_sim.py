#!/usr/bin/env python3
"""HKS 三种 dataflow 同构硬件软件仿真器 —— CLI 入口。

    python3 tools/hks_flow_sim/hks_sim.py \
        --config configs/p4_reproduction.yaml \
        --strategies dc,mp,oc \
        --output results/p4_repro

扫参：
    --sweep workload.ring_dimension=4096,8192,16384
    --sweep hardware.sram_capacity_bytes=1048576,4194304

每个扫描点都会完整跑一遍全部策略（方案 §4.2 的硬性要求）。
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from typing import Any, Dict, List, Sequence

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from hks_sim.config import ConfigError, SimConfig, sim_config_from_dict, with_overrides
from hks_sim.engine import (STRATEGIES, FairnessError, oc_tile_strategies,
                            parse_strategy, run_all)
from hks_sim.report import format_summary, write_results, write_summary


def _read(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:  # pragma: no cover
            raise SystemExit("读取 YAML 需要 pyyaml；或改用 .json 配置")
        return yaml.safe_load(text) or {}
    return json.loads(text)


def _deep_merge(base: Dict[str, Any], over: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _load_config(paths: Sequence[str]) -> SimConfig:
    """多个 --config 按顺序深合并。用法：base 配置 + keymult 档位配置。"""
    if not paths:
        return SimConfig()
    data: Dict[str, Any] = {}
    for p in paths:
        data = _deep_merge(data, _read(p))
    return sim_config_from_dict(data)


def _coerce(text: str) -> Any:
    low = text.strip().lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("none", "null"):
        return None
    for cast in (int, float):
        try:
            return cast(text)
        except ValueError:
            pass
    return text


def _parse_sweeps(specs: Sequence[str]) -> List[Dict[str, Any]]:
    """把多个 `key=v1,v2,...` 展开成笛卡尔积。"""
    if not specs:
        return [{}]
    axes = []
    for spec in specs:
        if "=" not in spec:
            raise SystemExit(f"--sweep 需要 key=v1,v2 形式，得到 {spec!r}")
        key, values = spec.split("=", 1)
        axes.append([(key.strip(), _coerce(v)) for v in values.split(",")])
    return [dict(combo) for combo in itertools.product(*axes)]


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", "--workload", dest="config", action="append",
                    default=[],
                    help="YAML/JSON 配置，可重复，按顺序深合并（base + keymult 档位）")
    ap.add_argument("--strategies", default=",".join(STRATEGIES),
                    help=f"逗号分隔，默认 {','.join(STRATEGIES)}；"
                         "OC 可写 oc-w1..oc-w<bconv_cols> 指定 output tile 宽度")
    ap.add_argument("--oc-tile-sweep", action="store_true",
                    help="把 oc-w1..oc-w<bconv_cols> 全部加入对比（M4 Pareto）")
    ap.add_argument("--output", help="输出目录；不给则只打印摘要")
    ap.add_argument("--sweep", action="append", default=[],
                    help="key=v1,v2,...；可重复，多个 --sweep 取笛卡尔积")
    ap.add_argument("--strict-p4", action="store_true",
                    help="只允许 P4 已 RTL 验证的形状，越界即报错")
    ap.add_argument("--enable-tiling", action="store_true")
    ap.add_argument("--emit-trace", action="store_true",
                    help="额外输出逐事件 CSV")
    ap.add_argument("--boundary", choices=("modup_only", "full_hks"))
    ap.add_argument("--invocation",
                    choices=("per_hks", "per_digit", "per_output_tower"))
    args = ap.parse_args(argv)

    strategies = [s.strip() for s in args.strategies.split(",") if s.strip()]
    for s in strategies:
        try:
            parse_strategy(s)
        except KeyError as exc:
            raise SystemExit(str(exc))

    try:
        base = _load_config(args.config)
    except ConfigError as exc:
        raise SystemExit(f"配置错误：{exc}")

    fixed: Dict[str, Any] = {}
    if args.strict_p4:
        fixed["hardware.strict_p4"] = True
    if args.enable_tiling:
        fixed["hardware.allow_tiling"] = True
    if args.boundary:
        fixed["boundary"] = args.boundary
    if args.invocation:
        fixed["invocation"] = args.invocation

    points = _parse_sweeps(args.sweep)
    exit_code = 0
    for i, point in enumerate(points):
        overrides = dict(fixed)
        overrides.update(point)
        try:
            cfg = with_overrides(base, **overrides) if overrides else base
            names = list(strategies)
            if args.oc_tile_sweep:
                # oc-w1..oc-w<bconv_cols>；bconv_cols 可能被 --sweep 改过，
                # 所以要在拿到本点的 cfg 之后再展开
                names = [s for s in names if s != "oc"]
                names += list(oc_tile_strategies(cfg.hardware.bconv_cols))
            results = run_all(cfg, names)
        except (ConfigError, FairnessError) as exc:
            print(f"[point {i}] 失败：{exc}", file=sys.stderr)
            exit_code = 1
            continue

        title = ", ".join(f"{k}={v}" for k, v in point.items()) or "base"
        print(f"\n### sweep point {i}: {title}")
        print(format_summary(results))

        if args.output:
            outdir = args.output if len(points) == 1 else os.path.join(
                args.output, f"point_{i:03d}")
            written = write_results(results, outdir, emit_trace=args.emit_trace)
            write_summary(results, outdir)
            print(f"\n输出：{outdir}  ({', '.join(sorted(written))}, summary.txt)")

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
