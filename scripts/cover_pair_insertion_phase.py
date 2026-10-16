#!/usr/bin/env python3
"""Replayable interval-dual cover away from the three duplicated-pair phases.

This is a partial first-order exclusion, not global optimality.  Its JSON
contains each fixed Gram dual, so verification never needs an SDP solver.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path

import certify_pair_insertion_box as cert
import test_pair_insertion as model


def format_decimal(value: Decimal) -> str:
    return format(value, "f")


def initial_intervals(radius: Decimal):
    root = json.loads((model.ROOT / "data/limiting_kkt_certificate.json").read_text())
    phases = sorted(abs(Decimal(root["variables"][i])) for i in (5, 6, 7))
    pi_upper = Decimal("3.1415926535897932384626433832795030")
    endpoints = [Decimal(0)]
    for phase in phases:
        endpoints.extend([phase - radius, phase + radius])
    endpoints.append(pi_upper)
    return [(endpoints[i], endpoints[i + 1]) for i in (0, 2, 4, 6)], phases


def select_dual(lo, hi, cap):
    center = float((lo + hi) / 2)
    record = model.solve_at_phase(center, 0, cap, False, 0, True)
    if "dual_multipliers" not in record:
        return None
    return {"dual_multipliers": record["dual_multipliers"]}


def verify_report(report):
    leaves = [tuple(map(Decimal, leaf["box"])) for leaf in report["certified_leaves"]]
    leaves += [tuple(map(Decimal, box)) for box in report["unresolved_leaves"]]
    leaves.sort()
    initial = [tuple(map(Decimal, box)) for box in report["initial_intervals"]]
    cursor = 0
    for lo, hi in initial:
        position = lo
        while cursor < len(leaves) and leaves[cursor][0] < hi:
            if leaves[cursor][0] != position or leaves[cursor][1] > hi:
                raise RuntimeError("leaf boxes have a gap or overlap")
            position = leaves[cursor][1]
            cursor += 1
        if position != hi:
            raise RuntimeError("leaf boxes do not cover an initial interval")
    if cursor != len(leaves):
        raise RuntimeError("leaf boxes extend outside the initial intervals")
    for leaf in report["certified_leaves"]:
        lo, hi = leaf["box"]
        upper = cert.bound_box(lo, hi, report["dual_selection_cap"], leaf["dual"], True)
        if upper >= 0:
            raise RuntimeError(f"dual no longer certifies box {lo},{hi}: {upper}")
    return True


def main():
    getcontext().prec = 80
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exclude-radius", default="0.02")
    parser.add_argument("--dual-selection-cap", type=float, default=100)
    parser.add_argument("--max-boxes", type=int, default=200)
    parser.add_argument("--out", type=Path, default=model.ROOT / "data/pair_insertion_partial_cover.json")
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    if args.replay:
        report = json.loads(args.replay.read_text())
        verify_report(report)
        print(json.dumps({"replay": "passed", "certified_leaves": len(report["certified_leaves"])}))
        return
    radius = Decimal(args.exclude_radius)
    intervals, phases = initial_intervals(radius)
    if args.resume:
        prior = json.loads(args.resume.read_text())
        if prior["dual_selection_cap"] != args.dual_selection_cap or prior["excluded_support_radius_radians"] != args.exclude_radius:
            raise RuntimeError("resume parameters differ from prior run")
        frontier = [(Decimal(lo), Decimal(hi), 0) for lo, hi in prior["unresolved_leaves"]]
        certified = prior["certified_leaves"]
        failures = prior["dual_selector_failures"]
        processed = prior["processed_boxes"]
    else:
        frontier = [(lo, hi, 0) for lo, hi in intervals]
        certified = []
        failures = 0
        processed = 0
    while frontier and processed < args.max_boxes:
        lo, hi, depth = frontier.pop()
        processed += 1
        selected = select_dual(lo, hi, args.dual_selection_cap)
        upper = None
        if selected is not None:
            upper = cert.bound_box(format_decimal(lo), format_decimal(hi), args.dual_selection_cap, selected, True)
        if upper is not None and upper < 0:
            certified.append({
                "box": [format_decimal(lo), format_decimal(hi)],
                "dual_upper": str(upper),
                "dual": selected,
            })
        else:
            if selected is None:
                failures += 1
            midpoint = (lo + hi) / 2
            frontier.extend([(midpoint, hi, depth + 1), (lo, midpoint, depth + 1)])
        if processed % 50 == 0:
            print(f"processed={processed} certified={len(certified)} frontier={len(frontier)}", flush=True)
    unresolved = [[format_decimal(lo), format_decimal(hi)] for lo, hi, _ in frontier]
    report = {
        "scope": "First-order fourth reflected-pair insertion at the certified six-class endpoint. Derivatives are unrestricted by exact KKT cancellation. Only listed boxes are excluded; three support-phase neighborhoods and unresolved boxes are not excluded.",
        "dual_selection_cap": args.dual_selection_cap,
        "unbounded_derivatives": True,
        "excluded_support_radius_radians": args.exclude_radius,
        "support_phase_centers": [format_decimal(x) for x in phases],
        "initial_intervals": [[format_decimal(lo), format_decimal(hi)] for lo, hi in intervals],
        "processed_boxes": processed,
        "dual_selector_failures": failures,
        "certified_leaves": certified,
        "unresolved_leaves": unresolved,
        "complete_away_from_support": not unresolved,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"processed": processed, "certified": len(certified), "unresolved": len(unresolved), "complete_away_from_support": not unresolved, "out": str(args.out)}))


if __name__ == "__main__":
    main()
