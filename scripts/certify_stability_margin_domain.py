#!/usr/bin/env python3
"""Prototype rigorous common-dual branch-and-bound for the stability margin.

All search boxes use exact rational endpoints in (u,v,wA,wB,wC), with
theta_i = pi*w_i.  Numerical conic solves only select a dual.  The selected
dual is rebuilt from exact decimal rationals, PSD matrices are exact Gram
matrices, and the normalization is imposed as an exact rational identity.
"""

from __future__ import annotations

import argparse
import heapq
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import cvxpy as cp
import mpmath as mp
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import solve_density_sdp as density_model  # noqa: E402


def frac_decimal(value: float, digits: int = 16, nonnegative: bool = False) -> Fraction:
    if nonnegative:
        value = max(0.0, float(value))
    return Fraction(format(float(value), f".{digits}g"))


def numerical_dual(p: np.ndarray, theta: np.ndarray, target: float, cap: float = 1e4):
    """Choose a finite feasible dual, even when the center primal is infeasible."""
    c = density_model.sector_coefficients(p, theta)
    degree = cp.Variable(3, nonneg=True)
    transverse = cp.Variable(3, nonneg=True)
    torque = cp.Variable(3)
    even = cp.Variable((2, 2), PSD=True)
    odd = cp.Variable((3, 3), PSD=True)
    support = cp.Variable(9, nonneg=True)
    constraints = [
        cp.sum(transverse) + cp.trace(even) + cp.trace(odd) == 1,
        degree <= cap,
        torque <= cap,
        torque >= -cap,
    ]
    for k in range(9):
        coefficient = (
            degree @ c["degree"][:, k]
            + torque @ c["torque"][:, k]
            + cp.trace(even @ c["even"][k])
            + cp.trace(odd @ c["odd"][k])
            + transverse @ c["transverse"][:, k]
        )
        constraints.append(support[k] >= coefficient)
    objective = degree @ (p - target) + transverse @ p + cp.sum(support)
    problem = cp.Problem(cp.Minimize(objective), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=2e-9,
        tol_gap_rel=2e-9,
        tol_feas=2e-9,
        max_iter=300,
    )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"dual selector failed: {problem.status}")
    return {
        "value": float(problem.value),
        "degree": np.asarray(degree.value, dtype=float),
        "transverse": np.asarray(transverse.value, dtype=float),
        "torque": np.asarray(torque.value, dtype=float),
        "even": np.asarray(even.value, dtype=float),
        "odd": np.asarray(odd.value, dtype=float),
    }


def psd_factor(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2)
    return eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))


def exact_dual(numerical):
    """Reconstruct exact rational cone variables with exact normalization."""
    degree_raw = [frac_decimal(x, nonnegative=True) for x in numerical["degree"]]
    torque_raw = [frac_decimal(x) for x in numerical["torque"]]
    transverse_raw = [
        frac_decimal(x, nonnegative=True) for x in numerical["transverse"]
    ]
    even_factor = [
        [frac_decimal(x) for x in row] for row in psd_factor(numerical["even"])
    ]
    odd_factor = [
        [frac_decimal(x) for x in row] for row in psd_factor(numerical["odd"])
    ]

    def gram(factor):
        return [
            [
                sum(
                    (factor[i][k] * factor[j][k] for k in range(len(factor[0]))),
                    Fraction(0),
                )
                for j in range(len(factor))
            ]
            for i in range(len(factor))
        ]

    even_raw, odd_raw = gram(even_factor), gram(odd_factor)
    normalization = (
        sum(transverse_raw, Fraction(0))
        + sum((even_raw[i][i] for i in range(2)), Fraction(0))
        + sum((odd_raw[i][i] for i in range(3)), Fraction(0))
    )
    if normalization <= 0:
        raise RuntimeError("zero cone normalization")
    # The dual variables form one ray.  Scaling only the PSD/transverse
    # variables would change that ray after rational Gram reconstruction, so
    # degree and torque multipliers must be scaled by the same exact factor.
    degree = [x / normalization for x in degree_raw]
    torque = [x / normalization for x in torque_raw]
    transverse = [x / normalization for x in transverse_raw]
    even = [[x / normalization for x in row] for row in even_raw]
    odd = [[x / normalization for x in row] for row in odd_raw]
    check = (
        sum(transverse, Fraction(0))
        + sum((even[i][i] for i in range(2)), Fraction(0))
        + sum((odd[i][i] for i in range(3)), Fraction(0))
    )
    assert check == 1
    return {
        "degree": degree,
        "torque": torque,
        "transverse": transverse,
        "even": even,
        "odd": odd,
    }


def iv_fraction(x: Fraction):
    return mp.iv.mpf(x.numerator) / mp.iv.mpf(x.denominator)


def iv_box(bounds):
    lo, hi = bounds
    lo_iv, hi_iv = iv_fraction(lo), iv_fraction(hi)
    # Each rational endpoint is itself evaluated outward; take the outer
    # endpoints again when forming the box interval.
    return mp.iv.mpf([mp.mpf(lo_iv._mpi_[0]), mp.mpf(hi_iv._mpi_[1])])


def iv_upper(x):
    return mp.mpf(x._mpi_[1])


def iv_positive_part(x):
    """Interval extension of max(0,x), retaining outward arithmetic."""
    lo, hi = mp.mpf(x._mpi_[0]), mp.mpf(x._mpi_[1])
    if hi <= 0:
        return mp.iv.mpf(0)
    if lo >= 0:
        return x
    return mp.iv.mpf([mp.mpf(0), hi])


def interval_bound(box, dual, target: Fraction):
    """Rigorous upper bound for one exact outer box and one fixed exact dual."""
    u, v, wa, wb, wc = [iv_box(b) for b in box]
    one = mp.iv.mpf(1)
    pa_c = (one - u) * (one - v) / 6
    pa_b = pa_c + (one - u) * v / 4
    pa_a = pa_b + u / 2
    p3 = [pa_a, pa_b, pa_c]
    masses = [pa_a, pa_a, pa_b, pa_b, pa_c, pa_c]
    angles = [mp.iv.pi * wa, mp.iv.pi * wb, mp.iv.pi * wc]
    phases = [angles[0], -angles[0], angles[1], -angles[1], angles[2], -angles[2]]

    # This ordered parameterization always has pA+pB bounded away from zero.
    ra, rb, rc = [mp.iv.sqrt(2 * p) for p in p3]
    scale = mp.iv.sqrt(ra**2 + rb**2)
    first = [rb / scale, -ra / scale, mp.iv.mpf(0)]
    second = [rc * ra / scale, rc * rb / scale, -scale]
    basis = [[first[i], second[i]] for i in range(3)]

    degree_cols, torque_cols, even_cols, odd_cols, transverse_cols = [], [], [], [], []
    for weights in density_model.DENSITY_BASIS:
        degree_cols.append(
            [
                sum(
                    (int(weights[2 * i, j]) * masses[j] for j in range(6)), mp.iv.mpf(0)
                )
                for i in range(3)
            ]
        )
        torque_cols.append(
            [
                sum(
                    (
                        masses[j]
                        * int(weights[2 * i, j])
                        * mp.iv.sin(phases[j] - phases[2 * i])
                        for j in range(6)
                    ),
                    mp.iv.mpf(0),
                )
                for i in range(3)
            ]
        )
        h = [[mp.iv.mpf(0) for _ in range(6)] for _ in range(6)]
        for i in range(6):
            for j in range(6):
                if i == j or weights[i, j] == 0:
                    continue
                cosine = mp.iv.cos(phases[i] - phases[j])
                h[i][i] += masses[j] * cosine
                h[i][j] = -mp.iv.sqrt(masses[i] * masses[j]) * cosine
        even = [
            [h[2 * i][2 * j] + h[2 * i][2 * j + 1] for j in range(3)] for i in range(3)
        ]
        odd = [
            [h[2 * i][2 * j] - h[2 * i][2 * j + 1] for j in range(3)] for i in range(3)
        ]
        even_cols.append(
            [
                [
                    sum(
                        (
                            basis[i][l] * even[i][j] * basis[j][r]
                            for i in range(3)
                            for j in range(3)
                        ),
                        mp.iv.mpf(0),
                    )
                    for r in range(2)
                ]
                for l in range(2)
            ]
        )
        odd_cols.append(odd)
        transverse_cols.append(
            [
                sum(
                    (
                        masses[j]
                        * int(weights[2 * i, j])
                        * mp.iv.cos(phases[2 * i] - phases[j])
                        for j in range(6)
                    ),
                    mp.iv.mpf(0),
                )
                for i in range(3)
            ]
        )

    degree = [iv_fraction(x) for x in dual["degree"]]
    torque = [iv_fraction(x) for x in dual["torque"]]
    transverse = [iv_fraction(x) for x in dual["transverse"]]
    even_dual = [[iv_fraction(x) for x in row] for row in dual["even"]]
    odd_dual = [[iv_fraction(x) for x in row] for row in dual["odd"]]
    target_iv = iv_fraction(target)

    def inner(left, right):
        return sum(
            (
                left[i][j] * right[i][j]
                for i in range(len(left))
                for j in range(len(left))
            ),
            mp.iv.mpf(0),
        )

    constant = sum(
        (degree[i] * (p3[i] - target_iv) + transverse[i] * p3[i] for i in range(3)),
        mp.iv.mpf(0),
    )
    # Keep the entire support-function sum in interval arithmetic.  Summing
    # scalar upper endpoints would not itself be directed-rounding safe.
    result = constant
    for k in range(9):
        coefficient = sum(
            (
                degree[i] * degree_cols[k][i]
                + torque[i] * torque_cols[k][i]
                + transverse[i] * transverse_cols[k][i]
                for i in range(3)
            ),
            mp.iv.mpf(0),
        )
        coefficient += inner(even_dual, even_cols[k]) + inner(odd_dual, odd_cols[k])
        result += iv_positive_part(coefficient)
    return iv_upper(result)


def midpoint(bounds):
    return (bounds[0] + bounds[1]) / 2


def center_outer(box):
    u, v, wa, wb, wc = [float(midpoint(b)) for b in box]
    c = (1 - u) * (1 - v) / 6
    b = c + (1 - u) * v / 4
    a = b + u / 2
    return np.array([a, b, c]), math.pi * np.array([wa, wb, wc])


def can_straddle(box):
    angle_boxes = box[2:]
    half = Fraction(1, 2)
    return min(b[0] for b in angle_boxes) <= half <= max(b[1] for b in angle_boxes)


def split_box(box):
    widths = [b[1] - b[0] for b in box]
    index = max(range(5), key=lambda i: widths[i])
    middle = midpoint(box[index])
    left, right = list(box), list(box)
    left[index] = (box[index][0], middle)
    right[index] = (middle, box[index][1])
    return tuple(left), tuple(right)


def fraction_json(value):
    return {
        "fraction": f"{value.numerator}/{value.denominator}",
        "decimal": mp.nstr(mp.mpf(value.numerator) / value.denominator, 20),
    }


def box_json(box):
    return [{"lower": fraction_json(lo), "upper": fraction_json(hi)} for lo, hi in box]


def box_volume(box):
    return math.prod((hi - lo for lo, hi in box), start=Fraction(1))


def run(args):
    mp.mp.dps = args.interval_dps
    mp.iv.dps = args.interval_dps
    certificate = json.loads((REPO / "data/limiting_kkt_certificate.json").read_text())
    target = Fraction(certificate["mu"]) + Fraction(args.target_offset)
    root = tuple((Fraction(0), Fraction(1)) for _ in range(5))
    heap = []
    counter = 0
    stats = {
        "processed": 0,
        "pruned_dual": 0,
        "rejected_semicircle": 0,
        "selector_failures": 0,
        "depth_limited": 0,
    }
    worst_pruned = None
    pruned_volume = Fraction(0)
    rejected_volume = Fraction(0)
    start = time.perf_counter()

    def evaluate(box, depth):
        nonlocal counter, worst_pruned, pruned_volume, rejected_volume
        if not can_straddle(box):
            stats["rejected_semicircle"] += 1
            rejected_volume += box_volume(box)
            return
        p, theta = center_outer(box)
        try:
            selected = numerical_dual(p, theta, float(target), args.dual_cap)
            exact = exact_dual(selected)
            upper = interval_bound(box, exact, target)
        except Exception:
            stats["selector_failures"] += 1
            upper = mp.inf
        if upper <= 0:
            stats["pruned_dual"] += 1
            pruned_volume += box_volume(box)
            if worst_pruned is None or upper > worst_pruned:
                worst_pruned = upper
            return
        counter += 1
        priority = -depth if args.strategy == "depth" else -float(upper)
        heapq.heappush(heap, (priority, counter, depth, box, str(upper)))

    evaluate(root, 0)
    frozen = []
    while heap and stats["processed"] < args.max_nodes:
        _, _, depth, box, _ = heapq.heappop(heap)
        stats["processed"] += 1
        if depth >= args.max_depth:
            stats["depth_limited"] += 1
            # Keep a depth-limited leaf in the final partition, but do not
            # reinsert it into the priority queue.  Reinsertion can repeatedly
            # select the same deepest box and consume the node budget without
            # examining any other branch.
            frozen.append((0.0, 0, depth, box, "depth_limit"))
            continue
        for child in split_box(box):
            evaluate(child, depth + 1)

    unresolved = sorted(heap + frozen, key=lambda item: item[0])
    unresolved_volume = sum((box_volume(item[3]) for item in unresolved), Fraction(0))
    accounted_volume = pruned_volume + rejected_volume + unresolved_volume
    assert accounted_volume == 1, (
        "branch partition lost or duplicated volume: "
        f"{accounted_volume.numerator}/{accounted_volume.denominator}"
    )
    elapsed = time.perf_counter() - start
    return {
        "scope": "Prototype rigorous dual pruning for target mu_star+offset; unresolved boxes are not certified.",
        "mu_star_decimal": certificate["mu"],
        "target_offset": args.target_offset,
        "target_exact_fraction": f"{target.numerator}/{target.denominator}",
        "certification_rule": "Only interval common-dual upper <= 0 is pruned; this proves no density with degree>=target and strictly positive common stability margin in that box.",
        "domain": "Exact rational u,v,wA,wB,wC in [0,1], theta=pi*w; boxes unable to straddle pi/2 are rejected.",
        "limits": {
            "max_nodes": args.max_nodes,
            "max_depth": args.max_depth,
            "interval_dps": args.interval_dps,
            "dual_cap": args.dual_cap,
            "strategy": args.strategy,
        },
        "stats": stats,
        "exact_volume_accounting": {
            "root": fraction_json(Fraction(1)),
            "dual_pruned": fraction_json(pruned_volume),
            "semicircle_rejected": fraction_json(rejected_volume),
            "unresolved": fraction_json(unresolved_volume),
            "sum": fraction_json(accounted_volume),
            "identity": "dual_pruned + semicircle_rejected + unresolved = root exactly over Fraction endpoints",
            "full_straddling_domain_volume": fraction_json(Fraction(3, 4)),
            "note": "Unresolved volume can still contain boxes later rejected by the semicircle test; dual-pruned/root is therefore only a conservative completed-coverage measure.",
        },
        "remaining_unresolved_boxes": len(unresolved),
        "maximum_depth_remaining": max((x[2] for x in unresolved), default=None),
        "worst_pruned_upper": (
            None if worst_pruned is None else mp.nstr(worst_pruned, 20)
        ),
        "elapsed_seconds": elapsed,
        "unresolved_preview": [
            {"depth": d, "selector_upper": u, "box": box_json(b)}
            for _, _, d, b, u in unresolved[: args.preview]
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-offset", default="1e-8")
    parser.add_argument("--max-nodes", type=int, default=300)
    parser.add_argument("--max-depth", type=int, default=12)
    parser.add_argument("--interval-dps", type=int, default=60)
    parser.add_argument("--dual-cap", type=float, default=1e4)
    parser.add_argument("--strategy", choices=("worst", "depth"), default="depth")
    parser.add_argument("--preview", type=int, default=12)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = run(args)
    text = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
