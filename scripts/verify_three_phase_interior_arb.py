#!/usr/bin/env python3
"""Independently verify the strict six-block point with Arb ball arithmetic.

Install the optional dependency with ``pip install -r requirements-arb.txt``.
This script uses only python-flint and the standard library; it does not import
the mpmath certificate or its reconstruction code. The three input strings
denote exact rational numbers in radians.

The torque identities follow algebraically from the reconstruction below.
Writing k=b/a, h=c/a, S=sin(alpha-beta), B=sin(beta+gamma),
A=sin(2 alpha), U=2 sin(alpha) cos(gamma), the three independent torques,
divided by a, are

    T_C/a = -2 cos(alpha) sin(gamma) - k B,
    T_B/a = y S - k sin(2 beta) - h B,
    T_A/a = -x A - k y S - h U.

The definitions of k, y, and x below set these expressions identically to
zero. The definition of h sets (d_A-d_C)/a=x+k*y+h-1-k to zero. Ball-valued
torque evaluations at the end are independent consistency checks, not proofs
of those exact identities.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import flint
from flint import arb, ctx

DEFAULT_ANGLES = (
    "-0.7296280156341436113668706822202972201",
    "-2.3587447703838763461050326004314526713",
    "0.95655249218679516481888457300326302366",
)
DEFAULT_DEGREE_THRESHOLD = "0.691553760567123749076110"


def excludes_zero(value: arb) -> bool:
    return bool(value > 0 or value < 0)


def contains_zero(value: arb) -> bool:
    return bool(value.contains(0))


def enclosure(value: arb) -> dict[str, object]:
    """Record a readable ball and exact dyadic outward endpoints."""
    lower_mantissa, lower_exponent = value.lower().man_exp()
    upper_mantissa, upper_exponent = value.upper().man_exp()
    return {
        "ball": value.str(36),
        "lower_dyadic": [str(lower_mantissa), int(lower_exponent)],
        "upper_dyadic": [str(upper_mantissa), int(upper_exponent)],
    }


def run(
    angle_strings: tuple[str, str, str], bits: int, degree_threshold: str
) -> dict[str, object]:
    ctx.prec = bits
    alpha, beta, gamma = (arb(text) for text in angle_strings)
    sin_ab = (alpha - beta).sin()
    sin_bg = (beta + gamma).sin()
    sin_2a = (2 * alpha).sin()
    if not all(map(excludes_zero, (sin_ab, sin_bg, sin_2a))):
        raise ValueError("a phase denominator contains zero")

    k = -2 * alpha.cos() * gamma.sin() / sin_bg
    y0 = k * (2 * beta).sin() / sin_ab
    y1 = sin_bg / sin_ab
    x0 = -k * k * (2 * beta).sin() / sin_2a
    x1 = -(k * sin_bg + 2 * alpha.sin() * gamma.cos()) / sin_2a
    degree_slope = x1 + k * y1 + 1
    if not excludes_zero(degree_slope):
        raise ValueError("the degree-elimination denominator contains zero")
    h = (1 + k - x0 - k * y0) / degree_slope
    mass_denominator = 2 * (1 + k + h)
    if not excludes_zero(mass_denominator):
        raise ValueError("the mass-normalization denominator contains zero")
    a = 1 / mass_denominator
    b, c = k * a, h * a
    x, y = x0 + h * x1, y0 + h * y1

    masses = [a, a, b, b, c, c]
    phases = [alpha, -alpha, beta, -beta, gamma, -gamma]
    one, zero = arb(1), arb(0)
    weights = [
        [one, x, y, zero, one, one],
        [x, one, zero, y, one, one],
        [y, zero, one, one, zero, one],
        [zero, y, one, one, one, zero],
        [one, one, zero, one, one, zero],
        [one, one, one, zero, zero, one],
    ]

    degrees = [
        sum((masses[j] * weights[i][j] for j in range(6)), arb(0)) for i in range(6)
    ]
    torques = [
        sum(
            (
                masses[j] * weights[i][j] * (phases[j] - phases[i]).sin()
                for j in range(6)
            ),
            arb(0),
        )
        for i in range(6)
    ]
    transverse = [
        sum(
            (
                masses[j] * weights[i][j] * (phases[j] - phases[i]).cos()
                for j in range(6)
            ),
            arb(0),
        )
        for i in range(6)
    ]

    quotient = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            coupling = weights[i][j] * (phases[j] - phases[i]).cos()
            quotient[i][i] = quotient[i][i] + masses[j] * coupling
            quotient[i][j] = -(masses[i] * masses[j]).sqrt() * coupling
    even = [
        [quotient[2 * i][2 * j] + quotient[2 * i][2 * j + 1] for j in range(3)]
        for i in range(3)
    ]
    odd = [
        [quotient[2 * i][2 * j] - quotient[2 * i][2 * j + 1] for j in range(3)]
        for i in range(3)
    ]
    even_trace = sum((even[i][i] for i in range(3)), arb(0))
    even_pair_sum = sum(
        (
            even[i][i] * even[j][j] - even[i][j] * even[j][i]
            for i, j in ((0, 1), (0, 2), (1, 2))
        ),
        arb(0),
    )
    odd_minors = [
        odd[0][0],
        odd[0][0] * odd[1][1] - odd[0][1] * odd[1][0],
        (
            odd[0][0] * (odd[1][1] * odd[2][2] - odd[1][2] * odd[2][1])
            - odd[0][1] * (odd[1][0] * odd[2][2] - odd[1][2] * odd[2][0])
            + odd[0][2] * (odd[1][0] * odd[2][1] - odd[1][1] * odd[2][0])
        ),
    ]
    degree_slacks = [value - arb("11/16") for value in degrees]
    threshold_slacks = [value - arb(degree_threshold) for value in degrees]
    rotation = [mass.sqrt() for mass in masses]
    rotation_residuals = [
        sum((quotient[i][j] * rotation[j] for j in range(6)), arb(0)) for i in range(6)
    ]
    checks = {
        "positive_masses": all(value > 0 for value in (a, b, c)),
        "valid_partial_densities": all(0 < value and value < 1 for value in (x, y)),
        "nonsynchronous": excludes_zero(sin_ab),
        "degree_above_11_over_16": all(value > 0 for value in degree_slacks),
        "degree_above_requested_lower_bound": all(
            value > 0 for value in threshold_slacks
        ),
        "middle_degree_exceeds_active_degrees": (
            degrees[2] > degrees[0] and degrees[2] > degrees[4]
        ),
        "even_nonrotation_positive": even_trace > 0 and even_pair_sum > 0,
        "odd_positive": all(value > 0 for value in odd_minors),
        "transverse_positive": all(value > 0 for value in transverse),
        "rotation_kernel_contains_zero": all(map(contains_zero, rotation_residuals)),
        "direct_torques_contain_zero": all(map(contains_zero, torques)),
        "active_degree_tie_contains_zero": contains_zero(degrees[0] - degrees[4]),
        "mass_sum_contains_one": contains_zero(sum(masses, arb(0)) - 1),
    }
    if not all(checks.values()):
        raise RuntimeError(f"Arb check failed: {checks}")

    return {
        "scope": (
            "Independent Arb inequality check for one exact six-block point. "
            "Exact torque and degree-tie identities follow from the displayed "
            "algebraic reconstruction; ball residuals only check consistency."
        ),
        "python_flint_version": flint.__version__,
        "precision_bits": bits,
        "exact_phase_radians": list(angle_strings),
        "strict_degree_lower_bound": degree_threshold,
        "pair_masses": {
            name: enclosure(value) for name, value in zip("abc", (a, b, c))
        },
        "partial_densities": {
            name: enclosure(value) for name, value in zip("xy", (x, y))
        },
        "class_degrees": {
            name: enclosure(degrees[i]) for name, i in (("A", 0), ("B", 2), ("C", 4))
        },
        "degree_slacks_above_11_over_16": {
            name: enclosure(degree_slacks[i])
            for name, i in (("A", 0), ("B", 2), ("C", 4))
        },
        "even_trace": enclosure(even_trace),
        "even_pair_sum": enclosure(even_pair_sum),
        "odd_leading_minors": [enclosure(value) for value in odd_minors],
        "transverse": {
            name: enclosure(transverse[i]) for name, i in (("A", 0), ("B", 2), ("C", 4))
        },
        "direct_torques": [enclosure(value) for value in torques],
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=512)
    parser.add_argument("--angles", nargs=3, default=DEFAULT_ANGLES)
    parser.add_argument("--degree-threshold", default=DEFAULT_DEGREE_THRESHOLD)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.bits < 128:
        parser.error("use at least 128 bits of precision")
    report = run(tuple(args.angles), args.bits, args.degree_threshold)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
