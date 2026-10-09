#!/usr/bin/env python3
"""Interval-check one exact strict six-block point above the 11/16 threshold.

The three signed phase angles are exact rational decimals in radians.  Masses
and the two free densities are exact real functions of those angles, defined
by the torque and active-degree formulas in ``optimize_three_phase_branch``.
This script evaluates those functions with outward interval arithmetic and
checks strict degree, density, quotient, and transverse inequalities.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ANGLES = (
    "-0.724131674725",
    "-2.359110639373",
    "0.963936529062",
)


def interval(value: str | int):
    return mp.iv.mpf(str(value))


def endpoint_fraction(endpoint) -> tuple[int, int]:
    sign, mantissa, exponent, _ = endpoint
    numerator = -mantissa if sign else mantissa
    if exponent >= 0:
        return numerator << exponent, 1
    return numerator, 1 << (-exponent)


def decimal_bound(endpoint, rounding: str, places: int = 24) -> str:
    numerator, denominator = endpoint_fraction(endpoint)
    with localcontext() as context:
        context.prec = 100
        context.rounding = rounding
        value = Decimal(numerator) / Decimal(denominator)
        return str(value.quantize(Decimal(1).scaleb(-places)))


def bounds(value) -> dict[str, str]:
    lower, upper = value._mpi_
    lower_num, lower_den = endpoint_fraction(lower)
    upper_num, upper_den = endpoint_fraction(upper)
    return {
        "lower": decimal_bound(lower, ROUND_FLOOR),
        "upper": decimal_bound(upper, ROUND_CEILING),
        "lower_exact_binary_fraction": f"{lower_num}/{lower_den}",
        "upper_exact_binary_fraction": f"{upper_num}/{upper_den}",
    }


def strictly_positive(value) -> bool:
    return mp.mpf(value._mpi_[0]) > 0


def strictly_below_one(value) -> bool:
    return mp.mpf(value._mpi_[1]) < 1


def verify_exact_elimination() -> bool:
    """Check six-class torques, degrees, Hessian sign, and eliminations."""
    alpha, beta, gamma = sp.symbols("alpha beta gamma", real=True)
    generic_k, generic_h = sp.symbols("k h", positive=True)
    generic_x, generic_y = sp.symbols("x y", real=True)
    generic_masses = (1, 1, generic_k, generic_k, generic_h, generic_h)
    generic_phases = (alpha, -alpha, beta, -beta, gamma, -gamma)
    generic_weights = sp.Matrix(
        [
            [1, generic_x, generic_y, 0, 1, 1],
            [generic_x, 1, 0, generic_y, 1, 1],
            [generic_y, 0, 1, 1, 0, 1],
            [0, generic_y, 1, 1, 1, 0],
            [1, 1, 0, 1, 1, 0],
            [1, 1, 1, 0, 0, 1],
        ]
    )
    direct_torques = tuple(
        sum(
            generic_masses[j]
            * generic_weights[i, j]
            * sp.sin(generic_phases[j] - generic_phases[i])
            for j in range(6)
        )
        for i in range(6)
    )
    torque_a = (
        -generic_x * sp.sin(2 * alpha)
        - generic_k * generic_y * sp.sin(alpha - beta)
        - 2 * generic_h * sp.cos(gamma) * sp.sin(alpha)
    )
    torque_b = (
        generic_y * sp.sin(alpha - beta)
        - generic_k * sp.sin(2 * beta)
        - generic_h * sp.sin(beta + gamma)
    )
    torque_c = -2 * sp.cos(alpha) * sp.sin(gamma) - generic_k * sp.sin(beta + gamma)
    expected_torques = (
        torque_a,
        -torque_a,
        torque_b,
        -torque_b,
        torque_c,
        -torque_c,
    )
    if not all(
        sp.simplify(sp.trigsimp(direct - expected)) == 0
        for direct, expected in zip(direct_torques, expected_torques)
    ):
        return False
    direct_degrees = tuple(
        sum(generic_masses[j] * generic_weights[i, j] for j in range(6))
        for i in (0, 2, 4)
    )
    expected_degrees = (
        1 + generic_x + generic_k * generic_y + 2 * generic_h,
        generic_y + 2 * generic_k + generic_h,
        2 + generic_k + generic_h,
    )
    if not all(
        sp.simplify(direct - expected) == 0
        for direct, expected in zip(direct_degrees, expected_degrees)
    ):
        return False
    generic_jacobian = sp.zeros(6)
    generic_hessian = sp.zeros(6)
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            coupling = generic_weights[i, j] * sp.cos(
                generic_phases[j] - generic_phases[i]
            )
            generic_jacobian[i, i] -= generic_masses[j] * coupling
            generic_jacobian[i, j] = generic_masses[j] * coupling
            generic_hessian[i, i] += generic_masses[j] * coupling
            generic_hessian[i, j] = (
                -sp.sqrt(generic_masses[i]) * sp.sqrt(generic_masses[j]) * coupling
            )
    if not all(
        sp.simplify(
            generic_hessian[i, j]
            + sp.sqrt(generic_masses[i])
            * generic_jacobian[i, j]
            / sp.sqrt(generic_masses[j])
        )
        == 0
        for i in range(6)
        for j in range(6)
    ):
        return False

    sin_ab = sp.sin(alpha - beta)
    sin_bg = sp.sin(beta + gamma)
    sin_2a = sp.sin(2 * alpha)
    k = -2 * sp.cos(alpha) * sp.sin(gamma) / sin_bg
    y_zero = k * sp.sin(2 * beta) / sin_ab
    y_slope = sin_bg / sin_ab
    x_zero = -k * k * sp.sin(2 * beta) / sin_2a
    x_slope = -(k * sin_bg + 2 * sp.cos(gamma) * sp.sin(alpha)) / sin_2a
    h = -(x_zero + k * y_zero - 1 - k) / (x_slope + k * y_slope + 1)
    x, y = x_zero + h * x_slope, y_zero + h * y_slope
    identities = (
        -2 * sp.cos(alpha) * sp.sin(gamma) - k * sin_bg,
        y * sin_ab - k * sp.sin(2 * beta) - h * sin_bg,
        -x * sin_2a - k * y * sin_ab - 2 * h * sp.cos(gamma) * sp.sin(alpha),
        x + k * y + h - 1 - k,
    )
    return all(sp.simplify(sp.cancel(value)) == 0 for value in identities)


def run(angles: tuple[str, str, str], precision: int) -> dict:
    exact_identity_check = verify_exact_elimination()
    if not exact_identity_check:
        raise RuntimeError("symbolic torque, degree, or Hessian identity failed")
    mp.mp.dps = precision
    mp.iv.dps = precision
    alpha, beta, gamma = map(interval, angles)
    sin_ab = mp.iv.sin(alpha - beta)
    sin_bg = mp.iv.sin(beta + gamma)
    sin_2a = mp.iv.sin(2 * alpha)
    if not all(
        mp.mpf(value._mpi_[0]) * mp.mpf(value._mpi_[1]) > 0
        for value in (sin_ab, sin_bg, sin_2a)
    ):
        raise RuntimeError("phase elimination has an interval containing zero")

    k = -2 * mp.iv.cos(alpha) * mp.iv.sin(gamma) / sin_bg
    y_zero = k * mp.iv.sin(2 * beta) / sin_ab
    y_slope = sin_bg / sin_ab
    x_zero = -(k * k * mp.iv.sin(2 * beta)) / sin_2a
    x_slope = -(k * sin_bg + 2 * mp.iv.cos(gamma) * mp.iv.sin(alpha)) / sin_2a
    degree_intercept = x_zero + k * y_zero - 1 - k
    degree_slope = x_slope + k * y_slope + 1
    if mp.mpf(degree_slope._mpi_[0]) * mp.mpf(degree_slope._mpi_[1]) <= 0:
        raise RuntimeError("active-degree elimination has an interval containing zero")
    h = -degree_intercept / degree_slope
    a = 1 / (2 * (1 + k + h))
    b, c = k * a, h * a
    x = x_zero + h * x_slope
    y = y_zero + h * y_slope
    masses = [a, a, b, b, c, c]
    phases = [alpha, -alpha, beta, -beta, gamma, -gamma]
    zero, one = interval(0), interval(1)
    weights = [
        [one, x, y, zero, one, one],
        [x, one, zero, y, one, one],
        [y, zero, one, one, zero, one],
        [zero, y, one, one, one, zero],
        [one, one, zero, one, one, zero],
        [one, one, one, zero, zero, one],
    ]

    degrees = [
        sum((masses[j] * weights[i][j] for j in range(6)), zero) for i in (0, 2, 4)
    ]
    transverse = [
        sum(
            (
                masses[j] * weights[i][j] * mp.iv.cos(phases[i] - phases[j])
                for j in range(6)
            ),
            zero,
        )
        for i in (0, 2, 4)
    ]

    hessian = [[zero for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            cosine = mp.iv.cos(phases[i] - phases[j])
            coupling = weights[i][j] * cosine
            hessian[i][i] += masses[j] * coupling
            hessian[i][j] = -mp.iv.sqrt(masses[i] * masses[j]) * coupling

    even = [
        [hessian[2 * i][2 * j] + hessian[2 * i][2 * j + 1] for j in range(3)]
        for i in range(3)
    ]
    odd = [
        [hessian[2 * i][2 * j] - hessian[2 * i][2 * j + 1] for j in range(3)]
        for i in range(3)
    ]
    even_trace = sum((even[i][i] for i in range(3)), zero)
    even_pair_sum = sum(
        (
            even[i][i] * even[j][j] - even[i][j] * even[j][i]
            for i, j in ((0, 1), (0, 2), (1, 2))
        ),
        zero,
    )
    odd_first = odd[0][0]
    odd_second = odd[0][0] * odd[1][1] - odd[0][1] * odd[1][0]
    odd_third = (
        odd[0][0] * (odd[1][1] * odd[2][2] - odd[1][2] * odd[2][1])
        - odd[0][1] * (odd[1][0] * odd[2][2] - odd[1][2] * odd[2][0])
        + odd[0][2] * (odd[1][0] * odd[2][1] - odd[1][1] * odd[2][0])
    )
    minimum_degree_above_threshold = [value - interval("0.6875") for value in degrees]
    checks = {
        "exact_six_torque_degree_and_hessian_identities": exact_identity_check,
        "positive_pair_masses": all(strictly_positive(v) for v in (a, b, c)),
        "nonsynchronous_phases": (
            mp.mpf(sin_ab._mpi_[0]) * mp.mpf(sin_ab._mpi_[1]) > 0
        ),
        "partial_densities_in_open_unit_interval": all(
            strictly_positive(v) and strictly_below_one(v) for v in (x, y)
        ),
        "all_degrees_above_11_over_16": all(
            strictly_positive(v) for v in minimum_degree_above_threshold
        ),
        "middle_degree_exceeds_active_degrees": all(
            strictly_positive(degrees[1] - degrees[i]) for i in (0, 2)
        ),
        "even_nonrotation_positive": strictly_positive(even_trace)
        and strictly_positive(even_pair_sum),
        "odd_positive": all(
            strictly_positive(v) for v in (odd_first, odd_second, odd_third)
        ),
        "transverse_positive": all(strictly_positive(v) for v in transverse),
    }
    if not all(checks.values()):
        raise RuntimeError(f"strict interior certificate failed: {checks}")
    return {
        "scope": (
            "One exact, strictly stable six-block step graphon above 11/16. "
            "The phases are rational decimals in radians; masses and two "
            "densities are defined exactly by torque and active-degree "
            "identities. This is a block-model certificate, not a global "
            "optimality result or a finite-graph certificate."
        ),
        "exact_phase_radians": list(angles),
        "interval_decimal_digits": precision,
        "pair_masses": {name: bounds(v) for name, v in zip("abc", (a, b, c))},
        "partial_densities": {"x": bounds(x), "y": bounds(y)},
        "class_degrees": {name: bounds(v) for name, v in zip("ABC", degrees)},
        "degree_minus_11_over_16": {
            name: bounds(v) for name, v in zip("ABC", minimum_degree_above_threshold)
        },
        "even_nonrotation_checks": {
            "trace": bounds(even_trace),
            "sum_of_pairwise_eigenvalue_products": bounds(even_pair_sum),
        },
        "odd_sylvester_checks": {
            "first": bounds(odd_first),
            "second": bounds(odd_second),
            "third": bounds(odd_third),
        },
        "transverse_values": {name: bounds(v) for name, v in zip("ABC", transverse)},
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--angles", nargs=3, default=DEFAULT_ANGLES)
    parser.add_argument("--precision", type=int, default=70)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.precision < 30:
        parser.error("use at least 30 decimal digits of interval precision")
    report = run(tuple(args.angles), args.precision)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
