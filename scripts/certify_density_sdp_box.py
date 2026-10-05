#!/usr/bin/env python3
"""Produce a rigorous common-dual upper bound on an outer parameter box.

For fixed pair masses and phases, the weak-stability density problem is a
convex SDP.  A feasible dual bounds its optimum from above.  This script solves
for a useful dual numerically at the box center, converts its nonnegative and
PSD variables into exact-decimal feasible objects, and evaluates the resulting
dual objective over the whole mass/phase box with outward-rounded interval
arithmetic.

The numerical SDP selects the certificate but is not trusted for validity.
The reported box bound follows from the reconstructed dual and interval
evaluation.  One box is not a global proof; this is the bounding primitive for
a future branch-and-bound calculation.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

import cvxpy as cp
import mpmath as mp
import numpy as np

import solve_density_sdp as model


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "data/limiting_kkt_certificate.json"


def numerical_dual(
    pair_masses: np.ndarray, angles: np.ndarray, stability_floor: float
) -> dict[str, Any]:
    coefficients = model.sector_coefficients(pair_masses, angles)
    degree_weights = cp.Variable(3, nonneg=True)
    transverse_weights = cp.Variable(3, nonneg=True)
    torque_weights = cp.Variable(3)
    even_dual = cp.Variable((2, 2), PSD=True)
    odd_dual = cp.Variable((3, 3), PSD=True)
    positive_parts = cp.Variable(9, nonneg=True)
    constraints = [cp.sum(degree_weights) == 1]
    for index in range(9):
        coefficient = (
            degree_weights @ coefficients["degree"][:, index]
            + torque_weights @ coefficients["torque"][:, index]
            + cp.trace(even_dual @ coefficients["even"][index])
            + cp.trace(odd_dual @ coefficients["odd"][index])
            + transverse_weights @ coefficients["transverse"][:, index]
        )
        constraints.append(positive_parts[index] >= coefficient)
    objective = (
        degree_weights @ pair_masses
        + transverse_weights @ (pair_masses - stability_floor)
        - stability_floor * (cp.trace(even_dual) + cp.trace(odd_dual))
        + cp.sum(positive_parts)
    )
    problem = cp.Problem(cp.Minimize(objective), constraints)
    problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-11,
        tol_gap_rel=1e-11,
        tol_feas=1e-11,
        max_iter=500,
    )
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"center dual solve failed: {problem.status}")
    return {
        "value": float(problem.value),
        "degree_weights": np.asarray(degree_weights.value, dtype=float),
        "transverse_weights": np.asarray(transverse_weights.value, dtype=float),
        "torque_weights": np.asarray(torque_weights.value, dtype=float),
        "even": np.asarray(even_dual.value, dtype=float),
        "odd": np.asarray(odd_dual.value, dtype=float),
    }


def psd_factor(matrix: np.ndarray) -> np.ndarray:
    """Return R with RR^T PSD; small negative numerical modes are discarded."""
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2.0)
    return eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0.0)))


def decimal_point(value: float, digits: int = 17):
    return mp.iv.mpf(format(float(value), f".{digits}g"))


def exact_decimal(value) -> Decimal:
    if isinstance(value, Decimal):
        return value
    if isinstance(value, str):
        return Decimal(value)
    return Decimal(format(float(value), ".17g"))


def interval(lower, upper):
    return mp.iv.mpf([str(exact_decimal(lower)), str(exact_decimal(upper))])


def lower(value) -> mp.mpf:
    return mp.mpf(value._mpi_[0])


def upper(value) -> mp.mpf:
    return mp.mpf(value._mpi_[1])


def outward_decimal(value: mp.mpf, digits: int, upward: bool) -> str:
    """Round an mp number outward to the requested significant digits."""
    if value == 0:
        return "0"
    exponent = int(mp.floor(mp.log10(abs(value))))
    unit = mp.power(10, exponent - digits + 1)
    scaled = value / unit
    rounded = mp.ceil(scaled) if upward else mp.floor(scaled)
    return mp.nstr(rounded * unit, digits)


def gram(factor):
    return [
        [
            sum(
                (factor[i][k] * factor[j][k] for k in range(len(factor[0]))),
                mp.iv.mpf(0),
            )
            for j in range(len(factor))
        ]
        for i in range(len(factor))
    ]


def matrix_inner(left, right):
    return sum(
        (
            left[i][j] * right[i][j]
            for i in range(len(left))
            for j in range(len(left))
        ),
        mp.iv.mpf(0),
    )


def certify_box(
    center_masses,
    center_angles,
    mass_radius,
    angle_radius,
    stability_floor,
) -> dict[str, Any]:
    getcontext().prec = 80
    center_mass_decimals = [exact_decimal(value) for value in center_masses]
    center_angle_decimals = [exact_decimal(value) for value in center_angles]
    mass_radius_decimal = exact_decimal(mass_radius)
    angle_radius_decimal = exact_decimal(angle_radius)
    stability_floor_decimal = exact_decimal(stability_floor)
    center_mass_floats = np.array(
        [float(value) for value in center_mass_decimals]
    )
    center_angle_floats = np.array(
        [float(value) for value in center_angle_decimals]
    )
    dual = numerical_dual(
        center_mass_floats,
        center_angle_floats,
        float(stability_floor_decimal),
    )
    mp.mp.dps = 70
    mp.iv.dps = 70

    # These are exact decimal point values.  Defining the final weight as the
    # complement makes the degree weights sum to exactly one.
    weight_zero = Decimal(
        format(max(0.0, dual["degree_weights"][0]), ".17g")
    )
    weight_one = Decimal(
        format(max(0.0, dual["degree_weights"][1]), ".17g")
    )
    weight_two = Decimal(1) - weight_zero - weight_one
    if weight_two < 0:
        raise RuntimeError("rounded degree weights are not dual feasible")
    degree_weights = [
        mp.iv.mpf(str(weight_zero)),
        mp.iv.mpf(str(weight_one)),
        mp.iv.mpf(str(weight_two)),
    ]
    transverse_weights = [
        decimal_point(max(0.0, value))
        for value in dual["transverse_weights"]
    ]
    torque_weights = [decimal_point(value) for value in dual["torque_weights"]]
    even_factor = [
        [decimal_point(value) for value in row]
        for row in psd_factor(dual["even"])
    ]
    odd_factor = [
        [decimal_point(value) for value in row]
        for row in psd_factor(dual["odd"])
    ]
    even_dual = gram(even_factor)
    odd_dual = gram(odd_factor)

    mass_a_lower = center_mass_decimals[0] - mass_radius_decimal
    mass_a_upper = center_mass_decimals[0] + mass_radius_decimal
    mass_b_lower = center_mass_decimals[1] - mass_radius_decimal
    mass_b_upper = center_mass_decimals[1] + mass_radius_decimal
    mass_c_lower = Decimal("0.5") - mass_a_upper - mass_b_upper
    mass_c_upper = Decimal("0.5") - mass_a_lower - mass_b_lower
    mass_a = interval(mass_a_lower, mass_a_upper)
    mass_b = interval(mass_b_lower, mass_b_upper)
    mass_c = mp.iv.mpf("0.5") - mass_a - mass_b
    pair_masses = [mass_a, mass_b, mass_c]
    if min(mp.mpf(value._mpi_[0]) for value in pair_masses) <= 0:
        raise ValueError("the mass box must stay strictly inside the simplex")
    masses = [mass_a, mass_a, mass_b, mass_b, mass_c, mass_c]
    angle_endpoints = [
        (
            center_angle_decimals[index] - angle_radius_decimal,
            center_angle_decimals[index] + angle_radius_decimal,
        )
        for index in range(3)
    ]
    angles = [interval(left, right) for left, right in angle_endpoints]
    phases = [
        angles[0],
        -angles[0],
        angles[1],
        -angles[1],
        angles[2],
        -angles[2],
    ]

    # Smooth mass-orthonormal basis perpendicular to the rotation vector.
    rotation_a, rotation_b, rotation_c = [
        mp.iv.sqrt(2 * value) for value in pair_masses
    ]
    scale = mp.iv.sqrt(rotation_a**2 + rotation_b**2)
    first = [rotation_b / scale, -rotation_a / scale, mp.iv.mpf(0)]
    second = [
        rotation_c * rotation_a / scale,
        rotation_c * rotation_b / scale,
        -scale,
    ]
    basis = [[first[index], second[index]] for index in range(3)]

    degree_columns = []
    torque_columns = []
    even_columns = []
    odd_columns = []
    transverse_columns = []
    for weights in model.DENSITY_BASIS:
        degree_columns.append(
            [
                sum(
                    (
                        mp.iv.mpf(int(weights[2 * i, j])) * masses[j]
                        for j in range(6)
                    ),
                    mp.iv.mpf(0),
                )
                for i in range(3)
            ]
        )
        torque_columns.append(
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
        hessian = [[mp.iv.mpf(0) for _ in range(6)] for _ in range(6)]
        for i in range(6):
            for j in range(6):
                if i == j or weights[i, j] == 0:
                    continue
                cosine = mp.iv.cos(phases[i] - phases[j])
                hessian[i][i] += masses[j] * cosine
                hessian[i][j] = -mp.iv.sqrt(masses[i] * masses[j]) * cosine
        even = [
            [
                hessian[2 * i][2 * j] + hessian[2 * i][2 * j + 1]
                for j in range(3)
            ]
            for i in range(3)
        ]
        odd = [
            [
                hessian[2 * i][2 * j] - hessian[2 * i][2 * j + 1]
                for j in range(3)
            ]
            for i in range(3)
        ]
        physical_even = [
            [
                sum(
                    (
                        basis[i][left]
                        * even[i][j]
                        * basis[j][right]
                        for i in range(3)
                        for j in range(3)
                    ),
                    mp.iv.mpf(0),
                )
                for right in range(2)
            ]
            for left in range(2)
        ]
        even_columns.append(physical_even)
        odd_columns.append(odd)
        transverse_columns.append(
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

    floor_interval = mp.iv.mpf(str(stability_floor_decimal))
    constant = sum(
        (
            degree_weights[i] * pair_masses[i]
            + transverse_weights[i] * (pair_masses[i] - floor_interval)
            for i in range(3)
        ),
        mp.iv.mpf(0),
    )
    constant -= floor_interval * (
        sum((even_dual[i][i] for i in range(2)), mp.iv.mpf(0))
        + sum((odd_dual[i][i] for i in range(3)), mp.iv.mpf(0))
    )
    result_interval = constant
    coefficient_upper_bounds = []
    for index in range(9):
        coefficient = sum(
            (
                degree_weights[i] * degree_columns[index][i]
                + torque_weights[i] * torque_columns[index][i]
                + transverse_weights[i] * transverse_columns[index][i]
                for i in range(3)
            ),
            mp.iv.mpf(0),
        )
        coefficient += matrix_inner(even_dual, even_columns[index])
        coefficient += matrix_inner(odd_dual, odd_columns[index])
        coefficient_lower = lower(coefficient)
        coefficient_upper = upper(coefficient)
        coefficient_upper_bounds.append(
            outward_decimal(coefficient_upper, 25, upward=True)
        )
        if coefficient_upper <= 0:
            support = mp.iv.mpf(0)
        elif coefficient_lower >= 0:
            support = coefficient
        else:
            support = mp.iv.mpf([0, coefficient_upper])
        result_interval += support
    result_upper = upper(result_interval)

    return {
        "scope": (
            "Rigorous common-dual upper bound for one mass/phase box in the "
            "weak-stability density SDP; not a global bound by itself."
        ),
        "center_pair_masses": [str(value) for value in center_mass_decimals],
        "center_angles": [str(value) for value in center_angle_decimals],
        "mass_radius": str(mass_radius_decimal),
        "angle_radius": str(angle_radius_decimal),
        "actual_pair_mass_box": [
            {"lower": str(left), "upper": str(right)}
            for left, right in (
                (mass_a_lower, mass_a_upper),
                (mass_b_lower, mass_b_upper),
                (mass_c_lower, mass_c_upper),
            )
        ],
        "actual_angle_box": [
            {"lower": str(left), "upper": str(right)}
            for left, right in angle_endpoints
        ],
        "stability_floor": str(stability_floor_decimal),
        "numerical_center_dual_value_not_used_as_certificate": dual["value"],
        "rigorous_box_upper_bound": outward_decimal(
            result_upper, 30, upward=True
        ),
        "density_coefficient_upper_bounds": coefficient_upper_bounds,
        "dual_feasibility_construction": (
            "Degree weights are nonnegative exact decimals summing to one; "
            "transverse weights are nonnegative; even and odd dual matrices "
            "are exact-decimal Gram matrices; torque weights are free."
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--mass-radius", default="1e-10")
    parser.add_argument("--angle-radius", default="1e-10")
    parser.add_argument("--stability-floor", default="0")
    parser.add_argument("--center-masses", nargs=3)
    parser.add_argument("--center-angles", nargs=3)
    parser.add_argument("--out", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if (args.center_masses is None) != (args.center_angles is None):
        raise ValueError("provide both --center-masses and --center-angles")
    if args.center_masses is None:
        certificate = json.loads(args.certificate.read_text())
        values = certificate["variables"]
        center_masses = values[:3]
        center_angles = values[5:8]
    else:
        center_masses = args.center_masses
        center_angles = args.center_angles
        if sum(exact_decimal(value) for value in center_masses) != Decimal("0.5"):
            raise ValueError("center pair masses must sum to 1/2")
    payload = certify_box(
        center_masses,
        center_angles,
        args.mass_radius,
        args.angle_radius,
        args.stability_floor,
    )
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
