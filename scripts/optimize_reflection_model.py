#!/usr/bin/env python3
"""Optimize the three-class-plus-reflections Kuramoto block model.

This is a local numerical optimizer, not a global or interval certificate.  It
uses continuation in a prescribed strict stability floor.  The three torque
equations are equality constraints; degree, even-sector, odd-sector, and
transverse stability are inequality constraints.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.optimize import minimize


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FLOORS = (1e-5, 1e-7, 1e-9, 1e-11)


def block_matrix(x: float, y: float) -> np.ndarray:
    return np.array(
        [
            [1.0, x, y, 0.0, 1.0, 1.0],
            [x, 1.0, 0.0, y, 1.0, 1.0],
            [y, 0.0, 1.0, 1.0, 0.0, 1.0],
            [0.0, y, 1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 0.0, 0.0, 1.0],
        ]
    )


def pair_masses(logit_a: float, logit_b: float) -> np.ndarray:
    """Return positive `(a,b,c)` with `a+b+c=1/2`."""
    logits = np.array([logit_a, logit_b, 0.0])
    exponentials = np.exp(logits - np.max(logits))
    return 0.5 * exponentials / np.sum(exponentials)


def physical_even_basis(pair_mass: np.ndarray) -> np.ndarray:
    """Orthonormal basis perpendicular to the even rotation vector."""
    rotation = np.sqrt(2.0 * pair_mass)
    first = np.array([rotation[1], -rotation[0], 0.0])
    first /= np.linalg.norm(first)
    second = np.cross(rotation, first)
    return np.column_stack([first, second])


def evaluate(vector: np.ndarray) -> dict[str, Any]:
    logit_a, logit_b, x, y, alpha, beta, gamma, mu = vector
    abc = pair_masses(logit_a, logit_b)
    a, b, c = abc
    masses = np.repeat(abc, 2)
    phases = np.array([alpha, -alpha, beta, -beta, gamma, -gamma])
    weights = block_matrix(x, y)

    degrees = weights @ masses
    torques = np.array(
        [
            np.sum(masses * weights[i] * np.sin(phases - phases[i]))
            for i in range(6)
        ]
    )

    hessian = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            cosine_weight = weights[i, j] * math.cos(phases[i] - phases[j])
            hessian[i, i] += masses[j] * cosine_weight
            hessian[i, j] = -math.sqrt(masses[i] * masses[j]) * cosine_weight

    even = np.array(
        [
            [hessian[2 * i, 2 * j] + hessian[2 * i, 2 * j + 1] for j in range(3)]
            for i in range(3)
        ]
    )
    odd = np.array(
        [
            [hessian[2 * i, 2 * j] - hessian[2 * i, 2 * j + 1] for j in range(3)]
            for i in range(3)
        ]
    )
    even_basis = physical_even_basis(abc)
    physical_even = even_basis.T @ even @ even_basis

    transverse = np.array(
        [
            np.sum(
                masses
                * weights[i]
                * np.cos(phases[i] - phases)
            )
            for i in range(6)
        ]
    )
    order_parameter = abs(
        2.0 * (a * math.cos(alpha) + b * math.cos(beta) + c * math.cos(gamma))
    )

    return {
        "pair_masses": abc,
        "densities": np.array([x, y]),
        "angles": np.array([alpha, beta, gamma]),
        "objective": float(mu),
        "degrees": degrees[[0, 2, 4]],
        "positive_class_torques": torques[[0, 2, 4]],
        "even_physical_eigenvalues": np.linalg.eigvalsh(physical_even),
        "odd_eigenvalues": np.linalg.eigvalsh(odd),
        "transverse_values": transverse[[0, 2, 4]],
        "order_parameter": float(order_parameter),
    }


def initial_vector() -> np.ndarray:
    candidate = json.loads((ROOT / "data/near_boundary_candidate.json").read_text())
    a, b, c = map(float, candidate["pair_class_masses"])
    x = float(candidate["partial_block_densities"]["x"])
    y = float(candidate["partial_block_densities"]["y"])
    alpha, beta, gamma = map(float, candidate["angles"])
    return np.array(
        [math.log(a / c), math.log(b / c), x, y, alpha, beta, gamma, 0.5 + a]
    )


def inequality_values(
    vector: np.ndarray, stability_floor: float, order_parameter_cap: float
) -> np.ndarray:
    state = evaluate(vector)
    return np.concatenate(
        [
            state["pair_masses"] - 1e-8,
            state["degrees"] - state["objective"],
            [state["even_physical_eigenvalues"][0] - stability_floor],
            [state["odd_eigenvalues"][0] - stability_floor],
            state["transverse_values"] - stability_floor,
            [order_parameter_cap - state["order_parameter"]],
        ]
    )


def solve(
    stability_floor: float,
    start: np.ndarray,
    order_parameter_cap: float = 0.95,
) -> tuple[np.ndarray, dict[str, Any]]:
    torque_constraint = {
        "type": "eq",
        "fun": lambda vector: evaluate(vector)["positive_class_torques"],
    }
    stability_constraint = {
        "type": "ineq",
        "fun": lambda vector: inequality_values(
            vector, stability_floor, order_parameter_cap
        ),
    }
    result = minimize(
        lambda vector: -vector[-1],
        start,
        method="SLSQP",
        bounds=[
            (-8.0, 8.0),
            (-8.0, 8.0),
            (0.0, 1.0),
            (0.0, 1.0),
            (-math.pi, math.pi),
            (-math.pi, math.pi),
            (-math.pi, math.pi),
            (0.0, 0.8),
        ],
        constraints=[torque_constraint, stability_constraint],
        options={"ftol": 1e-13, "maxiter": 1500, "disp": False},
    )
    state = evaluate(result.x)
    torque_error = float(np.max(np.abs(state["positive_class_torques"])))
    inequality_minimum = float(
        np.min(inequality_values(result.x, stability_floor, order_parameter_cap))
    )
    if not result.success:
        raise RuntimeError(f"SLSQP failed at floor {stability_floor}: {result.message}")
    if torque_error > 1e-10 or inequality_minimum < -1e-9:
        raise RuntimeError(
            f"invalid optimizer result at floor {stability_floor}: "
            f"torque={torque_error}, inequality={inequality_minimum}"
        )
    diagnostics = {
        "iterations": int(result.nit),
        "function_evaluations": int(result.nfev),
        "maximum_torque_error": torque_error,
        "minimum_constraint_slack": inequality_minimum,
        "message": str(result.message),
    }
    return result.x, diagnostics


def json_floats(values: np.ndarray) -> list[float]:
    return [float(value) for value in values]


def optimize_sequence(floors: Sequence[float]) -> dict[str, Any]:
    vector = initial_vector()
    results = []
    for floor in floors:
        vector, diagnostics = solve(floor, vector)
        state = evaluate(vector)
        results.append(
            {
                "stability_floor": float(floor),
                "mu": state["objective"],
                "pair_masses": json_floats(state["pair_masses"]),
                "partial_block_densities": json_floats(state["densities"]),
                "angles": json_floats(state["angles"]),
                "distinct_degrees": json_floats(state["degrees"]),
                "positive_class_torques": json_floats(
                    state["positive_class_torques"]
                ),
                "even_physical_eigenvalues": json_floats(
                    state["even_physical_eigenvalues"]
                ),
                "odd_eigenvalues": json_floats(state["odd_eigenvalues"]),
                "transverse_values": json_floats(state["transverse_values"]),
                "order_parameter": state["order_parameter"],
                "solver": diagnostics,
            }
        )

    penultimate, final = results[-2:]
    x1, y1 = penultimate["stability_floor"], penultimate["mu"]
    x2, y2 = final["stability_floor"], final["mu"]
    extrapolated = y2 - (y2 - y1) / (x2 - x1) * x2
    if any(
        results[index]["mu"] >= results[index + 1]["mu"]
        for index in range(len(results) - 1)
    ):
        raise RuntimeError("connectivity did not increase as the stability floor fell")
    for result in results:
        floor = result["stability_floor"]
        tolerance = max(1e-10, floor * 1e-3)
        if abs(result["even_physical_eigenvalues"][0] - floor) > tolerance:
            raise RuntimeError("weak even eigenvalue is not an active constraint")
        if abs(result["odd_eigenvalues"][0] - floor) > tolerance:
            raise RuntimeError("weak odd eigenvalue is not an active constraint")
    if not 0.6915 < extrapolated < 0.6916:
        raise RuntimeError("zero-floor extrapolation left the expected local branch")
    return {
        "schema": "three-pair-reflection-optimization-v1",
        "status": "local_numerical_continuation_not_global_certificate",
        "model": {
            "classes": ["A+", "A-", "B+", "B-", "C+", "C-"],
            "independent_pair_classes": 3,
            "free_masses": 2,
            "free_block_densities": 2,
            "free_phase_angles": 3,
            "torque_equalities": 3,
            "objective": "maximize min(dA,dB,dC)",
        },
        "results": results,
        "linear_zero_floor_extrapolation": float(extrapolated),
        "qualifications": [
            "SLSQP establishes only local numerical optima.",
            "The continuation starts from the packaged near-boundary architecture.",
            "No global upper bound or interval enclosure is claimed.",
            "The even rotation mode is removed before imposing stability.",
            "A nonsynchrony constraint excludes the synchronous optimum.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--floors",
        nargs="+",
        type=float,
        default=DEFAULT_FLOORS,
        help="Strict even/odd/transverse Hessian floors used for continuation.",
    )
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if len(args.floors) < 2 or any(floor <= 0 for floor in args.floors):
        raise SystemExit("provide at least two positive stability floors")
    result = optimize_sequence(args.floors)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
