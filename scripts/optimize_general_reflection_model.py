#!/usr/bin/env python3
"""Optimize the general three-pair reflection-symmetric Kuramoto step graphon.

The six blocks are ``A+, A-, B+, B-, C+, C-``.  Reflection symmetry leaves
nine independent off-diagonal block densities: three densities joining each
class to its own reflection and, for every pair of letters, one same-sign and
one opposite-sign density.  Diagonal block densities are fixed to one because
they improve degree and transverse stability without changing torque or the
block-constant Hessian sectors.

This is a deterministic multistart local search.  It is numerical evidence,
not a global-optimality or interval certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.optimize import minimize

import optimize_reflection_model as restricted


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FLOORS = (1e-5, 1e-7, 1e-9)
DENSITY_NAMES = (
    "ref_A",
    "ref_B",
    "ref_C",
    "same_AB",
    "opposite_AB",
    "same_AC",
    "opposite_AC",
    "same_BC",
    "opposite_BC",
)


def block_matrix(densities: np.ndarray) -> np.ndarray:
    """Construct the most general reflection-compatible symmetric matrix."""
    ref_a, ref_b, ref_c, s_ab, o_ab, s_ac, o_ac, s_bc, o_bc = densities
    matrix = np.eye(6)

    for pair, value in enumerate((ref_a, ref_b, ref_c)):
        plus = 2 * pair
        minus = plus + 1
        matrix[plus, minus] = matrix[minus, plus] = value

    for left, right, same, opposite in (
        (0, 1, s_ab, o_ab),
        (0, 2, s_ac, o_ac),
        (1, 2, s_bc, o_bc),
    ):
        lp, lm = 2 * left, 2 * left + 1
        rp, rm = 2 * right, 2 * right + 1
        matrix[lp, rp] = matrix[rp, lp] = same
        matrix[lm, rm] = matrix[rm, lm] = same
        matrix[lp, rm] = matrix[rm, lp] = opposite
        matrix[lm, rp] = matrix[rp, lm] = opposite
    return matrix


def evaluate(vector: np.ndarray) -> dict[str, Any]:
    logit_a, logit_b = vector[:2]
    densities = vector[2:11]
    alpha, beta, gamma, mu = vector[11:15]
    abc = restricted.pair_masses(logit_a, logit_b)
    masses = np.repeat(abc, 2)
    phases = np.array([alpha, -alpha, beta, -beta, gamma, -gamma])
    weights = block_matrix(densities)

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
    even_basis = restricted.physical_even_basis(abc)
    physical_even = even_basis.T @ even @ even_basis

    transverse = np.array(
        [
            np.sum(masses * weights[i] * np.cos(phases[i] - phases))
            for i in range(6)
        ]
    )
    order_parameter = abs(
        2.0
        * (
            abc[0] * math.cos(alpha)
            + abc[1] * math.cos(beta)
            + abc[2] * math.cos(gamma)
        )
    )

    return {
        "pair_masses": abc,
        "densities": densities,
        "angles": np.array([alpha, beta, gamma]),
        "objective": float(mu),
        "weights": weights,
        "masses": masses,
        "phases": phases,
        "all_degrees": degrees,
        "degrees": degrees[[0, 2, 4]],
        "all_torques": torques,
        "positive_class_torques": torques[[0, 2, 4]],
        "full_block_hessian": hessian,
        "full_block_eigenvalues": np.linalg.eigvalsh(hessian),
        "even_physical_eigenvalues": np.linalg.eigvalsh(physical_even),
        "odd_eigenvalues": np.linalg.eigvalsh(odd),
        "all_transverse_values": transverse,
        "transverse_values": transverse[[0, 2, 4]],
        "order_parameter": float(order_parameter),
    }


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


def restricted_seed(stability_floor: float) -> np.ndarray:
    """Embed the feasible two-density solution in the nine-density model."""
    old_vector, _ = restricted.solve(
        stability_floor, restricted.initial_vector()
    )
    logit_a, logit_b, x, y, alpha, beta, gamma, mu = old_vector
    densities = np.array([x, 1.0, 0.0, y, 0.0, 1.0, 1.0, 0.0, 1.0])
    return np.concatenate(
        [
            [logit_a, logit_b],
            densities,
            [alpha, beta, gamma, mu],
        ]
    )


def solve(
    stability_floor: float,
    start: np.ndarray,
    order_parameter_cap: float = 0.95,
    max_iterations: int = 4000,
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
            *[(0.0, 1.0)] * 9,
            (-math.pi, math.pi),
            (-math.pi, math.pi),
            (-math.pi, math.pi),
            (0.0, 0.9),
        ],
        constraints=[torque_constraint, stability_constraint],
        options={"ftol": 1e-12, "maxiter": max_iterations, "disp": False},
    )
    state = evaluate(result.x)
    torque_error = float(np.max(np.abs(state["positive_class_torques"])))
    inequality_minimum = float(
        np.min(inequality_values(result.x, stability_floor, order_parameter_cap))
    )
    feasibility_tolerance = max(1e-12, stability_floor * 1e-3)
    strictly_positive_spectra = bool(min(
        state["even_physical_eigenvalues"][0],
        state["odd_eigenvalues"][0],
        *state["transverse_values"],
    ) > 0.0)
    feasible = (
        torque_error <= 1e-9
        and inequality_minimum >= -feasibility_tolerance
        and strictly_positive_spectra
    )
    diagnostics = {
        "success": bool(result.success),
        "feasible": bool(feasible),
        "iterations": int(result.nit),
        "function_evaluations": int(result.nfev),
        "maximum_torque_error": torque_error,
        "minimum_constraint_slack": inequality_minimum,
        "feasibility_tolerance": feasibility_tolerance,
        "strictly_positive_stability_spectra": strictly_positive_spectra,
        "message": str(result.message),
    }
    return result.x, diagnostics


def perturb_start(
    base: np.ndarray, rng: np.random.Generator, scale: float
) -> np.ndarray:
    trial = base.copy()
    trial[:2] += rng.normal(0.0, scale, size=2)
    trial[2:11] = np.clip(
        trial[2:11] + rng.normal(0.0, scale, size=9), 0.0, 1.0
    )
    trial[11:14] = np.clip(
        trial[11:14] + rng.normal(0.0, scale, size=3), -math.pi, math.pi
    )
    trial[-1] = max(0.0, trial[-1] - 2.0 * scale)
    return trial


def json_floats(values: np.ndarray) -> list[float]:
    return [float(value) for value in values]


def optimize_sequence(
    floors: Sequence[float], starts: int, seed: int
) -> dict[str, Any]:
    if starts < 1:
        raise ValueError("starts must be positive")
    rng = np.random.default_rng(seed)
    vector = restricted_seed(float(floors[0]))
    records = []

    for floor_index, floor in enumerate(floors):
        candidates: list[tuple[np.ndarray, dict[str, Any]]] = []
        attempt_diagnostics: list[dict[str, Any]] = []
        starting_points = [vector]
        if floor_index == 0:
            scales = np.geomspace(0.01, 0.40, max(1, starts - 1))
            for index in range(starts - 1):
                scale = float(scales[index])
                starting_points.append(perturb_start(vector, rng, scale))

        for start in starting_points:
            candidate, diagnostics = solve(float(floor), start)
            attempt_diagnostics.append(
                {
                    **diagnostics,
                    "objective": float(candidate[-1]),
                }
            )
            if diagnostics["feasible"]:
                candidates.append((candidate, diagnostics))
        if not candidates:
            raise RuntimeError(f"no feasible optimizer result at floor {floor}")

        vector, diagnostics = max(candidates, key=lambda item: item[0][-1])
        state = evaluate(vector)
        records.append(
            {
                "stability_floor": float(floor),
                "mu": state["objective"],
                "pair_masses": json_floats(state["pair_masses"]),
                "densities": {
                    name: float(value)
                    for name, value in zip(DENSITY_NAMES, state["densities"])
                },
                "angles": json_floats(state["angles"]),
                "degrees": json_floats(state["degrees"]),
                "even_physical_eigenvalues": json_floats(
                    state["even_physical_eigenvalues"]
                ),
                "odd_eigenvalues": json_floats(state["odd_eigenvalues"]),
                "transverse_values": json_floats(state["transverse_values"]),
                "order_parameter": state["order_parameter"],
                "diagnostics": diagnostics,
                "search_summary": {
                    "attempted_starts": len(starting_points),
                    "feasible_starts": sum(
                        bool(item["feasible"]) for item in attempt_diagnostics
                    ),
                    "solver_successes": sum(
                        bool(item["success"]) for item in attempt_diagnostics
                    ),
                    "feasible_objectives": sorted(
                        float(item["objective"])
                        for item in attempt_diagnostics
                        if item["feasible"]
                    ),
                },
            }
        )

    objectives = [record["mu"] for record in records]
    if any(right + 1e-9 < left for left, right in zip(objectives, objectives[1:])):
        raise RuntimeError("objective decreased when the stability floor was relaxed")
    return {
        "model": "general reflection-symmetric three-pair step graphon",
        "density_order": list(DENSITY_NAMES),
        "random_seed": int(seed),
        "requested_multistarts": int(starts),
        "local_search_warning": (
            "Deterministic multistart SLSQP; this does not certify global optimality."
        ),
        "results": records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--floors",
        nargs="+",
        type=float,
        default=DEFAULT_FLOORS,
        help="strict stability floors used for continuation",
    )
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--out", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = optimize_sequence(args.floors, args.starts, args.seed)
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
