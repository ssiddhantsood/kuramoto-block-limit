#!/usr/bin/env python3
"""Solve the convex density layer numerically for fixed masses/phases.

Once the three pair masses and three reflection phases are fixed, torques,
degrees, all quotient Hessian entries, and transverse values are affine in the
nine reflection-compatible densities.  Maximizing minimum degree at fixed
masses/phases is therefore a small semidefinite program, not a nonlinear local
search.  CVXPY/Clarabel returns an approximate floating-point solution; this
script does not provide a certified SDP duality gap.

The optional outer differential-evolution search leaves only five nonlinear
variables (two mass logits and three phases).  It is global-search evidence,
not a proof that the five-dimensional search is exhaustive.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
from scipy.optimize import differential_evolution

import optimize_general_reflection_model as general


ROOT = Path(__file__).resolve().parents[1]


def density_basis() -> list[np.ndarray]:
    identity = np.eye(6)
    basis = []
    for index in range(9):
        densities = np.zeros(9)
        densities[index] = 1.0
        basis.append(general.block_matrix(densities) - identity)
    return basis


DENSITY_BASIS = density_basis()


def sector_coefficients(pair_masses: np.ndarray, angles: np.ndarray):
    masses = np.repeat(pair_masses, 2)
    phases = np.empty(6)
    phases[0::2] = angles
    phases[1::2] = -angles
    even_basis = general.restricted.physical_even_basis(pair_masses)
    degree_columns = []
    torque_columns = []
    even_columns = []
    odd_columns = []
    transverse_columns = []
    for weights in DENSITY_BASIS:
        degree_columns.append((weights @ masses)[0::2])
        torque_columns.append(
            np.array(
                [
                    np.sum(
                        masses
                        * weights[i]
                        * np.sin(phases - phases[i])
                    )
                    for i in range(6)
                ]
            )[0::2]
        )
        hessian = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                if i == j:
                    continue
                cosine_weight = weights[i, j] * math.cos(
                    phases[i] - phases[j]
                )
                hessian[i, i] += masses[j] * cosine_weight
                hessian[i, j] = (
                    -math.sqrt(masses[i] * masses[j]) * cosine_weight
                )
        even = np.array(
            [
                [
                    hessian[2 * i, 2 * j] + hessian[2 * i, 2 * j + 1]
                    for j in range(3)
                ]
                for i in range(3)
            ]
        )
        odd = np.array(
            [
                [
                    hessian[2 * i, 2 * j] - hessian[2 * i, 2 * j + 1]
                    for j in range(3)
                ]
                for i in range(3)
            ]
        )
        even_columns.append(even_basis.T @ even @ even_basis)
        odd_columns.append(odd)
        transverse_columns.append(
            np.array(
                [
                    np.sum(
                        masses
                        * weights[i]
                        * np.cos(phases[i] - phases)
                    )
                    for i in range(6)
                ]
            )[0::2]
        )
    return {
        "degree": np.asarray(degree_columns).T,
        "torque": np.asarray(torque_columns).T,
        "even": even_columns,
        "odd": odd_columns,
        "transverse": np.asarray(transverse_columns).T,
    }


def solve_density_sdp(
    pair_masses: np.ndarray,
    angles: np.ndarray,
    stability_floor: float,
    order_parameter_cap: float = 0.95,
) -> dict[str, Any]:
    order_parameter = abs(2.0 * np.sum(pair_masses * np.cos(angles)))
    if order_parameter > order_parameter_cap:
        return {
            "feasible": False,
            "status": "order_parameter_cap_violated",
            "order_parameter": float(order_parameter),
        }
    coefficients = sector_coefficients(pair_masses, angles)
    densities = cp.Variable(9)
    mu = cp.Variable()
    even = sum(densities[index] * coefficients["even"][index] for index in range(9))
    odd = sum(densities[index] * coefficients["odd"][index] for index in range(9))
    degrees = pair_masses + coefficients["degree"] @ densities
    transverse = pair_masses + coefficients["transverse"] @ densities
    constraints = [
        densities >= 0,
        densities <= 1,
        coefficients["torque"] @ densities == 0,
        degrees >= mu,
        even >> stability_floor * np.eye(2),
        odd >> stability_floor * np.eye(3),
        transverse >= stability_floor,
    ]
    problem = cp.Problem(cp.Maximize(mu), constraints)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            problem.solve(
                solver="CLARABEL",
                tol_gap_abs=1e-9,
                tol_gap_rel=1e-9,
                tol_feas=1e-9,
                max_iter=400,
            )
        except cp.error.SolverError:
            return {"feasible": False, "status": "solver_error"}
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        return {"feasible": False, "status": str(problem.status)}

    density_values = np.asarray(densities.value, dtype=float)
    mu_value = float(mu.value)
    vector = np.concatenate(
        (
            [
                math.log(pair_masses[0] / pair_masses[2]),
                math.log(pair_masses[1] / pair_masses[2]),
            ],
            density_values,
            angles,
            [mu_value],
        )
    )
    state = general.evaluate(vector)
    tolerance = (
        max(1e-12, stability_floor * 5e-3)
        if stability_floor > 0
        else 2e-8
    )
    torque_error = float(np.max(np.abs(state["positive_class_torques"])))
    minimum_stability = float(
        min(
            state["even_physical_eigenvalues"][0],
            state["odd_eigenvalues"][0],
            *state["transverse_values"],
        )
    )
    minimum_degree_slack = float(np.min(state["degrees"] - mu_value))
    density_slack = float(min(np.min(density_values), np.min(1.0 - density_values)))
    strictly_positive_when_requested = bool(
        stability_floor <= 0 or minimum_stability > 0.0
    )
    feasible = bool(
        torque_error <= tolerance
        and minimum_stability >= stability_floor - tolerance
        and strictly_positive_when_requested
        and minimum_degree_slack >= -tolerance
        and density_slack >= -tolerance
    )
    return {
        "feasible": feasible,
        "status": str(problem.status),
        "mu": mu_value,
        "pair_masses": [float(value) for value in pair_masses],
        "angles": [float(value) for value in angles],
        "densities": {
            name: float(value)
            for name, value in zip(general.DENSITY_NAMES, density_values)
        },
        "degrees": [float(value) for value in state["degrees"]],
        "even_physical_eigenvalues": [
            float(value) for value in state["even_physical_eigenvalues"]
        ],
        "odd_eigenvalues": [float(value) for value in state["odd_eigenvalues"]],
        "transverse_values": [float(value) for value in state["transverse_values"]],
        "order_parameter": state["order_parameter"],
        "maximum_torque_error": torque_error,
        "minimum_stability_value": minimum_stability,
        "minimum_degree_slack": minimum_degree_slack,
        "minimum_density_box_slack": density_slack,
        "validation_tolerance": tolerance,
        "strictly_positive_when_requested": strictly_positive_when_requested,
        "model_boundary": (
            "The sectors are complete for the ideal reflected step graphon. "
            "Finite partial-block realizations require separate centered-"
            "adjacency estimates."
        ),
    }


def pair_masses(logit_a: float, logit_b: float) -> np.ndarray:
    return general.restricted.pair_masses(logit_a, logit_b)


def known_outer_seed(stability_floor: float) -> np.ndarray:
    search = json.loads((ROOT / "data/general_reflection_multistart.json").read_text())
    records = search["branches"][0]["continuation"]
    record = min(
        records,
        key=lambda item: abs(float(item["stability_floor"]) - stability_floor),
    )
    a, b, c = map(float, record["pair_masses"])
    return np.array(
        [math.log(a / c), math.log(b / c), *map(float, record["angles"])]
    )


def outer_search(args: argparse.Namespace) -> dict[str, Any]:
    evaluations = 0
    feasible_evaluations = 0
    best: dict[str, Any] | None = None

    def objective(vector: np.ndarray) -> float:
        nonlocal evaluations, feasible_evaluations, best
        evaluations += 1
        state = solve_density_sdp(
            pair_masses(vector[0], vector[1]),
            vector[2:5],
            args.stability_floor,
            args.order_parameter_cap,
        )
        if not state["feasible"]:
            return 1.0
        feasible_evaluations += 1
        if best is None or state["mu"] > best["mu"]:
            best = state
        return -float(state["mu"])

    dimensions = 5
    population_size = max(5, args.population * dimensions)
    rng = np.random.default_rng(args.seed)
    initial_population = np.empty((population_size, dimensions))
    initial_population[:, :2] = rng.uniform(-3.0, 3.0, size=(population_size, 2))
    initial_population[:, 2:] = rng.uniform(
        -math.pi, math.pi, size=(population_size, 3)
    )
    initial_population[0] = known_outer_seed(args.stability_floor)
    initial_population[1] = np.array(
        [math.log(0.1 / 0.15), math.log(0.25 / 0.15), -math.pi / 4, -3 * math.pi / 4, math.pi / 4]
    )
    result = differential_evolution(
        objective,
        [(-3.0, 3.0), (-3.0, 3.0), *[(-math.pi, math.pi)] * 3],
        init=initial_population,
        maxiter=args.search_iterations,
        seed=args.seed,
        polish=False,
        updating="immediate",
        workers=1,
        atol=1e-8,
        tol=1e-7,
    )
    return {
        "scope": (
            "Differential evolution over five outer variables with an "
            "approximately solved convex density SDP at each evaluation. "
            "Neither the SDP values nor the outer search are interval-certified."
        ),
        "outer_bounds": {
            "mass_logits": [[-3.0, 3.0], [-3.0, 3.0]],
            "angles": [[-math.pi, math.pi]] * 3,
        },
        "stability_floor": args.stability_floor,
        "seed": args.seed,
        "iterations": args.search_iterations,
        "population_per_dimension": args.population,
        "evaluations": evaluations,
        "feasible_evaluations": feasible_evaluations,
        "optimizer_message": str(result.message),
        "best": best,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("known", "marginal"), default="known")
    parser.add_argument("--stability-floor", type=float, default=0.0)
    parser.add_argument("--order-parameter-cap", type=float, default=0.95)
    parser.add_argument("--search-iterations", type=int, default=0)
    parser.add_argument("--population", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--out", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.search_iterations:
        payload = outer_search(args)
    elif args.case == "known":
        certificate = json.loads((ROOT / "data/limiting_kkt_certificate.json").read_text())
        variables = list(map(float, certificate["variables"]))
        payload = solve_density_sdp(
            np.array(variables[:3]),
            np.array(variables[5:8]),
            args.stability_floor,
            args.order_parameter_cap,
        )
    else:
        payload = solve_density_sdp(
            np.array([0.1, 0.25, 0.15]),
            np.array([-math.pi / 4, -3 * math.pi / 4, math.pi / 4]),
            args.stability_floor,
            args.order_parameter_cap,
        )
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
