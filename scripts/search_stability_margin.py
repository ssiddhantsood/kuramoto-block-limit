#!/usr/bin/env python3
"""Search directly for a counterexample to the proposed three-pair bound.

For fixed pair masses and phases, this script maximizes the common stability
margin over all nine reflection-compatible densities while requiring every
class degree to exceed a target.  A positive optimum is a strictly stable
counterexample at that outer point.  Differential evolution searches the five
outer variables; the inner problem is a convex SDP.

This is a targeted numerical falsification attempt, not a global certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np
from scipy.optimize import differential_evolution

import optimize_general_reflection_model as general
import solve_density_sdp as density_sdp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "data/limiting_kkt_certificate.json"


def ordered_pair_masses(first: float, second: float) -> np.ndarray:
    """Map the unit square onto a>=b>=c>=0 with a+b+c=1/2."""
    largest_gap = first
    middle_weight = (1.0 - first) * second
    smallest_weight = (1.0 - first) * (1.0 - second)
    gap_ab = largest_gap / 2.0
    gap_bc = middle_weight / 4.0
    c = smallest_weight / 6.0
    b = c + gap_bc
    a = b + gap_ab
    return np.array([a, b, c])


def ordered_mass_coordinates(pair_masses: np.ndarray) -> tuple[float, float]:
    a, b, c = pair_masses
    first = 2.0 * (a - b)
    second = 4.0 * (b - c) / (1.0 - first)
    return float(first), float(second)


def phases_straddle_semicircle_boundary(angles: np.ndarray) -> bool:
    """Necessary condition for a connected nonsynchronous equilibrium.

    Independent swaps inside the three reflection pairs allow every angle to
    lie in [0,pi].  The six phases +/-theta fit in an open semicircle exactly
    when all three theta lie strictly on one side of pi/2.  A connected
    equilibrium in an open semicircle is synchronous by the maximum principle.
    """
    return bool(np.min(angles) <= math.pi / 2 <= np.max(angles))


def evaluate_density_solution(
    pair_masses: np.ndarray,
    angles: np.ndarray,
    densities: np.ndarray,
    target_degree: float,
) -> dict[str, Any]:
    if np.min(pair_masses) > 0:
        logits = [
            math.log(pair_masses[0] / pair_masses[2]),
            math.log(pair_masses[1] / pair_masses[2]),
        ]
        vector = np.concatenate((logits, densities, angles, [target_degree]))
        state = general.evaluate(vector)
        degrees = state["degrees"]
        torques = state["positive_class_torques"]
        even_eigenvalues = state["even_physical_eigenvalues"]
        odd_eigenvalues = state["odd_eigenvalues"]
        transverse_values = state["transverse_values"]
        validation_model = "independent direct six-block evaluator"
    else:
        # Boundary points of the ordered mass simplex cannot be represented by
        # finite logits.  Evaluate them directly from the affine sectors.
        coefficients = density_sdp.sector_coefficients(pair_masses, angles)
        degrees = pair_masses + coefficients["degree"] @ densities
        torques = coefficients["torque"] @ densities
        even = sum(
            densities[index] * coefficients["even"][index]
            for index in range(9)
        )
        odd = sum(
            densities[index] * coefficients["odd"][index]
            for index in range(9)
        )
        even_eigenvalues = np.linalg.eigvalsh(even)
        odd_eigenvalues = np.linalg.eigvalsh(odd)
        transverse_values = (
            pair_masses + coefficients["transverse"] @ densities
        )
        validation_model = "affine sector evaluator at zero-mass boundary"
    minimum_stability = float(
        min(
            even_eigenvalues[0],
            odd_eigenvalues[0],
            *transverse_values,
        )
    )
    return {
        "pair_masses": [float(value) for value in pair_masses],
        "angles": [float(value) for value in angles],
        "densities": {
            name: float(value)
            for name, value in zip(general.DENSITY_NAMES, densities)
        },
        "degrees": [float(value) for value in degrees],
        "even_physical_eigenvalues": [
            float(value) for value in even_eigenvalues
        ],
        "odd_eigenvalues": [float(value) for value in odd_eigenvalues],
        "transverse_values": [float(value) for value in transverse_values],
        "minimum_reconstructed_stability": minimum_stability,
        "maximum_torque_error": float(np.max(np.abs(torques))),
        "minimum_degree_slack": float(np.min(degrees - target_degree)),
        "minimum_density_box_slack": float(
            min(np.min(densities), np.min(1.0 - densities))
        ),
        "validation_model": validation_model,
    }


def solve_stability_margin(
    pair_masses: np.ndarray,
    angles: np.ndarray,
    target_degree: float,
) -> dict[str, Any]:
    coefficients = density_sdp.sector_coefficients(pair_masses, angles)
    densities = cp.Variable(9)
    margin = cp.Variable()
    even = sum(
        densities[index] * coefficients["even"][index] for index in range(9)
    )
    odd = sum(
        densities[index] * coefficients["odd"][index] for index in range(9)
    )
    constraints = [
        densities >= 0,
        densities <= 1,
        coefficients["torque"] @ densities == 0,
        pair_masses + coefficients["degree"] @ densities >= target_degree,
        even >> margin * np.eye(2),
        odd >> margin * np.eye(3),
        pair_masses + coefficients["transverse"] @ densities >= margin,
    ]
    problem = cp.Problem(cp.Maximize(margin), constraints)
    try:
        problem.solve(
            solver="CLARABEL",
            tol_gap_abs=2e-9,
            tol_gap_rel=2e-9,
            tol_feas=2e-9,
            max_iter=300,
        )
    except cp.error.SolverError:
        return {"feasible": False, "status": "solver_error"}
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        return {"feasible": False, "status": str(problem.status)}
    density_values = np.asarray(densities.value, dtype=float)
    reconstructed = evaluate_density_solution(
        pair_masses, angles, density_values, target_degree
    )
    return {
        "feasible": True,
        "status": str(problem.status),
        "solver_stability_margin": float(margin.value),
        **reconstructed,
    }


def seeds(certificate_path: Path) -> list[np.ndarray]:
    certificate = json.loads(certificate_path.read_text())
    values = np.array(list(map(float, certificate["variables"])))
    masses = values[:3]
    angles = np.abs(values[5:8])
    order = np.argsort(-masses)
    masses = masses[order]
    angles = angles[order]
    first, second = ordered_mass_coordinates(masses)
    result = [np.array([first, second, *angles])]
    for marginal_masses in (
        np.array([0.25, 0.15, 0.10]),
        np.array([0.25, 0.20, 0.05]),
    ):
        first, second = ordered_mass_coordinates(marginal_masses)
        result.append(
            np.array(
                [first, second, 3 * math.pi / 4, math.pi / 4, math.pi / 4]
            )
        )
    return result


def search(args: argparse.Namespace) -> dict[str, Any]:
    proposed_bound = args.target
    if proposed_bound is None:
        proposed_bound = float(
            json.loads(args.certificate.read_text())["mu"]
        )
    target_degree = proposed_bound + args.degree_excess
    rng = np.random.default_rng(args.seed)
    dimensions = 5
    population_size = max(5, args.population * dimensions)
    initial_population = rng.uniform(0.0, 1.0, size=(population_size, dimensions))
    initial_population[:, 2:] *= math.pi
    initial_seeds = seeds(args.certificate)
    for index, seed in enumerate(initial_seeds):
        initial_population[index] = seed

    evaluations = 0
    feasible_evaluations = 0
    semicircle_rejections = 0
    best_numerical: dict[str, Any] | None = None
    best_validated: dict[str, Any] | None = None

    def is_validated_counterexample(result: dict[str, Any]) -> bool:
        tolerance = args.counterexample_tolerance
        return bool(
            result["minimum_reconstructed_stability"] > tolerance
            and result["minimum_degree_slack"] >= 0.0
            and result["maximum_torque_error"] <= tolerance
            and result["minimum_density_box_slack"] >= -tolerance
        )

    def objective(vector: np.ndarray) -> float:
        nonlocal evaluations, feasible_evaluations, semicircle_rejections
        nonlocal best_numerical, best_validated
        evaluations += 1
        angles = vector[2:]
        if not phases_straddle_semicircle_boundary(angles):
            semicircle_rejections += 1
            return 1.0
        result = solve_stability_margin(
            ordered_pair_masses(vector[0], vector[1]),
            angles,
            target_degree,
        )
        if not result["feasible"]:
            return 1.0
        feasible_evaluations += 1
        candidate = {**result, "outer_coordinates": [float(x) for x in vector]}
        if (
            best_numerical is None
            or result["solver_stability_margin"]
            > best_numerical["solver_stability_margin"]
        ):
            best_numerical = candidate
        if is_validated_counterexample(candidate) and (
            best_validated is None
            or candidate["minimum_reconstructed_stability"]
            > best_validated["minimum_reconstructed_stability"]
        ):
            best_validated = candidate
        return -float(result["solver_stability_margin"])

    optimizer = differential_evolution(
        objective,
        [(0.0, 1.0), (0.0, 1.0), *[(0.0, math.pi)] * 3],
        init=initial_population,
        maxiter=args.iterations,
        seed=args.seed,
        polish=False,
        workers=1,
        updating="immediate",
        atol=1e-9,
        tol=1e-8,
    )
    return {
        "scope": (
            "Numerical counterexample search over the ordered-mass, reflected "
            "three-pair model. A positive validated margin would disprove the "
            "proposed bound; failure to find one is not a global proof."
        ),
        "proposed_bound": proposed_bound,
        "degree_excess_tested": args.degree_excess,
        "target_minimum_degree": target_degree,
        "outer_parameterization": (
            "Two unit-square coordinates cover a>=b>=c>=0, a+b+c=1/2; "
            "three angles lie in [0,pi] and must straddle pi/2."
        ),
        "seed": args.seed,
        "iterations": args.iterations,
        "population_per_dimension": args.population,
        "evaluations": evaluations,
        "sdp_solver_calls": evaluations - semicircle_rejections,
        "feasible_inner_sdps": feasible_evaluations,
        "semicircle_rejections": semicircle_rejections,
        "optimizer_message": str(optimizer.message),
        "counterexample_tolerance": args.counterexample_tolerance,
        "validated_strict_counterexample_found": best_validated is not None,
        "best_numerical": best_numerical,
        "best_validated_counterexample": best_validated,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--target", type=float)
    parser.add_argument("--degree-excess", type=float, default=1e-8)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--population", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--counterexample-tolerance", type=float, default=1e-7)
    parser.add_argument("--out", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = search(args)
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
