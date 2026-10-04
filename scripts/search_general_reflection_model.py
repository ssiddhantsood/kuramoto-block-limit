#!/usr/bin/env python3
"""Broad multistart search with the three torque equations eliminated.

For fixed masses, cross-pair densities, and phases, each torque equation is
linear in the density joining a class to its own reflection.  Solving those
three equations explicitly removes a difficult equality-constrained layer
from the general nine-density search.  The derived densities are then required
to lie in ``[0,1]``.

This improves exploration and finds multiple stable local branches.  It is not
an exhaustive or global optimization.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterator

import numpy as np
from scipy.optimize import minimize

import optimize_general_reflection_model as general


DEFAULT_FLOORS = (1e-6, 1e-8, 1e-10)


def lift(vector: np.ndarray) -> np.ndarray:
    """Insert the three reflection-pair densities determined by torque zero."""
    pair_masses = general.restricted.pair_masses(vector[0], vector[1])
    cross_densities = np.asarray(vector[2:8])
    angles = np.asarray(vector[8:11])
    s_ab, t_ab, s_ac, t_ac, s_bc, t_bc = cross_densities
    cross = {
        (0, 1): (s_ab, t_ab),
        (0, 2): (s_ac, t_ac),
        (1, 2): (s_bc, t_bc),
    }
    cross_torque = np.zeros(3)
    for (left, right), (same, opposite) in cross.items():
        cross_torque[left] += pair_masses[right] * (
            same * math.sin(angles[right] - angles[left])
            + opposite * math.sin(-angles[right] - angles[left])
        )
        cross_torque[right] += pair_masses[left] * (
            same * math.sin(angles[left] - angles[right])
            + opposite * math.sin(-angles[left] - angles[right])
        )
    denominator = pair_masses * np.sin(2.0 * angles)
    with np.errstate(divide="ignore", invalid="ignore"):
        reflected_densities = cross_torque / denominator
    densities = np.concatenate((reflected_densities, cross_densities))
    return np.concatenate(
        ([vector[0], vector[1]], densities, angles, [vector[-1]])
    )


def inequality_values(
    vector: np.ndarray, stability_floor: float, order_parameter_cap: float
) -> np.ndarray:
    full_vector = lift(vector)
    if not np.all(np.isfinite(full_vector)):
        return np.full(18, -1e100)
    reflected_densities = full_vector[2:5]
    return np.concatenate(
        (
            reflected_densities,
            1.0 - reflected_densities,
            general.inequality_values(
                full_vector, stability_floor, order_parameter_cap
            ),
        )
    )


def baseline(stability_floor: float) -> np.ndarray:
    full = general.restricted_seed(stability_floor)
    return np.concatenate((full[:2], full[5:11], full[11:15]))


BOUNDS = [
    (-8.0, 8.0),
    (-8.0, 8.0),
    *[(0.0, 1.0)] * 6,
    *[(-math.pi, math.pi)] * 3,
    (0.0, 0.9),
]


def solve(
    start: np.ndarray,
    stability_floor: float,
    order_parameter_cap: float,
    max_iterations: int = 2500,
):
    return minimize(
        lambda vector: -vector[-1],
        start,
        method="SLSQP",
        bounds=BOUNDS,
        constraints={
            "type": "ineq",
            "fun": lambda vector: inequality_values(
                vector, stability_floor, order_parameter_cap
            ),
        },
        options={"ftol": 2e-12, "maxiter": max_iterations, "disp": False},
    )


def pack_result(
    result, stability_floor: float, order_parameter_cap: float
) -> dict[str, Any]:
    full_vector = lift(result.x)
    finite = bool(np.all(np.isfinite(full_vector)))
    if not finite:
        return {
            "solver_success": bool(result.success),
            "feasible": False,
            "message": str(result.message),
        }
    state = general.evaluate(full_vector)
    inequalities = inequality_values(
        result.x, stability_floor, order_parameter_cap
    )
    tolerance = max(1e-12, stability_floor * 1e-3)
    torque_error = float(np.max(np.abs(state["positive_class_torques"])))
    strictly_positive = bool(min(
        state["even_physical_eigenvalues"][0],
        state["odd_eigenvalues"][0],
        *state["transverse_values"],
    ) > 0.0)
    feasible = bool(
        np.min(inequalities) >= -tolerance
        and torque_error <= 1e-9
        and strictly_positive
    )
    return {
        "solver_success": bool(result.success),
        "feasible": feasible,
        "message": str(result.message),
        "iterations": int(result.nit),
        "stability_floor": float(stability_floor),
        "mu": state["objective"],
        "pair_masses": [float(value) for value in state["pair_masses"]],
        "densities": {
            name: float(value)
            for name, value in zip(general.DENSITY_NAMES, state["densities"])
        },
        "angles": [float(value) for value in state["angles"]],
        "degrees": [float(value) for value in state["degrees"]],
        "even_physical_eigenvalues": [
            float(value) for value in state["even_physical_eigenvalues"]
        ],
        "odd_eigenvalues": [float(value) for value in state["odd_eigenvalues"]],
        "transverse_values": [float(value) for value in state["transverse_values"]],
        "order_parameter": state["order_parameter"],
        "maximum_torque_error": torque_error,
        "minimum_constraint_slack": float(np.min(inequalities)),
        "feasibility_tolerance": tolerance,
        "strictly_positive_stability_spectra": strictly_positive,
    }


def starting_points(
    stability_floor: float,
    near_starts: int,
    random_starts: int,
    seed: int,
) -> Iterator[tuple[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    base = baseline(stability_floor)
    yield "baseline", base
    for index in range(near_starts):
        vector = base.copy()
        vector[:2] += rng.normal(0.0, 0.7, size=2)
        vector[2:8] = np.clip(
            vector[2:8] + rng.normal(0.0, 0.3, size=6), 0.0, 1.0
        )
        vector[8:11] = np.clip(
            vector[8:11] + rng.normal(0.0, 0.35, size=3),
            -math.pi,
            math.pi,
        )
        vector[-1] = rng.uniform(0.50, 0.69)
        yield f"near_{index}", vector
    made = 0
    while made < random_starts:
        angles = rng.uniform(-math.pi, math.pi, size=3)
        if np.min(np.abs(np.sin(2.0 * angles))) < 0.12:
            continue
        vector = np.concatenate(
            (
                rng.uniform(-2.5, 2.5, size=2),
                rng.uniform(0.0, 1.0, size=6),
                angles,
                [rng.uniform(0.35, 0.69)],
            )
        )
        yield f"random_{made}", vector
        made += 1


def search(args: argparse.Namespace) -> dict[str, Any]:
    first_floor = float(args.floors[0])
    attempts = []
    feasible_endpoints: list[tuple[np.ndarray, dict[str, Any]]] = []
    for name, start in starting_points(
        first_floor, args.near_starts, args.random_starts, args.seed
    ):
        result = solve(start, first_floor, args.order_parameter_cap)
        packed = pack_result(result, first_floor, args.order_parameter_cap)
        attempts.append(
            {
                "start": name,
                "solver_success": packed["solver_success"],
                "feasible": packed["feasible"],
                "mu": packed.get("mu"),
                "message": packed["message"],
            }
        )
        if packed["feasible"]:
            feasible_endpoints.append((result.x, packed))
    if not feasible_endpoints:
        raise RuntimeError("the broad search found no feasible endpoints")

    feasible_endpoints.sort(key=lambda item: item[1]["mu"], reverse=True)
    representatives: list[tuple[np.ndarray, dict[str, Any]]] = []
    for vector, packed in feasible_endpoints:
        if all(
            abs(packed["mu"] - prior[1]["mu"]) > args.branch_tolerance
            for prior in representatives
        ):
            representatives.append((vector, packed))
        if len(representatives) >= args.max_branches:
            break

    branches = []
    for branch_index, (vector, first_record) in enumerate(representatives):
        continuation = [first_record]
        current = vector
        for floor in args.floors[1:]:
            result = solve(current, float(floor), args.order_parameter_cap)
            record = pack_result(result, float(floor), args.order_parameter_cap)
            continuation.append(record)
            if not record["feasible"]:
                break
            current = result.x
        branches.append(
            {
                "branch": branch_index,
                "continuation": continuation,
            }
        )

    return {
        "scope": "exploratory local multistart search; no global claim",
        "parameterization": (
            "Six cross-pair densities are free; three reflection-pair "
            "densities are solved exactly from the torque equations."
        ),
        "seed": args.seed,
        "floors": [float(value) for value in args.floors],
        "attempt_summary": {
            "attempts": len(attempts),
            "solver_successes": sum(item["solver_success"] for item in attempts),
            "feasible_endpoints": sum(item["feasible"] for item in attempts),
            "feasible_objectives": sorted(
                float(item["mu"])
                for item in attempts
                if item["feasible"] and item["mu"] is not None
            ),
        },
        "attempts": attempts,
        "branches": branches,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--floors", nargs="+", type=float, default=DEFAULT_FLOORS)
    parser.add_argument("--near-starts", type=int, default=12)
    parser.add_argument("--random-starts", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--max-branches", type=int, default=5)
    parser.add_argument("--branch-tolerance", type=float, default=1e-4)
    parser.add_argument("--order-parameter-cap", type=float, default=0.95)
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
