#!/usr/bin/env python3
"""Optimize the observed three-pair branch using only its three phase angles.

On this branch the three torque equations and the active-degree identity
``d_A=d_C`` determine the two mass ratios and both partial block densities.
The resulting optimization is local to this topology and regular phase chart;
it is not a global bound over other three-pair density faces.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize

import optimize_reflection_model as restricted

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/reduced_optimization_results.json"
DEFAULT_FLOORS = (1e-5, 1e-7, 1e-9)


class SingularPhaseChart(ValueError):
    """The torque/degree elimination is singular or leaves positive masses."""


def reconstruct(angles: np.ndarray) -> dict[str, Any]:
    """Recover masses and densities from three signed phases.

    Here ``k=b/a`` and ``h=c/a``.  The formula is exact wherever its three
    trigonometric and one linear denominators are nonzero.  It does not impose
    the density box or stability; callers check those inequalities separately.
    """
    alpha, beta, gamma = map(float, angles)
    sin_ab = math.sin(alpha - beta)
    sin_bg = math.sin(beta + gamma)
    sin_2a = math.sin(2.0 * alpha)
    for name, value in (
        ("sin(alpha-beta)", sin_ab),
        ("sin(beta+gamma)", sin_bg),
        ("sin(2*alpha)", sin_2a),
    ):
        if abs(value) < 1e-9:
            raise SingularPhaseChart(f"{name} is too small")

    k = -2.0 * math.cos(alpha) * math.sin(gamma) / sin_bg
    y_zero = k * math.sin(2.0 * beta) / sin_ab
    y_slope = sin_bg / sin_ab
    x_zero = -(k * k * math.sin(2.0 * beta)) / sin_2a
    x_slope = -(k * sin_bg + 2.0 * math.cos(gamma) * math.sin(alpha)) / sin_2a

    # d_A-d_C = a * (x+k*y+h-1-k).  Since x and y are affine in h, this
    # equality determines h without a nonlinear solve.
    degree_intercept = x_zero + k * y_zero - 1.0 - k
    degree_slope = x_slope + k * y_slope + 1.0
    if abs(degree_slope) < 1e-9:
        raise SingularPhaseChart("active-degree equation is singular")
    h = -degree_intercept / degree_slope
    if not (math.isfinite(k) and math.isfinite(h) and k > 0.0 and h > 0.0):
        raise SingularPhaseChart("reconstructed pair masses are not positive")

    a = 1.0 / (2.0 * (1.0 + k + h))
    b, c = k * a, h * a
    x = x_zero + h * x_slope
    y = y_zero + h * y_slope
    if not all(math.isfinite(value) for value in (a, b, c, x, y)):
        raise SingularPhaseChart("nonfinite reconstructed parameter")

    # The separate six-class evaluator supplies independent torque, degree,
    # and physical Hessian checks; the objective here is d_C=1/2+a.
    mu = 0.5 + a
    vector = np.array([math.log(a / c), math.log(b / c), x, y, alpha, beta, gamma, mu])
    direct = restricted.evaluate(vector)
    return {
        "pair_masses": np.array([a, b, c]),
        "partial_densities": np.array([x, y]),
        "angles": np.array([alpha, beta, gamma]),
        "mass_ratios": np.array([k, h]),
        "mu": mu,
        "direct": direct,
    }


def constraint_values(state: dict[str, Any], floor: float) -> np.ndarray:
    direct = state["direct"]
    x, y = state["partial_densities"]
    degrees = direct["degrees"]
    return np.array(
        [
            x,
            1.0 - x,
            y,
            1.0 - y,
            degrees[1] - state["mu"],
            direct["even_physical_eigenvalues"][0] - floor,
            direct["odd_eigenvalues"][0] - floor,
            *list(direct["transverse_values"] - floor),
            0.95 - direct["order_parameter"],
        ]
    )


def solve_floor(
    floor: float, start: np.ndarray, phase_radius: float = 0.2
) -> dict[str, Any]:
    """Run a local three-angle SLSQP optimization near the known branch."""
    if floor <= 0.0:
        raise ValueError("the stability floor must be positive")
    center = np.array([-0.729628, -2.358745, 0.956552], dtype=float)
    bounds = [(value - phase_radius, value + phase_radius) for value in center]
    calls = 0

    def safe_state(angles: np.ndarray) -> dict[str, Any] | None:
        nonlocal calls
        calls += 1
        try:
            return reconstruct(angles)
        except (SingularPhaseChart, OverflowError, ValueError):
            return None

    def objective(angles: np.ndarray) -> float:
        state = safe_state(angles)
        return -state["mu"] if state is not None else 1.0

    def inequalities(angles: np.ndarray) -> np.ndarray:
        state = safe_state(angles)
        return (
            constraint_values(state, floor) if state is not None else np.full(11, -1.0)
        )

    begun = time.perf_counter()
    result = minimize(
        objective,
        np.asarray(start, dtype=float),
        method="SLSQP",
        bounds=bounds,
        constraints=[{"type": "ineq", "fun": inequalities}],
        options={"ftol": 1e-13, "maxiter": 500, "disp": False},
    )
    elapsed = time.perf_counter() - begun
    state = reconstruct(result.x)
    direct = state["direct"]
    minimum_slack = float(np.min(constraint_values(state, floor)))
    torque_error = float(np.max(np.abs(direct["positive_class_torques"])))
    active_degree_error = float(abs(direct["degrees"][0] - direct["degrees"][2]))
    if (
        not result.success
        or minimum_slack < -2e-11
        or torque_error > 2e-12
        or active_degree_error > 2e-12
        or min(
            direct["even_physical_eigenvalues"][0],
            direct["odd_eigenvalues"][0],
            *direct["transverse_values"],
        )
        <= 0.0
    ):
        raise RuntimeError(
            f"invalid phase-only optimizer result: {result.message}; "
            f"slack={minimum_slack}, torque={torque_error}, "
            f"degree tie={active_degree_error}"
        )
    return {
        "stability_floor": floor,
        "mu": state["mu"],
        "pair_masses": state["pair_masses"].tolist(),
        "mass_ratios_b_over_a_c_over_a": state["mass_ratios"].tolist(),
        "partial_block_densities_x_y": state["partial_densities"].tolist(),
        "angles": state["angles"].tolist(),
        "degrees": direct["degrees"].tolist(),
        "even_physical_eigenvalues": direct["even_physical_eigenvalues"].tolist(),
        "odd_eigenvalues": direct["odd_eigenvalues"].tolist(),
        "transverse_values": direct["transverse_values"].tolist(),
        "maximum_torque_error": torque_error,
        "active_degree_error": active_degree_error,
        "minimum_constraint_slack": minimum_slack,
        "solver": {
            "success": bool(result.success),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_evaluations": int(result.nfev),
            "reconstruction_calls": calls,
            "elapsed_seconds": elapsed,
        },
    }


def run(floors: list[float], phase_radius: float) -> dict[str, Any]:
    reference = json.loads(REFERENCE.read_text())
    start = np.array(reference["results"][0]["angles"], dtype=float)
    records = []
    for floor in floors:
        record = solve_floor(floor, start, phase_radius)
        match = next(
            (
                item
                for item in reference["results"]
                if float(item["stability_floor"]) == floor
            ),
            None,
        )
        if match is not None:
            record["difference_from_previous_eight_variable_mu"] = record["mu"] - float(
                match["mu"]
            )
        records.append(record)
        start = np.array(record["angles"], dtype=float)
    return {
        "scope": (
            "Local optimization on the observed two-density, active-degree "
            "three-pair branch. Other density faces and singular phase charts "
            "are not covered. Numerical SLSQP values are not certificates."
        ),
        "derivation": (
            "T_C=0 gives k=b/a; T_B=0 gives y affine in h=c/a; "
            "T_A=0 gives x affine in h; d_A=d_C gives h. "
            "Thus only three signed phases are optimized and mu=d_C=1/2+a."
        ),
        "nonlinear_variables": 3,
        "eliminated_variables": ["b/a", "c/a", "x", "y"],
        "positive_stability_floors": floors,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--floors", nargs="+", type=float, default=DEFAULT_FLOORS)
    parser.add_argument("--phase-radius", type=float, default=0.2)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if not args.floors or any(floor <= 0.0 for floor in args.floors):
        parser.error("provide one or more positive stability floors")
    report = run(args.floors, args.phase_radius)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
