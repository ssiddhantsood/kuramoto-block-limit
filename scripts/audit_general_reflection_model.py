#!/usr/bin/env python3
"""Independent numerical consistency checks for the general reflected model."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import optimize_general_reflection_model as general


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data/general_reflection_optimization.json"


def reconstruct_vector(record: dict[str, Any]) -> np.ndarray:
    a, b, c = map(float, record["pair_masses"])
    densities = [float(record["densities"][name]) for name in general.DENSITY_NAMES]
    return np.array(
        [
            math.log(a / c),
            math.log(b / c),
            *densities,
            *map(float, record["angles"]),
            float(record["mu"]),
        ]
    )


def numerical_torque_jacobian(
    masses: np.ndarray, weights: np.ndarray, phases: np.ndarray, step: float
) -> np.ndarray:
    def torque(test_phases: np.ndarray) -> np.ndarray:
        return np.array(
            [
                np.sum(
                    masses
                    * weights[i]
                    * np.sin(test_phases - test_phases[i])
                )
                for i in range(6)
            ]
        )

    jacobian = np.zeros((6, 6))
    for column in range(6):
        offset = np.zeros(6)
        offset[column] = step
        jacobian[:, column] = (
            torque(phases + offset) - torque(phases - offset)
        ) / (2.0 * step)
    return jacobian


def audit_record(record: dict[str, Any], finite_difference_step: float) -> dict[str, Any]:
    vector = reconstruct_vector(record)
    state = general.evaluate(vector)
    masses = state["masses"]
    weights = state["weights"]
    phases = state["phases"]
    hessian = state["full_block_hessian"]

    reflection = np.zeros((6, 6))
    for pair in range(3):
        reflection[2 * pair, 2 * pair + 1] = 1.0
        reflection[2 * pair + 1, 2 * pair] = 1.0

    numerical_jacobian = numerical_torque_jacobian(
        masses, weights, phases, finite_difference_step
    )
    root_mass = np.diag(np.sqrt(masses))
    inverse_root_mass = np.diag(1.0 / np.sqrt(masses))
    finite_difference_hessian = (
        -root_mass @ numerical_jacobian @ inverse_root_mass
    )

    expected_sector_spectrum = np.sort(
        np.concatenate(
            (
                [0.0],
                state["even_physical_eigenvalues"],
                state["odd_eigenvalues"],
            )
        )
    )
    full_spectrum = np.sort(state["full_block_eigenvalues"])
    density_values = state["densities"]
    boundary_tolerance = 1e-7
    interior = {
        name: float(value)
        for name, value in zip(general.DENSITY_NAMES, density_values)
        if boundary_tolerance < value < 1.0 - boundary_tolerance
    }

    errors = {
        "weight_symmetry": float(np.max(np.abs(weights - weights.T))),
        "reflection_invariance": float(
            np.max(np.abs(weights - reflection @ weights @ reflection))
        ),
        "reflected_phases": float(
            np.max(np.abs(phases + reflection @ phases))
        ),
        "paired_degrees": float(
            np.max(np.abs(state["all_degrees"][0::2] - state["all_degrees"][1::2]))
        ),
        "opposite_torques": float(
            np.max(np.abs(state["all_torques"] + reflection @ state["all_torques"]))
        ),
        "paired_transverse_values": float(
            np.max(
                np.abs(
                    state["all_transverse_values"][0::2]
                    - state["all_transverse_values"][1::2]
                )
            )
        ),
        "rotation_kernel": float(
            np.max(np.abs(hessian @ np.sqrt(masses)))
        ),
        "sector_spectrum": float(
            np.max(np.abs(full_spectrum - expected_sector_spectrum))
        ),
        "finite_difference_hessian": float(
            np.max(np.abs(hessian - finite_difference_hessian))
        ),
        "positive_class_torque": float(
            np.max(np.abs(state["positive_class_torques"]))
        ),
    }
    thresholds = {
        **{name: 1e-10 for name in errors},
        "finite_difference_hessian": 2e-9,
    }
    failed = [name for name, value in errors.items() if value > thresholds[name]]
    floor = float(record["stability_floor"])
    feasibility_slacks = {
        "degree": float(np.min(state["degrees"] - state["objective"])),
        "even": float(state["even_physical_eigenvalues"][0] - floor),
        "odd": float(state["odd_eigenvalues"][0] - floor),
        "transverse": float(np.min(state["transverse_values"]) - floor),
        "nonsynchrony": float(0.95 - state["order_parameter"]),
    }
    feasibility_tolerance = max(1e-12, floor * 1e-3)
    strictly_positive_spectra = bool(min(
        state["even_physical_eigenvalues"][0],
        state["odd_eigenvalues"][0],
        *state["transverse_values"],
    ) > 0.0)
    if (
        min(feasibility_slacks.values()) < -feasibility_tolerance
        or not strictly_positive_spectra
    ):
        failed.append("feasibility")
    expected_zero = {"ref_A", "opposite_AB", "same_BC"}
    expected_one = {"ref_B", "same_AB", "same_AC", "opposite_AC"}
    observed_zero = {
        name
        for name, value in zip(general.DENSITY_NAMES, density_values)
        if value <= boundary_tolerance
    }
    observed_one = {
        name
        for name, value in zip(general.DENSITY_NAMES, density_values)
        if value >= 1.0 - boundary_tolerance
    }
    expected_interior = {"ref_C", "opposite_BC"}
    if (
        set(interior) != expected_interior
        or observed_zero != expected_zero
        or observed_one != expected_one
    ):
        failed.append("discovered_active_density_face")

    return {
        "stability_floor": floor,
        "mu": state["objective"],
        "passed": not failed,
        "failed_checks": failed,
        "maximum_identity_errors": errors,
        "feasibility_slacks": feasibility_slacks,
        "feasibility_tolerance": feasibility_tolerance,
        "strictly_positive_stability_spectra": strictly_positive_spectra,
        "interior_densities": interior,
        "boundary_density_counts": {
            "zero": int(np.sum(density_values <= boundary_tolerance)),
            "one": int(np.sum(density_values >= 1.0 - boundary_tolerance)),
            "interior": len(interior),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--finite-difference-step", type=float, default=1e-6)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = json.loads(args.input.read_text())
    audits = [
        audit_record(record, args.finite_difference_step)
        for record in source["results"]
    ]
    try:
        source_label = str(args.input.resolve().relative_to(ROOT))
    except ValueError:
        source_label = str(args.input)
    payload = {
        "source": source_label,
        "finite_difference_step": args.finite_difference_step,
        "passed": all(item["passed"] for item in audits),
        "audits": audits,
    }
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")
    if not payload["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
