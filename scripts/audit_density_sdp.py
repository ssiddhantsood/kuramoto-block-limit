#!/usr/bin/env python3
"""Audit the affine density decomposition used by ``solve_density_sdp.py``.

The test samples general three-pair reflected models and compares the affine
degree, torque, even-Hessian, odd-Hessian, and transverse reconstructions with
the direct six-block evaluator.  It tests the algebra behind the SDP
formulation; it does not turn the floating-point SDP or outer search into a
global certificate.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np

import optimize_general_reflection_model as general
import solve_density_sdp as density_sdp


def direct_sectors(state: dict) -> tuple[np.ndarray, np.ndarray]:
    hessian = state["full_block_hessian"]
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
    basis = general.restricted.physical_even_basis(state["pair_masses"])
    return basis.T @ even @ basis, odd


def run(samples: int, seed: int, tolerance: float) -> dict:
    rng = np.random.default_rng(seed)
    maximum_errors = {
        "degree": 0.0,
        "torque": 0.0,
        "even_matrix": 0.0,
        "odd_matrix": 0.0,
        "transverse": 0.0,
    }
    for _ in range(samples):
        logits = rng.uniform(-2.5, 2.5, size=2)
        pair_masses = general.restricted.pair_masses(*logits)
        angles = rng.uniform(-math.pi, math.pi, size=3)
        densities = rng.uniform(0.0, 1.0, size=9)
        vector = np.concatenate((logits, densities, angles, [0.0]))
        state = general.evaluate(vector)
        coefficients = density_sdp.sector_coefficients(pair_masses, angles)

        affine_degree = pair_masses + coefficients["degree"] @ densities
        affine_torque = coefficients["torque"] @ densities
        affine_even = sum(
            densities[index] * coefficients["even"][index]
            for index in range(9)
        )
        affine_odd = sum(
            densities[index] * coefficients["odd"][index]
            for index in range(9)
        )
        affine_transverse = (
            pair_masses + coefficients["transverse"] @ densities
        )
        direct_even, direct_odd = direct_sectors(state)
        comparisons = {
            "degree": (affine_degree, state["degrees"]),
            "torque": (affine_torque, state["positive_class_torques"]),
            "even_matrix": (affine_even, direct_even),
            "odd_matrix": (affine_odd, direct_odd),
            "transverse": (affine_transverse, state["transverse_values"]),
        }
        for name, (affine, direct) in comparisons.items():
            error = float(np.max(np.abs(affine - direct)))
            maximum_errors[name] = max(maximum_errors[name], error)

    passed = max(maximum_errors.values()) <= tolerance
    return {
        "scope": (
            "Randomized algebra audit of the affine density decomposition; "
            "not an optimization or global-optimality certificate."
        ),
        "samples": samples,
        "seed": seed,
        "tolerance": tolerance,
        "maximum_absolute_errors": maximum_errors,
        "passed": passed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    payload = run(args.samples, args.seed, args.tolerance)
    print(json.dumps(payload, indent=2))
    if not payload["passed"]:
        raise SystemExit("affine density audit failed")


if __name__ == "__main__":
    main()
