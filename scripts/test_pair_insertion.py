#!/usr/bin/env python3
"""First-order test: can a fourth reflected pair improve the six-class endpoint?

The old marginal six-class point has three quotient zero modes (rotation plus
one even and one odd). Give a new reflected pair mass epsilon per class and
phase +/-delta. At fixed delta, the first-order torque, degree, and Hessian
conditions are affine in the old parameter derivatives and six new cross-block
densities. A four-by-four semidefinite constraint opens *all* nonrotation zero
modes to first order. This is a directional discovery test, not a global bound:
failure at sampled deltas cannot exclude other phases or higher-order paths.
"""

from __future__ import annotations

import argparse
import json
import math
import warnings
from pathlib import Path

import cvxpy as cp
import numpy as np

import optimize_general_reflection_model as general

ROOT = Path(__file__).resolve().parents[1]


def base_point():
    saved = json.loads((ROOT / "data/limiting_kkt_certificate.json").read_text())
    values = np.array(list(map(float, saved["variables"])))
    pair_mass = values[:3]
    masses = np.repeat(pair_mass, 2)
    x, y = values[3:5]
    pair_angles = values[5:8]
    phases = np.repeat(pair_angles, 2) * np.tile([1.0, -1.0], 3)
    densities = np.array([x, 1.0, 0.0, y, 0.0, 1.0, 1.0, 0.0, 1.0])
    weights = general.block_matrix(densities)
    return pair_mass, masses, pair_angles, phases, densities, weights, values[8]


def old_density_basis() -> list[np.ndarray]:
    identity = np.eye(6)
    basis = []
    for index in range(9):
        coordinate = np.zeros(9)
        coordinate[index] = 1.0
        basis.append(general.block_matrix(coordinate) - identity)
    return basis


BASIS = old_density_basis()


def edge_laplacian(size: int, i: int, j: int) -> np.ndarray:
    vector = np.zeros(size)
    vector[i], vector[j] = 1.0, -1.0
    return np.outer(vector, vector)


def old_kernel(masses: np.ndarray, phases: np.ndarray, weights: np.ndarray):
    laplacian = np.zeros((8, 8))
    for i in range(6):
        for j in range(i + 1, 6):
            coefficient = (
                masses[i]
                * masses[j]
                * weights[i, j]
                * math.cos(phases[i] - phases[j])
            )
            laplacian += coefficient * edge_laplacian(8, i, j)
    eigenvalues = np.linalg.eigvalsh(laplacian[:6, :6])
    if max(abs(eigenvalues[index]) for index in range(3)) > 1e-10:
        raise RuntimeError("the starting six-class point does not have three zero modes")
    if eigenvalues[3] < 1e-6:
        raise RuntimeError("the starting point has an unexpected fourth zero mode")
    even = np.array(
        [[laplacian[2 * i, 2 * j] + laplacian[2 * i, 2 * j + 1] for j in range(3)] for i in range(3)]
    )
    odd = np.array(
        [[laplacian[2 * i, 2 * j] - laplacian[2 * i, 2 * j + 1] for j in range(3)] for i in range(3)]
    )
    even_vector = np.cross(np.ones(3), even[0])
    odd_vector = np.cross(odd[0], odd[1])
    if np.linalg.norm(even_vector) < 1e-10 or np.linalg.norm(odd_vector) < 1e-10:
        raise RuntimeError("cofactor kernel construction is singular")
    pair_masses = masses[::2]
    even_vector -= np.dot(pair_masses, even_vector) / np.sum(pair_masses)
    even_vector /= np.linalg.norm(even_vector) * math.sqrt(2)
    odd_vector /= np.linalg.norm(odd_vector) * math.sqrt(2)
    weak = np.zeros((8, 2))
    weak[:6, 0] = np.repeat(even_vector, 2)
    weak[:6, 1] = np.repeat(odd_vector, 2) * np.tile([1.0, -1.0], 3)
    if np.linalg.norm(laplacian @ weak, ord=np.inf) > 1e-10:
        raise RuntimeError("cofactor vectors are not old weak modes")
    # These columns, together with old rotation, span the entire five-dimensional
    # zero eigenspace. L1 annihilates global rotation, so PSD on this (nonorthogonal)
    # four-column basis is equivalent to PSD on the full first-order kernel.
    quotient_basis = weak.copy()
    quotient_basis = np.column_stack((quotient_basis, np.eye(8)[:, 6:8]))
    return quotient_basis, weak, eigenvalues.tolist()


def solve_at_phase(
    delta: float,
    sigma: float,
    derivative_cap: float,
    weak_only: bool = False,
    new_degree_slack: float = 1e-5,
    dual_upper: bool = False,
) -> dict:
    pair_mass, masses, pair_angles, phases, densities, weights, mu_star = base_point()
    quotient_basis, old_weak, old_eigenvalues = old_kernel(masses, phases, weights)
    new_phases = (delta, -delta)

    mass_derivative = cp.Variable(3)
    angle_derivative = cp.Variable(3)
    density_derivative = cp.Variable(9)
    new_same = cp.Variable(3)
    new_opposite = cp.Variable(3)
    degree_gain = cp.Variable()

    dm = [mass_derivative[i // 2] for i in range(6)]
    dt = [angle_derivative[i // 2] * (1 if i % 2 == 0 else -1) for i in range(6)]

    def dw(i: int, j: int):
        return sum(density_derivative[k] * BASIS[k][i, j] for k in range(9))

    def new_weight(i: int, new_index: int):
        same_sign = (i % 2) == (new_index - 6)
        return new_same[i // 2] if same_sign else new_opposite[i // 2]

    old_torque_derivatives = []
    old_degree_derivatives = []
    for i in (0, 2, 4):
        torque_terms = []
        degree_terms = []
        for j in range(6):
            phase_difference = phases[j] - phases[i]
            torque_terms.append(
                dm[j] * weights[i, j] * math.sin(phase_difference)
                + masses[j] * dw(i, j) * math.sin(phase_difference)
                + masses[j]
                * weights[i, j]
                * math.cos(phase_difference)
                * (dt[j] - dt[i])
            )
            degree_terms.append(dm[j] * weights[i, j] + masses[j] * dw(i, j))
        for j in (6, 7):
            torque_terms.append(
                new_weight(i, j) * math.sin(new_phases[j - 6] - phases[i])
            )
            degree_terms.append(new_weight(i, j))
        old_torque_derivatives.append(sum(torque_terms))
        old_degree_derivatives.append(sum(degree_terms))

    new_torque = sum(
        masses[j]
        * new_weight(j, 6)
        * math.sin(phases[j] - delta)
        for j in range(6)
    )
    new_degree = sum(masses[j] * new_weight(j, 6) for j in range(6))
    new_transverse = sum(
        masses[j]
        * new_weight(j, 6)
        * math.cos(phases[j] - delta)
        for j in range(6)
    )

    first_hessian = cp.Constant(np.zeros((8, 8)))
    for i in range(6):
        for j in range(i + 1, 6):
            phase_difference = phases[i] - phases[j]
            cosine, sine = math.cos(phase_difference), math.sin(phase_difference)
            coefficient = (
                (dm[i] * masses[j] + masses[i] * dm[j])
                * weights[i, j]
                * cosine
                + masses[i] * masses[j] * dw(i, j) * cosine
                - masses[i]
                * masses[j]
                * weights[i, j]
                * sine
                * (dt[i] - dt[j])
            )
            first_hessian += coefficient * edge_laplacian(8, i, j)
    for i in range(6):
        for j in (6, 7):
            coefficient = (
                masses[i]
                * new_weight(i, j)
                * math.cos(phases[i] - new_phases[j - 6])
            )
            first_hessian += coefficient * edge_laplacian(8, i, j)
    weak_opening = quotient_basis.T @ first_hessian @ quotient_basis
    necessary_old_weak_opening = old_weak.T @ first_hessian @ old_weak

    mass_constraint = cp.sum(mass_derivative) == -1
    torque_constraints = [value == 0 for value in old_torque_derivatives]
    new_torque_constraint = new_torque == 0
    new_degree_constraint = new_degree >= mu_star + new_degree_slack
    new_transverse_constraint = new_transverse >= 0
    active_degree_constraints = [
        old_degree_derivatives[0] >= degree_gain,
        old_degree_derivatives[2] >= degree_gain,
    ]
    psd_constraint = (
        necessary_old_weak_opening >> sigma * np.eye(2)
        if weak_only
        else weak_opening >> sigma * np.eye(4)
    )
    constraints = [
        mass_constraint,
        mass_derivative <= derivative_cap,
        mass_derivative >= -derivative_cap,
        angle_derivative <= derivative_cap,
        angle_derivative >= -derivative_cap,
        density_derivative <= derivative_cap,
        density_derivative >= -derivative_cap,
        new_same >= 0,
        new_same <= 1,
        new_opposite >= 0,
        new_opposite <= 1,
        *torque_constraints,
        new_torque_constraint,
        new_degree_constraint,
        new_transverse_constraint,
        *active_degree_constraints,
        psd_constraint,
    ]
    for index, value in enumerate(densities):
        if value < 1e-12:
            constraints.append(density_derivative[index] >= 0)
        elif value > 1 - 1e-12:
            constraints.append(density_derivative[index] <= 0)

    problem = cp.Problem(cp.Maximize(degree_gain), constraints)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            problem.solve(
                solver="CLARABEL",
                tol_gap_abs=1e-9,
                tol_gap_rel=1e-9,
                tol_feas=1e-9,
                max_iter=500,
            )
        except cp.error.SolverError:
            return {"delta": delta, "status": "solver_error"}
    record = {
        "delta": delta,
        "status": str(problem.status),
        "old_positive_eigenvalues": old_eigenvalues[3:],
    }
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        return record
    record.update(
        {
            "degree_gain_derivative": float(degree_gain.value),
            "new_pair_degree_at_zero_mass": float(new_degree.value),
            "new_pair_transverse_at_zero_mass": float(new_transverse.value),
            "weak_opening_eigenvalues": np.linalg.eigvalsh(weak_opening.value).tolist(),
            "old_weak_opening_eigenvalues": np.linalg.eigvalsh(
                necessary_old_weak_opening.value
            ).tolist(),
            "old_degree_derivatives": [float(v.value) for v in old_degree_derivatives],
            "old_torque_derivatives": [float(v.value) for v in old_torque_derivatives],
            "new_torque": float(new_torque.value),
            "mass_derivative": mass_derivative.value.tolist(),
            "angle_derivative": angle_derivative.value.tolist(),
            "density_derivative": density_derivative.value.tolist(),
            "new_same": new_same.value.tolist(),
            "new_opposite": new_opposite.value.tolist(),
            "first_hessian": first_hessian.value.tolist(),
        }
    )
    if dual_upper:
        # A Lagrangian upper bound is independently computable for any
        # nonnegative scalar multipliers and PSD matrix, even when the chosen
        # multipliers are only approximate. We retain explicit variable boxes
        # rather than relying on exact numerical stationarity.
        active_dual = np.array([c.dual_value for c in active_degree_constraints], dtype=float)
        active_dual = np.maximum(active_dual, 0.0)
        active_dual /= np.sum(active_dual)
        psd_raw = np.asarray(psd_constraint.dual_value, dtype=float)
        psd_eig, psd_vec = np.linalg.eigh((psd_raw + psd_raw.T) / 2)
        psd_factor = psd_vec * np.sqrt(np.maximum(psd_eig, 0.0))
        psd_dual = psd_factor @ psd_factor.T
        # CVXPY reports equality multipliers for the minimization of -gain;
        # the maximizing Lagrangian below uses their negatives.
        mass_dual = -float(mass_constraint.dual_value)
        torque_dual = [-float(c.dual_value) for c in torque_constraints]
        new_torque_dual = -float(new_torque_constraint.dual_value)
        lagrangian = (
            degree_gain
            + mass_dual * (cp.sum(mass_derivative) + 1)
            + sum(torque_dual[i] * old_torque_derivatives[i] for i in range(3))
            + new_torque_dual * new_torque
            + float(new_degree_constraint.dual_value) * (new_degree - mu_star - new_degree_slack)
            + float(new_transverse_constraint.dual_value) * new_transverse
            + sum(
                active_dual[k] * (old_degree_derivatives[i] - degree_gain)
                for k, i in enumerate((0, 2))
            )
            + cp.trace(psd_dual @ (
                necessary_old_weak_opening - sigma * np.eye(2)
                if weak_only else weak_opening - sigma * np.eye(4)
            ))
        )
        variables = [mass_derivative, angle_derivative, density_derivative, new_same, new_opposite]
        for variable in variables:
            variable.value = np.zeros(variable.size)
        degree_gain.value = 0.0
        constant = float(lagrangian.value)
        support = constant
        coefficients = []
        for variable in variables:
            for index in range(variable.size):
                vector = np.zeros(variable.size)
                vector[index] = 1.0
                variable.value = vector
                coefficient = float(lagrangian.value) - constant
                variable.value = np.zeros(variable.size)
                if variable is mass_derivative or variable is angle_derivative:
                    lo, hi = -derivative_cap, derivative_cap
                elif variable is density_derivative:
                    value = densities[index]
                    lo, hi = ((0.0, derivative_cap) if value < 1e-12 else
                              ((-derivative_cap, 0.0) if value > 1 - 1e-12 else
                               (-derivative_cap, derivative_cap)))
                else:
                    lo, hi = 0.0, 1.0
                support += max(coefficient * lo, coefficient * hi)
                coefficients.append(coefficient)
        record["dual_upper_float"] = float(support)
        record["dual_coefficient_max_abs"] = float(max(abs(x) for x in coefficients))
        record["dual_multipliers"] = {
            "mass": mass_dual,
            "old_torque": torque_dual,
            "new_torque": new_torque_dual,
            "new_degree": float(new_degree_constraint.dual_value),
            "new_transverse": float(new_transverse_constraint.dual_value),
            "active_degree": active_dual.tolist(),
            "psd": psd_dual.tolist(),
            "psd_factor": psd_factor.tolist(),
        }
    return record


def audit_finite_difference(record: dict, step: float = 1e-7) -> dict:
    """Independently assemble the eight-class system and check the linearization."""
    pair_mass, old_masses, pair_angles, _, densities, _, _ = base_point()
    delta = record["delta"]
    masses0 = np.r_[old_masses, 0.0, 0.0]
    def evaluate(epsilon: float):
        pair_masses = pair_mass + epsilon * np.array(record["mass_derivative"])
        masses = np.r_[np.repeat(pair_masses, 2), epsilon, epsilon]
        pair_phases = pair_angles + epsilon * np.array(record["angle_derivative"])
        phases = np.r_[np.repeat(pair_phases, 2) * np.tile([1.0, -1.0], 3), delta, -delta]
        weights = np.eye(8)
        weights[:6, :6] = general.block_matrix(
            densities + epsilon * np.array(record["density_derivative"])
        )
        for old in range(6):
            for new in (6, 7):
                coefficient = (
                    record["new_same"][old // 2]
                    if old % 2 == new - 6
                    else record["new_opposite"][old // 2]
                )
                weights[old, new] = weights[new, old] = coefficient
        degrees = weights @ masses
        torques = np.array(
            [np.sum(weights[i] * masses * np.sin(phases - phases[i])) for i in range(8)]
        )
        laplacian = np.zeros((8, 8))
        for i in range(8):
            for j in range(i + 1, 8):
                coefficient = (
                    masses[i]
                    * masses[j]
                    * weights[i, j]
                    * math.cos(phases[i] - phases[j])
                )
                laplacian += coefficient * edge_laplacian(8, i, j)
        return degrees, torques, laplacian

    degree0, torque0, laplacian0 = evaluate(0.0)
    degree1, torque1, laplacian1 = evaluate(step)
    fd_degree = (degree1[:6:2] - degree0[:6:2]) / step
    fd_torque = (torque1[:6:2] - torque0[:6:2]) / step
    fd_laplacian = (laplacian1 - laplacian0) / step
    spectrum = np.linalg.eigvalsh(laplacian1)
    return {
        "step": step,
        "degree_error": float(np.max(np.abs(fd_degree - record["old_degree_derivatives"]))),
        "old_torque_error": float(np.max(np.abs(fd_torque - record["old_torque_derivatives"]))),
        "new_torque_error": float(abs(torque0[6] - record["new_torque"])),
        "new_degree_error": float(abs(degree0[6] - record["new_pair_degree_at_zero_mass"])),
        "hessian_error": float(np.max(np.abs(fd_laplacian - record["first_hessian"]))),
        "mass_sum_error": float(abs(np.sum(masses0) - 1.0)),
        "rotation_error": float(np.linalg.norm(laplacian1 @ np.ones(8))),
        "lowest_nonrotation_eigenvalue": float(spectrum[0] if spectrum[0] < -1e-12 else spectrum[1]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=31)
    parser.add_argument("--sigma", type=float, default=1e-7)
    parser.add_argument("--derivative-cap", type=float, default=10.0)
    parser.add_argument("--weak-only", action="store_true")
    parser.add_argument("--new-degree-slack", type=float, default=1e-5)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if (
        args.grid < 2
        or args.sigma < 0
        or args.derivative_cap <= 0
        or args.new_degree_slack < 0
    ):
        parser.error("grid>=2, sigma>=0, derivative-cap>0, and slack>=0 are required")
    phases = np.linspace(0.0, math.pi, args.grid)
    records = [
        solve_at_phase(
            float(delta),
            args.sigma,
            args.derivative_cap,
            args.weak_only,
            args.new_degree_slack,
        )
        for delta in phases
    ]
    feasible = [record for record in records if "degree_gain_derivative" in record]
    best = max(feasible, key=lambda record: record["degree_gain_derivative"]) if feasible else None
    report = {
        "scope": (
            "First-order four-pair insertion test at a marginal six-class endpoint. "
            "A positive degree derivative and positive weak-opening eigenvalues "
            "give a candidate improving branch; negative results on a grid are "
            "not a global exclusion."
        ),
        "grid": args.grid,
        "sigma": args.sigma,
        "derivative_cap": args.derivative_cap,
        "weak_only": args.weak_only,
        "new_degree_slack": args.new_degree_slack,
        "best": best,
        "records": records,
    }
    if args.audit and best:
        report["finite_difference_audit"] = audit_finite_difference(best)
    rendered = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered)
        print(json.dumps({"best": best["degree_gain_derivative"], "audit": report.get("finite_difference_audit"), "out": str(args.out)}))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
