#!/usr/bin/env python3
"""Interval-enclose the limiting KKT point on the discovered density face.

The general nine-density search puts seven densities on box boundaries.  On
that face only ``x`` and ``y`` remain free, giving the topology used by
``optimize_reflection_model.py``.  At the limiting point two degree constraints
and the weak even/odd stability constraints are active.  This script:

1. solves the corresponding KKT system at high precision;
2. proves existence and uniqueness of a KKT root in a small box with a
   Krawczyk interval inclusion; and
3. reports multiplier signs, constraint rank, and reduced curvature.

This is a local certificate on the discovered active face.  It is not a global
upper bound over all reflection-symmetric graphons.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import root


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data/reduced_optimization_results.json"
POINT_DPS = 100
INTERVAL_DPS = 100


def symbolic_system():
    a, b, c, x, y, alpha, beta, gamma, mu = sp.symbols(
        "a b c x y alpha beta gamma mu", real=True
    )
    variables = (a, b, c, x, y, alpha, beta, gamma, mu)
    multipliers = sp.symbols("lambda0:8", real=True)
    masses = (a, a, b, b, c, c)
    phases = (alpha, -alpha, beta, -beta, gamma, -gamma)
    weights = sp.Matrix(
        [
            [1, x, y, 0, 1, 1],
            [x, 1, 0, y, 1, 1],
            [y, 0, 1, 1, 0, 1],
            [0, y, 1, 1, 1, 0],
            [1, 1, 0, 1, 1, 0],
            [1, 1, 1, 0, 0, 1],
        ]
    )

    torques = []
    degrees = []
    laplacian = sp.zeros(6, 6)
    for i in range(6):
        degrees.append(sum(masses[j] * weights[i, j] for j in range(6)))
        torques.append(
            sum(
                masses[j]
                * weights[i, j]
                * sp.sin(phases[j] - phases[i])
                for j in range(6)
            )
        )
        for j in range(6):
            if i == j:
                continue
            cosine_weight = masses[j] * weights[i, j] * sp.cos(
                phases[i] - phases[j]
            )
            laplacian[i, i] += cosine_weight
            laplacian[i, j] = -cosine_weight

    even = sp.Matrix(
        3,
        3,
        lambda i, j: laplacian[2 * i, 2 * j]
        + laplacian[2 * i, 2 * j + 1],
    )
    odd = sp.Matrix(
        3,
        3,
        lambda i, j: laplacian[2 * i, 2 * j]
        - laplacian[2 * i, 2 * j + 1],
    )
    even_second_coefficient = sum(
        even.extract(indices, indices).det()
        for indices in ((0, 1), (0, 2), (1, 2))
    )

    constraints = sp.Matrix(
        [
            a + b + c - sp.Rational(1, 2),
            torques[0],
            torques[2],
            torques[4],
            degrees[0] - mu,
            degrees[4] - mu,
            even_second_coefficient,
            odd.det(),
        ]
    )
    objective = -mu
    lagrangian = objective + sum(
        multipliers[index] * constraints[index] for index in range(8)
    )
    stationarity = sp.Matrix(
        [sp.diff(lagrangian, variable) for variable in variables]
    )
    kkt = stationarity.col_join(constraints)
    symbols = variables + multipliers
    jacobian = kkt.jacobian(symbols)
    constraint_jacobian = constraints.jacobian(variables)
    lagrangian_hessian = sp.hessian(lagrangian, variables)

    transverse = [
        sum(
            masses[j]
            * weights[i, j]
            * sp.cos(phases[i] - phases[j])
            for j in range(6)
        )
        for i in range(6)
    ]

    inactive = sp.Matrix(
        [
            degrees[2] - mu,
            # The other even eigenvalue is positive iff the trace is positive
            # when the second characteristic coefficient vanishes.
            sp.trace(even),
            # At det(odd)=0, these are the sum and product of the two nonzero
            # eigenvalues.
            sp.trace(odd),
            sum(
                odd.extract(indices, indices).det()
                for indices in ((0, 1), (0, 2), (1, 2))
            ),
            a,
            b,
            c,
            x,
            1 - x,
            y,
            1 - y,
            transverse[0],
            transverse[2],
            transverse[4],
            sp.Rational(19, 20) ** 2
            - (
                2
                * (
                    a * sp.cos(alpha)
                    + b * sp.cos(beta)
                    + c * sp.cos(gamma)
                )
            )
            ** 2,
            mu,
            sp.Rational(4, 5) - mu,
            sp.pi + alpha,
            sp.pi - alpha,
            sp.pi + beta,
            sp.pi - beta,
            sp.pi + gamma,
            sp.pi - gamma,
        ]
    )
    return {
        "variables": variables,
        "multipliers": multipliers,
        "symbols": symbols,
        "constraints": constraints,
        "kkt": kkt,
        "jacobian": jacobian,
        "constraint_jacobian": constraint_jacobian,
        "lagrangian_hessian": lagrangian_hessian,
        "inactive": inactive,
    }


def initial_variables(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text())
    final = payload["results"][-1]
    a, b, c = map(float, final["pair_masses"])
    x, y = map(float, final["partial_block_densities"])
    alpha, beta, gamma = map(float, final["angles"])
    mu = float(final["mu"])
    return np.array([a, b, c, x, y, alpha, beta, gamma, mu])


def point_functions(system):
    symbols = system["symbols"]
    return {
        name: sp.lambdify(symbols, system[name], modules="numpy", cse=True)
        for name in ("kkt", "jacobian")
    }


def solve_double(system, seed_variables: np.ndarray) -> np.ndarray:
    variables = system["variables"]
    constraint_jacobian = sp.lambdify(
        variables, system["constraint_jacobian"], modules="numpy", cse=True
    )
    jacobian_at_seed = np.asarray(
        constraint_jacobian(*seed_variables), dtype=float
    )
    objective_gradient = np.zeros(9)
    objective_gradient[-1] = -1.0
    seed_multipliers = np.linalg.lstsq(
        jacobian_at_seed.T, -objective_gradient, rcond=None
    )[0]
    seed = np.concatenate((seed_variables, seed_multipliers))
    functions = point_functions(system)

    result = root(
        lambda value: np.asarray(functions["kkt"](*value), dtype=float).reshape(-1),
        seed,
        jac=lambda value: np.asarray(functions["jacobian"](*value), dtype=float),
        method="lm",
        options={"ftol": 1e-13, "xtol": 1e-13, "gtol": 1e-13, "maxiter": 3000},
    )
    residual = float(
        np.max(
            np.abs(
                np.asarray(functions["kkt"](*result.x), dtype=float).reshape(-1)
            )
        )
    )
    if not result.success or residual > 1e-8:
        raise RuntimeError(f"double KKT solve failed: {result.message}; residual={residual}")
    return result.x


def high_precision_root(system, seed: Sequence[float]) -> list[mp.mpf]:
    mp.mp.dps = POINT_DPS
    functions = [
        sp.lambdify(system["symbols"], expression, modules="mpmath", cse=True)
        for expression in system["kkt"]
    ]
    root_values = mp.findroot(
        tuple(functions),
        tuple(mp.mpf(str(value)) for value in seed),
        solver="mdnewton",
        tol=mp.mpf("1e-85"),
        maxsteps=100,
    )
    return [mp.mpf(value) for value in root_values]


def interval_bounds(value) -> tuple[mp.mpf, mp.mpf]:
    if not hasattr(value, "_mpi_"):
        value = mp.iv.mpf(value)
    return mp.mpf(value._mpi_[0]), mp.mpf(value._mpi_[1])


def interval_lambdify(symbols, expression):
    if isinstance(expression, sp.MatrixBase):
        expression = expression.tolist()
    return sp.lambdify(
        symbols,
        expression,
        modules=[
            {
                "sin": mp.iv.sin,
                "cos": mp.iv.cos,
                "pi": mp.iv.pi,
                "mpf": mp.iv.mpf,
            },
            "mpmath",
        ],
        cse=True,
    )


def krawczyk_certificate(system, center: Sequence[mp.mpf], radius: mp.mpf):
    mp.iv.dps = INTERVAL_DPS
    point_jacobian = sp.lambdify(
        system["symbols"], system["jacobian"], modules="mpmath", cse=True
    )
    interval_kkt_function = interval_lambdify(
        system["symbols"], system["kkt"]
    )
    interval_jacobian_function = interval_lambdify(
        system["symbols"], system["jacobian"]
    )

    center_text = [mp.nstr(value, 80) for value in center]
    center_point = [mp.mpf(value) for value in center_text]
    center_interval = [mp.iv.mpf(value) for value in center_text]
    jacobian_center = mp.matrix(point_jacobian(*center_point))
    inverse = jacobian_center**-1
    preconditioner = [
        [mp.iv.mpf(mp.nstr(inverse[i, j], 80)) for j in range(len(center))]
        for i in range(len(center))
    ]
    box = [
        mp.iv.mpf(
            [
                mp.nstr(value - radius, 90),
                mp.nstr(value + radius, 90),
            ]
        )
        for value in center_point
    ]
    f_center = [row[0] for row in interval_kkt_function(*center_interval)]
    interval_jacobian = interval_jacobian_function(*box)
    dimension = len(center)
    remainder = [[mp.iv.mpf(0) for _ in range(dimension)] for _ in range(dimension)]
    for i in range(dimension):
        for j in range(dimension):
            value = mp.iv.mpf(1 if i == j else 0)
            for k in range(dimension):
                value -= preconditioner[i][k] * interval_jacobian[k][j]
            remainder[i][j] = value

    corrected_center = []
    for i in range(dimension):
        value = center_interval[i]
        for j in range(dimension):
            value -= preconditioner[i][j] * f_center[j]
        corrected_center.append(value)
    delta = [box[index] - center_interval[index] for index in range(dimension)]
    inclusions = []
    strict = True
    for i in range(dimension):
        value = corrected_center[i]
        for j in range(dimension):
            value += remainder[i][j] * delta[j]
        lower, upper = interval_bounds(value)
        box_lower, box_upper = interval_bounds(box[i])
        contained = box_lower < lower and upper < box_upper
        strict = strict and contained
        inclusions.append(
            {
                "image_lower": mp.nstr(lower, 35),
                "image_upper": mp.nstr(upper, 35),
                "box_lower": mp.nstr(box_lower, 35),
                "box_upper": mp.nstr(box_upper, 35),
                "strictly_contained": bool(contained),
            }
        )
    if not strict:
        raise RuntimeError("Krawczyk image was not strictly contained")
    return (
        {
            "radius": mp.nstr(radius, 20),
            "center_and_preconditioner_are_decimal_point_intervals": True,
            "center_residual_evaluated_with_interval_arithmetic": True,
            "strict_interior_inclusion": True,
            "unique_kkt_root_in_box": True,
            "global_uniqueness_claimed": False,
            "inclusions": inclusions,
        },
        box,
    )


def interval_excludes_zero(value) -> bool:
    lower, upper = interval_bounds(value)
    return upper < 0 or lower > 0


def interval_determinant(matrix: Sequence[Sequence[Any]]):
    """Interval Gaussian-elimination enclosure of a square determinant."""
    size = len(matrix)
    work = [[matrix[i][j] for j in range(size)] for i in range(size)]
    sign = 1
    for column in range(size):
        pivot_candidates = []
        for row in range(column, size):
            lower, upper = interval_bounds(work[row][column])
            distance = min(abs(lower), abs(upper)) if lower * upper > 0 else mp.mpf(0)
            pivot_candidates.append((distance, row))
        _, pivot_row = max(pivot_candidates)
        if not interval_excludes_zero(work[pivot_row][column]):
            raise RuntimeError(f"interval determinant pivot {column} contains zero")
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            factor = work[row][column] / pivot
            for entry in range(column + 1, size):
                work[row][entry] -= factor * work[column][entry]
            work[row][column] = mp.iv.mpf(0)
    determinant = mp.iv.mpf(sign)
    for index in range(size):
        determinant *= work[index][index]
    return determinant


def interval_linear_solve(
    matrix: Sequence[Sequence[Any]], right_hand_side: Sequence[Any]
):
    """Solve an interval linear system with certified nonzero pivots."""
    size = len(matrix)
    work = [
        [mp.iv.mpf(matrix[i][j]) for j in range(size)]
        + [mp.iv.mpf(right_hand_side[i])]
        for i in range(size)
    ]
    sign = 1
    for column in range(size):
        candidates = []
        for row in range(column, size):
            lower, upper = interval_bounds(work[row][column])
            distance = min(abs(lower), abs(upper)) if lower * upper > 0 else mp.mpf(0)
            candidates.append((distance, row))
        _, pivot_row = max(candidates)
        if not interval_excludes_zero(work[pivot_row][column]):
            raise RuntimeError(f"interval solve pivot {column} contains zero")
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        for row in range(column + 1, size):
            factor = work[row][column] / pivot
            for entry in range(column + 1, size + 1):
                work[row][entry] -= factor * work[column][entry]
            work[row][column] = mp.iv.mpf(0)
    solution = [mp.iv.mpf(0) for _ in range(size)]
    for row in range(size - 1, -1, -1):
        value = work[row][size]
        for column in range(row + 1, size):
            value -= work[row][column] * solution[column]
        solution[row] = value / work[row][row]
    determinant = mp.iv.mpf(sign)
    for index in range(size):
        determinant *= work[index][index]
    return solution, determinant


def interval_local_conditions(system, box: Sequence[Any]) -> dict[str, Any]:
    variable_box = box[:9]
    constraint_jacobian_function = interval_lambdify(
        system["variables"], system["constraint_jacobian"]
    )
    constraint_jacobian = constraint_jacobian_function(*variable_box)
    candidates = []
    for omitted in range(9):
        minor = [
            [row[column] for column in range(9) if column != omitted]
            for row in constraint_jacobian
        ]
        right_hand_side = [-row[omitted] for row in constraint_jacobian]
        try:
            solution, determinant = interval_linear_solve(minor, right_hand_side)
        except RuntimeError:
            continue
        distance = min(abs(bound) for bound in interval_bounds(determinant))
        candidates.append((distance, omitted, solution, determinant))
    if not candidates:
        raise RuntimeError("no active-gradient minor certifies LICQ")
    _, licq_index, tangent_solution, licq_value = max(candidates)
    tangent = []
    solution_index = 0
    for variable_index in range(9):
        if variable_index == licq_index:
            tangent.append(mp.iv.mpf(1))
        else:
            tangent.append(tangent_solution[solution_index])
            solution_index += 1

    hessian_function = interval_lambdify(
        system["symbols"], system["lagrangian_hessian"]
    )
    hessian = hessian_function(*box)
    reduced_curvature = mp.iv.mpf(0)
    for i in range(9):
        for j in range(9):
            reduced_curvature += tangent[i] * hessian[i][j] * tangent[j]
    curvature_lower, curvature_upper = interval_bounds(reduced_curvature)

    multiplier_intervals = box[13:17]
    multiplier_bounds = [interval_bounds(value) for value in multiplier_intervals]
    strict_multiplier_signs = all(upper < 0 for _, upper in multiplier_bounds)

    inactive_names = (
        "middle_degree_slack",
        "other_even_eigenvalue",
        "odd_nonzero_eigenvalue_sum",
        "odd_nonzero_eigenvalue_product",
        "mass_a",
        "mass_b",
        "mass_c",
        "density_x",
        "one_minus_density_x",
        "density_y",
        "one_minus_density_y",
        "transverse_A",
        "transverse_B",
        "transverse_C",
        "order_parameter_squared_slack",
        "mu_lower_bound",
            "mu_upper_bound_0.8",
        "alpha_lower_bound",
        "alpha_upper_bound",
        "beta_lower_bound",
        "beta_upper_bound",
        "gamma_lower_bound",
        "gamma_upper_bound",
    )
    inactive_function = interval_lambdify(
        system["variables"], system["inactive"]
    )
    inactive_values = [row[0] for row in inactive_function(*variable_box)]
    inactive_bounds = {
        name: {
            "lower": mp.nstr(interval_bounds(value)[0], 30),
            "upper": mp.nstr(interval_bounds(value)[1], 30),
        }
        for name, value in zip(inactive_names, inactive_values)
    }
    inactive_positive = all(interval_bounds(value)[0] > 0 for value in inactive_values)
    return {
        "licq_minor_omitted_variable_index": licq_index,
        "licq_minor_determinant": {
            "lower": mp.nstr(interval_bounds(licq_value)[0], 30),
            "upper": mp.nstr(interval_bounds(licq_value)[1], 30),
        },
        "licq_certified": interval_excludes_zero(licq_value),
        "active_inequality_multiplier_bounds": [
            {"lower": mp.nstr(lower, 30), "upper": mp.nstr(upper, 30)}
            for lower, upper in multiplier_bounds
        ],
        "strict_multiplier_signs_certified": strict_multiplier_signs,
        "cofactor_tangent_reduced_curvature": {
            "lower": mp.nstr(curvature_lower, 30),
            "upper": mp.nstr(curvature_upper, 30),
        },
        "positive_reduced_curvature_certified": curvature_lower > 0,
        "inactive_constraint_bounds": inactive_bounds,
        "inactive_constraints_positive_certified": inactive_positive,
    }


def diagnostics(system, root_values: Sequence[mp.mpf]) -> dict[str, Any]:
    root_float = np.array([float(value) for value in root_values])
    variables = root_float[:9]
    multipliers = root_float[9:]
    variable_symbols = system["variables"]
    all_symbols = system["symbols"]
    constraint_jacobian_function = sp.lambdify(
        variable_symbols,
        system["constraint_jacobian"],
        modules="numpy",
        cse=True,
    )
    hessian_function = sp.lambdify(
        all_symbols, system["lagrangian_hessian"], modules="numpy", cse=True
    )
    inactive_function = sp.lambdify(
        variable_symbols, system["inactive"], modules="numpy", cse=True
    )
    constraint_jacobian = np.asarray(
        constraint_jacobian_function(*variables), dtype=float
    )
    _, singular_values, right_vectors = np.linalg.svd(constraint_jacobian)
    tangent = right_vectors[-1]
    hessian = np.asarray(hessian_function(*root_float), dtype=float)
    reduced_curvature = float(tangent @ hessian @ tangent)
    inactive = np.asarray(inactive_function(*variables), dtype=float).reshape(-1)
    return {
        "constraint_jacobian_singular_values": [float(value) for value in singular_values],
        "constraint_jacobian_full_row_rank": bool(singular_values[-1] > 1e-8),
        "active_inequality_multipliers": [float(value) for value in multipliers[4:]],
        "multiplier_sign_convention": (
            "For minimizing -mu with constraints g>=0 and L=-mu+lambda*g, "
            "active inequality multipliers must be nonpositive."
        ),
        "strict_multiplier_signs": bool(np.max(multipliers[4:]) < -1e-8),
        "reduced_lagrangian_curvature": reduced_curvature,
        "positive_reduced_curvature": bool(reduced_curvature > 1e-8),
        "inactive_slacks": {
            "middle_degree": float(inactive[0]),
            "other_even_eigenvalue": float(inactive[1]),
            "odd_nonzero_eigenvalue_sum": float(inactive[2]),
            "odd_nonzero_eigenvalue_product": float(inactive[3]),
        },
        "inactive_slacks_positive": bool(np.min(inactive) > 1e-8),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--radius", default="1e-25")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    system = symbolic_system()
    double_root = solve_double(system, initial_variables(args.input))
    root_values = high_precision_root(system, double_root)
    certificate, root_box = krawczyk_certificate(
        system, root_values, mp.mpf(args.radius)
    )
    diagnostic_report = diagnostics(system, root_values)
    interval_conditions = interval_local_conditions(system, root_box)
    variables = root_values[:9]
    multiplier_values = root_values[9:]
    payload = {
        "scope": (
            "A unique limiting KKT root in the certified box on the density "
            "face discovered by the general nine-density search; no global "
            "uniqueness or optimality claim."
        ),
        "variable_order": [str(value) for value in system["variables"]],
        "variables": [mp.nstr(value, 60) for value in variables],
        "multiplier_order": [str(value) for value in system["multipliers"]],
        "multipliers": [mp.nstr(value, 60) for value in multiplier_values],
        "mu": mp.nstr(variables[-1], 60),
        "krawczyk": certificate,
        "interval_local_optimality_conditions": interval_conditions,
        "local_optimality_diagnostics": diagnostic_report,
    }
    checks = [
        certificate["strict_interior_inclusion"],
        interval_conditions["licq_certified"],
        interval_conditions["strict_multiplier_signs_certified"],
        interval_conditions["positive_reduced_curvature_certified"],
        interval_conditions["inactive_constraints_positive_certified"],
    ]
    payload["all_interval_local_certificate_checks_pass"] = all(checks)
    if not payload["all_interval_local_certificate_checks_pass"]:
        raise RuntimeError("one or more local optimality checks failed")
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
