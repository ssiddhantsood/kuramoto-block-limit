#!/usr/bin/env python3
"""Verify the exact marginal sheet and its second-jet tangent obstruction.

The degree-3/4 family lies on a larger six-dimensional fixed-phase sheet.
This script proves the sheet formulas symbolically and computes the exact
compatibility condition for a first critical jet to admit a second-flat jet.
It intentionally makes no all-orders or punctured-neighborhood claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


def zero_matrix(matrix: sp.Matrix) -> bool:
    return all(sp.simplify(value) == 0 for value in matrix)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    mass_a, mass_b = sp.symbols("mass_a mass_b", positive=True)
    mass_c = sp.Rational(1, 2) - mass_a - mass_b
    pair_masses = [mass_a, mass_b, mass_c]
    v_a, v_b, v_c, s_ab, t_ab, s_ac, t_ac, s_bc, t_bc = sp.symbols(
        "v_a v_b v_c s_ab t_ab s_ac t_ac s_bc t_bc"
    )
    theta_a, theta_b, theta_c = sp.symbols("theta_a theta_b theta_c", real=True)
    angles = [theta_a, theta_b, theta_c]
    reflected = [v_a, v_b, v_c]
    same = {(0, 1): s_ab, (0, 2): s_ac, (1, 2): s_bc}
    opposite = {(0, 1): t_ab, (0, 2): t_ac, (1, 2): t_bc}

    def pair_value(mapping, left, right):
        return mapping[(min(left, right), max(left, right))]

    even = sp.zeros(3)
    odd = sp.zeros(3)
    torques = []
    for left in range(3):
        torque = -pair_masses[left] * reflected[left] * sp.sin(2 * angles[left])
        for right in range(3):
            if left == right:
                continue
            same_value = pair_value(same, left, right)
            opposite_value = pair_value(opposite, left, right)
            even_weight = same_value * sp.cos(
                angles[left] - angles[right]
            ) + opposite_value * sp.cos(angles[left] + angles[right])
            odd_weight = same_value * sp.cos(
                angles[left] - angles[right]
            ) - opposite_value * sp.cos(angles[left] + angles[right])
            even[left, left] += pair_masses[right] * even_weight
            odd[left, left] += pair_masses[right] * even_weight
            even[left, right] -= (
                sp.sqrt(pair_masses[left] * pair_masses[right]) * even_weight
            )
            odd[left, right] -= (
                sp.sqrt(pair_masses[left] * pair_masses[right]) * odd_weight
            )
            torque += pair_masses[right] * (
                same_value * sp.sin(angles[right] - angles[left])
                - opposite_value * sp.sin(angles[right] + angles[left])
            )
        odd[left, left] += (
            2 * pair_masses[left] * reflected[left] * sp.cos(2 * angles[left])
        )
        torques.append(torque)

    # The exact fixed-phase marginal sheet.
    fixed_phases = {
        theta_a: -sp.pi / 4,
        theta_b: -3 * sp.pi / 4,
        theta_c: sp.pi / 4,
        t_ab: 0,
        s_bc: 0,
    }
    marginal_reflected = {
        v_a: (mass_b * s_ab - mass_c * s_ac) / mass_a,
        v_b: (mass_a * s_ab + mass_c * t_bc) / mass_b,
        v_c: (mass_b * t_bc - mass_a * s_ac) / mass_c,
    }
    marginal_substitution = {**fixed_phases, **marginal_reflected}
    expected_even = t_ac * sp.Matrix(
        [
            [mass_c, 0, -sp.sqrt(mass_a * mass_c)],
            [0, 0, 0],
            [-sp.sqrt(mass_a * mass_c), 0, mass_a],
        ]
    )
    expected_odd = t_ac * sp.Matrix(
        [
            [mass_c, 0, sp.sqrt(mass_a * mass_c)],
            [0, 0, 0],
            [sp.sqrt(mass_a * mass_c), 0, mass_a],
        ]
    )
    manifold_checks = {
        "marginal_torques_vanish": all(
            sp.simplify(value.subs(marginal_substitution)) == 0 for value in torques
        ),
        "marginal_even_formula": zero_matrix(
            even.subs(marginal_substitution) - expected_even
        ),
        "marginal_odd_formula": zero_matrix(
            odd.subs(marginal_substitution) - expected_odd
        ),
        "marginal_even_rank_one_identity": zero_matrix(
            expected_even**2 - t_ac * (mass_a + mass_c) * expected_even
        ),
        "marginal_odd_rank_one_identity": zero_matrix(
            expected_odd**2 - t_ac * (mass_a + mass_c) * expected_odd
        ),
    }

    # Work at a general point of the degree-3/4 subfamily.
    a = sp.symbols("a", positive=True)
    c = sp.Rational(1, 4) - a
    variables = [
        mass_a,
        mass_b,
        v_a,
        v_b,
        v_c,
        s_ab,
        t_ab,
        s_ac,
        t_ac,
        s_bc,
        t_bc,
        theta_a,
        theta_b,
        theta_c,
    ]
    base = {
        mass_a: a,
        mass_b: sp.Rational(1, 4),
        v_a: 1,
        v_b: 1,
        v_c: 1,
        s_ab: 1,
        t_ab: 0,
        s_ac: 1,
        t_ac: 1,
        s_bc: 0,
        t_bc: 1,
        theta_a: -sp.pi / 4,
        theta_b: -3 * sp.pi / 4,
        theta_c: sp.pi / 4,
    }
    rotation = sp.sqrt(2) * sp.Matrix([sp.sqrt(a), sp.Rational(1, 2), sp.sqrt(c)])
    weak_even = sp.sqrt(2) * sp.Matrix([sp.sqrt(a), -sp.Rational(1, 2), sp.sqrt(c)])
    strong_even = 2 * sp.Matrix([sp.sqrt(c), 0, -sp.sqrt(a)])
    odd_zero_one = sp.Matrix([0, 1, 0])
    odd_zero_two = 2 * sp.Matrix([sp.sqrt(a), 0, -sp.sqrt(c)])
    strong_odd = 2 * sp.Matrix([sp.sqrt(c), 0, sp.sqrt(a)])
    even_basis = sp.Matrix.hstack(rotation, weak_even, strong_even)
    odd_basis = sp.Matrix.hstack(odd_zero_one, odd_zero_two, strong_odd)

    parameter = sp.symbols("parameter", real=True)

    def directional_derivative(expression, direction, order):
        line = {
            variable: base[variable] + parameter * value
            for variable, value in zip(variables, direction)
        }
        return sp.simplify(
            sp.diff(expression.subs(line), parameter, order).subs(parameter, 0)
        )

    def matrix_derivative(matrix, direction, order):
        return matrix.applyfunc(
            lambda expression: directional_derivative(expression, direction, order)
        )

    # A generic second derivative, with the two missing densities flat.
    second_symbols = sp.symbols(
        "e_mass_a e_mass_b e_v_a e_v_b e_v_c "
        "e_s_ab e_t_ab e_s_ac e_t_ac e_s_bc e_t_bc "
        "e_theta_a e_theta_b e_theta_c"
    )
    second_direction = list(second_symbols)
    second_direction[6] = 0
    second_direction[9] = 0
    second_direction = sp.Matrix(second_direction)
    torque_first = sp.Matrix(
        [directional_derivative(value, second_direction, 1) for value in torques]
    )
    even_first = (
        even_basis.T * matrix_derivative(even, second_direction, 1) * even_basis
    )
    odd_first = odd_basis.T * matrix_derivative(odd, second_direction, 1) * odd_basis
    first_system = sp.Matrix(
        list(torque_first)
        + [even_first[0, 0], even_first[0, 1], even_first[1, 1]]
        + [odd_first[0, 0], odd_first[0, 1], odd_first[1, 1]]
    )
    retained_second_symbols = [
        second_symbols[index]
        for index in range(len(second_symbols))
        if index not in (6, 9)
    ]
    reduced_jacobian = sp.simplify(first_system.jacobian(retained_second_symbols))

    # Parameterize every first critical direction from the first-order trace,
    # torque, and compressed weak-matrix equations.
    (
        d_mass_a,
        d_mass_b,
        d_s_ab,
        d_s_ac,
        d_t_ac,
        d_t_bc,
        d_theta_a,
    ) = sp.symbols("d_mass_a d_mass_b d_s_ab d_s_ac d_t_ac d_t_bc d_theta_a")
    first_direction = [None] * len(variables)
    first_direction[0] = d_mass_a
    first_direction[1] = d_mass_b
    first_direction[5] = d_s_ab
    first_direction[6] = 0
    first_direction[7] = d_s_ac
    first_direction[8] = d_t_ac
    first_direction[9] = 0
    first_direction[10] = d_t_bc
    first_direction[11] = d_theta_a
    first_direction[12] = 0
    first_direction[13] = a / c * d_theta_a
    first_direction[2] = (
        2 * d_mass_b
        + sp.Rational(1, 4) * d_s_ab
        - c * d_s_ac
        + c * (d_theta_a + a / c * d_theta_a)
    ) / a
    first_direction[3] = -8 * d_mass_b + 4 * a * d_s_ab + 4 * c * d_t_bc
    first_direction[4] = (
        2 * d_mass_b
        - a * d_s_ac
        + sp.Rational(1, 4) * d_t_bc
        - a * (d_theta_a + a / c * d_theta_a)
    ) / c
    first_direction = sp.Matrix(first_direction)
    first_critical_check = zero_matrix(
        first_system.subs(dict(zip(second_symbols, first_direction)))
    )

    # If M(t)=M0+t M1+t^2 M2/2 and M0 has strong eigenvalue 1/4,
    # then the weak Schur complement has second derivative
    # A2-(2/(1/4)) b1 b1^T = A2-8 b1 b1^T.
    torque_second = sp.Matrix(
        [directional_derivative(value, first_direction, 2) for value in torques]
    )
    even_first_critical = (
        even_basis.T * matrix_derivative(even, first_direction, 1) * even_basis
    )
    odd_first_critical = (
        odd_basis.T * matrix_derivative(odd, first_direction, 1) * odd_basis
    )
    even_second = (
        even_basis.T * matrix_derivative(even, first_direction, 2) * even_basis
    )
    odd_second = odd_basis.T * matrix_derivative(odd, first_direction, 2) * odd_basis
    even_weak_second = (
        even_second[:2, :2]
        - 8 * even_first_critical[:2, 2:3] * even_first_critical[2:3, :2]
    )
    odd_weak_second = (
        odd_second[:2, :2]
        - 8 * odd_first_critical[:2, 2:3] * odd_first_critical[2:3, :2]
    )
    second_forcing = sp.Matrix(
        list(torque_second)
        + [
            even_weak_second[0, 0],
            even_weak_second[0, 1],
            even_weak_second[1, 1],
            odd_weak_second[0, 0],
            odd_weak_second[0, 1],
            odd_weak_second[1, 1],
        ]
    )
    left_nullspace = reduced_jacobian.T.nullspace()
    compatibility = [
        sp.factor((vector.T * second_forcing)[0]) for vector in left_nullspace
    ]
    expected_compatibility = [
        0,
        0,
        4 * a * d_theta_a**2 / (4 * a - 1),
        0,
    ]
    compatibility_check = len(compatibility) == len(expected_compatibility) and all(
        sp.simplify(actual - expected) == 0
        for actual, expected in zip(compatibility, expected_compatibility)
    )

    # With the phase derivative zero, the six remaining directions are
    # exactly the tangent space of the explicit marginal sheet.
    sheet_variables = [mass_a, mass_b, s_ab, s_ac, t_ac, t_bc]
    sheet_directions = [
        d_mass_a,
        d_mass_b,
        d_s_ab,
        d_s_ac,
        d_t_ac,
        d_t_bc,
    ]
    sheet_base = {
        mass_a: a,
        mass_b: sp.Rational(1, 4),
        s_ab: 1,
        s_ac: 1,
        t_ac: 1,
        t_bc: 1,
    }
    sheet_reflected_derivative = sp.Matrix(
        [
            sp.simplify(
                sum(
                    sp.diff(expression, variable) * direction
                    for variable, direction in zip(sheet_variables, sheet_directions)
                ).subs(sheet_base)
            )
            for expression in marginal_reflected.values()
        ]
    )
    tangent_check = zero_matrix(
        sheet_reflected_derivative - first_direction[2:5, :].subs(d_theta_a, 0)
    )

    # A non-uniform exact path on the sheet.  It shows why second-flat
    # directions are not restricted to mass variation plus uniform scaling.
    path_parameter = sp.symbols("path_parameter", real=True)
    explicit_path = {
        mass_a: a,
        mass_b: sp.Rational(1, 4),
        s_ab: 1 - path_parameter,
        s_ac: 1 - path_parameter,
        t_ac: 1 - path_parameter / 2,
        t_bc: 1 - path_parameter,
    }
    explicit_reflected = {
        key: sp.simplify(value.subs(explicit_path))
        for key, value in marginal_reflected.items()
    }
    explicit_path_check = all(
        sp.simplify(value - (1 - path_parameter)) == 0
        for value in explicit_reflected.values()
    )
    explicit_strong_eigenvalue = sp.simplify(t_ac * (mass_a + mass_c)).subs(
        explicit_path
    )
    explicit_cubic_flatness = (
        sp.diff(explicit_strong_eigenvalue, path_parameter, 3) == 0
        and explicit_path_check
    )

    checks = {
        **manifold_checks,
        "first_critical_parameterization": first_critical_check,
        "reduced_jacobian_rank_five": reduced_jacobian.rank() == 5,
        "compatibility_cokernel_dimension_four": len(left_nullspace) == 4,
        "second_jet_compatibility": compatibility_check,
        "compatible_space_equals_sheet_tangent": tangent_check,
        "explicit_nonuniform_path": explicit_path_check,
        "explicit_path_is_flat_through_cubic": explicit_cubic_flatness,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"marginal-manifold symbolic checks failed: {failed}")

    payload = {
        "scope": (
            "Exact marginal-sheet and second-jet calculation. No all-orders "
            "or punctured-neighborhood exclusion is claimed."
        ),
        "mass_domain": "A>0, B>0, C=1/2-A-B>0",
        "fixed_phases": ["-pi/4", "-3*pi/4", "pi/4"],
        "missing_densities": {"t_AB": 0, "s_BC": 0},
        "free_coordinates": ["A", "B", "s_AB", "s_AC", "t_AC", "t_BC"],
        "reflected_densities": {
            "v_A": "(B*s_AB-C*s_AC)/A",
            "v_B": "(A*s_AB+C*t_BC)/B",
            "v_C": "(B*t_BC-A*s_AC)/C",
        },
        "even_matrix": "t_AC*[[C,0,-sqrt(A*C)],[0,0,0],[-sqrt(A*C),0,A]]",
        "odd_matrix": "t_AC*[[C,0,+sqrt(A*C)],[0,0,0],[+sqrt(A*C),0,A]]",
        "strong_eigenvalue": "t_AC*(A+C)",
        "zero_eigenvalues_including_rotation": [0, 0, 0, 0],
        "extra_nonrotation_zero_eigenvalues": [0, 0, 0],
        "degree_three_quarters_base": {
            "masses": ["a", "1/4", "1/4-a"],
            "densities": [1, 1, 1, 1, 0, 1, 1, 0, 1],
        },
        "first_critical_phase_conditions": [
            "d(theta_B)=0",
            "(1/4-a)*d(theta_C)=a*d(theta_A)",
        ],
        "second_jet_compatibility": [
            "0",
            "0",
            "4*a*d(theta_A)^2/(4*a-1)",
            "0",
        ],
        "second_jet_conclusion": (
            "Since 0<a<1/4, a second-flat extension exists only when all "
            "three phase derivatives vanish. The six remaining first "
            "directions are exactly tangent to the explicit marginal sheet."
        ),
        "explicit_nonuniform_flat_path": {
            "mass_derivatives": [0, 0],
            "existing_density_derivatives": {
                "v_A": -1,
                "v_B": -1,
                "v_C": -1,
                "s_AB": -1,
                "s_AC": -1,
                "t_AC": "-1/2",
                "t_BC": -1,
            },
            "second_jet": 0,
            "third_jet": 0,
            "behavior": (
                "The path remains exactly marginal: every weak eigenvalue "
                "is identically zero and the strong eigenvalue is "
                "(1-t/2)/4."
            ),
        },
        "all_orders_status": (
            "The exact sheet and transverse quadratic compatibility suggest "
            "a Lyapunov-Schmidt or analytic-curve valuation proof, but the "
            "valuation step has not been completed here."
        ),
        "checks": checks,
    }
    text = json.dumps(payload, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
