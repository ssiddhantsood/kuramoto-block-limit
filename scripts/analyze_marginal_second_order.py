#!/usr/bin/env python3
"""Verify the second-order obstruction at the exact degree-3/4 family.

The existing first-order trace identity forces every potentially stable
tangent to have zero first-order opening of all three weak modes.  This script
computes the complete second spectral variation, including coupling to the
strong eigenvectors, and proves an exact cancellation.  Consequently no
one-sided C^2 equilibrium branch can open the modes at linear or quadratic
order.

This is a jet obstruction, not a punctured-neighborhood theorem: higher-order
or nonanalytic approaches are not excluded.
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
    mass_a, mass_b, family_parameter = sp.symbols(
        "mass_a mass_b family_parameter", positive=True
    )
    mass_c = sp.Rational(1, 2) - mass_a - mass_b
    pair_masses = [mass_a, mass_b, mass_c]
    v_a, v_b, v_c, s_ab, t_ab, s_ac, t_ac, s_bc, t_bc = sp.symbols(
        "v_a v_b v_c s_ab t_ab s_ac t_ac s_bc t_bc"
    )
    theta_a, theta_b, theta_c = sp.symbols("theta_a theta_b theta_c")
    angles = [theta_a, theta_b, theta_c]
    reflected = [v_a, v_b, v_c]
    same = {(0, 1): s_ab, (0, 2): s_ac, (1, 2): s_bc}
    opposite = {(0, 1): t_ab, (0, 2): t_ac, (1, 2): t_bc}

    def pair_value(mapping, left, right):
        return mapping[(min(left, right), max(left, right))]

    even_weights = {}
    odd_weights = {}
    for left in range(3):
        for right in range(left + 1, 3):
            same_value = pair_value(same, left, right)
            opposite_value = pair_value(opposite, left, right)
            even_weights[left, right] = (
                same_value * sp.cos(angles[left] - angles[right])
                + opposite_value * sp.cos(angles[left] + angles[right])
            )
            odd_weights[left, right] = (
                same_value * sp.cos(angles[left] - angles[right])
                - opposite_value * sp.cos(angles[left] + angles[right])
            )

    even = sp.zeros(3)
    odd = sp.zeros(3)
    torques = []
    for left in range(3):
        for right in range(3):
            if left == right:
                continue
            pair = (min(left, right), max(left, right))
            even_weight = even_weights[pair]
            odd_weight = odd_weights[pair]
            even[left, left] += pair_masses[right] * even_weight
            odd[left, left] += pair_masses[right] * even_weight
            even[left, right] -= sp.sqrt(
                pair_masses[left] * pair_masses[right]
            ) * even_weight
            odd[left, right] -= sp.sqrt(
                pair_masses[left] * pair_masses[right]
            ) * odd_weight
        odd[left, left] += (
            2
            * pair_masses[left]
            * reflected[left]
            * sp.cos(2 * angles[left])
        )
        torque = (
            -pair_masses[left]
            * reflected[left]
            * sp.sin(2 * angles[left])
        )
        for right in range(3):
            if left == right:
                continue
            torque += pair_masses[right] * (
                pair_value(same, left, right)
                * sp.sin(angles[right] - angles[left])
                - pair_value(opposite, left, right)
                * sp.sin(angles[right] + angles[left])
            )
        torques.append(torque)

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
    directions = sp.symbols(
        "d_mass_a d_mass_b d_v_a d_v_b d_v_c "
        "d_s_ab d_t_ab d_s_ac d_t_ac d_s_bc d_t_bc "
        "d_theta_a d_theta_b d_theta_c"
    )
    (
        d_mass_a,
        d_mass_b,
        d_v_a,
        d_v_b,
        d_v_c,
        d_s_ab,
        d_t_ab,
        d_s_ac,
        d_t_ac,
        d_s_bc,
        d_t_bc,
        d_theta_a,
        d_theta_b,
        d_theta_c,
    ) = directions
    a = family_parameter
    c = sp.Rational(1, 4) - a
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

    weak_even = sp.sqrt(2) * sp.Matrix(
        [sp.sqrt(a), -sp.Rational(1, 2), sp.sqrt(c)]
    )
    rotation = sp.sqrt(2) * sp.Matrix(
        [sp.sqrt(a), sp.Rational(1, 2), sp.sqrt(c)]
    )
    odd_zero_one = sp.Matrix([0, 1, 0])
    odd_zero_two = 2 * sp.Matrix([sp.sqrt(a), 0, -sp.sqrt(c)])
    strong_even = 2 * sp.Matrix([sp.sqrt(c), 0, -sp.sqrt(a)])
    strong_odd = 2 * sp.Matrix([sp.sqrt(c), 0, sp.sqrt(a)])
    even_base = even.subs(base)
    odd_base = odd.subs(base)
    strong_eigenvalue = sp.Rational(1, 4)

    vector_checks = {
        "weak_even_unit": sp.simplify(weak_even.dot(weak_even) - 1) == 0,
        "rotation_unit": sp.simplify(rotation.dot(rotation) - 1) == 0,
        "strong_even_unit": sp.simplify(strong_even.dot(strong_even) - 1) == 0,
        "odd_zero_basis_orthonormal": (
            sp.simplify(odd_zero_one.dot(odd_zero_one) - 1) == 0
            and sp.simplify(odd_zero_two.dot(odd_zero_two) - 1) == 0
            and sp.simplify(odd_zero_one.dot(odd_zero_two)) == 0
        ),
        "strong_odd_unit": sp.simplify(strong_odd.dot(strong_odd) - 1) == 0,
        "even_base_rank_one": zero_matrix(
            even_base - strong_eigenvalue * strong_even * strong_even.T
        ),
        "odd_base_rank_one": zero_matrix(
            odd_base - strong_eigenvalue * strong_odd * strong_odd.T
        ),
    }

    fixed_kernel_sum = sp.simplify(
        (weak_even.T * even * weak_even)[0]
        + (odd_zero_one.T * odd * odd_zero_one)[0]
        + (odd_zero_two.T * odd * odd_zero_two)[0]
    )
    gradient = sp.Matrix(
        [sp.simplify(sp.diff(fixed_kernel_sum, variable).subs(base)) for variable in variables]
    )
    expected_gradient = sp.zeros(len(variables), 1)
    expected_gradient[6] = -4 * a
    expected_gradient[9] = -4 * c
    first_order_trace_identity = zero_matrix(gradient - expected_gradient)

    def directional_derivative(matrix):
        return sum(
            (
                sp.diff(matrix, variable) * direction
                for variable, direction in zip(variables, directions)
            ),
            sp.zeros(*matrix.shape),
        ).subs(base)

    even_first = directional_derivative(even)
    odd_first = directional_derivative(odd)
    torque_first = sp.Matrix(
        [
            sp.simplify(
                sum(
                    sp.diff(torque, variable) * direction
                    for variable, direction in zip(variables, directions)
                ).subs(base)
            )
            for torque in torques
        ]
    )
    endpoint_flat = {d_t_ab: 0, d_s_bc: 0}
    expected_torque = sp.Matrix(
        [
            -2 * d_mass_b
            + a * d_v_a
            - sp.Rational(1, 4) * d_s_ab
            + c * d_s_ac
            - c * (d_theta_a + d_theta_c),
            -2 * d_mass_b
            - sp.Rational(1, 4) * d_v_b
            + a * d_s_ab
            + c * d_t_bc,
            2 * d_mass_b
            - c * d_v_c
            - a * d_s_ac
            + sp.Rational(1, 4) * d_t_bc
            - a * (d_theta_a + d_theta_c),
        ]
    )
    torque_linearization_check = zero_matrix(
        torque_first.subs(endpoint_flat) - expected_torque
    )
    torque_solution = {
        d_v_a: (
            2 * d_mass_b
            + sp.Rational(1, 4) * d_s_ab
            - c * d_s_ac
            + c * (d_theta_a + d_theta_c)
        )
        / a,
        d_v_b: -8 * d_mass_b + 4 * a * d_s_ab + 4 * c * d_t_bc,
        d_v_c: (
            2 * d_mass_b
            - a * d_s_ac
            + sp.Rational(1, 4) * d_t_bc
            - a * (d_theta_a + d_theta_c)
        )
        / c,
    }
    compressed_even = sp.simplify(
        (weak_even.T * even_first * weak_even)[0]
    ).subs(endpoint_flat).subs(torque_solution)
    compressed_odd = sp.Matrix(
        [
            [
                sp.simplify((left.T * odd_first * right)[0])
                .subs(endpoint_flat)
                .subs(torque_solution)
                for right in (odd_zero_one, odd_zero_two)
            ]
            for left in (odd_zero_one, odd_zero_two)
        ]
    )
    expected_compressed_even = -(
        4 * a * d_theta_a
        + 4 * a * d_theta_c
        - d_theta_b
        - d_theta_c
    ) / 2
    expected_compressed_odd = sp.Matrix(
        [
            [
                -(
                    4 * a * d_theta_a
                    + 4 * a * d_theta_c
                    + 3 * d_theta_b
                    - d_theta_c
                )
                / 4,
                (
                    4 * a * d_theta_a
                    + 4 * a * d_theta_c
                    - d_theta_b
                    - d_theta_c
                )
                / 4,
            ],
            [
                (
                    4 * a * d_theta_a
                    + 4 * a * d_theta_c
                    - d_theta_b
                    - d_theta_c
                )
                / 4,
                (
                    12 * a * d_theta_a
                    + 12 * a * d_theta_c
                    + d_theta_b
                    - 3 * d_theta_c
                )
                / 4,
            ],
        ]
    )
    compression_check = (
        sp.simplify(compressed_even - expected_compressed_even) == 0
        and zero_matrix(compressed_odd - expected_compressed_odd)
    )
    phase_nullspace = sp.Matrix(
        [expected_compressed_even, *list(expected_compressed_odd)]
    ).jacobian([d_theta_a, d_theta_b, d_theta_c]).nullspace()
    phase_reduction_check = (
        len(phase_nullspace) == 1
        and zero_matrix(
            phase_nullspace[0].cross(sp.Matrix([c, 0, a]))
        )
    )

    hessian = sp.hessian(fixed_kernel_sum, variables).subs(base)
    fixed_quadratic = sp.expand(
        (sp.Matrix(directions).T * hessian * sp.Matrix(directions))[0]
    )
    phase_solution = {
        d_theta_b: 0,
        d_theta_c: a / c * d_theta_a,
    }
    def restrict_to_critical(expression):
        return (
            expression.subs(endpoint_flat)
            .subs(torque_solution)
            .subs(phase_solution)
        )

    fixed_quadratic = sp.factor(
        sp.cancel(restrict_to_critical(fixed_quadratic))
    )
    expected_fixed_quadratic = (
        3 * (d_mass_a + 4 * a * d_mass_b) ** 2
        + 32 * a**2 * d_theta_a**2
    ) / (4 * a * (1 - 4 * a))
    fixed_quadratic_check = (
        sp.simplify(fixed_quadratic - expected_fixed_quadratic) == 0
    )

    even_coupling = sp.simplify(
        (weak_even.T * even_first * strong_even)[0]
    )
    even_coupling = restrict_to_critical(even_coupling)
    odd_coupling_one = sp.simplify(
        (odd_zero_one.T * odd_first * strong_odd)[0]
    )
    odd_coupling_one = restrict_to_critical(odd_coupling_one)
    odd_coupling_two = sp.simplify(
        (odd_zero_two.T * odd_first * strong_odd)[0]
    )
    odd_coupling_two = restrict_to_critical(odd_coupling_two)
    coupling = sp.simplify(
        2
        / strong_eigenvalue
        * (even_coupling**2 + odd_coupling_one**2 + odd_coupling_two**2)
    )
    spectral_cancellation = sp.simplify(
        fixed_quadratic - coupling
    ) == 0

    # The even zero cluster also contains rotation.  Its separate fixed-basis
    # second variation cancels its strong-mode coupling because the moving
    # rotation eigenvalue is identically zero.
    rotation_rayleigh = sp.simplify((rotation.T * even * rotation)[0])
    rotation_hessian = sp.hessian(rotation_rayleigh, variables).subs(base)
    rotation_fixed_quadratic = sp.expand(
        (sp.Matrix(directions).T * rotation_hessian * sp.Matrix(directions))[0]
    )
    rotation_fixed_quadratic = restrict_to_critical(rotation_fixed_quadratic)
    rotation_coupling = sp.simplify(
        (rotation.T * even_first * strong_even)[0]
    )
    rotation_coupling = restrict_to_critical(rotation_coupling)
    rotation_cancellation = sp.simplify(
        rotation_fixed_quadratic
        - 2 / strong_eigenvalue * rotation_coupling**2
    ) == 0

    checks = {
        **vector_checks,
        "first_order_trace_gradient": first_order_trace_identity,
        "linearized_torques": torque_linearization_check,
        "compressed_first_derivatives": compression_check,
        "critical_phase_nullspace": phase_reduction_check,
        "fixed_kernel_quadratic": fixed_quadratic_check,
        "strong_mode_coupling_cancels_quadratic": spectral_cancellation,
        "moving_rotation_contribution_cancels": rotation_cancellation,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"second-order symbolic checks failed: {failed}")

    payload = {
        "scope": (
            "Exact one-sided C2 jet calculation at the degree-3/4 marginal "
            "family; not a punctured-neighborhood exclusion."
        ),
        "parameter_range": "0 < a < 1/4",
        "first_order_spectral_sum": "-4*a*d(t_AB)-4*(1/4-a)*d(s_BC)",
        "critical_tangent_conditions": [
            "d(t_AB)=d(s_BC)=0",
            "d(theta_B)=0",
            "(1/4-a)*d(theta_C)=a*d(theta_A)",
            "linearized torque equations determine d(v_A), d(v_B), d(v_C)",
        ],
        "fixed_kernel_second_variation": (
            "[3*(dA+4*a*dB)^2+32*a^2*d(theta_A)^2]/"
            "[4*a*(1-4*a)]"
        ),
        "strong_mode_coupling": "exactly the same expression",
        "spectral_second_variation": (
            "-4*a*d2(t_AB)-4*(1/4-a)*d2(s_BC) <= 0"
        ),
        "conclusion": (
            "A one-sided C2 stable equilibrium branch cannot open the three "
            "weak modes at first or second order. Stability forces every "
            "individual linear and quadratic weak opening to vanish."
        ),
        "limitation": (
            "Cubic, higher-order, Puiseux/nonanalytic, and arbitrary-sequence "
            "approaches are not excluded."
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
