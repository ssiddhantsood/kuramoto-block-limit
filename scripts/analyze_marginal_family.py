#!/usr/bin/env python3
"""Symbolically verify an exact degree-3/4 marginal reflected family.

The family is an important control experiment for the generalized optimizer:
it has larger degree than the stable branch, but three additional quotient
zero modes.  It therefore cannot be reported as a strictly stable solution.
"""

from __future__ import annotations

import json

import sympy as sp


def main() -> None:
    a, spectral_variable = sp.symbols("a spectral_variable", real=True)
    b = sp.Rational(1, 4)
    c = sp.Rational(1, 4) - a
    masses = (a, a, b, b, c, c)
    phases = (
        -sp.pi / 4,
        sp.pi / 4,
        -3 * sp.pi / 4,
        3 * sp.pi / 4,
        sp.pi / 4,
        -sp.pi / 4,
    )
    weights = sp.Matrix(
        [
            [1, 1, 1, 0, 1, 1],
            [1, 1, 0, 1, 1, 1],
            [1, 0, 1, 1, 0, 1],
            [0, 1, 1, 1, 1, 0],
            [1, 1, 0, 1, 1, 1],
            [1, 1, 1, 0, 1, 1],
        ]
    )

    degrees = [
        sp.simplify(sum(masses[j] * weights[i, j] for j in range(6)))
        for i in range(6)
    ]
    torques = [
        sp.simplify(
            sum(
                masses[j]
                * weights[i, j]
                * sp.sin(phases[j] - phases[i])
                for j in range(6)
            )
        )
        for i in range(6)
    ]
    transverse = [
        sp.simplify(
            sum(
                masses[j]
                * weights[i, j]
                * sp.cos(phases[i] - phases[j])
                for j in range(6)
            )
        )
        for i in range(6)
    ]

    laplacian = sp.zeros(6, 6)
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            cosine_weight = (
                masses[j]
                * weights[i, j]
                * sp.cos(phases[i] - phases[j])
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
    even_characteristic = sp.factor(
        (spectral_variable * sp.eye(3) - even).det()
    )
    odd_characteristic = sp.factor(
        (spectral_variable * sp.eye(3) - odd).det()
    )

    # First-order obstruction in the fully general nine-density reflected
    # model.  The two zero densities are constrained to have nonnegative
    # one-sided derivatives, forcing the sum of the three weak-mode openings
    # to be nonpositive.
    pa, pb, pc = sp.symbols("pa pb pc", positive=True)
    ta, tb, tc = sp.symbols("theta_a theta_b theta_c", real=True)
    va, vb, vc = sp.symbols("v_a v_b v_c", real=True)
    sab, tab, sac, tac, sbc, tbc = sp.symbols(
        "s_ab t_ab s_ac t_ac s_bc t_bc", real=True
    )
    pair_masses = (pa, pb, pc)
    pair_phases = (ta, tb, tc)
    reflected = (va, vb, vc)
    same = {(0, 1): sab, (0, 2): sac, (1, 2): sbc}
    opposite = {(0, 1): tab, (0, 2): tac, (1, 2): tbc}

    def ordered(mapping, i, j):
        return mapping[(i, j) if i < j else (j, i)]

    even_general = sp.zeros(3, 3)
    odd_general = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            same_ij = ordered(same, i, j)
            opposite_ij = ordered(opposite, i, j)
            p_ij = same_ij * sp.cos(pair_phases[i] - pair_phases[j]) + opposite_ij * sp.cos(
                pair_phases[i] + pair_phases[j]
            )
            q_ij = same_ij * sp.cos(pair_phases[i] - pair_phases[j]) - opposite_ij * sp.cos(
                pair_phases[i] + pair_phases[j]
            )
            even_general[i, i] += pair_masses[j] * p_ij
            odd_general[i, i] += pair_masses[j] * p_ij
            even_general[i, j] = -sp.sqrt(
                pair_masses[i] * pair_masses[j]
            ) * p_ij
            odd_general[i, j] = -sp.sqrt(
                pair_masses[i] * pair_masses[j]
            ) * q_ij
        odd_general[i, i] += (
            2
            * pair_masses[i]
            * reflected[i]
            * sp.cos(2 * pair_phases[i])
        )

    parameters = (
        pa,
        pb,
        pc,
        ta,
        tb,
        tc,
        va,
        vb,
        vc,
        sab,
        tab,
        sac,
        tac,
        sbc,
        tbc,
    )
    dpa, dpb = sp.symbols("dpa dpb", real=True)
    directions = (
        dpa,
        dpb,
        -dpa - dpb,
        *sp.symbols("dtheta_a dtheta_b dtheta_c", real=True),
        *sp.symbols("dv_a dv_b dv_c", real=True),
        *sp.symbols("ds_ab dt_ab ds_ac dt_ac ds_bc dt_bc", real=True),
    )
    base_substitution = {
        pa: a,
        pb: sp.Rational(1, 4),
        pc: sp.Rational(1, 4) - a,
        ta: -sp.pi / 4,
        tb: -3 * sp.pi / 4,
        tc: sp.pi / 4,
        va: 1,
        vb: 1,
        vc: 1,
        sab: 1,
        tab: 0,
        sac: 1,
        tac: 1,
        sbc: 0,
        tbc: 1,
    }

    def directional_derivative(matrix):
        return matrix.applyfunc(
            lambda expression: sp.simplify(
                sum(
                    sp.diff(expression, parameter) * direction
                    for parameter, direction in zip(parameters, directions)
                ).subs(base_substitution)
            )
        )

    d_even = directional_derivative(even_general)
    d_odd = directional_derivative(odd_general)
    weak_even = sp.Matrix(
        [sp.sqrt(2 * a), -sp.sqrt(2) / 2, sp.sqrt(2 * (sp.Rational(1, 4) - a))]
    )
    odd_zero_one = sp.Matrix([0, 1, 0])
    odd_zero_two = sp.Matrix(
        [2 * sp.sqrt(a), 0, -2 * sp.sqrt(sp.Rational(1, 4) - a)]
    )
    rotation = sp.Matrix(
        [sp.sqrt(2 * a), sp.sqrt(sp.Rational(1, 2)), sp.sqrt(2 * (sp.Rational(1, 4) - a))]
    )
    even_at_base = even_general.subs(base_substitution)
    odd_at_base = odd_general.subs(base_substitution)
    weak_mode_trace = sp.simplify(
        (weak_even.T * d_even * weak_even)[0]
        + (odd_zero_one.T * d_odd * odd_zero_one)[0]
        + (odd_zero_two.T * d_odd * odd_zero_two)[0]
    )
    dt_ab = directions[10]
    ds_bc = directions[13]
    expected_trace = -4 * a * dt_ab - 4 * (sp.Rational(1, 4) - a) * ds_bc
    first_order_identity = sp.simplify(weak_mode_trace - expected_trace) == 0

    expected_characteristic = spectral_variable**2 * (
        4 * spectral_variable - 1
    ) / 4
    checks = {
        "all_degrees_are_three_quarters": all(
            sp.simplify(value - sp.Rational(3, 4)) == 0 for value in degrees
        ),
        "all_torques_are_zero": all(value == 0 for value in torques),
        "all_transverse_values_are_one_quarter": all(
            sp.simplify(value - sp.Rational(1, 4)) == 0
            for value in transverse
        ),
        "even_characteristic_polynomial": sp.simplify(
            even_characteristic - expected_characteristic
        )
        == 0,
        "odd_characteristic_polynomial": sp.simplify(
            odd_characteristic - expected_characteristic
        )
        == 0,
        "first_order_weak_mode_obstruction_identity": first_order_identity,
        "weak_even_vector_is_unit": sp.simplify((weak_even.T * weak_even)[0] - 1) == 0,
        "weak_even_vector_is_physical": sp.simplify((rotation.T * weak_even)[0]) == 0,
        "weak_even_vector_is_in_kernel": all(
            sp.simplify(value) == 0 for value in even_at_base * weak_even
        ),
        "odd_zero_vectors_are_orthonormal": (
            sp.simplify((odd_zero_one.T * odd_zero_one)[0] - 1) == 0
            and sp.simplify((odd_zero_two.T * odd_zero_two)[0] - 1) == 0
            and sp.simplify((odd_zero_one.T * odd_zero_two)[0]) == 0
        ),
        "odd_zero_vectors_are_in_kernel": all(
            sp.simplify(value) == 0
            for value in list(odd_at_base * odd_zero_one)
            + list(odd_at_base * odd_zero_two)
        ),
    }
    if not all(checks.values()):
        raise RuntimeError(f"symbolic marginal-family check failed: {checks}")

    payload = {
        "parameter_range": "0 < a < 1/4",
        "pair_masses": ["a", "1/4", "1/4-a"],
        "phases": ["-pi/4", "-3*pi/4", "pi/4"],
        "density_order": [
            "ref_A",
            "ref_B",
            "ref_C",
            "same_AB",
            "opposite_AB",
            "same_AC",
            "opposite_AC",
            "same_BC",
            "opposite_BC",
        ],
        "densities": [1, 1, 1, 1, 0, 1, 1, 0, 1],
        "degree": "3/4",
        "transverse_value": "1/4",
        "even_characteristic_polynomial": "lambda^2*(lambda-1/4)",
        "odd_characteristic_polynomial": "lambda^2*(lambda-1/4)",
        "extra_quotient_zero_modes": 3,
        "strictly_stable_modulo_rotation": False,
        "first_order_weak_mode_sum": "-4*a*d(t_AB)-4*(1/4-a)*d(s_BC)",
        "first_order_conclusion": (
            "Because t_AB=s_BC=0 and feasible one-sided density derivatives "
            "are nonnegative, no first-order perturbation can make all three "
            "extra zero modes strictly positive. Higher-order opening is not "
            "excluded by this identity."
        ),
        "checks": checks,
        "interpretation": (
            "This exact family is a marginal closure point, not a stable "
            "counterexample or a degree-3/4 construction."
        ),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
