#!/usr/bin/env python3
"""Check exact identities for the reflected three-pair fixed-angle slice.

At phases (-pi/4, -3*pi/4, pi/4), weak even stability forces two
cross-pair densities to vanish.  Torque then eliminates the three
reflection-pair densities.  The resulting degree slack identities prove
minimum degree at most 3/4, with equality only for the binary marginal
family.  The inequalities use positive masses and densities in [0,1];
SymPy checks their exact algebraic identities, while the accompanying
document spells out the nonnegative-term argument.

This makes no claim about nearby phases or a global three-pair bound.
"""

from __future__ import annotations

import sympy as sp


def require_zero(label: str, expression: sp.Expr | sp.MatrixBase) -> None:
    entries = (
        list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
    )
    for entry in entries:
        if sp.simplify(sp.trigsimp(entry)) != 0:
            raise AssertionError(f"{label}: nonzero expression {entry}")
    print(f"PASS {label}")


def main() -> None:
    a, b, c = sp.symbols("a b c", positive=True)
    v_a, v_b, v_c = sp.symbols("v_a v_b v_c", real=True)
    s_ab, t_ab, s_ac, t_ac, s_bc, t_bc = sp.symbols(
        "s_ab t_ab s_ac t_ac s_bc t_bc", real=True
    )
    phases = (-sp.pi / 4, -3 * sp.pi / 4, sp.pi / 4)
    six_masses = (a, a, b, b, c, c)
    six_phases = tuple(value for phase in phases for value in (phase, -phase))

    weights = sp.eye(6)
    for pair, value in enumerate((v_a, v_b, v_c)):
        plus, minus = 2 * pair, 2 * pair + 1
        weights[plus, minus] = weights[minus, plus] = value
    for left, right, same, opposite in (
        (0, 1, s_ab, t_ab),
        (0, 2, s_ac, t_ac),
        (1, 2, s_bc, t_bc),
    ):
        for sign in (0, 1):
            x, y = 2 * left + sign, 2 * right + sign
            weights[x, y] = weights[y, x] = same
            x, y = 2 * left + sign, 2 * right + 1 - sign
            weights[x, y] = weights[y, x] = opposite

    degrees = sp.Matrix(
        [sum(six_masses[j] * weights[i, j] for j in range(6)) for i in (0, 2, 4)]
    )
    torques = sp.Matrix(
        [
            sum(
                six_masses[j] * weights[i, j] * sp.sin(six_phases[j] - six_phases[i])
                for j in range(6)
            )
            for i in (0, 2, 4)
        ]
    )
    expected_torques = sp.Matrix(
        [
            a * v_a - b * s_ab + c * s_ac,
            a * s_ab - b * v_b + c * t_bc,
            -a * s_ac + b * t_bc - c * v_c,
        ]
    )
    require_zero("torque formulas", torques - expected_torques)

    hessian = sp.zeros(6)
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            cosine_weight = weights[i, j] * sp.cos(six_phases[i] - six_phases[j])
            hessian[i, i] += six_masses[j] * cosine_weight
            hessian[i, j] = -sp.sqrt(six_masses[i] * six_masses[j]) * cosine_weight
    even = sp.Matrix(
        3,
        3,
        lambda i, j: hessian[2 * i, 2 * j] + hessian[2 * i, 2 * j + 1],
    )
    odd = sp.Matrix(
        3,
        3,
        lambda i, j: hessian[2 * i, 2 * j] - hessian[2 * i, 2 * j + 1],
    )
    require_zero(
        "even rotation kernel",
        even * sp.Matrix([sp.sqrt(a), sp.sqrt(b), sp.sqrt(c)]),
    )
    require_zero("even B diagonal", even[1, 1] + a * t_ab + c * s_bc)

    # A PSD matrix has nonnegative diagonal.  For a,c>0 and
    # t_ab,s_bc>=0, the previous identity forces t_ab=s_bc=0.
    zero_missing = {t_ab: 0, s_bc: 0}
    reflection = {
        v_a: (b * s_ab - c * s_ac) / a,
        v_b: (a * s_ab + c * t_bc) / b,
        v_c: (b * t_bc - a * s_ac) / c,
    }
    for index, (variable, formula) in enumerate(reflection.items()):
        solved = sp.solve(torques[index].subs(zero_missing), variable)
        if len(solved) != 1:
            raise AssertionError(f"torque does not uniquely determine {variable}")
        require_zero(f"torque elimination for {variable}", solved[0] - formula)
    sheet = {**zero_missing, **reflection}
    require_zero("torque vanishes on sheet", torques.subs(sheet))

    expected_even = t_ac * sp.Matrix(
        [
            [c, 0, -sp.sqrt(a * c)],
            [0, 0, 0],
            [-sp.sqrt(a * c), 0, a],
        ]
    )
    expected_odd = t_ac * sp.Matrix(
        [
            [c, 0, sp.sqrt(a * c)],
            [0, 0, 0],
            [sp.sqrt(a * c), 0, a],
        ]
    )
    require_zero("even rank-one formula", even.subs(sheet) - expected_even)
    require_zero("odd rank-one formula", odd.subs(sheet) - expected_odd)
    require_zero(
        "even positive outer product",
        expected_even
        - t_ac
        * sp.Matrix([sp.sqrt(c), 0, -sp.sqrt(a)])
        * sp.Matrix([[sp.sqrt(c), 0, -sp.sqrt(a)]]),
    )
    require_zero(
        "odd positive outer product",
        expected_odd
        - t_ac
        * sp.Matrix([sp.sqrt(c), 0, sp.sqrt(a)])
        * sp.Matrix([[sp.sqrt(c), 0, sp.sqrt(a)]]),
    )
    require_zero(
        "even rank-one identity",
        expected_even**2 - t_ac * (a + c) * expected_even,
    )
    require_zero(
        "odd rank-one identity",
        expected_odd**2 - t_ac * (a + c) * expected_odd,
    )

    sheet_degrees = degrees.subs(sheet)
    expected_degrees = sp.Matrix(
        [
            a + 2 * b * s_ab + c * t_ac,
            b + 2 * a * s_ab + 2 * c * t_bc,
            c + a * t_ac + 2 * b * t_bc,
        ]
    )
    require_zero("sheet degree formulas", sheet_degrees - expected_degrees)
    require_zero(
        "B degree slack identity",
        3 * b - sheet_degrees[1] - 2 * b * (1 - reflection[v_b]),
    )
    require_zero(
        "A degree slack identity",
        3 * (a + c)
        - sheet_degrees[0]
        - 2 * a * (1 - reflection[v_a])
        - 2 * c * (1 - s_ac)
        - c * (1 - t_ac),
    )
    require_zero(
        "C degree slack identity",
        3 * (a + c)
        - sheet_degrees[2]
        - 2 * c * (1 - reflection[v_c])
        - 2 * a * (1 - s_ac)
        - a * (1 - t_ac),
    )
    require_zero(
        "mass-balance upper bound",
        (3 * b + 3 * (a + c)) / 2 - sp.Rational(3, 2) * (a + b + c),
    )

    # Equality in min{3b,3(a+c)}<=3/4 requires b=a+c=1/4.
    # The next identities, with the density box, successively force
    # s_ab=t_bc=1, then s_ac=1, and finally t_ac=1.
    equality_masses = {b: sp.Rational(1, 4), c: sp.Rational(1, 4) - a}
    require_zero(
        "equality forces two cross densities",
        (b * (1 - reflection[v_b]) - a * (1 - s_ab) - c * (1 - t_bc)).subs(
            equality_masses
        ),
    )
    require_zero(
        "equality forces s_ac from v_a box",
        (a * (1 - reflection[v_a]) + c * (1 - s_ac) - (a + c - b * s_ab)).subs(
            equality_masses
        ),
    )
    require_zero(
        "equality forces t_ac from d_A",
        (sp.Rational(3, 4) - sheet_degrees[0] - c * (1 - t_ac))
        .subs(equality_masses)
        .subs({s_ab: 1}),
    )
    binary = {
        **equality_masses,
        s_ab: 1,
        s_ac: 1,
        t_ac: 1,
        t_bc: 1,
    }
    require_zero(
        "binary equality densities",
        sp.Matrix([reflection[v] for v in (v_a, v_b, v_c)]).subs(binary)
        - sp.ones(3, 1),
    )
    require_zero(
        "binary equality degrees",
        sheet_degrees.subs(binary) - sp.Rational(3, 4) * sp.ones(3, 1),
    )
    transverse = sp.Matrix(
        [
            sum(
                six_masses[j] * weights[i, j] * sp.cos(six_phases[i] - six_phases[j])
                for j in range(6)
            )
            for i in (0, 2, 4)
        ]
    )
    require_zero(
        "binary transverse values",
        transverse.subs(sheet).subs(binary) - sp.Rational(1, 4) * sp.ones(3, 1),
    )
    print("All exact fixed-angle identities verified; no nearby-angle claim.")


if __name__ == "__main__":
    main()
