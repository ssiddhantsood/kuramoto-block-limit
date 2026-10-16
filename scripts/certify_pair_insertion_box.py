#!/usr/bin/env python3
"""Interval-dual exclusion for one phase box in the fourth-pair insertion test.

The numerical SDP only selects multipliers.  The reported upper bound is
recomputed from exact-decimal scalar multipliers, an exact-decimal Gram factor
for the PSD dual, an interval enclosure of the KKT root, and interval sine /
cosine on the entire requested phase box.  A negative bound excludes a
positive first-order degree gain on that box, with the derivative cap shown.
It does not exclude higher-order insertions or other graphon architectures.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

import test_pair_insertion as model


def iv(value):
    return mp.iv.mpf(str(value))


def iv_upper(value):
    return mp.mpf(value._mpi_[1])


def iv_lower(value):
    return mp.mpf(value._mpi_[0])


def block_matrix(q):
    one, zero = iv(1), iv(0)
    w = [[one if i == j else zero for j in range(6)] for i in range(6)]
    for pair in range(3):
        i = 2 * pair
        w[i][i + 1] = w[i + 1][i] = q[pair]
    for left, right, same, opposite in ((0, 1, 3, 4), (0, 2, 5, 6), (1, 2, 7, 8)):
        lp, lm, rp, rm = 2 * left, 2 * left + 1, 2 * right, 2 * right + 1
        w[lp][rp] = w[rp][lp] = w[lm][rm] = w[rm][lm] = q[same]
        w[lp][rm] = w[rm][lp] = w[lm][rp] = w[rp][lm] = q[opposite]
    return w


def root_data():
    saved = json.loads((model.ROOT / "data/limiting_kkt_certificate.json").read_text())
    if not saved["krawczyk"]["unique_kkt_root_in_box"]:
        raise RuntimeError("the endpoint KKT root is not certified")
    radius = iv(saved["krawczyk"]["radius"])
    root = [iv(x) + mp.iv.mpf([-iv_upper(radius), iv_upper(radius)]) for x in saved["variables"]]
    a, b, c, x, y, alpha, beta, gamma, mu = root
    p = [a, b, c]
    q = [x, iv(1), iv(0), y, iv(0), iv(1), iv(1), iv(0), iv(1)]
    theta = [alpha, beta, gamma]
    phases = [alpha, -alpha, beta, -beta, gamma, -gamma]
    masses = [a, a, b, b, c, c]
    return p, q, theta, phases, masses, block_matrix(q), mu


def weak_basis(masses, phases, weights):
    z = iv(0)
    lap = [[z for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            c = masses[i] * masses[j] * weights[i][j] * mp.iv.cos(phases[i] - phases[j])
            lap[i][i] += c
            lap[j][j] += c
            lap[i][j] -= c
            lap[j][i] -= c
    even = [[lap[2 * i][2 * j] + lap[2 * i][2 * j + 1] for j in range(3)] for i in range(3)]
    odd = [[lap[2 * i][2 * j] - lap[2 * i][2 * j + 1] for j in range(3)] for i in range(3)]
    e = [even[0][2] - even[0][1], even[0][0] - even[0][2], even[0][1] - even[0][0]]
    o = [odd[0][1] * odd[1][2] - odd[0][2] * odd[1][1],
         odd[0][2] * odd[1][0] - odd[0][0] * odd[1][2],
         odd[0][0] * odd[1][1] - odd[0][1] * odd[1][0]]
    pair_mass = [masses[0], masses[2], masses[4]]
    weighted_mean = sum(pair_mass[i] * e[i] for i in range(3)) / sum(pair_mass)
    e = [value - weighted_mean for value in e]
    en = mp.iv.sqrt(2 * sum(value * value for value in e))
    on = mp.iv.sqrt(2 * sum(value * value for value in o))
    e, o = [value / en for value in e], [value / on for value in o]
    basis = [[z for _ in range(4)] for _ in range(8)]
    for i in range(6):
        basis[i][0] = e[i // 2]
        basis[i][1] = o[i // 2] * (1 if i % 2 == 0 else -1)
    basis[6][2], basis[7][3] = iv(1), iv(1)
    return basis


def exact_decimal_dual(numerical, base, basis, unbounded=False):
    raw = numerical["dual_multipliers"]
    # All scalar inequality multipliers are clipped to exact nonnegative
    # decimals.  The two active-degree weights must sum to one exactly, so
    # derive the second from the first rather than separately rounding it.
    scalar = lambda x: iv(format(max(0.0, float(x)), ".16g"))
    root = json.loads((model.ROOT / "data/limiting_kkt_certificate.json").read_text())
    if not root["krawczyk"]["unique_kkt_root_in_box"]:
        raise RuntimeError("the endpoint KKT root is not certified")
    radius = iv(root["krawczyk"]["radius"])
    kappa = [iv(x) + mp.iv.mpf([-iv_upper(radius), iv_upper(radius)]) for x in root["multipliers"]]
    first = -kappa[4] if unbounded else scalar(raw["active_degree"][0])
    active = [first, iv(1) - first]
    if iv_lower(active[1]) < 0:
        raise RuntimeError("rounded active-degree weight is negative")
    factor = [[iv(format(float(x), ".16g")) for x in row] for row in raw["psd_factor"]]
    if unbounded:
        _, _, _, phases, masses, weights, _ = base
        H = [[iv(0) for _ in range(6)] for _ in range(6)]
        for i in range(6):
            for j in range(6):
                if i != j:
                    coefficient = masses[j] * weights[i][j] * mp.iv.cos(phases[i] - phases[j])
                    H[i][i] += coefficient
                    H[i][j] -= coefficient
        even_trace = sum(H[2 * i][2 * i] + H[2 * i][2 * i + 1] for i in range(3))
        odd = [[H[2 * i][2 * j] - H[2 * i][2 * j + 1] for j in range(3)] for i in range(3)]
        odd_product = sum(odd[i][i] * odd[j][j] - odd[i][j] * odd[j][i]
                          for i, j in ((0, 1), (0, 2), (1, 2)))
        weighted_norm_even = sum(masses[i] * basis[i][0] ** 2 for i in range(6))
        weighted_norm_odd = sum(masses[i] * basis[i][1] ** 2 for i in range(6))
        target_diagonal = [
            -kappa[6] * even_trace / weighted_norm_even,
            -kappa[7] * odd_product / weighted_norm_odd,
        ]
        for index in (0, 1):
            present = sum(value * value for value in factor[index])
            scale = mp.iv.sqrt(target_diagonal[index] / present)
            factor[index] = [value * scale for value in factor[index]]
    Y = [[sum(factor[i][k] * factor[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    return {
        "mass": -kappa[0] if unbounded else iv(format(raw["mass"], ".16g")),
        "old_torque": [-kappa[i + 1] for i in range(3)] if unbounded else [iv(format(x, ".16g")) for x in raw["old_torque"]],
        "new_torque": iv(format(raw["new_torque"], ".16g")),
        "new_degree": scalar(raw["new_degree"]),
        "new_transverse": scalar(raw["new_transverse"]),
        "active_degree": active,
        "Y": Y,
    }


def lagrangian(x, delta, base, basis, dual):
    pair_mass, q, theta, phases, masses, weights, mu = base
    u, dt, dq, same, opposite = x[:3], x[3:6], x[6:15], x[15:18], x[18:21]
    dm = [u[i // 2] for i in range(6)]
    dphase = [dt[i // 2] * (1 if i % 2 == 0 else -1) for i in range(6)]
    dw = block_matrix(dq)
    # A density derivative has zero diagonal, unlike the block matrix itself.
    for i in range(6):
        dw[i][i] = iv(0)
    new_phase = [delta, -delta]

    def nw(i, j):
        return same[i // 2] if i % 2 == j - 6 else opposite[i // 2]

    old_torque, old_degree = [], []
    for i in (0, 2, 4):
        torque = iv(0)
        degree = iv(0)
        for j in range(6):
            phi = phases[j] - phases[i]
            torque += dm[j] * weights[i][j] * mp.iv.sin(phi)
            torque += masses[j] * dw[i][j] * mp.iv.sin(phi)
            torque += masses[j] * weights[i][j] * mp.iv.cos(phi) * (dphase[j] - dphase[i])
            degree += dm[j] * weights[i][j] + masses[j] * dw[i][j]
        for j in (6, 7):
            torque += nw(i, j) * mp.iv.sin(new_phase[j - 6] - phases[i])
            degree += nw(i, j)
        old_torque.append(torque)
        old_degree.append(degree)
    new_torque = sum(masses[j] * nw(j, 6) * mp.iv.sin(phases[j] - delta) for j in range(6))
    new_degree = sum(masses[j] * nw(j, 6) for j in range(6))
    new_transverse = sum(masses[j] * nw(j, 6) * mp.iv.cos(phases[j] - delta) for j in range(6))

    def quadratic_on_basis(i, j):
        v = [basis[i][k] - basis[j][k] for k in range(4)]
        return sum(v[k] * dual["Y"][k][l] * v[l] for k in range(4) for l in range(4))

    psd_term = iv(0)
    for i in range(6):
        for j in range(i + 1, 6):
            phi = phases[i] - phases[j]
            coeff = (dm[i] * masses[j] + masses[i] * dm[j]) * weights[i][j] * mp.iv.cos(phi)
            coeff += masses[i] * masses[j] * dw[i][j] * mp.iv.cos(phi)
            coeff -= masses[i] * masses[j] * weights[i][j] * mp.iv.sin(phi) * (dphase[i] - dphase[j])
            psd_term += coeff * quadratic_on_basis(i, j)
    for i in range(6):
        for j in (6, 7):
            coeff = masses[i] * nw(i, j) * mp.iv.cos(phases[i] - new_phase[j - 6])
            psd_term += coeff * quadratic_on_basis(i, j)

    return (
        dual["mass"] * (sum(u) + 1)
        + sum(dual["old_torque"][k] * old_torque[k] for k in range(3))
        + dual["new_torque"] * new_torque
        + dual["new_degree"] * (new_degree - mu)
        + dual["new_transverse"] * new_transverse
        + dual["active_degree"][0] * old_degree[0]
        + dual["active_degree"][1] * old_degree[2]
        + psd_term
    )


def bound_box(lower, upper, derivative_cap, dual_record, unbounded=False):
    mp.mp.dps = 80
    mp.iv.dps = 80
    delta = mp.iv.mpf([lower, upper])
    base = root_data()
    basis = weak_basis(base[4], base[3], base[5])
    dual = exact_decimal_dual(dual_record, base, basis, unbounded)
    zero = [iv(0) for _ in range(21)]
    constant = lagrangian(zero, delta, base, basis, dual)
    total = constant
    for index in range(21):
        unit = zero.copy()
        unit[index] = iv(1)
        coeff = lagrangian(unit, delta, base, basis, dual) - constant
        if unbounded and index < 15:
            # For the eight free old coordinates, exact stationarity follows
            # from the certified KKT equations and the rank-one/rank-two
            # eigenvalue-derivative identity.  Interval arithmetic checks it
            # is consistent, but cannot symbolically cancel dependent balls.
            free = index < 6 or index in (6, 9)
            if free:
                if not (iv_lower(coeff) <= 0 <= iv_upper(coeff)):
                    raise RuntimeError(f"KKT cancellation failed at coordinate {index}")
            else:
                density = base[1][index - 6]
                if iv_upper(density) < mp.mpf("1e-12"):
                    if iv_upper(coeff) > 0:
                        raise RuntimeError(f"boundary density sign failed at coordinate {index}")
                elif iv_lower(density) > 1 - mp.mpf("1e-12"):
                    if iv_lower(coeff) < 0:
                        raise RuntimeError(f"boundary density sign failed at coordinate {index}")
                else:
                    raise RuntimeError("unexpected interior density")
            continue
        if index < 6:
            lo, hi = -derivative_cap, derivative_cap
        elif index < 15:
            density = (base[1])[index - 6]
            if iv_upper(density) < mp.mpf("1e-12"):
                lo, hi = 0, derivative_cap
            elif iv_lower(density) > 1 - mp.mpf("1e-12"):
                lo, hi = -derivative_cap, 0
            else:
                lo, hi = -derivative_cap, derivative_cap
        else:
            lo, hi = 0, 1
        support_upper = max(iv_upper(coeff * iv(lo)), iv_upper(coeff * iv(hi)))
        total += mp.iv.mpf([support_upper, support_upper])
    return iv_upper(total)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lower", default="1.5707963267948966")
    parser.add_argument("--upper", default="1.5707963267948967")
    parser.add_argument("--derivative-cap", type=float, default=100)
    parser.add_argument("--unbounded-derivatives", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    center = (float(args.lower) + float(args.upper)) / 2
    selected = model.solve_at_phase(center, 0, args.derivative_cap, False, 0, True)
    if selected.get("status") not in ("optimal", "optimal_inaccurate"):
        raise RuntimeError(f"numerical dual selection failed: {selected['status']}")
    upper_bound = bound_box(args.lower, args.upper, args.derivative_cap, selected, args.unbounded_derivatives)
    report = {
        "phase_box_radians": [args.lower, args.upper],
        "derivative_cap": args.derivative_cap,
        "unbounded_derivatives": args.unbounded_derivatives,
        "dual_selected_at": center,
        "numerical_primal_gain": selected["degree_gain_derivative"],
        "interval_dual_upper": mp.nstr(upper_bound, 30),
        "positive_first_order_gain_excluded": upper_bound < 0,
        "scope": ("Only reflected fourth-pair insertion at the certified six-class endpoint and stated phase interval; "
                  + ("derivatives are unbounded by exact KKT cancellation; " if args.unbounded_derivatives else "the derivative cap is required; ")
                  + "not a global optimality claim."),
    }
    rendered = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
