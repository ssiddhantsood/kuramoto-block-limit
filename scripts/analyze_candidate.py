#!/usr/bin/env python3
"""Dependency-free reconstruction of the six-block limiting candidate."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def jacobi_eigenvalues(matrix: list[list[float]], tolerance: float = 1e-15) -> list[float]:
    """Return eigenvalues of a small real symmetric matrix by Jacobi rotations."""
    a = [row[:] for row in matrix]
    n = len(a)
    for _ in range(100 * n * n):
        p, q = max(
            ((i, j) for i in range(n) for j in range(i + 1, n)),
            key=lambda ij: abs(a[ij[0]][ij[1]]),
        )
        if abs(a[p][q]) < tolerance:
            break
        tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
        t = math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau)) if tau else 1.0
        cosine = 1.0 / math.sqrt(1.0 + t * t)
        sine = t * cosine
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        a[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq
        a[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq
        a[p][q] = a[q][p] = 0.0
        for k in range(n):
            if k in (p, q):
                continue
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = cosine * akp - sine * akq
            a[k][q] = a[q][k] = sine * akp + cosine * akq
    return sorted(a[i][i] for i in range(n))


def block_matrix(x: float, y: float) -> list[list[float]]:
    return [
        [1.0, x, y, 0.0, 1.0, 1.0],
        [x, 1.0, 0.0, y, 1.0, 1.0],
        [y, 0.0, 1.0, 1.0, 0.0, 1.0],
        [0.0, y, 1.0, 1.0, 1.0, 0.0],
        [1.0, 1.0, 0.0, 1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 1.0],
    ]


def linear_boundary_extrapolation() -> float:
    sequence = json.loads((ROOT / "data/gap_sequence.json").read_text())
    first, second = sequence[-2:]
    x1, y1 = float(first["gap"]), float(first["mu"])
    x2, y2 = float(second["gap"]), float(second["mu"])
    slope = (y2 - y1) / (x2 - x1)
    return y2 - slope * x2


def main() -> None:
    candidate = json.loads((ROOT / "data/near_boundary_candidate.json").read_text())
    a, b, c = map(float, candidate["pair_class_masses"])
    x = float(candidate["partial_block_densities"]["x"])
    y = float(candidate["partial_block_densities"]["y"])
    alpha, beta, gamma = map(float, candidate["angles"])

    masses = [a, a, b, b, c, c]
    phases = [alpha, -alpha, beta, -beta, gamma, -gamma]
    weights = block_matrix(x, y)

    degrees = [sum(masses[j] * weights[i][j] for j in range(6)) for i in range(6)]
    torques = [
        sum(masses[j] * weights[i][j] * math.sin(phases[j] - phases[i]) for j in range(6))
        for i in range(6)
    ]

    hessian = [[0.0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            cosine_weight = weights[i][j] * math.cos(phases[i] - phases[j])
            hessian[i][i] += masses[j] * cosine_weight
            hessian[i][j] = -math.sqrt(masses[i] * masses[j]) * cosine_weight

    quotient_eigenvalues = jacobi_eigenvalues(hessian)
    transverse = [
        sum(
            masses[j] * weights[i][j] * math.cos(phases[i] - phases[j])
            for j in range(6)
        )
        for i in range(6)
    ]

    result = {
        "mass_sum": sum(masses),
        "distinct_degrees": [degrees[0], degrees[2], degrees[4]],
        "minimum_degree": min(degrees),
        "degree_identity_half_plus_a": 0.5 + a,
        "maximum_absolute_torque": max(map(abs, torques)),
        "quotient_hessian_eigenvalues": quotient_eigenvalues,
        "smallest_nonrotation_quotient_eigenvalue": quotient_eigenvalues[1],
        "distinct_transverse_multiplication_terms": [transverse[0], transverse[2], transverse[4]],
        "two_point_boundary_extrapolation": linear_boundary_extrapolation(),
    }

    assert abs(result["mass_sum"] - 1.0) < 1e-14
    assert result["maximum_absolute_torque"] < 1e-13
    assert min(transverse) > 0.0
    assert abs(result["minimum_degree"] - (0.5 + a)) < 1e-12
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

