#!/usr/bin/env python3
"""Exact degree arithmetic for the clique blow-ups of the N=80,002 witness."""

from fractions import Fraction
import json


N = 80_002
DELTA = 55_018


def record(multiplier: int) -> dict[str, object]:
    if multiplier < 1:
        raise ValueError("multiplier must be positive")
    vertices = multiplier * N
    minimum_degree = multiplier * (DELTA + 1) - 1
    ratio = Fraction(minimum_degree, vertices - 1)
    excess = 16 * minimum_degree - 11 * (vertices - 1)
    return {
        "multiplier": multiplier,
        "vertices": vertices,
        "minimum_degree": minimum_degree,
        "ratio": str(ratio),
        "ratio_decimal": float(ratio),
        "exact_excess_over_11_16": excess,
    }


def main() -> None:
    examples = [record(m) for m in (1, 2, 10, 1000)]
    limit = Fraction(DELTA + 1, N)
    assert all(item["exact_excess_over_11_16"] > 0 for item in examples)
    assert all(item["exact_excess_over_11_16"] == 282 * item["multiplier"] - 5 for item in examples)
    print(
        json.dumps(
            {
                "family": examples,
                "limiting_ratio": str(limit),
                "limiting_ratio_decimal": float(limit),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

