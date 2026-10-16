"""Regression tests for the fourth reflected-pair insertion calculation."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import certify_pair_insertion_box as certificate  # noqa: E402
import test_pair_insertion as model  # noqa: E402


class PairInsertionTests(unittest.TestCase):
    def test_direct_equations_match_first_derivatives(self):
        record = model.solve_at_phase(3 * math.pi / 4, 0, 100, False, 0)
        self.assertEqual(record["status"], "optimal")
        audit = model.audit_finite_difference(record, 1e-8)
        self.assertLess(audit["degree_error"], 2e-8)
        self.assertLess(audit["old_torque_error"], 2e-8)
        self.assertLess(audit["new_torque_error"], 1e-12)
        self.assertLess(audit["hessian_error"], 3e-8)
        self.assertLess(audit["rotation_error"], 1e-12)

    def test_old_weak_only_has_spurious_improvement(self):
        record = model.solve_at_phase(math.pi, 1e-7, 10, True, 0)
        self.assertGreater(record["degree_gain_derivative"], 3)
        self.assertGreater(min(record["old_weak_opening_eigenvalues"]), 0)
        self.assertLess(min(record["weak_opening_eigenvalues"]), -0.3)
        actual = model.audit_finite_difference(record, 1e-6)
        self.assertLess(actual["lowest_nonrotation_eigenvalue"], -3e-7)

    def test_unbounded_dual_excludes_interval(self):
        selected = model.solve_at_phase(1.57, 0, 100, False, 0, True)
        upper = certificate.bound_box("1.55", "1.59", 100, selected, True)
        self.assertLess(upper, -1)


if __name__ == "__main__":
    unittest.main()
