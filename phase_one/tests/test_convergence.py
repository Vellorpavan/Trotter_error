"""
tests/test_convergence.py — StepWise-Q Automated Tests
Verification of noiseless Trotter error scaling (Lie-Trotter O(1/n) and Suzuki O(1/n^2)).
"""

from __future__ import annotations

import unittest
import numpy as np
import scipy.linalg as la
from qiskit.quantum_info import Operator
from src.problem import tfim_hamiltonian
from src.routing import build_trotter_circuit


class TestConvergence(unittest.TestCase):
    """Test suite for noiseless Trotter convergence scaling."""

    def setUp(self):
        self.H = tfim_hamiltonian(4, J=1.0, g=1.0, boundary="open")
        self.H_mat = self.H.to_matrix()
        self.t = 0.5
        self.U_exact = la.expm(-1j * self.H_mat * self.t)
        self.steps = [2, 4, 8, 16]

    def test_order1_convergence_rate(self):
        """V3 gate: First-order Lie-Trotter operator norm error scales as O(1/n) (slope ≈ -1.0)."""
        errors = []
        for n in self.steps:
            qc = build_trotter_circuit(self.H, self.t, n, order=1).decompose()
            u_trotter = Operator(qc).data
            err = la.norm(u_trotter - self.U_exact, 2)
            errors.append(err)

        slope = float(np.polyfit(np.log(self.steps), np.log(errors), 1)[0])
        # Expected -1.0 with standard tolerance ±0.35
        self.assertAlmostEqual(slope, -1.0, delta=0.35,
                               msg=f"Order 1 convergence slope {slope:.3f} deviates from expected -1.0")

    def test_order2_convergence_rate(self):
        """V3 gate: Second-order Suzuki-Trotter operator norm error scales as O(1/n^2) (slope ≈ -2.0)."""
        errors = []
        for n in self.steps:
            qc = build_trotter_circuit(self.H, self.t, n, order=2).decompose()
            u_trotter = Operator(qc).data
            err = la.norm(u_trotter - self.U_exact, 2)
            errors.append(err)

        slope = float(np.polyfit(np.log(self.steps), np.log(errors), 1)[0])
        # Expected -2.0 with standard tolerance ±0.60
        self.assertAlmostEqual(slope, -2.0, delta=0.60,
                               msg=f"Order 2 convergence slope {slope:.3f} deviates from expected -2.0")

    def test_order2_beats_order1_in_asymptotic_regime(self):
        """In the asymptotic fine-step regime (n=16), Order 2 error is strictly smaller than Order 1."""
        qc1 = build_trotter_circuit(self.H, self.t, 16, order=1).decompose()
        qc2 = build_trotter_circuit(self.H, self.t, 16, order=2).decompose()
        err1 = la.norm(Operator(qc1).data - self.U_exact, 2)
        err2 = la.norm(Operator(qc2).data - self.U_exact, 2)
        self.assertLess(err2, err1, f"At n=16, Order 2 error ({err2:.2e}) should be < Order 1 ({err1:.2e})")


if __name__ == "__main__":
    unittest.main()

