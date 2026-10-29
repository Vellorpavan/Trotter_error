"""
tests/test_exact_reference.py — StepWise-Q Automated Tests
Verification of scipy.linalg.expm exact reference, unitarity, and observables.
"""

from __future__ import annotations

import unittest
import numpy as np
import scipy.linalg as la
from src.problem import tfim_hamiltonian, initial_state, exact_trajectory, compute_observable


class TestExactReference(unittest.TestCase):
    """Test suite for exact classical reference simulation."""

    def test_unitarity(self):
        """Verify U(t) = expm(-i H t) is unitary: U^dag U == I."""
        for g_val in [0.5, 1.0, 1.5]:
            H = tfim_hamiltonian(4, J=1.0, g=g_val, boundary="open")
            H_mat = H.to_matrix()
            for t_val in [0.5, 1.0, 2.0]:
                U = la.expm(-1j * H_mat * t_val)
                diff_u = np.max(np.abs(U.conj().T @ U - np.eye(16)))
                self.assertLess(diff_u, 1e-12, f"U(t={t_val}) at g={g_val} is not unitary.")

    def test_t0_observable(self):
        """V2 gate: <Z0> expectation value at t=0 matches state definition."""
        H = tfim_hamiltonian(4, J=1.0, g=1.0, boundary="open")
        sv_zeros = initial_state("zeros", 4)
        traj_zeros = exact_trajectory(H, sv_zeros, np.array([0.0]), ["Z0"])
        self.assertAlmostEqual(traj_zeros["Z0"][0], 1.0, places=9,
                               msg="<Z0> at t=0 for |0000> must be +1.0.")

    def test_energy_conservation(self):
        """Verify exact time evolution conserves energy <psi(t)|H|psi(t)>."""
        H = tfim_hamiltonian(4, J=1.0, g=1.0, boundary="open")
        H_mat = H.to_matrix()
        sv0 = initial_state("zeros", 4).data

        E0 = float(np.real(sv0.conj().T @ H_mat @ sv0))
        for t_val in [0.5, 1.0, 1.5, 2.0]:
            U = la.expm(-1j * H_mat * t_val)
            psi_t = U @ sv0
            E_t = float(np.real(psi_t.conj().T @ H_mat @ psi_t))
            self.assertAlmostEqual(E_t, E0, places=8,
                                   msg=f"Energy not conserved at t={t_val}: E(t)={E_t}, E0={E0}")


if __name__ == "__main__":
    unittest.main()

