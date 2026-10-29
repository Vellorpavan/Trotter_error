"""
tests/test_hamiltonian.py — StepWise-Q Automated Tests
Verification of Hamiltonian construction, Hermiticity, spectra, and term counting.
"""

from __future__ import annotations

import unittest
import numpy as np
import scipy.linalg as la
from src.problem import tfim_hamiltonian, N_SPINS


class TestHamiltonian(unittest.TestCase):
    """Test suite for TFIM Hamiltonian correctness and physical properties."""

    def test_2qubit_hand_calculated_matrix(self):
        """V1 gate: 2-qubit TFIM Hamiltonian matches analytical 4x4 matrix."""
        J = 1.0
        g = 1.0
        # H = -J Z0 Z1 - g (X0 + X1)
        # Z0 Z1 = diag(1, -1, -1, 1) in standard basis |00>, |01>, |10>, |11>
        # X0 + X1 has off-diagonal 1s between states differing by 1 bit flip
        H2 = tfim_hamiltonian(n_spins=2, J=J, g=g, boundary="open")
        mat = H2.to_matrix()

        # Hand-calculated analytical matrix
        expected = np.array([
            [-J, -g, -g,  0],
            [-g,  J,  0, -g],
            [-g,  0,  J, -g],
            [ 0, -g, -g, -J]
        ], dtype=complex)

        np.testing.assert_allclose(mat, expected, atol=1e-12,
                                   err_msg="2-qubit Hamiltonian matrix does not match analytical matrix.")

    def test_4qubit_hermiticity(self):
        """Verify H == H_dagger across all supported boundary conditions."""
        for bnd in ["open", "ring", "star", "all_to_all"]:
            for g_val in [0.5, 1.0, 1.5]:
                H = tfim_hamiltonian(n_spins=4, J=1.0, g=g_val, boundary=bnd)
                mat = H.to_matrix()
                herm_diff = np.max(np.abs(mat - mat.conj().T))
                self.assertLess(herm_diff, 1e-12, f"Hamiltonian {bnd} at g={g_val} is not Hermitian.")

    def test_pauli_term_counts(self):
        """Verify Pauli operator term count for all graph boundaries."""
        H_open = tfim_hamiltonian(4, 1.0, 1.0, "open")
        H_ring = tfim_hamiltonian(4, 1.0, 1.0, "ring")
        H_star = tfim_hamiltonian(4, 1.0, 1.0, "star")
        H_k4   = tfim_hamiltonian(4, 1.0, 1.0, "all_to_all")

        self.assertEqual(len(H_open), 7, "Open chain must have 3 ZZ + 4 X = 7 terms.")
        self.assertEqual(len(H_ring), 8, "Periodic ring must have 4 ZZ + 4 X = 8 terms.")
        self.assertEqual(len(H_star), 7, "Star graph must have 3 ZZ + 4 X = 7 terms.")
        self.assertEqual(len(H_k4), 10, "K4 complete graph must have 6 ZZ + 4 X = 10 terms.")

    def test_spectrum_ordering(self):
        """Verify ground state energy properties as transverse field g varies."""
        energies_gs = []
        for g_val in [0.5, 1.0, 1.5]:
            H = tfim_hamiltonian(4, 1.0, g_val, "open")
            mat = H.to_matrix()
            w = la.eigvalsh(mat)
            energies_gs.append(w[0])

        # As g increases, transverse kinetic energy lowers ground state energy
        self.assertTrue(energies_gs[0] > energies_gs[1] > energies_gs[2],
                        "Ground state energy must decrease monotonically with increasing g.")


if __name__ == "__main__":
    unittest.main()

