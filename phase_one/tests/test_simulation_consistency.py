"""
tests/test_simulation_consistency.py — StepWise-Q Automated Tests
Verification of DensityMatrix and Shot simulation consistency, zero-noise equivalence, and negative control.
"""

from __future__ import annotations

import unittest
from src.problem import tfim_hamiltonian
from src.processors import PROCESSOR_REGISTRY
from src.routing import build_trotter_circuit, transpile_circuit, simulate_density_matrix, simulate_shots


class TestSimulationConsistency(unittest.TestCase):
    """Test suite for simulation engines: density matrix vs finite-shot agreement."""

    def test_density_matrix_vs_shot_within_3se(self):
        """V8 gate: Density matrix expectation agrees with 4096-shot simulation within 3 SE."""
        proc = PROCESSOR_REGISTRY["A_star"]
        for g_val in [0.5, 1.0, 1.5]:
            H = tfim_hamiltonian(4, J=1.0, g=g_val, boundary="open")
            qc = build_trotter_circuit(H, 1.0, 2, order=1)
            qc_t = transpile_circuit(qc, proc)

            dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
            shot_res = simulate_shots(qc_t, proc, ["Z0"], shots=4096)

            val_dm = dm_res["Z0"]
            val_shot, se_shot = shot_res["Z0"]
            diff = abs(val_dm - val_shot)

            self.assertLessEqual(diff, 3.0 * se_shot + 1e-4,
                                 f"DM ({val_dm:.4f}) and shot ({val_shot:.4f}±{se_shot:.4f}) exceed 3 SE at g={g_val}")

    def test_negative_control_fails_3se(self):
        """V8 negative control: An intentionally corrupted expectation value fails the 3 SE check."""
        proc = PROCESSOR_REGISTRY["A_star"]
        H = tfim_hamiltonian(4, J=1.0, g=1.0, boundary="open")
        qc = build_trotter_circuit(H, 1.0, 2, order=1)
        qc_t = transpile_circuit(qc, proc)

        dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
        shot_res = simulate_shots(qc_t, proc, ["Z0"], shots=4096)

        val_corrupted = dm_res["Z0"] + 0.40  # Massive corruption
        val_shot, se_shot = shot_res["Z0"]

        diff = abs(val_corrupted - val_shot)
        self.assertGreater(diff, 3.0 * se_shot,
                           "Negative control failed: corrupted expectation should exceed 3 SE threshold.")

    def test_zero_noise_equals_noiseless(self):
        """V4 gate: Density matrix simulation with zero noise matches noiseless Trotter to machine precision."""
        from qiskit_aer import AerSimulator
        from qiskit.quantum_info import DensityMatrix as QiDM
        from src.problem import compute_observable

        H = tfim_hamiltonian(4, 1.0, 1.0, "open")
        qc = build_trotter_circuit(H, 2.0, 2, order=1)
        proc = PROCESSOR_REGISTRY["A0"]
        qc_t = transpile_circuit(qc, proc)

        # Noiseless via pure simulator
        sim_pure = AerSimulator(method="density_matrix")
        qc_pure = qc_t.copy()
        qc_pure.save_density_matrix()
        res_pure = sim_pure.run(qc_pure, shots=1, seed_simulator=42).result()
        dm_pure = QiDM(res_pure.data(0)["density_matrix"])
        val_pure = compute_observable(dm_pure, "Z0", 4)

        # Zero-noise via simulate_density_matrix with zero noise rates
        proc_zero = {**proc, "noise": {**proc["noise"], "p_1q": 0.0, "p_cx": 0.0, "p_readout": 0.0}}
        res_zero = simulate_density_matrix(qc_t, proc_zero, ["Z0"])
        val_zero = res_zero["Z0"]

        self.assertAlmostEqual(val_pure, val_zero, places=10,
                               msg="Zero-noise simulation does not match noiseless simulation to machine precision.")


if __name__ == "__main__":
    unittest.main()
