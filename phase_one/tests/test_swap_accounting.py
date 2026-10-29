"""
tests/test_swap_accounting.py — StepWise-Q Automated Tests
Verification of SWAP overhead, topological advantages, and parametrized CX algebraic identities.
"""

from __future__ import annotations

import unittest
from src.problem import tfim_hamiltonian
from src.processors import PROCESSOR_REGISTRY
from src.routing import build_trotter_circuit, transpile_circuit, count_swaps, circuit_stats


class TestSwapAccounting(unittest.TestCase):
    """Test suite for SWAP routing overhead and algebraic CX identities."""

    def test_open_chain_zero_swaps(self):
        """V5 gate: Open chain requires exactly 0 SWAPs across all processors."""
        H_open = tfim_hamiltonian(4, 1.0, 1.0, "open")
        qc_open = build_trotter_circuit(H_open, 1.0, 2, order=1)

        for pid in ["A0", "A1", "A_star", "B0", "B_star"]:
            proc = PROCESSOR_REGISTRY[pid]
            swaps = count_swaps(qc_open, proc)
            self.assertEqual(swaps, 0, f"Open chain on processor {pid} should have 0 SWAPs, got {swaps}.")

    def test_ring_topology_swap_advantage(self):
        """V5 gate: Ring requires >=2 SWAPs on linear A0, but 0 SWAPs on cycle-matched A* and B*."""
        H_ring = tfim_hamiltonian(4, 1.0, 1.0, "ring")
        qc_ring = build_trotter_circuit(H_ring, 1.0, 1, order=1)

        sw_a0 = count_swaps(qc_ring, PROCESSOR_REGISTRY["A0"])
        sw_astar = count_swaps(qc_ring, PROCESSOR_REGISTRY["A_star"])
        sw_bstar = count_swaps(qc_ring, PROCESSOR_REGISTRY["B_star"])

        self.assertGreaterEqual(sw_a0, 2, "Linear chain A0 must require >=2 SWAPs per step for ring.")
        self.assertEqual(sw_astar, 0, "Cycle-matched A* must require 0 SWAPs for ring.")
        self.assertEqual(sw_bstar, 0, "Heavy-hex / 4-cycle B* must require 0 SWAPs for ring.")

    def test_star_topology_swap_advantage(self):
        """Star graph topology: Tree A1 requires 0 SWAPs, while linear A0 requires >=1 SWAP."""
        H_star = tfim_hamiltonian(4, 1.0, 1.0, "star")
        qc_star = build_trotter_circuit(H_star, 1.0, 1, order=1)

        sw_a0 = count_swaps(qc_star, PROCESSOR_REGISTRY["A0"])
        sw_a1 = count_swaps(qc_star, PROCESSOR_REGISTRY["A1"])

        self.assertGreaterEqual(sw_a0, 1, "Linear chain A0 must require >= 1 SWAP for star topology.")
        self.assertEqual(sw_a1, 0, "Tree topology A1 must natively support star topology with 0 SWAPs.")

    def test_cx_algebraic_identity_across_processors(self):
        """V6 gate: cx_transpiled == cx_logical + 3 * swaps (modulo peephole cancellations <= 2)."""
        configs_tested = 0
        for pid in ["A0", "A1", "A_star", "B0", "B_star"]:
            proc = PROCESSOR_REGISTRY[pid]
            for boundary in ["open", "ring"]:
                H = tfim_hamiltonian(4, 1.0, 1.0, boundary)
                n_terms = len(H) - 4  # ZZ terms count
                for n_steps in [1, 2, 4]:
                    qc = build_trotter_circuit(H, 1.0, n_steps, order=1)
                    qc_t = transpile_circuit(qc, proc)
                    cx_actual = circuit_stats(qc_t)["cx"]
                    swaps = count_swaps(qc, proc)

                    cx_logical = 2 * n_terms * n_steps
                    cx_expected = cx_logical + 3 * swaps
                    diff = abs(cx_actual - cx_expected)
                    # Qiskit opt_level=1 can cancel up to 2 CX across adjacent decomposed SWAPs
                    self.assertLessEqual(diff, 2,
                                         f"Identity failed for {pid} {boundary} n={n_steps}: "
                                         f"actual={cx_actual}, expected={cx_expected}, diff={diff}")
                    configs_tested += 1

        self.assertEqual(configs_tested, 30, "Should have tested exactly 30 configurations.")


if __name__ == "__main__":
    unittest.main()

