"""
tests/test_transpilation.py — StepWise-Q Automated Tests
Verification of transpilation to basis gates, optimization levels, and circuit depth scaling.
"""

from __future__ import annotations

import unittest
from src.problem import tfim_hamiltonian
from src.processors import PROCESSOR_REGISTRY
from src.routing import build_trotter_circuit, transpile_circuit, circuit_stats


class TestTranspilation(unittest.TestCase):
    """Test suite for transpiler compliance and gate statistics."""

    def test_basis_gates_enforcement(self):
        """Transpiled circuits must only contain basis gates: rz, sx, x, cx."""
        H = tfim_hamiltonian(4, 1.0, 1.0, "open")
        allowed_gates = {"rz", "sx", "x", "cx", "barrier"}

        for pid in ["A0", "A1", "A_star", "B0", "B_star"]:
            proc = PROCESSOR_REGISTRY[pid]
            qc = build_trotter_circuit(H, 1.0, 2, order=1)
            qc_t = transpile_circuit(qc, proc)

            gate_ops = set(qc_t.count_ops().keys())
            unallowed = gate_ops - allowed_gates
            self.assertEqual(len(unallowed), 0,
                             f"Processor {pid} transpiled circuit contains unauthorized gates: {unallowed}")

    def test_linear_scaling_of_cx_with_steps(self):
        """Logical Trotter circuits without SWAP overhead scale CX count linearly with n."""
        H = tfim_hamiltonian(4, 1.0, 1.0, "open")
        proc_astar = PROCESSOR_REGISTRY["A_star"]

        cx_counts = {}
        for n in [1, 2, 4]:
            qc = build_trotter_circuit(H, 1.0, n, order=1)
            qc_t = transpile_circuit(qc, proc_astar)
            cx_counts[n] = circuit_stats(qc_t)["cx"]

        # Open chain on A* has 6 CX per step
        self.assertEqual(cx_counts[1], 6)
        self.assertEqual(cx_counts[2], 12)
        self.assertEqual(cx_counts[4], 24)

    def test_order2_depth_overhead(self):
        """Second-order Suzuki-Trotter circuit depth is strictly greater than first-order for same n."""
        H = tfim_hamiltonian(4, 1.0, 1.0, "open")
        proc_astar = PROCESSOR_REGISTRY["A_star"]

        qc1 = build_trotter_circuit(H, 1.0, 2, order=1)
        qc2 = build_trotter_circuit(H, 1.0, 2, order=2)

        stats1 = circuit_stats(transpile_circuit(qc1, proc_astar))
        stats2 = circuit_stats(transpile_circuit(qc2, proc_astar))

        self.assertGreater(stats2["cx"], stats1["cx"], "Order 2 must have greater CX count than Order 1.")
        self.assertGreater(stats2["depth"], stats1["depth"], "Order 2 must have greater depth than Order 1.")


if __name__ == "__main__":
    unittest.main()

