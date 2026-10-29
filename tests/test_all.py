"""
test_all.py - Complete Automated Verification Suite for QFF 2026 Phase 2 (PS-S3)
Can be executed with: python -m unittest discover tests or python tests/test_all.py
"""

import os
import sys
import json
import csv
import unittest
import numpy as np

# Ensure project root is in path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from qff_sim import (
    PROCESSOR_A,
    PROCESSOR_B,
    load_processor,
    validate_processor_d,
    tfim_hamiltonian,
    observable_ops,
    trotter_circuit,
    exact_evolution,
    circuit_metrics,
    noisy_state_fidelity,
    run_z_measurement,
    SEED,
    SHOTS,
)
from qiskit.quantum_info import Statevector, state_fidelity


class TestProcessorD(unittest.TestCase):
    """Audits Processor D configuration against official design rules."""

    def setUp(self):
        self.path = os.path.join(PROJECT_ROOT, "processors", "processor_D.json")
        self.processor = load_processor(self.path)

    def test_file_exists(self):
        self.assertTrue(os.path.exists(self.path), "processors/processor_D.json must exist")

    def test_official_validation(self):
        errors = validate_processor_d(self.processor, 4)
        self.assertEqual(errors, [], f"Processor D failed official validation: {errors}")

    def test_qubit_count(self):
        n = self.processor["num_qubits"]
        self.assertTrue(4 <= n <= 12, f"num_qubits must be between 4 and 12, got {n}")
        self.assertEqual(n, 9, "Processor D should have 9 qubits (3x3 grid)")

    def test_max_degree(self):
        deg = {i: 0 for i in range(self.processor["num_qubits"])}
        for u, v in self.processor["edges"]:
            deg[u] += 1
            deg[v] += 1
        max_deg = max(deg.values())
        self.assertLessEqual(max_deg, 4, f"Max node degree {max_deg} exceeds limit of 4")

    def test_basis_gates(self):
        allowed = [sorted(["rz", "sx", "x", "cx"]), sorted(["rz", "sx", "x", "cz"]), sorted(["rz", "sx", "x", "rzz"])]
        actual = sorted(self.processor["basis_gates"])
        self.assertIn(actual, allowed, "basis_gates must be an official allowed gate set")

    def test_error_rates(self):
        rates = self.processor["error_rates"]
        self.assertGreaterEqual(rates.get("p1", 0), 0.001)
        self.assertGreaterEqual(rates.get("p2", 0), 0.005)
        self.assertGreaterEqual(rates.get("p_readout", 0), 0.01)

    def test_justification_length(self):
        words = self.processor["justification"].split()
        self.assertGreaterEqual(len(words), 40, f"Justification must be >= 40 words, got {len(words)}")


class TestCSVOutputs(unittest.TestCase):
    """Audits all generated CSV tables against templates/results_columns.md schemas."""

    def setUp(self):
        self.tables_dir = os.path.join(PROJECT_ROOT, "results", "tables")
        self.expected_schemas = {
            "bell.csv": [
                "pair", "pair_distance", "condition", "two_qubit_gates",
                "XX", "XX_se", "YY", "YY_se", "ZZ", "ZZ_se", "fidelity", "fidelity_se"
            ],
            "resources.csv": [
                "ps", "circuit_label", "processor", "executable", "logical_qubits",
                "physical_qubits", "logical_2q", "two_q", "depth", "two_q_depth",
                "swaps", "routing_overhead_2q"
            ],
            "s3_scan.csv": [
                "processor", "noise_scale", "steps", "fidelity_noisy",
                "fidelity_trotter_ideal", "ms_noisy", "ms_noisy_se", "ms_exact",
                "logical_2q", "swap_2q", "two_q"
            ],
            "s3_summary.csv": [
                "processor", "noise_scale", "n_star", "fidelity_at_n_star"
            ],
            "mitigation.csv": [
                "method", "ms_raw", "ms_mitigated", "extra_shots"
            ]
        }

    def test_all_tables_exist_and_match_schema(self):
        for filename, expected_cols in self.expected_schemas.items():
            path = os.path.join(self.tables_dir, filename)
            self.assertTrue(os.path.exists(path), f"Required table {filename} does not exist")

            with open(path, mode="r", newline="") as f:
                reader = csv.DictReader(f)
                actual_cols = reader.fieldnames
                self.assertEqual(actual_cols, expected_cols,
                                 f"Column mismatch in {filename}.\nExpected: {expected_cols}\nActual: {actual_cols}")
                rows = list(reader)
                self.assertGreater(len(rows), 0, f"{filename} is empty")

                # Verify no NaN or empty values in any row
                for r_idx, row in enumerate(rows):
                    for col in expected_cols:
                        val = row[col]
                        self.assertIsNotNone(val, f"Null value in {filename} row {r_idx} col {col}")
                        self.assertNotEqual(val, "", f"Empty value in {filename} row {r_idx} col {col}")
                        self.assertNotEqual(val.lower(), "nan", f"NaN in {filename} row {r_idx} col {col}")

    def test_row_counts(self):
        # bell.csv: 2 pairs x 2 conditions = 4 rows
        with open(os.path.join(self.tables_dir, "bell.csv")) as f:
            self.assertEqual(len(list(csv.reader(f))) - 1, 4)

        # s3_summary.csv: A@1.0, B@1.0, D@0.5, D@1.0, D@2.0 = 5 rows
        with open(os.path.join(self.tables_dir, "s3_summary.csv")) as f:
            self.assertEqual(len(list(csv.reader(f))) - 1, 5)

        # s3_scan.csv: 5 configs x 12 steps = 60 rows
        with open(os.path.join(self.tables_dir, "s3_scan.csv")) as f:
            self.assertEqual(len(list(csv.reader(f))) - 1, 60)


class TestPhysicsAndSimulation(unittest.TestCase):
    """Validates mathematical correctness of quantum evolution, Trotterization, and n* identification."""

    def test_neel_state_preparation(self):
        qc = trotter_circuit(4, 1.0, 1.0, 0.0, 0)
        sv = Statevector(qc)
        # Qiskit '1010': qubit 1 and 3 in |1> -> integer index 10 in computational basis
        probs = sv.probabilities()
        self.assertAlmostEqual(probs[10], 1.0, places=6, msg="Initial state must be exact Neel 1010")

    def test_exact_staggered_magnetization(self):
        H_ring = tfim_hamiltonian(4, J=1.0, g=1.0, periodic=True)
        psi0 = Statevector(trotter_circuit(4, 1.0, 1.0, 0.0, 0))
        target_state = exact_evolution(H_ring, psi0, [2.0])[0]
        ops = observable_ops(4)
        ms_exact = float(np.real(target_state.expectation_value(ops["ms"])))
        self.assertAlmostEqual(ms_exact, 0.652444, places=4, msg="Exact Ms(T=2.0) must match official value ~0.652444")

    def test_optimal_step_counts(self):
        summary_path = os.path.join(PROJECT_ROOT, "results", "tables", "s3_summary.csv")
        with open(summary_path) as f:
            rows = list(csv.DictReader(f))

        results = {(r["processor"], float(r["noise_scale"])): int(r["n_star"]) for r in rows}
        self.assertEqual(results[("A", 1.0)], 6, "Processor A optimal step must be 6")
        self.assertEqual(results[("B", 1.0)], 7, "Processor B optimal step must be 7")
        self.assertEqual(results[("D", 1.0)], 7, "Processor D optimal step must be 7")
        self.assertEqual(results[("D", 0.5)], 8, "Processor D at 0.5x noise optimal step must be 8")
        self.assertEqual(results[("D", 2.0)], 6, "Processor D at 2.0x noise optimal step must be 6")

    def test_processor_d_zero_swap_overhead(self):
        res_path = os.path.join(PROJECT_ROOT, "results", "tables", "resources.csv")
        with open(res_path) as f:
            rows = list(csv.DictReader(f))

        d_rows = [r for r in rows if r["processor"] == "D"]
        for r in d_rows:
            self.assertEqual(int(r["swaps"]), 0, f"Processor D must have 0 SWAPs, got {r['swaps']} at {r['circuit_label']}")
            self.assertEqual(int(r["routing_overhead_2q"]), 0, "Processor D must have 0 routing overhead")


class TestBonusMitigation(unittest.TestCase):
    """Validates the Error Mitigation procedure and evidence."""

    def test_mitigation_recovery(self):
        mit_path = os.path.join(PROJECT_ROOT, "results", "tables", "mitigation.csv")
        with open(mit_path) as f:
            rows = list(csv.DictReader(f))

        target = 0.652444
        for r in rows:
            raw = float(r["ms_raw"])
            mit = float(r["ms_mitigated"])
            # Mitigated value should be closer to target than raw
            self.assertLess(abs(mit - target), abs(raw - target),
                            f"Method {r['method']} failed to improve accuracy: raw error {abs(raw - target)}, mit error {abs(mit - target)}")


if __name__ == "__main__":
    unittest.main()

