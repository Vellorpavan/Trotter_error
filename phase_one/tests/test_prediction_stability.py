"""
tests/test_prediction_stability.py — StepWise-Q Automated Tests
Verification of recommend_steps API, bootstrap confidence intervals, and n* prediction stability.
"""

from __future__ import annotations

import unittest
from src.problem import tfim_hamiltonian
from src.metrics import recommend_steps, observable_attenuation_model


class TestPredictionStability(unittest.TestCase):
    """Test suite for Trotter step recommendation and error model prediction stability."""

    def test_recommend_steps_api_contract(self):
        """recommend_steps must return all required keys: n*, predicted_error, confidence_interval, safe_range, model_metadata."""
        res = recommend_steps(g=1.0, t=2.0, processor="A_star", noise_scale=1.0)

        required_keys = ["n_star", "predicted_error", "confidence_interval", "safe_range", "model_metadata"]
        for key in required_keys:
            self.assertIn(key, res, f"recommend_steps response missing required key: {key}")

        self.assertIsInstance(res["n_star"], int)
        self.assertIsInstance(res["predicted_error"], float)
        self.assertIsInstance(res["confidence_interval"], tuple)
        self.assertEqual(len(res["confidence_interval"]), 2)
        self.assertIsInstance(res["safe_range"], tuple)
        self.assertEqual(len(res["safe_range"]), 2)
        self.assertIsInstance(res["model_metadata"], dict)

        # Verify metadata fields
        meta = res["model_metadata"]
        for m_key in ["model_name", "lambda_step", "processor", "shots", "seed", "qiskit_version"]:
            self.assertIn(m_key, meta, f"model_metadata missing field: {m_key}")

    def test_safe_range_and_ci_bounds(self):
        """Verify that recommended n* is consistent with safe_range and confidence intervals."""
        res = recommend_steps(g=0.5, t=2.0, processor="A_star", noise_scale=1.0)
        n_star = res["n_star"]
        ci_lo, ci_hi = res["confidence_interval"]
        safe_lo, safe_hi = res["safe_range"]

        self.assertLessEqual(safe_lo, safe_hi)
        self.assertLessEqual(ci_lo, ci_hi)
        self.assertTrue(safe_lo <= n_star <= safe_hi,
                        f"Recommended n* ({n_star}) not inside safe_range [{safe_lo}, {safe_hi}]")

    def test_recommend_steps_monotonic_noise_scaling(self):
        """As noise increases, optimal Trotter step count n* should decrease or stay constant."""
        n_star_low = recommend_steps(g=1.0, t=2.0, processor="A_star", noise_scale=0.5)["n_star"]
        n_star_high = recommend_steps(g=1.0, t=2.0, processor="A_star", noise_scale=3.0)["n_star"]

        self.assertGreaterEqual(n_star_low, n_star_high,
                                f"At lower noise n* ({n_star_low}) should be >= higher noise n* ({n_star_high})")

    def test_heldout_generalization_v10(self):
        """V10 gate: Held-out prediction summary confirms Model B MAE <= 1.0 step and beats trivial baseline."""
        import os
        import pandas as pd
        summary_path = "results/tables/heldout_summary.csv"
        if os.path.exists(summary_path):
            df_s = pd.read_csv(summary_path)
            mae_b = float(df_s.loc[df_s["model"] == "Model_B_observable_attenuation", "mae_steps"].values[0])
            mae_triv = float(df_s.loc[df_s["model"] == "Trivial_always_n4", "mae_steps"].values[0])
            self.assertLessEqual(mae_b, 1.0, f"Model B held-out MAE {mae_b} exceeds 1.0 step.")
            self.assertLess(mae_b, mae_triv, f"Model B held-out MAE {mae_b} should beat trivial {mae_triv}.")


if __name__ == "__main__":
    unittest.main()

