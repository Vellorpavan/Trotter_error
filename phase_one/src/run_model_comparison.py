"""
src/run_model_comparison.py — StepWise-Q
Phase 8 Extension: Model Comparison Study (Additive vs. Observable Attenuation).

Compares:
  Model A: Heuristic Additive Model:  epsilon(n) = c/n^k + (1 - Lambda^n)
  Model B: Physical Observable Model: epsilon(n) = | Lambda^n * Z_nl(n) - Z_exact |

Output:
  - results/tables/model_comparison.csv
"""

from __future__ import annotations

import pathlib
import numpy as np
import pandas as pd
import scipy.optimize

from src.metrics import trotter_error_model_fit, observable_attenuation_model


def run_model_comparison() -> pd.DataFrame:
    ext_path = pathlib.Path("results/tables/benchmark_extended.csv")
    if not ext_path.exists():
        raise FileNotFoundError("benchmark_extended.csv not found")

    df_ext = pd.read_csv(ext_path)
    out_dir = pathlib.Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    print("=== Comparing Error Models across g in {0.5, 1.0, 1.5} ===")

    for g in [0.5, 1.0, 1.5]:
        sub = df_ext[
            (df_ext["processor"] == "A_star") &
            (df_ext["boundary"] == "open") &
            (df_ext["order"] == 1) &
            (df_ext["g"] == g)
        ].sort_values("n_steps")

        n_arr = sub["n_steps"].values.astype(float)
        z_ex = sub["Z0_exact"].values[0]
        z_nl = sub["Z0_noiseless"].values
        z_ny = sub["Z0_noisy"].values
        err_meas = sub["error_noisy"].values
        lambda_phys = 0.985

        # 1. Model A: Additive
        fit_a = trotter_error_model_fit(n_arr, sub["error_noiseless"].values, err_meas, k_fixed=1)
        c_a = fit_a["c"]
        pred_err_a = c_a / n_arr + (1.0 - lambda_phys ** n_arr)
        mae_a = float(np.mean(np.abs(pred_err_a - err_meas)))

        # 2. Model B: Observable Attenuation
        obs_b = observable_attenuation_model(n_arr, z_ex, z_nl, z_ny, lambda_step=lambda_phys)
        pred_err_b = obs_b["err_pred"]
        mae_b = float(np.mean(np.abs(pred_err_b - err_meas)))

        for i, n in enumerate(n_arr):
            rows.append({
                "processor": "A_star",
                "g": g,
                "n_steps": int(n),
                "error_measured": round(float(err_meas[i]), 6),
                "model_A_additive_error": round(float(pred_err_a[i]), 6),
                "model_B_attenuation_error": round(float(pred_err_b[i]), 6),
                "discrepancy_A": round(float(abs(pred_err_a[i] - err_meas[i])), 6),
                "discrepancy_B": round(float(abs(pred_err_b[i] - err_meas[i])), 6),
            })

        print(f"  g={g:3.1f}: Measured n*={n_arr[np.argmin(err_meas)]} | "
              f"Model A n*={fit_a['n_star_pred']} (MAE={mae_a:.4f}) | "
              f"Model B n*={obs_b['n_star_obs_pred']} (MAE={mae_b:.4f})")

    df_res = pd.DataFrame(rows)
    csv_path = out_dir / "model_comparison.csv"
    df_res.to_csv(csv_path, index=False)
    print(f"Saved {len(df_res)} rows to {csv_path}")
    return df_res


if __name__ == "__main__":
    run_model_comparison()

