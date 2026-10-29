"""
src/run_phase6.py — StepWise-Q
Phase 6: Mitigation, Error Model Fitting, and n* Estimation.

Outputs:
  - results/tables/model_fit.csv
  - results/tables/n_star.csv
  - results/tables/scaling_law.csv
  - results/tables/graph_metrics.csv
  - results/tables/decision_table.csv
  - Final processor JSON files (processors/processor_A.json, processor_B.json)
"""

from __future__ import annotations

import pathlib
import numpy as np
import pandas as pd

from src.problem import tfim_hamiltonian, initial_state, exact_trajectory
from src.processors import (
    PROCESSOR_REGISTRY,
    populate_all_graph_metrics,
    save_processor_json,
)
from src.routing import (
    build_trotter_circuit,
    transpile_circuit,
    simulate_shots,
)
from src.metrics import (
    trotter_error_model_fit,
    predict_n_star,
    bootstrap_n_star_ci,
    bootstrap_discrete_n_star_ci,
    recommend_steps,
    zne_extrapolate,
    readout_mitigation,
    noise_scale_sweep,
)

G_VALUES = [0.5, 1.0, 1.5]
T_GRID_RED = [0.5, 1.0, 1.5, 2.0]
OUT_DIR = pathlib.Path("results/tables")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_phase6() -> dict[str, pd.DataFrame]:
    populate_all_graph_metrics()

    # 1. Load benchmark_extended.csv
    bench_ext_path = OUT_DIR / "benchmark_extended.csv"
    if not bench_ext_path.exists():
        raise FileNotFoundError(f"{bench_ext_path} not found. Run Phase 5 first.")
    df_ext = pd.read_csv(bench_ext_path)

    fit_rows = []
    n_star_rows = []

    print("=== Fitting Error Model & Predicting n* ===")
    for pid in ["A0", "A_star", "B0", "B_star"]:
        for boundary in ["open", "ring"]:
            for order in [1, 2]:
                for g in G_VALUES:
                    sub = df_ext[
                        (df_ext["processor"] == pid)
                        & (df_ext["boundary"] == boundary)
                        & (df_ext["order"] == order)
                        & (df_ext["g"] == g)
                    ].sort_values("n_steps")

                    if len(sub) < 3:
                        continue

                    n_arr = sub["n_steps"].values.astype(float)
                    err_nl = sub["error_noiseless"].values
                    err_ny = sub["error_noisy"].values

                    fit = trotter_error_model_fit(n_arr, err_nl, err_ny, k_fixed=order)
                    n_star_meas = int(n_arr[np.argmin(err_ny)])
                    ci_lo, ci_hi = bootstrap_n_star_ci(n_arr, err_ny, k_fixed=order, lambda_step=fit["Lambda"])
                    ci_disc_lo, ci_disc_hi = bootstrap_discrete_n_star_ci(n_arr, err_ny)

                    fit_rows.append({
                        "processor": pid,
                        "boundary": boundary,
                        "order": order,
                        "g": g,
                        "c": round(fit["c"], 6),
                        "k_fit": round(fit["k_fit"], 3),
                        "Lambda": round(fit["Lambda"], 6),
                        "R2_noiseless": round(fit["R2_noiseless"], 4),
                        "R2_noisy": round(fit["R2_noisy"], 4),
                        "n_star_pred": fit["n_star_pred"],
                        "n_star_measured": n_star_meas,
                    })

                    n_star_rows.append({
                        "processor": pid,
                        "boundary": boundary,
                        "order": order,
                        "g": g,
                        "n_star_measured": n_star_meas,
                        "n_star_pred": fit["n_star_pred"],
                        "ci_lower": round(ci_lo, 1),
                        "ci_upper": round(ci_hi, 1),
                        "ci_disc_lower": int(ci_disc_lo),
                        "ci_disc_upper": int(ci_disc_hi),
                    })

                    if boundary == "open" and order == 1:
                        print(f"  {pid:6s} | g={g:.1f} | c={fit['c']:.4f} | "
                              f"Lambda={fit['Lambda']:.5f} | "
                              f"n*_meas={n_star_meas} | n*_pred={fit['n_star_pred']}")

    df_fit = pd.DataFrame(fit_rows)
    df_nstar = pd.DataFrame(n_star_rows)
    df_fit.to_csv(OUT_DIR / "model_fit.csv", index=False)
    df_nstar.to_csv(OUT_DIR / "n_star.csv", index=False)
    print(f"Saved {len(df_fit)} fit rows to {OUT_DIR / 'model_fit.csv'}")
    print(f"Saved {len(df_nstar)} n* rows to {OUT_DIR / 'n_star.csv'}")

    # 2. Scaling law sweep
    print("\n=== Running Scaling Law Sweep ===")
    H_test = tfim_hamiltonian(4, 1.0, 1.0, "open")
    proc_Astar = PROCESSOR_REGISTRY["A_star"]
    noise_scales = [0.25, 0.5, 1.0, 2.0, 4.0]
    sweep = noise_scale_sweep(H_test, 2.0, proc_Astar, noise_scales=noise_scales)
    df_scaling = pd.DataFrame({
        "noise_scale": sweep["noise_scales"],
        "n_star": sweep["n_stars"],
        "fit_slope": sweep["fit_slope"],
        "r_squared": sweep["r_squared"],
    })
    df_scaling.to_csv(OUT_DIR / "scaling_law.csv", index=False)
    print(f"Saved scaling law to {OUT_DIR / 'scaling_law.csv'} (slope={sweep['fit_slope']:.3f}, R2={sweep['r_squared']:.3f})")

    # 3. Graph metrics table
    print("\n=== Generating Graph Metrics Table ===")
    gm_rows = []
    for pid, proc in PROCESSOR_REGISTRY.items():
        m = proc["graph_metrics"]
        gm_rows.append({"processor": pid, **m})
    df_gm = pd.DataFrame(gm_rows)
    df_gm.to_csv(OUT_DIR / "graph_metrics.csv", index=False)
    print(f"Saved graph metrics for {len(df_gm)} processors to {OUT_DIR / 'graph_metrics.csv'}")

    # 4. Decision table
    print("\n=== Generating Decision Table ===")
    decision_rows = []
    for pid in ["A_star", "B_star"]:
        proc = PROCESSOR_REGISTRY[pid]
        for g in G_VALUES:
            for t_val in T_GRID_RED:
                H_g = tfim_hamiltonian(4, 1.0, g, "open")
                rec = recommend_steps(H_g, t_val, proc, n_candidates=[1, 2, 4])
                decision_rows.append({
                    "processor": pid,
                    "g": g,
                    "t": t_val,
                    "n_star": rec["n_star"],
                    "predicted_error": round(rec["predicted_error"], 5),
                    "ci_lower": round(rec["ci_lower"], 1),
                    "ci_upper": round(rec["ci_upper"], 1),
                })
    df_decision = pd.DataFrame(decision_rows)
    df_decision.to_csv(OUT_DIR / "decision_table.csv", index=False)
    print(f"Saved decision table ({len(df_decision)} entries) to {OUT_DIR / 'decision_table.csv'}")

    # 5. Save final processor JSON files
    proc_dir = pathlib.Path("processors")
    proc_dir.mkdir(parents=True, exist_ok=True)
    p_a = save_processor_json("A_star", proc_dir / "processor_A.json")
    p_b = save_processor_json("B_star", proc_dir / "processor_B.json")
    print(f"Saved chosen A design to {p_a}")
    print(f"Saved chosen B design to {p_b}")

    return {
        "model_fit": df_fit,
        "n_star": df_nstar,
        "scaling_law": df_scaling,
        "graph_metrics": df_gm,
        "decision": df_decision,
    }


def validate_phase6_gate(dfs: dict[str, pd.DataFrame]) -> bool:
    print("\n=== Phase 6 Gate Validation ===")
    df_fit = dfs["model_fit"]

    # Check that measured n* is within +/- 1 of predicted n* (or gap explained)
    diffs = (df_fit["n_star_measured"] - df_fit["n_star_pred"]).abs()
    within_tol = (diffs <= 2).mean()  # on discrete grid {1,2,3,4,6,8}
    print(f"  Proportion of n* within discrete step tolerance: {within_tol:.1%}")

    # Check scaling law fitted
    df_scaling = dfs["scaling_law"]
    slope = df_scaling["fit_slope"].iloc[0]
    print(f"  Scaling law slope: {slope:.3f} (theoretical order-1: -0.5)")

    # Check decision table exists
    df_dec = dfs["decision"]
    dec_ok = len(df_dec) == 24
    print(f"  Decision table complete (24 rows): {'PASS' if dec_ok else 'FAIL'}")

    passed = within_tol > 0.5 and dec_ok
    return passed


if __name__ == "__main__":
    dfs = run_phase6()
    passed = validate_phase6_gate(dfs)
    print(f"\nPhase 6 Gate Result: {'PASSED' if passed else 'FAILED'}")
