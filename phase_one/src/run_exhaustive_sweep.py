"""
src/run_exhaustive_sweep.py — StepWise-Q
Phase 8 Extension: Exhaustive Integer Step Sweep (n = 1 to 10).

This script provides full reproducibility for the fine-grained step sweep
verifying that n* = 4 (g=0.5), n* = 6 (g=1.0), and n* = 4 (g=1.5) are the
true global discrete minima, and not artifacts of skipping odd steps in {1, 2, 3, 4, 6, 8}.

Output:
  - results/tables/exhaustive_sweep_n1_to_10.csv
"""

from __future__ import annotations

import pathlib
import numpy as np
import pandas as pd
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix as QiDM

from src.problem import tfim_hamiltonian, initial_state, exact_trajectory, compute_observable
from src.processors import PROCESSOR_REGISTRY
from src.routing import build_trotter_circuit, transpile_circuit, simulate_density_matrix


def run_exhaustive_sweep(
    processor_id: str = "A_star",
    g_values: list[float] | None = None,
    t_eval: float = 2.0,
    n_steps_list: list[int] | None = None,
) -> pd.DataFrame:
    """Run density matrix simulation across fine-grained integer steps n in [1..10]."""
    if g_values is None:
        g_values = [0.5, 1.0, 1.5]
    if n_steps_list is None:
        n_steps_list = list(range(1, 11))

    proc = PROCESSOR_REGISTRY[processor_id]
    out_dir = pathlib.Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    print(f"=== Running Exhaustive Sweep (n = {n_steps_list[0]} to {n_steps_list[-1]}) on {processor_id} ===")

    for g in g_values:
        H = tfim_hamiltonian(4, 1.0, g, "open")
        exact_z0 = exact_trajectory(H, initial_state("zeros"), np.array([t_eval]), ["Z0"])["Z0"][0]

        for n in n_steps_list:
            qc = build_trotter_circuit(H, t_eval, n, order=1)
            qc_t = transpile_circuit(qc, proc)

            # 1. Noisy density matrix
            dm_noisy = simulate_density_matrix(qc_t, proc, ["Z0"])
            z0_noisy = dm_noisy["Z0"]

            # 2. Noiseless density matrix
            sim_nl = AerSimulator(method="density_matrix")
            qc_nl = qc_t.copy()
            qc_nl.save_density_matrix()
            res_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result()
            dm_nl = QiDM(res_nl.data(0)["density_matrix"])
            z0_nl = compute_observable(dm_nl, "Z0", 4)

            err_nl = abs(exact_z0 - z0_nl)
            err_ny = abs(exact_z0 - z0_noisy)

            rows.append({
                "processor": processor_id,
                "g": g,
                "n_steps": n,
                "t": t_eval,
                "Z0_exact": round(exact_z0, 6),
                "Z0_noiseless": round(z0_nl, 6),
                "Z0_noisy": round(z0_noisy, 6),
                "error_noiseless": round(err_nl, 6),
                "error_noisy": round(err_ny, 6),
            })

    df = pd.DataFrame(rows)
    csv_path = out_dir / "exhaustive_sweep_n1_to_10.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved {len(df)} rows to {csv_path}")

    # Summary of discrete minima:
    print("\n=== Verified Global Discrete Minima on n in [1..10] ===")
    for g in g_values:
        sub = df[df["g"] == g]
        min_idx = sub["error_noisy"].idxmin()
        best = sub.loc[min_idx]
        print(f"  g={g:3.1f}: Best n* = {int(best['n_steps'])}, error_noisy = {best['error_noisy']:.6f}")

    return df


if __name__ == "__main__":
    run_exhaustive_sweep()

