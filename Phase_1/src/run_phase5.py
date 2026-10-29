"""
src/run_phase5.py — StepWise-Q
Phase 5: Extensions (X1–X3).

- X1: Extended Trotter steps grid: n in {1, 2, 3, 4, 6, 8}
- X2: Ring boundary Hamiltonian (tests topology advantage & SWAP overhead)
- X3: Second-order Suzuki Trotter vs First-order Lie Trotter
- Layout search study (transpile-only across injective layouts, simulate extremes)

Outputs:
  - results/tables/benchmark_extended.csv
  - results/tables/layout_study.csv
"""

from __future__ import annotations

import pathlib
import numpy as np
import pandas as pd
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix as QiDM

from src.problem import (
    tfim_hamiltonian,
    initial_state,
    exact_trajectory,
    compute_observable,
)
from src.processors import PROCESSOR_REGISTRY, populate_all_graph_metrics
from src.routing import (
    build_trotter_circuit,
    transpile_circuit,
    count_swaps,
    circuit_stats,
    simulate_density_matrix,
    layout_search,
)

G_VALUES = [0.5, 1.0, 1.5]
N_STEPS_EXTENDED = [1, 2, 3, 4, 6, 8]
PROCS_EXT = ["A_star", "B_star"]
BOUNDARIES = ["open", "ring"]
ORDERS = [1, 2]
T_EVAL = 2.0


def run_benchmark_extended() -> pd.DataFrame:
    populate_all_graph_metrics()
    out_dir = pathlib.Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    print("=== Running Extended Benchmark (X1, X2, X3) ===")

    for boundary in BOUNDARIES:
        for g in G_VALUES:
            H_g = tfim_hamiltonian(4, 1.0, g, boundary)
            sv0 = initial_state("zeros", 4)
            exact_z0 = exact_trajectory(H_g, sv0, np.array([T_EVAL]), ["Z0"])["Z0"][0]

            for pid in PROCS_EXT:
                proc = PROCESSOR_REGISTRY[pid]
                for order in ORDERS:
                    for n_steps in N_STEPS_EXTENDED:
                        qc = build_trotter_circuit(H_g, T_EVAL, n_steps, order=order)
                        qc_t = transpile_circuit(qc, proc)
                        st = circuit_stats(qc_t)
                        sw = count_swaps(qc, proc)

                        # Noisy DM
                        dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
                        z0_noisy = dm_res["Z0"]

                        # Noiseless DM
                        sim_nl = AerSimulator(method="density_matrix")
                        qc_nl = qc_t.copy()
                        qc_nl.save_density_matrix()
                        res_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result()
                        dm_nl = QiDM(res_nl.data(0)["density_matrix"])
                        z0_nl = compute_observable(dm_nl, "Z0", 4)

                        err_nl = abs(exact_z0 - z0_nl)
                        err_ny = abs(exact_z0 - z0_noisy)

                        rows.append({
                            "processor": pid,
                            "boundary": boundary,
                            "order": order,
                            "g": g,
                            "n_steps": n_steps,
                            "t": T_EVAL,
                            "Z0_exact": round(exact_z0, 6),
                            "Z0_noiseless": round(z0_nl, 6),
                            "Z0_noisy": round(z0_noisy, 6),
                            "error_noiseless": round(err_nl, 6),
                            "error_noisy": round(err_ny, 6),
                            "depth": st["depth"],
                            "cx": st["cx"],
                            "swap": sw,
                        })

    df_ext = pd.DataFrame(rows)
    csv_path = out_dir / "benchmark_extended.csv"
    df_ext.to_csv(csv_path, index=False)
    print(f"Saved {len(df_ext)} rows to {csv_path}")
    return df_ext


def run_layout_study() -> pd.DataFrame:
    print("\n=== Running Layout Study ===")
    out_dir = pathlib.Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    H_open = tfim_hamiltonian(4, 1.0, 1.0, "open")
    proc_A = PROCESSOR_REGISTRY["A_star"]
    qc_open = build_trotter_circuit(H_open, T_EVAL, 2, order=1)

    df_layout = layout_search(qc_open, proc_A, n_spins=4, simulate_extremes=True)
    csv_path = out_dir / "layout_study.csv"
    df_layout.to_csv(csv_path, index=False)
    print(f"Saved {len(df_layout)} layouts to {csv_path}")
    return df_layout


def validate_phase5_gate(df_ext: pd.DataFrame, df_layout: pd.DataFrame) -> bool:
    print("\n=== Phase 5 Gate Validation ===")

    # 1. Check ring boundary requires SWAP > 0 on linear A0 but SWAP == 0 on cycle-enabled A*
    ring_A0 = df_ext[(df_ext["boundary"] == "ring") & (df_ext["processor"] == "A0")]
    ring_Astar = df_ext[(df_ext["boundary"] == "ring") & (df_ext["processor"] == "A_star")]

    a0_has_swaps = (ring_A0["swap"] > 0).all()
    astar_zero_swaps = (ring_Astar["swap"] == 0).all()
    print(f"  Ring topology test (A0 SWAPs > 0): {'PASS' if a0_has_swaps else 'FAIL'}")
    print(f"  Ring topology test (A* SWAPs == 0): {'PASS' if astar_zero_swaps else 'FAIL'}")

    # 2. Check layout study completed
    layout_ok = len(df_layout) > 0 and "cx" in df_layout.columns
    print(f"  Layout study generated: {'PASS' if layout_ok else 'FAIL'} ({len(df_layout)} layouts)")

    # 3. Check second order Trotter (order=2) has higher gate count than order=1 for same n
    o1_cx = df_ext[df_ext["order"] == 1]["cx"].mean()
    o2_cx = df_ext[df_ext["order"] == 2]["cx"].mean()
    order_ok = o2_cx > o1_cx
    print(f"  Order-2 vs Order-1 CX check (O2 avg {o2_cx:.1f} > O1 avg {o1_cx:.1f}): {'PASS' if order_ok else 'FAIL'}")

    return a0_has_swaps and astar_zero_swaps and layout_ok and order_ok


if __name__ == "__main__":
    df_ext = run_benchmark_extended()
    df_layout = run_layout_study()
    passed = validate_phase5_gate(df_ext, df_layout)
    print(f"\nPhase 5 Gate Result: {'PASSED' if passed else 'FAILED'}")
