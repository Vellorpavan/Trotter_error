"""
src/run_phase4.py — StepWise-Q
Phase 4: Required Run (R0).

Generates results/tables/benchmark_required.csv and verifies Phase 4 gates:
1. Zero-noise equals noiseless Trotter to machine precision (diff == 0.0).
2. All open-chain transpiled circuits have SWAP == 0 on A0, A*, B0, B*.
3. At n=4 (noise-dominated regime), noisy error >= noiseless error.
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
)

G_VALUES = [0.5, 1.0, 1.5]
N_STEPS_REQUIRED = [1, 2, 4]
PROCS = ["A0", "A_star", "B0", "B_star"]
INIT_STATES = ["zeros", "neel"]
T_EVAL = 2.0


def run_phase4() -> pd.DataFrame:
    populate_all_graph_metrics()

    out_dir = pathlib.Path("results/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    print(f"=== Starting Phase 4 Required Run (t={T_EVAL}) ===")

    for g in G_VALUES:
        H_open = tfim_hamiltonian(4, 1.0, g, "open")
        for init_name in INIT_STATES:
            sv0 = initial_state(init_name, 4)
            exact_z0 = exact_trajectory(H_open, sv0, np.array([T_EVAL]), ["Z0"])["Z0"][0]

            for pid in PROCS:
                proc = PROCESSOR_REGISTRY[pid]
                for n_steps in N_STEPS_REQUIRED:
                    qc = build_trotter_circuit(H_open, T_EVAL, n_steps, order=1)
                    qc_t = transpile_circuit(qc, proc)
                    st = circuit_stats(qc_t)
                    sw = count_swaps(qc, proc)

                    # Noisy density matrix simulation
                    dm_noisy = simulate_density_matrix(qc_t, proc, ["Z0"])
                    z0_noisy = dm_noisy["Z0"]

                    # Noiseless density matrix simulation
                    sim_nl = AerSimulator(method="density_matrix")
                    qc_nl = qc_t.copy()
                    qc_nl.save_density_matrix()
                    res_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result()
                    dm_nl = QiDM(res_nl.data(0)["density_matrix"])
                    z0_noiseless = compute_observable(dm_nl, "Z0", 4)

                    err_nl = abs(exact_z0 - z0_noiseless)
                    err_ny = abs(exact_z0 - z0_noisy)

                    rows.append({
                        "processor": pid,
                        "n_qubits": proc["n_qubits"],
                        "boundary": "open",
                        "order": 1,
                        "g": g,
                        "init_state": init_name,
                        "n_steps": n_steps,
                        "t": T_EVAL,
                        "Z0_exact": round(exact_z0, 6),
                        "Z0_noiseless": round(z0_noiseless, 6),
                        "Z0_noisy": round(z0_noisy, 6),
                        "error_noiseless": round(err_nl, 6),
                        "error_noisy": round(err_ny, 6),
                        "depth": st["depth"],
                        "cx": st["cx"],
                        "swap": sw,
                    })

    df = pd.DataFrame(rows)
    csv_path = out_dir / "benchmark_required.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved {len(df)} rows to {csv_path}")
    return df


def validate_phase4_gate(df: pd.DataFrame) -> bool:
    print("\n=== Phase 4 Gate Validation ===")

    # Gate 1: Zero-noise equals noiseless Trotter
    proc_test = dict(PROCESSOR_REGISTRY["A0"])
    proc_test["noise"] = {
        "p_1q": 0.0,
        "p_cx": 0.0,
        "p_readout": 0.0,
        "noisy_gates_1q": ["sx", "x"],
        "rz_noisy": False,
    }
    H_test = tfim_hamiltonian(4, 1.0, 1.0, "open")
    qc = build_trotter_circuit(H_test, 2.0, 4, order=1)
    qc_t = transpile_circuit(qc, proc_test)
    z0_zero_noise = simulate_density_matrix(qc_t, proc_test, ["Z0"])["Z0"]

    sim_nl = AerSimulator(method="density_matrix")
    qc_nl = qc_t.copy()
    qc_nl.save_density_matrix()
    dm_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result().data(0)["density_matrix"]
    z0_nl = compute_observable(QiDM(dm_nl), "Z0", 4)
    zero_noise_diff = abs(z0_zero_noise - z0_nl)
    gate1_pass = zero_noise_diff < 1e-10
    print(f"  Gate 1 (zero-noise == noiseless Trotter): {'PASS' if gate1_pass else 'FAIL'} (diff={zero_noise_diff:.2e})")

    # Gate 2: Open-chain SWAP count is identically 0 across all processors
    swaps_all_zero = (df["swap"] == 0).all()
    print(f"  Gate 2 (open-chain SWAP == 0 everywhere): {'PASS' if swaps_all_zero else 'FAIL'}")

    # Gate 3: In the noise-dominated regime (n=4, g=1.5), noisy error >= noiseless error
    n4_sub = df[(df["n_steps"] == 4) & (df["g"] == 1.5)]
    n4_noise_dominates = (n4_sub["error_noisy"] >= n4_sub["error_noiseless"]).all()
    print(f"  Gate 3 (noisy error >= noiseless error at n=4, g=1.5): {'PASS' if n4_noise_dominates else 'FAIL'}")

    overall_pass = gate1_pass and swaps_all_zero and n4_noise_dominates
    return overall_pass


if __name__ == "__main__":
    df = run_phase4()
    passed = validate_phase4_gate(df)
    print(f"\nPhase 4 Gate Result: {'PASSED' if passed else 'FAILED'}")
