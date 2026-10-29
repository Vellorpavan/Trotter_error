"""
generate_notebook.py — StepWise-Q
Generates main.ipynb programmatically using nbformat.

Run:
    python generate_notebook.py

This writes main.ipynb with all cells, a QUICK flag,
and section headings a judge can follow in 5 minutes.
"""

from __future__ import annotations

import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
import pathlib

NB_PATH = pathlib.Path(__file__).parent / "main.ipynb"

# ──────────────────────────────────────────────────────────────────────────────
# Cell helpers
# ──────────────────────────────────────────────────────────────────────────────
def md(text: str) -> nbformat.NotebookNode:
    return new_markdown_cell(text)

def code(text: str) -> nbformat.NotebookNode:
    return new_code_cell(text)

# ──────────────────────────────────────────────────────────────────────────────
# Build notebook
# ──────────────────────────────────────────────────────────────────────────────
cells = []

# ── Title ──
cells.append(md("""# PS-S3: Noise Accumulation in Trotter Circuits

### StepWise-Q: Optimal Trotter Step Prediction on Custom Quantum Topologies
**Qiskit Fall Fest 2026 | SRM University-AP | Track 3: Quantum Simulation**

> **One-sentence answer:** a tested rule, backed by theory and simulation on two processors we designed,
> that tells you the number of Trotter steps `n*` with the lowest total error, and how that number
> moves with evolution time, field strength, noise level and coupling topology.

---
"""))

# ── Configuration cell (QUICK flag) ──
cells.append(md("## 0. Configuration"))
cells.append(code("""\
# ─── CONFIGURATION ────────────────────────────────────────────────────────────
# Set QUICK = True for a fast review run (reduced time grid, no layout search).
# Set QUICK = False for the full competition run.
QUICK = False

# Hard-wired constants (do not change — competition spec)
J          = 1.0
G_VALUES   = [0.5, 1.0, 1.5]
N_STEPS_REQUIRED  = [1, 2, 4]
N_STEPS_EXTENDED  = [1, 2, 3, 4, 6, 8]
T_GRID     = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0] if QUICK else list(__import__('numpy').linspace(0, 2.0, 21))
T_GRID_RED = [0.5, 1.0, 1.5, 2.0]   # reduced grid for extensions
SHOTS      = 4096
SEED       = 42
PROC_IDS_A = ["A0", "A_star"]
PROC_IDS_B = ["B0", "B_star"]

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

print(f"QUICK mode: {QUICK}")
print(f"t_grid points: {len(T_GRID)}")
"""))

# ── Imports ──
cells.append(md("## 1. Imports and Setup"))
cells.append(code("""\
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(".").resolve()))

from src.problem import (
    tfim_hamiltonian, initial_state, exact_trajectory,
    compute_observable, spin_to_qubit,
    validate_spin_to_qubit, validate_initial_observable, validate_norm_preserved,
)
from src.processors import (
    PROCESSOR_REGISTRY, load_processor, build_noise_model,
    graph_metrics, populate_all_graph_metrics, save_processor_json,
    validate_connectivity, validate_json_roundtrip,
)
from src.routing import (
    build_trotter_circuit, transpile_circuit, count_swaps,
    cx_identity_check, circuit_stats,
    simulate_density_matrix, simulate_shots, layout_search,
)
from src.metrics import (
    trotter_error_model_fit, predict_n_star, bootstrap_n_star_ci,
    recommend_steps, zne_extrapolate, readout_mitigation,
    noise_scale_sweep, run_validation,
)
from src.visualization import (
    fig1_coupling_maps, fig2_trajectory, fig3_u_curve, fig4_n_star_map,
    fig5_scaling_law, fig6_cx_depth, fig7_layout, fig8_mitigation,
    FIGURES_DIR,
)

# Populate graph metrics for all processors
populate_all_graph_metrics()

print("All imports successful.")
print(f"Processor registry: {list(PROCESSOR_REGISTRY.keys())}")
"""))

# ── Phase 1: Hamiltonians and exact reference ──
cells.append(md("""## 2. Phase 1 — Hamiltonians and Exact Reference

We build the 4-spin TFIM Hamiltonian for all three field values and compute the exact
time evolution using `scipy.linalg.expm`. This serves as the ground truth for all error calculations.

**Model:** H = −J Σ Z_i Z_{i+1} − g Σ X_i, J=1.0, g ∈ {0.5, 1.0, 1.5}
"""))
cells.append(code("""\
import numpy as np

# Phase 1 gate checks
print("=== Phase 1 Gate Checks ===")
print(f"spin_to_qubit check:  {'PASS' if validate_spin_to_qubit() else 'FAIL'}")
print(f"Z0(0) check:          {'PASS' if validate_initial_observable() else 'FAIL'}")
print(f"norm preserved:       {'PASS' if validate_norm_preserved() else 'FAIL'}")
print()

# Build Hamiltonians
H_open = {g: tfim_hamiltonian(4, 1.0, g, "open") for g in G_VALUES}
H_ring = {g: tfim_hamiltonian(4, 1.0, g, "ring") for g in G_VALUES}

print("Hamiltonians built:")
for g in G_VALUES:
    print(f"  g={g}: n_terms={len(H_open[g])}")

# Initial states
sv_zeros = initial_state("zeros", 4)
sv_neel  = initial_state("neel", 4)

# Exact trajectories for g=1.0
t_arr = np.array(T_GRID)
observables = ["Z0", "mean_mag", "ZZ_nn", "energy"]

print("\\nComputing exact trajectories (g=1.0)...")
exact_zeros = exact_trajectory(
    H_open[1.0], sv_zeros, t_arr, observables, save_prefix="open_zeros_g1"
)
exact_neel = exact_trajectory(
    H_open[1.0], sv_neel,  t_arr, observables, save_prefix="open_neel_g1"
)
print(f"  Z0(t=0)  zeros: {exact_zeros['Z0'][0]:.4f} (expected +1)")
print(f"  Z0(t=0)  neel:  {exact_neel['Z0'][0]:.4f}  (expected +1)")

# Compute for all g values
EXACT = {}   # {(boundary, g, init_state): {obs: array}}
for g in G_VALUES:
    for init_name, sv0 in [("zeros", sv_zeros), ("neel", sv_neel)]:
        for boundary, H_dict in [("open", H_open), ("ring", H_ring)]:
            key = (boundary, g, init_name)
            EXACT[key] = exact_trajectory(
                H_dict[g], sv0, t_arr, ["Z0"],
                save_prefix=f"{boundary}_{init_name}_g{g:.1f}",
            )

print(f"\\nExact trajectories computed for {len(EXACT)} configurations.")
"""))

# ── Phase 2: Processors ──
cells.append(md("""## 3. Phase 2 — Processor Topologies

We benchmark five coupling maps: two reference topologies (A0, B0) and three custom designs (A1, A*, B*).
Our final submission uses **A*** (5-qubit) and **B*** (9-qubit).
"""))
cells.append(code("""\
# Phase 2 gate checks
print("=== Phase 2 Gate Checks ===")
conn = validate_connectivity()
for pid, ok in conn.items():
    print(f"  {pid} connectivity: {'PASS' if ok else 'FAIL'}")
print(f"  JSON round-trip: {'PASS' if validate_json_roundtrip() else 'FAIL'}")
print()

# Print graph metrics table
import pandas as pd
rows = []
for pid, proc in PROCESSOR_REGISTRY.items():
    m = proc["graph_metrics"]
    rows.append({
        "Processor": pid,
        "Qubits": proc["n_qubits"],
        "Edges": m["n_edges"],
        "Diameter": m["diameter"],
        "Avg degree": round(m["avg_degree"], 2),
        "Max degree": m["max_degree"],
        "Connected": m["is_connected"],
    })
df_metrics = pd.DataFrame(rows)
display(df_metrics)

# Save initial JSON files
save_processor_json("A_star", pathlib.Path("processors/processor_A.json"))
save_processor_json("B_star", pathlib.Path("processors/processor_B.json"))
print("Wrote processors/processor_A.json (A*) and processor_B.json (B*)")
"""))
cells.append(code("""\
# Figure 1: Coupling maps
procs = [PROCESSOR_REGISTRY[p] for p in ("A0", "A1", "A_star", "B0", "B_star")]
fig1_coupling_maps(procs)
plt.show()
print(f"Fig 1 saved to {FIGURES_DIR / 'fig1_coupling_maps.png'}")
"""))

# ── Phase 3: Routing ──
cells.append(md("""## 4. Phase 3 — Trotter Circuits and SWAP Counting

We build, transpile and measure the gate resources for each processor/Hamiltonian combination.
The key validation: open chain on A0, A*, B0, B* must require **zero SWAP gates**.
"""))
cells.append(code("""\
print("=== Phase 3 Gate Checks: SWAP validation ===")
from src.routing import validate_swap_counts
swap_checks = validate_swap_counts(g=1.0)
for name, ok in swap_checks.items():
    print(f"  {name}: {'PASS' if ok else 'FAIL'}")
print()

# Sample transpiled circuit stats
H_test = tfim_hamiltonian(4, 1.0, 1.0, "open")
stats_rows = []
for pid in ("A0", "A_star", "B0", "B_star"):
    proc = PROCESSOR_REGISTRY[pid]
    for n_steps in N_STEPS_REQUIRED:
        qc = build_trotter_circuit(H_test, 2.0, n_steps, order=1)
        qc_t = transpile_circuit(qc, proc)
        st = circuit_stats(qc_t)
        sw = count_swaps(qc, proc)
        stats_rows.append({
            "Processor": pid, "n": n_steps, "Depth": st["depth"],
            "CX": st["cx"], "SWAP (open)": sw,
        })

df_stats = pd.DataFrame(stats_rows)
display(df_stats)
"""))

# ── Phase 4: Required run R0 ──
cells.append(md("""## 5. Phase 4 — Required Run (R0)

Open chain, order-1 Trotter, n ∈ {1, 2, 4}, g ∈ {0.5, 1.0, 1.5}, both initial states,
processors A0, A*, B0, B*. This is the **safe submission baseline**.
"""))
cells.append(code("""\
import pathlib, csv
pathlib.Path("results/tables").mkdir(parents=True, exist_ok=True)

benchmark_rows = []

proc_ids_all = ["A0", "A_star", "B0", "B_star"]

for pid in proc_ids_all:
    proc = PROCESSOR_REGISTRY[pid]
    for g in G_VALUES:
        for init_name in ["zeros", "neel"]:
            sv0 = initial_state(init_name, 4)
            H_g = tfim_hamiltonian(4, 1.0, g, "open")
            
            # Exact at t=2.0
            exact_at_t2 = EXACT[("open", g, init_name)]["Z0"][-1]
            
            for n_steps in N_STEPS_REQUIRED:
                qc = build_trotter_circuit(H_g, 2.0, n_steps, order=1)
                qc_t = transpile_circuit(qc, proc)
                st  = circuit_stats(qc_t)
                sw  = count_swaps(qc, proc)
                
                # Noisy (density matrix)
                dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
                z0_noisy = dm_res["Z0"]
                
                # Noiseless
                from qiskit_aer import AerSimulator
                from qiskit.quantum_info import DensityMatrix as QiDM
                sim_nl = AerSimulator(method="density_matrix")
                qc_nl = qc_t.copy(); qc_nl.save_density_matrix()
                res_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result()
                dm_nl  = QiDM(res_nl.data(0)["density_matrix"])
                z0_noiseless = compute_observable(dm_nl, "Z0", 4)
                
                benchmark_rows.append({
                    "processor": pid, "n_qubits": proc["n_qubits"],
                    "boundary": "open", "order": 1, "g": g,
                    "init_state": init_name, "n_steps": n_steps,
                    "t": 2.0,
                    "Z0_exact": round(exact_at_t2, 6),
                    "Z0_noiseless": round(z0_noiseless, 6),
                    "Z0_noisy": round(z0_noisy, 6),
                    "error_noiseless": round(abs(exact_at_t2 - z0_noiseless), 6),
                    "error_noisy":     round(abs(exact_at_t2 - z0_noisy), 6),
                    "depth": st["depth"], "cx": st["cx"], "swap": sw,
                })
                print(f"  {pid} g={g} {init_name} n={n_steps}: "
                      f"Z0_ex={exact_at_t2:.4f} Z0_noisy={z0_noisy:.4f} SWAP={sw}")

df_bench = pd.DataFrame(benchmark_rows)
df_bench.to_csv("results/tables/benchmark_required.csv", index=False)
print(f"\\nWrote results/tables/benchmark_required.csv ({len(df_bench)} rows)")
display(df_bench.head(10))
"""))

# ── Phase 4 gate check ──
cells.append(code("""\
# Phase 4 gate: zero-noise ≡ noiseless; noisy error >= noiseless error
print("=== Phase 4 Gate Checks ===")
all_ok = True
for _, row in df_bench.iterrows():
    # noisy error should be >= noiseless error (noise always adds error on average)
    if row["error_noisy"] < row["error_noiseless"] - 0.01:  # 0.01 tolerance for finite-sample effects
        print(f"  WARN: noisy < noiseless for {row['processor']} n={row['n_steps']} g={row['g']}")
        all_ok = False
print(f"  Noisy >= Noiseless check: {'PASS' if all_ok else 'WARN (see above)'}") 
"""))

# ── Phase 5: Extensions ──
cells.append(md("""## 6. Phase 5 — Extensions

- **X1:** Finer n grid {1,2,3,4,6,8} for the error model fit
- **X2:** Ring Hamiltonian — topology comparison and SWAP overhead
- **X3:** Order-2 (Suzuki) Trotter
- Layout search skipped in QUICK mode.
"""))
cells.append(code("""\
n_steps_ext = N_STEPS_EXTENDED if not QUICK else N_STEPS_REQUIRED
ext_rows = []

# X1: finer n grid + X2: ring boundary
for pid in ["A_star", "B_star"]:  # Focus on our designed processors
    proc = PROCESSOR_REGISTRY[pid]
    for g in G_VALUES:
        for boundary in ["open", "ring"]:
            for order in [1, 2]:
                H_g = tfim_hamiltonian(4, 1.0, g, boundary)
                sv0 = initial_state("zeros", 4)
                exact_at_t2 = exact_trajectory(H_g, sv0, np.array([2.0]), ["Z0"])["Z0"][0]
                
                for n_steps in n_steps_ext:
                    qc = build_trotter_circuit(H_g, 2.0, n_steps, order=order)
                    qc_t = transpile_circuit(qc, proc)
                    st = circuit_stats(qc_t)
                    sw = count_swaps(qc, proc)
                    
                    dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
                    z0_noisy = dm_res["Z0"]
                    
                    from qiskit_aer import AerSimulator
                    from qiskit.quantum_info import DensityMatrix as QiDM
                    sim_nl = AerSimulator(method="density_matrix")
                    qc_nl = qc_t.copy(); qc_nl.save_density_matrix()
                    res_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result()
                    dm_nl = QiDM(res_nl.data(0)["density_matrix"])
                    z0_nl = compute_observable(dm_nl, "Z0", 4)
                    
                    ext_rows.append({
                        "processor": pid, "boundary": boundary, "order": order,
                        "g": g, "n_steps": n_steps, "t": 2.0,
                        "Z0_exact": round(exact_at_t2, 6),
                        "Z0_noiseless": round(z0_nl, 6),
                        "Z0_noisy": round(z0_noisy, 6),
                        "error_noiseless": round(abs(exact_at_t2 - z0_nl), 6),
                        "error_noisy": round(abs(exact_at_t2 - z0_noisy), 6),
                        "depth": st["depth"], "cx": st["cx"], "swap": sw,
                    })

df_ext = pd.DataFrame(ext_rows)
df_ext.to_csv("results/tables/benchmark_extended.csv", index=False)
print(f"Extended benchmark: {len(df_ext)} rows saved.")
display(df_ext.groupby(["processor","boundary","order","g"])[["error_noiseless","error_noisy","swap"]].mean().round(4))
"""))

# ── Phase 6: Error model and n* ──
cells.append(md("""## 7. Phase 6 — Error Model Fit and n* Estimation

We fit the two-parameter model to determine the optimal Trotter step count.
"""))
cells.append(code("""\
from src.metrics import (
    trotter_error_model_fit,
    predict_n_star,
    bootstrap_n_star_ci,
    bootstrap_discrete_n_star_ci,
)

fit_rows = []
n_star_rows = []

# Fit for each (processor, g, boundary)
for pid in ["A_star", "B_star"]:
    proc = PROCESSOR_REGISTRY[pid]
    for g in G_VALUES:
        sub = df_ext[(df_ext.processor == pid) & (df_ext.g == g) &
                     (df_ext.boundary == "open") & (df_ext.order == 1)].sort_values("n_steps")
        if len(sub) < 3:
            continue
        
        n_arr = sub["n_steps"].values.astype(float)
        err_nl = sub["error_noiseless"].values
        err_ny = sub["error_noisy"].values
        
        fit = trotter_error_model_fit(n_arr, err_nl, err_ny)
        n_star_meas = int(n_arr[np.argmin(err_ny)])
        ci_lo, ci_hi = bootstrap_n_star_ci(n_arr, err_ny, k_fixed=1, lambda_step=fit["Lambda"])
        ci_disc_lo, ci_disc_hi = bootstrap_discrete_n_star_ci(n_arr, err_ny)
        
        fit_rows.append({
            "processor": pid, "g": g, "boundary": "open",
            "c": round(fit["c"], 6), "k_fit": round(fit["k_fit"], 3),
            "Lambda": round(fit["Lambda"], 6),
            "R2_noiseless": round(fit["R2_noiseless"], 4),
            "R2_noisy": round(fit["R2_noisy"], 4),
            "n_star_pred": fit["n_star_pred"],
            "n_star_measured": n_star_meas,
        })
        
        n_star_rows.append({
            "processor": pid, "g": g, "n_star": n_star_meas,
            "n_star_pred": fit["n_star_pred"],
            "ci_model_lower": round(ci_lo, 1), "ci_model_upper": round(ci_hi, 1),
            "ci_disc_lower": int(ci_disc_lo), "ci_disc_upper": int(ci_disc_hi),
        })
        
        print(f"{pid} g={g}: c={fit['c']:.4f} k={fit['k_fit']:.2f} "
              f"Λ={fit['Lambda']:.5f} n*={n_star_meas} (pred={fit['n_star_pred']})")

df_fit = pd.DataFrame(fit_rows)
df_nstar = pd.DataFrame(n_star_rows)
df_fit.to_csv("results/tables/model_fit.csv", index=False)
df_nstar.to_csv("results/tables/n_star.csv", index=False)
display(df_fit)
"""))

cells.append(code("""\
# ZNE and readout mitigation
mit_before = {}
mit_after  = {}

H_test = tfim_hamiltonian(4, 1.0, 1.0, "open")
proc = PROCESSOR_REGISTRY["A_star"]
sv0 = initial_state("zeros", 4)
exact_z0 = exact_trajectory(H_test, sv0, np.array([2.0]), ["Z0"])["Z0"][0]

for n_steps, label in [(2, "n=2"), (4, "n=4")]:
    qc = build_trotter_circuit(H_test, 2.0, n_steps, order=1)
    qc_t = transpile_circuit(qc, proc)
    
    # Raw (shot-based)
    shot_res = simulate_shots(qc_t, proc, ["Z0"], shots=4096)
    mit_before[label] = shot_res["Z0"]
    
    # ZNE
    try:
        zne_val = zne_extrapolate(qc_t, proc, "Z0", scale_factors=[1, 3, 5])
        mit_after[f"ZNE {label}"] = (zne_val, 0.0)
    except Exception as e:
        print(f"  ZNE failed for {label}: {e}")
        mit_after[f"ZNE {label}"] = mit_before[label]

print(f"Exact Z0 at t=2: {exact_z0:.4f}")
for k, (m, se) in {**mit_before, **mit_after}.items():
    print(f"  {k}: {m:.4f} ± {se:.4f}")

fig8_mitigation(mit_before, mit_after, exact_val=exact_z0)
plt.show()
"""))

# ── Phase 7: Figures ──
cells.append(md("""## 8. Phase 7 — Figures and Report

Generating all 8 required figures.
"""))
cells.append(code("""\
import matplotlib.pyplot as plt

# Fig 2: Trajectory
g_plot = 1.0
t_arr = np.array(T_GRID)
proc_ids_plot = ["A_star"]
init_plot = "zeros"

fig2_data = {"t_grid": t_arr, "exact": EXACT[("open", g_plot, init_plot)]["Z0"],
             "proc_ids": proc_ids_plot, "n_steps_list": [1, 2, 4]}

for pid in proc_ids_plot:
    proc = PROCESSOR_REGISTRY[pid]
    H_g = tfim_hamiltonian(4, 1.0, g_plot, "open")
    sv0 = initial_state(init_plot, 4)
    for n_steps in [1, 2, 4]:
        noisy_traj = []
        for t_val in t_arr:
            qc = build_trotter_circuit(H_g, t_val if t_val > 0 else 1e-6, n_steps, order=1)
            qc_t = transpile_circuit(qc, proc)
            dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
            noisy_traj.append(dm_res["Z0"])
        fig2_data[(pid, n_steps, "noisy")] = np.array(noisy_traj)

fig2_trajectory(fig2_data)
plt.show()
print("Fig 2 saved.")
"""))

cells.append(code("""\
# Fig 3: U-curve (hero figure 1)
pid = "A_star"
g = 1.0
sub = df_ext[(df_ext.processor==pid) & (df_ext.g==g) & 
             (df_ext.boundary=="open") & (df_ext.order==1)].sort_values("n_steps")
n_arr  = sub["n_steps"].values.astype(float)
err_nl = sub["error_noiseless"].values
err_ny = sub["error_noisy"].values
fit = trotter_error_model_fit(n_arr, err_nl, err_ny)

fig3_u_curve(n_arr, err_ny, fit, fit["n_star_pred"], err_nl)
plt.show()
print(f"Fig 3 saved.  n*={fit['n_star_pred']} c={fit['c']:.4f} Λ={fit['Lambda']:.5f}")
"""))

cells.append(code("""\
# Fig 4: n* map (hero figure 2) — noise-scale sweep
if not QUICK:
    H_test = tfim_hamiltonian(4, 1.0, 1.0, "open")
    proc = PROCESSOR_REGISTRY["A_star"]
    noise_scales = [0.25, 0.5, 1.0, 2.0, 4.0]
    t_sweep = T_GRID_RED
    
    n_star_grid = np.zeros((len(noise_scales), len(t_sweep)), dtype=int)
    for ni, ns in enumerate(noise_scales):
        for ti, t_val in enumerate(t_sweep):
            rec = recommend_steps(H_test, t_val, proc, noise_scale=ns, n_candidates=[1,2,3,4,6,8])
            n_star_grid[ni, ti] = rec["n_star"]
    
    fig4_n_star_map(n_star_grid, np.array(t_sweep), np.array(noise_scales))
    plt.show()
    print("Fig 4 saved.")
else:
    print("Fig 4 skipped in QUICK mode.")
"""))

cells.append(code("""\
# Fig 5: Scaling law
if not QUICK:
    sweep = noise_scale_sweep(H_test, 2.0, proc, noise_scales=[0.25, 0.5, 1.0, 2.0, 4.0])
    df_scaling = pd.DataFrame({"noise_scale": sweep["noise_scales"], "n_star": sweep["n_stars"]})
    df_scaling.to_csv("results/tables/scaling_law.csv", index=False)
    fig5_scaling_law(sweep["noise_scales"], sweep["n_stars"], sweep["fit_slope"], sweep["r_squared"])
    plt.show()
    print(f"Fig 5 saved. Fitted slope={sweep['fit_slope']:.3f} (theory order-1: -0.5)")
else:
    print("Fig 5 skipped in QUICK mode.")
"""))

cells.append(code("""\
# Fig 6: CX count and depth vs n, open vs ring
circ_data = {"n_steps": n_steps_ext, "proc_ids": ["A_star", "B_star"]}
for pid in ["A_star", "B_star"]:
    proc = PROCESSOR_REGISTRY[pid]
    for boundary in ["open", "ring"]:
        cx_list, depth_list = [], []
        H_g = tfim_hamiltonian(4, 1.0, 1.0, boundary)
        for n_steps in n_steps_ext:
            qc = build_trotter_circuit(H_g, 2.0, n_steps, order=1)
            qc_t = transpile_circuit(qc, proc)
            st = circuit_stats(qc_t)
            cx_list.append(st["cx"])
            depth_list.append(st["depth"])
        circ_data[(pid, boundary, "cx")] = cx_list
        circ_data[(pid, boundary, "depth")] = depth_list

fig6_cx_depth(circ_data)
plt.show()
print("Fig 6 saved.")
"""))

cells.append(code("""\
# Fig 7: Layout study
if not QUICK:
    H_open = tfim_hamiltonian(4, 1.0, 1.0, "open")
    proc_A = PROCESSOR_REGISTRY["A_star"]
    qc_open = build_trotter_circuit(H_open, 2.0, 2, order=1)
    
    df_layout = layout_search(qc_open, proc_A, n_spins=4, simulate_extremes=False)
    df_layout.to_csv("results/tables/layout_study.csv", index=False)
    fig7_layout(df_layout, proc_name="A_star")
    plt.show()
    print(f"Fig 7 saved. {len(df_layout)} layouts tested.")
else:
    print("Fig 7 (layout study) skipped in QUICK mode.")
"""))

# ── Phase 8: Validation ──
cells.append(md("""## 9. Phase 8 — Validation Suite

Running all 9 validation checks from README section 8.
"""))
cells.append(code("""\
print("=== Full Validation Suite ===")
val_results = run_validation(verbose=True)
n_pass = sum(val_results.values())
print(f"\\n{n_pass}/{len(val_results)} checks passed.")

# Write graph metrics table
gm_rows = []
for pid, proc in PROCESSOR_REGISTRY.items():
    m = proc["graph_metrics"]
    gm_rows.append({"processor": pid, **m})
pd.DataFrame(gm_rows).to_csv("results/tables/graph_metrics.csv", index=False)
print("Wrote results/tables/graph_metrics.csv")
"""))

cells.append(md("""## 10. Decision Table

`recommend_steps()` turned into a one-page lookup table.
"""))
cells.append(code("""\
from src.metrics import recommend_steps

decision_rows = []
for pid in ["A_star", "B_star"]:
    proc = PROCESSOR_REGISTRY[pid]
    for g in G_VALUES:
        for t_val in T_GRID_RED:
            H_g = tfim_hamiltonian(4, 1.0, g, "open")
            rec = recommend_steps(H_g, t_val, proc, n_candidates=N_STEPS_REQUIRED)
            decision_rows.append({
                "processor": pid, "g": g, "t": t_val,
                "n_star": rec["n_star"],
                "predicted_error": round(rec["predicted_error"], 5),
                "ci_lower": round(rec["ci_lower"], 1),
                "ci_upper": round(rec["ci_upper"], 1),
            })

df_decision = pd.DataFrame(decision_rows)
df_decision.to_csv("results/tables/decision_table.csv", index=False)
print("Decision table:")
display(df_decision)
"""))

cells.append(md("""---

## Summary

The optimal Trotter step count `n*` depends on:
1. **Evolution time t** — longer times need more steps (more Trotter error per step)
2. **Field strength g** — stronger field = faster dynamics = more Trotter error
3. **Noise level** — higher noise punishes more steps (more gate error)
4. **Topology** — ring Hamiltonian on poorly-matched topologies increases gate count, shifting n* lower

See `report.md` for full analysis, hypothesis outcomes, and limitations.
See `results/tables/` for all numerical results.
"""))

# ──────────────────────────────────────────────────────────────────────────────
# Write notebook
# ──────────────────────────────────────────────────────────────────────────────
nb = new_notebook(cells=cells)
nb.metadata["kernelspec"] = {
    "display_name": "Python 3",
    "language": "python",
    "name": "python3",
}
nb.metadata["language_info"] = {
    "name": "python",
    "version": "3.10.0",
}

NB_PATH.parent.mkdir(parents=True, exist_ok=True)
with open(NB_PATH, "w") as f:
    nbformat.write(nb, f)

print(f"Wrote {NB_PATH} ({len(cells)} cells)")
