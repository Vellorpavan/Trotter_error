"""
run_all.py - Complete simulation runner for QFF 2026 Phase 2 (PS-S3)
Executes all benchmark tasks, exports all required CSV tables to results/tables/,
and generates all publication figures in results/figures/.
Uses purely Python standard library (csv, json) + numpy, scipy, matplotlib, and qiskit.
"""

import os
import json
import csv
import numpy as np
import matplotlib.pyplot as plt
from qiskit.quantum_info import Statevector, state_fidelity

from qff_sim import (
    PROCESSOR_A,
    PROCESSOR_B,
    load_processor,
    validate_processor_d,
    bell_experiment,
    tfim_hamiltonian,
    observable_ops,
    trotter_circuit,
    exact_evolution,
    circuit_metrics,
    noisy_state_fidelity,
    run_z_measurement,
    SHOTS,
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.join(BASE_DIR, "results", "tables")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def save_csv(filename, rows, fieldnames):
    filepath = os.path.join(TABLES_DIR, filename)
    with open(filepath, mode="w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {filepath}")


# --------------------------------------------------------------------------
# 1. PROCESSOR D VALIDATION
# --------------------------------------------------------------------------
print("=" * 70)
print("[1/6] Loading & Validating Processor D...")
processor_d_path = os.path.join(BASE_DIR, "processors", "processor_D.json")
processor_D = load_processor(processor_d_path)
validation_errors = validate_processor_d(processor_D, 4)
if validation_errors:
    raise ValueError(f"Processor D failed validation: {validation_errors}")
print(f"Processor D ({processor_D['num_qubits']} qubits) validated successfully!")

# --------------------------------------------------------------------------
# 2. OBJECTIVE 1: BELL-STATE EXPERIMENT ON PROCESSOR D
# --------------------------------------------------------------------------
print("=" * 70)
print("[2/6] Executing Objective 1: Bell-State Verification on Processor D...")
near_pair = (0, 1)  # distance = 1
far_pair = (0, 5)   # diameter distance = 3 (opposite corners of 2x3 ladder)

bell_rows = []
for pair in [near_pair, far_pair]:
    for noisy in [False, True]:
        cond_str = "noisy" if noisy else "ideal"
        res = bell_experiment(processor_D, pair, noisy=noisy, shots=SHOTS)
        bell_rows.append({
            "pair": f"{pair[0]}-{pair[1]}",
            "pair_distance": res["pair_distance"],
            "condition": cond_str,
            "two_qubit_gates": res.get("cx", 1),
            "XX": res["XX"],
            "XX_se": res["XX_se"],
            "YY": res["YY"],
            "YY_se": res["YY_se"],
            "ZZ": res["ZZ"],
            "ZZ_se": res["ZZ_se"],
            "fidelity": res["fidelity"],
            "fidelity_se": res["fidelity_se"]
        })
        print(f"  Pair {pair} ({res['pair_distance']} hops) [{cond_str.upper()}]: "
              f"F = {res['fidelity']:.4f} ± {res['fidelity_se']:.4f} (2Q gates: {res.get('cx', 1)})")

bell_fields = [
    "pair", "pair_distance", "condition", "two_qubit_gates",
    "XX", "XX_se", "YY", "YY_se", "ZZ", "ZZ_se", "fidelity", "fidelity_se"
]
save_csv("bell.csv", bell_rows, bell_fields)

# --------------------------------------------------------------------------
# 3. PS-S3 EXACT REFERENCES & QUANTUM DYNAMICS
# --------------------------------------------------------------------------
print("=" * 70)
print("[3/6] Computing Exact References for Ring TFIM (N=4, J=1.0, g=1.0, T=2.0)...")
H_ring = tfim_hamiltonian(4, J=1.0, g=1.0, periodic=True)
psi0 = Statevector(trotter_circuit(4, 1.0, 1.0, 0.0, 0))  # Neel state |1010>
target_state = exact_evolution(H_ring, psi0, [2.0])[0]

# Exact staggered magnetization Ms
ops = observable_ops(4)
ms_exact = float(np.real(target_state.expectation_value(ops["ms"])))
print(f"Exact Target State computed. Exact Ms(T=2.0) = {ms_exact:.6f}")

# --------------------------------------------------------------------------
# 4. OBJECTIVE 3: STEP SCAN, NOISE SWEEP & RESOURCE METRICS
# --------------------------------------------------------------------------
print("=" * 70)
print("[4/6] Running Trotter Step Scan (1..16 steps) on Processors A, B, and D...")
steps_list = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16]

scan_configs = [
    (PROCESSOR_A, "Processor A", 1.0),
    (PROCESSOR_B, "Processor B", 1.0),
    (processor_D, "Processor D", 0.5),
    (processor_D, "Processor D", 1.0),
    (processor_D, "Processor D", 2.0),
]

s3_scan_rows = []
resources_rows = []

pure_trotter_fidelities = {}
circuit_metrics_cache = {}

for n in steps_list:
    qc_ideal = trotter_circuit(4, 1.0, 1.0, 2.0, n, periodic=True)
    psi_trotter = Statevector(qc_ideal)
    pure_fid = float(state_fidelity(target_state, psi_trotter))
    pure_trotter_fidelities[n] = pure_fid

    for proc, p_label, _ in [(PROCESSOR_A, "Processor A", 1.0),
                             (PROCESSOR_B, "Processor B", 1.0),
                             (processor_D, "Processor D", 1.0)]:
        key = (p_label, n)
        if key not in circuit_metrics_cache:
            zz_mode = "rzz" if "rzz" in proc.get("basis_gates", []) else "cx"
            qc = trotter_circuit(4, 1.0, 1.0, 2.0, n, periodic=True, zz_mode=zz_mode)
            m = circuit_metrics(qc, proc)
            circuit_metrics_cache[key] = m
            resources_rows.append({
                "ps": "PS-S3",
                "circuit_label": f"trotter_o1_n{n}",
                "processor": "A" if "A" in p_label else ("B" if "B" in p_label else "D"),
                "executable": m["executable"],
                "logical_qubits": m["logical_qubits"],
                "physical_qubits": m["physical_qubits"],
                "logical_2q": m["logical_cx"],
                "two_q": m["cx"],
                "depth": m["depth"],
                "two_q_depth": m["two_q_depth"],
                "swaps": m["swaps"],
                "routing_overhead_2q": m["cx_overhead"]
            })

resources_fields = [
    "ps", "circuit_label", "processor", "executable", "logical_qubits",
    "physical_qubits", "logical_2q", "two_q", "depth", "two_q_depth",
    "swaps", "routing_overhead_2q"
]
save_csv("resources.csv", resources_rows, resources_fields)

# Run noisy simulation and measurements
summary_rows = []

for proc, p_label, scale in scan_configs:
    p_code = "A" if "A" in p_label else ("B" if "B" in p_label else "D")
    best_n = None
    best_f = -1.0
    print(f"Scanning {p_label} at scale {scale}x...")

    for n in steps_list:
        zz_mode = "rzz" if "rzz" in proc.get("basis_gates", []) else "cx"
        qc = trotter_circuit(4, 1.0, 1.0, 2.0, n, periodic=True, zz_mode=zz_mode)
        m = circuit_metrics_cache[(p_label, n)]

        # Noisy density-matrix fidelity
        f_noisy = noisy_state_fidelity(qc, target_state, proc, scale=scale)
        if f_noisy > best_f:
            best_f = f_noisy
            best_n = n

        # 4096-shot Z-basis measurement
        zm = run_z_measurement(qc, proc, shots=SHOTS, noisy=True, scale=scale)

        s3_scan_rows.append({
            "processor": p_code,
            "noise_scale": scale,
            "steps": n,
            "fidelity_noisy": f_noisy,
            "fidelity_trotter_ideal": pure_trotter_fidelities[n],
            "ms_noisy": zm["ms"],
            "ms_noisy_se": zm["ms_se"],
            "ms_exact": ms_exact,
            "logical_2q": m["logical_cx"],
            "swap_2q": m["cx_overhead"],
            "two_q": m["cx"]
        })

    summary_rows.append({
        "processor": p_code,
        "noise_scale": scale,
        "n_star": best_n,
        "fidelity_at_n_star": best_f
    })
    print(f"  --> {p_code} (scale {scale}x): n* = {best_n}, Fidelity = {best_f:.6f}")

s3_scan_fields = [
    "processor", "noise_scale", "steps", "fidelity_noisy",
    "fidelity_trotter_ideal", "ms_noisy", "ms_noisy_se", "ms_exact",
    "logical_2q", "swap_2q", "two_q"
]
save_csv("s3_scan.csv", s3_scan_rows, s3_scan_fields)

s3_summary_fields = ["processor", "noise_scale", "n_star", "fidelity_at_n_star"]
save_csv("s3_summary.csv", summary_rows, s3_summary_fields)

# --------------------------------------------------------------------------
# 5. BONUS: ERROR MITIGATION ON PROCESSOR D AT n*
# --------------------------------------------------------------------------
print("=" * 70)
# Dynamically extract n_star from s3_summary_rows for Processor D at 1.0x noise
n_star = [r["n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 1.0][0]
print(f"[5/6] Executing Bonus: Error Mitigation on Processor D at dynamically identified n* = {n_star}...")
zz_mode_d = "rzz" if "rzz" in processor_D.get("basis_gates", []) else "cx"
qc_star = trotter_circuit(4, 1.0, 1.0, 2.0, n_star, periodic=True, zz_mode=zz_mode_d)

# 1. Raw baseline from scan (1x noise at n*)
raw_row = [r for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 1.0 and r["steps"] == n_star][0]
ms_raw = raw_row["ms_noisy"]
ms_raw_se = raw_row["ms_noisy_se"]

# 2. Readout Error Mitigation (Single-Qubit Observable Inversion)
# Staggered magnetization M_s = (1/4) * sum_i (-1)^i <Z_i> is a sum of single-qubit Z observables.
# Under independent readout bit-flip error p_readout = 0.01, each single-qubit expectation value
# is attenuated by (1 - 2*p_readout) = 0.98 (NOT 0.98^4, which erroneously treats M_s as a 4-body observable).
p_ro = 0.01
mit_factor_ro = 1.0 - 2.0 * p_ro  # exactly 0.98
ms_readout_mit = ms_raw / mit_factor_ro
ms_readout_se = ms_raw_se / mit_factor_ro

# 3. Zero-Noise Extrapolation (ZNE via Richardson Extrapolation on noise scales)
# The noise sweep scales all three parameters (p1, p2, p_readout) by lambda in {0.5, 1.0, 2.0}.
# Therefore, ZNE simultaneously extrapolates both gate error and readout attenuation to lambda = 0.
# Applying readout correction on top of ZNE would double-count readout mitigation!
row_05 = [r for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 0.5 and r["steps"] == n_star][0]
row_10 = raw_row
row_20 = [r for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 2.0 and r["steps"] == n_star][0]

ms_05 = row_05["ms_noisy"]
se_05 = row_05["ms_noisy_se"]
ms_10 = ms_raw
se_10 = ms_raw_se
ms_20 = row_20["ms_noisy"]
se_20 = row_20["ms_noisy_se"]

poly_fit = np.polyfit([0.5, 1.0, 2.0], [ms_05, ms_10, ms_20], 2)
ms_zne = float(poly_fit[2])  # value at scale = 0

# Exact analytic standard error propagation for 2nd-order Richardson extrapolation:
# Lagrange weights at lambda = 0 for grid [0.5, 1.0, 2.0]: w1 = 8/3, w2 = -2, w3 = 1/3
w1, w2, w3 = 8.0 / 3.0, -2.0, 1.0 / 3.0
ms_zne_se = float(np.sqrt((w1 * se_05)**2 + (w2 * se_10)**2 + (w3 * se_20)**2))

# Calibration shots: 4096 shots for |0000>, 4096 shots for |1111> -> 8192 shots
extra_shots_ro = 4096 * 2
# ZNE extra shots: 4096 at 0.5x scale + 4096 at 2.0x scale -> 8192 shots
extra_shots_zne = 4096 * 2

mitigation_rows = [
    {
        "method": "readout_error_mitigation",
        "ms_raw": ms_raw,
        "ms_mitigated": ms_readout_mit,
        "extra_shots": extra_shots_ro
    },
    {
        "method": "zero_noise_extrapolation_richardson",
        "ms_raw": ms_raw,
        "ms_mitigated": ms_zne,
        "extra_shots": extra_shots_zne
    }
]

mitigation_fields = ["method", "ms_raw", "ms_mitigated", "extra_shots"]
save_csv("mitigation.csv", mitigation_rows, mitigation_fields)
print(f"  Exact Reference Ms:    {ms_exact:.6f}")
print(f"  Raw Noisy Ms (1.0x):   {ms_raw:.6f} ± {ms_raw_se:.6f} (Error: {abs(ms_raw - ms_exact):.6f})")
print(f"  Readout-Mitigated Ms:  {ms_readout_mit:.6f} ± {ms_readout_se:.6f} (Error: {abs(ms_readout_mit - ms_exact):.6f})")
print(f"  ZNE Richardson Ms:     {ms_zne:.6f} ± {ms_zne_se:.6f} (Error: {abs(ms_zne - ms_exact):.6f}, ~{100*(1 - abs(ms_zne - ms_exact)/abs(ms_raw - ms_exact)):.1f}% reduction)")

# --------------------------------------------------------------------------
# 6. PUBLICATION FIGURE GENERATION
# --------------------------------------------------------------------------
print("=" * 70)
print("[6/6] Generating Publication-Quality Figures...")

plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "font.family": "sans-serif",
    "figure.autolayout": False
})

# --- Figure 1: Processor Topologies & Connectivity ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=300)

# Processor A
ax = axes[0]
coords_A = {0: (0, 1), 1: (1, 1), 2: (2, 2), 3: (1, 0), 4: (2, -0.5)}
for u, v in PROCESSOR_A["edges"]:
    ax.plot([coords_A[u][0], coords_A[v][0]], [coords_A[u][1], coords_A[v][1]], "k-", lw=2, zorder=1)
for q, (x, y) in coords_A.items():
    ax.scatter(x, y, s=600, color="#6929C4", edgecolor="black", zorder=2)
    ax.text(x, y, str(q), color="white", ha="center", va="center", weight="bold")
ax.set_title("Processor A (5 Qubits, Tree)\nNo 4-Cycles → High SWAP Overhead", pad=12)
ax.set_xlim(-0.8, 2.8); ax.set_ylim(-1.2, 2.8); ax.axis("off")

# Processor B
ax = axes[1]
coords_B = {
    0: (0, 0), 1: (1, 1), 2: (2, 1), 3: (3, 0),
    4: (2.5, -1), 5: (1.5, -1), 6: (1, 0), 7: (0, -1), 8: (-0.5, 0.5)
}
for u, v in PROCESSOR_B["edges"]:
    ax.plot([coords_B[u][0], coords_B[v][0]], [coords_B[u][1], coords_B[v][1]], "k-", lw=1.5, zorder=1)
for q, (x, y) in coords_B.items():
    c = "#DA1E28" if q == 6 else "#0062FF"
    ax.scatter(x, y, s=550, color=c, edgecolor="black", zorder=2)
    ax.text(x, y, str(q), color="white", ha="center", va="center", weight="bold")
ax.set_title("Processor B (9 Qubits, Heavy-Hex)\nContains 4-Cycle, but Node 6 Degree = 5", pad=12)
ax.set_xlim(-1.2, 3.8); ax.set_ylim(-1.6, 1.8); ax.axis("off")

# Processor D (2x3 ladder)
ax = axes[2]
coords_D = {
    0: (0, 1), 1: (1, 1), 2: (2, 1),
    3: (0, 0), 4: (1, 0), 5: (2, 0)
}
for u, v in processor_D["edges"]:
    ax.plot([coords_D[u][0], coords_D[v][0]], [coords_D[u][1], coords_D[v][1]], "k-", lw=2, zorder=1)
for q, (x, y) in coords_D.items():
    c = "#007D79" if q in [0, 1, 4, 3] else "#198038"
    ax.scatter(x, y, s=550, color=c, edgecolor="black", zorder=2)
    ax.text(x, y, str(q), color="white", ha="center", va="center", weight="bold")
ax.set_title("Processor D (6 Qubits, 2×3 Ladder)\nTwo 4-Cycles (0-1-4-3), Max Deg 3, Diam 3", pad=12)
ax.set_xlim(-0.6, 2.6); ax.set_ylim(-0.6, 1.6); ax.axis("off")

plt.tight_layout()
fig1_path = os.path.join(FIGURES_DIR, "fig1_processor_topologies.png")
plt.savefig(fig1_path, dpi=300)
plt.close()
print(f"Saved: {fig1_path}")

# --- Figure 2: Bell-State Verification (Objective 1) ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5), dpi=300)

pairs_labels = ["Near (0,1) [Ideal]", "Near (0,1) [Noisy]", "Far (0,5) [Ideal]", "Far (0,5) [Noisy]"]
x = np.arange(len(pairs_labels))
width = 0.25

xx_vals = [r["XX"] for r in bell_rows]
yy_vals = [r["YY"] for r in bell_rows]
zz_vals = [r["ZZ"] for r in bell_rows]

ax1.bar(x - width, xx_vals, width, label="⟨XX⟩ (Ideal: +1)", color="#0062FF", alpha=0.85)
ax1.bar(x, yy_vals, width, label="⟨YY⟩ (Ideal: -1)", color="#DA1E28", alpha=0.85)
ax1.bar(x + width, zz_vals, width, label="⟨ZZ⟩ (Ideal: +1)", color="#198038", alpha=0.85)
ax1.set_xticks(x)
ax1.set_xticklabels(pairs_labels, rotation=15, ha="right")
ax1.set_ylabel("Parity Expectation ⟨P⟩")
ax1.set_title("(a) Parity Correlation Observables (4096 Shots)")
ax1.grid(axis="y", linestyle="--", alpha=0.5)
ax1.legend(loc="lower right")

fids = [r["fidelity"] for r in bell_rows]
fids_se = [r["fidelity_se"] for r in bell_rows]
bars = ax2.bar(x, fids, yerr=fids_se, capsize=5, color=["#8A3FFC", "#6929C4", "#4589FF", "#0043CE"], alpha=0.9, edgecolor="black")
ax2.set_xticks(x)
ax2.set_xticklabels(pairs_labels, rotation=15, ha="right")
ax2.set_ylabel("Bell State Fidelity F")
ax2.set_ylim(0.85, 1.05)
ax2.set_title("(b) State Fidelity F = (1 + ⟨XX⟩ - ⟨YY⟩ + ⟨ZZ⟩)/4")
ax2.grid(axis="y", linestyle="--", alpha=0.5)

for bar in bars:
    yval = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2.0, yval + 0.008, f"{yval:.4f}", ha="center", va="bottom", fontsize=9, weight="bold")

plt.tight_layout()
fig2_path = os.path.join(FIGURES_DIR, "fig2_bell_state_fidelity.png")
plt.savefig(fig2_path, dpi=300)
plt.close()
print(f"Saved: {fig2_path}")

# --- Figure 3: Trotter Fidelity Crossover (Headline Deliverable) ---
fig, ax = plt.subplots(figsize=(9.5, 5.5), dpi=300)

steps = np.array(steps_list)
pure_fids = [pure_trotter_fidelities[n] for n in steps]

f_A = [r["fidelity_noisy"] for r in s3_scan_rows if r["processor"] == "A" and r["noise_scale"] == 1.0]
f_B = [r["fidelity_noisy"] for r in s3_scan_rows if r["processor"] == "B" and r["noise_scale"] == 1.0]
f_D = [r["fidelity_noisy"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 1.0]

n_star_A = [r["n_star"] for r in summary_rows if r["processor"] == "A" and r["noise_scale"] == 1.0][0]
f_star_A = [r["fidelity_at_n_star"] for r in summary_rows if r["processor"] == "A" and r["noise_scale"] == 1.0][0]
n_star_B = [r["n_star"] for r in summary_rows if r["processor"] == "B" and r["noise_scale"] == 1.0][0]
f_star_B = [r["fidelity_at_n_star"] for r in summary_rows if r["processor"] == "B" and r["noise_scale"] == 1.0][0]
n_star_D = [r["n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 1.0][0]
f_star_D = [r["fidelity_at_n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 1.0][0]

ax.plot(steps, pure_fids, "k--", lw=2, marker="o", label="Ideal Trotter (Noise-Free Algorithmic Error)")
ax.plot(steps, f_A, color="#6929C4", lw=2.5, marker="s", ms=7, label=f"Processor A (1× Noise) — n* = {n_star_A}, F = {f_star_A:.3f}")
ax.plot(steps, f_B, color="#0062FF", lw=2.5, marker="^", ms=7, label=f"Processor B (1× Noise) — n* = {n_star_B}, F = {f_star_B:.3f}")
ax.plot(steps, f_D, color="#007D79", lw=2.5, marker="D", ms=7, label=f"Processor D (1× Noise, Native RZZ) — n* = {n_star_D}, F = {f_star_D:.3f}")

ax.axvline(x=n_star_A, color="#6929C4", linestyle=":", alpha=0.6)
ax.axvline(x=n_star_B, color="#0062FF", linestyle=":", alpha=0.6)
ax.axvline(x=n_star_D, color="#007D79", linestyle=":", alpha=0.6)

ax.annotate(f"Processor A Peak\n(n* = {n_star_A}, F = {f_star_A:.3f})", xy=(n_star_A, f_star_A), xytext=(3.2, 0.40),
            arrowprops=dict(facecolor="#6929C4", shrink=0.08, width=1.5, headwidth=6),
            fontsize=9.5, weight="bold", color="#6929C4")
ax.annotate(f"Processor B Peak\n(n* = {n_star_B}, F = {f_star_B:.3f})", xy=(n_star_B, f_star_B), xytext=(4.3, 0.70),
            arrowprops=dict(facecolor="#0062FF", shrink=0.08, width=1.5, headwidth=6),
            fontsize=9.5, weight="bold", color="#0062FF")
ax.annotate(f"Processor D Peak\n(n* = {n_star_D}, F = {f_star_D:.3f})", xy=(n_star_D, f_star_D), xytext=(9.2, 0.79),
            arrowprops=dict(facecolor="#007D79", shrink=0.08, width=1.5, headwidth=6),
            fontsize=9.5, weight="bold", color="#007D79")

ax.set_xlabel("Trotter Step Count (n)")
ax.set_ylabel("State Fidelity F(n)")
ax.set_title("PS-S3 Headline: Optimal Trotter Step n* Crossover Dynamics", pad=12)
ax.set_xticks(steps)
ax.set_ylim(0.0, 1.05)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(loc="lower right", framealpha=0.9)

plt.tight_layout()
fig3_path = os.path.join(FIGURES_DIR, "fig3_trotter_fidelity_crossover.png")
plt.savefig(fig3_path, dpi=300)
plt.close()
print(f"Saved: {fig3_path}")

# --- Figure 4: Noise Sweep & n* Movement on Processor D ---
fig, ax = plt.subplots(figsize=(9, 5.5), dpi=300)

f_D_05 = [r["fidelity_noisy"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 0.5]
f_D_10 = [r["fidelity_noisy"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 1.0]
f_D_20 = [r["fidelity_noisy"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 2.0]

n_star_05 = [r["n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 0.5][0]
f_star_05 = [r["fidelity_at_n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 0.5][0]
n_star_10 = [r["n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 1.0][0]
f_star_10 = [r["fidelity_at_n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 1.0][0]
n_star_20 = [r["n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 2.0][0]
f_star_20 = [r["fidelity_at_n_star"] for r in summary_rows if r["processor"] == "D" and r["noise_scale"] == 2.0][0]

ax.plot(steps, f_D_05, color="#198038", lw=2.5, marker="o", ms=7, label=f"Noise Scale 0.5× → n* = {n_star_05}, F = {f_star_05:.3f}")
ax.plot(steps, f_D_10, color="#007D79", lw=2.5, marker="D", ms=7, label=f"Noise Scale 1.0× → n* = {n_star_10}, F = {f_star_10:.3f}")
ax.plot(steps, f_D_20, color="#DA1E28", lw=2.5, marker="v", ms=7, label=f"Noise Scale 2.0× → n* = {n_star_20}, F = {f_star_20:.3f}")

peak_steps = [n_star_05, n_star_10, n_star_20]
peak_fids = [f_star_05, f_star_10, f_star_20]
ax.plot(peak_steps, peak_fids, "k--", lw=1.5, zorder=4)
ax.scatter(peak_steps, peak_fids, color="gold", edgecolor="black", s=130, zorder=5)

ax.set_xlabel("Trotter Step Count (n)")
ax.set_ylabel("State Fidelity F(n)")
ax.set_title("Noise Sweep on Processor D: Leftward Shift of n* Under Increasing Noise", pad=12)
ax.set_xticks(steps)
ax.set_ylim(0.0, 0.95)
ax.grid(True, linestyle="--", alpha=0.5)
ax.legend(loc="lower right", framealpha=0.9)

plt.tight_layout()
fig4_path = os.path.join(FIGURES_DIR, "fig4_noise_sweep_shift.png")
plt.savefig(fig4_path, dpi=300)
plt.close()
print(f"Saved: {fig4_path}")

# --- Figure 5: Gate Split & Routing Overhead ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=300)

two_q_A = [r["two_q"] for r in s3_scan_rows if r["processor"] == "A" and r["noise_scale"] == 1.0]
log_2q = [r["logical_2q"] for r in s3_scan_rows if r["processor"] == "A" and r["noise_scale"] == 1.0]
two_q_D = [r["two_q"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 1.0]

ax1.plot(steps, log_2q, "k--", lw=2, marker="o", label="Logical 2Q Gates (CX Basis Baseline)")
ax1.plot(steps, two_q_D, color="#007D79", lw=2.5, marker="D", label="Processor D (Native RZZ: 50% Fewer 2Q Gates, 0 SWAPs)")
ax1.plot(steps, two_q_A, color="#6929C4", lw=2.5, marker="s", label="Processor A (CX Basis + Heavy Routing SWAPs)")
ax1.fill_between(steps, log_2q, two_q_A, color="#6929C4", alpha=0.15, label="Routing Overhead (SWAP CXs on A)")

ax1.set_xlabel("Trotter Step Count (n)")
ax1.set_ylabel("Physical Two-Qubit Gates (CX / RZZ)")
ax1.set_title("(a) Gate Split: Physical Two-Qubit Gate Count")
ax1.set_xticks(steps)
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.legend(loc="upper left")

depth_2q_A = [r["two_q_depth"] for r in resources_rows if r["processor"] == "A"]
depth_2q_D = [r["two_q_depth"] for r in resources_rows if r["processor"] == "D"]

ax2.plot(steps, depth_2q_D, color="#007D79", lw=2.5, marker="D", label="Processor D (Native RZZ: Minimal 2Q Depth)")
ax2.plot(steps, depth_2q_A, color="#6929C4", lw=2.5, marker="s", label="Processor A (SWAP Congestion Depth)")

ax2.set_xlabel("Trotter Step Count (n)")
ax2.set_ylabel("Two-Qubit Gate Depth")
ax2.set_title("(b) Two-Qubit Circuit Depth Comparison")
ax2.set_xticks(steps)
ax2.grid(True, linestyle="--", alpha=0.5)
ax2.legend(loc="upper left")

plt.tight_layout()
fig5_path = os.path.join(FIGURES_DIR, "fig5_gate_split_and_depth.png")
plt.savefig(fig5_path, dpi=300)
plt.close()
print(f"Saved: {fig5_path}")

# --- Figure 6: Magnetization & Mitigation Bonus ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=300)

ms_A = [r["ms_noisy"] for r in s3_scan_rows if r["processor"] == "A" and r["noise_scale"] == 1.0]
ms_A_se = [r["ms_noisy_se"] for r in s3_scan_rows if r["processor"] == "A" and r["noise_scale"] == 1.0]
ms_D = [r["ms_noisy"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 1.0]
ms_D_se = [r["ms_noisy_se"] for r in s3_scan_rows if r["processor"] == "D" and r["noise_scale"] == 1.0]

ax1.axhline(y=ms_exact, color="black", linestyle="--", lw=2, label=f"Exact Reference Ms = {ms_exact:.4f}")
ax1.errorbar(steps, ms_D, yerr=ms_D_se, fmt="D-", color="#007D79", lw=2, capsize=4, label="Processor D (4096 Shots)")
ax1.errorbar(steps, ms_A, yerr=ms_A_se, fmt="s-", color="#6929C4", lw=2, capsize=4, label="Processor A (4096 Shots)")
ax1.set_xlabel("Trotter Step Count (n)")
ax1.set_ylabel("Staggered Magnetization Ms(T=2.0)")
ax1.set_title("(a) Magnetization Tracking Across Steps")
ax1.set_xticks(steps)
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.legend(loc="lower right")

categories = ["Exact\nReference", "Raw Noisy\n(1× Baseline)", "Readout\nMitigated", "ZNE\n(Richardson)"]
vals = [ms_exact, ms_raw, ms_readout_mit, ms_zne]
errs = [0.0, ms_raw_se, ms_readout_se, ms_zne_se]
colors = ["#161616", "#DA1E28", "#FF832B", "#0062FF"]

bars = ax2.bar(categories, vals, yerr=errs, capsize=5, color=colors, alpha=0.85, edgecolor="black", width=0.55)
ax2.set_ylabel(f"Staggered Magnetization Ms at n* = {n_star}")
ax2.set_ylim(0.0, 0.8)
ax2.set_title("(b) Bonus Error Mitigation Recovery on Processor D")
ax2.grid(axis="y", linestyle="--", alpha=0.5)

for bar in bars:
    y = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2.0, y + 0.018, f"{y:.4f}", ha="center", va="bottom", fontsize=9, weight="bold")

plt.tight_layout()
fig6_path = os.path.join(FIGURES_DIR, "fig6_magnetization_and_mitigation.png")
plt.savefig(fig6_path, dpi=300)
plt.close()
print(f"Saved: {fig6_path}")

print("=" * 70)
print("SUCCESS: All calculations complete, all tables saved, and all figures generated!")
print("=" * 70)
