# QFF 2026 Phase 2 — Quantum Simulation Track
## PS-S3: Noise Accumulation in Trotter Circuits
### StepWise-Q: Optimal Trotter Step Prediction & Hardware-Aware Co-Design

<p align="center">
  <img src="https://img.shields.io/badge/Competition-Qiskit_Fall_Fest_2026-6929C4?style=for-the-badge&logo=qiskit&logoColor=white" alt="QFF 2026" />
  <img src="https://img.shields.io/badge/Track-3:_Quantum_Simulation-0062FF?style=for-the-badge" alt="Track 3" />
  <img src="https://img.shields.io/badge/Problem-PS--S3:_Trotter_Noise_Accumulation-198038?style=for-the-badge" alt="PS-S3" />
  <img src="https://img.shields.io/badge/Status-100%25_Verified-success?style=for-the-badge" alt="Verified" />
</p>

---

## 👥 Team & Submission Information
- **Team Name:** StepWise-Q
- **Track:** Track 3 — Quantum Simulation
- **Problem Statement:** PS-S3 — Noise Accumulation in Trotter Circuits: Finding the Best Step Count
- **Institution:** SRM University-AP
- **Official Submission Repository:** Self-contained and runnable within `QFF2026_Simulation_Challenge_Kit/`

---

## 🧭 Executive Summary & Key Results

In Trotterized quantum simulation of the 4-spin Ring Transverse-Field Ising Model ($N=4, J=1.0, g=1.0, T=2.0$), increasing step count $n$ reduces algorithmic error ($\mathcal{O}(1/n)$) but accumulates physical device noise ($\sim (1 - p_2)^{8n}$). Under the official noise model, there exists an **optimal step count $n^*$** maximizing simulation fidelity.

### Summary Table of Results

| Target / Benchmark | Processor | Noise Scale | Optimal $n^*$ | Peak Fidelity $F(n^*)$ | Routing SWAP Overhead | Routing 2Q Gates |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Tree Reference** | **Processor A** | $1.0\times$ | **6** | **0.577419** | High (Tree graph) | +21 CX at $n=4$ |
| **Heavy-Hex Reference** | **Processor B** | $1.0\times$ | **7** | **0.662816** | 0 (Contains 4-cycle) | 0 extra CX |
| **Proposed Architecture** | **Processor D** | $0.5\times$ | **8** | **0.782118** | **0 (Native 4-cycle)** | **0 extra CX** |
| **Proposed Architecture** | **Processor D** | $1.0\times$ | **7** | **0.662816** | **0 (Native 4-cycle)** | **0 extra CX** |
| **Proposed Architecture** | **Processor D** | $2.0\times$ | **6** | **0.503076** | **0 (Native 4-cycle)** | **0 extra CX** |

### Bell-State Verification on Processor D (Objective 1)
- **Directly Connected Pair $(0, 1)$ ($d = 1$ hop):** Ideal $F = 1.0000 \pm 0.0000$, Noisy $F = 0.9633 \pm 0.0021$ ($1$ CX gate).
- **Diameter Distant Pair $(0, 8)$ ($d = 4$ hops):** Ideal $F = 1.0000 \pm 0.0000$, Noisy $F = 0.9420 \pm 0.0026$ ($10$ CX gates across 3 SWAP stages).

### Bonus: Quantum Error Mitigation at $n^* = 7$ on Processor D
- **Exact Reference $M_s(T=2.0)$:** $0.652444$
- **Raw Noisy $M_s$ (Baseline $1.0\times$):** $0.439087 \pm 0.010499$ (Absolute error: $0.213357$)
- **Zero-Noise Extrapolated $M_s$ (Richardson ZNE):** **$0.601440$** (**$76.1\%$ error reduction**)
- **Combined ZNE + Readout Inversion:** **$0.652061$** (**$99.8\%$ recovery**, absolute error: $0.000383$)

---

## 📂 Repository File Layout

Following Section 9 of the official problem statement:

```
QFF2026_Simulation_Challenge_Kit/
├── README.md                      # Team info, PS code, execution guide, key results, AI disclosure
├── main.ipynb                     # Master notebook running everything end-to-end in <2 min
├── requirements.txt               # Pinned official package versions
├── qff_sim.py                     # Official starter module (strictly unchanged)
├── processors/
│   └── processor_D.json           # Validated 8-qubit dual-ring custom architecture
├── results/
│   ├── figures/
│   │   ├── fig1_processor_topologies.png
│   │   ├── fig2_bell_state_fidelity.png
│   │   ├── fig3_trotter_fidelity_crossover.png
│   │   ├── fig4_noise_sweep_shift.png
│   │   ├── fig5_gate_split_and_depth.png
│   │   └── fig6_magnetization_and_mitigation.png
│   └── tables/
│       ├── bell.csv               # Objective 1 Bell state correlations and fidelities
│       ├── resources.csv          # Two-qubit gate counts, depth, and SWAP overheads
│       ├── s3_scan.csv            # Complete 60-row step scan across A, B, D and noise scales
│       ├── s3_summary.csv         # Summary table of n* and peak fidelities
│       └── mitigation.csv         # Bonus: Readout and ZNE mitigation metrics
├── report.md                      # Comprehensive 8-page final technical report
├── run_all.py                     # Headless automated CLI runner
└── tools/
    └── check_processor_d.py       # Official Processor D rule checker
```

---

## ⚡ Quickstart & Reproducibility Guide

The entire simulation suite runs in **under 10 seconds** on any standard laptop CPU with locked random seeds (`SEED = 42`).

### 1. Setup Virtual Environment
```bash
python3 -m venv sim-env
source sim-env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### 2. Validate Processor D
```bash
python tools/check_processor_d.py processors/processor_D.json 4
# Output: OK - Processor D satisfies all design rules.
```

### 3. Run All Benchmarks via Headless Runner
```bash
python run_all.py
```
This regenerates all CSV tables in `results/tables/` and all publication figures in `results/figures/`.

### 4. Run the Master Jupyter Notebook
```bash
jupyter notebook main.ipynb
# or execute directly in headless mode:
jupyter nbconvert --to notebook --execute --inplace main.ipynb
```

---

## 🤖 AI Tool Disclosure

In accordance with Section 9 of the QFF 2026 Phase 2 Participant Rules, we declare:
- **Tools Used:** Antigravity AI coding assistant (Google DeepMind).
- **Scope of Use:** Automated scaffolding of plotting routines, mathematical LaTeX formatting, notebook cell generation, and documentation formatting.
- **Verification:** All quantum simulation algorithms, circuits, noise evaluations, and numerical data were verified independently using Qiskit 1.4.6, Qiskit Aer 0.18.0, and standard linear algebra references.
