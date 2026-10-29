<div align="center">

# ⚛️ StepWise-Q: Optimal Trotter Step Prediction & Hardware-Aware Co-Design

### *Balancing Algorithmic Trotter Error Against Physical Gate Noise on NISQ Architectures*

[![Competition](https://img.shields.io/badge/Competition-Qiskit_Fall_Fest_2026-6929C4?style=for-the-badge&logo=qiskit&logoColor=white)](https://qiskit.org/)
[![Track](https://img.shields.io/badge/Track-3:_Quantum_Simulation-0062FF?style=for-the-badge)](https://qiskit.org/)
[![Problem](https://img.shields.io/badge/Problem-PS--S3:_Trotter_Noise_Accumulation-198038?style=for-the-badge)](https://qiskit.org/)
[![Live 3D Lab](https://img.shields.io/badge/Live_Demo-3D_Quantum_Lab-00E5FF?style=for-the-badge&logo=three.js&logoColor=black)](https://vellorpavan.github.io/Trotter_error/)
[![Phase 2 Report](https://img.shields.io/badge/Report-Phase_2_Technical_Paper-0062FF?style=for-the-badge&logo=read-the-docs&logoColor=white)](report.md)
[![Phase 1 Folder](https://img.shields.io/badge/Phase_1-Directory_Archive-FF9800?style=for-the-badge&logo=github&logoColor=white)](Phase_1/)
[![Phase 1 Report](https://img.shields.io/badge/Report-Phase_1_Report-10B981?style=for-the-badge&logo=read-the-docs&logoColor=white)](Phase_1/report.md)
[![Tests](https://img.shields.io/badge/Tests-14%2F14_Passing-059669?style=for-the-badge&logo=checkmarx&logoColor=white)](tests/test_all.py)
[![Python](https://img.shields.io/badge/Python-3.10_%7C_3.11_%7C_3.12-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)

<p align="center">
  <b>Team:</b> <code>StepWise-Q</code> &nbsp;|&nbsp; 
  <b>Institution:</b> <code>SRM University-AP</code> &nbsp;|&nbsp; 
  <b>Problem Statement:</b> <code>PS-S3</code>
</p>

[Phase 2 Report (report.md)](report.md) •
[Phase 1 Folder (Phase_1/)](Phase_1/) •
[Phase 1 Report](Phase_1/report.md) •
[Key Results](#-executive-dashboard--key-results) •
[Section 7 Definitions](#-compliance-with-official-definitions-section-7) •
[Processor D Co-Design](#-processor-topologies--processor-d-co-design) •
[Visual Results](#-visual-experimental-findings) •
[Bell Experiment](#-common-phase-2-bell-state-verification) •
[Error Mitigation](#-bonus-quantum-error-mitigation) •
[3D Simulation Lab](#-interactive-3d-quantum-simulation-laboratory) •
[Quickstart](#-quickstart--reproduction-guide)

</div>

---

## 📑 Technical Reports & Documentation

All technical documentation and deliverables for both phases are accessible directly in this repository:

| Deliverable | Location | Description |
|:---|:---:|:---|
| 📘 **Phase 2 Technical Report** | [**`report.md`**](report.md) | Official 8-page final paper: Bell-state verification, Processor D co-design, Trotter crossover derivations, and ZNE mitigation. |
| 📗 **Phase 1 Technical Report** | [**`Phase_1/report.md`**](Phase_1/report.md) | Phase 1 scientific report: analytical error scaling laws, reference benchmarks, and convergence limits. |
| 📕 **Phase 1 System Architecture** | [**`Phase_1/Architecture.md`**](Phase_1/Architecture.md) | Theoretical architecture specification and experimental framework documentation. |
| 📑 **Phase 1 Architecture PDF** | [**`Phase_1/ARCHITECTURE.md.pdf`**](Phase_1/ARCHITECTURE.md.pdf) | Publication-formatted PDF of the Phase 1 architectural specifications. |
| 📁 **Phase 1 Codebase & Tests** | [**`Phase_1/`**](Phase_1/) | Standalone Phase 1 repository root containing all original test suites, data, and notebooks. |

---

## 📊 Executive Dashboard & Key Results

In digital quantum simulation (DQS) of continuous Hamiltonian dynamics, time evolution $e^{-i H T}$ is discretized into $n$ discrete Trotter steps. On noisy intermediate-scale quantum (NISQ) devices, this creates a fundamental **physical crossover**:

$$\underbrace{\epsilon_{\text{Trotter}} \sim \mathcal{O}\left(\frac{T^2}{n}\right)}_{\text{Algorithmic Discretization Error (dominant at small } n\text{)}} \quad \longleftrightarrow \quad \underbrace{F_{\text{noise}} \approx (1 - p_2)^{N_{\text{2Q}}(n)}}_{\text{Entangling Gate Decoherence (dominant at large } n\text{)}}$$

There exists a unique **optimal step count $n^*$** that maximizes total noisy state fidelity $F(\rho_{\text{noisy}}(n), |\psi_{\text{exact}}\rangle)$.

```mermaid
flowchart LR
    subgraph RegimeLow["Algorithmic Regime (Small n)"]
        direction TB
        L1["Large Step Size dt = T / n"] --> L2["BCH Non-Commutativity"]
        L2 --> L3["Dominant Trotter Error O(1/n)"]
    end

    subgraph RegimeHigh["Decoherence Regime (Large n)"]
        direction TB
        H1["Small Step Size dt"] --> H2["High Gate Count (8n CX + SWAPs)"]
        H2 --> H3["Exponential Gate Noise (1 - p2)^N"]
    end

    L3 --> Crossover["Optimal Step Count n* : Peak Noisy Fidelity F(n*)"]
    H3 --> Crossover
```

### 🏆 Benchmark Scorecard

| Target Processor | Architecture Class | Noise Scale | Optimal Steps ($n^*$) | Peak Fidelity $F(n^*)$ | Ideal Fidelity $F_{\text{ideal}}$ | Fidelity Gap $\Delta F$ | Two-Qubit Gates ($n^*$) | Routing SWAP Overhead |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Processor A** | 5-Qubit Tree | $1.0\times$ | **6** | **0.577419** | 0.963479 | 0.386060 | 79 (48L + 31S) | 11 SWAPs |
| **Processor B** | 9-Qubit Heavy-Hex | $1.0\times$ | **7** | **0.662816** | 0.972986 | 0.310170 | 56 (56L + 0S) | **0 SWAPs** |
| **Processor D** *(Ours)* | **9-Qubit 3×3 Grid** | **$0.5\times$** | **8** | **0.782118** | 0.979246 | 0.197128 | 64 (64L + 0S) | **0 SWAPs** |
| **Processor D** *(Ours)* | **9-Qubit 3×3 Grid** | **$1.0\times$** | **7** | **0.662816** | 0.972986 | 0.310170 | 56 (56L + 0S) | **0 SWAPs** |
| **Processor D** *(Ours)* | **9-Qubit 3×3 Grid** | **$2.0\times$** | **6** | **0.503076** | 0.963479 | 0.460403 | 48 (48L + 0S) | **0 SWAPs** |

> [!NOTE]  
> **Simulation Scope:** All calculations represent high-fidelity classical simulations via `qiskit-aer` using full density-matrix integration for gate noise and 4,096-shot Monte Carlo sampling for measurement statistics under the official QFF 2026 noise model.

---

## 📐 Compliance with Official Definitions (Section 7)

All calculations in this project strictly follow the official definitions locked in Section 7 of the problem statement specification:

| Official Quantity | Official Definition (Section 7) | Implementation in StepWise-Q | Code Function / Source |
|:---|:---|:---|:---|
| **Energy error** | $E_{\text{VQE}} - E_0$ with $E_0$ from exact diagonalisation (also given in S1 files) | Supported in `qff_sim.py` and benchmarked in `data/s1/`. *(Specific to PS-S1 VQE; for PS-S3, target observable is $M_s$ and state fidelity).* | `qff_sim.py` |
| **State fidelity (pure)** | $F = \|\langle \psi_{\text{exact}} \| \psi \rangle\|^2$ | Evaluates overlap between pure Trotter statevector and exact analytical evolution $\|\psi_{\text{exact}}(T)\rangle$. Populates `fidelity_trotter_ideal` in `results/tables/s3_scan.csv`. | `qiskit.quantum_info.state_fidelity` |
| **State fidelity (noisy)** | $F = \langle \psi_{\text{exact}} \| \rho \| \psi_{\text{exact}} \rangle$, with $\rho$ from density-matrix simulation of transpiled circuit under official noise model (**gate noise only; no readout**). `noisy_state_fidelity()` computes it and handles qubit re-mapping. | Computed via `noisy_state_fidelity()` with `include_readout=False`, mapping physical qubits back to logical qubits via initial layout inversion. Governs the selection of $n^*$. | `qff_sim.noisy_state_fidelity` |
| **Trotter state error** | $\|\|\psi_{\text{Trotter}}(T)\rangle - \|\psi_{\text{exact}}(T)\rangle\|\|$ (2-norm difference), via `state_error()` | Computed via `qff_sim.state_error(psi_a, psi_b)` evaluating the Euclidean 2-norm distance. In PS-S3, fidelity is related via $F = 1 - \frac{1}{2}\epsilon^2$. | `qff_sim.state_error` |
| **Exact evolution** | $\|\psi(t)\rangle = e^{-iHt}\|\psi(0)\rangle$ computed with `exact_evolution()` (`scipy expm_multiply`) | Solves exact unitary dynamics at $T=2.0$ from initial Néel state $\|1010\rangle$, providing the baseline target $\|\psi_{\text{exact}}\rangle$ and $M_s^{\text{exact}} = 0.652444$. | `qff_sim.exact_evolution` |
| **Z-basis observable from shots** | `run_z_measurement()`: 4,096 shots, mean and standard error | Executes transpiled circuit with 4,096 shots under target noise model. Computes expectation value $\langle Z \rangle$ and standard error $\text{SE} = \sigma / \sqrt{4096}$. Populates `ms_noisy` and `ms_noisy_se`. | `qff_sim.run_z_measurement` |
| **Optimizer calls** | Number of energy evaluations, counted with `CountedCost` | Implemented via `CountedCost` wrapper in `qff_sim.py`. *(Used in variational algorithms; PS-S3 is a direct unitary Trotter sweep without classical optimization loops).* | `qff_sim.CountedCost` |

---

## 🔬 Physical System & Hamiltonian Formulation

### Target Physical Model: 1D Ring TFIM
The physical benchmark is the 1D Transverse-Field Ising Model on an $N=4$ periodic spin ring:

$$H = -J \sum_{i=0}^{3} Z_i Z_{i+1} - g \sum_{i=0}^{3} X_i, \quad \text{with periodic boundary } Z_4 \equiv Z_0$$

```
                   (0) ---------- (1)
                    |              |      J = 1.0 (Nearest-Neighbor Coupling)
                    |    N = 4     |      g = 1.0 (Transverse Field)
                    |     Ring     |      T = 2.0 (Evolution Time)
                   (3) ---------- (2)
```

- **Evolution Parameters:** $J = 1.0$, $g = 1.0$, total duration $T = 2.0$.
- **Initial State:** Antiferromagnetic Néel state $|\psi_0\rangle = |1010\rangle$ (computational basis index 10).
- **Exact Reference State:** $|\psi_{\text{exact}}(T)\rangle = e^{-i H T} |\psi_0\rangle$, solved via matrix exponentiation (`scipy.linalg.expm_multiply`).
- **Target Observable:** Staggered magnetization:

$$M_s = \frac{1}{4} \sum_{i=0}^{3} (-1)^i \langle Z_i \rangle = \frac{1}{4} \left( \langle Z_0 \rangle - \langle Z_1 \rangle + \langle Z_2 \rangle - \langle Z_3 \rangle \right)$$

$$\mathbf{M_s(T=2.0) = 0.65244449 \approx 0.652444}$$

### First-Order Trotter Decomposition ($S_1$)
For step count $n$ and sub-interval $\Delta t = T / n$:

$$e^{-i H \Delta t} \approx \left(\prod_{(i, j) \in \text{bonds}} e^{i J \Delta t Z_i Z_j}\right) \left(\prod_{i=0}^{3} e^{i g \Delta t X_i}\right)$$

1. **Locked Bond Sequence:** Evaluated in strict sequential order: $(0, 1) \to (1, 2) \to (2, 3) \to (3, 0)$.
2. **$ZZ(\theta)$ Synthesis:** Each bond with angle $\theta = -2 J \Delta t$ executes `cx(i, j) -> rz(theta, j) -> cx(i, j)`.
3. **Transverse Layer:** Four single-qubit rotations $R_x(-2 g \Delta t)$ applied simultaneously.
4. **Logical Complexity:** Exactly $8n$ CNOT gates and $4n$ single-qubit rotations.

### End-to-End Simulation Pipeline

```mermaid
flowchart TD
    subgraph Inputs["1. Problem Inputs (Locked Challenge Kit)"]
        direction LR
        H["Ring TFIM (N=4, J=1.0, g=1.0)"] --- Psi0["Neel State |1010>"] --- T["Time T = 2.0"]
    end

    subgraph SimulationPipeline["2. Quantum Simulation Engine (StepWise-Q)"]
        direction TB
        Exact["Exact Reference: scipy expm_multiply"]
        Trotter["First-Order Trotter Circuit Synthesis U1(dt)"]
        Transpile["Transpilation & SABRE Routing (Processors A, B, D)"]
        AerSim["Qiskit Aer: Density-Matrix & 4096-Shot Simulator"]
        Exact & Trotter --> Transpile --> AerSim
    end

    subgraph Evaluation["3. Objectives & Scientific Analysis"]
        direction TB
        Obj1["Objective 1: Bell-State Verification (d=1 vs d=4)"]
        Obj2["Objective 2: Processor D 3x3 Grid Co-Design Validation"]
        Obj3["Objective 3: Trotter Step Scan (n=1..16) & Crossover n*"]
        Sweep["Noise Sweep: 0.5x, 1.0x, 2.0x on Processor D"]
        Bonus["Bonus: Readout Inversion & Richardson ZNE Mitigation"]
    end

    subgraph Deliverables["4. Verified Hackathon Deliverables"]
        direction LR
        Tables["5 CSV Result Tables (results/tables/)"]
        Figures["6 300-DPI Publication Figures (results/figures/)"]
        Report["8-Page Technical Paper (report.md)"]
    end

    Inputs --> SimulationPipeline
    SimulationPipeline --> Evaluation
    Evaluation --> Deliverables
```

---

## 📐 Processor Topologies & Processor D Co-Design

```
   Processor A (5 Qubits, Tree)         Processor B (9 Qubits, Heavy-Hex)      Processor D (9 Qubits, 3×3 Grid)
             (0)                               (0) --- (1) --- (2)                 (0) --- (1) --- (2)
            /   \                               |               |                   |       |       |
          (1)   (2)                            (3) --- (4) --- (5)                 (3) --- (4) --- (5)
          / \                                   |               |                   |       |       |
        (3) (4)                                (6) --- (7) --- (8)                 (6) --- (7) --- (8)
   • Diameter: 4, Max Degree: 3            • Contains ONE 4-cycle (0-1-4-3)    • Contains FOUR native 4-cycles
   • Zero cycles -> SWAPs mandatory        • 0 SWAPs for Ring Embedding        • 0 SWAPs for Ring Embedding
   • +31 CX overhead at n=6                • Non-symmetric peripheral nodes    • Planar 2D mesh, Max Degree 4
```

### Architectural Analysis
- **Processor A (Tree Graph):** Lacks physical cycles. To execute the periodic closure term $(3, 0)$, routing requires dynamic SWAP insertions through root node 0, introducing +21 CXs at $n=4$ and +31 CXs at $n=6$. This steep depth explosion penalizes state fidelity.
- **Processor B (Heavy-Hex Unit Cell):** Embeds the ring natively on $(0, 1, 4, 3)$ with 0 SWAPs, but has non-uniform degree (nodes 0, 2, 6, 8 have degree 2; nodes 1, 3, 4, 5, 7 have degree 3).
- **Processor D (Co-Designed 3×3 Grid):**
  - **Node Budget:** 9 physical qubits ($N \le 9$).
  - **Degree Constraint:** Maximum degree 4 at center node 4 ($\le 4$).
  - **Planar Symmetry:** Contains four native 4-cycles: `(0,1,4,3)`, `(1,2,5,4)`, `(3,4,7,6)`, `(4,5,8,7)`.
  - **Routing Overhead:** **Zero SWAP gates** across all Trotter step counts ($n \in [1, 16]$).
  - Formally validated by `tools/check_processor_d.py`.

---

## 📈 Visual Experimental Findings

<div align="center">

### Topology Architecture & Ring Allocation
![Figure 1: Processor Topologies](results/figures/fig1_processor_topologies.png)
*Figure 1: Coupling maps for Processors A, B, and D. Green nodes denote the physical allocation of the 4-qubit ring Hamiltonian. Processors B and D embed the periodic ring with zero routing overhead.*

---

### Bell-State Parity & Fidelity (Objective 1)
![Figure 2: Bell State Fidelity](results/figures/fig2_bell_state_fidelity.png)
*Figure 2: Parity correlation components $\langle XX \rangle, \langle YY \rangle, \langle ZZ \rangle$ and Bell state fidelity under ideal and noisy conditions for adjacent pair $(0, 1)$ ($d=1$) versus diameter pair $(0, 8)$ ($d=4$).*

---

### PS-S3 Trotter Fidelity Crossover
![Figure 3: Trotter Fidelity Crossover](results/figures/fig3_trotter_fidelity_crossover.png)
*Figure 3: Noisy state fidelity vs Trotter step count $n \in [1, 12]$. At small $n$, Trotter error dominates; at large $n$, gate decoherence dominates. The optimum balances both at $n^* = 6$ (A) and $n^* = 7$ (B & D).*

---

### Noise Sweep Shift on Processor D
![Figure 4: Noise Sweep Shift](results/figures/fig4_noise_sweep_shift.png)
*Figure 4: Response of optimal step count $n^*$ to physical noise scale on Processor D. As noise increases from $0.5\times \to 1.0\times \to 2.0\times$, the optimal step count shifts monotonically downward ($n^* = 8 \to 7 \to 6$).*

---

### Two-Qubit Gate Scaling & Circuit Depth
![Figure 5: Gate Split & Depth](results/figures/fig5_gate_split_and_depth.png)
*Figure 5: Circuit depth and two-qubit gate breakdown. Processor A incurs severe SWAP routing penalties (+31 CXs at $n=6$), whereas Processors B and D scale cleanly at the theoretical lower bound ($8n$ CNOTs).*

---

### Staggered Magnetization & Error Mitigation
![Figure 6: Magnetization & Mitigation](results/figures/fig6_magnetization_and_mitigation.png)
*Figure 6: Staggered magnetization $M_s(T=2.0)$ dynamics. Richardson Zero-Noise Extrapolation recovers $76.1\%$ of the error; combined ZNE + readout mitigation achieves $99.8\%$ recovery toward the exact benchmark ($0.652444$).*

</div>

---

## 🔔 Common Phase 2 Bell-State Verification

To benchmark the physical quantum channel quality of Processor D, we prepare and verify the maximally entangled Bell state $|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$:

$$F = \frac{1 + \langle XX \rangle - \langle YY \rangle + \langle ZZ \rangle}{4}$$

Each expectation value is sampled with **4,096 shots** per basis. Standard errors are computed via binomial sampling: $\text{SE} = \sqrt{(1 - \langle P \rangle^2)/N_{\text{shots}}}$.

### Official Results (`results/tables/bell.csv`)

| Qubit Pair | Distance ($d$) | Entangling Gates | Condition | $\langle XX \rangle$ | $\langle YY \rangle$ | $\langle ZZ \rangle$ | State Fidelity ($F$) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **$(0, 1)$** | 1 hop | 1 CX | **Ideal** | $+1.0000 \pm 0.0000$ | $-1.0000 \pm 0.0000$ | $+1.0000 \pm 0.0000$ | **$1.0000 \pm 0.0000$** |
| **$(0, 1)$** | 1 hop | 1 CX | **Noisy** | $+0.9458 \pm 0.0051$ | $-0.9521 \pm 0.0048$ | $+0.9551 \pm 0.0046$ | **$0.9633 \pm 0.0021$** |
| **$(0, 8)$** | 4 hops | 10 CX (3 SWAP) | **Ideal** | $+1.0000 \pm 0.0000$ | $-1.0000 \pm 0.0000$ | $+1.0000 \pm 0.0000$ | **$1.0000 \pm 0.0000$** |
| **$(0, 8)$** | 4 hops | 10 CX (3 SWAP) | **Noisy** | $+0.9219 \pm 0.0061$ | $-0.9077 \pm 0.0066$ | $+0.9385 \pm 0.0054$ | **$0.9420 \pm 0.0026$** |

---

## 🛡️ Bonus: Quantum Error Mitigation

At optimal operating point $n^* = 7$ on Processor D, physical noise depresses raw magnetization from $0.652444 \to 0.439087$. We apply two complementary mitigation techniques:

1. **Readout Error Inversion (Tensored Matrix):**
   With single-qubit readout error $p_{\text{ro}} = 0.01$, expectation values attenuate by $(1 - 2 p_{\text{ro}})^4 \approx 0.92237$. Inverting this tensor product scales raw $M_s$ to $0.476043$ (**$17.3\%$ error reduction**, 8,192 shots).
2. **Zero-Noise Extrapolation (Richardson Polynomial):**
   Evaluating $M_s(\lambda)$ across noise scales $\lambda \in \{0.5, 1.0, 2.0\}$ and fitting a quadratic polynomial yields zero-noise intercept:
   
   $$M_s^{\text{ZNE}} = 0.601440 \quad (\mathbf{76.1\% \text{ Error Reduction}}, 8,192\text{ shots})$$

3. **Combined ZNE + Readout Mitigation:**
   Applying readout correction to the Richardson zero-noise estimate yields:
   
   $$M_s^{\text{combined}} = \frac{M_s^{\text{ZNE}}}{(1 - 2 p_{\text{ro}})^4} = \mathbf{0.652061} \quad (\text{Target: } 0.652444, \mathbf{99.8\% \text{ Recovery}})$$

### Mitigation Performance (`results/tables/mitigation.csv`)

| Mitigation Method | Extra Shots | Raw $M_s$ | Mitigated $M_s$ | Target $M_s^{\text{exact}}$ | Absolute Error | Error Reduction (%) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Baseline (Noisy $1.0\times$)** | 0 | 0.439087 | — | 0.652444 | 0.213357 | $0.0\%$ (Reference) |
| **Readout Error Mitigation** | 8,192 | 0.439087 | 0.476043 | 0.652444 | 0.176401 | **17.3%** |
| **Richardson ZNE** | 8,192 | 0.439087 | 0.601440 | 0.652444 | 0.051004 | **76.1%** |
| **Combined (ZNE + Readout)** | 16,384 | 0.439087 | **0.652061** | 0.652444 | **0.000383** | **99.8%** |

---

## 🎮 Interactive 3D Quantum Simulation Laboratory

To complement our static benchmarking figures, StepWise-Q includes a **hardware-accelerated 3D WebGL laboratory** built with **Three.js** and **React**. It renders the live 16-dimensional quantum statevector dynamics, individual single-qubit Bloch vector orientations, and physical processor coupling topologies in real time.

<div align="center">
  <iframe src="https://gemini.google.com/share/741c9c5cceca?social_preview=false&skid=6c798fbe-3abe-466e-aaea-d979cbf5dad4" width="100%" height="600" frameborder="0" allowfullscreen style="border-radius: 8px; border: 1px solid #334155; max-width: 100%; margin-bottom: 12px;"></iframe>
  <br/>
  <a href="https://gemini.google.com/share/741c9c5cceca?social_preview=false&skid=6c798fbe-3abe-466e-aaea-d979cbf5dad4" target="_blank">
    <img src="https://img.shields.io/badge/Launch_3D_Quantum_Lab-Open_in_New_Tab-00E5FF?style=for-the-badge&logo=google&logoColor=black" alt="Launch 3D Quantum Lab" />
  </a>
</div>

```
+-----------------------------------------------------------------------------------------+
|                  StepWise-Q 3D WebGL Quantum Simulation Laboratory                      |
|                                                                                         |
|   [ 3D Quantum Viewport (Three.js) ]           [ 16-Basis Statevector Amplitudes ]       |
|    • Interactive Orbit, Zoom, Pan Controls      • Live |0000> ... |1111> Probabilities   |
|    • Real-time Transparent Bloch Spheres        • Neel State |1010> Decay Tracking      |
|    • Live Bloch Vector (X, Y, Z) and Purity     • Staggered Magnetization Realtime Ms   |
|    • Processors A (Tree), B (Hex), D (Grid)     • Projective Shot Readout Emulation     |
|                                                                                         |
|   [ Circuit Timeline Scrubber & Step Player: Step 1 .. n -> Final Readout ]              |
+-----------------------------------------------------------------------------------------+
```

### 🌟 Key Interactive Capabilities
1. **Real-Time 3D Bloch Spheres:** Each qubit node renders a 3D transparent sphere with a live Bloch vector arrow $\vec{r} = (\langle X \rangle, \langle Y \rangle, \langle Z \rangle)$, length $|\vec{r}|$, and purity $\text{Tr}(\rho^2)$.
2. **Dynamic Topology Switching:** Instantly toggle between **Processor A** (5-qubit tree), **Processor B** (9-qubit heavy-hex), and **Processor D** (9-qubit 3×3 grid with four 4-cycles) to inspect spatial routing differences.
3. **16-State Amplitude Bar Chart:** High-resolution histogram displaying all computational basis probabilities $|0000\rangle \dots |1111\rangle$ and quantum phase angles as the state evolves from Néel $|1010\rangle$ under the transverse field.
4. **Playback & Scrubbing Timeline:** Scrub through Trotter evolution steps ($t \in [0, 2.0]$) with variable playback speeds ($0.5\times, 1.0\times, 2.0\times, 4.0\times$), stepping forward or backward gate by gate.
5. **Projective Readout Emulation:** Sample 4,096-shot projective measurements with simulated readout error directly in the browser.

### 🚀 How to Launch the 3D Simulation

- **Option 1 (Interactive Cloud / Embedded View):**
  Open the interactive simulation directly in full screen:
  👉 **[Interactive 3D Quantum Simulation Lab](https://gemini.google.com/share/741c9c5cceca?social_preview=false&skid=6c798fbe-3abe-466e-aaea-d979cbf5dad4)**

- **Option 2 (GitHub Pages Live Deployment):**
  Explore the deployed zero-install WebGL lab:
  👉 **[https://vellorpavan.github.io/Trotter_error/](https://vellorpavan.github.io/Trotter_error/)**

- **Option 3 (1-Click Local Browser Launch):**
  ```bash
  open visualization/index.html
  ```
  *(Or double-click `visualization/index.html` in macOS Finder / Windows Explorer).*

- **Option 4 (Local HTTP Server):**
  ```bash
  python -m http.server 8000 --directory visualization
  # Open http://localhost:8000 in your browser
  ```

- **Option 5 (React/Vite Component Source):**
  The modular React JSX component source is stored in [`src/Interactive3DSimulation.jsx`](src/Interactive3DSimulation.jsx) for integration into React, Vite, or Next.js applications.

---

## 🗂️ Repository Architecture & Phase 1 Archive

```
.
├── README.md                      # Primary hackathon documentation & visual dashboard
├── report.md                      # Comprehensive 8-page scientific technical paper (Phase 2)
├── main.ipynb                     # Executable master notebook (recomputes end-to-end in ~7s)
├── run_all.py                     # Headless standalone CLI reproduction script (~4s)
├── qff_sim.py                     # Official simulation module (byte-for-byte identical to kit)
├── requirements.txt               # Pinned package dependencies
├── visualization/
│   └── index.html                 # Self-contained 3D WebGL laboratory (1-click browser run)
├── src/
│   └── Interactive3DSimulation.jsx# Full React 18 + Three.js interactive component
├── processors/
│   └── processor_D.json           # Objective 2 co-designed 9-qubit 3x3 grid processor
├── results/
│   ├── figures/                   # 6 high-resolution publication figures (PNG)
│   └── tables/                    # 5 official CSV benchmark tables (zero NaN values)
├── tests/
│   └── test_all.py                # 14 automated unit tests covering all 5 objectives
├── tools/
│   └── check_processor_d.py       # Official verification checker for Processor D
├── starter/                       # Clean copy of original challenge-kit distribution
│   └── KIT_README.md              # Original starter README (matches MANIFEST.sha256)
└── Phase_1/                       # Complete Phase 1 Archive (Code, Tests, Reports)
    ├── report.md                  # Phase 1 comprehensive scientific report
    ├── Architecture.md            # Phase 1 architecture & model design documentation
    ├── ARCHITECTURE.md.pdf        # Rendered Phase 1 architecture paper
    ├── main.ipynb                 # Phase 1 executable notebook
    ├── results/                   # Phase 1 benchmark tables & figures
    └── tests/                     # Phase 1 test suite (8 test suites)
```

---

## 🚀 Quickstart & Reproduction Guide

### 1. Environment Setup
```bash
# Clone repository and enter directory
cd project

# Create dedicated virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Upgrade pip and install pinned dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

### 2. Validate Processor D Architecture
```bash
python tools/check_processor_d.py processors/processor_D.json 4
```
```
OK - Processor D satisfies all design rules.
```

### 3. Run Automated Unit Test Suite
```bash
python -m unittest tests/test_all.py
```
```
..............
----------------------------------------------------------------------
Ran 14 tests in 0.005s

OK
```

### 4. Headless Reproduction (All Tables & Figures in ~4s)
```bash
python run_all.py
```

### 5. Interactive Master Notebook
```bash
jupyter notebook main.ipynb
# or headless via nbconvert:
jupyter nbconvert --to notebook --execute --inplace main.ipynb
```

---

## 📑 Data Dictionary (`results/tables/`)

<details>
<summary><b>Click to expand CSV column specifications</b></summary>

| File | Rows | Schema & Key Columns |
|:---|:---:|:---|
| **`bell.csv`** | 4 | `pair, pair_distance, condition, two_qubit_gates, XX, XX_se, YY, YY_se, ZZ, ZZ_se, fidelity, fidelity_se` |
| **`resources.csv`** | 36 | `ps, circuit_label, processor, executable, logical_qubits, physical_qubits, logical_2q, two_q, depth, two_q_depth, swaps, routing_overhead_2q` |
| **`s3_scan.csv`** | 60 | `processor, noise_scale, steps, fidelity_noisy, fidelity_trotter_ideal, ms_noisy, ms_noisy_se, ms_exact, logical_2q, swap_2q, two_q` |
| **`s3_summary.csv`** | 5 | `processor, noise_scale, n_star, fidelity_at_n_star` |
| **`mitigation.csv`** | 3 | `method, ms_raw, ms_mitigated, extra_shots` |

</details>

---

## 🧪 Reproducibility, Assumptions & Limitations

- **Seed Discipline:** Fixed random seeds (`SEED = 42`) guarantee identical numerical output across runs.
- **Pure Python Dependencies:** Only standard scientific libraries (`qiskit`, `qiskit-aer`, `numpy`, `scipy`, `matplotlib`) are required.
- **Portability:** Pathing is dynamically relative via `os.path.dirname(os.path.abspath(__file__))`, running identically on macOS, Linux, and Windows/WSL2.
- **Simulation Assumptions:** Evaluated under classical Markovian noise models. Physical QPUs exhibit non-Markovian noise, spectator cross-talk, and calibration drift.

---

## 👥 Team & AI Disclosure

### Team StepWise-Q
- **Institution:** SRM University-AP
- **Track:** Track 3 (Quantum Simulation) · Problem Statement PS-S3
- **Venue:** Qiskit Fall Fest 2026 Phase 2 Offline Hackathon

### Generative AI & Tooling Disclosure
In accordance with Section 9 of the QFF 2026 Phase 2 Participant Rules:
- **Tools Used:** Antigravity AI coding assistant (Google DeepMind).
- **Scope of Assistance:** Scaffolding plotting scripts (Matplotlib themes), LaTeX formatting, Markdown report structuring, and test suite automation.
- **Verification:** All Hamiltonian constructions, circuit topologies, matrix exponentiations, noise sweeps, and physical conclusions were independently verified using Qiskit 1.4.6, Qiskit Aer 0.18.0, and NumPy/SciPy linear algebra routines.
