# PS-S3: Noise Accumulation in Trotter Circuits
### Final Technical Report — StepWise-Q
**Qiskit Fall Fest 2026 | Track 3: Quantum Simulation | SRM University-AP**
**Authors:** Jagadesh ([@jagadesh-jpg](https://github.com/jagadesh-jpg)) & Vellor Pavan ([@Vellorpavan](https://github.com/Vellorpavan))

---

## 🧭 Executive Summary

In Trotterized quantum simulation of interacting many-body quantum systems, algorithmic discretization error scales as $\mathcal{O}(\Delta t) = \mathcal{O}(T/n)$ for first-order Trotter-Suzuki decompositions. While increasing the step count $n$ systematically suppresses this discretization error, each additional step inserts a discrete layer of non-commuting unitary gates. On noisy intermediate-scale quantum (NISQ) devices, physical two-qubit gate errors ($p_2 = 0.005$) and single-qubit errors ($p_1 = 0.001$) compound exponentially with circuit depth. Consequently, the physical simulation fidelity $F(n) = \langle \psi_{\text{exact}} | \rho(n) | \psi_{\text{exact}} \rangle$ exhibits a strict convex peak at an **optimal step count** $n^*$.

This report details **StepWise-Q**'s comprehensive theoretical modeling, hardware co-design, and empirical benchmarking for **Problem Statement PS-S3 (Ring TFIM, $N=4, J=1.0, g=1.0, T=2.0$)** across Processors A, B, and our custom-designed **Processor D**:

1. **Objective 1 (Bell-State Verification):** Evaluated $|\Phi^+\rangle$ on Processor D. Direct nearest-neighbor pair $(0, 1)$ ($d=1$) achieves noisy state fidelity $F = 0.9633 \pm 0.0021$ ($1$ CX gate), where physical loss is predominantly driven by readout attenuation ($F_{\text{readout}} \approx 0.9703$). The diameter-separated pair $(0, 5)$ across $d=3$ hops achieves $F = 0.9471 \pm 0.0025$ ($7$ CX gates: 1 entangling + 2 SWAPs), following the exact $1 + 3(d-1)$ compilation pattern and demonstrating the fidelity penalty of non-local routing.
2. **Objective 2 (Processor D Proposal):** Co-designed a 6-qubit $2 \times 3$ ladder architecture ($\mathcal{E}_D = \{(0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)\}$) with maximum vertex degree 3 and lattice diameter 3. It natively embeds two symmetric 4-cycles ($\{0, 1, 4, 3\}$ and $\{1, 2, 5, 4\}$), eliminating all routing SWAPs for the periodic TFIM ring. While S3 metrics tie Processor B by construction (since any native 4-cycle yields identical S3 depth and gates), Processor D delivers decisive hardware advantages: bounded degree $\le 3$ (eliminating Processor B's degree-5 crosstalk hub) and a shorter diameter of 3 (requiring only 7 CXs for far-pair entanglement vs. 10 CXs on B).
3. **Objective 3 (Headline Step Scan & Crossover):**
   - **Processor A (Tree):** $n^* = 6$ with peak fidelity $F(n^*) = 0.577419$. The absence of cycles forces 11 SWAPs ($33$ extra CX gates, totaling 81 CXs) at $n=6$, shifting $n^*$ to lower depth to avoid fatal gate noise.
   - **Processor B (Heavy-Hex):** $n^* = 7$ with peak fidelity $F(n^*) = 0.662816$ ($56$ CX gates, $0$ SWAPs). Contains an asymmetric degree-5 hub prone to crosstalk.
   - **Processor D (6-Qubit Ladder):** $n^* = 7$ with peak fidelity $F(n^*) = 0.662816$ ($56$ CX gates, $0$ SWAPs). Matches the theoretical limit for 4-cycle embedding while maintaining degree $\le 3$.
4. **Noise Sweep Dynamics:** Across noise scaling factors $\lambda \in \{0.5, 1.0, 2.0\}$ on Processor D, the optimal step count shifts monotonically leftward:
   $$n^*(\lambda = 0.5) = 8 \; (F = 0.7821) \;\longrightarrow\; n^*(\lambda = 1.0) = 7 \; (F = 0.6628) \;\longrightarrow\; n^*(\lambda = 2.0) = 6 \; (F = 0.5031)$$
   confirming our theoretical model: higher physical error rates truncate allowable algorithmic depth.
5. **Bonus (Error Mitigation):** At $n^*=7$ on Processor D, raw noisy staggered magnetization $M_s = 0.4288 \pm 0.0106$ was mitigated via Readout Inversion ($M_s \to 0.4376 \pm 0.0108$, using the linear $0.98$ factor for single-qubit $Z$ sums) and Richardson Zero-Noise Extrapolation ($M_s \to 0.6205 \pm 0.0352$, evaluated with 8,192 calibration shots), achieving an honest **$85.7\%$ error reduction** toward the exact reference $M_s^{\text{exact}} = 0.652444$ without double-counting.
6. **Phase 1 Foundations & Continuity:** This Phase 2 report builds upon our Phase 1 analytical benchmark archive (located in directory [`Phase_1/`](Phase_1/) with dedicated report [`Phase_1/report.md`](Phase_1/report.md)), scaling theoretical prediction laws to full quantum density-matrix evolution, circuit transpilation, and hardware co-design.
7. **Interactive 3D WebGL Visualization:** An interactive hardware-accelerated 3D WebGL laboratory is provided at `visualization/index.html` (and live online at [vellorpavan.github.io/Trotter_error](https://vellorpavan.github.io/Trotter_error/)), enabling real-time 3D Bloch sphere inspection and dynamic topology exploration.

---

## 1. Mathematical Formulation & Theoretical Error Model

### 1.1 The Transverse-Field Ising Model (TFIM) on a Ring

The physical Hamiltonian under investigation is the 1D Transverse-Field Ising Model with periodic boundary conditions (ring topology):
$$H = -J \sum_{i=0}^{N-1} Z_i Z_{i+1} - g \sum_{i=0}^{N-1} X_i, \quad Z_N \equiv Z_0$$
with locked competition parameters $N = 4$, $J = 1.0$, $g = 1.0$, and total evolution time $T = 2.0$. The initial state is the Néel antiferromagnetic state:
$$|\psi(0)\rangle = |1010\rangle \quad \left(\text{Qiskit qubit ordering: } q_1, q_3 \in |1\rangle \implies Z = (+1, -1, +1, -1)\right)$$

The exact quantum state evolution $|\psi_{\text{exact}}(T)\rangle = e^{-i H T} |\psi(0)\rangle$ is evaluated numerically via Scipy's sparse matrix exponential (`scipy.sparse.linalg.expm_multiply`), yielding the exact reference staggered magnetization:
$$M_s(T) = \frac{1}{N} \sum_{i=0}^{N-1} (-1)^i \langle Z_i \rangle(T) \implies M_s^{\text{exact}}(T=2.0) = 0.65244449$$

### 1.2 First-Order Trotter Discretization

The continuous time evolution operator $U(T) = e^{-i H T}$ is approximated by $n$ discrete Trotter steps of duration $\Delta t = T/n$:
$$U_1(\Delta t) = \exp(-i H_X \Delta t) \exp(-i H_{ZZ} \Delta t)$$
Because the ring interaction terms do not all commute, the interaction Hamiltonian $H_{ZZ} = -J \sum_{i=0}^3 Z_i Z_{i+1}$ is decomposed into two non-overlapping, parallelizable bond layers:
1. **Even Bonds:** $(0, 1)$ and $(2, 3)$
2. **Odd Bonds:** $(1, 2)$ and $(3, 0)$

Each two-qubit bond evolution $\exp(i J \Delta t Z_a Z_b)$ compiles under the standard CX basis (`["rz", "sx", "x", "cx"]`) as:
$$R_{ZZ}(\theta) = \text{CX}(a, b) \cdot R_Z(b, \theta) \cdot \text{CX}(a, b), \quad \text{where } \theta = -2 J \Delta t$$
Consequently, each Trotter step requires exactly $4 \text{ bonds} \times 2 \text{ CX/bond} = 8 \text{ logical CX gates}$. For $n$ Trotter steps, the baseline logical two-qubit gate count is strictly:
$$N_{\text{2Q}}^{\text{logical}}(n) = 8n$$

### 1.3 The Two-Component Error Model for $n^*$

The total simulation error $\epsilon(n) = 1 - F(n)$ arises from two competing physical mechanisms:

$$\epsilon_{\text{total}}(n) \approx \epsilon_{\text{Trotter}}(n) + \epsilon_{\text{noise}}(n)$$

1. **Coherent Algorithmic Error:** From the Baker-Campbell-Hausdorff (BCH) expansion, the commutator $[H_{ZZ}, H_X] \neq 0$ induces a state error scaling as:
   $$\epsilon_{\text{Trotter}}(n) \approx \frac{T^2}{2n} \left\| [H_{ZZ}, H_X] \right\| = \frac{c_1}{n}$$
2. **Incoherent Gate Depolarizing Noise:** In a circuit with $N_{1Q}(n)$ single-qubit gates and $N_{2Q}(n)$ two-qubit gates, under depolarizing probabilities $p_1$ and $p_2$, the survival fidelity of the quantum state decays as:
   $$F_{\text{noise}}(n) \approx \prod_{g \in \text{gates}} (1 - p_g) \approx (1 - p_1)^{N_{1Q}(n)} (1 - p_2)^{N_{2Q}(n)} \approx \exp\left( - \lambda [N_{1Q}(n) p_1 + N_{2Q}(n) p_2] \right)$$
   For processors with topological routing overhead, $N_{2Q}(n) = 8n + 3 \cdot N_{\text{SWAP}}(n)$.

Taking the derivative $\frac{d F(n)}{dn} = 0$ yields the optimal step count $n^*$:
$$n^* \propto \sqrt{\frac{c_1}{\lambda \cdot p_2 \cdot (8 + 3 \alpha_{\text{SWAP}})}}$$
This analytical relationship predicts two central physical behaviors verified in our experiments:
- **Topology Dependency:** Increasing routing overhead ($\alpha_{\text{SWAP}} > 0$) increases the effective noise per step, forcing $n^*$ to a lower value.
- **Noise Scale Dependency:** Increasing the noise scale factor $\lambda$ shifts $n^*$ monotonically to smaller step counts ($n^* \propto 1/\sqrt{\lambda}$).

### 1.4 Compliance with Section 7 Official Definitions

In accordance with Section 7 of the QFF 2026 Phase 2 Problem Statements (*Definitions used by all simulation problem statements*), all quantitative metrics in this project strictly follow the locked official implementations:

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

## 2. Objective 2 — Processor D Proposal & Architectural Co-Design

### 2.1 Hardware Bottleneck Analysis on Reference Processors

![Figure 1: Processor Topologies](results/figures/fig1_processor_topologies.png)

1. **Processor A (5 Qubits, Tree Topology):**
   - Coupling edges: $(0,1), (1,2), (1,3), (3,4)$.
   - Maximum degree: 3 (node 1); **Diameter: 3** (longest shortest path is $0-1-3-4$).
   - **Fatal Bottleneck:** Processor A is a tree graph containing **zero cycles**. When compiling a periodic 4-spin ring TFIM Hamiltonian, the closing bond $(3, 0)$ cannot be executed natively. SABRE routing must insert SWAP gates (at least 1 to 2 SWAPs per step). At $n=4$, Processor A incurs 7 SWAPs ($21$ extra CX gates), inflating total two-qubit gates from 32 to 53 ($+65.6\%$ overhead) and total circuit depth from 45 to 66 ($+46.7\%$). At $n=6$, it requires 11 SWAPs ($33$ extra CXs, reaching 81 CXs vs 48 logical). This severe overhead accelerates noise accumulation and degrades peak fidelity.
2. **Processor B (9 Qubits, Heavy-Hex Inspired):**
   - Coupling edges: 12 edges connecting 9 qubits.
   - **Bottlenecks:** While Processor B contains a 4-cycle (qubits $\{2, 3, 4, 5\}$), node 6 has a high degree of **5** (connected to 1, 2, 5, 7, 8). In planar superconducting transmon architectures, high vertex degrees induce severe frequency collision (crowding) and parasitic $ZZ$ crosstalk. Furthermore, its lattice diameter is 4.

### 2.2 The Processor D Proposal: 6-Qubit $2 \times 3$ Ladder Architecture

Rather than relying on generic square grids, Processor D is deliberately co-designed as a **6-qubit $2 \times 3$ ladder architecture**, formalized in `processors/processor_D.json`:
- **Qubit Count:** 6 physical qubits (strictly within the allowed $[4, 12]$ range).
- **Topology:** $2 \times 3$ ladder mesh with 7 undirected coupling edges:
  $$\mathcal{E}_D = \{(0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)\}$$
- **Maximum Degree:** Strictly bounded to **3** (nodes 1 and 4 have degree 3; boundary nodes 0, 2, 3, 5 have degree 2). This eliminates high-degree hubs, drastically suppressing frequency crowding and parasitic $ZZ$ crosstalk.
- **Dedicated TFIM 4-Cycles:** The ladder contains two planar, symmetric 4-cycles: $\{0, 1, 4, 3\}$ and $\{1, 2, 5, 4\}$. The 4-spin ring TFIM maps natively onto either 4-cycle with **zero SWAPs and zero routing overhead**.
- **Lattice Diameter:** Strictly bounded to **3** (distance between opposite corners 0 and 5, or 2 and 3). This enables non-local 2-qubit operations across the chip diameter in just 7 CX gates ($1 + 3(d-1) = 7$), compared to 10 CX gates on the 9-qubit grid.
- **Basis Gates & Noise Model:** Pinned strictly to the official baseline: `["rz", "sx", "x", "cx"]`, $p_1 = 0.001$, $p_2 = 0.005$, $p_{\text{ro}} = 0.01$.

`processors/processor_D.json` was validated using `tools/check_processor_d.py`, passing all official design rules with zero errors.

### 2.3 Architectural Target vs. Measured Experimental Benchmarks

```
                    Processor D Graph Drawing (2×3 Ladder)
                           (0) ----------- (1) ----------- (2)
                            |               |               |
                            |   4-Cycle 1   |   4-Cycle 2   |
                            |               |               |
                           (3) ----------- (4) ----------- (5)
                    • Vertices: 6 qubits   • Edges: 7 undirected edges
                    • Max Degree: 3        • Graph Diameter: 3
                    • Planar ring embedding: (0-1-4-3) with 0 SWAPs
```

The table below contrasts our architectural targets with measured compiler and simulator telemetry:

| Benchmark Metric | Reference Tree (Proc A) | Architectural Target / Design Specification | Proposed Ladder (Proc D) Measured | Empirical Architectural Impact |
|:---|:---:|:---:|:---:|:---|
| **Routing SWAPs ($n=6$)** | 11 SWAPs | **0 SWAPs** | **0 SWAPs** | **100% elimination** of routing gates |
| **Two-Qubit Gates ($n=6$)** | 81 CX (48L + 33S) | 48 CX (48L + 0S) | **48 CX** (48L + 0S) | **$-40.7\%$ fewer 2Q gates** ($33$ CXs saved) |
| **Two-Qubit Depth ($n=4$)** | 37 layers | 16 layers | **16 layers** | **$-56.8\%$ shallower 2Q depth** |
| **Total Circuit Depth ($n=4$)** | 66 layers | 45 layers | **45 layers** | **$-31.8\%$ circuit depth reduction** |
| **Peak Noisy Fidelity $F(n^*)$** | $0.577419$ ($n^*=6$) | $> 0.650000$ | **$0.662816$ ($n^*=7$)** | **$+14.8\%$ relative fidelity improvement** |
| **Optimal Step Count ($n^*$)** | $n^* = 6$ | $n^* \ge 7$ | **$n^* = 7$** | Eliminating SWAPs shifts $n^*$ deeper |

> [!NOTE]
> **Open Architectural Clarification on S3 Metrics:** On the ring TFIM S3 metrics, Processor D ties Processor B by construction ($n^* = 7$, $F = 0.662816$, 56 CXs). This is mathematically expected: any processor providing a native, isolated 4-cycle compiles the ring Hamiltonian with 0 routing SWAPs and identical unitary depth. Processor D's decisive advantages over Processor B reside in its hardware graph metrics: maximum vertex degree $\le 3$ (vs B's 5), mitigating parasitic $ZZ$ crosstalk, and a shorter graph diameter of 3 (vs B's 4), saving 3 CX gates on long-range entangling operations.

---

## 3. Objective 1 — Bell-State Verification on Processor D

### 3.1 Experimental Protocol

Under Objective 1, we prepare the maximally entangled Bell state $|\Phi^+\rangle = \frac{|00\rangle + |11\rangle}{\sqrt{2}}$ on Processor D using a Hadamard gate followed by a CX gate ($H_0 \cdot \text{CX}_{0,1}$). We benchmark:
- **Near Pair:** Physical qubits $(0, 1)$ — distance $d = 1$ hop (directly connected).
- **Far Pair:** Physical qubits $(0, 5)$ — distance $d = 3$ hops (opposite corners of the $2 \times 3$ ladder).

For each pair, we execute 12 independent quantum circuits (2 pairs $\times$ 2 conditions $\times$ 3 measurement bases: $XX, YY, ZZ$) with $N = 4,096$ shots each. The measurement bases are implemented via pre-measurement rotations:
- $ZZ$ Basis: Direct computational measurement.
- $XX$ Basis: Hadamard $H$ on both qubits prior to measurement.
- $YY$ Basis: $S^\dagger$ followed by $H$ on both qubits prior to measurement.

The two-qubit parity expectation value $\langle P \rangle$ and its standard error are computed as:
$$\langle P \rangle = \frac{N_{\text{even}} - N_{\text{odd}}}{N}, \quad \text{SE}(\langle P \rangle) = \sqrt{\frac{1 - \langle P \rangle^2}{N}}$$
The state fidelity is extracted using the exact analytical identity:
$$F = \frac{1 + \langle XX \rangle - \langle YY \rangle + \langle ZZ \rangle}{4}, \quad \text{SE}(F) = \frac{1}{4} \sqrt{\text{SE}_{XX}^2 + \text{SE}_{YY}^2 + \text{SE}_{ZZ}^2}$$

### 3.2 Quantitative Results

![Figure 2: Bell-State Verification](results/figures/fig2_bell_state_fidelity.png)

The full experimental dataset stored in `results/tables/bell.csv` is summarized below:

| Physical Pair | Graph Distance | Condition | 2Q Gates (CX) | $\langle XX \rangle \pm \text{SE}$ | $\langle YY \rangle \pm \text{SE}$ | $\langle ZZ \rangle \pm \text{SE}$ | Fidelity $F \pm \text{SE}$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **$(0, 1)$** | 1 hop | **Ideal** | 1 | $+1.0000 \pm 0.0000$ | $-1.0000 \pm 0.0000$ | $+1.0000 \pm 0.0000$ | **$1.0000 \pm 0.0000$** |
| **$(0, 1)$** | 1 hop | **Noisy** | 1 | $+0.9458 \pm 0.0051$ | $-0.9521 \pm 0.0048$ | $+0.9551 \pm 0.0046$ | **$0.9633 \pm 0.0021$** |
| **$(0, 5)$** | 3 hops | **Ideal** | 7 | $+1.0000 \pm 0.0000$ | $-1.0000 \pm 0.0000$ | $+1.0000 \pm 0.0000$ | **$1.0000 \pm 0.0000$** |
| **$(0, 5)$** | 3 hops | **Noisy** | 7 | $+0.9302 \pm 0.0057$ | $-0.9199 \pm 0.0061$ | $+0.9385 \pm 0.0054$ | **$0.9471 \pm 0.0025$** |

### 3.3 Physical Analysis of Fidelity Degradation

1. **Why Noisy Fidelity Drops on Connected Pair $(0, 1)$ (Readout Dominance):**
   Under ideal execution, parity correlations are exactly $\langle XX \rangle = +1$, $\langle YY \rangle = -1$, $\langle ZZ \rangle = +1$, giving $F = 1.0000$. Under the official noise model, state fidelity drops to $F = 0.9633 \pm 0.0021$.
   A rigorous physical error breakdown reveals that **readout attenuation, not gate noise, is the dominant loss mechanism**:
   - Each single-qubit measurement has bit-flip probability $p_{\text{ro}} = 0.01$. The two-qubit parity expectation value attenuates by:
     $$\langle P \rangle_{\text{noisy}} = (1 - 2 p_{\text{ro}})^2 \langle P \rangle_{\text{ideal}} = (0.98)^2 = 0.9604$$
   - Under readout noise alone, the extracted fidelity would be:
     $$F_{\text{readout}} = \frac{1 + 3 \times 0.9604}{4} = 0.9703$$
   - The measured fidelity is $0.9633 \pm 0.0021$. Thus, readout error accounts for **$73.8\%$** of the fidelity loss ($1.0000 \to 0.9703$), while the single entangling CX gate ($p_2 = 0.005$) and single-qubit gate ($p_1 = 0.001$) account for the remaining $\sim 0.0070$. Readout error is the primary fidelity bottleneck on directly connected transmon qubits.
2. **Why Distant Pair $(0, 5)$ Suffers Further Degradation (Routing Overhead):**
   Physical qubits 0 and 5 are separated by $d=3$ edges across the ladder diagonal ($0 \leftrightarrow 1 \leftrightarrow 2 \leftrightarrow 5$ or $0 \leftrightarrow 3 \leftrightarrow 4 \leftrightarrow 5$). Entangling non-adjacent qubits requires routing through intermediate nodes via SWAP gates.
   Following the optimal routing pattern:
   $$N_{\text{CX}}(d) = 1 + 3(d - 1) \implies N_{\text{CX}}(3) = 1 + 3(2) = 7 \text{ CX gates}$$
   Each additional CX gate applies an independent two-qubit depolarizing channel $\mathcal{E}(\rho) = (1 - p_2)\rho + \frac{p_2}{4} I$. The 7-gate sequence compounds depolarizing errors, lowering parity correlations ($\langle XX \rangle = 0.9302$, $\langle YY \rangle = -0.9199$) and dropping fidelity to $F = 0.9471 \pm 0.0025$. This verifies that distance-dependent routing adds substantial gate noise, but Processor D's diameter of 3 saves 3 CX gates compared to 9-qubit grid topologies ($10$ CXs, $F = 0.9420$).

---

## 4. Objective 3 — PS-S3 Headline Benchmark: Optimal Step Count $n^*$

### 4.1 Step Scan Methodology

We benchmarked the 4-spin ring TFIM simulation across the locked step grid:
$$n \in \{1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16\}$$
For each processor (A, B, and D) and step count $n$:
1. We computed the exact density-matrix fidelity $F(n) = \langle \psi_{\text{exact}}| \rho_n | \psi_{\text{exact}} \rangle$ using `noisy_state_fidelity()`, taking into account transpiler initial layout and final routing permutations.
2. We measured the staggered magnetization $M_s(n)$ via $4,096$ noisy measurement shots using `run_z_measurement()`.
3. We mapped out circuit depth, two-qubit depth, logical CX gates, and routing SWAP counts using `circuit_metrics()`.

### 4.2 Crossover Dynamics & Headline Findings

![Figure 3: Trotter Fidelity Crossover](results/figures/fig3_trotter_fidelity_crossover.png)

The complete step scan data (`results/tables/s3_scan.csv` and `results/tables/s3_summary.csv`) yields the definitive performance comparison:

| Processor | Noise Scale | Optimal $n^*$ | Peak State Fidelity $F(n^*)$ | $M_s(n^*) \pm \text{SE}$ | Total 2Q Gates at $n^*$ | Routing SWAP Gates at $n^*$ |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Processor A** | $1.0\times$ | **6** | **0.577419** | $0.3698 \pm 0.0104$ | 81 | 11 |
| **Processor B** | $1.0\times$ | **7** | **0.662816** | $0.4409 \pm 0.0105$ | 56 | **0** |
| **Processor D (Ours)** | $1.0\times$ | **7** | **0.662816** | $0.4288 \pm 0.0106$ | 56 | **0** |

#### Why $F(n)$ Rises Then Falls (The Crossover Explained)
- **Ascending Regime ($n = 1 \to 6$):** At low step counts, algorithmic Trotter error dominates. For $n=1$, the time step $\Delta t = 2.0$ generates massive coherent phase errors; pure Trotter fidelity is only $0.1150$. As $n$ increases to 6, $\Delta t$ shrinks to $0.333$, suppressing algorithmic discretization error by over $85\%$. The reduction in Trotter error far outweighs the gate noise added by the modest number of gates, driving $F(n)$ rapidly upward.
- **Peak Regime ($n = n^*$):** At the optimal balance point ($n^* = 6$ on A; $n^* = 7$ on B and D), the marginal gain in algorithmic accuracy from an additional step exactly balances the marginal loss of purity from the 8 additional CX gates.
- **Descending Regime ($n > n^*$):** For large $n$, pure Trotter fidelity asymptotically approaches $1.0$ (reaching $0.9856$ at $n=16$). However, physical gate noise dominates. Total two-qubit gates scale linearly ($N_{\text{2Q}} = 8n$), and physical depolarizing noise exponentially degrades state purity:
  $$F(n) \approx F_{\text{Trotter}}(n) \cdot (1 - p_2)^{8n}$$
  At $n=16$, the circuit contains 128 CX gates, driving noisy fidelity down to $0.4819$ on Processor D and $0.3452$ on Processor A.

#### Why Processor A Peaks at $n^* = 6$ with Lower Fidelity
Processor A cannot natively embed a 4-cycle. Routing the ring Hamiltonian on A requires continuous SWAP insertions. At $n=6$, Processor A executes **81 CX gates** (48 logical + 33 routing, from 11 SWAPs). The $68.8\%$ excess in two-qubit gates injects severe depolarizing noise at each step, causing the noise penalty to overwhelm Trotter improvements earlier. Consequently:
- $n^*$ shifts leftward to $n^* = 6$.
- Peak fidelity is capped at $0.5774$, an **$8.54$ percentage point deficit** compared to Processor D ($F = 0.6628$).

---

## 5. Objective 3 — Noise Sweep Analysis on Processor D

### 5.1 Noise Sweep Behavior

![Figure 4: Noise Sweep Dynamics](results/figures/fig4_noise_sweep_shift.png)

We investigated the sensitivity of the optimal step count $n^*$ to physical device quality by sweeping the noise scaling parameter $\lambda \in \{0.5, 1.0, 2.0\}$ on Processor D. The noise model scales all physical error rates proportionally ($p_1 = 0.001\lambda, p_2 = 0.005\lambda, p_{\text{ro}} = 0.01\lambda$):

| Noise Scale $\lambda$ | $p_1$ | $p_2$ | $p_{\text{ro}}$ | Optimal $n^*$ | Peak Fidelity $F(n^*)$ | Noisy $M_s(n^*)$ |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **$0.5\times$ (Low Noise)** | $0.0005$ | $0.0025$ | $0.005$ | **8** | **0.782118** | $0.5157 \pm 0.0105$ |
| **$1.0\times$ (Baseline)** | $0.0010$ | $0.0050$ | $0.010$ | **7** | **0.662816** | $0.4288 \pm 0.0106$ |
| **$2.0\times$ (High Noise)** | $0.0020$ | $0.0100$ | $0.020$ | **6** | **0.503076** | $0.3149 \pm 0.0103$ |

### 5.2 Physical Interpretation of $n^*(\lambda)$ Movement

As noise strength $\lambda$ increases, the optimal step count shifts strictly to smaller values:
$$n^*(\lambda = 0.5) = 8 \;\longrightarrow\; n^*(\lambda = 1.0) = 7 \;\longrightarrow\; n^*(\lambda = 2.0) = 6$$

This shift confirms our mathematical model from §1.3:
1. **Low Noise Regime ($\lambda = 0.5$):** Because gate error rates are halved ($p_2 = 0.0025$), the quantum state retains its purity through deeper circuits. The simulator can execute more Trotter steps ($n^* = 8$, requiring 64 CX gates) to suppress algorithmic Trotter error before noise becomes prohibitive, attaining a peak fidelity of $F = 0.7821$.
2. **High Noise Regime ($\lambda = 2.0$):** At doubled error rates ($p_2 = 0.010$), each Trotter step incurs double the depolarizing penalty. The noise threshold is reached earlier, truncating allowable circuit depth to $n^* = 6$ (48 CX gates). Attempting to run deeper circuits at $\lambda = 2.0$ results in severe fidelity decay ($F(8) = 0.4586$).

---

## 6. Objective 3 — Gate Split & Topological Routing Analysis

### 6.1 Gate Decomposition Across Topologies

![Figure 5: Gate Split & Depth Scaling](results/figures/fig5_gate_split_and_depth.png)

The complete resource metrics compiled in `results/tables/resources.csv` reveal the structural overhead imposed by processor coupling graphs:

| Step Count $n$ | Logical 2Q Gates | Processor D 2Q Gates | Processor D SWAPs | Processor A 2Q Gates | Processor A SWAPs | Routing Overhead on A ($\Delta \text{CX}$) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 1 | 8 | 8 | 0 | 11 | 1 | +3 CX (+37.5%) |
| 2 | 16 | 16 | 0 | 25 | 3 | +9 CX (+56.3%) |
| 4 | 32 | 32 | 0 | 53 | 7 | +21 CX (+65.6%) |
| 6 | 48 | 48 | 0 | 81 | 11 | +33 CX (+68.8%) |
| 8 | 64 | 64 | 0 | 109 | 15 | +45 CX (+70.3%) |
| 16 | 128 | 128 | 0 | 221 | 31 | +93 CX (+72.7%) |

### 6.2 Topological Routing Breakdown

1. **Native 4-Cycle Advantage (Processor D & B):**
   Because the ring TFIM requires closed periodic interactions $(0,1), (1,2), (2,3), (3,0)$, processors embedding a native 4-cycle can place the 4 logical qubits directly onto physical hardware nodes with matching nearest-neighbor connectivity. Consequently, routing SWAPs are identically **zero** across all step counts.
2. **Linear Tree Penalty (Processor A):**
   Processor A's coupling map is a tree graph. Embedding the ring forces the transpiler to route non-adjacent interactions via SWAPs. At step count $n=4$, Processor A has total depth 66 (vs 45 on D) and 2Q depth 37 (vs 16 on D). At $n=6$, Processor A incurs 11 SWAPs (33 extra CXs), yielding total depth 100 (vs 67 on D) and 2Q depth 57 (vs 24 on D), accelerating decoherence and entropy production.

---

## 7. Bonus — Quantum Error Mitigation on Processor D

![Figure 6: Magnetization & Mitigation](results/figures/fig6_magnetization_and_mitigation.png)

### 7.1 Mitigation Methodology at $n^* = 7$

At the optimal step count $n^* = 7$ on Processor D, the exact analytical staggered magnetization is $M_s^{\text{exact}} = 0.652444$. However, raw 4,096-shot measurement under the $1.0\times$ noise model yields:
$$M_s^{\text{raw}} = 0.428833 \pm 0.010605 \quad (\text{Absolute Error } |\Delta M_s| = 0.223611)$$

We implemented two rigorous, mathematically sound error mitigation techniques:

1. **Readout Error Inversion (Single-Qubit Observable Inversion):**
   The staggered magnetization $M_s = \frac{1}{4} \sum_{i=0}^{3} (-1)^i \langle Z_i \rangle$ is a linear combination of single-qubit expectation values $\langle Z_i \rangle$. Under symmetric single-qubit readout bit-flip error $p_{\text{ro}} = 0.01$, each individual expectation value attenuates as:
   $$\langle Z_i \rangle_{\text{noisy}} = (1 - 2 p_{\text{ro}}) \langle Z_i \rangle_{\text{ideal}} = 0.98 \, \langle Z_i \rangle_{\text{ideal}}$$
   Therefore, the unbiased expectation value is recovered by dividing by the single-qubit attenuation factor $0.98$:
   $$M_s^{\text{RO}} = \frac{M_s^{\text{raw}}}{1 - 2 p_{\text{ro}}} = \frac{0.428833}{0.98} = 0.437585 \pm 0.010822$$
   *(Note: Applying $(0.98)^4$ would be mathematically invalid because $M_s$ is a sum of single-qubit observables, not a 4-qubit parity product. The $0.98$ factor correctly inverts the physical single-qubit readout channel).*
   *Calibration Cost:* 8,192 calibration shots.

2. **Zero-Noise Extrapolation (ZNE via Richardson Extrapolation):**
   We utilize empirical measurements at scaled noise factors $\lambda \in \{0.5, 1.0, 2.0\}$ on Processor D. Because the noise scaling scales all physical noise channels simultaneously (gate depolarizing and readout errors), ZNE simultaneously extrapolates away both gate and measurement errors.
   Using second-order Richardson extrapolation (Lagrange polynomial evaluated at $\lambda = 0$):
   $$M_s^{\text{ZNE}} = \frac{8}{3} M_s(0.5) - 2 M_s(1.0) + \frac{1}{3} M_s(2.0) = 0.620483 \pm 0.035246$$
   *Calibration Cost:* 8,192 extra shots ($4,096$ shots at $\lambda = 0.5$ and $\lambda = 2.0$).
   *(Note: Because ZNE already scales the readout error rate alongside gate noise, applying readout inversion on top of ZNE would constitute double-counting. We report ZNE independently to ensure rigorous, non-inflated scientific conclusions).*

### 7.2 Mitigation Performance Summary

The mitigation dataset saved in `results/tables/mitigation.csv` is summarized below:

| Method | Extra Shots | Raw $M_s \pm \text{SE}$ | Mitigated $M_s \pm \text{SE}$ | Target $M_s^{\text{exact}}$ | Absolute Error $|\Delta M_s|$ | Error Reduction (%) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| **Raw Baseline (Noisy $1.0\times$)** | 0 | $0.428833 \pm 0.010605$ | — | $0.652444$ | $0.223611$ | $0.0\%$ (Reference) |
| **Readout Error Inversion** | 8,192 | $0.428833 \pm 0.010605$ | $0.437585 \pm 0.010822$ | $0.652444$ | $0.214859$ | **3.9%** |
| **Richardson ZNE** | 8,192 | $0.428833 \pm 0.010605$ | **$0.620483 \pm 0.035246$** | $0.652444$ | **0.031961** | **85.7%** |

Zero-Noise Extrapolation reduces observable bias by **$85.7\%$** ($0.2236 \to 0.0320$), bringing the measured staggered magnetization within $1 \sigma$ of the exact analytical value.

---

## 8. Limitations & Physical Hardware Roadmap

### 8.1 Methodological Limitations
1. **Uncorrelated Depolarizing Approximation:** The Aer noise model implements independent single- and two-qubit depolarizing channels. Real superconducting quantum processors exhibit non-Markovian noise, coherent phase drift, and state-dependent $T_1 / T_2$ relaxation.
2. **Fixed First-Order Trotter:** While first-order Trotter-Suzuki provides predictable $\mathcal{O}(\Delta t)$ scaling, higher-order decompositions (such as second-order symmetric Suzuki $U_2(\Delta t)$) cancel odd-order commutator errors and may attain higher peak fidelities at lower step counts.

### 8.2 Roadmap to Physical IBM Quantum Processors
To deploy StepWise-Q on physical IBM Quantum QPUs (e.g., Eagle or Heron architectures):
1. **Heavy-Hex Subgraph Embedding:** Physical QPUs utilize heavy-hexagonal coupling lattices with degree $\le 3$. Embedding our 4-spin ring requires identifying a hexagonal 6-cycle and applying a single bridge SWAP to emulate the periodic boundary condition.
2. **Dynamical Decoupling (DD):** During idle periods in Trotter layers, idle qubits suffer dephasing. Applying XY4 or CPMG periodic pulse sequences will suppress environmental phase noise.
3. **Hardware-Native RZZ Compilation:** In superconducting architectures with direct cross-resonance or flux-tunable couplers, implementing native $R_{ZZ}(\theta)$ gates replaces the 2-CX decomposition with a single calibrated pulse, cutting two-qubit gate counts in half and extending $n^*$ deeper into the low-discretization regime.

---

## 9. Conclusion

The StepWise-Q investigation for QFF 2026 Phase 2 establishes clear, quantitatively verified conclusions:
1. **$n^*$ Crossover is Real and Governed by Topology:** The optimal Trotter step count $n^*$ represents an unavoidable physical trade-off between algorithmic error ($\propto 1/n$) and physical gate noise ($\propto n$). Processor A is limited to $n^* = 6$ ($F = 0.5774$) due to tree routing overhead, whereas Processor D reaches $n^* = 7$ ($F = 0.6628$) via native 4-cycle co-design.
2. **Noise Strength Shifts $n^*$ Predictably:** Doubling physical noise shifts $n^*$ from 7 to 6, while halving noise allows scaling to $n^* = 8$ ($F = 0.7821$).
3. **Hardware Co-Design Outperforms Pure Compilation:** Providing native cycles matching problem Hamiltonian symmetries eliminates routing SWAPs entirely, providing greater fidelity gains than transpiler optimization alone.
4. **Error Mitigation is Highly Effective:** Richardson Zero-Noise Extrapolation reduces observable bias by **$85.7\%$**, recovering the physical signal lost to device attenuation.
