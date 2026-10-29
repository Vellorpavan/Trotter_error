# Qiskit Fall Fest 2026 — Official Submission Form
### Team StepWise-Q | Problem Statement PS-S3: Noise Accumulation in Trotter Circuits

Use the exact text and fields below when completing the submission form on the hackathon portal:
👉 **Submission Portal:** [https://qffsrmap2026.com/learning/hackathon](https://qffsrmap2026.com/learning/hackathon)

---

## 📋 Portal Form Fields (Copy & Paste Ready)

### 1. Project Title
```
StepWise-Q: Optimal Trotter Step Prediction & Hardware-Aware Co-Design for Ring TFIM Quantum Simulation
```

### 2. Team Information
- **Team Name:** `StepWise-Q`
- **Institution:** `SRM University-AP`
- **Track:** `Track 3: Quantum Simulation`
- **Problem Statement:** `PS-S3: Noise Accumulation in Trotter Circuits`
- **Team Members:**
  1. **Jagadesh** (`Lead`) — GitHub: [`@jagadesh-jpg`](https://github.com/jagadesh-jpg) | Email: `jagadeshai.svpcet@gmail.com`
  2. **Vellor Pavan** (`Lead Developer`) — GitHub: [`@Vellorpavan`](https://github.com/Vellorpavan) | Email: `91930057+Vellorpavan@users.noreply.github.com`

### 3. Repository & Live Demo URLs
- **Primary GitHub Repository:** `https://github.com/jagadesh-jpg/qiskit-noise-lab`
- **Secondary GitHub Repository:** `https://github.com/Vellorpavan/Trotter_error`
- **Live 3D WebGL Simulation Lab:** `https://vellorpavan.github.io/Trotter_error/`
- **Interactive Cloud View:** `https://gemini.google.com/share/741c9c5cceca?social_preview=false&skid=6c798fbe-3abe-466e-aaea-d979cbf5dad4`

---

## 📝 Project Description / Abstract (Form Text Box)

```
In digital quantum simulation (DQS) of continuous Hamiltonian dynamics, time evolution e^(-iHT) is approximated by discrete Trotter-Suzuki steps. This creates a fundamental physical crossover on NISQ hardware: increasing the step count n suppresses coherent algorithmic discretization error as O(T^2/n), but each step inserts non-commuting entangling gates whose physical depolarizing noise (p2 = 0.005) compounds exponentially as (1 - p2)^(8n). The physical state fidelity F(n) therefore reaches a strict convex peak at an optimal step count n*.

Focusing on Problem Statement PS-S3 (4-spin periodic Ring Transverse-Field Ising Model, J=1.0, g=1.0, T=2.0, initial Neel state |1010>), StepWise-Q demonstrates that processor coupling topology directly dictates simulation accuracy. On Processor A (5-qubit tree), the absence of cycles forces dynamic SWAP routing (+33 CX gates at n=6, an overhead of +68.8%), prematurely capping peak fidelity at F = 0.5774 (n* = 6). 

To overcome this, we co-designed Processor D: a 6-qubit 2x3 ladder architecture (edges: (0,1), (1,2), (3,4), (4,5), (0,3), (1,4), (2,5)) featuring native RZZ basis gates (["rz", "sx", "x", "rzz"]) satisfying all competition design constraints (node count <= 12, max degree <= 3, diameter 3). Processor D natively embeds two 4-cycles, eliminating 100% of ring routing SWAPs, while its native RZZ pulses compile each ZZ interaction directly, cutting two-qubit entangling gates per step from 8 CX to 4 RZZ. This cuts entangling gate noise in half, unlocking deeper optimal depth (n* = 8) and boosting peak state fidelity to F = 0.7528 (beating Processor B by +9.0 percentage points and Processor A by +17.5 percentage points). Bounding the degree to <= 3 completely avoids the high-degree (degree 5) frequency-crowding crosstalk hub found in Processor B, while the diameter of 3 allows far-pair Bell state generation in only 7 CX gates (vs 10 CX on 9-qubit grids).

Finally, we address observable fidelity degradation at n* = 8. We prove that staggered magnetization Ms is a linear combination of single-qubit Z observables, meaning readout attenuation scales as (1 - 2*p_ro) = 0.98 (rather than the invalid (0.98)^4 product). Applying 2nd-order Richardson Zero-Noise Extrapolation (ZNE) across noise scales lambda in {0.5, 1.0, 2.0} recovers 80.9% of the lost physical signal toward the exact reference (Ms = 0.652444), bringing the mitigated value (Ms = 0.6224 +/- 0.0352) within 1 sigma. All benchmarks execute end-to-end in under 5 seconds (well within the 30-minute limit) and are verified by a 14-test suite and a live 3D WebGL laboratory.
```

---

## 📊 Summary of Headline Benchmark Results

| Metric | Reference Tree (Processor A) | Heavy-Hex (Processor B) | Co-Designed Ladder (Processor D) | Impact / Physical Meaning |
|:---|:---:|:---:|:---:|:---|
| **Physical Qubits** | 5 qubits | 9 qubits | **6 qubits** | Compact planar footprint |
| **Max Vertex Degree** | 3 (node 1) | 5 (node 6) | **3 (nodes 1, 4)** | Eliminates parasitic $ZZ$ crosstalk |
| **Graph Diameter** | 3 hops | 4 hops | **3 hops** | 7 CX far-pair Bell generation |
| **Optimal Step Count ($n^*$)** | **$n^* = 6$** | **$n^* = 7$** | **$n^* = 8$** | Co-design enables deeper Trotter circuits |
| **Peak State Fidelity $F(n^*)$** | $0.577419$ | $0.662816$ | **$0.752817$** | **$+13.6\%$ relative gain** vs B (**$+30.4\%$** vs A) |
| **Two-Qubit Gates at $n=6$** | 81 CX (48L + 33S) | 48 CX (48L + 0S) | **24 RZZ** (24L + 0S) | **$-70.4\%$ fewer entangling gates** vs A ($-50\%$ vs B) |
| **Routing SWAP Overhead** | 11 SWAPs | **0 SWAPs** | **0 SWAPs** | 100% elimination of routing gates |
| **Circuit Depth at $n=4$** | 66 layers | 45 layers | **29 layers** | **$-56.1\%$ circuit depth reduction** vs A |
| **Bell Pair Fidelity (Near: 1 hop)** | — | — | **$0.9598 \pm 0.0022$** (1 CX) | Readout is dominant loss ($F_{\text{ro}} \approx 0.970$) |
| **Bell Pair Fidelity (Far: 3 hops)** | — | — | **$0.9384 \pm 0.0027$** (7 CX) | Matches $1 + 3(d-1)$ compilation |
| **Raw Magnetization $M_s(n^*)$** | $0.3698 \pm 0.0104$ | $0.4409 \pm 0.0105$ | **$0.4949 \pm 0.0105$** | Target $M_s^{\text{exact}} = 0.652444$ |
| **ZNE Mitigated $M_s$** | — | — | **$0.6224 \pm 0.0352$** | **$80.9\%$ observable error recovery** |

---

## 🛠️ Verification & Reproduction Instructions for Judges

To replicate our entire dataset from scratch on any machine:

```bash
# 1. Clone repository
git clone https://github.com/jagadesh-jpg/qiskit-noise-lab.git
cd qiskit-noise-lab

# 2. Set up virtual environment and install dependencies
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 3. Validate Processor D architectural rules
python tools/check_processor_d.py processors/processor_D.json 4
# Output: OK - Processor D satisfies all design rules.

# 4. Run automated 14-test verification suite (< 0.02s)
python -m unittest tests/test_all.py
# Output: Ran 14 tests in 0.009s - OK

# 5. Headless end-to-end reproduction (~4s)
python run_all.py

# 6. Execute master notebook (~4.5s)
jupyter nbconvert --to notebook --execute --inplace main.ipynb
```
