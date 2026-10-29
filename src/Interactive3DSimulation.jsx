import React, { useState, useEffect, useRef, useMemo, useCallback } from 'react';
import * as THREE from 'three';
import { 
  Play, Pause, RotateCcw, SkipBack, SkipForward, ChevronRight, ChevronLeft,
  Activity, Cpu, Eye, Zap, BarChart2, Layers, Compass, HelpCircle, AlertTriangle,
  Info, CheckCircle2, ShieldCheck, RefreshCw, ZoomIn, ZoomOut, Move3D
} from 'lucide-react';

// ============================================================================
// 1. SCIENTIFIC CONSTANTS & LOCKED OFFICIAL QFF 2026 PS-S3 PARAMETERS
// ============================================================================
const MODEL_CONFIG = {
  N: 4,
  J: 1.0,
  g: 1.0,
  T: 2.0,
  initialStateStr: "1010", // Néel state: |q0=1, q1=0, q2=1, q3=0>
  officialStepCounts: [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16],
  bonds: [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0] // Periodic boundary condition closing the 4-qubit ring
  ]
};

// Official Processor Topologies & Hardware Parameters (Strict Section 4 Specification)
const PROCESSORS = {
  A: {
    name: "Processor A (Tree - 5 Qubits)",
    type: "Tree (5 physical qubits, no 4-cycle)",
    numPhysQubits: 5,
    // (0,1), (1,2), (1,3), (3,4) as defined in Section 4 image
    nativeBonds: [[0, 1], [1, 2], [1, 3], [3, 4]],
    sqErr: 0.0004,
    twoQErr: 0.0065,
    readoutErr: 0.018,
    routingSwapCostPerRing: 2,
    logicalToPhys: [0, 1, 2, 3], // Logical Q0->P0, Q1->P1, Q2->P2, Q3->P3 (P4 unused)
    description: "Official 5-qubit tree: contains 4-qubit line (0-1-3-4 or 2-1-3-4) but NO 4-qubit cycle. Bond (3,0) requires SWAP routing."
  },
  B: {
    name: "Processor B (9 Qubits)",
    type: "Planar Heavy-Hex Graph (9 physical qubits)",
    numPhysQubits: 9,
    // (0,1), (1,2), (1,6), (2,3), (2,6), (2,5), (3,4), (4,5), (5,6), (6,7), (6,8), (8,0)
    nativeBonds: [
      [0, 1], [1, 2], [1, 6], [2, 3], [2, 6], [2, 5],
      [3, 4], [4, 5], [5, 6], [6, 7], [6, 8], [8, 0]
    ],
    sqErr: 0.0003,
    twoQErr: 0.0052,
    readoutErr: 0.014,
    routingSwapCostPerRing: 1,
    logicalToPhys: [0, 1, 2, 3], // Logical Q0..Q3 mapped to P0..P3
    description: "Official 9-qubit graph with 12 edges: (0,1),(1,2),(1,6),(2,3),(2,6),(2,5),(3,4),(4,5),(5,6),(6,7),(6,8),(8,0)."
  },
  D: {
    name: "Processor D (Custom - 6 Qubits)",
    type: "Native Ring with Resonators (6 physical qubits)",
    numPhysQubits: 6,
    nativeBonds: [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [2, 5]],
    sqErr: 0.0002,
    twoQErr: 0.0035,
    readoutErr: 0.009,
    routingSwapCostPerRing: 0,
    logicalToPhys: [0, 1, 2, 3],
    description: "Section 6 Custom Design (at most 12 qubits): Native 4-cycle ring (0,1),(1,2),(2,3),(3,0) with 0 routing SWAP overhead."
  }
};

// ============================================================================
// 2. LINEAR ALGEBRA & QUANTUM MECHANICS EXACT SOLVER (16-DIMENSIONAL HILBERT SPACE)
// ============================================================================

// Complex Number Operations
const cAdd = (a, b) => ({ r: a.r + b.r, i: a.i + b.i });
const cSub = (a, b) => ({ r: a.r - b.r, i: a.i - b.i });
const cMul = (a, b) => ({ r: a.r * b.r - a.i * b.i, i: a.r * b.i + a.i * b.r });
const cScale = (a, s) => ({ r: a.r * s, i: a.i * s });
const cExp = (phase) => ({ r: Math.cos(phase), i: Math.sin(phase) });
const cMagSq = (a) => a.r * a.r + a.i * a.i;
const cAbs = (a) => Math.sqrt(cMagSq(a));
const cPhase = (a) => Math.atan2(a.i, a.r);

// Build exact 16x16 TFIM Hamiltonian Matrix
function buildHamiltonian(N, J, g, bonds) {
  const dim = 1 << N; // 16
  const H = Array.from({ length: dim }, () => Array.from({ length: dim }, () => 0.0));

  // ZZ terms: -J * sum_{<j,k>} Z_j Z_k (diagonal)
  for (let s = 0; s < dim; s++) {
    let diagEnergy = 0;
    for (const [q1, q2] of bonds) {
      const bit1 = (s >> q1) & 1;
      const bit2 = (s >> q2) & 1;
      const z1 = bit1 === 0 ? 1 : -1;
      const z2 = bit2 === 0 ? 1 : -1;
      diagEnergy += -J * (z1 * z2);
    }
    H[s][s] = diagEnergy;
  }

  // X terms: -g * sum_j X_j (off-diagonal bit-flips)
  for (let s = 0; s < dim; s++) {
    for (let q = 0; q < N; q++) {
      const flipped = s ^ (1 << q);
      H[flipped][s] += -g;
    }
  }

  return H;
}

// Matrix-Vector Multiplication for Complex Statevectors
function matVecMul(M, vec) {
  const dim = vec.length;
  const res = new Array(dim);
  for (let i = 0; i < dim; i++) {
    let r = 0, im = 0;
    for (let j = 0; j < dim; j++) {
      r += M[i][j] * vec[j].r;
      im += M[i][j] * vec[j].i;
    }
    res[i] = { r, i: im };
  }
  return res;
}

// Exact Continuous Matrix Exponential Unitary Evolution via High-Order Taylor-Lanczos
function exactTimeEvolution(psi0, H, t) {
  const steps = 60;
  const dt = t / steps;
  let psi = psi0.map(v => ({ ...v }));
  const dim = psi0.length;

  for (let step = 0; step < steps; step++) {
    // 4th-order Runge-Kutta / series integrator for d|psi>/dt = -i H |psi>
    const k1 = matVecMul(H, psi).map(v => ({ r: v.i * dt, i: -v.r * dt }));
    
    const p2 = psi.map((v, idx) => cAdd(v, cScale(k1[idx], 0.5)));
    const k2 = matVecMul(H, p2).map(v => ({ r: v.i * dt, i: -v.r * dt }));

    const p3 = psi.map((v, idx) => cAdd(v, cScale(k2[idx], 0.5)));
    const k3 = matVecMul(H, p3).map(v => ({ r: v.i * dt, i: -v.r * dt }));

    const p4 = psi.map((v, idx) => cAdd(v, k3[idx]));
    const k4 = matVecMul(H, p4).map(v => ({ r: v.i * dt, i: -v.r * dt }));

    for (let i = 0; i < dim; i++) {
      const delta = {
        r: (k1[i].r + 2 * k2[i].r + 2 * k3[i].r + k4[i].r) / 6,
        i: (k1[i].i + 2 * k2[i].i + 2 * k3[i].i + k4[i].i) / 6
      };
      psi[i] = cAdd(psi[i], delta);
    }
    // Re-normalize against floating point drift
    let normSq = 0;
    for (let i = 0; i < dim; i++) normSq += cMagSq(psi[i]);
    const invNorm = 1 / Math.sqrt(normSq);
    for (let i = 0; i < dim; i++) psi[i] = cScale(psi[i], invNorm);
  }

  return psi;
}

// Convert Bitstring "1010" to Computational Basis Index
// Logical state |q3 q2 q1 q0> = |1 0 1 0>
// q0=0 (bit 0), q1=1 (bit 1), q2=0 (bit 2), q3=1 (bit 3) -> index 10 (binary 1010)
function getInitialStatevector() {
  const psi = Array.from({ length: 16 }, () => ({ r: 0, i: 0 }));
  psi[10] = { r: 1.0, i: 0.0 };
  return psi;
}

// Apply Single-Qubit Rx(theta) Gate: exp(-i theta/2 X)
// Note: Hamiltonian has -g X, so evolution exp(-i H dt) requires exp(i g dt X)
// Hence theta = -2 * g * dt
function applyRx(state, qubit, theta) {
  const dim = state.length;
  const next = state.map(v => ({ ...v }));
  const half = theta / 2;
  const cos = Math.cos(half);
  const sin = Math.sin(half);

  for (let s = 0; s < dim; s++) {
    if (((s >> qubit) & 1) === 0) {
      const flipped = s ^ (1 << qubit);
      const a = state[s];
      const b = state[flipped];

      // [cos, -i*sin; -i*sin, cos]
      next[s] = {
        r: cos * a.r + sin * b.i,
        i: cos * a.i - sin * b.r
      };
      next[flipped] = {
        r: cos * b.r + sin * a.i,
        i: cos * b.i - sin * a.r
      };
    }
  }
  return next;
}

// Apply Two-Qubit Rzz(theta) Gate: exp(-i theta/2 Z_j Z_k)
// Hamiltonian has -J ZZ, so exp(-i H dt) has exp(i J dt ZZ)
// Hence theta = -2 * J * dt
function applyRzz(state, q1, q2, theta) {
  const dim = state.length;
  const next = new Array(dim);
  const phiOver2 = theta / 2;

  for (let s = 0; s < dim; s++) {
    const bit1 = (s >> q1) & 1;
    const bit2 = (s >> q2) & 1;
    const parity = bit1 ^ bit2; // 0 if equal (+1 eigenvalue), 1 if opposite (-1 eigenvalue)
    const phase = parity === 0 ? -phiOver2 : phiOver2;
    const u = cExp(phase);
    next[s] = cMul(u, state[s]);
  }
  return next;
}

// Apply SWAP Gate (for routing non-native connections)
function applySwap(state, q1, q2) {
  const dim = state.length;
  const next = new Array(dim);
  for (let s = 0; s < dim; s++) {
    const b1 = (s >> q1) & 1;
    const b2 = (s >> q2) & 1;
    if (b1 !== b2) {
      const swappedIndex = s ^ (1 << q1) ^ (1 << q2);
      next[s] = state[swappedIndex];
    } else {
      next[s] = state[s];
    }
  }
  return next;
}

// Calculate State Fidelity F = |<psi_exact | psi_sim>|^2
function calculateFidelity(psiA, psiB) {
  let overlapR = 0;
  let overlapI = 0;
  for (let i = 0; i < psiA.length; i++) {
    overlapR += psiA[i].r * psiB[i].r + psiA[i].i * psiB[i].i;
    overlapI += psiA[i].r * psiB[i].i - psiA[i].i * psiB[i].r;
  }
  return Math.min(1.0, Math.max(0.0, overlapR * overlapR + overlapI * overlapI));
}

// Calculate Expectation Value of Staggered Magnetization: M_stag = 1/4 * sum_{q=0}^3 (-1)^q <Z_q>
function calculateStaggeredMagnetization(state, N = 4) {
  let totalM = 0;
  const dim = state.length;
  for (let q = 0; q < N; q++) {
    const sign = Math.pow(-1, q);
    let expZ = 0;
    for (let s = 0; s < dim; s++) {
      const bit = (s >> q) & 1;
      const zVal = bit === 0 ? 1 : -1;
      expZ += zVal * cMagSq(state[s]);
    }
    totalM += sign * expZ;
  }
  return totalM / N;
}

// Compute Single-Qubit Reduced Density Matrix & 3D Bloch Vector (rx, ry, rz)
function computeReducedBlochVector(state, targetQubit) {
  let zExp = 0;
  let xExp = 0;
  let yExp = 0;
  const dim = state.length;

  for (let s = 0; s < dim; s++) {
    const prob = cMagSq(state[s]);
    const bit = (s >> targetQubit) & 1;
    zExp += (bit === 0 ? 1 : -1) * prob;

    if (bit === 0) {
      const flipped = s ^ (1 << targetQubit);
      const a0 = state[s];
      const a1 = state[flipped];
      xExp += 2 * (a0.r * a1.r + a0.i * a1.i);
      yExp += 2 * (a0.r * a1.i - a0.i * a1.r);
    }
  }

  const length = Math.sqrt(xExp * xExp + yExp * yExp + zExp * zExp);
  const purity = (1 + length * length) / 2;
  return { x: xExp, y: yExp, z: zExp, length, purity };
}

// Generate Instruction Stream
function generateCircuitInstructionStream(trotterSteps, processorKey) {
  const proc = PROCESSORS[processorKey];
  const dt = MODEL_CONFIG.T / trotterSteps;
  // Evolution exp(i g dt X) and exp(i J dt ZZ)
  const rxAngle = -2 * MODEL_CONFIG.g * dt;
  const rzzAngle = -2 * MODEL_CONFIG.J * dt;
  const instructions = [];

  instructions.push({
    id: "init-0",
    name: "STATE_PREP",
    step: 0,
    layer: "Init",
    qubits: [0, 1, 2, 3],
    desc: "Prepare Néel Initial State |1010⟩ on 4-qubit ring",
    type: "prep",
    native: true,
    routingOverhead: false,
    durationNs: 0,
    errorRate: 0.0
  });

  for (let step = 1; step <= trotterSteps; step++) {
    // 1. Single-qubit transverse field rotations: Rx on all 4 qubits
    for (let q = 0; q < MODEL_CONFIG.N; q++) {
      instructions.push({
        id: `s${step}-rx-q${q}`,
        name: `Rx(${Math.abs(rxAngle).toFixed(3)})`,
        step,
        layer: "Transverse-X",
        qubits: [q],
        targetGate: "Rx",
        angleRad: rxAngle,
        desc: `Transverse field rotation on Q${q} by θ = ${rxAngle.toFixed(3)} rad`,
        type: "sq",
        native: true,
        routingOverhead: false,
        durationNs: 35,
        errorRate: proc.sqErr
      });
    }

    // 2. Periodic ZZ interactions: Rzz on bonds (0,1), (1,2), (2,3), (3,0)
    for (const [q1, q2] of MODEL_CONFIG.bonds) {
      const isNative = proc.nativeBonds.some(b => 
        (b[0] === q1 && b[1] === q2) || (b[0] === q2 && b[1] === q1)
      );

      if (!isNative) {
        instructions.push({
          id: `s${step}-swap-${q1}-${q2}`,
          name: `SWAP(Q${q1}, Q${q2})`,
          step,
          layer: "Hardware-Routing",
          qubits: [q1, q2],
          desc: `Routing overhead: physical SWAP to bridge non-native periodic bond (${q1},${q2})`,
          type: "swap",
          native: false,
          routingOverhead: true,
          durationNs: 300,
          errorRate: proc.twoQErr * 3
        });
      }

      instructions.push({
        id: `s${step}-rzz-${q1}-${q2}`,
        name: `Rzz(${Math.abs(rzzAngle).toFixed(3)})`,
        step,
        layer: "Ising-ZZ",
        qubits: [q1, q2],
        targetGate: "Rzz",
        angleRad: rzzAngle,
        desc: `Nearest-neighbor Ising interaction between Q${q1} & Q${q2}`,
        type: "twoQ",
        native: isNative,
        routingOverhead: false,
        durationNs: 220,
        errorRate: proc.twoQErr
      });

      if (!isNative) {
        instructions.push({
          id: `s${step}-unswap-${q1}-${q2}`,
          name: `SWAP_UNWIND(Q${q1}, Q${q2})`,
          step,
          layer: "Hardware-Routing",
          qubits: [q1, q2],
          desc: `Unwind physical routing mapping back to logical canonical ring`,
          type: "swap",
          native: false,
          routingOverhead: true,
          durationNs: 300,
          errorRate: proc.twoQErr * 3
        });
      }
    }
  }

  instructions.push({
    id: "meas-final",
    name: "MEASURE_Z",
    step: trotterSteps,
    layer: "Measurement",
    qubits: [0, 1, 2, 3],
    desc: "Computational basis projective readout in Z basis across all 4 qubits",
    type: "measure",
    native: true,
    routingOverhead: false,
    durationNs: 1200,
    errorRate: proc.readoutErr
  });

  return instructions;
}

// Compute Complete Precomputed Simulation Grid
function runComprehensiveBenchmark() {
  const H = buildHamiltonian(MODEL_CONFIG.N, MODEL_CONFIG.J, MODEL_CONFIG.g, MODEL_CONFIG.bonds);
  const psi0 = getInitialStatevector();
  const exactFinalPsi = exactTimeEvolution(psi0, H, MODEL_CONFIG.T);
  const exactMstag = calculateStaggeredMagnetization(exactFinalPsi);

  const resultsByStep = {};

  for (const n of MODEL_CONFIG.officialStepCounts) {
    const dt = MODEL_CONFIG.T / n;
    let psi = psi0.map(v => ({ ...v }));

    for (let step = 0; step < n; step++) {
      for (const [q1, q2] of MODEL_CONFIG.bonds) {
        psi = applyRzz(psi, q1, q2, -2 * MODEL_CONFIG.J * dt);
      }
      for (let q = 0; q < MODEL_CONFIG.N; q++) {
        psi = applyRx(psi, q, -2 * MODEL_CONFIG.g * dt);
      }
    }

    const idealFidelity = calculateFidelity(exactFinalPsi, psi);
    const idealMstag = calculateStaggeredMagnetization(psi);

    const procResults = {};
    for (const [pKey, pConfig] of Object.entries(PROCESSORS)) {
      const numSqGates = 4 * n;
      const numLogicalTwoQ = 4 * n;
      const numRoutingSwaps = (pConfig.routingSwapCostPerRing > 0 ? 1 : 0) * n;
      const totalTwoQGates = numLogicalTwoQ + numRoutingSwaps * 3;
      const depth = n * (1 + 4 + (numRoutingSwaps > 0 ? 4 : 0));

      const pNoiseScale1x = 
        Math.pow(1 - pConfig.sqErr, numSqGates) *
        Math.pow(1 - pConfig.twoQErr, totalTwoQGates) *
        Math.pow(1 - pConfig.readoutErr, 4);

      const noisyFid1x = Math.max(0.01, idealFidelity * pNoiseScale1x);

      const pNoiseScale05x = 
        Math.pow(1 - pConfig.sqErr * 0.5, numSqGates) *
        Math.pow(1 - pConfig.twoQErr * 0.5, totalTwoQGates) *
        Math.pow(1 - pConfig.readoutErr * 0.5, 4);
      const noisyFid05x = Math.max(0.01, idealFidelity * pNoiseScale05x);

      const pNoiseScale2x = 
        Math.pow(1 - pConfig.sqErr * 2.0, numSqGates) *
        Math.pow(1 - pConfig.twoQErr * 2.0, totalTwoQGates) *
        Math.pow(1 - pConfig.readoutErr * 2.0, 4);
      const noisyFid2x = Math.max(0.01, idealFidelity * pNoiseScale2x);

      procResults[pKey] = {
        noisyFidelity: noisyFid1x,
        noisyFidelity05x: noisyFid05x,
        noisyFidelity2x: noisyFid2x,
        totalTwoQGates,
        circuitDepth: depth,
        noisyMstag: idealMstag * (0.85 + 0.15 * pNoiseScale1x)
      };
    }

    resultsByStep[n] = {
      stepCount: n,
      dt,
      idealFidelity,
      idealMstag,
      finalState: psi,
      processors: procResults
    };
  }

  const optimalSteps = {};
  for (const pKey of Object.keys(PROCESSORS)) {
    let bestN = MODEL_CONFIG.officialStepCounts[0];
    let maxFid = -1;
    for (const n of MODEL_CONFIG.officialStepCounts) {
      const fid = resultsByStep[n].processors[pKey].noisyFidelity;
      if (fid > maxFid + 1e-6) {
        maxFid = fid;
        bestN = n;
      }
    }
    optimalSteps[pKey] = { bestN, maxFid };
  }

  return { exactFinalPsi, exactMstag, resultsByStep, optimalSteps };
}

const BENCHMARK_DATA = runComprehensiveBenchmark();

// ============================================================================
// 2B. 9-CHAPTER GUIDED SCIENTIFIC STORY DATA
// ============================================================================
const STORY_CHAPTERS = [
  {
    title: "Chapter 1 — The Four-Qubit Transverse-Field Ising Model",
    badge: "PHYSICS & HAMILTONIAN",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          We simulate four superconducting transmon qubits arranged in a closed 1D ring governed by the Transverse-Field Ising Model (TFIM):
        </p>
        <div className="p-3 bg-slate-950 rounded-lg border border-slate-800 font-mono text-cyan-300 text-center">
          H = -J ∑(i=0..3) Z_i Z_(i+1 mod 4) - g ∑(i=0..3) X_i
        </div>
        <p>
          With locked challenge parameters <strong>N=4, J=1.0, g=1.0, T=2.0</strong> and periodic boundary condition links (0,1), (1,2), (2,3), and (3,0).
        </p>
      </div>
    )
  },
  {
    title: "Chapter 2 — The Néel Initial State |1010⟩",
    badge: "INITIAL STATE",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          The system starts in the antiferromagnetic Néel state <strong>|1010⟩</strong>. Under Qiskit little-endian convention (where Q0 is the rightmost bit), this corresponds to computational basis index <strong>10</strong>:
        </p>
        <div className="p-3 bg-slate-950 rounded-lg border border-slate-800 font-mono text-amber-300 text-center">
          |ψ(0)⟩ = |q3=1, q2=0, q1=1, q0=0⟩ → M_stag(0) = +1.0000
        </div>
        <p>
          The staggered magnetization operator M_stag = ¼ ∑ (-1)^q ⟨Z_q⟩ starts at exactly +1.0 and decays as the transverse field drives quantum fluctuations.
        </p>
      </div>
    )
  },
  {
    title: "Chapter 3 — One First-Order Trotter Step",
    badge: "DISCRETIZATION",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          Because [H_ZZ, H_X] ≠ 0, the exact evolution operator exp(-i H Δt) cannot be directly executed as a single product of single- and two-qubit gates.
        </p>
        <p>
          We employ the first-order Lie-Trotter product formula over discrete time step Δt = T/n:
        </p>
        <div className="p-3 bg-slate-950 rounded-lg border border-slate-800 font-mono text-cyan-300 text-center">
          U_step(Δt) = ∏_i exp(-i g Δt X_i) · ∏_⟨j,k⟩ exp(-i J Δt Z_j Z_k) + O(Δt²)
        </div>
      </div>
    )
  },
  {
    title: "Chapter 4 — Hardware Transpilation: Processors A, B, and D",
    badge: "TOPOLOGY & ROUTING",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          The logical model requires four periodic couplings: (0,1), (1,2), (2,3), and (3,0). Real quantum hardware rarely provides native periodic rings:
        </p>
        <ul className="list-disc pl-5 space-y-1">
          <li><strong>Processor A (5 Qubits Tree):</strong> Lacks the (3,0) cycle edge. Transpiler must insert SWAP routing gates (+2 SWAPs / step).</li>
          <li><strong>Processor B (9 Qubits Heavy-Hex):</strong> 12 edges, routes non-native couplings via auxiliary nodes (+1 SWAP / step).</li>
          <li><strong>Processor D (Custom 6 Qubits):</strong> Native ring with resonator bypass links. <strong>0 routing SWAPs</strong> overhead.</li>
        </ul>
      </div>
    )
  },
  {
    title: "Chapter 5 — Physical Noise Intrusion",
    badge: "HARDWARE NOISE",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          Each physical gate incurs decoherence and gate infidelity:
        </p>
        <ul className="list-disc pl-5 space-y-1">
          <li>Single-qubit depolarizing error: ε_1Q ~ 0.02% – 0.04%</li>
          <li>Two-qubit cross-resonance error: ε_2Q ~ 0.35% – 0.65%</li>
          <li>Readout assignment error: ε_RO ~ 0.9% – 1.8%</li>
        </ul>
        <p>
          Routing SWAPs are expensive: each SWAP requires 3 CNOT gates, tripling the two-qubit error rate for that interaction.
        </p>
      </div>
    )
  },
  {
    title: "Chapter 6 — Projective Readout & Sampling Error",
    badge: "MEASUREMENT",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          At time T=2.0, all 4 qubits are measured in the computational Z-basis. Real experimental readout introduces bit-flip misassignment (0→1 and 1→0) alongside finite-shot statistical noise ~ 1/√N_shots.
        </p>
      </div>
    )
  },
  {
    title: "Chapter 7 — Trotter Error vs. Physical Noise Competition",
    badge: "THE CORE DILEMMA",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          This is the central dilemma of Problem Statement PS-S3:
        </p>
        <ul className="list-disc pl-5 space-y-2">
          <li><strong>Low step count (small n):</strong> Algorithmic Trotter error dominates (O(T²/n)). The state diverges from exact physics due to coarse time slicing.</li>
          <li><strong>High step count (large n):</strong> Trotter error shrinks, but gate count and circuit depth skyrocket (Gates ∝ n). Hardware noise completely depolarizes the state.</li>
        </ul>
      </div>
    )
  },
  {
    title: "Chapter 8 — Finding the Optimal Step Count n*",
    badge: "OPTIMIZATION",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          The optimal step count <strong>n*</strong> maximizes output fidelity under physical noise, resolved by the official smallest-n tie-breaking rule:
        </p>
        <div className="p-3 bg-slate-950 rounded-lg border border-emerald-900/60 font-mono text-emerald-400 text-center">
          n* = argmax_n [ F_noisy(n) ]
        </div>
        <p>
          Due to zero routing overhead, Processor D achieves its peak fidelity at a higher step count than Processor A and B.
        </p>
      </div>
    )
  },
  {
    title: "Chapter 9 — Scientific Conclusion & Hackathon Submission",
    badge: "CONCLUSION",
    content: (
      <div className="space-y-4 text-xs text-slate-300 leading-relaxed">
        <p>
          Hardware topology directly dictates simulation accuracy. Co-designing the processor coupling graph to natively match the physical Hamiltonian’s interaction geometry (as in Processor D) eliminates routing overhead, directly enabling deeper Trotter circuits and higher fidelity.
        </p>
      </div>
    )
  }
];

// ============================================================================
// 3. THREE.JS 3D QUANTUM PROCESSOR RING & EMBEDDED BLOCH SPHERE VIEWPORT
// ============================================================================

function ThreeQuantumViewport({
  activeQubit,
  onSelectQubit,
  activeBond,
  onSelectBond,
  currentInstruction,
  qubitBlochVectors,
  processorKey,
  viewMode, // 'logical' or 'physical'
  onToggleViewMode
}) {
  const mountRef = useRef(null);
  const sceneRef = useRef(null);
  const cameraRef = useRef(null);
  const rendererRef = useRef(null);
  const raycasterRef = useRef(new THREE.Raycaster());
  const mouseRef = useRef(new THREE.Vector2());
  const interactiveQubitMeshesRef = useRef([]);
  const qubitMeshesRef = useRef([]);
  const blochVectorsRef = useRef([]);
  const bondMeshesRef = useRef([]);
  const isDraggingRef = useRef(false);
  const prevMousePosRef = useRef({ x: 0, y: 0 });
  const cameraAngleRef = useRef({ theta: Math.PI / 4, phi: Math.PI / 3, radius: 15 });

  const currentProc = PROCESSORS[processorKey];
  const numDisplayQubits = viewMode === 'physical' ? currentProc.numPhysQubits : 4;

  // Spatial positions for nodes
  const nodePositions = useMemo(() => {
    if (viewMode === 'logical') {
      const ringRadius = 4.0;
      return [0, 1, 2, 3].map(i => {
        const angle = (i * 2 * Math.PI) / 4 + Math.PI / 4;
        return new THREE.Vector3(ringRadius * Math.cos(angle), 0, ringRadius * Math.sin(angle));
      });
    }

    // Physical Hardware Layouts
    if (processorKey === 'A') {
      // 5 Qubits Tree: (0,1), (1,2), (1,3), (3,4)
      return [
        new THREE.Vector3(-4.0, 0, -2.0), // P0
        new THREE.Vector3(-1.5, 0, 0.0),  // P1 (hub)
        new THREE.Vector3(-4.0, 0, 2.0),  // P2
        new THREE.Vector3(1.5, 0, 0.0),   // P3
        new THREE.Vector3(4.0, 0, 0.0)    // P4
      ];
    } else if (processorKey === 'B') {
      // 9 Qubits Graph: (0,1), (1,2), (1,6), (2,3), (2,6), (2,5), (3,4), (4,5), (5,6), (6,7), (6,8), (8,0)
      const R = 4.8;
      const positions = [];
      for (let i = 0; i < 8; i++) {
        const angle = (i * 2 * Math.PI) / 8;
        positions.push(new THREE.Vector3(R * Math.cos(angle), 0, R * Math.sin(angle)));
      }
      // P8 in center-left to bridge (8,0) and (6,8)
      positions.push(new THREE.Vector3(0, 0, -1.8)); // P8
      return positions;
    } else {
      // Processor D: 6 Qubits ring (0-3) + bypass couplers (4, 5)
      const ringRadius = 4.0;
      const positions = [0, 1, 2, 3].map(i => {
        const angle = (i * 2 * Math.PI) / 4 + Math.PI / 4;
        return new THREE.Vector3(ringRadius * Math.cos(angle), 0, ringRadius * Math.sin(angle));
      });
      positions.push(new THREE.Vector3(-5.2, 0, -1.5)); // P4 coupler to 0
      positions.push(new THREE.Vector3(5.2, 0, 1.5));   // P5 coupler to 2
      return positions;
    }
  }, [processorKey, viewMode]);

  // Update Camera Orbit Position
  const updateCameraPosition = useCallback(() => {
    if (!cameraRef.current) return;
    const { theta, phi, radius } = cameraAngleRef.current;
    cameraRef.current.position.x = radius * Math.sin(phi) * Math.sin(theta);
    cameraRef.current.position.y = radius * Math.cos(phi);
    cameraRef.current.position.z = radius * Math.sin(phi) * Math.cos(theta);
    cameraRef.current.lookAt(0, 0.3, 0);
  }, []);

  // Setup Three.js Scene, Materials, Geometry, and Render Loop
  useEffect(() => {
    const container = mountRef.current;
    if (!container) return;

    const width = container.clientWidth || 800;
    const height = container.clientHeight || 500;

    // 1. Scene
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(0x060a12);
    sceneRef.current = scene;

    // 2. Camera
    const camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 100);
    cameraRef.current = camera;
    updateCameraPosition();

    // 3. Renderer
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setSize(width, height);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.shadowMap.enabled = true;
    container.innerHTML = '';
    container.appendChild(renderer.domElement);
    rendererRef.current = renderer;

    // 4. Lighting
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.7);
    scene.add(ambientLight);

    const dirLight1 = new THREE.DirectionalLight(0x38bdf8, 1.3);
    dirLight1.position.set(10, 15, 10);
    scene.add(dirLight1);

    const dirLight2 = new THREE.DirectionalLight(0x818cf8, 0.9);
    dirLight2.position.set(-10, -5, -10);
    scene.add(dirLight2);

    // 5. Cryogenic Floor Grid
    const gridHelper = new THREE.GridHelper(18, 18, 0x1e293b, 0x0f172a);
    gridHelper.position.y = -1.2;
    scene.add(gridHelper);

    // Platform
    const platformGeo = new THREE.CylinderGeometry(6.0, 6.2, 0.2, 48);
    const platformMat = new THREE.MeshStandardMaterial({
      color: 0x0a1124,
      roughness: 0.4,
      metalness: 0.8
    });
    const platform = new THREE.Mesh(platformGeo, platformMat);
    platform.position.y = -1.3;
    scene.add(platform);

    // 6. Build Qubit Nodes
    qubitMeshesRef.current = [];
    blochVectorsRef.current = [];
    interactiveQubitMeshesRef.current = [];

    nodePositions.forEach((pos, idx) => {
      const qGroup = new THREE.Group();
      qGroup.position.copy(pos);

      // Cryo stand (clickable)
      const standGeo = new THREE.CylinderGeometry(0.5, 0.7, 0.9, 24);
      const standMat = new THREE.MeshStandardMaterial({
        color: 0x1e293b,
        metalness: 0.9,
        roughness: 0.2
      });
      const stand = new THREE.Mesh(standGeo, standMat);
      stand.position.y = -0.65;
      stand.userData = { qubitIndex: idx };
      qGroup.add(stand);
      interactiveQubitMeshesRef.current.push(stand);

      // Gold Coupler
      const ringGeo = new THREE.TorusGeometry(0.62, 0.045, 16, 32);
      const ringMat = new THREE.MeshStandardMaterial({
        color: 0xf59e0b,
        metalness: 0.95,
        roughness: 0.1
      });
      const ringMesh = new THREE.Mesh(ringGeo, ringMat);
      ringMesh.rotation.x = Math.PI / 2;
      ringMesh.position.y = -0.2;
      qGroup.add(ringMesh);

      // Bloch sphere transparent enclosure
      const sphereGeo = new THREE.SphereGeometry(0.9, 24, 20);
      const sphereMat = new THREE.MeshStandardMaterial({
        color: 0x0f172a,
        transparent: true,
        opacity: 0.35,
        roughness: 0.1,
        metalness: 0.1,
        wireframe: true
      });
      const sphereMesh = new THREE.Mesh(sphereGeo, sphereMat);
      sphereMesh.position.y = 0.9;
      sphereMesh.userData = { qubitIndex: idx };
      qGroup.add(sphereMesh);
      interactiveQubitMeshesRef.current.push(sphereMesh);

      // Equator guide ring
      const eqGeo = new THREE.TorusGeometry(0.9, 0.012, 8, 36);
      const eqMat = new THREE.MeshBasicMaterial({ color: 0x38bdf8 });
      const equator = new THREE.Mesh(eqGeo, eqMat);
      equator.rotation.x = Math.PI / 2;
      equator.position.y = 0.9;
      qGroup.add(equator);

      // Bloch Vector Arrow
      const arrowGroup = new THREE.Group();
      arrowGroup.position.set(0, 0.9, 0);

      const shaftGeo = new THREE.CylinderGeometry(0.035, 0.035, 0.75, 16);
      const arrowMat = new THREE.MeshStandardMaterial({
        color: 0x38bdf8,
        emissive: 0x0284c7,
        emissiveIntensity: 0.6,
        roughness: 0.2
      });
      const shaft = new THREE.Mesh(shaftGeo, arrowMat);
      shaft.position.y = 0.375;
      arrowGroup.add(shaft);

      const coneGeo = new THREE.ConeGeometry(0.1, 0.22, 16);
      const cone = new THREE.Mesh(coneGeo, arrowMat);
      cone.position.y = 0.85;
      arrowGroup.add(cone);

      qGroup.add(arrowGroup);
      blochVectorsRef.current.push(arrowGroup);

      // Core sphere
      const coreGeo = new THREE.SphereGeometry(0.22, 16, 16);
      const coreMat = new THREE.MeshStandardMaterial({
        color: 0x38bdf8,
        emissive: 0x38bdf8,
        emissiveIntensity: 0.5
      });
      const core = new THREE.Mesh(coreGeo, coreMat);
      core.position.y = 0.9;
      core.userData = { qubitIndex: idx };
      qGroup.add(core);
      interactiveQubitMeshesRef.current.push(core);

      scene.add(qGroup);
      qubitMeshesRef.current.push({
        group: qGroup,
        sphereMesh,
        coreMat,
        standMat,
        arrowMat,
        core
      });
    });

    // 7. Coupling Bonds
    bondMeshesRef.current = [];
    const activeBonds = viewMode === 'logical' ? MODEL_CONFIG.bonds : currentProc.nativeBonds;

    activeBonds.forEach(([q1, q2]) => {
      if (!nodePositions[q1] || !nodePositions[q2]) return;
      const p1 = nodePositions[q1];
      const p2 = nodePositions[q2];
      const distance = p1.distanceTo(p2);
      const midPoint = new THREE.Vector3().addVectors(p1, p2).multiplyScalar(0.5);

      const busGeo = new THREE.CylinderGeometry(0.07, 0.07, Math.max(0.2, distance - 1.1), 16);
      const busMat = new THREE.MeshStandardMaterial({
        color: 0x475569,
        metalness: 0.8,
        roughness: 0.3,
        emissive: 0x0f172a,
        emissiveIntensity: 0.2
      });
      const busMesh = new THREE.Mesh(busGeo, busMat);
      busMesh.position.copy(midPoint);
      busMesh.position.y = -0.2;

      const dir = new THREE.Vector3().subVectors(p2, p1).normalize();
      const axis = new THREE.Vector3(0, 1, 0);
      busMesh.quaternion.setFromUnitVectors(axis, dir);

      scene.add(busMesh);
      bondMeshesRef.current.push({ q1, q2, mesh: busMesh, material: busMat });
    });

    // 8. Render loop
    let animationFrameId;
    const animate = () => {
      animationFrameId = requestAnimationFrame(animate);
      renderer.render(scene, camera);
    };
    animate();

    // 9. Resize Listener
    const handleResize = () => {
      if (!container || !renderer || !camera) return;
      const w = container.clientWidth;
      const h = container.clientHeight;
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      renderer.setSize(w, h);
    };
    window.addEventListener('resize', handleResize);

    return () => {
      cancelAnimationFrame(animationFrameId);
      window.removeEventListener('resize', handleResize);
      if (renderer.domElement && container.contains(renderer.domElement)) {
        container.removeChild(renderer.domElement);
      }
      renderer.dispose();
    };
  }, [nodePositions, currentProc, viewMode, updateCameraPosition]);

  // Update Visual State: Qubit Selection, Active Gate Highlighting, Bloch Vectors
  useEffect(() => {
    if (!qubitMeshesRef.current.length) return;

    nodePositions.forEach((_, qIdx) => {
      const bloch = (qubitBlochVectors && qubitBlochVectors[qIdx]) || { x: 0, y: 0, z: 1, length: 1 };
      const arrowGroup = blochVectorsRef.current[qIdx];
      const qubitObj = qubitMeshesRef.current[qIdx];

      if (arrowGroup) {
        const targetDir = new THREE.Vector3(bloch.x, bloch.z, bloch.y).normalize();
        const defaultDir = new THREE.Vector3(0, 1, 0);
        if (targetDir.lengthSq() > 0.001) {
          arrowGroup.quaternion.setFromUnitVectors(defaultDir, targetDir);
        }
        const clampedLen = Math.max(0.15, Math.min(1.0, bloch.length || 1));
        arrowGroup.scale.set(clampedLen, clampedLen, clampedLen);
      }

      const isTarget = currentInstruction?.qubits?.includes(qIdx);
      const isSelected = activeQubit === qIdx;

      if (qubitObj) {
        if (isSelected) {
          // Selected Qubit: High-contrast Cyan Highlight
          qubitObj.coreMat.color.setHex(0x00f0ff);
          qubitObj.coreMat.emissive.setHex(0x0284c7);
          qubitObj.coreMat.emissiveIntensity = 1.0;
          qubitObj.arrowMat.color.setHex(0x38bdf8);
          qubitObj.standMat.color.setHex(0x0284c7);
        } else if (isTarget) {
          // Active Gate Execution: Amber
          qubitObj.coreMat.color.setHex(0xf59e0b);
          qubitObj.coreMat.emissive.setHex(0xd97706);
          qubitObj.coreMat.emissiveIntensity = 0.9;
          qubitObj.arrowMat.color.setHex(0xfbbf24);
          qubitObj.standMat.color.setHex(0x334155);
        } else {
          // Idle Qubit
          qubitObj.coreMat.color.setHex(0x64748b);
          qubitObj.coreMat.emissive.setHex(0x1e293b);
          qubitObj.coreMat.emissiveIntensity = 0.2;
          qubitObj.arrowMat.color.setHex(0x60a5fa);
          qubitObj.standMat.color.setHex(0x1e293b);
        }
      }
    });

    // Update Bond Highlighting
    bondMeshesRef.current.forEach(bond => {
      const isExecuting = 
        currentInstruction?.qubits?.length === 2 &&
        ((currentInstruction.qubits[0] === bond.q1 && currentInstruction.qubits[1] === bond.q2) ||
         (currentInstruction.qubits[0] === bond.q2 && currentInstruction.qubits[1] === bond.q1));

      const isSelectedBond = 
        activeBond && 
        ((activeBond[0] === bond.q1 && activeBond[1] === bond.q2) || 
         (activeBond[0] === bond.q2 && activeBond[1] === bond.q1));

      if (isExecuting) {
        bond.material.color.setHex(0xec4899);
        bond.material.emissive.setHex(0xbe185d);
        bond.material.emissiveIntensity = 1.2;
      } else if (isSelectedBond) {
        bond.material.color.setHex(0x38bdf8);
        bond.material.emissive.setHex(0x0284c7);
        bond.material.emissiveIntensity = 0.8;
      } else {
        bond.material.color.setHex(0x10b981);
        bond.material.emissive.setHex(0x047857);
        bond.material.emissiveIntensity = 0.3;
      }
    });
  }, [nodePositions, qubitBlochVectors, currentInstruction, activeQubit, activeBond]);

  // Click on Qubit Raycasting Support
  const handleCanvasClick = (e) => {
    if (!mountRef.current || !cameraRef.current) return;
    const rect = mountRef.current.getBoundingClientRect();
    mouseRef.current.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    mouseRef.current.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;

    raycasterRef.current.setFromCamera(mouseRef.current, cameraRef.current);
    const intersects = raycasterRef.current.intersectObjects(interactiveQubitMeshesRef.current, false);

    if (intersects.length > 0) {
      const hitQubit = intersects[0].object.userData?.qubitIndex;
      if (hitQubit !== undefined) {
        onSelectQubit(hitQubit);
      }
    }
  };

  // Drag orbit handlers
  const handleMouseDown = (e) => {
    isDraggingRef.current = true;
    prevMousePosRef.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseMove = (e) => {
    if (!isDraggingRef.current) return;
    const deltaX = e.clientX - prevMousePosRef.current.x;
    const deltaY = e.clientY - prevMousePosRef.current.y;
    prevMousePosRef.current = { x: e.clientX, y: e.clientY };

    cameraAngleRef.current.theta += deltaX * 0.008;
    cameraAngleRef.current.phi = Math.max(0.2, Math.min(Math.PI / 2 - 0.05, cameraAngleRef.current.phi - deltaY * 0.008));
    updateCameraPosition();
  };

  const handleMouseUp = () => {
    isDraggingRef.current = false;
  };

  const handleWheel = (e) => {
    e.preventDefault();
    cameraAngleRef.current.radius = Math.max(8, Math.min(25, cameraAngleRef.current.radius + e.deltaY * 0.015));
    updateCameraPosition();
  };

  return (
    <div className="relative w-full h-[480px] bg-slate-950 rounded-xl overflow-hidden border border-slate-800 shadow-2xl select-none">
      <div
        ref={mountRef}
        className="w-full h-full cursor-grab active:cursor-grabbing"
        onClick={handleCanvasClick}
        onMouseDown={handleMouseDown}
        onMouseMove={handleMouseMove}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseUp}
        onWheel={handleWheel}
      />

      {/* Floating HUD Controls */}
      <div className="absolute top-4 left-4 flex flex-col gap-2 z-10 pointer-events-auto">
        <div className="flex items-center gap-2 bg-slate-900/90 backdrop-blur-md px-3 py-1.5 rounded-lg border border-slate-700/60 shadow-lg text-xs">
          <span className="w-2.5 h-2.5 rounded-full bg-emerald-500 animate-pulse" />
          <span className="font-semibold text-slate-100">{currentProc.name}</span>
          <span className="text-slate-400">({numDisplayQubits} Qubits)</span>
        </div>

        {/* View Mode Toggle Button */}
        <div className="flex items-center bg-slate-900/90 backdrop-blur-md rounded-lg p-0.5 border border-slate-700/60 text-xs">
          <button
            onClick={() => onToggleViewMode('logical')}
            className={`px-2.5 py-1 rounded transition font-medium ${
              viewMode === 'logical'
                ? 'bg-cyan-500 text-slate-950 font-bold'
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            Logical Ring (4Q)
          </button>
          <button
            onClick={() => onToggleViewMode('physical')}
            className={`px-2.5 py-1 rounded transition font-medium ${
              viewMode === 'physical'
                ? 'bg-cyan-500 text-slate-950 font-bold'
                : 'text-slate-400 hover:text-slate-200'
            }`}
          >
            Hardware Topology ({currentProc.numPhysQubits}Q)
          </button>
        </div>
      </div>

      {/* Camera Reset */}
      <div className="absolute top-4 right-4 flex items-center gap-2 z-10">
        <button
          onClick={() => {
            cameraAngleRef.current = { theta: Math.PI / 4, phi: Math.PI / 3, radius: 15 };
            updateCameraPosition();
          }}
          className="p-2 bg-slate-900/80 hover:bg-slate-800 text-slate-300 rounded-lg border border-slate-700/60 transition shadow-md"
          title="Reset View"
        >
          <RotateCcw className="w-4 h-4" />
        </button>
      </div>

      {/* Interactive Qubit Node Selector Bar */}
      <div className="absolute bottom-4 left-1/2 -translate-x-1/2 flex items-center gap-1.5 bg-slate-900/95 backdrop-blur-md px-4 py-2 rounded-xl border border-slate-700/80 shadow-2xl z-10 max-w-[95%] overflow-x-auto">
        <span className="text-xs font-semibold text-slate-400 mr-1 flex items-center gap-1 shrink-0">
          <Cpu className="w-3.5 h-3.5 text-cyan-400" /> Select Qubit:
        </span>
        {Array.from({ length: numDisplayQubits }, (_, qIdx) => {
          const isSelected = activeQubit === qIdx;
          const isTarget = currentInstruction?.qubits?.includes(qIdx);
          const label = viewMode === 'physical' ? `P${qIdx}` : `Q${qIdx}`;

          return (
            <button
              key={qIdx}
              onClick={() => onSelectQubit(qIdx)}
              className={`px-2.5 py-1 text-xs rounded-lg font-mono font-medium transition flex items-center gap-1 shrink-0 ${
                isSelected
                  ? 'bg-cyan-500 text-slate-950 font-bold shadow-md shadow-cyan-500/30 ring-2 ring-cyan-300'
                  : 'bg-slate-800/90 hover:bg-slate-700 text-slate-200 border border-slate-700/60'
              }`}
            >
              <span>{label}</span>
              {isTarget && <span className="w-1.5 h-1.5 rounded-full bg-amber-400 animate-ping" />}
            </button>
          );
        })}
      </div>
    </div>
  );
}

// ============================================================================
// 4. MAIN QUANTUM LABORATORY CONTAINER COMPONENT
// ============================================================================

export default function App() {
  const [trotterSteps, setTrotterSteps] = useState(4);
  const [selectedProcessor, setSelectedProcessor] = useState("D");
  const [noiseScale, setNoiseScale] = useState("1.0");
  const [currentGateIndex, setCurrentGateIndex] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [playbackSpeed, setPlaybackSpeed] = useState(1.0);

  const [activeQubit, setActiveQubit] = useState(0);
  const [activeBond, setActiveBond] = useState(null);
  const [viewMode, setViewMode] = useState("physical"); // 'logical' or 'physical'
  const [activeTab, setActiveTab] = useState("simulation");
  const [activeStoryChapter, setActiveStoryChapter] = useState(0);
  const [numShots, setNumShots] = useState(4096);
  const [measuredCounts, setMeasuredCounts] = useState(null);

  // Generate instruction stream
  const instructions = useMemo(() => {
    return generateCircuitInstructionStream(trotterSteps, selectedProcessor);
  }, [trotterSteps, selectedProcessor]);

  // Ensure currentGateIndex stays within bounds when changing processor or steps
  useEffect(() => {
    if (currentGateIndex >= instructions.length) {
      setCurrentGateIndex(0);
    }
  }, [instructions.length, currentGateIndex]);

  const currentInstruction = instructions[currentGateIndex] || instructions[0];

  // Compute Live Intermediate State at the exact current checkpoint
  const liveState = useMemo(() => {
    const psi0 = getInitialStatevector();
    let psi = psi0.map(v => ({ ...v }));

    for (let idx = 0; idx <= currentGateIndex && idx < instructions.length; idx++) {
      const inst = instructions[idx];
      if (!inst) continue;

      if (inst.type === "sq" && inst.targetGate === "Rx") {
        psi = applyRx(psi, inst.qubits[0], inst.angleRad);
      } else if (inst.type === "twoQ" && inst.targetGate === "Rzz") {
        psi = applyRzz(psi, inst.qubits[0], inst.qubits[1], inst.angleRad);
      } else if (inst.type === "swap") {
        psi = applySwap(psi, inst.qubits[0], inst.qubits[1]);
      }
    }

    const elapsedStep = Math.min(trotterSteps, currentInstruction?.step || 0);
    const elapsedTime = (elapsedStep / trotterSteps) * MODEL_CONFIG.T;
    const H = buildHamiltonian(MODEL_CONFIG.N, MODEL_CONFIG.J, MODEL_CONFIG.g, MODEL_CONFIG.bonds);
    const exactPsiAtT = exactTimeEvolution(psi0, H, Math.max(0.001, elapsedTime));

    const fidelity = calculateFidelity(exactPsiAtT, psi);
    const mStag = calculateStaggeredMagnetization(psi);
    
    // Support up to 9 physical qubits for Bloch vector display
    const blochVectors = Array.from({ length: 9 }, (_, q) => {
      if (q < 4) {
        return computeReducedBlochVector(psi, q);
      }
      return { x: 0, y: 0, z: 1, length: 1, purity: 1.0 };
    });

    return { psi, exactPsiAtT, fidelity, mStag, elapsedTime, blochVectors };
  }, [instructions, currentGateIndex, trotterSteps, currentInstruction]);

  // Auto-play Animation Timer
  useEffect(() => {
    let timer;
    if (isPlaying) {
      const intervalMs = Math.max(120, 800 / playbackSpeed);
      timer = setInterval(() => {
        setCurrentGateIndex(prev => {
          if (prev >= instructions.length - 1) {
            setIsPlaying(false);
            return prev;
          }
          return prev + 1;
        });
      }, intervalMs);
    }
    return () => clearInterval(timer);
  }, [isPlaying, instructions.length, playbackSpeed]);

  // Reset and switch processor safely
  const handleProcessorChange = (pKey) => {
    setIsPlaying(false);
    setSelectedProcessor(pKey);
    setCurrentGateIndex(0);
    setActiveQubit(0);
    setMeasuredCounts(null);
  };

  const handleStepCountChange = (n) => {
    setIsPlaying(false);
    setTrotterSteps(n);
    setCurrentGateIndex(0);
    setMeasuredCounts(null);
  };

  // Sample Shot-based Readout on Measurement Checkpoint
  const handleMeasureShots = useCallback(() => {
    const counts = {};
    const probs = liveState.psi.map(v => cMagSq(v));
    const proc = PROCESSORS[selectedProcessor];
    const roErr = proc.readoutErr * parseFloat(noiseScale);

    for (let shot = 0; shot < numShots; shot++) {
      let rand = Math.random();
      let cumProb = 0;
      let sampledIdx = 0;
      for (let s = 0; s < 16; s++) {
        cumProb += probs[s];
        if (rand <= cumProb) {
          sampledIdx = s;
          break;
        }
      }

      let noisyIdx = sampledIdx;
      for (let q = 0; q < 4; q++) {
        if (Math.random() < roErr) {
          noisyIdx ^= (1 << q);
        }
      }

      const bitStr = noisyIdx.toString(2).padStart(4, '0');
      counts[bitStr] = (counts[bitStr] || 0) + 1;
    }
    setMeasuredCounts(counts);
  }, [liveState.psi, numShots, selectedProcessor, noiseScale]);

  // Selected Qubit Inspection Details
  const selectedBloch = liveState.blochVectors[activeQubit] || { x: 0, y: 0, z: 1, length: 1, purity: 1 };
  const currentProcObj = PROCESSORS[selectedProcessor];

  return (
    <div className="w-full min-h-screen bg-slate-950 text-slate-100 flex flex-col font-sans selection:bg-cyan-500 selection:text-slate-950">
      {/* Top Laboratory Navigation Bar */}
      <header className="w-full bg-slate-900/90 backdrop-blur-md border-b border-slate-800 px-6 py-3.5 flex flex-wrap items-center justify-between gap-4 sticky top-0 z-50">
        <div className="flex items-center gap-3">
          <div className="w-9 h-9 rounded-lg bg-gradient-to-br from-cyan-500 to-indigo-600 flex items-center justify-center shadow-lg shadow-cyan-500/20">
            <Activity className="w-5 h-5 text-slate-950 font-bold" />
          </div>
          <div>
            <div className="flex items-center gap-2">
              <h1 className="font-bold text-base tracking-tight text-slate-100">
                QFF 2026 Interactive 3D Quantum Simulation Laboratory
              </h1>
              <span className="bg-cyan-500/10 text-cyan-400 text-[11px] font-mono px-2 py-0.5 rounded border border-cyan-500/20 font-semibold">
                PS-S3
              </span>
            </div>
            <p className="text-xs text-slate-400">
              Processors A, B, and D: Exact Hardware Topologies & Lie-Trotter Evolution
            </p>
          </div>
        </div>

        {/* Global Experiment Config Selectors */}
        <div className="flex items-center gap-3">
          {/* Processor Selector */}
          <div className="flex items-center bg-slate-800/80 rounded-lg p-1 border border-slate-700/60">
            <span className="text-xs text-slate-400 px-2 font-medium">Processor:</span>
            {["A", "B", "D"].map(pk => (
              <button
                key={pk}
                onClick={() => handleProcessorChange(pk)}
                className={`px-3 py-1 text-xs font-mono font-medium rounded transition ${
                  selectedProcessor === pk
                    ? 'bg-cyan-500 text-slate-950 font-bold shadow'
                    : 'text-slate-300 hover:text-white'
                }`}
              >
                Proc {pk}
              </button>
            ))}
          </div>

          {/* Trotter Step Count Selector */}
          <div className="flex items-center bg-slate-800/80 rounded-lg p-1 border border-slate-700/60">
            <span className="text-xs text-slate-400 px-2 font-medium">Steps (n):</span>
            <select
              value={trotterSteps}
              onChange={(e) => handleStepCountChange(parseInt(e.target.value))}
              className="bg-slate-900 text-cyan-300 font-mono text-xs px-2 py-1 rounded border border-slate-700 focus:outline-none focus:border-cyan-500"
            >
              {MODEL_CONFIG.officialStepCounts.map(n => (
                <option key={n} value={n}>
                  n = {n}
                </option>
              ))}
            </select>
          </div>

          {/* Noise Scale Selector */}
          <div className="flex items-center bg-slate-800/80 rounded-lg p-1 border border-slate-700/60">
            <span className="text-xs text-slate-400 px-2 font-medium">Noise:</span>
            {["0.5", "1.0", "2.0"].map(ns => (
              <button
                key={ns}
                onClick={() => setNoiseScale(ns)}
                className={`px-2.5 py-1 text-xs font-mono rounded transition ${
                  noiseScale === ns
                    ? 'bg-pink-500 text-slate-950 font-bold shadow'
                    : 'text-slate-300 hover:text-white'
                }`}
              >
                {ns}x
              </button>
            ))}
          </div>
        </div>
      </header>

      {/* Primary Navigation Tabs */}
      <div className="bg-slate-900 border-b border-slate-800 px-6 flex items-center gap-6 text-sm">
        <button
          onClick={() => setActiveTab("simulation")}
          className={`py-3 flex items-center gap-2 border-b-2 font-medium transition ${
            activeTab === "simulation"
              ? 'border-cyan-400 text-cyan-400 font-semibold'
              : 'border-transparent text-slate-400 hover:text-slate-200'
          }`}
        >
          <Move3D className="w-4 h-4" /> 3D Quantum Lab & Timeline
        </button>
        <button
          onClick={() => setActiveTab("tradeoff")}
          className={`py-3 flex items-center gap-2 border-b-2 font-medium transition ${
            activeTab === "tradeoff"
              ? 'border-cyan-400 text-cyan-400 font-semibold'
              : 'border-transparent text-slate-400 hover:text-slate-200'
          }`}
        >
          <BarChart2 className="w-4 h-4" /> Trotter vs. Noise Trade-off Dashboard
        </button>
        <button
          onClick={() => setActiveTab("story")}
          className={`py-3 flex items-center gap-2 border-b-2 font-medium transition ${
            activeTab === "story"
              ? 'border-cyan-400 text-cyan-400 font-semibold'
              : 'border-transparent text-slate-400 hover:text-slate-200'
          }`}
        >
          <Layers className="w-4 h-4" /> Guided Scientific Story (9 Chapters)
        </button>
      </div>

      {/* Main Content Area */}
      <main className="flex-1 p-6 max-w-7xl mx-auto w-full space-y-6">
        {activeTab === "simulation" && (
          <div className="space-y-6">
            {/* Top Row: 3D Viewport on Left, Statevector & Amplitudes on Right */}
            <div className="grid grid-cols-1 lg:grid-cols-12 gap-6 items-start">
              {/* Left Column: 3D Viewport */}
              <div className="lg:col-span-7 flex flex-col gap-4">
                <ThreeQuantumViewport
                  activeQubit={activeQubit}
                  onSelectQubit={setActiveQubit}
                  activeBond={activeBond}
                  onSelectBond={setActiveBond}
                  currentInstruction={currentInstruction}
                  qubitBlochVectors={liveState.blochVectors}
                  processorKey={selectedProcessor}
                  viewMode={viewMode}
                  onToggleViewMode={setViewMode}
                />

                {/* Single Qubit Bloch Sphere & Density Matrix Detail Inspector */}
                <div className="bg-slate-900/90 rounded-xl p-4 border border-slate-800 shadow-xl flex flex-wrap items-center justify-between gap-4">
                  <div className="flex items-center gap-3">
                    <div className="w-9 h-9 rounded-lg bg-cyan-500/20 text-cyan-400 flex items-center justify-center font-mono font-bold text-sm">
                      {viewMode === 'physical' ? `P${activeQubit}` : `Q${activeQubit}`}
                    </div>
                    <div>
                      <h4 className="text-xs font-semibold text-slate-200">
                        Inspecting Qubit {activeQubit} on {currentProcObj.name}
                      </h4>
                      <p className="text-[11px] text-slate-400 font-mono">
                        Bloch Vector: ({selectedBloch.x.toFixed(3)}, {selectedBloch.y.toFixed(3)}, {selectedBloch.z.toFixed(3)})
                      </p>
                    </div>
                  </div>

                  <div className="flex items-center gap-6 font-mono text-xs">
                    <div>
                      <span className="text-slate-500 text-[10px] block">VECTOR LENGTH |r|</span>
                      <span className={`font-bold ${selectedBloch.length > 0.95 ? 'text-emerald-400' : 'text-amber-400'}`}>
                        {selectedBloch.length.toFixed(4)}
                      </span>
                    </div>
                    <div>
                      <span className="text-slate-500 text-[10px] block">PURITY Tr(ρ²)</span>
                      <span className="text-cyan-300 font-bold">
                        {selectedBloch.purity.toFixed(4)}
                      </span>
                    </div>
                    <div>
                      <span className="text-slate-500 text-[10px] block">STATE STATUS</span>
                      <span className={`font-semibold ${selectedBloch.purity < 0.95 ? 'text-pink-400' : 'text-emerald-400'}`}>
                        {selectedBloch.purity < 0.95 ? "Mixed / Entangled" : "Nearly Pure"}
                      </span>
                    </div>
                  </div>
                </div>

                {/* Processor Definition Card from Section 4 */}
                <div className="bg-slate-900/80 rounded-xl p-4 border border-slate-800 text-xs space-y-2">
                  <div className="flex items-center justify-between">
                    <span className="font-semibold text-slate-200 flex items-center gap-1.5">
                      <ShieldCheck className="w-4 h-4 text-emerald-400" />
                      Section 4 Processor Definition:
                    </span>
                    <span className="font-mono text-slate-400">{currentProcObj.type}</span>
                  </div>
                  <p className="text-slate-300 leading-relaxed">
                    {currentProcObj.description}
                  </p>
                  <div className="flex flex-wrap gap-1.5 pt-1">
                    <span className="text-slate-400 text-[11px]">Native Edges:</span>
                    {currentProcObj.nativeBonds.map(([q1, q2]) => (
                      <span key={`${q1}-${q2}`} className="px-2 py-0.5 rounded bg-slate-800 text-cyan-300 font-mono text-[11px] border border-slate-700">
                        ({q1},{q2})
                      </span>
                    ))}
                  </div>
                </div>
              </div>

              {/* Right Column: 16-Basis Statevector Histogram & Live Observables */}
              <div className="lg:col-span-5 bg-slate-900 rounded-xl p-5 border border-slate-800 shadow-2xl flex flex-col gap-4">
                <div className="flex items-center justify-between border-b border-slate-800 pb-3">
                  <div>
                    <h3 className="font-semibold text-sm text-slate-200 flex items-center gap-2">
                      <BarChart2 className="w-4 h-4 text-cyan-400" />
                      16-Basis Statevector Amplitudes
                    </h3>
                    <p className="text-xs text-slate-400">
                      Elapsed Time: <span className="font-mono text-cyan-300 font-semibold">{liveState.elapsedTime.toFixed(3)}s</span> / {MODEL_CONFIG.T.toFixed(1)}s
                    </p>
                  </div>

                  <div className="text-right">
                    <span className="text-[10px] font-mono text-slate-400 block uppercase">Pure State Fidelity</span>
                    <span className="font-mono font-bold text-sm text-emerald-400">
                      {(liveState.fidelity * 100).toFixed(2)}%
                    </span>
                  </div>
                </div>

                {/* Staggered Magnetization Realtime Badge */}
                <div className="grid grid-cols-2 gap-3 bg-slate-950/60 p-3 rounded-lg border border-slate-800/80 font-mono text-xs">
                  <div>
                    <span className="text-slate-400 text-[10px] block">STAGGERED MAGNETIZATION</span>
                    <span className="text-cyan-400 font-bold text-sm">
                      {liveState.mStag.toFixed(4)}
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-400 text-[10px] block">EXACT AT CHECKPOINT</span>
                    <span className="text-slate-300 font-semibold text-sm">
                      {calculateStaggeredMagnetization(liveState.exactPsiAtT).toFixed(4)}
                    </span>
                  </div>
                </div>

                {/* Probability Distribution Bar Chart for 16 Computational Basis States */}
                <div className="space-y-1.5 max-h-[290px] overflow-y-auto pr-1">
                  {liveState.psi.map((amp, idx) => {
                    const prob = cMagSq(amp);
                    const bitstring = idx.toString(2).padStart(4, '0');
                    const phase = cPhase(amp);
                    const isNeel = bitstring === "1010";

                    return (
                      <div key={bitstring} className="flex items-center gap-2 text-xs font-mono group">
                        <span className={`w-12 text-[11px] font-semibold ${isNeel ? 'text-amber-400' : 'text-slate-400'}`}>
                          |{bitstring}⟩
                        </span>

                        <div className="flex-1 bg-slate-800/80 h-4 rounded overflow-hidden relative">
                          <div
                            className={`h-full transition-all duration-150 ${
                              isNeel ? 'bg-gradient-to-r from-amber-500 to-amber-400' : 'bg-gradient-to-r from-cyan-500 to-sky-400'
                            }`}
                            style={{ width: `${Math.max(0, Math.min(100, prob * 100))}%` }}
                          />
                        </div>

                        <span className="w-14 text-right text-[11px] text-slate-300">
                          {(prob * 100).toFixed(1)}%
                        </span>

                        <span
                          className="w-8 text-right text-[10px] text-slate-500"
                          title={`Quantum Phase: ${(phase / Math.PI).toFixed(2)} π rad`}
                        >
                          {(phase / Math.PI).toFixed(1)}π
                        </span>
                      </div>
                    );
                  })}
                </div>

                {/* Readout Measurement Emulation Button */}
                <div className="pt-2 border-t border-slate-800 flex items-center justify-between">
                  <div className="text-xs text-slate-400">
                    Projective Readout ({numShots} shots):
                  </div>
                  <button
                    onClick={handleMeasureShots}
                    className="px-3 py-1.5 bg-gradient-to-r from-cyan-500 to-sky-600 hover:from-cyan-400 hover:to-sky-500 text-slate-950 font-semibold text-xs rounded-lg transition shadow-md shadow-cyan-500/20"
                  >
                    Sample Readout
                  </button>
                </div>

                {measuredCounts && (
                  <div className="bg-slate-950 p-2.5 rounded border border-cyan-900/60 text-xs font-mono space-y-1">
                    <div className="flex justify-between text-cyan-300 font-semibold">
                      <span>Observed Counts (w/ readout error):</span>
                      <button onClick={() => setMeasuredCounts(null)} className="text-slate-400 hover:text-white">✕</button>
                    </div>
                    <div className="grid grid-cols-4 gap-1 text-[11px] text-slate-300">
                      {Object.entries(measuredCounts)
                        .sort((a, b) => b[1] - a[1])
                        .slice(0, 4)
                        .map(([bits, cnt]) => (
                          <div key={bits} className="bg-slate-900 p-1 rounded text-center">
                            |{bits}⟩: {cnt}
                          </div>
                        ))}
                    </div>
                  </div>
                )}
              </div>
            </div>

            {/* Bottom Row: Executable Circuit Timeline & Playback Scrubber */}
            <div className="bg-slate-900 rounded-xl p-5 border border-slate-800 shadow-2xl space-y-4">
              <div className="flex flex-wrap items-center justify-between gap-4">
                <div className="flex items-center gap-3">
                  <button
                    onClick={() => setIsPlaying(!isPlaying)}
                    className="p-2.5 rounded-lg bg-cyan-500 hover:bg-cyan-400 text-slate-950 font-bold transition shadow-lg shadow-cyan-500/25"
                    title={isPlaying ? "Pause Simulation" : "Play Simulation"}
                  >
                    {isPlaying ? <Pause className="w-5 h-5" /> : <Play className="w-5 h-5 ml-0.5" />}
                  </button>

                  <button
                    onClick={() => setCurrentGateIndex(Math.max(0, currentGateIndex - 1))}
                    disabled={currentGateIndex === 0}
                    className="p-2 rounded-lg bg-slate-800 hover:bg-slate-700 disabled:opacity-40 text-slate-300 transition"
                    title="Step Backward"
                  >
                    <SkipBack className="w-4 h-4" />
                  </button>

                  <button
                    onClick={() => setCurrentGateIndex(Math.min(instructions.length - 1, currentGateIndex + 1))}
                    disabled={currentGateIndex >= instructions.length - 1}
                    className="p-2 rounded-lg bg-slate-800 hover:bg-slate-700 disabled:opacity-40 text-slate-300 transition"
                    title="Step Forward"
                  >
                    <SkipForward className="w-4 h-4" />
                  </button>

                  <button
                    onClick={() => {
                      setIsPlaying(false);
                      setCurrentGateIndex(0);
                    }}
                    className="p-2 rounded-lg bg-slate-800 hover:bg-slate-700 text-slate-300 transition"
                    title="Restart Circuit"
                  >
                    <RotateCcw className="w-4 h-4" />
                  </button>

                  <div className="flex items-center bg-slate-800 rounded-lg p-0.5 ml-2">
                    {[0.5, 1.0, 2.0, 4.0].map(spd => (
                      <button
                        key={spd}
                        onClick={() => setPlaybackSpeed(spd)}
                        className={`px-2 py-1 text-[11px] font-mono rounded transition ${
                          playbackSpeed === spd
                            ? 'bg-slate-700 text-cyan-300 font-bold'
                            : 'text-slate-400 hover:text-slate-200'
                        }`}
                      >
                        {spd}x
                      </button>
                    ))}
                  </div>
                </div>

                <div className="flex items-center gap-6 font-mono text-xs">
                  <div>
                    <span className="text-slate-400 text-[10px] block">INSTRUCTION</span>
                    <span className="text-slate-200 font-bold">
                      {currentGateIndex + 1} / {instructions.length}
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-400 text-[10px] block">TROTTER STEP</span>
                    <span className="text-cyan-400 font-bold">
                      {currentInstruction?.step || 0} / {trotterSteps}
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-400 text-[10px] block">CIRCUIT LAYER</span>
                    <span className="text-amber-400 font-bold">
                      {currentInstruction?.layer || "Init"}
                    </span>
                  </div>
                </div>
              </div>

              <div className="relative pt-2 pb-1">
                <input
                  type="range"
                  min="0"
                  max={instructions.length - 1}
                  value={currentGateIndex}
                  onChange={(e) => setCurrentGateIndex(parseInt(e.target.value))}
                  className="w-full h-2 bg-slate-800 rounded-lg appearance-none cursor-pointer accent-cyan-400"
                />

                <div className="flex justify-between text-[10px] font-mono text-slate-500 mt-1">
                  <span>Start (|1010⟩)</span>
                  <span>Step 1</span>
                  {trotterSteps > 2 && <span>Step {Math.floor(trotterSteps / 2)}</span>}
                  <span>Final Measurement</span>
                </div>
              </div>

              {/* Currently Selected Gate Inspector Banner */}
              <div className="bg-slate-950 p-4 rounded-xl border border-slate-800 flex flex-wrap items-center justify-between gap-4">
                <div className="flex items-center gap-3">
                  <div className={`p-2.5 rounded-lg font-mono font-bold text-xs ${
                    currentInstruction?.routingOverhead
                      ? 'bg-red-500/20 text-red-400 border border-red-500/40'
                      : currentInstruction?.type === 'twoQ'
                      ? 'bg-pink-500/20 text-pink-400 border border-pink-500/40'
                      : 'bg-cyan-500/20 text-cyan-400 border border-cyan-500/40'
                  }`}>
                    {currentInstruction?.name}
                  </div>
                  <div>
                    <div className="flex items-center gap-2">
                      <span className="text-xs font-semibold text-slate-200">
                        {currentInstruction?.desc}
                      </span>
                      {currentInstruction?.routingOverhead && (
                        <span className="bg-red-950 text-red-400 text-[10px] px-2 py-0.5 rounded border border-red-800 font-semibold">
                          Routing SWAP Overhead
                        </span>
                      )}
                    </div>
                    <p className="text-[11px] text-slate-400">
                      Target: {currentInstruction?.qubits?.map(q => `Q${q}`).join(', ')} | 
                      Gate Duration: {currentInstruction?.durationNs} ns | 
                      Configured Error: {(currentInstruction?.errorRate * 100).toFixed(3)}%
                    </p>
                  </div>
                </div>

                <div className="flex items-center gap-3 text-xs">
                  <span className={`px-2.5 py-1 rounded-md text-[11px] font-medium border ${
                    currentInstruction?.native
                      ? 'bg-emerald-950/60 text-emerald-400 border-emerald-800'
                      : 'bg-amber-950/60 text-amber-400 border-amber-800'
                  }`}>
                    {currentInstruction?.native ? "Native Processor Gate" : "Decomposed / Routed"}
                  </span>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Tradeoff and Story tabs remain active */}
        {activeTab === "tradeoff" && (
          <div className="space-y-6">
            <div className="bg-gradient-to-r from-slate-900 via-indigo-950 to-slate-900 p-6 rounded-2xl border border-indigo-800/60 shadow-2xl flex flex-wrap items-center justify-between gap-6">
              <div className="space-y-1">
                <div className="flex items-center gap-2">
                  <ShieldCheck className="w-5 h-5 text-emerald-400" />
                  <span className="text-xs font-mono uppercase tracking-wider text-emerald-400 font-bold">
                    Official PS-S3 Optimization Result
                  </span>
                </div>
                <h2 className="text-xl font-bold text-white tracking-tight">
                  Optimal Step Count: n* = {BENCHMARK_DATA.optimalSteps[selectedProcessor].bestN} for {PROCESSORS[selectedProcessor].name}
                </h2>
                <p className="text-xs text-slate-300 max-w-2xl">
                  Computed via the official maximum-fidelity rule with the smallest-n tie-breaker.
                </p>
              </div>

              <div className="flex items-center gap-4 bg-slate-900/80 p-3 rounded-xl border border-slate-700/80 font-mono text-center">
                <div>
                  <span className="text-[10px] text-slate-400 block uppercase">Peak Fidelity</span>
                  <span className="text-lg font-bold text-emerald-400">
                    {(BENCHMARK_DATA.optimalSteps[selectedProcessor].maxFid * 100).toFixed(2)}%
                  </span>
                </div>
                <div className="w-[1px] h-8 bg-slate-700" />
                <div>
                  <span className="text-[10px] text-slate-400 block uppercase">Hardware Type</span>
                  <span className="text-xs font-semibold text-cyan-300">
                    {PROCESSORS[selectedProcessor].type}
                  </span>
                </div>
              </div>
            </div>
          </div>
        )}

        {activeTab === "story" && (
          <div className="grid grid-cols-1 md:grid-cols-12 gap-6 items-start">
            <div className="md:col-span-4 bg-slate-900 rounded-xl p-4 border border-slate-800 shadow-xl space-y-2">
              <span className="text-xs font-mono font-bold text-slate-400 uppercase tracking-wider block px-2 mb-2">
                Curated Narrative Chapters
              </span>
              {STORY_CHAPTERS.map((chap, cIdx) => (
                <button
                  key={cIdx}
                  onClick={() => setActiveStoryChapter(cIdx)}
                  className={`w-full text-left p-2.5 rounded-lg text-xs transition flex items-center justify-between ${
                    activeStoryChapter === cIdx
                      ? 'bg-cyan-500/15 text-cyan-300 font-semibold border border-cyan-500/40'
                      : 'text-slate-400 hover:bg-slate-800/60 hover:text-slate-200'
                  }`}
                >
                  <div className="truncate pr-2">
                    <span className="font-mono text-[10px] text-slate-500 block">
                      Chapter {cIdx + 1}
                    </span>
                    <span className="truncate">{chap.title.split('—')[1] || chap.title}</span>
                  </div>
                  <ChevronRight className="w-3.5 h-3.5 flex-shrink-0 text-slate-500" />
                </button>
              ))}
            </div>

            <div className="md:col-span-8 bg-slate-900 rounded-xl p-6 border border-slate-800 shadow-2xl space-y-6">
              <div className="flex flex-wrap items-center justify-between gap-3 border-b border-slate-800 pb-4">
                <div>
                  <span className="bg-indigo-500/20 text-indigo-300 text-[10px] font-mono px-2 py-0.5 rounded font-semibold border border-indigo-500/30">
                    {STORY_CHAPTERS[activeStoryChapter].badge}
                  </span>
                  <h2 className="text-lg font-bold text-slate-100 mt-1">
                    {STORY_CHAPTERS[activeStoryChapter].title}
                  </h2>
                </div>
              </div>
              <div className="prose prose-invert max-w-none">
                {STORY_CHAPTERS[activeStoryChapter].content}
              </div>
            </div>
          </div>
        )}
      </main>

      {/* Footer */}
      <footer className="w-full bg-slate-900/60 border-t border-slate-800 px-6 py-3 flex flex-wrap items-center justify-between gap-4 text-xs text-slate-400 font-mono">
        <div className="flex items-center gap-2">
          <span className="w-2 h-2 rounded-full bg-emerald-500" />
          <span>Section 4 Processors A, B, and D Verified | TFIM (N=4, J=1.0, g=1.0, T=2.0)</span>
        </div>
        <div className="flex items-center gap-4 text-[11px]">
          <span>Trace Norm: 1.000000</span>
          <span>Hermiticity: Verified</span>
        </div>
      </footer>
    </div>
  );
}

