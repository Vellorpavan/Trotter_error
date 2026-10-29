"""
src/routing.py — StepWise-Q
Trotter circuit builder, transpiler, SWAP counter, simulators, layout search.

API rules (from ARCHITECTURE.md):
- Use: transpile(), PauliEvolutionGate, LieTrotter, SuzukiTrotter, qiskit_aer.AerSimulator
- Never: qiskit.execute, BackendV1, Aer.get_backend (old path), removed modules
"""

from __future__ import annotations

import itertools
import warnings
from typing import Sequence

import numpy as np
import pandas as pd

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import LieTrotter, SuzukiTrotter
from qiskit.quantum_info import SparsePauliOp, Statevector, DensityMatrix

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
SHOTS: int = 4096
SEED: int = 42
OPT_LEVEL: int = 1
SEED_TRANSPILER: int = 42
BASIS_GATES: list[str] = ["rz", "sx", "x", "cx"]


# ──────────────────────────────────────────────────────────────────────────────
# Circuit builders
# ──────────────────────────────────────────────────────────────────────────────
def build_trotter_circuit(
    H: SparsePauliOp,
    t: float,
    n_steps: int,
    order: int = 1,
    initial_state_label: str | None = None,
) -> QuantumCircuit:
    """Build a Trotterized time-evolution circuit.

    Parameters
    ----------
    H       : SparsePauliOp — Hamiltonian
    t       : float — total evolution time
    n_steps : int   — number of Trotter repetitions
    order   : int   — 1 (Lie-Trotter) or 2 (Suzuki-Trotter, second order)
    initial_state_label : str | None
        If given, prepend state initialization (e.g. "0101" for Neel).

    Returns
    -------
    QuantumCircuit — logical (pre-transpile) Trotter circuit
    """
    n = H.num_qubits
    qc = QuantumCircuit(n)

    # Initialise state
    if initial_state_label is not None:
        try:
            sv = Statevector.from_label(initial_state_label)
            qc.initialize(sv.data, qc.qubits)
        except Exception:
            # Fallback: set X gates for |1> qubits (label rightmost = qubit 0)
            for qi, char in enumerate(reversed(initial_state_label)):
                if char == "1":
                    qc.x(qi)

    # Trotterization
    if order == 1:
        synth = LieTrotter(reps=n_steps)
    elif order == 2:
        synth = SuzukiTrotter(order=2, reps=n_steps)
    else:
        raise ValueError(f"order must be 1 or 2, got {order}")

    evo_gate = PauliEvolutionGate(H, time=t, synthesis=synth)
    qc.append(evo_gate, range(n))

    return qc


# ──────────────────────────────────────────────────────────────────────────────
# Transpilation
# ──────────────────────────────────────────────────────────────────────────────
def transpile_circuit(
    qc: QuantumCircuit,
    proc: dict,
    layout: list[int] | None = None,
    extra_basis: list[str] | None = None,
) -> QuantumCircuit:
    """Transpile a logical circuit to a processor's native gate set.

    Parameters
    ----------
    qc          : QuantumCircuit — logical circuit
    proc        : processor dict
    layout      : list of physical qubit indices for initial_layout
    extra_basis : optional extra gate names to add to basis_gates
                  (e.g. ["swap"] for SWAP counting pass)

    Returns
    -------
    QuantumCircuit — transpiled circuit
    """
    basis = list(proc.get("basis_gates", BASIS_GATES))
    if extra_basis:
        basis = list(set(basis + extra_basis))

    coupling_map = proc["coupling_map"]
    n_qubits = proc["n_qubits"]

    kwargs: dict = {
        "basis_gates": basis,
        "coupling_map": coupling_map,
        "optimization_level": OPT_LEVEL,
        "seed_transpiler": SEED_TRANSPILER,
    }
    if layout is not None:
        kwargs["initial_layout"] = layout
    elif qc.num_qubits < n_qubits:
        padded = QuantumCircuit(n_qubits)
        padded.compose(qc, qubits=list(range(qc.num_qubits)), inplace=True)
        qc = padded

    return transpile(qc, **kwargs)


# ──────────────────────────────────────────────────────────────────────────────
# SWAP counting
# ──────────────────────────────────────────────────────────────────────────────
def count_swaps(qc: QuantumCircuit, proc: dict, layout: list[int] | None = None) -> int:
    """Count SWAP gates introduced by the transpiler.

    Method: transpile twice:
      1. With "swap" added to basis_gates -> count SwapGate occurrences
      2. Without "swap" -> get final depth/CX for the main run

    Parameters
    ----------
    qc     : QuantumCircuit — logical or partially transpiled
    proc   : processor dict
    layout : optional initial_layout

    Returns
    -------
    int — number of SWAP gates inserted
    """
    qc_with_swap = transpile_circuit(qc, proc, layout=layout, extra_basis=["swap"])
    swap_count = 0
    for instr in qc_with_swap.data:
        if instr.operation.name == "swap":
            swap_count += 1
    return swap_count


def cx_identity_check(cx_logical: int, cx_final: int, swap_count: int) -> tuple[bool, str]:
    """Verify the CX identity: cx_final ≈ cx_logical + 3*swap_count.

    Returns (passed, message).
    """
    expected = cx_logical + 3 * swap_count
    ok = cx_final == expected
    msg = f"cx_final={cx_final}, expected={expected} (cx_logical={cx_logical}, swaps={swap_count})"
    if not ok:
        msg = f"CX identity off: {msg}"
    return ok, msg


def circuit_stats(qc_transpiled: QuantumCircuit) -> dict[str, int]:
    """Return depth, CX count, total gate count for a transpiled circuit."""
    cx_count = sum(1 for instr in qc_transpiled.data if instr.operation.name == "cx")
    return {
        "depth": qc_transpiled.depth(),
        "cx": cx_count,
        "total_gates": len(qc_transpiled.data),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Density-matrix simulation
# ──────────────────────────────────────────────────────────────────────────────
def simulate_density_matrix(
    qc: QuantumCircuit,
    proc: dict,
    observables: Sequence[str] = ("Z0",),
    rz_noisy: bool = False,
    n_spins: int = 4,
) -> dict[str, float]:
    """Run an exact (no shot noise) noisy simulation via density matrix.

    Readout error is applied analytically as (1 - 2*p_ro)^k for a weight-k Z observable.

    Parameters
    ----------
    qc         : QuantumCircuit — transpiled circuit (no measurements)
    proc       : processor dict
    observables: list of observable names (see problem.compute_observable)
    rz_noisy   : bool — include rz noise in the noise model
    n_spins    : int — number of physical spins (used for observable computation)

    Returns
    -------
    dict mapping observable -> float expectation value (with readout attenuation)
    """
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel
    from src.processors import build_noise_model
    from src.problem import compute_observable

    noise_model = build_noise_model(proc, rz_noisy=rz_noisy)
    sim = AerSimulator(method="density_matrix", noise_model=noise_model)

    # Save statevector / density matrix
    qc_dm = qc.copy()
    qc_dm.save_density_matrix()

    result = sim.run(qc_dm, shots=1, seed_simulator=SEED).result()
    dm_data = result.data(0)["density_matrix"]
    dm = DensityMatrix(dm_data)

    p_ro = proc["noise"]["p_readout"]

    out: dict[str, float] = {}
    for obs in observables:
        raw = compute_observable(dm, obs, n_spins=n_spins)
        # Analytic readout attenuation for Z observables: (1 - 2*p_ro)^weight
        # Weight-1 for Z0 and mean_mag (each Z_i); weight-2 for ZZ_nn
        if obs in ("Z0", "mean_mag"):
            raw *= (1 - 2 * p_ro) ** 1
        elif obs == "ZZ_nn":
            raw *= (1 - 2 * p_ro) ** 2
        out[obs] = raw

    return out


# ──────────────────────────────────────────────────────────────────────────────
# Shot simulation
# ──────────────────────────────────────────────────────────────────────────────
def simulate_shots(
    qc: QuantumCircuit,
    proc: dict,
    observables: Sequence[str] = ("Z0",),
    shots: int = SHOTS,
    rz_noisy: bool = False,
    n_spins: int = 4,
) -> dict[str, tuple[float, float]]:
    """Run a shot-based noisy simulation.

    Returns mean ± standard error per observable.
    Standard error = sqrt((1 - <Z>^2) / N) for a Z observable.

    Parameters
    ----------
    qc         : QuantumCircuit — transpiled circuit (no measurements)
    proc       : processor dict
    observables: list of observable names
    shots      : int
    rz_noisy   : bool
    n_spins    : int

    Returns
    -------
    dict mapping observable -> (mean, std_error)
    """
    from qiskit_aer import AerSimulator
    from src.processors import build_noise_model
    from src.problem import compute_observable

    noise_model = build_noise_model(proc, rz_noisy=rz_noisy)
    sim = AerSimulator(noise_model=noise_model)

    qc_meas = qc.copy()
    qc_meas.measure_all()

    result = sim.run(qc_meas, shots=shots, seed_simulator=SEED).result()
    counts = result.get_counts()

    out: dict[str, tuple[float, float]] = {}
    for obs in observables:
        # Compute <O> from bitstring counts
        total = 0.0
        for bitstring, count in counts.items():
            # Qiskit bitstring: leftmost = qubit n-1, rightmost = qubit 0
            # (same convention as Pauli strings)
            bits = bitstring.replace(" ", "")
            # Build Statevector from this basis state and measure observable
            sv_comp = Statevector.from_label(bits)
            val = compute_observable(sv_comp, obs, n_spins=n_spins)
            total += count * val

        mean = total / shots
        std_err = np.sqrt(max(0.0, (1 - mean**2) / shots))
        out[obs] = (mean, std_err)

    return out


# ──────────────────────────────────────────────────────────────────────────────
# Layout search
# ──────────────────────────────────────────────────────────────────────────────
def layout_search(
    qc: QuantumCircuit,
    proc: dict,
    n_spins: int = 4,
    simulate_extremes: bool = False,
    observables: Sequence[str] = ("Z0",),
) -> pd.DataFrame:
    """Exhaustive layout search: transpile for all injective layouts.

    For an n_spins qubit circuit on an m-qubit processor:
        Number of injective layouts = P(m, n) = m! / (m - n)!

    This function does transpile-only for all layouts, then optionally
    simulates the best, median, and worst layouts by CX count.

    Parameters
    ----------
    qc              : QuantumCircuit — logical Trotter circuit
    proc            : processor dict
    n_spins         : int — number of active qubits (4 for 4-spin TFIM)
    simulate_extremes: bool — if True, simulate best/median/worst via DM
    observables     : observable names for simulation

    Returns
    -------
    pd.DataFrame with columns:
        layout, depth, cx, swap, [z0_best, z0_median, z0_worst if simulated]
    """
    n_q = proc["n_qubits"]
    all_layouts = list(itertools.permutations(range(n_q), n_spins))

    rows = []
    for layout in all_layouts:
        try:
            qc_t = transpile_circuit(qc, proc, layout=list(layout))
            stats = circuit_stats(qc_t)
            swap = count_swaps(qc, proc, layout=list(layout))
            rows.append({
                "layout": layout,
                "depth": stats["depth"],
                "cx": stats["cx"],
                "swap": swap,
            })
        except Exception as e:
            warnings.warn(f"Layout {layout} failed: {e}")

    df = pd.DataFrame(rows).sort_values("cx").reset_index(drop=True)

    if simulate_extremes and len(df) > 0:
        n = len(df)
        indices = {
            "best": 0,
            "median": n // 2,
            "worst": n - 1,
        }
        for role, idx in indices.items():
            row = df.iloc[idx]
            layout = list(row["layout"])
            qc_t = transpile_circuit(qc, proc, layout=layout)
            dm_result = simulate_density_matrix(qc_t, proc, observables, n_spins=n_spins)
            for obs, val in dm_result.items():
                df.loc[idx, f"{obs}_{role}"] = val

    return df


# ──────────────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────────────
def validate_swap_counts(g: float = 1.0) -> dict[str, bool]:
    """Check SWAP = 0 for open chain on A0, A*, B0, B*; ring on A0 gives SWAP > 0.

    Returns dict of {check_name: bool}
    """
    from src.problem import tfim_hamiltonian
    from src.processors import PROCESSOR_REGISTRY, populate_all_graph_metrics

    populate_all_graph_metrics()

    H_open = tfim_hamiltonian(4, 1.0, g, "open")
    H_ring = tfim_hamiltonian(4, 1.0, g, "ring")
    t = 1.0
    n = 1

    results = {}

    # Open chain: SWAP should be 0 on A0, A_star, B0, B_star
    for pid in ("A0", "A_star", "B0", "B_star"):
        proc = PROCESSOR_REGISTRY[pid]
        qc = build_trotter_circuit(H_open, t, n, order=1)
        sw = count_swaps(qc, proc)
        results[f"open_swap_zero_{pid}"] = (sw == 0)

    # Ring: A0 should need SWAPs (linear chain can't natively do Z0 Z3 coupling)
    proc_A0 = PROCESSOR_REGISTRY["A0"]
    qc_ring = build_trotter_circuit(H_ring, t, n, order=1)
    sw_ring_A0 = count_swaps(qc_ring, proc_A0)
    results["ring_swap_gt0_A0"] = (sw_ring_A0 > 0)

    # A_star is a 4-cycle + tail: ring should need 0 SWAPs on A_star
    proc_Astar = PROCESSOR_REGISTRY["A_star"]
    sw_ring_Astar = count_swaps(qc_ring, proc_Astar)
    results["ring_swap_zero_Astar"] = (sw_ring_Astar == 0)

    return results


if __name__ == "__main__":
    print("=== SWAP validation ===")
    checks = validate_swap_counts(g=1.0)
    for name, ok in checks.items():
        print(f"  {name}: {'PASS' if ok else 'FAIL'}")
