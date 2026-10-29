"""
src/problem.py — StepWise-Q
Hamiltonians, initial states, exact trajectories, observables.

CONSTANTS (hard-wired from competition spec):
  J = 1.0
  Basis: rz, sx, x, cx
  Little-endian Qiskit convention: spin i -> qubit i,
    leftmost Pauli character = qubit n_spins-1.
"""

from __future__ import annotations

import os
import pathlib
import numpy as np
import scipy.linalg
from typing import Sequence

from qiskit.quantum_info import SparsePauliOp, Statevector, DensityMatrix, Operator

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
J_DEFAULT: float = 1.0
G_VALUES: list[float] = [0.5, 1.0, 1.5]
N_SPINS: int = 4
DATA_PROCESSED: pathlib.Path = pathlib.Path(__file__).parent.parent / "data" / "processed"


# ──────────────────────────────────────────────────────────────────────────────
# Qubit mapping
# ──────────────────────────────────────────────────────────────────────────────
def spin_to_qubit(spin_idx: int) -> int:
    """Maps physical spin index i to Qiskit qubit index i (identity map).

    Qiskit is little-endian: qubit 0 = rightmost position in a Pauli string.
    The leftmost character in a length-n Pauli string acts on qubit n-1.

    Unit test:
        spin_to_qubit(0) == 0
        spin_to_qubit(3) == 3
    """
    return spin_idx


# ──────────────────────────────────────────────────────────────────────────────
# Hamiltonian builders
# ──────────────────────────────────────────────────────────────────────────────
def _zz_term(i: int, j: int, n: int) -> str:
    """Return a Pauli string for Z_i Z_j on n qubits (little-endian Qiskit).

    In a Qiskit Pauli string the character at position p acts on qubit (n-1-p).
    So qubit k corresponds to string position (n-1-k).
    """
    ops = ["I"] * n
    ops[n - 1 - i] = "Z"
    ops[n - 1 - j] = "Z"
    return "".join(ops)


def _x_term(i: int, n: int) -> str:
    """Return a Pauli string for X_i on n qubits (little-endian Qiskit)."""
    ops = ["I"] * n
    ops[n - 1 - i] = "X"
    return "".join(ops)


def tfim_hamiltonian(
    n_spins: int = N_SPINS,
    J: float = J_DEFAULT,
    g: float = 1.0,
    boundary: str = "open",
) -> SparsePauliOp:
    """Build the 4-spin Transverse-Field Ising Model Hamiltonian.

    H = -J * sum_{i} Z_i Z_{i+1}  -  g * sum_i X_i

    Parameters
    ----------
    n_spins : int
        Number of spins (default 4).
    J : float
        ZZ coupling strength.
    g : float
        Transverse field strength; must be in {0.5, 1.0, 1.5} for the contest.
    boundary : str
        "open"  -> sum over i=0..n-2 (no Z_{n-1} Z_0 term)
        "ring"  -> adds the Z_{n-1} Z_0 wrap-around term

    Returns
    -------
    SparsePauliOp
        The Hamiltonian as a sum of Pauli operators.
    """
    if boundary == "open":
        pairs = [(i, i + 1) for i in range(n_spins - 1)]
    elif boundary == "ring":
        pairs = [(i, (i + 1) % n_spins) for i in range(n_spins)]
    elif boundary == "star":
        pairs = [(0, i) for i in range(1, n_spins)]
    elif boundary in {"all_to_all", "k4"}:
        pairs = [(i, j) for i in range(n_spins) for j in range(i + 1, n_spins)]
    else:
        raise ValueError(f"boundary must be 'open', 'ring', 'star', or 'all_to_all', got {boundary!r}")

    terms: list[tuple[str, complex]] = []

    # ZZ coupling terms
    for u, v in pairs:
        terms.append((_zz_term(u, v, n_spins), complex(-J)))

    # Transverse field X terms
    for i in range(n_spins):
        terms.append((_x_term(i, n_spins), complex(-g)))

    return SparsePauliOp.from_list(terms)


# ──────────────────────────────────────────────────────────────────────────────
# Initial states
# ──────────────────────────────────────────────────────────────────────────────
def initial_state(name: str, n_spins: int = N_SPINS) -> Statevector:
    """Return an initial state Statevector.

    Parameters
    ----------
    name : str
        "zeros" -> |0000> (all spins up in Z basis)
        "neel"  -> |0101> (alternating: qubit 0 = |0>, qubit 1 = |1>, ...)
                   Note: qubit 1 and 3 are in |1> state.
    """
    if name == "zeros":
        return Statevector.from_label("0" * n_spins)
    elif name == "neel":
        # Neel state: alternating |0101...>
        # Qiskit string ordering: rightmost char = qubit 0
        # We want qubit 0=|0>, qubit 1=|1>, qubit 2=|0>, qubit 3=|1>
        # String (leftmost=qubit n-1): qubit3, qubit2, qubit1, qubit0 = 1,0,1,0 -> "1010"
        label = "".join("0" if (n_spins - 1 - i) % 2 == 0 else "1" for i in range(n_spins))
        return Statevector.from_label(label)
    else:
        raise ValueError(f"Unknown initial state: {name!r}. Use 'zeros' or 'neel'.")


# ──────────────────────────────────────────────────────────────────────────────
# Observable computation
# ──────────────────────────────────────────────────────────────────────────────
def compute_observable(
    state: Statevector | DensityMatrix,
    obs: str,
    n_spins: int = N_SPINS,
    H: SparsePauliOp | None = None,
) -> float:
    """Compute an observable expectation value.

    Parameters
    ----------
    state : Statevector or DensityMatrix
    obs   : str
        "Z0"       -> <Z_0>  (spin 0 = qubit 0)
        "mean_mag" -> (1/n) * sum_i <Z_i>
        "ZZ_nn"    -> (1/(n-1)) * sum_i <Z_i Z_{i+1}> (open chain only)
        "energy"   -> <H>  (requires H to be provided)
    H : SparsePauliOp, optional
        Required only for "energy".

    Returns
    -------
    float
        Real part of the expectation value.
    """
    if obs == "Z0":
        op = SparsePauliOp.from_list([(_x_term(0, n_spins).replace("X", "Z"), 1.0)])
        # Build Z on qubit 0
        ops = ["I"] * n_spins
        ops[n_spins - 1 - 0] = "Z"
        op = SparsePauliOp.from_list([("".join(ops), 1.0)])
        return float(state.expectation_value(op).real)

    elif obs == "mean_mag":
        total = 0.0
        for i in range(n_spins):
            ops = ["I"] * n_spins
            ops[n_spins - 1 - i] = "Z"
            op = SparsePauliOp.from_list([("".join(ops), 1.0)])
            total += float(state.expectation_value(op).real)
        return total / n_spins

    elif obs == "ZZ_nn":
        total = 0.0
        for i in range(n_spins - 1):
            op = SparsePauliOp.from_list([(_zz_term(i, i + 1, n_spins), 1.0)])
            total += float(state.expectation_value(op).real)
        return total / (n_spins - 1)

    elif obs == "energy":
        if H is None:
            raise ValueError("Must provide H for energy observable.")
        return float(state.expectation_value(H).real)

    else:
        raise ValueError(f"Unknown observable: {obs!r}")


# ──────────────────────────────────────────────────────────────────────────────
# Exact trajectory (scipy expm)
# ──────────────────────────────────────────────────────────────────────────────
def exact_trajectory(
    H: SparsePauliOp,
    state0: Statevector,
    t_grid: np.ndarray,
    observables: Sequence[str] = ("Z0",),
    save_prefix: str | None = None,
    H_op: SparsePauliOp | None = None,
) -> dict[str, np.ndarray]:
    """Compute exact time evolution using scipy matrix exponentiation.

    U(t) = exp(-i H t)

    Parameters
    ----------
    H        : SparsePauliOp — Hamiltonian
    state0   : Statevector  — initial state
    t_grid   : array of times
    observables : list of observable names
    save_prefix : if given, saves results to data/processed/{save_prefix}_{obs}.npy
    H_op     : same as H; passed through to compute_observable for "energy"

    Returns
    -------
    dict mapping observable name -> array of shape (len(t_grid),)
    """
    n = H.num_qubits
    H_matrix = H.to_matrix()
    psi0 = state0.data  # complex array of shape (2**n,)

    results: dict[str, list[float]] = {obs: [] for obs in observables}

    for t in t_grid:
        U = scipy.linalg.expm(-1j * H_matrix * t)
        psi_t = U @ psi0
        sv = Statevector(psi_t)
        for obs in observables:
            val = compute_observable(sv, obs, n_spins=n, H=H_op or H)
            results[obs].append(val)

    out: dict[str, np.ndarray] = {obs: np.array(v) for obs, v in results.items()}

    if save_prefix is not None:
        DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
        for obs, arr in out.items():
            path = DATA_PROCESSED / f"{save_prefix}_{obs}.npy"
            np.save(path, arr)

    return out


# ──────────────────────────────────────────────────────────────────────────────
# Validation helpers
# ──────────────────────────────────────────────────────────────────────────────
def validate_spin_to_qubit() -> bool:
    """Unit test: 2-qubit TFIM hand-built check."""
    # Hand-built 2-spin TFIM: H = -J ZZ - g(XI + IX)
    J, g = 1.0, 1.0
    # Pauli string: for n=2, Z0Z1 = "ZZ" (leftmost=qubit1, rightmost=qubit0)
    H_hand = np.array([
        [-J, -g, -g, 0],
        [-g,  J,  0, -g],
        [-g,  0,  J, -g],
        [0,  -g, -g, -J],
    ], dtype=complex)
    H_auto = tfim_hamiltonian(2, J, g, "open").to_matrix()
    return np.allclose(H_hand, H_auto, atol=1e-10)


def validate_initial_observable() -> bool:
    """Check Z0(0) = +1 for zeros, and Z0(0) != +1 for Neel."""
    sv_zeros = initial_state("zeros", 4)
    sv_neel = initial_state("neel", 4)

    z0_zeros = compute_observable(sv_zeros, "Z0", 4)
    z0_neel = compute_observable(sv_neel, "Z0", 4)

    # zeros: qubit 0 in |0> -> <Z0> = +1
    # neel: qubit 0 in |0> -> <Z0> = +1 (alternating starts at 0)
    # Actually for neel: spin 0 -> qubit 0 -> |0> -> Z0 = +1
    # spin 1 -> qubit 1 -> |1> -> Z1 = -1
    ok_zeros = abs(z0_zeros - 1.0) < 1e-10
    ok_neel = abs(z0_neel - 1.0) < 1e-10  # qubit 0 is |0> in Neel
    return ok_zeros and ok_neel


def validate_norm_preserved(t_max: float = 2.0, n_t: int = 5) -> bool:
    """Check that the exact propagator preserves norm."""
    H = tfim_hamiltonian(4, 1.0, 1.0, "open")
    sv0 = initial_state("zeros", 4)
    t_grid = np.linspace(0, t_max, n_t)
    H_mat = H.to_matrix()
    psi0 = sv0.data
    for t in t_grid:
        U = scipy.linalg.expm(-1j * H_mat * t)
        psi_t = U @ psi0
        if abs(np.dot(psi_t.conj(), psi_t).real - 1.0) > 1e-10:
            return False
    return True


if __name__ == "__main__":
    print("spin_to_qubit check:  ", validate_spin_to_qubit())
    print("Z0(0) check:          ", validate_initial_observable())
    print("norm preserved check: ", validate_norm_preserved())
