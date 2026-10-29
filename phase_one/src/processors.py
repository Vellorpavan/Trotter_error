"""
src/processors.py — StepWise-Q
Processor registry, noise model builder, graph metrics.

All 5 processor topologies are hard-coded exactly as in README section 3.
Edge lists are bidirectional (both directions stored in coupling_map).
"""

from __future__ import annotations

import json
import pathlib
from typing import Any

import networkx as nx
import numpy as np

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
P_1Q: float = 0.001
P_CX: float = 0.005
P_READOUT: float = 0.01
BASIS_GATES: list[str] = ["rz", "sx", "x", "cx"]
OPT_LEVEL: int = 1
SEED_TRANSPILER: int = 42
SHOTS: int = 4096

PROCESSORS_DIR: pathlib.Path = pathlib.Path(__file__).parent.parent / "processors"


def _bidir(edges: list[tuple[int, int]]) -> list[list[int]]:
    """Convert undirected edge list to bidirectional coupling_map."""
    result = []
    for u, v in edges:
        result.append([u, v])
        result.append([v, u])
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Processor registry
# ──────────────────────────────────────────────────────────────────────────────
PROCESSOR_REGISTRY: dict[str, dict[str, Any]] = {
    "A0": {
        "name": "A0",
        "n_qubits": 5,
        "coupling_map": _bidir([(0, 1), (1, 2), (2, 3), (3, 4)]),
        "basis_gates": BASIS_GATES,
        "noise": {
            "p_1q": P_1Q,
            "p_cx": P_CX,
            "p_readout": P_READOUT,
            "noisy_gates_1q": ["sx", "x"],
            "rz_noisy": False,
        },
        "transpile": {
            "optimization_level": OPT_LEVEL,
            "seed_transpiler": SEED_TRANSPILER,
            "shots": SHOTS,
        },
        "graph_metrics": {},
        "design_rationale": "Reference linear chain from the competition guide.",
    },
    "A1": {
        "name": "A1",
        "n_qubits": 5,
        "coupling_map": _bidir([(0, 1), (1, 2), (1, 3), (3, 4)]),
        "basis_gates": BASIS_GATES,
        "noise": {
            "p_1q": P_1Q,
            "p_cx": P_CX,
            "p_readout": P_READOUT,
            "noisy_gates_1q": ["sx", "x"],
            "rz_noisy": False,
        },
        "transpile": {
            "optimization_level": OPT_LEVEL,
            "seed_transpiler": SEED_TRANSPILER,
            "shots": SHOTS,
        },
        "graph_metrics": {},
        "design_rationale": "Reference tree topology from the competition guide.",
    },
    "A_star": {
        "name": "A_star",
        "n_qubits": 5,
        "coupling_map": _bidir([(0, 1), (1, 2), (2, 3), (3, 0), (3, 4)]),
        "basis_gates": BASIS_GATES,
        "noise": {
            "p_1q": P_1Q,
            "p_cx": P_CX,
            "p_readout": P_READOUT,
            "noisy_gates_1q": ["sx", "x"],
            "rz_noisy": False,
        },
        "transpile": {
            "optimization_level": OPT_LEVEL,
            "seed_transpiler": SEED_TRANSPILER,
            "shots": SHOTS,
        },
        "graph_metrics": {},
        "design_rationale": (
            "Designed: 4-cycle (qubits 0-1-2-3) plus tail (3-4). "
            "Natively supports the ring Hamiltonian on the 4-cycle without SWAPs. "
            "Max degree 3 (qubit 3 has degree 3). 5 edges."
        ),
    },
    "B0": {
        "name": "B0",
        "n_qubits": 9,
        "coupling_map": _bidir([
            (0, 1), (1, 2), (1, 6), (2, 3), (2, 6), (2, 5),
            (3, 4), (4, 5), (5, 6), (6, 7), (6, 8), (8, 0)
        ]),
        "basis_gates": BASIS_GATES,
        "noise": {
            "p_1q": P_1Q,
            "p_cx": P_CX,
            "p_readout": P_READOUT,
            "noisy_gates_1q": ["sx", "x"],
            "rz_noisy": False,
        },
        "transpile": {
            "optimization_level": OPT_LEVEL,
            "seed_transpiler": SEED_TRANSPILER,
            "shots": SHOTS,
        },
        "graph_metrics": {},
        "design_rationale": (
            "Reference heavy-hex-inspired topology from the competition guide. "
            "Qubit 6 has degree 5 (exempt from max-degree rule as a reference graph)."
        ),
    },
    "B_star": {
        "name": "B_star",
        "n_qubits": 9,
        "coupling_map": _bidir([
            (0, 1), (1, 2), (2, 3), (3, 0),   # 4-cycle A
            (4, 5), (5, 6), (6, 7), (7, 4),   # 4-cycle B
            (3, 4),                             # bridge
            (6, 8),                             # pendant
        ]),
        "basis_gates": BASIS_GATES,
        "noise": {
            "p_1q": P_1Q,
            "p_cx": P_CX,
            "p_readout": P_READOUT,
            "noisy_gates_1q": ["sx", "x"],
            "rz_noisy": False,
        },
        "transpile": {
            "optimization_level": OPT_LEVEL,
            "seed_transpiler": SEED_TRANSPILER,
            "shots": SHOTS,
        },
        "graph_metrics": {},
        "design_rationale": (
            "Designed: two 4-cycles connected by a bridge, plus one pendant qubit. "
            "10 edges. Max degree 3. Both 4-cycles can host the ring Hamiltonian. "
            "Bridge enables flexible layout for 4-spin chains."
        ),
    },
}


# ──────────────────────────────────────────────────────────────────────────────
# Graph metrics
# ──────────────────────────────────────────────────────────────────────────────
def graph_metrics(coupling_map: list[list[int]], n_qubits: int) -> dict[str, Any]:
    """Compute graph-theoretic metrics for a processor coupling map.

    Parameters
    ----------
    coupling_map : list of [u, v] pairs (bidirectional)
    n_qubits     : total number of qubits

    Returns
    -------
    dict with keys:
        n_qubits, n_edges, diameter, avg_degree, max_degree,
        clustering_avg, is_connected, is_bipartite
    """
    G = nx.Graph()
    G.add_nodes_from(range(n_qubits))
    for u, v in coupling_map:
        G.add_edge(u, v)

    degrees = [d for _, d in G.degree()]
    n_edges = G.number_of_edges()

    try:
        diam = nx.diameter(G)
    except nx.NetworkXError:
        diam = -1  # disconnected

    return {
        "n_qubits": n_qubits,
        "n_edges": n_edges,
        "diameter": diam,
        "avg_degree": float(np.mean(degrees)),
        "max_degree": int(np.max(degrees)),
        "clustering_avg": float(nx.average_clustering(G)),
        "is_connected": nx.is_connected(G),
        "is_bipartite": nx.is_bipartite(G),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Noise model builder
# ──────────────────────────────────────────────────────────────────────────────
def build_noise_model(proc: dict, rz_noisy: bool = False):
    """Build a Qiskit Aer NoiseModel from a processor specification.

    Parameters
    ----------
    proc     : processor dict (from PROCESSOR_REGISTRY or load_processor)
    rz_noisy : bool — if True, adds depolarizing noise to rz gates too

    Returns
    -------
    qiskit_aer.noise.NoiseModel
    """
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError

    nm = proc["noise"]
    p_1q = nm["p_1q"]
    p_cx = nm["p_cx"]
    p_ro = nm["p_readout"]
    gates_1q = nm["noisy_gates_1q"]  # ["sx", "x"]

    noise_model = NoiseModel()

    # 1-qubit depolarizing noise
    err_1q = depolarizing_error(p_1q, 1)
    for gate in gates_1q:
        noise_model.add_all_qubit_quantum_error(err_1q, gate)

    if rz_noisy:
        noise_model.add_all_qubit_quantum_error(err_1q, "rz")

    # 2-qubit CX depolarizing noise
    err_cx = depolarizing_error(p_cx, 2)
    noise_model.add_all_qubit_quantum_error(err_cx, "cx")

    # Symmetric readout error
    # ReadoutError([[p(0|0), p(1|0)], [p(0|1), p(1|1)]])
    # = [[1-p_ro, p_ro], [p_ro, 1-p_ro]]
    ro_err = ReadoutError([[1 - p_ro, p_ro], [p_ro, 1 - p_ro]])
    for q in range(proc["n_qubits"]):
        noise_model.add_readout_error(ro_err, [q])

    return noise_model


# ──────────────────────────────────────────────────────────────────────────────
# Registry access and JSON I/O
# ──────────────────────────────────────────────────────────────────────────────
def load_processor(path_or_id: str) -> dict[str, Any]:
    """Load a processor by registry ID or JSON file path.

    Parameters
    ----------
    path_or_id : str
        A key in PROCESSOR_REGISTRY (e.g. "A_star") or a path to a JSON file.

    Returns
    -------
    dict with the processor specification.
    """
    if path_or_id in PROCESSOR_REGISTRY:
        return PROCESSOR_REGISTRY[path_or_id]
    p = pathlib.Path(path_or_id)
    if not p.exists():
        raise FileNotFoundError(f"Processor not found: {path_or_id!r}")
    with open(p) as f:
        return json.load(f)


def save_processor_json(proc_id: str, path: pathlib.Path | str | None = None) -> pathlib.Path:
    """Save a processor from the registry to a JSON file.

    Parameters
    ----------
    proc_id : str — registry key
    path    : output path; defaults to processors/{proc_id}.json

    Returns
    -------
    pathlib.Path — where the file was written
    """
    proc = PROCESSOR_REGISTRY[proc_id]
    # Compute metrics before saving
    metrics = graph_metrics(proc["coupling_map"], proc["n_qubits"])
    proc = {**proc, "graph_metrics": metrics}

    if path is None:
        PROCESSORS_DIR.mkdir(parents=True, exist_ok=True)
        path = PROCESSORS_DIR / f"{proc_id}.json"

    with open(path, "w") as f:
        json.dump(proc, f, indent=2)

    return pathlib.Path(path)


def populate_all_graph_metrics() -> None:
    """Fill graph_metrics for every entry in PROCESSOR_REGISTRY (in-place)."""
    for pid, proc in PROCESSOR_REGISTRY.items():
        proc["graph_metrics"] = graph_metrics(proc["coupling_map"], proc["n_qubits"])


# ──────────────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────────────
def validate_connectivity() -> dict[str, bool]:
    """Check that all processors are connected and have the expected edge counts.

    Expected edge counts (undirected):
        A0: 4, A1: 4, A*: 5, B0: 12, B*: 10
    """
    expected_edges = {"A0": 4, "A1": 4, "A_star": 5, "B0": 12, "B_star": 10}
    results = {}
    for pid, proc in PROCESSOR_REGISTRY.items():
        metrics = graph_metrics(proc["coupling_map"], proc["n_qubits"])
        connected = metrics["is_connected"]
        edge_ok = metrics["n_edges"] == expected_edges.get(pid, metrics["n_edges"])
        results[pid] = connected and edge_ok
    return results


def validate_json_roundtrip() -> bool:
    """Write and re-read processor_A.json and processor_B.json; check equality."""
    import tempfile

    for pid in ("A_star", "B_star"):
        proc = PROCESSOR_REGISTRY[pid]
        with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as f:
            json.dump(proc, f)
            fname = f.name
        with open(fname) as f:
            loaded = json.load(f)
        # Compare coupling_map (lists may serialize differently)
        if loaded["coupling_map"] != proc["coupling_map"]:
            return False
        if loaded["n_qubits"] != proc["n_qubits"]:
            return False
    return True


if __name__ == "__main__":
    populate_all_graph_metrics()

    print("=== Connectivity checks ===")
    for pid, ok in validate_connectivity().items():
        print(f"  {pid}: {'PASS' if ok else 'FAIL'}")

    print("\n=== Graph metrics ===")
    for pid, proc in PROCESSOR_REGISTRY.items():
        m = proc["graph_metrics"]
        print(f"  {pid}: {m}")

    print("\n=== JSON round-trip ===")
    print(f"  {'PASS' if validate_json_roundtrip() else 'FAIL'}")

    print("\n=== Saving JSON files ===")
    # Write the initial choices: A* and B* as final processors
    p_a = save_processor_json("A_star", PROCESSORS_DIR / "processor_A.json")
    p_b = save_processor_json("B_star", PROCESSORS_DIR / "processor_B.json")
    print(f"  Wrote {p_a}")
    print(f"  Wrote {p_b}")
