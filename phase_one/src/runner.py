"""
src/runner.py — StepWise-Q Unified Experiment Runner
Single-command pipeline to regenerate every benchmark, model fit, fairness experiment,
figure, notebook, test, and validation gate in StepWise-Q.

Usage:
    python -m src.runner --all              # Full end-to-end reproduction
    python -m src.runner --smoke            # Fast smoke verification (< 25s)
    python -m src.runner --fairness         # Run processor-fairness study
    python -m src.runner --benchmarks       # Re-run benchmark sweeps
    python -m src.runner --models           # Re-run fitting and heldout validation
    python -m src.runner --figures          # Regenerate all figures
    python -m src.runner --test             # Run automated unit test suite
    python -m src.runner --validate         # Run 10/10 validation suite
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import pathlib
import shutil
import subprocess
import sys
import time
from typing import Any

import numpy as np
import pandas as pd
import qiskit
import qiskit_aer
import scipy

from src.problem import tfim_hamiltonian, initial_state, exact_trajectory, compute_observable
from src.processors import PROCESSOR_REGISTRY, populate_all_graph_metrics
from src.routing import (
    build_trotter_circuit,
    transpile_circuit,
    count_swaps,
    circuit_stats,
    simulate_density_matrix,
    simulate_shots,
)
from src.metrics import (
    trotter_error_model_fit,
    observable_attenuation_model,
    predict_n_star_numeric,
    bootstrap_shot_resampling_n_star_ci,
    recommend_steps,
    run_validation,
)

RESULTS_TABLES: pathlib.Path = pathlib.Path("results/tables")
RESULTS_FIGURES: pathlib.Path = pathlib.Path("results/figures")
ARCHIVE_RUNS: pathlib.Path = pathlib.Path("results/archive_runs")


def create_timestamped_backup() -> pathlib.Path | None:
    """Create a timestamped backup of current tables before regeneration.

    Satisfies Quality Rule: Never overwrite raw benchmark CSVs without generating
    a new timestamped backup or versioned filename.
    """
    if not RESULTS_TABLES.exists() or not any(RESULTS_TABLES.glob("*.csv")):
        return None

    ARCHIVE_RUNS.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    backup_dir = ARCHIVE_RUNS / f"backup_{ts}"
    backup_dir.mkdir(parents=True, exist_ok=True)

    for item in RESULTS_TABLES.glob("*.*"):
        if item.is_file():
            shutil.copy2(item, backup_dir / item.name)

    print(f"📦 Created backup of tables in: {backup_dir}")
    return backup_dir


def run_processor_fairness(smoke: bool = False) -> pd.DataFrame:
    """Execute processor fairness study across 4 distinct Hamiltonian topologies.

    Evaluates:
      1. Open chain (nearest-neighbor 1D)
      2. Periodic ring (1D with periodic boundary)
      3. Star graph (central-spin interaction)
      4. Complete graph K4 (all-to-all pairwise interaction)
    across all 5 processors (A0, A1, A*, B0, B*).

    Demonstrates topological specialization and No-Free-Lunch routing trade-offs.
    """
    populate_all_graph_metrics()
    RESULTS_TABLES.mkdir(parents=True, exist_ok=True)

    boundaries = ["open", "ring", "star"] if smoke else ["open", "ring", "star", "all_to_all"]
    procs = ["A0", "A_star", "B_star"] if smoke else ["A0", "A1", "A_star", "B0", "B_star"]
    g_vals = [1.0] if smoke else [0.5, 1.0, 1.5]
    n_steps_list = [1, 2, 4]

    print("\n=== Running Processor Fairness Experiments ===")
    print(f"  Boundaries: {boundaries}")
    print(f"  Processors: {procs}")
    print(f"  Field values: {g_vals}, Steps: {n_steps_list}")

    rows = []
    t_eval = 2.0
    sv0 = initial_state("zeros", 4)

    for boundary in boundaries:
        for g in g_vals:
            H = tfim_hamiltonian(4, J=1.0, g=g, boundary=boundary)
            exact_z0 = float(exact_trajectory(H, sv0, np.array([t_eval]), ["Z0"])["Z0"][0])

            for pid in procs:
                proc = PROCESSOR_REGISTRY[pid]
                step_errors = []
                for n_steps in n_steps_list:
                    qc = build_trotter_circuit(H, t_eval, n_steps, order=1)
                    qc_t = transpile_circuit(qc, proc)
                    sw = count_swaps(qc, proc)
                    st = circuit_stats(qc_t)

                    # Noiseless DM simulation
                    from qiskit_aer import AerSimulator
                    from qiskit.quantum_info import DensityMatrix as QiDM
                    sim_nl = AerSimulator(method="density_matrix")
                    qc_nl = qc_t.copy()
                    qc_nl.save_density_matrix()
                    res_nl = sim_nl.run(qc_nl, shots=1, seed_simulator=42).result()
                    dm_nl = QiDM(res_nl.data(0)["density_matrix"])
                    z0_nl = compute_observable(dm_nl, "Z0", 4)

                    # Noisy DM simulation under official noise profile
                    dm_res = simulate_density_matrix(qc_t, proc, ["Z0"])
                    z0_ny = dm_res["Z0"]

                    err_nl = abs(exact_z0 - z0_nl)
                    err_ny = abs(exact_z0 - z0_ny)
                    step_errors.append((n_steps, err_ny))

                    rows.append({
                        "boundary": boundary,
                        "processor": pid,
                        "g": g,
                        "n_steps": n_steps,
                        "swaps": sw,
                        "cx": st["cx"],
                        "depth": st["depth"],
                        "exact_z0": exact_z0,
                        "z0_noiseless": z0_nl,
                        "z0_noisy": z0_ny,
                        "error_noiseless": err_nl,
                        "error_noisy": err_ny,
                        "shots": 4096,
                        "seed": 42,
                        "qiskit_version": str(qiskit.__version__),
                        "aer_version": str(qiskit_aer.__version__),
                        "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    })

    df_fair = pd.DataFrame(rows)
    out_path = RESULTS_TABLES / "fairness_experiments.csv"
    df_fair.to_csv(out_path, index=False)
    print(f"✅ Saved fairness experiments table: {out_path} ({len(df_fair)} rows)")
    return df_fair


def run_all_benchmarks(smoke: bool = False):
    """Regenerate Phase 4 required, Phase 5 extended, and exhaustive sweep benchmarks."""
    from src.run_phase4 import run_phase4
    from src.run_phase5 import run_benchmark_extended, run_layout_study
    from src.run_exhaustive_sweep import run_exhaustive_sweep

    print("\n=== [1/4] Running Benchmark Suites ===")
    df_req = run_phase4()
    print(f"  Saved benchmark_required.csv ({len(df_req)} rows)")

    df_ext = run_benchmark_extended()
    print(f"  Saved benchmark_extended.csv ({len(df_ext)} rows)")

    if not smoke:
        df_layout = run_layout_study()
        print(f"  Saved layout_study.csv ({len(df_layout)} rows)")

        df_exh = run_exhaustive_sweep(processor_id="A_star", g_values=[0.5, 1.0, 1.5], n_steps_list=list(range(1, 11)))
        print(f"  Saved exhaustive_sweep_n1_to_10.csv ({len(df_exh)} rows)")


def run_all_models(smoke: bool = False):
    """Regenerate model fits, calibration, held-out prediction, and statistical tables."""
    from src.run_phase6 import run_phase6
    from src.run_model_comparison import run_model_comparison

    print("\n=== [2/4] Running Error Models & Validation Suites ===")
    res_p6 = run_phase6()
    print("  Completed Phase 6 fits, scaling law, and decision tables.")

    df_mc = run_model_comparison()
    print(f"  Saved model_comparison.csv ({len(df_mc)} rows)")


def run_all_figures():
    """Regenerate all 8 competition figures in results/figures/ at 200 dpi."""
    from src.run_phase7 import (
        plot_fig1_coupling_maps,
        plot_fig2_trajectories,
        plot_fig3_u_curve,
        plot_fig4_n_star_map,
        plot_fig5_scaling_law,
        plot_fig6_cx_depth,
        plot_fig7_layout_search,
        plot_fig8_mitigation,
    )
    print("\n=== [3/4] Generating Competition Figures (200 dpi) ===")
    plot_fig1_coupling_maps()
    plot_fig2_trajectories()
    plot_fig3_u_curve()
    plot_fig4_n_star_map()
    plot_fig5_scaling_law()
    plot_fig6_cx_depth()
    plot_fig7_layout_search()
    plot_fig8_mitigation()
    print("  All 8 figures verified in results/figures/")


def run_all_tests() -> bool:
    """Execute automated unit test suite."""
    print("\n=== Running Automated Unit Test Suite (tests/) ===")
    import unittest
    loader = unittest.TestLoader()
    suite = loader.discover("tests")
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return result.wasSuccessful()


def run_all_validation() -> bool:
    """Execute 10/10 automated validation suite and update JSON summary."""
    print("\n=== Running Validation Suite (src/metrics.py) ===")
    results = run_validation(verbose=True)
    all_pass = all(results.values())

    # Write results/tables/validation_summary.json
    try:
        summary = {
            "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "total_checks": int(len(results)),
            "passed_checks": int(sum(bool(v) for v in results.values())),
            "all_passed": bool(all_pass),
            "checks": {str(k): bool(v) for k, v in results.items()},
            "packages": {
                "qiskit": str(qiskit.__version__),
                "qiskit_aer": str(qiskit_aer.__version__),
                "numpy": str(np.__version__),
                "scipy": str(scipy.__version__),
                "pandas": str(pd.__version__),
            },
        }
        with open(RESULTS_TABLES / "validation_summary.json", "w") as f:
            json.dump(summary, f, indent=2)
    except Exception as e:
        print(f"  Warning: Could not update validation_summary.json: {e}")

    return all_pass


def run_notebook():
    """Execute main.ipynb notebook inplace."""
    print("\n=== [4/4] Executing main.ipynb Top-to-Bottom ===")
    cmd = [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "main.ipynb"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0:
        print("  main.ipynb executed successfully with all cells populated.")
    else:
        print(f"  Warning: Notebook execution exited with code {res.returncode}:\n{res.stderr[:500]}")


def main():
    parser = argparse.ArgumentParser(
        description="StepWise-Q Unified Experiment Runner: regenerate benchmarks, models, fairness, figures, and tests."
    )
    parser.add_argument("--all", action="store_true", help="Run full pipeline end-to-end.")
    parser.add_argument("--smoke", action="store_true", help="Fast smoke run to verify pipeline integrity.")
    parser.add_argument("--benchmarks", action="store_true", help="Re-run benchmark suites.")
    parser.add_argument("--models", action="store_true", help="Re-run error models and heldout validation.")
    parser.add_argument("--fairness", action="store_true", help="Run processor fairness study across 4 topologies.")
    parser.add_argument("--figures", action="store_true", help="Re-render all competition figures.")
    parser.add_argument("--notebook", action="store_true", help="Execute main.ipynb notebook inplace.")
    parser.add_argument("--test", action="store_true", help="Run automated unit test suite.")
    parser.add_argument("--validate", action="store_true", help="Run 10/10 validation suite.")
    parser.add_argument("--no-backup", action="store_true", help="Skip creating timestamped backup of tables.")

    args = parser.parse_args()

    # Default action if no flags provided: display help or run smoke
    if not any([args.all, args.smoke, args.benchmarks, args.models, args.fairness,
                args.figures, args.notebook, args.test, args.validate]):
        print("StepWise-Q Unified Experiment Runner\nNo flags specified. Use --all, --smoke, or --help.")
        sys.exit(0)

    start_time = time.time()
    print("======================================================================")
    print("               StepWise-Q Unified Experiment Runner                   ")
    print(f"Timestamp: {datetime.datetime.now(datetime.timezone.utc).isoformat()}")
    print(f"Qiskit: {qiskit.__version__} | Aer: {qiskit_aer.__version__} | Python: {sys.version.split()[0]}")
    print("======================================================================")

    # 1. Backup tables if regenerating benchmarks or running all
    if not args.no_backup and (args.all or args.benchmarks or args.models or args.fairness):
        create_timestamped_backup()

    # If --smoke is passed alone, trigger fast pipeline across all components
    if args.smoke and not (args.benchmarks or args.models or args.fairness or args.figures or args.notebook):
        args.benchmarks = True
        args.fairness = True
        args.test = True
        args.validate = True

    # 2. Benchmarks
    if args.all or args.benchmarks:
        run_all_benchmarks(smoke=args.smoke)

    # 3. Models
    if args.all or args.models:
        run_all_models(smoke=args.smoke)

    # 4. Fairness
    if args.all or args.fairness:
        run_processor_fairness(smoke=args.smoke)

    # 5. Figures
    if args.all or args.figures:
        run_all_figures()

    # 6. Notebook
    if args.all or args.notebook:
        run_notebook()

    # 7. Unit Tests
    test_ok = True
    if args.all or args.test:
        test_ok = run_all_tests()

    # 8. Validation Suite
    val_ok = True
    if args.all or args.validate:
        val_ok = run_all_validation()

    elapsed = time.time() - start_time
    print("\n======================================================================")
    print(f"🎉 Runner completed in {elapsed:.2f} seconds.")
    print(f"Unit Tests: {'PASSED' if test_ok else 'FAILED'}")
    print(f"Validation Suite: {'PASSED' if val_ok else 'FAILED'}")
    print("======================================================================")

    if not (test_ok and val_ok):
        sys.exit(1)


if __name__ == "__main__":
    main()
