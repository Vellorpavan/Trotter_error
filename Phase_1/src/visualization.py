"""
src/visualization.py — StepWise-Q
All 8 figures required by the competition, plus coupling-map diagrams.

All figures are saved at 200 dpi as PNG.
"""

from __future__ import annotations

import pathlib
from typing import Sequence

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import networkx as nx
import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# Defaults
# ──────────────────────────────────────────────────────────────────────────────
DPI: int = 200
FIGURES_DIR: pathlib.Path = pathlib.Path(__file__).parent.parent / "results" / "figures"
STYLE_KW = dict(fontfamily="DejaVu Sans")

PROC_COLORS = {
    "A0": "#2196F3",    # blue
    "A1": "#9C27B0",    # purple
    "A_star": "#4CAF50",  # green
    "B0": "#FF5722",    # deep orange
    "B_star": "#009688",  # teal
}

N_COLORS = {1: "#E53935", 2: "#FB8C00", 4: "#43A047", 6: "#1E88E5", 8: "#8E24AA"}


def _ensure_dir(path: str | pathlib.Path) -> pathlib.Path:
    p = pathlib.Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# ──────────────────────────────────────────────────────────────────────────────
# Figure 1: Coupling maps
# ──────────────────────────────────────────────────────────────────────────────
def fig1_coupling_maps(
    processors: list[dict],
    save_path: str | pathlib.Path | None = None,
    highlight_qubits: dict[str, list[int]] | None = None,
) -> plt.Figure:
    """Plot all processor topologies with highlighted active-qubit embeddings.

    Parameters
    ----------
    processors       : list of processor dicts
    save_path        : output path (default: results/figures/fig1_coupling_maps.png)
    highlight_qubits : dict mapping processor ID -> list of qubit indices to highlight
                       (default: first 4 qubits = the TFIM spins)
    """
    n_procs = len(processors)
    ncols = min(3, n_procs)
    nrows = (n_procs + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    fig.patch.set_facecolor("#0f1117")

    if n_procs == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)

    for idx, proc in enumerate(processors):
        ax = axes[idx // ncols][idx % ncols]
        ax.set_facecolor("#1a1d27")

        pid = proc["name"]
        n_q = proc["n_qubits"]
        coupling_map = proc["coupling_map"]

        G = nx.Graph()
        G.add_nodes_from(range(n_q))
        for u, v in coupling_map:
            G.add_edge(u, v)

        # Layout
        pos = nx.spring_layout(G, seed=42, k=1.5)

        # Active qubits (TFIM spins 0-3)
        active = (highlight_qubits or {}).get(pid, list(range(4)))
        node_colors = [
            PROC_COLORS.get(pid, "#4fc3f7") if i in active else "#546e7a"
            for i in range(n_q)
        ]
        edge_colors = []
        for u, v in G.edges():
            if u in active and v in active:
                edge_colors.append(PROC_COLORS.get(pid, "#4fc3f7"))
            else:
                edge_colors.append("#546e7a")

        nx.draw(
            G,
            pos=pos,
            ax=ax,
            node_color=node_colors,
            edge_color=edge_colors,
            node_size=500,
            width=2.5,
            with_labels=True,
            font_color="white",
            font_size=10,
            font_weight="bold",
        )

        metrics = proc.get("graph_metrics", {})
        n_edges = metrics.get("n_edges", len(G.edges()))
        diam = metrics.get("diameter", "?")
        ax.set_title(
            f"{pid}  ·  {n_q}q  ·  {n_edges} edges  ·  diam={diam}",
            color="white", fontsize=11, pad=8,
        )
        ax.set_axis_off()

    # Hide unused subplots
    for idx in range(n_procs, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(
        "Processor Coupling Maps  (highlighted = active TFIM qubits)",
        color="white", fontsize=13, y=1.01,
    )
    plt.tight_layout()

    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig1_coupling_maps.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 2: Trajectory — exact vs noisy
# ──────────────────────────────────────────────────────────────────────────────
def fig2_trajectory(
    results: dict,
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """Plot ⟨Z₀(t)⟩ exact vs noisy for n=1,2,4 and two processors.

    Parameters
    ----------
    results : dict with keys like ("A_star", 1, "Z0_exact"), ("A_star", 1, "Z0_noisy"), etc.
        Expected keys:
            t_grid           : np.ndarray
            exact            : np.ndarray  (exact trajectory)
            (proc_id, n_steps, "noisy") : np.ndarray
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    fig.patch.set_facecolor("#0f1117")

    t_grid = results["t_grid"]
    exact = results["exact"]
    proc_ids = results.get("proc_ids", ["A_star", "B_star"])
    n_steps_list = results.get("n_steps_list", [1, 2, 4])

    for ax_idx, n_steps in enumerate(n_steps_list):
        ax = axes[ax_idx]
        ax.set_facecolor("#1a1d27")

        # Exact
        ax.plot(t_grid, exact, color="white", lw=2, ls="--", label="Exact", zorder=5)

        # Noisy per processor
        for pid in proc_ids:
            key = (pid, n_steps, "noisy")
            if key in results:
                arr = results[key]
                color = PROC_COLORS.get(pid, "#4fc3f7")
                ax.plot(t_grid, arr, color=color, lw=1.8, label=f"{pid} noisy")

        ax.set_xlabel("t", color="#b0bec5")
        ax.set_title(f"n = {n_steps}", color="white", fontsize=12)
        ax.tick_params(colors="#b0bec5")
        ax.spines[:].set_color("#37474f")
        for spine in ax.spines.values():
            spine.set_edgecolor("#37474f")

        if ax_idx == 0:
            ax.set_ylabel("⟨Z₀(t)⟩", color="#b0bec5")
        if ax_idx == 2:
            ax.legend(facecolor="#1a1d27", edgecolor="#37474f",
                      labelcolor="white", fontsize=9)

    fig.suptitle(
        "⟨Z₀(t)⟩: Exact vs Noisy Trotter  (open chain, g=1.0, |0000⟩)",
        color="white", fontsize=13,
    )
    plt.tight_layout()

    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig2_trajectory.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 3: U-curve ε(n) with fitted model (hero figure 1)
# ──────────────────────────────────────────────────────────────────────────────
def fig3_u_curve(
    n_vals: np.ndarray,
    errors_noisy: np.ndarray,
    model_fit: dict,
    n_star: int,
    errors_noiseless: np.ndarray | None = None,
    std_errors: np.ndarray | None = None,
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """U-curve: total error ε(n) with fitted model and n* marked.

    This is the primary hero figure.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    n_vals = np.asarray(n_vals, dtype=float)
    errors_noisy = np.asarray(errors_noisy, dtype=float)

    # Data points with error bars
    if std_errors is not None:
        ax.errorbar(n_vals, errors_noisy, yerr=std_errors, fmt="o",
                    color="#4fc3f7", ms=8, capsize=5, lw=2, zorder=4, label="Noisy error (data)")
    else:
        ax.plot(n_vals, errors_noisy, "o-", color="#4fc3f7", ms=8, lw=2,
                zorder=4, label="Noisy error (data)")

    if errors_noiseless is not None:
        ax.plot(n_vals, np.asarray(errors_noiseless), "s--", color="#81c784", ms=7, lw=1.5,
                zorder=3, label="Noiseless Trotter error")

    # Fitted model curve
    c = model_fit["c"]
    k = model_fit["k_fit"]
    Lambda = model_fit["Lambda"]
    n_fine = np.linspace(0.5, max(n_vals) * 1.3, 300)

    def model_curve(n):
        return c / n**k + (1 - Lambda**n)

    ax.plot(n_fine, model_curve(n_fine), color="#ff9800", lw=2.5, ls="-",
            label=f"Fitted model  (c={c:.3f}, Λ={Lambda:.4f})", zorder=3)

    # n* vertical marker
    ax.axvline(n_star, color="#ef5350", lw=2, ls="--", zorder=5)
    ax.text(n_star + 0.05 * max(n_vals), ax.get_ylim()[1] * 0.05 if len(errors_noisy) else 0.5,
            f"n* = {n_star}", color="#ef5350", fontsize=12, fontweight="bold")

    ax.set_xlabel("Trotter steps n", color="#b0bec5", fontsize=12)
    ax.set_ylabel("Total error  ε(n) = |⟨Z₀⟩_ex − ⟨Z₀⟩_noisy|", color="#b0bec5", fontsize=12)
    ax.set_title(
        "U-curve: Trotter Error vs Step Count\n"
        f"(c/n^k = Trotter term; 1−Λ^n = noise term)",
        color="white", fontsize=12,
    )
    ax.tick_params(colors="#b0bec5")
    for spine in ax.spines.values():
        spine.set_edgecolor("#37474f")
    ax.legend(facecolor="#1a1d27", edgecolor="#37474f", labelcolor="white", fontsize=10)

    plt.tight_layout()
    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig3_u_curve.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 4: n* map over (t, noise scale) (hero figure 2)
# ──────────────────────────────────────────────────────────────────────────────
def fig4_n_star_map(
    n_star_grid: np.ndarray,
    t_grid: np.ndarray,
    noise_grid: np.ndarray,
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """n* map: 2D heatmap over (evolution time t, noise scale)."""
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    im = ax.imshow(
        n_star_grid,
        extent=[t_grid[0], t_grid[-1], noise_grid[0], noise_grid[-1]],
        origin="lower",
        aspect="auto",
        cmap="plasma",
        vmin=1,
    )

    # Contour lines
    try:
        ct = ax.contour(
            t_grid, noise_grid, n_star_grid,
            levels=sorted(set(n_star_grid.flatten())),
            colors="white", linewidths=0.8, alpha=0.5,
        )
        ax.clabel(ct, inline=True, fontsize=8, colors="white")
    except Exception:
        pass

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Optimal Trotter steps n*", color="white", fontsize=11)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    ax.set_xlabel("Evolution time t", color="#b0bec5", fontsize=12)
    ax.set_ylabel("Noise scale (× official rates)", color="#b0bec5", fontsize=12)
    ax.set_title("n* Map: Optimal Trotter Steps vs (t, Noise Scale)", color="white", fontsize=13)
    ax.tick_params(colors="#b0bec5")
    for spine in ax.spines.values():
        spine.set_edgecolor("#37474f")

    plt.tight_layout()
    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig4_n_star_map.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 5: Scaling law log-log
# ──────────────────────────────────────────────────────────────────────────────
def fig5_scaling_law(
    noise_scales: list[float],
    n_stars: list[int],
    fit_slope: float,
    r_squared: float = 0.0,
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """Scaling law: n* ∝ p^α on a log-log plot with fitted slope."""
    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    noise_arr = np.array(noise_scales)
    n_arr = np.array(n_stars, dtype=float)

    ax.loglog(noise_arr, n_arr, "o", color="#4fc3f7", ms=10, zorder=4, label="Measured n*")

    # Fitted line
    n_fit = noise_arr[0] ** fit_slope * np.exp(np.log(n_arr).mean() - fit_slope * np.log(noise_arr).mean())
    ax.loglog(noise_arr, n_fit * noise_arr ** fit_slope, color="#ff9800", lw=2.5,
              label=f"Fit: n* ∝ p^{fit_slope:.3f}  (R²={r_squared:.3f})")

    # Theory lines
    for k, label, color in [
        (1, "Theory order-1: slope = -0.5", "#81c784"),
        (2, "Theory order-2: slope = -0.33", "#ce93d8"),
    ]:
        slope_th = -1 / (k + 1)
        y_th = n_arr.mean() * (noise_arr / noise_arr.mean()) ** slope_th
        ax.loglog(noise_arr, y_th, ls="--", color=color, lw=1.5, alpha=0.7, label=label)

    ax.set_xlabel("Noise scale (× official)", color="#b0bec5", fontsize=12)
    ax.set_ylabel("Optimal step count n*", color="#b0bec5", fontsize=12)
    ax.set_title("Scaling Law: n* vs Noise Level", color="white", fontsize=13)
    ax.tick_params(colors="#b0bec5")
    for spine in ax.spines.values():
        spine.set_edgecolor("#37474f")
    ax.legend(facecolor="#1a1d27", edgecolor="#37474f", labelcolor="white", fontsize=9)
    ax.grid(True, which="both", color="#2a2d3e", ls="--", alpha=0.6)

    plt.tight_layout()
    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig5_scaling_law.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 6: CX count and depth vs n, open vs ring
# ──────────────────────────────────────────────────────────────────────────────
def fig6_cx_depth(
    results: dict,
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """Plot CX count and circuit depth vs Trotter steps, open vs ring Hamiltonian.

    Parameters
    ----------
    results : dict with keys:
        n_steps    : list of n values
        (proc_id, boundary, "cx")    : list of CX counts
        (proc_id, boundary, "depth") : list of depths
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor("#0f1117")

    for ax in (ax1, ax2):
        ax.set_facecolor("#1a1d27")
        ax.tick_params(colors="#b0bec5")
        for spine in ax.spines.values():
            spine.set_edgecolor("#37474f")
        ax.grid(True, color="#2a2d3e", ls="--", alpha=0.6)

    n_steps = results.get("n_steps", [1, 2, 4])
    proc_ids = results.get("proc_ids", ["A_star", "B_star"])
    boundaries = ["open", "ring"]

    ls_map = {"open": "-", "ring": "--"}
    for pid in proc_ids:
        for boundary in boundaries:
            cx_key = (pid, boundary, "cx")
            depth_key = (pid, boundary, "depth")
            if cx_key not in results:
                continue
            color = PROC_COLORS.get(pid, "#4fc3f7")
            ls = ls_map[boundary]
            label = f"{pid} {boundary}"
            ax1.plot(n_steps, results[cx_key], "o" + ls, color=color, lw=2, ms=7, label=label)
            ax2.plot(n_steps, results[depth_key], "s" + ls, color=color, lw=2, ms=7, label=label)

    ax1.set_xlabel("Trotter steps n", color="#b0bec5")
    ax1.set_ylabel("CX gate count", color="#b0bec5")
    ax1.set_title("CX Count vs n", color="white", fontsize=12)
    ax1.legend(facecolor="#1a1d27", edgecolor="#37474f", labelcolor="white", fontsize=9)

    ax2.set_xlabel("Trotter steps n", color="#b0bec5")
    ax2.set_ylabel("Circuit depth", color="#b0bec5")
    ax2.set_title("Circuit Depth vs n", color="white", fontsize=12)

    fig.suptitle("CX Count and Circuit Depth: Open vs Ring Hamiltonian", color="white", fontsize=13)
    plt.tight_layout()

    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig6_cx_depth.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 7: Layout study — SWAP overhead and box plots
# ──────────────────────────────────────────────────────────────────────────────
def fig7_layout(
    layout_df: pd.DataFrame,
    proc_name: str = "A_star",
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """SWAP overhead and CX count distribution across all layouts.

    Parameters
    ----------
    layout_df : pd.DataFrame with columns: layout, depth, cx, swap
    proc_name : str — processor name for the title
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor("#0f1117")

    for ax in (ax1, ax2):
        ax.set_facecolor("#1a1d27")
        ax.tick_params(colors="#b0bec5")
        for spine in ax.spines.values():
            spine.set_edgecolor("#37474f")

    # Box plot of CX counts
    bp = ax1.boxplot(
        layout_df["cx"].values,
        patch_artist=True,
        notch=False,
        boxprops=dict(facecolor="#1565c0", color="#4fc3f7"),
        whiskerprops=dict(color="#4fc3f7"),
        capprops=dict(color="#4fc3f7"),
        medianprops=dict(color="#ef5350", lw=2),
        flierprops=dict(marker=".", color="#4fc3f7", ms=4),
    )
    ax1.set_xticks([1])
    ax1.set_xticklabels(["All layouts"])
    ax1.set_ylabel("CX count", color="#b0bec5")
    ax1.set_title(f"CX Distribution ({proc_name})", color="white", fontsize=11)

    # Best / median / worst highlighted
    sorted_df = layout_df.sort_values("cx").reset_index(drop=True)
    n = len(sorted_df)
    special = {
        "best":   (sorted_df.iloc[0]["cx"], "#43A047"),
        "median": (sorted_df.iloc[n // 2]["cx"], "#FB8C00"),
        "worst":  (sorted_df.iloc[-1]["cx"], "#E53935"),
    }
    for label, (val, color) in special.items():
        ax1.axhline(val, color=color, lw=1.5, ls="--", label=f"{label} = {val:.0f}")
    ax1.legend(facecolor="#1a1d27", edgecolor="#37474f", labelcolor="white", fontsize=9)

    # Scatter: SWAP count vs CX count
    ax2.scatter(layout_df["swap"], layout_df["cx"], c="#4fc3f7", alpha=0.5, s=20)
    ax2.set_xlabel("SWAP count", color="#b0bec5")
    ax2.set_ylabel("CX count", color="#b0bec5")
    ax2.set_title(f"SWAP vs CX ({proc_name})", color="white", fontsize=11)
    ax2.grid(True, color="#2a2d3e", ls="--", alpha=0.6)

    fig.suptitle(f"Layout Study — {proc_name}", color="white", fontsize=13)
    plt.tight_layout()

    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig7_layout.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Figure 8: Mitigation and extrapolation — before/after bars
# ──────────────────────────────────────────────────────────────────────────────
def fig8_mitigation(
    before: dict[str, tuple[float, float]],
    after: dict[str, tuple[float, float]],
    exact_val: float | None = None,
    save_path: str | pathlib.Path | None = None,
) -> plt.Figure:
    """Before/after bars for ZNE and readout mitigation.

    Parameters
    ----------
    before : dict mapping method_name -> (mean, std_err)
    after  : dict mapping method_name -> (mean, std_err)
    exact_val : float | None — exact reference value
    """
    fig, ax = plt.subplots(figsize=(9, 5))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#1a1d27")

    methods = sorted(set(list(before.keys()) + list(after.keys())))
    x = np.arange(len(methods))
    width = 0.35

    before_means = [before.get(m, (0, 0))[0] for m in methods]
    before_errs = [before.get(m, (0, 0))[1] for m in methods]
    after_means = [after.get(m, (0, 0))[0] for m in methods]
    after_errs = [after.get(m, (0, 0))[1] for m in methods]

    bars1 = ax.bar(x - width / 2, before_means, width, yerr=before_errs,
                   color="#1565c0", capsize=5, label="Before mitigation",
                   edgecolor="#4fc3f7", linewidth=1.2, error_kw=dict(ecolor="#4fc3f7"))
    bars2 = ax.bar(x + width / 2, after_means, width, yerr=after_errs,
                   color="#1b5e20", capsize=5, label="After mitigation",
                   edgecolor="#81c784", linewidth=1.2, error_kw=dict(ecolor="#81c784"))

    if exact_val is not None:
        ax.axhline(exact_val, color="#ef5350", lw=2, ls="--", label=f"Exact = {exact_val:.4f}")

    ax.set_xticks(x)
    ax.set_xticklabels(methods, color="#b0bec5", fontsize=10)
    ax.set_ylabel("⟨Z₀⟩", color="#b0bec5", fontsize=12)
    ax.set_title("Error Mitigation: ZNE and Readout Correction", color="white", fontsize=13)
    ax.tick_params(colors="#b0bec5")
    for spine in ax.spines.values():
        spine.set_edgecolor("#37474f")
    ax.legend(facecolor="#1a1d27", edgecolor="#37474f", labelcolor="white", fontsize=10)
    ax.grid(True, axis="y", color="#2a2d3e", ls="--", alpha=0.6)

    plt.tight_layout()
    if save_path is None:
        _ensure_dir(FIGURES_DIR)
        save_path = FIGURES_DIR / "fig8_mitigation.png"
    fig.savefig(save_path, dpi=DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# Utility: generate all placeholder/stub figures
# ──────────────────────────────────────────────────────────────────────────────
def generate_stub_figures() -> None:
    """Generate stub figures with placeholder data for testing the pipeline."""
    _ensure_dir(FIGURES_DIR)

    # Fig 1: coupling maps (need processors)
    from src.processors import PROCESSOR_REGISTRY, populate_all_graph_metrics
    populate_all_graph_metrics()
    procs = [PROCESSOR_REGISTRY[p] for p in ("A0", "A1", "A_star", "B0", "B_star")]
    fig1_coupling_maps(procs)
    print(f"  Fig 1 saved: {FIGURES_DIR / 'fig1_coupling_maps.png'}")

    # Fig 2: placeholder trajectory
    t_grid = np.linspace(0, 2.0, 21)
    exact = np.cos(t_grid)
    dummy_results = {
        "t_grid": t_grid, "exact": exact,
        "proc_ids": ["A_star"],
        "n_steps_list": [1, 2, 4],
        ("A_star", 1, "noisy"): exact * 0.85 + 0.05 * np.random.randn(21),
        ("A_star", 2, "noisy"): exact * 0.92 + 0.03 * np.random.randn(21),
        ("A_star", 4, "noisy"): exact * 0.95 + 0.02 * np.random.randn(21),
    }
    fig2_trajectory(dummy_results)
    print(f"  Fig 2 saved: {FIGURES_DIR / 'fig2_trajectory.png'}")

    # Fig 3: U-curve
    n_vals = np.array([1, 2, 3, 4, 6, 8], dtype=float)
    errors_ny = np.array([0.35, 0.18, 0.12, 0.11, 0.16, 0.28])
    errors_nl = np.array([0.30, 0.15, 0.10, 0.075, 0.05, 0.037])
    from src.metrics import trotter_error_model_fit
    fit = trotter_error_model_fit(n_vals, errors_nl, errors_ny)
    fig3_u_curve(n_vals, errors_ny, fit, fit["n_star_pred"], errors_nl)
    print(f"  Fig 3 saved: {FIGURES_DIR / 'fig3_u_curve.png'}")

    # Fig 4: n* map
    t_grid_2 = np.linspace(0.5, 2.0, 8)
    noise_grid = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
    n_star_grid = np.array([
        [4, 3, 2, 2, 1],
        [5, 4, 3, 2, 1],
        [6, 5, 4, 3, 2],
        [8, 6, 4, 3, 2],
        [8, 7, 5, 4, 3],
        [8, 8, 6, 4, 3],
        [8, 8, 7, 5, 4],
        [8, 8, 8, 6, 4],
    ]).T  # shape (len(noise), len(t))
    fig4_n_star_map(n_star_grid, t_grid_2, noise_grid)
    print(f"  Fig 4 saved: {FIGURES_DIR / 'fig4_n_star_map.png'}")

    # Fig 5: scaling law
    fig5_scaling_law([0.25, 0.5, 1.0, 2.0, 4.0], [6, 5, 4, 3, 2],
                     fit_slope=-0.48, r_squared=0.92)
    print(f"  Fig 5 saved: {FIGURES_DIR / 'fig5_scaling_law.png'}")

    # Fig 6: CX depth
    dummy_circ = {
        "n_steps": [1, 2, 4],
        "proc_ids": ["A_star"],
        ("A_star", "open", "cx"): [12, 24, 48],
        ("A_star", "open", "depth"): [10, 20, 40],
        ("A_star", "ring", "cx"): [18, 36, 72],
        ("A_star", "ring", "depth"): [15, 30, 60],
    }
    fig6_cx_depth(dummy_circ)
    print(f"  Fig 6 saved: {FIGURES_DIR / 'fig6_cx_depth.png'}")

    # Fig 7: layout study
    dummy_df = pd.DataFrame({
        "layout": [tuple(range(4))] * 20,
        "cx": np.random.randint(10, 50, 20),
        "depth": np.random.randint(8, 40, 20),
        "swap": np.random.randint(0, 5, 20),
    })
    fig7_layout(dummy_df)
    print(f"  Fig 7 saved: {FIGURES_DIR / 'fig7_layout.png'}")

    # Fig 8: mitigation
    fig8_mitigation(
        before={"Raw (n=2)": (0.65, 0.02), "Raw (n=4)": (0.58, 0.02)},
        after={"ZNE (n=2)": (0.72, 0.03), "Readout (n=4)": (0.64, 0.02)},
        exact_val=0.75,
    )
    print(f"  Fig 8 saved: {FIGURES_DIR / 'fig8_mitigation.png'}")


if __name__ == "__main__":
    print("Generating stub figures...")
    generate_stub_figures()
    print("Done.")
