"""
generate_3d_gif.py
Generates a publication-grade, hardware-accelerated animated GIF of the 
StepWise-Q 3D Quantum Simulation Laboratory for Processor D (6-Qubit 2x3 Ladder).
Embeddable natively in GitHub README.md without iframe security blocks.
"""

import os
import io
import numpy as np
import scipy.sparse.linalg as sla
import scipy.linalg as la
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image

def get_paulis():
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return I, X, Y, Z

def kron_n(ops):
    res = ops[0]
    for op in ops[1:]:
        res = np.kron(res, op)
    return res

def build_hamiltonian(N=4, J=1.0, g=1.0):
    I, X, Y, Z = get_paulis()
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    
    # ZZ terms on periodic ring: (0,1), (1,2), (2,3), (3,0)
    bonds = [(0, 1), (1, 2), (2, 3), (3, 0)]
    for q1, q2 in bonds:
        op_list = [I] * N
        op_list[q1] = Z
        op_list[q2] = Z
        H += -J * kron_n(op_list)
        
    # X terms
    for q in range(N):
        op_list = [I] * N
        op_list[q] = X
        H += -g * kron_n(op_list)
        
    return H

def get_bloch_vectors(psi, N=4):
    I, X, Y, Z = get_paulis()
    dim = 2**N
    rho = np.outer(psi, np.conj(psi))
    bloch_vecs = []
    
    for q in range(N):
        op_x = kron_n([X if i == q else I for i in range(N)])
        op_y = kron_n([Y if i == q else I for i in range(N)])
        op_z = kron_n([Z if i == q else I for i in range(N)])
        
        rx = float(np.real(np.trace(rho @ op_x)))
        ry = float(np.real(np.trace(rho @ op_y)))
        rz = float(np.real(np.trace(rho @ op_z)))
        bloch_vecs.append((rx, ry, rz))
        
    return bloch_vecs

def main():
    print("Computing quantum dynamics for 3D simulation animation...")
    N = 4
    J = 1.0
    g = 1.0
    T = 2.0
    H = build_hamiltonian(N, J, g)
    
    # Neel state |1010> -> binary index 10
    psi0 = np.zeros(2**N, dtype=complex)
    psi0[10] = 1.0
    
    # Generate 36 animation frames across t in [0, T]
    num_frames = 36
    t_vals = np.linspace(0, T, num_frames)
    
    # Precompute states
    states = [sla.expm_multiply(-1j * H * t, psi0) for t in t_vals]
    ms_vals = []
    for psi in states:
        bloch = get_bloch_vectors(psi, N)
        ms = 0.25 * (bloch[0][2] - bloch[1][2] + bloch[2][2] - bloch[3][2])
        ms_vals.append(ms)
        
    # 2x3 Ladder layout for Processor D (6 Qubits)
    # Positions: (x, y, z)
    # Row 0 (top): Q0(-1.5, 1.2, 0), Q1(0, 1.2, 0), Q2(1.5, 1.2, 0)
    # Row 1 (bot): Q3(-1.5, -1.2, 0), Q4(0, -1.2, 0), Q5(1.5, -1.2, 0)
    node_coords = {
        0: np.array([-1.8, 1.1, 0.0]),
        1: np.array([0.0, 1.1, 0.0]),
        2: np.array([1.8, 1.1, 0.0]),
        3: np.array([-1.8, -1.1, 0.0]),
        4: np.array([0.0, -1.1, 0.0]),
        5: np.array([1.8, -1.1, 0.0]),
    }
    
    # Logical to Physical mapping for 4-qubit ring:
    # Logical Q0 -> Phys 0
    # Logical Q1 -> Phys 1
    # Logical Q2 -> Phys 4
    # Logical Q3 -> Phys 3
    # Ring bonds: (0,1), (1,4), (4,3), (3,0) -> Native 4-cycle on Ladder!
    log_to_phys = {0: 0, 1: 1, 2: 4, 3: 3}
    edges = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
    ring_edges = [(0, 1), (1, 4), (4, 3), (3, 0)]
    
    # Sphere mesh for Bloch spheres
    u = np.linspace(0, 2 * np.pi, 16)
    v = np.linspace(0, np.pi, 10)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones(np.size(u)), np.cos(v))
    sphere_radius = 0.55
    
    frames = []
    print("Rendering 3D frames...")
    
    for idx, (t, psi, ms) in enumerate(zip(t_vals, states, ms_vals)):
        fig = plt.figure(figsize=(12, 6.2), facecolor='#08090a')
        
        # 3D Subplot for Architecture & Bloch Spheres
        ax3d = fig.add_axes([0.02, 0.05, 0.58, 0.90], projection='3d', facecolor='#08090a')
        ax3d.set_box_aspect((1.6, 1.0, 0.8))
        
        # Rotating camera angle
        azimuth = -60 + (idx * 360.0 / num_frames) * 0.4
        elevation = 28 + 6 * np.sin(idx * 2 * np.pi / num_frames)
        ax3d.view_init(elev=elevation, azim=azimuth)
        
        # Draw ladder edges
        for u_n, v_n in edges:
            p1 = node_coords[u_n]
            p2 = node_coords[v_n]
            is_ring = (u_n, v_n) in ring_edges or (v_n, u_n) in ring_edges
            color = '#00f0ff' if is_ring else '#334155'
            lw = 3.0 if is_ring else 1.5
            alpha = 0.95 if is_ring else 0.45
            ax3d.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 
                      color=color, linewidth=lw, alpha=alpha, zorder=2)
            
        bloch_vecs = get_bloch_vectors(psi, N)
        
        # Draw nodes & Bloch spheres
        for phys_idx in range(6):
            coord = node_coords[phys_idx]
            is_active = phys_idx in log_to_phys.values()
            
            if is_active:
                log_idx = [k for k, v in log_to_phys.items() if v == phys_idx][0]
                rx, ry, rz = bloch_vecs[log_idx]
                
                # Draw transparent Bloch sphere wireframe
                ax3d.plot_wireframe(coord[0] + sx * sphere_radius, 
                                    coord[1] + sy * sphere_radius, 
                                    coord[2] + sz * sphere_radius, 
                                    color='#00f0ff', alpha=0.15, linewidth=0.6)
                
                # Draw central sphere node
                ax3d.scatter([coord[0]], [coord[1]], [coord[2]], 
                             color='#00e5ff', s=120, edgecolors='white', linewidths=1.2, alpha=0.9, zorder=5)
                
                # Draw 3D Bloch vector arrow
                arrow_len = sphere_radius * 1.15
                ax3d.quiver(coord[0], coord[1], coord[2], 
                            rx * arrow_len, ry * arrow_len, rz * arrow_len,
                            color='#ff007f', linewidth=2.4, arrow_length_ratio=0.32, zorder=6)
                
                # Label
                ax3d.text(coord[0], coord[1], coord[2] + sphere_radius + 0.25, 
                          f"Q{log_idx} (P{phys_idx})", color='white', 
                          fontsize=8, fontweight='bold', ha='center')
            else:
                # Auxiliary / Spectator qubit on ladder (P2 or P5)
                ax3d.scatter([coord[0]], [coord[1]], [coord[2]], 
                             color='#475569', s=70, edgecolors='#94a3b8', linewidths=0.8, alpha=0.6)
                ax3d.text(coord[0], coord[1], coord[2] + 0.35, 
                          f"P{phys_idx} (Idle)", color='#94a3b8', fontsize=7, ha='center')

        # Clean 3D view styling
        ax3d.grid(False)
        ax3d.set_xticks([])
        ax3d.set_yticks([])
        ax3d.set_zticks([])
        ax3d.xaxis.pane.fill = False
        ax3d.yaxis.pane.fill = False
        ax3d.zaxis.pane.fill = False
        ax3d.xaxis.pane.set_edgecolor('#08090a')
        ax3d.yaxis.pane.set_edgecolor('#08090a')
        ax3d.zaxis.pane.set_edgecolor('#08090a')
        ax3d.set_xlim([-2.8, 2.8])
        ax3d.set_ylim([-2.2, 2.2])
        ax3d.set_zlim([-1.5, 1.5])
        
        ax3d.set_title("PROCESSOR D (6-QUBIT 2×3 LADDER) — LIVE 3D DYNAMICS\nNative 4-Cycle Ring Embedding (Q0-Q1-Q2-Q3) with Zero SWAPs", 
                       color='#00e5ff', fontsize=10, fontweight='bold', pad=-15)
        
        # -------------------------------------------------------------
        # Right Side 2D Panels: Magnetization Curve & State Probabilities
        # -------------------------------------------------------------
        # 1. Staggered Magnetization curve
        ax_mag = fig.add_axes([0.65, 0.56, 0.31, 0.34], facecolor='#0f172a')
        ax_mag.plot(t_vals, ms_vals, color='#38bdf8', lw=2.0, label=r'Exact $M_s(t)$')
        ax_mag.axhline(0.652444, color='#10b981', ls='--', lw=1.2, alpha=0.8, label=r'Target $M_s(T)=0.6524$')
        # Current progress point
        ax_mag.scatter([t], [ms], color='#ff007f', s=80, edgecolors='white', zorder=5)
        ax_mag.axvline(t, color='#ff007f', ls=':', lw=1.0, alpha=0.7)
        ax_mag.set_xlim([0, T])
        ax_mag.set_ylim([-0.8, 0.85])
        ax_mag.set_ylabel(r"Staggered $M_s$", color='white', fontsize=8)
        ax_mag.tick_params(colors='#94a3b8', labelsize=7)
        ax_mag.spines['top'].set_visible(False)
        ax_mag.spines['right'].set_visible(False)
        ax_mag.spines['left'].set_color('#334155')
        ax_mag.spines['bottom'].set_color('#334155')
        ax_mag.grid(True, color='#1e293b', ls=':', lw=0.6)
        ax_mag.legend(loc='lower left', fontsize=7, facecolor='#0f172a', edgecolor='#334155', labelcolor='white')
        ax_mag.set_title(f"TIME EVOLUTION: t = {t:.2f} / {T:.1f} s | Ms = {ms:+.4f}", 
                         color='white', fontsize=8, fontweight='bold')

        # 2. 16-State Amplitude Bar Chart
        ax_bar = fig.add_axes([0.65, 0.12, 0.31, 0.32], facecolor='#0f172a')
        probs = np.abs(psi)**2
        indices = np.arange(16)
        bar_colors = ['#10b981' if i == 10 else '#0ea5e9' for i in indices]
        ax_bar.bar(indices, probs, color=bar_colors, edgecolor='#0284c7', lw=0.6, alpha=0.85)
        ax_bar.set_xticks(indices)
        ax_bar.set_xticklabels([f"{i:04b}" for i in indices], rotation=90, fontsize=6, color='#94a3b8')
        ax_bar.set_xlim([-0.8, 15.8])
        ax_bar.set_ylim([0, 1.05])
        ax_bar.set_ylabel("Probability", color='white', fontsize=8)
        ax_bar.tick_params(colors='#94a3b8', labelsize=7)
        ax_bar.spines['top'].set_visible(False)
        ax_bar.spines['right'].set_visible(False)
        ax_bar.spines['left'].set_color('#334155')
        ax_bar.spines['bottom'].set_color('#334155')
        ax_bar.grid(True, axis='y', color='#1e293b', ls=':', lw=0.6)
        ax_bar.set_title("16-BASIS COMPUTATIONAL PROBABILITIES (|1010> Green)", 
                         color='white', fontsize=8, fontweight='bold')
        
        # Save frame to buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100, facecolor=fig.get_facecolor(), edgecolor='none')
        plt.close(fig)
        buf.seek(0)
        img = Image.open(buf)
        frames.append(img.convert('RGB'))
        print(f"  Frame {idx+1}/{num_frames} rendered.")

    out_gif = "results/figures/simulation_3d_demo.gif"
    print(f"Compiling animated GIF to {out_gif}...")
    
    # Optimize palette and save looping GIF
    frames[0].save(
        out_gif,
        save_all=True,
        append_images=frames[1:],
        duration=90,  # 90 ms per frame -> ~3.2 second loop
        loop=0,
        optimize=True
    )
    
    size_mb = os.path.getsize(out_gif) / (1024 * 1024)
    print(f"Successfully generated {out_gif} ({size_mb:.2f} MB)")

if __name__ == "__main__":
    main()
