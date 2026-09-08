#!/usr/bin/env python3
"""Visualize WaveQLab3D-Q's curvilinear two-block mesh for TPV36.

The formulas and defaults reproduce inputfile/TPV36.in and the
``analytical_tpv36`` branch in src/grid.f90.  Coordinates are in km.

The SCEC description uses (x, y, z) = (along strike, vertical depth,
horizontal dip direction).  WaveQLab's mesh instead uses
(X, Y, Z) = (horizontal dip direction, vertical depth, along strike),
with an X translation of -15 km at the surface fault trace.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import plotly.graph_objects as go
except ImportError:  # Report the optional dependency clearly at runtime.
    go = None


DIP_DEG = 15.0
DIP = np.deg2rad(DIP_DEG)
X_MIN, X_MAX = -30.0, 30.0
Y_MIN, Y_MAX = 0.0, 15.0
Z_MIN, Z_MAX = -17.5, 17.5
FAULT_DOWN_DIP = 28.0
FAULT_HALF_STRIKE = 15.0
HYPO_DOWN_DIP = 18.0
BEND_DEPTH = 8.0


def interface_x(y: np.ndarray | float) -> np.ndarray:
    """Exact analytical_tpv36 shared-interface formula from grid.f90."""
    y = np.asarray(y, dtype=float)
    planar = y / np.tan(DIP) - 15.0
    bent = (
        BEND_DEPTH / np.tan(DIP)
        - 15.0
        + 0.25 * np.arctan(4.0 * (y - BEND_DEPTH))
        # grid.f90 contains exp(-0*(y-8)), which is identically one.
    )
    return np.where(y <= BEND_DEPTH, planar, bent)


def block_x(q: np.ndarray | float, y: np.ndarray | float, block: int) -> np.ndarray:
    """X coordinate for a constant computational-q mesh line.

    For this z-invariant geometry, the transfinite interpolation reduces in
    an X-Y cross-section to interpolation between each block's two X faces.
    q runs from 0 to 1 within a block.
    """
    q = np.asarray(q, dtype=float)
    xf = interface_x(y)
    if block == 1:
        return X_MIN + q * (xf - X_MIN)
    if block == 2:
        return xf + q * (X_MAX - xf)
    raise ValueError("block must be 1 or 2")


def draw_cross_section(ax: plt.Axes, nq_lines: int, ny_lines: int) -> None:
    y_dense = np.linspace(Y_MIN, Y_MAX, 600)
    y_rows = np.linspace(Y_MIN, Y_MAX, ny_lines)
    colors = {1: "#277da1", 2: "#f8961e"}

    for block in (1, 2):
        for q in np.linspace(0.0, 1.0, nq_lines):
            ax.plot(block_x(q, y_dense, block), y_dense, color=colors[block], lw=0.55)
        for y in y_rows:
            q = np.linspace(0.0, 1.0, 150)
            ax.plot(block_x(q, y, block), np.full_like(q, y), color=colors[block], lw=0.55)

    fault_y_end = FAULT_DOWN_DIP * np.sin(DIP)
    fault_y = np.linspace(0.0, fault_y_end, 200)
    extension_y = np.linspace(fault_y_end, Y_MAX, 200)
    ax.plot(interface_x(fault_y), fault_y, color="crimson", lw=3.0, label="active 28 km fault")
    ax.plot(
        interface_x(extension_y), extension_y, color="black", lw=2.0,
        ls="--", label="locked interface extension",
    )

    hypo_y = HYPO_DOWN_DIP * np.sin(DIP)
    ax.scatter(interface_x(hypo_y), hypo_y, marker="*", s=150, color="gold",
               edgecolor="black", zorder=6, label="hypocenter")
    ax.axhline(BEND_DEPTH, color="0.35", ls=":", lw=1.2)
    ax.annotate("atan bend begins (Y = 8 km)", (interface_x(BEND_DEPTH), BEND_DEPTH),
                xytext=(-105, 12), textcoords="offset points", fontsize=8,
                arrowprops={"arrowstyle": "->", "lw": 0.8})
    ax.set(xlabel="WaveQLab X: horizontal dip direction (km)",
           ylabel="WaveQLab Y: depth (km)",
           title="X-Y cross-section: fitted two-block mesh")
    ax.set_xlim(X_MIN - 1, X_MAX + 1)
    ax.set_ylim(Y_MAX + 0.4, Y_MIN - 0.4)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    ax.legend(loc="lower left", fontsize=8)


def draw_mapping(ax: plt.Axes) -> None:
    w = np.linspace(0.0, FAULT_DOWN_DIP, 250)
    y = w * np.sin(DIP)
    pdf_horizontal = w * np.cos(DIP)
    waveqlab_horizontal = interface_x(y) + 15.0
    ax.plot(w, pdf_horizontal, lw=3, color="crimson", label=r"PDF planar fault: $w\cos15^\circ$")
    ax.plot(w, waveqlab_horizontal, "--", lw=2, color="black",
            label="WaveQLab active-interface formula + 15 km")
    ax.axvline(HYPO_DOWN_DIP, color="goldenrod", ls=":", label="hypocenter, w=18 km")
    ax.set(xlabel="distance down dip, w (km)", ylabel="horizontal distance from trace (km)",
           title="Benchmark geometry vs. implemented active fault")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    error_m = 1000.0 * np.max(np.abs(pdf_horizontal - waveqlab_horizontal))
    ax.text(0.03, 0.08, f"maximum difference on active fault: {error_m:.2e} m",
            transform=ax.transAxes, fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"})


def make_interactive_3d(nq_lines: int, ny_lines: int):
    """Return a Plotly figure with independently toggleable mesh layers."""
    if go is None:
        raise SystemExit("Plotly is required for 3-D output: python3 -m pip install plotly")

    fig = go.Figure()
    fault_y_end = FAULT_DOWN_DIP * np.sin(DIP)
    z = np.linspace(-FAULT_HALF_STRIKE, FAULT_HALF_STRIKE, 90)
    y = np.linspace(0.0, fault_y_end, 65)
    zz, yy = np.meshgrid(z, y)
    xx = interface_x(yy)
    fig.add_trace(go.Surface(
        x=xx, y=zz, z=yy,
        name="active 28 km fault",
        colorscale=[[0, "crimson"], [1, "crimson"]],
        opacity=0.62,
        showscale=False,
        hovertemplate="X=%{x:.2f} km<br>Z=%{y:.2f} km<br>depth=%{z:.2f} km<extra>active fault</extra>",
    ))

    # Volume grid lines on three along-strike slices. None separators let each
    # block be one toggleable Plotly trace instead of dozens of traces.
    y_full = np.linspace(Y_MIN, Y_MAX, 180)
    for block, color in ((1, "#277da1"), (2, "#f8961e")):
        xs, zs, depths = [], [], []
        for z0 in (Z_MIN, 0.0, Z_MAX):
            for q in np.linspace(0.0, 1.0, nq_lines):
                xs.extend(block_x(q, y_full, block).tolist() + [None])
                zs.extend(np.full_like(y_full, z0).tolist() + [None])
                depths.extend(y_full.tolist() + [None])
            for depth in np.linspace(Y_MIN, Y_MAX, ny_lines):
                q = np.linspace(0.0, 1.0, 100)
                xs.extend(block_x(q, depth, block).tolist() + [None])
                zs.extend(np.full_like(q, z0).tolist() + [None])
                depths.extend(np.full_like(q, depth).tolist() + [None])
        fig.add_trace(go.Scatter3d(
            x=xs, y=zs, z=depths, mode="lines", name=f"block {block} mesh slices",
            line={"color": color, "width": 2}, opacity=0.72,
            hoverinfo="skip",
        ))

    # Show the full shared interface, including the bent locked extension.
    y_extension = np.linspace(fault_y_end, Y_MAX, 180)
    for z0, show_legend in ((Z_MIN, True), (0.0, False), (Z_MAX, False)):
        fig.add_trace(go.Scatter3d(
            x=interface_x(y_extension), y=np.full_like(y_extension, z0), z=y_extension,
            mode="lines", name="locked interface extension", legendgroup="extension",
            showlegend=show_legend, line={"color": "black", "width": 5, "dash": "dash"},
            hovertemplate="X=%{x:.2f} km<br>Z=%{y:.2f} km<br>depth=%{z:.2f} km<extra>locked extension</extra>",
        ))

    hypo_y = HYPO_DOWN_DIP * np.sin(DIP)
    fig.add_trace(go.Scatter3d(
        x=[float(interface_x(hypo_y))], y=[0.0], z=[hypo_y], mode="markers",
        name="hypocenter", marker={"size": 8, "color": "gold", "symbol": "diamond",
                                    "line": {"color": "black", "width": 2}},
        hovertemplate="w=18 km<br>X=%{x:.3f} km<br>Z=0 km<br>depth=%{z:.3f} km<extra>hypocenter</extra>",
    ))
    fig.update_layout(
        title="TPV36: interactive curvilinear two-block mesh and active fault",
        template="plotly_white",
        legend={"x": 0.01, "y": 0.99},
        margin={"l": 0, "r": 0, "t": 55, "b": 0},
        scene={
            "xaxis": {"title": "X: dip direction (km)", "range": [X_MIN, X_MAX]},
            "yaxis": {"title": "Z: along strike (km)", "range": [Z_MIN, Z_MAX]},
            "zaxis": {"title": "Y: depth (km)", "range": [Y_MAX, Y_MIN]},
            "aspectmode": "manual",
            "aspectratio": {"x": 1.5, "y": 0.9, "z": 0.55},
            "camera": {"eye": {"x": 1.55, "y": -1.6, "z": 0.9}},
        },
    )
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="directory for all three outputs (default: current directory)")
    parser.add_argument("--q-lines", type=int, default=13,
                        help="constant-q lines per block")
    parser.add_argument("--y-lines", type=int, default=16,
                        help="constant-depth lines in the cross-section")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--show", action="store_true", help="also open an interactive window")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.q_lines < 2 or args.y_lines < 2:
        raise SystemExit("--q-lines and --y-lines must be at least 2")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    cross_path = args.output_dir / "tpv36_mesh_cross_section.png"
    fig_cross, ax_cross = plt.subplots(figsize=(12, 7), constrained_layout=True)
    draw_cross_section(ax_cross, args.q_lines, args.y_lines)
    fig_cross.savefig(cross_path, dpi=args.dpi, bbox_inches="tight")

    comparison_path = args.output_dir / "tpv36_geometry_comparison.png"
    fig_comparison, ax_comparison = plt.subplots(figsize=(11, 7), constrained_layout=True)
    draw_mapping(ax_comparison)
    fig_comparison.savefig(comparison_path, dpi=args.dpi, bbox_inches="tight")

    interactive_path = args.output_dir / "tpv36_mesh_3d_interactive.html"
    interactive = make_interactive_3d(max(5, args.q_lines // 2), max(5, args.y_lines // 2))
    interactive.write_html(interactive_path, include_plotlyjs=True, full_html=True)

    for path in (cross_path, comparison_path, interactive_path):
        print(f"wrote {path.resolve()}")
    if args.show:
        interactive.show()
        plt.show()


if __name__ == "__main__":
    main()
