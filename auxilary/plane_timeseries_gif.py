#!/usr/bin/env python3
"""Create a GIF from all .plane slices described by a domain file."""

from __future__ import annotations

import argparse
import glob
import os
import re
import struct
import warnings
from typing import Dict, List, TypedDict

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

warnings.filterwarnings("ignore")

try:
    from PIL import Image

    HAS_PIL = True
except ImportError:
    HAS_PIL = False

STRUCT_HEADER = "<iiiidi"
HEADER_SIZE = struct.calcsize(STRUCT_HEADER)
BLOCK_PATTERN = re.compile(r"_block(\d+)\.plane$")

AXIS_MAP = {
    1: {"name": "yz", "dim1": "y", "dim2": "z"},
    2: {"name": "xz", "dim1": "x", "dim2": "z"},
    3: {"name": "xy", "dim1": "x", "dim2": "y"},
}


class PlaneSummary(TypedDict):
    axis_name: str
    dim1: str
    dim2: str
    coord: float
    fixed_index: int
    n1: int
    n2: int
    times: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray


class DomainDesc(TypedDict):
    global_axes: Dict[str, tuple[float, float]]
    blocks: Dict[int, Dict[str, tuple[float, float]]]


def parse_domain(path: str) -> DomainDesc:
    domain: DomainDesc = {"global_axes": {}, "blocks": {}}
    current_block: int | None = None

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            key = parts[0].lower()
            if key == "block":
                current_block = int(parts[1])
                domain["blocks"].setdefault(current_block, {})
                continue
            values = [float(v) for v in parts[1:]]
            if not values:
                continue
            axis_range = (values[0], values[-1]) if len(values) > 1 else (values[0], values[0])
            if current_block is None:
                domain["global_axes"][key] = axis_range
            else:
                domain["blocks"][current_block][key] = axis_range

    return domain


def linear_coords(axis_range: tuple[float, float], n: int) -> np.ndarray:
    start, end = axis_range
    if n <= 1:
        return np.array([start])
    return np.linspace(start, end, n)


def read_plane_file(path: str) -> PlaneSummary:
    with open(path, "rb") as f:
        header = f.read(HEADER_SIZE)
        if len(header) != HEADER_SIZE:
            raise IOError("Incomplete header")
        _, axis_id, n1, n2, coord, fixed_index = struct.unpack(STRUCT_HEADER, header)
        if axis_id not in AXIS_MAP:
            raise ValueError(f"Unsupported axis_id {axis_id}")
        raw = np.fromfile(f, dtype="<f8")

    per_step = 1 + 3 * n1 * n2
    if raw.size % per_step != 0:
        raise ValueError("Plane payload size mismatch")

    nt = raw.size // per_step
    records = raw.reshape((nt, per_step))
    times = records[:, 0].copy()

    vx = np.empty((nt, n1, n2), dtype=float)
    vy = np.empty((nt, n1, n2), dtype=float)
    vz = np.empty((nt, n1, n2), dtype=float)

    for step in range(nt):
        snapshot = records[step, 1:]
        for comp, arr in enumerate((vx, vy, vz)):
            offset = comp * n1 * n2
            block = snapshot[offset : offset + n1 * n2]
            arr[step] = block.reshape((n1, n2), order="F")

    axis_info = AXIS_MAP[axis_id]
    return {
        "axis_name": axis_info["name"],
        "dim1": axis_info["dim1"],
        "dim2": axis_info["dim2"],
        "coord": float(coord),
        "fixed_index": int(fixed_index),
        "n1": n1,
        "n2": n2,
        "times": times,
        "vx": vx,
        "vy": vy,
        "vz": vz,
    }


def gather_plane_summaries(folder: str) -> Dict[int, PlaneSummary]:
    plane_paths = sorted(glob.glob(os.path.join(folder, "*.plane")))
    if not plane_paths:
        raise FileNotFoundError(f"No .plane files found in {folder}")

    summaries: Dict[int, PlaneSummary] = {}
    for path in plane_paths:
        match = BLOCK_PATTERN.search(os.path.basename(path))
        if not match:
            raise ValueError(f"Cannot determine block number from {path}")
        block_num = int(match.group(1))
        summaries[block_num] = read_plane_file(path)
    return summaries


def assemble_grid(summaries: Dict[int, PlaneSummary], domain: DomainDesc) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    block_nums = sorted(summaries)
    first = summaries[block_nums[0]]
    dim1 = first["dim1"]
    dim2 = first["dim2"]
    times = first["times"]
    nt = times.size

    for summary in summaries.values():
        if summary["times"].size != nt or not np.allclose(summary["times"], times):
            raise ValueError("Plane files do not share the same time axis")

    coords_dim1_segments: List[np.ndarray] = []
    total_n1 = 0
    data_vx = []
    data_vy = []
    data_vz = []

    for block in block_nums:
        summary = summaries[block]
        axis_range = domain["blocks"].get(block, {}).get(dim1) or domain["global_axes"].get(dim1)
        if axis_range is None:
            raise KeyError(f"Missing {dim1} range for block {block}")
        coords = linear_coords(axis_range, summary["n1"])
        coords_dim1_segments.append(coords)
        total_n1 += summary["n1"]
        data_vx.append(summary["vx"])
        data_vy.append(summary["vy"])
        data_vz.append(summary["vz"])

    coords_dim1 = np.concatenate(coords_dim1_segments)
    dim2_range = domain["global_axes"].get(dim2) or domain["blocks"][block_nums[0]].get(dim2)
    if dim2_range is None:
        raise KeyError(f"Missing {dim2} range in domain description")
    coords_dim2 = linear_coords(dim2_range, first["n2"])

    global_vx = np.full((nt, total_n1, first["n2"]), np.nan)
    global_vy = np.full_like(global_vx, np.nan)
    global_vz = np.full_like(global_vx, np.nan)

    offset = 0
    for vx_block, vy_block, vz_block in zip(data_vx, data_vy, data_vz):
        n1 = vx_block.shape[1]
        global_vx[:, offset : offset + n1, :] = vx_block
        global_vy[:, offset : offset + n1, :] = vy_block
        global_vz[:, offset : offset + n1, :] = vz_block
        offset += n1

    return coords_dim1, coords_dim2, times, global_vx, global_vy, global_vz


def build_frames(
    coords_x: np.ndarray,
    coords_z: np.ndarray,
    times: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    stride: int,
    labels: tuple[str, str],
    layout: str = "vertical",
    shared_clim: bool = False,
    coord_value: float = 0.0,
    plane_axis: str = "xz",
    title_fontsize: int = 18,
) -> list[np.ndarray]:
    nt = times.size
    grids = [vx, vy, vz]
    titles = ["vx", "vy", "vz"]
    
    stats = ((np.nanmin(vx), np.nanmax(vx)), (np.nanmin(vy), np.nanmax(vy)), (np.nanmin(vz), np.nanmax(vz)))
    vmin = np.array([s[0] for s in stats])
    vmax = np.array([s[1] for s in stats])
    
    extent = (coords_x[0], coords_x[-1], coords_z[0], coords_z[-1])
    
    cmap = plt.get_cmap("seismic").copy()
    cmap.set_bad(color=(0.9, 0.9, 0.9))
    
    comps = [0, 1, 2]
    
    if layout == "horizontal":
        nrows, ncols = 1, len(comps)
        figsize = (25, 10)
    else:
        nrows, ncols = len(comps), 1
        figsize = (10, 25)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
    axes_list = list(np.ravel(np.array(axes)))
    
    use_shared = shared_clim and len(comps) > 1
    if use_shared:
        shared_vmin = float(np.min(vmin[comps]))
        shared_vmax = float(np.max(vmax[comps]))
        shared_norm = mcolors.SymLogNorm(
            linthresh=1e-6,
            vmin=shared_vmin,
            vmax=shared_vmax,
        )
    else:
        shared_norm = None
    
    images = []
    for ax, comp in zip(axes_list, comps):
        norm = shared_norm
        if norm is None:
            norm = mcolors.SymLogNorm(
                linthresh=1e-6,
                vmin=float(vmin[comp]),
                vmax=float(vmax[comp]),
            )
        im = ax.imshow(
            np.zeros((coords_z.size, coords_x.size)),
            extent=extent,
            origin="lower",
            cmap=cmap,
            aspect="equal",
            norm=norm,
        )
        ax.set_title(titles[comp])
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        images.append((comp, im))
    
    if use_shared:
        fig.colorbar(images[0][1], ax=axes_list, shrink=0.8)
    else:
        for _comp, im in images:
            plt.colorbar(im, ax=im.axes, shrink=0.8)
    
    def _header(time_value: float) -> str:
        fixed_axis_map = {"yz": "x", "xz": "y", "xy": "z"}
        fixed_axis = fixed_axis_map.get(str(plane_axis).lower(), "?")
        return (
            f"Particle velocities on {plane_axis}-plane at {fixed_axis} = {coord_value:.3g} km "
            f"at time = {time_value:.6g} s"
        )
    
    time_text = fig.suptitle(
        _header(times[0]),
        fontsize=int(title_fontsize),
        fontweight="bold",
    )

    frames: List[np.ndarray] = []
    iterator = range(0, nt, stride)
    if HAS_PIL:
        iterator = tqdm(iterator, desc="Rendering frames", unit="frame")
    
    for idx in iterator:
        for comp, im in images:
            im.set_data(grids[comp][idx].T)
        time_text.set_text(_header(times[idx]))
        fig.canvas.draw()
        buf = fig.canvas.buffer_rgba()
        frame = np.frombuffer(buf, dtype=np.uint8).reshape(fig.canvas.get_width_height()[::-1] + (4,))
        frames.append(frame[:, :, :3].copy())

    plt.close(fig)
    return frames


def save_animation(frames: list[np.ndarray], output: str) -> None:
    if not frames:
        raise ValueError("No frames to save")

    if HAS_PIL:
        images = [Image.fromarray(frame) for frame in frames]
        images[0].save(output, save_all=True, append_images=images[1:], duration=100, loop=0)
        print(f"Saved GIF to {output}")
    else:
        tmp_dir = os.path.join("/tmp", "plane_frames")
        os.makedirs(tmp_dir, exist_ok=True)
        for idx, frame in enumerate(frames):
            Image.fromarray(frame).save(os.path.join(tmp_dir, f"frame_{idx:04d}.png"))
        print("PIL missing; exported PNG frames to", tmp_dir)
        print("Convert with ffmpeg: ffmpeg -i frame_%04d.png -vf 'fps=10' output.gif")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Animate plane outputs with domain metadata")
    parser.add_argument("plane_dir", help="Directory containing .plane files")
    parser.add_argument("-d", "--domain", default="plane/domain.txt", help="Path to domain description file")
    parser.add_argument("-fs", "--frame-stride", type=int, default=1, help="Frame stride")
    parser.add_argument("-o", "--output", default="plane_animation.gif", help="Output GIF path")
    parser.add_argument(
        "--layout",
        type=str,
        default="vertical",
        choices=["vertical", "horizontal", "separate"],
        help="Plot layout: vertical (3x1), horizontal (1x3), or separate (three GIFs)",
    )
    parser.add_argument(
        "--shared-clim",
        action="store_true",
        help="Use a shared color scale for vx/vy/vz and draw a single colorbar when stacked",
    )
    parser.add_argument(
        "--title-fontsize",
        type=int,
        default=18,
        help="Header title font size",
    )
    return parser.parse_args()


def build_separate_gifs(
    coords_x: np.ndarray,
    coords_z: np.ndarray,
    times: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    stride: int,
    labels: tuple[str, str],
    output_path: str,
    coord_value: float,
    plane_axis: str,
    title_fontsize: int,
) -> None:
    """Build three separate GIF files, one for each velocity component."""
    output_base = output_path.replace(".gif", "")
    grids = [(vx, "vx"), (vy, "vy"), (vz, "vz")]
    
    for grid, name in grids:
        out_path = f"{output_base}_{name}.gif"
        frames = build_frames(
            coords_x, coords_z, times,
            grid, grid, grid,  # Pass same grid 3 times but only use first
            stride, labels,
            layout="vertical",
            shared_clim=False,
            coord_value=coord_value,
            plane_axis=plane_axis,
            title_fontsize=title_fontsize,
        )
        save_animation(frames, out_path)


def main() -> None:
    args = parse_args()
    domain_file = args.domain
    if not os.path.isfile(domain_file):
        raise FileNotFoundError(f"domain definition not found: {domain_file}")

    domain = parse_domain(domain_file)
    summaries = gather_plane_summaries(args.plane_dir)
    coords_x, coords_z, times, vx, vy, vz = assemble_grid(summaries, domain)

    first_summary = summaries[next(iter(summaries))]
    axis_name = first_summary["axis_name"]
    coord_value = first_summary["coord"]
    
    print(f"Plane axis: {axis_name}, nt={times.size}")
    print(f"dim1 range: {coords_x[0]:.3f} .. {coords_x[-1]:.3f}")
    print(f"dim2 range: {coords_z[0]:.3f} .. {coords_z[-1]:.3f}")

    if args.layout == "separate":
        build_separate_gifs(
            coords_x, coords_z, times, vx, vy, vz,
            args.frame_stride, ("x (km)", "z (km)"),
            args.output, coord_value, axis_name, args.title_fontsize
        )
    else:
        frames = build_frames(
            coords_x, coords_z, times, vx, vy, vz,
            args.frame_stride, ("x (km)", "z (km)"),
            layout=args.layout,
            shared_clim=args.shared_clim,
            coord_value=coord_value,
            plane_axis=axis_name,
            title_fontsize=args.title_fontsize,
        )
        save_animation(frames, args.output)


if __name__ == "__main__":
    main()
