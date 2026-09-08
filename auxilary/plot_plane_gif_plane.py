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


def build_frames(coords_x: np.ndarray, coords_z: np.ndarray, times: np.ndarray, vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, stride: int, labels: tuple[str, str]) -> list[np.ndarray]:
    nt = times.size
    stats = ((np.nanmin(vx), np.nanmax(vx)), (np.nanmin(vy), np.nanmax(vy)), (np.nanmin(vz), np.nanmax(vz)))
    norms = [
        mcolors.SymLogNorm(linthresh=1e-6, vmin=minv, vmax=maxv)
        for minv, maxv in stats
    ]

    fig, axes = plt.subplots(3, 1, figsize=(10, 25))
    extent = (coords_x[0], coords_x[-1], coords_z[0], coords_z[-1])

    images = []
    for ax, norm, title in zip(axes, norms, ("vx", "vy", "vz")):
        im = ax.imshow(
            np.zeros((coords_z.size, coords_x.size)),
            extent=extent,
            origin="lower",
            cmap="seismic",
            aspect="equal",
            norm=norm,
        )
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_title(title)
        plt.colorbar(im, ax=ax, shrink=0.8)
        images.append(im)

    frames: List[np.ndarray] = []
    for idx in tqdm(range(0, nt, stride), desc="Rendering frames"):
        for im, plane in zip(images, (vx, vy, vz)):
            im.set_data(plane[idx].T)
        fig.suptitle(f"t={times[idx]:.3f}s")
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    domain_file = args.domain
    if not os.path.isfile(domain_file):
        raise FileNotFoundError(f"domain definition not found: {domain_file}")

    domain = parse_domain(domain_file)
    summaries = gather_plane_summaries(args.plane_dir)
    coords_x, coords_z, times, vx, vy, vz = assemble_grid(summaries, domain)

    axis_name = summaries[next(iter(summaries))]["axis_name"]
    print(f"Plane axis: {axis_name}, nt={times.size}")
    print(f"dim1 range: {coords_x[0]:.3f} .. {coords_x[-1]:.3f}")
    print(f"dim2 range: {coords_z[0]:.3f} .. {coords_z[-1]:.3f}")

    frames = build_frames(coords_x, coords_z, times, vx, vy, vz, args.frame_stride, ("x (km)", "z (km)"))
    save_animation(frames, args.output)


if __name__ == "__main__":
    main()
