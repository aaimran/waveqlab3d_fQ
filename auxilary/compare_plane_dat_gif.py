#!/usr/bin/env python3
"""Compare plane (.plane) and station (.dat) time series on the same plane.

Usage
-----
python3 compare_plane_dat_gif.py plane/ seismogram/

Produces a 2x3 GIF:
  Row 1: plane
  Row 2: dat (stations mapped onto the same x/z grid)
  Cols: vx, vy, vz
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import struct
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

warnings.filterwarnings("ignore")


STRUCT_HEADER = "<iiiidi"
HEADER_SIZE = struct.calcsize(STRUCT_HEADER)
BLOCK_PATTERN = re.compile(r"_block(\d+)\.plane$")

AXIS_MAP = {
	1: {"name": "yz", "dim1": "y", "dim2": "z"},
	2: {"name": "xz", "dim1": "x", "dim2": "z"},
	3: {"name": "xy", "dim1": "x", "dim2": "y"},
}


def parse_args() -> argparse.Namespace:
	ap = argparse.ArgumentParser(description="Compare plane vs dat time series on the same plane")
	ap.add_argument("plane_dir", type=str, help="Directory containing *.plane files")
	ap.add_argument("seisdir", type=str, help="Seismogram base directory (contains block*/)")
	ap.add_argument("-d", "--domain", default="plane/domain.txt", help="Path to domain description file")
	ap.add_argument("-o", "--output", default="plane_vs_dat_compare.gif", help="Output GIF")
	ap.add_argument("--time-stride", type=int, default=1, help="Use every Nth time step")
	ap.add_argument("--fps", type=int, default=12, help="GIF frames per second")
	ap.add_argument("--dpi", type=int, default=120, help="GIF DPI")
	ap.add_argument("--linthresh", type=float, default=1e-6, help="SymLogNorm linthresh")
	ap.add_argument("--y-tol", type=float, default=1e-6, help="Tolerance for selecting station y==plane_y")
	ap.add_argument("--title-fontsize", type=int, default=18, help="Main title font size")
	ap.add_argument("--subtitle-fontsize", type=int, default=14, help="Subtitle font size")
	return ap.parse_args()


def parse_domain(path: str) -> Dict[str, Dict]:
	domain: Dict[str, Dict] = {"global_axes": {}, "blocks": {}}
	current_block: Optional[int] = None
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


def linear_coords(axis_range: Tuple[float, float], n: int) -> np.ndarray:
	start, end = axis_range
	if n <= 1:
		return np.array([start], dtype=float)
	return np.linspace(start, end, n)


def read_plane_file(path: str) -> Dict:
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
		"n1": int(n1),
		"n2": int(n2),
		"times": times,
		"vx": vx,
		"vy": vy,
		"vz": vz,
	}


def gather_plane_summaries(folder: str) -> Dict[int, Dict]:
	plane_paths = sorted(glob.glob(os.path.join(folder, "*.plane")))
	if not plane_paths:
		raise FileNotFoundError(f"No .plane files found in {folder}")

	summaries: Dict[int, Dict] = {}
	for path in plane_paths:
		match = BLOCK_PATTERN.search(os.path.basename(path))
		if not match:
			raise ValueError(f"Cannot determine block number from {path}")
		block_num = int(match.group(1))
		summaries[block_num] = read_plane_file(path)
	return summaries


def assemble_plane_grid(
	summaries: Dict[int, Dict],
	domain: Dict[str, Dict],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, str, float]:
	block_nums = sorted(summaries)
	first = summaries[block_nums[0]]
	dim1 = first["dim1"]
	dim2 = first["dim2"]
	axis_name = str(first["axis_name"])
	coord_value = float(first["coord"])
	times = first["times"]
	nt = times.size

	for summary in summaries.values():
		if summary["times"].size != nt or not np.allclose(summary["times"], times):
			raise ValueError("Plane files do not share the same time axis")
		if summary["axis_name"] != axis_name:
			raise ValueError("Plane files do not share the same plane axis")

	coords_dim1_segments: List[np.ndarray] = []
	total_n1 = 0
	vx_blocks = []
	vy_blocks = []
	vz_blocks = []

	for block in block_nums:
		summary = summaries[block]
		axis_range = domain["blocks"].get(block, {}).get(dim1) or domain["global_axes"].get(dim1)
		if axis_range is None:
			raise KeyError(f"Missing {dim1} range for block {block}")
		coords = linear_coords(axis_range, int(summary["n1"]))
		coords_dim1_segments.append(coords)
		total_n1 += int(summary["n1"])
		vx_blocks.append(summary["vx"])
		vy_blocks.append(summary["vy"])
		vz_blocks.append(summary["vz"])

	coords_x = np.concatenate(coords_dim1_segments)
	dim2_range = domain["global_axes"].get(dim2) or domain["blocks"][block_nums[0]].get(dim2)
	if dim2_range is None:
		raise KeyError(f"Missing {dim2} range in domain description")
	coords_z = linear_coords(dim2_range, int(first["n2"]))

	plane_vx = np.full((nt, total_n1, int(first["n2"])), np.nan, dtype=float)
	plane_vy = np.full_like(plane_vx, np.nan)
	plane_vz = np.full_like(plane_vx, np.nan)

	offset = 0
	for vx_b, vy_b, vz_b in zip(vx_blocks, vy_blocks, vz_blocks):
		n1 = vx_b.shape[1]
		plane_vx[:, offset : offset + n1, :] = vx_b
		plane_vy[:, offset : offset + n1, :] = vy_b
		plane_vz[:, offset : offset + n1, :] = vz_b
		offset += n1

	return coords_x, coords_z, times, plane_vx, plane_vy, plane_vz, axis_name, coord_value


def parse_xyz_from_filename(path: Path) -> Optional[Tuple[float, float, float]]:
	stem = path.stem
	if "_exact_" in stem:
		return None
	parts = stem.split("_")
	if len(parts) < 4:
		return None
	tail = parts[-3:]
	try:
		x = float(tail[0])
		y = float(tail[1])
		z = float(tail[2])
		return x, y, z
	except ValueError:
		return None


def load_station_series(path: Path) -> Optional[Tuple[np.ndarray, np.ndarray]]:
	try:
		data = np.loadtxt(path, dtype=float)
	except Exception:
		return None

	if data.ndim == 1:
		if data.shape[0] < 4:
			return None
		data = data.reshape(1, -1)
	if data.shape[1] < 4:
		return None
	t = data[:, 0]
	v = data[:, 1:4]
	return t, v


def _min_positive_spacing(coords: np.ndarray) -> Optional[float]:
	coords = np.asarray(coords, dtype=float)
	if coords.size < 2:
		return None
	diffs = np.diff(np.sort(np.unique(coords)))
	diffs = diffs[diffs > 0]
	if diffs.size == 0:
		return None
	return float(diffs.min())


def _nearest_index(sorted_coords: np.ndarray, value: float) -> int:
	idx = int(np.searchsorted(sorted_coords, value))
	if idx <= 0:
		return 0
	if idx >= sorted_coords.size:
		return int(sorted_coords.size - 1)
	left = sorted_coords[idx - 1]
	right = sorted_coords[idx]
	return idx - 1 if abs(value - left) <= abs(value - right) else idx


def build_dat_grids_on_plane(
	seisdir: Path,
	plane_times: np.ndarray,
	plane_x: np.ndarray,
	plane_z: np.ndarray,
	plane_axis: str,
	plane_coord_value: float,
	y_tol: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
	if str(plane_axis).lower() != "xz":
		raise ValueError(f"This comparer currently expects an xz-plane; got {plane_axis}")

	block_dirs = sorted(
		[p for p in seisdir.iterdir() if p.is_dir() and p.name.startswith("block")],
		key=lambda p: int(re.sub(r"\D+", "", p.name) or 0),
	)
	if not block_dirs:
		raise FileNotFoundError(f"No block* directories found under {seisdir}")

	# Prefer earlier blocks for duplicate xyz.
	by_xyz: Dict[Tuple[float, float, float], Path] = {}
	for bdir in block_dirs:
		for p in sorted(bdir.glob("*.dat")):
			xyz = parse_xyz_from_filename(p)
			if xyz is None:
				continue
			if xyz not in by_xyz:
				by_xyz[xyz] = p

	if not by_xyz:
		raise RuntimeError("No xyz-named station .dat files found.")

	# Filter to y == plane y
	selected = [(xyz, p) for (xyz, p) in by_xyz.items() if abs(float(xyz[1]) - plane_coord_value) <= float(y_tol)]
	if not selected:
		raise RuntimeError(
			f"No station files matched y={plane_coord_value} within tol={y_tol}. "
			f"(Found {len(by_xyz)} total xyz-named station files.)"
		)

	# Prepare target grids: (nt, nz, nx)
	nt = int(plane_times.size)
	nx = int(plane_x.size)
	nz = int(plane_z.size)
	vx = np.full((nt, nz, nx), np.nan, dtype=np.float32)
	vy = np.full_like(vx, np.nan)
	vz = np.full_like(vx, np.nan)

	dx = _min_positive_spacing(plane_x)
	dz = _min_positive_spacing(plane_z)
	x_tol = 0.55 * dx if dx is not None else 1e-6
	z_tol = 0.55 * dz if dz is not None else 1e-6

	# Ensure sorted coordinate arrays for nearest lookup
	x_sorted = np.asarray(plane_x, dtype=float)
	z_sorted = np.asarray(plane_z, dtype=float)

	import matplotlib

	matplotlib.use("Agg")
	from tqdm import tqdm

	for (x, _y, z), path in tqdm(selected, desc="Loading station .dat", unit="file"):
		loaded = load_station_series(path)
		if loaded is None:
			continue
		t, v = loaded

		# Interpolate to plane time axis if needed
		if len(t) != nt or np.max(np.abs(t - plane_times)) > 1e-9:
			v_interp = np.empty((nt, 3), dtype=float)
			for ci in range(3):
				v_interp[:, ci] = np.interp(plane_times, t, v[:, ci])
			v = v_interp

		xi = _nearest_index(x_sorted, float(x))
		zi = _nearest_index(z_sorted, float(z))
		if abs(float(x) - float(x_sorted[xi])) > x_tol:
			continue
		if abs(float(z) - float(z_sorted[zi])) > z_tol:
			continue

		vx[:, zi, xi] = v[:, 0]
		vy[:, zi, xi] = v[:, 1]
		vz[:, zi, xi] = v[:, 2]

	return vx, vy, vz


def main() -> None:
	args = parse_args()

	plane_dir = Path(args.plane_dir)
	seisdir = Path(args.seisdir)
	domain_file = Path(args.domain)
	if not domain_file.is_file():
		raise FileNotFoundError(f"domain definition not found: {domain_file}")

	domain = parse_domain(str(domain_file))
	summaries = gather_plane_summaries(str(plane_dir))
	coords_x, coords_z, times, plane_vx, plane_vy, plane_vz, axis_name, coord_value = assemble_plane_grid(
		summaries, domain
	)

	if axis_name.lower() != "xz":
		raise SystemExit(f"Expected xz-plane for this comparison; got {axis_name}")

	dat_vx, dat_vy, dat_vz = build_dat_grids_on_plane(
		seisdir=seisdir,
		plane_times=times,
		plane_x=coords_x,
		plane_z=coords_z,
		plane_axis=axis_name,
		plane_coord_value=coord_value,
		y_tol=float(args.y_tol),
	)

	# Plot / animate
	import matplotlib

	matplotlib.use("Agg")
	import matplotlib.colors as mcolors
	import matplotlib.pyplot as plt
	from matplotlib.animation import PillowWriter
	from tqdm import tqdm

	extent = (float(coords_x[0]), float(coords_x[-1]), float(coords_z[0]), float(coords_z[-1]))

	plane_grids = [plane_vx, plane_vy, plane_vz]
	dat_grids = [dat_vx, dat_vy, dat_vz]
	titles = ["vx", "vy", "vz"]

	# Shared per-component norm across (plane, dat)
	norms = []
	for comp in range(3):
		pmin = float(np.nanmin(plane_grids[comp]))
		pmax = float(np.nanmax(plane_grids[comp]))
		dmin = float(np.nanmin(dat_grids[comp]))
		dmax = float(np.nanmax(dat_grids[comp]))
		vmax_abs = max(abs(pmin), abs(pmax), abs(dmin), abs(dmax))
		if not np.isfinite(vmax_abs) or vmax_abs == 0.0:
			vmax_abs = 1.0
		norms.append(
			mcolors.SymLogNorm(
				linthresh=max(1e-30, float(args.linthresh)),
				vmin=-vmax_abs,
				vmax=vmax_abs,
			)
		)

	cmap = plt.get_cmap("seismic").copy()
	cmap.set_bad(color=(0.9, 0.9, 0.9))

	fig, axes = plt.subplots(2, 3, figsize=(25, 12), constrained_layout=True)

	# Row labels
	fig.text(0.01, 0.72, "plane", rotation=90, va="center", ha="left", fontsize=14, fontweight="bold")
	fig.text(0.01, 0.28, "dat", rotation=90, va="center", ha="left", fontsize=14, fontweight="bold")

	main_title = f"timeseries comparison of particle velocities at xz-plane at y={coord_value:g} km (plane vs dat)"
	fig.suptitle(main_title, fontsize=int(args.title_fontsize), fontweight="bold")
	subtitle = fig.text(0.5, 0.92, "time = 0 s", ha="center", va="center", fontsize=int(args.subtitle_fontsize))

	ims_plane = []
	ims_dat = []
	for j in range(3):
		axp = axes[0, j]
		axd = axes[1, j]

		im_p = axp.imshow(
			np.full((coords_z.size, coords_x.size), np.nan, dtype=float),
			origin="lower",
			aspect="equal",
			interpolation="nearest",
			extent=extent,
			cmap=cmap,
			norm=norms[j],
		)
		im_d = axd.imshow(
			np.full((coords_z.size, coords_x.size), np.nan, dtype=float),
			origin="lower",
			aspect="equal",
			interpolation="nearest",
			extent=extent,
			cmap=cmap,
			norm=norms[j],
		)
		ims_plane.append(im_p)
		ims_dat.append(im_d)

		axp.set_title(titles[j])
		axd.set_xlabel("x (km)")
		if j == 0:
			axp.set_ylabel("z (km)")
			axd.set_ylabel("z (km)")

		# One shared colorbar per column (shared across the two rows)
		fig.colorbar(im_p, ax=[axp, axd], shrink=0.85)

	frame_indices = list(range(0, int(times.size), max(1, int(args.time_stride))))
	writer = PillowWriter(fps=max(1, int(args.fps)))

	out_path = Path(args.output)
	out_path.parent.mkdir(parents=True, exist_ok=True)

	with writer.saving(fig, str(out_path), dpi=int(args.dpi)):
		for ti in tqdm(frame_indices, desc="Writing GIF frames", unit="frame"):
			# plane grids are (nt, nx, nz) and need transpose to (nz, nx)
			for j in range(3):
				ims_plane[j].set_data(plane_grids[j][ti].T)
				ims_dat[j].set_data(dat_grids[j][ti])
			subtitle.set_text(f"time = {float(times[ti]):g} s")
			fig.canvas.draw()
			writer.grab_frame()

	plt.close(fig)
	print(f"Saved comparison GIF to {out_path}")


if __name__ == "__main__":
	main()

