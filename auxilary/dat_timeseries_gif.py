#!/usr/bin/env python3
"""Generate a VX/VY/VZ time-series GIF from seismogram station outputs.

Assumptions about station files:
- Each station file is plain text with columns: time vx vy vz
- Files live under a per-block directory, e.g. seismogram/block1/*.dat

"Common stations" means stations that exist in BOTH block1 and block2.
This script will (by default) select that intersection and then use the
block1 values to build the animation.

Examples
--------
# From simulation/Test directory
python3 ../../python/make_station_timeseries_gif.py \
  --seisdir seismogram \
  --out vxvyvz_block1_common.gif

# Limit frames and stations for quick iteration
python3 ../../python/make_station_timeseries_gif.py \
  --seisdir seismogram --out quick.gif --time-stride 10 --max-stations 300

# If you only have block1 files (no common filtering)
python3 ../../python/make_station_timeseries_gif.py \
  --block1 seismogram/block1 --no-common
"""

from __future__ import annotations

import argparse
import math
import os
import re
import multiprocessing as mp
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

try:
    from tqdm import tqdm

    HAS_TQDM = True
except Exception:
    HAS_TQDM = False


@dataclass(frozen=True)
class StationSeries:
    key: str
    path: Path
    x: float
    z: float
    t: np.ndarray  # shape (nt,)
    v: np.ndarray  # shape (nt, 3) for vx,vy,vz


@dataclass(frozen=True)
class BlockStats:
    name: str
    total_dat_files: int
    station_files: int
    x_range: Tuple[float, float]
    y_range: Tuple[float, float]
    z_range: Tuple[float, float]
    dx: Optional[float]
    dy: Optional[float]
    dz: Optional[float]


_FLOAT_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$")


def _is_float(token: str) -> bool:
    return bool(_FLOAT_RE.match(token.strip()))


def _read_first_data_line(path: Path) -> Optional[List[str]]:
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                # ignore potential markers
                if s.startswith("!---"):
                    continue
                parts = s.split()
                return parts
    except OSError:
        return None
    return None


def is_station_timeseries_file(path: Path) -> bool:
    """Heuristic: station time series has ~4 numeric columns per line."""
    parts = _read_first_data_line(path)
    if not parts:
        return False
    if not (4 <= len(parts) <= 6):
        return False
    return all(_is_float(p) for p in parts[:4])


def station_key_from_filename(path: Path) -> Optional[str]:
    stem = path.stem
    # Skip exact-moment outputs for now (requested vx/vy/vz)
    if "_exact_" in stem:
        return None
    # Normalize away trailing _blockN
    stem = re.sub(r"_block\d+$", "", stem)
    return stem


def parse_xz_from_filename(path: Path, key: str) -> Optional[Tuple[float, float]]:
    """Try to recover (x,z) from xyz-based naming, else fall back to index naming."""
    stem = path.stem

    # xyz naming: <name>_<x>_<y>_<z>.dat  (block subdir avoids collisions)
    parts = stem.split("_")
    if len(parts) >= 4:
        tail = parts[-3:]
        try:
            x = float(tail[0])
            _y = float(tail[1])
            z = float(tail[2])
            return x, z
        except ValueError:
            pass

    # index naming: <name>_<i>_<j>_<k>_blockN.dat -> key is without _blockN
    parts = key.split("_")
    if len(parts) >= 3:
        try:
            i = float(int(parts[-3]))
            _j = float(int(parts[-2]))
            k = float(int(parts[-1]))
            return i, k
        except ValueError:
            return None
    return None


def parse_xyz_from_filename(path: Path) -> Optional[Tuple[float, float, float]]:
    """Parse x,y,z from xyz-style seismogram filenames.

    Expected filename: <name>_<x>_<y>_<z>.dat (no block suffix when using per-block dirs)
    """
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
        # single line
        if data.shape[0] < 4:
            return None
        data = data.reshape(1, -1)

    if data.shape[1] < 4:
        return None

    t = data[:, 0]
    v = data[:, 1:4]
    return t, v


def _cpu_available() -> Tuple[int, Dict[str, str]]:
    n = os.cpu_count() or 1
    env_keys = [
        "SLURM_CPUS_PER_TASK",
        "SLURM_CPUS_ON_NODE",
        "SLURM_JOB_CPUS_PER_NODE",
        "OMP_NUM_THREADS",
    ]
    env = {k: os.environ.get(k, "") for k in env_keys if os.environ.get(k)}
    return n, env


def _resolution(values: np.ndarray) -> Optional[float]:
    values = np.unique(values)
    if values.size < 2:
        return None
    values = np.sort(values)

    # For a uniform grid with (N+1) points, spacing should be (max-min)/N.
    span = float(values[-1] - values[0])
    n = int(values.size)
    if n <= 1:
        return None
    expected = span / float(n - 1)

    diffs = np.diff(values)
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return None

    # If diffs are (almost) constant, prefer expected spacing.
    tol = max(1e-12, 1e-6 * abs(expected))
    if np.all(np.abs(diffs - expected) <= tol):
        return float(np.round(expected, 12))
    return float(np.round(float(diffs.min()), 12))


def _collect_station_paths_by_xyz_fast(block_dir: Path) -> Dict[Tuple[float, float, float], Path]:
    """Fast discovery of station series files based only on filename xyz pattern.

    Avoids opening files (important for 10k+ files). This is reliable for this
    codebase because station outputs are named <name>_<x>_<y>_<z>.dat.
    """
    if not block_dir.exists():
        return {}
    out: Dict[Tuple[float, float, float], Path] = {}
    for p in sorted(block_dir.glob("*.dat")):
        xyz = parse_xyz_from_filename(p)
        if xyz is None:
            continue
        out[xyz] = p
    return out


def _block_stats(block_name: str, block_dir: Path) -> BlockStats:
    total_dat = len(list(block_dir.glob("*.dat"))) if block_dir.exists() else 0
    xyz_map = _collect_station_paths_by_xyz_fast(block_dir)
    if not xyz_map:
        return BlockStats(
            name=block_name,
            total_dat_files=total_dat,
            station_files=0,
            x_range=(math.nan, math.nan),
            y_range=(math.nan, math.nan),
            z_range=(math.nan, math.nan),
            dx=None,
            dy=None,
            dz=None,
        )

    xs = np.array([k[0] for k in xyz_map.keys()], dtype=float)
    ys = np.array([k[1] for k in xyz_map.keys()], dtype=float)
    zs = np.array([k[2] for k in xyz_map.keys()], dtype=float)
    return BlockStats(
        name=block_name,
        total_dat_files=total_dat,
        station_files=len(xyz_map),
        x_range=(float(xs.min()), float(xs.max())),
        y_range=(float(ys.min()), float(ys.max())),
        z_range=(float(zs.min()), float(zs.max())),
        dx=_resolution(xs),
        dy=_resolution(ys),
        dz=_resolution(zs),
    )


def _load_station_worker(path_str: str):
    path = Path(path_str)
    loaded = load_station_series(path)
    if loaded is None:
        return path_str, None
    t, v = loaded
    return path_str, (t, v)


def collect_station_files(block_dir: Path) -> Dict[str, Path]:
    if not block_dir.exists():
        return {}
    out: Dict[str, Path] = {}
    # On some filesystems, repeated stat() calls are slow; glob("*.dat") already
    # yields directory entries we can attempt to parse.
    for p in sorted(block_dir.glob("*.dat")):
        k = station_key_from_filename(p)
        if not k:
            continue
        if not is_station_timeseries_file(p):
            continue
        out[k] = p
    return out


def collect_station_files_by_xyz(block_dir: Path) -> Dict[Tuple[float, float, float], Path]:
    """Collect station series files keyed by (x,y,z) parsed from filename.

    This is used for stitching blocks: filenames already encode physical coordinates.
    """
    if not block_dir.exists():
        return {}

    out: Dict[Tuple[float, float, float], Path] = {}
    for p in sorted(block_dir.glob("*.dat")):
        xyz = parse_xyz_from_filename(p)
        if xyz is None:
            continue
        out[xyz] = p
    return out


def select_common_keys(block1: Dict[str, Path], block2: Dict[str, Path]) -> List[str]:
    return sorted(set(block1.keys()) & set(block2.keys()))


def build_series(
    block1_dir: Path,
    block2_dir: Optional[Path],
    common: bool,
    max_stations: Optional[int],
    jobs: int,
) -> List[StationSeries]:
    b1 = collect_station_files(block1_dir)
    if not b1:
        raise RuntimeError(f"No station time-series .dat files found in: {block1_dir}")

    keys: List[str]
    if common and block2_dir is not None and block2_dir.exists():
        b2 = collect_station_files(block2_dir)
        keys = select_common_keys(b1, b2)
        if not keys:
            # Some configurations (e.g. non-overlapping station sets across blocks)
            # may legitimately have no intersection.
            print(
                f"Warning: no common stations found between {block1_dir} and {block2_dir}; "
                "using all block1 stations instead."
            )
            keys = sorted(b1.keys())
    else:
        keys = sorted(b1.keys())

    if max_stations is not None and len(keys) > max_stations:
        # deterministic downsample
        step = max(1, len(keys) // max_stations)
        keys = keys[::step][:max_stations]

    series: List[StationSeries] = []
    base_t: Optional[np.ndarray] = None

    metas: List[Tuple[str, str, Optional[Tuple[float, float]]]] = []
    for k in keys:
        p = b1[k]
        metas.append((k, str(p), parse_xz_from_filename(p, k)))

    path_to_meta = {p_str: (k, xz) for (k, p_str, xz) in metas}
    paths = [p_str for (_k, p_str, _xz) in metas]

    if jobs <= 0:
        jobs = os.cpu_count() or 1
    jobs = max(1, jobs)

    def _iter_loaded():
        if jobs == 1 or len(paths) < 50:
            it: Iterable[str] = paths
            if HAS_TQDM:
                it = tqdm(paths, total=len(paths), desc="Parsing block1", unit="file")
            for p in it:
                loaded = load_station_series(Path(p))
                yield p, loaded
        else:
            ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
            with ctx.Pool(processes=jobs) as pool:
                it = pool.imap_unordered(_load_station_worker, paths, chunksize=200)
                if HAS_TQDM:
                    it = tqdm(it, total=len(paths), desc="Parsing block1", unit="file")
                for p_str, loaded in it:
                    yield p_str, loaded

    for p_str, loaded in _iter_loaded():
        if loaded is None:
            continue
        meta = path_to_meta.get(p_str)
        if meta is None:
            continue
        k, xz = meta
        if xz is None:
            continue
        p = Path(p_str)
        t, v = loaded

        if base_t is None:
            base_t = t
        else:
            # If time grids differ, interpolate to base_t
            if len(t) != len(base_t) or np.max(np.abs(t - base_t)) > 1e-9:
                v_interp = np.empty((len(base_t), 3), dtype=float)
                for ci in range(3):
                    v_interp[:, ci] = np.interp(base_t, t, v[:, ci])
                t = base_t
                v = v_interp

        series.append(StationSeries(key=k, path=p, x=xz[0], z=xz[1], t=t, v=v))

    if not series:
        raise RuntimeError("No usable station series were parsed (could not infer x/z or read data).")

    return series


def build_series_stitched(
    block_dirs: List[Path],
    max_stations: Optional[int],
    jobs: int,
) -> List[StationSeries]:
    """Stitch station outputs from multiple block directories into one set.

    Stations are matched by (x,y,z) encoded in filenames.
    """
    # Collect per-block so we can show per-block progress during parsing.
    by_block: List[Tuple[str, Dict[Tuple[float, float, float], Path]]] = []
    for bdir in block_dirs:
        name = bdir.name
        found = collect_station_files_by_xyz(bdir)
        by_block.append((name, found))

    # Union with block order preference (earlier blocks win if duplicates).
    by_xyz: Dict[Tuple[float, float, float], Path] = {}
    for _name, found in by_block:
        for xyz, p in found.items():
            if xyz not in by_xyz:
                by_xyz[xyz] = p

    if not by_xyz:
        raise RuntimeError("No station time-series .dat files found for stitching.")

    items = sorted(by_xyz.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][1]))

    if max_stations is not None and len(items) > max_stations:
        step = max(1, len(items) // max_stations)
        items = items[::step][:max_stations]

    paths = [str(p) for _xyz, p in items]

    if jobs <= 0:
        jobs = os.cpu_count() or 1
    jobs = max(1, jobs)

    # Load one reference file to define the canonical time axis.
    ref_loaded = load_station_series(Path(paths[0]))
    if ref_loaded is None:
        raise RuntimeError(f"Failed to load reference station file: {paths[0]}")
    base_t, _ = ref_loaded

    # Parallel load per-block, with a progress bar per block.
    ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context("spawn")
    results: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    # Build reverse lookup: path -> xyz for final items
    selected_paths = set(paths)
    for block_name, found in by_block:
        block_paths = [str(p) for xyz, p in found.items() if str(p) in selected_paths]
        if not block_paths:
            continue

        if jobs == 1 or len(block_paths) < 50:
            it: Iterable[str] = block_paths
            if HAS_TQDM:
                it = tqdm(block_paths, total=len(block_paths), desc=f"Parsing {block_name}", unit="file")
            for p in it:
                loaded = load_station_series(Path(p))
                if loaded is None:
                    continue
                results[p] = loaded
        else:
            with ctx.Pool(processes=jobs) as pool:
                it = pool.imap_unordered(_load_station_worker, block_paths, chunksize=200)
                if HAS_TQDM:
                    it = tqdm(it, total=len(block_paths), desc=f"Parsing {block_name}", unit="file")
                for p_str, loaded in it:
                    if loaded is None:
                        continue
                    results[p_str] = loaded

    series: List[StationSeries] = []
    for (x, y, z), path in items:
        p_str = str(path)
        if p_str not in results:
            continue
        t, v = results[p_str]
        # Interpolate if needed
        if len(t) != len(base_t) or np.max(np.abs(t - base_t)) > 1e-9:
            v_interp = np.empty((len(base_t), 3), dtype=float)
            for ci in range(3):
                v_interp[:, ci] = np.interp(base_t, t, v[:, ci])
            t = base_t
            v = v_interp

        # Keep a key that is stable and readable
        key = f"{x:.3f}_{y:.3f}_{z:.3f}"
        series.append(StationSeries(key=key, path=path, x=float(x), z=float(z), t=t, v=v))

    if not series:
        raise RuntimeError("No usable station series were parsed for stitched blocks.")
    return series


def make_gif(
    series: List[StationSeries],
    out_path: Path,
    time_stride: int,
    fps: int,
    symmetric_clim: bool,
    dpi: int,
    linthresh: float,
    layout: str,
    shared_clim: bool,
    plane_y_km: float,
    title_fontsize: int,
) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.animation import PillowWriter

    t = series[0].t
    nt = len(t)

    xs = np.array([s.x for s in series], dtype=float)
    zs = np.array([s.z for s in series], dtype=float)

    x_unique = np.unique(xs)
    z_unique = np.unique(zs)
    x_unique.sort()
    z_unique.sort()

    # For imshow extent: treat values as cell centers
    x_min, x_max = float(x_unique.min()), float(x_unique.max())
    z_min, z_max = float(z_unique.min()), float(z_unique.max())
    if len(x_unique) > 1:
        dx = float(np.min(np.diff(x_unique)))
    else:
        dx = 1.0
    if len(z_unique) > 1:
        dz = float(np.min(np.diff(z_unique)))
    else:
        dz = 1.0
    extent = (x_min - 0.5 * dx, x_max + 0.5 * dx, z_min - 0.5 * dz, z_max + 0.5 * dz)

    x_index = {float(x): i for i, x in enumerate(x_unique)}
    z_index = {float(z): i for i, z in enumerate(z_unique)}

    nx = len(x_unique)
    nz = len(z_unique)

    # Precompute full grids once: vx[nt,nz,nx], vy[nt,nz,nx], vz[nt,nz,nx]
    # This makes "all stations + all timesteps" runs much faster.
    zi = np.array([z_index[float(z)] for z in zs], dtype=int)
    xi = np.array([x_index[float(x)] for x in xs], dtype=int)

    vx_grid = np.full((nt, nz, nx), np.nan, dtype=np.float32)
    vy_grid = np.full_like(vx_grid, np.nan)
    vz_grid = np.full_like(vx_grid, np.nan)

    v_stack = np.stack([s.v.astype(np.float32, copy=False) for s in series], axis=1)  # (nt, nstations, 3)
    vx_grid[:, zi, xi] = v_stack[:, :, 0]
    vy_grid[:, zi, xi] = v_stack[:, :, 1]
    vz_grid[:, zi, xi] = v_stack[:, :, 2]

    # Compute color limits from grids (ignoring NaNs)
    vmin = np.array([
        float(np.nanmin(vx_grid)),
        float(np.nanmin(vy_grid)),
        float(np.nanmin(vz_grid)),
    ])
    vmax = np.array([
        float(np.nanmax(vx_grid)),
        float(np.nanmax(vy_grid)),
        float(np.nanmax(vz_grid)),
    ])
    if symmetric_clim:
        vmax_abs = np.maximum(np.abs(vmin), np.abs(vmax))
        vmin = -vmax_abs
        vmax = vmax_abs

    titles = ["vx", "vy", "vz"]
    grids = [vx_grid, vy_grid, vz_grid]

    cmap = plt.get_cmap("seismic").copy()
    cmap.set_bad(color=(0.9, 0.9, 0.9))

    frame_indices = list(range(0, nt, max(1, time_stride)))
    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _make_one_gif(path: Path, comps: List[int], layout_mode: str) -> None:
        if layout_mode == "horizontal":
            nrows, ncols = 1, len(comps)
            figsize = (25, 10) if len(comps) == 3 else (10, 8)
        else:
            # Default to vertical stacking
            nrows, ncols = len(comps), 1
            figsize = (10, 25) if len(comps) == 3 else (10, 8)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
        axes_list = list(np.ravel(np.array(axes)))

        def _header(time_value: float) -> str:
            return f"Particle velocities on xz-plane at y = {plane_y_km:g} km at time = {time_value:.6g}"

        # If requested and plotting multiple components, share color scale + a single colorbar.
        use_shared = bool(shared_clim) and len(comps) > 1
        if use_shared:
            shared_vmin = float(np.min(vmin[comps]))
            shared_vmax = float(np.max(vmax[comps]))
            shared_norm = mcolors.SymLogNorm(
                linthresh=max(1e-30, float(linthresh)),
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
                    linthresh=max(1e-30, float(linthresh)),
                    vmin=float(vmin[comp]),
                    vmax=float(vmax[comp]),
                )
            im = ax.imshow(
                np.full((nz, nx), np.nan, dtype=float),
                origin="lower",
                aspect="equal",
                interpolation="nearest",
                extent=extent,
                cmap=cmap,
                norm=norm,
            )
            ax.set_title(titles[comp])
            ax.set_xlabel("x (km)")
            ax.set_ylabel("z (km)")
            images.append((comp, im))

        if use_shared:
            # One colorbar for all subplots
            fig.colorbar(images[0][1], ax=axes_list, shrink=0.8)
        else:
            # One colorbar per subplot
            for _comp, im in images:
                fig.colorbar(im, ax=im.axes, shrink=0.8)

        time_text = fig.suptitle(
            _header(0.0),
            fontsize=int(title_fontsize),
            fontweight="bold",
        )

        writer = PillowWriter(fps=max(1, fps))
        frame_iter: Iterable[int] = frame_indices
        if HAS_TQDM:
            if len(comps) == 1:
                desc = f"Writing {titles[comps[0]]} GIF frames"
            else:
                desc = "Writing GIF frames"
            frame_iter = tqdm(frame_indices, desc=desc, unit="frame")

        with writer.saving(fig, str(path), dpi=dpi):
            for ti in frame_iter:
                for comp, im in images:
                    im.set_data(grids[comp][ti])
                time_text.set_text(_header(float(t[ti])))
                fig.canvas.draw()
                writer.grab_frame()
        plt.close(fig)

    if layout == "separate":
        # If out_path ends with .gif, use it as a stem, else append suffix.
        if out_path.suffix.lower() == ".gif":
            base = out_path.with_suffix("")
            vx_path = Path(str(base) + "_vx.gif")
            vy_path = Path(str(base) + "_vy.gif")
            vz_path = Path(str(base) + "_vz.gif")
        else:
            vx_path = Path(str(out_path) + "_vx.gif")
            vy_path = Path(str(out_path) + "_vy.gif")
            vz_path = Path(str(out_path) + "_vz.gif")

        _make_one_gif(vx_path, [0], "vertical")
        _make_one_gif(vy_path, [1], "vertical")
        _make_one_gif(vz_path, [2], "vertical")
        return [vx_path, vy_path, vz_path]

    layout_mode = "horizontal" if layout == "horizontal" else "vertical"
    _make_one_gif(out_path, [0, 1, 2], layout_mode)
    return [out_path]


def main() -> None:
    ap = argparse.ArgumentParser(description="Make vx/vy/vz station time-series GIF (block1 common stations).")
    ap.add_argument("--seisdir", type=str, default=None, help="Base seismogram directory (contains block1/, block2/)")
    ap.add_argument("--block1", type=str, default=None, help="Path to block1 output directory")
    ap.add_argument("--block2", type=str, default=None, help="Path to block2 output directory")
    ap.add_argument("--stitch", action="store_true", help="Stitch all block* directories under --seisdir (e.g., block1+block2)")
    ap.add_argument("--blocks", type=str, default=None, help="Comma-separated block numbers to stitch (e.g. 1,2). Default: all block* dirs")
    ap.add_argument("--jobs", type=int, default=0, help="Parallel processes for reading station files (0=all available)")
    ap.add_argument("--no-common", action="store_true", help="Do not restrict to stations common with block2")
    ap.add_argument("--max-stations", type=int, default=None, help="Downsample stations to at most this many")
    ap.add_argument("--time-stride", type=int, default=5, help="Use every Nth time step")
    ap.add_argument("--fps", type=int, default=12, help="GIF frames per second")
    ap.add_argument("--linthresh", type=float, default=1e-6, help="SymLogNorm linthresh (like plane GIF)")
    ap.add_argument("--no-symmetric-clim", action="store_true", help="Use min/max instead of symmetric color limits")
    ap.add_argument(
        "--shared-clim",
        action="store_true",
        help="Use a shared color scale for vx/vy/vz and draw a single colorbar when stacked",
    )
    ap.add_argument(
        "--plane-y-km",
        type=float,
        default=0.0,
        help="Plane y coordinate to show in the title (km)",
    )
    ap.add_argument(
        "--title-fontsize",
        type=int,
        default=18,
        help="Header title font size",
    )
    ap.add_argument(
        "--layout",
        type=str,
        default="vertical",
        choices=["vertical", "horizontal", "separate"],
        help="Plot layout: vertical (3x1), horizontal (1x3), or separate (three GIFs)",
    )
    ap.add_argument("--dpi", type=int, default=120, help="GIF DPI")
    ap.add_argument("--out", type=str, default="vxvyvz.gif", help="Output GIF path")

    args = ap.parse_args()

    if args.seisdir is not None:
        seisdir = Path(args.seisdir)
        block1_dir = seisdir / "block1"
        block2_dir = seisdir / "block2"
    else:
        if args.block1 is None:
            raise SystemExit("Provide --seisdir or --block1")
        block1_dir = Path(args.block1)
        block2_dir = Path(args.block2) if args.block2 is not None else None

    out_path = Path(args.out)

    # --- Show processors available
    cpu_n, cpu_env = _cpu_available()
    jobs = cpu_n if args.jobs == 0 else max(1, int(args.jobs))
    print(f"Processors available (os.cpu_count): {cpu_n}")
    if cpu_env:
        print("Scheduler/env CPU hints:")
        for k, v in cpu_env.items():
            print(f"  {k}={v}")
    print(f"Parallel processes used for parsing: {jobs}")

    # --- Per-block station file counts and ranges
    block1_stats = _block_stats("block1", block1_dir)
    block2_stats = _block_stats("block2", block2_dir) if block2_dir is not None else None
    print("Files per block:")
    print(
        f"  block1: total .dat={block1_stats.total_dat_files}, station .dat (xyz-named)={block1_stats.station_files}"
    )
    if block2_stats is not None:
        print(
            f"  block2: total .dat={block2_stats.total_dat_files}, station .dat (xyz-named)={block2_stats.station_files}"
        )

    def _print_ranges(label: str, st: BlockStats):
        if st.station_files == 0:
            print(f"{label}: (no xyz-named station files found)")
            return
        print(f"{label} ranges:")
        print(f"  x: {st.x_range[0]} .. {st.x_range[1]}")
        print(f"  y: {st.y_range[0]} .. {st.y_range[1]}")
        print(f"  z: {st.z_range[0]} .. {st.z_range[1]}")
        print(f"{label} resolution:")
        print(f"  dx={st.dx}  dy={st.dy}  dz={st.dz}")

    _print_ranges("block1", block1_stats)
    if block2_stats is not None:
        _print_ranges("block2", block2_stats)

    # --- Global ranges/resolution + common points between block1 and block2
    if block2_dir is not None and block2_dir.exists():
        b1_xyz = set(_collect_station_paths_by_xyz_fast(block1_dir).keys())
        b2_xyz = set(_collect_station_paths_by_xyz_fast(block2_dir).keys())
        common_xyz = b1_xyz & b2_xyz
        print(f"Common points between blocks (x,y,z match): {len(common_xyz)}")

        all_xyz = list(b1_xyz | b2_xyz)
        if all_xyz:
            xs = np.array([p[0] for p in all_xyz], dtype=float)
            ys = np.array([p[1] for p in all_xyz], dtype=float)
            zs = np.array([p[2] for p in all_xyz], dtype=float)
            print("Global ranges (block1+block2):")
            print(f"  x: {float(xs.min())} .. {float(xs.max())}")
            print(f"  y: {float(ys.min())} .. {float(ys.max())}")
            print(f"  z: {float(zs.min())} .. {float(zs.max())}")
            print("Global resolution (block1+block2):")
            print(f"  dx={_resolution(xs)}  dy={_resolution(ys)}  dz={_resolution(zs)}")

    if args.stitch:
        if args.seisdir is None:
            raise SystemExit("--stitch requires --seisdir")

        blocks: Optional[List[int]] = None
        if args.blocks:
            blocks = [int(x.strip()) for x in args.blocks.split(",") if x.strip()]

        if blocks is None:
            block_dirs = sorted(
                [p for p in Path(args.seisdir).iterdir() if p.is_dir() and p.name.startswith("block")],
                key=lambda p: int(re.sub(r"\D+", "", p.name) or 0),
            )
        else:
            block_dirs = [Path(args.seisdir) / f"block{b}" for b in blocks]

        series = build_series_stitched(
            block_dirs=block_dirs,
            max_stations=args.max_stations,
            jobs=jobs,
        )
    else:
        series = build_series(
            block1_dir=block1_dir,
            block2_dir=block2_dir,
            common=(not args.no_common),
            max_stations=args.max_stations,
            jobs=jobs,
        )

    written = make_gif(
        series=series,
        out_path=out_path,
        time_stride=max(1, args.time_stride),
        fps=max(1, args.fps),
        symmetric_clim=(not args.no_symmetric_clim),
        dpi=args.dpi,
        linthresh=args.linthresh,
        layout=args.layout,
        shared_clim=args.shared_clim,
        plane_y_km=float(args.plane_y_km),
        title_fontsize=int(args.title_fontsize),
    )

    if len(written) == 1:
        print(f"Wrote GIF: {written[0]}")
    else:
        for p in written:
            print(f"Wrote GIF: {p}")


if __name__ == "__main__":
    main()
