#!/usr/bin/env python3
"""Script to check dimensions and fields in .dat files."""

from __future__ import annotations

import glob
import os
import re
import sys
from collections import defaultdict

import numpy as np


BLOCK_PATTERN = re.compile(r"_block(\w+)$")


def parse_dat_file(filepath: str) -> dict[str, float | int] | None:
    """Return summary stats for a single .dat file."""
    data = np.loadtxt(filepath)
    if data.size == 0:
        return None
    if data.ndim == 1:
        data = data.reshape(1, -1)

    times = data[:, 0]
    vx = data[:, 1]
    vy = data[:, 2]
    vz = data[:, 3]
    return {
        "nt": data.shape[0],
        "time_min": float(times.min()),
        "time_max": float(times.max()),
        "vx_min": float(np.min(vx)),
        "vx_max": float(np.max(vx)),
        "vy_min": float(np.min(vy)),
        "vy_max": float(np.max(vy)),
        "vz_min": float(np.min(vz)),
        "vz_max": float(np.max(vz)),
    }


def block_key_from_path(filepath: str) -> str:
    stem = os.path.splitext(os.path.basename(filepath))[0]
    match = BLOCK_PATTERN.search(stem)
    return match.group(0) if match else stem


def summary_for_block(stats_list: list[dict[str, float | int]]) -> dict[str, float | int]:
    agg = {
        "files": len(stats_list),
        "nt_total": 0,
        "nt_min": float("inf"),
        "nt_max": 0,
        "time_min": float("inf"),
        "time_max": -float("inf"),
        "vx_min": float("inf"),
        "vx_max": -float("inf"),
        "vy_min": float("inf"),
        "vy_max": -float("inf"),
        "vz_min": float("inf"),
        "vz_max": -float("inf"),
    }

    for stats in stats_list:
        nt = int(stats["nt"])
        agg["nt_total"] += nt
        agg["nt_min"] = min(agg["nt_min"], nt)
        agg["nt_max"] = max(agg["nt_max"], nt)
        agg["time_min"] = min(agg["time_min"], float(stats["time_min"]))
        agg["time_max"] = max(agg["time_max"], float(stats["time_max"]))
        agg["vx_min"] = min(agg["vx_min"], float(stats["vx_min"]))
        agg["vx_max"] = max(agg["vx_max"], float(stats["vx_max"]))
        agg["vy_min"] = min(agg["vy_min"], float(stats["vy_min"]))
        agg["vy_max"] = max(agg["vy_max"], float(stats["vy_max"]))
        agg["vz_min"] = min(agg["vz_min"], float(stats["vz_min"]))
        agg["vz_max"] = max(agg["vz_max"], float(stats["vz_max"]))

    return agg


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python check_dat.py <dat_folder>")
        sys.exit(1)

    dat_folder = sys.argv[1]
    dat_files = sorted(glob.glob(os.path.join(dat_folder, "*.dat")))

    if not dat_files:
        print(f"No .dat files found in {dat_folder}")
        return

    per_block: dict[str, list[dict[str, float | int]]] = defaultdict(list)

    for filepath in dat_files:
        stats = parse_dat_file(filepath)
        if stats is None:
            continue
        block_key = block_key_from_path(filepath)
        per_block[block_key].append(stats)

    for block_key in sorted(per_block):
        block_stats = summary_for_block(per_block[block_key])
        files = block_stats["files"]
        block_label = f"Block {block_key}"
        print(f"\n{block_label} ({files} file{'s' if files != 1 else ''})")
        if files == 0:
            print("  (no data)")
            continue

        print(
            f"  nt: {block_stats['nt_min']} .. {block_stats['nt_max']} (total {block_stats['nt_total']})"
        )
        print(
            f"  time range: {block_stats['time_min']:.6f} .. {block_stats['time_max']:.6f}"
        )
        print(f"  vx range: {block_stats['vx_min']:.6e} .. {block_stats['vx_max']:.6e}")
        print(f"  vy range: {block_stats['vy_min']:.6e} .. {block_stats['vy_max']:.6e}")
        print(f"  vz range: {block_stats['vz_min']:.6e} .. {block_stats['vz_max']:.6e}")


if __name__ == "__main__":
    main()
