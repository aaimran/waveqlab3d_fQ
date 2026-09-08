#!/usr/bin/env python3
"""Run WaveQLab3D input preflight without initializing or evolving a domain."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys


def executable_path(script_dir: Path) -> Path:
    override = os.environ.get("WAVEQLAB3D_EXE")
    candidates = [
        Path(override).expanduser() if override else None,
        script_dir / "build" / "waveqlab3d",
        script_dir / "bin" / "waveqlab3d",
    ]
    for candidate in candidates:
        if candidate and candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()

    from_path = shutil.which("waveqlab3d")
    if from_path:
        return Path(from_path).resolve()

    raise FileNotFoundError(
        "waveqlab3d executable not found; build it or set WAVEQLAB3D_EXE"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run all WaveQLab3D preflight checks and report the resolved MPI "
            "decomposition and time-step parameters without starting the simulation."
        )
    )
    parser.add_argument("input_file", type=Path, metavar="fname.in")
    parser.add_argument("np", type=int, help="number of MPI processes")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_file = args.input_file.expanduser().resolve()
    if not input_file.is_file():
        print(f"error: input file not found: {input_file}", file=sys.stderr)
        return 2
    if args.np < 1:
        print("error: np must be a positive integer", file=sys.stderr)
        return 2

    script_dir = Path(__file__).resolve().parent
    try:
        executable = executable_path(script_dir)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    mpirun = os.environ.get("MPIEXEC") or shutil.which("mpirun")
    if not mpirun:
        print("error: mpirun not found; set MPIEXEC to the MPI launcher", file=sys.stderr)
        return 2

    environment = os.environ.copy()
    environment["WAVEQLAB3D_PREFLIGHT_ONLY"] = "1"
    command = [mpirun, "-np", str(args.np), str(executable), str(input_file)]

    print(f"Input:      {input_file}")
    print(f"Processes:  {args.np}")
    print(f"Executable: {executable}")
    print(flush=True)
    try:
        return subprocess.run(command, env=environment, check=False).returncode
    except OSError as exc:
        print(f"error: failed to launch MPI preflight: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
