#!/usr/bin/env python3
"""
Script to check dimensions and fields in .plane files.
"""

import os
import glob
import numpy as np


def parse_plane_file(filepath):
    """Parse a .plane file, return nx, nz, list of (time, vx, vy, vz)"""
    results = []
    with open(filepath, 'rb') as f:
        # Read metadata
        version = np.frombuffer(f.read(4), dtype='<i4')[0]
        axis = np.frombuffer(f.read(4), dtype='<i4')[0]
        nx = np.frombuffer(f.read(4), dtype='<i4')[0]
        nz = np.frombuffer(f.read(4), dtype='<i4')[0]
        fixed = np.frombuffer(f.read(4), dtype='<i4')[0]
        # Now read data records
        try:
            while True:
                time = np.frombuffer(f.read(8), dtype='<f8')[0]
                vx = np.frombuffer(f.read(nx * nz * 8), dtype='<f8').reshape((nx, nz))
                vy = np.frombuffer(f.read(nx * nz * 8), dtype='<f8').reshape((nx, nz))
                vz = np.frombuffer(f.read(nx * nz * 8), dtype='<f8').reshape((nx, nz))
                results.append((time, vx, vy, vz))
        except Exception:
            pass  # EOF or partial read
    return nx, nz, results


def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: python check_plane.py <plane_folder>")
        sys.exit(1)

    plane_folder = sys.argv[1]
    plane_files = glob.glob(os.path.join(plane_folder, '*.plane'))

    print(f"Found {len(plane_files)} .plane files")

    for filepath in plane_files:
        basename = os.path.basename(filepath)
        nx, nz, records = parse_plane_file(filepath)
        nt = len(records)
        print(f"\nFile: {basename}")
        print(f"  nx: {nx}, nz: {nz}, nt: {nt}")
        if records:
            times = [r[0] for r in records]
            print(f"  Time range: {min(times):.6f} to {max(times):.6f}")
            vx_all = np.concatenate([r[1].flatten() for r in records])
            vy_all = np.concatenate([r[2].flatten() for r in records])
            vz_all = np.concatenate([r[3].flatten() for r in records])
            print(f"  vx range: {np.nanmin(vx_all):.6e} to {np.nanmax(vx_all):.6e}")
            print(f"  vy range: {np.nanmin(vy_all):.6e} to {np.nanmax(vy_all):.6e}")
            print(f"  vz range: {np.nanmin(vz_all):.6e} to {np.nanmax(vz_all):.6e}")


if __name__ == '__main__':
    main()
