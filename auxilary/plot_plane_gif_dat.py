#!/usr/bin/env python3
"""
Script to parse seismogram .dat files and generate .gif animations for vx, vy, vz time series on the xz plane at y=0.

Default domain settings (can be modified):
- x = [0, 30] km, y = 0 km, z = [0, 30] km
- Block1: x=[0, 15] km, nx=76
- Block2: x=[15, 30] km, nx=76
- z=[0, 30] km, nz=151
- Stations at 1 km resolution on the plane

Usage: python plot_plane_gif.py [options] <dat_folder>

Options:
  -np <num_processes>  Number of parallel processes (default: 1)
  -fs <frame_stride>   Stride for frames (default: 1, use every frame)
  -o <output_gif>      Output gif filename (default: plane_animation.gif)
  -h                   Show help

Console output includes:
- Total .dat files
- .dat files per block
- Common points between blocks
- Total unique points
- Progress bar during parsing
"""

import os
import sys
import glob
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from multiprocessing import Pool
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

# Default domain parameters (can be changed)
X_MIN, X_MAX = 0.0, 30.0  # km
Y_FIXED = 0.0  # km
Z_MIN, Z_MAX = 0.0, 30.0  # km
BLOCK1_X_MIN, BLOCK1_X_MAX = 0.0, 15.0  # km
BLOCK2_X_MIN, BLOCK2_X_MAX = 15.0, 30.0  # km
NX_BLOCK1 = 76
NX_BLOCK2 = 76
NZ = 151
STATION_RES = 1.0  # km

# Grid for plotting (31x31 for 0 to 30 at 1 km)
GRID_X = np.arange(X_MIN, X_MAX + STATION_RES, STATION_RES)
GRID_Z = np.arange(Z_MIN, Z_MAX + STATION_RES, STATION_RES)
NX_GRID = len(GRID_X)
NZ_GRID = len(GRID_Z)

def parse_filename(filename):
    """Parse filename to extract block_num, i, j, k (logical indices)"""
    basename = os.path.basename(filename)
    match = re.match(r'.*_(\d+)_(\d+)_(\d+)_block(\d+)\.dat', basename)
    if not match:
        raise ValueError(f"Cannot parse filename: {basename}")
    i, j, k, block_num = map(int, match.groups())
    return block_num, i, j, k

def logical_to_physical(block_num, i, j, k):
    """Convert logical indices to physical coordinates (km)"""
    if block_num == 1:
        x_min, x_max, nx = BLOCK1_X_MIN, BLOCK1_X_MAX, NX_BLOCK1
        y_min, y_max, ny = 0.0, 20.0, 101  # from input
        z_min, z_max, nz = Z_MIN, Z_MAX, NZ
    elif block_num == 2:
        x_min, x_max, nx = BLOCK2_X_MIN, BLOCK2_X_MAX, NX_BLOCK2
        y_min, y_max, ny = 0.0, 20.0, 101
        z_min, z_max, nz = Z_MIN, Z_MAX, NZ
    else:
        raise ValueError(f"Unknown block_num: {block_num}")

    # Assuming uniform spacing, indices start from 1 (interior)
    hx = (x_max - x_min) / (nx - 1)
    hy = (y_max - y_min) / (ny - 1)
    hz = (z_max - z_min) / (nz - 1)

    x = x_min + (i - 1) * hx
    y = y_min + (j - 1) * hy
    z = z_min + (k - 1) * hz

    return x, y, z

def parse_dat_file(filepath):
    """Parse a single .dat file, return (block_num, x, z, times, vx, vy, vz)"""
    block_num, i, j, k = parse_filename(filepath)
    x, y, z = logical_to_physical(block_num, i, j, k)

    # Read data
    data = np.loadtxt(filepath)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    times = data[:, 0]
    vx = data[:, 1]
    vy = data[:, 2]
    vz = data[:, 3]

    return block_num, x, z, times, vx, vy, vz

def process_files_parallel(filepaths, num_processes):
    """Process files in parallel, return list of results"""
    with Pool(num_processes) as pool:
        results = list(tqdm(pool.imap(parse_dat_file, filepaths), total=len(filepaths), desc="Parsing files"))
    return results

def main():
    parser = argparse.ArgumentParser(description="Generate .gif from seismogram .dat files")
    parser.add_argument('dat_folder', help='Folder containing .dat files')
    parser.add_argument('-np', '--num_processes', type=int, default=1, help='Number of parallel processes')
    parser.add_argument('-fs', '--frame_stride', type=int, default=1, help='Stride for frames (1 = every frame)')
    parser.add_argument('-o', '--output', default='plane_animation.gif', help='Output gif filename')
    args = parser.parse_args()

    dat_folder = args.dat_folder
    num_processes = args.num_processes
    frame_stride = args.frame_stride
    output_gif = args.output

    # Find all .dat files
    dat_files = glob.glob(os.path.join(dat_folder, '*.dat'))
    total_files = len(dat_files)
    print(f"Total .dat files: {total_files}")

    # Group by block
    block1_files = [f for f in dat_files if '_block1.dat' in f]
    block2_files = [f for f in dat_files if '_block2.dat' in f]
    print(f".dat files in block1: {len(block1_files)}")
    print(f".dat files in block2: {len(block2_files)}")

    # Parse files
    print("Parsing block1 files...")
    block1_data = process_files_parallel(block1_files, num_processes)
    print("Parsing block2 files...")
    block2_data = process_files_parallel(block2_files, num_processes)

    # Collect unique points
    points_block1 = {(x, z) for _, x, z, _, _, _, _ in block1_data}
    points_block2 = {(x, z) for _, x, z, _, _, _, _ in block2_data}
    common_points = points_block1 & points_block2
    unique_points = points_block1 | points_block2
    print(f"Common points between blocks: {len(common_points)}")
    print(f"Total unique points: {len(unique_points)}")

    # Prepare data structures
    # Assume all files have same time steps
    if block1_data:
        times = block1_data[0][3]
    elif block2_data:
        times = block2_data[0][3]
    else:
        raise ValueError("No data found")

    n_times = len(times)
    vx_grid = np.full((n_times, NX_GRID, NZ_GRID), np.nan)
    vy_grid = np.full((n_times, NX_GRID, NZ_GRID), np.nan)
    vz_grid = np.full((n_times, NX_GRID, NZ_GRID), np.nan)

    # Fill grids, preferring block1 for common points
    for block_data in [block1_data, block2_data]:
        for block_num, x, z, _, vx, vy, vz in block_data:
            if block_num == 2 and (x, z) in points_block1:
                continue  # Skip block2 common points
            ix = int(round((x - X_MIN) / STATION_RES))
            iz = int(round((z - Z_MIN) / STATION_RES))
            if 0 <= ix < NX_GRID and 0 <= iz < NZ_GRID:
                vx_grid[:, ix, iz] = vx
                vy_grid[:, ix, iz] = vy
                vz_grid[:, ix, iz] = vz

    # Compute global min-max for fixed colorbars
    vmin_vx, vmax_vx = np.nanmin(vx_grid), np.nanmax(vx_grid)
    vmin_vy, vmax_vy = np.nanmin(vy_grid), np.nanmax(vy_grid)
    vmin_vz, vmax_vz = np.nanmin(vz_grid), np.nanmax(vz_grid)

    # Create symmetric log norms for better contrast
    norm_vx = mcolors.SymLogNorm(linthresh=1e-6, vmin=vmin_vx, vmax=vmax_vx)
    norm_vy = mcolors.SymLogNorm(linthresh=1e-6, vmin=vmin_vy, vmax=vmax_vy)
    norm_vz = mcolors.SymLogNorm(linthresh=1e-6, vmin=vmin_vz, vmax=vmax_vz)

    # Generate frames
    frames = []
    fig, axes = plt.subplots(3, 1, figsize=(10, 25))
    # Create dummy images and colorbars to avoid accumulation
    ims = []
    cbs = []
    for ax, norm in zip(axes, [norm_vx, norm_vy, norm_vz]):
        im = ax.imshow(np.zeros((NZ_GRID, NX_GRID)), extent=[X_MIN, X_MAX, Z_MIN, Z_MAX], origin='lower', cmap='seismic', aspect='equal', norm=norm)
        ims.append(im)
        cb = plt.colorbar(im, ax=ax, shrink=0.8)
        cbs.append(cb)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('z (km)')
    # No need for set_clim since norm handles it
    for t_idx in tqdm(range(0, n_times, frame_stride), desc="Generating frames"):
        vx_data = vx_grid[t_idx]
        vy_data = vy_grid[t_idx]
        vz_data = vz_grid[t_idx]
        for im, data, title in zip(ims, [vx_data, vy_data, vz_data], ['vx', 'vy', 'vz']):
            im.set_data(data.T)
            im.axes.set_title(f'{title} at t={times[t_idx]:.2f} s')

        fig.suptitle(f'Plane at y=0 km, t={times[t_idx]:.2f} s')
        fig.canvas.draw()
        buf = fig.canvas.buffer_rgba()
        image = np.frombuffer(buf, dtype=np.uint8).reshape(fig.canvas.get_width_height()[::-1] + (4,))
        image = image[:, :, :3]  # RGB, drop alpha
        if HAS_PIL:
            frames.append(Image.fromarray(image))
        else:
            frames.append(image)

    plt.close(fig)

    # Save gif
    if HAS_PIL:
        if frames:
            frames[0].save(output_gif, save_all=True, append_images=frames[1:], duration=100, loop=0)
        print(f"Saved animation to {output_gif}")
    else:
        print("PIL not available, cannot create GIF. Frames saved as PNGs in /tmp/frames/")
        os.makedirs('/tmp/frames', exist_ok=True)
        for i, frame in enumerate(frames):
            Image.fromarray(frame).save(f'/tmp/frames/frame_{i:04d}.png')
        print("Run: ffmpeg -i /tmp/frames/frame_%04d.png -vf 'fps=10' output.gif")

if __name__ == '__main__':
    main()