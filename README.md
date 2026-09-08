# README #

WaveQLab3D is a code for 3D seismic wave propagation and earthquake rupture dynamics. It solves the elastic wave equation in curvilinear coordinates (i.e., complex geometries) with a possibly nonplanar frictional fault interface. The current version supports off-fault viscoplasticity, spatially variable elastic properties, and several friction laws (including rate-and-state and slip-weakening). The code is under development and is available under the MIT license. Authors include Kenneth Duru, Sam Bydlon, Eric Dunham, and Kyle Withers with parallelization by Hari Radhakrishnan.

Supported attenuation response options currently include `anelastic`, `anelastic-Q`, `anelastic-Q8`, `anelastic-cQ8-b2`, `anelastic-cQ`, `anelastic-Qf`, `constant-Q-4M`, `constant-Q-8M`, `frequency-Q-4M`, and `frequency-Q-8M`.

For the fixed eight-mechanism constant-Q response, prefer explicit P- and
S-wave quality factors in anelastic_Q8_list:

    &problem_list
      response = 'anelastic-Q8'
    /

    &anelastic_Q8_list
      Qs0  = 50.0
      Qp0  = 50.0
      fref = 1.0
    /

Qs0 and Qp0 must be supplied together and must be positive. The stored eight-mechanism weights are normalized
spectral-shape coefficients; the RHS scales them once by the local inverse Q.

For a two-block model with a different constant Q pair in each block, use the
additive `anelastic-cQ8-b2` response. It requires `nblocks=2`; array element 1
applies to block 1 and element 2 applies to block 2:

    &problem_list
      response = 'anelastic-cQ8-b2'
      nblocks = 2
    /

    &anelastic_cQ8_b2_list
      Qs0 = 40.0, 100.0
      Qp0 = 80.0, 180.0
      fref = 1.0
      fmin = 0.05
      fmax = 20.0
      weight_method = 'fixed-q50'
    /

The reference frequency, approximation band, and eight-mechanism weight method
are shared by the blocks. Existing `anelastic-Q8` inputs and behavior are unchanged.

The independent fitted constant-Q response is selected with `anelastic-cQ`.
It requires two blocks, accepts 4 through 8 relaxation mechanisms, and keeps
its configuration and runtime memory separate from the fixed Q8 responses:

    &problem_list
      response = 'anelastic-cQ'
      nblocks = 2
    /

    &anelastic_cQ_list
      Qs0 = 40.0, 100.0
      Qp0 = 80.0, 180.0
      fref = 1.0
      fmin = 0.05
      fmax = 20.0
      n_mechanisms = 6
      coefficient_policy = 'nnls-block-ps'
      nnls_samples = 256
      nnls_objective = 'relative-q'
      nnls_tolerance = 1.0e-10
      max_fit_error = 0.10
    /

Supported coefficient policies are `nnls-shared`, `nnls-block`,
`nnls-block-ps`, and `fixed-q50`. The fixed table requires eight mechanisms;
the NNLS policies support 4, 5, 6, 7, or 8. `max_fit_error` is the maximum
allowed relative error over the requested frequency band and should be chosen
to match the accuracy required by the simulation.
The response supports `fd_type='traditional'`, `fd_type='upwind'`, and
`fd_type='upwind_drp'` with the order combinations accepted by preflight.

Station output columns default to `t vx vy vz`. Their order can be changed in
the `output_list` namelist; for example:

    &output_list
      output_seismograms = T,
      station_output_order = 't vz vx vy'
    /

`station_output_order` is case-insensitive, accepts spaces or commas between
names, and must contain each of `t`, `vx`, `vy`, and `vz` exactly once.

A leading integer station number can optionally be included on every station
list row and used in the output filename:

    &output_list
      output_seismograms = T,
      station_number_in_list = T,
      station_number_in_filename = T
    /

    !---begin:station_list---
    1   0.693d0   0.000d0   0.000d0
    2   5.543d0   0.000d0   0.000d0
    3  10.392d0   0.000d0   0.000d0
    !---end:station_list---

This produces names such as `fname_station-1.dat`. Both options default to
false. `station_number_in_filename = T` requires
`station_number_in_list = T`.

Station files can be written directly in `station_file_directory` instead of
its `block1` and `block2` subdirectories, and stations can be restricted to one
block:

    &output_list
      output_seismograms = T,
      station_use_block_subdirectories = F,
      common_stations_blocks = 'both'
    /

`common_stations_blocks` accepts `block1`, `block2`, or `both` and defaults to `both`.
`station_use_block_subdirectories` defaults to true. When `both` is selected,
station files that use physical coordinates or station numbers in their names
receive a `_block1` or `_block2` suffix, preventing common-plane stations from
overwriting one another.

Optional commented headers and station metadata can be written at the start of
each station `.dat` file:

    &output_list
      station_add_header = T,
      station_add_metadata = T
    /

For `station_output_order = 't vz vx vy'`, the beginning of a numbered station
file is:

    # station_number: 1
    # x y z:  6.9300000000000000E-001  0.0000000000000000E+000  0.0000000000000000E+000
    # grid_i j k: 58 1 51
    # grid_x y z:  7.0000000000000000E-001  0.0000000000000000E+000  0.0000000000000000E+000
    # mapping_distance:  7.0000000000000000E-003
    # t vz vx vy

Requested and mapped physical coordinates, grid indices, and mapping distance
are included when metadata is enabled. The station number line is included when
`station_number_in_list = T`. Both preamble options default to false.

The station-to-grid mapping printed during startup is controlled separately
from the boxed station summary:

    &output_list
      output_station_mapping = F
    /

It defaults to true. When enabled, every matched station is printed on one
clearly delimited line ending in a semicolon, for example:

    station 1: distance= 7.000000E-003, indices=(58 1 51), grid_xyz=(...), requested_xyz=(...);

`output_station_info` only controls the boxed configuration summary;
`output_station_mapping` controls these individual station mapping lines.
