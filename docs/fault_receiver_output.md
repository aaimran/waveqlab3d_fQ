# Pointwise fault-receiver output

`waveqlab3d_Q` can write benchmark time histories at selected fault receivers
without writing full fault-plane arrays. The receiver output is independent of
`w_fault`; set `w_fault = F` when full-plane binary output is not needed.

Add this namelist and receiver list to the input file:

```fortran
&fault_receiver_output_list
 enable_fault_receivers = T,
 number_fault_receivers = 2,
 fault_receiver_stride = 1,
 fault_receiver_coordinate_mode = 'strike_down_dip',
 fault_receiver_directory = 'output/fault_receivers',
 fault_dip_degrees = 15.0 /

!---begin:fault_receiver_list---
 fault_receiver_1   0.0   7.5
 fault_receiver_2   3.0  12.0
!---end:fault_receiver_list---
```

For `strike_down_dip`, each record is:

```text
name  strike_km  down_dip_km
```

For a planar fault, the target Cartesian point is computed as
`(x,y,z) = (down_dip*cos(dip), down_dip*sin(dip), strike)`. Alternatively,
set `fault_receiver_coordinate_mode = 'xyz'`; each record then has
`name x_km y_km z_km`.

Each target is mapped once at startup to the nearest physical node on the
block-1 side of the fault. The file header records the requested coordinate,
mapped coordinate, and mapping distance in km. A large mapping distance means
the receiver or mesh should be checked. This first implementation deliberately
uses nearest-node sampling; it does not interpolate.

One benchmark-compatible text file is written per receiver in
`fault_receiver_directory`. Its un-commented field-list line and eight columns
are:

```text
t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate  v-shear-stress  n-stress
```

The units are `s, m, m/s, MPa, m, m/s, MPa, MPa`. Horizontal is the local
along-strike (`l`) component and vertical is the local down-dip (`m`) component.
Normal stress follows the solver sign convention, so compression is negative.
`fault_receiver_stride = N` writes every Nth full time step, beginning at
`t = 0`.
