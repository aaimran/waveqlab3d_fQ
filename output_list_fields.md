# `&output_list` Namelist Fields

| Field | Type | Default | Purpose |
|---|---|---|---|
| `output_exact_moment` | logical | F | Write exact moment tensor info |
| `output_seismograms` | logical | F | Write station seismograms |
| `output_station_info` | logical | F | Write station metadata |
| `output_station_mapping` | logical | F | Write station-to-grid mapping |
| `output_fault_topo` | logical | F | Write fault/topography data |
| `output_fields_block1` | logical | F | Write full 3D field snapshots (block 1) |
| `output_fields_block2` | logical | F | Write full 3D field snapshots (block 2) |
| `stride_fields` | integer | 1 | Temporal stride for field output |
| `station_xyz_index` | logical | F | Use xyz coordinates vs grid indices for stations |
| `station_list` | string | `'infile'` | Station source: `'infile'` (from input file) or `'extfile'` (external file) |
| `station_list_file` | string | `''` | External station file path (required when `station_list='extfile'`) |
| `station_file_directory` | string | `'seismogram'` | Output directory for seismogram files |
| `station_output_order` | string | | Column ordering in station output |
| `station_number_in_list` | logical | F | Station rows include a station number column |
| `station_number_in_filename` | logical | F | Include station number in output filenames |
| `station_use_block_subdirectories` | logical | F | Separate output dirs per block |
| `interface_stations` | string | `'block_1_2'` | Which block(s) output interface stations: `block_1`, `block_2`, `block_1_2` |
| `append_block` | logical | F | Append `_blockN` to output filenames. Forced T when `interface_stations='block_1_2'` |
| `station_add_header` | logical | F | Add header line to seismogram files |
| `station_add_metadata` | logical | F | Add metadata to seismogram files |

## Station List Format

When `output_seismograms=T`, stations are defined between markers in the input file (or external file):

```
!---begin:station_listU---
x1  y1  z1
x2  y2  z2
!---end:station_listU---
```

If `station_number_in_list=T`, each row has 4 fields: `station_number x y z`.

For 2-block setups, use `station_listU` (block 1) and `station_listV` (block 2).
