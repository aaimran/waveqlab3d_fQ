# WaveQLab3D Response Types

## Response Summary

| `response` | Namelist | Mechanisms | Blocks | Physics |
|---|---|---|---|---|
| `elastic` | — | 0 | 1–2 | Pure elastic wave propagation, no attenuation |
| `plastic` | — | 0 | 1–2 | Elastoplastic (Drucker-Prager yield) |
| `anelastic` | `&anelastic_list` | 3 | 1–2 | Legacy low-pass viscoelastic (deprecated) |
| `low-pass` | `&anelastic_list` | 3 | 1–2 | Alias for `anelastic` |
| `anelastic-Q4` | `&anelastic_Q4_list` | 4 | 1–2 | Constant-Q attenuation, 4 relaxation mechanisms |
| `anelastic-Q8` | `&anelastic_Q8_list` | 8 | 1–2 | Constant-Q attenuation, 8 relaxation mechanisms |
| `anelastic-cQ8-b2` | `&anelastic_cQ8_b2_list` | 8 | 2 only | Per-block constant-Q, 8 mechanisms |
| `anelastic-cQ` | `&anelastic_cQ_list` | 4–8 | 2 only | Per-block constant-Q, variable mechanisms, NNLS fitting |
| `anelastic-Qf` | `&anelastic_Qf_list` | 4 | 1–2 | Frequency-dependent Q(f) = Q₀·f^γ, 4 mechanisms |
| `constant-Q-4M` | `&constant_Q_4M_list` | 4 | 1–2 | Constant-Q with manual weight control, 4 mechanisms |
| `constant-Q-8M` | `&constant_Q_8M_list` | 8 | 1–2 | Constant-Q with manual weight control, 8 mechanisms |
| `anelastic-fQ8` | `&anelastic_fQ8_list` | 8 | 1–2 | Frequency-dependent Q(f) = Q₀·f^γ, 8 mechanisms (Withers 2015) |

### Deprecated Aliases

| Input | Normalized to |
|---|---|
| `anelastic-Q` | `anelastic-Q4` |
| `frequency-Q-8M` | `anelastic-fQ8` |
| `frequency-Q-4M` | `anelastic-Qf` (same init routine) |

---

## Parameter Reference by Namelist

### `&anelastic_list` (responses: `anelastic`, `low-pass`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `c` | real | — | Direct relaxation coefficient (legacy, not Q-based) |
| `weight_exp` | real | — | Weighting exponent for mechanism distribution |
| `fref` | real | — | Reference frequency (Hz) where unrelaxed modulus is defined |

**Deprecated.** Use `anelastic-Q4` or `anelastic-Q8` instead. Specifies attenuation through direct coefficient `c` rather than target Q values.

---

### `&anelastic_Q4_list` (response: `anelastic-Q4`)

| Parameter | Type | Default | Options | Description |
|---|---|---|---|---|
| `weight_method` | string | `'fixed-q50'` | `'fixed-q50'` | Coefficient fitting method. `fixed-q50`: pre-computed weights optimized for Q≈50 over standard band |
| `Qs0` | real | — (required) | > 0 | Target shear-wave quality factor Q_S |
| `Qp0` | real | — (required) | > 0 | Target compressional-wave quality factor Q_P |
| `fref` | real | 1.0 | > 0 | Reference frequency (Hz). Elastic moduli are defined at this frequency |
| `fmin` | real | 0.05 | > 0 | Lower bound of fitting band (Hz). Must be 0.05 for `fixed-q50` |
| `fmax` | real | 20.0 | > fmin | Upper bound of fitting band (Hz). Must be 20.0 for `fixed-q50` |

4 relaxation mechanisms. Constant Q(f) ≈ Q₀ over [fmin, fmax]. Weights pre-computed for Q=50 and rescaled to target Q.

---

### `&anelastic_Q8_list` (response: `anelastic-Q8`)

| Parameter | Type | Default | Options | Description |
|---|---|---|---|---|
| `weight_method` | string | `'fixed-q50'` | `'fixed-q50'` | Coefficient fitting method. Same as Q4 but with 8 mechanisms |
| `Qs0` | real | — (required) | > 0 | Target Q_S |
| `Qp0` | real | — (required) | > 0 | Target Q_P |
| `fref` | real | 1.0 | > 0 | Reference frequency (Hz) |
| `fmin` | real | 0.05 | > 0 | Lower fitting band (Hz). Must be 0.05 for `fixed-q50` |
| `fmax` | real | 20.0 | > fmin | Upper fitting band (Hz). Must be 20.0 for `fixed-q50` |

8 relaxation mechanisms. Better Q-approximation flatness than Q4 across the frequency band.

---

### `&anelastic_cQ8_b2_list` (response: `anelastic-cQ8-b2`)

Requires `nblocks = 2`.

| Parameter | Type | Default | Options | Description |
|---|---|---|---|---|
| `weight_method` | string | `'fixed-q50'` | `'fixed-q50'` | Coefficient fitting method |
| `Qs0(2)` | real(2) | — (required) | > 0 each | Target Q_S per block: `Qs0(1)` for block 1, `Qs0(2)` for block 2 |
| `Qp0(2)` | real(2) | — (required) | > 0 each | Target Q_P per block |
| `fref` | real | 1.0 | > 0 | Reference frequency (Hz) |
| `fmin` | real | 0.05 | > 0 | Lower fitting band (Hz) |
| `fmax` | real | 20.0 | > fmin | Upper fitting band (Hz) |

Two-block contrast model with independent Q in each block. 8 mechanisms per block. For layered models with different attenuation in each layer (e.g., sediment over bedrock).

---

### `&anelastic_cQ_list` (response: `anelastic-cQ`)

Requires `nblocks = 2`.

| Parameter | Type | Default | Options | Description |
|---|---|---|---|---|
| `Qs0(2)` | real(2) | — (required) | > 0 each | Target Q_S per block |
| `Qp0(2)` | real(2) | — (required) | > 0 each | Target Q_P per block |
| `fref` | real | 1.0 | > 0 | Reference frequency (Hz) |
| `fmin` | real | 0.05 | > 0 | Lower fitting band (Hz) |
| `fmax` | real | 20.0 | > fmin | Upper fitting band (Hz) |
| `n_mechanisms` | integer | 8 | 4, 5, 6, 7, 8 | Number of relaxation mechanisms |
| `nnls_samples` | integer | 256 | ≥ n_mechanisms | Number of frequency samples for NNLS fitting |
| `coefficient_policy` | string | `'nnls-block-ps'` | `'fixed-q50'`, `'nnls-shared'`, `'nnls-block'`, `'nnls-block-ps'` | How relaxation weights are computed (see below) |
| `nnls_objective` | string | `'relative-q'` | `'relative-q'` | NNLS objective function |
| `nnls_tolerance` | real | 1e-10 | > 0 | NNLS convergence tolerance |
| `max_fit_error` | real | 0.10 | > 0 | Maximum allowed relative Q fitting error (fraction) |

**Coefficient policies:**
- `fixed-q50`: Pre-computed weights for Q=50, requires n_mechanisms=8, fmin=0.05, fmax=20.0
- `nnls-shared`: NNLS-fitted weights shared between blocks
- `nnls-block`: Independent NNLS fit per block
- `nnls-block-ps`: Independent NNLS fit per block, separate P and S weights

Most flexible constant-Q model. Per-block Q with tunable mechanism count and fitting strategy.

---

### `&anelastic_Qf_list` (responses: `anelastic-Qf`, `frequency-Q-4M`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `c` | real | — | Direct relaxation coefficient |
| `gamma` | real | — | Frequency exponent: Q(f) = Q₀·f^γ |
| `f_trans` | real | — | Transition frequency (Hz) |
| `fref` | real | — | Reference frequency (Hz) |

4-mechanism frequency-dependent Q. Power-law Q(f) model.

---

### `&constant_Q_4M_list` (response: `constant-Q-4M`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `c` | real | — | Direct relaxation coefficient |
| `fmin` | real | — | Lower fitting band (Hz) |
| `fmax` | real | — | Upper fitting band (Hz) |
| `target_Q` | real | — | Target constant Q value |
| `weight_method` | string | — | Weight computation method |
| `manual_weights` | real array | — | User-specified mechanism weights (optional) |
| `fref` | real | — | Reference frequency (Hz) |

4-mechanism constant-Q with explicit weight control. Allows manual weight specification.

---

### `&constant_Q_8M_list` (response: `constant-Q-8M`)

Same parameters as `&constant_Q_4M_list` but with 8 mechanisms.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `c` | real | — | Direct relaxation coefficient |
| `fmin` | real | — | Lower fitting band (Hz) |
| `fmax` | real | — | Upper fitting band (Hz) |
| `target_Q` | real | — | Target constant Q value |
| `weight_method` | string | — | Weight computation method |
| `manual_weights` | real array | — | User-specified mechanism weights (optional) |
| `fref` | real | — | Reference frequency (Hz) |

---

### `&anelastic_fQ8_list` (response: `anelastic-fQ8`)

| Parameter | Type | Default | Options | Description |
|---|---|---|---|---|
| `coefficient_method` | string | `'conventional-nnls'` | `'conventional-nnls'`, `'withers-2015'` | How relaxation times and weights are computed (see below) |
| `weight_policy` | string | `'table-exact'` | `'table-exact'`, `'nonnegative-refit'` | Post-processing of weights. `nonnegative-refit` only valid with `withers-2015` |
| `coarse_grain` | integer | -1 | 0, 2 | Coarse-graining level. 0: no coarse-graining (use with `conventional-nnls`). 2: coarse-grain 8→4 mechanism pairs (use with `withers-2015`) |
| `Qs0` | real | — (required) | ≥ 15 | Target shear-wave quality factor Q_S |
| `Qp0` | real | — (required) | ≥ 15 | Target compressional-wave quality factor Q_P |
| `gamma` | real | — (required) | ≥ 0 | Frequency exponent: Q(f) = Q₀·f^γ. 0 = constant Q. Typical: 0.0–1.0 |
| `f_transition` | real | — (required) | > 0 | Transition frequency (Hz) separating low/high frequency behavior |
| `fref` | real | — (required) | > 0 | Reference frequency (Hz). Elastic moduli defined here |

**Coefficient methods:**
- `conventional-nnls`: Standard NNLS optimization. Log-spaced relaxation times, NNLS-fit weights to match target Q(f). Use with `coarse_grain=0`.
- `withers-2015`: Withers, Olsen & Day (2015, BSSA) Table 1 coefficients. Pre-tabulated relaxation times and strength ratios for accurate power-law Q(f). Use with `coarse_grain=2`.

**Weight policies:**
- `table-exact`: Use tabulated weights directly (default).
- `nonnegative-refit`: Refit negative weights to nearest non-negative solution. Only with `withers-2015`.

8-mechanism frequency-dependent Q model. The primary method for Withers et al. (2015) benchmarks. Supports both constant Q (γ=0) and power-law frequency-dependent Q (γ>0).

---

## Typical Configurations

### Elastic (no attenuation)
```fortran
&problem_list response = 'elastic' /
```

### Constant Q = 50 (8 mechanisms)
```fortran
&problem_list response = 'anelastic-Q8' /
&anelastic_Q8_list Qs0 = 50d0, Qp0 = 100d0, fref = 1d0 /
```

### Power-law Q(f) = 50·f^0.6 (Withers 2015)
```fortran
&problem_list response = 'anelastic-fQ8' /
&anelastic_fQ8_list
 coefficient_method = 'withers-2015', weight_policy = 'table-exact', coarse_grain = 2,
 Qs0 = 50d0, Qp0 = 100d0, gamma = 0.6d0, f_transition = 1d0, fref = 1d0
/
```

### Layered model — different Q per block
```fortran
&problem_list response = 'anelastic-cQ8-b2', nblocks = 2 /
&anelastic_cQ8_b2_list
 Qs0 = 20d0, 210d0, Qp0 = 40d0, 420d0, fref = 1d0
/
```

### Constant Q via NNLS fit (conventional)
```fortran
&problem_list response = 'anelastic-fQ8' /
&anelastic_fQ8_list
 coefficient_method = 'conventional-nnls', weight_policy = 'table-exact', coarse_grain = 0,
 Qs0 = 50d0, Qp0 = 100d0, gamma = 0d0, f_transition = 1d0, fref = 1d0
/
```

---

## Source Files

| Response | Init routine | Source file |
|---|---|---|
| `elastic` | (none — elastic is baseline) | `src/block.f90` |
| `plastic` | `init_plastic_material` | `src/block.f90` |
| `anelastic` / `low-pass` | `init_anelastic_properties` | `src/material.f90:22` |
| `anelastic-Q4` | `init_anelastic_Q4_properties` | `src/anelastic_q4_model.f90` |
| `anelastic-Q8` | `init_anelastic_Q8_properties` | `src/anelastic_q8_model.f90` |
| `anelastic-cQ8-b2` | `init_anelastic_Q8_properties` (per block) | `src/anelastic_cq8_b2_model.f90` |
| `anelastic-cQ` | `init_anelastic_cq_properties` | `src/anelastic_cq_model.f90` |
| `anelastic-Qf` / `frequency-Q-4M` | `init_anelastic_Qf_properties` | `src/material.f90:517` |
| `constant-Q-4M` | `init_const_Q_4M_properties` | `src/material.f90:1265` |
| `constant-Q-8M` | `init_const_Q_8M_properties` | `src/material.f90:1688` |
| `anelastic-fQ8` | `init_anelastic_Qf8_properties` | `src/anelastic_fq8_model.f90` |

## RHS Dispatch Flags

Each response sets a boolean flag on the material type `M`. The RHS kernel (`src/RHS_Interior.f90`) checks these flags at each grid point:

| Flag | Responses | Dispatch routine |
|---|---|---|
| `M%anelastic` | `anelastic`, `low-pass` | `apply_anelastic_point_dispatch` |
| `M%anelastic_Q` | `anelastic-Q4` | `apply_anelastic_Q_point_dispatch` |
| `M%anelastic_Q8` | `anelastic-Q8`, `anelastic-cQ8-b2` | `apply_anelastic_Q8_point_dispatch` |
| `M%anelastic_cQ` | `anelastic-cQ` | `apply_anelastic_cQ_point_dispatch` |
| `M%anelastic_Qf` | `anelastic-Qf`, `frequency-Q-4M` | `apply_anelastic_Qf_point_dispatch` |
| `M%anelastic_const_Q_4M` | `constant-Q-4M` | `apply_const_Q_4M_point_dispatch` |
| `M%anelastic_const_Q_8M` | `constant-Q-8M` | (8M variant) |
| `M%anelastic_Qf8` | `anelastic-fQ8` | `apply_anelastic_Qf_point_dispatch` (collocated) |
