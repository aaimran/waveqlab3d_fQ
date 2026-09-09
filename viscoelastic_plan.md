# Plan: Unified `viscoelastic` Response

## Motivation

The current codebase has **7 separate anelastic response types** with nearly identical physics but different:
- Memory variable arrays (`eta4Q`, `eta4Q8`, `eta4cQ`, `eta4Qf`, `eta4Qf8`, ...)
- RHS kernels (122 references in `RHS_Interior.f90` alone, 3578 lines)
- Time-stepping blocks (6 duplicated blocks in `fields.f90`)
- Parameter types (`q4_parameters`, `q8_parameters`, `cq_parameters`, `fq8_parameters`, ...)
- Init routines across 5+ source files

All implement the **same generalized standard linear solid (SLS)** / Maxwell body:
```
dη_k/dt = (σ_target - η_k) / τ_k     for k = 1..N
σ_total = σ_elastic - Σ η_k
```

The only differences: number of mechanisms N, how τ_k and weights are computed, and whether Q varies per block.

## Design: `response = 'viscoelastic'`

### Namelist: `&viscoelastic_list`

```fortran
&viscoelastic_list
 attenuation = 'constant-Q',      ! or 'frequency-Q'
 Qs0 = 50d0,                      ! scalar or array(2) for per-block
 Qp0 = 100d0,                     ! scalar or array(2) for per-block
 gamma = 0d0,                     ! power-law exponent (frequency-Q only)
 fmin = 0.05d0,                   ! lower fitting band (Hz)
 fmax = 20d0,                     ! upper fitting band (Hz)
 fref = 1d0,                      ! reference frequency (Hz)
 f_transition = 1d0,              ! transition frequency (frequency-Q only)
 n_mechanisms = 8,                ! 3, 4, 5, 6, 7, or 8
 weight_method = 'nnls',          ! coefficient fitting strategy
 nnls_samples = 256,              ! frequency samples for NNLS
 nnls_tolerance = 1d-10,          ! convergence tolerance
 max_fit_error = 0.10d0           ! max allowed relative Q error
/
```

### Parameter Definitions

| Parameter | Type | Default | Valid | Description |
|---|---|---|---|---|
| `attenuation` | string | `'constant-Q'` | `'constant-Q'`, `'frequency-Q'` | Attenuation model. `constant-Q`: Q(f) = Q₀ (frequency-independent). `frequency-Q`: Q(f) = Q₀·(f/f_transition)^γ for f ≥ f_transition, Q₀ below |
| `Qs0` | real or real(2) | — (required) | ≥ 15 | Target shear quality factor. Scalar: same for all blocks. Array(2): per-block [block1, block2] |
| `Qp0` | real or real(2) | — (required) | ≥ 15 | Target compressional quality factor. Same scalar/array rules as Qs0 |
| `gamma` | real | 0.0 | 0.0–1.0 | Power-law exponent. Only used when `attenuation = 'frequency-Q'`. Ignored for `constant-Q`. 0.0 = constant Q (within frequency-Q framework). Typical: 0.1–0.8 |
| `fmin` | real | 0.05 | > 0 | Lower bound of frequency band (Hz) where Q approximation is valid |
| `fmax` | real | 20.0 | > fmin | Upper bound of frequency band (Hz) |
| `fref` | real | 1.0 | fmin ≤ fref ≤ fmax | Reference frequency (Hz). Elastic moduli (Vp, Vs) are exact at fref. Physical velocity dispersion shifts phase velocity away from fref |
| `f_transition` | real | 1.0 | > 0 | Transition frequency (Hz) for `frequency-Q`. Below f_transition: Q = Q₀. Above: Q = Q₀·(f/f_transition)^γ. Ignored for `constant-Q` |
| `n_mechanisms` | integer | 8 | 3, 4, 5, 6, 7, 8 | Number of relaxation mechanisms per stress component. More = flatter Q(f) over [fmin, fmax] but more memory (6 × N 3D arrays per block) and stricter dt limit |
| `weight_method` | string | `'nnls'` | see below | How relaxation times τ_k and weights w_k are computed |
| `nnls_samples` | integer | 256 | ≥ n_mechanisms | Number of log-spaced frequency samples for NNLS fitting. Only used by `nnls` and `withers-method` |
| `nnls_tolerance` | real | 1e-10 | > 0 | NNLS convergence tolerance (max weight change per sweep) |
| `max_fit_error` | real | 0.10 | > 0 | Maximum allowed relative Q fitting error. Preflight warns if exceeded, error if > 2× |

### `weight_method` Options

| Value | Description |
|---|---|
| `'nnls'` | **General purpose.** Standard NNLS optimization. Log-spaced τ_k from 1/(2π·fmax) to 1/(2π·fmin). Fits non-negative weights to minimize Q(f) misfit over [fmin, fmax]. Works for any `attenuation`, any N, any Q, any frequency band |
| `'withers-table'` | **Paper replication.** Direct lookup from Withers, Olsen & Day (2015, BSSA) Tables 1–2. Uses tabulated τ_min/τ_max for Withers eq. (15) spacing, tabulated W_k (high Q: w_k = W_k/Q) and a_k, b_k (low Q: w_k = a_k/Q² + b_k/Q). Instant — no NNLS solve. **Constraints:** `n_mechanisms = 8`, `attenuation = 'frequency-Q'`, γ ∈ [0.0, 0.9], Q ≥ 15 |
| `'withers-method'` | **Generalized Withers.** Same algorithm the paper authors used to *generate* their tables, but run at runtime for arbitrary parameters. Uses Withers eq. (15) τ_k spacing (γ-dependent τ_min/τ_max, interpolated or extrapolated), then NNLS-fits weights with the Withers objective (Q-misfit with non-negative constraint). **No constraints on N, γ, or frequency band** — extends the methodology beyond the tabulated range |
| `'fixed-q50'` | **Legacy.** Pre-computed weights optimized for Q = 50 over [0.05, 20] Hz, linearly rescaled for other Q values. Fast (no solve). **Constraints:** `attenuation = 'constant-Q'`, `n_mechanisms ∈ {4, 8}`, fmin = 0.05, fmax = 20.0 |

#### How the three methods differ

| | τ_k spacing | Weight computation | Speed | Accuracy guarantee |
|---|---|---|---|---|
| `nnls` | Log-spaced in [1/(2π·fmax), 1/(2π·fmin)] | NNLS minimize Σ\|Q_realized(f_i) - Q_target(f_i)\| | ~20k sweeps | Bounded by `max_fit_error` |
| `withers-table` | Withers eq. (15) with tabulated τ_min(γ), τ_max(γ) | Table 1 (Q>200): w_k = W_k/Q. Table 2 (15≤Q≤200): w_k = a_k/Q² + b_k/Q | Instant | Validated in paper for tabulated γ |
| `withers-method` | Withers eq. (15) with γ-dependent τ_min, τ_max (interpolated/extrapolated for any γ, N) | NNLS with Withers objective + non-negative constraint | ~20k sweeps | Bounded by `max_fit_error` |

**Key distinction:** `withers-table` uses the *precomputed result* from the paper. `withers-method` uses the *algorithm* from the paper to compute new results at runtime. `nnls` uses standard log-spaced relaxation times (not γ-dependent Withers spacing).

#### `withers-method` τ_k computation for arbitrary N

For N ≠ 8, Withers eq. (15) generalizes naturally:
```
τ_k = exp( ln(τ_min) + (2k-1)/(2N) · ln(τ_max/τ_min) )    k = 1..N
```

τ_min(γ) and τ_max(γ) are interpolated from the Withers table for tabulated γ values (0.0–0.9) and extrapolated (clamped) beyond. This gives the same logarithmic spacing philosophy but with fewer or more mechanisms.

### Interaction Rules

| Condition | Behavior |
|---|---|
| `attenuation = 'constant-Q'` | `gamma` and `f_transition` ignored (internally gamma=0, f_transition=fref) |
| `attenuation = 'frequency-Q'` and `gamma = 0` | Valid. Constant Q via frequency-Q machinery (Withers τ spacing, not log-spaced) |
| `nblocks = 1` and `Qs0` is array(2) | Error: per-block Q requires nblocks = 2 |
| `nblocks = 2` and `Qs0` is scalar | Both blocks get same Q (uniform model, still 2-block geometry) |
| `weight_method = 'withers-table'` and `n_mechanisms ≠ 8` | Error |
| `weight_method = 'withers-table'` and `attenuation = 'constant-Q'` | Error (use `nnls` or `fixed-q50`) |
| `weight_method = 'withers-table'` and `gamma > 0.9` | Error: outside tabulated range |
| `weight_method = 'withers-method'` and `attenuation = 'constant-Q'` | Valid. Uses Withers τ spacing with gamma=0 internally |
| `weight_method = 'fixed-q50'` and `n_mechanisms ∉ {4, 8}` | Error |
| `weight_method = 'fixed-q50'` and `attenuation = 'frequency-Q'` | Error |
| `weight_method = 'fixed-q50'` and `fmin ≠ 0.05` or `fmax ≠ 20.0` | Error |
| `Qs0 < 15` or `Qp0 < 15` | Error: below SLS validity range |

---

## Implementation Plan

### Phase 1: New Module `src/viscoelastic_model.f90`

Self-contained module, no coupling to legacy code.

```fortran
module viscoelastic_model
  use common, only : wp
  implicit none
  private

  integer, parameter, public :: VE_MAX_MECHANISMS = 8
  integer, parameter, public :: VE_MIN_MECHANISMS = 3

  type, public :: viscoelastic_parameters
     character(len=32) :: attenuation = 'constant-Q'
     real(wp) :: Qs0(2) = -1.0_wp
     real(wp) :: Qp0(2) = -1.0_wp
     real(wp) :: gamma = 0.0_wp
     real(wp) :: fmin = 0.05_wp
     real(wp) :: fmax = 20.0_wp
     real(wp) :: fref = 1.0_wp
     real(wp) :: f_transition = 1.0_wp
     integer  :: n_mechanisms = 8
     character(len=32) :: weight_method = 'nnls'
     integer  :: nnls_samples = 256
     real(wp) :: nnls_tolerance = 1.0e-10_wp
     real(wp) :: max_fit_error = 0.10_wp
     integer  :: nblocks = 1
  end type

  type, public :: viscoelastic_coefficients
     integer :: n
     real(wp) :: tau(VE_MAX_MECHANISMS)
     real(wp) :: w_s(VE_MAX_MECHANISMS)
     real(wp) :: w_p(VE_MAX_MECHANISMS)
  end type

  public :: read_viscoelastic_parameters
  public :: validate_viscoelastic_parameters
  public :: build_viscoelastic_coefficients
  public :: viscoelastic_relaxation_dt_limit
  public :: viscoelastic_max_relative_error
  public :: viscoelastic_realized_q
end module
```

**Key functions:**

- `read_viscoelastic_parameters(infile, params, nblocks, status, message)` — reads `&viscoelastic_list`, stores in `viscoelastic_parameters`
- `validate_viscoelastic_parameters(params, issues)` — preflight: interaction rules, band, Q range
- `build_viscoelastic_coefficients(params, block_id, coeffs)` — dispatches based on `weight_method`:
  - `'nnls'` → `build_nnls_coefficients` (log-spaced τ, NNLS fit)
  - `'withers-table'` → `build_withers_table_coefficients` (table lookup from `withers_tables.f90`)
  - `'withers-method'` → `build_withers_method_coefficients` (Withers τ spacing, NNLS fit)
  - `'fixed-q50'` → `build_fixed_q50_coefficients` (legacy pre-computed)
- `viscoelastic_relaxation_dt_limit(coeffs)` — `min(2·τ_k)` stability limit
- `viscoelastic_max_relative_error(target_Q, coeffs, fmin, fmax)` — max |ΔQ/Q| over band
- `viscoelastic_realized_q(coeffs, f)` — Q(f) from given coefficients (diagnostic)

#### Internal coefficient builders

**`build_nnls_coefficients(N, Q_s, Q_p, gamma, f_transition, fmin, fmax, nnls_samples, nnls_tolerance, coeffs)`**
1. Compute τ_k: log-spaced from `1/(2π·fmax)` to `1/(2π·fmin)`, N mechanisms
2. Build target Q(f) curve at `nnls_samples` log-spaced frequencies over [fmin, fmax]:
   - `constant-Q`: Q_target(f) = Q₀
   - `frequency-Q`: Q_target(f) = Q₀ for f < f_transition, Q₀·(f/f_transition)^γ above
3. NNLS solve for non-negative weights w_k that minimize Σ|Q_realized(f_i) - Q_target(f_i)|
4. Separate fits for S and P weights

**`build_withers_table_coefficients(Q_s, Q_p, gamma, f_transition, coeffs)`**
1. `call get_relaxation_times(gamma, tau)` → 8 τ_k from Withers eq. (15) with tabulated τ_min/τ_max
2. Scale: `tau = tau / f_transition`
3. `call get_withers_weights(gamma, Q_s, w_s)` → table lookup:
   - Q > 200: `w_k = W_k(γ) / Q` (Table 1)
   - 15 ≤ Q ≤ 200: `w_k = a_k(γ)/Q² + b_k(γ)/Q` (Table 2)
4. Same for Q_p → w_p

**`build_withers_method_coefficients(N, Q_s, Q_p, gamma, f_transition, fmin, fmax, nnls_samples, nnls_tolerance, coeffs)`**
1. Compute γ-dependent τ_min, τ_max:
   - For γ ∈ [0.0, 0.9]: interpolate from Withers TAU_MIN/TAU_MAX tables
   - For γ > 0.9: extrapolate (clamped to γ=0.9 values)
2. Compute τ_k via generalized Withers eq. (15): `τ_k = exp(ln(τ_min) + (2k-1)/(2N) · ln(τ_max/τ_min))`
3. Scale: `tau = tau / f_transition`
4. NNLS solve for weights (same solver as `nnls`, but with Withers-spaced τ_k instead of log-spaced)

**`build_fixed_q50_coefficients(N, Q_s, Q_p, coeffs)`**
1. Pre-stored τ_k and W_k for N=4 or N=8, optimized for Q=50 over [0.05, 20] Hz
2. Scale: `w_s = W_k / Q_s`, `w_p = W_k / Q_p`

### Phase 2: Unified Memory Variables in `datatypes.f90`

Replace 6 separate sets of `eta4X/Deta4X` (X = Q, Q8, cQ, Qf, Qf8, constQ) with **one**:

```fortran
! In block_material:
logical :: viscoelastic = .false.
integer :: n_mechanism_ve = 0
real(wp), allocatable :: tau_ve(:)               ! (n_mech)
real(wp), allocatable :: weight_s_ve(:)          ! (n_mech)
real(wp), allocatable :: weight_p_ve(:)          ! (n_mech)
real(wp), allocatable :: Qs_inv_ve(:,:,:)        ! (nx,ny,nz) — 1/Q field
real(wp), allocatable :: Qp_inv_ve(:,:,:)        ! (nx,ny,nz)
real(wp), allocatable :: eta4_ve(:,:,:,:)        ! (nx,ny,nz,n_mech) — σ_xx memory
real(wp), allocatable :: eta5_ve(:,:,:,:)        ! σ_yy
real(wp), allocatable :: eta6_ve(:,:,:,:)        ! σ_zz
real(wp), allocatable :: eta7_ve(:,:,:,:)        ! σ_xy
real(wp), allocatable :: eta8_ve(:,:,:,:)        ! σ_xz
real(wp), allocatable :: eta9_ve(:,:,:,:)        ! σ_yz
real(wp), allocatable :: Deta4_ve(:,:,:,:)       ! rates (same shapes)
real(wp), allocatable :: Deta5_ve(:,:,:,:)
real(wp), allocatable :: Deta6_ve(:,:,:,:)
real(wp), allocatable :: Deta7_ve(:,:,:,:)
real(wp), allocatable :: Deta8_ve(:,:,:,:)
real(wp), allocatable :: Deta9_ve(:,:,:,:)
```

**Memory per grid point:** `12 × N × 8 bytes` (eta + Deta, double precision).
- N=3: 288 bytes/point
- N=4: 384 bytes/point
- N=8: 768 bytes/point

### Phase 3: Single RHS Kernel

One `apply_viscoelastic_point` replaces 6+ current routines:

```fortran
subroutine apply_viscoelastic_point(M, x, y, z, Ux, Uy, Uz, DU)
  type(block_material), intent(inout) :: M
  integer, intent(in) :: x, y, z
  real(wp), intent(in) :: Ux(:), Uy(:), Uz(:)
  real(wp), intent(inout) :: DU(:)
  integer :: k, n
  real(wp) :: tr, mu2, sm, pm, bulk

  n = M%n_mechanism_ve

  DU(4) = DU(4) - sum(M%eta4_ve(x,y,z,1:n))
  DU(5) = DU(5) - sum(M%eta5_ve(x,y,z,1:n))
  DU(6) = DU(6) - sum(M%eta6_ve(x,y,z,1:n))
  DU(7) = DU(7) - sum(M%eta7_ve(x,y,z,1:n))
  DU(8) = DU(8) - sum(M%eta8_ve(x,y,z,1:n))
  DU(9) = DU(9) - sum(M%eta9_ve(x,y,z,1:n))

  tr = Ux(1) + Uy(2) + Uz(3)
  mu2 = 2.0_wp * M%M(x,y,z,2)

  do k = 1, n
     sm = M%weight_s_ve(k) * M%Qs_inv_ve(x,y,z)
     pm = M%weight_p_ve(k) * M%Qp_inv_ve(x,y,z)
     bulk = (M%M(x,y,z,1) + mu2) * pm - mu2 * sm

     M%Deta4_ve(x,y,z,k) = M%Deta4_ve(x,y,z,k) + &
          (mu2*sm*Ux(1) + bulk*tr - M%eta4_ve(x,y,z,k)) / M%tau_ve(k)
     M%Deta5_ve(x,y,z,k) = M%Deta5_ve(x,y,z,k) + &
          (mu2*sm*Uy(2) + bulk*tr - M%eta5_ve(x,y,z,k)) / M%tau_ve(k)
     M%Deta6_ve(x,y,z,k) = M%Deta6_ve(x,y,z,k) + &
          (mu2*sm*Uz(3) + bulk*tr - M%eta6_ve(x,y,z,k)) / M%tau_ve(k)
     M%Deta7_ve(x,y,z,k) = M%Deta7_ve(x,y,z,k) + &
          (M%M(x,y,z,2)*sm*(Uy(1)+Ux(2)) - M%eta7_ve(x,y,z,k)) / M%tau_ve(k)
     M%Deta8_ve(x,y,z,k) = M%Deta8_ve(x,y,z,k) + &
          (M%M(x,y,z,2)*sm*(Uz(1)+Ux(3)) - M%eta8_ve(x,y,z,k)) / M%tau_ve(k)
     M%Deta9_ve(x,y,z,k) = M%Deta9_ve(x,y,z,k) + &
          (M%M(x,y,z,2)*sm*(Uz(2)+Uy(3)) - M%eta9_ve(x,y,z,k)) / M%tau_ve(k)
  end do
end subroutine
```

Plus `apply_viscoelastic_point_pml` variant (same structure, PML damping added to strains).

### Phase 4: Time Integration in `fields.f90`

Replace 6 if-blocks with one:

```fortran
! In zero_rates:
if (F%M%viscoelastic) then
   F%M%Deta4_ve = A * F%M%Deta4_ve
   F%M%Deta5_ve = A * F%M%Deta5_ve
   F%M%Deta6_ve = A * F%M%Deta6_ve
   F%M%Deta7_ve = A * F%M%Deta7_ve
   F%M%Deta8_ve = A * F%M%Deta8_ve
   F%M%Deta9_ve = A * F%M%Deta9_ve
end if

! In update_fields:
if (F%M%viscoelastic) then
   F%M%eta4_ve = F%M%eta4_ve + dt * F%M%Deta4_ve
   F%M%eta5_ve = F%M%eta5_ve + dt * F%M%Deta5_ve
   F%M%eta6_ve = F%M%eta6_ve + dt * F%M%Deta6_ve
   F%M%eta7_ve = F%M%eta7_ve + dt * F%M%Deta7_ve
   F%M%eta8_ve = F%M%eta8_ve + dt * F%M%Deta8_ve
   F%M%eta9_ve = F%M%eta9_ve + dt * F%M%Deta9_ve
end if
```

### Phase 5: Integration Points

| File | Change |
|---|---|
| `src/viscoelastic_model.f90` | **New.** Parameter type, reader, validator, coefficient builders |
| `src/datatypes.f90` | Add `viscoelastic` flag + unified `_ve` arrays. Keep legacy arrays for backward compat initially |
| `src/simulation_config.f90` | Add `viscoelastic_parameters` to `simulation_config_t` |
| `src/input_preflight.f90` | Add `read_viscoelastic_parameters` call for `response = 'viscoelastic'`. Validation in `validate_viscoelastic_parameters` |
| `src/block.f90` | Add `init_viscoelastic_properties` call in `init_block`. Allocates `_ve` arrays, builds coefficients, sets `M%viscoelastic = .true.` |
| `src/domain.f90` | Add `'viscoelastic'` to response dispatch: dt limit, init, finiteness check |
| `src/RHS_Interior.f90` | Add `if (M%viscoelastic) call apply_viscoelastic_point_dispatch(...)` at all integration paths. One dispatch handles interior and PML |
| `src/fields.f90` | Add `viscoelastic` block in `zero_rates` and `update_fields_interior` |
| `src/material.f90` | Add `init_viscoelastic_properties(M, G, coeffs)` — allocates arrays, fills tau/weight/Q_inv fields |

### Phase 6: Preflight Diagnostics

| Code | Condition |
|---|---|
| `CFG-VE-001` | Cannot parse `&viscoelastic_list` |
| `CFG-VE-002` | Invalid `attenuation` value |
| `CFG-VE-003` | `n_mechanisms` out of range (3–8) |
| `CFG-VE-004` | `weight_method` / `attenuation` / `n_mechanisms` incompatible combination |
| `CFG-VE-005` | Q₀ < 15 (below SLS validity) |
| `CFG-VE-006` | Band invalid (fmin ≤ 0, fmax ≤ fmin, fref outside band) |
| `CFG-VE-007` | Per-block Q with nblocks = 1 |
| `CFG-VE-008` | NNLS fit error exceeds `max_fit_error` |
| `CFG-VE-009` | `withers-table` with γ > 0.9 (outside tabulated range) |

Preflight prints Q(f) diagnostic:
```
viscoelastic: n_mechanisms = 8, weight_method = nnls
  block 1: Qs0 = 50.0, Qp0 = 100.0
  Q_S fit: max relative error = 1.2% over [0.05, 20.0] Hz
  Q_P fit: max relative error = 0.8% over [0.05, 20.0] Hz
  relaxation dt limit = 0.00234 s
```

---

## Typical Input Configurations

### Constant Q, simple
```fortran
&problem_list response = 'viscoelastic' /
&viscoelastic_list
 attenuation = 'constant-Q',
 Qs0 = 50d0, Qp0 = 100d0,
 fref = 1d0,
 n_mechanisms = 8,
 weight_method = 'nnls'
/
```

### Frequency-dependent Q — exact Withers paper replication
```fortran
&problem_list response = 'viscoelastic' /
&viscoelastic_list
 attenuation = 'frequency-Q',
 Qs0 = 50d0, Qp0 = 100d0,
 gamma = 0.6d0,
 f_transition = 1d0, fref = 1d0,
 n_mechanisms = 8,
 weight_method = 'withers-table'
/
```

### Frequency-dependent Q — Withers methodology, 5 mechanisms
```fortran
&problem_list response = 'viscoelastic' /
&viscoelastic_list
 attenuation = 'frequency-Q',
 Qs0 = 50d0, Qp0 = 100d0,
 gamma = 0.6d0,
 f_transition = 1d0, fref = 1d0,
 n_mechanisms = 5,
 weight_method = 'withers-method'
/
```

### Layered model (per-block Q)
```fortran
&problem_list response = 'viscoelastic', nblocks = 2 /
&viscoelastic_list
 attenuation = 'constant-Q',
 Qs0 = 20d0, 210d0,
 Qp0 = 40d0, 420d0,
 fref = 1d0,
 n_mechanisms = 6,
 weight_method = 'nnls'
/
```

### Low-cost 3-mechanism (memory-constrained runs)
```fortran
&problem_list response = 'viscoelastic' /
&viscoelastic_list
 attenuation = 'constant-Q',
 Qs0 = 100d0, Qp0 = 200d0,
 fmin = 0.1d0, fmax = 5d0, fref = 1d0,
 n_mechanisms = 3,
 weight_method = 'nnls'
/
```

---

## Backward Compatibility

- All existing `response` values (`anelastic-Q4`, `anelastic-Q8`, `anelastic-fQ8`, etc.) continue to work unchanged
- Legacy code paths remain functional — no deletion in this phase
- `viscoelastic` is a new, parallel code path validated against existing responses
- Migration: once validated, legacy responses emit deprecation warnings pointing to `viscoelastic`

## Verification Plan

1. **Coefficient identity**: for each legacy response, construct equivalent `viscoelastic` parameters and verify `build_viscoelastic_coefficients` produces identical τ_k and w_k
2. **Seismogram identity**: run each Withers benchmark case with both legacy and `viscoelastic` response, compare seismograms (should be bit-identical)
3. **Variable mechanisms**: run constant-Q with N = 3, 4, 5, 6, 7, 8 and verify Q(f) fit accuracy monotonically improves with N
4. **Per-block Q**: run layered benchmark with `viscoelastic` and compare to `anelastic-cQ8-b2`
5. **Memory**: verify N=3 allocates 18 3D fields vs N=8 allocates 48
6. **Preflight**: test all error codes with invalid inputs
7. **withers-method vs withers-table**: for N=8 and tabulated γ, verify `withers-method` produces coefficients matching `withers-table` within NNLS tolerance

## Files Created/Modified (Summary)

| Action | File | Lines (est.) |
|---|---|---|
| **New** | `src/viscoelastic_model.f90` | ~500 |
| Modify | `src/datatypes.f90` | +20 |
| Modify | `src/simulation_config.f90` | +5 |
| Modify | `src/input_preflight.f90` | +30 |
| Modify | `src/block.f90` | +15 |
| Modify | `src/material.f90` | +60 |
| Modify | `src/domain.f90` | +20 |
| Modify | `src/RHS_Interior.f90` | +80 |
| Modify | `src/fields.f90` | +15 |

---

## Future: Coarse-Graining (Not In Scope)

Withers (2015) §3.2 describes a memory-saving coarse-grained layout where each grid node stores only **one** mechanism instead of all N. In a period-two (2×2×2) layout, 8 mechanisms cycle across 8 nodes. Memory drops from 6N to 6 fields per block — an 8× reduction for N=8.

This is orthogonal to the `weight_method` / `attenuation` / `n_mechanisms` design and can be added later as an optional parameter:

```fortran
coarse_grain = 0    ! default: all N mechanisms at every grid point (collocated)
coarse_grain = 2    ! period-two layout, N=8 only, one mechanism per node
```

Implementation requires:
- Per-node mechanism assignment map based on grid index parity
- Modified RHS kernel: each point updates only its assigned mechanism
- Weight scaling: w_k × N to compensate for reduced update frequency
- Constraint: `n_mechanisms = 8` and grid dimensions divisible by 2

This is a pure performance optimization with no change to the physics or user-facing Q behavior. It can be validated by comparing coarse-grained seismograms against the collocated baseline.
