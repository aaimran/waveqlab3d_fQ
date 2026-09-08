# Instructions: Match WaveQLab3D-fQ to the Withers et al. (2015) f-k Benchmark

## Context

We need the `anelastic-fQ8` response (full-mechanism, collocated layout,
`coefficient_method='conventional-nnls'`, `coarse_grain=0`) to correctly
reproduce the Withers, Olsen & Day (2015, BSSA) f-k benchmark. A prior audit
of `src/` (see summary below) found the coefficient-fitting side is already
implemented correctly. The main open item is verifying the RHS time-stepping
kernel actually applies the collocated (full 8-mechanism-per-node) update
consistently, and adding a dedicated validation test against a freshly
computed f-k reference (not digitized paper figures).

## Background: audit summary (as of this session)

**Correct, already implemented:**
- `src/withers_tables.f90` — verbatim Table 1 (high-Q, `Q>200`) and Table 2
  (low-Q coefficients `a_k`, `b_k`, `15<Q<200`) from the paper, plus
  `get_relaxation_times` (Eq. 15 log-spacing). These are valid **only** for
  the coarse-grained, one-mechanism-per-node, period-2 layout (`w_k = N*λ_k`).
- `src/anelastic_fq8_model.f90` — implements two coefficient methods:
  - `'withers-2015'`: uses the raw published table via `get_withers_weights`.
    Only legal with `coarse_grain=2` (the code already rejects
    `coarse_grain=0` + `withers-2015` as invalid — see
    `read_fq8_parameters` / `build_fq8_coefficients`).
  - `'conventional-nnls'`: fits λ_k directly against the **exact** relation
    `Im(M) - Re(M)/Q = 0` (not the low-loss approximation) via nonnegative
    coordinate-descent least squares (`fit_conventional_strengths`), summed
    over all 8 mechanisms with no coarse-cell scaling. This is the correct
    method for `coarse_grain=0` (full mechanism / collocated) and is
    mathematically equivalent to the iterative exact fit the original paper
    used to build its own f-k reference.
- `src/material.f90` → `init_anelastic_Qf8_properties` — branches correctly
  on `M%coarse_grained_Qf8` (`= parameters%coarse_grain == 2`). For the
  `.false.` (collocated) case, the unrelaxed-modulus correction loop sums
  all 8 mechanisms per node with `denom = (ω_ref τ_k)² + 1`, matching a
  true full-mechanism collocated scheme (no harmonic/coarse-cell averaging).
  Startup diagnostics print `"layout = full 8 mechanisms at every grid
  point"` and `"modulus normalization = collocated full-mechanism
  reference"` for this path, and `fq8_max_relative_error` checks realized-Q
  error over `[0.1, 10] x f_transition`.

**Known-bad / do not use for benchmarking:**
- `src/material.f90` → `init_anelastic_QF_properties` (the older, separate
  `anelastic-Qf`, N=4 response) calls
  `withers_tables.f90:get_withers_weights_Qf`, which just **truncates the
  8-mechanism table to its first 4 entries**
  (`weights_Qf(1:4) = weights_full(1:4)`). This is not a valid N=4 fit to
  Q(f) and will not reproduce the benchmark at any Q. Treat `anelastic-Qf`
  as deprecated for this work; use `anelastic-fQ8` instead.

## Tasks

### 1. Verify the RHS kernel applies the collocated update correctly

Find where `RHS_Interior.f90` (and any dispatch in `time_step.f90` /
`CouplingForcing.f90`) consumes the `Qf8` state
(`M%tau_Qf8`, `M%strength_s_Qf8`, `M%strength_p_Qf8`, `M%eta{4..9}Qf8`,
`M%Deta{4..9}Qf8`, `M%coarse_grained_Qf8`).

Confirm:
- When `M%coarse_grained_Qf8 == .false.`, the memory-variable update loop
  advances **all 8 mechanisms at every grid point** (no `MODULO(i+j+k,2)`-
  style coarse/period-2 indexing inherited from the `anelastic_Q8`
  coarse-grained kernel).
- When `M%coarse_grained_Qf8 == .true.`, it correctly uses the period-2,
  one-mechanism-per-node coarse layout instead.
- The stress-update equation matches Eq. (17)/(19) of Withers et al.
  (2015): per-mechanism exponential-integrator decay
  `exp(-Δt/τ_k)` and forcing `strength_k * (1 - exp(-Δt/τ_k)) * ε(t)`,
  with independent S and P strength sets.

If the collocated case is silently reusing the coarse-grained update path,
fix it to iterate all 8 mechanisms per node, gated on
`M%coarse_grained_Qf8`. Do not modify the coarse-grained (`anelastic-Q8`,
`anelastic-cQ8-b2`) code paths — this response must remain fully additive
(see `roust_layered_Q.md`'s non-interference contract, which the `fQ8`
work should also honor).

### 2. Add a from-scratch f-k reference generator (do not use digitized figures)

Create a small standalone tool (Python or Fortran, under `tools/` or
`python/`) that:
- Takes the same `(τ_k, strength_s, strength_p)` produced by
  `build_fq8_coefficients` (conventional-nnls, coarse_grain=0) for a given
  `(Qs0, Qp0, gamma, f_transition, fref)`.
- Builds the complex modulus `M(ω) = M_u * (1 - Σ λ_k/(1+iωτ_k))` (Eq. 6),
  and uses it to compute a complex wavenumber / dispersive wavespeed for an
  f-k (frequency-wavenumber) solution — following the iterative approach
  described in Withers et al. (2015) (Numerical Tests section): initialize
  λ_k = 0 in the exact Eq. (7) denominator, solve the resulting linear
  system by least squares, iterate to convergence (~2-3% tolerance).
- Emits synthetic seismograms (or at minimum the Q(f) and dispersion
  curves) that WaveQLab3D-fQ output can be diffed against directly, so the
  benchmark check is self-consistent rather than relying on pixel-extracted
  paper figures.

### 3. Test cases to reproduce (Withers et al. 2015)

Using `anelastic-fQ8`, `coefficient_method='conventional-nnls'`,
`coarse_grain=0`:

1. Elastic half-space (Q=∞) — sanity check zero-attenuation path against
   the f-k reference (Fig. 3 config).
2. Constant-Q half-space: `gamma=0.0`, `Qs0=Qp0=50` (Fig. 4).
3. Power-law half-space: `gamma=0.6`, `Qs0=Qp0=50` (Fig. 5).
4. Layered model, `gamma=0.6`, `Qs0=Qp0=20` (shallow layer) /
   `Qs0=Qp0=210` (half-space) — Table 3 / Fig. 6. Tests low-Q accuracy and
   a sharp Q discontinuity.

Acceptance criteria (per the paper): envelope misfit (EM) and phase misfit
(PM), computed the way Kristekova et al. (2009) define them, should be a
few percent for the half-space cases and up to ~8-10% EM / ~2% PM for the
later surface-wave arrivals in the layered case. Add these as automated
tests (parallel to the existing `fq8_effective_response_test` /
`cq_coefficients_test` pattern in `CMakeLists.txt`), not just one-off runs.

### 4. Non-goals

- Do not touch `anelastic-Q8`, `anelastic-cQ8-b2`, or `anelastic-cQ`
  numerical paths.
- Do not attempt to "fix" `anelastic-Qf` (N=4) — it's superseded by
  `anelastic-fQ8` with `n_mechanisms` effectively fixed at 8; if a smaller
  mechanism count is wanted later, that should be a new, explicitly-fitted
  variant of `fQ8`, not a patch to the truncation bug.
