# Withers et al. (2015) f-k Benchmark Implementation

## Task 1: RHS Kernel Verification (COMPLETE)

Verified `src/RHS_Interior.f90` lines 3460-3576. Collocated path (`coarse_grained_Qf8 == .false.`) correctly updates all 8 mechanisms at every grid node. Coarse-grained path uses period-2 indexing via `1+mod(x-1,2)+2*mod(y-1,2)+4*mod(z-1,2)`. No code changes needed.

Key dispatch at line 3431-3457 routes via `allocated(M%eta4Qf8)`. Stress forcing matches Withers Eq. 17/19: `strength_s*2μ*ε_xx + (strength_p*M_p - strength_s*2μ)*tr` for normals, `strength_s*μ*γ_ij` for shears.

## Task 2: f-k Reference Generator (COMPLETE)

**File:** `python/fk_reference.py`

### Coefficient Fitting
- `get_relaxation_times`: Withers Eq. 15 log-spacing, scaled by f_transition
- `fit_conventional_strengths`: NNLS coordinate-descent matching `anelastic_fq8_model.f90`
- `build_fq8_coefficients`: (tau, strength_s, strength_p) for conventional-nnls, coarse_grain=0

### Complex Modulus & Dispersion
- `complex_modulus_ratio`: M(f)/M_u = 1 - Σ λ_k/(1+iωτ_k) [Withers Eq. 6]
- `modulus_scale`: Unrelaxed→physical correction so c(fref) = c_ref
- `realized_q`: Q(f) = Re(M)/Im(M)
- `complex_velocity`, `phase_velocity`: Dispersive velocities with modulus-scale normalization

### Propagation
- `plane_wave_seismogram`: 1D IFFT[S(ω)·exp(-ik(ω)·x)] with complex wavenumber
- `farfield_halfspace_seismograms`: Point source with 1/R spreading, free-surface factor 2, frequency-dependent complex velocities

### Kristekova (2009) Misfit Metrics
- `envelope_phase_misfit(ref, test, dt)`: Returns (EM%, PM%) — envelope and phase misfits via analytic signal / Hilbert transform
- `compare_seismograms`: File-based wrapper for EM/PM comparison

### 4 Benchmark Configurations
- `make_elastic_halfspace`: Q=∞, rho=2700, vs=3464, vp=6000 (Fig. 3)
- `make_constant_q_halfspace`: gamma=0, Qs=Qp=50 (Fig. 4)
- `make_powerlaw_halfspace`: gamma=0.6, Qs=Qp=50 (Fig. 5)
- `make_layered_model`: Q=20/210, gamma=0.6 (Table 3 / Fig. 6)

### CLI
```
python fk_reference.py --selftest          # 15 verification checks
python fk_reference.py --case elastic      # single benchmark
python fk_reference.py --all               # all 4 benchmarks
python fk_reference.py --diagnostics       # Q(f) and dispersion curves
python fk_reference.py --compare REF TEST  # EM/PM comparison
```

### Selftest Checks (15 total)
1. Constant-Q (gamma=0, Q=50): P/S symmetry, sum(strength)<1, Q accuracy <0.2%, modulus scale
2. Power-law (gamma=0.6, Q=50): P/S symmetry, Q accuracy <10%
3. Elastic plane-wave arrival time (within 10ms of R/c_ref)
4. Q=50 amplitude reduction vs elastic
5. Low-Q (Q=20, gamma=0.6): Q accuracy <10%, sum(strength)<1
6. High-Q (Q=210, gamma=0.6): Q accuracy <10%
7. Kristekova EM/PM consistency: identical→0%, shifted→nonzero PM, 80% scale→~20% EM

## Task 3: Automated CTest Tests (COMPLETE)

### Fortran Unit Test
**File:** `tests/fq8_withers_benchmark_test.f90`

Tests coefficient fitting for all 4 benchmark cases:
- Elastic (Q=1e10): strengths ≈ 0, modulus scale = 1.0
- Constant-Q (Q=50, gamma=0): Q accuracy <0.2%, non-negative strengths, P/S symmetry
- Power-law (Q=50, gamma=0.6): Q accuracy <10%, non-negative, sum<1
- Layered Q=20 (gamma=0.6): Q accuracy <10%, non-negative, sum<1
- Layered Q=210 (gamma=0.6): Q accuracy <10%, non-negative

**CTest name:** `fq8_withers_benchmark_unit`

### Python Selftest
**CTest name:** `fk_reference_selftest` (requires Python3)

### MPI Benchmark Tests
Using `run_q8_dynamic_regression.cmake` (decomposition-consistency: 1-rank vs 2-rank):
- `fq8_benchmark_elastic` — `inputfile/test_fq8_benchmark_elastic.in`
- `fq8_benchmark_constant_q` — `inputfile/test_fq8_benchmark_constant_q.in`
- `fq8_benchmark_powerlaw` — `inputfile/test_fq8_benchmark_powerlaw.in`

### Comparison Infrastructure
**File:** `cmake/run_fq8_benchmark.cmake` — template for full solver-vs-reference EM/PM comparison (runs solver, generates reference, compares seismograms)

### Acceptance Criteria (from instructions)
- Half-space cases: EM and PM should be a few percent
- Layered case: up to ~8-10% EM / ~2% PM for later surface-wave arrivals

## Files Created/Modified

| File | Action |
|------|--------|
| `python/fk_reference.py` | Created — f-k reference generator |
| `tests/fq8_withers_benchmark_test.f90` | Created — Fortran coefficient unit test |
| `cmake/run_fq8_benchmark.cmake` | Created — benchmark comparison cmake script |
| `inputfile/test_fq8_benchmark_elastic.in` | Created — elastic benchmark input |
| `inputfile/test_fq8_benchmark_constant_q.in` | Created — constant-Q benchmark input |
| `inputfile/test_fq8_benchmark_powerlaw.in` | Created — power-law benchmark input |
| `src/CMakeLists.txt` | Modified — added 6 new CTest entries |

## Task 4: Non-goals (honored)
- Did not touch `anelastic-Q8`, `anelastic-cQ8-b2`, or `anelastic-cQ`
- Did not attempt to fix `anelastic-Qf` (N=4)
