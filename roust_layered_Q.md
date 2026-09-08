# Robust Layered-Q Implementation Plan

## Implementation status

The additive `anelastic-cQ` response described here is implemented. It has a
separate namelist, configuration object, fitted coefficients, material state,
memory arrays, modulus correction, constitutive routines, diagnostics, and
cleanup path. Mechanism counts 4 through 8 are supported and covered by
parameterized coefficient and dynamic tests. Traditional, upwind, and
upwind-DRP derivative kernels dispatch the independent cQ update. The
established `anelastic-Q8` and `anelastic-cQ8-b2` numerical routines remain
unchanged.

## Purpose

Add a new, independent response for layered models that fits coefficient sets
to each block's `Qs0` and `Qp0`. The design must not extend, reinterpret, or
modify the numerical paths used by `anelastic-Q8` or `anelastic-cQ8-b2`.

Use a distinct response name throughout this plan:

```fortran
response = 'anelastic-cQ'
```

Here, `cQ` identifies the future robust constant-Q formulation. It remains
different from all existing response names and supports a configurable number
of relaxation mechanisms.

The current implementation remains the baseline:

- Each block has independent `Qs0` and `Qp0`.
- Both blocks share `fref`, `fmin`, `fmax`, and `weight_method`.
- Each block owns its Q8 material and memory-variable arrays.
- One normalized eight-element weight vector is used for both P and S terms in
  a block.
- `weight_method='fixed-q50'` is the only accepted method.

This baseline is appropriate for initial testing. The work below should be
implemented only after those simulations establish the required accuracy and
frequency range.

### Non-interference contract

The future implementation must obey these rules:

- Do not add new methods to `&anelastic_Q8_list` or
  `&anelastic_cQ8_b2_list`.
- Do not change the defaults, validation, coefficient construction, material
  initialization, RHS dispatch, PML kernel, diagnostics, or cleanup behavior
  of either current response.
- Do not reuse current Q8 state arrays for the robust response.
- Do not route the new response through `M%anelastic_Q8`.
- Do not refactor current NNLS or Q8 routines as part of the initial robust
  implementation. Copy required numerical logic into new focused modules,
  with provenance comments and independent tests.
- Existing response tests must pass without updated expected values.
- Deleting the new response's source files and registrations should leave the
  current implementation behaviorally identical.

## Target capability

The new `anelastic-cQ` response may support these coefficient policies:

1. `fixed-q50`: the present shared spectral shape, retained exactly.
2. `nnls-shared`: one fitted normalized shape shared by all blocks and by P/S.
3. `nnls-block`: one fitted shape per block, shared by P and S within that block.
4. `nnls-block-ps`: independent P and S coefficient sets in every block.

`fixed-q50` is valid only when `n_mechanisms=8`, because the frozen table has
eight entries. All NNLS policies must support `n_mechanisms=4,5,6,7,8`.

The most robust final representation for `N = n_mechanisms` is:

```text
block 1: tau_S1, weight_S1, tau_P1, weight_P1
block 2: tau_S2, weight_S2, tau_P2, weight_P2
```

If all materials use the same approximation band, the relaxation times should
normally remain common:

```text
tau_S1 = tau_P1 = tau_S2 = tau_P2
```

Only the weights then need to differ. Separate relaxation times should be
allowed only when blocks explicitly request different frequency bands, because
they introduce multiple stiffness scales and a more restrictive global time
step.

## Input design

Keep the current namelist valid and unchanged. It is shown only as a frozen
compatibility baseline and must not receive new fields:

```fortran
&anelastic_cQ8_b2_list
  Qs0 = 40.0, 100.0
  Qp0 = 80.0, 180.0
  fref = 1.0
  fmin = 0.05
  fmax = 20.0
  weight_method = 'fixed-q50'
/
```

Create a separate namelist for the robust response:

```fortran
&anelastic_cQ_list
  Qs0 = 40.0, 100.0
  Qp0 = 80.0, 180.0
  fref = 1.0
  fmin = 0.05
  fmax = 20.0
  n_mechanisms = 8
  coefficient_policy = 'nnls-block-ps'
  nnls_samples = 256
  nnls_objective = 'relative-q'
  nnls_tolerance = 1.0e-10
  max_fit_error = 0.05
/
```

Recommended validation rules:

- Require exactly two finite, positive `Qs0` and `Qp0` values.
- Require integer `n_mechanisms` in the inclusive range 4 through 8. Supported
  values are exactly `4`, `5`, `6`, `7`, and `8`.
- Require `0 < fmin < fmax` and `fmin <= fref <= fmax`.
- Require `nnls_samples >= n_mechanisms` and recommend at least 128 samples.
- Reject an unsupported weight method rather than silently falling back.
- Require all fitted weights to be finite and nonnegative.
- Reject a singular or nonconverged fit.
- Reject a modulus correction whose denominator is zero or negative.
- Make `max_fit_error` an explicit acceptance threshold. Do not merely warn
  when the user has requested a fitted method.

If different bands are later needed, introduce explicit two-element arrays
such as `fmin_block(2)` and `fmax_block(2)`. Do not overload scalar fields with
implicit block-dependent behavior.

For the initial additive response, require `nblocks=2` because the namelist
contains two Q pairs. The shorter response name does not imply arbitrary block
count. Generalization to `nblocks > 2` should later replace fixed-size Q arrays
with allocatable configuration data rather than changing the meaning of the
two-block input.

## Configurable mechanism count

Set `n_mechanisms=8` as the new response's default, while accepting every
integer from 4 through 8. The selected value controls:

- The number of logarithmically distributed relaxation times
- The length of every P and S weight vector
- The mechanism dimension of every cQ memory and derivative array
- NNLS matrix column count
- All initialization, RHS, PML, update, diagnostic, and cleanup loops
- The relaxation stability calculation

No robust implementation loop should retain a literal upper bound of 8.
Allocate arrays after preflight validation and iterate to
`M%cQ%n_mechanisms`. Generate relaxation times with the selected `N` in the
same log-spacing formula, so changing `N` does not change the requested
frequency band.

Startup output must report the resolved mechanism count and print exactly `N`
relaxation times and weights for each active spectrum.

## Constitutive representation

The existing Q8 members of `block_material`, including the following, are
reserved for the current responses and must not be read or written by the new
response:

```fortran
weight_Q8(8)
```

Add a separate nested state type, preferably in a new
`anelastic_cq_types.f90` module. Coefficient and memory arrays must be
allocatable because their last dimension is selected at runtime:

```fortran
type :: anelastic_cq_state
  logical :: active = .false.
  integer :: n_mechanisms = 0
  real(wp), allocatable :: weight_s(:), weight_p(:)
  real(wp), allocatable :: tau_s(:), tau_p(:)
  real(wp), allocatable :: qs_inv(:,:,:), qp_inv(:,:,:)
  ! Independent cQ memory variables. The mechanism dimension has extent
  ! n_mechanisms; exact leading dimensions follow the P/S derivation.
  real(wp), allocatable :: eta_s(:,:,:,:,:), deta_s(:,:,:,:,:)
  real(wp), allocatable :: eta_p(:,:,:,:,:), deta_p(:,:,:,:,:)
end type anelastic_cq_state
```

Attach this as a new member of `block_material`, for example `M%cQ`, without
changing existing Q8 members or their meaning.

The block itself provides layer separation. `D%B(i)%M%cQ` should own all robust
coefficients and memory variables for block `i`. The new response must have its
own allocation, initialization, finite-state checks, and destruction routines.

In this document, “number of memory variables” means the number of relaxation
mechanisms per stress-memory family. For three-dimensional stress there are
multiple component families, each allocated with a final extent of
`n_mechanisms`. The input should therefore use the physically explicit name
`n_mechanisms` while documentation can describe it as the requested memory-
variable count.

## NNLS formulation

Create a new `anelastic_cq_nnls.f90` coefficient module. Initially copy
and adapt the proven NNLS algorithm from `material.f90`; do not move or alter
the existing routine because that could change established responses. Record
the source routine and copy date in comments so later fixes can be audited.
The new module should expose a side-effect-free interface
similar to:

```fortran
subroutine fit_constant_q_weights(target_q, tau, fmin, fmax, nsamples, &
                                  objective, weights, status, message)
```

The fit should:

1. Sample frequency logarithmically over `[fmin, fmax]`.
2. Construct the linearized generalized-standard-linear-solid response matrix.
3. Solve the nonnegative least-squares problem.
4. Evaluate the realized Q using the full constitutive response, not only the
   linearized fitting residual.
5. Report maximum and RMS relative-Q errors.
6. Evaluate phase-velocity/modulus error at `fref` after correction.

For `nnls-block-ps`, run four independent fits:

```text
fit(Qs0(1)) -> block 1 S weights
fit(Qp0(1)) -> block 1 P weights
fit(Qs0(2)) -> block 2 S weights
fit(Qp0(2)) -> block 2 P weights
```

For `nnls-block`, define the compromise objective explicitly. A recommended
choice is a stacked relative-error system containing both the block's `Qs0`
and `Qp0` targets with equal weighting. Do not quietly fit only `Qs0` and reuse
that result for P.

## RHS changes

Do not edit the current `apply_anelastic_Q8_*` routines. Add independent robust
kernels:

```text
apply_anelastic_cQ_point
apply_anelastic_cQ_point_pml
apply_anelastic_cQ_point_dispatch
```

The new kernels split the contributions carefully:

- Shear and deviatoric terms use `M%cQ%weight_s` and `M%cQ%qs_inv`.
- The P-wave/volumetric modulus term uses `M%cQ%weight_p` and
  `M%cQ%qp_inv`.
- The bulk contribution must remain the difference between the P and shear
  relaxation effects; it cannot be converted by replacing every existing
  current-Q8 weight reference independently.

Derive the updated normal-stress memory equations on paper and add them to the
code documentation before implementation. A new-response unit test should
verify that equal robust P/S weights reduce algebraically to the established
Q8 equations, but the production implementation must not call the established
Q8 kernel.

Implement both new non-PML and new PML kernels together. No current Q8 kernel
should be modified.

## Reference-modulus correction

Add a new `init_anelastic_cQ_properties` routine rather than modifying
`init_anelastic_Q8_properties`. It corrects the unrelaxed P and S moduli so
input velocities are recovered at `fref`:

- Compute the S correction using `M%cQ%weight_s` and `M%cQ%tau_s`.
- Compute the P correction using `M%cQ%weight_p` and `M%cQ%tau_p`.
- Validate both correction denominators independently.
- Store the requested and corrected moduli for diagnostic testing.

Add a unit test that reconstructs phase velocities at `fref` and verifies they
match the input velocities within a documented tolerance.

## Time-step stability

Compute the relaxation stability limit from every active relaxation spectrum:

```text
minimum over blocks, P spectra, S spectra, and mechanisms
```

If P and S share `tau`, this is identical to the current result. If later
extensions permit different bands, print which block and spectrum limits the
time step.

## Diagnostics

Print a deterministic summary on rank zero:

```text
block 1 S: target Q, max error, RMS error, weights
block 1 P: target Q, max error, RMS error, weights
block 2 S: target Q, max error, RMS error, weights
block 2 P: target Q, max error, RMS error, weights
```

Also print:

- Weight policy and objective
- Frequency band and sample count
- Reference frequency
- Minimum and maximum relaxation times
- Selected relaxation time-step limit
- Whether coefficient sets are shared or independent

The summary should be generated from the resolved configuration before MPI
block ownership divides the model, ensuring all four fits are visible.

## Additive module and integration layout

Create new implementation files:

```text
anelastic_cq_types.f90       configuration and runtime-state types
anelastic_cq_model.f90       namelist validation and fit orchestration
anelastic_cq_nnls.f90        independent NNLS implementation
anelastic_cq_material.f90    allocation, modulus correction, destruction
anelastic_cq_rhs.f90         non-PML and PML constitutive kernels
```

Changes to shared orchestration files must be limited to additive registration:

- Accept the new response name in preflight.
- Parse and broadcast the new configuration object.
- Select the new global relaxation limit.
- Call the new initializer only for `anelastic-cQ`.
- Dispatch the new RHS only when `M%cQ%active` is true.
- Scale/update only the new derivative arrays in field operations.
- Run the new finite-state check and destructor only for the new response.
- Register new source files and tests in CMake.

Do not combine response conditions such as
`anelastic-Q8 .or. anelastic-cQ` around a shared numerical routine. The
selector may be adjacent, but each branch must call a response-specific entry
point. Shared low-level utilities are acceptable only if they are newly added,
pure, independently tested, and cannot change an existing code path.

All new diagnostics should use separate identifiers, for example
`CFG-CQ-*` and `RUN-CQ-*`, so failures cannot be confused with current
Q8 behavior.

## Implementation stages

### Stage 1: coefficient module and unit tests

- Copy/generalize the existing NNLS solver into the new robust module.
- Keep all runtime material and RHS behavior unchanged.
- Test known systems, nonnegative constraints, convergence failure, and
  deterministic results.
- Compare the extracted solver against the existing constant-Q 8M output.

### Stage 2: offline fit diagnostics

- Add `anelastic-cQ`, `&anelastic_cQ_list`, and a separate
  `anelastic_cq_config` object.
- Add `nnls-shared` and `nnls-block` parsing only to the new namelist.
- Compute and report coefficients and realized-Q errors.
- Do not initialize any current Q8 state.
- Compare fixed and fitted spectra for representative Q values.

### Stage 3: block-specific shared P/S weights

- Populate each block's new `M%cQ` state independently.
- Add new response dispatch and robust shared-P/S RHS kernels.
- Validate that different blocks actually receive different coefficients.
- Establish whether this compromise is accurate enough before adding P/S
  arrays.

### Stage 4: independent P/S weights

- Populate the independent P/S fields already reserved in `M%cQ`.
- Update only the robust modulus correction, robust non-PML RHS, and robust PML
  RHS.
- Preserve the robust response's shared policy by assigning identical P/S
  vectors inside `M%cQ`.
- Require algebraic and numerical compatibility tests.

### Stage 5: optional block-specific bands

- Add per-block frequency controls only if simulations require them.
- Extend time-step selection and diagnostics.
- Test interfaces where adjacent blocks have different relaxation spectra.

## Test matrix

### Unit tests

- NNLS returns finite, nonnegative coefficients.
- Repeated fits are deterministic.
- Parameterized tests cover `n_mechanisms=4,5,6,7,8` and verify array extents,
  loop bounds, relaxation-time ordering, and finite fitted results.
- Values below 4, above 8, non-integers at the input boundary, and incompatible
  `fixed-q50`/mechanism-count combinations are rejected collectively.
- Realized Q meets the configured maximum-error threshold.
- Each of the four targets is evaluated independently.
- Modulus correction recovers P- and S-wave velocities at `fref`.
- Invalid bands, sample counts, Q values, and solver failures are rejected.

### Compatibility tests

- Existing `anelastic-Q8` output is bitwise unchanged.
- Existing `anelastic-cQ8-b2` output is bitwise unchanged.
- Existing response source files need no numerical edits and existing expected
  values need no updates.
- With matched parameters, the separate robust response agrees with current
  Q8 within a stated tolerance; it does not achieve this by sharing state or
  calling the current kernel.
- Identical robust P/S coefficient vectors reproduce the established Q8
  equations in an isolated unit test.
- Identical Q pairs in both blocks produce identical coefficients and material
  corrections.

### Dynamic tests

- Run at least one finite, nonzero-memory dynamic case for each supported
  mechanism count `4,5,6,7,8`.
- One-rank shared two-block and two-rank distributed runs agree.
- Asymmetric three-rank decomposition agrees with the supported reference.
- PML and non-PML paths remain finite and evolve nonzero memory variables.
- A source dominated by shear energy responds to the selected `Qs` fit.
- A source dominated by volumetric energy responds to the selected `Qp` fit.
- Swapping block Q pairs swaps the measured block attenuation without changing
  interface stability.

### Physical verification

- Measure amplitude decay and phase shift for monochromatic P and S waves at
  several logarithmically spaced frequencies.
- Recover apparent Q from the simulated decay.
- Compare apparent Q against each block's target and the analytical realized-Q
  curve.
- Include frequencies near `fmin`, `fref`, and `fmax`, where fitting and phase
  errors are most informative.

## Acceptance criteria

The robust implementation is ready when:

- Existing fixed-Q8 and cQ8-b2 results are bitwise unchanged.
- The new response has separate configuration, coefficient, state,
  initialization, RHS/PML, diagnostics, and cleanup paths.
- No current Q8 state flag or memory array is active during an `anelastic-cQ`
  run.
- Every fitted spectrum satisfies its configured maximum relative-Q error.
- All five supported mechanism counts pass unit, non-PML, and representative
  PML coverage; invalid counts fail during preflight.
- P- and S-wave velocities at `fref` match the requested material velocities.
- Serial and MPI decompositions give matching diagnostics and fields.
- PML and non-PML simulations remain finite.
- No unsupported method silently falls back to another coefficient policy.
- Startup output identifies the coefficients and errors for every block and
  wave type.
- Documentation states whether weights and relaxation times are shared or
  independent for every supported policy.

## Recommended decision after current testing

Collect the following from the current fixed-weight simulations before choosing
the next stage:

- The dominant resolved frequency range in each block
- The four realized-Q maximum errors for `Qs0(1)`, `Qp0(1)`, `Qs0(2)`, and
  `Qp0(2)`
- Sensitivity of amplitudes and arrival phases to those errors
- Whether P or S attenuation accuracy dominates the scientific objective
- Runtime and memory constraints for adding separate P/S coefficient arrays

If all four errors are acceptable, retain `fixed-q50`. If errors differ mainly
between blocks, implement `nnls-block`. If P and S errors within a block are
both important and substantially different, proceed directly to
`nnls-block-ps`.
