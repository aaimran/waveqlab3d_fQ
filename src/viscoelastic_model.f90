module viscoelastic_model
  use common, only : wp
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: VE_MAX_MECHANISMS = 8
  integer, parameter, public :: VE_MIN_MECHANISMS = 3
  real(wp), parameter :: pi = 3.141592653589793238462643383279_wp

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
     logical  :: has_ve = .false.
  end type

  type, public :: viscoelastic_coefficients
     integer :: n = 0
     real(wp) :: tau(VE_MAX_MECHANISMS) = 0.0_wp
     real(wp) :: w_s(VE_MAX_MECHANISMS) = 0.0_wp
     real(wp) :: w_p(VE_MAX_MECHANISMS) = 0.0_wp
  end type

  public :: read_viscoelastic_parameters
  public :: build_viscoelastic_coefficients
  public :: ve_relaxation_dt_limit
  public :: ve_params_dt_limit
  public :: ve_realized_q
  public :: ve_max_relative_error

contains

  subroutine read_viscoelastic_parameters(infile, params, nblocks, status, message)
    integer, intent(in) :: infile, nblocks
    type(viscoelastic_parameters), intent(out) :: params
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=32) :: attenuation, weight_method
    real(wp) :: Qs0(2), Qp0(2), gamma, fmin, fmax, fref, f_transition
    real(wp) :: nnls_tolerance, max_fit_error
    integer :: n_mechanisms, nnls_samples, stat

    namelist /viscoelastic_list/ attenuation, Qs0, Qp0, gamma, fmin, fmax, fref, &
         f_transition, n_mechanisms, weight_method, nnls_samples, nnls_tolerance, &
         max_fit_error

    attenuation = params%attenuation
    Qs0 = params%Qs0; Qp0 = params%Qp0
    gamma = params%gamma
    fmin = params%fmin; fmax = params%fmax; fref = params%fref
    f_transition = params%f_transition
    n_mechanisms = params%n_mechanisms
    weight_method = params%weight_method
    nnls_samples = params%nnls_samples
    nnls_tolerance = params%nnls_tolerance
    max_fit_error = params%max_fit_error
    status = 0; message = ''

    rewind(infile)
    read(infile, nml=viscoelastic_list, iostat=stat)
    if (stat /= 0) then
       status = 1
       message = 'viscoelastic requires a valid &viscoelastic_list namelist'
       return
    end if

    attenuation = trim(adjustl(attenuation))
    weight_method = trim(adjustl(weight_method))

    if (attenuation /= 'constant-Q' .and. attenuation /= 'frequency-Q') then
       status = 1; message = 'viscoelastic attenuation must be constant-Q or frequency-Q'; return
    end if

    if (n_mechanisms < VE_MIN_MECHANISMS .or. n_mechanisms > VE_MAX_MECHANISMS) then
       status = 1; message = 'viscoelastic n_mechanisms must be 3-8'; return
    end if

    if (any(.not.ieee_is_finite(Qs0)) .or. any(.not.ieee_is_finite(Qp0))) then
       status = 1; message = 'viscoelastic Qs0/Qp0 must be finite'; return
    end if

    if (nblocks == 1) then
       if (Qs0(1) <= 0.0_wp .or. Qp0(1) <= 0.0_wp) then
          status = 1; message = 'viscoelastic requires positive Qs0 and Qp0'; return
       end if
       if (Qs0(1) < 15.0_wp .or. Qp0(1) < 15.0_wp) then
          status = 1; message = 'viscoelastic Q < 15 is below SLS validity range'; return
       end if
       Qs0(2) = Qs0(1); Qp0(2) = Qp0(1)
    else
       if (Qs0(2) < 0.0_wp) then
          Qs0(2) = Qs0(1); Qp0(2) = Qp0(1)
       end if
       if (any(Qs0(1:2) <= 0.0_wp) .or. any(Qp0(1:2) <= 0.0_wp)) then
          status = 1; message = 'viscoelastic requires positive Qs0 and Qp0 per block'; return
       end if
       if (any(Qs0(1:2) < 15.0_wp) .or. any(Qp0(1:2) < 15.0_wp)) then
          status = 1; message = 'viscoelastic Q < 15 is below SLS validity range'; return
       end if
    end if

    if (.not.ieee_is_finite(fmin) .or. .not.ieee_is_finite(fmax) .or. &
        .not.ieee_is_finite(fref) .or. fmin <= 0.0_wp .or. fmax <= fmin .or. &
        fref < fmin .or. fref > fmax) then
       status = 1; message = 'viscoelastic requires 0 < fmin <= fref <= fmax'; return
    end if

    if (attenuation == 'constant-Q') then
       gamma = 0.0_wp
       f_transition = fref
    else
       if (gamma < 0.0_wp .or. gamma > 1.0_wp) then
          status = 1; message = 'viscoelastic gamma must be in [0, 1]'; return
       end if
       if (.not.ieee_is_finite(f_transition) .or. f_transition <= 0.0_wp) then
          status = 1; message = 'viscoelastic f_transition must be positive'; return
       end if
    end if

    select case (weight_method)
    case ('nnls')
       ! no extra constraints
    case ('withers-table')
       if (n_mechanisms /= 8) then
          status = 1; message = 'withers-table requires n_mechanisms = 8'; return
       end if
       if (attenuation /= 'frequency-Q') then
          status = 1; message = 'withers-table requires attenuation = frequency-Q'; return
       end if
       if (gamma > 0.9_wp) then
          status = 1; message = 'withers-table requires gamma <= 0.9'; return
       end if
    case ('withers-method')
       ! works for any n_mechanisms and attenuation
    case ('fixed-q50')
       if (attenuation /= 'constant-Q') then
          status = 1; message = 'fixed-q50 requires attenuation = constant-Q'; return
       end if
       if (n_mechanisms /= 4 .and. n_mechanisms /= 8) then
          status = 1; message = 'fixed-q50 requires n_mechanisms = 4 or 8'; return
       end if
       if (abs(fmin - 0.05_wp) > 100.0_wp*epsilon(1.0_wp) .or. &
           abs(fmax - 20.0_wp) > 100.0_wp*epsilon(1.0_wp)) then
          status = 1; message = 'fixed-q50 requires fmin=0.05 and fmax=20.0'; return
       end if
    case default
       status = 1
       message = 'viscoelastic weight_method must be nnls, withers-table, withers-method, or fixed-q50'
       return
    end select

    if (nnls_samples < n_mechanisms) then
       status = 1; message = 'viscoelastic nnls_samples must be >= n_mechanisms'; return
    end if

    params%attenuation = attenuation
    params%Qs0 = Qs0; params%Qp0 = Qp0
    params%gamma = gamma; params%fmin = fmin; params%fmax = fmax
    params%fref = fref; params%f_transition = f_transition
    params%n_mechanisms = n_mechanisms; params%weight_method = weight_method
    params%nnls_samples = nnls_samples; params%nnls_tolerance = nnls_tolerance
    params%max_fit_error = max_fit_error; params%nblocks = nblocks
    params%has_ve = .true.
  end subroutine read_viscoelastic_parameters


  subroutine build_viscoelastic_coefficients(params, block_id, coeffs, status, message)
    type(viscoelastic_parameters), intent(in) :: params
    integer, intent(in) :: block_id
    type(viscoelastic_coefficients), intent(out) :: coeffs
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(wp) :: qs, qp

    status = 0; message = ''
    coeffs%n = params%n_mechanisms
    qs = params%Qs0(block_id)
    qp = params%Qp0(block_id)

    select case (trim(params%weight_method))
    case ('nnls')
       call build_nnls(params%n_mechanisms, qs, qp, params%gamma, &
            params%f_transition, params%fmin, params%fmax, &
            params%nnls_samples, params%nnls_tolerance, coeffs)
    case ('withers-table')
       call build_withers_table(qs, qp, params%gamma, params%f_transition, coeffs)
    case ('withers-method')
       call build_withers_method(params%n_mechanisms, qs, qp, params%gamma, &
            params%f_transition, params%fmin, params%fmax, &
            params%nnls_samples, params%nnls_tolerance, coeffs)
    case ('fixed-q50')
       call build_fixed_q50(params%n_mechanisms, qs, qp, params%fmin, params%fmax, coeffs)
    end select

    if (minval(coeffs%tau(1:coeffs%n)) <= 0.0_wp) then
       status = 1; message = 'viscoelastic coefficient build produced non-positive tau'; return
    end if
  end subroutine build_viscoelastic_coefficients


  ! --- τ spacing routines ---

  pure subroutine log_spaced_tau(n, fmin, fmax, tau)
    integer, intent(in) :: n
    real(wp), intent(in) :: fmin, fmax
    real(wp), intent(out) :: tau(VE_MAX_MECHANISMS)
    real(wp) :: taumin, taumax
    integer :: k
    taumin = 1.0_wp / (2.0_wp * pi * fmax)
    taumax = 1.0_wp / (2.0_wp * pi * fmin)
    do k = 1, n
       tau(k) = exp(log(taumin) + (2.0_wp*k - 1.0_wp) / (2.0_wp*n) * log(taumax/taumin))
    end do
  end subroutine log_spaced_tau

  subroutine withers_tau(n, gamma, f_transition, tau)
    use withers_tables, only : get_relaxation_times, N_MECH
    integer, intent(in) :: n
    real(wp), intent(in) :: gamma, f_transition
    real(wp), intent(out) :: tau(VE_MAX_MECHANISMS)
    real(wp) :: tau8(N_MECH), taumin, taumax
    integer :: k

    call get_relaxation_times(gamma, tau8)
    taumin = tau8(1) * exp(-log(tau8(N_MECH)/tau8(1)) / (2.0_wp * N_MECH))
    taumax = tau8(N_MECH) * exp(log(tau8(N_MECH)/tau8(1)) / (2.0_wp * N_MECH))

    if (n == N_MECH) then
       tau(1:n) = tau8(1:n) / f_transition
    else
       do k = 1, n
          tau(k) = exp(log(taumin) + (2.0_wp*k - 1.0_wp) / (2.0_wp*n) * &
               log(taumax/taumin)) / f_transition
       end do
    end if
  end subroutine withers_tau


  ! --- NNLS solver (generalized for any N) ---

  pure subroutine nnls_fit_weights(n, q0, gamma, f_transition, tau, nfreq, tolerance, weights)
    integer, intent(in) :: n, nfreq
    real(wp), intent(in) :: q0, gamma, f_transition, tau(VE_MAX_MECHANISMS), tolerance
    real(wp), intent(out) :: weights(VE_MAX_MECHANISMS)
    integer, parameter :: max_sweeps = 20000
    real(wp) :: a(nfreq, VE_MAX_MECHANISMS), b(nfreq), residual(nfreq)
    real(wp) :: f, target_q, x, denom, old_w, new_w, delta, max_delta, ata
    integer :: i, j, sweep

    weights = 0.0_wp
    do i = 1, nfreq
       f = 0.1_wp * f_transition * 100.0_wp**(real(i-1,wp)/real(nfreq-1,wp))
       if (f < f_transition) then
          target_q = q0
       else
          target_q = q0 * (f / f_transition)**gamma
       end if
       b(i) = 1.0_wp
       do j = 1, n
          x = 2.0_wp * pi * f * tau(j)
          denom = 1.0_wp + x*x
          a(i,j) = (target_q * x + 1.0_wp) / denom
       end do
    end do

    residual(1:nfreq) = b(1:nfreq)
    do sweep = 1, max_sweeps
       max_delta = 0.0_wp
       do j = 1, n
          old_w = weights(j)
          ata = dot_product(a(1:nfreq,j), a(1:nfreq,j))
          new_w = max(0.0_wp, old_w + dot_product(a(1:nfreq,j), residual(1:nfreq)) / &
               max(ata, tiny(1.0_wp)))
          delta = new_w - old_w
          if (delta /= 0.0_wp) residual(1:nfreq) = residual(1:nfreq) - delta * a(1:nfreq,j)
          weights(j) = new_w
          max_delta = max(max_delta, abs(delta))
       end do
       if (max_delta < tolerance) exit
    end do
  end subroutine nnls_fit_weights


  ! --- Coefficient builders ---

  subroutine build_nnls(n, qs, qp, gamma, f_transition, fmin, fmax, nfreq, tol, coeffs)
    integer, intent(in) :: n, nfreq
    real(wp), intent(in) :: qs, qp, gamma, f_transition, fmin, fmax, tol
    type(viscoelastic_coefficients), intent(inout) :: coeffs

    call log_spaced_tau(n, fmin, fmax, coeffs%tau)
    call nnls_fit_weights(n, qs, gamma, f_transition, coeffs%tau, nfreq, tol, coeffs%w_s)
    call nnls_fit_weights(n, qp, gamma, f_transition, coeffs%tau, nfreq, tol, coeffs%w_p)
  end subroutine build_nnls

  subroutine build_withers_table(qs, qp, gamma, f_transition, coeffs)
    use withers_tables, only : get_relaxation_times, get_withers_weights
    real(wp), intent(in) :: qs, qp, gamma, f_transition
    type(viscoelastic_coefficients), intent(inout) :: coeffs
    real(wp) :: tau8(8), ws8(8), wp8(8)

    call get_relaxation_times(gamma, tau8)
    coeffs%tau(1:8) = tau8 / f_transition
    call get_withers_weights(gamma, qs, ws8)
    coeffs%w_s(1:8) = ws8
    call get_withers_weights(gamma, qp, wp8)
    coeffs%w_p(1:8) = wp8
  end subroutine build_withers_table

  subroutine build_withers_method(n, qs, qp, gamma, f_transition, fmin, fmax, nfreq, tol, coeffs)
    integer, intent(in) :: n, nfreq
    real(wp), intent(in) :: qs, qp, gamma, f_transition, fmin, fmax, tol
    type(viscoelastic_coefficients), intent(inout) :: coeffs

    call withers_tau(n, gamma, f_transition, coeffs%tau)
    call nnls_fit_weights(n, qs, gamma, f_transition, coeffs%tau, nfreq, tol, coeffs%w_s)
    call nnls_fit_weights(n, qp, gamma, f_transition, coeffs%tau, nfreq, tol, coeffs%w_p)
  end subroutine build_withers_method

  pure subroutine build_fixed_q50(n, qs, qp, fmin, fmax, coeffs)
    integer, intent(in) :: n
    real(wp), intent(in) :: qs, qp, fmin, fmax
    type(viscoelastic_coefficients), intent(inout) :: coeffs
    real(wp) :: taumin, taumax
    integer :: k

    real(wp), parameter :: w4(4) = [1.549360_wp, 0.804277_wp, 0.887718_wp, 1.464160_wp]
    real(wp), parameter :: w8(8) = [1.50589707_wp, 0.0_wp, 0.52793567_wp, 0.53065494_wp, &
         0.32862132_wp, 0.64375916_wp, 0.0_wp, 1.32751442_wp]

    taumin = 1.0_wp / (2.0_wp * pi * fmax)
    taumax = 1.0_wp / (2.0_wp * pi * fmin)
    do k = 1, n
       coeffs%tau(k) = exp(log(taumin) + (2.0_wp*k - 1.0_wp) / (2.0_wp*n) * log(taumax/taumin))
    end do

    if (n == 4) then
       coeffs%w_s(1:4) = w4 / qs
       coeffs%w_p(1:4) = w4 / qp
    else
       coeffs%w_s(1:8) = w8 / qs
       coeffs%w_p(1:8) = w8 / qp
    end if
  end subroutine build_fixed_q50


  ! --- Diagnostics ---

  pure real(wp) function ve_realized_q(frequency, n, tau, strength) result(q)
    real(wp), intent(in) :: frequency, tau(:), strength(:)
    integer, intent(in) :: n
    real(wp) :: x, re, im
    integer :: l
    re = 1.0_wp; im = 0.0_wp
    do l = 1, n
       x = 2.0_wp * pi * frequency * tau(l)
       re = re - strength(l) / (1.0_wp + x*x)
       im = im + strength(l) * x / (1.0_wp + x*x)
    end do
    if (im <= tiny(1.0_wp)) then
       q = huge(1.0_wp)
    else
       q = abs(re / im)
    end if
  end function ve_realized_q

  pure subroutine ve_max_relative_error(q0, gamma, f_transition, n, tau, strength, &
       fmin, fmax, max_error)
    real(wp), intent(in) :: q0, gamma, f_transition, tau(:), strength(:), fmin, fmax
    integer, intent(in) :: n
    real(wp), intent(out) :: max_error
    integer, parameter :: nfreq = 256
    real(wp) :: f, target, realized
    integer :: i
    max_error = 0.0_wp
    do i = 1, nfreq
       f = fmin * (fmax/fmin)**(real(i-1,wp)/real(nfreq-1,wp))
       if (f < f_transition) then
          target = q0
       else
          target = q0 * (f / f_transition)**gamma
       end if
       realized = ve_realized_q(f, n, tau, strength)
       max_error = max(max_error, abs(realized/target - 1.0_wp))
    end do
  end subroutine ve_max_relative_error

  real(wp) function ve_params_dt_limit(params) result(dt_limit)
    type(viscoelastic_parameters), intent(in) :: params
    type(viscoelastic_coefficients) :: coeffs
    integer :: stat, b
    character(len=256) :: msg
    dt_limit = huge(1.0_wp)
    do b = 1, params%nblocks
       call build_viscoelastic_coefficients(params, b, coeffs, stat, msg)
       if (stat /= 0) then; dt_limit = 0.0_wp; return; end if
       dt_limit = min(dt_limit, ve_relaxation_dt_limit(coeffs))
    end do
  end function ve_params_dt_limit

  pure real(wp) function ve_relaxation_dt_limit(coeffs) result(dt_limit)
    type(viscoelastic_coefficients), intent(in) :: coeffs
    if (coeffs%n > 0) then
       dt_limit = 2.0_wp * minval(coeffs%tau(1:coeffs%n))
    else
       dt_limit = huge(1.0_wp)
    end if
  end function ve_relaxation_dt_limit

end module viscoelastic_model
