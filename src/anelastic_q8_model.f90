module anelastic_q8_model

  use common, only : wp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: q8_nmechanisms = 8
  real(wp), parameter, public :: q8_standard_fmin = 0.05_wp
  real(wp), parameter, public :: q8_standard_fmax = 20.0_wp

  type, public :: q8_parameters
     character(len=32) :: weight_method = 'fixed-q50'
     real(wp) :: Qs0 = -1.0_wp
     real(wp) :: Qp0 = -1.0_wp
     real(wp) :: fref = 1.0_wp
     real(wp) :: fmin = q8_standard_fmin
     real(wp) :: fmax = q8_standard_fmax
  end type q8_parameters

  public :: read_q8_parameters, build_q8_fixed_coefficients
  public :: q8_max_relative_error
  public :: q8_relaxation_dt_limit

contains

  subroutine read_q8_parameters(infile, parameters, status, message)
    integer, intent(in) :: infile
    type(q8_parameters), intent(out) :: parameters
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=32) :: weight_method
    real(wp) :: Qs0, Qp0, fref, fmin, fmax
    integer :: stat
    real(wp), parameter :: band_tolerance = 100.0_wp * epsilon(1.0_wp)
    namelist /anelastic_Q8_list/ weight_method, Qs0, Qp0, fref, fmin, fmax

    weight_method = parameters%weight_method
    Qs0 = parameters%Qs0
    Qp0 = parameters%Qp0
    fref = parameters%fref
    fmin = parameters%fmin
    fmax = parameters%fmax
    status = 0
    message = ''

    rewind(infile)
    read(infile, nml=anelastic_Q8_list, iostat=stat)
    if (stat /= 0) then
       status = 1
       message = 'response anelastic-Q8 requires a valid &anelastic_Q8_list namelist'
       return
    end if

    weight_method = trim(adjustl(weight_method))

    if (.not.ieee_is_finite(Qs0) .or. .not.ieee_is_finite(Qp0) .or. &
        Qs0 <= 0.0_wp .or. Qp0 <= 0.0_wp) then
       status = 1
       message = 'anelastic-Q8 requires finite, positive Qs0 and Qp0'
       return
    end if

    if (.not.ieee_is_finite(fref) .or. fref <= 0.0_wp) then
       status = 1
       message = 'anelastic-Q8 fref must be finite and positive'
       return
    end if
    if (.not.ieee_is_finite(fmin) .or. .not.ieee_is_finite(fmax) .or. &
        fmin <= 0.0_wp .or. fmax <= fmin) then
       status = 1
       message = 'anelastic-Q8 requires finite frequencies with 0 < fmin < fmax'
       return
    end if

    select case (weight_method)
    case ('fixed-q50')
       if (abs(fmin-q8_standard_fmin) > band_tolerance .or. &
           abs(fmax-q8_standard_fmax) > band_tolerance) then
          status = 1
          message = 'fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz'
          return
       end if
    case default
       status = 1
       message = 'unsupported anelastic-Q8 weight_method in this increment: '//trim(weight_method)
       return
    end select

    parameters%weight_method = weight_method
    parameters%Qs0 = Qs0
    parameters%Qp0 = Qp0
    parameters%fref = fref
    parameters%fmin = fmin
    parameters%fmax = fmax
  end subroutine read_q8_parameters


  pure subroutine build_q8_fixed_coefficients(parameters, tau, weight)
    type(q8_parameters), intent(in) :: parameters
    real(wp), intent(out) :: tau(q8_nmechanisms), weight(q8_nmechanisms)
    real(wp) :: taumin, taumax
    real(wp), parameter :: pi = 3.141592653589793_wp
    integer :: k

    taumin = 1.0_wp / (2.0_wp*pi*parameters%fmax)
    taumax = 1.0_wp / (2.0_wp*pi*parameters%fmin)
    do k = 1, q8_nmechanisms
       tau(k) = exp(log(taumin) + (2.0_wp*k-1.0_wp) / &
                (2.0_wp*q8_nmechanisms) * log(taumax/taumin))
    end do

    weight = [1.50589707_wp, 0.0_wp, 0.52793567_wp, 0.53065494_wp, &
              0.32862132_wp, 0.64375916_wp, 0.0_wp, 1.32751442_wp]
  end subroutine build_q8_fixed_coefficients


  pure subroutine q8_max_relative_error(target_Q, tau, weight, fmin, fmax, max_error)
    real(wp), intent(in) :: target_Q, tau(:), weight(:), fmin, fmax
    real(wp), intent(out) :: max_error
    integer, parameter :: nfreq = 256
    real(wp), parameter :: pi = 3.141592653589793_wp
    real(wp) :: frequency, omega, x, modulus_real, modulus_imag, realized_Q
    integer :: i, l

    max_error = 0.0_wp
    do i = 1, nfreq
       frequency = fmin * (fmax/fmin)**(real(i-1,wp)/real(nfreq-1,wp))
       omega = 2.0_wp*pi*frequency
       modulus_real = 1.0_wp
       modulus_imag = 0.0_wp
       do l = 1, size(tau)
          x = omega*tau(l)
          modulus_real = modulus_real - (weight(l)/target_Q)/(1.0_wp+x*x)
          modulus_imag = modulus_imag + (weight(l)/target_Q)*x/(1.0_wp+x*x)
       end do
       if (modulus_imag <= tiny(1.0_wp)) then
          max_error = huge(1.0_wp)
          return
       end if
       realized_Q = abs(modulus_real/modulus_imag)
       max_error = max(max_error, abs(realized_Q/target_Q-1.0_wp))
    end do
  end subroutine q8_max_relative_error


  pure real(wp) function q8_relaxation_dt_limit(parameters) result(dt_limit)
    type(q8_parameters), intent(in) :: parameters
    real(wp) :: tau(q8_nmechanisms), weight(q8_nmechanisms)

    call build_q8_fixed_coefficients(parameters, tau, weight)
    ! Conservative negative-real-axis bound.  The configured five-stage RK4
    ! scheme has a wider stability interval, so 2*tau_min retains margin.
    dt_limit = 2.0_wp*minval(tau)
  end function q8_relaxation_dt_limit

end module anelastic_q8_model
