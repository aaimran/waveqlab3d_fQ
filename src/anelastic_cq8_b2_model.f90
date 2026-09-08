module anelastic_cq8_b2_model

  use common, only : wp
  use anelastic_q8_model, only : q8_parameters, q8_standard_fmin, q8_standard_fmax, &
       q8_relaxation_dt_limit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  type, public :: cq8_b2_parameters
     type(q8_parameters) :: block(2)
  end type cq8_b2_parameters

  public :: read_cq8_b2_parameters, cq8_b2_relaxation_dt_limit

contains

  subroutine read_cq8_b2_parameters(infile, parameters, status, message)
    integer, intent(in) :: infile
    type(cq8_b2_parameters), intent(out) :: parameters
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    character(len=32) :: weight_method
    real(wp) :: Qs0(2), Qp0(2), fref, fmin, fmax
    real(wp), parameter :: band_tolerance = 100.0_wp * epsilon(1.0_wp)
    integer :: stat, i
    namelist /anelastic_cQ8_b2_list/ weight_method, Qs0, Qp0, fref, fmin, fmax

    weight_method = 'fixed-q50'
    Qs0 = -1.0_wp
    Qp0 = -1.0_wp
    fref = 1.0_wp
    fmin = q8_standard_fmin
    fmax = q8_standard_fmax
    status = 0
    message = ''

    rewind(infile)
    read(infile, nml=anelastic_cQ8_b2_list, iostat=stat)
    if (stat /= 0) then
       status = 1
       message = 'response anelastic-cQ8-b2 requires a valid &anelastic_cQ8_b2_list namelist'
       return
    end if

    weight_method = trim(adjustl(weight_method))
    if (any(.not.ieee_is_finite(Qs0)) .or. any(.not.ieee_is_finite(Qp0)) .or. &
        any(Qs0 <= 0.0_wp) .or. any(Qp0 <= 0.0_wp)) then
       status = 1
       message = 'anelastic-cQ8-b2 requires two finite, positive Qs0 and Qp0 values'
       return
    end if
    if (.not.ieee_is_finite(fref) .or. fref <= 0.0_wp) then
       status = 1
       message = 'anelastic-cQ8-b2 fref must be finite and positive'
       return
    end if
    if (.not.ieee_is_finite(fmin) .or. .not.ieee_is_finite(fmax) .or. &
        fmin <= 0.0_wp .or. fmax <= fmin) then
       status = 1
       message = 'anelastic-cQ8-b2 requires finite frequencies with 0 < fmin < fmax'
       return
    end if
    if (weight_method /= 'fixed-q50') then
       status = 1
       message = 'unsupported anelastic-cQ8-b2 weight_method: '//trim(weight_method)
       return
    end if
    if (abs(fmin-q8_standard_fmin) > band_tolerance .or. &
        abs(fmax-q8_standard_fmax) > band_tolerance) then
       status = 1
       message = 'anelastic-cQ8-b2 fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz'
       return
    end if

    do i = 1, 2
       parameters%block(i)%weight_method = weight_method
       parameters%block(i)%Qs0 = Qs0(i)
       parameters%block(i)%Qp0 = Qp0(i)
       parameters%block(i)%fref = fref
       parameters%block(i)%fmin = fmin
       parameters%block(i)%fmax = fmax
    end do
  end subroutine read_cq8_b2_parameters


  pure real(wp) function cq8_b2_relaxation_dt_limit(parameters) result(dt_limit)
    type(cq8_b2_parameters), intent(in) :: parameters
    dt_limit = min(q8_relaxation_dt_limit(parameters%block(1)), &
                   q8_relaxation_dt_limit(parameters%block(2)))
  end function cq8_b2_relaxation_dt_limit

end module anelastic_cq8_b2_model
