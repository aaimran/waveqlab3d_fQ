module anelastic_fq8_model

  use common, only : wp
  use withers_tables, only : get_relaxation_times, get_withers_weights
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: fq8_nmechanisms = 8
  real(wp), parameter, public :: fq8_minimum_q = 15.0_wp

  type, public :: fq8_parameters
     character(len=32) :: coefficient_method = 'conventional-nnls'
     character(len=32) :: weight_policy = 'table-exact'
     integer :: coarse_grain = -1
     real(wp) :: Qs0 = -1.0_wp
     real(wp) :: Qp0 = -1.0_wp
     real(wp) :: gamma = -1.0_wp
     real(wp) :: f_transition = -1.0_wp
     real(wp) :: fref = -1.0_wp
  end type fq8_parameters

  public :: read_fq8_parameters, build_fq8_coefficients
  public :: fq8_realized_q, fq8_max_relative_error
  public :: fq8_relaxation_dt_limit
  public :: fq8_harmonic_response, fq8_effective_q
  public :: fq8_common_modulus_scale, fq8_phase_velocity_ratio

contains

  subroutine read_fq8_parameters(infile, parameters, status, message)
    integer, intent(in) :: infile
    type(fq8_parameters), intent(out) :: parameters
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    character(len=32) :: coefficient_method, weight_policy
    integer :: coarse_grain
    real(wp) :: Qs0, Qp0, gamma, f_transition, fref
    integer :: stat
    namelist /anelastic_fQ8_list/ coefficient_method, weight_policy, coarse_grain, &
         Qs0, Qp0, gamma, f_transition, fref

    coefficient_method=parameters%coefficient_method
    weight_policy=parameters%weight_policy
    coarse_grain=parameters%coarse_grain
    Qs0=parameters%Qs0; Qp0=parameters%Qp0; gamma=parameters%gamma
    f_transition=parameters%f_transition; fref=parameters%fref
    status=0; message=''
    rewind(infile)
    read(infile,nml=anelastic_fQ8_list,iostat=stat)
    if (stat /= 0) then
       status=1
       message='anelastic-fQ8 requires a valid &anelastic_fQ8_list (legacy c is unsupported)'
       return
    end if
    coefficient_method=trim(adjustl(coefficient_method))
    weight_policy=trim(adjustl(weight_policy))
    if (coefficient_method /= 'withers-2015' .and. &
        coefficient_method /= 'conventional-nnls') then
       status=1
       message='anelastic-fQ8 coefficient_method must be conventional-nnls or withers-2015'
       return
    end if
    if (weight_policy /= 'table-exact' .and. weight_policy /= 'nonnegative-refit') then
       status=1
       message='anelastic-fQ8 weight_policy must be table-exact or nonnegative-refit'
       return
    end if
    if (coefficient_method /= 'withers-2015' .and. weight_policy /= 'table-exact') then
       status=1
       message='anelastic-fQ8 nonnegative-refit applies only to coefficient_method=withers-2015'
       return
    end if
    if (coarse_grain == -1) then
       if (coefficient_method == 'withers-2015') then
          coarse_grain=2
       else
          coarse_grain=0
       end if
    end if
    if (coarse_grain /= 0 .and. coarse_grain /= 2) then
       status=1
       message='anelastic-fQ8 coarse_grain must be 0 or 2'
       return
    end if
    if (coarse_grain == 0 .and. coefficient_method == 'withers-2015') then
       status=1
       message='anelastic-fQ8 coarse_grain=0 requires coefficient_method=conventional-nnls; '// &
            'raw withers-2015 strengths are for the 2x2x2 coarse layout'
       return
    end if
    if (.not.ieee_is_finite(Qs0) .or. .not.ieee_is_finite(Qp0) .or. &
        Qs0 < fq8_minimum_q .or. Qp0 < fq8_minimum_q) then
       status=1; message='anelastic-fQ8 requires finite Qs0 and Qp0 >= 15'; return
    end if
    if (.not.ieee_is_finite(gamma) .or. gamma < 0.0_wp .or. gamma > 0.9_wp) then
       status=1; message='anelastic-fQ8 gamma must be finite and in [0,0.9]'; return
    end if
    if (.not.ieee_is_finite(f_transition) .or. f_transition <= 0.0_wp) then
       status=1; message='anelastic-fQ8 f_transition must be finite and positive'; return
    end if
    if (.not.ieee_is_finite(fref) .or. fref <= 0.0_wp) then
       status=1; message='anelastic-fQ8 fref must be finite and positive'; return
    end if
    parameters%coefficient_method=coefficient_method
    parameters%weight_policy=weight_policy
    parameters%coarse_grain=coarse_grain
    parameters%Qs0=Qs0; parameters%Qp0=Qp0; parameters%gamma=gamma
    parameters%f_transition=f_transition; parameters%fref=fref
  end subroutine read_fq8_parameters

  subroutine build_fq8_coefficients(parameters,tau,strength_s,strength_p,status,message)
    type(fq8_parameters), intent(in) :: parameters
    real(wp), intent(out) :: tau(fq8_nmechanisms), strength_s(fq8_nmechanisms), &
         strength_p(fq8_nmechanisms)
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    if (parameters%coarse_grain == 0 .and. parameters%coefficient_method == 'withers-2015') then
       status=1
       message='anelastic-fQ8 coarse_grain=0 requires coefficient_method=conventional-nnls; '// &
            'raw withers-2015 strengths are for the 2x2x2 coarse layout'
       return
    end if
    call get_relaxation_times(parameters%gamma,tau)
    tau=tau/parameters%f_transition
    if (parameters%coefficient_method == 'conventional-nnls') then
       call fit_conventional_strengths(parameters%Qs0,parameters%gamma, &
            parameters%f_transition,tau,strength_s)
       call fit_conventional_strengths(parameters%Qp0,parameters%gamma, &
            parameters%f_transition,tau,strength_p)
    else
       call get_withers_weights(parameters%gamma,parameters%Qs0,strength_s)
       call get_withers_weights(parameters%gamma,parameters%Qp0,strength_p)
       ! Published weights are w_k=N*lambda_k and are used directly with one
       ! mechanism per node in the deterministic period-two coarse layout.
       if (parameters%weight_policy == 'nonnegative-refit') then
          call refit_nonnegative_at_reference(parameters%fref,tau,strength_s,status,message)
          if (status /= 0) return
          call refit_nonnegative_at_reference(parameters%fref,tau,strength_p,status,message)
          if (status /= 0) return
       end if
    end if
    status=0; message=''
    if (.not.all(ieee_is_finite(tau)) .or. any(tau <= 0.0_wp) .or. &
        .not.all(ieee_is_finite(strength_s)) .or. &
        .not.all(ieee_is_finite(strength_p))) then
       status=1; message='anelastic-fQ8 coefficient construction produced invalid values'; return
    end if
    if (sum(strength_s) >= 1.0_wp .or. sum(strength_p) >= 1.0_wp) then
       status=1; message='anelastic-fQ8 coefficients leave a non-positive relaxed modulus'; return
    end if
  end subroutine build_fq8_coefficients

  pure subroutine refit_nonnegative_at_reference(fref,tau,strength,status,message)
    real(wp), intent(in) :: fref,tau(fq8_nmechanisms)
    real(wp), intent(inout) :: strength(fq8_nmechanisms)
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(wp) :: raw(fq8_nmechanisms),target_q,lo,hi,mid,qmid,max_factor
    integer :: iteration

    status=0; message=''
    if (all(strength >= 0.0_wp)) return
    raw=strength
    target_q=fq8_effective_q(fref,tau,raw)
    strength=max(strength,0.0_wp)
    if (sum(strength) <= tiny(1.0_wp)) then
       status=1; message='nonnegative-refit has no positive strengths to rescale'; return
    end if
    max_factor=(1.0_wp-100.0_wp*epsilon(1.0_wp))/sum(strength)
    lo=0.0_wp; hi=min(2.0_wp,max_factor)
    if (fq8_effective_q(fref,tau,hi*strength) > target_q) then
       status=1; message='nonnegative-refit cannot bracket reference effective Q'; return
    end if
    do iteration=1,100
       mid=0.5_wp*(lo+hi)
       qmid=fq8_effective_q(fref,tau,mid*strength)
       if (qmid > target_q) then
          lo=mid
       else
          hi=mid
       end if
    end do
    strength=0.5_wp*(lo+hi)*strength
  end subroutine refit_nonnegative_at_reference

  pure subroutine fit_conventional_strengths(q0,gamma,f_transition,tau,strength)
    real(wp), intent(in) :: q0,gamma,f_transition,tau(fq8_nmechanisms)
    real(wp), intent(out) :: strength(fq8_nmechanisms)
    integer, parameter :: nfreq=256, max_sweeps=20000
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: a(nfreq,fq8_nmechanisms),b(nfreq),residual(nfreq)
    real(wp) :: f,target_q,x,denom,old,new,delta,max_delta
    integer :: i,j,sweep

    do i=1,nfreq
       f=0.1_wp*f_transition*100.0_wp**(real(i-1,wp)/real(nfreq-1,wp))
       if (f < f_transition) then
          target_q=q0
       else
          target_q=q0*(f/f_transition)**gamma
       end if
       b(i)=1.0_wp
       do j=1,fq8_nmechanisms
          x=2.0_wp*pi*f*tau(j); denom=1.0_wp+x*x
          ! Im(M)-Re(M)/Q=0 is linear in the co-located strengths.
          a(i,j)=(target_q*x+1.0_wp)/denom
       end do
    end do
    strength=0.0_wp; residual=b
    do sweep=1,max_sweeps
       max_delta=0.0_wp
       do j=1,fq8_nmechanisms
          old=strength(j)
          new=max(0.0_wp,old+dot_product(a(:,j),residual)/ &
                  max(dot_product(a(:,j),a(:,j)),tiny(1.0_wp)))
          delta=new-old
          if (delta /= 0.0_wp) residual=residual-delta*a(:,j)
          strength(j)=new; max_delta=max(max_delta,abs(delta))
       end do
       if (max_delta < 1.0e-13_wp) exit
    end do
  end subroutine fit_conventional_strengths

  pure real(wp) function fq8_realized_q(frequency,tau,strength) result(q)
    real(wp), intent(in) :: frequency, tau(:), strength(:)
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: x, re, im
    integer :: l
    re=1.0_wp; im=0.0_wp
    do l=1,size(tau)
       x=2.0_wp*pi*frequency*tau(l)
       re=re-strength(l)/(1.0_wp+x*x)
       im=im+strength(l)*x/(1.0_wp+x*x)
    end do
    if (im <= tiny(1.0_wp)) then
       q=huge(1.0_wp)
    else
       q=abs(re/im)
    end if
  end function fq8_realized_q

  pure complex(wp) function fq8_harmonic_response(frequency,tau,strength) result(response)
    real(wp), intent(in) :: frequency, tau(:), strength(:)
    real(wp), parameter :: pi=3.141592653589793_wp
    complex(wp) :: local_response, inverse_sum
    integer :: l

    inverse_sum=cmplx(0.0_wp,0.0_wp,kind=wp)
    do l=1,size(tau)
       local_response=1.0_wp-strength(l)/ &
            cmplx(1.0_wp,2.0_wp*pi*frequency*tau(l),kind=wp)
       inverse_sum=inverse_sum+1.0_wp/local_response
    end do
    response=real(size(tau),wp)/inverse_sum
  end function fq8_harmonic_response

  pure real(wp) function fq8_effective_q(frequency,tau,strength) result(q)
    real(wp), intent(in) :: frequency,tau(:),strength(:)
    complex(wp) :: response
    response=fq8_harmonic_response(frequency,tau,strength)
    if (aimag(response) <= tiny(1.0_wp)) then
       q=huge(1.0_wp)
    else
       q=real(response,wp)/aimag(response)
    end if
  end function fq8_effective_q

  pure real(wp) function fq8_common_modulus_scale(fref,tau,strength) result(scale)
    real(wp), intent(in) :: fref,tau(:),strength(:)
    complex(wp) :: response
    response=fq8_harmonic_response(fref,tau,strength)
    scale=real(1.0_wp/sqrt(response),wp)**2
  end function fq8_common_modulus_scale

  pure real(wp) function fq8_phase_velocity_ratio(frequency,fref,tau,strength) result(ratio)
    real(wp), intent(in) :: frequency,fref,tau(:),strength(:)
    real(wp) :: scale
    complex(wp) :: response
    scale=fq8_common_modulus_scale(fref,tau,strength)
    response=scale*fq8_harmonic_response(frequency,tau,strength)
    ratio=1.0_wp/real(1.0_wp/sqrt(response),wp)
  end function fq8_phase_velocity_ratio

  pure subroutine fq8_max_relative_error(q0,gamma,f_transition,tau,strength, &
       fmin,fmax,max_error)
    real(wp), intent(in) :: q0,gamma,f_transition,tau(:),strength(:),fmin,fmax
    real(wp), intent(out) :: max_error
    integer, parameter :: nfreq=256
    real(wp) :: f,target,realized
    integer :: i
    max_error=0.0_wp
    do i=1,nfreq
       f=fmin*(fmax/fmin)**(real(i-1,wp)/real(nfreq-1,wp))
       if (f < f_transition) then
          target=q0
       else
          target=q0*(f/f_transition)**gamma
       end if
       realized=fq8_realized_q(f,tau,strength)
       max_error=max(max_error,abs(realized/target-1.0_wp))
    end do
  end subroutine fq8_max_relative_error

  real(wp) function fq8_relaxation_dt_limit(parameters) result(limit)
    type(fq8_parameters), intent(in) :: parameters
    real(wp) :: tau(8),ss(8),sp(8)
    integer :: status
    character(len=128) :: message
    call build_fq8_coefficients(parameters,tau,ss,sp,status,message)
    if (status == 0) then
       limit=2.0_wp*minval(tau)
    else
       limit=0.0_wp
    end if
  end function fq8_relaxation_dt_limit

end module anelastic_fq8_model
