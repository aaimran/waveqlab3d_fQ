module anelastic_cq_model

  use common, only : wp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: cq_min_mechanisms = 4
  integer, parameter, public :: cq_max_mechanisms = 8

  type, public :: cq_parameters
     real(wp) :: Qs0(2) = -1.0_wp
     real(wp) :: Qp0(2) = -1.0_wp
     real(wp) :: fref = 1.0_wp
     real(wp) :: fmin = 0.05_wp
     real(wp) :: fmax = 20.0_wp
     integer :: n_mechanisms = 8
     integer :: nnls_samples = 256
     character(len=32) :: coefficient_policy = 'nnls-block-ps'
     character(len=32) :: nnls_objective = 'relative-q'
     real(wp) :: nnls_tolerance = 1.0e-10_wp
     real(wp) :: max_fit_error = 0.10_wp
  end type cq_parameters

  public :: read_cq_parameters, build_cq_coefficients
  public :: cq_relaxation_dt_limit, cq_max_relative_error

contains

  subroutine read_cq_parameters(infile, parameters, status, message)
    integer, intent(in) :: infile
    type(cq_parameters), intent(out) :: parameters
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(wp) :: Qs0(2), Qp0(2), fref, fmin, fmax, nnls_tolerance, max_fit_error
    integer :: n_mechanisms, nnls_samples, stat
    character(len=32) :: coefficient_policy, nnls_objective
    namelist /anelastic_cQ_list/ Qs0, Qp0, fref, fmin, fmax, n_mechanisms, &
         nnls_samples, coefficient_policy, nnls_objective, nnls_tolerance, max_fit_error

    Qs0=parameters%Qs0; Qp0=parameters%Qp0; fref=parameters%fref
    fmin=parameters%fmin; fmax=parameters%fmax
    n_mechanisms=parameters%n_mechanisms; nnls_samples=parameters%nnls_samples
    coefficient_policy=parameters%coefficient_policy
    nnls_objective=parameters%nnls_objective
    nnls_tolerance=parameters%nnls_tolerance; max_fit_error=parameters%max_fit_error
    status=0; message=''

    rewind(infile)
    read(infile,nml=anelastic_cQ_list,iostat=stat)
    if (stat /= 0) then
       status=1; message='response anelastic-cQ requires a valid &anelastic_cQ_list namelist'; return
    end if
    coefficient_policy=trim(adjustl(coefficient_policy))
    nnls_objective=trim(adjustl(nnls_objective))
    if (any(.not.ieee_is_finite(Qs0)) .or. any(.not.ieee_is_finite(Qp0)) .or. &
        any(Qs0 <= 0.0_wp) .or. any(Qp0 <= 0.0_wp)) then
       status=1; message='anelastic-cQ requires two finite, positive Qs0 and Qp0 values'; return
    end if
    if (.not.ieee_is_finite(fref) .or. .not.ieee_is_finite(fmin) .or. &
        .not.ieee_is_finite(fmax) .or. fmin <= 0.0_wp .or. fmax <= fmin .or. &
        fref < fmin .or. fref > fmax) then
       status=1; message='anelastic-cQ requires 0 < fmin <= fref <= fmax'; return
    end if
    if (n_mechanisms < cq_min_mechanisms .or. n_mechanisms > cq_max_mechanisms) then
       status=1; message='anelastic-cQ n_mechanisms must be one of 4, 5, 6, 7, or 8'; return
    end if
    if (nnls_samples < n_mechanisms) then
       status=1; message='anelastic-cQ nnls_samples must be at least n_mechanisms'; return
    end if
    if (nnls_objective /= 'relative-q') then
       status=1; message='anelastic-cQ nnls_objective must be relative-q'; return
    end if
    select case(coefficient_policy)
    case('fixed-q50')
       if (n_mechanisms /= 8) then
          status=1; message='anelastic-cQ fixed-q50 requires n_mechanisms=8'; return
       end if
       if (abs(fmin-0.05_wp) > 100.0_wp*epsilon(1.0_wp) .or. &
           abs(fmax-20.0_wp) > 100.0_wp*epsilon(1.0_wp)) then
          status=1; message='anelastic-cQ fixed-q50 requires fmin=0.05 Hz and fmax=20 Hz'; return
       end if
    case('nnls-shared','nnls-block','nnls-block-ps')
    case default
       status=1; message='unsupported anelastic-cQ coefficient_policy: '//trim(coefficient_policy); return
    end select
    if (.not.ieee_is_finite(nnls_tolerance) .or. nnls_tolerance <= 0.0_wp .or. &
        .not.ieee_is_finite(max_fit_error) .or. max_fit_error <= 0.0_wp) then
       status=1; message='anelastic-cQ NNLS tolerances must be finite and positive'; return
    end if

    parameters%Qs0=Qs0; parameters%Qp0=Qp0; parameters%fref=fref
    parameters%fmin=fmin; parameters%fmax=fmax
    parameters%n_mechanisms=n_mechanisms; parameters%nnls_samples=nnls_samples
    parameters%coefficient_policy=coefficient_policy
    parameters%nnls_objective=nnls_objective
    parameters%nnls_tolerance=nnls_tolerance; parameters%max_fit_error=max_fit_error
  end subroutine read_cq_parameters


  subroutine build_cq_coefficients(parameters, block_id, tau, weight_s, weight_p)
    type(cq_parameters), intent(in) :: parameters
    integer, intent(in) :: block_id
    real(wp), intent(out) :: tau(:), weight_s(:), weight_p(:)
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: taumin, taumax, target
    integer :: k, n
    n=parameters%n_mechanisms
    if (size(tau) /= n .or. size(weight_s) /= n .or. size(weight_p) /= n) &
         error stop 'invalid anelastic-cQ coefficient array extent'
    taumin=1.0_wp/(2.0_wp*pi*parameters%fmax)
    taumax=1.0_wp/(2.0_wp*pi*parameters%fmin)
    if (trim(parameters%coefficient_policy) == 'fixed-q50') then
       ! Preserve the relaxation times associated with the legacy fixed table.
       do k=1,n
          tau(k)=exp(log(taumin)+(2.0_wp*k-1.0_wp)/(2.0_wp*n)*log(taumax/taumin))
       end do
    else
       ! Include both requested band edges in an NNLS fit.  Placing every
       ! mechanism strictly inside the band leaves the edge response
       ! under-resolved, producing about 12% error even with eight mechanisms.
       do k=1,n
          tau(k)=exp(log(taumin)+real(k-1,wp)/real(n-1,wp)*log(taumax/taumin))
       end do
    end if
    select case(trim(parameters%coefficient_policy))
    case('fixed-q50')
       weight_s=[1.50589707_wp,0.0_wp,0.52793567_wp,0.53065494_wp, &
            0.32862132_wp,0.64375916_wp,0.0_wp,1.32751442_wp]
       weight_p=weight_s
    case('nnls-shared')
       target=0.25_wp*sum(parameters%Qs0+parameters%Qp0)
       call fit_cq_nnls(target,tau,parameters,weight_s)
       weight_p=weight_s
    case('nnls-block')
       target=0.5_wp*(parameters%Qs0(block_id)+parameters%Qp0(block_id))
       call fit_cq_nnls(target,tau,parameters,weight_s)
       weight_p=weight_s
    case('nnls-block-ps')
       call fit_cq_nnls(parameters%Qs0(block_id),tau,parameters,weight_s)
       call fit_cq_nnls(parameters%Qp0(block_id),tau,parameters,weight_p)
    end select
  end subroutine build_cq_coefficients


  subroutine fit_cq_nnls(target_q,tau,parameters,weights)
    real(wp), intent(in) :: target_q,tau(:)
    type(cq_parameters), intent(in) :: parameters
    real(wp), intent(out) :: weights(:)
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp), allocatable :: a(:,:), gradient(:), trial(:), residual(:)
    real(wp) :: f,omega,x,step,lipschitz,change
    integer :: i,k,iter,n,ns
    n=size(tau); ns=parameters%nnls_samples
    allocate(a(ns,n),gradient(n),trial(n),residual(ns))
    do i=1,ns
       f=parameters%fmin*(parameters%fmax/parameters%fmin)**(real(i-1,wp)/real(ns-1,wp))
       omega=2.0_wp*pi*f
       do k=1,n
          x=omega*tau(k)
          a(i,k)=(x+1.0_wp/target_q)/(1.0_wp+x*x)
       end do
    end do
    ! Projected gradient for min ||A*w-1||_2 with w >= 0.  The Frobenius
    ! norm is a conservative Lipschitz bound for A^T*A.
    lipschitz=max(sum(a*a),tiny(1.0_wp)); step=0.95_wp/lipschitz
    weights=1.0_wp/real(n,wp)
    do iter=1,50000
       residual=matmul(a,weights)-1.0_wp
       gradient=matmul(transpose(a),residual)
       trial=max(weights-step*gradient,0.0_wp)
       change=maxval(abs(trial-weights)); weights=trial
       if (change <= parameters%nnls_tolerance) exit
    end do
    deallocate(a,gradient,trial,residual)
  end subroutine fit_cq_nnls


  pure subroutine cq_max_relative_error(target_q,tau,weight,fmin,fmax,max_error)
    real(wp), intent(in) :: target_q,tau(:),weight(:),fmin,fmax
    real(wp), intent(out) :: max_error
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: f,omega,x,mr,mi,q
    integer :: i,k
    max_error=0.0_wp
    do i=1,256
       f=fmin*(fmax/fmin)**(real(i-1,wp)/255.0_wp); omega=2.0_wp*pi*f
       mr=1.0_wp; mi=0.0_wp
       do k=1,size(tau)
          x=omega*tau(k); mr=mr-(weight(k)/target_q)/(1.0_wp+x*x)
          mi=mi+(weight(k)/target_q)*x/(1.0_wp+x*x)
       end do
       if (mi <= tiny(1.0_wp)) then; max_error=huge(1.0_wp); return; end if
       q=abs(mr/mi); max_error=max(max_error,abs(q/target_q-1.0_wp))
    end do
  end subroutine cq_max_relative_error


  pure real(wp) function cq_relaxation_dt_limit(parameters) result(limit)
    type(cq_parameters), intent(in) :: parameters
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: taumin,taumax
    taumin=1.0_wp/(2.0_wp*pi*parameters%fmax)
    taumax=1.0_wp/(2.0_wp*pi*parameters%fmin)
    if (trim(parameters%coefficient_policy) == 'fixed-q50') then
       limit=2.0_wp*exp(log(taumin)+1.0_wp/(2.0_wp*parameters%n_mechanisms)*log(taumax/taumin))
    else
       limit=2.0_wp*taumin
    end if
  end function cq_relaxation_dt_limit

end module anelastic_cq_model
