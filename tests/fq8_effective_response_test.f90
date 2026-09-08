program fq8_effective_response_test
  use common, only : wp
  use anelastic_fq8_model, only : fq8_parameters, build_fq8_coefficients, &
       fq8_effective_q, fq8_common_modulus_scale, fq8_phase_velocity_ratio, &
       read_fq8_parameters
  implicit none
  type(fq8_parameters) :: parameters
  real(wp) :: tau(8),strength_s(8),strength_p(8),q,scale,ratio
  real(wp) :: tau_refit(8),strength_s_refit(8),strength_p_refit(8)
  real(wp) :: frequency,max_q_delta,max_velocity_delta
  real(wp) :: rho,vs_ref,vp_ref,mu_nodes(8),lambda_nodes(8)
  integer :: status
  character(len=256) :: message

  call test_fq8_parameter_parsing()

  parameters%coefficient_method='withers-2015'
  parameters%coarse_grain=2
  parameters%Qs0=50.0_wp; parameters%Qp0=50.0_wp
  parameters%gamma=0.6_wp; parameters%f_transition=1.0_wp
  parameters%fref=1.0_wp
  call build_fq8_coefficients(parameters,tau,strength_s,strength_p,status,message)
  if (status /= 0) then
       write(*,'(A)') trim(message)
       stop 1
  end if

  q=fq8_effective_q(1.0_wp,tau,strength_s)
  scale=fq8_common_modulus_scale(1.0_wp,tau,strength_s)
  ratio=fq8_phase_velocity_ratio(1.0_wp,1.0_wp,tau,strength_s)
  call require_close('effective Q at 1 Hz',q,54.8324966843126_wp,1.0e-11_wp)
  call require_close('common modulus scale',scale,1.01727165422423_wp,1.0e-11_wp)
  call require_close('reference phase velocity ratio',ratio,1.0_wp,1.0e-13_wp)
  call require_close('P/S equal-Q coefficient identity',maxval(abs(strength_s-strength_p)), &
       0.0_wp,1.0e-15_wp)
  call require_close('published negative strength',strength_s(1),-0.000688_wp,1.0e-14_wp)

  ! The coarse Withers cell must have one common elastic modulus.  Mechanism
  ! identity changes only the memory-variable coefficients, never impedance.
  rho=2.7_wp; vs_ref=3.464_wp; vp_ref=6.0_wp
  mu_nodes=rho*vs_ref**2*scale
  lambda_nodes=rho*vp_ref**2*scale-2.0_wp*mu_nodes
  call require_close('coarse-cell mu spread',maxval(mu_nodes)-minval(mu_nodes), &
       0.0_wp,1.0e-14_wp)
  call require_close('coarse-cell lambda spread',maxval(lambda_nodes)-minval(lambda_nodes), &
       0.0_wp,1.0e-14_wp)
  call require_close('S reference velocity', &
       sqrt(mu_nodes(1)/(rho*scale))*ratio,vs_ref,1.0e-13_wp)
  call require_close('P reference velocity', &
       sqrt((lambda_nodes(1)+2.0_wp*mu_nodes(1))/(rho*scale))*ratio, &
       vp_ref,1.0e-13_wp)

  parameters%weight_policy='nonnegative-refit'
  call build_fq8_coefficients(parameters,tau_refit,strength_s_refit, &
       strength_p_refit,status,message)
  if (status /= 0) then
     write(*,'(A)') trim(message)
     stop 1
  end if
  if (any(strength_s_refit < 0.0_wp) .or. any(strength_p_refit < 0.0_wp)) &
       error stop 'nonnegative-refit retained a negative strength'
  call require_close('refit preserves Qeff at fref', &
       fq8_effective_q(1.0_wp,tau_refit,strength_s_refit),q,1.0e-11_wp)
  call require_close('refit P/S coefficient identity', &
       maxval(abs(strength_s_refit-strength_p_refit)),0.0_wp,1.0e-15_wp)
  max_q_delta=0.0_wp; max_velocity_delta=0.0_wp
  do status=0,511
     frequency=0.05_wp*(20.0_wp/0.05_wp)**(real(status,wp)/511.0_wp)
     max_q_delta=max(max_q_delta,abs(fq8_effective_q(frequency,tau_refit,strength_s_refit)/ &
          fq8_effective_q(frequency,tau,strength_s)-1.0_wp))
     max_velocity_delta=max(max_velocity_delta,abs( &
          fq8_phase_velocity_ratio(frequency,1.0_wp,tau_refit,strength_s_refit)/ &
          fq8_phase_velocity_ratio(frequency,1.0_wp,tau,strength_s)-1.0_wp))
  end do
  if (max_q_delta > 0.0128_wp) error stop 'nonnegative-refit Q delta exceeds 1.28 percent'
  if (max_velocity_delta > 1.3e-5_wp) &
       error stop 'nonnegative-refit phase-velocity delta exceeds 1.3e-5'

  parameters%coefficient_method='withers-2015'
  parameters%coarse_grain=0
  call build_fq8_coefficients(parameters,tau,strength_s,strength_p,status,message)
  if (status == 0) error stop 'coarse_grain=0 accepted raw withers-2015 coefficients'

  write(*,'(A)') 'fq8 effective-response test passed'

contains
  subroutine test_fq8_parameter_parsing()
    type(fq8_parameters) :: parsed
    integer :: unit,status
    character(len=256) :: message

    open(newunit=unit,status='scratch',action='readwrite')
    write(unit,'(A)') '&anelastic_fQ8_list'
    write(unit,'(A)') " coefficient_method='withers-2015',"
    write(unit,'(A)') " weight_policy='table-exact',"
    write(unit,'(A)') ' Qs0=50d0, Qp0=50d0, gamma=0d0, f_transition=1d0, fref=1d0 /'
    rewind(unit)
      call read_fq8_parameters(unit,parsed,status,message)
      close(unit)
      if (status /= 0) then
           write(*,'(A)') trim(message)
           stop 1
      end if
    if (parsed%coarse_grain /= 2) error stop 'missing coarse_grain did not default to 2'

    open(newunit=unit,status='scratch',action='readwrite')
    write(unit,'(A)') '&anelastic_fQ8_list'
    write(unit,'(A)') " coefficient_method='conventional-nnls',"
    write(unit,'(A)') ' Qs0=50d0, Qp0=50d0, gamma=0d0, f_transition=1d0, fref=1d0 /'
    rewind(unit)
      call read_fq8_parameters(unit,parsed,status,message)
      close(unit)
      if (status /= 0) then
           write(*,'(A)') trim(message)
           stop 1
      end if
    if (parsed%coarse_grain /= 0) &
         error stop 'missing coarse_grain did not preserve conventional full mode'

    open(newunit=unit,status='scratch',action='readwrite')
    write(unit,'(A)') '&anelastic_fQ8_list'
    write(unit,'(A)') " coefficient_method='conventional-nnls',"
    write(unit,'(A)') ' coarse_grain=0,'
    write(unit,'(A)') ' Qs0=50d0, Qp0=50d0, gamma=0d0, f_transition=1d0, fref=1d0 /'
    rewind(unit)
      call read_fq8_parameters(unit,parsed,status,message)
      close(unit)
      if (status /= 0) then
           write(*,'(A)') trim(message)
           stop 1
      end if
    if (parsed%coarse_grain /= 0) error stop 'coarse_grain=0 did not parse'

    open(newunit=unit,status='scratch',action='readwrite')
    write(unit,'(A)') '&anelastic_fQ8_list'
    write(unit,'(A)') " coefficient_method='withers-2015',"
    write(unit,'(A)') ' coarse_grain=0,'
    write(unit,'(A)') ' Qs0=50d0, Qp0=50d0, gamma=0d0, f_transition=1d0, fref=1d0 /'
    rewind(unit)
    call read_fq8_parameters(unit,parsed,status,message)
    close(unit)
    if (status == 0) error stop 'parser accepted coarse_grain=0 with withers-2015'

    open(newunit=unit,status='scratch',action='readwrite')
    write(unit,'(A)') '&anelastic_fQ8_list'
    write(unit,'(A)') " coefficient_method='conventional-nnls',"
    write(unit,'(A)') ' coarse_grain=1,'
    write(unit,'(A)') ' Qs0=50d0, Qp0=50d0, gamma=0d0, f_transition=1d0, fref=1d0 /'
    rewind(unit)
    call read_fq8_parameters(unit,parsed,status,message)
    close(unit)
    if (status == 0) error stop 'parser accepted invalid coarse_grain=1'
  end subroutine test_fq8_parameter_parsing

  subroutine require_close(label,value,expected,tolerance)
    character(*), intent(in) :: label
    real(wp), intent(in) :: value,expected,tolerance
    if (abs(value-expected) > tolerance) then
       write(*,'(A,2ES24.15,A,ES12.4)') trim(label)//': ',value,expected, &
            ' tolerance=',tolerance
       error stop 1
    end if
  end subroutine require_close
end program fq8_effective_response_test
