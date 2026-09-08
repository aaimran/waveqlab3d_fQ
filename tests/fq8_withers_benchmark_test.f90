program fq8_withers_benchmark_test
  use common, only : wp
  use anelastic_fq8_model, only : fq8_parameters, build_fq8_coefficients, &
       fq8_effective_q, fq8_common_modulus_scale, fq8_phase_velocity_ratio
  implicit none

  call test_elastic_halfspace()
  call test_constant_q_halfspace()
  call test_powerlaw_halfspace()
  call test_layered_model_Q20()
  call test_layered_model_Q210()

  write(*,'(A)') 'fq8 Withers benchmark coefficient test passed'

contains

  subroutine test_elastic_halfspace()
    real(wp) :: tau(8), strength_s(8), strength_p(8)
    type(fq8_parameters) :: p
    integer :: status
    character(len=256) :: msg

    p%coefficient_method = 'conventional-nnls'
    p%coarse_grain = 0
    p%Qs0 = 1.0e10_wp; p%Qp0 = 1.0e10_wp
    p%gamma = 0.0_wp; p%f_transition = 1.0_wp; p%fref = 1.0_wp
    call build_fq8_coefficients(p, tau, strength_s, strength_p, status, msg)
    if (status /= 0) then
       write(*,'(A)') 'elastic: '//trim(msg); error stop 1
    end if

    call require_close('elastic: sum(strength_s)', sum(strength_s), 0.0_wp, 1.0e-10_wp)
    call require_close('elastic: sum(strength_p)', sum(strength_p), 0.0_wp, 1.0e-10_wp)
    call require_close('elastic: modulus scale', &
         fq8_common_modulus_scale(1.0_wp, tau, strength_s), 1.0_wp, 1.0e-12_wp)
    write(*,'(A)') '  elastic half-space: OK'
  end subroutine

  subroutine test_constant_q_halfspace()
    real(wp) :: tau(8), strength_s(8), strength_p(8)
    real(wp) :: freq, q_s, q_target, max_err, scale
    type(fq8_parameters) :: p
    integer :: i, status
    character(len=256) :: msg

    p%coefficient_method = 'conventional-nnls'
    p%coarse_grain = 0
    p%Qs0 = 50.0_wp; p%Qp0 = 50.0_wp
    p%gamma = 0.0_wp; p%f_transition = 1.0_wp; p%fref = 1.0_wp
    call build_fq8_coefficients(p, tau, strength_s, strength_p, status, msg)
    if (status /= 0) then
       write(*,'(A)') 'constant-Q: '//trim(msg); error stop 1
    end if

    call require_close('constant-Q: P/S symmetry', &
         maxval(abs(strength_s - strength_p)), 0.0_wp, 1.0e-14_wp)
    if (any(strength_s < 0.0_wp)) error stop 'constant-Q: negative strength (nnls)'
    if (sum(strength_s) >= 1.0_wp) error stop 'constant-Q: sum(strength) >= 1'

    max_err = 0.0_wp
    do i = 0, 255
       freq = 0.1_wp * (100.0_wp ** (real(i,wp)/255.0_wp))
       q_s = fq8_effective_q(freq, tau, strength_s)
       q_target = 50.0_wp
       max_err = max(max_err, abs(q_s/q_target - 1.0_wp))
    end do
    if (max_err > 0.002_wp) then
       write(*,'(A,ES12.4)') 'constant-Q: max Q error = ', max_err
       error stop 'constant-Q: Q accuracy exceeds 0.2%'
    end if

    scale = fq8_common_modulus_scale(1.0_wp, tau, strength_s)
    if (scale < 1.0_wp .or. scale > 1.10_wp) then
       write(*,'(A,ES15.8)') 'constant-Q: modulus scale = ', scale
       error stop 'constant-Q: modulus scale out of range [1.0, 1.10]'
    end if
    write(*,'(A,ES12.4)') '  constant-Q half-space: OK, max Q err = ', max_err
  end subroutine

  subroutine test_powerlaw_halfspace()
    real(wp) :: tau(8), strength_s(8), strength_p(8)
    real(wp) :: freq, q_s, q_target, max_err, scale
    type(fq8_parameters) :: p
    integer :: i, status
    character(len=256) :: msg

    p%coefficient_method = 'conventional-nnls'
    p%coarse_grain = 0
    p%Qs0 = 50.0_wp; p%Qp0 = 50.0_wp
    p%gamma = 0.6_wp; p%f_transition = 1.0_wp; p%fref = 1.0_wp
    call build_fq8_coefficients(p, tau, strength_s, strength_p, status, msg)
    if (status /= 0) then
       write(*,'(A)') 'powerlaw: '//trim(msg); error stop 1
    end if

    call require_close('powerlaw: P/S symmetry', &
         maxval(abs(strength_s - strength_p)), 0.0_wp, 1.0e-14_wp)
    if (any(strength_s < 0.0_wp)) error stop 'powerlaw: negative strength (nnls)'
    if (sum(strength_s) >= 1.0_wp) error stop 'powerlaw: sum(strength) >= 1'

    max_err = 0.0_wp
    do i = 0, 255
       freq = 0.1_wp * (100.0_wp ** (real(i,wp)/255.0_wp))
       q_target = 50.0_wp
       if (freq >= 1.0_wp) q_target = 50.0_wp * (freq / 1.0_wp) ** 0.6_wp
       q_s = fq8_effective_q(freq, tau, strength_s)
       max_err = max(max_err, abs(q_s/q_target - 1.0_wp))
    end do
    if (max_err > 0.10_wp) then
       write(*,'(A,ES12.4)') 'powerlaw: max Q error = ', max_err
       error stop 'powerlaw: Q accuracy exceeds 10%'
    end if

    scale = fq8_common_modulus_scale(1.0_wp, tau, strength_s)
    if (scale < 1.0_wp .or. scale > 1.10_wp) then
       write(*,'(A,ES15.8)') 'powerlaw: modulus scale = ', scale
       error stop 'powerlaw: modulus scale out of range'
    end if
    write(*,'(A,ES12.4)') '  powerlaw half-space: OK, max Q err = ', max_err
  end subroutine

  subroutine test_layered_model_Q20()
    real(wp) :: tau(8), strength_s(8), strength_p(8)
    real(wp) :: freq, q_s, q_target, max_err
    type(fq8_parameters) :: p
    integer :: i, status
    character(len=256) :: msg

    p%coefficient_method = 'conventional-nnls'
    p%coarse_grain = 0
    p%Qs0 = 20.0_wp; p%Qp0 = 20.0_wp
    p%gamma = 0.6_wp; p%f_transition = 1.0_wp; p%fref = 1.0_wp
    call build_fq8_coefficients(p, tau, strength_s, strength_p, status, msg)
    if (status /= 0) then
       write(*,'(A)') 'Q=20: '//trim(msg); error stop 1
    end if

    if (any(strength_s < 0.0_wp)) error stop 'Q=20: negative strength (nnls)'
    if (sum(strength_s) >= 1.0_wp) error stop 'Q=20: sum(strength) >= 1'

    max_err = 0.0_wp
    do i = 0, 255
       freq = 0.1_wp * (100.0_wp ** (real(i,wp)/255.0_wp))
       q_target = 20.0_wp
       if (freq >= 1.0_wp) q_target = 20.0_wp * (freq / 1.0_wp) ** 0.6_wp
       q_s = fq8_effective_q(freq, tau, strength_s)
       max_err = max(max_err, abs(q_s/q_target - 1.0_wp))
    end do
    if (max_err > 0.10_wp) then
       write(*,'(A,ES12.4)') 'Q=20: max Q error = ', max_err
       error stop 'Q=20: Q accuracy exceeds 10%'
    end if
    write(*,'(A,ES12.4)') '  layered Q=20 (gamma=0.6): OK, max Q err = ', max_err
  end subroutine

  subroutine test_layered_model_Q210()
    real(wp) :: tau(8), strength_s(8), strength_p(8)
    real(wp) :: freq, q_s, q_target, max_err
    type(fq8_parameters) :: p
    integer :: i, status
    character(len=256) :: msg

    p%coefficient_method = 'conventional-nnls'
    p%coarse_grain = 0
    p%Qs0 = 210.0_wp; p%Qp0 = 210.0_wp
    p%gamma = 0.6_wp; p%f_transition = 1.0_wp; p%fref = 1.0_wp
    call build_fq8_coefficients(p, tau, strength_s, strength_p, status, msg)
    if (status /= 0) then
       write(*,'(A)') 'Q=210: '//trim(msg); error stop 1
    end if

    if (any(strength_s < 0.0_wp)) error stop 'Q=210: negative strength (nnls)'

    max_err = 0.0_wp
    do i = 0, 255
       freq = 0.1_wp * (100.0_wp ** (real(i,wp)/255.0_wp))
       q_target = 210.0_wp
       if (freq >= 1.0_wp) q_target = 210.0_wp * (freq / 1.0_wp) ** 0.6_wp
       q_s = fq8_effective_q(freq, tau, strength_s)
       max_err = max(max_err, abs(q_s/q_target - 1.0_wp))
    end do
    if (max_err > 0.10_wp) then
       write(*,'(A,ES12.4)') 'Q=210: max Q error = ', max_err
       error stop 'Q=210: Q accuracy exceeds 10%'
    end if
    write(*,'(A,ES12.4)') '  layered Q=210 (gamma=0.6): OK, max Q err = ', max_err
  end subroutine

  subroutine require_close(label, value, expected, tolerance)
    character(*), intent(in) :: label
    real(wp), intent(in) :: value, expected, tolerance
    if (abs(value - expected) > tolerance) then
       write(*,'(A,2ES24.15,A,ES12.4)') trim(label)//': ', value, expected, &
            ' tolerance=', tolerance
       error stop 1
    end if
  end subroutine require_close

end program fq8_withers_benchmark_test
