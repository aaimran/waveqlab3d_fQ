program cq_coefficients_test
  use common, only : wp
  use anelastic_cq_model, only : cq_parameters, build_cq_coefficients, cq_max_relative_error
  implicit none
  type(cq_parameters) :: p
  real(wp), allocatable :: tau(:),ws(:),wpv(:),tau2(:),ws2(:),wp2(:)
  real(wp) :: es,ep
  integer :: n
  p%Qs0=[40.0_wp,100.0_wp]; p%Qp0=[80.0_wp,180.0_wp]
  p%coefficient_policy='nnls-block-ps'; p%nnls_samples=256
  do n=4,8
     p%n_mechanisms=n; allocate(tau(n),ws(n),wpv(n))
     call build_cq_coefficients(p,1,tau,ws,wpv)
     if (any(tau <= 0.0_wp) .or. any(ws < 0.0_wp) .or. any(wpv < 0.0_wp)) error stop 1
     if (any(tau(2:n) <= tau(1:n-1))) error stop 2
     if (maxval(abs(ws-wpv)) <= 100.0_wp*epsilon(1.0_wp)) error stop 3
     call cq_max_relative_error(p%Qs0(1),tau,ws,p%fmin,p%fmax,es)
     call cq_max_relative_error(p%Qp0(1),tau,wpv,p%fmin,p%fmax,ep)
     if (es >= huge(1.0_wp) .or. ep >= huge(1.0_wp)) error stop 4
     deallocate(tau,ws,wpv)
  end do
  p%Qs0=[40.0_wp,120.0_wp]; p%Qp0=[69.3_wp,155.9_wp]
  p%n_mechanisms=8; n=8
  allocate(tau(n),ws(n),wpv(n))
  do n=1,2
     call build_cq_coefficients(p,n,tau,ws,wpv)
     call cq_max_relative_error(p%Qs0(n),tau,ws,p%fmin,p%fmax,es)
     call cq_max_relative_error(p%Qp0(n),tau,wpv,p%fmin,p%fmax,ep)
     write(*,'(A,I0,A,F10.6,A,F10.6)') 'production block ',n,': S error=',es,', P error=',ep
     if (max(es,ep) > 0.10_wp) error stop 8
  end do
  deallocate(tau,ws,wpv)
  p%Qs0=[40.0_wp,100.0_wp]; p%Qp0=[80.0_wp,180.0_wp]
  n=8; p%n_mechanisms=n
  allocate(tau(n),ws(n),wpv(n),tau2(n),ws2(n),wp2(n))
  p%coefficient_policy='nnls-shared'
  call build_cq_coefficients(p,1,tau,ws,wpv)
  call build_cq_coefficients(p,2,tau2,ws2,wp2)
  if (maxval(abs(ws-ws2)) > 100.0_wp*epsilon(1.0_wp) .or. &
      maxval(abs(ws-wpv)) > 100.0_wp*epsilon(1.0_wp)) error stop 5
  p%coefficient_policy='nnls-block'
  call build_cq_coefficients(p,1,tau,ws,wpv)
  call build_cq_coefficients(p,2,tau2,ws2,wp2)
  if (maxval(abs(ws-wpv)) > 100.0_wp*epsilon(1.0_wp) .or. &
      maxval(abs(ws-ws2)) <= 100.0_wp*epsilon(1.0_wp)) error stop 6
  p%coefficient_policy='fixed-q50'
  call build_cq_coefficients(p,1,tau,ws,wpv)
  if (maxval(abs(ws-wpv)) > 100.0_wp*epsilon(1.0_wp)) error stop 7
  deallocate(tau,ws,wpv,tau2,ws2,wp2)
  write(*,'(A)') 'anelastic-cQ coefficient tests passed for N=4..8'
end program cq_coefficients_test
