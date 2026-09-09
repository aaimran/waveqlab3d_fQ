module viscoelastic_material
  use common, only : wp
  implicit none
  private

  public :: init_viscoelastic_properties
  public :: destroy_viscoelastic_properties

contains

  subroutine init_viscoelastic_properties(M, G, params, block_id)
    use mpi3dcomm, only : allocate_array_body
    use mpi3dbasic, only : error, rank
    use datatypes, only : block_material
    use grid, only : block_grid_t
    use viscoelastic_model, only : viscoelastic_parameters, viscoelastic_coefficients, &
         build_viscoelastic_coefficients, ve_max_relative_error, ve_relaxation_dt_limit
    type(block_material), intent(inout) :: M
    type(block_grid_t), intent(in) :: G
    type(viscoelastic_parameters), intent(in) :: params
    integer, intent(in) :: block_id

    real(wp), parameter :: pi = 3.141592653589793238462643383279_wp
    type(viscoelastic_coefficients) :: coeffs
    integer :: stat, n, i, j, k, l
    character(len=256) :: msg
    real(wp) :: wref, val_s, val_p, vs, vp, mu_s, mu_p, max_s, max_p

    if (M%viscoelastic .or. allocated(M%eta4_ve)) &
         call error('viscoelastic material already initialized', 'init_viscoelastic_properties')

    call build_viscoelastic_coefficients(params, block_id, coeffs, stat, msg)
    if (stat /= 0) call error(trim(msg), 'init_viscoelastic_properties')

    n = coeffs%n
    M%viscoelastic = .true.
    M%n_mechanism_ve = n

    allocate(M%tau_ve(n), M%weight_s_ve(n), M%weight_p_ve(n))
    M%tau_ve = coeffs%tau(1:n)
    M%weight_s_ve = coeffs%w_s(1:n)
    M%weight_p_ve = coeffs%w_p(1:n)

    if (any(M%tau_ve <= 0.0_wp) .or. any(M%weight_s_ve < 0.0_wp) .or. &
        any(M%weight_p_ve < 0.0_wp)) &
         call error('viscoelastic produced invalid coefficients', 'init_viscoelastic_properties')

    call ve_max_relative_error(params%Qs0(block_id), params%gamma, params%f_transition, &
         n, M%tau_ve, M%weight_s_ve, params%fmin, params%fmax, max_s)
    call ve_max_relative_error(params%Qp0(block_id), params%gamma, params%f_transition, &
         n, M%tau_ve, M%weight_p_ve, params%fmin, params%fmax, max_p)
    if (max(max_s, max_p) > params%max_fit_error) then
       write(msg, '(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'viscoelastic fitted response exceeds max_fit_error: S=', max_s, &
            ', P=', max_p, ', limit=', params%max_fit_error
       call error(trim(msg), 'init_viscoelastic_properties')
    end if

    call allocate_array_body(M%Qs_inv_ve, G%C, ghost_nodes=.true.)
    call allocate_array_body(M%Qp_inv_ve, G%C, ghost_nodes=.true.)
    M%Qs_inv_ve = 1.0_wp / params%Qs0(block_id)
    M%Qp_inv_ve = 1.0_wp / params%Qp0(block_id)

    call allocate_array_body(M%eta4_ve, G%C, n, ghost_nodes=.true.); M%eta4_ve = 0.0_wp
    call allocate_array_body(M%Deta4_ve, G%C, n, ghost_nodes=.true.); M%Deta4_ve = 0.0_wp
    call allocate_array_body(M%eta5_ve, G%C, n, ghost_nodes=.true.); M%eta5_ve = 0.0_wp
    call allocate_array_body(M%Deta5_ve, G%C, n, ghost_nodes=.true.); M%Deta5_ve = 0.0_wp
    call allocate_array_body(M%eta6_ve, G%C, n, ghost_nodes=.true.); M%eta6_ve = 0.0_wp
    call allocate_array_body(M%Deta6_ve, G%C, n, ghost_nodes=.true.); M%Deta6_ve = 0.0_wp
    call allocate_array_body(M%eta7_ve, G%C, n, ghost_nodes=.true.); M%eta7_ve = 0.0_wp
    call allocate_array_body(M%Deta7_ve, G%C, n, ghost_nodes=.true.); M%Deta7_ve = 0.0_wp
    call allocate_array_body(M%eta8_ve, G%C, n, ghost_nodes=.true.); M%eta8_ve = 0.0_wp
    call allocate_array_body(M%Deta8_ve, G%C, n, ghost_nodes=.true.); M%Deta8_ve = 0.0_wp
    call allocate_array_body(M%eta9_ve, G%C, n, ghost_nodes=.true.); M%eta9_ve = 0.0_wp
    call allocate_array_body(M%Deta9_ve, G%C, n, ghost_nodes=.true.); M%Deta9_ve = 0.0_wp

    wref = 2.0_wp * pi * params%fref
    do i = G%C%mq, G%C%pq; do j = G%C%mr, G%C%pr; do k = G%C%ms, G%C%ps
       if (M%M(i,j,k,2) <= 0.0_wp .or. M%M(i,j,k,3) <= 0.0_wp) &
            call error('mu and density must be positive for viscoelastic', 'init_viscoelastic_properties')
       val_s = 0.0_wp; val_p = 0.0_wp
       do l = 1, n
          val_s = val_s + M%weight_s_ve(l) / ((wref*wref*M%tau_ve(l)**2 + 1.0_wp) * params%Qs0(block_id))
          val_p = val_p + M%weight_p_ve(l) / ((wref*wref*M%tau_ve(l)**2 + 1.0_wp) * params%Qp0(block_id))
       end do
       if (val_s >= 1.0_wp .or. val_p >= 1.0_wp) &
            call error('invalid viscoelastic modulus correction', 'init_viscoelastic_properties')
       vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
       vp = sqrt((M%M(i,j,k,1) + 2.0_wp * M%M(i,j,k,2)) / M%M(i,j,k,3))
       mu_s = M%M(i,j,k,3) * vs * vs / (1.0_wp - val_s)
       mu_p = M%M(i,j,k,3) * vp * vp / (1.0_wp - val_p)
       M%M(i,j,k,2) = mu_s
       M%M(i,j,k,1) = mu_p - 2.0_wp * mu_s
    end do; end do; end do

    if (rank == 0) then
       write(*, '(A,I0,A,I0)') 'viscoelastic block ', block_id, ': mechanisms=', n
       write(*, '(A,ES12.4,A,ES12.4)') '  Qs0=', params%Qs0(block_id), ', Qp0=', params%Qp0(block_id)
       write(*, '(A,A)') '  weight_method=', trim(params%weight_method)
       write(*, '(A,F8.3,A,F8.3,A)') '  max relative Q error: S=', 100.0_wp*max_s, &
            ' %, P=', 100.0_wp*max_p, ' %'
       write(*, '(A,ES12.4)') '  relaxation dt limit=', ve_relaxation_dt_limit(coeffs)
    end if
  end subroutine init_viscoelastic_properties

  subroutine destroy_viscoelastic_properties(M)
    use datatypes, only : block_material
    type(block_material), intent(inout) :: M
    if (allocated(M%tau_ve)) deallocate(M%tau_ve, M%weight_s_ve, M%weight_p_ve)
    if (allocated(M%Qs_inv_ve)) deallocate(M%Qs_inv_ve, M%Qp_inv_ve)
    if (allocated(M%eta4_ve)) then
       deallocate(M%eta4_ve, M%eta5_ve, M%eta6_ve, M%eta7_ve, M%eta8_ve, M%eta9_ve)
       deallocate(M%Deta4_ve, M%Deta5_ve, M%Deta6_ve, M%Deta7_ve, M%Deta8_ve, M%Deta9_ve)
    end if
    M%viscoelastic = .false.; M%n_mechanism_ve = 0
  end subroutine destroy_viscoelastic_properties

end module viscoelastic_material
