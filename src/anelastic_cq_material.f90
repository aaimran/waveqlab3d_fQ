module anelastic_cq_material

  use common, only : wp
  use datatypes, only : block_material, block_grid_t
  use anelastic_cq_model, only : cq_parameters, build_cq_coefficients, cq_max_relative_error
  implicit none
  private
  public :: init_anelastic_cq_properties, destroy_anelastic_cq_properties
  public :: apply_anelastic_cq_strain

contains

  subroutine apply_anelastic_cq_strain(M,x,y,z,Dx,Dy,Dz,rate)
    type(block_material), intent(inout) :: M
    integer, intent(in) :: x,y,z
    real(wp), intent(in) :: Dx(:),Dy(:),Dz(:)
    real(wp), intent(inout) :: rate(:)
    integer :: i,n
    real(wp) :: tr,mu2,pm,sm,bulk
    n=M%n_mechanism_cQ
    rate(4)=rate(4)-sum(M%eta4cQ(x,y,z,1:n)); rate(5)=rate(5)-sum(M%eta5cQ(x,y,z,1:n))
    rate(6)=rate(6)-sum(M%eta6cQ(x,y,z,1:n)); rate(7)=rate(7)-sum(M%eta7cQ(x,y,z,1:n))
    rate(8)=rate(8)-sum(M%eta8cQ(x,y,z,1:n)); rate(9)=rate(9)-sum(M%eta9cQ(x,y,z,1:n))
    tr=Dx(1)+Dy(2)+Dz(3); mu2=2.0_wp*M%M(x,y,z,2)
    do i=1,n
       sm=M%weight_s_cQ(i)*M%Qs_inv_cQ(x,y,z)
       pm=M%weight_p_cQ(i)*M%Qp_inv_cQ(x,y,z)
       bulk=(M%M(x,y,z,1)+mu2)*pm-mu2*sm
       M%Deta4cQ(x,y,z,i)=M%Deta4cQ(x,y,z,i)+(mu2*sm*Dx(1)+bulk*tr-M%eta4cQ(x,y,z,i))/M%tau_cQ(i)
       M%Deta5cQ(x,y,z,i)=M%Deta5cQ(x,y,z,i)+(mu2*sm*Dy(2)+bulk*tr-M%eta5cQ(x,y,z,i))/M%tau_cQ(i)
       M%Deta6cQ(x,y,z,i)=M%Deta6cQ(x,y,z,i)+(mu2*sm*Dz(3)+bulk*tr-M%eta6cQ(x,y,z,i))/M%tau_cQ(i)
       M%Deta7cQ(x,y,z,i)=M%Deta7cQ(x,y,z,i)+(M%M(x,y,z,2)*sm*(Dy(1)+Dx(2))-M%eta7cQ(x,y,z,i))/M%tau_cQ(i)
       M%Deta8cQ(x,y,z,i)=M%Deta8cQ(x,y,z,i)+(M%M(x,y,z,2)*sm*(Dz(1)+Dx(3))-M%eta8cQ(x,y,z,i))/M%tau_cQ(i)
       M%Deta9cQ(x,y,z,i)=M%Deta9cQ(x,y,z,i)+(M%M(x,y,z,2)*sm*(Dz(2)+Dy(3))-M%eta9cQ(x,y,z,i))/M%tau_cQ(i)
    end do
  end subroutine apply_anelastic_cq_strain

  subroutine init_anelastic_cq_properties(M,G,parameters,block_id)
    use mpi3dcomm, only : allocate_array_body
    use mpi3dbasic, only : error, rank
    type(block_material), intent(inout) :: M
    type(block_grid_t), intent(in) :: G
    type(cq_parameters), intent(in) :: parameters
    integer, intent(in) :: block_id
    real(wp), parameter :: pi=3.141592653589793_wp
    real(wp) :: wref,val_s,val_p,vs,vp,mu_s,mu_p,max_s,max_p
    integer :: i,j,k,l,n
    character(len=256) :: fit_message

    if (M%anelastic_cQ .or. allocated(M%eta4cQ)) &
         call error('anelastic-cQ material is already initialized','init_anelastic_cq_properties')
    n=parameters%n_mechanisms
    M%anelastic_cQ=.true.; M%n_mechanism_cQ=n
    M%fref_cQ=parameters%fref; M%fmin_cQ=parameters%fmin; M%fmax_cQ=parameters%fmax
    M%Qs0_cQ=parameters%Qs0(block_id); M%Qp0_cQ=parameters%Qp0(block_id)
    M%coefficient_policy_cQ=parameters%coefficient_policy
    allocate(M%tau_cQ(n),M%weight_s_cQ(n),M%weight_p_cQ(n))
    call build_cq_coefficients(parameters,block_id,M%tau_cQ,M%weight_s_cQ,M%weight_p_cQ)
    if (any(M%tau_cQ <= 0.0_wp) .or. any(M%weight_s_cQ < 0.0_wp) .or. &
        any(M%weight_p_cQ < 0.0_wp)) &
         call error('anelastic-cQ produced invalid coefficients','init_anelastic_cq_properties')

    call cq_max_relative_error(M%Qs0_cQ,M%tau_cQ,M%weight_s_cQ, &
         parameters%fmin,parameters%fmax,max_s)
    call cq_max_relative_error(M%Qp0_cQ,M%tau_cQ,M%weight_p_cQ, &
         parameters%fmin,parameters%fmax,max_p)
    if (max(max_s,max_p) > parameters%max_fit_error) then
       write(fit_message,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'anelastic-cQ fitted response exceeds max_fit_error: S=',max_s, &
            ', P=',max_p,', limit=',parameters%max_fit_error
       call error(trim(fit_message),'init_anelastic_cq_properties')
    end if

    call allocate_array_body(M%Qp_inv_cQ,G%C,ghost_nodes=.true.)
    call allocate_array_body(M%Qs_inv_cQ,G%C,ghost_nodes=.true.)
    M%Qp_inv_cQ=1.0_wp/M%Qp0_cQ; M%Qs_inv_cQ=1.0_wp/M%Qs0_cQ
    call allocate_array_body(M%eta4cQ,G%C,n,ghost_nodes=.true.); M%eta4cQ=0.0_wp
    call allocate_array_body(M%Deta4cQ,G%C,n,ghost_nodes=.true.); M%Deta4cQ=0.0_wp
    call allocate_array_body(M%eta5cQ,G%C,n,ghost_nodes=.true.); M%eta5cQ=0.0_wp
    call allocate_array_body(M%Deta5cQ,G%C,n,ghost_nodes=.true.); M%Deta5cQ=0.0_wp
    call allocate_array_body(M%eta6cQ,G%C,n,ghost_nodes=.true.); M%eta6cQ=0.0_wp
    call allocate_array_body(M%Deta6cQ,G%C,n,ghost_nodes=.true.); M%Deta6cQ=0.0_wp
    call allocate_array_body(M%eta7cQ,G%C,n,ghost_nodes=.true.); M%eta7cQ=0.0_wp
    call allocate_array_body(M%Deta7cQ,G%C,n,ghost_nodes=.true.); M%Deta7cQ=0.0_wp
    call allocate_array_body(M%eta8cQ,G%C,n,ghost_nodes=.true.); M%eta8cQ=0.0_wp
    call allocate_array_body(M%Deta8cQ,G%C,n,ghost_nodes=.true.); M%Deta8cQ=0.0_wp
    call allocate_array_body(M%eta9cQ,G%C,n,ghost_nodes=.true.); M%eta9cQ=0.0_wp
    call allocate_array_body(M%Deta9cQ,G%C,n,ghost_nodes=.true.); M%Deta9cQ=0.0_wp

    wref=2.0_wp*pi*parameters%fref
    do i=G%C%mq,G%C%pq; do j=G%C%mr,G%C%pr; do k=G%C%ms,G%C%ps
       if (M%M(i,j,k,2) <= 0.0_wp .or. M%M(i,j,k,3) <= 0.0_wp) &
            call error('mu and density must be positive for anelastic-cQ','init_anelastic_cq_properties')
       val_s=0.0_wp; val_p=0.0_wp
       do l=1,n
          val_s=val_s+M%weight_s_cQ(l)/((wref*wref*M%tau_cQ(l)**2+1.0_wp)*M%Qs0_cQ)
          val_p=val_p+M%weight_p_cQ(l)/((wref*wref*M%tau_cQ(l)**2+1.0_wp)*M%Qp0_cQ)
       end do
       if (val_s >= 1.0_wp .or. val_p >= 1.0_wp) &
            call error('invalid anelastic-cQ modulus correction','init_anelastic_cq_properties')
       vs=sqrt(M%M(i,j,k,2)/M%M(i,j,k,3))
       vp=sqrt((M%M(i,j,k,1)+2.0_wp*M%M(i,j,k,2))/M%M(i,j,k,3))
       mu_s=M%M(i,j,k,3)*vs*vs/(1.0_wp-val_s)
       mu_p=M%M(i,j,k,3)*vp*vp/(1.0_wp-val_p)
       M%M(i,j,k,2)=mu_s
       M%M(i,j,k,1)=mu_p-2.0_wp*mu_s
    end do; end do; end do

    if (rank == 0) then
       write(*,'(A,I0,A,I0)') 'anelastic-cQ block ',block_id,': mechanisms=',n
       write(*,'(A,ES12.4,A,ES12.4)') '  Qs0=',M%Qs0_cQ,', Qp0=',M%Qp0_cQ
       write(*,'(A,A)') '  coefficient policy=',trim(parameters%coefficient_policy)
       write(*,'(A,F8.3,A,F8.3,A)') '  max relative Q error: S=',100.0_wp*max_s, &
            ' %, P=',100.0_wp*max_p,' %'
    end if
  end subroutine init_anelastic_cq_properties


  subroutine destroy_anelastic_cq_properties(M)
    type(block_material), intent(inout) :: M
    if (allocated(M%tau_cQ)) deallocate(M%tau_cQ,M%weight_s_cQ,M%weight_p_cQ)
    if (allocated(M%Qp_inv_cQ)) deallocate(M%Qp_inv_cQ,M%Qs_inv_cQ)
    if (allocated(M%eta4cQ)) then
       deallocate(M%eta4cQ,M%eta5cQ,M%eta6cQ,M%eta7cQ,M%eta8cQ,M%eta9cQ)
       deallocate(M%Deta4cQ,M%Deta5cQ,M%Deta6cQ,M%Deta7cQ,M%Deta8cQ,M%Deta9cQ)
    end if
    M%anelastic_cQ=.false.; M%n_mechanism_cQ=0
  end subroutine destroy_anelastic_cq_properties

end module anelastic_cq_material
