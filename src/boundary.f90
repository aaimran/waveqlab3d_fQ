module boundary

  !> boundary module handles fields and properties on block boundaries
  !> and has routines used to enforce boundary conditions

  use common, only : wp
  use datatypes, only : block_boundary

  implicit none

contains


  subroutine init_boundary(B, G, physics)

    !> @brief initialize block_boundary type

    use mpi3dcomm, only: cartesian3d_t, allocate_array_boundary
    use grid, only : block_grid_t
    use unit_normals, only : local_orth_vectors_q,local_orth_vectors_r,local_orth_vectors_s

    implicit none

    type(block_boundary),intent(out) :: B(6)
    type(block_grid_t),intent(in) :: G
    character(*),intent(in) :: physics

    integer :: side,i
    character(1) :: direction
    type(cartesian3d_t) :: C

    ! loop over all sides of the block

    C = G%C

    do side = 1,6

       ! keep track of directions in computational coordinates q,r,s instead of x,y,z

       select case(side)
       case(1,2)
          ! side 1 = negative q
          ! side 2 = positive q
          direction = 'q'
       case(3,4)
          ! side 3 = negative r
          ! side 4 = positive r
          direction = 'r'
       case(5,6)
          ! side 5 = negative s
          ! side 6 = positive s
          direction = 's'
       end select

       ! set unit normal vectors on boundary
       ! (formulas below are incorrect, but illustrate use of metric derivative arrays)

       call allocate_array_boundary(B(side)%n_l, C, 3, direction, ghost_nodes=.true.)
       call allocate_array_boundary(B(side)%n_m, C, 3, direction, ghost_nodes=.true.)
       call allocate_array_boundary(B(side)%n_n, C, 3, direction, ghost_nodes=.true.)

       select case(side)
       
       case(1)
          i = C%mq
          call local_orth_vectors_q(i, G, B(side))
       case(2)
          i = C%pq
          call local_orth_vectors_q(i, G, B(side))
          
       case(3)
          i = C%mr
          call local_orth_vectors_r(i, G, B(side))
       case(4)
          i = C%pr
          call local_orth_vectors_r(i, G, B(side))
       case(5)
          i = C%ms
          call local_orth_vectors_s(i, G, B(side))
       case(6)
          i = C%ps
          call local_orth_vectors_s(i, G, B(side))
          
       end select
       
       ! initialize arrays that are physics-specific
       
       select case(physics)
          
       case default

          stop 'invalid block physics in init_boundary'

       case('elastic')

          ! allocate arrays and set initial values to zero
          call allocate_array_boundary(B(side)%X, C, 2, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%M, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%Mopp, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%Fopp, C, 9, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%F, C, 9, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%DF, C, 9, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%U, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%DU, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)

       case('acoustic')

          ! allocate arrays and set initial values to zero
          call allocate_array_boundary(B(side)%X, C, 2, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%M, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%Mopp, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%Fopp, C, 9, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%F, C, 9, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%DF, C, 9, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%U, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)
          call allocate_array_boundary(B(side)%DU, C, 3, direction, ghost_nodes = .true., Fval = 0.0_wp)

       end select

    end do

  end subroutine init_boundary

  subroutine exchange_fields_across_interface(B, C, I)

    use mpi3dcomm, only : cartesian3d_t
    use mpi3d_interface, only : interface3d, exchange_interface_neighbors

    type(block_boundary), intent(inout) :: B
    type(cartesian3d_t), intent(in) :: C
    type(interface3d), intent(in) :: I
    integer :: nf

    do nf = 1, 9
      call exchange_interface_neighbors(B%F(:,:,nf), B%Fopp(:,:,nf), C, I)
    end do

  end subroutine exchange_fields_across_interface

  subroutine exchange_materials_across_interface(B, C, I)

    use mpi3dcomm, only : cartesian3d_t
    use mpi3d_interface, only : interface3d, exchange_interface_neighbors

    type(block_boundary), intent(inout) :: B
    type(cartesian3d_t), intent(in) :: C
    type(interface3d), intent(in) :: I
    integer :: nf

    do nf = 1, 3
      call exchange_interface_neighbors(B%M(:,:,nf), B%Mopp(:,:,nf), C, I)
    end do

  end subroutine exchange_materials_across_interface

end module boundary
