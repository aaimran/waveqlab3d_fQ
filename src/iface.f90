module iface

  !> iface module defines interface type and has routines related to fields on
  !> and enforcement of interface conditions; note that enforcing interface
  !> will also require fields on block boundaries adjacent to interface

  use common, only : wp
  use datatypes, only : iface_type
  implicit none

contains

  subroutine init_iface(I,Bm, II)

!!! initialize interface, based on information stored in neighboring blocks

    use common, only : wp
    use block, only : block_type
    use mpi3d_interface, only : interface3d
    use mpi3dcomm, only : allocate_array_boundary

    implicit none

    type(iface_type),intent(inout) :: I
    type(block_type),intent(in) :: Bm ! blocks on either side of interface
    type(interface3d), intent(in) :: II

!!! check that block sizes are compatible

!  @todo   if (.not.compatible) stop 'block sizes not compatible, cannot create interface'

!!! initialize arrays on interface

    call allocate_array_boundary(I%T ,Bm%G%C,3, I%direction, ghost_nodes=.true.)
    call allocate_array_boundary(I%V ,Bm%G%C,3, I%direction, ghost_nodes=.true.)
    call allocate_array_boundary(I%DV ,Bm%G%C,3, I%direction, ghost_nodes=.true.)
    call allocate_array_boundary(I%S ,Bm%G%C,4, I%direction, ghost_nodes=.true., Fval = 0.0_wp)
    call allocate_array_boundary(I%DS,Bm%G%C,4, I%direction, ghost_nodes=.true., Fval = 1.0e40_wp)
    call allocate_array_boundary(I%W ,Bm%G%C,1, I%direction, ghost_nodes=.true., Fval = 0.0_wp)
    call allocate_array_boundary(I%DW ,Bm%G%C,1, I%direction, ghost_nodes=.true., Fval = 1.0e40_wp)
    call allocate_array_boundary(I%trup ,Bm%G%C,1, I%direction, ghost_nodes=.true., Fval = 1.0e9_wp)
    call allocate_array_boundary(I%Svel ,Bm%G%C,3, I%direction, ghost_nodes=.true., Fval = 0.0_wp)

    !> Set the interface ID
!     I%id = Bm%id + Bp%id
    I%II = II

  end subroutine init_iface


  subroutine scale_rates_iface(I,A)

!!! multiply rates by RK coefficient A

    implicit none

    type(iface_type),intent(inout) :: I
    real(kind = wp),intent(in) :: A

    I%DS = A*I%DS
    I%DW = A*I%DW

  end subroutine scale_rates_iface


  subroutine update_fields_iface(I,dt)

!!! update fields using rates

    implicit none

    type(iface_type),intent(inout) :: I
    real(kind = wp),intent(in) :: dt

    I%S = I%S + dt*I%DS
    I%W = I%W + dt*I%DW

  end subroutine update_fields_iface

end module iface
