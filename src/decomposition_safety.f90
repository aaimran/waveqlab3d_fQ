module decomposition_safety

  use common, only : wp
  implicit none
  private

  integer, parameter, public :: CONSERVATIVE_MINIMUM_OWNED = 20

  type, public :: stencil_requirements_t
     integer :: halo_width = 0
     integer :: boundary_width = 0
     integer :: minimum_global_points = 0
     integer :: minimum_owned_points = 0
     logical :: supported = .false.
     character(len=64) :: operator_name = ''
  end type stencil_requirements_t

  public :: get_stencil_requirements, topology_fits, select_single_block_topology

contains

  pure function get_stencil_requirements(fd_type, order) result(req)
    character(*), intent(in) :: fd_type
    integer, intent(in) :: order
    type(stencil_requirements_t) :: req

    select case (trim(adjustl(fd_type)))
    case ('traditional')
       req%halo_width = 3
       req%boundary_width = 6
       req%operator_name = 'traditional-6'
       req%supported = .true.
    case ('upwind')
       select case (order)
       case (2,3)
          req%halo_width = 3; req%boundary_width = 2
       case (4,5)
          req%halo_width = 3; req%boundary_width = 4
       case (6,7)
          req%halo_width = 4; req%boundary_width = 6
       case (8,9)
          req%halo_width = 5; req%boundary_width = 8
       case default
          return
       end select
       write(req%operator_name,'(A,I0)') 'upwind-', order
       req%supported = .true.
    case ('upwind_drp')
       select case (order)
       case (3)
          req%halo_width = 3; req%boundary_width = 4
       case (4,5,66)
          req%halo_width = 4; req%boundary_width = 6
       case (6,7,679)
          req%halo_width = 5; req%boundary_width = 8
       case default
          return
       end select
       write(req%operator_name,'(A,I0)') 'upwind-drp-', order
       req%supported = .true.
    case default
       return
    end select

    req%minimum_global_points = 2*req%boundary_width
    req%minimum_owned_points = max(CONSERVATIVE_MINIMUM_OWNED, &
         2*req%halo_width+1, 2*req%boundary_width)
  end function get_stencil_requirements


  pure logical function topology_fits(points, dims, required_owned) result(fits)
    integer, intent(in) :: points(3), dims(3), required_owned
    integer :: d

    fits = .true.
    do d = 1, 3
       if (dims(d) < 1 .or. dims(d) > points(d)) then
          fits = .false.
          return
       end if
       if (dims(d) > 1 .and. points(d)/dims(d) < required_owned) then
          fits = .false.
          return
       end if
    end do
  end function topology_fits


  subroutine select_single_block_topology(points, nranks, required_owned, dims, found)
    integer, intent(in) :: points(3), nranks, required_owned
    integer, intent(out) :: dims(3)
    logical, intent(out) :: found
    integer :: p, q, r
    real(wp) :: score, best_score, local(3)

    dims = 0
    found = .false.
    best_score = huge(1.0_wp)
    do p = 1, nranks
       if (mod(nranks,p) /= 0) cycle
       do q = 1, nranks/p
          if (mod(nranks/p,q) /= 0) cycle
          r = nranks/(p*q)
          if (.not.topology_fits(points, [p,q,r], required_owned)) cycle
          local = real(points,wp)/real([p,q,r],wp)
          ! Approximate communication area per rank, normalized by volume.
          score = 1.0_wp/local(1) + 1.0_wp/local(2) + 1.0_wp/local(3) + &
               100.0_wp*epsilon(1.0_wp)*real((q-1)+2*(r-1),wp)
          if (score < best_score) then
             best_score = score
             dims = [p,q,r]
             found = .true.
          end if
       end do
    end do
  end subroutine select_single_block_topology

end module decomposition_safety
