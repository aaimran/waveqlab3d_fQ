module mpi3d_interface

  use common, only : wp
  use diagnostics, only : fatal_local

  implicit none

  public :: is_on_interface2d, get_interface_coord2d, &
       get_dimensions2d, interface_communicator2d, interface2d

  ! n             : number of processes on the interface
  ! coord         : global, cartesian coordinates of process
  ! icoord        : interface coordinate
  ! normal        : outward pointing unit normal, with respect to this
  ! side of the interface
  ! comm          : interface communicator (all processes on the interface)
  ! rank          : rank of process
  ! rank_neighbor : rank of process on opposite side of interface
  ! on_interface  : true, if this process is on the interface
  type interface2d
     integer :: n,coord(2),icoord,normal(2),comm,rank,rank_neighbor
     logical :: on_interface
  end type interface2d

  type interface3d
     integer :: n(2),coord(3),icoord(2),normal(3),comm,rank,rank_neighbor
     logical :: on_interface
     integer :: face_qr,face_qs,face_rs
  end type interface3d

  interface new_interface
     module procedure new_interface2d,new_interface3d
  end interface new_interface

contains

  ! Construct a new interface between two blocks
  ! processes on the interface are determined by the normal.
  ! The normal must be outward pointing with respect to a block side
  ! This means that n_1 = -n_2, where n_1 is the normal on one side
  ! of the interface, and n_2 is the normal on the opposite side.
  ! All data needed for future communication are contained in I
  subroutine new_interface2d(coord,cart_size,normal,comm,I)

    use mpi

    implicit none

    integer,intent(in) :: coord(2),cart_size(2),normal(2),comm
    type(interface2d),intent(out) :: I

    integer :: ierr

    call interface_communicator2d(coord,cart_size,normal,comm,&
         I%comm,I%on_interface)

    if(.not. I%on_interface) return

    call MPI_Comm_rank(I%comm,I%rank,ierr)

    call find_neighbor2d(coord,cart_size,normal,I%comm,&
         I%rank_neighbor)


    I%n = get_dimensions2d(cart_size,normal)
    I%coord = coord
    I%icoord = get_interface_coord2d(coord,normal)
    I%normal = normal

  end subroutine new_interface2d

  subroutine new_interface3d(coord,cart_size,normal,comm,C,I)

    use mpi
    use mpi3dcomm, only : cartesian3d_t

    implicit none

    integer,intent(in) :: coord(3),cart_size(3),normal(3),comm
    type(cartesian3d_t),intent(in) :: C
    type(interface3d),intent(out) :: I

    integer :: ierr, comm_size, expected_size

    call interface_communicator3d(coord,cart_size,normal,comm,&
         I%comm,I%on_interface)

    if(.not. I%on_interface) return

    call MPI_Comm_rank(I%comm,I%rank,ierr)

    ! A conforming two-sided interface has exactly two ranks for every
    ! tangential Cartesian coordinate.  Checking this before neighbor lookup
    ! prevents a malformed communicator from producing an arbitrary pairing.
    call MPI_Comm_size(I%comm,comm_size,ierr)
    I%n = get_dimensions3d(cart_size,normal)
    expected_size = 2*I%n(1)*I%n(2)
    if (comm_size /= expected_size .and. comm_size /= 1) then
       call fatal_local('RUN-IFACE-001', &
            'Interface communicator size does not match the two-sided tangential topology.', &
            'new_interface3d')
    end if

    call find_neighbor3d(coord,cart_size,normal,I%comm,&
         I%rank_neighbor)

    I%coord = coord
    I%icoord = get_interface_coord3d(coord,normal)
    I%normal = normal

    call define_faces(C,I)

  end subroutine new_interface3d


  ! Check that a process is on the interface.
  function is_on_interface2d(coord,cart_size,normal) &
       result(is_interface_node)

    use mpi3dcomm, only : cartesian3d_t
    use mpi, only : MPI_PROC_NULL

    implicit none

    integer :: coord(2),cart_size(2),normal(2)
    logical :: is_interface_node

    ! Check all four interface sides defined by the outward unit
    ! normal with respect to each side
    if (normal(1) == 1 .and. normal(2) == 0 ) then
       is_interface_node = (coord(1) == cart_size(1) - 1)
    elseif (normal(1) == 0 .and. normal(2) == 1 ) then
       is_interface_node = (coord(2) == cart_size(2) - 1)
    elseif (normal(1) == -1 .and. normal(2) == 0 ) then
       is_interface_node = (coord(1) == 0)
    elseif (normal(1) ==  0 .and. normal(2) == -1 ) then
       is_interface_node = (coord(2) == 0)
    endif

  end function is_on_interface2d

  ! Check that a process is on the interface.
  function is_on_interface3d(coord,cart_size,normal) &
       result(is_interface_node)

    use mpi3dcomm, only : cartesian3d_t
    use mpi, only : MPI_PROC_NULL

    implicit none

    integer :: coord(3),cart_size(3),normal(3)
    logical :: is_interface_node

    ! Check all four interface sides defined by the outward unit
    ! normal with respect to each side
    if (normal(1) == 1 .and. normal(2) == 0 .and. normal(3) == 0) then
       is_interface_node = (coord(1) == cart_size(1) - 1)
    elseif (normal(1) == 0 .and. normal(2) == 1 .and. normal(3) == 0) then
       is_interface_node = (coord(2) == cart_size(2) - 1)
    elseif (normal(1) == -1 .and. normal(2) == 0 .and. normal(3) == 0 ) then
       is_interface_node = (coord(1) == 0)
    elseif (normal(1) ==  0 .and. normal(2) == -1 .and. normal(3) == 0 ) then
       is_interface_node = (coord(2) == 0)
    elseif (normal(1) ==  0 .and. normal(2) == 0 .and. normal(3) ==  1 ) then
       is_interface_node = (coord(3) == cart_size(3) - 1)
    elseif (normal(1) ==  0 .and. normal(2) == 0 .and. normal(3) == -1 ) then
       is_interface_node = (coord(3) == 0)
    endif

  end function is_on_interface3d

  ! Determine the number of processes on the interface
  function get_dimensions2d(cart_size,normal) &
       result(interface_size)

    implicit none

    integer,intent(in) :: cart_size(2),normal(2)
    integer :: interface_size

    if (abs(normal(1)) == 1 .and. normal(2) == 0 ) then
       interface_size = cart_size(2)
    elseif (normal(1) == 0 .and. abs(normal(2)) == 1 ) then
       interface_size = cart_size(1)
    endif

  end function get_dimensions2d

  function get_dimensions3d(cart_size,normal) &
       result(interface_size)

    implicit none

    integer,intent(in) :: cart_size(3),normal(3)
    integer :: interface_size(2)

    ! Plane convention
    ! x : (y,z)
    ! y : (x,z)
    ! z : (x,y)

    if (abs(normal(1)) == 1 .and. normal(2) == 0  .and. normal(3) == 0) then
       interface_size = (/ cart_size(2), cart_size(3) /)
    elseif (normal(1) == 0 .and. abs(normal(2)) == 1 .and. normal(3) == 0) then
       interface_size = (/ cart_size(1), cart_size(3) /)
    elseif (normal(1) == 0 .and. normal(2) == 0 .and. abs(normal(3)) ==  1 ) then
       interface_size = (/ cart_size(1), cart_size(2) /)
    endif

  end function get_dimensions3d
  ! Convert between the global block coordinate system
  ! to a local, interface coordinate system
  function get_interface_coord2d(coord,normal) &
       result(interface_coord)

    implicit none

    integer,intent(in) :: coord(2), normal(2)
    integer :: interface_coord, tangent(2)

    ! Tangent vector with positive components
    tangent = (/ abs(normal(2)), abs(normal(1)) /)

    interface_coord = tangent(1) * coord(1) + tangent(2) * coord(2)

  end function get_interface_coord2d


  function get_interface_coord3d(coord,normal) &
       result(interface_coord)

    implicit none

    integer,intent(in) :: coord(3), normal(3)
    integer :: interface_coord(2)

    if (abs(normal(1)) == 1 .and. normal(2) == 0  .and. normal(3) == 0) then
       interface_coord = (/ coord(2), coord(3) /)
    elseif (normal(1) == 0 .and. abs(normal(2)) == 1 .and. normal(3) == 0) then
       interface_coord = (/ coord(1), coord(3) /)
    elseif (normal(1) == 0 .and. normal(2) == 0 .and. abs(normal(3)) ==  1 ) then
       interface_coord = (/ coord(1), coord(2) /)
    endif

  end function get_interface_coord3d

  ! Construct an interface communicator,
  ! only processes on the interface with opposing normals
  ! and in 'comm' are included in 'comm_interface'
  !
  ! Example: Block1:  (0, 1) (1, 1)
  !                   (0, 0) (1, 0)
  !          Block2:  (0, 1) (1, 1)
  !                   (0, 0) (1, 0)
  ! The interface between blocks are defined with normal (0,-1), for Block1
  ! and normal (0,1), for Block2.
  subroutine interface_communicator2d(coord,cart_size,normal,comm, &
       comm_interface,on_interface)

    use mpi
    use mpi3dbasic, only : new_communicator

    implicit none

    integer,intent(in) :: coord(2),cart_size(2),normal(2),comm
    integer,intent(out) :: comm_interface
    logical,intent(out) :: on_interface

    on_interface = .false.

    ! Determine processes belonging to the interface
    if (is_on_interface2d(coord,cart_size,normal)) on_interface = .true.

    call new_communicator(on_interface,comm_interface,comm)

  end subroutine interface_communicator2d

  subroutine interface_communicator3d(coord,cart_size,normal,comm, &
       comm_interface,on_interface)

    use mpi
    use mpi3dbasic, only : new_communicator

    implicit none

    integer,intent(in) :: coord(3),cart_size(3),normal(3),comm
    integer,intent(out) :: comm_interface
    logical,intent(out) :: on_interface

    on_interface = .false.

    ! Determine processes belonging to the interface
    if (is_on_interface3d(coord,cart_size,normal)) on_interface = .true.

    call new_communicator(on_interface,comm_interface,comm)

  end subroutine interface_communicator3d

  ! Find neighbors on opposite side of interface.
  ! This subroutine must only be called by processes in interface_comm
  subroutine find_neighbor2d(coord,cart_size,normal,comm,rank_neighbor)

    use mpi

    implicit none

    integer,intent(in) :: coord(2), cart_size(2), normal(2),comm
    integer :: n,rank,interface_coord,ierr,i
    integer,dimension(:),allocatable :: ranks,x
    integer,intent(out) :: rank_neighbor

    ! Ensure that only processes on the interface participate
    if (.not. is_on_interface2d(coord,cart_size,normal)) return

    n = get_dimensions2d(cart_size,normal)

    call MPI_Comm_rank(comm,rank,ierr)
    interface_coord = get_interface_coord2d(coord,normal)

    ! There are 2*n processes on the interface (n on each side)
    allocate(ranks(2*n),x(2*n))

    call MPI_Allgather(rank,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(interface_coord,1,MPI_INTEGER,x,1,MPI_INTEGER,comm, ierr)

    ! Search for process with the same interface coordinates
    ! as self
    do i=1,2*n
       if (x(i) == interface_coord .and. ranks(i) /= rank) exit
    end do
    ! @todo should check that x(i) == interface_coord, if not true,
    ! then the block sides are not adjacent
    if (i > 2*n) i = 2*n

    rank_neighbor = ranks(i)

    deallocate(ranks,x)

  end subroutine find_neighbor2d

  subroutine find_neighbor3d(coord,cart_size,normal,comm,rank_neighbor)

    use mpi
    use mpi3dbasic, only : nprocs
    implicit none

    integer,intent(in) :: coord(3), cart_size(3), normal(3),comm
    integer :: n(2),rank,interface_coord(2),ierr,i,nx,ny,comm_size
    integer,dimension(:),allocatable :: ranks,x,y
    integer,intent(out) :: rank_neighbor

    ! Ensure that only processes on the interface participate
    if (.not. is_on_interface3d(coord,cart_size,normal)) return

    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,comm_size,ierr)

    if (nprocs == 1) then
      rank_neighbor = rank
      return
    end if

    n = get_dimensions3d(cart_size,normal)
    nx = n(1)
    ny = n(2)
    interface_coord = get_interface_coord3d(coord,normal)

    ! Allocate based on actual communicator size to handle both 1-block and 2-block cases
    allocate(ranks(comm_size),x(comm_size),y(comm_size))

    call MPI_Allgather(rank,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm,ierr)
    call MPI_Allgather(interface_coord(1),1,MPI_INTEGER,x,1,MPI_INTEGER,comm, ierr)
    call MPI_Allgather(interface_coord(2),1,MPI_INTEGER,y,1,MPI_INTEGER,comm, ierr)

    ! Search for process with the same interface coordinates as self
    ! (on the opposite side of the interface, i.e., different rank)
    rank_neighbor = -1  ! Default value: not found
    do i=1,comm_size
       if (x(i) == interface_coord(1) .and. &
            y(i) == interface_coord(2) .and. ranks(i) /= rank) then
          rank_neighbor = ranks(i)
          exit
       end if
    end do

    if (rank_neighbor == -1) then
       call fatal_local('RUN-IFACE-002', &
            'No opposite-side rank has matching tangential Cartesian coordinates.', &
            'find_neighbor3d')
    end if

    deallocate(ranks,x,y)

  end subroutine find_neighbor3d

  subroutine exchange_interface_neighbors2d(C,F,I)

    use mpi
    use mpi3dcomm, only : cartesian3d_t

    implicit none

    type(cartesian3d_t),intent(in) :: C
    real(kind = wp),dimension(C%mbq:C%pbq,C%mbr:C%pbr),intent(inout) :: F
    type(interface2d),intent(in) :: I
    integer,parameter :: tag1 = 1,tag2 = 2,tag3 = 3,tag4 = 4
    integer :: ierr

    if(.not. I%on_interface) return

    if (I%normal(1) == 1 .and. I%normal(2) == 0) then
       call MPI_SendRecv( &
            F(C%pq,C%mr),1,C%line_r,I%rank_neighbor,tag1, &
            F(C%pq+1,C%mr),1,C%line_r,I%rank_neighbor,tag2, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    elseif (I%normal(1) == -1 .and. I%normal(2) == 0) then
       call MPI_SendRecv( &
            F(C%mq,C%mr),1,C%line_r,I%rank_neighbor,tag2, &
            F(C%mq-1,C%mr),1,C%line_r,I%rank_neighbor,tag1, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    elseif (I%normal(1) == 0 .and. I%normal(2) == 1) then
       call MPI_SendRecv( &
            F(C%mq,C%pr),1,C%line_q,I%rank_neighbor,tag3, &
            F(C%mq,C%pr+1),1,C%line_q,I%rank_neighbor,tag4, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    elseif (I%normal(1) == 0 .and. I%normal(2) == -1) then
       call MPI_SendRecv( &
            F(C%mq,C%mr),1,C%line_q,I%rank_neighbor,tag4, &
            F(C%mq,C%mr-1),1,C%line_q,I%rank_neighbor,tag3, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

  end subroutine exchange_interface_neighbors2d

  subroutine exchange_interface_neighbors3d(C,F,I)

    use mpi
    use mpi3dcomm, only : cartesian3d_t

    implicit none

    type(cartesian3d_t),intent(in) :: C
    real(kind = wp),dimension(C%mbq:C%pbq,C%mbr:C%pbr,C%mbs:C%pbs),intent(inout) :: F
    type(interface3d),intent(in) :: I
    integer,parameter :: tag1 = 1,tag2 = 2,tag3 = 3,tag4 = 4, &
         tag5 = 5,tag6 = 6
    integer :: ierr,nx,ny,nz

    if(.not. I%on_interface) return

    nx = I%normal(1)
    ny = I%normal(2)
    nz = I%normal(3)


    if (nx == 1 .and. ny ==  0 .and. nz == 0) then
       call MPI_SendRecv( &
            F(C%pq,C%mr,C%ms),1,I%face_rs,I%rank_neighbor,tag1, &
            F(C%pq+1,C%mr,C%ms),1,I%face_rs,I%rank_neighbor,tag2, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

    if (nx == -1 .and. ny ==  0 .and. nz == 0) then
       call MPI_SendRecv( &
            F(C%mq,C%mr,C%ms),1,I%face_rs,I%rank_neighbor,tag2, &
            F(C%mq-1,C%mr,C%ms),1,I%face_rs,I%rank_neighbor,tag1, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

    if (nx == 0 .and. ny ==  1 .and. nz == 0) then
       call MPI_SendRecv( &
            F(C%mq,C%pr,C%ms),1,I%face_qs,I%rank_neighbor,tag3, &
            F(C%mq,C%pr+1,C%ms),1,I%face_qs,I%rank_neighbor,tag4, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

    if (nx == 0 .and. ny == -1 .and. nz == 0) then
       call MPI_SendRecv( &
            F(C%mq,C%mr,C%ms),1,I%face_qs,I%rank_neighbor,tag4, &
            F(C%mq,C%mr-1,C%ms),1,I%face_qs,I%rank_neighbor,tag3, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

    if (nx == 0 .and. ny ==  0 .and. nz == 1) then
       call MPI_SendRecv( &
            F(C%mq,C%mr,C%ps),1,I%face_qr,I%rank_neighbor,tag5, &
            F(C%mq,C%mr,C%ps+1),1,I%face_qr,I%rank_neighbor,tag6, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

    if (nx == 0 .and. ny ==  0 .and. nz == -1) then
       call MPI_SendRecv( &
            F(C%mq,C%mr,C%ms),1,I%face_qr,I%rank_neighbor,tag6, &
            F(C%mq,C%mr,C%ms-1),1,I%face_qr,I%rank_neighbor,tag5, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

  end subroutine exchange_interface_neighbors3d

  subroutine exchange_interface_neighbors(F, Fopp, C, I)

    use mpi
    use mpi3dbasic, only : MPI_REAL_PW, rank
    use mpi3dcomm, only : cartesian3d_t

    implicit none

    type(cartesian3d_t),intent(in) :: C
    real(kind = wp),dimension(C%mbr:C%pbr,C%mbs:C%pbs),intent(in) :: F
    real(kind = wp),dimension(C%mbr:C%pbr,C%mbs:C%pbs),intent(out) :: Fopp
    type(interface3d),intent(in) :: I
    integer,parameter :: tag1 = 1,tag2 = 2,tag3 = 3,tag4 = 4, &
         tag5 = 5,tag6 = 6
    integer :: ierr,nx,ny,nz, npts

    if(.not. I%on_interface) return

    nx = I%normal(1)
    ny = I%normal(2)
    nz = I%normal(3)

    npts = (C%pbr-C%mbr + 1)*(C%pbs - C%mbs + 1)

    if (nx == 1 .and. ny ==  0 .and. nz == 0) then
       call MPI_SendRecv( &
            F,npts,MPI_REAL_PW,I%rank_neighbor,tag1, &
            Fopp,npts,MPI_REAL_PW,I%rank_neighbor,tag2, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

    if (nx == -1 .and. ny ==  0 .and. nz == 0) then
       call MPI_SendRecv( &
            F,npts,MPI_REAL_PW,I%rank_neighbor,tag2, &
            Fopp,npts,MPI_REAL_PW,I%rank_neighbor,tag1, &
            I%comm,MPI_STATUS_IGNORE,ierr)
    endif

  end subroutine exchange_interface_neighbors

  ! Consider moving this to Cartesian
  subroutine define_faces(C,I)

    use mpi3dbasic, only : MPI_REAL_PW
    use mpi3dcomm, only : cartesian3d_t
    use mpi

    implicit none

    type(cartesian3d_t),intent(in) :: C
    type(interface3d),intent(inout) :: I

    integer :: count_qr,length_qr,stride_qr, &
         count_rs,length_rs,stride_rs, &
         count_qs,length_qs,stride_qs

    integer :: ierr,pw

    ! Define the face_qr, given by the plane with normal: (0 0 a)
    count_qr = C%lnr
    stride_qr =  (C%lnq + 2*C%nb)
    length_qr = C%lnq
    call MPI_Type_vector(&
         count_qr,length_qr,stride_qr,MPI_REAL_PW,I%face_qr,ierr)
    call MPI_Type_commit(I%face_qr,ierr)

    ! Define the face_rs given by the plane with normal: (a 0 0)
    call MPI_Type_size(MPI_REAL_PW,pw,ierr)
    count_rs = C%lns
    stride_rs = (C%lnq + 2*C%nb)*(C%lnr + 2*C%nb)*pw
    length_rs = 1

    ! Note the use of hvector, which specifies the stride in bytes
    call MPI_Type_hvector(&
         count_rs,length_rs,stride_rs,C%line_r,I%face_rs,ierr)
    call MPI_Type_commit(I%face_rs,ierr)

    ! Define the face_qs given by the plane with normal: (0 a 0)
    count_qs = C%lns
    stride_qs = (C%lnq + 2*C%nb)*(C%lnr + 2*C%nb)*pw
    length_qs = 1

    call MPI_Type_hvector(&
         count_qs,length_qs,stride_qs,C%line_q,I%face_qs,ierr)
    call MPI_Type_commit(I%face_qs,ierr)
  end subroutine define_faces

end module mpi3d_interface
