module mpi3dcomm

  use common, only : wp

  implicit none

  ! Cartesian communicator type, including
  ! COMM = communicator
  ! RANK = rank of process
  ! SIZE = size of communicator
  ! COORD = integer coordinates of process
  ! RANK_M = rank of neighbor in minus direction
  ! RANK_P = rank of neighbor in plus  direction

  type :: cartesian3d_t
     integer :: nb, &
          nq,mq,pq,mbq,pbq,lnq, &
          nr,mr,pr,mbr,pbr,lnr, &
          ns,ms,ps,mbs,pbs,lns
     integer :: comm,rank,size_q,size_r,size_s,coord(3), & ! maybe turn size_q, etc. into vector: size(3)
          rank_mq,rank_pq,rank_mr,rank_pr,rank_ms,rank_ps ! maybe turn rank_mq, etc. into vectors: rank_m(3), rank_p(3)
     integer :: line_q,line_r,line_s, &
                block_q,block_r,block_s, &
                block3d_qr,block3d_qs,block3d_rs ! MPI types for communication, again maybe organized in vectors?
     integer :: array_w,array_s ! MPI types for I/O (leave here, Eric will use in I/O module)
  end type cartesian3d_t

  ! naming convention on indices, example for x direction:
  ! nb = number of additional boundary points (ghost nodes)
  ! nq = total number of points (on all processes, excluding ghost nodes)
  ! mq = starting (minus) index on process
  ! pq = ending   (plus ) index on process
  ! mbq = starting index on process, including ghost nodes
  ! pbq = ending   index on process, including ghost nodes
  ! lnq = total number of points (on process)

  ! module procedures for allocating arrays
  interface allocate_array_body
     module procedure allocate_array_body_2d,allocate_array_body_2d_nf,allocate_array_body_3d,allocate_array_body_3d_nF
  end interface allocate_array_body

  interface allocate_array_boundary
     module procedure allocate_array_boundary_3d_nF
  end interface allocate_array_boundary

contains


  ! 3D decomposition subroutine

  subroutine decompose3d(C,nq,nr,ns,nb,periodic_q,periodic_r,periodic_s,comm, &
       method,size_q_in,size_r_in,size_s_in)

    use mpi3dbasic, only : is_master,decompose1d,error,MPI_REAL_PW,pw
    use mpi

    implicit none

    type(cartesian3d_t),intent(out) :: C
    integer,intent(in) :: nq,nr,ns,nb,comm,size_q_in,size_r_in,size_s_in
    logical,intent(in) :: periodic_q,periodic_r,periodic_s
    character(*),intent(inout) :: method

    integer,parameter :: dim=3
    integer :: ierr,size,index,shift,size_qrs(3),l,count_qr,stride_qr,length_qr
    logical :: periodic(3),reorder
    integer,dimension(:),allocatable :: blocklengths,types
    integer,dimension(:),allocatable :: displacements_mpi1 ! MPI-1
    !integer(MPI_ADDRESS_KIND),dimension(:),allocatable :: displacements ! MPI-2

    ! domain size

    C%nq = nq
    C%nr = nr
    C%ns = ns
    C%nb = nb
!    C%nb = 4
    ! properties of original communicator

    call MPI_Comm_size(comm,size,ierr)

    ! processor layout in Cartesian topology

    if (method=='manual') then

       ! manual decomposition, user specifies number of processes in x and y directions

       if (size/=size_q_in*size_r_in*size_s_in.and.is_master(comm)) &
            call error('Error: Incorrect number of processors for manual decomposition','decompose3d')
       size_qrs = (/ size_q_in,size_r_in,size_s_in /)

    else

       select case(method)
       case default ! 3D: MPI_Dims_create selects number of processes in all directions
          size_qrs = 0
       case('1Dx') ! 1D, x: MPI_Dims_create selects number of processes in x direction only
          size_qrs(1) = 0
          size_qrs(2) = 1
          size_qrs(3) = 1
       case('1Dy') ! 1D, y: MPI_Dims_create selects number of processes in y direction only
          size_qrs(1) = 1
          size_qrs(2) = 0
          size_qrs(3) = 1
       case('1Dz') ! 1D, z: MPI_Dims_create selects number of processes in z direction only
          size_qrs(1) = 1
          size_qrs(2) = 1
          size_qrs(3) = 0
       case('2Dxy') ! 2D, xy: MPI_Dims_create selects number of processes in x and y directions only
          size_qrs(1) = 0
          size_qrs(2) = 0
          size_qrs(3) = 1
       case('2Dxz') ! 2D, xz: MPI_Dims_create selects number of processes in x and z directions only
          size_qrs(1) = 0
          size_qrs(2) = 1
          size_qrs(3) = 0
       case('2Dyz') ! 2D, yz: MPI_Dims_create selects number of processes in y and z directions only
          size_qrs(1) = 1
          size_qrs(2) = 0
          size_qrs(3) = 0
       end select

       ! special cases
       if (C%nq==1) size_qrs(1) = 1
       if (C%nr==1) size_qrs(2) = 1
       if (C%ns==1) size_qrs(3) = 1

       call MPI_Dims_create(size,dim,size_qrs,ierr)

    end if

    C%size_q = size_qrs(1) ! number of processes in x direction
    C%size_r = size_qrs(2) ! number of processes in y direction
    C%size_s = size_qrs(3) ! number of processes in z direction

    ! 3D Cartesian communicator and coordinates

    periodic = (/ periodic_q,periodic_r,periodic_s /)
    reorder = .true.
    call MPI_Cart_create(comm,dim,size_qrs,periodic,reorder,C%comm,ierr)
    call MPI_Comm_rank(C%comm,C%rank,ierr)
    call MPI_Cart_coords(C%comm,C%rank,dim,C%coord,ierr)

    ! nearest neighbors in x-direction

    index = 0; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mq,C%rank_pq,ierr)

    ! nearest neighbors in y-direction

    index = 1; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mr,C%rank_pr,ierr)

    ! nearest neighbors in z-direction

    index = 2; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_ms,C%rank_ps,ierr)


    ! initial data distribution on processors

    call decompose1d(C%nq,C%size_q,C%coord(1),C%mq,C%pq,C%lnq)
    call decompose1d(C%nr,C%size_r,C%coord(2),C%mr,C%pr,C%lnr)
    call decompose1d(C%ns,C%size_s,C%coord(3),C%ms,C%ps,C%lns)
    C%mbq = C%mq-C%nb
    C%pbq = C%pq+C%nb
    C%mbr = C%mr-C%nb
    C%pbr = C%pr+C%nb
    C%mbs = C%ms-C%nb
    C%pbs = C%ps+C%nb

    ! MPI types containing lines of constant x, y, and z

    ! variable x, constant y, constant z
    call MPI_Type_vector(C%lnq,1,1,MPI_REAL_PW,C%line_q,ierr)
    call MPI_Type_commit(C%line_q,ierr)

    ! variable y, constant x, constant z
    call MPI_Type_vector(C%lnr,1,C%lnq+2*C%nb,MPI_REAL_PW,C%line_r,ierr)
    call MPI_Type_commit(C%line_r,ierr)

    ! variable z, constant x, constant y
    call MPI_Type_vector(C%lns,1,(C%lnq+2*C%nb)*(C%lnr+2*C%nb),MPI_REAL_PW,C%line_s,ierr)
    call MPI_Type_commit(C%line_s,ierr)


    ! MPI types containing 2D boundary blocks

    allocate(blocklengths(C%nb),types(C%nb))
    blocklengths = 1
    allocate(displacements_mpi1(C%nb)) ! MPI-1
    !allocate(displacements(C%nb)) ! MPI-2

    !types = C%line_r
    !! MPI-1
    !displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnq+2*C%nb)*(C%lnr+2*C%nb)
    !call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_s,ierr)
    !call MPI_Type_commit(C%block_s,ierr)
    !! MPI-2
    !!displacements = (/ (l,l=0,C%nb-1) /)*(C%lnq+2*C%nb)
    !!call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_r,ierr)

    count_qr = C%lnr
    stride_qr =  (C%lnq + 2*C%nb)
    length_qr = C%lnq
    call MPI_Type_vector(&
         count_qr,length_qr,stride_qr,MPI_REAL_PW,C%block_s,ierr)
    call MPI_Type_commit(C%block_s,ierr)

    types = C%line_r
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_q,ierr)
    call MPI_Type_commit(C%block_q,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_q,ierr)

    types = C%line_q
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnq+2*C%nb)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_r,ierr)
    call MPI_Type_commit(C%block_r,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)*(C%lnq+2*C%nb)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_r,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2


    ! MPI Types containing 3D boundary blocks


    ! 3D MPI block covering faces in x-y plane

    allocate(blocklengths(C%lnq),types(C%lnq))
    blocklengths = 1
    allocate(displacements_mpi1(C%lnq)) ! MPI-1
    !allocate(displacements(C%nb)) ! MPI-2

    types = C%block_s
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%lnq-1) /)*(C%lnq+2*C%nb)*(C%lnr+2*C%nb)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,&
         C%block3d_qr,ierr)
    call MPI_Type_commit(C%block3d_qr,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%lnq-1) /)
    !call MPI_Type_create_struct(C%lnq,blocklengths,displacements_mpi1*pw,types,&
    !   C%block3d_qr,ierr)

    !types = C%block_s
    !! MPI-1
    !displacements_mpi1 = (/ (l,l=0,C%lnq-1) /)
    !call MPI_Type_struct(C%lnq,blocklengths,displacements_mpi1*pw,types,&
    !    C%block3d_qr,ierr)
    !call MPI_Type_commit(C%block3d_qr,ierr)
    !! MPI-2
    !!displacements = (/ (l,l=0,C%lnq-1) /)
    !!call MPI_Type_create_struct(C%lnq,blocklengths,displacements_mpi1*pw,types,&
    !!   C%block3d_qr,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2


    ! 3D MPI block covering faces in x-z plane

    allocate(blocklengths(C%lns),types(C%lns))
    blocklengths = 1
    allocate(displacements_mpi1(C%lns)) ! MPI-1
    !allocate(displacements(C%lns)) ! MPI-2

    types = C%block_r
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%lns-1) /)*(C%lnq+2*C%nb)*(C%lnr+2*C%nb)
    call MPI_Type_struct(C%lns,blocklengths,displacements_mpi1*pw,types,&
         C%block3d_qs,ierr)
    call MPI_Type_commit(C%block3d_qs,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%lns-1) /)*(C%lnq+2*C%nb)*(C%lnr+2*C%nb)
    !call MPI_Type_create_struct(C%lns,blocklengths,displacements*pw,types,&
    !   C%block3d_qs,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2


    ! 3D MPI block covering faces in x-y plane

    allocate(blocklengths(C%lns),types(C%lns))
    blocklengths = 1
    allocate(displacements_mpi1(C%lns)) ! MPI-1
    !allocate(displacements(C%lns)) ! MPI-2

    types = C%block_q
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%lns-1) /)*(C%lnq+2*C%nb)*(C%lnr+2*C%nb)
    call MPI_Type_struct(C%lns,blocklengths,displacements_mpi1*pw,types,&
         C%block3d_rs,ierr)
    call MPI_Type_commit(C%block3d_rs,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%lns-1) /)*(C%lnq+2*C%nb)*(C%lnr+2*C%nb)
    !call MPI_Type_create_struct(C%lns,blocklengths,displacements*pw,types,&
    !   C%block3d_rs,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2


    ! TO DO: ADD NEW DATA TYPE THAT PACKS ALL NF FIELDS TOGETHER FOR MORE EFFICIENT COMMUNICATION

  end subroutine decompose3d

  !2D Decomposition subroutine

  subroutine decompose2d(C,nq,nr,nb,periodic_q,periodic_r,comm, &
       method,size_q_in,size_r_in)

    use mpi3dbasic, only : is_master,decompose1d,error,MPI_REAL_PW,pw
    use mpi

    implicit none

    type(cartesian3d_t),intent(out) :: C
    integer,intent(in) :: nq,nr,nb,comm,size_q_in,size_r_in
    logical,intent(in) :: periodic_q,periodic_r
    character(*),intent(inout) :: method

    integer,parameter :: dim=2
    integer :: ierr,size,index,shift,size_qr(2),l
    logical :: periodic(2),reorder
    integer,dimension(:),allocatable :: blocklengths,types
    integer,dimension(:),allocatable :: displacements_mpi1 ! MPI-1
    !integer(MPI_ADDRESS_KIND),dimension(:),allocatable :: displacements ! MPI-2

    ! domain size

    C%nq = nq
    C%nr = nr
    C%nb = nb

    ! properties of original communicator

    call MPI_Comm_size(comm,size,ierr)

    ! processor layout in Cartesian topology

    if (method=='manual') then

       ! manual decomposition, user specifies number of processes in x and y directions

       if (size/=size_q_in*size_r_in.and.is_master(comm)) &
            call error('Error: Incorrect number of processors for manual decomposition','decompose2d')
       size_qr = (/ size_q_in,size_r_in /)

    else

       select case(method)
       case default ! 2D: MPI_Dims_create selects number of processes in both directions
          size_qr = 0
       case('1Dx') ! 1D, x: MPI_Dims_create selects number of processes in x direction only
          size_qr(1) = 0
          size_qr(2) = 1
       case('1Dy') ! 1D, y: MPI_Dims_create selects number of processes in y direction only
          size_qr(1) = 1
          size_qr(2) = 0
       end select

       ! special cases
       if (C%nq==1) size_qr(1) = 1
       if (C%nr==1) size_qr(2) = 1

       call MPI_Dims_create(size,dim,size_qr,ierr)

    end if

    C%size_q = size_qr(1) ! number of processes in x direction
    C%size_r = size_qr(2) ! number of processes in y direction

    ! 2D Cartesian communicator and coordinates

    periodic = (/ periodic_q,periodic_r /)
    reorder = .true.
    call MPI_Cart_create(comm,dim,size_qr,periodic,reorder,C%comm,ierr)
    call MPI_Comm_rank(C%comm,C%rank,ierr)
    call MPI_Cart_coords(C%comm,C%rank,dim,C%coord,ierr)

    ! nearest neighbors in x-direction

    index = 0; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mq,C%rank_pq,ierr)

    ! nearest neighbors in y-direction

    index = 1; shift = 1
    call MPI_Cart_shift(C%comm,index,shift,C%rank_mr,C%rank_pr,ierr)

    ! initial data distribution on processors

    call decompose1d(C%nq,C%size_q,C%coord(1),C%mq,C%pq,C%lnq)
    call decompose1d(C%nr,C%size_r,C%coord(2),C%mr,C%pr,C%lnr)
    C%mbq = C%mq-C%nb
    C%pbq = C%pq+C%nb
    C%mbr = C%mr-C%nb
    C%pbr = C%pr+C%nb

    ! MPI types containing lines of constant x and y

    ! variable x, constant y
    call MPI_Type_vector(C%lnq,1,1,MPI_REAL_PW,C%line_q,ierr)
    call MPI_Type_commit(C%line_q,ierr)

    ! variable y, constant x
    call MPI_Type_vector(C%lnr,1,C%lnq+2*C%nb,MPI_REAL_PW,C%line_r,ierr)
    call MPI_Type_commit(C%line_r,ierr)

    ! MPI types containing boundary blocks

    allocate(blocklengths(C%nb),types(C%nb))
    blocklengths = 1
    allocate(displacements_mpi1(C%nb)) ! MPI-1
    !allocate(displacements(C%nb)) ! MPI-2


    types = C%line_r
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_q,ierr)
    call MPI_Type_commit(C%block_q,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_q,ierr)


    types = C%line_q
    ! MPI-1
    displacements_mpi1 = (/ (l,l=0,C%nb-1) /)*(C%lnq+2*C%nb)
    call MPI_Type_struct(C%nb,blocklengths,displacements_mpi1*pw,types,C%block_r,ierr)
    call MPI_Type_commit(C%block_r,ierr)
    ! MPI-2
    !displacements = (/ (l,l=0,C%nb-1) /)*(C%lnq+2*C%nb)
    !call MPI_Type_create_struct(C%nb,blocklengths,displacements*pw,types,C%block_r,ierr)

    deallocate(blocklengths,types)
    deallocate(displacements_mpi1) ! MPI-1
    !deallocate(displacements) ! MPI-2

    ! TO DO: ADD NEW DATA TYPE THAT PACKS ALL NF FIELDS TOGETHER FOR MORE EFFICIENT COMMUNICATION

  end subroutine decompose2d


  subroutine exchange_all_neighbors(C,F)

    use mpi3dbasic, only: rank

    implicit none

    type(cartesian3d_t),intent(in) :: C
    real(kind = wp),dimension(C%mbq:C%pbq,C%mbr:C%pbr,C%mbs:C%pbs),intent(inout) :: F

    call exchange_neighbors(C,F,'xm')
    call exchange_neighbors(C,F,'xp')
    call exchange_neighbors(C,F,'ym')
    call exchange_neighbors(C,F,'yp')
    call exchange_neighbors(C,F,'zm')
    call exchange_neighbors(C,F,'zp')

  end subroutine exchange_all_neighbors


  subroutine exchange_neighbors(C,F,side)

    use mpi

    implicit none

    type(cartesian3d_t),intent(in) :: C
    real(kind = wp),dimension(C%mbq:C%pbq,C%mbr:C%pbr,C%mbs:C%pbs),intent(inout) :: F
    character(*),intent(in) :: side

    integer :: ierr
    integer,parameter :: tagxm=1, tagxp=2, tagym=3, tagyp=4, tagzm=5, tagzp=6

    ! share data between processes using blocking sendrecv with 2D communicator


    select case(side)
    case('xm') ! pq --> mbq
       call MPI_SendRecv( &
            F(C%pq-C%nb+1,C%mr,C%ms),1,C%block3d_rs,C%rank_pq,tagxm, &
            F(C%mq-C%nb  ,C%mr,C%ms),1,C%block3d_rs,C%rank_mq,tagxm, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('xp') ! pbq <-- mq
       call MPI_SendRecv( &
            F(C%mq  ,C%mr,C%ms),1,C%block3d_rs,C%rank_mq,tagxp, &
            F(C%pq+1,C%mr,C%ms),1,C%block3d_rs,C%rank_pq,tagxp, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('ym') ! pr --> mbr
       call MPI_SendRecv( &
            F(C%mq,C%pr-C%nb+1,C%ms),1,C%block3d_qs,C%rank_pr,tagym, &
            F(C%mq,C%mr-C%nb  ,C%ms),1,C%block3d_qs,C%rank_mr,tagym, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('yp') ! pbr <-- mr
       call MPI_SendRecv( &
            F(C%mq,C%mr  ,C%ms),1,C%block3d_qs,C%rank_mr,tagyp, &
            F(C%mq,C%pr+1,C%ms),1,C%block3d_qs,C%rank_pr,tagyp, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('zm') ! ps --> mbs
       call MPI_SendRecv( &
            F(C%mq,C%mr,C%ps-C%nb+1),1,C%block3d_qr,C%rank_ps,tagzm, &
            F(C%mq,C%mr,C%ms-C%nb  ),1,C%block3d_qr,C%rank_ms,tagzm, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    case('zp') ! ps --> mbs
       call MPI_SendRecv( &
            F(C%mq,C%mr,C%ms  ),1,C%block3d_qr,C%rank_ms,tagzp, &
            F(C%mq,C%mr,C%ps+1),1,C%block3d_qr,C%rank_ps,tagzp, &
            C%comm,MPI_STATUS_IGNORE,ierr)
    end select

  end subroutine exchange_neighbors


  subroutine test_exchange(C,F)

    use common, only : wp
    implicit none

    type(cartesian3d_t),intent(in) :: C
    real(kind = wp),dimension(C%mbq:C%pbq,C%mbr:C%pbr,C%mbs:C%pbs),intent(inout) :: F

    integer :: i,j,k
    real(kind = wp) :: G

    ! initialize array


    do j = C%mr,C%pr
       do i = C%mq,C%pq
          do k = C%ms,C%ps
             F(i,j,k) = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(k, wp)
          end do
       end do
    end do

    ! exchange

    call exchange_all_neighbors(C,F)

    ! check

    ! left
    if (C%mq/=1) then
       do j = C%mr,C%pr
          do i = C%mbq,C%mq-1
             do k = C%ms,C%ps
                G = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(K, wp)
                if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if

    ! right
    if (C%pq/=C%nq) then
       do j = C%mr,C%pr
          do i = C%pq+1,C%pbq
             do k = C%ms,C%ps
                G = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(k, wp)
                if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if

    ! bottom
    if (C%mr/=1) then
       do j = C%mbr,C%mr-1
          do i = C%mq,C%pq
             do k = C%ms,C%ps
                G = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(k, wp)
                if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if

    ! top
    if (C%pr/=C%nr) then
       do j = C%pr+1,C%pbr
          do i = C%mq,C%pq
             do k = C%ms,C%ps
                G = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(k, wp)
                if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if

    ! up (z-direction)
    if (C%mr/=1) then
       do j = C%mr,C%pr
          do i = C%mq,C%pq
             do k = C%mbs,C%ms-1
                G = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(k, wp)
                if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if

    ! down (z-direction)
    if (C%pr/=C%nr) then
       do j = C%mr,C%pr
          do i = C%mq,C%pq
             do k = C%ps+1,C%pbs
                G = real(i, wp)+1000.0_wp*real(j, wp)+1000000.0_wp*real(k, wp)
                if (F(i,j,k)/=G) print *, i,j,k,G,F(i,j,k)
             end do
          end do
       end do
    end if
  end subroutine test_exchange

  subroutine allocate_array_body_2d(F,C,ghost_nodes,Fval)

    use common, only : wp
    implicit none

    real(kind = wp),dimension(:,:),allocatable,intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    logical,intent(in) :: ghost_nodes
    real(kind = wp),intent(in),optional :: Fval

!     if (associated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbq:C%pbq,C%mbr:C%pbr))
    else
       allocate(F(C%mq :C%pq ,C%mr :C%pr))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1.0e20_wp
    end if

  end subroutine allocate_array_body_2d


  subroutine allocate_array_body_2d_nF(F,C,nF,ghost_nodes,Fval)

    use common, only : wp
    implicit none

    real(kind = wp),dimension(:,:,:),allocatable,intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    integer,intent(in) :: nF
    logical,intent(in) :: ghost_nodes
    real(kind = wp),intent(in),optional :: Fval

!     if (associated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbq:C%pbq,C%mbr:C%pbr,nF))
    else
       allocate(F(C%mq :C%pq ,C%mr :C%pr ,nF))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1.0e20_wp
    end if

  end subroutine allocate_array_body_2d_nF



  subroutine allocate_array_body_3d(F,C,ghost_nodes,Fval)

    use common, only : wp
    implicit none

    real(kind = wp),dimension(:,:,:),allocatable,intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    logical,intent(in) :: ghost_nodes
    real(kind = wp),intent(in),optional :: Fval

!     if (associated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbq:C%pbq,C%mbr:C%pbr,C%mbs:C%pbs))
    else
       allocate(F(C%mq :C%pq ,C%mr :C%pr ,C%ms :C%ps ))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1.0e20_wp
    end if

  end subroutine allocate_array_body_3d


  subroutine allocate_array_body_3d_nF(F,C,nF,ghost_nodes,Fval)

    use common, only : wp
    implicit none

    real(kind = wp),dimension(:,:,:,:),allocatable,intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    integer,intent(in) :: nF
    logical,intent(in) :: ghost_nodes
    real(kind = wp),intent(in),optional :: Fval

!     if (associated(F)) return

    if (ghost_nodes) then
       allocate(F(C%mbq:C%pbq,C%mbr:C%pbr,C%mbs:C%pbs,nF))
    else
       allocate(F(C%mq :C%pq ,C%mr :C%pr ,C%ms :C%ps ,nF))
    end if
    if (present(Fval)) then
       F = Fval
    else
       F = 1.0e20_wp
    end if

  end subroutine allocate_array_body_3d_nF

  subroutine allocate_array_boundary_3d_nF(F, C, nF, direction, ghost_nodes, Fval)
    use common, only : wp
    implicit none

    real(kind = wp),dimension(:,:,:),allocatable,intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    integer,intent(in) :: nF
    logical,intent(in) :: ghost_nodes
    real(kind = wp),intent(in),optional :: Fval
    character(1), intent(in) :: direction


!     if (associated(F)) return
    select case(direction)
    case default
       stop 'invalid direction in allocate_array_boundary'
    case('q')
       if (ghost_nodes) then
          allocate(F(C%mbr:C%pbr,C%mbs:C%pbs,nF))
       else
          allocate(F(C%mr :C%pr ,C%ms :C%ps ,nF))
       end if
    case('r')
       if (ghost_nodes) then
          allocate(F(C%mbq:C%pbq, C%mbs:C%pbs,nF))
       else
          allocate(F(C%mq :C%pq, C%ms :C%ps ,nF))
       end if
    case('s')
       if (ghost_nodes) then
          allocate(F(C%mbq:C%pbq,C%mbr:C%pbr, nF))
       else
          allocate(F(C%mq :C%pq ,C%mr :C%pr,nF))
       end if
    end select

    if (present(Fval)) then
       F = Fval
    else
       F = 1.0e20_wp
    end if
  end subroutine allocate_array_boundary_3d_nF


end module mpi3dcomm
