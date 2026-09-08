program pre_wql3d

    use mpi
    use common, only : wp
    use datatypes, only : block_grid_t, block_temp_parameters
    use grid, only : cartesian_grid_3d, Curve_Grid_3D3

    implicit none

    type(block_grid_t) :: G, Gt
    real(kind = wp) :: starts(3), ends(3), lc, rc
    real(kind = wp) :: CFL, t_final, topo
    character(64) :: coupling, mesh_source, type_of_mesh, profile_type, profile_path
    character(64) :: topography_type,topography_path
    logical :: use_mms, use_topography
    character(len=256) :: name, problem, response, ifname
    integer :: nblocks
    logical :: w_fault
    integer :: w_stride, ny, nz
    integer :: i
    integer :: infile, stat, len_ifname
    integer :: rank, nprocs, ierr
    type(block_temp_parameters) :: btp(2)

    namelist /problem_list/ name, problem, nblocks, t_final, CFL, coupling, topo, mesh_source, type_of_mesh, w_stride, w_fault
    namelist /block_list/ btp

    ! Read the command-line to get the input filename
    call get_command_argument(1, ifname, length=len_ifname, status=stat)
    if (stat /= 0) stop ':: Problem reading input filename from commandline.'

    ! Read the input file
    open(newunit=infile, file=ifname, iostat=stat, status='old')
    if (stat/=0) stop ':: Cannot open input file'

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

    read(infile, nml = problem_list, iostat = stat)
    read(infile, nml = block_list, iostat = stat)

    use_mms = .false.

    do i = 1, 2
        starts = btp(i)%aqrs
        ends = btp(i)%bqrs
        lc = btp(i)%lc
        rc = btp(i)%rc
        profile_type = btp(i)%profile_type
        profile_path = btp(i)%profile_path
        G%C%nq = btp(i)%nqrs(1)
        G%C%nr = btp(i)%nqrs(2)
        G%C%ns = btp(i)%nqrs(3)
        G%C%mq = 1
        G%C%pq = btp(i)%nqrs(1)
        G%C%mr = 1
        G%C%pr = btp(i)%nqrs(2)
        G%C%ms = 1
        G%C%ps = btp(i)%nqrs(3)

        allocate(G%X(G%C%nq, G%C%nr, G%C%ns, 3))

        if (type_of_mesh == "cartesian") then
            call cartesian_grid_3d(G,starts, ends)
        else if (type_of_mesh == "curvilinear") then
            call Curve_Grid_3D3(G, starts(1), ends(1), starts(2), ends(2), &
                 starts(3), ends(3), lc, rc, use_mms, profile_type, profile_path, &
                 use_topography, topography_type,topography_path, topo, ny, nz)
        end if

        call write_mesh_serial(i, G)
        deallocate(G%X)
    end do

    call MPI_Finalize(ierr)

contains
    subroutine write_mesh_serial(i, G)
        integer, intent(in) :: i
        type(block_grid_t), intent(in) :: G
        character(len = 1) :: is, ext(3)
        integer :: k, fh(3)
        character(len=256) :: fname

        ext = ['X', 'Y', 'Z']

        ! Save the block id as string
        write(is, '(i1)') i

        do k = 1, 3
            fname = "block_grid_"//is//"."//ext(k)
            call MPI_File_Open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                               MPI_INFO_NULL, fh(k), ierr)
            call MPI_File_write(fh(k), real(G%X(:,:,:,k), 8), size(G%X(:,:,:,k)), &
                                MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(fh(k), ierr)
        end do

    end subroutine write_mesh_serial

    subroutine read_mesh_serial(i, G)
        integer, intent(in) :: i
        type(block_grid_t), intent(inout) :: G
        character(len = 1) :: is, ext(3)
        integer :: k, fh(3)
        character(len=256) :: fname

        ext = ['X', 'Y', 'Z']

        ! Save the block id as string
        write(is, '(i1)') i

        do k = 1, 3
            fname = "block_grid_"//is//"."//ext(k)
            call MPI_File_Open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, &
                               MPI_INFO_NULL, fh(k), ierr)
            call MPI_File_read(fh(k), real(G%X(:,:,:,k), 4), size(G%X(:,:,:,k)), &
                                MPI_REAL, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(fh(k), ierr)
        end do

    end subroutine read_mesh_serial
end program pre_wql3d
