module plane_output

  use mpi
  use common, only : wp
  use mpi3dbasic, only : MPI_REAL_PW, error
  use datatypes, only : plane_output_type, plane_output_plane, block_grid_t

  implicit none

  integer, parameter :: max_planes = 64
  integer, parameter :: tag_meta = 31001
  integer, parameter :: tag_vx   = 31002
  integer, parameter :: tag_vy   = 31003
  integer, parameter :: tag_vz   = 31004

contains

  subroutine init_plane_output(input, domain_name, planes, G, block_num)

    implicit none

    integer, intent(in) :: input
    character(*), intent(in) :: domain_name
    type(plane_output_type), intent(inout) :: planes
    type(block_grid_t), intent(in) :: G
    integer, intent(in) :: block_num

    logical :: enable_plane_output
    integer :: nplanes
    integer :: plane_stride
    character(len=2) :: plane_axis(max_planes)
    real(kind=wp) :: plane_coord(max_planes)
    character(len=64) :: plane_name(max_planes)
    character(len=256) :: plane_file_directory

    namelist /plane_output_list/ enable_plane_output, nplanes, plane_axis, plane_coord, plane_name, plane_stride, plane_file_directory

    integer :: stat, p

    enable_plane_output = .false.
    nplanes = 0
    plane_stride = 1
    plane_axis = ''
    plane_coord = 0.0_wp
    plane_name = ''
    plane_file_directory = '.'

    rewind(input)
    read(input, nml=plane_output_list, iostat=stat)

    if (stat /= 0) then
      planes%enabled = .false.
      planes%nplanes = 0
      if (allocated(planes%P)) deallocate(planes%P)
      return
    end if

    if (.not. enable_plane_output .or. nplanes < 1) then
      planes%enabled = .false.
      planes%nplanes = 0
      if (allocated(planes%P)) deallocate(planes%P)
      return
    end if

    if (nplanes > max_planes) call error('nplanes exceeds max_planes in &plane_output_list', 'init_plane_output')

    planes%enabled = .true.
    planes%nplanes = nplanes
    if (allocated(planes%P)) deallocate(planes%P)
    allocate(planes%P(nplanes))

    do p = 1, nplanes
      planes%P(p)%enabled = .true.
      planes%P(p)%axis = adjustl(plane_axis(p))
      planes%P(p)%coord = plane_coord(p)
      planes%P(p)%name = adjustl(plane_name(p))
      planes%P(p)%file_directory = adjustl(plane_file_directory)
      planes%P(p)%stride = max(1, plane_stride)
      planes%P(p)%step_counter = 0

      call init_one_plane(domain_name, planes%P(p), G, block_num)
    end do

  end subroutine init_plane_output


  subroutine init_one_plane(domain_name, P, G, block_num)

    implicit none

    character(*), intent(in) :: domain_name
    type(plane_output_plane), intent(inout) :: P
    type(block_grid_t), intent(in) :: G
    integer, intent(in) :: block_num

    integer :: comm, ierr
    real(kind=wp) :: local_min, local_max, global_min, global_max
    real(kind=wp) :: best_dist_local, best_dist_global, tol_dist
    integer :: best_idx_local, best_idx_candidate, best_idx_global
    integer :: axis_id
    integer :: normal_comp
    integer :: logical_max
    logical :: participates

    comm = G%C%comm

    call axis_to_params(P%axis, axis_id, P%fixed_dir, normal_comp, P%n1, P%n2, P%use_logical_coord, G%C)

    if (P%fixed_dir == 0) then
      P%active = .false.
      return
    end if

    if (P%use_logical_coord) then
      logical_max = fixed_dir_extent(G%C, P%fixed_dir)
      best_idx_global = nint(P%coord)

      if (abs(P%coord - real(best_idx_global, wp)) > 1.0e-6_wp) then
        P%active = .false.
        return
      end if

      if (best_idx_global < 1 .or. best_idx_global > logical_max) then
        P%active = .false.
        return
      end if
    else
      call block_minmax(G, normal_comp, local_min, local_max)

      call MPI_Allreduce(local_min, global_min, 1, MPI_REAL_PW, MPI_MIN, comm, ierr)
      call MPI_Allreduce(local_max, global_max, 1, MPI_REAL_PW, MPI_MAX, comm, ierr)

      if (P%coord < global_min .or. P%coord > global_max) then
        P%active = .false.
        return
      end if

      call find_best_fixed_index(G, normal_comp, P%coord, best_dist_local, best_idx_local)

      call MPI_Allreduce(best_dist_local, best_dist_global, 1, MPI_REAL_PW, MPI_MIN, comm, ierr)

      tol_dist = max(1.0e-10_wp, 1.0e-6_wp * max(1.0_wp, best_dist_global))

      if (abs(best_dist_local - best_dist_global) <= tol_dist) then
        best_idx_candidate = best_idx_local
      else
        best_idx_candidate = huge(1)
      end if

      call MPI_Allreduce(best_idx_candidate, best_idx_global, 1, MPI_INTEGER, MPI_MIN, comm, ierr)

      if (best_idx_global == huge(1)) then
        P%active = .false.
        return
      end if
    end if

    P%fixed_index_global = best_idx_global

    participates = owns_fixed_index(G%C, P%fixed_dir, P%fixed_index_global)

    if (participates) then
      call MPI_Comm_split(comm, 1, G%C%rank, P%plane_comm, ierr)
      call MPI_Comm_rank(P%plane_comm, P%plane_rank, ierr)
      call MPI_Comm_size(P%plane_comm, P%plane_size, ierr)
    else
      call MPI_Comm_split(comm, MPI_UNDEFINED, G%C%rank, P%plane_comm, ierr)
      P%plane_rank = -1
      P%plane_size = 0
    end if

    P%active = .true.
    ! Note: active means plane is valid for this block; participation is via plane_comm.

    if (participates .and. P%plane_rank == 0) then
      call open_plane_file(domain_name, P, block_num, axis_id)
    end if

  end subroutine init_one_plane


  subroutine open_plane_file(domain_name, P, block_num, axis_id)

    use, intrinsic :: iso_fortran_env, only : int32

    implicit none

    character(*), intent(in) :: domain_name
    type(plane_output_plane), intent(inout) :: P
    integer, intent(in) :: block_num
    integer, intent(in) :: axis_id

    integer(int32) :: version32, axis32, n132, n232, fixed32

    character(len=256) :: filename
    character(len=512) :: filepath

    if (len_trim(P%name) > 0) then
      write(filename,'(a,a,a,i0,a)') trim(adjustl(domain_name)) // '_', trim(adjustl(P%name)), '_block', block_num, '.plane'
    else
      write(filename,'(a,a,i0,a)') trim(adjustl(domain_name)) // '_plane_block', block_num, '.plane'
    end if

    if (len_trim(P%file_directory) > 0 .and. trim(adjustl(P%file_directory)) /= '.') then
      call execute_command_line('mkdir -p "' // trim(adjustl(P%file_directory)) // '"')
      filepath = trim(adjustl(P%file_directory)) // '/' // trim(filename)
    else
      filepath = trim(filename)
    end if

    open(newunit=P%file_unit, file=filepath, access='stream', form='unformatted', status='replace', action='write')

    version32 = 1_int32
    axis32 = int(axis_id, int32)
    n132 = int(P%n1, int32)
    n232 = int(P%n2, int32)
    fixed32 = int(P%fixed_index_global, int32)

    write(P%file_unit) version32
    write(P%file_unit) axis32
    write(P%file_unit) n132
    write(P%file_unit) n232
    write(P%file_unit) P%coord
    write(P%file_unit) fixed32

    P%header_written = .true.

  end subroutine open_plane_file


  subroutine write_plane_output(planes, t, F, G)

    implicit none

    type(plane_output_type), intent(inout) :: planes
    real(kind=wp), intent(in) :: t
    real(kind=wp), dimension(:,:,:,:), allocatable, intent(in) :: F
    type(block_grid_t), intent(in) :: G

    integer :: p

    if (.not. planes%enabled) return
    if (.not. allocated(planes%P)) return

    do p = 1, planes%nplanes
      call write_one_plane(planes%P(p), t, F, G)
    end do

  end subroutine write_plane_output


  subroutine write_one_plane(P, t, F, G)

    implicit none

    type(plane_output_plane), intent(inout) :: P
    real(kind=wp), intent(in) :: t
    real(kind=wp), dimension(:,:,:,:), allocatable, intent(in) :: F
    type(block_grid_t), intent(in) :: G

    integer :: ierr
    integer :: m1, p1, m2, p2
    integer :: nloc1, nloc2
    integer :: src, meta(4)
    integer :: sender
    integer :: status(MPI_STATUS_SIZE)

    real(kind=wp), allocatable :: vx_loc(:,:), vy_loc(:,:), vz_loc(:,:)
    real(kind=wp), allocatable :: tmp(:,:)

    if (.not. P%enabled) return
    if (.not. P%active) return

    P%step_counter = P%step_counter + 1
    if (mod(P%step_counter - 1, P%stride) /= 0) return

    if (P%plane_comm == MPI_COMM_NULL) return

    call local_plane_bounds(G%C, P%fixed_dir, P%fixed_index_global, m1, p1, m2, p2)

    nloc1 = p1 - m1 + 1
    nloc2 = p2 - m2 + 1

    allocate(vx_loc(nloc1, nloc2), vy_loc(nloc1, nloc2), vz_loc(nloc1, nloc2))

    call extract_local_plane(F, G%C, P%fixed_dir, P%fixed_index_global, m1, p1, m2, p2, vx_loc, vy_loc, vz_loc)

    if (P%plane_rank == 0) then

      if (.not. allocated(P%vx)) then
        allocate(P%vx(P%n1, P%n2), P%vy(P%n1, P%n2), P%vz(P%n1, P%n2))
      end if

      call place_patch(P%fixed_dir, P%vx, vx_loc, m1, p1, m2, p2)
      call place_patch(P%fixed_dir, P%vy, vy_loc, m1, p1, m2, p2)
      call place_patch(P%fixed_dir, P%vz, vz_loc, m1, p1, m2, p2)

      do src = 1, P%plane_size - 1
        call MPI_Recv(meta, 4, MPI_INTEGER, MPI_ANY_SOURCE, tag_meta, P%plane_comm, status, ierr)
        sender = status(MPI_SOURCE)

        nloc1 = meta(2) - meta(1) + 1
        nloc2 = meta(4) - meta(3) + 1

        allocate(tmp(nloc1, nloc2))
        call MPI_Recv(tmp, nloc1*nloc2, MPI_REAL_PW, sender, tag_vx, P%plane_comm, MPI_STATUS_IGNORE, ierr)
        call place_patch(P%fixed_dir, P%vx, tmp, meta(1), meta(2), meta(3), meta(4))

        call MPI_Recv(tmp, nloc1*nloc2, MPI_REAL_PW, sender, tag_vy, P%plane_comm, MPI_STATUS_IGNORE, ierr)
        call place_patch(P%fixed_dir, P%vy, tmp, meta(1), meta(2), meta(3), meta(4))

        call MPI_Recv(tmp, nloc1*nloc2, MPI_REAL_PW, sender, tag_vz, P%plane_comm, MPI_STATUS_IGNORE, ierr)
        call place_patch(P%fixed_dir, P%vz, tmp, meta(1), meta(2), meta(3), meta(4))

        deallocate(tmp)
      end do

      write(P%file_unit) t
      write(P%file_unit) P%vx
      write(P%file_unit) P%vy
      write(P%file_unit) P%vz

    else

      meta = [m1, p1, m2, p2]
      call MPI_Send(meta, 4, MPI_INTEGER, 0, tag_meta, P%plane_comm, ierr)
      call MPI_Send(vx_loc, size(vx_loc), MPI_REAL_PW, 0, tag_vx, P%plane_comm, ierr)
      call MPI_Send(vy_loc, size(vy_loc), MPI_REAL_PW, 0, tag_vy, P%plane_comm, ierr)
      call MPI_Send(vz_loc, size(vz_loc), MPI_REAL_PW, 0, tag_vz, P%plane_comm, ierr)

    end if

    deallocate(vx_loc, vy_loc, vz_loc)

  end subroutine write_one_plane


  subroutine end_plane_output(planes)

    implicit none

    type(plane_output_type), intent(inout) :: planes

    integer :: p, ierr

    if (.not. allocated(planes%P)) return

    do p = 1, planes%nplanes

      if (planes%P(p)%plane_comm /= MPI_COMM_NULL) then
        call MPI_Comm_free(planes%P(p)%plane_comm, ierr)
        planes%P(p)%plane_comm = MPI_COMM_NULL
      end if

      if (planes%P(p)%file_unit > 0 .and. planes%P(p)%plane_rank == 0) then
        close(planes%P(p)%file_unit)
      end if

      if (allocated(planes%P(p)%vx)) deallocate(planes%P(p)%vx)
      if (allocated(planes%P(p)%vy)) deallocate(planes%P(p)%vy)
      if (allocated(planes%P(p)%vz)) deallocate(planes%P(p)%vz)

    end do

  end subroutine end_plane_output


  subroutine axis_to_params(axis, axis_id, fixed_dir, normal_comp, n1, n2, use_logical_coord, C)

    use mpi3dcomm, only : cartesian3d_t

    implicit none

    character(len=*), intent(in) :: axis
    integer, intent(out) :: axis_id
    integer, intent(out) :: fixed_dir
    integer, intent(out) :: normal_comp
    integer, intent(out) :: n1, n2
    logical, intent(out) :: use_logical_coord
    type(cartesian3d_t), intent(in) :: C

    character(len=2) :: a

    a = adjustl(axis)
    use_logical_coord = .false.

    select case (a)
    case ('yz')
      axis_id = 1
      fixed_dir = 1
      normal_comp = 1
      n1 = C%nr
      n2 = C%ns
    case ('xz')
      axis_id = 2
      fixed_dir = 2
      normal_comp = 2
      n1 = C%nq
      n2 = C%ns
    case ('xy')
      axis_id = 3
      fixed_dir = 3
      normal_comp = 3
      n1 = C%nq
      n2 = C%nr
    case ('rs')
      axis_id = 1
      fixed_dir = 1
      normal_comp = 1
      n1 = C%nr
      n2 = C%ns
      use_logical_coord = .true.
    case ('qs')
      axis_id = 2
      fixed_dir = 2
      normal_comp = 2
      n1 = C%nq
      n2 = C%ns
      use_logical_coord = .true.
    case ('qr')
      axis_id = 3
      fixed_dir = 3
      normal_comp = 3
      n1 = C%nq
      n2 = C%nr
      use_logical_coord = .true.
    case default
      axis_id = 0
      fixed_dir = 0
      normal_comp = 0
      n1 = 0
      n2 = 0
    end select

  end subroutine axis_to_params


  integer function fixed_dir_extent(C, fixed_dir)

    use mpi3dcomm, only : cartesian3d_t

    implicit none

    type(cartesian3d_t), intent(in) :: C
    integer, intent(in) :: fixed_dir

    select case (fixed_dir)
    case (1)
      fixed_dir_extent = C%nq
    case (2)
      fixed_dir_extent = C%nr
    case (3)
      fixed_dir_extent = C%ns
    case default
      fixed_dir_extent = 0
    end select

  end function fixed_dir_extent


  subroutine block_minmax(G, comp, local_min, local_max)

    implicit none

    type(block_grid_t), intent(in) :: G
    integer, intent(in) :: comp
    real(kind=wp), intent(out) :: local_min, local_max

    integer :: i, j, k

    local_min = huge(1.0_wp)
    local_max = -huge(1.0_wp)

    do k = G%C%ms, G%C%ps
      do j = G%C%mr, G%C%pr
        do i = G%C%mq, G%C%pq
          local_min = min(local_min, G%X(i,j,k,comp))
          local_max = max(local_max, G%X(i,j,k,comp))
        end do
      end do
    end do

  end subroutine block_minmax


  subroutine find_best_fixed_index(G, comp, coord, best_dist, best_idx)

    implicit none

    type(block_grid_t), intent(in) :: G
    integer, intent(in) :: comp
    real(kind=wp), intent(in) :: coord
    real(kind=wp), intent(out) :: best_dist
    integer, intent(out) :: best_idx

    integer :: i, j, k
    real(kind=wp) :: d

    best_dist = huge(1.0_wp)
    best_idx = huge(1)

    do k = G%C%ms, G%C%ps
      do j = G%C%mr, G%C%pr
        do i = G%C%mq, G%C%pq
          d = abs(G%X(i,j,k,comp) - coord)
          if (d < best_dist) then
            best_dist = d
            select case (comp)
            case (1)
              best_idx = i
            case (2)
              best_idx = j
            case (3)
              best_idx = k
            end select
          end if
        end do
      end do
    end do

  end subroutine find_best_fixed_index


  logical function owns_fixed_index(C, fixed_dir, fixed_index)

    use mpi3dcomm, only : cartesian3d_t

    implicit none

    type(cartesian3d_t), intent(in) :: C
    integer, intent(in) :: fixed_dir
    integer, intent(in) :: fixed_index

    select case (fixed_dir)
    case (1)
      owns_fixed_index = (C%mq <= fixed_index .and. fixed_index <= C%pq)
    case (2)
      owns_fixed_index = (C%mr <= fixed_index .and. fixed_index <= C%pr)
    case (3)
      owns_fixed_index = (C%ms <= fixed_index .and. fixed_index <= C%ps)
    case default
      owns_fixed_index = .false.
    end select

  end function owns_fixed_index


  subroutine local_plane_bounds(C, fixed_dir, fixed_index, m1, p1, m2, p2)

    use mpi3dcomm, only : cartesian3d_t

    implicit none

    type(cartesian3d_t), intent(in) :: C
    integer, intent(in) :: fixed_dir
    integer, intent(in) :: fixed_index
    integer, intent(out) :: m1, p1, m2, p2

    select case (fixed_dir)
    case (2) ! fixed r -> xz plane: vary q and s
      m1 = C%mq
      p1 = C%pq
      m2 = C%ms
      p2 = C%ps
    case (1) ! fixed q -> yz plane: vary r and s
      m1 = C%mr
      p1 = C%pr
      m2 = C%ms
      p2 = C%ps
    case (3) ! fixed s -> xy plane: vary q and r
      m1 = C%mq
      p1 = C%pq
      m2 = C%mr
      p2 = C%pr
    case default
      m1 = 1
      p1 = 0
      m2 = 1
      p2 = 0
    end select

  end subroutine local_plane_bounds


  subroutine extract_local_plane(F, C, fixed_dir, fixed_index, m1, p1, m2, p2, vx, vy, vz)

    use mpi3dcomm, only : cartesian3d_t

    implicit none

    real(kind=wp), dimension(:,:,:,:), allocatable, intent(in) :: F
    type(cartesian3d_t), intent(in) :: C
    integer, intent(in) :: fixed_dir
    integer, intent(in) :: fixed_index
    integer, intent(in) :: m1, p1, m2, p2
    real(kind=wp), intent(out) :: vx(:,:), vy(:,:), vz(:,:)

    integer :: a, b
    integer :: i, j, k

    select case (fixed_dir)
    case (2) ! xz: fixed j
      j = fixed_index
      do b = 1, size(vx,2)
        k = m2 + (b-1)
        do a = 1, size(vx,1)
          i = m1 + (a-1)
          vx(a,b) = F(i,j,k,1)
          vy(a,b) = F(i,j,k,2)
          vz(a,b) = F(i,j,k,3)
        end do
      end do

    case (1) ! yz: fixed i
      i = fixed_index
      do b = 1, size(vx,2)
        k = m2 + (b-1)
        do a = 1, size(vx,1)
          j = m1 + (a-1)
          vx(a,b) = F(i,j,k,1)
          vy(a,b) = F(i,j,k,2)
          vz(a,b) = F(i,j,k,3)
        end do
      end do

    case (3) ! xy: fixed k
      k = fixed_index
      do b = 1, size(vx,2)
        j = m2 + (b-1)
        do a = 1, size(vx,1)
          i = m1 + (a-1)
          vx(a,b) = F(i,j,k,1)
          vy(a,b) = F(i,j,k,2)
          vz(a,b) = F(i,j,k,3)
        end do
      end do

    case default
      vx = 0.0_wp
      vy = 0.0_wp
      vz = 0.0_wp
    end select

  end subroutine extract_local_plane


  subroutine place_patch(fixed_dir, global, patch, m1, p1, m2, p2)

    implicit none

    integer, intent(in) :: fixed_dir
    real(kind=wp), intent(inout) :: global(:,:)
    real(kind=wp), intent(in) :: patch(:,:)
    integer, intent(in) :: m1, p1, m2, p2

    global(m1:p1, m2:p2) = patch

  end subroutine place_patch

end module plane_output
