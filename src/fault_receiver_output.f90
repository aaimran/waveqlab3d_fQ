module fault_receiver_output

  use mpi
  use common, only : wp
  use mpi3dbasic, only : MPI_REAL_PW, error
  use datatypes, only : fault_receiver_output_type, fault_type, iface_type, block_grid_t

  implicit none

contains

  subroutine init_fault_receiver_output(input, name, R, I, F, G, comm)
    integer, intent(in) :: input, comm
    character(*), intent(in) :: name
    type(fault_receiver_output_type), intent(inout) :: R
    type(iface_type), intent(in) :: I
    type(fault_type), intent(in) :: F
    type(block_grid_t), intent(in) :: G

    logical :: enable_fault_receivers
    integer :: number_fault_receivers, fault_receiver_stride
    real(kind=wp) :: fault_dip_degrees
    character(len=32) :: fault_receiver_coordinate_mode
    character(len=256) :: fault_receiver_directory
    namelist /fault_receiver_output_list/ enable_fault_receivers, &
         number_fault_receivers, fault_receiver_stride, &
         fault_receiver_coordinate_mode, fault_receiver_directory, &
         fault_dip_degrees

    integer :: stat, n, ierr, winner, candidate, j, k, best_j, best_k
    integer :: pq, unit
    real(kind=wp) :: a, b, c, dip, d2, local_best, global_best
    character(len=64) :: receiver_name
    character(len=512) :: line, path
    logical :: found

    enable_fault_receivers = .false.
    number_fault_receivers = 0
    fault_receiver_stride = 1
    fault_receiver_coordinate_mode = 'strike_down_dip'
    fault_receiver_directory = 'fault_receivers'
    fault_dip_degrees = 15.0_wp

    rewind(input)
    read(input, nml=fault_receiver_output_list, iostat=stat)
    if (stat > 0) call error('invalid &fault_receiver_output_list', &
                             'init_fault_receiver_output')
    if (.not.enable_fault_receivers) return
    if (number_fault_receivers < 1) call error( &
         'number_fault_receivers must be positive', 'init_fault_receiver_output')

    R%enabled = .true.
    R%nreceivers = number_fault_receivers
    R%stride = max(1, fault_receiver_stride)
    R%directory = trim(adjustl(fault_receiver_directory))
    if (len_trim(R%directory) == 0) R%directory = '.'
    R%comm = comm
    call MPI_Comm_rank(comm, R%comm_rank, ierr)

    allocate(R%name(R%nreceivers), R%target_xyz(3,R%nreceivers), &
             R%actual_xyz(3,R%nreceivers), R%distance(R%nreceivers), &
             R%owner(R%nreceivers), R%j(R%nreceivers), R%k(R%nreceivers), &
             R%file_unit(R%nreceivers))
    R%file_unit = -1

    call seek_receiver_list(input, found)
    if (.not.found) call error('fault receiver list markers not found', &
                               'init_fault_receiver_output')
    dip = fault_dip_degrees*acos(-1.0_wp)/180.0_wp
    do n = 1,R%nreceivers
       read(input,'(a)',iostat=stat) line
       if (stat /= 0 .or. trim(adjustl(line)) == '!---end:fault_receiver_list---') &
            call error('fault receiver list is shorter than number_fault_receivers', &
                       'init_fault_receiver_output')
       select case(trim(adjustl(fault_receiver_coordinate_mode)))
       case('strike_down_dip')
          read(line,*,iostat=stat) receiver_name, a, b
          if (stat /= 0) call error('expected: name strike_km down_dip_km', &
                                    'init_fault_receiver_output')
          R%target_xyz(:,n) = (/ b*cos(dip), b*sin(dip), a /)
       case('xyz')
          read(line,*,iostat=stat) receiver_name, a, b, c
          if (stat /= 0) call error('expected: name x_km y_km z_km', &
                                    'init_fault_receiver_output')
          R%target_xyz(:,n) = (/ a, b, c /)
       case default
          call error('fault_receiver_coordinate_mode must be strike_down_dip or xyz', &
                     'init_fault_receiver_output')
       end select
       R%name(n) = trim(adjustl(receiver_name))
    end do

    if (R%comm_rank == 0 .and. trim(R%directory) /= '.') &
         call execute_command_line('mkdir -p "' // trim(R%directory) // '"')
    call MPI_Barrier(comm, ierr)

    pq = G%C%pq
    do n = 1,R%nreceivers
       local_best = huge(1.0_wp)
       best_j = G%C%mr
       best_k = G%C%ms
       do k = G%C%ms,G%C%ps
          do j = G%C%mr,G%C%pr
             d2 = sum((G%x(pq,j,k,1:3)-R%target_xyz(:,n))**2)
             if (d2 < local_best) then
                local_best = d2
                best_j = j
                best_k = k
             end if
          end do
       end do
       call MPI_Allreduce(local_best, global_best, 1, MPI_REAL_PW, MPI_MIN, comm, ierr)
       candidate = huge(1)
       if (abs(local_best-global_best) <= &
           max(1.0e-20_wp,1.0e-12_wp*max(1.0_wp,global_best))) candidate = R%comm_rank
       call MPI_Allreduce(candidate, winner, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
       R%owner(n) = winner
       R%j(n) = best_j
       R%k(n) = best_k
       R%distance(n) = sqrt(global_best)
       R%actual_xyz(:,n) = 0.0_wp
       if (R%comm_rank == winner) R%actual_xyz(:,n) = G%x(pq,best_j,best_k,1:3)
       call MPI_Bcast(R%actual_xyz(:,n), 3, MPI_REAL_PW, winner, comm, ierr)

       if (R%comm_rank == winner) then
          path = trim(R%name(n)) // '.txt'
          if (trim(R%directory) /= '.') path = trim(R%directory) // '/' // trim(path)
          open(newunit=unit, file=trim(path), status='replace', action='write', iostat=stat)
          if (stat /= 0) call error('cannot open fault receiver output ' // trim(path), &
                                    'init_fault_receiver_output')
          R%file_unit(n) = unit
          write(unit,'(a)') '# code=waveqlab3d_Q'
          write(unit,'(a,a)') '# simulation=', trim(name)
          write(unit,'(a,a)') '# receiver=', trim(R%name(n))
          write(unit,'(a,3(es24.16,1x))') '# target_xyz_km=', R%target_xyz(:,n)
          write(unit,'(a,3(es24.16,1x))') '# mapped_xyz_km=', R%actual_xyz(:,n)
          write(unit,'(a,es24.16)') '# mapping_distance_km=', R%distance(n)
          write(unit,'(a)') '# units: s m m/s MPa m m/s MPa MPa'
          write(unit,'(a)') 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress'
       end if
    end do
  end subroutine init_fault_receiver_output


  subroutine seek_receiver_list(input, found)
    integer, intent(in) :: input
    logical, intent(out) :: found
    character(len=512) :: line
    integer :: stat
    found = .false.
    rewind(input)
    do
       read(input,'(a)',iostat=stat) line
       if (stat /= 0) exit
       if (trim(adjustl(line)) == '!---begin:fault_receiver_list---') then
          found = .true.
          exit
       end if
    end do
  end subroutine seek_receiver_list


  subroutine write_fault_receiver_output(R, t, I, F)
    type(fault_receiver_output_type), intent(inout) :: R
    real(kind=wp), intent(in) :: t
    type(iface_type), intent(in) :: I
    type(fault_type), intent(in) :: F
    integer :: n, j, k, u
    real(kind=wp) :: row(8)
    if (.not.R%enabled) return
    if (mod(R%step_counter,R%stride) == 0) then
       do n = 1,R%nreceivers
          if (R%comm_rank /= R%owner(n)) cycle
          j = R%j(n); k = R%k(n); u = R%file_unit(n)
          row = (/ t, &
               I%S(j,k,3), &
               I%Svel(j,k,3), &
               F%Uhat_pluspres(j,k,6), &
               I%S(j,k,2), &
               I%Svel(j,k,2), &
               F%Uhat_pluspres(j,k,5), &
               F%Uhat_pluspres(j,k,4) /)
          write(u,'(8(es24.16,1x))') row
       end do
    end if
    R%step_counter = R%step_counter+1
  end subroutine write_fault_receiver_output


  subroutine end_fault_receiver_output(R)
    type(fault_receiver_output_type), intent(inout) :: R
    integer :: n
    if (.not.R%enabled) return
    do n = 1,R%nreceivers
       if (R%comm_rank == R%owner(n) .and. R%file_unit(n) /= -1) close(R%file_unit(n))
    end do
    R%enabled = .false.
  end subroutine end_fault_receiver_output

end module fault_receiver_output
