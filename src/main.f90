!----------------------------------------------------------
! WAVEQLAB3D
!----------------------------------------------------------
!
!> @author Eric Dunham
!> @author Kenneth Duru
!> @author Hari Radhakrishnan
!> @copyright Stanford University

!> @file main.f90

!> @fn main
!> @brief This is the main program.
!> @details Its purpose is to initialize a domain
!> (a set of blocks coupled across interfaces, with fields, grid,
!> material properties, etc. defined both within blocks and on interfaces)
!> and then evolve the fields according to some partial differential equation
!> (e.g., elastic or acoustic wave equation) with appropriate boundary and
!> interface conditions.

program main

  use common, only : wp
  use mpi3dbasic, only: start_mpi, finish_mpi, is_master
  !use hdf5_output
  use domain, only : domain_type,init_domain, close_domain, &
                     eval_mms, norm_fields, write_fault_output
  use simulation_config, only : simulation_config_t
  use input_preflight, only : preflight_input, print_preflight_time_parameters
  use time_step, only : RK_type,init_RK,time_step_RK
  implicit none

  ! domain_type contains all properties of a domain, including fields
  ! defined at a particular time t; evolving the solution in time is
  ! an operation done to the domain

  type(domain_type) :: D
  type(simulation_config_t) :: resolved_config

  ! time stepping parameters: these will be specific to the time-stepping
  ! method being used. we typically use a low-storage Runge-Kutta scheme
  ! (either 3rd or 4th order accurate) but the idea is that any time-stepping
  ! method could be used to advance the solution defined in the domain_type

  type(RK_type) :: RK ! Runge-Kutta coefficients
  !integer,parameter :: nt = 11 ! number of time steps
  real(kind = wp),parameter :: dt = 0.1_wp ! time step

  integer :: n, i

  ! Variables for parsing the input file
  character(len=256) :: ifname !< Input filename
  character(len=32) :: preflight_env
  integer :: infile, len_ifname, stat
  logical :: preflight_only
  real(kind = wp) :: t1, t2     ! timing information
  real(kind = wp) :: step_s, progress_pct




  !---------------------------------------------------------------------------
  !                       END OF VARIABLES
  !---------------------------------------------------------------------------

  ! initialize

  call start_mpi ! Initialize MPI for parallel processing
  !call start_hdf_output() ! Initialize the parallel I/O using HF5

  ! Read the command-line to get the input filename
  call get_command_argument(1, ifname, length=len_ifname, status=stat)
  if (stat /= 0) stop ':: Problem reading input filename from commandline.'

  call preflight_input(ifname, resolved_config)

  preflight_env = ''
  call get_environment_variable('WAVEQLAB3D_PREFLIGHT_ONLY', preflight_env, status=stat)
  preflight_only = stat == 0 .and. &
       any(trim(adjustl(preflight_env)) == [character(len=5) :: '1', 'true', 'TRUE', 'yes', 'YES'])
  if (preflight_only) then
     if (is_master()) call print_preflight_time_parameters(resolved_config)
     call finish_mpi(.false.)
     stop
  end if

  ! Read the input file
  open(newunit=infile, file=ifname, iostat=stat, status='old')
  if (stat/=0) stop ':: Cannot open input file'

  call init_RK(RK) ! Runge-Kutta coefficients
  call init_domain(D, infile, ifname, resolved_config) ! domain

  !call create_dataset("seismogram1", 0, [0])
  !call create_dataset("seismogram2", 0, [0])

  ! loop over all time steps, advancing solution according to wave equation
  ! and boundary and interface conditions; output is done within the time-
  ! stepping routine
  do n = 1, D%nt
    call cpu_time(t1)
    call time_step_RK(D,D%dt,RK, n)
    
    if (is_master()) then 
      call cpu_time(t2)
      step_s = t2 - t1
      progress_pct = 100.0_wp * real(n, kind=wp) / real(D%nt, kind=wp)

      if (n == 1) then
        write(*,'(A6,1X,A20,1X,A20,1X,A6)') 'n', 't', 'T(s)', '%'
      end if
      write(*,'(I6,1X,ES20.12E3,1X,F20.16,1X,F6.1)') n, D%t, step_s, progress_pct
    end if

     if (D%mms_vars%use_mms) then

        ! Analytical solution
        call eval_mms(D)

        ! norm of error
        call norm_fields(D)

        !         if (my_id == 0) then
        do i = 1, D%nblocks
          ! On multi-rank runs only the locally owned block has field storage.
          if (.not.allocated(D%B(i)%F%F)) cycle
          print *, 'MMS error in block ', i, ': ', &
               D%B(i)%sum*sqrt(D%B(i)%G%hq * D%B(i)%G%hr * D%B(i)%G%hs), ' at time = ', D%t
        end do
        !            write(30, *) t, sum1*sqrt(hq1*hr1*hs1), sum2*sqrt(hq2*hr2*hs2)
        !         end if

     else
        ! timing information
        ! call cpu_time(t2)

        ! status update to stdout

        !         if (my_id == 0) then
        !    print *, 'n = ',n,' of ', nt,' time steps, t = ',t,', wall clock time / time step = ',t2-t1
        !print*,t,handles%Uhat_pluspres(26,26,3),V(1,26,26,3)-U(mx1,26,26,3) !,W(26,201)
        !            print*,t,handles%Uhat_pluspres(11,11,3), Theta(11,11)
        !      print *, n, D%t, D%B(1)%G%mr, D%B(1)%G%pr, D%B(1)%G%ms, D%B(1)%G%ps, &
        !         D%handles%Uhat_pluspres(1,(D%B(1)%G%mr+D%B(1)%G%pr)/2,(D%B(1)%G%ms+D%B(1)%G%ps)/2,3)
        !             print *, n, D%t, D%fault%Uhat_pluspres(1,6,6,3),  D%fault%Uhat_pluspres(1,6,16,3), &
        !               D%fault%Uhat_pluspres(1,11,11,3), D%fault%Uhat_pluspres(1,16,6,3),  D%fault%Uhat_pluspres(1,16,16,3)
        !         end if
        !    print *, 'coordinates : ', X1(mx1, 76, 201, :), X2(1, 76, 201, :)
     end if
  end do

  call close_domain(D)
  !call finish_hdf_output()

  call finish_mpi

end program main
