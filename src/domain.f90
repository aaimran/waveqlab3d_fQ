!> @brief domain module holds multiple blocks and interfaces, and is used by time-stepping routines

!> @todo Hari, I've never been careful about public vs private in modules, maybe you can help with that?

!> domain = set of blocks and interfaces, with fields defined at time t
!> blocks and interfaces are allocatable to allow for any number of them,
!> but initialization is specific to two blocks and one interface (which I think
!> should be the initial focus since that is sufficient to do many interesting problems),
!> but we should generalize this to an arbitrary number of blocks and interfaces
!> with an unstructured connectivity pattern (presumably then the domain_type would also
!> contain some type that would handle the block/interface layout)

module domain

  use mpi
  use common, only : wp
  use datatypes, only : domain_type, block_type, iface_type, mms_type,fault_type,moment_tensor
  use block, only : block_temp_parameters, eval_block_mms, norm_fields_block
  use slice_output, only : init_slice_output, end_slice_output
  use seismogram
  use plane_output, only : init_plane_output, write_plane_output, end_plane_output
  use fault_receiver_output, only : init_fault_receiver_output, &
       write_fault_receiver_output, end_fault_receiver_output

  implicit none

  logical :: in_block_comm(2), in_fault_comm(2), in_hslice_comm(2), in_vslice_comm(2)

contains


  subroutine init_domain(D, infile, input_filename, config)

    !> initialize the domain at t=0

    use block, only : init_block, block_time_step
    use iface, only : init_iface
    use fault_output, only : init_fault_output
    use mms, only : init_mms
    use interface_condition, only : init_vel_state
    use mpi3dbasic, only: is_master, rank, nprocs, new_communicator, warning, print_run_info
    use mpi3d_interface
    use simulation_config, only : simulation_config_t
    use anelastic_q4_model, only : q4_relaxation_dt_limit
    use anelastic_q8_model, only : q8_relaxation_dt_limit
    use anelastic_q8_model, only : q8_parameters
    use anelastic_cq8_b2_model, only : cq8_b2_relaxation_dt_limit
    use anelastic_cq_model, only : cq_relaxation_dt_limit
    use anelastic_fq8_model, only : fq8_relaxation_dt_limit

    implicit none

    type(domain_type),intent(out) :: D
    integer, intent(in) :: infile
    character(*), intent(in), optional :: input_filename
    type(simulation_config_t), intent(in) :: config


    integer :: i, cart_size(3)
    integer :: n, ierr, cart_rank, comm_cart, coord(3), normal(3), &
               block_comms(2), fault_comms(2), hslice_comms(2), &
               vslice_comms(2), my_block_comm
    logical :: periodic(3) = (/ .false., .false., .false. /), reorder=.true.
    type(interface3d) :: II

    ! Variables for parsing the input data file
    integer :: stat
    character(256) :: name, problem,response,plastic_model,fd_type
    character(256) :: response_norm
    logical :: invalid_response
    integer :: nt, nblocks
    logical :: output_fault_topo, w_fault, interpol, use_topography, mollify_source
    integer :: w_stride, ny, nz, order
    real(kind = wp) :: CFL, t_final, topo
    ! interface conditions
    character(64) :: coupling !< locked, slip-weakening_friction, linear_friction
    character(64) :: mesh_source, type_of_mesh, material_source !< cartesian or curvilinear

    real(kind = wp) :: dt, dtmin, dt1, dt2, spatialsetep1(3), spatialsetep2(3)
    real(kind = wp) :: dt_i, spat(3), elastic_dt_limit, relaxation_dt_limit
    
    real(kind = wp) :: meshvolume1, meshvolume2, totalvolume, ratio1, ratio2
    integer :: nprocs_1, nprocs_2

    type(block_temp_parameters) :: btp(2)
    type(q8_parameters) :: selected_q8


    !---------------------------------------------------------------------------
    !                       END OF VARIABLES
    !---------------------------------------------------------------------------
    ! Consume the root-parsed, validated, and broadcast configuration.
    name = config%problem%name
    problem = config%problem%problem
    response = config%problem%response
    response_norm = config%problem%response
    plastic_model = config%problem%plastic_model
    nblocks = config%problem%nblocks
    nt = config%problem%nt
    CFL = config%problem%CFL
    coupling = config%problem%coupling
    t_final = config%problem%t_final
    w_fault = config%problem%w_fault
    output_fault_topo = .true.
    mesh_source = config%problem%mesh_source
    type_of_mesh = config%problem%type_of_mesh
    material_source = config%problem%material_source
    w_stride = config%problem%w_stride
    interpol = config%problem%interpol
    use_topography = config%problem%use_topography
    mollify_source = config%problem%mollify_source
    topo = config%problem%topo
    fd_type = config%problem%fd_type
    order = config%problem%order
    btp(1:nblocks) = config%blocks

    if (is_master() .and. present(input_filename)) then
      call print_run_info(input_filename, name, problem)
      write(*,'(/,/)')
    end if

    D%name = name
    D%problem = problem
    D%response = response_norm
    D%plastic_model = plastic_model
    D%coupling = coupling
    D%w_fault = w_fault
    D%output_fault_topo = output_fault_topo
    D%t_final = t_final
    D%type_of_mesh = type_of_mesh
    D%mesh_source = mesh_source
    D%material_source = material_source
    D%CFL = CFL
    D%w_stride = w_stride
    D%interpol = interpol
    D%use_topography = use_topography
    D%topo = topo
    D%mollify_source = mollify_source 
    D%fd_type = fd_type
    D%order = order
    D%debug = .false.  ! Default off; enable for diagnostics when needed
    ! Defensive fallback only: a one-block input need not specify btp(2), but
    ! legacy setup code that inspects it still sees a valid, initialized value.
    ! Runtime field code must nevertheless use only the allocated D%B(1).
    if (nblocks == 1) btp(2) = btp(1)


    
       

    !> set initial time

    D%t = 0.0_wp

    ! compute global time-step using only the blocks specified by `nblocks`
     dtmin = huge(1.0_wp)
     do i = 1, nblocks
       spat(1) = (btp(i)%bqrs(1)-btp(i)%aqrs(1))/real(btp(i)%nqrs(1)-1)
       spat(2) = (btp(i)%bqrs(2)-btp(i)%aqrs(2))/real(btp(i)%nqrs(2)-1)
       spat(3) = (btp(i)%bqrs(3)-btp(i)%aqrs(3))/real(btp(i)%nqrs(3)-1)
       dt_i = block_time_step(spat, CFL, btp(i)%rho_s_p)
       dtmin = min(dtmin, dt_i)
     end do

     elastic_dt_limit = dtmin
     relaxation_dt_limit = huge(1.0_wp)
     if (trim(response_norm) == 'anelastic-Q4') then
       relaxation_dt_limit = q4_relaxation_dt_limit(config%q4)
       dtmin = min(dtmin, relaxation_dt_limit)
     else if (trim(response_norm) == 'anelastic-Q8') then
       relaxation_dt_limit = q8_relaxation_dt_limit(config%q8)
       dtmin = min(dtmin, relaxation_dt_limit)
     else if (trim(response_norm) == 'anelastic-cQ8-b2') then
       relaxation_dt_limit = cq8_b2_relaxation_dt_limit(config%cq8_b2)
       dtmin = min(dtmin, relaxation_dt_limit)
     else if (trim(response_norm) == 'anelastic-cQ') then
       relaxation_dt_limit = cq_relaxation_dt_limit(config%cq)
       dtmin = min(dtmin, relaxation_dt_limit)
     else if (trim(response_norm) == 'anelastic-fQ8') then
       relaxation_dt_limit = fq8_relaxation_dt_limit(config%fq8)
       dtmin = min(dtmin, relaxation_dt_limit)
     end if
     D%dt = dtmin

        if (is_master()) call warning('warning: current method for setting time step does not use material properties; ' // &
          'mesh information from files is ignored (only scalars from btp%... in input file are used)', 'init_domain')
    
    !> set number of blocks and interfaces
    !> (now supports 1 or 2 blocks; read from input file)
    
    D%nblocks = nblocks
    D%nifaces = max(0, nblocks - 1)  ! 0 for 1-block, 1 for 2-blocks

    ! MPI distribution resolved during preflight.
    if (D%nblocks == 1) then
      in_block_comm(1) = .true.
      in_block_comm(2) = .false.
    else if (config%serial_shared_blocks) then
      in_block_comm = .true.
    else
      in_block_comm(1) = rank >= config%rank_begin(1) .and. rank <= config%rank_end(1)
      in_block_comm(2) = rank >= config%rank_begin(2) .and. rank <= config%rank_end(2)
    end if

    call new_communicator(in_block_comm(1), block_comms(1))
    call new_communicator(in_block_comm(2), block_comms(2))

    if (in_block_comm(1)) my_block_comm = block_comms(1)
    if (in_block_comm(2)) my_block_comm = block_comms(2)

    allocate(D%B(D%nblocks))
    allocate(D%I(D%nifaces))
    allocate(D%seismometers(D%nblocks))
    allocate(D%plane_outputs(D%nblocks))

    !> initialize blocks
    do i = 1,D%nblocks
       ny = btp(i)%nqrs(2)
       nz = btp(i)%nqrs(3)

    if((btp(i)%profile_type == 'read_from_memomry_fractal') .and. &
         (btp(i)%topography_type == 'read_topo_from_memory')) then

       ny = btp(i)%faultsize(1)
       nz = btp(i)%faultsize(2)

    end if
       
      if (.not.in_block_comm(i)) cycle
      selected_q8 = config%q8
      if (trim(response_norm) == 'anelastic-cQ8-b2') selected_q8 = config%cq8_b2%block(i)
      call init_block(D%mesh_source, D%type_of_mesh, D%material_source,&
           D%response, D%fd_type,  D%order, D%interpol, D%use_topography, topo, D%B(i), &
         problem, btp(i), block_comms(i),infile,i, ny, nz, config%q4, selected_q8, config%cq, config%fq8, &
         config%process_dims(i,:), D%debug)

      cart_size = [D%B(i)%G%C%size_q,D%B(i)%G%C%size_r,D%B(i)%G%C%size_s]
      coord = D%B(i)%G%C%coord

      
    end do


    D%nt = floor(D%t_final/dtmin)
   
    if (is_master()) then
      write (*,*) name  
      write (*,*) "Domain time parameters :"
      write (*,*) "        Number of time steps = ", D%nt
      write (*,*) "        Initial dt = ", D%dt
      write (*,*) "        Elastic CFL limit = ", elastic_dt_limit
      if (trim(response_norm) == 'anelastic-Q4') &
        write (*,*) "        Q4 relaxation limit = ", relaxation_dt_limit
      if (trim(response_norm) == 'anelastic-Q8') &
        write (*,*) "        Q8 relaxation limit = ", relaxation_dt_limit
      if (trim(response_norm) == 'anelastic-cQ8-b2') &
        write (*,*) "        cQ8-b2 relaxation limit = ", relaxation_dt_limit
      if (trim(response_norm) == 'anelastic-cQ') &
        write (*,*) "        cQ relaxation limit = ", relaxation_dt_limit
      write (*,*) "        Final time = ", D%t_final
    end if

    ! Setup communicator for interface partitioning.
    ! Only create interface communicators for 2-block mode (when nifaces > 0)

    if (D%nifaces > 0) then
      if (in_block_comm(1)) then
        normal = [1, 0, 0]
        call new_interface(coord, cart_size, normal, MPI_COMM_WORLD, D%B(1)%G%C, II)
      end if
      if (in_block_comm(2).and. .not.in_block_comm(1)) then
        normal = [-1, 0, 0]
        call new_interface(coord, cart_size, normal, MPI_COMM_WORLD, D%B(2)%G%C, II)
      end if
    end if

    call validate_initialized_geometry(D, II, config%serial_shared_blocks)

    !> initialize interfaces (only if nifaces > 0, i.e., for 2-block mode)

    do i = 1, D%nifaces
       !> create an interface that joins two blocks
       !> first define which blocks will be coupled: block 1 to block 2
       D%I(i)%im = 1
       D%I(i)%ip = 2
       !> and in which direction they will be coupled: q
       D%I(i)%direction = 'q'
       !> then initialize interface
       if (in_block_comm(1)) call init_iface(D%I(i),D%B(1), II)
       if (in_block_comm(2) .and. .not.in_block_comm(1)) call init_iface(D%I(i),D%B(2), II)
       if (in_block_comm(1)) call init_vel_state(D%problem, D%I(i), D%B(1), -1.0_wp, 2)
       if (in_block_comm(2) .and. .not.in_block_comm(1)) call init_vel_state(D%problem, D%I(i), D%B(2), 1.0_wp, 1)
    end do

    ! Exchange materials across interface (only for 2-block mode)
    if (D%nifaces > 0) call exchange_materials_interface(D)

    ! Fault output initialization (only for 2-block mode with interface)
    if (D%nifaces > 0) then
      in_fault_comm(1) = (in_block_comm(1) .and. D%I(1)%II%on_interface)
      in_fault_comm(2) = (in_block_comm(2) .and. D%I(1)%II%on_interface)
      call new_communicator(in_fault_comm(1), fault_comms(1))
      call new_communicator(in_fault_comm(2), fault_comms(2))

      if (in_fault_comm(1)) call init_fault_output(infile, D%w_fault,D%name, D%fault, D%B(1)%G%C, fault_comms(1))
      if (in_fault_comm(2) .and. .not.in_fault_comm(1)) call init_fault_output(infile, D%w_fault,D%name, D%fault, D%B(2)%G%C, fault_comms(2))
      if (in_fault_comm(1)) call init_fault_receiver_output(infile, D%name, &
           D%fault_receivers, D%I(1), D%fault, D%B(1)%G, fault_comms(1))
    else
      ! 1-block mode: no interface communicators
      in_fault_comm(1) = .false.
      in_fault_comm(2) = .false.
    end if

    call init_mms(infile, D%mms_vars)

    ! Initialize seismometers and plane outputs for each block (generalized for 1 or 2 blocks)
    do i = 1, D%nblocks
      D%seismometers(i)%block_num = i
      if (in_block_comm(i)) call init_seismogram(infile, D%seismometers(i), D%name, D%B(i)%G)
      if (in_block_comm(i)) call init_plane_output(infile, D%name, D%plane_outputs(i), D%B(i)%G, i)
    end do
    
    !> above assumes that blocks contain fields defined on a structured mesh,
    !> which is a cube in the computational domain; the computational coordinates
    !> are q,r,s, and there is a mapping to a curvilinear mesh in physical coordinates x,y,z
    !> I recommend that we use q,r,s in all references to directions, sides, etc., because
    !> it is possible that a boundary in q is mapped to a boundary in y or z rather than x
    !> (this will be a different naming convention than in original version of 3D code)

  end subroutine init_domain


  subroutine validate_initialized_geometry(D, II, serial_shared_blocks)
    use mpi3d_interface, only : interface3d
    use diagnostics, only : fatal_local
    use mpi3dbasic, only : rank
    implicit none

    type(domain_type), intent(in), target :: D
    type(interface3d), intent(in) :: II
    logical, intent(in) :: serial_shared_blocks
    integer :: ierr, i, j, k, c, qface, count, position
    integer :: local_shape(2), peer_shape(2)
    real(wp) :: local_mismatch, global_mismatch, local_scale, global_scale
    real(wp) :: local_min_j, global_min_j, tolerance
    real(wp), allocatable :: send_face(:), recv_face(:)
    type(block_type), pointer :: active
    logical :: on_interface

    local_min_j = huge(1.0_wp)
    do i = 1, D%nblocks
       if (.not.allocated(D%B(i)%G%J)) cycle
       local_min_j = min(local_min_j, minval(D%B(i)%G%J( &
            D%B(i)%G%C%mq:D%B(i)%G%C%pq, D%B(i)%G%C%mr:D%B(i)%G%C%pr, &
            D%B(i)%G%C%ms:D%B(i)%G%C%ps)))
    end do
    call MPI_Allreduce(local_min_j, global_min_j, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, MPI_COMM_WORLD, ierr)
    if (global_min_j <= 0.0_wp .and. rank == 0) call fatal_local('RUN-GRID-001', &
         'Initialized grid contains a nonpositive Jacobian.', 'validate_initialized_geometry')

    if (D%nifaces == 0) return
    local_mismatch = 0.0_wp
    local_scale = 1.0_wp

    if (serial_shared_blocks) then
       do c = 1, 3
          local_mismatch = max(local_mismatch, maxval(abs( &
               D%B(1)%G%x(D%B(1)%G%C%pq,D%B(1)%G%C%mr:D%B(1)%G%C%pr, &
                          D%B(1)%G%C%ms:D%B(1)%G%C%ps,c) - &
               D%B(2)%G%x(D%B(2)%G%C%mq,D%B(2)%G%C%mr:D%B(2)%G%C%pr, &
                          D%B(2)%G%C%ms:D%B(2)%G%C%ps,c))))
          local_scale = max(local_scale, maxval(abs(D%B(1)%G%x( &
               D%B(1)%G%C%pq,D%B(1)%G%C%mr:D%B(1)%G%C%pr, &
               D%B(1)%G%C%ms:D%B(1)%G%C%ps,c))))
       end do
    else
       nullify(active)
       on_interface = .false.
       if (in_block_comm(1)) then
          active => D%B(1)
          on_interface = active%G%C%coord(1) == active%G%C%size_q-1
          qface = active%G%C%pq
       else if (in_block_comm(2)) then
          active => D%B(2)
          on_interface = active%G%C%coord(1) == 0
          qface = active%G%C%mq
       end if

       if (on_interface) then
          local_shape = [active%G%C%lnr, active%G%C%lns]
          call MPI_Sendrecv(local_shape, 2, MPI_INTEGER, II%rank_neighbor, 31, &
               peer_shape, 2, MPI_INTEGER, II%rank_neighbor, 31, II%comm, &
               MPI_STATUS_IGNORE, ierr)
          if (any(local_shape /= peer_shape)) call fatal_local('RUN-IFACE-003', &
               'Paired interface ranks have different owned face shapes.', &
               'validate_initialized_geometry')
          count = 3*product(local_shape)
          allocate(send_face(count), recv_face(count))
          position = 0
          do c = 1, 3
             do k = active%G%C%ms, active%G%C%ps
                do j = active%G%C%mr, active%G%C%pr
                   position = position + 1
                   send_face(position) = active%G%x(qface,j,k,c)
                end do
             end do
          end do
          call MPI_Sendrecv(send_face, count, MPI_DOUBLE_PRECISION, II%rank_neighbor, 32, &
               recv_face, count, MPI_DOUBLE_PRECISION, II%rank_neighbor, 32, II%comm, &
               MPI_STATUS_IGNORE, ierr)
          local_mismatch = maxval(abs(send_face-recv_face))
          local_scale = max(local_scale, maxval(abs(send_face)), maxval(abs(recv_face)))
          deallocate(send_face, recv_face)
       end if
    end if

    call MPI_Allreduce(local_mismatch, global_mismatch, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_scale, global_scale, 1, MPI_DOUBLE_PRECISION, &
         MPI_MAX, MPI_COMM_WORLD, ierr)
    tolerance = 1000.0_wp*epsilon(1.0_wp)*global_scale
    if (global_mismatch > tolerance .and. rank == 0) call fatal_local('RUN-IFACE-004', &
         'Post-grid interface coordinates differ pointwise.', &
         'validate_initialized_geometry')
  end subroutine validate_initialized_geometry


  subroutine close_domain(D)

    use fault_output, only: destroy_fault
    use material, only : destroy_anelastic_Q4_properties, destroy_anelastic_Q8_properties, &
         destroy_anelastic_Qf8_properties
    use anelastic_cq_material, only : destroy_anelastic_cq_properties
    use diagnostics, only : fatal_local
    use mpi3dbasic, only : rank
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    type(domain_type),intent(inout) :: D
    real(wp) :: local_eta, global_eta, local_field, global_field
    integer :: i, ierr
    logical :: local_finite, global_finite

    if (trim(D%response) == 'anelastic-Q4' .or. trim(D%response) == 'anelastic-Q8' .or. &
        trim(D%response) == 'anelastic-cQ8-b2' .or. &
        trim(D%response) == 'anelastic-cQ' .or. &
        trim(D%response) == 'anelastic-fQ8') then
       local_eta = 0.0_wp
       local_field = 0.0_wp
       local_finite = .true.
       do i = 1, D%nblocks
          if (.not.allocated(D%B(i)%F%F)) cycle
          local_field = max(local_field, maxval(abs(D%B(i)%F%F)))
          if (trim(D%response) == 'anelastic-Q4' .and. allocated(D%B(i)%M%eta4Q)) then
             local_eta = max(local_eta, maxval(abs(D%B(i)%M%eta4Q)), &
                  maxval(abs(D%B(i)%M%eta5Q)), maxval(abs(D%B(i)%M%eta6Q)), &
                  maxval(abs(D%B(i)%M%eta7Q)), maxval(abs(D%B(i)%M%eta8Q)), &
                  maxval(abs(D%B(i)%M%eta9Q)))
          else if ((trim(D%response) == 'anelastic-Q8' .or. &
                    trim(D%response) == 'anelastic-cQ8-b2') .and. &
                   allocated(D%B(i)%M%eta4Q8)) then
             local_eta = max(local_eta, maxval(abs(D%B(i)%M%eta4Q8)), &
                  maxval(abs(D%B(i)%M%eta5Q8)), maxval(abs(D%B(i)%M%eta6Q8)), &
                  maxval(abs(D%B(i)%M%eta7Q8)), maxval(abs(D%B(i)%M%eta8Q8)), &
                  maxval(abs(D%B(i)%M%eta9Q8)))
          else if (trim(D%response) == 'anelastic-cQ' .and. allocated(D%B(i)%M%eta4cQ)) then
             local_eta=max(local_eta,maxval(abs(D%B(i)%M%eta4cQ)), &
                  maxval(abs(D%B(i)%M%eta5cQ)),maxval(abs(D%B(i)%M%eta6cQ)), &
                  maxval(abs(D%B(i)%M%eta7cQ)),maxval(abs(D%B(i)%M%eta8cQ)), &
                  maxval(abs(D%B(i)%M%eta9cQ)))
          else if (allocated(D%B(i)%M%eta4Qf8)) then
             local_eta = max(local_eta,maxval(abs(D%B(i)%M%eta4Qf8)), &
                  maxval(abs(D%B(i)%M%eta5Qf8)),maxval(abs(D%B(i)%M%eta6Qf8)), &
                  maxval(abs(D%B(i)%M%eta7Qf8)),maxval(abs(D%B(i)%M%eta8Qf8)), &
                  maxval(abs(D%B(i)%M%eta9Qf8)))
          end if
          local_finite = local_finite .and. ieee_is_finite(local_field) .and. &
               ieee_is_finite(local_eta)
       end do
       call MPI_Allreduce(local_eta, global_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
            MPI_COMM_WORLD, ierr)
       call MPI_Allreduce(local_field, global_field, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
            MPI_COMM_WORLD, ierr)
       call MPI_Allreduce(local_finite, global_finite, 1, MPI_LOGICAL, MPI_LAND, &
            MPI_COMM_WORLD, ierr)
       if (.not.global_finite .and. trim(D%response) == 'anelastic-Q4') &
            call fatal_local('RUN-Q4-001', &
            'Non-finite Q4 field or memory state detected at shutdown.', 'close_domain')
       if (.not.global_finite .and. trim(D%response) == 'anelastic-Q8') &
            call fatal_local('RUN-Q8-001', &
            'Non-finite Q8 field or memory state detected at shutdown.', 'close_domain')
       if (.not.global_finite .and. trim(D%response) == 'anelastic-cQ8-b2') &
            call fatal_local('RUN-CQ8-B2-001', &
            'Non-finite cQ8-b2 field or memory state detected at shutdown.', 'close_domain')
       if (.not.global_finite .and. trim(D%response) == 'anelastic-cQ') &
            call fatal_local('RUN-CQ-001', &
            'Non-finite cQ field or memory state detected at shutdown.','close_domain')
       if (.not.global_finite .and. trim(D%response) == 'anelastic-fQ8') &
            call fatal_local('RUN-FQ8-001', &
            'Non-finite fQ8 field or memory state detected at shutdown.', 'close_domain')
       if (rank == 0 .and. trim(D%response) == 'anelastic-Q4') &
            write(*,'(A,ES12.4,A,ES12.4)') &
            'Q4 final state: max|field|=', global_field, ', max|memory|=', global_eta
       if (rank == 0 .and. trim(D%response) == 'anelastic-Q8') &
            write(*,'(A,ES12.4,A,ES12.4)') &
            'Q8 final state: max|field|=', global_field, ', max|memory|=', global_eta
       if (rank == 0 .and. trim(D%response) == 'anelastic-cQ8-b2') &
            write(*,'(A,ES12.4,A,ES12.4)') &
            'cQ8-b2 final state: max|field|=', global_field, ', max|memory|=', global_eta
       if (rank == 0 .and. trim(D%response) == 'anelastic-cQ') &
            write(*,'(A,ES12.4,A,ES12.4)') &
            'cQ final state: max|field|=',global_field,', max|memory|=',global_eta
       if (rank == 0 .and. trim(D%response) == 'anelastic-fQ8') &
            write(*,'(A,ES12.4,A,ES12.4)') &
            'fQ8 final state: max|field|=', global_field, ', max|memory|=', global_eta
    end if

    if ( D%w_fault .eqv.  .true.) then
       if (in_fault_comm(1)) call destroy_fault(D%fault)
       if (in_fault_comm(2) .and. .not.in_fault_comm(1)) call destroy_fault(D%fault)
    end if
    if (in_fault_comm(1)) call end_fault_receiver_output(D%fault_receivers)
    !call end_slice_output(D%slicer)
    if (in_block_comm(1)) call destroy_seismogram(D%seismometers(1))
    if (in_block_comm(2)) call destroy_seismogram(D%seismometers(2))

    if (in_block_comm(1)) call end_plane_output(D%plane_outputs(1))
    if (in_block_comm(2)) call end_plane_output(D%plane_outputs(2))

    do i = 1, D%nblocks
       if (D%B(i)%M%anelastic_Q) call destroy_anelastic_Q4_properties(D%B(i)%M)
       if (D%B(i)%M%anelastic_Q8) call destroy_anelastic_Q8_properties(D%B(i)%M)
       if (D%B(i)%M%anelastic_cQ) call destroy_anelastic_cq_properties(D%B(i)%M)
       if (D%B(i)%M%anelastic_Qf8) call destroy_anelastic_Qf8_properties(D%B(i)%M)
    end do

  end subroutine close_domain

  subroutine enforce_bound_iface_conditions(D, stage)

    !> enforce boundary and interface conditions

    use block, only : enforce_bound_conditions

    implicit none

    type(domain_type),intent(inout) :: D
    integer, intent(in) :: stage

    !> enforce boundary conditions on external sides of each block
    !> enforce interface conditions to couple blocks (only for 2-block mode)

    if (in_block_comm(1)) then
      call enforce_bound_conditions(D%B(1), D%mms_vars, D%t)
      if (D%nifaces > 0) then
        call enforce_iface_conditions(D%problem, D%coupling, D%I(1), &
          D%B(1),2,D%t, stage, D%mms_vars, D%fault)
      end if
    end if
    if (in_block_comm(2)) then
      call enforce_bound_conditions(D%B(2), D%mms_vars, D%t)
      if (D%nifaces > 0) then
        call enforce_iface_conditions(D%problem, D%coupling, D%I(1), &
          D%B(2),1,D%t, stage, D%mms_vars, D%fault)
      end if
    end if

  end subroutine enforce_bound_iface_conditions

  subroutine enforce_iface_conditions(problem, coupling, I,B, ib, t, stage, mms_vars, handles)

    use datatypes, only : block_type, iface_type, mms_type, fault_type, block_boundary
    use RHS_Interior, only : Impose_Interface_Condition
    implicit none

    character(*), intent(in) :: problem, coupling
    type(iface_type),intent(inout) :: I
    type(block_type),intent(inout) :: B
    type(mms_type), intent(inout) :: mms_vars
    type(fault_type), intent(inout) :: handles
    real(kind = wp),intent(in) :: t
    integer, intent(in) :: stage, ib

    select case(I%direction)

    case('q')

!!! coupling: side 2 of block 1, B(2) <==> side 1 of block 2, Bp(1)

!!! this solves interface conditions for hat variables, constructs SAT forcing terms,
!!! and adds SAT forcing terms to rates

      if (B%boundary_vars%Rx == 0 .or. B%boundary_vars%Lx == 0) then
        call Impose_Interface_Condition(problem, coupling, I, B, ib, &
                              t, stage, mms_vars, handles)
      end if

    case('r','s')

       stop 'interfaces in r and s direction not implemented in enforce_iface_conditions'

    end select

  end subroutine enforce_iface_conditions

  subroutine exchange_fields(D)

    use block, only : exchange_fields_block

    implicit none

    type(domain_type),intent(inout) :: D
    integer :: i

    if (in_block_comm(1)) call exchange_fields_block(D%B(1))
    if (in_block_comm(2)) call exchange_fields_block(D%B(2))

  end subroutine exchange_fields

  subroutine exchange_fields_interface(D)

    use mpi3dbasic, only : nprocs
    use block, only : copy_fields_to_boundary
    use boundary, only : exchange_fields_across_interface

    implicit none
    type(domain_type),intent(inout) :: D

    ! Skip interface field exchange if no interfaces exist (1-block mode)
    if (D%nifaces == 0) return

    if (nprocs == 1) then
      call copy_fields_to_boundary(D%B(1))
      if (D%nblocks >= 2) then
        call copy_fields_to_boundary(D%B(2))
        D%B(1)%B(2)%Fopp(:,:,:) = D%B(2)%B(1)%F(:,:,:)
        D%B(2)%B(1)%Fopp(:,:,:) = D%B(1)%B(2)%F(:,:,:)
      end if
    else
      if (in_block_comm(1)) then
        call copy_fields_to_boundary(D%B(1))
        call exchange_fields_across_interface(D%B(1)%B(2), D%B(1)%G%C, D%I(1)%II)
      end if
      if (in_block_comm(2)) then
        call copy_fields_to_boundary(D%B(2))
        call exchange_fields_across_interface(D%B(2)%B(1), D%B(2)%G%C, D%I(1)%II)
      end if
    end if

  end subroutine exchange_fields_interface

  subroutine exchange_materials_interface(D)

    use mpi3dbasic, only : nprocs
    use block, only : copy_fields_to_boundary
    use boundary, only : exchange_materials_across_interface

    implicit none
    type(domain_type),intent(inout) :: D

    ! Skip interface material exchange if no interfaces exist (1-block mode)
    if (D%nifaces == 0) return

    if (nprocs == 1) then
      if (D%nblocks >= 2) then
        D%B(1)%B(2)%Mopp = D%B(2)%B(1)%M
        D%B(2)%B(1)%Mopp = D%B(1)%B(2)%M
      end if
    else
      if (in_block_comm(1)) then
        call exchange_materials_across_interface(D%B(1)%B(2), D%B(1)%G%C, D%I(1)%II)
      end if
      if (in_block_comm(2)) then
        call exchange_materials_across_interface(D%B(2)%B(1), D%B(2)%G%C, D%I(1)%II)
      end if
    end if

  end subroutine exchange_materials_interface

  subroutine write_output(D)

    !> @brief write fields (and, in some cases, rates) at time t
    !> note that the way it is written exposes the details of how fields are stored
    !> in blocks and on interfaces; a potentially better way is to use subroutines
    !> that retrieve certain fields and return them in an output array

    use,intrinsic :: iso_fortran_env, only : output_unit

    implicit none

    type(domain_type),intent(inout) :: D

    if ( D%w_fault .eqv.  .true.) then
       call write_fault_output(D)
       !call write_slice_output(D)
       call write_hat_output(D)
    end if

     if (in_block_comm(1)) call write_plane_output(D%plane_outputs(1), D%t, D%B(1)%F%F, D%B(1)%G)
     if (in_block_comm(2)) call write_plane_output(D%plane_outputs(2), D%t, D%B(2)%F%F, D%B(2)%G)
    
    call write_seismogram_output(D)

  end subroutine write_output

  subroutine write_fault_receivers_domain(D)
    implicit none
    type(domain_type), intent(inout) :: D
    if (D%nifaces < 1) return
    if (in_fault_comm(1)) call write_fault_receiver_output( &
         D%fault_receivers, D%t, D%I(1), D%fault)
  end subroutine write_fault_receivers_domain

  subroutine write_fault_output(D)

    use fault_output
    use mpi3dbasic, only : rank
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1, mq2, mr2, ms2, pq2, pr2, ps2

    if ( D%w_fault .eqv.  .true.) then
       if (in_fault_comm(1)) then
          
          mq1 = D%B(1)%G%C%mq
          mr1 = D%B(1)%G%C%mr
          ms1 = D%B(1)%G%C%ms
          pq1 = D%B(1)%G%C%pq
          pr1 = D%B(1)%G%C%pr
          ps1 = D%B(1)%G%C%ps
          
          call write_fault(D%B(1)%F%F(pq1,mr1:pr1,ms1:ps1,:),  &
               D%I(1)%S(mr1:pr1,ms1:ps1,:), &
               D%I(1)%W(mr1:pr1,ms1:ps1,:), D%fault)
          
       end if
       if (in_fault_comm(2)) then
          
          mq2 = D%B(2)%G%C%mq
          mr2 = D%B(2)%G%C%mr
          ms2 = D%B(2)%G%C%ms
          pq2 = D%B(2)%G%C%pq
          pr2 = D%B(2)%G%C%pr
          ps2 = D%B(2)%G%C%ps
          
          if (D%fault%output_enabled(2)) call write_file_distributed( &
               D%fault%handles(2), D%B(2)%F%F(mq2,mr2:pr2,ms2:ps2,:))
       end if
    end if


  end subroutine write_fault_output

  subroutine write_hat_output(D)
    use fault_output
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1

    if ( D%w_fault .eqv.  .true.) then

       if (.not.in_fault_comm(1)) return
       
       mq1 = D%B(1)%G%C%mq
       mr1 = D%B(1)%G%C%mr
       ms1 = D%B(1)%G%C%ms
       pq1 = D%B(1)%G%C%pq
       pr1 = D%B(1)%G%C%pr
       ps1 = D%B(1)%G%C%ps
       
       call write_hats(D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,:), &
            D%fault%Vhat_pluspres(mr1:pr1,ms1:ps1,:), &
            D%fault%Uhat_pluspres(mr1:pr1,ms1:ps1,1:3), &
            D%fault%time_rup(mr1:pr1,ms1:ps1,1),D%fault)
       
    end if
       
     end subroutine write_hat_output

  subroutine write_slice_output(D)

    use slice_output
    implicit none

    type(domain_type), intent(in) :: D

    if (in_block_comm(1)) call write_slice(D%B(1)%F%F,D%B(1)%G%C, D%slicer)

  end subroutine write_slice_output

  subroutine write_seismogram_output(D)

    use slice_output
    implicit none

    type(domain_type), intent(in) :: D
    integer :: mq1, mr1, ms1, pq1, pr1, ps1, mq2, mr2, ms2, pq2, pr2, ps2

    if(in_block_comm(1)) then
        mq1 = D%B(1)%G%C%mq
        mr1 = D%B(1)%G%C%mr
        ms1 = D%B(1)%G%C%ms
        pq1 = D%B(1)%G%C%pq
        pr1 = D%B(1)%G%C%pr
        ps1 = D%B(1)%G%C%ps

        call write_seismogram(D%seismometers(1), D%t, D%B(1)%F%F)

    end if

    if(in_block_comm(2)) then
        mq2 = D%B(2)%G%C%mq
        mr2 = D%B(2)%G%C%mr
        ms2 = D%B(2)%G%C%ms
        pq2 = D%B(2)%G%C%pq
        pr2 = D%B(2)%G%C%pr
        ps2 = D%B(2)%G%C%ps

        call write_seismogram(D%seismometers(2), D%t, D%B(2)%F%F) 

    end if


  end subroutine write_seismogram_output

  subroutine eval_mms(D)

    implicit none
    type(domain_type),intent(inout) :: D

    integer :: i

    do i = 1, D%nblocks
      if (.not.in_block_comm(i)) cycle
       call eval_block_mms(D%B(i), D%t, D%mms_vars)
    end do

  end subroutine eval_mms

  subroutine norm_fields(D)

    implicit none
    type(domain_type),intent(inout) :: D

    integer :: i

     do i = 1, D%nblocks
       if (in_block_comm(i)) call norm_fields_block(D%B(i))
     end do

  end subroutine norm_fields

  subroutine scale_rates(D,A)

    !> @brief multiply all rates by RK coefficient A

    use block, only : scale_rates_block
    use iface, only : scale_rates_iface

    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: A

    integer :: i

    ! first within blocks and on their boundaries

      do i = 1, D%nblocks
        if (in_block_comm(i)) call scale_rates_block(D%B(i),A)
      end do


    ! and then on interfaces

    do i = 1,D%nifaces
       call scale_rates_iface(D%I(i),A)
    end do

  end subroutine scale_rates

  subroutine set_rates(D)

    !> @brief set rates using the PDE
    !> @details set rates using the PDE (in a low storage RK method, the new rates
    !> at the current stage are added to the old rates, instead of overwriting
    !> the rates array)

    use block, only : set_rates_block
    use moment_tensor, only : set_moment_tensor, moment_tensor_body_force, set_moment_tensor_smooth
    use mpi3dbasic, only : rank

    implicit none

    type(domain_type),intent(inout) :: D
    integer :: i

    do i = 1, D%nblocks
      if (.not.in_block_comm(i)) cycle
      call set_rates_block(D%B(i), D%type_of_mesh)

      if (D%B(i)%MT%use_moment_tensor) then
        if (D%mollify_source) then
          call set_moment_tensor_smooth(D%B(i),D%t)
        else
          call set_moment_tensor(D%B(i),D%t)
        end if
      end if
    end do
  end subroutine set_rates


  subroutine update_fields(D,dt,stage,RKstage)

    !> @brief use rates to update fields

    use block, only : update_fields_block
    use iface, only : update_fields_iface
    use plastic, only : update_fields_plastic
    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: dt
    integer,intent(in) :: stage,RKstage

    integer :: i

    ! first within blocks and on their boundaries

      if (in_block_comm(1)) call update_fields_block(D%B(1),dt)
      if (in_block_comm(2)) call update_fields_block(D%B(2),dt)


    ! and then on interfaces

    do i = 1,D%nifaces
       call update_fields_iface(D%I(i),dt)
    end do

    if (stage==RKstage) then
       !if (in_block_comm(1)) call update_fields_plastic(D%B(1),D%B(1)%P,D%B(1)%G,D%B(1)%M,D%dt,D%t,D%problem,D%response)
       !if (in_block_comm(2)) call update_fields_plastic(D%B(2),D%B(2)%P,D%B(2)%G,D%B(2)%M,D%dt,D%t,D%problem,D%response)
       if (in_block_comm(1)) call update_fields_plastic(D%B(1),D%B(1)%P,D%B(1)%G,D%B(1)%M,D%dt,D%t,D%problem,D%response,&
            D%plastic_model)
       if (in_block_comm(2)) call update_fields_plastic(D%B(2),D%B(2)%P,D%B(2)%G,D%B(2)%M,D%dt,D%t,D%problem,D%response,&
            D%plastic_model)
    end if
  end subroutine update_fields

end module domain
