!> @brief block module defines blocks and subroutines that perform operations on blocks
!> @details this defines the properties of a block
!> most properties are obvious, but the block_boundary type deserves more discussion:
!> enforcing boundary and interface conditions will require use of material properties
!> and the coordinate transform (i.e., unit normal vector) on the boundary; the
!> block_boundary type has lower-rank arrays to store these for convenience in later use;
!> it also has fields that are defined on the boundary

module block

   use mpi
   use common, only : wp
   use datatypes, only : block_type, block_temp_parameters, block_pml
   use grid, only : block_grid_t
  use material, only : block_material, init_material, init_material_from_file, &
                init_anelastic_properties, init_anelastic_Q4_properties, &
                init_anelastic_Q8_properties, init_anelastic_Qf_properties, &
                init_const_Q_4M_properties, init_const_Q_8M_properties, &
                init_anelastic_Qf8_properties
   use plastic_material, only : block_plastic
   use fields, only : block_fields
   use boundary, only : block_boundary
   use BoundaryConditions, only : boundary_type
   use mpi3dio
   use mpi3dcomm
 
   implicit none
 
 contains
 
     subroutine init_block(mesh_source, type_of_mesh,material_source,response,fd_type, order_fd, interpol, &
       use_topography, topo, B, problem, btp, block_comm,infile,id,ny,nz,q4_config,q8_config,cq_config,fq8_config,process_dims,debug)
 
     !> @brief initialize a block
 
     use common, only : wp
     use mpi3dbasic, only : rank
     use grid, only : init_grid_cartesian, init_grid_curve, init_grid_from_file
    use material, only : init_material, init_material_from_file, init_anelastic_properties, &
         init_anelastic_Q4_properties, init_anelastic_Q8_properties, init_anelastic_Qf_properties
     use plastic_material, only : init_plastic_material
     use pml, only : init_pml
     use fields, only : init_fields
     use boundary, only : init_boundary
     use BoundaryConditions, only : init_boundaries
     use moment_tensor, only : init_moment_tensor
     use  inter_material, only : interpolatematerials
     use anelastic_q4_model, only : q4_parameters
     use anelastic_q8_model, only : q8_parameters
     use anelastic_cq_model, only : cq_parameters
     use anelastic_cq_material, only : init_anelastic_cq_properties
     use anelastic_fq8_model, only : fq8_parameters
     use decomposition_safety, only : stencil_requirements_t, get_stencil_requirements
     implicit none
 
     integer,intent(in) :: id, ny, nz, order_fd, process_dims(3)
     character(64), intent(in) :: type_of_mesh,problem, mesh_source, material_source, response, fd_type
     type(block_temp_parameters), intent(in) :: btp
     type(q4_parameters), intent(in) :: q4_config
     type(q8_parameters), intent(in) :: q8_config
     type(cq_parameters), intent(in) :: cq_config
     type(fq8_parameters), intent(in) :: fq8_config
     type(block_type),intent(out) :: B
     integer, intent(in) :: block_comm,infile
     logical, intent(in) :: interpol, use_topography
     real(kind = wp), intent(in) :: topo 
    logical, intent(in), optional :: debug
 
     integer :: size_in(3), npml, local_lqrs(3), local_rqrs(3)
     integer :: stat,order,i,j,k
     logical :: periodic(3)
     character(len=256) :: method
     logical :: use_moment_tensor
     type(stencil_requirements_t) :: stencil
 
 
     namelist /moment_list/ use_moment_tensor,order
 
     ! Set the block ID
     ! The blocks are numbered as 1,2,4,...2^N
     
     B%MT%block_id = id
     B%id = id
     stencil = get_stencil_requirements(fd_type, order_fd)
     if (.not.stencil%supported) error stop 'unsupported FD type/order reached init_block'
     B%nb = stencil%halo_width
    B%debug = .false.
    if (present(debug)) B%debug = debug
 
     ! Prep for initializing moment tensor source, if applicable
     use_moment_tensor = .FALSE.
     order = 2
     rewind(infile)
     read(infile,nml=moment_list,iostat=stat)
     if (stat>0) stop 'error reading namelist moment_list' 
     B%MT%use_moment_tensor = use_moment_tensor
     B%MT%order = order
 
     ! set parameters (really, read these from input file)
     ! set block physics (elastic, acoustic, etc.)
 
     B%physics = 'elastic'
 
     ! set number of grid points
 
     B%I%nq = btp%nqrs(1)
     B%I%nr = btp%nqrs(2)
     B%I%ns = btp%nqrs(3)
 
     B%rho_s_p = btp%rho_s_p
 
     !print*,  btp%pml_lqrs,  btp%pml_rqrs
     B%fd_type =  fd_type
     B%order =  order_fd
     method = 'manual'
     periodic = (/ .false., .false., .false. /)
     size_in = process_dims
     
     call partition_block(B, method, periodic, size_in, block_comm)
     
     ! initialize grid, material properties, fields
     ! (allocate arrays and set initial conditions)
     if (mesh_source == 'file') then
        call init_grid_from_file(B%id, B%G, btp%aqrs, btp%bqrs, B%nb)
     else if (mesh_source == 'compute') then
        if (type_of_mesh == 'cartesian') then
          call init_grid_cartesian(B%G, btp%aqrs, btp%bqrs,B%nb, B%debug)
        else if (type_of_mesh == 'curvilinear') then
           call init_grid_curve(B%G, btp%aqrs, btp%bqrs, &
                btp%lc, btp%rc, btp%profile_type, btp%profile_path, &
             btp%topography_type, btp%topography_path, .false., use_topography, topo, ny, nz,B%nb, B%debug)
        end if
     end if
 
     ! Initialize materials properties, either from file or based on problem type
     if (material_source == 'file') then
         call init_material_from_file(B%id, B%nb, B%M, B%G,btp%material_path)
     else
         call init_material(B%M,B%G,B%I,B%physics,problem, btp%rho_s_p, B%nb)
     end if
 
     
     if (interpol .eqv. .true.) call interpolatematerials(B,btp)

      if (trim(response) == 'low-pass' .or. trim(response) == 'anelastic') then
        call init_anelastic_properties(B%M, B%G, infile)
      else
        B%M%anelastic = .false.
      end if

      if (trim(response) == 'anelastic-Q4') then
        call init_anelastic_Q4_properties(B%M, B%G, infile, q4_config)
      else
        B%M%anelastic_Q = .false.
      end if

      if (trim(response) == 'anelastic-Q8' .or. trim(response) == 'anelastic-cQ8-b2') then
        call init_anelastic_Q8_properties(B%M, B%G, infile, q8_config)
      else
        B%M%anelastic_Q8 = .false.
      end if

      if (trim(response) == 'anelastic-cQ') then
        call init_anelastic_cq_properties(B%M,B%G,cq_config,id)
      else
        B%M%anelastic_cQ=.false.
      end if

      if (trim(response) == 'anelastic-Qf') then
        call init_anelastic_Qf_properties(B%M, B%G, infile)
      else
        B%M%anelastic_Qf = .false.
      end if

      if (trim(response) == 'constant-Q-4M') then
        call init_const_Q_4M_properties(B%M, B%G, infile)
      else
        B%M%anelastic_const_Q_4M = .false.
      end if
      
      if (trim(response) == 'constant-Q-8M') then
        call init_const_Q_8M_properties(B%M, B%G, infile)
      else
        B%M%anelastic_const_Q_8M = .false.
      end if

      if (trim(response) == 'frequency-Q-4M') then
        call init_anelastic_Qf_properties(B%M, B%G, infile)
      else
        B%M%anelastic_Qf = .false.
      end if

      if (trim(response) == 'anelastic-fQ8') then
        call init_anelastic_Qf8_properties(B%M, B%G, fq8_config)
      else
        B%M%anelastic_Qf8 = .false.
      end if
     
     if(response == 'plastic') call init_plastic_material(B%P,B%G,B%I,problem,btp%mu_beta_eta)
 
     if (B%MT%use_moment_tensor) then
         call init_moment_tensor(B%G,B%MT,infile)
      end if
 
 
     B%PMLB(1)%pml =  btp%pml_lqrs(1)
     B%PMLB(3)%pml =  btp%pml_lqrs(2)
     B%PMLB(5)%pml =  btp%pml_lqrs(3)
 
     B%PMLB(2)%pml =  btp%pml_rqrs(1)
     B%PMLB(4)%pml =  btp%pml_rqrs(2)
     B%PMLB(6)%pml =  btp%pml_rqrs(3)
 
     npml = btp%npml
     
     call init_pml(B%PMLB(1), B%G, 1, npml)
     call init_pml(B%PMLB(2), B%G, 2, npml)
     call init_pml(B%PMLB(3), B%G, 3, npml)
     call init_pml(B%PMLB(4), B%G, 4, npml)
     call init_pml(B%PMLB(5), B%G, 5, npml)
     call init_pml(B%PMLB(6), B%G, 6, npml)
     
     
     call init_fields(B%F,B%G%C,B%physics)
     call exchange_fields_block(B)
 
     
 
     call init_boundary(B%B, B%G, B%physics)
     call copy_materials_to_boundary(B)
     
     ! Physical boundary conditions apply only on the global exterior.  A
     ! Cartesian subdomain face with an MPI neighbor is an internal halo face,
     ! even though every rank received the same block-level boundary codes.
     local_lqrs = btp%lqrs
     local_rqrs = btp%rqrs
     if (B%G%C%rank_mq /= MPI_PROC_NULL) local_lqrs(1) = 0
     if (B%G%C%rank_pq /= MPI_PROC_NULL) local_rqrs(1) = 0
     if (B%G%C%rank_mr /= MPI_PROC_NULL) local_lqrs(2) = 0
     if (B%G%C%rank_pr /= MPI_PROC_NULL) local_rqrs(2) = 0
     if (B%G%C%rank_ms /= MPI_PROC_NULL) local_lqrs(3) = 0
     if (B%G%C%rank_ps /= MPI_PROC_NULL) local_rqrs(3) = 0
     call init_boundaries(B%boundary_vars, B%G%C, local_lqrs, local_rqrs)
     
     B%tau0 = 1.0_wp/(0.315949074074074_wp)
     
     B%fd_type =  fd_type
     B%order =  order_fd

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 2) then
       B%tau0 =  4.0_wp
     end if 

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 3) then
       B%tau0 =  12.0_wp/5.0_wp
     end if 

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 4) then
       B%tau0 =  144.0_wp/49.0_wp
     end if 

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 5) then
       B%tau0 =  720.0_wp/251.0_wp
     end if 

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 6) then
       B%tau0 =  42300.0_wp/13613.0_wp
     end if

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 7) then
       B%tau0 =  60480.0_wp/19087.0_wp  
     end if

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 8) then
       B%tau0 =  25401600.0_wp/7489399.0_wp  
     end if

     if (B%fd_type .eq. 'upwind' .and. B%order .eq. 9) then
       B%tau0 = 3628800.0_wp/1070017.0_wp
     end if

     !drp schemes

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 3) then
       B%tau0 =  2.4559485579992_wp
     end if

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 4) then
       B%tau0 =  3.16870508082731_wp
     end if

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 5) then
       B%tau0 =  3.14387257518354_wp
     end if

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 6) then
       B%tau0 =  3.39644507697204_wp
     end if

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 7) then
       B%tau0 =  3.39097726640256_wp
     end if

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 66) then
       B%tau0 =  1198.0_wp/377.0_wp
     end if

     if (B%fd_type .eq. 'upwind_drp' .and. B%order .eq. 679) then
       B%tau0 =  3.39644507697204_wp
     end if
    
     ! note that the grid, material properties, etc., are used
     ! when initializing the block_boundary type
 

!      print *,B%fd_type, B%order, B%tau0, B%nb 

!      STOP
     
 
     
     
   end subroutine init_block
   
   subroutine partition_block(B, method, periodic, size_in, block_comm)
     type(block_type), intent(inout) :: B
     character(256), intent(inout) :: method
     logical, dimension(3), intent(in) :: periodic
     integer, dimension(3), intent(in) :: size_in
     integer, intent(in) :: block_comm
 
     integer :: nb, nF
 
     nF = 9
     nb = B%nb
     
     
     call decompose3d(B%G%C, B%I%nq, B%I%nr, B%I%ns, B%nb, &
          periodic(1), periodic(2), periodic(3), &
          block_comm, method, &
          size_in(1), size_in(2), size_in(3))
 
     call subarray(B%G%C%nq,B%G%C%nr,B%G%C%ns,B%G%C%mq,B%G%C%pq,&
       B%G%C%mr,B%G%C%pr,B%G%C%ms,B%G%C%ps,MPI_DOUBLE_PRECISION,B%G%C%array_w)
     call subarray(B%G%C%nq,B%G%C%nr,B%G%C%ns,B%G%C%mq,B%G%C%pq,&
       B%G%C%mr,B%G%C%pr,B%G%C%ms,B%G%C%ps,MPI_REAL,B%G%C%array_s)
 
   end subroutine partition_block
 
   subroutine copy_fields_to_boundary(B)
 
     type(block_type), intent(inout) :: B
 
     B%B(1)%F = B%F%F(B%G%C%mq,:,:,:)
     B%B(2)%F = B%F%F(B%G%C%pq,:,:,:)
 !     B%B(3)%F = B%F%F(:,B%G%C%mr,:,:)
 !     B%B(4)%F = B%F%F(:,B%G%C%pr,:,:)
 !     B%B(5)%F = B%F%F(:,:,B%G%C%ms,:)
 !     B%B(6)%F = B%F%F(:,:,B%G%C%ps,:)
 
   end subroutine copy_fields_to_boundary
 
   subroutine copy_materials_to_boundary(B)
 
     type(block_type), intent(inout) :: B
 
     B%B(1)%M = B%M%M(B%G%C%mq,:,:,:)
     B%B(2)%M = B%M%M(B%G%C%pq,:,:,:)
 !     B%B(3)%M = B%M%M(:,B%G%C%mr,:,:)
 !     B%B(4)%M = B%M%M(:,B%G%C%pr,:,:)
 !     B%B(5)%M = B%M%M(:,:,B%G%C%ms,:)
 !     B%B(6)%M = B%M%M(:,:,B%G%C%ps,:)
 
   end subroutine copy_materials_to_boundary
 
   subroutine scale_rates_block(B,A)
 
     !> @brief scale rates by RK coefficient A
 
     use fields, only : scale_rates_interior
 
     implicit none
 
     type(block_type),intent(inout) :: B
     real(kind = wp),intent(in) :: A
 
     call scale_rates_interior(B, A)
 
   end subroutine scale_rates_block
 
 
   subroutine set_rates_block(B, type_of_mesh)
 
     !> @brief set rates using PDE
 
     use elastic, only : set_rates_elastic
 
     implicit none
 
     type(block_type),intent(inout) :: B
     character(len=64) :: type_of_mesh
 
     ! separate routines are needed to handle different physics
     ! since governing PDE will be different
 
     select case(B%physics)
     case default
        stop 'invalid block physics in set_rates_block'
     case('elastic')
        ! arguments are representative of actual case
        call set_rates_elastic(B, B%G, B%M, type_of_mesh)
     case('acoustic')
        stop 'acoustic physics not yet implemented in set_rates_block'
     end select
 
 
     !> @brief add SAT penalty terms
     !> @todo (eventually put this in separate subroutine since expressions will become much longer)
     !> @details note how the penalty terms are added to the rates array, but only along the boundaries,
     !> and that the penalty terms can be evaluated using only arrays stored in block_boundary type
     !> (that's part of the rationale for defining the block_boundary type)
 
   end subroutine set_rates_block
 
 
   subroutine update_fields_block(B,dt)
 
     !> @brief use rates to update fields
 
     use fields, only : update_fields_interior
 
     implicit none
 
     type(block_type),intent(inout) :: B
     real(kind = wp),intent(in) :: dt
 
     call update_fields_interior(B,dt)
 
   end subroutine update_fields_block
 
   subroutine exchange_fields_block(B)
 
     use fields, only : exchange_all_fields
 
     implicit none
 
     type(block_type),intent(inout) :: B
 
     call exchange_all_fields(B%F, B%G%C)
 
   end subroutine exchange_fields_block
 
 
   subroutine enforce_bound_conditions(B, mms_vars, t)
 
     !> @brief enforce boundary conditions on block boundaries
     !> @details the general idea is to
 
     use datatypes, only : mms_type
     use RHS_Interior, only: impose_boundary_condition
 
     implicit none
 
     type(block_type),intent(inout) :: B
     real(kind = wp),intent(in) :: t
 
     type(mms_type), intent(inout) :: mms_vars
 
     ! first copy fields from section of interior array to boundary array
     ! (this isn't really necessary, but the idea is to put everything that is
     ! required for boundary and interface conditions into the block_boundary type)
 
     ! it might be desirable to separate out this "preparation" step
     ! from enforcing the actual boundary conditions
 
     ! now enforce boundary conditions: the end result of this is setting the
     ! hat variables, which are fields defined on the boundaries that exactly satisfy
     ! the boundary or interface conditions
 
     ! right now this is done for all boundaries, but we would want some variable
     ! that could be used to identify when a block boundary is part of an interface,
     ! in which case an interface condition, rather than a boundary condition, will
     ! be enforced (in our 2D code, we use a character string to specify the boundary condition,
     ! like a free surface or rigid wall, and it can be set to 'none' for an interface)
 
     call impose_boundary_condition(B, mms_vars, t)
 
   end subroutine enforce_bound_conditions
 
 
   real(kind = wp) function block_time_step(spatialsetep, CFL, rho_s_p)
 
     implicit none
 
     real(kind = wp), intent(in) :: spatialsetep(3)
     real(kind = wp), intent(in) ::  rho_s_p(3)
     real(kind = wp), intent(in) :: CFL
     real(kind = wp) :: h
 
     ! calculate time step
 
     h = min(spatialsetep(1), spatialsetep(2), spatialsetep(3))
 
     block_time_step = cfl*h/sqrt(rho_s_p(2)**2 + rho_s_p(3)**2)
 
   end function block_time_step
   
 
   subroutine eval_block_mms(B, t, mms_vars)
 
     use mms, only : mms_type, eval_mms
 
     implicit none
 
     type(block_type), intent(inout) :: B
     real(kind = wp), intent(in) :: t
     type(mms_type), intent(in) :: mms_vars
 
 
     call eval_mms(B%F%DF, t, B%G%X, mms_vars, B%G%C%mq, B%G%C%mr, B%G%C%ms, B%G%C%pq, B%G%C%pr, B%G%C%ps)
 
 
   end subroutine eval_block_mms
 
   subroutine norm_fields_block(B)
 
     use fields, only : norm_fields
 
     implicit none
 
     type(block_type), intent(inout) :: B
 
 
     call norm_fields(B%F%F, B%G%C%mq, B%G%C%mr, B%G%C%ms, B%G%C%pq, B%G%C%pr, B%G%C%ps, 9, B%sum)
 
 
   end subroutine norm_fields_block
 
 end module block
 
