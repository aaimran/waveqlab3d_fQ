
!> This module contains all the derived datatypes used in WaveQLab.
!> Initially, these were distributed amongst the different source
!> files but that was causing circular dependencies. To maintain
!> source compatability with KD3D while killing the circular
!> dependencies, this approach was chosen.
module datatypes

   use common, only : wp
   use mpi3dcomm, only : cartesian3d_t
   use mpi3d_interface, only : interface3d
   use mpi3dio, only : file_distributed
 
   implicit none
 
   !> block_indices datatype originally created by Eric.
   !> Its functionality is redundant and should be deleted at some point.
   type :: block_indices
      integer :: nq,nr,ns !< number of grid points in each direction in computational coordinates q,r,s
      real(kind = wp) :: hq, hr, hs !< grid spacing along the 3 computational axes
   end type block_indices
 
   !> block_grid datatype that contains the grid positions, derivatives, and jacobians
   type :: block_grid_t
      type(cartesian3d_t) :: C
      real(kind = wp) :: hq, hr, hs, bhq, bhr, bhs      !< Grid spacing in unit cube (h*) and block (bh*)
      real(kind = wp),dimension(:,:,:,:), allocatable :: x !< coordinates
      real(kind = wp),dimension(:,:,:,:), allocatable :: metricx , &
                                         metricy , &
                                         metricz !< metric derivative
      real(kind = wp),dimension(:,:,:), allocatable :: J !> Jacobian
   end type block_grid_t
 
   !> block_material datatype to hold the material properties
   type :: block_material
      real(kind = wp),dimension(:,:,:,:), allocatable :: M !< material properties

      ! --- Anelastic (Q) / attenuation state (used only when response == 'anelastic')
      logical :: anelastic = .false.
      real(kind = wp) :: dt = 0.0_wp
      real(kind = wp) :: c = 1.0_wp
      real(kind = wp) :: weight_exp = 0.0_wp
      real(kind = wp) :: fref = 1.0_wp

      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv, Qs_inv
      real(kind = wp), dimension(4) :: tau = 0.0_wp
      real(kind = wp), dimension(4) :: weight = 0.0_wp

      ! Six stress-component memory variables (4 mechanisms each)
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4,eta5,eta6,eta7,eta8,eta9
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4,Deta5,Deta6,Deta7,Deta8,Deta9

      ! --- anelastic-Q4: four-mechanism broadband attenuation
      logical :: anelastic_Q = .false.
      integer :: n_mechanism_Q = 4
      real(kind = wp) :: fref_Q = 1.0_wp
      real(kind = wp) :: Qs0_Q = -1.0_wp
      real(kind = wp) :: Qp0_Q = -1.0_wp
      real(kind = wp) :: fmin_Q = -1.0_wp
      real(kind = wp) :: fmax_Q = -1.0_wp
      character(len=32) :: weight_method_Q = ''
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_Q, Qs_inv_Q
      real(kind = wp), dimension(4) :: tau_Q = 0.0_wp
      real(kind = wp), dimension(4) :: weight_Q = 0.0_wp
      ! Six stress-component memory variables for anelastic-Q4 (4 mechanisms each)
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4Q,eta5Q,eta6Q,eta7Q,eta8Q,eta9Q
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4Q,Deta5Q,Deta6Q,Deta7Q,Deta8Q,Deta9Q

      ! --- anelastic-Q8: 8-mechanism broadband attenuation (response == 'anelastic-Q8')
      logical :: anelastic_Q8 = .false.
      integer :: n_mechanism_Q8 = 8
      real(kind = wp) :: fref_Q8 = 1.0_wp
      real(kind = wp) :: Qs0_Q8 = -1.0_wp
      real(kind = wp) :: Qp0_Q8 = -1.0_wp
      real(kind = wp) :: fmin_Q8 = -1.0_wp
      real(kind = wp) :: fmax_Q8 = -1.0_wp
      character(len=32) :: weight_method_Q8 = ''
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_Q8, Qs_inv_Q8
      real(kind = wp), dimension(8) :: tau_Q8    = 0.0_wp
      real(kind = wp), dimension(8) :: weight_Q8 = 0.0_wp
      ! Six stress-component memory variables for anelastic-Q8 (8 mechanisms each)
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4Q8,eta5Q8,eta6Q8,eta7Q8,eta8Q8,eta9Q8
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4Q8,Deta5Q8,Deta6Q8,Deta7Q8,Deta8Q8,Deta9Q8

      ! --- independent configurable constant-Q response (response == 'anelastic-cQ')
      logical :: anelastic_cQ = .false.
      integer :: n_mechanism_cQ = 0
      real(kind = wp) :: fref_cQ = 1.0_wp
      real(kind = wp) :: Qs0_cQ = -1.0_wp, Qp0_cQ = -1.0_wp
      real(kind = wp) :: fmin_cQ = -1.0_wp, fmax_cQ = -1.0_wp
      character(len=32) :: coefficient_policy_cQ = ''
      real(kind = wp), dimension(:), allocatable :: tau_cQ, weight_s_cQ, weight_p_cQ
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_cQ, Qs_inv_cQ
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4cQ,eta5cQ,eta6cQ,eta7cQ,eta8cQ,eta9cQ
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4cQ,Deta5cQ,Deta6cQ,Deta7cQ,Deta8cQ,Deta9cQ

      ! --- anelastic-Qf: frequency-dependent Q (response == 'anelastic-Qf')
      logical :: anelastic_Qf = .false.
      integer :: n_mechanism_Qf = 4
      real(kind = wp) :: fref_Qf = 1.0_wp           ! Reference frequency for unrelaxed modulus correction (Hz)
      real(kind = wp) :: gamma_Qf = 0.0_wp          ! Power-law exponent (0.0-0.9)
      real(kind = wp) :: f_trans_Qf = 1.0_wp        ! Transition frequency (Hz)
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_Qf, Qs_inv_Qf
      real(kind = wp), dimension(4) :: tau_Qf = 0.0_wp
      real(kind = wp), dimension(4) :: weight_Qf = 0.0_wp
      ! Six stress-component memory variables for anelastic-Qf (4 mechanisms each)
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4Qf,eta5Qf,eta6Qf,eta7Qf,eta8Qf,eta9Qf
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4Qf,Deta5Qf,Deta6Qf,Deta7Qf,Deta8Qf,Deta9Qf

      ! --- constant-Q-4M: configurable 4-mechanism constant-Q attenuation (response == 'constant-Q-4M')
      logical :: anelastic_const_Q_4M = .false.
      integer :: n_mechanism_const_Q_4M = 4
      real(kind = wp) :: fref_const_Q_4M = 1.0_wp           ! Reference frequency for unrelaxed modulus correction (Hz)
      real(kind = wp) :: fmin_const_Q_4M = 0.05_wp          ! Minimum frequency (Hz) for tau-band
      real(kind = wp) :: fmax_const_Q_4M = 20.0_wp          ! Maximum frequency (Hz) for tau-band
      real(kind = wp) :: target_Q_const_Q_4M = 50.0_wp      ! Target constant-Q value
      character(16) :: weight_method_const_Q_4M = 'nnls'    ! 'nnls' | 'withers' | 'lookup'
      real(kind = wp), dimension(4) :: manual_weights_const_Q_4M = 0.0_wp  ! Manual weight override (if all > 0)
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_const_Q_4M, Qs_inv_const_Q_4M
      real(kind = wp), dimension(4) :: tau_const_Q_4M = 0.0_wp
      real(kind = wp), dimension(4) :: weight_const_Q_4M = 0.0_wp
      ! Six stress-component memory variables for constant-Q-4M (4 mechanisms each)
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4_4M,eta5_4M,eta6_4M,eta7_4M,eta8_4M,eta9_4M
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4_4M,Deta5_4M,Deta6_4M,Deta7_4M,Deta8_4M,Deta9_4M
      
      ! --- constant-Q-8M: configurable 8-mechanism constant-Q attenuation (response == 'constant-Q-8M')
      logical :: anelastic_const_Q_8M = .false.
      integer :: n_mechanism_const_Q_8M = 8
      real(kind = wp) :: fref_const_Q_8M = 1.0_wp
      real(kind = wp) :: fmin_const_Q_8M = 0.05_wp
      real(kind = wp) :: fmax_const_Q_8M = 20.0_wp
      real(kind = wp) :: target_Q_const_Q_8M = 50.0_wp
      character(16) :: weight_method_const_Q_8M = 'nnls'
      real(kind = wp), dimension(8) :: manual_weights_const_Q_8M = 0.0_wp
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_const_Q_8M, Qs_inv_const_Q_8M
      real(kind = wp), dimension(8) :: tau_const_Q_8M = 0.0_wp
      real(kind = wp), dimension(8) :: weight_const_Q_8M = 0.0_wp
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4_8M,eta5_8M,eta6_8M,eta7_8M,eta8_8M,eta9_8M
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4_8M,Deta5_8M,Deta6_8M,Deta7_8M,Deta8_8M,Deta9_8M

      ! --- frequency-dependent Q: canonical response anelastic-fQ8 (8 mechanisms)
      logical :: anelastic_Qf8 = .false.
      logical :: coarse_grained_Qf8 = .false.
      integer :: n_mechanism_Qf8 = 8
      real(kind = wp) :: fref_Qf8 = 1.0_wp
      real(kind = wp) :: gamma_Qf8 = 0.0_wp
      real(kind = wp) :: f_trans_Qf8 = 1.0_wp
      real(kind = wp), dimension(:,:,:), allocatable :: Qp_inv_Qf8, Qs_inv_Qf8
      real(kind = wp), dimension(8) :: tau_Qf8 = 0.0_wp
      real(kind = wp), dimension(8) :: strength_s_Qf8 = 0.0_wp
      real(kind = wp), dimension(8) :: strength_p_Qf8 = 0.0_wp
      real(kind = wp), dimension(:,:,:,:), allocatable :: eta4Qf8,eta5Qf8,eta6Qf8,eta7Qf8,eta8Qf8,eta9Qf8
      real(kind = wp), dimension(:,:,:,:), allocatable :: Deta4Qf8,Deta5Qf8,Deta6Qf8,Deta7Qf8,Deta8Qf8,Deta9Qf8
   end type block_material
 
   !> block_plastic datatype to hold Drucker-Prager plasticity variables
   type :: block_plastic
      real(kind = wp) :: mu_beta_eta(3) !< internal friction: mu, plastic dilatancy: beta and viscosity: eta
      real(kind = wp),dimension(:,:,:,:), allocatable :: P !< plastic array (plaststic strain, and strain rate) 
   end type block_plastic
 
   !> block_fields datatype holds the different velocity and stress fields,
   !> and their temporal derivatives
   type :: block_fields
      real(kind = wp),dimension(:,:,:,:), allocatable :: F !< fields array (velocities, stresses)
      real(kind = wp),dimension(:,:,:,:), allocatable :: DF !< rates array (stores time derivative of fields, dF/dt)     
   end type block_fields
   
   !> block_boundary datatype holds the different components of the 6 2D block boundaries.
   type :: block_boundary
      real(kind = wp),dimension(:,:,:), allocatable :: X !< grid values
      real(kind = wp),dimension(:,:,:), allocatable :: M !< material properties
      real(kind = wp),dimension(:,:,:), allocatable :: n_l , &
           n_m , &
           n_n !< unit normal vector
      real(kind = wp),dimension(:,:,:), allocatable :: F !< fields on boundary (copied from block_fields)
      real(kind = wp),dimension(:,:,:), allocatable :: DF !< rates on boundary (copied to block_rates)
      real(kind = wp),dimension(:,:,:), allocatable :: Fopp !< field values from other side
      real(kind = wp),dimension(:,:,:), allocatable :: Mopp !< material values from other side
      real(kind = wp),dimension(:,:,:), allocatable :: U , &
           DU !< displacements (and associated rate array)
   end type block_boundary
   
   !> boundary_type datatype created by Kenneth. Its functionality should be
   !> folded into the block_boundary above
   type :: boundary_type
      integer :: block_num
      integer :: Rx, Ry, Rz, Lx, Ly, Lz
   end type boundary_type
 
   type :: moment_tensor
       logical :: use_moment_tensor
       character(64), dimension(:), allocatable :: source_type
       real(kind = wp), dimension(:), allocatable :: mXX,mXY,mXZ,mYY,mYZ,mZZ
       real(kind = wp), dimension(:), allocatable :: location_x,location_y, location_z
       integer, dimension(:), allocatable :: near_x, near_y, near_z, alpha
       real(kind = wp), dimension(:), allocatable :: near_phys_x, near_phys_y, near_phys_z
       real(kind = wp), dimension(:,:,:,:), allocatable :: exact
       real(kind = wp), dimension(:), allocatable :: duration,t_init
       integer :: order,block_id,num_tensor
   end type moment_tensor  
 
   !> block_pml datatype holds the different components of the pml auxiliary variables and rates.
   type :: block_pml
      logical :: pml !< logical variable (true, false) denoting if PML needs to implemented
      integer :: N_pml !< number of grid points inside the PML
      type(cartesian3d_t) :: C
      real(kind = wp),dimension(:,:,:,:), allocatable :: Q !< PML auxiliary fields 
      real(kind = wp),dimension(:,:,:,:), allocatable :: DQ !< PML auxiliary rates 
   end type block_pml
   
 
   !> block datatype holds all the different pieces above.
   type :: block_type
      integer :: id !< ID of the block, is equal to 2^N where N is the block number
      character(64) :: physics !< physics (e.g., elastic or acoustic medium)
      character(64) :: fd_type
      integer :: order, nb
      type(block_indices) :: I !< indices and such for defining fields on a structured mesh
      type(block_grid_t) :: G !< spatial grid (physical coordinates and metric derivatives of mapping)
      type(block_material) :: M !< material properties
      type(block_fields) :: F !< fields (on mesh within block)
      type(moment_tensor) :: MT   !<moment tensor sources within block
      type(block_plastic) :: P !< plasticity variables 
      type(block_boundary) :: B(6) !< six sides on a block in 3D, with various arrays defined on them
      type(boundary_type) :: boundary_vars
      type(block_pml) :: PMLB(6) !< six sides on a block in 3D, with various PML auxiliary arrays defined on them
      real(kind = wp) :: tau0
      real(kind = wp) :: rho_s_p(3)
      real(kind = wp) :: sum
      logical :: debug = .false.
   end type block_type
 
   !> iface datatype, which includes various fields defined on the interface as well as
   !> information like indices of blocks on either side of interface
 
   type :: iface_type
      integer :: id !< ID of the interface, is equal to the sum of the two block IDs
      integer :: im,ip !< indices of minus and plus blocks
      character(1) :: direction !< direction of interface (q,r,s)
      real(kind = wp),dimension(:,:,:), allocatable :: V , &
                                       DV !< discontinuity in velocity vector across interface
      real(kind = wp),dimension(:,:,:), allocatable :: T !< traction vector on interface
      real(kind = wp),dimension(:,:,:), allocatable :: S , &
                                       DS !< displacement discontinuity across interface and its rates array
      real(kind = wp),dimension(:,:,:), allocatable :: W , &
                                       DW !< state variable and its rate
      real(kind = wp),dimension(:,:,:), allocatable :: Svel !< slip velocity across the interface
      real(kind = wp),dimension(:,:,:), allocatable :: trup !< slip velocity across the interface
      type(interface3d) :: II !< Interface communicator
   end type iface_type
 
   !> fault datatype for outputting the different values at the interface
   type fault_type
      type(file_distributed), dimension(8) :: handles
      logical, dimension(8) :: output_enabled = .true.
      character(len=256) :: output_directory = '.'
      integer :: array_s
      real(kind = wp),  allocatable, dimension(:, :, :) :: time_rup
      real(kind = wp),  allocatable, dimension(:, :, :) :: W
      real(kind = wp),  allocatable, dimension(:, :, :) :: slip , &
                                           Svel
      real(kind = wp),  allocatable, dimension(:, :, :) :: U_pluspres , &
                                           V_pluspres , &
                                           Uhat_pluspres , &
                                           Vhat_pluspres
   end type fault_type

   type :: fault_receiver_output_type
      logical :: enabled = .false.
      integer :: nreceivers = 0
      integer :: stride = 1
      integer :: step_counter = 0
      integer :: comm = -1
      integer :: comm_rank = -1
      character(len=256) :: directory = 'fault_receivers'
      character(len=64), allocatable :: name(:)
      real(kind=wp), allocatable :: target_xyz(:,:), actual_xyz(:,:)
      real(kind=wp), allocatable :: distance(:)
      integer, allocatable :: owner(:), j(:), k(:), file_unit(:)
   end type fault_receiver_output_type
 
   !> mms datatype
   type :: mms_type
      logical :: use_mms
      real(kind = wp) :: nx, ny, nz, nt ! spatial wavenumbers and frequency
   end type mms_type
 
   !> seismogram datatype
   type :: seismogram_type
      logical :: output_exact_moment
      logical :: output_seismograms,output_fields_block1,output_fields_block2
        logical :: output_station_info
      logical :: station_xyz_index
      logical :: station_number_in_list, station_number_in_filename
      logical :: station_use_block_subdirectories, station_both_blocks
      integer :: stride_fields,file_unit_block1(9),file_unit_block2(9)
      integer :: nstations,block_num
      integer :: station_output_order(4) = [1, 2, 3, 4]
      integer,dimension(:),allocatable :: i,j,k,file_unit,station_number
      real(kind = wp), dimension(:), allocatable :: i_phys,j_phys,k_phys
   end type seismogram_type
 
   !> slice datatype
   type slice_type
      type(file_distributed) :: Hslice,Vslice
      logical :: horz_slice,vert_slice
      integer :: locationH,locationV
      integer :: array_s
   end type slice_type

   !> plane output datatype (per block)
   type :: plane_output_plane
      logical :: enabled = .false.
      logical :: active = .false.
      logical :: use_logical_coord = .false.
      character(len=2) :: axis = ''
      real(kind = wp) :: coord = 0.0_wp
      character(len=64) :: name = ''
      character(len=256) :: file_directory = '.'
      integer :: stride = 1
      integer :: step_counter = 0

      integer :: fixed_dir = 0
      integer :: fixed_index_global = 0

      integer :: plane_comm = -1
      integer :: plane_rank = -1
      integer :: plane_size = 0

      integer :: n1 = 0
      integer :: n2 = 0

      integer :: file_unit = -1
      logical :: header_written = .false.

      real(kind = wp), allocatable :: vx(:,:), vy(:,:), vz(:,:)
   end type plane_output_plane

   type :: plane_output_type
      logical :: enabled = .false.
      integer :: nplanes = 0
      type(plane_output_plane), allocatable :: P(:)
   end type plane_output_type
 
 
   !> domain datatype contains the different blocks and interfaces.
   !> For now, the output datatypes are also held here but that
   !> should be cleaned out.
 
   type :: domain_type
      integer :: nblocks,nifaces !< number of blocks and interfaces
      integer :: w_stride !< output stride
      real(kind = wp) :: dt, t, t_final !< time
      real(kind = wp) :: CFL !< CFL
      real(kind = wp) :: topo !< topography prefactor
      integer :: nt !< number of time steps
      type(block_type),allocatable :: B(:) !< blocks
      type(iface_type),allocatable :: I(:) !< interfaces
      type(fault_type) :: fault
      type(fault_receiver_output_type) :: fault_receivers
      type(mms_type) :: mms_vars
      character(256) :: name, problem,response,plastic_model
      character(64) :: coupling, fd_type
      character(64) :: type_of_mesh, mesh_source,material_source
      logical :: w_fault, output_fault_topo, interpol,use_topography,mollify_source
      integer :: order !< order of accuracy of fd operator
      type(seismogram_type), allocatable :: seismometers(:)   !< seismograms
      type(slice_type) :: slicer
      type(plane_output_type), allocatable :: plane_outputs(:)
      logical :: debug !< enable diagnostic output
   end type domain_type
 
   type:: block_temp_parameters
      ! Defaults make omitted/inactive namelist entries deterministic.  Active
      ! blocks are validated in init_domain before any of these values are used.
      integer :: nqrs(3) = 0
      real(kind = wp) :: aqrs(3) = 0.0_wp, bqrs(3) = 0.0_wp
      real(kind = wp) :: lc = 0.0_wp, rc = 0.0_wp
      real(kind = wp) :: rho_s_p(3) = 0.0_wp
      real(kind = wp) :: mu_beta_eta(3) = 0.0_wp
      integer :: lqrs(3) = 0, rqrs(3) = 0
      character(len=256) :: profile_type = '', profile_path = '', material_path(3) = ''
      character(len=256) :: topography_type = '', topography_path = ''
      logical :: pml_lqrs(3) = .false., pml_rqrs(3) = .false.
      integer :: npml = 0, faultsize(2) = 0
   end type block_temp_parameters
 
 end module datatypes
 
