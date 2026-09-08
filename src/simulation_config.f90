module simulation_config

  use common, only : wp
  use datatypes, only : block_temp_parameters
  use anelastic_q4_model, only : q4_parameters
  use anelastic_q8_model, only : q8_parameters
  use anelastic_cq8_b2_model, only : cq8_b2_parameters
  use anelastic_cq_model, only : cq_parameters
  use anelastic_fq8_model, only : fq8_parameters
  implicit none
  private

  type, public :: problem_config_t
     character(len=256) :: name = 'default'
     character(len=256) :: problem = 'TPV5'
     character(len=256) :: response = 'elastic'
     character(len=256) :: plastic_model = 'default'
     character(len=64) :: coupling = 'locked'
     character(len=64) :: fd_type = 'traditional'
     character(len=64) :: mesh_source = 'compute'
     character(len=64) :: type_of_mesh = 'cartesian'
     character(len=64) :: material_source = 'hardcode'
     integer :: nblocks = 2
     integer :: nt = 0
     integer :: order = 5
     integer :: w_stride = 1
     real(wp) :: CFL = 0.5_wp
     real(wp) :: t_final = 0.0_wp
     real(wp) :: topo = 1.0_wp
     logical :: w_fault = .true.
     logical :: interpol = .false.
     logical :: use_topography = .false.
     logical :: mollify_source = .false.
  end type problem_config_t

  type, public :: simulation_config_t
     type(problem_config_t) :: problem
     type(block_temp_parameters), allocatable :: blocks(:)
     type(q4_parameters) :: q4
     type(q8_parameters) :: q8
     type(cq8_b2_parameters) :: cq8_b2
     type(cq_parameters) :: cq
     type(fq8_parameters) :: fq8
     logical :: has_q4 = .false.
     logical :: has_q8 = .false.
     logical :: has_cq8_b2 = .false.
     logical :: has_cq = .false.
     logical :: has_fq8 = .false.
     integer :: process_dims(2,3) = 1
     integer :: block_sizes(2) = 0
     integer :: rank_begin(2) = 0
     integer :: rank_end(2) = -1
     logical :: serial_shared_blocks = .false.
  end type simulation_config_t

end module simulation_config
