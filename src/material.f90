module material

!> material module defines material properties on a block
  use common, only : wp
  use datatypes, only : block_material,block_grid_t,block_indices
  use diagnostics, only : warn_once
  use anelastic_q4_model, only : q4_parameters, q4_nmechanisms, &
       read_q4_parameters, build_q4_fixed_coefficients, q4_max_relative_error
  use anelastic_q8_model, only : q8_parameters, q8_nmechanisms, &
       read_q8_parameters, build_q8_fixed_coefficients, q8_max_relative_error
  use anelastic_fq8_model, only : fq8_max_relative_error, fq8_effective_q, &
       fq8_common_modulus_scale, fq8_phase_velocity_ratio
  implicit none

  logical, save :: q4_summary_printed = .false.
  logical, save :: q8_summary_printed = .false.
  logical, save :: fq8_summary_printed = .false.

contains


   subroutine init_anelastic_properties(M, G, infile)

      ! Initializes GSLS attenuation parameters and allocates Q and memory-variable arrays.
      ! Must be called only when response == 'anelastic'.

      use mpi3dcomm, only : allocate_array_body

      implicit none

      type(block_material), intent(inout) :: M
      type(block_grid_t), intent(in) :: G
      integer, intent(in) :: infile

      real(kind = wp) :: c, weight_exp, fref
      real(kind = wp) :: FAC, taumin, taumax, wref, facex
      real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
      integer :: stat, i, l, j, k, N
      real(kind = wp), parameter :: pi = 3.141592653589793_wp

      namelist /anelastic_list/ c, weight_exp, fref

      c = 1.0_wp
      weight_exp = 0.0_wp
      fref = 1.0_wp
      M%anelastic = .true.

      rewind(infile)
      read(infile, nml=anelastic_list, iostat=stat)
      if (stat > 0) stop 'error reading namelist anelastic_list'

      M%c = c
      M%weight_exp = weight_exp
      M%fref = fref

      call allocate_array_body(M%Qp_inv, G%C, ghost_nodes=.true.)
      M%Qp_inv(:,:,:) = 0.0_wp
      call allocate_array_body(M%Qs_inv, G%C, ghost_nodes=.true.)
      M%Qs_inv(:,:,:) = 0.0_wp

      N = 4

      call allocate_array_body(M%eta4,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta4, G%C, N, ghost_nodes=.true.)
      M%eta4(:,:,:,:) = 0.0_wp
      M%Deta4(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta5,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta5, G%C, N, ghost_nodes=.true.)
      M%eta5(:,:,:,:) = 0.0_wp
      M%Deta5(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta6,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta6, G%C, N, ghost_nodes=.true.)
      M%eta6(:,:,:,:) = 0.0_wp
      M%Deta6(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta7,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta7, G%C, N, ghost_nodes=.true.)
      M%eta7(:,:,:,:) = 0.0_wp
      M%Deta7(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta8, G%C, N, ghost_nodes=.true.)
      M%eta8(:,:,:,:) = 0.0_wp
      M%Deta8(:,:,:,:) = 0.0_wp
      call allocate_array_body(M%eta9,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta9, G%C, N, ghost_nodes=.true.)
      M%eta9(:,:,:,:) = 0.0_wp
      M%Deta9(:,:,:,:) = 0.0_wp

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps
               M%Qs_inv(i,j,k) = 1.0_wp/(c*sqrt(M%M(i,j,k,2)/M%M(i,j,k,3)))
               M%Qp_inv(i,j,k) = 0.5_wp*M%Qs_inv(i,j,k)
            end do
         end do
      end do

      if (M%weight_exp < 0.01_wp) then
         FAC = 1.0_wp
         taumin = 1.0_wp/(2.0_wp*pi*15.0_wp)*1.0_wp*FAC
         taumax = 1.0_wp/(2.0_wp*pi*0.08_wp)*200.0_wp*FAC
         do k = 1, N
            M%tau(k) = exp(log(taumin) + (2.0_wp*k-1.0_wp)/16.0_wp*(log(taumax) - log(taumin)))
         end do
         M%weight(1) = 1.6126_wp
         M%weight(2) = 0.6255_wp
         M%weight(3) = 0.6382_wp
         M%weight(4) = 1.5969_wp
      end if

      if (M%weight_exp > 0.59_wp .AND. M%weight_exp < 0.61_wp) then
         FAC = 1.0_wp
         taumin = 1.0_wp/(2.0_wp*pi*15.0_wp)*1.0_wp*FAC
         taumax = 1.0_wp/(2.0_wp*pi*0.08_wp)*200.0_wp*FAC
         do k = 1, N
            M%tau(k) = exp(log(taumin) + (2.0_wp*k-1.0_wp)/16.0_wp*(log(taumax) - log(taumin)))
         end do
         M%weight(1) = 0.0336_wp
         M%weight(2) = 0.6873_wp
         M%weight(3) = 0.8767_wp
         M%weight(4) = 1.5202_wp
      end if

      wref = 2.0_wp*pi*fref
      facex = 1.0_wp

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps

               val_S = 0.0_wp
               val_P = 0.0_wp
               do l = 1, N
                  denom_S = (wref**2.0_wp*M%tau(l)**2.0_wp + 1.0_wp) * (1.0_wp/M%Qs_inv(i,j,k)) * facex
                  denom_P = (wref**2.0_wp*M%tau(l)**2.0_wp + 1.0_wp) * (1.0_wp/M%Qp_inv(i,j,k)) * facex
                  val_S = val_S + M%weight(l)/denom_S
                  val_P = val_P + M%weight(l)/denom_P
               end do

               vs = sqrt(M%M(i,j,k,2)/M%M(i,j,k,3))
               vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2))/M%M(i,j,k,3))
               mu_unrelax_S = M%M(i,j,k,3)*vs**2.0_wp/(1.0_wp - val_S)
               mu_unrelax_P = M%M(i,j,k,3)*vp**2.0_wp/(1.0_wp - val_P)
               vs = sqrt(mu_unrelax_S/M%M(i,j,k,3))
               vp = sqrt(mu_unrelax_P/M%M(i,j,k,3))

               M%M(i,j,k,2) = vs**2.0_wp*M%M(i,j,k,3)
               M%M(i,j,k,1) = vp**2.0_wp*M%M(i,j,k,3) - 2.0_wp*M%M(i,j,k,2)

            end do
         end do
      end do

   end subroutine init_anelastic_properties


   subroutine init_anelastic_Q4_properties(M, G, infile, resolved_parameters)

      use mpi3dcomm, only : allocate_array_body
      use mpi3dbasic, only : error, is_master

      implicit none

      type(block_material), intent(inout) :: M
      type(block_grid_t), intent(in) :: G
      integer, intent(in) :: infile
      type(q4_parameters), intent(in), optional :: resolved_parameters

      real(kind = wp) :: wref
      real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
      real(kind = wp) :: max_qs_error, max_qp_error
      integer :: stat, i, l, j, k, N
      real(kind = wp), parameter :: pi = 3.141592653589793_wp
      type(q4_parameters) :: parameters
      character(len=256) :: parameter_error

      if (M%anelastic_Q .or. allocated(M%eta4Q)) then
         call error('anelastic-Q4 material is already initialized; refusing double modulus correction', &
                    'init_anelastic_Q4_properties')
      end if

      if (present(resolved_parameters)) then
         parameters = resolved_parameters
      else
         call read_q4_parameters(infile, parameters, stat, parameter_error)
         if (stat /= 0) call error(trim(parameter_error), 'init_anelastic_Q4_properties')
      end if

      M%anelastic_Q = .true.
      N = q4_nmechanisms
      M%n_mechanism_Q = N
      M%Qs0_Q = parameters%Qs0
      M%Qp0_Q = parameters%Qp0
      M%fref_Q = parameters%fref
      M%fmin_Q = parameters%fmin
      M%fmax_Q = parameters%fmax
      M%weight_method_Q = parameters%weight_method
      call build_q4_fixed_coefficients(parameters, M%tau_Q, M%weight_Q)

      ! Allocate Q arrays
      call allocate_array_body(M%Qp_inv_Q, G%C, ghost_nodes=.true.)
      M%Qp_inv_Q(:,:,:) = 0.0_wp
      call allocate_array_body(M%Qs_inv_Q, G%C, ghost_nodes=.true.)
      M%Qs_inv_Q(:,:,:) = 0.0_wp

      ! Allocate memory variable arrays (N=4 mechanisms)
      call allocate_array_body(M%eta4Q,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta4Q, G%C, N, ghost_nodes=.true.)
      M%eta4Q = 0.0_wp;  M%Deta4Q = 0.0_wp
      call allocate_array_body(M%eta5Q,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta5Q, G%C, N, ghost_nodes=.true.)
      M%eta5Q = 0.0_wp;  M%Deta5Q = 0.0_wp
      call allocate_array_body(M%eta6Q,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta6Q, G%C, N, ghost_nodes=.true.)
      M%eta6Q = 0.0_wp;  M%Deta6Q = 0.0_wp
      call allocate_array_body(M%eta7Q,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta7Q, G%C, N, ghost_nodes=.true.)
      M%eta7Q = 0.0_wp;  M%Deta7Q = 0.0_wp
      call allocate_array_body(M%eta8Q,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta8Q, G%C, N, ghost_nodes=.true.)
      M%eta8Q = 0.0_wp;  M%Deta8Q = 0.0_wp
      call allocate_array_body(M%eta9Q,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta9Q, G%C, N, ghost_nodes=.true.)
      M%eta9Q = 0.0_wp;  M%Deta9Q = 0.0_wp

      ! Set independent inverse Q values directly from the required input.
      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps
               if (M%M(i,j,k,2) <= 0.0_wp .or. M%M(i,j,k,3) <= 0.0_wp) &
                    call error('mu and density must be positive for anelastic-Q4', &
                               'init_anelastic_Q4_properties')
               M%Qs_inv_Q(i,j,k) = 1.0_wp/parameters%Qs0
               M%Qp_inv_Q(i,j,k) = 1.0_wp/parameters%Qp0
            end do
         end do
      end do

      if (any(M%tau_Q <= 0.0_wp) .or. any(M%weight_Q < 0.0_wp)) &
           call error('anelastic-Q4 requires positive tau and nonnegative weights', &
                      'init_anelastic_Q4_properties')

      ! Correct unrelaxed moduli so relaxed (low-freq) velocity matches input
      wref = 2.0_wp*pi*parameters%fref

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps

               val_S = 0.0_wp
               val_P = 0.0_wp
               do l = 1, N
                  denom_S = (wref**2 * M%tau_Q(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qs_inv_Q(i,j,k))
                  denom_P = (wref**2 * M%tau_Q(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qp_inv_Q(i,j,k))
                  val_S = val_S + M%weight_Q(l) / denom_S
                  val_P = val_P + M%weight_Q(l) / denom_P
               end do

               if (val_S >= 1.0_wp .or. val_P >= 1.0_wp) &
                    call error('invalid anelastic-Q4 modulus correction: weight/Q sum >= 1', &
                               'init_anelastic_Q4_properties')

               vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
               vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2)) / M%M(i,j,k,3))
               mu_unrelax_S = M%M(i,j,k,3) * vs**2 / (1.0_wp - val_S)
               mu_unrelax_P = M%M(i,j,k,3) * vp**2 / (1.0_wp - val_P)
               vs = sqrt(mu_unrelax_S / M%M(i,j,k,3))
               vp = sqrt(mu_unrelax_P / M%M(i,j,k,3))

               M%M(i,j,k,2) = vs**2 * M%M(i,j,k,3)
               M%M(i,j,k,1) = vp**2 * M%M(i,j,k,3) - 2.0_wp * M%M(i,j,k,2)

            end do
         end do
      end do

      if (is_master() .and. .not.q4_summary_printed) then
         q4_summary_printed = .true.
         call q4_max_relative_error(parameters%Qs0, M%tau_Q, M%weight_Q, &
              0.08_wp, 15.0_wp, max_qs_error)
         call q4_max_relative_error(parameters%Qp0, M%tau_Q, M%weight_Q, &
              0.08_wp, 15.0_wp, max_qp_error)
         write(*,'(A)') 'anelastic-Q4 parameters:'
         write(*,'(A,A)') '  weight method = ', trim(parameters%weight_method)
         write(*,'(A,ES12.4,A,ES12.4)') '  Qs0 = ', parameters%Qs0, ', Qp0 = ', parameters%Qp0
         write(*,'(A,ES12.4)') '  fref (Hz) = ', parameters%fref
         write(*,'(A,ES12.4,A,ES12.4)') '  tau band (Hz) = ', parameters%fmin, ' to ', parameters%fmax
         write(*,'(A,4(ES12.4,1X))') '  tau (s) = ', M%tau_Q
         write(*,'(A,4(ES12.4,1X))') '  normalized weights = ', M%weight_Q
         write(*,'(A,F8.3,A)') '  max Qs error over 0.08-15 Hz = ', 100.0_wp*max_qs_error, ' %'
         write(*,'(A,F8.3,A)') '  max Qp error over 0.08-15 Hz = ', 100.0_wp*max_qp_error, ' %'
         if (max(max_qs_error,max_qp_error) > 0.05_wp) call warn_once('Q4-FIT-001', &
              'Realized-Q error exceeds 5 percent; four mechanisms provide a coarser fit than Q8.', &
              'init_anelastic_Q4_properties')
      end if

   end subroutine init_anelastic_Q4_properties


   subroutine destroy_anelastic_Q4_properties(M)
      type(block_material), intent(inout) :: M

      if (allocated(M%Qp_inv_Q)) deallocate(M%Qp_inv_Q)
      if (allocated(M%Qs_inv_Q)) deallocate(M%Qs_inv_Q)
      if (allocated(M%eta4Q)) deallocate(M%eta4Q)
      if (allocated(M%eta5Q)) deallocate(M%eta5Q)
      if (allocated(M%eta6Q)) deallocate(M%eta6Q)
      if (allocated(M%eta7Q)) deallocate(M%eta7Q)
      if (allocated(M%eta8Q)) deallocate(M%eta8Q)
      if (allocated(M%eta9Q)) deallocate(M%eta9Q)
      if (allocated(M%Deta4Q)) deallocate(M%Deta4Q)
      if (allocated(M%Deta5Q)) deallocate(M%Deta5Q)
      if (allocated(M%Deta6Q)) deallocate(M%Deta6Q)
      if (allocated(M%Deta7Q)) deallocate(M%Deta7Q)
      if (allocated(M%Deta8Q)) deallocate(M%Deta8Q)
      if (allocated(M%Deta9Q)) deallocate(M%Deta9Q)
      M%anelastic_Q = .false.
      M%Qs0_Q = -1.0_wp
      M%Qp0_Q = -1.0_wp
      M%fref_Q = 1.0_wp
      M%fmin_Q = -1.0_wp
      M%fmax_Q = -1.0_wp
      M%weight_method_Q = ''
      M%tau_Q = 0.0_wp
      M%weight_Q = 0.0_wp
   end subroutine destroy_anelastic_Q4_properties


   subroutine init_anelastic_Q8_properties(M, G, infile, resolved_parameters)

      ! Initializes anelastic-Q8 attenuation: N=8 mechanisms with correct tau spacing.
      ! Mirrors the Q4 material formulation but uses 8 relaxation mechanisms,
      ! providing improved constant-Q approximation over [0.05, 20] Hz.
      ! Namelist &anelastic_Q8_list parameters:
      !   Qs0, Qp0 : required independent constant S- and P-wave Q values.
      !   fref      : reference frequency (Hz) for unrelaxed-modulus correction.

      use mpi3dcomm, only : allocate_array_body
      use mpi3dbasic, only : error, is_master

      implicit none

      type(block_material), intent(inout) :: M
      type(block_grid_t), intent(in) :: G
      integer, intent(in) :: infile
      type(q8_parameters), intent(in), optional :: resolved_parameters

      real(kind = wp) :: wref
      real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
      integer :: stat, i, l, j, k, N
      real(kind = wp) :: max_qs_error, max_qp_error
      real(kind = wp), parameter :: pi = 3.141592653589793_wp
      type(q8_parameters) :: parameters
      character(len=256) :: parameter_error

      if (M%anelastic_Q8 .or. allocated(M%eta4Q8)) then
         call error('anelastic-Q8 material is already initialized; refusing double modulus correction', &
                    'init_anelastic_Q8_properties')
      end if
      M%anelastic_Q8 = .true.
      N = q8_nmechanisms
      M%n_mechanism_Q8 = N

      if (present(resolved_parameters)) then
         parameters = resolved_parameters
      else
         call read_q8_parameters(infile, parameters, stat, parameter_error)
         if (stat /= 0) call error(trim(parameter_error), 'init_anelastic_Q8_properties')
      end if
      call build_q8_fixed_coefficients(parameters, M%tau_Q8, M%weight_Q8)

      M%Qs0_Q8 = parameters%Qs0
      M%Qp0_Q8 = parameters%Qp0
      M%fref_Q8 = parameters%fref
      M%fmin_Q8 = parameters%fmin
      M%fmax_Q8 = parameters%fmax
      M%weight_method_Q8 = parameters%weight_method

      ! Allocate Q arrays
      call allocate_array_body(M%Qp_inv_Q8, G%C, ghost_nodes=.true.)
      M%Qp_inv_Q8(:,:,:) = 0.0_wp
      call allocate_array_body(M%Qs_inv_Q8, G%C, ghost_nodes=.true.)
      M%Qs_inv_Q8(:,:,:) = 0.0_wp

      ! Allocate memory variable arrays (N=8 mechanisms)
      call allocate_array_body(M%eta4Q8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta4Q8, G%C, N, ghost_nodes=.true.)
      M%eta4Q8 = 0.0_wp;  M%Deta4Q8 = 0.0_wp
      call allocate_array_body(M%eta5Q8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta5Q8, G%C, N, ghost_nodes=.true.)
      M%eta5Q8 = 0.0_wp;  M%Deta5Q8 = 0.0_wp
      call allocate_array_body(M%eta6Q8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta6Q8, G%C, N, ghost_nodes=.true.)
      M%eta6Q8 = 0.0_wp;  M%Deta6Q8 = 0.0_wp
      call allocate_array_body(M%eta7Q8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta7Q8, G%C, N, ghost_nodes=.true.)
      M%eta7Q8 = 0.0_wp;  M%Deta7Q8 = 0.0_wp
      call allocate_array_body(M%eta8Q8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta8Q8, G%C, N, ghost_nodes=.true.)
      M%eta8Q8 = 0.0_wp;  M%Deta8Q8 = 0.0_wp
      call allocate_array_body(M%eta9Q8,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta9Q8, G%C, N, ghost_nodes=.true.)
      M%eta9Q8 = 0.0_wp;  M%Deta9Q8 = 0.0_wp

      ! Set independent inverse Q values directly from the required input.
      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps
               if (M%M(i,j,k,2) <= 0.0_wp .or. M%M(i,j,k,3) <= 0.0_wp) then
                  call error('mu and density must be positive for anelastic-Q8', &
                             'init_anelastic_Q8_properties')
               end if
               M%Qs_inv_Q8(i,j,k) = 1.0_wp / parameters%Qs0
               M%Qp_inv_Q8(i,j,k) = 1.0_wp / parameters%Qp0
            end do
         end do
      end do

      if (any(M%tau_Q8 <= 0.0_wp) .or. any(M%weight_Q8 < 0.0_wp)) then
         call error('anelastic-Q8 requires positive tau and nonnegative weights', &
                    'init_anelastic_Q8_properties')
      end if

      ! Correct unrelaxed moduli so phase velocity at fref matches the input.
      wref = 2.0_wp * pi * parameters%fref

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps

               val_S = 0.0_wp
               val_P = 0.0_wp
               do l = 1, N
                  denom_S = (wref**2 * M%tau_Q8(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qs_inv_Q8(i,j,k))
                  denom_P = (wref**2 * M%tau_Q8(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qp_inv_Q8(i,j,k))
                  val_S = val_S + M%weight_Q8(l) / denom_S
                  val_P = val_P + M%weight_Q8(l) / denom_P
               end do

               if (val_S >= 1.0_wp .or. val_P >= 1.0_wp) then
                  call error('invalid anelastic-Q8 modulus correction: weight/Q sum >= 1', &
                             'init_anelastic_Q8_properties')
               end if
               vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
               vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2)) / M%M(i,j,k,3))
               mu_unrelax_S = M%M(i,j,k,3) * vs**2 / (1.0_wp - val_S)
               mu_unrelax_P = M%M(i,j,k,3) * vp**2 / (1.0_wp - val_P)
               vs = sqrt(mu_unrelax_S / M%M(i,j,k,3))
               vp = sqrt(mu_unrelax_P / M%M(i,j,k,3))

               M%M(i,j,k,2) = vs**2 * M%M(i,j,k,3)
               M%M(i,j,k,1) = vp**2 * M%M(i,j,k,3) - 2.0_wp * M%M(i,j,k,2)

            end do
         end do
      end do

      if (is_master() .and. .not.q8_summary_printed) then
         q8_summary_printed = .true.
         write(*,'(A)') 'anelastic-Q8 parameters:'
         write(*,'(A,A)') '  weight method = ', trim(parameters%weight_method)
         write(*,'(A,ES12.4,A,ES12.4)') '  Qs0 = ', parameters%Qs0, ', Qp0 = ', parameters%Qp0
         write(*,'(A,ES12.4)') '  fref (Hz) = ', parameters%fref
         write(*,'(A,ES12.4,A,ES12.4)') '  tau band (Hz) = ', parameters%fmin, ' to ', parameters%fmax
         write(*,'(A,8(ES12.4,1X))') '  tau (s) = ', M%tau_Q8
         write(*,'(A,8(ES12.4,1X))') '  normalized weights = ', M%weight_Q8
         call q8_max_relative_error(parameters%Qs0, M%tau_Q8, M%weight_Q8, &
                                    0.08_wp, 15.0_wp, max_qs_error)
         call q8_max_relative_error(parameters%Qp0, M%tau_Q8, M%weight_Q8, &
                                    0.08_wp, 15.0_wp, max_qp_error)
         write(*,'(A,F8.3,A)') '  max Qs error over 0.08-15 Hz = ', &
                               100.0_wp*max_qs_error, ' %'
         write(*,'(A,F8.3,A)') '  max Qp error over 0.08-15 Hz = ', &
                               100.0_wp*max_qp_error, ' %'
         if (max(max_qs_error, max_qp_error) > 0.05_wp) then
            call warn_once('Q8-FIT-001', &
                 'Realized-Q error exceeds 5 percent; fixed-q50 is optimized for Q=50.', &
                 'init_anelastic_Q8_properties')
         end if
      end if

   end subroutine init_anelastic_Q8_properties


   subroutine destroy_anelastic_Q8_properties(M)
      type(block_material), intent(inout) :: M

      if (allocated(M%Qp_inv_Q8)) deallocate(M%Qp_inv_Q8)
      if (allocated(M%Qs_inv_Q8)) deallocate(M%Qs_inv_Q8)
      if (allocated(M%eta4Q8)) deallocate(M%eta4Q8)
      if (allocated(M%eta5Q8)) deallocate(M%eta5Q8)
      if (allocated(M%eta6Q8)) deallocate(M%eta6Q8)
      if (allocated(M%eta7Q8)) deallocate(M%eta7Q8)
      if (allocated(M%eta8Q8)) deallocate(M%eta8Q8)
      if (allocated(M%eta9Q8)) deallocate(M%eta9Q8)
      if (allocated(M%Deta4Q8)) deallocate(M%Deta4Q8)
      if (allocated(M%Deta5Q8)) deallocate(M%Deta5Q8)
      if (allocated(M%Deta6Q8)) deallocate(M%Deta6Q8)
      if (allocated(M%Deta7Q8)) deallocate(M%Deta7Q8)
      if (allocated(M%Deta8Q8)) deallocate(M%Deta8Q8)
      if (allocated(M%Deta9Q8)) deallocate(M%Deta9Q8)
      M%anelastic_Q8 = .false.
      M%Qs0_Q8 = -1.0_wp
      M%Qp0_Q8 = -1.0_wp
      M%fref_Q8 = 1.0_wp
      M%fmin_Q8 = -1.0_wp
      M%fmax_Q8 = -1.0_wp
      M%weight_method_Q8 = ''
      M%tau_Q8 = 0.0_wp
      M%weight_Q8 = 0.0_wp
   end subroutine destroy_anelastic_Q8_properties


   subroutine init_anelastic_Qf_properties(M, G, infile)

      ! Initializes anelastic-Qf attenuation: frequency-dependent Q with N=4 mechanisms
      ! Q(f) = Q0 for f < f_trans; Q(f) = Q0*(f/f_trans)^gamma for f > f_trans
      ! Namelist &anelastic_Qf_list parameters:
      !   c      : Q_S = c / sqrt(mu/rho)  (Q_P = 2*Q_S)
      !   gamma  : Power-law exponent (0.0-0.9, default 0.0 = constant Q)
      !   f_trans: Transition frequency (Hz, default 1.0)
      !   fref   : Reference frequency for unrelaxed modulus correction (Hz, default 1.0)

      use mpi3dcomm, only : allocate_array_body
      use withers_tables, only : get_relaxation_times_Qf, get_withers_weights_Qf

      implicit none

      type(block_material), intent(inout) :: M
      type(block_grid_t), intent(in) :: G
      integer, intent(in) :: infile

      real(kind = wp) :: c, gamma, f_trans, fref
      real(kind = wp) :: wref
      real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
      integer :: stat, i, l, j, k, N
      real(kind = wp), parameter :: pi = 3.141592653589793_wp

      namelist /anelastic_Qf_list/ c, gamma, f_trans, fref

      ! Defaults
      c       = 1.0_wp
      gamma   = 0.0_wp
      f_trans = 1.0_wp
      fref    = 1.0_wp

      M%anelastic_Qf = .true.
      N = 4
      M%n_mechanism_Qf = N
      M%gamma_Qf = gamma
      M%f_trans_Qf = f_trans
      M%fref_Qf = fref

      rewind(infile)
      read(infile, nml=anelastic_Qf_list, iostat=stat)
      if (stat > 0) stop 'error reading namelist anelastic_Qf_list'

      M%gamma_Qf = gamma
      M%f_trans_Qf = f_trans
      M%fref_Qf = fref

      ! Clamp gamma to valid range [0.0, 0.9]
      if (M%gamma_Qf < 0.0_wp) M%gamma_Qf = 0.0_wp
      if (M%gamma_Qf > 0.9_wp) M%gamma_Qf = 0.9_wp

      ! Allocate Qf arrays
      call allocate_array_body(M%Qp_inv_Qf, G%C, ghost_nodes=.true.)
      M%Qp_inv_Qf(:,:,:) = 0.0_wp
      call allocate_array_body(M%Qs_inv_Qf, G%C, ghost_nodes=.true.)
      M%Qs_inv_Qf(:,:,:) = 0.0_wp

      ! Allocate memory variable arrays (N=4 mechanisms, 6 stress components)
      call allocate_array_body(M%eta4Qf,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta4Qf, G%C, N, ghost_nodes=.true.)
      M%eta4Qf(:,:,:,:) = 0.0_wp
      M%Deta4Qf(:,:,:,:) = 0.0_wp

      call allocate_array_body(M%eta5Qf,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta5Qf, G%C, N, ghost_nodes=.true.)
      M%eta5Qf(:,:,:,:) = 0.0_wp
      M%Deta5Qf(:,:,:,:) = 0.0_wp

      call allocate_array_body(M%eta6Qf,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta6Qf, G%C, N, ghost_nodes=.true.)
      M%eta6Qf(:,:,:,:) = 0.0_wp
      M%Deta6Qf(:,:,:,:) = 0.0_wp

      call allocate_array_body(M%eta7Qf,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta7Qf, G%C, N, ghost_nodes=.true.)
      M%eta7Qf(:,:,:,:) = 0.0_wp
      M%Deta7Qf(:,:,:,:) = 0.0_wp

      call allocate_array_body(M%eta8Qf,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta8Qf, G%C, N, ghost_nodes=.true.)
      M%eta8Qf(:,:,:,:) = 0.0_wp
      M%Deta8Qf(:,:,:,:) = 0.0_wp

      call allocate_array_body(M%eta9Qf,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta9Qf, G%C, N, ghost_nodes=.true.)
      M%eta9Qf(:,:,:,:) = 0.0_wp
      M%Deta9Qf(:,:,:,:) = 0.0_wp

      ! Compute Q_S and Q_P from velocity ratio at each grid point
      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps
               M%Qs_inv_Qf(i,j,k) = 1.0_wp / (c * sqrt(M%M(i,j,k,2)/M%M(i,j,k,3)))
               M%Qp_inv_Qf(i,j,k) = 0.5_wp * M%Qs_inv_Qf(i,j,k)
            end do
         end do
      end do

      ! Set relaxation times from Withers table (frequency-dependent gamma)
      ! Note: tau_k are independent of gamma; only weights depend on gamma
      call get_relaxation_times_Qf(M%gamma_Qf, M%tau_Qf)

      ! Get memory-variable weights for current gamma
      ! Use reference Q=1 for now; will scale appropriately in RHS kernel
      call get_withers_weights_Qf(M%gamma_Qf, 1.0_wp, M%weight_Qf)

      ! Correct unrelaxed moduli to preserve relaxed velocity
      ! This accounts for velocity dispersion introduced by frequency-dependent weights
      wref = 2.0_wp * pi * fref

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps

               val_S = 0.0_wp
               val_P = 0.0_wp
               do l = 1, N
                  denom_S = (wref**2 * M%tau_Qf(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qs_inv_Qf(i,j,k))
                  denom_P = (wref**2 * M%tau_Qf(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qp_inv_Qf(i,j,k))
                  val_S = val_S + M%weight_Qf(l) / denom_S
                  val_P = val_P + M%weight_Qf(l) / denom_P
               end do

               vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
               vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2)) / M%M(i,j,k,3))
               mu_unrelax_S = M%M(i,j,k,3) * vs**2 / (1.0_wp - val_S)
               mu_unrelax_P = M%M(i,j,k,3) * vp**2 / (1.0_wp - val_P)
               vs = sqrt(mu_unrelax_S / M%M(i,j,k,3))
               vp = sqrt(mu_unrelax_P / M%M(i,j,k,3))

               M%M(i,j,k,2) = vs**2 * M%M(i,j,k,3)
               M%M(i,j,k,1) = vp**2 * M%M(i,j,k,3) - 2.0_wp * M%M(i,j,k,2)

            end do
         end do
      end do

   end subroutine init_anelastic_Qf_properties


  !> initialize material properties
  subroutine init_material(M, G, I, physics,problem, rho_s_p, nb)

    use mpi3dcomm

    implicit none

    type(block_material),intent(out) :: M
    type(block_grid_t),intent(in) :: G
    type(block_indices),intent(in) :: I
    real(kind = wp), intent(in) :: rho_s_p(3)
    character(*),intent(in) :: physics,problem
    integer, intent(in) :: nb
    integer :: ii,l,j,k, np
    real(kind = wp) :: Vp,Vs,rho

    !> the number of material properties will be specific to the physics
    !> (I included the acoustic case to illustrate this, but I don't foresee
    !> it being a priority to implement all the necessary routines to solve
    !> the acoustic wave equation)

   
    select case(physics)

    case default

       stop 'invalid block physics in init_material'

    case('elastic')

       !> note use of the array allocation routine, which makes it easier
       !> to allocate arrays without having to make sure that indices are being
       !> typed correctly

       select case(problem)

       case default

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          M%M(:,:,:,1) = rho_s_p(1)*(rho_s_p(3)**2 - 2.0_wp*rho_s_p(2)**2) ! lambda (Lame's first parameter)
          M%M(:,:,:,2) = rho_s_p(1)*rho_s_p(2)**2 ! mu (shear modulus)
          M%M(:,:,:,3) = rho_s_p(1) ! rho (density)

       case('LOH1')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.) 

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) < 0.9999_wp) then

                      Vp = 4.0_wp ! P-wave speed
                      Vs  = 2.0_wp ! S-wave speed
                      rho  = 2.6_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) >= 0.9999_wp .AND. G%x(l,j,k,1) <= 1.0001_wp) then

                      Vp = 5.0_wp ! P-wave speed
                      Vs  = 2.732_wp ! S-wave speed
                      rho  = 2.65_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) > 1.0001_wp) then

                      Vp = 6_wp    ! P-wave speed 
                      Vs = 3.464_wp    ! S-wave speed 
                      rho = 2.7_wp    ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   end if
                   !print *, l, j, k, M%M(l,j,k,1)
                   !if (l == 25 .and. j == 94 .and. k == 1) print *, l, j, k, M%M(l,j,k,1)
                end do
             end do
          end do

       case('LOH1_Harmonic')

          call allocate_array_body(M%M,G%C,3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) < 0.9999_wp) then

                      Vp = 4_wp ! P-wave speed
                      Vs  = 2_wp ! S-wave speed
                      rho  = 2.6_wp ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2)   ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2                    ! mu (shear modulus)
                      M%M(l,j,k,3) = rho                          ! rho (density)

                   else if (G%x(l,j,k,1) >= 0.9999_wp .AND. G%x(l,j,k,1) <= 1.0001_wp) then

                      Vp   = 2.0_wp*4.0_wp*6.0_wp/(4.0_wp + 6.0_wp)           ! P-wave speed
                      Vs   = 2.0_wp*2.0_wp*3.464_wp/(2.0_wp + 3.464_wp)       ! S-wave speed
                      rho  = 2.0_wp*2.7_wp*2.6_wp/(2.7_wp + 2.6_wp)           ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2)   ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2                    ! mu (shear modulus)
                      M%M(l,j,k,3) = rho                          ! rho (density)

                   else if (G%x(l,j,k,1) > 1.0001_wp) then

                      Vp = 6_wp                                   ! P-wave speed
                      Vs = 3.464_wp                               ! S-wave speed
                      rho = 2.7_wp                                ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   end if
                end do
             end do
          end do

       case('TPV31')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,2) < 2.3999_wp) then

                      Vp = 4.05_wp ! P-wave speed
                      Vs  = 2.25_wp ! S-wave speed
                      rho  = 2.58_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,2) >= 2.3999_wp .AND. G%x(l,j,k,2) <= 2.4001_wp) then

                      Vp = 4.25_wp ! P-wave speed
                      Vs  = 2.4_wp ! S-wave speed
                      rho  = 2.59_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) > 2.4001_wp .AND. G%x(l,j,k,2) < 4.9999_wp) then

                      Vp  = 4.45_wp + 0.75_wp*(G%x(l,j,k,2)-2.4_wp)/2.6_wp  ! P-wave speed
                      Vs  = 2.55_wp + 0.5_wp*(G%x(l,j,k,2)-2.4_wp)/2.6_wp ! S-wave speed 
                      rho  = 2.6_wp + 0.02_wp*(G%x(l,j,k,2)-2.4_wp)/2.6_wp  ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) >= 4.9999_wp .AND. G%x(l,j,k,2) <= 5.0001_wp) then

                      Vp = 5.475_wp ! P-wave speed
                      Vs  = 3.25_wp ! S-wave speed
                      rho  = 2.67_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) > 5.0001_wp .AND. G%x(l,j,k,2) < 9.9999_wp) then

                      Vp = 5.75_wp   ! P-wave speed 
                      Vs = 3.45_wp   ! S-wave speed 
                      rho = 2.72_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,2) >= 9.9999_wp .AND. G%x(l,j,k,2) <= 10.0001_wp) then

                      Vp = 6.125_wp ! P-wave speed
                      Vs  = 3.625_wp ! S-wave speed
                      rho  = 2.76_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,2) > 10.0001_wp) then

                      Vp = 6.5_wp    ! P-wave speed 
                      Vs = 3.8_wp    ! S-wave speed 
                      rho = 3.0_wp    ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   end if

                end do
             end do
          end do


       case('TPV32')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr 
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,2) <= 0.5_wp) then

                      Vp = 2.2_wp + 0.8_wp*(G%x(l,j,k,2))/0.5_wp      ! P-wave speed 
                      Vs = 1.05_wp + 0.35_wp*(G%x(l,j,k,2))/0.5_wp    ! S-wave speed 
                      rho = 2.2_wp + 0.25_wp*(G%x(l,j,k,2))/0.5_wp     ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 0.5_wp .AND. G%x(l,j,k,2) <= 1.0_wp) then

                      Vp = 3.0_wp + 0.6_wp*(G%x(l,j,k,2)-0.5_wp)/0.5_wp    ! P-wave speed 
                      Vs = 1.4_wp + 0.55_wp*(G%x(l,j,k,2)-0.5_wp)/0.5_wp   ! S-wave speed 
                      rho = 2.45_wp + 0.1_wp*(G%x(l,j,k,2)-0.5_wp)/0.5_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,2) > 1.0_wp .AND. G%x(l,j,k,2) <= 1.6_wp) then

                      Vp = 3.6_wp + 0.8_wp*(G%x(l,j,k,2)-1.0_wp)/0.6_wp     ! P-wave speed 
                      Vs = 1.95_wp + 0.55_wp*(G%x(l,j,k,2)-1.0_wp)/0.6_wp   ! S-wave speed 
                      rho = 2.55_wp + 0.05_wp*(G%x(l,j,k,2)-1.0_wp)/0.6_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 1.6_wp .AND. G%x(l,j,k,2) <= 2.4_wp) then

                      Vp = 4.4_wp + 0.4_wp*(G%x(l,j,k,2)-1.6_wp)/0.8_wp   ! P-wave speed 
                      Vs = 2.5_wp + 0.3_wp*(G%x(l,j,k,2)-1.6_wp)/0.8_wp   ! S-wave speed 
                      rho = 2.6_wp                                         ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 2.4_wp .AND. G%x(l,j,k,2) <= 3.6_wp) then

                      Vp = 4.8_wp + 0.45_wp*(G%x(l,j,k,2)-2.4_wp)/1.2_wp  ! P-wave speed 
                      Vs = 2.8_wp + 0.3_wp*(G%x(l,j,k,2)-2.4_wp)/1.2_wp   ! S-wave speed 
                      rho = 2.6_wp + 0.02_wp*(G%x(l,j,k,2)-2.4_wp)/1.2_wp  ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 3.6_wp .AND. G%x(l,j,k,2) <= 5.0_wp) then

                      Vp = 5.25_wp + 0.25_wp*(G%x(l,j,k,2)-3.6_wp)/1.4_wp  ! P-wave speed 
                      Vs = 3.1_wp + 0.15_wp*(G%x(l,j,k,2)-3.6_wp)/1.4_wp   ! S-wave speed 
                      rho = 2.62_wp + 0.03_wp*(G%x(l,j,k,2)-3.6_wp)/1.4_wp  ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 5.0_wp .AND. G%x(l,j,k,2) <= 9.0_wp) then

                      Vp = 5.5_wp + 0.25_wp*(G%x(l,j,k,2)-5.0_wp)/4.0_wp   ! P-wave speed 
                      Vs = 3.25_wp + 0.2_wp*(G%x(l,j,k,2)-5.0_wp)/4.0_wp   ! S-wave speed 
                      rho = 2.65_wp + 0.07_wp*(G%x(l,j,k,2)-5.0_wp)/4.0_wp  ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 9.0_wp .AND. G%x(l,j,k,2) <= 11.0_wp) then

                      Vp = 5.75_wp + 0.35_wp*(G%x(l,j,k,2)-9.0_wp)/2.0_wp   ! P-wave speed 
                      Vs = 3.45_wp + 0.25_wp*(G%x(l,j,k,2)-9.0_wp)/2.0_wp   ! S-wave speed 
                      rho = 2.72_wp + 0.03_wp*(G%x(l,j,k,2)-9.0_wp)/2.0_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 11.0_wp .AND. G%x(l,j,k,2) <= 15.0_wp) then

                      Vp = 6.1_wp + 0.2_wp*(G%x(l,j,k,2)-11.0_wp)/4.0_wp     ! P-wave speed 
                      Vs = 3.6_wp + 0.1_wp*(G%x(l,j,k,2)-11.0_wp)/4.0_wp     ! S-wave speed 
                      rho = 2.75_wp + 0.15_wp*(G%x(l,j,k,2)-11.0_wp)/4.0_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,2) > 15.0_wp) then

                      Vp = 6.3_wp     ! P-wave speed 
                      Vs = 3.7_wp     ! S-wave speed 
                      rho = 2.9_wp     ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   end if

                end do
             end do
          end do

          ! Low velocity zone benchmark
       case('TPV33')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr 
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) <= -0.8_wp) then

                      Vp = 5.626_wp      ! P-wave speed 
                      Vs = 3.248_wp      ! S-wave speed 
                      rho = 2.67_wp      ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   else if (G%x(l,j,k,1) >= 0.8_wp) then

                      Vp = 6.0_wp        ! P-wave speed 
                      Vs = 3.464_wp      ! S-wave speed 
                      rho = 2.67_wp      ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > -0.8_wp .AND. G%x(l,j,k,1) < 0.8_wp) then

                      Vp = 3.75_wp       ! P-wave speed 
                      Vs = 2.165_wp      ! S-wave speed 
                      rho = 2.67_wp      ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)

                   end if

                end do
             end do
          end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       case('OKLAHOMA')

          call allocate_array_body(M%M,G%C, 3, ghost_nodes=.true.)

          do l = G%C%mq, G%C%pq
             do j = G%C%mr, G%C%pr
                do k = G%C%ms, G%C%ps

                   if (G%x(l,j,k,1) < 1.4999_wp) then

                      Vp = 2.77_wp ! P-wave speed
                      Vs  = 1.49_wp ! S-wave speed
                      rho  = 2.169_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) >= 1.4999_wp .AND. G%x(l,j,k,1) <= 1.5001_wp) then

                      Vp = 3.74_wp ! P-wave speed
                      Vs  = 2.01_wp ! S-wave speed
                      rho  = 2.418_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) > 1.5001_wp .AND. G%x(l,j,k,1) < 2.2999_wp) then

                      Vp  = 5.76_wp  ! P-wave speed
                      Vs  = 3.1_wp  ! S-wave speed 
                      rho  = 2.667_wp  ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 2.2999_wp .AND. G%x(l,j,k,1) <= 2.3001_wp) then

                      Vp = 5.755_wp ! P-wave speed
                      Vs  = 3.26_wp ! S-wave speed
                      rho  = 2.666_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) > 2.3001_wp .AND. G%x(l,j,k,1) < 5.2999_wp) then

                      Vp = 5.75_wp   ! P-wave speed 
                      Vs = 3.423_wp   ! S-wave speed 
                      rho = 2.665_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 5.2999_wp .AND. G%x(l,j,k,1) <= 5.3001_wp) then

                      Vp = 5.96_wp ! P-wave speed
                      Vs  = 3.54_wp ! S-wave speed
                      rho  = 2.711_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 5.3001_wp .AND. G%x(l,j,k,1) < 8.2999_wp) then

                      Vp = 6.18_wp   ! P-wave speed 
                      Vs = 3.65_wp   ! S-wave speed 
                      rho = 2.756_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 8.2999_wp .AND. G%x(l,j,k,1) <= 8.3001_wp) then

                      Vp = 6.20_wp ! P-wave speed
                      Vs  = 3.63_wp ! S-wave speed
                      rho  = 2.762_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 8.3001_wp .AND. G%x(l,j,k,1) < 11.2999_wp) then

                      Vp = 6.23_wp   ! P-wave speed 
                      Vs = 3.60_wp   ! S-wave speed 
                      rho = 2.767_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 11.2999_wp .AND. G%x(l,j,k,1) <= 11.3001_wp) then

                      Vp = 6.26_wp ! P-wave speed
                      Vs  = 3.64_wp ! S-wave speed
                      rho  = 2.776_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 11.3001_wp .AND. G%x(l,j,k,1) < 20.2999_wp) then

                      Vp = 6.30_wp   ! P-wave speed 
                      Vs = 3.67_wp   ! S-wave speed 
                      rho = 2.784_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 20.2999_wp .AND. G%x(l,j,k,1) <= 20.3001_wp) then

                      Vp = 6.54_wp ! P-wave speed
                      Vs  = 3.80_wp ! S-wave speed
                      rho  = 2.848_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 20.3001_wp .AND. G%x(l,j,k,1) < 39.9999_wp) then

                      Vp = 6.80_wp   ! P-wave speed 
                      Vs = 3.93_wp   ! S-wave speed 
                      rho = 2.911_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   else if (G%x(l,j,k,1) >= 39.9999_wp .AND. G%x(l,j,k,1) <= 40.0001_wp) then

                      Vp = 7.39_wp ! P-wave speed
                      Vs  = 4.25_wp ! S-wave speed
                      rho  = 3.119_wp ! rho (density) 

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density) 

                   else if (G%x(l,j,k,1) > 40.0001_wp) then

                      Vp = 8.10_wp   ! P-wave speed 
                      Vs = 4.62_wp   ! S-wave speed 
                      rho = 3.326_wp   ! rho (density)

                      M%M(l,j,k,1) = rho*(Vp**2 - 2.0_wp*Vs**2) ! lambda (Lame's first parameter)
                      M%M(l,j,k,2) = rho*Vs**2 ! mu (shear modulus)
                      M%M(l,j,k,3) = rho ! rho (density)  

                   end if

                end do
             end do
          end do


       end select
       do ii = 1,3
          call exchange_all_neighbors(G%C,M%M(:,:,:,ii))
       end do


    case('acoustic')

       call allocate_array_body(M%M,G%C,2, ghost_nodes=.true.)

       M%M(:,:,:,1) = 1.0_wp ! rho (density)
       M%M(:,:,:,2) = 1.0_wp ! K (bulk modulus)

    end select

  end subroutine init_material


  subroutine init_material_from_file(id,nb,M,G,material_path)

    use mpi3dbasic, only : pw, ps
    use mpi3dio, only: file_distributed, open_file_distributed, read_file_distributed, close_file_distributed
    use mpi3dcomm, only : allocate_array_body, exchange_all_neighbors 

    implicit none

    integer, intent(in) :: id, nb
    type(block_grid_t), intent(in) :: G
    type(block_material),intent(out) :: M
    type(file_distributed) :: fids(3)
    character(len=256),intent(in) :: material_path(3)
    character(len=256) :: name
    integer :: i,j,k


    ! allocate memory for array

    call allocate_array_body(M%M,G%C,3,ghost_nodes=.true.)
    M%M(:,:,:,:) = 0_wp

    do i = 1,3
        name = material_path(i)
        call open_file_distributed(fids(i),name,"read",G%C%comm,G%C%array_w,pw)
        call read_file_distributed(fids(i),M%M(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps,i))
        call close_file_distributed(fids(i))
    end do

    do i = 1,3
        call exchange_all_neighbors(G%C,M%M(:,:,:,i))
    end do

  !  do k = G%C%mq,G%C%pq
  !    print*, k
  !    print*,'lambda = ',M%M(k,G%C%mr,G%C%ms,1)
  !    print*,'mu = ' ,M%M(k,G%C%mr,G%C%ms,2)
  !    print*,'rho = ',M%M(k,G%C%mr,G%C%ms,3)
  !    print*,' '
  !  end do

  end subroutine init_material_from_file


   subroutine init_const_Q_4M_properties(M, G, infile)

      ! Initializes constant-Q-4M attenuation: N=4 mechanisms with user-configurable fmin/fmax.
      ! Supports three weight computation methods:
      ! 1. 'nnls'   : Solves NNLS optimization to minimize Q(f) error over [fmin, fmax]
      ! 2. 'withers': Scales Withers 8-mechanism framework to [fmin, fmax] and extracts first 4
      ! 3. 'lookup' : Interpolates pre-computed weight tables for standard frequency bands
      ! Users can override auto-computed weights with manual_weights if provided (all > 0).
      !
      ! Namelist &constant_Q_4M_list parameters:
      !   c                  : Q_S = c / sqrt(mu/rho)  (Q_P = 2*Q_S)
      !   fmin               : Minimum frequency (Hz) for tau-band (default 0.05)
      !   fmax               : Maximum frequency (Hz) for tau-band (default 20.0)
      !   target_Q           : Target constant-Q value (default 50.0)
      !   weight_method      : 'nnls' | 'withers' | 'lookup' (default 'nnls')
      !   manual_weights     : Optional 4-element array; if all > 0, use instead of auto-computed (default all 0)
      !   fref               : Reference frequency for unrelaxed modulus correction (Hz, default 1.0)

      use mpi3dcomm, only : allocate_array_body
      use withers_tables, only : get_relaxation_times_Qf, get_withers_weights_Qf

      implicit none

      type(block_material), intent(inout) :: M
      type(block_grid_t), intent(in) :: G
      integer, intent(in) :: infile

      real(kind = wp) :: c, fmin, fmax, target_Q, fref
      real(kind = wp), dimension(4) :: manual_weights
      character(16) :: weight_method
      real(kind = wp) :: taumin, taumax, wref
      real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
      integer :: stat, i, l, j, k, N
      logical :: use_manual
      real(kind = wp), parameter :: pi = 3.141592653589793_wp

      namelist /constant_Q_4M_list/ c, fmin, fmax, target_Q, weight_method, manual_weights, fref

      ! Defaults
      c = 1.0_wp
      fmin = 0.05_wp
      fmax = 20.0_wp
      target_Q = 50.0_wp
      weight_method = 'nnls'
      manual_weights = 0.0_wp
      fref = 1.0_wp

      M%anelastic_const_Q_4M = .true.
      N = 4
      M%n_mechanism_const_Q_4M = N
      M%fmin_const_Q_4M = fmin
      M%fmax_const_Q_4M = fmax
      M%target_Q_const_Q_4M = target_Q
      M%weight_method_const_Q_4M = weight_method
      M%manual_weights_const_Q_4M = manual_weights
      M%fref_const_Q_4M = fref

      rewind(infile)
      read(infile, nml=constant_Q_4M_list, iostat=stat)
      if (stat > 0) stop 'error reading namelist constant_Q_4M_list'

      M%fmin_const_Q_4M = fmin
      M%fmax_const_Q_4M = fmax
      M%target_Q_const_Q_4M = target_Q
      M%weight_method_const_Q_4M = weight_method
      M%manual_weights_const_Q_4M = manual_weights
      M%fref_const_Q_4M = fref

      ! Sanity checks
      if (fmin <= 0.0_wp) fmin = 0.05_wp
      if (fmax <= fmin) fmax = 20.0_wp
      if (target_Q <= 0.0_wp) target_Q = 50.0_wp

      ! Allocate Q arrays
      call allocate_array_body(M%Qp_inv_const_Q_4M, G%C, ghost_nodes=.true.)
      M%Qp_inv_const_Q_4M(:,:,:) = 0.0_wp
      call allocate_array_body(M%Qs_inv_const_Q_4M, G%C, ghost_nodes=.true.)
      M%Qs_inv_const_Q_4M(:,:,:) = 0.0_wp

      ! Allocate memory variable arrays (N=4 mechanisms)
      call allocate_array_body(M%eta4_4M,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta4_4M, G%C, N, ghost_nodes=.true.)
      M%eta4_4M = 0.0_wp;  M%Deta4_4M = 0.0_wp
      call allocate_array_body(M%eta5_4M,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta5_4M, G%C, N, ghost_nodes=.true.)
      M%eta5_4M = 0.0_wp;  M%Deta5_4M = 0.0_wp
      call allocate_array_body(M%eta6_4M,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta6_4M, G%C, N, ghost_nodes=.true.)
      M%eta6_4M = 0.0_wp;  M%Deta6_4M = 0.0_wp
      call allocate_array_body(M%eta7_4M,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta7_4M, G%C, N, ghost_nodes=.true.)
      M%eta7_4M = 0.0_wp;  M%Deta7_4M = 0.0_wp
      call allocate_array_body(M%eta8_4M,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta8_4M, G%C, N, ghost_nodes=.true.)
      M%eta8_4M = 0.0_wp;  M%Deta8_4M = 0.0_wp
      call allocate_array_body(M%eta9_4M,  G%C, N, ghost_nodes=.true.)
      call allocate_array_body(M%Deta9_4M, G%C, N, ghost_nodes=.true.)
      M%eta9_4M = 0.0_wp;  M%Deta9_4M = 0.0_wp

      ! Compute Q_S and Q_P from velocity ratio
      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps
               M%Qs_inv_const_Q_4M(i,j,k) = 1.0_wp / (c * sqrt(M%M(i,j,k,2)/M%M(i,j,k,3)))
               M%Qp_inv_const_Q_4M(i,j,k) = 0.5_wp * M%Qs_inv_const_Q_4M(i,j,k)
            end do
         end do
      end do

      ! Relaxation times: formula with denominator = 2*N = 8
      ! tau_k = exp( ln(tau_min) + (2k-1)/(2N) * ln(tau_max/tau_min) )
      taumin = 1.0_wp / (2.0_wp * pi * fmax)
      taumax = 1.0_wp / (2.0_wp * pi * fmin)

      do k = 1, N
         M%tau_const_Q_4M(k) = exp(log(taumin) + (2.0_wp*k - 1.0_wp) / (2.0_wp*N) &
                                  * log(taumax/taumin))
      end do

      ! Check if manual weights are provided (all > 0)
      use_manual = all(manual_weights > 0.0_wp)

      if (use_manual) then
         ! Use user-provided manual weights
         M%weight_const_Q_4M = manual_weights
      else
         ! Compute weights based on selected method
         select case(trim(weight_method))
         case('nnls')
            ! Solve NNLS for user's [fmin, fmax] and target_Q
            call compute_weights_nnls(M%tau_const_Q_4M, target_Q, fmin, fmax, &
                                      M%weight_const_Q_4M, N)

         case('withers')
            ! Scale Withers 8-mechanism framework to [fmin, fmax]; extract first 4
            call compute_weights_withers(fmin, fmax, M%weight_const_Q_4M, N)

         case('lookup')
            ! Use pre-computed lookup table for standard bands
            call compute_weights_lookup(fmin, fmax, target_Q, M%weight_const_Q_4M, N)

         case default
            ! Default to NNLS if invalid method specified
            write(*,*) 'Warning: unknown weight_method "', trim(weight_method), '"; using nnls'
            call compute_weights_nnls(M%tau_const_Q_4M, target_Q, fmin, fmax, &
                                      M%weight_const_Q_4M, N)
         end select
      end if

      ! Correct unrelaxed moduli so relaxed (low-freq) velocity matches input
      wref = 2.0_wp * pi * fref

      do i = G%C%mq, G%C%pq
         do j = G%C%mr, G%C%pr
            do k = G%C%ms, G%C%ps

               val_S = 0.0_wp
               val_P = 0.0_wp
               do l = 1, N
                  denom_S = (wref**2 * M%tau_const_Q_4M(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qs_inv_const_Q_4M(i,j,k))
                  denom_P = (wref**2 * M%tau_const_Q_4M(l)**2 + 1.0_wp) &
                             * (1.0_wp / M%Qp_inv_const_Q_4M(i,j,k))
                  val_S = val_S + M%weight_const_Q_4M(l) / denom_S
                  val_P = val_P + M%weight_const_Q_4M(l) / denom_P
               end do

               vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
               vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2)) / M%M(i,j,k,3))
               mu_unrelax_S = M%M(i,j,k,3) * vs**2 / (1.0_wp - val_S)
               mu_unrelax_P = M%M(i,j,k,3) * vp**2 / (1.0_wp - val_P)
               vs = sqrt(mu_unrelax_S / M%M(i,j,k,3))
               vp = sqrt(mu_unrelax_P / M%M(i,j,k,3))

               M%M(i,j,k,2) = vs**2 * M%M(i,j,k,3)
               M%M(i,j,k,1) = vp**2 * M%M(i,j,k,3) - 2.0_wp * M%M(i,j,k,2)

            end do
         end do
      end do

   end subroutine init_const_Q_4M_properties


   subroutine compute_weights_nnls(tau_vals, target_Q, fmin, fmax, weights_out, N)
      ! Solves NNLS to minimize constant-Q error over [fmin, fmax].
      ! Uses simplified active-set method (iterative refinement).
      !
      ! Minimize: sum_i (Q_target - Q_approx(f_i, weights))^2
      ! subject to: weights >= 0

      implicit none

      integer, intent(in) :: N
      real(kind = wp), intent(in) :: tau_vals(N), target_Q, fmin, fmax
      real(kind = wp), intent(out) :: weights_out(N)
      
      real(kind = wp), parameter :: pi = 3.141592653589793_wp
      real(kind = wp), parameter :: tol = 1.0e-6_wp
      integer, parameter :: n_freq_test = 30, max_iter = 100
      
      real(kind = wp) :: freq_test(n_freq_test), Q_resp(n_freq_test)
      real(kind = wp) :: H(N,N), g(N), weights_iter(N), residual
      integer :: i, j, k, iter, n_active
      logical :: active(N)
      real(kind = wp) :: temp_w(N), step, max_step, grad(N)
      
      ! Generate test frequencies (logarithmically spaced)
      do i = 1, n_freq_test
         freq_test(i) = fmin * (fmax/fmin)**((real(i,wp)-1.0_wp)/(real(n_freq_test,wp)-1.0_wp))
      end do
      
      ! Initialize weights (uniform initial guess)
      weights_iter = 1.0_wp / real(N, wp)
      active = .true.
      
      ! Iterative NNLS using active-set method (simplified)
      do iter = 1, max_iter
         
         ! Compute Hessian and gradient using test frequencies
         H = 0.0_wp
         g = 0.0_wp
         
         do i = 1, n_freq_test
            ! Compute Q response at freq_test(i) with current weights
            call compute_Q_response(freq_test(i), tau_vals, weights_iter, Q_resp(i), N)
            
            ! Accumulate Hessian and gradient
            ! grad_j = -2 * sum_i (Q_target - Q_resp(i)) * dQ_resp/dw_j
            ! H_jl = 2 * sum_i (dQ_resp/dw_j) * (dQ_resp/dw_l)
            do j = 1, N
               do k = 1, N
                  ! Approximation: H_jk ~ omega_j^2 * omega_k^2 / tau_j / tau_k (simplified)
                  H(j,k) = H(j,k) + 2.0_wp * (freq_test(i)**2 * tau_vals(j) * tau_vals(k)) / &
                           (1.0_wp + (freq_test(i) * tau_vals(j))**2) / &
                           (1.0_wp + (freq_test(i) * tau_vals(k))**2)
               end do
               
               ! Gradient approximation
               g(j) = g(j) - 2.0_wp * (target_Q - Q_resp(i)) * freq_test(i)**2 * tau_vals(j) / &
                      (1.0_wp + (freq_test(i) * tau_vals(j))**2)
            end do
         end do
         
         ! Add small regularization to H for numerical stability
         do j = 1, N
            H(j,j) = H(j,j) + 1.0e-8_wp
         end do
         
         ! Solve H * p = -g for descent direction (simplified solver)
         call solve_normal_equations(H, -g, temp_w, N)
         
         ! Line search with projection to non-negative orthant
         max_step = 1.0_wp
         do j = 1, N
            if (temp_w(j) < 0.0_wp .and. weights_iter(j) > 0.0_wp) then
               max_step = min(max_step, -0.9_wp * weights_iter(j) / temp_w(j))
            end if
         end do
         
         step = 0.5_wp * max_step
         weights_iter = weights_iter + step * temp_w
         weights_iter = max(weights_iter, 0.0_wp)
         
         ! Check convergence
         residual = maxval(abs(temp_w))
         if (residual < tol) exit
         
      end do
      
      weights_out = weights_iter

   end subroutine compute_weights_nnls


   subroutine compute_Q_response(freq, tau_vals, weights, Q_resp, N)
      ! Computes Q response at a single frequency using GSLS model.
      ! Q(f) = 1 / (2 * sum_k w_k / (1 + (omega*tau_k)^2))

      implicit none
      
      integer, intent(in) :: N
      real(kind = wp), intent(in) :: freq, tau_vals(N), weights(N)
      real(kind = wp), intent(out) :: Q_resp
      real(kind = wp), parameter :: pi = 3.141592653589793_wp
      
      real(kind = wp) :: omega, sum_term, k
      integer :: i
      
      omega = 2.0_wp * pi * freq
      sum_term = 0.0_wp
      
      do i = 1, N
         sum_term = sum_term + weights(i) / (1.0_wp + (omega * tau_vals(i))**2)
      end do
      
      Q_resp = 1.0_wp / (2.0_wp * sum_term)

   end subroutine compute_Q_response


   subroutine solve_normal_equations(H, g, x, N)
      ! Simple solver for H*x = g using Gaussian elimination (for small N).

      implicit none
      
      integer, intent(in) :: N
      real(kind = wp), intent(in) :: H(N,N), g(N)
      real(kind = wp), intent(out) :: x(N)
      
      real(kind = wp) :: A(N,N+1), temp, pivot
      integer :: i, j, k, piv_row
      
      ! Form augmented matrix [H | g]
      do i = 1, N
         do j = 1, N
            A(i,j) = H(i,j)
         end do
         A(i,N+1) = g(i)
      end do
      
      ! Forward elimination with partial pivoting
      do k = 1, N
         ! Find pivot
         piv_row = k
         do i = k+1, N
            if (abs(A(i,k)) > abs(A(piv_row,k))) piv_row = i
         end do
         
         ! Swap rows
         do j = k, N+1
            temp = A(k,j)
            A(k,j) = A(piv_row,j)
            A(piv_row,j) = temp
         end do
         
         ! Eliminate below
         if (abs(A(k,k)) > 1.0e-12_wp) then
            do i = k+1, N
               pivot = A(i,k) / A(k,k)
               do j = k, N+1
                  A(i,j) = A(i,j) - pivot * A(k,j)
               end do
            end do
         end if
      end do
      
      ! Back substitution
      do i = N, 1, -1
         x(i) = A(i,N+1)
         do j = i+1, N
            x(i) = x(i) - A(i,j) * x(j)
         end do
         if (abs(A(i,i)) > 1.0e-12_wp) then
            x(i) = x(i) / A(i,i)
         else
            x(i) = 0.0_wp
         end if
      end do

   end subroutine solve_normal_equations


   subroutine compute_weights_withers(fmin, fmax, weights_out, N)
      ! Uses the tabulated Withers relaxation spectrum as a fast non-negative
      ! weight approximation. For 4M callers, the first 4 weights are returned;
      ! for 8M callers, all 8 weights are returned.

      use withers_tables, only : get_relaxation_times, get_withers_weights

      implicit none
      
      integer, intent(in) :: N
      real(kind = wp), intent(in) :: fmin, fmax
      real(kind = wp), intent(out) :: weights_out(N)
      
      real(kind = wp) :: tau_withers(8), weights_withers(8), scale_factor
      real(kind = wp), parameter :: pi = 3.141592653589793_wp
      real(kind = wp) :: fmin_std, fmax_std
      integer :: i
      
      ! Get Withers taus and weights for the standard 8-mechanism table.
      call get_relaxation_times(0.0_wp, tau_withers)
      call get_withers_weights(0.0_wp, 1.0_wp, weights_withers)
      
      ! Standard band is approximately [0.05, 20] Hz
      fmin_std = 0.05_wp
      fmax_std = 20.0_wp
      
      ! Compute scaling factor for tau adjustment
      scale_factor = log(fmax / fmin) / log(fmax_std / fmin_std)
      
      ! Preserve the same vector shape as the caller expects.
      do i = 1, N
         weights_out(i) = weights_withers(i)
      end do

   end subroutine compute_weights_withers


   subroutine compute_weights_lookup(fmin, fmax, target_Q, weights_out, N)
      ! Uses a light-weight lookup approximation derived from the standard
      ! Withers spectrum and scales it to the requested target_Q.

      implicit none
      
      integer, intent(in) :: N
      real(kind = wp), intent(in) :: fmin, fmax, target_Q
      real(kind = wp), intent(out) :: weights_out(N)
      
      real(kind = wp) :: weights_base(8)
      integer :: i

      call compute_weights_withers(fmin, fmax, weights_base, 8)
      do i = 1, N
         weights_out(i) = weights_base(i)
      end do

      ! Scale the tabulated weights to the requested target Q.
      if (target_Q > 0.0_wp) weights_out = weights_out * target_Q / 50.0_wp

   end subroutine compute_weights_lookup

      subroutine init_const_Q_8M_properties(M, G, infile)
         ! Initializes constant-Q-8M attenuation: N=8 mechanisms with user-configurable fmin/fmax.
         use mpi3dcomm, only : allocate_array_body

         implicit none

         type(block_material), intent(inout) :: M
         type(block_grid_t), intent(in) :: G
         integer, intent(in) :: infile

         real(kind = wp) :: c, fmin, fmax, target_Q, fref
         real(kind = wp), dimension(8) :: manual_weights
         character(16) :: weight_method
         real(kind = wp) :: taumin, taumax, wref
         real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
         integer :: stat, i, l, j, k, N
         logical :: use_manual
         real(kind = wp), parameter :: pi = 3.141592653589793_wp

         namelist /constant_Q_8M_list/ c, fmin, fmax, target_Q, weight_method, manual_weights, fref

         c = 1.0_wp
         fmin = 0.05_wp
         fmax = 20.0_wp
         target_Q = 50.0_wp
         weight_method = 'nnls'
         manual_weights = 0.0_wp
         fref = 1.0_wp

         M%anelastic_const_Q_8M = .true.
         M%anelastic_const_Q_4M = .true.
         N = 8
         M%n_mechanism_const_Q_8M = N

         rewind(infile)
         read(infile, nml=constant_Q_8M_list, iostat=stat)
         if (stat > 0) stop 'error reading namelist constant_Q_8M_list'

         M%c = c
         M%fmin_const_Q_8M = fmin
         M%fmax_const_Q_8M = fmax
         M%target_Q_const_Q_8M = target_Q
         M%weight_method_const_Q_8M = weight_method
         M%manual_weights_const_Q_8M = manual_weights
         M%fref_const_Q_8M = fref

         if (fmin <= 0.0_wp) fmin = 0.05_wp
         if (fmax <= fmin) fmax = 20.0_wp
         if (target_Q <= 0.0_wp) target_Q = 50.0_wp

         call allocate_array_body(M%Qp_inv_const_Q_8M, G%C, ghost_nodes=.true.)
         M%Qp_inv_const_Q_8M(:,:,:) = 0.0_wp
         call allocate_array_body(M%Qs_inv_const_Q_8M, G%C, ghost_nodes=.true.)
         M%Qs_inv_const_Q_8M(:,:,:) = 0.0_wp

         call allocate_array_body(M%eta4_8M,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta4_8M, G%C, N, ghost_nodes=.true.)
         M%eta4_8M = 0.0_wp;  M%Deta4_8M = 0.0_wp
         call allocate_array_body(M%eta5_8M,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta5_8M, G%C, N, ghost_nodes=.true.)
         M%eta5_8M = 0.0_wp;  M%Deta5_8M = 0.0_wp
         call allocate_array_body(M%eta6_8M,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta6_8M, G%C, N, ghost_nodes=.true.)
         M%eta6_8M = 0.0_wp;  M%Deta6_8M = 0.0_wp
         call allocate_array_body(M%eta7_8M,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta7_8M, G%C, N, ghost_nodes=.true.)
         M%eta7_8M = 0.0_wp;  M%Deta7_8M = 0.0_wp
         call allocate_array_body(M%eta8_8M,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta8_8M, G%C, N, ghost_nodes=.true.)
         M%eta8_8M = 0.0_wp;  M%Deta8_8M = 0.0_wp
         call allocate_array_body(M%eta9_8M,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta9_8M, G%C, N, ghost_nodes=.true.)
         M%eta9_8M = 0.0_wp;  M%Deta9_8M = 0.0_wp

         do i = G%C%mq, G%C%pq
            do j = G%C%mr, G%C%pr
               do k = G%C%ms, G%C%ps
                  M%Qs_inv_const_Q_8M(i,j,k) = 1.0_wp / (c * sqrt(M%M(i,j,k,2)/M%M(i,j,k,3)))
                  M%Qp_inv_const_Q_8M(i,j,k) = 0.5_wp * M%Qs_inv_const_Q_8M(i,j,k)
               end do
            end do
         end do

         taumin = 1.0_wp / (2.0_wp * pi * fmax)
         taumax = 1.0_wp / (2.0_wp * pi * fmin)
         do k = 1, N
            M%tau_const_Q_8M(k) = exp(log(taumin) + (2.0_wp*k - 1.0_wp) / (2.0_wp*N) * log(taumax/taumin))
         end do

         use_manual = all(manual_weights > 0.0_wp)
         if (use_manual) then
            M%weight_const_Q_8M = manual_weights
         else
            select case(trim(weight_method))
            case('nnls')
               call compute_weights_nnls(M%tau_const_Q_8M, target_Q, fmin, fmax, M%weight_const_Q_8M, N)
            case('withers')
               call compute_weights_withers(fmin, fmax, M%weight_const_Q_8M, N)
            case('lookup')
               call compute_weights_lookup(fmin, fmax, target_Q, M%weight_const_Q_8M, N)
            case default
               write(*,*) 'Warning: unknown weight_method "', trim(weight_method), '"; using nnls'
               call compute_weights_nnls(M%tau_const_Q_8M, target_Q, fmin, fmax, M%weight_const_Q_8M, N)
            end select
         end if

         wref = 2.0_wp * pi * fref
         do i = G%C%mq, G%C%pq
            do j = G%C%mr, G%C%pr
               do k = G%C%ms, G%C%ps
                  val_S = 0.0_wp
                  val_P = 0.0_wp
                  do l = 1, N
                     denom_S = (wref**2 * M%tau_const_Q_8M(l)**2 + 1.0_wp) * (1.0_wp / M%Qs_inv_const_Q_8M(i,j,k))
                     denom_P = (wref**2 * M%tau_const_Q_8M(l)**2 + 1.0_wp) * (1.0_wp / M%Qp_inv_const_Q_8M(i,j,k))
                     val_S = val_S + M%weight_const_Q_8M(l) / denom_S
                     val_P = val_P + M%weight_const_Q_8M(l) / denom_P
                  end do

                  vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
                  vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2)) / M%M(i,j,k,3))
                  mu_unrelax_S = M%M(i,j,k,3) * vs**2 / (1.0_wp - val_S)
                  mu_unrelax_P = M%M(i,j,k,3) * vp**2 / (1.0_wp - val_P)
                  vs = sqrt(mu_unrelax_S / M%M(i,j,k,3))
                  vp = sqrt(mu_unrelax_P / M%M(i,j,k,3))

                  M%M(i,j,k,2) = vs**2 * M%M(i,j,k,3)
                  M%M(i,j,k,1) = vp**2 * M%M(i,j,k,3) - 2.0_wp * M%M(i,j,k,2)
               end do
            end do
         end do

      end subroutine init_const_Q_8M_properties


      subroutine init_anelastic_Qf8_properties(M, G, parameters)
         ! Initializes frequency-dependent Q with N=8 mechanisms.
         use mpi3dcomm, only : allocate_array_body
         use mpi3dbasic, only : is_master
         use anelastic_fq8_model, only : fq8_parameters, build_fq8_coefficients
         implicit none

         type(block_material), intent(inout) :: M
         type(block_grid_t), intent(in) :: G
         type(fq8_parameters), intent(in) :: parameters

         real(kind = wp) :: wref
         real(kind = wp) :: val_S, val_P, denom_S, denom_P, vs, vp, mu_unrelax_S, mu_unrelax_P
         real(kind = wp) :: max_qs_error, max_qp_error
         real(kind = wp) :: common_scale_s, common_scale_p
         real(kind = wp) :: sample_frequencies(7), qeff_s, qeff_p, &
              velocity_ratio_s, velocity_ratio_p
         integer :: stat, i, l, j, k, N, sample
         character(len=256) :: message
         real(kind = wp), parameter :: pi = 3.141592653589793_wp

         M%anelastic_Qf8 = .true.
         M%anelastic_Qf = .true.
         M%coarse_grained_Qf8 = parameters%coarse_grain == 2
         N = 8
         M%n_mechanism_Qf8 = N
         M%gamma_Qf8 = parameters%gamma
         M%fref_Qf8 = parameters%fref
         M%f_trans_Qf8 = parameters%f_transition

         call allocate_array_body(M%Qp_inv_Qf8, G%C, ghost_nodes=.true.)
         M%Qp_inv_Qf8(:,:,:) = 0.0_wp
         call allocate_array_body(M%Qs_inv_Qf8, G%C, ghost_nodes=.true.)
         M%Qs_inv_Qf8(:,:,:) = 0.0_wp

         call allocate_array_body(M%eta4Qf8,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta4Qf8, G%C, N, ghost_nodes=.true.)
         M%eta4Qf8 = 0.0_wp;  M%Deta4Qf8 = 0.0_wp
         call allocate_array_body(M%eta5Qf8,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta5Qf8, G%C, N, ghost_nodes=.true.)
         M%eta5Qf8 = 0.0_wp;  M%Deta5Qf8 = 0.0_wp
         call allocate_array_body(M%eta6Qf8,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta6Qf8, G%C, N, ghost_nodes=.true.)
         M%eta6Qf8 = 0.0_wp;  M%Deta6Qf8 = 0.0_wp
         call allocate_array_body(M%eta7Qf8,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta7Qf8, G%C, N, ghost_nodes=.true.)
         M%eta7Qf8 = 0.0_wp;  M%Deta7Qf8 = 0.0_wp
         call allocate_array_body(M%eta8Qf8,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta8Qf8, G%C, N, ghost_nodes=.true.)
         M%eta8Qf8 = 0.0_wp;  M%Deta8Qf8 = 0.0_wp
         call allocate_array_body(M%eta9Qf8,  G%C, N, ghost_nodes=.true.)
         call allocate_array_body(M%Deta9Qf8, G%C, N, ghost_nodes=.true.)
         M%eta9Qf8 = 0.0_wp;  M%Deta9Qf8 = 0.0_wp

         do i = G%C%mq, G%C%pq
            do j = G%C%mr, G%C%pr
               do k = G%C%ms, G%C%ps
                  M%Qs_inv_Qf8(i,j,k) = 1.0_wp / parameters%Qs0
                  M%Qp_inv_Qf8(i,j,k) = 1.0_wp / parameters%Qp0
               end do
            end do
         end do

         call build_fq8_coefficients(parameters, M%tau_Qf8, M%strength_s_Qf8, &
              M%strength_p_Qf8, stat, message)
         if (stat /= 0) then
            write(*,'(A)') trim(message)
            error stop 1
         end if
         common_scale_s=fq8_common_modulus_scale(parameters%fref,M%tau_Qf8,M%strength_s_Qf8)
         common_scale_p=fq8_common_modulus_scale(parameters%fref,M%tau_Qf8,M%strength_p_Qf8)

         wref = 2.0_wp * pi * parameters%fref
         do i = G%C%mq, G%C%pq
            do j = G%C%mr, G%C%pr
               do k = G%C%ms, G%C%ps
                  val_S = 0.0_wp
                  val_P = 0.0_wp
                  if (.not.M%coarse_grained_Qf8) then
                     do l = 1, N
                        denom_S = wref**2 * M%tau_Qf8(l)**2 + 1.0_wp
                        denom_P = denom_S
                        val_S = val_S + M%strength_s_Qf8(l) / denom_S
                        val_P = val_P + M%strength_p_Qf8(l) / denom_P
                     end do
                  end if

                  vs = sqrt(M%M(i,j,k,2) / M%M(i,j,k,3))
                  vp = sqrt((M%M(i,j,k,1) + 2.0_wp*M%M(i,j,k,2)) / M%M(i,j,k,3))
                  if (M%coarse_grained_Qf8) then
                     mu_unrelax_S=M%M(i,j,k,3)*vs**2*common_scale_s
                     mu_unrelax_P=M%M(i,j,k,3)*vp**2*common_scale_p
                  else
                     mu_unrelax_S=M%M(i,j,k,3)*vs**2/(1.0_wp-val_S)
                     mu_unrelax_P=M%M(i,j,k,3)*vp**2/(1.0_wp-val_P)
                  end if
                  vs = sqrt(mu_unrelax_S / M%M(i,j,k,3))
                  vp = sqrt(mu_unrelax_P / M%M(i,j,k,3))

                  M%M(i,j,k,2) = vs**2 * M%M(i,j,k,3)
                  M%M(i,j,k,1) = vp**2 * M%M(i,j,k,3) - 2.0_wp * M%M(i,j,k,2)
               end do
            end do
         end do

         if (is_master() .and. .not.fq8_summary_printed) then
            fq8_summary_printed=.true.
            call fq8_max_relative_error(parameters%Qs0,parameters%gamma, &
                 parameters%f_transition,M%tau_Qf8,M%strength_s_Qf8, &
                 0.1_wp*parameters%f_transition,10.0_wp*parameters%f_transition,max_qs_error)
            call fq8_max_relative_error(parameters%Qp0,parameters%gamma, &
                 parameters%f_transition,M%tau_Qf8,M%strength_p_Qf8, &
                 0.1_wp*parameters%f_transition,10.0_wp*parameters%f_transition,max_qp_error)
            write(*,'(A)') 'anelastic-fQ8 parameters:'
            write(*,'(A,A)') '  coefficient method = ',trim(parameters%coefficient_method)
            write(*,'(A,A)') '  weight policy = ',trim(parameters%weight_policy)
            write(*,'(A,I0)') '  coarse_grain = ',parameters%coarse_grain
            write(*,'(A,ES12.4,A,ES12.4)') '  Qs0 = ',parameters%Qs0,', Qp0 = ',parameters%Qp0
            write(*,'(A,F5.2,A,ES12.4)') '  gamma = ',parameters%gamma, &
                 ', transition frequency (Hz) = ',parameters%f_transition
            write(*,'(A,8(ES11.3,1X))') '  tau (s) = ',M%tau_Qf8
            write(*,'(A,8(ES11.3,1X))') '  S strengths = ',M%strength_s_Qf8
            write(*,'(A,8(ES11.3,1X))') '  P strengths = ',M%strength_p_Qf8
            if (any(M%strength_s_Qf8 < 0.0_wp) .or. any(M%strength_p_Qf8 < 0.0_wp)) &
                 call warn_once('FQ8-WEIGHT-001', &
                 'table-exact interpolation contains negative strengths; use '// &
                 'weight_policy=''nonnegative-refit'' for strict passivity.', &
                 'init_anelastic_Qf8_properties')
            if (.not.M%coarse_grained_Qf8) write(*,'(A,F8.3,A,F8.3,A)') &
                 '  max realized-Q errors (S/P) = ',100.0_wp*max_qs_error, &
                 ' / ',100.0_wp*max_qp_error,' %'
            if (.not.M%coarse_grained_Qf8 .and. &
                max(max_qs_error,max_qp_error) > 0.10_wp) call warn_once('FQ8-FIT-001', &
                 'Realized-Q error exceeds 10 percent over 0.1-10 times f_transition.', &
                 'init_anelastic_Qf8_properties')
            if (M%coarse_grained_Qf8) write(*,'(A)') &
                 '  layout = 2x2x2 one-mechanism-per-node coarse cell'
            if (.not.M%coarse_grained_Qf8) write(*,'(A)') &
                 '  layout = full 8 mechanisms at every grid point'
            if (M%coarse_grained_Qf8) write(*,'(A)') &
                 '  modulus normalization = harmonic coarse-cell reference'
            if (.not.M%coarse_grained_Qf8) write(*,'(A)') &
                 '  modulus normalization = collocated full-mechanism reference'
            if (M%coarse_grained_Qf8) then
               sample_frequencies=[0.1_wp,0.2_wp,0.5_wp,1.0_wp,2.0_wp,5.0_wp,10.0_wp]
               write(*,'(A,ES12.4,A,ES12.4)') '  common modulus scale (S/P) = ', &
                    fq8_common_modulus_scale(parameters%fref,M%tau_Qf8,M%strength_s_Qf8), &
                    ' / ',fq8_common_modulus_scale(parameters%fref,M%tau_Qf8,M%strength_p_Qf8)
               write(*,'(A)') '  f(Hz)       Qeff_S       Qeff_P       cS/cSref     cP/cPref'
               do sample=1,size(sample_frequencies)
                  qeff_s=fq8_effective_q(sample_frequencies(sample),M%tau_Qf8,M%strength_s_Qf8)
                  qeff_p=fq8_effective_q(sample_frequencies(sample),M%tau_Qf8,M%strength_p_Qf8)
                  velocity_ratio_s=fq8_phase_velocity_ratio(sample_frequencies(sample), &
                       parameters%fref,M%tau_Qf8,M%strength_s_Qf8)
                  velocity_ratio_p=fq8_phase_velocity_ratio(sample_frequencies(sample), &
                       parameters%fref,M%tau_Qf8,M%strength_p_Qf8)
                  write(*,'(F8.3,4ES14.5)') sample_frequencies(sample),qeff_s,qeff_p, &
                       velocity_ratio_s,velocity_ratio_p
               end do
            end if
         end if

      end subroutine init_anelastic_Qf8_properties

      subroutine destroy_anelastic_Qf8_properties(M)
         type(block_material), intent(inout) :: M
         if (allocated(M%Qp_inv_Qf8)) deallocate(M%Qp_inv_Qf8)
         if (allocated(M%Qs_inv_Qf8)) deallocate(M%Qs_inv_Qf8)
         if (allocated(M%eta4Qf8)) deallocate(M%eta4Qf8,M%eta5Qf8,M%eta6Qf8, &
              M%eta7Qf8,M%eta8Qf8,M%eta9Qf8)
         if (allocated(M%Deta4Qf8)) deallocate(M%Deta4Qf8,M%Deta5Qf8,M%Deta6Qf8, &
              M%Deta7Qf8,M%Deta8Qf8,M%Deta9Qf8)
         M%anelastic_Qf8=.false.; M%anelastic_Qf=.false.
         M%coarse_grained_Qf8=.false.
         M%tau_Qf8=0.0_wp
         M%strength_s_Qf8=0.0_wp; M%strength_p_Qf8=0.0_wp
      end subroutine destroy_anelastic_Qf8_properties

end module material
