module withers_tables
   !> @brief Withers et al. (2015) Q(f) tables for memory-variable weights
   !> @details Contains Table 1 (high Q) and Table 2 (low Q) from:
   !>          "Memory-Efficient Simulation of Frequency-Dependent Q"
   !>          Bull. Seismol. Soc. Am., Vol. 105, No. 6, pp. 3129-3142
   
   use common, only : wp
   implicit none
   private
   
   ! Public interfaces
   public :: get_withers_weights, get_relaxation_times
   public :: interpolate_gamma, N_GAMMA, N_MECH
   public :: get_withers_weights_Qf, get_relaxation_times_Qf
   
   !> Number of tabulated gamma values (0.0 to 0.9 in steps of 0.1)
   integer, parameter :: N_GAMMA = 10
   
   !> Number of relaxation mechanisms (Withers uses N=8)
   integer, parameter :: N_MECH = 8
   
   !> Gamma values: 0.0, 0.1, 0.2, ..., 0.9
   real(kind=wp), parameter, dimension(N_GAMMA) :: GAMMA_VALUES = &
      [0.0_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp, &
       0.5_wp, 0.6_wp, 0.7_wp, 0.8_wp, 0.9_wp]
   
   !> Minimum relaxation times (seconds) for each gamma
   !> Bandwidth is approximately 0.01-2.5 s (0.4-100 Hz) for gamma < 0.7
   real(kind=wp), parameter, dimension(N_GAMMA) :: TAU_MIN = &
      [0.0032_wp, 0.0032_wp, 0.0032_wp, 0.0032_wp, 0.0032_wp, &
       0.0032_wp, 0.0032_wp, 0.0066_wp, 0.0066_wp, 0.0085_wp]
   
   !> Maximum relaxation times (seconds) for each gamma
   real(kind=wp), parameter, dimension(N_GAMMA) :: TAU_MAX = &
      [15.9155_wp, 15.9155_wp, 15.9155_wp, 15.9155_wp, 15.9155_wp, &
       15.9155_wp, 15.9155_wp,  3.9789_wp,  3.9789_wp,  3.9789_wp]
   
   !> Table 1: Weights for HIGH Q (Q > 200), unscaled (Q* = 1)
   !> To use: w_k(Q) = W_HIGH_Q(k, gamma_index) / Q
   !> Shape: [mechanism 1-8, gamma index 1-10]
   real(kind=wp), parameter, dimension(N_MECH, N_GAMMA) :: W_HIGH_Q = reshape([ &
      ! gamma = 0.0
      0.8867_wp, 0.8323_wp, 0.5615_wp, 0.8110_wp, &
      0.4641_wp, 1.0440_wp, 0.0423_wp, 1.7275_wp, &
      ! gamma = 0.1
      0.3273_wp, 0.8478_wp, 0.3690_wp, 0.9393_wp, &
      0.4474_wp, 1.0434_wp, 0.0440_wp, 1.7268_wp, &
      ! gamma = 0.2
      0.0010_wp, 0.8040_wp, 0.2005_wp, 1.0407_wp, &
      0.4452_wp, 1.0349_wp, 0.0497_wp, 1.7245_wp, &
      ! gamma = 0.3
      0.0010_wp, 0.6143_wp, 0.0918_wp, 1.1003_wp, &
      0.4659_wp, 1.0135_wp, 0.0621_wp, 1.7198_wp, &
      ! gamma = 0.4
      0.0010_wp, 0.4639_wp, 0.0010_wp, 1.1275_wp, &
      0.5090_wp, 0.9782_wp, 0.0820_wp, 1.7122_wp, &
      ! gamma = 0.5
      0.2073_wp, 0.1872_wp, 0.0010_wp, 1.0810_wp, &
      0.6016_wp, 0.9120_wp, 0.1186_wp, 1.6984_wp, &
      ! gamma = 0.6
      0.3112_wp, 0.0010_wp, 0.0010_wp, 1.0117_wp, &
      0.7123_wp, 0.8339_wp, 0.1616_wp, 1.6821_wp, &
      ! gamma = 0.7
      0.1219_wp, 0.0010_wp, 0.0010_wp, 0.2999_wp, &
      1.3635_wp, 0.0010_wp, 0.5084_wp, 1.2197_wp, &
      ! gamma = 0.8
      0.0462_wp, 0.0010_wp, 0.0010_wp, 0.1585_wp, &
      1.4986_wp, 0.0010_wp, 0.4157_wp, 1.3005_wp, &
      ! gamma = 0.9
      0.0010_wp, 0.0010_wp, 0.0010_wp, 0.1935_wp, &
      1.5297_wp, 0.0010_wp, 0.1342_wp, 1.5755_wp &
   ], shape(W_HIGH_Q), order=[1,2])
   
   !> Table 2: Coefficients a_k for LOW Q (15 < Q < 200)
   !> Weight formula: w_k = a_k/Q^2 + b_k/Q
   !> Shape: [mechanism 1-8, gamma index 1-10]
   real(kind=wp), parameter, dimension(N_MECH, N_GAMMA) :: A_COEF = reshape([ &
      ! gamma = 0.0
      -27.50_wp, -34.10_wp,  -1.62_wp, -27.70_wp, &
       14.60_wp, -52.20_wp,  72.00_wp, -82.80_wp, &
      ! gamma = 0.1
        7.37_wp, -37.60_wp,  13.10_wp, -36.10_wp, &
       12.30_wp, -51.40_wp,  69.00_wp, -83.10_wp, &
      ! gamma = 0.2
       31.80_wp, -42.00_wp,  25.70_wp, -40.80_wp, &
        7.02_wp, -49.20_wp,  65.40_wp, -83.20_wp, &
      ! gamma = 0.3
       43.70_wp, -43.40_wp,  34.30_wp, -41.40_wp, &
       -2.87_wp, -45.30_wp,  60.90_wp, -83.10_wp, &
      ! gamma = 0.4
       41.60_wp, -41.10_wp,  38.00_wp, -43.20_wp, &
        5.63_wp, -73.00_wp, 103.00_wp,-164.00_wp, &
      ! gamma = 0.5
       20.00_wp, -23.07_wp,  31.40_wp, -25.10_wp, &
      -45.20_wp, -27.80_wp,  45.90_wp, -81.60_wp, &
      ! gamma = 0.6
        8.08_wp, -13.00_wp,  25.40_wp, -10.40_wp, &
      -75.90_wp, -13.20_wp,  35.70_wp, -79.90_wp, &
      ! gamma = 0.7
        1.99_wp,  -2.70_wp,   0.00_wp,  41.30_wp, &
      -88.80_wp,   0.00_wp,  40.70_wp, -76.60_wp, &
      ! gamma = 0.8
        5.16_wp,  -8.20_wp,   0.00_wp,  58.90_wp, &
     -108.60_wp,  15.02_wp,  -5.88_wp, -46.50_wp, &
      ! gamma = 0.9
       -0.811_wp,  0.00_wp,   0.00_wp,  56.03_wp, &
     -116.90_wp,  22.00_wp,   0.03_wp, -61.90_wp &
   ], shape(A_COEF), order=[1,2])
   
   !> Table 2: Coefficients b_k for LOW Q (15 < Q < 200)
   !> Weight formula: w_k = a_k/Q^2 + b_k/Q
   !> Shape: [mechanism 1-8, gamma index 1-10]
   real(kind=wp), parameter, dimension(N_MECH, N_GAMMA) :: B_COEF = reshape([ &
      ! gamma = 0.0
       7.410_wp,  6.020_wp,  4.680_wp,  6.280_wp, &
       3.880_wp,  8.170_wp,  0.529_wp, 13.190_wp, &
      ! gamma = 0.1
       4.165_wp,  5.520_wp,  3.470_wp,  7.210_wp, &
       3.610_wp,  8.193_wp,  0.498_wp, 13.130_wp, &
      ! gamma = 0.2
       1.612_wp,  5.080_wp,  2.280_wp,  7.931_wp, &
       3.460_wp,  8.150_wp,  0.511_wp, 13.070_wp, &
      ! gamma = 0.3
      -0.109_wp,  4.580_wp,  1.190_wp,  8.390_wp, &
       3.530_wp,  8.020_wp,  0.592_wp, 13.000_wp, &
      ! gamma = 0.4
      -0.734_wp,  3.820_wp,  0.393_wp,  8.670_wp, &
       3.320_wp,  8.580_wp, -0.419_wp, 14.900_wp, &
      ! gamma = 0.5
      -0.435_wp,  2.670_wp, -0.043_wp,  8.245_wp, &
       4.847_wp,  7.190_wp,  1.150_wp, 12.800_wp, &
      ! gamma = 0.6
      -0.196_wp,  1.810_wp, -0.394_wp,  7.657_wp, &
       6.170_wp,  6.360_wp,  1.680_wp, 12.700_wp, &
      ! gamma = 0.7
       0.418_wp,  0.590_wp,  0.000_wp,  2.180_wp, &
      11.000_wp,  0.000_wp,  1.950_wp, 11.300_wp, &
      ! gamma = 0.8
       0.212_wp,  0.345_wp,  0.000_wp,  0.813_wp, &
      12.400_wp, -0.283_wp,  1.420_wp, 11.700_wp, &
      ! gamma = 0.9
       0.162_wp,  0.000_wp,  0.000_wp,  0.797_wp, &
      13.020_wp, -0.402_wp, -0.001_wp, 12.500_wp &
   ], shape(B_COEF), order=[1,2])
   
contains

   !> @brief Get relaxation times for a given gamma
   !> @param[in] gamma Power-law exponent (0.0 to 0.9)
   !> @param[out] tau Relaxation times (seconds), dimension(N_MECH)
   subroutine get_relaxation_times(gamma, tau)
      real(kind=wp), intent(in) :: gamma
      real(kind=wp), dimension(N_MECH), intent(out) :: tau
      
      integer :: k, idx_low, idx_high
      real(kind=wp) :: alpha, taumin, taumax
      
      ! Find gamma indices for interpolation
      call find_gamma_indices(gamma, idx_low, idx_high, alpha)
      
      ! Interpolate tau_min and tau_max
      taumin = (1.0_wp - alpha) * TAU_MIN(idx_low) + alpha * TAU_MIN(idx_high)
      taumax = (1.0_wp - alpha) * TAU_MAX(idx_low) + alpha * TAU_MAX(idx_high)
      
      ! Withers equation (15): logarithmically spaced
      do k = 1, N_MECH
         tau(k) = exp(log(taumin) + real(2*k - 1, wp) / 16.0_wp * &
                      (log(taumax) - log(taumin)))
      end do
      
   end subroutine get_relaxation_times
   
   
   !> @brief Get memory-variable weights for given gamma and Q
   !> @param[in] gamma Power-law exponent (0.0 to 0.9)
   !> @param[in] Q Quality factor
   !> @param[out] weights Memory-variable weights, dimension(N_MECH)
   subroutine get_withers_weights(gamma, Q, weights)
      real(kind=wp), intent(in) :: gamma, Q
      real(kind=wp), dimension(N_MECH), intent(out) :: weights
      
      integer :: k, idx_low, idx_high
      real(kind=wp) :: alpha
      
      ! Find gamma indices for interpolation
      call find_gamma_indices(gamma, idx_low, idx_high, alpha)
      
      if (Q > 200.0_wp) then
         ! High Q: Use Table 1 with linear scaling
         call compute_high_q_weights(idx_low, idx_high, alpha, Q, weights)
      else if (Q >= 15.0_wp) then
         ! Low Q: Use Table 2 interpolation formula
         call compute_low_q_weights(idx_low, idx_high, alpha, Q, weights)
      else
         ! Q < 15: Use Q=15 weights (lower bound of validity)
         call compute_low_q_weights(idx_low, idx_high, alpha, 15.0_wp, weights)
      end if
      
   end subroutine get_withers_weights
   
   
   !> @brief Find gamma indices and interpolation weight
   !> @param[in] gamma Target gamma value
   !> @param[out] idx_low Lower index in GAMMA_VALUES
   !> @param[out] idx_high Upper index in GAMMA_VALUES
   !> @param[out] alpha Interpolation weight (0 to 1)
   subroutine find_gamma_indices(gamma, idx_low, idx_high, alpha)
      real(kind=wp), intent(in) :: gamma
      integer, intent(out) :: idx_low, idx_high
      real(kind=wp), intent(out) :: alpha
      
      real(kind=wp) :: gamma_clamped
      integer :: i
      
      ! Clamp gamma to valid range [0.0, 0.9]
      gamma_clamped = max(0.0_wp, min(0.9_wp, gamma))
      
      ! Find bracketing indices
      idx_low = 1
      do i = 1, N_GAMMA - 1
         if (gamma_clamped >= GAMMA_VALUES(i) .and. &
             gamma_clamped <= GAMMA_VALUES(i+1)) then
            idx_low = i
            idx_high = i + 1
            exit
         end if
      end do
      
      ! Handle edge cases
      if (gamma_clamped <= GAMMA_VALUES(1)) then
         idx_low = 1
         idx_high = 1
         alpha = 0.0_wp
      else if (gamma_clamped >= GAMMA_VALUES(N_GAMMA)) then
         idx_low = N_GAMMA
         idx_high = N_GAMMA
         alpha = 0.0_wp
      else
         ! Linear interpolation weight
         alpha = (gamma_clamped - GAMMA_VALUES(idx_low)) / &
                 (GAMMA_VALUES(idx_high) - GAMMA_VALUES(idx_low))
      end if
      
   end subroutine find_gamma_indices
   
   
   !> @brief Compute weights for high Q using Table 1
   subroutine compute_high_q_weights(idx_low, idx_high, alpha, Q, weights)
      integer, intent(in) :: idx_low, idx_high
      real(kind=wp), intent(in) :: alpha, Q
      real(kind=wp), dimension(N_MECH), intent(out) :: weights
      
      integer :: k
      real(kind=wp) :: w_low, w_high
      
      do k = 1, N_MECH
         w_low = W_HIGH_Q(k, idx_low)
         w_high = W_HIGH_Q(k, idx_high)
         
         ! Linear interpolation in gamma
         weights(k) = (1.0_wp - alpha) * w_low + alpha * w_high
         
         ! Scale by 1/Q (Withers scaling property)
         weights(k) = weights(k) / Q
      end do
      
   end subroutine compute_high_q_weights
   
   
   !> @brief Compute weights for low Q using Table 2 formula
   subroutine compute_low_q_weights(idx_low, idx_high, alpha, Q, weights)
      integer, intent(in) :: idx_low, idx_high
      real(kind=wp), intent(in) :: alpha, Q
      real(kind=wp), dimension(N_MECH), intent(out) :: weights
      
      integer :: k
      real(kind=wp) :: ak_low, ak_high, bk_low, bk_high, ak, bk
      
      do k = 1, N_MECH
         ! Interpolate a_k coefficients
         ak_low = A_COEF(k, idx_low)
         ak_high = A_COEF(k, idx_high)
         ak = (1.0_wp - alpha) * ak_low + alpha * ak_high
         
         ! Interpolate b_k coefficients
         bk_low = B_COEF(k, idx_low)
         bk_high = B_COEF(k, idx_high)
         bk = (1.0_wp - alpha) * bk_low + alpha * bk_high
         
         ! Withers equation (16): w_k = a_k/Q^2 + b_k/Q
         weights(k) = ak / (Q**2) + bk / Q
      end do
      
   end subroutine compute_low_q_weights
   
   
   !> @brief Interpolate gamma value from spatial arrays (utility function)
   !> @param[in] gamma_array Spatial array of gamma values
   !> @param[in] i,j,k Grid indices
   !> @param[out] gamma Interpolated gamma value
   subroutine interpolate_gamma(gamma_array, i, j, k, gamma)
      real(kind=wp), dimension(:,:,:), intent(in) :: gamma_array
      integer, intent(in) :: i, j, k
      real(kind=wp), intent(out) :: gamma
      
      ! Simple: just take the value at this point
      ! Could be extended to do spatial interpolation if needed
      gamma = gamma_array(i, j, k)
      
      ! Clamp to valid range
      gamma = max(0.0_wp, min(0.9_wp, gamma))
      
   end subroutine interpolate_gamma


   !> @brief Get relaxation times for 4-mechanism Qf variant
   !> @details Extracts first 4 mechanisms from Withers Table 1 (N_MECH=8)
   !> @param[in] gamma Power-law exponent (0.0 to 0.9)
   !> @param[out] tau_Qf Relaxation times (seconds), dimension(4)
   subroutine get_relaxation_times_Qf(gamma, tau_Qf)
      real(kind=wp), intent(in) :: gamma
      real(kind=wp), dimension(4), intent(out) :: tau_Qf
      
      real(kind=wp), dimension(N_MECH) :: tau_full
      
      ! Get full 8-mechanism relaxation times
      call get_relaxation_times(gamma, tau_full)
      
      ! Extract first 4 mechanisms
      tau_Qf(1:4) = tau_full(1:4)
      
   end subroutine get_relaxation_times_Qf


   !> @brief Get memory-variable weights for 4-mechanism Qf variant
   !> @details Extracts first 4 mechanisms from Withers Tables 1 & 2 (N_MECH=8)
   !> @param[in] gamma Power-law exponent (0.0 to 0.9)
   !> @param[in] Q Quality factor
   !> @param[out] weights_Qf Memory-variable weights, dimension(4)
   subroutine get_withers_weights_Qf(gamma, Q, weights_Qf)
      real(kind=wp), intent(in) :: gamma, Q
      real(kind=wp), dimension(4), intent(out) :: weights_Qf
      
      real(kind=wp), dimension(N_MECH) :: weights_full
      
      ! Get full 8-mechanism weights
      call get_withers_weights(gamma, Q, weights_full)
      
      ! Extract first 4 mechanisms (Note: this uses the first 4 from the 8-mechanism table)
      ! For anelastic-Qf (m4), we use mechanisms 1-4 which are optimized for that case
      weights_Qf(1:4) = weights_full(1:4)
      
   end subroutine get_withers_weights_Qf

end module withers_tables
