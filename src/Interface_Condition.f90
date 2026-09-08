module Interface_Condition

  use common, only : wp
  use  initial_stress_condition, only : initial_stress_tensor
  implicit none

contains

  !> This subroutine solves for the hat variables u_out and v_out using the grid variables u_in and v_in
  !> All variables are in the locally rotated othorgonal cordinates

  subroutine Interface_RHS(problem, coupling, t, stage, handles, u_in, v_in, &
                       u_out, v_out, G, mat_block, LambdaMu1, LambdaMu2, n_l, n_m, n_n, I, ib)

    use mpi3dbasic, only : rank
    use datatypes, only : block_grid_t, block_boundary, iface_type, fault_type, block_material
    implicit none

    character(256), intent(in) :: problem
    real(kind = wp), intent(in) :: t                                        !< time
    character(*),intent(in) :: coupling
    integer, intent(in) :: stage, ib
    type(fault_type), intent(inout) :: handles
    real(kind = wp), dimension(:,:,:), allocatable, intent(in) :: u_in, v_in           !< grid variables
    real(kind = wp), dimension(:,:,:), allocatable, intent(inout) :: u_out, v_out      !< hat variables
    real(kind = wp), dimension(:,:,:), allocatable, intent(in) :: LambdaMu1, LambdaMu2           !< material properties
    real(kind = wp), dimension(:,:,:), allocatable, intent(in) :: n_l, n_m, n_n          !< normals

    type(block_grid_t), intent(in) :: G
    type(iface_type), intent(inout) :: I
    type(block_material), intent(in) :: mat_block


    integer :: mx1, my1, mz1, px1, py1, pz1 !< grid bounds

    ! all variables below are point-wise on the fault plane
    real(kind = wp) :: c1_p, c1_s, c2_p, c2_s, rho1, rho2                   !< local velocities and densities
    real(kind = wp) :: lambda1, lambda2, mu1, mu2                           !< local lambda and mu
    real(kind = wp) :: Z1_s, Z2_s, Z1_p, Z2_p                               !< local p-wave impedance and s-wave impedance
    real(kind = wp) :: u_n, u_m, u_l, Tu_n, Tu_m, Tu_l                      !< local velocities and stresses for the left block
    real(kind = wp) :: v_n, v_m, v_l, Tv_n, Tv_m, Tv_l                      !< local velocities and stresses for the right block

    integer :: y, z, x                                             !< local grid indices
    real(kind = wp) :: alpha                                                !< friction coeff. for the linearized friction law
    real(kind = wp) :: slip, slip_rate                                      !< local slip and slip rate
    real(kind = wp), dimension(1:6) ::  U_loc, V_loc, P_loc, Q_loc          !< work arrays for local grid functions and hat-functions
    real(kind = wp), dimension(1:3) :: x_ij                                 !< local grid points on the fault plane
    real(kind = wp) :: T0_n, T0_m, T0_l                                     !< tractions acting at point on the fault plane
    real(kind = wp), dimension(1:3) :: l, m, n                              !< local unit vectors
    real(kind = wp) :: psi                                                  !< state variable
    real(kind = wp) ::  a, V0, Fss, L0, b, f0                               !< friction parameters
    real(kind = wp) :: Vs, rho, u, Mat(3)                                   !< shear modulus
    real(kind = wp) :: Vel, V_w

    mx1 = G%C%mq
    px1 = G%C%pq
    my1 = G%C%mr
    py1 = G%C%pr
    mz1 = G%C%ms
    pz1 = G%C%ps

    do z = mz1, pz1
       do y = my1, py1

          !==================================================================================
          ! Set local material parameters and compute local velocities
          !==================================================================================

          lambda1 = LambdaMu1(y, z, 1)                                                    ! Lambda1
          mu1 = LambdaMu1(y, z, 2)                                                        ! Mu1
          rho1 = LambdaMu1(y, z, 3)                                                       ! density1
          c1_p = sqrt((2.0_wp*mu1 + lambda1)/rho1)                                                 ! P-wave speed1
          c1_s = sqrt(mu1/rho1)                                                                 ! S-wave speed1

          lambda2 = LambdaMu2(y, z, 1)                                                     ! Lambda2
          mu2 = LambdaMu2(y, z, 2)                                                         ! Mu2
          rho2 = LambdaMu2(y, z, 3)                                                        ! density2
          c2_p = sqrt((2.0_wp*mu2 + lambda2)/rho2)                                                 ! P-wave speed2
          c2_s = sqrt(mu2/rho2)                                                                 ! S-wave speed2

          if (ib == 1) then
             x_ij(:) = G%X(mx1, y, z, 1:3)

              !if(my1 == 1) x_ij(2) = x_ij(2)- G%X(mx1, 1, z, 2)
             
             !print*, 'x = ',  x_ij(1), 'y = ',  x_ij(2), 'z = ',  x_ij(3)
             !print*, 'x = ', G%X(mx1, y, z, 1), 'y = ', G%X(mx1, y, z, 2), 'z = ', G%X(mx1, y, z, 3)
          else
             x_ij(:) = G%X(px1, y, z, 1:3)

             
             !if(my1 == 1) x_ij(2) = x_ij(2) - G%X(px1, 1, z, 2)

             
             !print*, 'x = ',  x_ij(1), 'y = ',  x_ij(2), 'z = ',  x_ij(3)
          end if


         !if (z == 1) print*,  x_ij(2)

          select case(problem)

          case default

              Mat(2) = 1.0e10_wp

          case('LOH1')

             if (G%x(mx1,y,z,1) < 0.9999_wp) then

              Vs  = 2_wp ! S-wave speed
              Mat(3)  = 2.6_wp ! rho (density)       

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)

             else if (G%x(mx1,y,z,1) >= 0.9999_wp .AND. G%x(mx1,y,z,1) <= 1.0001_wp) then

              Vs  = 2.732_wp ! S-wave speed
              Mat(3)  = 2.65_wp ! rho (density)       

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)  

             else if (G%x(mx1,y,z,1) > 1.0001_wp) then

               Vs = 3.464_wp    ! S-wave speed 
               Mat(3) = 2.7_wp    ! rho (density)

               Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             end if


          case('TPV31')

             if (G%x(mx1,y,z,2) < 2.3999_wp) then

              Vs  = 2.25_wp ! S-wave speed
              Mat(3)  = 2.58_wp ! rho (density)       

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)

             else if (G%x(mx1,y,z,2) >= 2.3999_wp .AND. G%x(mx1,y,z,2) <= 2.4001_wp) then

              Vs  = 2.4_wp ! S-wave speed
              Mat(3)  = 2.59_wp ! rho (density)       

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 2.4001_wp .AND. G%x(mx1,y,z,2) < 4.9999_wp) then
 
               Vs  = 2.55_wp + 0.5_wp*(G%x(mx1,y,z,2)-2.4_wp)/2.6_wp ! S-wave speed 
               Mat(3)  = 2.6_wp + 0.02_wp*(G%x(mx1,y,z,2)-2.4_wp)/2.6_wp  ! rho (density) 

               Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)

             else if (G%x(mx1,y,z,2) >= 4.9999_wp .AND. G%x(mx1,y,z,2) <= 5.0001_wp) then

              Vs  = 3.25_wp ! S-wave speed
              Mat(3)  = 2.67_wp ! rho (density)       

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)

             else if (G%x(mx1,y,z,2) > 5.0001_wp .AND. G%x(mx1,y,z,2) < 9.9999_wp) then

               Vs = 3.45_wp   ! S-wave speed 
               Mat(3) = 2.72_wp   ! rho (density)

               Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)

             else if (G%x(mx1,y,z,2) >= 9.9999_wp .AND. G%x(mx1,y,z,2) <= 10.0001_wp) then

              Vs  = 3.625_wp ! S-wave speed
              Mat(3)  = 2.76_wp ! rho (density)       

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus)

             else if (G%x(mx1,y,z,2) > 10.0001_wp) then

               Vs = 3.8_wp    ! S-wave speed 
               Mat(3) = 3.0_wp    ! rho (density)

               Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             end if  


         case('TPV32')

             if (G%x(mx1,y,z,2) <= 0.5_wp) then

              Vs = 1.05_wp + 0.35_wp*(G%x(mx1,y,z,2))/0.5_wp    ! S-wave speed 
              Mat(3) = 2.2_wp + 0.25_wp*(G%x(mx1,y,z,2))/0.5_wp     ! rho (density)   

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 0.5_wp .AND. G%x(mx1,y,z,2) <= 1.0_wp) then
 
              Vs = 1.4_wp + 0.55_wp*(G%x(mx1,y,z,2)-0.5_wp)/0.5_wp   ! S-wave speed 
              Mat(3) = 2.45_wp + 0.1_wp*(G%x(mx1,y,z,2)-0.5_wp)/0.5_wp   ! rho (density) 

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 1.0_wp .AND. G%x(mx1,y,z,2) <= 1.6_wp) then

              Vs = 1.95_wp + 0.55_wp*(G%x(mx1,y,z,2)-1.0_wp)/0.6_wp   ! S-wave speed 
              Mat(3) = 2.55_wp + 0.05_wp*(G%x(mx1,y,z,2)-1.0_wp)/0.6_wp   ! rho (density)               

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 1.6_wp .AND. G%x(mx1,y,z,2) <= 2.4_wp) then

              Vs = 2.5_wp + 0.3_wp*(G%x(mx1,y,z,2)-1.6_wp)/0.8_wp   ! S-wave speed 
              Mat(3) = 2.6_wp                                         ! rho (density)

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 
               
             else if (G%x(mx1,y,z,2) > 2.4_wp .AND. G%x(mx1,y,z,2) <= 3.6_wp) then

              Vs = 2.8_wp + 0.3_wp*(G%x(mx1,y,z,2)-2.4_wp)/1.2_wp   ! S-wave speed 
              Mat(3) = 2.6_wp + 0.02_wp*(G%x(mx1,y,z,2)-2.4_wp)/1.2_wp  ! rho (density) 

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 3.6_wp .AND. G%x(mx1,y,z,2) <= 5.0_wp) then

              Vs = 3.1_wp + 0.15_wp*(G%x(mx1,y,z,2)-3.6_wp)/1.4_wp   ! S-wave speed 
              Mat(3) = 2.62_wp + 0.03_wp*(G%x(mx1,y,z,2)-3.6_wp)/1.4_wp  ! rho (density)

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 5.0_wp .AND. G%x(mx1,y,z,2) <= 9.0_wp) then

              Vs = 3.25_wp + 0.2_wp*(G%x(mx1,y,z,2)-5.0_wp)/4.0_wp   ! S-wave speed 
              Mat(3) = 2.65_wp + 0.07_wp*(G%x(mx1,y,z,2)-5.0_wp)/4.0_wp  ! rho (density)

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 9.0_wp .AND. G%x(mx1,y,z,2) <= 11.0_wp) then

              Vs = 3.45_wp + 0.25_wp*(G%x(mx1,y,z,2)-9.0_wp)/2.0_wp   ! S-wave speed 
              Mat(3) = 2.72_wp + 0.03_wp*(G%x(mx1,y,z,2)-9.0_wp)/2.0_wp   ! rho (density) 

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 11.0_wp .AND. G%x(mx1,y,z,2) <= 15.0_wp) then

              Vs = 3.6_wp + 0.1_wp*(G%x(mx1,y,z,2)-11.0_wp)/4.0_wp     ! S-wave speed 
              Mat(3) = 2.75_wp + 0.15_wp*(G%x(mx1,y,z,2)-11.0_wp)/4.0_wp   ! rho (density)

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             else if (G%x(mx1,y,z,2) > 15.0_wp) then

              Vs = 3.7_wp     ! S-wave speed 
              Mat(3) = 2.9_wp     ! rho (density) 

              Mat(2) = Mat(3)*Vs**2 ! mu (shear modulus) 

             end if 

         case('TPV34')

            if (G%x(mx1,y,z,1) > -0.01_wp .AND. G%x(mx1,y,z,1) < 0.01_wp) then
              Mat(1) = mat_block%M(mx1,y,z,1)
              Mat(2) = mat_block%M(mx1,y,z,2)
              Mat(3) = mat_block%M(mx1,y,z,3)

            else                   
              Mat(1) = mat_block%M(px1,y,z,1)
              Mat(2) = mat_block%M(px1,y,z,2)
              Mat(3) = mat_block%M(px1,y,z,3) 

            end if

         case('SCITS2016')

            if (G%x(mx1,y,z,1) > 19.99_wp .AND. G%x(mx1,y,z,1) < 20.01_wp) then
              Mat(1) = mat_block%M(mx1,y,z,1)
              Mat(2) = mat_block%M(mx1,y,z,2)
              Mat(3) = mat_block%M(mx1,y,z,3)

            else                   
              Mat(1) = mat_block%M(px1,y,z,1)
              Mat(2) = mat_block%M(px1,y,z,2)
              Mat(3) = mat_block%M(px1,y,z,3) 

            end if  


         end select


          n(:) = n_n(y,z,:)
          m(:) = n_m(y,z,:)
          l(:) = n_l(y,z,:)

          !==================================================================================
          ! Set local grid-functions
          !==================================================================================

          ! velocities in the left elastic block
          u_n = u_in(y, z, 1)
          u_m = u_in(y, z, 2)
          u_l = u_in(y, z, 3)

          ! velocities in the right elastic block
          v_n = v_in(y, z, 1)
          v_m = v_in(y, z, 2)
          v_l = v_in(y, z, 3)

          ! stress fields in the left elastic block
          Tu_n = u_in(y, z, 4)
          Tu_m = u_in(y, z, 5)
          Tu_l = u_in(y, z, 6)
          ! stress fields in the right elastic block
          Tv_n = v_in(y, z, 4)
          Tv_m = v_in(y, z, 5)
          Tv_l = v_in(y, z, 6)

          ! initialize local work array, U_loc = (/ u1,u2,u3, T1_nn,.../)

          U_loc(1) = u_n
          U_loc(2) = u_m
          U_loc(3) = u_l

          U_loc(4) = Tu_n
          U_loc(5) = Tu_m
          U_loc(6) = Tu_l

          V_loc(1) = v_n
          V_loc(2) = v_m
          V_loc(3) = v_l

          V_loc(4) = Tv_n
          V_loc(5) = Tv_m
          V_loc(6) = Tv_l

          ! s-wave impedance
          Z1_s = rho1*c1_s
          Z2_s = rho2*c2_s

          ! p-wave impedance
          Z1_p = rho1*c1_p
          Z2_p = rho2*c2_p

          ! ==================================================================================
          ! COMPUTE HAT-VARIABLES USING INTERFACE CONDITIONS
          ! ==================================================================================

          select case(coupling)

          case default

             stop 'invalid coupling condition'

          case('locked') ! locked or welded interface

             call Welded_Interface(U_loc,V_loc, P_loc, Q_loc, Z1_s, Z2_s, Z1_p, Z2_p)

          case ('linear_friction') ! linear friction law: tau = alpha*v
             ! There are two limiting values of alpha:
             ! Frictionless interface: alpha = 0
             ! Locked interface: alpha = very large

             alpha = 1.0e14_wp
             call Linear_Friction(U_loc,V_loc, P_loc, Q_loc, Z1_s, Z2_s, Z1_p, Z2_p, alpha)

          case('slip-weakening_friction') ! slip-weakning friction law

             slip = I%S(y,z,4)

             call Slip_Weakening_Friction(problem,U_loc, V_loc, P_loc, Q_loc, Vel, Z1_s, Z2_s, Z1_p, Z2_p, &
                  slip, x_ij, t, l, m, n,coupling,Mat)

          case('rate-and-state_friction')
             
             ! set friction parameters
             ! L0 = 0.2572_wp
             ! f0 = 0.7_wp
             ! b =  0.02_wp
             ! V0 = 1.0e-6_wp
             ! a = 0.016_wp
             
             select case(problem)
                
             case('TPV101')
                L0 = 0.02_wp
                f0 = 0.6_wp
                b = 0.012_wp
                V0 = 1.0e-6_wp
                a = 0.008_wp+0.008_wp*(1.0_wp-smooth_boxcar(x_ij(2),15.0_wp,3.0_wp)&
                     *smooth_boxcar(x_ij(3),15.0_wp,3.0_wp))
                
             case('TPV102')
                L0 = 0.02_wp
                f0 = 0.6_wp
                b = 0.012_wp
                V0 = 1.0e-6_wp                
                a = 0.008_wp+0.008_wp*(1.0_wp-smooth_boxcar(x_ij(2),15.0_wp,3.0_wp)&
                     *smooth_boxcar(x_ij(3),15.0_wp,3.0_wp))
                
             case('KDN')
                L0 = 0.4_wp
                f0 = 0.3_wp
                b =  0.012_wp
                V0 = 1.0e-6_wp
                a = 0.008_wp+0.008_wp*(1.0_wp-smooth_boxcar(x_ij(2),15.0_wp,3.0_wp)&
                     *smooth_boxcar(x_ij(3),15.0_wp,3.0_wp))
                
             case('rate-weakening')
                L0 = 2.0_wp*0.2572_wp
                f0 = 0.7_wp
                b =  0.02_wp
                V0 = 1.0e-6_wp
                a = 0.016_wp+0.016_wp*(1.0_wp-smooth_boxcar(x_ij(2),20.0_wp,2.5_wp)&
                     *smooth_boxcar(x_ij(3),20.0_wp,2.5_wp))
             end select
             
             
             psi = I%W(y,z,1)  ! state variable 
             
             ! compute hat-variables from the friction law
             call rate_state_friction(problem,U_loc,V_loc, P_loc,Q_loc, Vel, Z1_s, Z2_s, Z1_p, Z2_p, &
                  psi,a,V0, x_ij,t, l,m,n, coupling)
             
             slip_rate = sqrt((Q_loc(2)-P_loc(2))**2 + (Q_loc(3)-P_loc(3))**2) ! hat variable slip velocity
                          
             ! set the right-handside of the rate of state to be used in the Runge-Kutta integration for state
             select case(problem)
                
             case('TPV101')
                !I%DW(y,z,1) = I%DW(y,z,1) + 1.0_wp-(1.0_wp/L0)*slip_rate*I%W(y,z,1)
                I%DW(y,z,1) = I%DW(y,z,1) + b*V0/L0*exp(-(I%W(y,z,1)-f0)/b) - Vel*b/L0
             case('TPV102')
                !I%DW(y,z,1) = I%DW(y,z,1) + 1.0_wp-(1.0_wp/L0)*slip_rate*I%W(y,z,1)
                I%DW(y,z,1) = I%DW(y,z,1) + b*V0/L0*exp(-(I%W(y,z,1)-f0)/b) - Vel*b/L0
                             
             case('rate-weakening')
                V_w  = 0.17_wp + 5.0_wp*(1.0_wp-smooth_boxcar(x_ij(2),20.0_wp,2.5_wp)&
                     *smooth_boxcar(x_ij(3),20.0_wp,2.5_wp))
                Fss = f0 -(b-a)*log(Vel/V0)
                Fss = 0.13_wp + (Fss - 0.13_wp)/((1.0_wp + (Vel/V_w)**8)**0.125_wp)
                I%DW(y,z,1) = I%DW(y,z,1) - &
                     Vel/L0*(a*asinh(Vel/(2.0_wp*V0)*exp(I%W(y,z,1)/a))-Fss)
                
             case('KDN')
                Fss = f0 -(b-a)*log(Vel/V0)
                I%DW(y,z,1) = I%DW(y,z,1) - &
                     (1.0_wp/L0)*slip_rate*(asinh(Vel/(2.0_wp*V0)*exp(I%W(y,z,1)/a))-Fss)
                
             end select
          end select

          u_out(y, z, 1) = P_loc(1)
          u_out(y, z, 2) = P_loc(2)
          u_out(y, z, 3) = P_loc(3)

          v_out(y, z, 1) = Q_loc(1)
          v_out(y, z, 2) = Q_loc(2)
          v_out(y, z, 3) = Q_loc(3)

          u_out(y, z, 4) = P_loc(4)
          u_out(y, z, 5) = P_loc(5)
          u_out(y, z, 6) = P_loc(6)

          v_out(y, z, 4) = Q_loc(4)
          v_out(y, z, 5) = Q_loc(5)
          v_out(y, z, 6) = Q_loc(6)

          slip_rate = sqrt((Q_loc(2)-P_loc(2))**2 + (Q_loc(3)-P_loc(3))**2) ! hat variable slip velocity?
          Vel = slip_rate

          ! used in Runge-Kutta time integration for total slip and components of slip

          I%DS(y,z,1:3) = I%DS(y,z,1:3) + Q_loc(1:3) - P_loc(1:3)
          I%DS(y,z,4) = I%DS(y,z,4) + Vel

          ! Update hat-variables on the fault for output

          if (stage == 1) then


             if (coupling == 'locked') then

                T0_n = 0d0
                T0_m = 0d0
                T0_l = 0d0

             else
                call prestress(T0_n, T0_m, T0_l, x_ij, t, problem, l, m, n, Mat)

             end if

             handles%Uhat_pluspres(y,z,1:3) = P_loc(1:3)
             handles%Vhat_pluspres(y,z,1:3) = Q_loc(1:3)

             handles%Uhat_pluspres(y,z,4) = P_loc(4) + T0_n
             handles%Vhat_pluspres(y,z,4) = Q_loc(4) + T0_n

             handles%Uhat_pluspres(y,z,5) = P_loc(5) + T0_m
             handles%Vhat_pluspres(y,z,5) = Q_loc(5) + T0_m

             handles%Uhat_pluspres(y,z,6) = P_loc(6) + T0_l
             handles%Vhat_pluspres(y,z,6) = Q_loc(6) + T0_l


             

          end if
       end do
    end do

  end subroutine Interface_RHS


  subroutine Welded_Interface(u_in,v_in, u_out,v_out, Z1_s, Z2_s, Z1_p, Z2_p)

    ! This subroutine solves for the hat-variables governed by a welded interface
    ! All variables are in the locally rotated orthogonal coordinates

    implicit none
    real(kind = wp), dimension(:), intent(in) :: u_in, v_in
    real(kind = wp), dimension(:), intent(out) :: u_out, v_out
    real(kind = wp), intent(in) :: Z1_s, Z2_s, Z1_p, Z2_p                      ! Local velocities and densities

    real(kind = wp) :: u_n, u_m, u_l, Tu_n, Tu_m, Tu_l                                 ! Local velocities and stresses for the left block
    real(kind = wp) :: v_n, v_m, v_l, Tv_n, Tv_m, Tv_l                                 ! Local velocities and stresses for the right block
    real(kind = wp) :: u_n_hat, u_m_hat, u_l_hat, &
            Tu_n_hat, Tu_m_hat, Tu_l_hat ! Local velocities_hat and stresses_hat for the left block
    real(kind = wp) :: v_n_hat, v_m_hat, v_l_hat, &
         Tv_n_hat, Tv_m_hat, Tv_l_hat ! Local velocities_hat and stresses_hat for the right block

    real(kind = wp) :: p_n, p_m, p_l, q_n, q_m, q_l
    real(kind = wp) :: phi_n, phi_m, phi_l,  eta_s,  eta_p

    ! velocities in the left elastic block
    u_n = u_in(1)
    u_m = u_in(2)
    u_l = u_in(3)

    ! velocities in the right elastic block
    v_n = v_in(1)
    v_m = v_in(2)
    v_l = v_in(3)

    ! stress fields in the left elastic block
    Tu_n = u_in(4)
    Tu_m = u_in(5)
    Tu_l = u_in(6)

    ! stress fields in the right elastic block
    Tv_n = v_in(4)
    Tv_m = v_in(5)
    Tv_l = v_in(6)


     ! compute characteristics
    p_n = Z1_p*u_n - Tu_n
    p_m = Z1_s*u_m - Tu_m
    p_l = Z1_s*u_l - Tu_l

    q_n = Z2_p*v_n + Tv_n
    q_m = Z2_s*v_m + Tv_m
    q_l = Z2_s*v_l + Tv_l
    
    eta_p = 1.0_wp/(1.0_wp/Z1_p + 1.0_wp/Z2_p)                                      ! half of the harmonic mean of Z1_p, Z2_p
    eta_s = 1.0_wp/(1.0_wp/Z1_s + 1.0_wp/Z2_s)                                      ! half of the harmonic mean of Z1_s, Z2_s
    
    ! stress transfer function
    phi_n = eta_p*(q_n/Z2_p - p_n/Z1_p)
    phi_m = eta_s*(q_m/Z2_s - p_m/Z1_s)
    phi_l = eta_s*(q_l/Z2_s - p_l/Z1_s)

    ! Compute Hat-Variables
    u_n_hat = (q_n - phi_n)/Z2_p                                        ! continuity of normal velocity (vn-un = 0, no opening)
    Tu_n_hat = phi_n                                                    ! continuity of normal stress

    v_n_hat = (p_n + phi_n)/Z1_p                                        ! continuity of normal velocity (vn-un = 0, no opening)
    Tv_n_hat = phi_n                                                    ! continuity of normal stress

     ! Compute Hat-Variables
    u_l_hat = (q_l - phi_l)/Z2_s                                        ! continuity of normal velocity (vl-ul = 0, no slip)
    Tu_l_hat = phi_l                                                    ! continuity of shear stress

    v_l_hat = (p_l + phi_l)/Z1_s                                        ! continuity of normal velocity (vl-ul = 0, no slip)
    Tv_l_hat = phi_l                                                    ! continuity of shear stress


    ! Compute Hat-Variables
    u_m_hat = (q_m - phi_m)/Z2_s                                        ! continuity of normal velocity (vm-um = 0, no slip)
    Tu_m_hat = phi_m                                                    ! continuity of shear stress
    
    v_m_hat = (p_m + phi_m)/Z1_s                                        ! continuity of normal velocity (vm-um = 0, no slip)
    Tv_m_hat = phi_m                                                    ! continuity of shear stress
    

!!$    ! Continuity of normal traction and normal velocities
!!$    u_n_hat = v_n
!!$    Tu_n_hat = Tv_n
!!$
!!$    v_n_hat = u_n
!!$    Tv_n_hat = Tu_n
!!$
!!$    ! Continuity of shear traction and shear velocities
!!$
!!$    u_m_hat = v_m
!!$    Tu_m_hat = Tv_m
!!$
!!$    v_m_hat = u_m
!!$    Tv_m_hat = Tu_m
!!$
!!$    u_l_hat = v_l
!!$    Tu_l_hat = Tv_l
!!$
!!$    v_l_hat = u_l
!!$    Tv_l_hat = Tu_l
    
    ! Save all hat-variables in the work-array  u_out, v_out

    u_out(1) = u_n_hat
    u_out(2) = u_m_hat
    u_out(3) = u_l_hat

    u_out(4) = Tu_n_hat
    u_out(5) = Tu_m_hat
    u_out(6) = Tu_l_hat

    v_out(1) = v_n_hat
    v_out(2) = v_m_hat
    v_out(3) = v_l_hat

    v_out(4) = Tv_n_hat
    v_out(5) = Tv_m_hat
    v_out(6) = Tv_l_hat

    !==================================================================================
  end subroutine Welded_Interface

  subroutine Linear_Friction(u_in,v_in, u_out,v_out, Z1_s, Z2_s, Z1_p, Z2_p, alpha)

    ! This subroutine solves for the hat-variables governed by a linear friction law: tau = alpha*v
    ! All variables are in the locally rotated othorgonal cordinates

    implicit none

    real(kind = wp), dimension(:), intent(in) :: u_in, v_in
    real(kind = wp), dimension(:), intent(out) :: u_out, v_out
    real(kind = wp), intent(in) :: Z1_s, Z2_s, Z1_p, Z2_p, alpha                       ! Local velocities and densities

    real(kind = wp) :: u_n, u_m, u_l, Tu_n, Tu_m, Tu_l   ! Local velocities and stresses for the left block
    real(kind = wp) :: v_n, v_m, v_l, Tv_n, Tv_m, Tv_l   ! Local velocities and stresses for the right block
    real(kind = wp) :: u_n_hat, u_m_hat, u_l_hat, &
            Tu_n_hat, Tu_m_hat, Tu_l_hat ! Local velocities_hat and stresses_hat for the left block
    real(kind = wp) :: v_n_hat, v_m_hat, v_l_hat, &
            Tv_n_hat, Tv_m_hat, Tv_l_hat ! Local velocities_hat and stresses_hat for the right block
    real(kind = wp) :: p_n, p_m, p_l, q_n, q_m, q_l
    real(kind = wp) :: phi_n, phi_m, phi_l,  eta_s,  eta_p

    if((Z1_s .le. 1.0e-9_wp) .or. (Z2_s .le. 1.0e-9_wp)) STOP 'shear impedance MUST be NONZERO'

    ! velocities in the left elastic block
    u_n = u_in(1)
    u_m = u_in(2)
    u_l = u_in(3)

    ! velocities in the right elastic block
    v_n = v_in(1)
    v_m = v_in(2)
    v_l = v_in(3)

    ! stress fields in the left elastic block
    Tu_n = u_in(4)
    Tu_m = u_in(5)
    Tu_l = u_in(6)

    ! stress fields in the right elastic block
    Tv_n = v_in(4)
    Tv_m = v_in(5)
    Tv_l = v_in(6)

    ! compute characteristics
    p_n = Z1_p*u_n - Tu_n
    p_m = Z1_s*u_m - Tu_m
    p_l = Z1_s*u_l - Tu_l

    q_n = Z2_p*v_n + Tv_n
    q_m = Z2_s*v_m + Tv_m
    q_l = Z2_s*v_l + Tv_l
    
    eta_p = 1.0_wp/(1.0_wp/Z1_p + 1.0_wp/Z2_p)                                      ! half of the harmonic mean of Z1_p, Z2_p
    eta_s = 1.0_wp/(1.0_wp/Z1_s + 1.0_wp/Z2_s)                                      ! half of the harmonic mean of Z1_s, Z2_s
    
    ! stress transfer function
    phi_n = eta_p*(q_n/Z2_p - p_n/Z1_p)
    phi_m = eta_s*(q_m/Z2_s - p_m/Z1_s)
    phi_l = eta_s*(q_l/Z2_s - p_l/Z1_s)

    ! Compute Hat-Variables
    u_n_hat = (q_n - phi_n)/Z2_p                                        ! continuity of normal velocity (v1-u1 = 0, no opening)
    Tu_n_hat = phi_n                                                    ! continuity of normal stress

    v_n_hat = (p_n + phi_n)/Z1_p                                        ! continuity of normal velocity (v1-u1 = 0, no opening)
    Tv_n_hat = phi_n                                                    ! continuity of normal stress

    ! Solve the equations: tau1 - alpha*v1 = 0; tau1 + eta*v1 = phi_1; tau2 - alpha*v2 = 0; tau2 + eta*v2 = phi_2

    Tu_m_hat = alpha/(alpha + eta_s)*phi_m                                ! = tau1
    Tu_l_hat = alpha/(alpha + eta_s)*phi_l                                ! = tau2
    
    Tv_m_hat = alpha/(alpha + eta_s)*phi_m                                ! = tau1
    Tv_l_hat = alpha/(alpha + eta_s)*phi_l                                ! = tau2
    
    v_m_hat = (p_m + Tu_m_hat)/(Z1_s) + 1.0_wp/(alpha + eta_s)*phi_m          !
    v_l_hat = (p_l + Tu_l_hat)/(Z1_s) + 1.0_wp/(alpha + eta_s)*phi_l

    u_m_hat = (q_m - Tv_m_hat)/(Z2_s) - 1.0_wp/(alpha + eta_s)*phi_m
    u_l_hat = (q_l - Tv_l_hat)/(Z2_s) - 1.0_wp/(alpha + eta_s)*phi_l

    ! Save all hat-variables in the work-array  u_out, v_out

    u_out(1) = u_n_hat                                                      ! continuity of normal velocity

    ! constrain tangential velocities by the continuity of the left going shear characteristics
    u_out(3) = u_l_hat
    u_out(4) = Tu_n_hat                                                   ! continuity of normal stress

    ! shear stress given by friction law
    u_out(5) = Tu_m_hat
    u_out(6) = Tu_l_hat

    v_out(1) = v_n_hat                                                    ! continuity of normal velocity

    v_out(2) = v_m_hat
    v_out(3) = v_l_hat

    v_out(4) = Tv_n_hat                                                   ! continuity of normal stress

    ! shear stress given by friction law
    v_out(5) = Tv_m_hat
    v_out(6) = Tv_l_hat

  end subroutine Linear_Friction


  subroutine  Slip_Weakening_Friction(problem,u_in, v_in, u_out, v_out, Vel, Z1_s, Z2_s, Z1_p, Z2_p, &
       slip, x_ij, t, l, m, n, coupling, Mat)
    
    ! This subroutine solves for the hat variables governed by a slip weakening friction law
    ! All variables are in the locally rotated orthogonal coordinates
    
    use datatypes, only : block_material
    
    implicit none
    
    character(256), intent(in) :: problem
    real(kind = wp), dimension(1:6), intent(in) :: u_in, v_in                                             ! Hat-variables
    real(kind = wp), dimension(1:6), intent(out) :: u_out, v_out                                          ! Hat-variables
    real(kind = wp), intent(out) :: Vel                                                 ! slip-rate
    real(kind = wp), intent(in) :: Z1_s, Z2_s, Z1_p, Z2_p                                               ! Local impedance
    real(kind = wp), intent(in) :: slip                                                                    ! Slip
    real(kind = wp), intent(in) :: t                                                                    ! Time
    real(kind = wp), intent(in) :: Mat(3)                                                                    ! Shear modulus
    real(kind = wp), dimension(1:3), intent(in) :: x_ij
    real(kind = wp), dimension(1:3), intent(in) :: n, l, m
    character(*), intent(in) :: coupling
    
    
    real(kind = wp) :: T0_n, T0_m, T0_l                                                     ! Background stress
    
    real(kind = wp) :: u_n, u_m, u_l, Tu_n, Tu_m, Tu_l  ! Local velocities and stresses for the left block
    real(kind = wp) :: v_n, v_m, v_l, Tv_n, Tv_m, Tv_l  ! Local velocities and stresses for the right block
    
    real(kind = wp) :: u_n_hat, u_m_hat, u_l_hat, &
         Tu_n_hat, Tu_m_hat, Tu_l_hat   ! Local velocities_hat and stresses_hat for the left block
    real(kind = wp) :: v_n_hat, v_m_hat, v_l_hat, &
            Tv_n_hat, Tv_m_hat, Tv_l_hat   ! Local velocities_hat and stresses_hat for the right block
    real(kind = wp) :: p_n, p_m, p_l, q_n, q_m, q_l                                                     ! Local characteristic variables
    real(kind = wp) :: phi_l, phi_m, phi_n                                                              ! Stress tranfer functions
    real(kind = wp) :: eta_s, eta_p                                                                     ! Radiation damping coefficients
    real(kind = wp) :: vv_m, vv_l                                                                       ! Slip velocities
    real(kind = wp) :: tau_str, tau_m, tau_l, sigma_n                                    ! Friction parameters
    integer :: No_it = 200                                                        ! Maximum number of Newton iterations

    real(kind = wp) :: tau_lock
    
    if((Z1_s .le. 1.0e-9_wp) .or. (Z2_s .le. 1.0e-9_wp)) STOP 'shear impedance must be nonzero'
    
    ! velocities in the left elastic block
    u_n = u_in(1)
    u_m = u_in(2)
    u_l = u_in(3)
    
    ! velocities in the right elastic block
    v_n = v_in(1)
    v_m = v_in(2)
    v_l = v_in(3)

    ! stress fields in the left elastic block
    Tu_n = u_in(4)
    Tu_m = u_in(5)
    Tu_l = u_in(6)
    
    ! stress fields in the right elastic block
    Tv_n = v_in(4)
    Tv_m = v_in(5)
    Tv_l = v_in(6)
    
    ! compute characteristics
    p_n = Z1_p*u_n - Tu_n
    p_m = Z1_s*u_m - Tu_m
    p_l = Z1_s*u_l - Tu_l
    
    q_n = Z2_p*v_n + Tv_n
    q_m = Z2_s*v_m + Tv_m
    q_l = Z2_s*v_l + Tv_l

    eta_s = 1.0_wp/(1.0_wp/Z1_s + 1.0_wp/Z2_s)                                    ! half of the harmonic mean of Z1_s, Z2_s
    eta_p = 1.0_wp/(1.0_wp/Z1_p + 1.0_wp/Z2_p)                                    ! half of the harmonic mean of Z1_p, Z2_p
    
    ! get prestress (where normal traction is effective normal traction)

    call prestress(T0_n, T0_m, T0_l, x_ij, t, problem, l, m, n, Mat)

    ! stress transfer function, without prestress
    
    phi_n = eta_p*(q_n/Z2_p - p_n/Z1_p)
    phi_m = eta_s*(q_m/Z2_s - p_m/Z1_s)
    phi_l = eta_s*(q_l/Z2_s - p_l/Z1_s)
    
    ! Compute Hat-Variables for normal fields
    u_n_hat = (q_n - phi_n)/Z2_p                                                                   ! continuity of normal velocity
    Tu_n_hat = phi_n                                                                              ! continuity of normal stress
    v_n_hat = (p_n + phi_n)/Z1_p                                                                   ! continuity of normal velocity
    Tv_n_hat = phi_n                                                                              ! continuity of normal stress
    
    ! Compute Hat-Variables for tangential fields

    ! **note that tau below includes prestress**
    
    ! first calculate magnitude of shear stress, if fault is constrained against further slip
    ! (including prestress)
    
    
    tau_lock = sqrt((T0_l + phi_l)**2 + (T0_m + phi_m)**2)
    
    ! effective normal stress and shear strength
    sigma_n = max(0.0_wp, -(T0_n + phi_n)) ! including prestress
    
    call Tau_strength(problem, tau_str, sigma_n, slip, x_ij, t)
    
    ! now compare shear stress if fault is locked to frictional strength
    if (tau_lock > tau_str) then ! fault slipping

       
       ! solve for slip velocity and shear stress
       call solve(problem,vv_m,vv_l, Vel, tau_m,tau_l,T0_m + phi_m,T0_l + phi_l,eta_s,&
            tau_str,No_it,coupling,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp)
       
       ! set all hat-variables (subtracting prestress)
       Tu_m_hat = tau_m - T0_m
       Tu_l_hat = tau_l - T0_l

       Tv_m_hat = tau_m - T0_m
       Tv_l_hat = tau_l - T0_l
       
    else ! fault locked
       
       Tu_m_hat = phi_m
       Tu_l_hat = phi_l
       
       Tv_m_hat = phi_m
       Tv_l_hat = phi_l
       
       vv_m = 0.0_wp
       vv_l = 0.0_wp

       Vel = 0.0_wp
    end if
    
    v_m_hat = (Tu_m_hat + p_m)/(Z1_s) + vv_m
    v_l_hat = (Tu_l_hat + p_l)/(Z1_s) + vv_l
    
    u_m_hat = (q_m - Tv_m_hat)/(Z2_s) - vv_m
    u_l_hat = (q_l - Tv_l_hat)/(Z2_s) - vv_l
    
    !print *, v_m_hat, v_l_hat, u_m_hat, u_l_hat
    !stop

    ! Save all hat-variables in the work-array  u_out, v_out
    
    ! continuity of normal velocity
    
    u_out(1) = u_n_hat
    v_out(1) = v_n_hat
    
    ! constrain tangential velocities by the continuity of the left- and right-going shear characteristics
    
    u_out(2) = u_m_hat
    u_out(3) = u_l_hat
    
    v_out(2) = v_m_hat
    v_out(3) = v_l_hat

    ! continuity of normal stress
    
    u_out(4) = Tu_n_hat
    v_out(4) = Tv_n_hat

    ! shear stress given by friction law
    u_out(5) = Tu_m_hat
    u_out(6) = Tu_l_hat
    
    v_out(5) = Tv_m_hat
    v_out(6) = Tv_l_hat
    
  end subroutine Slip_Weakening_Friction
  
  subroutine  rate_state_friction(problem,u_in, v_in, u_out, v_out, Vel, Z1_s, Z2_s, Z1_p, Z2_p, &
       psi,a,V0, x_ij, t, l, m, n,coupling)
    
    ! This subroutine solves for the hat variables governed by a rate-and-state friction law
    ! All variables are in the locally rotated orthogonal coordinates
    
    implicit none
    
    character(256), intent(in) :: problem
    real(kind = wp), dimension(:), intent(in) :: u_in, v_in           ! local grid-functions to be used to compute hat-functions
    real(kind = wp), dimension(:), intent(out) :: u_out, v_out        ! array used to store hat-variables for output
    real(kind = wp), intent(out) :: Vel                               ! slip-rate
    real(kind = wp), intent(in) :: Z1_s, Z2_s, Z1_p, Z2_p             ! local impedance
    
    real(kind = wp), intent(in) :: psi                ! state variable
    real(kind = wp), intent(in) :: a                  ! direct effect parameter
    real(kind = wp), intent(in) :: V0                 ! reference velocity

    real(kind = wp), intent(in) :: t                                  ! time
    real(kind = wp), dimension(:), intent(in) :: x_ij                 ! spatial coordinates
    real(kind = wp), dimension(:), intent(in) :: n, l, m              ! unit vectors
    character(*), intent(in) :: coupling                   ! specifies the type of interface condition used
    
    real(kind = wp) :: T0_n, T0_m, T0_l                                                    ! background tractions
    real(kind = wp) :: u_n, u_m, u_l, Tu_n, Tu_m, Tu_l                                     ! local velocities and stresses for the left block
    real(kind = wp) :: v_n, v_m, v_l, Tv_n, Tv_m, Tv_l                                     ! local velocities and stresses for the right block
    real(kind = wp) :: u_n_hat, u_m_hat, u_l_hat, &
    Tu_n_hat, Tu_m_hat, Tu_l_hat             ! local velocities_hat and stresses_hat for the left block
    real(kind = wp) :: v_n_hat, v_m_hat, v_l_hat, &
         Tv_n_hat, Tv_m_hat, Tv_l_hat             ! local velocities_hat and stresses_hat for the right block

    real(kind = wp) :: sigma_n                                                             ! compressive normal stress
    
    
    real(kind = wp) :: p_n, p_m, p_l, q_n, q_m, q_l                                                     ! local characteristic variables
    real(kind = wp) :: phi_l, phi_m, phi_n                                                              ! stress tranfer functions
    real(kind = wp) :: eta_s, eta_p                                                                     ! radiation damping coefficients
    real(kind = wp) :: vv_m, vv_l                                                       ! slip velocities (hat-variables)
    real(kind = wp) :: tau_m, tau_l                                                     ! shear tractions (hat-variables)
    integer :: No_it = 200                                                   ! maximum number of Newton iterations
    real(kind = wp) :: tau_str = 1.0_wp
    real(kind = wp) :: r
    
    if((Z1_s .le. 1.0e-9_wp) .or. (Z2_s .le. 1.0e-9_wp)) STOP 'shear impedance MUST be NONZERO'
    
    ! velocities in the left elastic block
    u_n = u_in(1)
    u_m = u_in(2)
    u_l = u_in(3)
    
    ! velocities in the right elastic block
    v_n = v_in(1)
    v_m = v_in(2)
    v_l = v_in(3)

    ! stress fields in the left elastic block
    Tu_n = u_in(4)
    Tu_m = u_in(5)
    Tu_l = u_in(6)
    
    Tv_n = v_in(4)
    Tv_m = v_in(5)
    Tv_l = v_in(6)
    
    ! compute characteristics

    ! positive going characteristics
    p_n = Z1_p*u_n - Tu_n
    p_m = Z1_s*u_m - Tu_m
    p_l = Z1_s*u_l - Tu_l
    
    ! negative going characteristics
    q_n = Z2_p*v_n + Tv_n
    q_m = Z2_s*v_m + Tv_m
    q_l = Z2_s*v_l + Tv_l
    
    eta_s = 1.0_wp/(1.0_wp/Z1_s + 1.0_wp/Z2_s)                                    ! half of the harmonic mean of Z1_s, Z2_s
    eta_p = 1.0_wp/(1.0_wp/Z1_p + 1.0_wp/Z2_p)                                    ! half of the harmonic mean of Z1_p, Z2_p
    
    ! get prestress (where normal traction is effective normal traction)
    call prestress(T0_n, T0_m, T0_l, x_ij, t, problem, l, m, n)
    
    ! stress transfer function, without prestress
    phi_n = eta_p*(q_n/Z2_p - p_n/Z1_p)
    phi_m = eta_s*(q_m/Z2_s - p_m/Z1_s)
    phi_l = eta_s*(q_l/Z2_s - p_l/Z1_s)
    
    ! Compute Hat-Variables for normal fields
    u_n_hat = (q_n - phi_n)/Z2_p                                                                  ! continuity of normal velocity
    Tu_n_hat = phi_n                                                                              ! continuity of normal stress
    
    v_n_hat = (p_n + phi_n)/Z1_p                                                                  ! continuity of normal velocity
    Tv_n_hat = phi_n                                                                              ! continuity of normal stress
    
    ! Compute Hat-Variables for tangential fields
    ! **note that tau below includes prestress**
    ! first calculate magnitude of shear stress, if fault is constrained against further slip
    
    ! effective normal stress and shear strength
    sigma_n = max(0.0_wp, -(T0_n + phi_n)) ! including prestress
    r = sqrt((x_ij(2)-7.5_wp)**2 + x_ij(3)**2)
    

    tau_m = 0.0_wp
    tau_l = 0.0_wp

    vv_m = v_in(2)-u_in(2)
    vv_l = v_in(3)-u_in(3)
    
    ! solve for slip velocity and shear stress
    
    call solve(problem, vv_m,vv_l, Vel, tau_m,tau_l, T0_m+phi_m,T0_l+phi_l, eta_s,&
         tau_str,No_it,coupling,sigma_n,psi,V0,a,r)
    
    ! set all hat-variables (subtracting prestress)
    Tu_m_hat = tau_m - T0_m
    Tu_l_hat = tau_l - T0_l
    
    Tv_m_hat = tau_m - T0_m
    Tv_l_hat = tau_l - T0_l
    
    v_m_hat = (Tu_m_hat + p_m)/(Z1_s) + vv_m
    v_l_hat = (Tu_l_hat + p_l)/(Z1_s) + vv_l
    
    u_m_hat = (q_m - Tv_m_hat)/(Z2_s) - vv_m
    u_l_hat = (q_l - Tv_l_hat)/(Z2_s) - vv_l
        
    ! continuity of normal velocity
    u_out(1) = u_n_hat
    v_out(1) = v_n_hat
    
    ! constrain tangential velocities by the continuity of the left- and right-going shear characteristics
    u_out(2) = u_m_hat
    u_out(3) = u_l_hat
    
    v_out(2) = v_m_hat
    v_out(3) = v_l_hat
    
    ! continuity of normal stress
    u_out(4) = Tu_n_hat
    v_out(4) = Tv_n_hat

    ! shear stress given by friction law
    u_out(5) = Tu_m_hat
    u_out(6) = Tu_l_hat
    
    v_out(5) = Tv_m_hat
    v_out(6) = Tv_l_hat
    
  end subroutine rate_state_friction
  
  subroutine test_nonlinear_solver()
    
    implicit none
    
    real(kind = wp) :: v1,v2,tau1,tau2,phi1,phi2
    real(kind = wp),parameter :: eta=1.0_wp, tau_str=0.0_wp
    integer,parameter :: N=50
    character(*),parameter ::coupling = 'slip-weakening_friction'
    character(256),parameter ::problem = 'TPV5'
    real(kind = wp) :: Vel
    phi1 = 1.0_wp
    phi2 = 0.0_wp
    
    v1 = 1.0_wp
    v2 = 1.0_wp
    tau1 = 1.0_wp
    tau2 = 1.0_wp
    
    call solve(problem,v1,v2, Vel, tau1,tau2, phi1,phi2,eta,tau_str,N,coupling,1.0_wp,1.0_wp,1.0_wp,1.0_wp,1.0_wp)
    
    print *, v1, v2, tau1, tau2, sqrt(tau1**2+tau2**2)
    stop
    
  end subroutine test_nonlinear_solver

  
  subroutine solve(problem,v1, v2, Vel, tau1, tau2, phi_1, phi_2, eta, tau_str, N,coupling,sigma_n,psi,V0,a,r)
    
    implicit none
    
    character(256),intent(in)::problem
    real(kind = wp), intent(in) :: phi_1, phi_2, eta, tau_str
    integer, intent(in) :: N
    real(kind = wp), intent(out) :: Vel
    real(kind = wp), intent(inout) :: v1, v2, tau1, tau2
    character(*),intent(in)::coupling
    real(kind = wp), intent(in) :: sigma_n            ! compressive normal stress
    real(kind = wp), intent(in) :: psi                ! state variable
    real(kind = wp), intent(in) :: a                  ! direct effect parameter
    real(kind = wp), intent(in) :: V0                 ! reference velocity
    real(kind = wp), intent(in) :: r
    
    real(kind = wp) ::  Jf1, Jf2, err
    integer :: k
    real(kind = wp),parameter :: tol=1.0e-9_wp
    real(kind = wp) ::  fv, V, Phi, fpv 
    err = 1.0_wp
    k = 1
    
    select case(coupling)
    
    case('rate-and-state_friction')
       ! solve a nonlinear problem for slip rate
       Phi = sqrt(phi_1**2 + phi_2**2)                 ! stress-transfer functional
       
       ! initialize V
       V = sqrt(v1**2 + v2**2)
       if (V > Phi) V = 0.5_wp*Phi/eta  
       
       ! solve a nonlinear problem for slip-rate: V
       call Regula_Falsi(problem,V,Phi,eta,sigma_n,psi,V0,a)
       
       if (V < 0.0_wp) then 
          print*, 'negative slip-rate.'
          stop
       end if
       
       Vel = V
       
       ! compute slip velocities
       fv = sigma_n*a*asinh(0.5_wp*V/V0*exp(psi/a))
       v1 = phi_1/(eta+fv/V)
       v2 = phi_2/(eta+fv/V)
       
       ! compute shear stress on the fault
       tau1 = phi_1 - eta*v1 
       tau2 = phi_2 - eta*v2 
       
    case('slip-weakening_friction')
       !solve a linear problem for slip-rate: V  
       Phi = sqrt(phi_1**2 + phi_2**2)        ! stress-transfer functional
       V = (Phi - tau_str)/eta                ! slip-rate

       if (V < 0.0_wp) then 
          print*, 'negative slip-rate.'
          stop
       end if
       
       Vel = V
       !  compute slip velocities
       v1 = phi_1/(eta+tau_str/V)
       v2 = phi_2/(eta+tau_str/V)
       
       !  compute shear stress on the fault
       tau1 = phi_1 - eta*v1
       tau2 = phi_2 - eta*v2
       
    end select
    
  end subroutine solve
  
  subroutine Regula_Falsi(problem,V,Phi,eta,sigma_n,psi,V0,a)
    ! a nonlinear solver for slip-rate using Regula-Falsi
    ! solve: V + sigma_n/eta*f(V,psi) - Phi/eta = 0                                               
    
    implicit none
    character(256),intent(in)::problem
    real(kind = wp), intent(inout) :: V          ! magnitude of slip-velocity (slip-rate)                                   
    real(kind = wp), intent(in) :: Phi           ! stress transfer functions                                                    
    real(kind = wp), intent(in) :: eta           ! half of harmonic average of shear impedance                                  
    real(kind = wp), intent(in) :: sigma_n       ! compressive normal stress                                                    
    real(kind = wp), intent(in) :: psi           ! state variable                                                               
    real(kind = wp), intent(in) :: a             ! direct effect parameter                                                      
    real(kind = wp), intent(in) :: V0            ! reference velocity                                                           
    real(kind = wp) :: fv, fl, fr, Vl, Vr, V1    ! function values
    real(kind = wp) :: err, tol 
    integer ::  k, maxit
    
    err = 1.0_wp
    tol = 1.0e-12_wp
    
    Vl = 0.0_wp             ! lower bound
    Vr = Phi/eta            ! upper bound
    
    k = 1
    maxit = 5000
    
    call f(problem,fv,V,Phi,eta,sigma_n,psi,V0,a)
    call f(problem,fl,Vl,Phi,eta,sigma_n,psi,V0,a)
    call f(problem,fr,Vr,Phi,eta,sigma_n,psi,V0,a)
    
    if (abs(fv) .le. tol) then ! check if the guess is a solution 
       V = V 
       return
    end if
    
    if (abs(fl) .le. tol) then ! check if the lower bound is a solution 
       V = Vl 
       return
    end if
    
    if (abs(fr) .le. tol) then  ! check if the upper bound is a solution 
       V = Vr 
       return
    end if
    
    ! else solve for V iteratively using regula-falsi
    do while((k .le. maxit) .and. (err .ge. tol))
       
       V1 = V                                                             
       
       if (fv*fl > tol) then
          Vl = V
          V = Vr - (Vr-Vl)*fr/(fr-fl)
          fl = fv
       elseif (fv*fr > tol) then 
          Vr = V
          V = Vr -(Vr-Vl)*fr/(fr-fl)
          fr = fv
       end if
       
       call f(problem,fv,V,Phi,eta,sigma_n,psi,V0,a)
       
       err = abs(V-V1)
       k = k + 1
       
    end do
    
    if (k .ge. maxit) print*, err, k, 'maximum iteration reached. solution didnt convergence!'
    
  end subroutine Regula_Falsi
  
  subroutine f(problem,fv,V,Phi,eta,sigma_n,psi,V0,a)
    !                                                                   
    ! This subroutine evaluates the function:                                                                  
    ! f:= V + sigma_n/eta*a*f(V,psi) - Phi/eta = 0
    ! where V = sqrt(v1**2 + v2**2)
    
    implicit none
    character(256),intent(in)::problem
    real(kind = wp), intent(in) :: V             ! magnitude of slip-velocity (slip-rate) 
    real(kind = wp), intent(in) :: Phi           ! stress transfer functions                                                
    real(kind = wp), intent(in) :: eta           ! half of harmonic average of shear impedance                              
    real(kind = wp), intent(in) :: sigma_n       ! compressive normal stress                                                
    real(kind = wp), intent(in) :: psi           ! state variable                                                           
    real(kind = wp), intent(in) :: a             ! direct effect parameter                                                  
    real(kind = wp), intent(in) :: V0            ! reference velocity                                                       
    real(kind = wp), intent(inout) ::fv          ! function values                                              
    
    !select case(problem)
       
    !case('rate-weakening')
       
       fv = V + a*sigma_n/eta*asinh(0.5_wp*V/V0*exp(psi/a)) - Phi/eta
       
    !end select
  end subroutine f
  
  
  subroutine Tau_strength(problem, tau_str, sigma_n, D, x, t)
    
    implicit none
    
    character(256), intent(in) :: problem
    real(kind = wp), intent(in) :: D,x(:),t,sigma_n
    real(kind = wp), intent(out) :: tau_str
    
    real(kind = wp) :: fs, fd, Dc, f, f1, f2, t0, r
    real(kind = wp) :: time_of_forced_rup, C0, Vs
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp) :: ang
    
    select case(problem)
       
    case('TPV5')

       ! interior region

       fs = 0.677_wp
       fd = 0.525_wp
       Dc = 0.4_wp

       ! add high-strength barrier outside fault
       ! Compute static frictional coefficient
       fs = fs + 1.0e10_wp*(1.0_wp-boxcar(x(2),15.0_wp)*boxcar(x(3),15.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = sigma_n*f
       
    case('MG01', 'MG2')

       ! interior region                        
       fs = 0.6778_wp
       fd = 0.525_wp
       Dc = 0.4_wp

       ! add high-strength barrier outside fault
       ! Compute static frictional coefficient 

       fs = fs + 1.0e10_wp*(1.0_wp-boxcar(x(2),15.0_wp)*boxcar(x(3),15.0_wp))

       ! Compute frictional coefficient                                      
       f = fs-(fs-fd)*min(D,Dc)/Dc                            ! linear slip-weakening

       ! Compute tau_strength                                                                                    
       tau_str = sigma_n*f

    case('mg_a1', 'mg_b1a')

       ! interior region                                                                                                                             
       fs = 0.677_wp
       fd = 0.373_wp
       Dc = 0.4_wp
       
       ! add high-strength barrier outside fault                                                                                                     
       ! Compute static frictional coefficient                                                                                                       

       fs = fs + (1.0e9_wp-fs)*(1.0_wp-boxcar(x(2),15.4_wp)*boxcar(x(3),15.4_wp))
       fd = fd + (1.0e9_wp-fd)*(1.0_wp-boxcar(x(2),15.4_wp)*boxcar(x(3),15.4_wp))
       Dc = Dc + (1.0e3_wp-Dc)*(1.0_wp-boxcar(x(2),15.4_wp)*boxcar(x(3),15.4_wp))

        ! Compute frictional cohesion
       call Frictional_Cohesion(C0, x(2), problem) 

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc                  ! linear slip-weakening
       
       ! Compute tau_strength
       tau_str = C0 + sigma_n*f  
       ! Compute frictional coefficient                                                                                                              


    case('TPV5_2D')

       fs = 0.677_wp
       fd = 0.525_wp
       Dc = 0.4_wp

       ! Compute static frictional coefficient
       fs = fs + 1.0e10_wp*(1.0_wp-boxcar(x(3),15.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening
       ! Compute tau_strength
       tau_str = sigma_n*f

    case('TPV5_fractal')

       ! interior region
       fs = 0.677_wp
       fd = 0.373_wp
       Dc = 0.4_wp

       ! add high-strength barrier outside fault
       ! Compute static frictional coefficient
       fs = fs + 1.0e10_wp*(1.0_wp-boxcar(x(2),15.0_wp)*boxcar(x(3),18.0_wp))

       ! Compute frictional coefficient

       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = sigma_n*f

    case('TPV26')

       fs = 0.18_wp
       fd = 0.12_wp
       Dc = 0.3_wp
       t0 = 0.5_wp
       Vs = 3.464_wp
       !f2 = 0.0_wp

       r = sqrt((x(2)-10.0_wp)**2 + (x(3)+5.0_wp)**2)

       ! Compute frictional cohesion
       call Frictional_Cohesion(C0, x(2), 'TPV26')

       ! Compute time of forced rupture
       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, 'TPV26')

       f1 = min(D/Dc, 1.0_wp)

       if (t < time_of_forced_rup) then

          f2 = 0.0_wp

       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then

          f2 = (t-time_of_forced_rup)/t0

       else if(t .ge. time_of_forced_rup+t0) then

          f2 = 1.0_wp

       end if

       fs = fs + 1.0e10_wp*(1.0_wp-boxcar(x(2),20.0_wp)*boxcar(x(3),20.0_wp))

       ! Compute frictional coefficient

       f = fs-(fs-fd)*max(f1,f2)

       ! Compute tau_strength
       tau_str = C0 + f*sigma_n

    case('TPV27')
       
       fs = 0.18_wp
       fd = 0.12_wp
       Dc = 0.3_wp
       t0 = 0.5_wp
       Vs = 3.464_wp
       !f2 = 0.0_wp                                  
       
       r = sqrt((x(2)-10.0_wp)**2 + (x(3)+5.0_wp)**2)
       ! Compute frictional cohesion                                                                                                                          
       call Frictional_Cohesion(C0, x(2), 'TPV27')
       
       ! Compute time of forced rupture                                                                                                                       
       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, 'TPV27')

       f1 = min(D/Dc, 1.0_wp)

       if (t < time_of_forced_rup) then

          f2 = 0.0_wp

       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then

          f2 = (t-time_of_forced_rup)/t0

       else if(t .ge. time_of_forced_rup+t0) then

          f2 = 1.0_wp

       end if

       fs = fs + 1.0e10*(1.0_wp-boxcar(x(2),20.0_wp)*boxcar(x(3),20.0_wp))
       
       ! Compute frictional coefficient                                                                                                                       

       f = fs-(fs-fd)*max(f1,f2)

       ! Compute tau_strength                                                                                                                                 
       tau_str = C0 + f*sigma_n

    case('TPV28')

       ! interior region
       fs = 0.677_wp
       fd = 0.373_wp
       Dc = 0.4_wp

       ! add high-strength barrier outside fault
       ! Compute static frictional coefficient
       fs = fs + 1.0e10_wp*(1.0_wp-boxcar(x(2),15.0_wp)*boxcar(x(3),18.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = sigma_n*f

      
    case('TPV29','TPV30')
       
       fs = 0.18_wp
       fd = 0.12_wp
       Dc = 0.3_wp
       t0 = 0.5_wp
       Vs = 3.464_wp

       ! x(1) = depth, x(2) = along strike
       
       r = sqrt((x(2)-10.0_wp)**2 + (x(3)+5.0_wp)**2)
       
       ! Compute frictional cohesion                                                                        
       call Frictional_Cohesion(C0, x(2), problem)

       ! lock fault by increasing cohesion (and hence strength) outside of slipping region
       C0 = C0 + 1.0e10*(1.0_wp-boxcar(x(2),20.0_wp)*boxcar(x(3),20.0_wp))
       
       ! Compute time of forced rupture                                                                     
       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, problem)
       
       f1 = min(D/Dc, 1.0_wp)

       if (t < time_of_forced_rup) then
          
          f2 = 0.0_wp
       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then
          
          f2 = (t-time_of_forced_rup)/t0
          
       else if(t .ge. time_of_forced_rup+t0) then
          
          f2 = 1.0_wp
          
       end if
       
       ! Compute frictional coefficient                                                                     
       f = fs-(fs-fd)*max(f1,f2)

       ! Compute tau_strength                                                                               
       tau_str = C0 + f*sigma_n


    case('TPV31','TPV32')

       fs = 0.58_wp
       fd = 0.45_wp
       Dc = 0.18_wp

       ! Compute frictional cohesion
       call Frictional_Cohesion(C0, x(2), 'TPV31') 

       ! Compute static frictional coefficient
       fs = fs + 1d10*(1_wp-boxcar(x(2),15.0_wp)*boxcar(x(3),15.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = C0 + max(sigma_n,0.0_wp)*f 

   case('TPV33')

       fs = 0.55_wp
       fd = 0.45_wp
       Dc = 0.18_wp

       C0 = 0.0_wp

       ! Compute static frictional coefficient
       fs = fs + 1d10*(1_wp-boxcar(x(1),10.0_wp)*boxcar(x(2)+4.0_wp,8.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = C0 + max(sigma_n,0.0_wp)*f 

   case('TPV34')

       fs = 0.58_wp
       fd = 0.45_wp
       Dc = 0.18_wp

       ! Compute frictional cohesion
       call Frictional_Cohesion(C0, x(1), 'TPV34') 

       ! Compute static frictional coefficient
       fs = fs + 1d10*(1.0_wp-boxcar(x(1),15.0_wp)*boxcar(x(2),15.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = C0 + max(sigma_n,0.0_wp)*f 

   case('SCITS2016')

       fs = 0.58_wp
       fd = 0.35_wp
       Dc = 0.09_wp

       ! Compute frictional cohesion
       call Frictional_Cohesion(C0, x(2), 'SCITS2016') 

       ! Compute static frictional coefficient
       fs = fs + 1d10*(1_wp-boxcar(x(2)-5.0_wp,2.0_wp)*boxcar(x(3)-20.0_wp,3.0_wp))

       ! Compute frictional coefficient
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening

       ! Compute tau_strength
       tau_str = C0 + max(sigma_n,0.0_wp)*f   

    case('TPV10')

       fs = 0.760_wp
       fd = 0.448_wp
       Dc = 0.50_wp

       C0 = 0.2_wp

       ! Compute static frictional coefficient                                                                                                                                                   
       fs = fs + (1.0e4_wp-0.760_wp)*(1_wp-boxcar(x(2)*(2.0_wp/(sqrt(3.0_wp))),15.0_wp)*boxcar(x(3)-0.0_wp,15.0_wp))

       ! Compute frictional coefficient                                                                                                                                                          
       f = fs-(fs-fd)*min(D,Dc)/Dc ! linear slip-weakening                                                                                                                                       
       C0 = C0 + 1.0e3_wp*(1.0_wp-boxcar(x(2)*(2.0_wp/(sqrt(3.0_wp))),15.0_wp)*boxcar(x(3)-0.0_wp,15.0_wp))
       ! Compute tau_strength                                                                                                                                                                    
       tau_str = C0 + max(sigma_n,0.0_wp)*f


    case('TPV36')

       fs = 0.575_wp
       fd = 0.450_wp
       Dc = 0.180_wp
       t0 = 0.5_wp
       Vs = 3.464_wp
       
       ! x(2) = depth, x(3) = along strike
       r = sqrt((x(2)/sin(15.0_wp/180.0_wp*pi)-18.0_wp)**2 + x(3)**2)

       ! Compute frictional cohesion
       call Frictional_Cohesion(C0, x(2)/sin(15.0_wp/180.0_wp*pi), problem)

       ! lock fault by increasing cohesion (and hence strength) outside of slipping region
       C0 = C0 + 1.0e10*(1.0_wp-boxcar(x(2)/sin(15.0_wp/180.0_wp*pi),28.0_wp)*boxcar(x(3),15.0_wp))

       ! Compute time of forced rupture
       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, problem)

       f1 = min(D/Dc, 1.0_wp)
       if (t < time_of_forced_rup) then

          f2 = 0.0_wp
       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then

          f2 = (t-time_of_forced_rup)/t0

       else if(t .ge. time_of_forced_rup+t0) then

          f2 = 1.0_wp

       end if
       
       f = fs-(fs-fd)*max(f1,f2)

       ! Compute tau_strength
       tau_str = C0 + f*max(sigma_n,0.0_wp)


       case('TPV37')

       fs = 0.575_wp
       fd = 0.450_wp
       Dc = 0.180_wp
       t0 = 0.5_wp
       Vs = 3.464_wp

       ! x(2) = depth, x(3) = along strike                                                                                                                                                                 
       r = sqrt((x(2)/sin(15.0_wp/180.0_wp*pi)-18.0_wp)**2 + x(3)**2)

       ! Compute frictional cohesion                                                                                                                                                                       
       call Frictional_Cohesion(C0, x(2)/sin(15.0_wp/180.0_wp*pi), problem)

       ! lock fault by increasing cohesion (and hence strength) outside of slipping region                                                                                                                 
       C0 = C0 + 1.0e10*(1.0_wp-boxcar(x(2)/sin(15.0_wp/180.0_wp*pi),28.0_wp)*boxcar(x(3),15.0_wp))

       ! Compute time of forced rupture                                                                                                                                                                    
       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, problem)

       f1 = min(D/Dc, 1.0_wp)
       if (t < time_of_forced_rup) then

          f2 = 0.0_wp
       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then

          f2 = (t-time_of_forced_rup)/t0

       else if(t .ge. time_of_forced_rup+t0) then

          f2 = 1.0_wp

       end if

       f = fs-(fs-fd)*max(f1,f2)

       ! Compute tau_strength                                                                                                                                                                              
       tau_str = C0 + f*max(sigma_n,0.0_wp)
       
    case('wasatch100')

       ! interior region                                                                                                                                                                                
       fs = 0.6_wp
!       fd = 0.45_wp                                                                                                                                                                                    

       if ((x(2)+2.149d0) .le. (4.0_wp)) then
          fd = 0.44_wp + (0.6_wp-0.44_wp)*0.25_wp*(4.0_wp - (x(2)+2.149d0))
       elseif (((x(2)+2.149d0) .gt. 4.0_wp) .and. ((x(2)+2.149d0) .le. 17.0_wp)) then
          fd = 0.44_wp
       elseif (((x(2)+2.149d0) .gt. 17.0_wp) .and. (x(2)+2.149d0) .le. 18.0_wp) then
          fd = 0.44_wp+(0.6_wp-0.44_wp)*((x(2)+2.149d0)-17.0_wp)
       elseif (((x(2)+2.149d0) .gt. 18.0_wp)) then
             fd = 0.6_wp
       end if


       if ((x(2)+2.149d0) .le. 4.0_wp) then
          Dc = 0.3_wp + 0.2_wp*(4.0_wp - (x(2)+2.149d0))
       elseif ((x(2)+2.149d0) .ge. 4.0_wp) then
          Dc = 0.3_wp
       end if

!       Dc = 0.3_wp                                                                                                                                                                                     

       t0 = 0.5_wp !added in for tapered circle                                                                                                                                                         
       Vs = 3.464_wp !added in for tapered circle                                                                                                                                                                                                                                                                                           
       r = sqrt((1.3054_wp*((x(2)+2.149d0)-14.0_wp))**2 + (x(3)-38.0_wp)**2) !for ellipse, i.e. equiavelent circle
       ! add high-strength barrier outside fault                                                                                                                                                        
       ! Compute static frictional coefficient                                                                                                                                                          
       fs = fs + 1.0e10*(1.0_wp-boxcar((x(2)+2.149d0),18.0_wp)*boxcar(x(3)-38.0_wp,18.6_wp))
       ! Compute frictional cohesion                                                                                                                                                                    
       call Frictional_Cohesion(C0, (x(2)+2.149d0), problem)

       ! lock fault by increasing cohesion (and hence strength) outside of slipping region                                                                                                              
       C0 = C0 + 1.0e10*(1.0_wp-boxcar((x(2)+2.149d0),18.0_wp)*boxcar(x(3)-38.0_wp,18.6_wp))
       ! Compute frictional coefficient                                                                                                                                                                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                                                   
       ! Compute time of forced rupture                                                                 \                                                                                               

       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, problem)
       f1 = min(D/Dc, 1.0_wp)
       if (t < time_of_forced_rup) then
          f2 = 0.0_wp
       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then
          f2 = (t-time_of_forced_rup)/t0
       else if(t .ge. time_of_forced_rup+t0) then
          f2 = 1.0_wp
       end if
       f = fs-(fs-fd)*max(f1,f2) !for time of forced rupture weakening                                     
       tau_str = C0 + max(sigma_n,0.0_wp)*f

!!$       case('wasatch100')
!!$
!!$       ! interior region                                                                                                                                                                                
!!$       fs = 0.6_wp
!!$!       fd = 0.45_wp                                                                                                                                                                                    
!!$
!!$       if ((x(2)+0d0*2.781d0) .le. (4.0_wp)) then
!!$          fd = 0.44_wp + (0.6_wp-0.44_wp)*0.25_wp*(4.0_wp - (x(2)+0d0*2.781d0))
!!$       elseif (((x(2)+0d0*2.781d0) .gt. 4.0_wp) .and. ((x(2)+0d0*2.781d0) .le. 17.0_wp)) then
!!$          fd = 0.44_wp
!!$       elseif (((x(2)+0d0*2.781d0) .gt. 17.0_wp) .and. (x(2)+0d0*2.781d0) .le. 18.0_wp) then
!!$          fd = 0.44_wp+(0.6_wp-0.44_wp)*((x(2)+0d0*2.781d0)-17.0_wp)
!!$       elseif (((x(2)+0d0*2.781d0) .gt. 18.0_wp)) then
!!$             fd = 0.6_wp
!!$       end if
!!$
!!$
!!$       if ((x(2)+0d0*2.781d0) .le. 4.0_wp) then
!!$          Dc = 0.3_wp + 0.2_wp*(4.0_wp - (x(2)+0d0*2.781d0))
!!$       elseif ((x(2)+0d0*2.781d0) .ge. 4.0_wp) then
!!$          Dc = 0.3_wp
!!$       end if
!!$
!!$!       Dc = 0.3_wp                                                                                                                                                                                     
!!$
!!$       t0 = 0.5_wp !added in for tapered circle                                                                                                                                                         
!!$       Vs = 3.464_wp !added in for tapered circle                                                                                                                                                                                                                                                                                           
!!$       r = sqrt((1.3054_wp*((x(2)+0d0*2.781d0)-14.0_wp))**2 + (x(3)-38.0_wp)**2) !for ellipse, i.e. equiavelent circle
!!$       ! add high-strength barrier outside fault                                                                                                                                                        
!!$       ! Compute static frictional coefficient                                                                                                                                                          
!!$       fs = fs + 1.0e10*(1.0_wp-boxcar((x(2)+0d0*2.781d0),18.0_wp)*boxcar(x(3)-38.0_wp,20.0_wp))
!!$       ! Compute frictional cohesion                                                                                                                                                                    
!!$       call Frictional_Cohesion(C0, (x(2)+0d0*2.781d0), problem)
!!$
!!$       ! lock fault by increasing cohesion (and hence strength) outside of slipping region                                                                                                              
!!$       C0 = C0 + 1.0e10*(1.0_wp-boxcar((x(2)+0d0*2.781d0),18.0_wp)*boxcar(x(3)-38.0_wp,20.0_wp))
!!$       ! Compute frictional coefficient                                                                                                                                                                 
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                                                   
!!$       ! Compute time of forced rupture                                                                 \                                                                                               
!!$
!!$       call Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, problem)
!!$       f1 = min(D/Dc, 1.0_wp)
!!$       if (t < time_of_forced_rup) then
!!$          f2 = 0.0_wp
!!$       else if((t .ge. time_of_forced_rup) .and. (t < time_of_forced_rup+t0)) then
!!$          f2 = (t-time_of_forced_rup)/t0
!!$       else if(t .ge. time_of_forced_rup+t0) then
!!$          f2 = 1.0_wp
!!$       end if
!!$       f = fs-(fs-fd)*max(f1,f2) !for time of forced rupture weakening                                     
!!$       tau_str = C0 + max(sigma_n,0.0_wp)*f
    case default

       stop 'invalid problem in slip-weakening friction'

    end select

  end subroutine Tau_strength

  


  subroutine prestress(T0_n, T0_m, T0_l, x, t, problem, l, m, n, Mat)

    use datatypes, only : block_material

    implicit none

    real(kind = wp), intent(in) :: x(:),t
    real(kind = wp), intent(out) :: T0_n, T0_m, T0_l
    real(kind = wp), dimension(:), intent(in) :: n, m, l
    character(*),intent(in) :: problem
    real(kind = wp) :: Omega, b11, b33, b13, sigma22, Pf
    real(kind = wp) :: sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, u0
    real(kind = wp) :: Tx, Ty, Tz
    real(kind = wp) :: r, tau_nuke, C0
    real(kind = wp) :: F, G,tau_pulse,r0,s0(6),y0,z0
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp), optional, intent(in) :: Mat(3)

    select case(problem)

    case('TPV5','LOH1','OKLAHOMA')

       T0_n = -120.0_wp
       T0_m = 0.0_wp
       T0_l = 70.0_wp

       T0_l = T0_l+8.0_wp   *boxcar(x(2)-7.5_wp,1.5_wp)*boxcar(x(3)+7.5,1.5_wp)
       T0_l = T0_l+11.6_wp*boxcar(x(2)-7.5_wp,1.5_wp)*boxcar(x(3),1.5_wp)
       T0_l = T0_l-8.0_wp   *boxcar(x(2)-7.5_wp,1.5_wp)*boxcar(x(3)-7.5,1.5_wp)


    case('MG01')
       
       T0_n = -120.0_wp
       T0_m = 0.0_wp
       T0_l = 79.666667_wp

       T0_l = T0_l+1.675_wp*boxcar(x(2)-7.5_wp,0.734_wp)*boxcar(x(3),0.734_wp)

    case('MG2')
       
       T0_n = -120.0_wp
       T0_m = 0.0_wp
       T0_l = 69.111_wp
       
       T0_l = T0_l+12.283_wp*boxcar(x(2)-7.5_wp,1.4_wp)*boxcar(x(3),1.4_wp)
       

    case('TPV5_fractal')

       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       r = sqrt((x(1)-7.5_wp)**2 + x(2)**2)

       if (r <= 2.0_wp) then
          tau_nuke = 1.6_wp * 11.60_wp
       elseif(r >= 2.0_wp .and. r<=3.0_wp) then
          tau_nuke = 1.6_wp * 5.8_wp * (1.0_wp + cos(pi*(r-2.0_wp)/1.0_wp))
       else
          tau_nuke = 0.0_wp
       end if

       T0_l = T0_l + tau_nuke

    case('TPV10')
       
       sigma_xx = -7.378_wp*x(2)*(2.0_wp/sqrt(3.0_wp))
       sigma_yy = 0.0_wp
       sigma_zz = -7.378_wp*x(2)*(2.0_wp/sqrt(3.0_wp))
       sigma_xy = 0.0_wp
       sigma_xz = 0.0_wp
       sigma_yz = -0.55_wp*sigma_xx

       !call  initial_stress_tensor(s0,x,problem)
       !sigma_xx = s0(1)
       !sigma_yy = s0(2)
       !sigma_zz = s0(3)
       !sigma_xy = s0(4)
       !sigma_xz = s0(5)
       !sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                                                                                            
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz
       
       ! Compute tractions in local n, m, l coordinate system                                                                                                                                    
       T0_n = -7.378_wp*x(2)*(2.0_wp/(sqrt(3.0_wp)))
       T0_m = -0.55_wp*T0_n
       T0_l = 0.0_wp

       !(2.0_wp/sqrt(3.0_wp))
       ! Compute tractions in local n, m, l coordinate system                                                                                                                                  
       !T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       !T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       !T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       T0_m = T0_m + (-(0.760_wp+0.0057_wp)*T0_n+0.2_wp-T0_m)*boxcar(x(2)*(2.0_wp/sqrt(3.0_wp))-12.0_wp,1.5_wp)&
            *boxcar(x(3)-0.0_wp,1.5_wp)

    case('TPV36','TPV37')

       T0_n = -4.24_wp*x(2)/sin(15.0_wp/180.0_wp*pi)
       T0_m = -2.12_wp*x(2)/sin(15.0_wp/180.0_wp*pi)
       T0_l = 0.0_wp

          

    case('wasatch100')

!       call  initial_stress_tensor(s0,x,problem)                                                                                                                                                       
!       sigma_xx = s0(1)                                                                                                                                                                                
!       sigma_yy = s0(2)                                                                                                                                                                                
!       sigma_zz = s0(3)                                                                                                                                                                                
!       sigma_xy = s0(4)                                                                                                                                                                                
!       sigma_xz = s0(5)                                                                                                                                                                                
!       sigma_yz = s0(6)                                                                                                                                                                                

       ! Compute tractions in x,y,z coordinate system                                                                            \                                                                      

!       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz                                                                                                                                              
!       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz                                                                                                                                              
!       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz                                                                                                                                              

       ! Compute tractions in local n, m, l coordinate system                                                                    \                                                                      

!       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz                                                                                                                                                              
!       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz                                                                                                                                                              
!       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz                                                                                                                                                              

!TPV10 works                                                                                                                                                                                            
!       T0_n = -7.378_wp*1.154701_wp*x(2)                                                                                                                                                               
!       T0_m = 0.55_wp*7.378_wp*1.154701_wp*x(2)                                                                                                                                                        
!       T0_l = 0.0_wp                                                                                                                                                                                   
!TPV10 opposite                                                                                                                                                                                         
!       T0_n = -7.378_wp*1.154701_wp*x(2)                                                                                                                                                               
!       T0_m = -0.55_wp*7.378_wp*1.154701_wp*x(2)                                                                                                                                                       
!       T0_l = 0.0_wp                                                                                                                                                                                   

!       T0_n = -7.378_wp*x(2)                                                                                                                                                                           
!       T0_m = 0.55_wp*7.378_wp*x(2)                                                                                                                                                                    
!       T0_l = 0.0_wp                                                                                                                                                                                   

       !call Frictional_Cohesion(C0, x(2)+2.149_wp, problem)
       !call Frictional_Cohesion(C0, x(2)+2.781_wp, problem)

       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)

       ! Compute tractions in x,y,z coordinate system                                                                                  \                                                                

       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                          \                                                                

       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = (m(1)*Tx + m(2)*Ty + m(3)*Tz)
       T0_l = (l(1)*Tx + l(2)*Ty + l(3)*Tz)


       !T0_n = sigma_xx
       !T0_m = sigma_xy
       !T0_l = sigma_xz

       !T0_m = 0d0
!TPV10 works                                                                                                                                                                                            
       !T0_m=  T0_m+boxcar(x(2)+2.149_wp-10.392_wp,1.299_wp)*boxcar(x(3)-38.0_wp,1.5_wp)&
       !     *((C0+(0.760_wp+0.0057_wp)*7.378_wp*1.154701_wp*(x(2)+2.149_wp))-T0_m)                                !oppostie                             
       !T0_m=  T0_m+boxcar(x(2)-10.392_wp,1.299_wp)*boxcar(x(3),1.5_wp)*((C0+(0.760_wp+0.0057_wp)*-7.378_wp*1.154701_wp*x(2))-T0_m)                                                                     


       
    case('TPV5_2D')
       
       T0_n = -120.0_wp
       T0_m = 0.0_wp
       T0_l = 70.0_wp

       T0_l = T0_l+8.0_wp  * boxcar(x(2)+7.5,1.5_wp)
       T0_l = T0_l+11.6_wp * boxcar(x(2),1.5_wp)
       T0_l = T0_l-8.0_wp  * boxcar(x(2)-7.5,1.5_wp)


    case('TPV26')

       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)

       ! Compute tractions in x,y,z coordinate system                                                                                                         
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                                                 
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz
       

    case('TPV27')
          
       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                                                         
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                                                 
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz
       
    case('TPV28')
       
       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz
       
       ! Compute tractions in local n, m, l coordinate system
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       r = sqrt((x(2)-7.5_wp)**2 + x(3)**2)

       if (r <= 1.4_wp) then
          tau_nuke = 11.60_wp
       elseif(r >= 1.4_wp .and. r<=2.0_wp) then
          tau_nuke = 5.8_wp*(1.0_wp + cos(pi*(r-1.4_wp)/0.6_wp))
       else
          tau_nuke = 0.0_wp
       end if

       T0_l = T0_l + tau_nuke

    case('TPV29','TPV30')

       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)

       ! Compute tractions in x,y,z coordinate system                                                                                                                         
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                                                                 
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

    
    case('TPV31','TPV32')

       call  initial_stress_tensor(s0,x,problem,Mat)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       u0 = 32.03812032_wp

       ! Compute tractions in x,y,z coordinate system
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       r = sqrt((x(2)-7.5d0)**2 + x(3)**2)

       if (r <= 1.4_wp) then
          tau_nuke = 4.95_wp*Mat(2)/u0
       elseif(r >= 1.4_wp .and. r<=2_wp) then
          tau_nuke = 2.475_wp*(1_wp + cos(pi*(r-1.4_wp)/0.6_wp))*Mat(2)/u0
       else
          tau_nuke = 0_wp
       end if

       T0_l = T0_l + tau_nuke 

    case('TPV33')

       call  initial_stress_tensor(s0,x,problem,Mat)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)

       ! Compute tractions in x,y,z coordinate system
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz 

       r = sqrt((x(2)-6.0d0)**2 + (x(3)+6.0_wp)**2)

       if (r <= 0.55_wp) then
          tau_nuke = 3.15_wp
       elseif(r >= 0.55_wp .and. r<=0.8_wp) then
          tau_nuke = 1.575_wp*(1.0_wp + cos(pi*(r-0.55_wp)/0.25_wp))
       else
          tau_nuke = 0_wp
       end if

       T0_l = T0_l + tau_nuke 


    case('TPV34')

       call  initial_stress_tensor(s0,x,problem,Mat)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)

       ! Compute tractions in x,y,z coordinate system
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz 

       r = sqrt((x(2)-7.5d0)**2 + x(3)**2)

       u0 = 32.03812032_wp

       if (r <= 1.4_wp) then
          tau_nuke = 4.95_wp*Mat(2)/u0
       elseif(r >= 1.4_wp .and. r<=2.0_wp) then
          tau_nuke = 2.475_wp*(1.0_wp + cos(pi*(r-1.4_wp)/0.6_wp))*Mat(2)/u0

          !print*,'strike = ',x(2),'depth = ',x(1)
          !print*,'r_out = ',r,'u = ',Mat(3)*1e3_wp

       else
          tau_nuke = 0_wp
       end if

       T0_l = T0_l + tau_nuke  

    case('SCITS2016')

       call  initial_stress_tensor(s0,x,problem,Mat)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6) 

       ! Compute tractions in x,y,z coordinate system
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz 

       r = sqrt((x(2)-4.5_wp)**2.0_wp + (x(3)-20.0_wp)**2.0_wp)

       u0 = 32.03812032_wp

       if (r <= 1.4_wp) then
          tau_nuke = 14.95_wp*Mat(2)/u0
       elseif(r >= 1.4_wp .and. r<=2.0_wp) then
          tau_nuke = 7.475_wp*(1.0_wp + cos(pi*(r-1.4_wp)/0.6_wp))*Mat(2)/u0

          !print*,'strike = ',x(2),'depth = ',x(1)
          !print*,'r_out = ',r,'u = ',Mat(3)*1e3_wp

       else
          tau_nuke = 0.0_wp
       end if

       T0_l = T0_l + tau_nuke

    case('mg_a1')
       
       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                              
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                      
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       y0 = 7.5_wp
       z0 = -12.0_wp
       
       r = sqrt((x(2)-y0)**2.0_wp + (x(3) - z0)**2.0_wp)

       if (r <= 2.0_wp) then
          tau_nuke = 5.8_wp
       elseif(r >= 2.0_wp .and. r<=3.0_wp) then

          tau_nuke = 5.8_wp*0.5_wp*(1.0_wp + cos(pi*(r-2.0_wp)))

       else
          
          tau_nuke = 0.0_wp

       end if

       T0_l = T0_l + tau_nuke


        case('mg_b1a')
       
       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                              
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                      
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       y0 = 7.5_wp
       z0 = -12.0_wp
       
       r = sqrt((x(2)-y0)**2.0_wp + (x(3) - z0)**2.0_wp)

       if (r <= 2.0_wp) then
          tau_nuke = 7.5_wp
       elseif(r >= 2.0_wp .and. r<=3.0_wp) then

          tau_nuke = 14.0_wp*0.5_wp*(1.0_wp + cos(pi*(r-2.0_wp)))

       else
          
          tau_nuke = 0.0_wp

       end if

       T0_l = T0_l + tau_nuke

    case('TPV101')

       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                              
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                      
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz

       r = sqrt((x(2)-7.5_wp)**2 + x(3)**2)

       F = 0.0_wp
       if(r < 3.0_wp) F = exp(r**2/(r**2 - 9.0_wp))

       G = 0.0_wp
       if (t > 0.0_wp .and. t < 1.0_wp) G = exp((t-1.0_wp)**2/(t*(t-2.0_wp)))
       if (t .ge. 1.0_wp) G = 1.0_wp

       T0_l = T0_l + 25.0_wp*F*G

    case('TPV102')

       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                              
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz

       ! Compute tractions in local n, m, l coordinate system                                                                      
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz
       
       r = sqrt((x(2)-7.5_wp)**2 + x(3)**2)
       
       F = 0.0_wp
       if(r < 3.0_wp) F = exp(r**2/(r**2 - 9.0_wp))
       
       G = 0.0_wp
       if (t > 0.0_wp .and. t < 1.0_wp) G = exp((t-1.0_wp)**2/(t*(t-2.0_wp)))
       if (t .ge. 1.0_wp) G = 1.0_wp

       T0_l = T0_l + 25.0_wp*F*G


    case('rate-weakening')
       
       !call prestresstensor(sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, x(1), problem)
       call  initial_stress_tensor(s0,x,problem)
       sigma_xx = s0(1)
       sigma_yy = s0(2)
       sigma_zz = s0(3)
       sigma_xy = s0(4)
       sigma_xz = s0(5)
       sigma_yz = s0(6)
       
       ! Compute tractions in x,y,z coordinate system                                                                                     
       Tx = n(1)*sigma_xx + n(2)*sigma_xy + n(3)*sigma_xz
       Ty = n(1)*sigma_xy + n(2)*sigma_yy + n(3)*sigma_yz
       Tz = n(1)*sigma_xz + n(2)*sigma_yz + n(3)*sigma_zz
       
       ! Compute tractions in local n, m, l coordinate system                                                                             
       T0_n = n(1)*Tx + n(2)*Ty + n(3)*Tz
       T0_m = m(1)*Tx + m(2)*Ty + m(3)*Tz
       T0_l = l(1)*Tx + l(2)*Ty + l(3)*Tz
       
       tau_pulse = -0.2429_wp*T0_n
       
       r0 = 3.0_wp*1.5_wp*0.3_wp
       !r0 = 3.0_wp*0.3_wp

       y0 = 8.0_wp
       z0 = 0.0_wp
       
       r = sqrt((x(2) - y0)**2 + (x(3)+ z0)**2)
       
       F = exp(-((r/r0)**2)/2.0_wp)
       
       G = 0.0_wp
       
       if (t > 0.0_wp) G = 1.0_wp
       
       
       T0_l = T0_l + 2.0_wp*tau_pulse*F*G
       
       
    case default
       
       stop 'invalid problem in prestress'

    end select


  end subroutine Prestress


  elemental function boxcar(x,w) result(f)

    implicit none

    ! f(x) is boxcar of unit amplitude in (x-w,x+w)

    real(kind = wp),intent(in) :: x, w
    real(kind = wp) :: f

    real(kind = wp),parameter :: tol = 100.0_wp*epsilon(x)

    if (-w+tol<x.and.x<w-tol) then ! inside
       f = 1.0_wp
    elseif (abs(-w-x)<=tol.or.abs(x-w)<=tol) then ! boundary
       f = 0.5_wp
    else ! outside
       f = 0.0_wp
    end if

  end function boxcar


  elemental function  smooth_boxcar(x, w, wr) result(f)

    ! f(x) is a boxcar of unit amplitude in (x-w,x+w) with                                                       
    ! a smooth boundary layer of with wr: w<abs(x)<w+wr

    real(kind = wp),intent(in) :: x, w, wr
    real(kind = wp) :: f

    if (-w<=x .and. x<=w)  then ! inside
       f = 1.0_wp
    elseif (abs(x)>w .and. abs(x)< w+wr)  then ! boundary 
       f = 0.5_wp*(1.0_wp + tanh(wr/(abs(x)-w-wr) + wr/(abs(x)-w)))
    else
       f = 0.0_wp;
    end if

  end function smooth_boxcar

  elemental function arcsinh(x) result(f)

    real(kind = wp),intent(in) :: x
    real(kind = wp) :: f

    real(kind = wp) :: y

    y = abs(x)
    f = log(y+sqrt(y**2+1.0_wp))
    f = sign(f,x)

  end function arcsinh

  subroutine Frictional_Cohesion(C0, depth, problem)

    implicit none

    real(kind = wp), intent(in) :: depth
    character(*), intent(in) :: problem
    real(kind = wp), intent(out) :: C0


    select case(problem)

    case('TPV26')

       ! converted units from meter to kilometer

       if (depth .le. 5.0_wp) then

          C0 = 0.4_wp + 0.72_wp*(5.0_wp - depth)

       elseif (depth .ge. 5.0_wp) then

          C0 = 0.4_wp

       end if
       
       case('TPV27')

       ! converted units from meter to kilometer                                                                                                              

       if (depth .le. 5.0_wp) then

          C0 = 0.4_wp + 0.72_wp*(5.0_wp - depth)

       elseif (depth .ge. 5.0_wp) then

          C0 = 0.4_wp

       end if
       
    case('TPV29','TPV30')
       
       ! converted units from meter to kilometer                                                                        
       
       if (depth .le. 4.0_wp) then

          C0 = 0.4_wp + 0.2_wp*(4.0_wp - depth)
          
       elseif (depth .ge. 4.0_wp) then
          
          C0 = 0.4_wp

       end if
       
    
    case('TPV31','TPV32','TPV34','SCITS2016')

      ! converted units from meters to kilometers

      if (depth .le. 2.4_wp) then

          C0 = 0.425_wp*(2.4_wp - depth)

      elseif (depth .ge. 2.4_wp) then

          C0 = 0.0_wp

       endif


    case('mg_a1', 'mg_b1a')

       ! converted units from meter to kilometer                                                                        
       
       if (depth .le. 4.0_wp) then

          C0 = 0.4_wp + 0.5_wp*(4.0_wp - depth)
          
       elseif (depth .ge. 4.0_wp) then
          
          C0 = 0.4_wp

       end if

       case('wasatch100')

!          C0=0.2_wp                                                                                                                                                                                    
!added in frictional tapering of cohesion (larger near surface)                                                                                                                                         
       if (depth .le. 4.0_wp) then

          C0 = 0.4_wp + 0.2_wp*(4.0_wp - depth)

       elseif (depth .ge. 4.0_wp) then

          C0 = 0.4_wp

       end if

       case('TPV36')

       ! converted units from meter to kilometer                                                                                                                                                                                             

       if (depth .le. 8.0_wp) then

          C0 =  0.5_wp*(8.0_wp - depth)

       elseif (depth .ge. 8.0_wp) then

          C0 = 0.0_wp

       end if
       
    case('TPV37')

       ! converted units from meter to kilometer                                                                                      
                                                                                                                                        
       if (depth .le. 8.0_wp) then

          C0 =  1.875_wp*(8.0_wp - depth)

       elseif (depth .ge. 8.0_wp) then

          C0 = 0.0_wp

       end if   

    case default

       stop 'invalid problem in frictional cohesion'

    end select

  end subroutine Frictional_Cohesion



  subroutine Time_of_Forced_Rupture(time_of_forced_rup, r, Vs, problem)

    implicit none

    real(kind = wp), intent(in) :: r, Vs
    character(*), intent(in) :: problem
    real(kind = wp), intent(out) :: time_of_forced_rup
    real(kind = wp) :: r_crit


    select case(problem)

    case('TPV26')

       ! converted units from meter to kilometer and Kg/m^3 to g/cm^3
       r_crit = 4.0_wp

       if (r < r_crit-(1.0e-9_wp)) then

          time_of_forced_rup = r/(0.7_wp*Vs) + 0.081_wp*r_crit/(0.7_wp*Vs)*(1.0_wp/(1.0_wp - (r/r_crit)**2) - 1.0_wp)

       elseif (r .ge. r_crit-(1.0e-9_wp)) then

          time_of_forced_rup = 1.0e9_wp

       end if

       case('TPV27')

       ! converted units from meter to kilometer and Kg/m^3 to g/cm^3                                                                                         
       r_crit = 4.0_wp

       if (r < r_crit-(1d-9)) then

          time_of_forced_rup = r/(0.7_wp*Vs) + 0.081_wp*r_crit/(0.7_wp*Vs)*(1.0_wp/(1.0_wp - (r/r_crit)**2) - 1.0_wp)

       elseif (r .ge. r_crit-(1d-9)) then

          time_of_forced_rup = 1d9

       end if

       case('TPV29', 'TPV30')

       ! converted units from meter to kilometer and Kg/m^3 to g/cm^3                                                   
       r_crit = 4.0_wp

       if (r < r_crit-(1d-9)) then

          time_of_forced_rup = r/(0.7_wp*Vs) + 0.081_wp*r_crit/(0.7_wp*Vs)*(1.0_wp/(1.0_wp - (r/r_crit)**2) - 1.0_wp)
          
       elseif (r .ge. r_crit-(1.0e-9)) then
          
          time_of_forced_rup = 1.0e9
          
       end if
       
       case('TPV36', 'TPV37')

       ! converted units from meter to kilometer and Kg/m^3 to g/cm^3                                                                                                                                                                          
       r_crit = 4.0_wp

       if (r < r_crit-(10d-10)) then

          time_of_forced_rup = r/(0.7_wp*Vs) + 0.081_wp*r_crit/(0.7_wp*Vs)*(1.0_wp/(1.0_wp - (r/r_crit)**2) - 1.0_wp)

       elseif (r .ge. r_crit-(10d-10)) then

          time_of_forced_rup = 10d10

       end if

        case('wasatch100')

       ! converted units from meter to kilometer and Kg/m^3 to g/cm^3                                                                                                                                   
       r_crit = 4.0_wp

       if (r < r_crit-(1d-9)) then

          time_of_forced_rup = r/(0.7_wp*Vs) + 0.081_wp*r_crit/(0.7_wp*Vs)*(1.0_wp/(1.0_wp - (r/r_crit)**2) - 1.0_wp)

       elseif (r .ge. r_crit-(1.0e-9)) then

          time_of_forced_rup = 1.0e9

       end if
       
    
    case default

       stop 'invalid problem in time of forced rupture'

    end select

  end subroutine Time_of_Forced_Rupture
  
  subroutine Init_Vel_State(problem, II, BB, aa, ib)

    use datatypes, only : domain_type, block_type, iface_type
    
    character(256), intent(in) :: problem
    type(iface_type),intent(inout) :: II
    type(block_type),intent(inout) :: BB
    real(kind = wp), intent(in) :: aa
    integer, intent(in) :: ib
    
    integer :: mx, my, mz, px, py, pz !< grid lower and upper bounds
    
    real(kind = wp) :: T0_n, T0_m, T0_l                                ! tractions at a point lying on the fault plane
    real(kind = wp) :: x_jk(3)                                         ! coordinates of a point lying on the fault plane
    integer :: j,k                                        ! the grid indices
    real(kind = wp),dimension(1:3)::l,m,n                              ! unit vectors at a point lying on the fault plane
    real(kind = wp)::L0,f0,b ,V0,a,sigma_n,theta                             ! friction parameters
    real(kind = wp) :: Vin                                             ! initial background slip velocity on the fault
    integer :: xend
    
    mx = BB%G%C%mq
    my = BB%G%C%mr
    mz = BB%G%C%ms
    px = BB%G%C%pq
    py = BB%G%C%pr
    pz = BB%G%C%ps
    
    
    if (ib == 2) then
      xend = px
   else
      xend = mx
   end if

   select case(problem)
    case('TPV101')
       
       ! initialize particle velocities in the medium
       Vin = 1.0e-12_wp
      BB%F%F(mx:px, my:py, mz:pz, 3) = sign(0.5_wp*Vin, aa)
      
      ! initialize state variable on the fault
       do k=mz, pz
          do j= my, py
             
             x_jk(:) = BB%G%X(xend,j,k,:)
             
             n(:) = BB%B(ib)%n_n(j,k,:)
             m(:) = BB%B(ib)%n_m(j,k,:)
             l(:) = BB%B(ib)%n_l(j,k,:)

             call prestress(T0_n, T0_m, T0_l, x_jk, 0.0_wp, problem, l, m, n)
             L0 = 0.02_wp
             f0 = 0.6_wp
             b = 0.012_wp
             V0 = 1.0e-6
             a = 0.008_wp+0.008_wp*(1.0_wp-smooth_boxcar(x_jk(2),15.0_wp,3.0_wp)&
                  *smooth_boxcar(x_jk(3),15.0_wp,3.0_wp))
             sigma_n = -T0_n
             
             theta =L0/V0*exp(((a*log(2.0_wp*sinh(T0_l/(a*sigma_n))))-f0-a*log(Vin/V0))/b)
             II%W(j,k,1) = f0 + b*log(V0/L0*theta)
             
          end do
       end do
       
    case('TPV102')
       ! initialize particle velocities in the medium
       
       Vin = 1.0e-12_wp
       BB%F%F(mx:px, my:py, mz:pz, 3) = sign(0.5_wp*Vin, aa)
       
       ! initialize state variable on the fault
       do k= mz, pz
          do j= my, py
             
             x_jk(:) = BB%G%X(xend,j,k,:)
                          
             n(:) = BB%B(ib)%n_n(j,k,:)
             m(:) = BB%B(ib)%n_m(j,k,:)
             l(:) = BB%B(ib)%n_l(j,k,:)
             !         
             call prestress(T0_n, T0_m, T0_l, x_jk, 0.0_wp, problem, l, m, n)
             L0 = 0.02_wp
             f0 = 0.6_wp
             b = 0.012_wp
             V0 = 1.0e-6
             
             a = 0.008_wp+0.008_wp*(1.0_wp-smooth_boxcar(x_jk(2),15.0_wp,3.0_wp) &
                  *smooth_boxcar(x_jk(3),15.0_wp,3.0_wp))
             sigma_n = -T0_n
             
             theta =L0/V0*exp(((a*log(2.0_wp*sinh(T0_l/(a*sigma_n))))-f0-a*log(Vin/V0))/b)
             II%W(j,k,1) = f0 + b*log(V0/L0*theta)
             
          end do
       end do
       
    case('rate-weakening')
       ! initialize particle velocities in the medium                            
       
       Vin = 0.0_wp !1.0e-12_wp
       BB%F%F(mx:px, my:py, mz:pz, 3) = sign(0.5_wp*Vin, aa)
       
       ! initialize state variable on the fault    
       do k= mz, pz
          do j= my, py
             
             x_jk(:) = BB%G%X(xend,j,k,:)
                          
             n(:) = BB%B(ib)%n_n(j,k,:)
             m(:) = BB%B(ib)%n_m(j,k,:)
             l(:) = BB%B(ib)%n_l(j,k,:)
             
             call prestress(T0_n, T0_m, T0_l, x_jk, 0.0_wp, problem, l, m, n)
             L0 = 2.0_wp*0.2572_wp
             f0 = 0.7_wp
             b =  0.02_wp
             V0 = 1.0e-6_wp
             
             theta = 0.4367_wp
             II%W(j,k,1) = theta ! f0 + b*log(V0/L0*theta)
             
          end do
       end do
       
    end select
    
  end subroutine Init_Vel_State
  

  subroutine Set_Theta(W, Theta, problem)
    ! Sets the state variable Theta                                                                                                             
    real(kind = wp),dimension(:,:),intent(in):: W                   ! the state variable                                                                   
    real(kind = wp),dimension(:,:),intent(inout):: Theta            ! the state variable                                                                   
    character(256),intent(in)::problem                   ! the problem type                                                                     
    real(kind = wp)::L0,f0,b,V0                             ! friction parameters                                                                          

    select case(problem)
       
    case ('TPV101')
       L0 = 0.02_wp
       f0 = 0.6_wp
       b = 0.012_wp
       V0 = 1.0e-6_wp
       Theta = L0/V0*exp((W-f0)/b)
       
    case ('TPV102')
       
       L0 = 0.02_wp
       f0 = 0.6_wp
       b = 0.012_wp
       V0 = 1.0e-6_wp
       Theta= L0/V0*exp((W-f0)/b)

    case('rate-weakening')
       L0 = 2.0_wp*0.2572_wp
       f0 = 0.7_wp
       b =  0.02_wp
       V0 = 1.0e-6_wp
       Theta= L0/V0*exp((W-f0)/b)
       
    end select
  end subroutine Set_Theta

end module Interface_Condition


