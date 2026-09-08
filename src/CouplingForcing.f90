module CouplingForcing

  use common, only : wp
  implicit none

contains


  subroutine Couple_Interface_x(u_bc, v_bc, problem, coupling, I, G, M, U, V, &
                            LambdaMu1, LambdaMu2, n_l, n_m, n_n, ib, &
                            t, stage, mms_vars, handles)

    use common, only : wp
    use datatypes, only : block_grid_t, block_boundary, iface_type, mms_type,block_material
    use mms, only : mms_type,evaluate_mms
    use Interface_Condition, only :  Interface_RHS
    use fault_output, only : fault_type, write_hats
    use mpi3dbasic, only:  rank
    
    ! compute hat variables and SAT forcing terms, but does not add SAT terms to rates
  
    implicit none
    real(kind = wp),dimension(:,:,:), allocatable, intent(inout) :: u_bc, v_bc                           ! SAT-forcing
    character(256), intent(in) :: problem
    character(*),intent(in) :: coupling
    integer, intent(in) :: stage, ib
    type(iface_type),intent(inout) :: I
    type(block_grid_t), intent(in) :: G
    type(block_material), intent(in) :: M
    type(mms_type), intent(inout) :: mms_vars
    type(fault_type), intent(inout) :: handles
    real(kind = wp),intent(in) :: t

    integer :: mbx1, mby1, mbz1, pbx1, pby1, pbz1
    real(kind = wp),dimension(:,:,:), allocatable, intent(in) :: u                  ! Field variables
    real(kind = wp),dimension(:,:,:), allocatable, intent(in) :: v                  ! Field variables
    real(kind = wp),dimension(:,:,:), allocatable, intent(in) :: LambdaMu1                    ! Material parameters
    real(kind = wp),dimension(:,:,:), allocatable, intent(in) :: LambdaMu2                    ! Material parameters
    real(kind = wp),dimension(:,:,:), allocatable, intent(in) :: n_l, n_m, n_n           !< normal vectors

    integer :: mx1, my1, mz1, px1, py1, pz1

    real(kind = wp), dimension(:, :, :), allocatable, save :: u_rotated, v_rotated   ! Local grid varibles in rotated cordinates
    real(kind = wp), dimension(:, :, :), allocatable, save :: u_hat, v_hat           ! Local hat variables

    ! Local characteristics in local rotated grid variables
    real(kind = wp) :: p1_p, p2_p, p3_p
    real(kind = wp) :: q1_m, q2_m, q3_m

    ! Local characteristics in local rotated hat variables
    real(kind = wp) :: p1_phat, p2_phat, p3_phat
    real(kind = wp) :: q1_mhat, q2_mhat, q3_mhat

    real(kind = wp) :: c1_p, c1_s, c2_p, c2_s                                     ! Local velocites
    real(kind = wp) :: Tu_x, Tu_y, Tu_z                                           ! Local stresses in the left block
    real(kind = wp) :: Tv_x, Tv_y, Tv_z                                           ! Local stresses in the right block

    real(kind = wp) :: u_x, u_y, u_z                                              ! Local particle velocities in the left block
    real(kind = wp) :: v_x, v_y, v_z                                              ! Local particle velocites in the right block

    integer :: y, z, xend                                                    ! Local grid indices

    real(kind = wp), dimension(:), allocatable, save :: fu(:), fv(:) ! exact MMS solutions

    real(kind = wp) :: norm1, norm2, dir
    real(kind = wp) :: q1_x, q1_y, q1_z, q2_x, q2_y, q2_z
    real(kind = wp) :: Z1_p, Z1_s, Z2_p, Z2_s, rho1, rho2
    real(kind = wp) :: lambda1, lambda2, mu1, mu2
    !real :: Z_p, Z_s
    real(kind = wp) :: n_1, n_2, n_3, m_1, m_2, m_3, l_1, l_2, l_3
    real(kind = wp) :: m11, m12, m13, m21, m22, m23, m31, m32, m33
    real(kind = wp) :: p11, p12, p13, p21, p22, p23, p31, p32, p33
    real(kind = wp) :: n11, n12, n13, n21, n22, n23, n31, n32, n33
    real(kind = wp) :: q11, q12, q13, q21, q22, q23, q31, q32, q33
    !======================================================


    mx1 = G%C%mq
    px1 = G%C%pq
    my1 = G%C%mr
    py1 = G%C%pr
    mz1 = G%C%ms
    pz1 = G%C%ps

    mbx1 = G%C%mbq
    pbx1 = G%C%pbq
    mby1 = G%C%mbr
    pby1 = G%C%pbr
    mbz1 = G%C%mbs
    pbz1 = G%C%pbs

    if (.not. allocated(fu)) allocate(fu(9))
    if (.not. allocated(fv)) allocate(fv(9))
    
    if (.not. allocated(u_rotated)) allocate(u_rotated(my1:py1, mz1:pz1, 6), v_rotated(my1:py1, mz1:pz1, 6))
    if (.not. allocated(u_hat)) allocate(u_hat(my1:py1, mz1:pz1, 6), v_hat(my1:py1, mz1:pz1, 6))

    !> Block Boundary Fields
    if (ib == 2) then
      xend = px1
    else
      xend = mx1
    end if

    dir = 1.0_wp*I%II%normal(1)


    do z = mz1, pz1
       do y = my1, py1

          !==================================================================================
          ! compute tractions  and particle velocities (components in x,y,z coordinates)
          !==================================================================================

          Tu_x =  n_n(y, z, 1)*(u(y, z, 4)) &
               + n_n(y, z, 2)*(u(y, z, 7))&
               + n_n(y, z, 3)*(u(y, z, 8))

          Tu_y =  n_n(y, z, 1)*(u(y, z, 7))&
               + n_n(y, z, 2)*(u(y, z, 5)) &
               + n_n(y, z, 3)*(u(y, z, 9))

          Tu_z =  n_n(y, z, 1)*(u(y, z, 8)) &
               + n_n(y, z, 2)*(u(y, z, 9)) &
               + n_n(y, z, 3)*(u(y, z, 6))

          Tv_x =  n_n(y, z, 1)*(v(y, z, 4))&
               + n_n(y, z, 2)*(v(y, z, 7))&
               + n_n(y, z, 3)*(v(y, z, 8))

          Tv_y =  n_n(y, z, 1)*(v(y, z, 7))&
               + n_n(y, z, 2)*(v(y, z, 5))&
               + n_n(y, z, 3)*(v(y, z, 9))

          Tv_z =  n_n(y, z, 1)*(v(y, z, 8))&
               + n_n(y, z, 2)*(v(y, z, 9))&
               + n_n(y, z, 3)*(v(y, z, 6))

          u_x = u(y, z, 1)
          u_y = u(y, z, 2)
          u_z = u(y, z, 3)

          v_x = v(y, z, 1)
          v_y = v(y, z, 2)
          v_z = v(y, z, 3)

          if (mms_vars%use_mms) then ! forcing

            !> @TODO fix pointers to grids
             !fu = evaluate_mms(t, XX1(px1, y, z, :), mms_vars)
             !fv = evaluate_mms(t, XX2(mx1, y, z, :), mms_vars)

             Tu_x =  Tu_x - (n_n(y, z, 1)*fu(4)&
                  + n_n(y, z, 2)*fu(7)&
                  + n_n(y, z, 3)*fu(8))

             Tu_y =  Tu_y - (n_n(y, z, 1)*fu(7)&
                  + n_n(y, z, 2)*fu(5)&
                  + n_n(y, z, 3)*fu(9))

             Tu_z =  Tu_z - (n_n(y, z, 1)*fu(8)&
                  + n_n(y, z, 2)*fu(9)&
                  + n_n(y, z, 3)*fu(6))

             Tv_x = Tv_x - (n_n(y, z, 1)*fv(4)&
                  + n_n(y, z, 2)*fv(7)&
                  + n_n(y, z, 3)*fv(8))

             Tv_y = Tv_y - (n_n(y, z, 1)*fv(7)&
                  + n_n(y, z, 2)*fv(5)&
                  + n_n(y, z, 3)*fv(9))

             Tv_z = Tv_z - (n_n(y, z, 1)*fv(8)&
                  + n_n(y, z, 2)*fv(9)&
                  + n_n(y, z, 3)*fv(6))

             u_x =  u_x - fu(1)
             u_y =  u_y - fu(2)
             u_z =  u_z - fu(3)

             v_x =  v_x - fv(1)
             v_y =  v_y - fv(2)
             v_z =  v_z - fv(3)

          end if

          ! Initialize variables to be rotated below
          u_rotated(y, z, 1) = u_x
          u_rotated(y, z, 2) = u_y
          u_rotated(y, z, 3) = u_z
          u_rotated(y, z, 4) = Tu_x
          u_rotated(y, z, 5) = Tu_y
          u_rotated(y, z, 6) = Tu_z

          v_rotated(y, z, 1) = v_x
          v_rotated(y, z, 2) = v_y
          v_rotated(y, z, 3) = v_z
          v_rotated(y, z, 4) = Tv_x
          v_rotated(y, z, 5) = Tv_y
          v_rotated(y, z, 6) = Tv_z

       end do
    end do

    ! Rotate in local (n,m,l) orthogonal coordinates
    call Rotate_in_Local_Cordinates(u_rotated, v_rotated, n_l, n_m, n_n, my1, mz1, py1, pz1)

    !> Compute hat-variables: u_hat, v_hat (in n,m,l coordinates)

    call Interface_RHS(problem, coupling, t, stage, handles, u_rotated, v_rotated, &
                       u_hat, v_hat, G, M, LambdaMu1, LambdaMu2, n_l, n_m, n_n, I, ib)


    !> @todo need an explanation from KD on what is going on here.

    if(stage == 1) then

       ! only used for output

       I%Svel(my1:py1,mz1:pz1,1:3) = v_hat(my1:py1,mz1:pz1,1:3) - u_hat(my1:py1,mz1:pz1,1:3)

       call rupture_front(my1, py1, mz1, pz1, I%Svel, I%trup, t)
       handles%slip(my1:py1,mz1:pz1,:) = I%S(my1:py1,mz1:pz1,:)
       handles%time_rup(my1:py1,mz1:pz1,1) = I%trup(my1:py1,mz1:pz1,1)
       handles%Uhat_pluspres(my1:py1,mz1:pz1,1:3) = v_hat(my1:py1,mz1:pz1,1:3) - u_hat(my1:py1,mz1:pz1,1:3)
       handles%Vhat_pluspres(my1:py1,mz1:pz1,1:3) = v_hat(my1:py1,mz1:pz1,1:3) - u_hat(my1:py1,mz1:pz1,1:3)

    end if


    !if (write_interface) call write_hats(u_hat,v_hat,I%Svel,I%trup,handles)
    ! Rotate back to x,y,z coordinates
    call Rotate_Back_Physical_Cordinates(u_rotated, v_rotated, n_l, n_m, n_n, my1, mz1, py1, pz1)
    call Rotate_Back_Physical_Cordinates(u_hat, v_hat, n_l, n_m, n_n, my1, mz1, py1, pz1)

    ! enforce interface conditions, in local orthogonal coordinate system, and set hat variables

    do z =mz1, pz1
       do y = my1, py1

          ! something related to unit normal?
          q1_x = G%metricx(xend, y, z, 1)
          q1_y = G%metricy(xend, y, z, 1)
          q1_z = G%metricz(xend, y, z, 1)

          q2_x = G%metricx(xend, y, z, 1)
          q2_y = G%metricy(xend, y, z, 1)
          q2_z = G%metricz(xend, y, z, 1)

          norm1 = sqrt(q1_x**2 + q1_y**2 + q1_z**2)
          norm2 = sqrt(q2_x**2 + q2_y**2 + q2_z**2)

          !==================================================================================
          ! Left Block
          !==================================================================================

          lambda1 = LambdaMu1(y, z, 1)                                                         ! bulk modulus
          mu1     = LambdaMu1(y, z, 2)                                                         ! shear modulus
          rho1    = LambdaMu1(y, z, 3)                                                         ! density

          lambda2 = LambdaMu2(y, z, 1)                                                           ! bulk modulus
          mu2     = LambdaMu2(y, z, 2)                                                           ! shear modulus
          rho2    = LambdaMu2(y, z, 3)                                                           ! density

          c1_p = sqrt((2.0_wp*mu1 + lambda1)/rho1)                                                     ! P-wave speed
          c1_s = sqrt(mu1/rho1)                                                                     ! S-wave speed
          Z1_p = rho1*c1_p                                                                          ! P-wave impedance
          Z1_s = rho1*c1_s                                                                          ! S-wave impedance

          c2_p = sqrt((2.0_wp*mu2 + lambda2)/rho2)                                                     ! P-wave speed
          c2_s = sqrt(mu2/rho2)                                                                     ! S-wave speed
          Z2_p = rho2*c2_p                                                                          ! P-wave impedance
          Z2_s = rho2*c2_s                                                                          ! S-wave impedance


          ! Take the average of the impedances on the interface (to avoid division by zero)
          !if (min(Z1_s,Z2_s)/max(Z1_s,Z2_s) .le. 1d-5) then
             ! Take the arithmetic average
           !  Z_p = 0.5_wp*(Z1_p + Z2_p)
           !  Z_s = 0.5_wp*(Z1_s + Z2_s)  ! DO NOT CHANGE THIS

          !else
             ! Take the harmonic average
             !Z_p = 2.0_wp*Z1_p*Z2_p/(Z1_p + Z2_p)
             !Z_s = 2.0_wp*Z1_s*Z2_s/(Z1_s + Z2_s)

          !end if

          !==================================================================================
          ! Left Block
          !==================================================================================


          !==================================================================================
          ! Construct Characteristics
          !==================================================================================

          ! Unit normals
          n_1 = n_n(y, z, 1)
          n_2 = n_n(y, z, 2)
          n_3 = n_n(y, z, 3)

          m_1 = n_m(y, z, 1)
          m_2 = n_m(y, z, 2)
          m_3 = n_m(y, z, 3)

          l_1 = n_l(y, z, 1)
          l_2 = n_l(y, z, 2)
          l_3 = n_l(y, z, 3)

          ! Compute M = S^TZ1S
          m11 = n_1*n_1*Z1_p + n_2*n_2*Z1_s + n_3*n_3*Z1_s
          m12 = n_1*m_1*Z1_p + n_2*m_2*Z1_s + n_3*m_3*Z1_s
          m13 = n_1*l_1*Z1_p + n_2*l_2*Z1_s + n_3*l_3*Z1_s

          m21 = n_1*m_1*Z1_p + n_2*m_2*Z1_s + n_3*m_3*Z1_s
          m22 = m_1*m_1*Z1_p + m_2*m_2*Z1_s + m_3*m_3*Z1_s
          m23 = m_1*l_1*Z1_p + m_2*l_2*Z1_s + m_3*l_3*Z1_s

          m31 = l_1*n_1*Z1_p + l_2*n_2*Z1_s + l_3*n_3*Z1_s
          m32 = l_1*m_1*Z1_p + l_2*m_2*Z1_s + l_3*m_3*Z1_s
          m33 = l_1*l_1*Z1_p + l_2*l_2*Z1_s + l_3*l_3*Z1_s

          ! Compute P = S^TZ1^(-1)S
          p11 = n_1*n_1*1.0_wp/Z1_p + n_2*n_2*1.0_wp/Z1_s + n_3*n_3*1.0_wp/Z1_s
          p12 = n_1*m_1*1.0_wp/Z1_p + n_2*m_2*1.0_wp/Z1_s + n_3*m_3*1.0_wp/Z1_s
          p13 = n_1*l_1*1.0_wp/Z1_p + n_2*l_2*1.0_wp/Z1_s + n_3*l_3*1.0_wp/Z1_s

          p21 = n_1*m_1*1.0_wp/Z1_p + n_2*m_2*1.0_wp/Z1_s + n_3*m_3*1.0_wp/Z1_s
          p22 = m_1*m_1*1.0_wp/Z1_p + m_2*m_2*1.0_wp/Z1_s + m_3*m_3*1.0_wp/Z1_s
          p23 = m_1*l_1*1.0_wp/Z1_p + m_2*l_2*1.0_wp/Z1_s + m_3*l_3*1.0_wp/Z1_s
          p31 = l_1*n_1*1.0_wp/Z1_p + l_2*n_2*1.0_wp/Z1_s + l_3*n_3*1.0_wp/Z1_s
          p32 = l_1*m_1*1.0_wp/Z1_p + l_2*m_2*1.0_wp/Z1_s + l_3*m_3*1.0_wp/Z1_s
          p33 = l_1*l_1*1.0_wp/Z1_p + l_2*l_2*1.0_wp/Z1_s + l_3*l_3*1.0_wp/Z1_s

          ! The Hat-variable characteristics for the left block
          p1_phat = 0.5_wp*(m11*u_hat(y, z, 1) + m12*u_hat(y, z, 2) &
               + m13*u_hat(y, z, 3) + u_hat(y, z, 4))
          p2_phat = 0.5_wp*(m21*u_hat(y, z, 1) + m22*u_hat(y, z, 2) &
               + m23*u_hat(y, z, 3) + u_hat(y, z, 5))
          p3_phat = 0.5_wp*(m31*u_hat(y, z, 1) + m32*u_hat(y, z, 2) &
               + m33*u_hat(y, z, 3) + u_hat(y, z, 6))

          ! The Grid-function characteristics for the left block
          p1_p = 0.5_wp*(m11*u_rotated(y, z, 1) + m12*u_rotated(y, z, 2) &
               + m13*u_rotated(y, z, 3) + u_rotated(y, z, 4))
          p2_p = 0.5_wp*(m21*u_rotated(y, z, 1) + m22*u_rotated(y, z, 2) &
               + m23*u_rotated(y, z, 3) + u_rotated(y, z, 5))
          p3_p = 0.5_wp*(m31*u_rotated(y, z, 1) + m32*u_rotated(y, z, 2) &
               + m33*u_rotated(y, z, 3) + u_rotated(y, z, 6))



          !==================================================================================
          ! Right block
          ! =================================================================================


          ! Compute N = S^TZ2S
          n11 = n_1*n_1*Z2_p + n_2*n_2*Z2_s + n_3*n_3*Z2_s
          n12 = n_1*m_1*Z2_p + n_2*m_2*Z2_s + n_3*m_3*Z2_s
          n13 = n_1*l_1*Z2_p + n_2*l_2*Z2_s + n_3*l_3*Z2_s

          n21 = n_1*m_1*Z2_p + n_2*m_2*Z2_s + n_3*m_3*Z2_s
          n22 = m_1*m_1*Z2_p + m_2*m_2*Z2_s + m_3*m_3*Z2_s
          n23 = m_1*l_1*Z2_p + m_2*l_2*Z2_s + m_3*l_3*Z2_s

          n31 = l_1*n_1*Z2_p + l_2*n_2*Z2_s + l_3*n_3*Z2_s
          n32 = l_1*m_1*Z2_p + l_2*m_2*Z2_s + l_3*m_3*Z2_s
          n33 = l_1*l_1*Z2_p + l_2*l_2*Z2_s + l_3*l_3*Z2_s

          ! Compute Q = S^TZ2^(-1)S
          q11 = n_1*n_1*1.0_wp/Z2_p + n_2*n_2*1.0_wp/Z2_s + n_3*n_3*1.0_wp/Z2_s
          q12 = n_1*m_1*1.0_wp/Z2_p + n_2*m_2*1.0_wp/Z2_s + n_3*m_3*1.0_wp/Z2_s
          q13 = n_1*l_1*1.0_wp/Z2_p + n_2*l_2*1.0_wp/Z2_s + n_3*l_3*1.0_wp/Z2_s

          q21 = n_1*m_1*1.0_wp/Z2_p + n_2*m_2*1.0_wp/Z2_s + n_3*m_3*1.0_wp/Z2_s
          q22 = m_1*m_1*1.0_wp/Z2_p + m_2*m_2*1.0_wp/Z2_s + m_3*m_3*1.0_wp/Z2_s
          q23 = m_1*l_1*1.0_wp/Z2_p + m_2*l_2*1.0_wp/Z2_s + m_3*l_3*1.0_wp/Z2_s

          q31 = l_1*n_1*1.0_wp/Z2_p + l_2*n_2*1.0_wp/Z2_s + l_3*n_3*1.0_wp/Z2_s
          q32 = l_1*m_1*1.0_wp/Z2_p + l_2*m_2*1.0_wp/Z2_s + l_3*m_3*1.0_wp/Z2_s
          q33 = l_1*l_1*1.0_wp/Z2_p + l_2*l_2*1.0_wp/Z2_s + l_3*l_3*1.0_wp/Z2_s

          !==================================================================================
          ! Construct Characteristics
          !==================================================================================

          ! The Hat-variables characteristics for the right block
          q1_mhat = 0.5_wp*(n11*v_hat(y, z, 1) + n12*v_hat(y, z, 2) &
               + n13*v_hat(y, z, 3) - v_hat(y, z, 4))
          q2_mhat = 0.5_wp*(n21*v_hat(y, z, 1) + n22*v_hat(y, z, 2) &
               + n23*v_hat(y, z, 3) - v_hat(y, z, 5))
          q3_mhat = 0.5_wp*(n31*v_hat(y, z, 1) + n32*v_hat(y, z, 2) &
               + n33*v_hat(y, z, 3) - v_hat(y, z, 6))

          ! The Grid-function characteristics for the right block
          q1_m = 0.5_wp*(n11*v_rotated(y, z, 1) + n12*v_rotated(y, z, 2) &
               + n13*v_rotated(y, z, 3) - v_rotated(y, z, 4))
          q2_m = 0.5_wp*(n21*v_rotated(y, z, 1) + n22*v_rotated(y, z, 2) &
               + n23*v_rotated(y, z, 3) - v_rotated(y, z, 5))
          q3_m = 0.5_wp*(n31*v_rotated(y, z, 1) + n32*v_rotated(y, z, 2) &
               + n33*v_rotated(y, z, 3) - v_rotated(y, z, 6))


          ! ================================================================
          ! Construct SAT forcing for the left block
          ! ================================================================

          u_bc(y, z, 1) = norm1/rho1*(p1_p - p1_phat)
          u_bc(y, z, 2) = norm1/rho1*(p2_p - p2_phat)
          u_bc(y, z, 3) = norm1/rho1*(p3_p - p3_phat)

          u_bc(y, z, 4) = (2.0_wp*mu1+lambda1)*q1_x*(p11*(p1_p - p1_phat)+p12*(p2_p - p2_phat)+p13*(p3_p - p3_phat)) &
                + lambda1*q1_y*(p21*(p1_p - p1_phat)+p22*(p2_p - p2_phat)+p23*(p3_p - p3_phat))&
                + lambda1*q1_z*(p31*(p1_p - p1_phat)+p32*(p2_p - p2_phat)+p33*(p3_p - p3_phat))

          u_bc(y, z, 5) = (2.0_wp*mu1+lambda1)*q1_y*(p21*(p1_p - p1_phat)+p22*(p2_p - p2_phat)+p23*(p3_p - p3_phat)) &
               + lambda1*q1_z*(p31*(p1_p - p1_phat)+p32*(p2_p - p2_phat)+p33*(p3_p - p3_phat))&
               + lambda1*q1_x*(p11*(p1_p - p1_phat)+p12*(p2_p - p2_phat)+p13*(p3_p - p3_phat))

          u_bc(y, z, 6) = (2.0_wp*mu1+lambda1)*q1_z*(p31*(p1_p - p1_phat)+p32*(p2_p - p2_phat)+p33*(p3_p - p3_phat)) &
               + lambda1*q1_y*(p21*(p1_p - p1_phat)+p22*(p2_p - p2_phat)+p23*(p3_p - p3_phat))&
               + lambda1*q1_x*(p11*(p1_p - p1_phat)+p12*(p2_p - p2_phat)+p13*(p3_p - p3_phat))

          u_bc(y, z, 7) = mu1*(q1_y*(p11*(p1_p - p1_phat)+p12*(p2_p - p2_phat)+p13*(p3_p - p3_phat))&
               + q1_x*(p21*(p1_p - p1_phat)+p22*(p2_p - p2_phat)+p23*(p3_p - p3_phat)))

          u_bc(y, z, 8) = mu1*(q1_z*(p11*(p1_p - p1_phat)+p12*(p2_p - p2_phat)+p13*(p3_p - p3_phat))&
               + q1_x*(p31*(p1_p - p1_phat)+p32*(p2_p - p2_phat)+p33*(p3_p - p3_phat)))

          u_bc(y, z, 9) = mu1*(q1_z*(p21*(p1_p - p1_phat)+p22*(p2_p - p2_phat)+p23*(p3_p - p3_phat))&
               + q1_y*(p31*(p1_p - p1_phat)+p32*(p2_p - p2_phat)+p33*(p3_p - p3_phat)))


          ! ================================================================
          ! Construct SAT forcing for the right block
          ! ================================================================

          v_bc(y, z, 1) = norm2/rho2*(q1_m - q1_mhat)
          v_bc(y, z, 2) = norm2/rho2*(q2_m - q2_mhat)
          v_bc(y, z, 3) = norm2/rho2*(q3_m - q3_mhat)

          v_bc(y, z, 4) = -(2.0_wp*mu2+lambda2)*q2_x*(q11*(q1_m - q1_mhat)+q12*(q2_m - q2_mhat)+q13*(q3_m - q3_mhat)) &
                - lambda2*q2_y*(q21*(q1_m - q1_mhat)+q22*(q2_m - q2_mhat)+q23*(q3_m - q3_mhat))&
                - lambda2*q2_z*(q31*(q1_m - q1_mhat)+q32*(q2_m - q2_mhat)+q33*(q3_m - q3_mhat))

          v_bc(y, z, 5) = -(2.0_wp*mu2+lambda2)*q2_y*(q21*(q1_m - q1_mhat)+q22*(q2_m - q2_mhat)+q23*(q3_m - q3_mhat)) &
               - lambda2*q2_z*(q31*(q1_m - q1_mhat)+q32*(q2_m - q2_mhat)+q33*(q3_m - q3_mhat))&
               - lambda2*q2_x*(q11*(q1_m - q1_mhat)+q12*(q2_m - q2_mhat)+q13*(q3_m - q3_mhat))

          v_bc(y, z, 6) = -(2.0_wp*mu2+lambda2)*q2_z*(q31*(q1_m - q1_mhat)+q32*(q2_m - q2_mhat)+q33*(q3_m - q3_mhat)) &
               - lambda2*q2_y*(q21*(q1_m - q1_mhat)+q22*(q2_m - q2_mhat)+q23*(q3_m - q3_mhat))&
               - lambda2*q2_x*(q11*(q1_m - q1_mhat)+q12*(q2_m - q2_mhat)+q13*(q3_m - q3_mhat))

          v_bc(y, z, 7) = -mu2*(q2_y*(q11*(q1_m - q1_mhat)+q12*(q2_m - q2_mhat)+q13*(q3_m - q3_mhat))&
               + q2_x*(q21*(q1_m - q1_mhat)+q22*(q2_m - q2_mhat)+q23*(q3_m - q3_mhat)))

          v_bc(y, z, 8) = -mu2*(q2_z*(q11*(q1_m - q1_mhat)+q12*(q2_m - q2_mhat)+q13*(q3_m - q3_mhat))&
               + q2_x*(q31*(q1_m - q1_mhat)+q32*(q2_m - q2_mhat)+q33*(q3_m - q3_mhat)))

          v_bc(y, z, 9) = -mu2*(q2_z*(q21*(q1_m - q1_mhat)+q22*(q2_m - q2_mhat)+q23*(q3_m - q3_mhat))&
               + q2_y*(q31*(q1_m - q1_mhat)+q32*(q2_m - q2_mhat)+q33*(q3_m - q3_mhat)))

       end do
    end do
  end subroutine Couple_Interface_x


  subroutine Rotate_in_Local_Cordinates(u, v, X_l1, X_m1, X_n1, my1, mz1, py1, pz1)

    implicit none

    ! Rotate field variables in local basis vectors

    integer, intent(in) :: my1, mz1, py1, pz1 !< Grid length
    real(kind = wp), dimension(:,:,:), allocatable, intent(in) :: X_l1, X_m1, X_n1 !< Local basis vectors
    real(kind = wp), dimension(:,:,:), allocatable, intent(inout) :: u, v !< Field variables

    real(kind = wp) :: u1, u2, u3, v1, v2, v3, Tu_x, Tu_y, Tu_z, Tv_x, Tv_y, Tv_z !< Local work variables
    integer :: y, z !< Local grid indices

    do z = mz1, pz1
       do y = my1, py1

          ! Set local grid functions

          u1 = u(y, z, 1)
          u2 = u(y, z, 2)
          u3 = u(y, z, 3)

          v1 = v(y, z, 1)
          v2 = v(y, z, 2)
          v3 = v(y, z, 3)

          Tu_x = u(y, z, 4)
          Tu_y = u(y, z, 5)
          Tu_z = u(y, z, 6)

          Tv_x = v(y, z, 4)
          Tv_y = v(y, z, 5)
          Tv_z = v(y, z, 6)

          ! Rotate in local basis vectors

          u(y, z, 1) = X_n1(y, z, 1)*u1 + X_n1(y, z, 2)*u2 + X_n1(y, z, 3)*u3
          u(y, z, 2) = X_m1(y, z, 1)*u1 + X_m1(y, z, 2)*u2 + X_m1(y, z, 3)*u3
          u(y, z, 3) = X_l1(y, z, 1)*u1 + X_l1(y, z, 2)*u2 + X_l1(y, z, 3)*u3

          u(y, z, 4) = X_n1(y, z, 1)*Tu_x + X_n1(y, z, 2)*Tu_y + X_n1(y, z, 3)*Tu_z
          u(y, z, 5) = X_m1(y, z, 1)*Tu_x + X_m1(y, z, 2)*Tu_y + X_m1(y, z, 3)*Tu_z
          u(y, z, 6) = X_l1(y, z, 1)*Tu_x + X_l1(y, z, 2)*Tu_y + X_l1(y, z, 3)*Tu_z

          v(y, z, 1) = X_n1(y, z, 1)*v1 + X_n1(y, z, 2)*v2 + X_n1(y, z, 3)*v3
          v(y, z, 2) = X_m1(y, z, 1)*v1 + X_m1(y, z, 2)*v2 + X_m1(y, z, 3)*v3
          v(y, z, 3) = X_l1(y, z, 1)*v1 + X_l1(y, z, 2)*v2 + X_l1(y, z, 3)*v3

          v(y, z, 4) = X_n1(y, z, 1)*Tv_x + X_n1(y, z, 2)*Tv_y + X_n1(y, z, 3)*Tv_z
          v(y, z, 5) = X_m1(y, z, 1)*Tv_x + X_m1(y, z, 2)*Tv_y + X_m1(y, z, 3)*Tv_z
          v(y, z, 6) = X_l1(y, z, 1)*Tv_x + X_l1(y, z, 2)*Tv_y + X_l1(y, z, 3)*Tv_z

       end do
    end do

  end subroutine Rotate_in_Local_Cordinates


  subroutine Rotate_Back_Physical_Cordinates(u, v, X_l1, X_m1, X_n1, my1, mz1, py1, pz1)

    implicit none

    ! Rotate field variables in local basis vectors
    integer, intent(in) :: my1, mz1, py1, pz1 ! Grid length
    real(kind = wp), dimension(:,:,:), allocatable, intent(in) :: X_l1, X_m1, X_n1                               ! Local basis vectors
    real(kind = wp), dimension(:,:,:), allocatable, intent(inout) :: u, v                                        ! Field variables

    real(kind = wp) :: u1, u2, u3, v1, v2, v3, Tu_x, Tu_y, Tu_z, Tv_x, Tv_y, Tv_z               ! Local work variables
    integer :: y, z                                                                        ! Local grid indices

    do z = mz1, pz1
       do y = my1, py1

          ! Set local grid functions
          u1 = u(y, z, 1)
          u2 = u(y, z, 2)
          u3 = u(y, z, 3)

          v1 = v(y, z, 1)
          v2 = v(y, z, 2)
          v3 = v(y, z, 3)

          Tu_x = u(y, z, 4)
          Tu_y = u(y, z, 5)
          Tu_z = u(y, z, 6)

          Tv_x = v(y, z, 4)
          Tv_y = v(y, z, 5)
          Tv_z = v(y, z, 6)

          ! Rotate back to Cartensian basis vectors

          u(y, z, 1) = X_n1(y, z, 1)*u1 + X_m1(y, z, 1)*u2 + X_l1(y, z, 1)*u3
          u(y, z, 2) = X_n1(y, z, 2)*u1 + X_m1(y, z, 2)*u2 + X_l1(y, z, 2)*u3
          u(y, z, 3) = X_n1(y, z, 3)*u1 + X_m1(y, z, 3)*u2 + X_l1(y, z, 3)*u3

          u(y, z, 4) = X_n1(y, z, 1)*Tu_x + X_m1(y, z, 1)*Tu_y + X_l1(y, z, 1)*Tu_z
          u(y, z, 5) = X_n1(y, z, 2)*Tu_x + X_m1(y, z, 2)*Tu_y + X_l1(y, z, 2)*Tu_z
          u(y, z, 6) = X_n1(y, z, 3)*Tu_x + X_m1(y, z, 3)*Tu_y + X_l1(y, z, 3)*Tu_z

          v(y, z, 1) = X_n1(y, z, 1)*v1 + X_m1(y, z, 1)*v2 + X_l1(y, z, 1)*v3
          v(y, z, 2) = X_n1(y, z, 2)*v1 + X_m1(y, z, 2)*v2 + X_l1(y, z, 2)*v3
          v(y, z, 3) = X_n1(y, z, 3)*v1 + X_m1(y, z, 3)*v2 + X_l1(y, z, 3)*v3

          v(y, z, 4) = X_n1(y, z, 1)*Tv_x + X_m1(y, z, 1)*Tv_y + X_l1(y, z, 1)*Tv_z
          v(y, z, 5) = X_n1(y, z, 2)*Tv_x + X_m1(y, z, 2)*Tv_y + X_l1(y, z, 2)*Tv_z
          v(y, z, 6) = X_n1(y, z, 3)*Tv_x + X_m1(y, z, 3)*Tv_y + X_l1(y, z, 3)*Tv_z

       end do
    end do

  end subroutine Rotate_Back_Physical_Cordinates


  subroutine rupture_front(my1, py1, mz1, pz1, Svel,trup,t)

    implicit none

    integer :: my1, py1, mz1, pz1
    real(kind = wp) :: threshold
    real(kind = wp), dimension(:,:,:), allocatable, intent(in) :: Svel
    real(kind = wp), dimension(:,:,:), allocatable, intent(inout) :: trup
    real(kind = wp), intent(in) :: t
    integer :: i,j

    threshold = 1d-3

    do i = my1, py1
      do j = mz1, pz1
        if (sqrt(Svel(i,j,2)**2+Svel(i,j,3)**2) > threshold) then
          trup(i,j,1) = min(trup(i,j,1),t)
        end if
      end do
    end do

  end subroutine rupture_front

end module CouplingForcing
