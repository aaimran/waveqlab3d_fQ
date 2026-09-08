module BoundaryConditions

  use common, only : wp
  use datatypes, only : boundary_type
  use mms, only : mms_type,evaluate_mms

  implicit none

contains


  ! Construct boundary procedures for the six faces of a unit cube: Lx, Rx, Ly, Ry, Lz, Rz
  ! Here only two types of BCs can be implemented: a free-surface BC or a characteristic BC

  subroutine init_boundaries(B, C, lqrs, rqrs)

    use mpi3dbasic, only : rank
    use mpi3dcomm, only : cartesian3d_t

    implicit none

    integer, intent(in) :: lqrs(3), rqrs(3)
    type(cartesian3d_t), intent(in) :: C
    type(boundary_type), intent(inout) :: B

    ! Set default parameters
    ! Parameters to set boundary conditions
    ! Characteristics: Set type_of_bc = 1
    ! Free_Surface: Set type_of_bc = 2
    ! Material Interface: Set type_of_bc = 0     

    B%Lx = -1
    B%Ly = -1
    B%Lz = -1

    B%Rx = -1
    B%Ry = -1
    B%Rz = -1


    if (C%pq == C%nq) B%Rx = rqrs(1)
    if (C%pr == C%nr) B%Ry = rqrs(2)
    if (C%ps == C%ns) B%Rz = rqrs(3)
    if (C%mq == 1) B%Lx = lqrs(1)
    if (C%mr == 1) B%Ly = lqrs(2)
    if (C%ms == 1) B%Lz = lqrs(3)

  end subroutine init_boundaries

  function BC_Lx(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, nf
    real(kind = wp) :: c_p, c_s, rho, mu,lambda
    integer :: y, z
    real(kind = wp) :: q_x, q_y, q_z, norm, fu(9), ubc(9)
    real(kind = wp) :: switch


    real(kind = wp) :: l(3), m(3), n(3)
    real(kind = wp) :: u(9)
    real(kind = wp) :: Zp, Zs, r, theta, alpha


    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    nf = size(B%F%F, 4)
   
    if (.not. allocated(u_bc)) allocate(u_bc(my:py,mz:pz,nf))

     ! replace this as in above routines
    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if
    ! end replace


    do z = mz, pz
       do y = my, py

          ! evaluate the mms solutions on the boundary
          fu = evaluate_mms(t, B%G%X(mx, y, z, :), mms_vars)

          ! material properties
          rho = B%M%M(mx, y, z, 3)
          mu = B%M%M(mx, y, z, 2)
          lambda =  B%M%M(mx, y, z, 1)

          ! wave speeds
          c_p = sqrt((2.0_wp*mu + lambda)/rho)                  ! P-wave speed
          c_s = sqrt(mu/rho)                                    ! S-wave speed

          ! impedances
          Zp = rho*c_p
          Zs = rho*c_s

          ! normal metric derivatives
          q_x = B%G%metricx(mx, y, z, 1)
          q_y = B%G%metricy(mx, y, z, 1)
          q_z = B%G%metricz(mx, y, z, 1)

          norm = sqrt(q_x**2 + q_y**2 + q_z**2)

          ! local basis vectors
          l(:) = B%B(1)%n_l(y,z,:)
          m(:) = B%B(1)%n_m(y,z,:)
          n(:) = B%B(1)%n_n(y,z,:)

          ! extract boundary fields and modify with mms forcing if needed
          u(:) = B%F%F(mx,y,z,:) - switch*fu

          ! compute the SAT boundary forcing for various boundary consitions
          if (type_of_bc .eq. 1) then

             r = 0.0_wp
             alpha = 2.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(y, z,:) = ubc(:)

          end if

          ! free-surface boundary conditions 
          if (type_of_bc .eq. 2) then

             r = 1.0_wp
             alpha = 1.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(y, z,:) = ubc(:)

          end if

          ! clamped wall boundary conditions 
          if (type_of_bc .eq. 3) then
             
             r = -1.0_wp
             alpha = 0.0_wp
             theta = -1.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

             u_bc(y, z,:) = ubc(:)
          end if

          
       end do
    end do



  end function BC_Lx


  function BC_Rx(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable  :: u_bc

    integer :: mx, my, mz, px, py, pz, nf
    real(kind = wp) :: c_p, c_s, rho, mu,lambda
    integer :: y, z
    real(kind = wp) :: norm, q_x, q_y, q_z, fu(9), ubc(9)
    real(kind = wp) :: switch
    real(kind = wp) :: l(3), m(3), n(3)
    real(kind = wp) :: u(9)
    real(kind = wp) :: Zp, Zs, r, theta, alpha

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    nf = size(B%F%F, 4)

    if (.not.allocated(u_bc)) allocate(u_bc(my:py,mz:pz,nf))

    ! replace this as in above routines
    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if
    ! end replace

    do z = mz, pz
       do y = my, py

          ! mms forcing (used only if mms switch is on)
          fu = evaluate_mms(t, B%G%X(px, y, z, :), mms_vars)

          ! material parameters
          rho = B%M%M(px, y, z, 3)
          mu = B%M%M(px, y, z, 2)
          lambda =  B%M%M(px, y, z, 1)

          ! wave speeds
          c_p = sqrt((2.0_wp*mu + lambda)/rho)                  
          c_s = sqrt(mu/rho)                                   

          ! impedances
          Zp = rho*c_p
          Zs = rho*c_s

          ! normal metric derivatives
          q_x = B%G%metricx(px, y, z, 1)
          q_y = B%G%metricy(px, y, z, 1)
          q_z = B%G%metricz(px, y, z, 1)

          norm = sqrt(q_x**2 + q_y**2 + q_z**2)

          ! orthonormal basis vectors
          l(:) = B%B(2)%n_l(y,z,:)
          m(:) = B%B(2)%n_m(y,z,:)
          n(:) = B%B(2)%n_n(y,z,:)

          ! extract the boundary fields
          !(subtracting the mms solution if mms swich is on)
          u(:) = B%F%F(px,y,z,:) - switch*fu

           ! compute the SAT boundary forcing for various boundary conditions
           if (type_of_bc .eq. 1) then

             r = 0.0_wp
             alpha = 2.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(y, z,:) = ubc(:)

          end if

          ! free-surface boundary conditions 
          if (type_of_bc .eq. 2) then

             r = 1.0_wp
             alpha = 1.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(y, z,:) = ubc(:)

          end if

          ! clamped wall boundary conditions 
          if (type_of_bc .eq. 3) then
             
             r = -1.0_wp
             alpha = 0.0_wp
             theta = 1.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

             u_bc(y, z,:) = ubc(:)
          end if
          
          
       end do
    end do


  end function BC_Rx


  function BC_Ly(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, nf
    real(kind = wp) :: c_p, c_s, rho, mu,lambda
    integer :: x, z
    real(kind = wp) :: norm, r_x, r_y, r_z, fu(9), ubc(9)
    real(kind = wp) :: switch

    real(kind = wp) :: l(3), m(3), n(3)
    real(kind = wp) :: u(9)
    real(kind = wp) :: Zp, Zs, r, theta, alpha

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    nf = size(B%F%F, 4)


    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,mz:pz,nf))

    ! replace this as in above routines
    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if
    ! end replace

    do z = mz, pz
       do x = mx, px

          ! evaluate the mms solutions on the boundary
          fu = evaluate_mms(t, B%G%X(x, my, z, :), mms_vars)

          ! material properties
          rho = B%M%M(x, my, z, 3)
          mu = B%M%M(x, my, z, 2)
          lambda =  B%M%M(x, my, z, 1)
          
          ! wave speeds
          c_p = sqrt((2.0_wp*mu + lambda)/rho)                  
          c_s = sqrt(mu/rho)                                   

          ! impedances
          Zp = rho*c_p
          Zs = rho*c_s

          ! normal metric derivatives
          r_x = B%G%metricx(x, my, z, 2)
          r_y = B%G%metricy(x, my, z, 2)
          r_z = B%G%metricz(x, my, z, 2)

          norm = sqrt(r_x**2 + r_y**2 + r_z**2)

          ! orthonormal basis vectors
          l(:) = B%B(3)%n_l(x,z,:)
          m(:) = B%B(3)%n_m(x,z,:)
          n(:) = B%B(3)%n_n(x,z,:)

          ! extract the boundary fields
          !(subtracting the mms solution if mms swich is on)
          u(:) = B%F%F(x,my,z,:) - switch*fu

          ! compute the SAT boundary forcing for various boundary consitions
          if (type_of_bc .eq. 1) then
             
             r = 0.0_wp
             alpha = 2.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, z,:) = ubc(:)

          end if

          ! free-surface boundary conditions 
          if (type_of_bc .eq. 2) then

             r = 1.0_wp
             alpha = 1.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, z,:) = ubc(:)

          end if

          ! clamped wall boundary conditions 
          if (type_of_bc .eq. 3) then
             
             r = -1.0_wp
             alpha = 0.0_wp
             theta = -1.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

             u_bc(x, z,:) = ubc(:)
          end if

       end do
    end do
    
    
  end function BC_Ly


  function BC_Ry(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, nf
    real(kind = wp) :: c_p, c_s, rho, mu,lambda
    integer :: x, z
    real(kind = wp) :: norm, r_x, r_y, r_z, fu(9), ubc(9)
    real(kind = wp) :: switch

    real(kind = wp) :: l(3), m(3), n(3)
    real(kind = wp) :: u(9)
    real(kind = wp) :: Zp, Zs, r, theta, alpha

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    nf = size(B%F%F, 4)

    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,mz:pz,nf))

    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if

    do z = mz, pz
       do x = mx, px

          ! evaluate the mms solutions on the boundary
          fu = evaluate_mms(t, B%G%X(x, py, z,:), mms_vars)

          ! material properties
          rho = B%M%M(x, py, z, 3)
          mu = B%M%M(x, py, z, 2)
          lambda =  B%M%M(x, py, z, 1)
                  
          ! wave speeds
          c_p = sqrt((2.0_wp*mu + lambda)/rho)                  
          c_s = sqrt(mu/rho)                                   

          ! impedances
          Zp = rho*c_p
          Zs = rho*c_s

          ! normal metric derivatives
          r_x = B%G%metricx(x, py, z, 2)
          r_y = B%G%metricy(x, py, z, 2)
          r_z = B%G%metricz(x, py, z, 2)

          norm = sqrt(r_x**2 + r_y**2 + r_z**2)

          ! orthonormal basis vectors
          l(:) = B%B(4)%n_l(x,z,:)
          m(:) = B%B(4)%n_m(x,z,:)
          n(:) = B%B(4)%n_n(x,z,:)
          
          ! extract the boundary fields
          !(subtracting the mms solution if mms swich is on)
          u(:) = B%F%F(x,py,z,:) - switch*fu

          ! compute the SAT boundary forcing for various boundary consitions
          if (type_of_bc .eq. 1) then

             r = 0.0_wp
             alpha = 2.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, z,:) = ubc(:)

          end if

          ! free-surface boundary conditions 
          if (type_of_bc .eq. 2) then

             r = 1.0_wp
             alpha = 1.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, z,:) = ubc(:)

          end if

          ! clamped wall boundary conditions 
          if (type_of_bc .eq. 3) then
             
             r = -1.0_wp
             alpha = 0.0_wp
             theta = 1.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

             u_bc(x, z,:) = ubc(:)
          end if
          
       end do
    end do

  end function BC_Ry


  function BC_Lz(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, nf
    real(kind = wp) :: c_p, c_s, rho, mu, lambda
    integer :: x, y
    real(kind = wp) :: norm, s_x, s_y, s_z, fu(9), ubc(9)
    real(kind = wp) :: switch

    real(kind = wp) :: l(3), m(3), n(3)
    real(kind = wp) :: u(9)
    real(kind = wp) :: Zp, Zs, r, theta, alpha

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    nf = size(B%F%F, 4)

    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,my:py,nf))

    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if

    do y = my, py
       do x = mx, px

          ! evaluate the mms solutions on the boundary
          fu = evaluate_mms(t, B%G%X(x, y, mz,:), mms_vars)
          
          ! material properties
          lambda = B%M%M(x, y, mz, 1)
          mu = B%M%M(x, y, mz, 2)
          rho = B%M%M(x, y, mz, 3)
          
          ! wave speeds
          c_p = sqrt((2.0_wp*mu + lambda)/rho)                  
          c_s = sqrt(mu/rho)                                   

          ! impedances
          Zp = rho*c_p
          Zs = rho*c_s

          ! normal metric derivatives
          s_x = B%G%metricx(x, y, mz, 3)
          s_y = B%G%metricy(x, y, mz, 3)
          s_z = B%G%metricz(x, y, mz, 3)

          norm = sqrt(s_x**2 + s_y**2 + s_z**2)

          ! orthonormal basis vectors
          l(:) = B%B(5)%n_l(x,y,:)
          m(:) = B%B(5)%n_m(x,y,:)
          n(:) = B%B(5)%n_n(x,y,:)

          ! extract the boundary fields
          !(subtracting the mms solution if mms swich is on)
          u(:) = B%F%F(x,y,mz,:) - switch*fu

          ! compute the SAT boundary forcing for various boundary consitions
          if (type_of_bc .eq. 1) then

             r = 0.0_wp
             alpha = 2.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, y,:) = ubc(:)

          end if

          ! free-surface boundary conditions 
          if (type_of_bc .eq. 2) then

             r = 1.0_wp
             alpha = 1.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, y,:) = ubc(:)

          end if

          ! clamped wall boundary conditions 
          if (type_of_bc .eq. 3) then
             
             r = -1.0_wp
             alpha = 0.0_wp
             theta = -1.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

             u_bc(x, y,:) = ubc(:)
          end if

       end do
    end do


  end function BC_Lz


  function BC_Rz(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, nf

    real(kind = wp) :: c_p, c_s, rho, mu, lambda
    integer :: x, y

    real(kind = wp) :: norm, s_x, s_y, s_z, fu(9), ubc(9)
    real(kind = wp) :: switch

    real(kind = wp) :: l(3), m(3), n(3)
    real(kind = wp) :: u(9)
    real(kind = wp) :: Zp, Zs, r, theta, alpha

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    nf = size(B%F%F, 4)


    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,my:py,nf))

    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if

    do y = my, py
       do x = mx, px

          ! evaluate the mms solutions on the boundary
          fu = evaluate_mms(t, B%G%X(x, y, pz, :), mms_vars)

          ! material properties
          lambda = B%M%M(x, y, pz, 1)
          mu = B%M%M(x, y, pz, 2)
          rho = B%M%M(x, y, pz, 3)
          
          ! wave speeds
          c_p = sqrt((2.0_wp*mu + lambda)/rho)                  
          c_s = sqrt(mu/rho)                                   

          ! impedances
          Zp = rho*c_p
          Zs = rho*c_s

          ! normal metric derivatives
          s_x = B%G%metricx(x, y, pz, 3)
          s_y = B%G%metricy(x, y, pz, 3)
          s_z = B%G%metricz(x, y, pz, 3)

          norm = sqrt(s_x**2 + s_y**2 + s_z**2) 

          ! orthonormal basis vectors
          l(:) = B%B(6)%n_l(x,y,:)
          m(:) = B%B(6)%n_m(x,y,:)
          n(:) = B%B(6)%n_n(x,y,:)

          ! extract the boundary fields
          !(subtracting the mms solution if mms swich is on)
          u(:) = B%F%F(x,y,pz,:) - switch*fu

          ! compute the SAT boundary forcing for various boundary consitions
          if (type_of_bc .eq. 1) then

             r = 0.0_wp
             alpha = 2.0_wp
             theta = 0.0_wp
             
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, y,:) = ubc(:)

          end if

          ! free-surface boundary conditions 
          if (type_of_bc .eq. 2) then

             r = 1.0_wp
             alpha = 1.0_wp
             theta = 0.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)
             
             u_bc(x, y,:) = ubc(:)

          end if

          ! clamped wall boundary conditions 
          if (type_of_bc .eq. 3) then
             
             r = -1.0_wp
             alpha = 0.0_wp
             theta = 1.0_wp
             ! set characterististics bc with the reflection coeffients r
             call BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

             u_bc(x, y,:) = ubc(:)
          end if

       end do
    end do


  end function BC_Rz

  subroutine BCp(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

    real(kind = wp), intent(in) :: u(9), l(3), m(3), n(3)
    real(kind = wp), intent(in) :: lambda, mu, rho, norm, Zp, Zs, r, alpha, theta
    real(kind = wp), intent(inout) :: ubc(9)
    
    real(kind = wp) :: v(3), T(3), lambda2mu
    real(kind = wp) :: Tl, Tm, Tn, vl, vm, vn
    real(kind = wp) :: pl, pm, pn, ql, qm, qn
    real(kind = wp) :: Fx, Fy, Fz, F_x, F_y, F_z
    
    ! extract particle velocities
    v(:) = u(1:3)                              

    ! extract tractions
    T(1) =  n(1)*u(4) + n(2)*u(7) + n(3)*u(8)
    T(2) =  n(1)*u(7) + n(2)*u(5) + n(3)*u(9)
    T(3) =  n(1)*u(8) + n(2)*u(9) + n(3)*u(6)

    ! rotate into local orthogonal coordinates: l, m, n
    vl =  l(1)*v(1) + l(2)*v(2) + l(3)*v(3)
    vm =  m(1)*v(1) + m(2)*v(2) + m(3)*v(3)
    vn =  n(1)*v(1) + n(2)*v(2) + n(3)*v(3)

    Tl =  l(1)*T(1) + l(2)*T(2) + l(3)*T(3)
    Tm =  m(1)*T(1) + m(2)*T(2) + m(3)*T(3)
    Tn =  n(1)*T(1) + n(2)*T(2) + n(3)*T(3)

    ! set characteristics bc with the refelction coefficient r
    pl = 0.5_wp*((1.0_wp-r)*Zs*vl + (1.0_wp + r)*Tl)
    pm = 0.5_wp*((1.0_wp-r)*Zs*vm + (1.0_wp + r)*Tm)
    pn = 0.5_wp*((1.0_wp-r)*Zp*vn + (1.0_wp + r)*Tn)

    ql = 0.5_wp*((1.0_wp-r)*vl + (1.0_wp + r)/Zs*Tl)
    qm = 0.5_wp*((1.0_wp-r)*vm + (1.0_wp + r)/Zs*Tm)
    qn = 0.5_wp*((1.0_wp-r)*vn + (1.0_wp + r)/Zp*Tn)


    ! rotate the bc back to x,y,z coordinates
    Fx = l(1)*pl + m(1)*pm + n(1)*pn
    Fy = l(2)*pl + m(2)*pm + n(2)*pn
    Fz = l(3)*pl + m(3)*pm + n(3)*pn

    F_x = l(1)*ql + m(1)*qm + n(1)*qn
    F_y = l(2)*ql + m(2)*qm + n(2)*qn
    F_z = l(3)*ql + m(3)*qm + n(3)*qn

    lambda2mu = lambda + 2.0_wp*mu
    
    ! set the boundary forcing to be penalized on the boundary
    ubc(1) = alpha*norm/rho*Fx
    ubc(2) = alpha*norm/rho*Fy
    ubc(3) = alpha*norm/rho*Fz
    ubc(4) = theta*norm*(lambda2mu*n(1)*F_x + lambda*(n(2)*F_y + n(3)*F_z))
    ubc(5) = theta*norm*(lambda2mu*n(2)*F_y + lambda*(n(1)*F_x + n(3)*F_z))
    ubc(6) = theta*norm*(lambda2mu*n(3)*F_z + lambda*(n(1)*F_x + n(2)*F_y))
    ubc(7) = theta*mu*norm*(n(2)*F_x + n(1)*F_y)
    ubc(8) = theta*mu*norm*(n(3)*F_x + n(1)*F_z)
    ubc(9) = theta*mu*norm*(n(3)*F_y + n(2)*F_z)

  end subroutine BCp


   subroutine BCm(ubc, u, rho, mu, lambda, l, m, n, norm, Zs, Zp, r, alpha, theta)

    real(kind = wp), intent(in) :: u(9), l(3), m(3), n(3)
    real(kind = wp), intent(in) :: lambda, mu, rho, norm, Zp, Zs, r, alpha, theta
    real(kind = wp), intent(inout) :: ubc(9)
    
    real(kind = wp) :: v(3), T(3), lambda2mu
    real(kind = wp) :: Tl, Tm, Tn, vl, vm, vn
    real(kind = wp) :: pl, pm, pn, ql, qm, qn
    real(kind = wp) :: Fx, Fy, Fz, F_x, F_y, F_z
    
    ! extract particle velocities
    v(:) = u(1:3)                              

    ! extract tractions
    T(1) =  n(1)*u(4) + n(2)*u(7) + n(3)*u(8)
    T(2) =  n(1)*u(7) + n(2)*u(5) + n(3)*u(9)
    T(3) =  n(1)*u(8) + n(2)*u(9) + n(3)*u(6)

    ! rotate into local orthogonal coordinates: l, m, n
    vl =  l(1)*v(1) + l(2)*v(2) + l(3)*v(3)
    vm =  m(1)*v(1) + m(2)*v(2) + m(3)*v(3)
    vn =  n(1)*v(1) + n(2)*v(2) + n(3)*v(3)

    Tl =  l(1)*T(1) + l(2)*T(2) + l(3)*T(3)
    Tm =  m(1)*T(1) + m(2)*T(2) + m(3)*T(3)
    Tn =  n(1)*T(1) + n(2)*T(2) + n(3)*T(3)

    ! set characteristics bc with the refelction coefficient r
    pl = 0.5_wp*((1.0_wp-r)*Zs*vl - (1.0_wp + r)*Tl)
    pm = 0.5_wp*((1.0_wp-r)*Zs*vm - (1.0_wp + r)*Tm)
    pn = 0.5_wp*((1.0_wp-r)*Zp*vn - (1.0_wp + r)*Tn)

    ql = 0.5_wp*((1.0_wp-r)*vl - (1.0_wp + r)/Zs*Tl)
    qm = 0.5_wp*((1.0_wp-r)*vm - (1.0_wp + r)/Zs*Tm)
    qn = 0.5_wp*((1.0_wp-r)*vn - (1.0_wp + r)/Zp*Tn)


    ! rotate the bc back to x,y,z coordinates
    Fx = l(1)*pl + m(1)*pm + n(1)*pn
    Fy = l(2)*pl + m(2)*pm + n(2)*pn
    Fz = l(3)*pl + m(3)*pm + n(3)*pn

    F_x = l(1)*ql + m(1)*qm + n(1)*qn
    F_y = l(2)*ql + m(2)*qm + n(2)*qn
    F_z = l(3)*ql + m(3)*qm + n(3)*qn

    lambda2mu = lambda + 2.0_wp*mu
    
    ! set the boundary forcing to be penalized on the boundary
    ubc(1) = alpha*norm/rho*Fx
    ubc(2) = alpha*norm/rho*Fy
    ubc(3) = alpha*norm/rho*Fz
    ubc(4) = theta*norm*(lambda2mu*n(1)*F_x + lambda*(n(2)*F_y + n(3)*F_z))
    ubc(5) = theta*norm*(lambda2mu*n(2)*F_y + lambda*(n(1)*F_x + n(3)*F_z))
    ubc(6) = theta*norm*(lambda2mu*n(3)*F_z + lambda*(n(1)*F_x + n(2)*F_y))
    ubc(7) = theta*mu*norm*(n(2)*F_x + n(1)*F_y)
    ubc(8) = theta*mu*norm*(n(3)*F_x + n(1)*F_z)
    ubc(9) = theta*mu*norm*(n(3)*F_y + n(2)*F_z)
    
  end subroutine BCm

  
end module BoundaryConditions
