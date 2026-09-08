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

    integer :: mx, my, mz, px, py, pz, n

    real(kind = wp), dimension(:,:), allocatable, save :: w1_p, w2_p, w3_p, w1_m, w2_m, w3_m

    real(kind = wp) :: c_p, c_s, rho
    integer :: y, z

    real(kind = wp) :: q_x, q_y, q_z, norm, fu(9)

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%F, 4)

    if (.not. allocated(w1_m)) allocate(w1_m(my:py, mz:pz), w2_m(my:py, mz:pz), &
             w3_m(my:py, mz:pz), w1_p(my:py, mz:pz), &
            w2_p(my:py, mz:pz), w3_p(my:py, mz:pz))
    if (.not. allocated(u_bc)) allocate(u_bc(my:py,mz:pz,n))

    do z = mz, pz
       do y = my, py

          rho = B%M%M(mx, y, z, 3)
          c_p = sqrt((2.0_wp*B%M%M(mx, y, z, 2) &
               & + B%M%M(mx, y, z, 1))/B%M%M(mx, y, z, 3))                                     ! P-wave speed
          c_s = sqrt(B%M%M(mx, y, z, 2)/B%M%M(mx, y, z, 3))                                    ! S-wave speed

          q_x = B%G%metricx(mx, y, z, 1)
          q_y = B%G%metricy(mx, y, z, 1)
          q_z = B%G%metricz(mx, y, z, 1)

          norm = sqrt(q_x**2 + q_y**2 + q_z**2)

          w1_p(y, z) = 1.0_wp/rho*(norm*rho*c_p*B%F%F(mx, y, z, 1) &
               + (q_x*B%F%F(mx, y, z, 4) + q_y*B%F%F(mx, y, z, 7) + q_z*B%F%F(mx, y, z, 8)))

          w1_m(y, z) = 1.0_wp/rho*(norm*rho*c_p*B%F%F(mx, y, z, 1) &
               - (q_x*B%F%F(mx, y, z, 4) + q_y*B%F%F(mx, y, z, 7) + q_z*B%F%F(mx, y, z, 8)))

          w2_p(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(mx, y, z, 2) &
               + (q_x*B%F%F(mx, y, z, 7) + q_y*B%F%F(mx, y, z, 5) + q_z*B%F%F(mx, y, z, 9)))

          w2_m(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(mx, y, z, 2) &
               - (q_x*B%F%F(mx, y, z, 7) + q_y*B%F%F(mx, y, z, 5) + q_z*B%F%F(mx, y, z, 9)))

          w3_p(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(mx, y, z, 3) &
               + (q_x*B%F%F(mx, y, z, 8) + q_y*B%F%F(mx, y, z, 9) + q_z*B%F%F(mx, y, z, 6)))

          w3_m(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(mx, y, z, 3) &
               - (q_x*B%F%F(mx, y, z, 8) + q_y*B%F%F(mx, y, z, 9) + q_z*B%F%F(mx, y, z, 6)))

          if (mms_vars%use_mms) then

             fu = evaluate_mms(t, B%G%X(mx, y, z, :), mms_vars)

             w1_p(y, z) = w1_p(y, z) + 1.0_wp/rho*( &
                  - (norm*rho*c_p*fu(1) &
                  + (q_x*fu(4) + q_y*fu(7) + q_z*fu(8))))

             w1_m(y, z) = w1_m(y, z) + 1.0_wp/rho*( &
                  - (norm*rho*c_p*fu(1) &
                  - (q_x*fu(4) + q_y*fu(7) + q_z*fu(8))))

             w2_p(y, z) = w2_p(y, z) + 1.0_wp/rho*( &
                  - (norm*rho*c_s*fu(2) &
                  + (q_x*fu(7) + q_y*fu(5) + q_z*fu(9))))

             w2_m(y, z) = w2_m(y, z) + 1.0_wp/rho*( &
                  - (norm*rho*c_s*fu(2) &
                  - (q_x*fu(7) + q_y*fu(5) + q_z*fu(9))))

             w3_p(y, z) = w3_p(y, z) + 1.0_wp/rho*( &
                  - (norm*rho*c_s*fu(3) &
                  + (q_x*fu(8) + q_y*fu(9) + q_z*fu(6))))

             w3_m(y, z) = w3_m(y, z) + 1.0_wp/rho*( &
                  - (norm*rho*c_s*fu(3) &
                  - (q_x*fu(8) + q_y*fu(9) + q_z*fu(6))))

          end if

       end do
    end do

    select case(type_of_bc)

    case(1) ! absorbing boundary condition

       u_bc(:, :, 1) = w1_m
       u_bc(:, :, 2) = w2_m
       u_bc(:, :, 3) = w3_m
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    case(2) ! free-surface boundary condition

       u_bc(:, :, 1) = 0.5_wp*(w1_m - w1_p)
       u_bc(:, :, 2) = 0.5_wp*(w2_m - w2_p)
       u_bc(:, :, 3) = 0.5_wp*(w3_m - w3_p)
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end select

  end function BC_Lx


  function BC_Rx(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable  :: u_bc

    integer :: mx, my, mz, px, py, pz, n

    real(kind = wp), dimension(:,:), allocatable, save :: w1_p, w2_p, w3_p, w1_m, w2_m, w3_m

    real(kind = wp) :: c_p, c_s, rho
    integer :: y, z

    real(kind = wp) :: norm, q_x, q_y, q_z, fu(9)

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%F, 4)

    if (.not. allocated(w1_m)) allocate(w1_m(my:py, mz:pz), w2_m(my:py, mz:pz), &
             w3_m(my:py, mz:pz), w1_p(my:py, mz:pz), &
            w2_p(my:py, mz:pz), w3_p(my:py, mz:pz))
    if (.not.allocated(u_bc)) allocate(u_bc(my:py,mz:pz,n))

    do z = mz, pz
       do y = my, py

          rho = B%M%M(px, y, z, 3) 

          q_x = B%G%metricx(px, y, z, 1)
          q_y = B%G%metricy(px, y, z, 1)
          q_z = B%G%metricz(px, y, z, 1)

          c_p = sqrt((2.0_wp*B%M%M(px, y, z, 2) &
               & + B%M%M(px, y, z, 1))/B%M%M(px, y, z, 3))                                     ! P-wave speed
          c_s = sqrt(B%M%M(px, y, z, 2)/B%M%M(px, y, z, 3))                                    ! S-wave speed

          norm = sqrt(q_x**2 + q_y**2 + q_z**2)

          w1_p(y, z) = 1.0_wp/rho*(norm*rho*c_p*B%F%F(px, y, z, 1) &
               & + (q_x*B%F%F(px, y, z, 4) &
               & + q_y*B%F%F(px, y, z, 7) &
               & + q_z*B%F%F(px, y, z, 8)))

          w1_m(y, z) = 1.0_wp/rho*(norm*rho*c_p*B%F%F(px, y, z, 1) &
               & - (q_x*B%F%F(px, y, z, 4) &
               & + q_y*B%F%F(px, y, z, 7) &
               & + q_z*B%F%F(px, y, z, 8)))

          w2_p(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(px, y, z, 2) &
               & + (q_x*B%F%F(px, y, z, 7) &
               & + q_y*B%F%F(px, y, z, 5) &
               & + q_z*B%F%F(px, y, z, 9)))

          w2_m(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(px, y, z, 2) &
               & - (q_x*B%F%F(px, y, z, 7) &
               & + q_y*B%F%F(px, y, z, 5) &
               & + q_z*B%F%F(px, y, z, 9)))

          w3_p(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(px, y, z, 3) &
               & + (q_x*B%F%F(px, y, z, 8) &
               & + q_y*B%F%F(px, y, z, 9) &
               & + q_z*B%F%F(px, y, z, 6)))

          w3_m(y, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(px, y, z, 3) &
               & - (q_x*B%F%F(px, y, z, 8) &
               & + q_y*B%F%F(px, y, z, 9) &
               & + q_z*B%F%F(px, y, z, 6)))

          if (mms_vars%use_mms) then

             fu = evaluate_mms(t, B%G%X(px, y, z, :), mms_vars)

             w1_p(y, z) = w1_p(y, z) + 1.0_wp/rho*( &
                  & - (norm*rho*c_p*fu(1) &
                  & + (q_x*fu(4) &
                  & + q_y*fu(7) &
                  & + q_z*fu(8))))

             w1_m(y, z) = w1_m(y, z) + 1.0_wp/rho*( &
                  & - (norm*rho*c_p*fu(1) &
                  & - (q_x*fu(4) &
                  & + q_y*fu(7) &
                  & + q_z*fu(8))))

             w2_p(y, z) = w2_p(y, z) + 1.0_wp/rho*( &
                  & - (norm*rho*c_s*fu(2) &
                  & + (q_x*fu(7) &
                  & + q_y*fu(5) &
                  & + q_z*fu(9))))

             w2_m(y, z) = w2_m(y, z) + 1.0_wp/rho*( &
                  & - (norm*rho*c_s*fu(2) &
                  & - (q_x*fu(7) &
                  & + q_y*fu(5) &
                  & + q_z*fu(9))))

             w3_p(y, z) = w3_p(y, z) + 1.0_wp/rho*( &
                  & - (norm*rho*c_s*fu(3) &
                  & + (q_x*fu(8) &
                  & + q_y*fu(9) &
                  & + q_z*fu(6))))

             w3_m(y, z) = w3_m(y, z) + 1.0_wp/rho*( &
                  & - (norm*rho*c_s*fu(3) &
                  & - (q_x*fu(8) &
                  & + q_y*fu(9) &
                  & + q_z*fu(6))))

          end if

       end do
    end do

    select case(type_of_bc)

    case(1) ! absorbing boundary condition

       u_bc(:, :, 1) = w1_p
       u_bc(:, :, 2) = w2_p
       u_bc(:, :, 3) = w3_p
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    case(2) ! free-surface boundary condition

       u_bc(:, :, 1) = 0.5_wp*(w1_p - w1_m)
       u_bc(:, :, 2) = 0.5_wp*(w2_p - w2_m)
       u_bc(:, :, 3) = 0.5_wp*(w3_p - w3_m)
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end select

  end function BC_Rx


  function BC_Ly(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, n

    real(kind = wp), dimension(:,:), allocatable, save :: w1_p, w2_p, w3_p, w1_m, w2_m, w3_m

    real(kind = wp) :: c_p, c_s, rho
    integer :: x, z

    real(kind = wp) :: norm, r_x, r_y, r_z, fu(9)
    real(kind = wp) :: switch

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%F, 4)

    if (.not. allocated(w1_m)) allocate(w1_m(mx:px, mz:pz), w2_m(mx:px, mz:pz), &
             w3_m(mx:px, mz:pz), w1_p(mx:px, mz:pz), &
            w2_p(mx:px, mz:pz), w3_p(mx:px, mz:pz))
    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,mz:pz,n))

    ! replace this as in above routines
    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if
    ! end replace

    do z = mz, pz
       do x = mx, px

          fu = evaluate_mms(t, B%G%X(x, my, z, :), mms_vars)

          rho = B%M%M(x, my, z, 3)
          c_p = sqrt((2.0_wp*B%M%M(x, my, z, 2)  + B%M%M(x, my, z, 1))/rho)                         ! P-wave speed
          c_s = sqrt(B%M%M(x, my, z, 2)/rho)                                                         ! S-wave speed

          r_x = B%G%metricx(x, my, z, 2)
          r_y = B%G%metricy(x, my, z, 2)
          r_z = B%G%metricz(x, my, z, 2)

          norm = sqrt(r_x**2 + r_y**2 + r_z**2)

          w1_p(x, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, my, z, 1) &
               & + (r_x*B%F%F(x, my, z, 4) &
               & + r_y*B%F%F(x, my, z, 7) &
               & + r_z*B%F%F(x, my, z, 8)) &
               & - switch*(norm*rho*c_s*fu(1) &
               & + (r_x*fu(4) &
               & + r_y*fu(7) &
               & + r_z*fu(8))))

          w1_m(x, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, my, z, 1) &
               & - (r_x*B%F%F(x, my, z, 4) &
               & + r_y*B%F%F(x, my, z, 7) &
               & + r_z*B%F%F(x, my, z, 8)) &
               & - switch*(norm*rho*c_s*fu(1) &
               & - (r_x*fu(4) &
               & + r_y*fu(7) &
               & + r_z*fu(8))))

          w2_p(x, z) = 1.0_wp/rho*(norm*rho*c_p*B%F%F(x, my, z, 2) &
               & + (r_x*B%F%F(x, my, z, 7) &
               & + r_y*B%F%F(x, my, z, 5) &
               & + r_z*B%F%F(x, my, z, 9)) &
               & - switch*(norm*rho*c_p*fu(2) &
               & + (r_x*fu(7) &
               & + r_y*fu(5) &
               & + r_z*fu(9))))

          w2_m(x, z) = 1.0_wp/rho*(norm*rho*c_p*B%F%F(x, my, z, 2) &
               & - (r_x*B%F%F(x, my, z, 7) &
               & + r_y*B%F%F(x, my, z, 5) &
               & + r_z*B%F%F(x, my, z, 9)) &
               & - switch*(norm*rho*c_p*fu(2) &
               & - (r_x*fu(7) &
               & + r_y*fu(5) &
               & + r_z*fu(9))))

          w3_p(x, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, my, z, 3) &
               & + (r_x*B%F%F(x, my, z, 8) &
               & + r_y*B%F%F(x, my, z, 9) &
               & + r_z*B%F%F(x, my, z, 6))&
               & - switch*(norm*rho*c_s*fu(3) &
               & + (r_x*fu(8) &
               & + r_y*fu(9) &
               & + r_z*fu(6))))

          w3_m(x, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, my, z, 3) &
               & - (r_x*B%F%F(x, my, z, 8) &
               & + r_y*B%F%F(x, my, z, 9) &
               & + r_z*B%F%F(x, my, z, 6))&
               & - switch*(norm*rho*c_s*fu(3) &
               & - (r_x*fu(8) &
               & + r_y*fu(9) &
               & + r_z*fu(6))))

       end do
    end do

    ! absorbing boundary conditions

    if (type_of_bc .eq. 1) then

       u_bc(:, :, 1) = w1_m
       u_bc(:, :, 2) = w2_m
       u_bc(:, :, 3) = w3_m
       u_bc(:, :, 4) = 0.0_wp*w3_m
       u_bc(:, :, 5) = 0.0_wp*w3_m
       u_bc(:, :, 6) = 0.0_wp*w3_m
       u_bc(:, :, 7) = 0.0_wp*w3_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m
    end if

    ! free-surface boundary conditions

    if (type_of_bc .eq. 2) then

       u_bc(:, :, 1) = 0.5_wp*(w1_m - w1_p)
       u_bc(:, :, 2) = 0.5_wp*(w2_m - w2_p)
       u_bc(:, :, 3) = 0.5_wp*(w3_m - w3_p)
       u_bc(:, :, 4) = 0.0_wp*w3_m
       u_bc(:, :, 5) = 0.0_wp*w3_m
       u_bc(:, :, 6) = 0.0_wp*w3_m
       u_bc(:, :, 7) = 0.0_wp*w3_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end if

  end function BC_Ly


  function BC_Ry(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, n

    real(kind = wp), dimension(:,:), allocatable, save :: w1_p, w2_p, w3_p, w1_m, w2_m, w3_m

    real(kind = wp) :: c_p, c_s, rho
    integer :: x, z

    real(kind = wp) :: norm, r_x, r_y, r_z, fu(9)
    real(kind = wp) :: switch

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%F, 4)

    if (.not. allocated(w1_m)) allocate(w1_m(mx:px, mz:pz), w2_m(mx:px, mz:pz), &
             w3_m(mx:px, mz:pz), w1_p(mx:px, mz:pz), &
            w2_p(mx:px, mz:pz), w3_p(mx:px, mz:pz))
    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,mz:pz,n))

    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if
    
    do z = mz, pz
       do x = mx, px

          fu = evaluate_mms(t, B%G%X(x, py, z,:), mms_vars)

          rho = B%M%M(x, py, z, 3)

          c_p = sqrt((2.0_wp*B%M%M(x, py, z, 2) &
               & + B%M%M(x, py, z, 1))/B%M%M(x, py, z, 3))                                     ! P-wave speed
          c_s = sqrt(B%M%M(x, py, z, 2)/B%M%M(x, py, z, 3))                        ! S-wave speed

          r_x = B%G%metricx(x, py, z, 2)
          r_y = B%G%metricy(x, py, z, 2)
          r_z = B%G%metricz(x, py, z, 2)

          norm = sqrt(r_x**2 + r_y**2 + r_z**2)


          w1_p(x, z) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, py, z, 1) &
               & + (r_x*B%F%F(x, py, z, 4) &
               & + r_y*B%F%F(x, py, z, 7) &
               & + r_z*B%F%F(x, py, z, 8)) &
               & - switch*(norm*rho*c_s*fu(1) &
               & + (r_x*fu(4) &
               & + r_y*fu(7) &
               & + r_z*fu(8))))

          w1_m(x, z) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, py, z, 1) &
               & - (r_x*B%F%F(x, py, z, 4) &
               & + r_y*B%F%F(x, py, z, 7) &
               & + r_z*B%F%F(x, py, z, 8)) &
               & - switch*(norm*rho*c_s*fu(1) &
               & - (r_x*fu(4) &
               & + r_y*fu(7) &
               & + r_z*fu(8))))

          w2_p(x, z) =  1.0_wp/rho*(norm*rho*c_p*B%F%F(x, py, z, 2) &
               & + (r_x*B%F%F(x, py, z, 7) &
               & + r_y*B%F%F(x, py, z, 5) &
               & + r_z*B%F%F(x, py, z, 9)) &
               & - switch*(norm*rho*c_p*fu(2) &
               & + (r_x*fu(7) &
               & + r_y*fu(5) &
               & + r_z*fu(9))))

          w2_m(x, z) =  1.0_wp/rho*(norm*rho*c_p*B%F%F(x, py, z, 2) &
               & - (r_x*B%F%F(x, py, z, 7) &
               & + r_y*B%F%F(x, py, z, 5) &
               & + r_z*B%F%F(x, py, z, 9)) &
               & - switch*(norm*rho*c_p*fu(2) &
               & - (r_x*fu(7) &
               & + r_y*fu(5) &
               & + r_z*fu(9))))

          w3_p(x, z) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, py, z, 3) &
               & + (r_x*B%F%F(x, py, z, 8) &
               & + r_y*B%F%F(x, py, z, 9) &
               & + r_z*B%F%F(x, py, z, 6)) &
               & - switch*(norm*rho*c_s*fu(3) &
               & + (r_x*fu(8) &
               & + r_y*fu(9) &
               & + r_z*fu(6))))

          w3_m(x, z) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, py, z, 3) &
               & - (r_x*B%F%F(x, py, z, 8) &
               & + r_y*B%F%F(x, py, z, 9) &
               & + r_z*B%F%F(x, py, z, 6)) &
               & - switch*(norm*rho*c_s*fu(3) &
               & - (r_x*fu(8) &
               & + r_y*fu(9) &
               & + r_z*fu(6))))

       end do
    end do

    ! absorbing boundary condition
    if (type_of_bc .eq. 1) then

       u_bc(:, :, 1) = w1_p
       u_bc(:, :, 2) = w2_p
       u_bc(:, :, 3) = w3_p
       u_bc(:, :, 4) = 0.0_wp*w3_m
       u_bc(:, :, 5) = 0.0_wp*w3_m
       u_bc(:, :, 6) = 0.0_wp*w3_m
       u_bc(:, :, 7) = 0.0_wp*w3_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m
    end if

    ! free-surface boundary condition
    if (type_of_bc .eq. 2) then

       u_bc(:, :, 1) = 0.5_wp*(w1_p - w1_m)
       u_bc(:, :, 2) = 0.5_wp*(w2_p - w2_m)
       u_bc(:, :, 3) = 0.5_wp*(w3_p - w3_m)
       u_bc(:, :, 4) = 0.0_wp*w3_m
       u_bc(:, :, 5) = 0.0_wp*w3_m
       u_bc(:, :, 6) = 0.0_wp*w3_m
       u_bc(:, :, 7) = 0.0_wp*w3_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end if

  end function BC_Ry


  function BC_Lz(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, n

    real(kind = wp), dimension(:,:), allocatable, save :: w1_p, w2_p, w3_p, w1_m, w2_m, w3_m

    real(kind = wp) :: c_p, c_s, rho
    integer :: x, y

    real(kind = wp) :: norm, s_x, s_y, s_z, fu(9)
    real(kind = wp) :: switch

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%F, 4)

    if (.not. allocated(w1_m)) allocate(w1_m(mx:px, my:py), w2_m(mx:px, my:py), &
             w3_m(mx:px, my:py), w1_p(mx:px, my:py), &
            w2_p(mx:px, my:py), w3_p(mx:px, my:py))
    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,my:py,n))

    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if

    do y = my, py
       do x = mx, px

          fu = evaluate_mms(t, B%G%X(x, y, mz,:), mms_vars)

          rho = B%M%M(x, y, mz, 3)

          c_p = sqrt((2.0_wp*B%M%M(x, y, mz, 2) + B%M%M(x, y, mz, 1))/rho)                       ! P-wave speed
          c_s = sqrt(B%M%M(x, y, mz, 2)/rho)                                                      ! S-wave speed

          s_x = B%G%metricx(x, y, mz, 3)
          s_y = B%G%metricy(x, y, mz, 3)
          s_z = B%G%metricz(x, y, mz, 3)

          norm = sqrt(s_x**2 + s_y**2 + s_z**2)

          w1_p(x, y) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, mz, 1) &
               & + (s_x*B%F%F(x, y, mz, 4) &
               & + s_y*B%F%F(x, y, mz, 7) &
               & + s_z*B%F%F(x, y, mz, 8)) &
               & - switch*(norm*c_s*rho*fu(1) &
               & + (s_x*fu(4) &
               & + s_y*fu(7) &
               & + s_z*fu(8))))

          w1_m(x, y) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, mz, 1) &
               & - (s_x*B%F%F(x, y, mz, 4) &
               & + s_y*B%F%F(x, y, mz, 7) &
               & + s_z*B%F%F(x, y, mz, 8)) &
               & - switch*(norm*c_s*rho*fu(1) &
               & - (s_x*fu(4) &
               & + s_y*fu(7) &
               & + s_z*fu(8))))

          w2_p(x, y) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, mz, 2) &
               & + (s_x*B%F%F(x, y, mz, 7) &
               & + s_y*B%F%F(x, y, mz, 5) &
               & + s_z*B%F%F(x, y, mz, 9)) &
               & - switch*(norm*rho*c_s*fu(2) &
               & + (s_x*fu(7) &
               & + s_y*fu(5) &
               & + s_z*fu(9))))

          w2_m(x, y) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, mz, 2) &
               & - (s_x*B%F%F(x, y, mz, 7) &
               & + s_y*B%F%F(x, y, mz, 5) &
               & + s_z*B%F%F(x, y, mz, 9)) &
               & - switch*(norm*rho*c_s*fu(2) &
               & - (s_x*fu(7) &
               & + s_y*fu(5) &
               & + s_z*fu(9))))

          w3_p(x, y) =  1.0_wp/rho*(norm*rho*c_p*B%F%F(x, y, mz, 3) &
               & + (s_x*B%F%F(x, y, mz, 8) &
               & + s_y*B%F%F(x, y, mz, 9) &
               & + s_z*B%F%F(x, y, mz, 6)) &
               & - switch*(norm*rho*c_p*fu(3) &
               & + (s_x*fu(8) &
               & + s_y*fu(9) &
               & + s_z*fu(6))))

          w3_m(x, y) =  1.0_wp/rho*(norm*rho*c_p*B%F%F(x, y, mz, 3) &
               & - (s_x*B%F%F(x, y, mz, 8) &
               & + s_y*B%F%F(x, y, mz, 9) &
               & + s_z*B%F%F(x, y, mz, 6)) &
               & - switch*(norm*rho*c_p*fu(3) &
               & - (s_x*fu(8) &
               & + s_y*fu(9) &
               & + s_z*fu(6))))


       end do
    end do



    ! absorbing boundary conditions
    if (type_of_bc .eq. 1) then

       u_bc(:, :, 1) = w1_m
       u_bc(:, :, 2) = w2_m
       u_bc(:, :, 3) = w3_m
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end if
    ! free-surface boundary conditions
    if (type_of_bc .eq. 2) then

       u_bc(:, :, 1) = 0.5_wp*(w1_m - w1_p)
       u_bc(:, :, 2) = 0.5_wp*(w2_m - w2_p)
       u_bc(:, :, 3) = 0.5_wp*(w3_m - w3_p)
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end if

  end function BC_Lz


  function BC_Rz(B, type_of_bc, mms_vars, t) result(u_bc)

    use datatypes, only : block_type, mms_type

    implicit none

    type(block_type), intent(in) :: B
    integer, intent(in) :: type_of_bc
    type(mms_type),intent(in) :: mms_vars 
    real(kind = wp), intent(in) :: t
    real(kind = wp), dimension(:,:,:), allocatable :: u_bc

    integer :: mx, my, mz, px, py, pz, n

    real(kind = wp), dimension(:,:), allocatable, save :: w1_p, w2_p, w3_p, w1_m, w2_m, w3_m

    real(kind = wp) :: c_p, c_s, rho
    integer :: x, y

    real(kind = wp) :: norm, s_x, s_y, s_z, fu(9)
    real(kind = wp) :: switch

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%F, 4)

    if (.not. allocated(w1_m)) allocate(w1_m(mx:px, my:py), w2_m(mx:px, my:py), &
             w3_m(mx:px, my:py), w1_p(mx:px, my:py), &
            w2_p(mx:px, my:py), w3_p(mx:px, my:py))
    if (.not.allocated(u_bc)) allocate(u_bc(mx:px,my:py,n))

    if (mms_vars%use_mms) then
       switch = 1.0_wp
    else
       switch = 0.0_wp
    end if

    do y = my, py
       do x = mx, px

          fu = evaluate_mms(t, B%G%X(x, y, pz, :), mms_vars)

          rho = B%M%M(x, y, pz, 3)
          c_p = sqrt((2.0_wp*B%M%M(x, y, pz, 2) + B%M%M(x, y, pz, 1))/rho)                         ! P-wave speed
          c_s = sqrt(B%M%M(x, y, pz, 2)/rho)                                                           ! S-wave speed

          s_x = B%G%metricx(x, y, pz, 3)
          s_y = B%G%metricy(x, y, pz, 3)
          s_z = B%G%metricz(x, y, pz, 3)

          norm = sqrt(s_x**2 + s_y**2 + s_z**2)   

          w1_p(x, y) = 1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, pz, 1) &
               & + (s_x*B%F%F(x, y, pz, 4) &
               & + s_y*B%F%F(x, y, pz, 7) &
               & + s_z*B%F%F(x, y, pz, 8)) &
               & - switch*(norm*rho*c_s*fu(1) &
               & + (s_x*fu(4) &
               & + s_y*fu(7) &
               & + s_z*fu(8))))

          w1_m(x, y) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, pz, 1) &
               & - (s_x*B%F%F(x, y, pz, 4) &
               & + s_y*B%F%F(x, y, pz, 7) &
               & + s_z*B%F%F(x, y, pz, 8)) &
               & - switch*(norm*rho*c_s*fu(1) &
               & - (s_x*fu(4) &
               & + s_y*fu(7) &
               & + s_z*fu(8))))

          w2_p(x, y) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, pz, 2) &
               & + (s_x*B%F%F(x, y, pz, 7) &
               & + s_y*B%F%F(x, y, pz, 5) &
               & + s_z*B%F%F(x, y, pz, 9)) &
               & - switch*(norm*rho*c_s*fu(2) &
               & + (s_x*fu(7) &
               & + s_y*fu(5) &
               & + s_z*fu(9))))

          w2_m(x, y) =  1.0_wp/rho*(norm*rho*c_s*B%F%F(x, y, pz, 2) &
               & - (s_x*B%F%F(x, y, pz, 7) &
               & + s_y*B%F%F(x, y, pz, 5) &
               & + s_z*B%F%F(x, y, pz, 9)) &
               & - switch*(norm*rho*c_s*fu(2) &
               & - (s_x*fu(7) &
               & + s_y*fu(5) &
               & + s_z*fu(9))))


          w3_p(x, y) =  1.0_wp/rho*(norm*rho*c_p*B%F%F(x, y, pz, 3) &
               & + (s_x*B%F%F(x, y, pz, 8) &
               & + s_y*B%F%F(x, y, pz, 9) &
               & + s_z*B%F%F(x, y, pz, 6)) &
               & - switch*(norm*rho*c_p*fu(3) &
               & + (s_x*fu(8) &
               & + s_y*fu(9) &
               & + s_z*fu(6))))

          w3_m(x, y) =  1.0_wp/rho*(norm*rho*c_p*B%F%F(x, y, pz, 3) &
               & - (s_x*B%F%F(x, y, pz, 8) &
               & + s_y*B%F%F(x, y, pz, 9) &
               & + s_z*B%F%F(x, y, pz, 6)) &
               & - switch*(norm*rho*c_p*fu(3) &
               & - (s_x*fu(8) &
               & + s_y*fu(9) &
               & + s_z*fu(6))))

       end do
    end do

    ! absorbing boundary conditions

    if (type_of_bc .eq. 1) then

       u_bc(:, :, 1) = w1_p
       u_bc(:, :, 2) = w2_p
       u_bc(:, :, 3) = w3_p
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end if

    ! free-surface boundary conditions 
    if (type_of_bc .eq. 2) then

       u_bc(:, :, 1) = 0.5_wp*(w1_p - w1_m)
       u_bc(:, :, 2) = 0.5_wp*(w2_p - w2_m)
       u_bc(:, :, 3) = 0.5_wp*(w3_p - w3_m)
       u_bc(:, :, 4) = 0.0_wp*w1_m
       u_bc(:, :, 5) = 0.0_wp*w1_m
       u_bc(:, :, 6) = 0.0_wp*w1_m
       u_bc(:, :, 7) = 0.0_wp*w2_m
       u_bc(:, :, 8) = 0.0_wp*w3_m
       u_bc(:, :, 9) = 0.0_wp*w3_m

    end if

  end function BC_Rz


end module BoundaryConditions
