module mms

  ! MMS = method of manufactured solutions
  ! pure sinusoidal solutions for velocity and stress fields

  use common, only : wp
  use datatypes, only : mms_type

  implicit none

  real(kind = wp),parameter :: pi = 3.141592653589793_wp


contains


  subroutine init_mms(input,M)

    implicit none

    integer,intent(in) :: input
    type(mms_type),intent(out) :: M

    logical :: use_mms
    real(kind = wp) :: nt, nx, ny, nz
    integer :: stat

    namelist /mms_list/ use_mms, nt, nx, ny, nz

    ! default parameters

    use_mms = .false.                                   ! true for MMS
    nt = 2.0_wp                                   ! temporal frequency                               
    nx = 1.0_wp                               ! spatial wavenumber in the x-direction
    ny = 1.0_wp                               ! spatial wavenumber in the y-direction
    nz = 1.0_wp                               ! spatial wavenumber in the z-direction

    rewind(input)
    read(input,nml=mms_list,iostat=stat)
    if (stat>0) stop 'error reading namelist mms_list'

    M%use_mms = use_mms
    M%nx = nx
    M%ny = ny
    M%nz = nz
    M%nt = nt

  end subroutine init_mms


  subroutine eval_mms(u, t, XX, mms_vars, mx, my, mz, px, py, pz)

    ! evaluate the manufactured solution at time t

    implicit none

    integer, intent(in) :: mx, my, mz, px, py, pz                                                 ! spatial indices
    real(kind = wp), intent(in) :: t                                                          ! time
    integer :: x, y, z                                                             ! spatial indices
    type(mms_type),intent(in) :: mms_vars
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: XX                                     ! grid   
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(inout) :: u                                   ! manufactured fields

    if (.not.mms_vars%use_mms) return ! return if no forcing

    do z = mz, pz
       do y = my, py
          do x = mx, px
             u(x,y,z,:) =  evaluate_mms(t, XX(x,y,z,:), mms_vars)
          end do
       end do
    end do

  end subroutine eval_mms


  subroutine mms_forcing(DU, t, XX, LambdaMu, mx, my, mz, px, py, pz, mms_vars)

    ! evaluate mms source terms in interior

    implicit none

    real(kind = wp), intent(in) :: t                                                 ! time
    integer, intent(in) :: mx, my, mz, px, py, pz                                     ! array dimensions
    type(mms_type),intent(in) :: mms_vars
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: XX                            ! grid   
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: LambdaMu                      ! material lambda and mu
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(inout) :: DU                         ! rates array

    real(kind = wp), dimension(9) :: u_t, u_x, u_y, u_z, f                           ! work arrays for derivatives and source term
    integer :: x,y,z

    if (.not.mms_vars%use_mms) return                                     ! return if no forcing

    do z = mz, pz
       do y = my, py
          do x = mx, px

             ! compute all derivatives at a single point
             u_t = evaluate_mms_dt(t, XX(x, y, z, :), mms_vars)
             u_x = evaluate_mms_dx(t, XX(x, y, z, :), mms_vars)
             u_y = evaluate_mms_dy(t, XX(x, y, z, :), mms_vars)
             u_z = evaluate_mms_dz(t, XX(x, y, z, :), mms_vars)

             ! construct:   F(x, y, z, t) = AU_x + BU_y + CU_z - U_t

             f(1) = (1.0_wp/LambdaMu(x, y, z, 3))*(u_x(4) + u_y(7) + u_z(8))
             f(2) = (1.0_wp/LambdaMu(x, y, z, 3))*(u_x(7) + u_y(5) + u_z(9))
             f(3) = (1.0_wp/LambdaMu(x, y, z, 3))*(u_x(8) + u_y(9) + u_z(6))

             f(4) = (2.0_wp*LambdaMu(x, y, z, 2) + LambdaMu(x, y, z, 1))*u_x(1) &
                  + LambdaMu(x, y, z, 1)*u_y(2)  + LambdaMu(x, y, z, 1)*u_z(3)
             f(5) = LambdaMu(x, y, z, 1)*u_x(1) &
                  + (2.0_wp*LambdaMu(x, y, z, 2) + LambdaMu(x, y, z, 1))*u_y(2) &
                  + LambdaMu(x, y, z, 1)*u_z(3)
             f(6) = LambdaMu(x, y, z, 1)*u_x(1) + LambdaMu(x, y, z, 1)*u_y(2) &
                  + (2.0_wp*LambdaMu(x, y, z, 2) + LambdaMu(x, y, z, 1))*u_z(3)

             f(7) = LambdaMu(x, y, z, 2)*(u_y(1) + u_x(2))
             f(8) = LambdaMu(x, y, z, 2)*(u_z(1) + u_x(3))
             f(9) = LambdaMu(x, y, z, 2)*(u_z(2) + u_y(3))

             f = f-u_t

             ! add source terms to rates

             DU(x,y,z,:) = DU(x,y,z,:) - f

          end do
       end do
    end do

  end subroutine mms_forcing


  pure function evaluate_mms(t, XX, M) result(u)

    implicit none

    real(kind = wp), dimension(9) :: u !< manufactured fields  
    real(kind = wp), intent(in) :: t !< time
    real(kind = wp), dimension(:), intent(in) :: XX !< grid                            
    type(mms_type),intent(in) :: M

    ! velocities
    u(1:3) = cos(M%nt*pi*t)*sin(M%nx*pi*XX(1))*sin(M%ny*pi*XX(2))*sin(M%nz*pi*XX(3))

    ! stresses
    u(4:9) = cos(M%nt*pi*t)*cos(M%nx*pi*XX(1))*cos(M%ny*pi*XX(2))*cos(M%nz*pi*XX(3))

  end function evaluate_mms


  pure function evaluate_mms_dt(t, XX, M) result(u)

    implicit none

    real(kind = wp), dimension(9) :: u ! manufactured fields  
    real(kind = wp), intent(in) :: t                                                          ! time
    real(kind = wp), dimension(:), intent(in) :: XX !< grid                            
    type(mms_type),intent(in) :: M

    ! Temporal Derivatives of solutions du/dt

    ! temporal derivatives of velocities
    u(1:3) = -M%nt*pi*sin(M%nt*pi*t)*sin(M%nx*pi*XX(1))*sin(M%ny*pi*XX(2))*sin(M%nz*pi*XX(3))

    ! temporal derivatives of stresses
    u(4:9) = -M%nt*pi*sin(M%nt*pi*t)*cos(M%nx*pi*XX(1))*cos(M%ny*pi*XX(2))*cos(M%nz*pi*XX(3))

  end function evaluate_mms_dt


  pure function evaluate_mms_dx(t, XX, M) result(u)

    implicit none

    real(kind = wp), dimension(9) :: u ! manufactured fields  
    real(kind = wp), intent(in) :: t                                                          ! time
    type(mms_type),intent(in) :: M
    real(kind = wp), dimension(:), intent(in) :: XX !< grid                            

    ! spatial derivatives of solutions with respect to x:  du/dx

    !  velocities: du/dx
    u(1:3) = M%nx*pi*cos(M%nt*pi*t)*cos(M%nx*pi*XX(1))*sin(M%ny*pi*XX(2))*sin(M%nz*pi*XX(3))

    !  stresses: du/dx
    u(4:9) = -M%nx*pi*cos(M%nt*pi*t)*sin(M%nx*pi*XX(1))*cos(M%ny*pi*XX(2))*cos(M%nz*pi*XX(3))

  end function evaluate_mms_dx



  pure function evaluate_mms_dy(t, XX, M) result(u)

    implicit none

    real(kind = wp), dimension(9) :: u ! manufactured fields  
    real(kind = wp), intent(in) :: t                                                          ! time
    real(kind = wp), dimension(:), intent(in) :: XX !< grid                            
    type(mms_type),intent(in) :: M

    ! spatial derivatives of solutions with respect to y:  du/dy

    !  velocities: du/dy
    u(1:3) = M%ny*pi*cos(M%nt*pi*t)*sin(M%nx*pi*XX(1))*cos(M%ny*pi*XX(2))*sin(M%nz*pi*XX(3))

    !  stresses: du/dy
    u(4:9) = -M%ny*pi*cos(M%nt*pi*t)*cos(M%nx*pi*XX(1))*sin(M%ny*pi*XX(2))*cos(M%nz*pi*XX(3))

  end function evaluate_mms_dy


  pure function evaluate_mms_dz(t, XX, M) result(u)

    implicit none

    real(kind = wp), dimension(9) :: u ! manufactured fields  
    real(kind = wp), intent(in) :: t                                                          ! time
    type(mms_type),intent(in) :: M
    real(kind = wp), dimension(:), intent(in) :: XX !< grid                            

    ! spatial derivatives of solutions with respect to z:  du/dz

    !  velocities: du/dz
    u(1:3) = M%nz*pi*cos(M%nt*pi*t)*sin(M%nx*pi*XX(1))*sin(M%ny*pi*XX(2))*cos(M%nz*pi*XX(3))

    !  stresses: du/dz
    u(4:9) = -M%nz*pi*cos(M%nt*pi*t)*cos(M%nx*pi*XX(1))*cos(M%ny*pi*XX(2))*sin(M%nz*pi*XX(3))

  end function evaluate_mms_dz


end module mms
