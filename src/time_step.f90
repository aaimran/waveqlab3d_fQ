!> @brief time_step module
!> @details this module contains derived types and routines that are used to
!> advance the solution defined on a domain in time according to some
!> partial differential equation and boundary and interface conditions

module time_step

  use common, only : wp

  implicit none

  !> @brief low storage Runge-Kutta method parameters
  type :: RK_type
     integer :: nstage ! number of stages
     real(kind = wp),dimension(:),allocatable :: A,B,C ! coefficients
  end type RK_type

contains

  subroutine init_RK(RK)

    !> @brief initializes RK_type that stores Runge-Kutta coefficients
    use common, only : wp
    implicit none

    integer,parameter :: order = 4 !< RK order of accuracy (read from input file)

    type(RK_type),intent(out) :: RK

    ! set RK coefficients

    select case(order)
    case default
       stop 'Invalid Runge-Kutta order in init_RK'
    case(1) ! explicit Euler
       RK%nstage = 1
       allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
       RK%A = (/ 0.0_wp /)
       RK%B = (/ 1.0_wp /)
       RK%C = (/ 0.0_wp,1.0_wp /)
    case(2) ! Heun
       RK%nstage = 2
       allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
       RK%A = (/ 0.0_wp,-1.0_wp /)
       RK%B = (/ 1.0_wp,1.0_wp/2.0_wp /)
       RK%C = (/ 0.0_wp,1.0_wp,1.0_wp /)
    case(3) ! Williamson (3,3)
       RK%nstage = 3
       allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
       RK%A = (/ 0.0_wp, -5.0_wp/9.0_wp, -153.0_wp/128.0_wp /)
       RK%B = (/ 1.0_wp/3.0_wp, 15.0_wp/16.0_wp, 8.0_wp/15.0_wp /)
       RK%C = (/ 0.0_wp, 1.0_wp/3.0_wp, 3.0_wp/4.0_wp, 1.0_wp /)
    case(4) ! Kennedy-Carpenter (5,4)
       RK%nstage = 5
       allocate(RK%A(RK%nstage),RK%B(RK%nstage),RK%C(0:RK%nstage))
       RK%A = (/ 0.0_wp, -567301805773.0_wp/1357537059087.0_wp, -2404267990393.0_wp/2016746695238.0_wp &
            , -3550918686646.0_wp/2091501179385.0_wp, -1275806237668.0_wp/842570457699.0_wp /)
       RK%B = (/ 1432997174477.0_wp/9575080441755.0_wp, 5161836677717.0_wp/13612068292357.0_wp &
            , 1720146321549.0_wp/2090206949498.0_wp, 3134564353537.0_wp/4481467310338.0_wp &
            , 2277821191437.0_wp/14882151754819.0_wp /)
       RK%C = (/ 0.0_wp, 1432997174477.0_wp/9575080441755.0_wp, 2526269341429.0_wp/6820363962896.0_wp &
            , 2006345519317.0_wp/3224310063776.0_wp, 2802321613138.0_wp/2924317926251.0_wp, 1.0_wp /)
    end select

  end subroutine init_RK


  subroutine time_step_RK(D,dt,RK, n)

    !> @brief advance solution in domain by one time step dt using low storage Runge-Kutta method

    use domain, only : domain_type,write_output,enforce_bound_iface_conditions, &
         exchange_fields, scale_rates, set_rates, update_fields, exchange_fields_interface, &
         write_fault_receivers_domain
    !use plastic, only : update_fields_plastic
    !use moment_tensor, only : exact_moment_tensor

    implicit none

    type(domain_type),intent(inout) :: D
    real(kind = wp),intent(in) :: dt
    type(RK_type),intent(in) :: RK
    integer :: stage,stat
    integer,intent(in) :: n
    real(kind = wp) :: t0
    integer :: nqU,nrU,nsU,nqV,nrV,nsV
    real(kind = wp), allocatable, dimension(:,:,:,:) :: u,v
    real(kind = wp), allocatable, dimension(:,:) :: M_fin

    ! initial time

    t0 = D%t

    ! Get moment tensor exact solution 

  !  if(D%B(1)%MT%use_moment_tensor) then

  !    call exact_moment_tensor(D%B(1),D%B(1)%MT,D%t) 

  !  end if

  !  if(D%B(2)%MT%use_moment_tensor) then

  !    call exact_moment_tensor(D%B(2),D%B(2)%MT,D%t) 

  !  end if

    ! loop over RK stages

    do stage = 1,RK%nstage

       call exchange_fields(D)

       ! multiply rates by RK coefficient A

       call scale_rates(D,RK%A(stage))

       ! set rates (only within blocks using elastic wave equation,
       ! does not include SAT forcing, does not set hat variables)

       call set_rates(D)

       ! prepare to enforce interface and boundary conditions

       call exchange_fields_interface(D)

       ! enforce boundary and interface conditions: set hat variables, computes SAT forcing,
       ! and adds SAT forcing to rates
       
       call enforce_bound_iface_conditions(D, stage)

       ! output fields at time t0
       ! it might seem that this could be done externally,
       ! in the main program time-step loop, before/after advancing the solution
       ! by one time step, but this would not provide "rates" (like slip velocity
       ! fault tractions) at time t0, which are set above during first stage of RK method

       if (stage==1) call write_fault_receivers_domain(D)
       if ((stage==1) .and. (mod(n-1,D%w_stride)==0)) call write_output(D)

       ! update fields

       D%t = t0 + RK%C(stage)*dt
       call update_fields(D, RK%B(stage)*dt,stage,RK%nstage)


       !if (stage == 4) stop

    end do


  end subroutine time_step_RK


end module time_step
