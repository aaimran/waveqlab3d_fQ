module plastic

  use common, only : wp
  

  implicit none

contains

  !subroutine update_fields_plastic(U,lambda,gammap,X,LambdaMu,mu_beta_eta,mx,my,mz,t,dt,response,problem)
  subroutine update_fields_plastic(F, P, G, M, dt, t,problem,response,plastic_model)

    use datatypes, only: block_type, block_fields, block_plastic, block_grid_t, block_material
    use initial_stress_condition, only : initial_stress_tensor
    
    implicit none

    !type(block_fields), intent(inout):: F
    type(block_type), intent(inout):: F
    type(block_plastic), intent(inout):: P
    type(block_grid_t), intent(in) :: G
    type(block_material), intent(in) :: M

    real(kind = wp),intent(in) :: t,dt                                  ! time t, time-step dt
    !real(kind = wp),dimension(:),intent(in):: mu_beta_eta              ! mu,beta,eta
    integer :: mx,my,mz,px,py,pz                                        ! number of grid-points in x,y,z
    integer :: ix,iy,iz,fx,fy,fz                                        ! number of grid-points in x,y,z
    character(*),intent(in) :: response,problem,plastic_model           ! response and problem types

    !real(kind = wp),dimension(:,:,:),intent(inout) :: lambda           ! plastic strain rate
    !real(kind = wp),dimension(:,:,:),intent(inout) :: gammap           ! plastic strain

    integer :: i,j,k                                                    ! grid-indices
    real(kind = wp) :: Gel,Kel                                          ! local shear and bulk modulus
    real(kind = wp) :: s0(6)                                            ! local initial stress tensor
    real(kind = wp) :: mu, beta, eta
    integer :: n, nx, ny, nz, npml

    nx = G%C%nq
    ny = G%C%nr
    nz = G%C%ns

    mx = G%C%mq
    my = G%C%mr
    mz = G%C%ms
    px = G%C%pq
    py = G%C%pr
    pz = G%C%ps

    if (response/='plastic') return

    mu = P%mu_beta_eta(1)
    beta = P%mu_beta_eta(2)
    eta = P%mu_beta_eta(3)
    !
    !mu = 5.734623443633283e-1_wp 
    !beta = mu/2.0_wp
    !eta = 2.774826789838337e-1_wp !.2775_wp 
    !
    !print *, P%mu_beta_eta

    !STOP
    !print *,  mu,  beta, eta
    !STOP
    !
    ! PML can trigger plasticity and should be avoided
    iz = mz
    fz = pz
    iy = my
    fy = py
    ix = mx
    fx = px
    if (F%PMLB(5)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(5)%N_pml
       if ((mz == 1) .or. (pz .le. npml)) iz = npml + 1
                 
    end if
    
    if (F%PMLB(6)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(6)%N_pml
       if ((pz == nz) .or. (mz .ge. (nz-(npml-1)))) fz = pz - npml
       
    end if
    
    if (F%PMLB(3)%pml .EQV. .TRUE.) then

       npml = F%PMLB(3)%N_pml
       if ((my == 1) .or. (py .le. npml)) iy = npml + 1

    end if

    if (F%PMLB(4)%pml .EQV. .TRUE.) then

       npml = F%PMLB(4)%N_pml
       if ((py == ny) .or. (my .ge. (ny-(npml-1)))) fy = py - npml
       
    end if

    if (F%PMLB(1)%pml .EQV. .TRUE.) then

       npml = F%PMLB(1)%N_pml
       if ((mx == 1) .or. (px .le. npml)) ix = npml + 1

    end if

    if (F%PMLB(2)%pml .EQV. .TRUE.) then

       npml = F%PMLB(2)%N_pml
       if ((px == nx) .or. (mx .ge. (nx-(npml-1)))) fx = px - npml
       
    end if

    if (plastic_model .eq. 'scec') then
        
       do k=iz,fz
          do j=iy,fy
             do i=ix,fx
                !
                call initial_stress_tensor(s0,G%X(i,j,k,1:3),problem)
                !
                call plastic_flow2(F%F%F(i,j,k,4:9),s0,G%X(i,j,k,2),problem,dt)
                !
             end do
          end do
       end do

    else

       do k=iz,fz
          do j=iy,fy
             do i=ix,fx
                !                                                                                                                                                                   
                Gel = M%M(i,j,k,2)
                !                                                                                                                                                                   
                Kel = M%M(i,j,k,1) + 2.0_wp/3.0_wp*M%M(i,j,k,2)
                !                                                                                                                                                                   
                call initial_stress_tensor(s0,G%X(i,j,k,1:3),problem)
                !                                                                                                                                                                   
                call plastic_flow(F%F%F(i,j,k,4:9),P%P(i,j,k,2),P%P(i,j,k,1),mu,beta,eta,Gel,Kel,s0,G%X(i,j,k,2),problem,dt)
                !
             end do
          end do
       end do
       
    end if

       

  end subroutine update_fields_plastic

  subroutine plastic_flow2(s,s0,y_ijk,problem,dt)

    ! The return map algorithm for computation of plastic flow                                                                                                                 
    implicit none

    ! lambda = plastic strain rate in current time step (or RK stage)  
    ! s0 = initial stress            
    ! s = stress change (added to s0 to get absolute stress)
    ! sa = absolute stress tensor
    ! sd = deviatoric stress tensor  

    real(kind = wp),intent(inout) :: s(6)
    real(kind = wp),intent(in) :: s0(6),dt,y_ijk
    character(*),intent(in) :: problem             ! problem type
    real(kind = wp) :: Yf, Ys
    real(kind = wp) :: c, nu, Tv
    real(kind = wp) :: Pf,r
    real(kind = wp) :: tau,sigma
    real(kind = wp),dimension(6) :: sa,sd
    real(kind = wp),dimension(6),parameter :: delta = (/1d0,1d0,1d0,0d0,0d0,0d0/) ! Kronecker delta 

    c = 1.36d0                        ! cohesion
    nu = 0.1934d0                     ! bulk friction
    Tv = 0.03d0                       ! viscoplastic relaxation time

    if (problem == 'TPV30') then
       c = 1.18d0                        ! cohesion                                                                                                                           
       nu = 0.1680d0                     ! bulk friction                                                                                                                        
       Tv = 0.05d0                       ! viscoplastic relaxation time  
    end if

    Pf = 9.8d0*y_ijk                  ! pore pressure

    ! compute components of the stress tensor
    call stress(s,s0,sa,tau,sigma)

    ! compute trial Drucker-Prager yield stress
    Ys = max(0d0, c*cos(atan(nu))-(sigma)*sin(atan(nu)))

    ! compute trial Drucker-Prager yield function
    Yf = tau - Ys

    ! check if trial update violates yield condition 
    if (Yf<=0d0) then ! no plastic flow
       return
    end if

    ! plastic flow:
    r = exp(-dt/Tv) + (1d0-exp(-dt/Tv))*Ys/tau

    ! compute deviatoric stress (trial)
    sd = sa-sigma*delta

    ! update absolute stress sa
    sa = sigma*delta + r*sd

    ! extract s from s = sa-s0
    s(1) = sa(1)-s0(1) ! sxx 
    s(2) = sa(2)-s0(2) ! syy                         
    s(3) = sa(3)-s0(3) ! szz
    s(4) = sa(4)-s0(4) ! sxy
    s(5) = sa(5)-s0(5) ! sxz              
    s(6) = sa(6)-s0(6) ! syz

  end subroutine plastic_flow2


  !subroutine plastic_flow(lambda,gammap,mu,beta,eta,G,K,s,s0,dt)
  subroutine plastic_flow(s,gammap,lambda,mu,beta,eta,G,K,s0,y_ijk,problem,dt)
    !subroutine plastic_flow(s,mu,beta,eta,G,K,s0,y_ijk,problem,dt)
    ! A general algorithm for computation of plastic flow

    implicit none


    ! lambda = plastic strain rate in current time step (or RK stage)                                                                   
    ! s0 = initial stress
    ! s = stress change (added to s0 to get absolute stress)                                                                            

    !real(kind = wp),intent(inout) :: lambda,gammap,s(6)
    real(kind = wp),intent(inout) :: s(6)
    real(kind = wp),intent(inout):: lambda,gammap
    real(kind = wp),intent(in) :: G,K,s0(6),dt,y_ijk
    real(kind = wp),intent(in) :: mu, beta, eta
    character(*),intent(in) :: problem             ! problem type  

    ! sa = absolute stress tensor
    ! sd = deviatoric stress tensor                                                                                                     

    real(kind = wp) :: tau,sigma,Y
    real(kind = wp),dimension(6) :: sa,sd
    real(kind = wp),dimension(6),parameter :: delta = (/1d0,1d0,1d0,0d0,0d0,0d0/) ! Kronecker delta                                     

    ! compute components of the trial stress tensor, sa, tau, sigma
    call stress(s,s0,sa,tau,sigma)

    ! compute the yield function Y
    call yield(tau,sigma,mu,Y)

    ! check if trial update violates yield condition   
    if (Y<=0d0) then ! no plastic flow
       lambda = 0d0
       return
    end if

    ! plastic flow:     
    ! compute deviatoric stress (trial)                                                                                                                                   
    sd = sa-sigma*delta

    ! implicit solution for lambda,s to    
    ! eta*lambda = Y(s)

    ! rate-independent limit is eta = 0
    ! solution for lambda in closed-form:                                                           

    ! update lambda and gammap
    lambda = (tau+mu*sigma)/(eta+dt*(G+mu*beta*K))
    gammap = gammap+dt*lambda

    ! update tau and sigma
    tau = tau-dt*lambda*G
    sigma = sigma-dt*lambda*beta*K

    ! update deviatoric stress
    sd = sd*tau/(tau+dt*lambda*G)

    ! update absolute stress
    sa = sigma*delta+sd

    ! extract stress changes from s = sa-s0    
    s(1) = sa(1)-s0(1) ! sxx
    s(2) = sa(2)-s0(2) ! syy
    s(3) = sa(3)-s0(3) ! szz
    s(4) = sa(4)-s0(4) ! sxy
    s(5) = sa(5)-s0(5) ! sxz
    s(6) = sa(6)-s0(6) ! syz

  end subroutine plastic_flow


  subroutine invariants(sxx,syy,szz,sxy,sxz,syz,tau,sigma)

    implicit none

    real(kind = wp),intent(in) :: sxx,syy,szz,sxy,sxz,syz
    real(kind = wp),intent(out) :: tau,sigma

    real(kind = wp) :: J2

    sigma = (sxx+syy+szz)/3d0                          
    J2 = ((sxx-syy)**2+(syy-szz)**2+(szz-sxx)**2)/6d0+ &
         sxy**2+sxz**2+syz**2
    tau = sqrt(J2)

  end subroutine invariants

  subroutine stress(s,s0,sa,tau,sigma)

    ! s: stress perturbation
    ! s0: prestress
    ! sa: absolute stress
    ! sigma: mean stress                                             
    ! tau: second invariant of the deviatoric stress tensor 

    implicit none

    real(kind = wp),intent(in) :: s(6),s0(6)
    real(kind = wp),intent(out) :: sa(6),tau,sigma

    ! s0 = (/sxx0,syy0,szz0,sxy0,sxz0,syz0 /)
    !           1    2    3    4    5    6

    ! compute sa
    sa(1) = s(1)+s0(1)
    sa(2) = s(2)+s0(2)
    sa(3) = s(3)+s0(3)
    sa(4) = s(4)+s0(4)
    sa(5) = s(5)+s0(5)
    sa(6) = s(6)+s0(6)

    ! compute tau,sigma
    call invariants(sa(1),sa(2),sa(3),sa(4),sa(5),sa(6),tau,sigma)

  end subroutine stress

  subroutine yield(tau,sigma,mu,Y)

    ! Y:  yield function
    ! sigma: mean stress
    ! mu: ?
    ! tau: second invariant of the deviatoric stress tensor

    implicit none

    real(kind = wp),intent(in) :: tau,sigma,mu
    real(kind = wp),intent(out) :: Y

    Y = tau+mu*sigma

  end subroutine yield



end module plastic
