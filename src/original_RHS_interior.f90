module RHS_Interior

  use common, only : wp

contains


    subroutine RHS_Center(F, G, M, type_of_mesh)

    use datatypes, only: block_type,block_fields, block_grid_t, block_material
    use JU_xJU_yJU_z6, only : JJU_x6_interior

    implicit none

    !type(block_fields), intent(inout):: F
    type(block_type), intent(inout):: F
    type(block_grid_t), intent(in) :: G
    type(block_material), intent(in) :: M

    character(len=64), intent(in) :: type_of_mesh

   
    ! compute all spatial derivatives and add interior rates to rates array, no forcing
    call JJU_x6_interior(F, G, M, type_of_mesh)

  end subroutine RHS_center


  subroutine RHS_near_boundaries(F, G, M)

    use common, only : wp
    use datatypes, only: block_type, block_fields, block_grid_t, block_material,block_pml
    use JU_xJU_yJU_z6, only : JJU_x6

    implicit none

    !type(block_fields), intent(inout):: F
    type(block_type), intent(inout):: F
    type(block_grid_t), intent(in) :: G
    type(block_material), intent(in) :: M

    integer :: mbx, mby, mbz, pbx, pby, pbz, &
               mx,my,mz, px, py, pz, nx, ny, nz, n ! number of grid points and fields
    real(kind = wp) :: hq, hr, hs                                          ! spatial steps

    ! work arrays
    real(kind = wp), dimension(:), allocatable :: Ux, Uy, Uz ! to hold derivatives
    real(kind = wp), dimension(:), allocatable :: DU  ! to hold rates
    integer :: x, y, z
    real(kind = wp) :: rhoJ_inv, lambda2mu
    integer :: ix, iy, iz,  fx, fy, fz
    
    
    ! compute all spatial derivatives and add interior rates to rates array, no forcing
    !
    ! this routine can be used for any point in the interior, both near boundaries and away from them
    !
    ! here, it is used only for points near boundaries and an optimized version is used above
    ! for points away from boundaries, where FD stencil is identical at every point

    mbx = G%C%mbq
    mby = G%C%mbr
    mbz = G%C%mbs
    pbx = G%C%pbq
    pby = G%C%pbr
    pbz = G%C%pbs

    mx = G%C%mq
    my = G%C%mr
    mz = G%C%ms
    px = G%C%pq
    py = G%C%pr
    pz = G%C%ps
    nx = G%C%nq
    ny = G%C%nr
    nz = G%C%ns

    n = size(F%F%F, 4)
    
    hq = G%hq
    hr = G%hr
    hs = G%hs
    
    allocate(Ux(n), Uy(n), Uz(n), DU(n))
    
    ix = 6
    iy = 6
    iz = 6

    fx = 6
    fy = 6
    fz = 6

    if (F%PMLB(5)%pml .EQV. .TRUE.) iz = max(F%PMLB(5)%N_pml, 6)
    if (F%PMLB(6)%pml .EQV. .TRUE.) fz = max(F%PMLB(6)%N_pml, 6)
    
    if (F%PMLB(3)%pml .EQV. .TRUE.) iy = max(F%PMLB(3)%N_pml, 6)
    if (F%PMLB(4)%pml .EQV. .TRUE.) fy = max(F%PMLB(4)%N_pml, 6)
    
    if (F%PMLB(1)%pml .EQV. .TRUE.) ix = max(F%PMLB(1)%N_pml, 6)
    if (F%PMLB(2)%pml .EQV. .TRUE.) fx = max(F%PMLB(2)%N_pml, 6)
    
    
    do z = mz, pz
       do y = my, py
          do x = mx, px
             
             ! cycle if an interior point
             if (((iz+1 <= z) .and. (z <= nz - fz)) .and. &
                  ((iy+1 <= y) .and. (y <= ny - fy)) .and. &
                  ((ix+1 <= x) .and. (x <= nx - fx))) cycle
             
             ! compute spatial derivatives, 6th order accuracy
             call JJU_x6(x, y, z, F%F, G, G%metricx, Ux)
             call JJU_x6(x, y, z, F%F, G, G%metricy, Uy)
             call JJU_x6(x, y, z, F%F, G, G%metricz, Uz)
             
             rhoJ_inv = 1.0_wp/(M%M(x,y,z,3)*G%J(x,y,z)) ! 1/(rho*J)
             lambda2mu = M%M(x,y,z,1) + 2.0_wp*M%M(x,y,z,2) ! lambda+2*mu
             
            
             
             ! compute rates
             DU(1) = (Ux(4) + Uy(7) + Uz(8))*rhoJ_inv
             DU(2) = (Ux(7) + Uy(5) + Uz(9))*rhoJ_inv
             DU(3) = (Ux(8) + Uy(9) + Uz(6))*rhoJ_inv

             DU(4) = lambda2mu*Ux(1) + M%M(x,y,z,1)*(Uy(2) + Uz(3))
             DU(5) = lambda2mu*Uy(2) + M%M(x,y,z,1)*(Ux(1) + Uz(3))
             DU(6) = lambda2mu*Uz(3) + M%M(x,y,z,1)*(Ux(1) + Uy(2))
             
             DU(7) = M%M(x,y,z,2)*(Uy(1) + Ux(2))
             DU(8) = M%M(x,y,z,2)*(Uz(1) + Ux(3))
             DU(9) = M%M(x,y,z,2)*(Uz(2) + Uy(3))
             
             ! compute pml rhs to be used in runge-kutta time-step
             ! and append pml auxiliary functions to elastic-rates DU
             call PML(F, M, G, DU, Ux, Uy, Uz, x,y,z,n,rhoJ_inv,lambda2mu)
             
             ! add new rates to rates arrays
             F%F%DF(x,y,z,1) = F%F%DF(x,y,z,1) + DU(1)
             F%F%DF(x,y,z,2) = F%F%DF(x,y,z,2) + DU(2)
             F%F%DF(x,y,z,3) = F%F%DF(x,y,z,3) + DU(3)
             F%F%DF(x,y,z,4) = F%F%DF(x,y,z,4) + DU(4)
             F%F%DF(x,y,z,5) = F%F%DF(x,y,z,5) + DU(5)
             F%F%DF(x,y,z,6) = F%F%DF(x,y,z,6) + DU(6)
             F%F%DF(x,y,z,7) = F%F%DF(x,y,z,7) + DU(7)
             F%F%DF(x,y,z,8) = F%F%DF(x,y,z,8) + DU(8)
             F%F%DF(x,y,z,9) = F%F%DF(x,y,z,9) + DU(9)
             
          end do
       end do
    end do

  end subroutine RHS_near_boundaries

  subroutine impose_boundary_condition(B, mms_vars, t)

    use datatypes, only : block_type, mms_type, boundary_type
    use BoundaryConditions, only : BC_Lx, BC_Ly, BC_Lz, BC_Rx, BC_Ry, BC_Rz

    implicit none

    type(block_type), intent(inout) :: B
    type(mms_type), intent(in) :: mms_vars                        ! parameters for MMS
    real(kind = wp), intent(in) :: t                                         ! time

    integer :: mx, my, mz, px, py, pz, n
    type(boundary_type) :: boundary_vars
    real(kind = wp) :: tau0, hx, hy, hz                          ! penalty the in x-direction, y-direction, z-direction

    real(kind = wp), dimension(:,:,:), allocatable, save :: U_x, U_y, U_z     ! to hold boundary forcing

    integer :: x, y, z

    boundary_vars = B%boundary_vars
    tau0 = B%tau0

    hx = B%G%hq
    hy = B%G%hr
    hz = B%G%hs

    mx = B%G%C%mq
    my = B%G%C%mr
    mz = B%G%C%ms
    px = B%G%C%pq
    py = B%G%C%pr
    pz = B%G%C%ps

    n = size(B%F%DF,4)

    ! initialize work array

    if (.not.allocated(U_x)) allocate(U_x(my:py,mz:pz,n))
    if (.not.allocated(U_y)) allocate(U_y(mx:px,mz:pz,n))
    if (.not.allocated(U_z)) allocate(U_z(mx:px,my:py,n))

    ! construct boundary forcing
    if (boundary_vars%Lx > 0) then
      U_x = BC_Lx(B, boundary_vars%Lx, mms_vars, t)
       do z = mz, pz
          do y = my, py
             B%F%DF(mx, y, z, 1:n) = B%F%DF(mx, y, z, 1:n) - tau0/(hx)*U_x(y, z, 1:n)
             if(B%PMLB(1)%pml .EQV. .TRUE.) then
                B%PMLB(1)%DQ(mx, y, z, 1:n) = B%PMLB(1)%DQ(mx, y, z, 1:n) - tau0/(hx)*U_x(y, z, 1:n)
                
             end if
          end do
       end do
    end if

    if (boundary_vars%Rx > 0) then
       U_x = BC_Rx(B, boundary_vars%Rx, mms_vars, t)
       do z = mz, pz
          do y = my, py
             B%F%DF(px, y, z, 1:n) = B%F%DF(px, y, z, 1:n) - tau0/(hx)*U_x(y, z, 1:n)
             if(B%PMLB(2)%pml .EQV. .TRUE.) then
                B%PMLB(2)%DQ(px, y, z, 1:n) = B%PMLB(2)%DQ(px, y, z, 1:n) - tau0/(hx)*U_x(y, z, 1:n)
              
             end if
          end do
       end do
    end if


    if (boundary_vars%Ly > 0) then
       U_y = BC_Ly(B, boundary_vars%Ly, mms_vars, t)
       do z = mz, pz
          do x = mx, px
             B%F%DF(x, my, z, 1:n) = B%F%DF(x, my, z, 1:n) - tau0/(hy)*U_y(x, z, 1:n)
             if(B%PMLB(3)%pml .EQV. .TRUE.) then
                B%PMLB(3)%DQ(x, my, z, 1:n) =  B%PMLB(3)%DQ(x, my, z, 1:n) - tau0/(hy)*U_y(x, z, 1:n)
              end if
          end do
       end do
    end if

    if (boundary_vars%Ry > 0) then
       U_y = BC_Ry(B, boundary_vars%Ry, mms_vars, t)
       do z = mz, pz
          do x = mx, px
             B%F%DF(x, py, z, 1:n) = B%F%DF(x, py, z, 1:n) - (tau0/(hy))*U_y(x, z, 1:n)
             if(B%PMLB(4)%pml .EQV. .TRUE.) then
                B%PMLB(4)%DQ(x, py, z, 1:n) =  B%PMLB(4)%DQ(x, py, z, 1:n) - tau0/(hy)*U_y(x, z, 1:n)
              end if
          end do
       end do
    end if


    if (boundary_vars%Lz > 0) then
       U_z = BC_Lz(B, boundary_vars%Lz, mms_vars, t)
       do y = my, py
          do x = mx, px
             B%F%DF(x, y, mz, 1:n) = B%F%DF(x, y, mz, 1:n) - (tau0/(hz))*U_z(x, y, 1:n)
             if(B%PMLB(5)%pml .EQV. .TRUE.) then
                B%PMLB(5)%DQ(x, y, mz, 1:n) = B%PMLB(5)%DQ(x, y, mz, 1:n) - (tau0/(hz))*U_z(x, y, 1:n)
             end if
          end do
       end do
    end if

    if (boundary_vars%Rz > 0) then
       U_z = BC_Rz(B, boundary_vars%Rz, mms_vars, t)
       do y = my, py
          do x = mx, px
             B%F%DF(x, y, pz, 1:n) = B%F%DF(x, y, pz, 1:n) - (tau0/(hz))*U_z(x, y, 1:n)
             if(B%PMLB(6)%pml .EQV. .TRUE.) then
                B%PMLB(6)%DQ(x, y, pz, 1:n) = B%PMLB(6)%DQ(x, y, pz, 1:n) - (tau0/hz)*U_z(x, y, 1:n)
             end if
          end do
       end do
    end if

  end subroutine Impose_Boundary_Condition


  subroutine Impose_Interface_Condition(problem, coupling, I, B, ib, &
                                        t, stage, mms_vars, handles)

    use mpi3dbasic, only : rank
    use datatypes, only : block_type, iface_type, mms_type
    use CouplingForcing, only :  Couple_Interface_x
    use fault_output, only : fault_type

    ! compute hat variables, SAT forcing terms, and add SAT forcing to rates
    
    implicit none

    character(256), intent(in) :: problem
    character(*),intent(in) :: coupling
    integer, intent(in) :: stage, ib
    type(iface_type),intent(inout) :: I
    type(block_type),intent(inout) :: B
    type(mms_type), intent(inout) :: mms_vars
    type(fault_type), intent(inout) :: handles
    real(kind = wp),intent(in) :: t
    integer :: mx, my, mz, px, py, pz, n             ! number of grid points and field variables

    real(kind = wp) :: tau0, hx
    real(kind = wp), dimension(:,:,:), allocatable, save :: F_x, G_x  ! to hold boundary forcing in the negative directions
    integer :: y, z

    mx = B%G%C%mq
    px = B%G%C%pq
    my = B%G%C%mr
    py = B%G%C%pr
    mz = B%G%C%ms
    pz = B%G%C%ps

    hx = B%G%hq

    tau0 = B%tau0

    n = size(B%F%DF, 4)
    if (.not.allocated(F_x)) allocate(F_x(my:py, mz:pz, n))
    if (.not.allocated(G_x)) allocate(G_x(my:py, mz:pz, n))

    ! Couple_Interface_x computes hat variables (by solving friction law or interface conditions)
    ! and also computes SAT forcing terms, but it does not add SAT forcing to rates
    
    if (ib == 2) then
      call Couple_Interface_x(F_x, G_x, problem, coupling, I, B%G, B%M,&
                            B%B(ib)%F, B%B(ib)%Fopp, B%B(ib)%M, B%B(ib)%Mopp, &
                            B%B(ib)%n_l, B%B(ib)%n_m, B%B(ib)%n_n, &
                            ib, t, stage, mms_vars, handles)
     ! then add SAT forcing to rates
        do z = mz, pz
           do y = my, py
            B%F%DF(px, y, z, 1:n) = B%F%DF(px, y, z, 1:n) - tau0/hx*F_x(y, z, 1:n)
            end do
        end do

      else if (ib == 1) then
       call Couple_Interface_x(F_x, G_x, problem, coupling, I, B%G, B%M, &
                            B%B(ib)%Fopp, B%B(ib)%F, B%B(ib)%Mopp, B%B(ib)%M, &
                            B%B(ib)%n_l, B%B(ib)%n_m, B%B(ib)%n_n, &
                            ib, t, stage, mms_vars, handles)
     ! then add SAT forcing to rates
      do z = mz, pz
         do y = my, py
          B%F%DF(mx, y, z, 1:n) = B%F%DF(mx, y, z, 1:n) - tau0/hx*G_x(y, z, 1:n)
          end do
      end do

     end if

!      print *, ib, G_x(11,11,:)
     

   end subroutine Impose_Interface_Condition


   subroutine PML(F, M, G, DU, Ux, Uy, Uz, x,y,z, n, rhoJ_inv,lambda2mu)

     use common, only : wp
     use datatypes, only: block_type, block_grid_t, block_material, block_pml
     use JU_xJU_yJU_z6, only : JJU_x6
     
     implicit none
          
     type(block_type), intent(inout):: F
     type(block_grid_t), intent(in) :: G
     type(block_material), intent(in) :: M
     real(kind = wp),intent(in) :: rhoJ_inv, lambda2mu
     integer, intent(in) :: n                      ! number of  fields
     integer, intent(in) :: x, y, z
     
     ! work arrays
     real(kind = wp), intent(in) :: Ux(n), Uy(n), Uz(n) ! to hold derivatives
     real(kind = wp),intent(inout) :: DU(n)
    
    
     integer :: nx, ny, nz                          ! number of grid points
     real(kind = wp) :: hq, hr, hs                  ! spatial steps

     real(kind = wp) :: lam, mu, rho
   
   
   
    
    ! pml parameters
    integer :: npml
    real(kind = wp) :: cfs, dx, dy, dz, d0x, d0y, d0z, Rf                

    nx = G%C%nq
    ny = G%C%nr
    nz = G%C%ns

    !n = size(F%F%F, 4)

    hq = G%hq
    hr = G%hr
    hs = G%hs
    
    
    cfs = 0.1_wp                     !complex frequency shift
    Rf = 1.0e-3_wp                   !relative pml error


    lam = M%M(x,y,z,1)
    mu  = M%M(x,y,z,2)
    rho = M%M(x,y,z,3)

    d0x = 6.0_wp*4.0_wp/(2.0_wp*1.0_wp)*log(1.0_wp/Rf)      !damping strength
    d0y = 6.0_wp*4.0_wp/(2.0_wp*1.0_wp)*log(1.0_wp/Rf)      !damping strength
    d0z = 6.0_wp*4.0_wp/(2.0_wp*1.0_wp)*log(1.0_wp/Rf)      !damping strength

   
    
    ! z-dependent layer
    if(F%PMLB(5)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(5)%N_pml
       
       if (z .le. npml) then
          
          dz = d0z*(abs(real(z-npml)/real(npml)))**3

          DU(1) = DU(1) - dz*G%J(x,y,z)*F%PMLB(5)%Q(x,y,z,1)*rhoJ_inv
          DU(2) = DU(2) - dz*G%J(x,y,z)*F%PMLB(5)%Q(x,y,z,2)*rhoJ_inv
          DU(3) = DU(3) - dz*G%J(x,y,z)*F%PMLB(5)%Q(x,y,z,3)*rhoJ_inv
          
          DU(4) = DU(4) - dz*(lambda2mu*F%PMLB(5)%Q(x,y,z,4) &
               + lam*F%PMLB(5)%Q(x,y,z,5) &
               + lam*F%PMLB(5)%Q(x,y,z,6))
          
          DU(5) = DU(5) - dz*(lam*F%PMLB(5)%Q(x,y,z,4) &
               + lambda2mu*F%PMLB(5)%Q(x,y,z,5) &
               + lam*F%PMLB(5)%Q(x,y,z,6))
          
          DU(6) = DU(6) - dz*(lam*F%PMLB(5)%Q(x,y,z,4) &
               + lam*F%PMLB(5)%Q(x,y,z,5) &
               + lambda2mu*F%PMLB(5)%Q(x,y,z,6))
          
          
          DU(7) = DU(7) - mu*dz*F%PMLB(5)%Q(x,y,z,7)
          DU(8) = DU(8) - mu*dz*F%PMLB(5)%Q(x,y,z,8)
          DU(9) = DU(9) - mu*dz*F%PMLB(5)%Q(x,y,z,9)
          
          F%PMLB(5)%DQ(x,y,z,1) = F%PMLB(5)%DQ(x,y,z,1) + Uz(8)/G%J(x,y,z) &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,1)
          F%PMLB(5)%DQ(x,y,z,2) = F%PMLB(5)%DQ(x,y,z,2) + Uz(9)/G%J(x,y,z) &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,2)
          F%PMLB(5)%DQ(x,y,z,3) = F%PMLB(5)%DQ(x,y,z,3) + Uz(6)/G%J(x,y,z) &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,3)
          
          F%PMLB(5)%DQ(x,y,z,4) = F%PMLB(5)%DQ(x,y,z,4) + 0d0 &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,4)
          
          F%PMLB(5)%DQ(x,y,z,5) = F%PMLB(5)%DQ(x,y,z,5) + 0d0 &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,5)

          F%PMLB(5)%DQ(x,y,z,6) = F%PMLB(5)%DQ(x,y,z,6) + Uz(3) &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,6)

          F%PMLB(5)%DQ(x,y,z,7) = F%PMLB(5)%DQ(x,y,z,7) + 0d0 &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,7)

          F%PMLB(5)%DQ(x,y,z,8) = F%PMLB(5)%DQ(x,y,z,8) + Uz(1) &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,8)
          
          F%PMLB(5)%DQ(x,y,z,9) = F%PMLB(5)%DQ(x,y,z,9) + Uz(2) &
               - (dz+cfs)*F%PMLB(5)%Q(x,y,z,9)
          
       end if
    end if
    
    
    if(F%PMLB(6)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(6)%N_pml
       
       if (z .ge. nz-npml+1) then
          
          dz = d0z*(abs(real(z-(nz-npml))/real(npml)))**3

          DU(1) = DU(1) - dz*G%J(x,y,z)*F%PMLB(6)%Q(x,y,z,1)*rhoJ_inv
          DU(2) = DU(2) - dz*G%J(x,y,z)*F%PMLB(6)%Q(x,y,z,2)*rhoJ_inv
          DU(3) = DU(3) - dz*G%J(x,y,z)*F%PMLB(6)%Q(x,y,z,3)*rhoJ_inv
          
          DU(4) = DU(4) - dz*(lambda2mu*F%PMLB(6)%Q(x,y,z,4) &
               + lam*F%PMLB(6)%Q(x,y,z,5) &
               + lam*F%PMLB(6)%Q(x,y,z,6))
          
          DU(5) = DU(5) - dz*(lam*F%PMLB(6)%Q(x,y,z,4) &
               + lambda2mu*F%PMLB(6)%Q(x,y,z,5) &
               + lam*F%PMLB(6)%Q(x,y,z,6))
          
          DU(6) = DU(6) - dz*(lam*F%PMLB(6)%Q(x,y,z,4) &
               + lam*F%PMLB(6)%Q(x,y,z,5) &
               + lambda2mu*F%PMLB(6)%Q(x,y,z,6))
          
          
          DU(7) = DU(7) - mu*dz*F%PMLB(6)%Q(x,y,z,7)
          DU(8) = DU(8) - mu*dz*F%PMLB(6)%Q(x,y,z,8)
          DU(9) = DU(9) - mu*dz*F%PMLB(6)%Q(x,y,z,9)
          
          F%PMLB(6)%DQ(x,y,z,1) = F%PMLB(6)%DQ(x,y,z,1) + Uz(8)/G%J(x,y,z) &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,1)
          F%PMLB(6)%DQ(x,y,z,2) = F%PMLB(6)%DQ(x,y,z,2) + Uz(9)/G%J(x,y,z) &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,2)
          F%PMLB(6)%DQ(x,y,z,3) = F%PMLB(6)%DQ(x,y,z,3) + Uz(6)/G%J(x,y,z) &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,3)
          
          F%PMLB(6)%DQ(x,y,z,4) = F%PMLB(6)%DQ(x,y,z,4) + 0d0 &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,4)
          
          F%PMLB(6)%DQ(x,y,z,5) = F%PMLB(6)%DQ(x,y,z,5) + 0d0 &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,5)

          F%PMLB(6)%DQ(x,y,z,6) = F%PMLB(6)%DQ(x,y,z,6) + Uz(3) &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,6)

          F%PMLB(6)%DQ(x,y,z,7) = F%PMLB(6)%DQ(x,y,z,7) + 0d0 &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,7)

          F%PMLB(6)%DQ(x,y,z,8) = F%PMLB(6)%DQ(x,y,z,8) + Uz(1) &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,8)
          
          F%PMLB(6)%DQ(x,y,z,9) = F%PMLB(6)%DQ(x,y,z,9) + Uz(2) &
               - (dz+cfs)*F%PMLB(6)%Q(x,y,z,9)
          
          
       end if
    end if
    !=====================================================================================================
    
    ! y-dependent layer
    if(F%PMLB(3)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(3)%N_pml
       
       if (y .le. npml) then
          
          dy = d0y*(abs(real(y-(1+npml))/real(npml)))**3

          DU(1) = DU(1) - dy*G%J(x,y,z)*F%PMLB(3)%Q(x,y,z,1)*rhoJ_inv
          DU(2) = DU(2) - dy*G%J(x,y,z)*F%PMLB(3)%Q(x,y,z,2)*rhoJ_inv
          DU(3) = DU(3) - dy*G%J(x,y,z)*F%PMLB(3)%Q(x,y,z,3)*rhoJ_inv
          
          DU(4) = DU(4) - dy*(lambda2mu*F%PMLB(3)%Q(x,y,z,4) &
               + lam*F%PMLB(3)%Q(x,y,z,5) &
               + lam*F%PMLB(3)%Q(x,y,z,6))
          
          DU(5) = DU(5) - dy*(lam*F%PMLB(3)%Q(x,y,z,4) &
               + lambda2mu*F%PMLB(3)%Q(x,y,z,5) &
               + lam*F%PMLB(3)%Q(x,y,z,6))
          
          DU(6) = DU(6) - dy*(lam*F%PMLB(3)%Q(x,y,z,4) &
               + lam*F%PMLB(3)%Q(x,y,z,5) &
               + lambda2mu*F%PMLB(3)%Q(x,y,z,6))
          
          
          DU(7) = DU(7) - mu*dy*F%PMLB(3)%Q(x,y,z,7)
          DU(8) = DU(8) - mu*dy*F%PMLB(3)%Q(x,y,z,8)
          DU(9) = DU(9) - mu*dy*F%PMLB(3)%Q(x,y,z,9)
          
          F%PMLB(3)%DQ(x,y,z,1) = F%PMLB(3)%DQ(x,y,z,1) + Uy(7)/G%J(x,y,z) &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,1)
          F%PMLB(3)%DQ(x,y,z,2) = F%PMLB(3)%DQ(x,y,z,2) + Uy(5)/G%J(x,y,z) &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,2)
          F%PMLB(3)%DQ(x,y,z,3) = F%PMLB(3)%DQ(x,y,z,3) + Uy(9)/G%J(x,y,z) &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,3)
          
          F%PMLB(3)%DQ(x,y,z,4) = F%PMLB(3)%DQ(x,y,z,4) + 0d0 &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,4)
          
          F%PMLB(3)%DQ(x,y,z,5) = F%PMLB(3)%DQ(x,y,z,5) + Uy(2) &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,5)

          F%PMLB(3)%DQ(x,y,z,6) = F%PMLB(3)%DQ(x,y,z,6) + 0d0 &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,6)

          F%PMLB(3)%DQ(x,y,z,7) = F%PMLB(3)%DQ(x,y,z,7) + Uy(1) &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,7)

          F%PMLB(3)%DQ(x,y,z,8) = F%PMLB(3)%DQ(x,y,z,8) + 0d0 &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,8)
          
          F%PMLB(3)%DQ(x,y,z,9) = F%PMLB(3)%DQ(x,y,z,9) + Uy(3) &
               - (dy+cfs)*F%PMLB(3)%Q(x,y,z,9)
          
          
          
       end if
    end if
    
    
    if(F%PMLB(4)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(4)%N_pml
       
       if (y .ge. ny-npml+1) then
          
          dy = d0y*(abs(real(y-(ny-npml))/real(npml)))**3

          DU(1) = DU(1) - dy*G%J(x,y,z)*F%PMLB(4)%Q(x,y,z,1)*rhoJ_inv
          DU(2) = DU(2) - dy*G%J(x,y,z)*F%PMLB(4)%Q(x,y,z,2)*rhoJ_inv
          DU(3) = DU(3) - dy*G%J(x,y,z)*F%PMLB(4)%Q(x,y,z,3)*rhoJ_inv
          
          DU(4) = DU(4) - dy*(lambda2mu*F%PMLB(4)%Q(x,y,z,4) &
               + lam*F%PMLB(4)%Q(x,y,z,5) &
               + lam*F%PMLB(4)%Q(x,y,z,6))
          
          DU(5) = DU(5) - dy*(lam*F%PMLB(4)%Q(x,y,z,4) &
               + lambda2mu*F%PMLB(4)%Q(x,y,z,5) &
               + lam*F%PMLB(4)%Q(x,y,z,6))
          
          DU(6) = DU(6) - dy*(lam*F%PMLB(4)%Q(x,y,z,4) &
               + lam*F%PMLB(4)%Q(x,y,z,5) &
               + lambda2mu*F%PMLB(4)%Q(x,y,z,6))
          
          
          DU(7) = DU(7) - mu*dy*F%PMLB(4)%Q(x,y,z,7)
          DU(8) = DU(8) - mu*dy*F%PMLB(4)%Q(x,y,z,8)
          DU(9) = DU(9) - mu*dy*F%PMLB(4)%Q(x,y,z,9)
          
          F%PMLB(4)%DQ(x,y,z,1) = F%PMLB(4)%DQ(x,y,z,1) + Uy(7)/G%J(x,y,z) &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,1)
          F%PMLB(4)%DQ(x,y,z,2) = F%PMLB(4)%DQ(x,y,z,2) + Uy(5)/G%J(x,y,z) &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,2)
          F%PMLB(4)%DQ(x,y,z,3) = F%PMLB(4)%DQ(x,y,z,3) + Uy(9)/G%J(x,y,z) &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,3)
          
          F%PMLB(4)%DQ(x,y,z,4) = F%PMLB(4)%DQ(x,y,z,4) + 0d0 &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,4)
          
          F%PMLB(4)%DQ(x,y,z,5) = F%PMLB(4)%DQ(x,y,z,5) + Uy(2) &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,5)

          F%PMLB(4)%DQ(x,y,z,6) = F%PMLB(4)%DQ(x,y,z,6) + 0d0 &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,6)

          F%PMLB(4)%DQ(x,y,z,7) = F%PMLB(4)%DQ(x,y,z,7) + Uy(1) &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,7)

          F%PMLB(4)%DQ(x,y,z,8) = F%PMLB(4)%DQ(x,y,z,8) + 0d0 &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,8)
          
          F%PMLB(4)%DQ(x,y,z,9) = F%PMLB(4)%DQ(x,y,z,9) + Uy(3) &
               - (dy+cfs)*F%PMLB(4)%Q(x,y,z,9)
          
       end if
    end if
             
    
    ! x-dependent layer
    if(F%PMLB(1)%pml .EQV. .TRUE.) then
       
       npml = F%PMLB(1)%N_pml
       
       if (x .le. npml) then
          
          dx = d0x*(abs(real(x-(1+npml))/real(npml)))**3
          
          DU(1) = DU(1) - dx*G%J(x,y,z)*F%PMLB(1)%Q(x,y,z,1)*rhoJ_inv
          DU(2) = DU(2) - dx*G%J(x,y,z)*F%PMLB(1)%Q(x,y,z,2)*rhoJ_inv
          DU(3) = DU(3) - dx*G%J(x,y,z)*F%PMLB(1)%Q(x,y,z,3)*rhoJ_inv
          
          DU(4) = DU(4) - dx*(lambda2mu*F%PMLB(1)%Q(x,y,z,4) &
               + lam*F%PMLB(1)%Q(x,y,z,5) &
               + lam*F%PMLB(1)%Q(x,y,z,6))
          
          DU(5) = DU(5) - dx*(lam*F%PMLB(1)%Q(x,y,z,4) &
               + lambda2mu*F%PMLB(1)%Q(x,y,z,5) &
               + lam*F%PMLB(1)%Q(x,y,z,6))
          
          DU(6) = DU(6) - dx*(lam*F%PMLB(1)%Q(x,y,z,4) &
               + lam*F%PMLB(1)%Q(x,y,z,5) &
               + lambda2mu*F%PMLB(1)%Q(x,y,z,6))
          
          
          DU(7) = DU(7) - mu*dx*F%PMLB(1)%Q(x,y,z,7)
          DU(8) = DU(8) - mu*dx*F%PMLB(1)%Q(x,y,z,8)
          DU(9) = DU(9) - mu*dx*F%PMLB(1)%Q(x,y,z,9)
          
          F%PMLB(1)%DQ(x,y,z,1) = F%PMLB(1)%DQ(x,y,z,1) + Ux(4)/G%J(x,y,z) &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,1)
          F%PMLB(1)%DQ(x,y,z,2) = F%PMLB(1)%DQ(x,y,z,2) + Ux(7)/G%J(x,y,z) &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,2)
          F%PMLB(1)%DQ(x,y,z,3) = F%PMLB(1)%DQ(x,y,z,3) + Ux(8)/G%J(x,y,z) &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,3)
          
          F%PMLB(1)%DQ(x,y,z,4) = F%PMLB(1)%DQ(x,y,z,4) + Ux(1) &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,4)
          
          F%PMLB(1)%DQ(x,y,z,5) = F%PMLB(1)%DQ(x,y,z,5) + 0d0 &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,5)

          F%PMLB(1)%DQ(x,y,z,6) = F%PMLB(1)%DQ(x,y,z,6) + 0d0 &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,6)

          F%PMLB(1)%DQ(x,y,z,7) = F%PMLB(1)%DQ(x,y,z,7) + Ux(2) &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,7)

          F%PMLB(1)%DQ(x,y,z,8) = F%PMLB(1)%DQ(x,y,z,8) + Ux(3) &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,8)
          
          F%PMLB(1)%DQ(x,y,z,9) = F%PMLB(1)%DQ(x,y,z,9) + 0d0 &
               - (dx+cfs)*F%PMLB(1)%Q(x,y,z,9)
          
       end if
    end if

    
    if(F%PMLB(2)%pml .EQV. .TRUE.) then

       npml = F%PMLB(2)%N_pml
       
       if (x .ge. nx-npml+1) then
          
          dx = d0x*(abs(real(x-(nx-npml))/real(npml)))**3
                   
          DU(1) = DU(1) - dx*G%J(x,y,z)*F%PMLB(2)%Q(x,y,z,1)*rhoJ_inv
          DU(2) = DU(2) - dx*G%J(x,y,z)*F%PMLB(2)%Q(x,y,z,2)*rhoJ_inv
          DU(3) = DU(3) - dx*G%J(x,y,z)*F%PMLB(2)%Q(x,y,z,3)*rhoJ_inv
          
          DU(4) = DU(4) - dx*(lambda2mu*F%PMLB(2)%Q(x,y,z,4) &
               + lam*F%PMLB(2)%Q(x,y,z,5) &
               + lam*F%PMLB(2)%Q(x,y,z,6))
          
          DU(5) = DU(5) - dx*(lam*F%PMLB(2)%Q(x,y,z,4) &
               + lambda2mu*F%PMLB(2)%Q(x,y,z,5) &
               + lam*F%PMLB(2)%Q(x,y,z,6))
          
          DU(6) = DU(6) - dx*(lam*F%PMLB(2)%Q(x,y,z,4) &
               + lam*F%PMLB(2)%Q(x,y,z,5) &
               + lambda2mu*F%PMLB(2)%Q(x,y,z,6))
          
          
          DU(7) = DU(7) - mu*dx*F%PMLB(2)%Q(x,y,z,7)
          DU(8) = DU(8) - mu*dx*F%PMLB(2)%Q(x,y,z,8)
          DU(9) = DU(9) - mu*dx*F%PMLB(2)%Q(x,y,z,9)
          
          F%PMLB(2)%DQ(x,y,z,1) = F%PMLB(2)%DQ(x,y,z,1) + Ux(4)/G%J(x,y,z) &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,1)
          F%PMLB(2)%DQ(x,y,z,2) = F%PMLB(2)%DQ(x,y,z,2) + Ux(7)/G%J(x,y,z) &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,2)
          F%PMLB(2)%DQ(x,y,z,3) = F%PMLB(2)%DQ(x,y,z,3) + Ux(8)/G%J(x,y,z) &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,3)
          
          F%PMLB(2)%DQ(x,y,z,4) = F%PMLB(2)%DQ(x,y,z,4) + Ux(1) &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,4)
          
          F%PMLB(2)%DQ(x,y,z,5) = F%PMLB(2)%DQ(x,y,z,5) + 0d0 &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,5)

          F%PMLB(2)%DQ(x,y,z,6) = F%PMLB(2)%DQ(x,y,z,6) + 0d0 &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,6)

          F%PMLB(2)%DQ(x,y,z,7) = F%PMLB(2)%DQ(x,y,z,7) + Ux(2) &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,7)

          F%PMLB(2)%DQ(x,y,z,8) = F%PMLB(2)%DQ(x,y,z,8) + Ux(3) &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,8)
          
          F%PMLB(2)%DQ(x,y,z,9) = F%PMLB(2)%DQ(x,y,z,9) + 0d0 &
               - (dx+cfs)*F%PMLB(2)%Q(x,y,z,9)
          
       end if
    end if    
    
  end subroutine PML
  

end module RHS_Interior
