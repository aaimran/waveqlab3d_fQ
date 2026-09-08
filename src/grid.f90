module grid

  !> grid module defines spatial grid and coordinate transform on a block

  use common, only : wp
  use mpi3dcomm, only : cartesian3d_t
  use datatypes, only : block_grid_t

  implicit none

contains

   subroutine init_grid_cartesian(G, aqrs, bqrs, nb, debug)

!!! initialize the grid

    use mpi3dcomm, only : allocate_array_body, exchange_all_neighbors

    implicit none

    type(block_grid_t),intent(inout) :: G
    type(cartesian3d_t) :: C
    real(kind = wp), intent(in) :: aqrs(3), bqrs(3)
    integer, intent(in) :: nb
   logical, intent(in), optional :: debug
    integer :: i

    C = G%C

    ! allocate memory for arrays

    call allocate_array_body(G%x, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricx, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricy, C, 3, ghost_nodes=.true.)
    call allocate_array_body(G%metricz, C, 3, ghost_nodes=.true.)
    call allocate_array_body(G%J, C, ghost_nodes=.true.)

    call cartesian_grid_3d(G, aqrs, bqrs)

    do i = 1,3
       call exchange_all_neighbors(C, G%x(:,:,:,i))
    end do

!!! compute metric derivatives (these are partial derivatives of the mapping)
!!! these will be computed using finite difference routines applied to the x array

    call derivatives(G)

    call exchange_all_neighbors(C, G%J(:, :, :))

    do i = 1,3
       call exchange_all_neighbors(C, G%metricx(:, :, :, i))
       call exchange_all_neighbors(C, G%metricy(:, :, :, i))
       call exchange_all_neighbors(C, G%metricz(:, :, :, i))
    end do

  end subroutine init_grid_cartesian

  subroutine init_grid_curve(G, aqrs, bqrs, lc, rc, profile_type, profile_path, &
     topography_type,topography_path, use_mms, use_topography, topo, ny, nz, nb, debug)

!!! initialize the grid

    use mpi3dcomm, only : allocate_array_body, exchange_all_neighbors

    implicit none

    type(block_grid_t),intent(inout) :: G
    type(cartesian3d_t) :: C
    real(kind = wp), intent(in) :: aqrs(3), bqrs(3), lc, rc, topo
    integer, intent(in) :: ny, nz, nb
    logical, intent(in) :: use_mms, use_topography
   logical, intent(in), optional :: debug
    character(64), intent(in) :: profile_type,profile_path
    character(64), intent(in) :: topography_type,topography_path




    integer :: i,j,k
   real(kind=wp) :: minJ, maxJ, maxMetric
    

    C = G%C

    ! allocate memory for arrays

    call allocate_array_body(G%x, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricx, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricy, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricz, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%J, C, ghost_nodes=.true.)

!!! discretize unit cube with uniform Cartesian mesh
!!! (this would be replaced with a mesh generator for the specific block geometry)

    G%bhq = (bqrs(1) - aqrs(1))/real(G%C%nq-1, wp)
    G%bhr = (bqrs(2) - aqrs(2))/real(G%C%nr-1, wp)
    G%bhs = (bqrs(3) - aqrs(3))/real(G%C%ns-1, wp)

    !call Curve_Grid_3D3(G, aqrs(1), bqrs(1), aqrs(2), bqrs(2), &
    !     aqrs(3), bqrs(3), lc, rc, use_mms, profile_type, profile_path)

    call Curve_Grid_3D3(G, aqrs(1), bqrs(1), aqrs(2), bqrs(2), &
         aqrs(3), bqrs(3), lc, rc, use_mms, profile_type, profile_path, &
         use_topography, topography_type,topography_path, topo, ny, nz)

    !call stretched_grid_3d(G,aqrs,bqrs,lc,rc)

    G%hq = 1.0_wp/(real(G%C%nq, wp)-1.0_wp)
    G%hr = 1.0_wp/(real(G%C%nr, wp)-1.0_wp)
    G%hs = 1.0_wp/(real(G%C%ns, wp)-1.0_wp)

    !print*, 'okay'
    do i = 1,3
       call exchange_all_neighbors(C, G%x(:,:,:,i))
    end do

!!! compute metric derivatives (these are partial derivatives of the mapping)
!!! these will be computed using finite difference routines applied to the x array

    call derivatives(G)

    call exchange_all_neighbors(C, G%J(:, :, :))

    do i = 1,3
       call exchange_all_neighbors(C, G%x(:, :, :, i))
       call exchange_all_neighbors(C, G%metricx(:, :, :, i))
       call exchange_all_neighbors(C, G%metricy(:, :, :, i))
       call exchange_all_neighbors(C, G%metricz(:, :, :, i))
    end do

    if (present(debug)) then
       if (debug) then
          minJ = minval(G%J(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps))
          maxJ = maxval(G%J(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps))
          maxMetric = max(maxval(abs(G%metricx(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, :))), &
                          maxval(abs(G%metricy(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, :))), &
                          maxval(abs(G%metricz(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, :))))
          write(*,'(A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') 'debug grid diagnostics: minJ=', minJ, &
               'maxJ=', maxJ, 'max|metric|=', maxMetric
       end if
    end if

  end subroutine init_grid_curve

   subroutine init_grid_from_file(id, G, aqrs, bqrs, nb)

!!! initialize the grid

    use mpi3dbasic, only : pw, ps
    use mpi3dio, only: file_distributed, open_file_distributed, read_file_distributed, close_file_distributed
    use mpi3dcomm, only : allocate_array_body, exchange_all_neighbors

    implicit none

    integer, intent(in) :: id, nb
    type(block_grid_t),intent(inout) :: G
    type(cartesian3d_t) :: C
    type(file_distributed) :: fids(3)
    real(kind = wp), intent(in) :: aqrs(3), bqrs(3)
    integer :: i, k
    character(len=1) :: ext(3), id_string
    character(len=256) :: name

    C = G%C

    ext = ['X', 'Y', 'Z']

    ! allocate memory for arrays

    call allocate_array_body(G%x, C, 3, ghost_nodes=.true.)
    call allocate_array_body(G%metricx, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricy, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%metricz, C,3, ghost_nodes=.true.)
    call allocate_array_body(G%J, C, ghost_nodes=.true.)

    G%bhq = (bqrs(1) - aqrs(1))/real(G%C%nq-1, wp)
    G%bhr = (bqrs(2) - aqrs(2))/real(G%C%nr-1, wp)
    G%bhs = (bqrs(3) - aqrs(3))/real(G%C%ns-1, wp)

    write(id_string, '(i1)') id

    do i = 1, 3
       name = "block_grid_"//id_string//"."//ext(i)
       call open_file_distributed(fids(i), name, "read", C%comm, C%array_w, pw)
       call read_file_distributed(fids(i), G%X(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps,i))
       call close_file_distributed(fids(i))
    end do

    G%hq = 1.0_wp/real(G%C%nq-1, wp)
    G%hr = 1.0_wp/real(G%C%nr-1, wp)
    G%hs = 1.0_wp/real(G%C%ns-1, wp)

    do i = 1,3
       call exchange_all_neighbors(C, G%x(:,:,:,i))
    end do

!!! compute metric derivatives (these are partial derivatives of the mapping)
!!! these will be computed using finite difference routines applied to the x array

    call derivatives(G)

    call exchange_all_neighbors(C, G%J(:, :, :))

    do i = 1,3
       call exchange_all_neighbors(C, G%metricx(:,:,:,i))
       call exchange_all_neighbors(C, G%metricy(:,:,:,i))
       call exchange_all_neighbors(C, G%metricz(:,:,:,i))
    end do

  end subroutine init_grid_from_file

  subroutine cartesian_grid_3d(G, aqrs, bqrs)
    real(kind = wp), intent(in) :: aqrs(3), bqrs(3)
    type(block_grid_t), intent(inout) :: G
    integer :: i,j,k

    !> discretize unit cube with uniform Cartesian mesh
    !> (this would be replaced with a mesh generator for the specific block geometry)

    G%bhq = (bqrs(1) - aqrs(1))/real(G%C%nq-1, wp)
    G%bhr = (bqrs(2) - aqrs(2))/real(G%C%nr-1, wp)
    G%bhs = (bqrs(3) - aqrs(3))/real(G%C%ns-1, wp)


    do k = G%C%ms, G%C%ps
       do j = G%C%mr, G%C%pr
          do i = G%C%mq, G%C%pq
             G%x(i,j,k,1) = aqrs(1) + real(i-1, wp)*G%bhq ! x
             G%x(i,j,k,2) = aqrs(2) + real(j-1, wp)*G%bhr ! y
             G%x(i,j,k,3) = aqrs(3) + real(k-1, wp)*G%bhs ! z
          end do
       end do
    end do

    G%hq = 1.0_wp/real(G%C%nq-1, wp)
    G%hr = 1.0_wp/real(G%C%nr-1, wp)
    G%hs = 1.0_wp/real(G%C%ns-1, wp)

  end subroutine cartesian_grid_3d

  subroutine stretched_grid_3d(G,aqrs,bqrs,lc,rc)
    real(kind = wp), intent(in) :: aqrs(3), bqrs(3)
    real(kind = wp), intent(in) :: lc, rc
    type(block_grid_t), intent(inout) :: G
    integer :: i,j,k, mx,my,mz,nnz
    real(kind = wp) :: q, r, s, gamma, ax,bx,ay,by,az,bz,a0z,b0z,anz,bnz
    real(kind = wp) :: alpha,hfx,hfy,hfz,x_0,y_0,z_0,hnz
    real(kind = wp),dimension(:),allocatable :: x,y,z,x0,xn,z0,zn                                                                                                                                     

    mx = G%C%nq
    my = G%C%nr
    mz = G%C%ns

    nnz = (mz-1)/2+1

    allocate(x(mx),y(my),z(mz),x0(mx),xn(mx),z0(nnz),zn(nnz))


    !> discretize unit cube with uniform Cartesian mesh                                                                                                                                                
    !> (this would be replaced with a mesh generator for the specific block geometry)

    G%bhq = (bqrs(1) - aqrs(1))/real(G%C%nq-1, wp)
    G%bhr = (bqrs(2) - aqrs(2))/real(G%C%nr-1, wp)
    G%bhs = (bqrs(3) - aqrs(3))/real(G%C%ns-1, wp)

    G%hq = 1.0_wp/real(G%C%nq-1, wp)
    G%hr = 1.0_wp/real(G%C%nr-1, wp)
    G%hs = 1.0_wp/real(G%C%ns-1, wp)

    hnz = 1.0_wp/real(nnz-1, wp)

    ax = aqrs(1)
    bx = bqrs(1)
    ay = aqrs(2)
    by = bqrs(2)
    az = aqrs(3)
    bz = bqrs(3)

    x_0 = 10.0_wp
    y_0 = 10.0_wp
    z_0 = 10.0_wp

    alpha = 0.0225_wp


    nnz = (mz-1)/2+1 
    hnz = 1.0_wp/real(nnz-1, wp)

    a0z = az 
    b0z = (az+bz)/2.0_wp

    anz = b0z
    bnz = bz

    hfx = 13.0_wp*G%hq
    hfy = 13.0_wp*G%hr
    hfz = 13.0_wp*hnz

    do k = 1,nnz
       s = real(k-1, wp)*hnz
       if (s <= hnz/hfz*z_0) then
          zn(k) = (s*hfz/hnz)/(z_0 + hfz/alpha*(exp(alpha*(1.0_wp/hnz-z_0/hfz)) - 1.0_wp))
       else
          zn(k) = (z_0 + hfz/alpha*(exp(alpha*(s/hnz-z_0/hfz)) - 1.0_wp))/&
               (z_0 + hfz/alpha*(exp(alpha*(1.0_wp/hnz-z_0/hfz)) - 1.0_wp))
       end if
    end do


    do k = 1,nnz    
       j = nnz -(k-1)
       z0(k) = 1.0_wp - zn(j)
    end do

    do k = 1,nnz
       z0(k) = a0z + (b0z-a0z)*z0(k)
       zn(k) = anz + (bnz-anz)*zn(k);
    end do

    do k = 1,mz
       if (k .le. nnz) z(k) = z0(k)
       if (k > nnz) z(k) = zn(k-(nnz-1))
       !if (k == nnz) print *, z(k)
    end do

    do j = 1,my
       r = real(j-1, wp)*G%hr       
       if (r <= G%hr/hfy*y_0) then
          y(j) = (r*hfy/G%hr)/(y_0 + hfy/alpha*(exp(alpha*(1.0_wp/G%hr-y_0/hfy)) - 1.0_wp))
       else
          y(j) = (y_0 + hfy/alpha*(exp(alpha*(r/G%hr-y_0/hfy)) - 1.0_wp))/&
               (y_0 + hfy/alpha*(exp(alpha*(1.0_wp/G%hr-y_0/hfy)) - 1.0_wp))
       end if
    end do

    do j = 1,my
       y(j) = ay + (by-ay)*y(j);
    end do


    do i = 1,mx
       q = real(i-1, wp)*G%hq
       if (q <= G%hq/hfx*x_0) then
          xn(i) = (q*hfx/G%hq)/(x_0 + hfx/alpha*(exp(alpha*(1.0_wp/G%hq-x_0/hfx)) - 1.0_wp))
       else

          xn(i) = (x_0 + hfx/alpha*(exp(alpha*(q/G%hq-x_0/hfx)) - 1.0_wp))/&
               (x_0 + hfx/alpha*(exp(alpha*(1.0_wp/G%hq-x_0/hfx)) - 1.0_wp))
       end if
    end do

    do i = 1,mx
       j = mx -(i-1)
       x0(i) = 1.0_wp - xn(j)
       !print *, i, x0(i), xn(i)
    end do

    do i = 1,mx
       x0(i) = ax + (bx-ax)*x0(i)
       xn(i) = ax + (bx-ax)*xn(i);
    end do



    do k = G%C%ms, G%C%ps
       do j = G%C%mr, G%C%pr
          do i = G%C%mq, G%C%pq

             ! x
             if (rc == 1 .and. lc == 0) then
                G%x(i,j,k,1) = xn(i)
                !print *, i, x0(i)
             elseif (rc == 0 .and. lc == 1) then 
                G%x(i,j,k,1) = x0(i)
                !print *, i, x0(i)
             else
                ! invalid case, stop with error message
                stop 'invalid case for block interface surface (try setting lc,rc in btp)'
             end if

             G%x(i,j,k,2) = y(j) ! y    
             G%x(i,j,k,3) = z(k) ! z                                                                                                                                                                        
             !print *, G%x(i,j,k,1), G%x(i,j,k,2), G%x(i,j,k,3)         
          end do
       end do
    end do

    G%hq = 1.0_wp/real(G%C%nq-1, wp)
    G%hr = 1.0_wp/real(G%C%nr-1, wp)
    G%hs = 1.0_wp/real(G%C%ns-1, wp)

    deallocate(x0, xn, x, y, z, z0, zn)                                                                                                                                                                                   

  end subroutine stretched_grid_3d

  subroutine derivatives(G)
    ! Compute the metric coefficients q_x,q_y,q_z, r_x,r_y,r_z, s_x,s_y,s_z
    ! and the Jacobian J

    use metrics, only : D6x, D6y, D6z, Dfx6, Dfy6, Dfz6
    use mpi3dcomm

    implicit none
    type(block_grid_t), intent(inout) :: G

    real(kind = wp),dimension(3) :: X_q, X_r, X_s                                 ! work arrays
    integer :: q, r, s                                                 ! indices in the tranformed cordinate
    real(kind = wp) :: hq, hr, hs                                                 ! spatial steps discertizing the unit cube

    hq = 1.0_wp/real(G%C%nq-1, wp)
    hr = 1.0_wp/real(G%C%nr-1, wp)
    hs = 1.0_wp/real(G%C%ns-1, wp)

    ! compute the spatial steps

    do s = G%C%ms,G%C%ps
       do r = G%C%mr,G%C%pr
          do q = G%C%mq,G%C%pq

             ! compute x_q, y_q, z_q
             X_q = D6x(hq, G%x, q, r, s, G%C%mbq, G%C%pbq, G%C%mbr, G%C%pbr, G%C%mbs, G%C%pbs, G%C%nq)

             ! compute x_r, y_r, z_r
             X_r = D6y(hr, G%x, q, r, s, G%C%mbq, G%C%pbq, G%C%mbr, G%C%pbr, G%C%mbs, G%C%pbs, G%C%nr)

             ! compute x_s, y_s, z_s
             X_s = D6z(hs, G%x, q, r, s, G%C%mbq, G%C%pbq, G%C%mbr, G%C%pbr, G%C%mbs, G%C%pbs, G%C%ns)

             ! compute x_q, y_q, z_q
             !X_q = Dfx6(hq, G%x, q, r, s, G%C%mbq, G%C%pbq, G%C%mbr, G%C%pbr, G%C%mbs, G%C%pbs, G%C%nq, 1, 3)

             ! compute x_r, y_r, z_r
             !X_r = Dfy6(hr, G%x, q, r, s, G%C%mbq, G%C%pbq, G%C%mbr, G%C%pbr, G%C%mbs, G%C%pbs, G%C%nr, 1, 3)

             ! compute x_s, y_s, z_s
             !X_s = Dfz6(hs, G%x, q, r, s, G%C%mbq, G%C%pbq, G%C%mbr, G%C%pbr, G%C%mbs, G%C%pbs, G%C%ns, 1, 3)

             !============================================
            !print *, X_s

             ! compute Jacobian
             G%J(q, r, s) = &
                  X_q(1)*(X_r(2)*X_s(3) - X_r(3)*X_s(2)) &
                  & - X_q(2)*(X_r(1)*X_s(3) - X_r(3)*X_s(1)) &
                  & + X_q(3)*(X_r(1)*X_s(2) - X_r(2)*X_s(1))
             !============================================

             ! compute metric coeffs: q_x, r_x, s_x
             G%metricx(q, r, s, 1) = (1.0_wp/G%J(q, r, s))*(X_r(2)*X_s(3) - X_r(3)*X_s(2))
             G%metricx(q, r, s, 2) = (1.0_wp/G%J(q, r, s))*(X_q(3)*X_s(2) - X_q(2)*X_s(3))
             G%metricx(q, r, s, 3) = (1.0_wp/G%J(q, r, s))*(X_q(2)*X_r(3) - X_q(3)*X_r(2))

             ! compute metric coeffs: q_y, r_y, s_y
             G%metricy(q, r, s, 1) = (1.0_wp/G%J(q, r, s))*(X_r(3)*X_s(1) - X_r(1)*X_s(3))
             G%metricy(q, r, s, 2) = (1.0_wp/G%J(q, r, s))*(X_q(1)*X_s(3) - X_q(3)*X_s(1))
             G%metricy(q, r, s, 3) = (1.0_wp/G%J(q, r, s))*(X_q(3)*X_r(1) - X_q(1)*X_r(3))

             ! compute metric coeffs: q_z, r_z, s_z
             G%metricz(q, r, s, 1) = (1.0_wp/G%J(q, r, s))*(X_r(1)*X_s(2) - X_r(2)*X_s(1))
             G%metricz(q, r, s, 2) = (1.0_wp/G%J(q, r, s))*(X_q(2)*X_s(1) - X_q(1)*X_s(2))
             G%metricz(q, r, s, 3) = (1.0_wp/G%J(q, r, s))*(X_q(1)*X_r(2) - X_q(2)*X_r(1))

          end do
       end do
    end do

!STOP

  end subroutine derivatives

  subroutine Curve_Grid_3D3(G,ax,bx,ay,by,az,bz,lc,rc,use_mms,profile_type,profile_path,&
       use_topography,topography_type,topography_path, topo, ny,nz)
    !> Generate a 3D curvilinear mesh using the transfinite interpolation grid generator
    !> implemented in Gen_Curve_Mesh below. Any user defined smooth surfaces can be used

    type(block_grid_t),intent(inout) :: G
    integer :: mx,my,mz,px,py,pz,nq,nr,ns            ! number of grid points
    real(kind = wp), intent(in) :: ax,bx, ay,by, az,bz          ! dimensions of a 3D computational elastic block
    real(kind = wp), intent(in) :: lc, rc, topo                       ! lc, rc  add boundary curves to the x-boundaries
    logical,intent(in) :: use_mms, use_topography
    integer,intent(in) :: ny,nz
    character(64), intent(in) :: profile_type,profile_path, topography_type, topography_path

    !real(kind = wp), dimension(:,:,:), allocatable:: X, Y, Z      ! x y z components of the grid

    ! Six boundary surfaces for a 3D physical domain
    real(kind = wp),dimension(:,:  ),allocatable :: Xleft,Xright,Xtop,Xbottom,Xfront,Xback
    real(kind = wp),dimension(:,:  ),allocatable :: Yleft,Yright,Ytop,Ybottom,Yfront,Yback
    real(kind = wp),dimension(:,:  ),allocatable :: Zleft,Zright,Ztop,Zbottom,Zfront,Zback

    ! Seventy two boundary edges needed to generate 2D boundary curves
    real(kind = wp),dimension(:),allocatable :: X0xL,X0xR,Y0xL,Y0xR,Z0xL,Z0xR
    real(kind = wp),dimension(:),allocatable :: X0xB,X0xT,Y0xB,Y0xT,Z0xB,Z0xT
    real(kind = wp),dimension(:),allocatable :: X0yL,X0yR,Y0yL,Y0yR,Z0yL,Z0yR
    real(kind = wp),dimension(:),allocatable :: X0yB,X0yT,Y0yB,Y0yT,Z0yB,Z0yT
    real(kind = wp),dimension(:),allocatable :: X0zL,X0zR,Y0zL,Y0zR,Z0zL,Z0zR
    real(kind = wp),dimension(:),allocatable :: X0zB,X0zT,Y0zB,Y0zT,Z0zB,Z0zT
    real(kind = wp),dimension(:),allocatable :: XnxL,XnxR,YnxL,YnxR,ZnxL,ZnxR
    real(kind = wp),dimension(:),allocatable :: XnxB,XnxT,YnxB,YnxT,ZnxB,ZnxT
    real(kind = wp),dimension(:),allocatable :: XnyL,XnyR,YnyL,YnyR,ZnyL,ZnyR
    real(kind = wp),dimension(:),allocatable :: XnyB,XnyT,YnyB,YnyT,ZnyB,ZnyT
    real(kind = wp),dimension(:),allocatable :: XnzL,XnzR,YnzL,YnzR,ZnzL,ZnzR
    real(kind = wp),dimension(:),allocatable :: XnzB,XnzT,YnzB,YnzT,ZnzB,ZnzT

    real(kind = wp),dimension(:),allocatable :: X0avg, Xmavg

    real(kind = wp),dimension(:,:  ),allocatable :: Fault_Geometry
    real(kind = wp),dimension(:,:  ),allocatable :: X1,Y1,Z1
    real(kind = wp),dimension(:,:  ),allocatable :: XX1,YY1,ZZ1

    real(kind = wp),dimension(:,:  ),allocatable :: Topography
    real(kind = wp),dimension(:,:  ),allocatable :: XT,YT,ZT
    real(kind = wp),dimension(:,:  ),allocatable :: X0T,Y0T,Z0T

    real(kind = wp) :: hx, hy, hz
    integer :: i, j, k
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp) :: r1, r2, f, ang, tpv36_shift
    real(kind = wp) :: r0y, rny, r0z, rnz, q, r, s

    real(kind = wp) :: y, z, dz, dy, gg, delta, stretch
    integer :: jj, kk, nmy, npy, nmz, npz, stenc, ndp


    mx = G%C%nq
    my = G%C%nr
    mz = G%C%ns

    nq = G%C%nq
    nr = G%C%nr
    ns = G%C%ns


    allocate(Xfront(mx,my),Xback(mx,my),Yfront(mx,my),&
         Yback(mx,my),Zfront(mx,my),Zback(mx,my))

    allocate(Xtop(mx,mz),Xbottom(mx,mz),Ytop(mx,mz),&
         Ybottom(mx,mz),Ztop(mx,mz),Zbottom(mx,mz))

    allocate(Xleft(my,mz),Xright(my,mz),Yleft(my,mz),&
         Yright(my,mz),Zleft(my,mz),Zright(my,mz))
    !allocate(Fault_Geometry(my*mz, 3))
    allocate(Fault_Geometry(ny*nz, 3))
    allocate(Topography(mx*mz, 3))
    !allocate(X1(my,mz),Y1(my,mz),Z1(my,mz))
    allocate(XT(mz,mx),YT(mz,mx),ZT(mz,mx))
    !allocate(X0T(mx,mx),Y0T(mz,mx),Z0T(mz,mx))

    allocate(X1(ny,nz),Y1(ny,nz),Z1(ny,nz))
    allocate(XX1(my,mz),YY1(my,mz),ZZ1(my,mz))

    allocate(Y0xL(mz),Y0xR(mz),Z0xL(mz),Z0xR(mz))
    allocate(Y0xB(my),Y0xT(my),Z0xB(my),Z0xT(my))

    allocate(X0yL(mz),X0yR(mz),Z0yL(mz),Z0yR(mz))
    allocate(X0yB(mx),X0yT(mx),Z0yB(mx),Z0yT(mx))

    allocate(X0zL(my),X0zR(my),Y0zL(my),Y0zR(my))
    allocate(X0zB(mx),X0zT(mx),Y0zB(mx),Y0zT(mx))

    allocate(YnxL(mz),YnxR(mz),ZnxL(mz),ZnxR(mz))
    allocate(YnxB(my),YnxT(my),ZnxB(my),ZnxT(my))

    allocate(XnyL(mz),XnyR(mz),ZnyL(mz),ZnyR(mz))
    allocate(XnyB(mx),XnyT(mx),ZnyB(mx),ZnyT(mx))

    allocate(XnzL(my),XnzR(my),YnzL(my),YnzR(my))
    allocate(XnzB(mx),XnzT(mx),YnzB(mx),YnzT(mx))

    allocate(X0avg(mz), Xmavg(mz))




    hx = (bx-ax)/real(mx-1, wp)
    hy = (by-ay)/real(my-1, wp)
    hz = (bz-az)/real(mz-1, wp)

    r0y = 1.0_wp
    r0z = 1.0_wp
    rny = 1.0_wp
    rnz = 1.0_wp



    if (use_mms .NEQV.  .TRUE.) then

       ! x-face

       ! default

       if (use_topography .neqv.  .true.) then


          hy = (by-ay)/real(my-1, wp)
          hz = (bz-az)/real(mz-1, wp)

          do k = 1,mz
             do j = 1,my

                Yleft(j, k) = ay+real(j-1, wp)*hy;
                Zleft(j, k) = az+real(k-1, wp)*hz;
                Yright(j, k) = ay+real(j-1, wp)*hy;
                Zright(j, k) = az+real(k-1, wp)*hz;
                Xleft(j, k) = ax
                Xright(j, k) = bx

                ! Geometry for TPV28
                f = 0.0_wp
                r1 = sqrt((Zleft(j, k)+10.5_wp)**2 + (Yleft(j, k)-7.5_wp)**2)
                r2 = sqrt((Zleft(j, k)-10.5_wp)**2 + (Yleft(j, k)-7.5_wp)**2)

                if (r1 .le. 3.0_wp) f = -0.3_wp*(1.0_wp + cos(pi*r1/3.0_wp))
                if (r2 .le. 3.0_wp) f = -0.3_wp*(1.0_wp + cos(pi*r2/3.0_wp))

                Xleft(j, k) =  Xleft(j, k) + 0d0*rc*f
                Xright(j, k) =  Xright(j, k) + 0d0*lc*f

                !==============================================

             end do
          end do


          select case(profile_type)

          case('analytical_tpv28')

             hy = (by-ay)/real(my-1, wp)
             hz = (bz-az)/real(mz-1, wp)

             do k = 1,mz
                do j = 1,my

                   Yleft(j, k) = ay+real(j-1, wp)*hy;
                   Zleft(j, k) = az+real(k-1, wp)*hz;
                   Yright(j, k) = ay+real(j-1, wp)*hy;
                   Zright(j, k) = az+real(k-1, wp)*hz;
                   Xleft(j, k) = ax
                   Xright(j, k) = bx

                   ! Geometry for TPV28
                   f = 0.0_wp
                   r1 = sqrt((Zleft(j, k)+10.5_wp)**2 + (Yleft(j, k)-7.5_wp)**2)
                   r2 = sqrt((Zleft(j, k)-10.5_wp)**2 + (Yleft(j, k)-7.5_wp)**2)

                   if (r1 .le. 3.0_wp) f = -0.3_wp*(1.0_wp + cos(pi*r1/3.0_wp))
                   if (r2 .le. 3.0_wp) f = -0.3_wp*(1.0_wp + cos(pi*r2/3.0_wp))

                   Xleft(j, k) =  Xleft(j, k) + rc*f
                   Xright(j, k) =  Xright(j, k) + lc*f

                   !==============================================

                end do
             end do

          case('analytical_tpv10')

             hy = (by-ay)/real(my-1, wp)
             hz = (bz-az)/real(mz-1, wp)

             ang = 60.0_wp/180.0_wp*pi

             do k = 1,mz
                do j = 1,my

                   Yleft(j, k) = ay+real(j-1, wp)*hy;
                   Zleft(j, k) = az+real(k-1, wp)*hz;
                   Yright(j, k) = ay+real(j-1, wp)*hy;
                   Zright(j, k) = az+real(k-1, wp)*hz;
                   Xleft(j, k) = ax
                   Xright(j, k) = bx



                   ! Geometry for TPV10
                   f = 0.0_wp

                   if (rc > 0.0_wp) then
                      if (Yleft(j,k) .le. 15.0_wp) then

                         f = (1.0_wp/tan(ang))*(Yleft(j,k) - 0.5_wp*(ay + by))

                      else

                         f = ((1.0_wp/tan(ang))*(15.0_wp - 0.5_wp*(ay + by)) + &
                              0.5_wp*atan(4.0_wp*(Yleft(j,k)-15.0_wp))*exp(-5.0_wp*(Yleft(j,k)-15.0_wp)))

                      end if

                   end if

                   if (lc > 0.0_wp) then
                      if (Yright(j,k) .le. 15.0_wp) then

                         f = (1.0_wp/tan(ang))*(Yright(j,k) - 0.5_wp*(ay + by))

                      else

                         f = ((1.0_wp/tan(ang))*(15.0_wp - 0.5_wp*(ay + by)) + &
                              0.5_wp*atan(4.0_wp*(Yright(j,k)-15.0_wp))*exp(-5.0_wp*(Yright(j,k)-15.0_wp)))

                      end if
                   end if


                   Xleft(j, k) =  Xleft(j, k) + rc*f
                   Xright(j, k) =  Xright(j, k) + lc*f
                   !==============================================

                end do
             end do

          case('analytical_tpv36', 'analytical_tpv36_unshifted', &
               'analytical_tpv36_ver_a', 'analytical_tpv36_ver_b', &
               'analytical_tpv36_ver_c', 'analytical_tpv36_ver_d', &
               'analytical_tpv36_a', 'analytical_tpv36_b', &
               'analytical_tpv36_c', 'analytical_tpv36_d', &
               'analytical_tpv36_e')

             hy = (by-ay)/real(my-1, wp)
             hz = (bz-az)/real(mz-1, wp)

             ang = 15.0_wp/180.0_wp*pi
             tpv36_shift = 15.0_wp
             if (trim(profile_type) /= 'analytical_tpv36') &
                  tpv36_shift = 0.0_wp

             do k = 1,mz
                do j = 1,my

                   Yleft(j, k) = ay+real(j-1, wp)*hy;
                   Zleft(j, k) = az+real(k-1, wp)*hz;
                   Yright(j, k) = ay+real(j-1, wp)*hy;
                   Zright(j, k) = az+real(k-1, wp)*hz;
                   Xleft(j, k) = ax
                   Xright(j, k) = bx



                   ! Geometry for TPV36                                                                                                                                                                                                                                             
                   f = 0.0_wp

                   if (rc > 0.0_wp) then
                      if (trim(profile_type) == 'analytical_tpv36_ver_a' .or. &
                          trim(profile_type) == 'analytical_tpv36_ver_b' .or. &
                          trim(profile_type) == 'analytical_tpv36_ver_c' .or. &
                          trim(profile_type) == 'analytical_tpv36_a' .or. &
                          trim(profile_type) == 'analytical_tpv36_b' .or. &
                          trim(profile_type) == 'analytical_tpv36_c' .or. &
                          trim(profile_type) == 'analytical_tpv36_e') then

                         f = tpv36_ver_interface(Yleft(j,k))

                      else if (trim(profile_type) == 'analytical_tpv36_ver_d' .or. &
                               trim(profile_type) == 'analytical_tpv36_d') then

                         f = Yleft(j,k)/tan(ang)

                      else if (Yleft(j,k) .le. 8.0_wp) then

                         f = 1.0_wp/tan(ang)*Yleft(j,k) - tpv36_shift
                         
                      else

                         f = 1.0_wp/tan(ang)*8.0_wp - tpv36_shift + &
                              0.25_wp*atan(4.0_wp*(Yleft(j,k)-8.0_wp))*exp(-0.0_wp*(Yleft(j,k)-8.0_wp))
                              !1.0_wp*atan(1.0_wp*(Yleft(j,k)-8.5_wp))*exp(-2.0_wp*(Yleft(j,k)-8.5_wp))

                      end if

                   end if

                   if (lc > 0.0_wp) then
                      if (trim(profile_type) == 'analytical_tpv36_ver_a' .or. &
                          trim(profile_type) == 'analytical_tpv36_ver_b' .or. &
                          trim(profile_type) == 'analytical_tpv36_ver_c' .or. &
                          trim(profile_type) == 'analytical_tpv36_a' .or. &
                          trim(profile_type) == 'analytical_tpv36_b' .or. &
                          trim(profile_type) == 'analytical_tpv36_c' .or. &
                          trim(profile_type) == 'analytical_tpv36_e') then

                         f = tpv36_ver_interface(Yright(j,k))

                      else if (trim(profile_type) == 'analytical_tpv36_ver_d' .or. &
                               trim(profile_type) == 'analytical_tpv36_d') then

                         f = Yright(j,k)/tan(ang)

                      else if (Yright(j,k) .le. 8.0_wp) then

                         f = 1.0_wp/tan(ang)*Yright(j,k) - tpv36_shift

                      else

                         f = 1.0_wp/tan(ang)*8.0_wp - tpv36_shift + &
                               0.25_wp*atan(4.0_wp*(Yright(j,k)-8.0_wp))*exp(-0.0_wp*(Yright(j,k)-8.0_wp))
                              !1.0_wp*atan(1.0_wp*(Yright(j,k)-8.5_wp))*exp(-2.0_wp*(Yright(j,k)-8.5_wp))

                      end if
                   end if

                   
                   Xleft(j, k) =  Xleft(j, k) + rc*f
                   Xright(j, k) =  Xright(j, k) + lc*f
                   !==============================================                                                                                                                                                                                                                  

                end do
             end do

   
          case('analytical_test_problem')

             hy = (by-ay)/real(my-1, wp)
             hz = (bz-az)/real(mz-1, wp)

             do j = 1,my
                Y0xT(j) = ay+real(j-1, wp)*hy
                Y0xB(j) = ay+real(j-1, wp)*hy
                Z0xT(j) = bz
                Z0xB(j) = az
             end do

             do k = 1,mz
                Z0xL(k) = az+real(k-1, wp)*hz
                Z0xR(k) = az+real(k-1, wp)*hz
                Y0xL(k) = ay
                Y0xR(k) = by
             end do


             do j = 1,my
                YnxT(j) = ay+real(j-1, wp)*hy
                YnxB(j) = ay+real(j-1, wp)*hy
                ZnxT(j) = bz
                ZnxB(j) = az
             end do

             do k = 1,mz
                ZnxL(k) = az+real(k-1, wp)*hz
                ZnxR(k) = az+real(k-1, wp)*hz
                YnxL(k) = ay
                YnxR(k) = by
             end do

             call Gen_Surf_Mesh(my,mz,Y0xL,Y0xR,Y0xB,Y0xT,Z0xL,Z0xR,Z0xB,Z0xT,Yleft,Zleft)
             call Gen_Surf_Mesh(my,mz,YnxL,YnxR,YnxB,YnxT,ZnxL,ZnxR,ZnxB,ZnxT,Yright,Zright)

             do k = 1,mz
                do j = 1,my

                   Xleft(j, k) = ax&
                        + 0.1_wp*rc*exp(-((1.0_wp*Yleft(j, k)-0.5_wp*(ay+by))**2 &
                        + 1.0_wp*(Zleft(j, k)-0.5_wp*(az+bz))**2)/0.025_wp)&
                        + 0.1_wp*rc*sin(2.0_wp/(by-ay)*pi*Yleft(j, k))&
                        *sin(2.0_wp/(bz-az)*pi*Zleft(j, k))

                   Xright(j, k) = bx&
                        + 0.1_wp*lc*exp(-(1.0_wp*(Yright(j, k)-0.5_wp*(ay+by))**2 &
                        + 1.0_wp*(Zright(j, k)-0.5_wp*(az+bz))**2)/0.025_wp)&
                        + 0.1_wp*lc*sin(2.0_wp/(by-ay)*pi*Yright(j, k))&
                        *sin(2.0_wp/(bz-az)*pi*Zright(j, k))


                   !==============================================
                end do
             end do

          case('read_from_memomry_fractal')

             !read fault geometry from memory
             call read_2darray(Fault_Geometry, my*mz,3,profile_path)
             call reshape_2darray(X1, Fault_Geometry(:,1), my,mz)

             call reshape_2darray(Y1, Fault_Geometry(:,2), my,mz)
             call reshape_2darray(Z1, Fault_Geometry(:,3), my,mz)

             hy = (by-ay)/real(my-1, wp)
             hz = (bz-az)/real(mz-1, wp)


             !                                                   left-face: i = 1
             !
             !                                                       Top
             !                                          !-----------------------------!
             !                                          !                             ! 
             !                                     Left !z                            ! Right
             !                                          !                             !
             !                                          !             y               !
             !                                          !-----------------------------!
             !                                                        Bottom

             !L
             do j = 1, my
                Y0xT(j) = ay+real(j-1, wp)*hy
                Y0xB(j) = ay+real(j-1, wp)*hy
                Z0xT(j) = bz
                Z0xB(j) = az
             end do

             do k = 1, mz
                Z0xL(k) = az+real(k-1, wp)*hz
                Z0xR(k) = az+real(k-1, wp)*hz
                Y0xL(k) = ay
                Y0xR(k) = by
             end do

             call Gen_Surf_Mesh(my,mz,Y0xL,Y0xR,Y0xB,Y0xT,Z0xL,Z0xR,Z0xB,Z0xT,Yleft,Zleft)

             do k = 1,mz
                do j = 1,my

                   Xleft(j, k) = ax
                   if (rc /= 0.0)  Xleft(j, k) = Xleft(j, k)+rc*X1(j,k)
                end do
             end do








             !                                                   right-face: i = mx
             !
             !                                                       Top
             !                                          !-----------------------------!
             !                                          !                             ! 
             !                                     Left !z                            ! Right
             !                                          !                             !
             !                                          !             y               !
             !                                          !-----------------------------!
             !                                                        Bottom

             do j = 1,my
                YnxT(j) = ay+real(j-1, wp)*hy
                YnxB(j) = ay+real(j-1, wp)*hy
                ZnxT(j) = bz
                ZnxB(j) = az
             end do

             do k = 1,mz
                ZnxL(k) = az+real(k-1, wp)*hz
                ZnxR(k) = az+real(k-1, wp)*hz
                YnxL(k) = ay
                YnxR(k) = by
             end do

            
             call Gen_Surf_Mesh(my,mz,YnxL,YnxR,YnxB,YnxT,ZnxL,ZnxR,ZnxB,ZnxT,Yright,Zright)

             do k = 1,mz
                do j = 1,my

                   Xright(j, k) = bx

                   if (lc /=0.0)  Xright(j, k) = Xright(j, k)+lc*X1(j,k)

                end do
             end do


          end select

          
          !                                                        z-face
          !
          !
          !                                                   front-face: k = 1
          !
          !                                                       Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !y                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom   


          do j = 1,my

             Y0zL(j) = Yleft(j,1)  
             Y0zR(j) = Yright(j,1)

             X0zL(j) = Xleft(j,1)
             X0zR(j) = Xright(j,1)


          end do


          do i = 1,mx

             hx = (Xright(1,1)-Xleft(1,1))/real(mx-1, wp)
             X0zB(i) = Xleft(1,1)+real(i-1, wp)*hx

             hx = (Xright(my,1)-Xleft(my,1))/real(mx-1, wp)
             X0zT(i) = Xleft(my,1)+real(i-1, wp)*hx


             hy = (Yright(1,1)-Yleft(1,1))/real(mx-1, wp)
             Y0zB(i) = Yleft(1,1)+ real(i-1, wp)*hy

             hy = (Yright(my,1)-Yleft(my,1))/real(mx-1, wp)
             Y0zT(i) = Yleft(my,1)+ real(i-1, wp)*hy
          end do


          call Gen_Surf_Mesh(mx,my,X0zL,X0zR,X0zB,X0zT,Y0zL,Y0zR,Y0zB,Y0zT,Xfront,Yfront)

          do j = 1,my
             do i = 1,mx

                Zfront(i, j) = az

             end do
          end do




          !                                                   back-face: k = mz
          !
          !                                                       Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !y                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom         

          do i = 1,mx             

             hx = (Xright(1,mz)-Xleft(1,mz))/real(mx-1, wp)
             XnzB(i) = Xleft(1,mz)+real(i-1, wp)*hx

             hx = (Xright(my,mz)-Xleft(my,mz))/real(mx-1, wp)
             XnzT(i) = Xleft(my,mz)+real(i-1, wp)*hx




             hy = (Yright(1,mz)-Yleft(1,mz))/real(mx-1, wp)
             YnzB(i) = Yleft(1,mz)+real(i-1, wp)*hy

             hy = (Yright(my,mz)-Yleft(my,mz))/real(mx-1, wp)
             YnzT(i) = Yleft(my,mz)+real(i-1, wp)*hy

          end do




          do j = 1,my


             YnzL(j) = Yleft(j,mz) 
             YnzR(j) = Yright(j,mz)

             XnzL(j) = Xleft(j,mz)
             XnzR(j) = Xright(j,mz) 

          end do


          call Gen_Surf_Mesh(mx,my,XnzL,XnzR,XnzB,XnzT,YnzL,YnzR,YnzB,YnzT,Xback,Yback)


          do j = 1,my
             do i = 1,mx

                Zback(i, j) = bz

             end do
          end do


          ! y-face:
          ! boundary curves at the top and bottom surfaces

          !                                                   top-face: j = 1
          !
          !                                                       Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !z                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom

          do i = 1,mx

             Z0yT(i) = Zback(i,1)
             Z0yB(i) = Zfront(i,1)

          end do

          do i = 1,mx


             hx = (Xright(1,1)-Xleft(1,1))/real(mx-1, wp)
             X0yB(i) = Xleft(1,1)+real(i-1, wp)*hx

             hx = (Xright(1,mz)-Xleft(1,mz))/real(mx-1, wp)
             X0yT(i) = Xleft(1,mz)+real(i-1, wp)*hx

          end do

          do k = 1,mz

             X0yL(k) = Xleft(1,k)
             X0yR(k) = Xright(1,k)

             Z0yL(k) = Zleft(1,k)
             Z0yR(k) = Zright(1,k)


          end do


          call Gen_Surf_Mesh(mx,mz,X0yL,X0yR,X0yB,X0yT,Z0yL,Z0yR,Z0yB,Z0yT,Xtop,Ztop)


          do k = 1, mz
             do i = 1,mx

                Ytop(i, k) = ay
             end do
          end do





          !                                                   bottom-face: j = my
          !
          !                                                       Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !z                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom

          do i = 1,mx

             ZnyT(i) = Zback(i,my)
             ZnyB(i) = Zfront(i,my)

          end do

          do i = 1,mx


             hx = (Xright(my,1)-Xleft(my,1))/real(mx-1, wp)
             XnyB(i) = Xleft(my,1)+real(i-1, wp)*hx

             hx = (Xright(my,mz)-Xleft(my,mz))/real(mx-1, wp)
             XnyT(i) = Xleft(my,mz)+real(i-1, wp)*hx

          end do

          do k = 1,mz

             XnyL(k) = Xleft(my,k)
             XnyR(k) = Xright(my,k)

             ZnyL(k) = Zleft(my,k)
             ZnyR(k) = Zright(my,k)


          end do



          call Gen_Surf_Mesh(mx,mz,XnyL,XnyR,XnyB,XnyT,ZnyL,ZnyR,ZnyB,ZnyT,Xbottom,Zbottom)


          do k = 1, mz
             do i = 1,mx

                Ybottom(i, k) = by
             end do
          end do


!!$             do j = 1,my
!!$                Y0xT(j) = ay+real(j-1, wp)*hy
!!$                Y0xB(j) = ay+real(j-1, wp)*hy
!!$                Z0xT(j) = bz
!!$                Z0xB(j) = az
!!$             end do
!!$
!!$             do k = 1,mz
!!$                Z0xL(k) = az+real(k-1, wp)*hz
!!$                Z0xR(k) = az+real(k-1, wp)*hz
!!$                Y0xL(k) = by
!!$                Y0xR(k) = ay
!!$             end do
!!$
!!$
!!$             do j = 1,my
!!$                YnxT(j) = ay+real(j-1, wp)*hy
!!$                YnxB(j) = ay+real(j-1, wp)*hy
!!$                ZnxT(j) = bz
!!$                ZnxB(j) = az
!!$             end do
!!$
!!$             do k = 1,mz
!!$                ZnxL(k) = az+real(k-1, wp)*hz
!!$                ZnxR(k) = az+real(k-1, wp)*hz
!!$                YnxL(k) = by
!!$                YnxR(k) = ay
!!$             end do
!!$
!!$             call Gen_Surf_Mesh(my,mz,Y0xL,Y0xR,Y0xB,Y0xT,Z0xL,Z0xR,Z0xB,Z0xT,Yleft,Zleft)
!!$             call Gen_Surf_Mesh(my,mz,YnxL,YnxR,YnxB,YnxT,ZnxL,ZnxR,ZnxB,ZnxT,Yright,Zright)
!!$
!!$
!!$             do k = 1,mz
!!$                do j = 1,my
!!$
!!$
!!$                   Xleft(j, k) = ax
!!$                   Xright(j, k) = bx
!!$
!!$                   if (rc /=0.0_wp) then
!!$                      Xleft(j, k) = Xleft(j, k)+rc*X1(j,k)
!!$                      Yleft(j, k) = ay+real(j-1, wp)*hy !Y1(j,k)
!!$                      Zleft(j, k) = az+real(k-1, wp)*hz !Z1(j,k)
!!$                   end if
!!$
!!$                   if (lc /=0.0_wp) then
!!$                      Xright(j, k) = Xright(j, k)+lc*X1(j,k)
!!$                      Yright(j, k) = ay+real(j-1, wp)*hy !Y1(j,k)
!!$                      Zright(j, k) = az+real(k-1, wp)*hz !Z1(j,k)
!!$                   end if
!!$
!!$                   !==============================================
!!$                end do
!!$             end do

!!$       end select
!!$
!!$          ! z-face
!!$          do i = 1,mx
!!$
!!$             hx = (Xright(1,1)-Xleft(1,1))/real(mx-1, wp)
!!$             X0zT(i) = Xleft(1,1)+real(i-1, wp)*hx
!!$
!!$             hx = (Xright(my,1)-Xleft(my,1))/real(mx-1, wp)
!!$             X0zB(i) = Xleft(my,1)+real(i-1, wp)*hx
!!$
!!$             hx = (Xright(1,mz)-Xleft(1,mz))/real(mx-1, wp)
!!$             XnzT(i) = Xleft(1,mz)+real(i-1, wp)*hx
!!$
!!$             hx = (Xright(my,mz)-Xleft(my,mz))/real(mx-1, wp)
!!$             XnzB(i) = Xleft(my,mz)+real(i-1, wp)*hx
!!$
!!$
!!$             hy = (Yright(1,1)-Yleft(1,1))/real(mx-1, wp)
!!$             Y0zT(i) = Yleft(1,1)+real(i-1, wp)*hy
!!$
!!$             hy = (Yright(my,1)-Yleft(my,1))/real(mx-1, wp)
!!$             Y0zB(i) = Yleft(my,1)+real(i-1, wp)*hy
!!$
!!$             hy = (Yright(1,mz)-Yleft(1,mz))/real(mx-1, wp)
!!$             YnzT(i) = Yleft(1,mz)+real(i-1, wp)*hy
!!$
!!$             hy = (Yright(my,mz)-Yleft(my,mz))/real(mx-1, wp)
!!$             YnzB(i) = Yleft(my,mz)+real(i-1, wp)*hy
!!$
!!$          end do
!!$
!!$          !stop
!!$
!!$          do j = 1,my
!!$
!!$             Y0zL(j) = Yleft(j,1)   !
!!$             Y0zR(j) = Yright(j,1)  !
!!$
!!$             X0zL(j) = Xleft(j,1)
!!$             X0zR(j) = Xright(j,1)  !
!!$
!!$             YnzL(j) = Yleft(j,mz)  !
!!$             YnzR(j) = Yright(j,mz) !
!!$
!!$             XnzL(j) = Xleft(j,mz)
!!$             XnzR(j) = Xright(j,mz) !
!!$
!!$          end do
!!$
!!$          !stop
!!$
!!$          call Gen_Surf_Mesh(mx,my,X0zL,X0zR,X0zB,X0zT,Y0zL,Y0zR,Y0zB,Y0zT,Xfront,Yfront)
!!$          call Gen_Surf_Mesh(mx,my,XnzL,XnzR,XnzB,XnzT,YnzL,YnzR,YnzB,YnzT,Xback,Yback)
!!$
!!$          do j = 1,my
!!$             do i = 1,mx
!!$
!!$                Zfront(i, j) = az
!!$                Zback(i, j) = bz
!!$
!!$             end do
!!$          end do
!!$
!!$
!!$          ! y-face:
!!$          ! boundary curves at the top and bottom surfaces
!!$
!!$          do i = 1,mx
!!$
!!$             Z0yT(i) = Zback(i,my)
!!$             Z0yB(i) = Zfront(i,my)
!!$             ZnyT(i) = Zback(i,1)
!!$             ZnyB(i) = Zfront(i,1)
!!$
!!$          end do
!!$
!!$
!!$
!!$          do k = 1,mz
!!$
!!$             X0yL(k) = Xleft(my,k)
!!$             X0yR(k) = Xright(my,k)
!!$
!!$             Z0yL(k) = Zleft(my,k)
!!$             Z0yR(k) = Zright(my,k)
!!$
!!$             XnyL(k) = Xleft(1,k)
!!$             XnyR(k) = Xright(1,k)
!!$
!!$             ZnyL(k) = Zleft(1,k)
!!$             ZnyR(k) = Zright(1,k)
!!$
!!$          end do
!!$          !stop
!!$
!!$          do i = 1,mx
!!$
!!$             hx = (Xright(1,1)-Xleft(1,1))/real(mx-1, wp)
!!$             XnyB(i) = Xleft(1,1)+real(i-1, wp)*hx
!!$
!!$             hx = (Xright(my,1)-Xleft(my,1))/real(mx-1, wp)
!!$             X0yB(i) = Xleft(my,1)+real(i-1, wp)*hx
!!$
!!$             hx = (Xright(1,mz)-Xleft(1,mz))/real(mx-1, wp)
!!$             XnyT(i) = Xleft(1,mz)+real(i-1, wp)*hx
!!$
!!$             hx = (Xright(my,mz)-Xleft(my,mz))/real(mx-1, wp)
!!$             X0yT(i) = Xleft(my,mz)+real(i-1, wp)*hx
!!$
!!$
!!$             hz = (Zright(1,1)-Zleft(1,1))/real(mx-1, wp)
!!$             ZnyB(i) = Zleft(1,1)+real(i-1, wp)*hz
!!$
!!$             hz = (Zright(my,1)-Zleft(my,1))/real(mx-1, wp)
!!$             Z0yB(i) = Zleft(my,1)+real(i-1, wp)*hz
!!$
!!$             hz = (Zright(1,mz)-Zleft(1,mz))/real(mx-1, wp)
!!$             ZnyT(i) = Zleft(1,mz)+real(i-1, wp)*hz
!!$
!!$             hz = (Zright(my,mz)-Zleft(my,mz))/real(mx-1, wp)
!!$             Z0yT(i) = Zleft(my,mz)+real(i-1, wp)*hz
!!$
!!$          end do
!!$
!!$          call Gen_Surf_Mesh(mx,mz,X0yL,X0yR,X0yB,X0yT,Z0yL,Z0yR,Z0yB,Z0yT,Xbottom,Zbottom)
!!$          call Gen_Surf_Mesh(mx,mz,XnyL,XnyR,XnyB,XnyT,ZnyL,ZnyR,ZnyB,ZnyT,Xtop,Ztop)
!!$
!!$          do k = 1, mz
!!$             do i = 1,mx
!!$
!!$                Ybottom(i, k) = by
!!$                Ytop(i, k) = ay 
!!$             end do
!!$          end do

       end if

       if (use_topography .eqv.  .true.) then

          ! To implement topography, we modify the top surface  Ytop(i, k)

          do k = 1, mz
             do i = 1,mx

                Ybottom(i, k) = by
                Ytop(i, k) = ay

             end do
          end do

          select case(topography_type)

          case('analytical_topography')

             ! Implement the topography analytically

             hx = (bx-ax)/real(mx-1, wp)
             hz = (bz-az)/real(mz-1, wp)

             do k = 1, mz
                do i = 1,mx

                   Xtop(i, k) = ax + real(i-1, wp)*hx
                   Ztop(i, k) = az + real(k-1, wp)*hz

                   ! Analytical topogragraphy
                   f = 0.0_wp

                   r1 = sqrt((Xtop(i, k)-10.0_wp)**2 + (Ztop(i, k)-10.0_wp)**2)
                   r2 = sqrt((Xtop(i, k)-10.0_wp)**2 + (Ztop(i, k)-10.0_wp)**2)

                   if (r1 .le. 3.0_wp) then
                      f = -0.3_wp*(1.0_wp + cos(pi*r1/3.0_wp))

                   end if


                   Ytop(i, k) = ay - topo*1.5_wp*f 

                end do
             end do

          case('read_topo_from_memory')

             !read topography from memory exactly the same way fault geometry was treated
             call read_2darray(Topography, mx*mz, 3, topography_path)

             call reshape_2darray(XT, Topography(:,1), mz, mx)
             call reshape_2darray(YT, Topography(:,2), mz, mx)
             call reshape_2darray(ZT, Topography(:,3), mz, mx)

             !                                                   top-face: j = 1
             !
             !                                                        Top
             !                                          !-----------------------------!
             !                                          !                             ! 
             !                                     Left !z                            ! Right
             !                                          !                             !
             !                                          !             x               !
             !                                          !-----------------------------!
             !                                                        Bottom

             do k = 1, mz
                do i = 1,mx


                   ! User generated topography
                   !--------------------------
                   Xtop(i, k) = XT(k, i)
                   Ytop(i, k) = ay + topo*YT(k, i)
                   Ztop(i, k) = ZT(k, i)


                end do
             end do

          case default

             stop 'specify topography profile'

          end select

          hx = 1.0_wp/real(mx-1, wp)
          hy = 1.0_wp/real(my-1, wp)
          hz = 1.0_wp/real(mz-1, wp)

          !delta = 3.0_wp


          !                                                   left-face: i = 1

          !                                                        Top
          !-----------------------------!
          !                             ! 
          !                                     Left !z                            ! Right
          !                             !
          !             y               !
          !-----------------------------!
          !                                                        Bottom

          ! x-face
          do j = 1,my

             r = real(j-1, wp)*hy

             !call stretch_hyp(r, delta,stretch)

             stretch = r

             Y0xT(j) = Ytop(1,mz)+(by-Ytop(1, mz))*stretch
             Y0xB(j) = Ytop(1,1)+(by-Ytop(1, 1))*stretch

             Z0xT(j) = bz
             Z0xB(j) = az



          end do


          do k = 1,mz

             s = real(k-1, wp)*hz

             Z0xL(k) = az+(bz-az)*s
             Z0xR(k) = az+(bz-az)*s

             Y0xR(k) = by
             Y0xL(k) = Ytop(1, k) !ay

          end do

          !                                                   right-face: i = mx
          !
          !                                                        Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !z                            ! Right
          !                                          !                             !
          !                                          !             y               !
          !                                          !-----------------------------!
          !                                                        Bottom

          hx = 1.0_wp/real(mx-1, wp)
          hy = 1.0_wp/real(my-1, wp)
          hz = 1.0_wp/real(mz-1, wp)

          do j = 1,my

             r = real(j-1, wp)*hy

             !call stretch_hyp(r, delta, stretch)

             stretch = r

             YnxT(j) = Ytop(mx,mz)+(by-Ytop(mx, mz))*stretch
             YnxB(j) = Ytop(mx,1) +(by-Ytop(mx, 1) )*stretch

             ZnxT(j) = bz
             ZnxB(j) = az
          end do

          do k = 1,mz

             s = real(k-1, wp)*hz

             ZnxL(k) = az+(bz-az)*s
             ZnxR(k) = az+(bz-az)*s

             YnxR(k) = by
             YnxL(k) = Ytop(mx, k) !ay

          end do




          call Gen_Surf_Mesh(my,mz,Y0xL,Y0xR,Y0xB,Y0xT,Z0xL,Z0xR,Z0xB,Z0xT,Yleft,Zleft)
          call Gen_Surf_Mesh(my,mz,YnxL,YnxR,YnxB,YnxT,ZnxL,ZnxR,ZnxB,ZnxT,Yright,Zright)


          do k = 1,mz
             do j = 1,my


                Xleft(j, k)  = ax
                Xright(j, k) = bx

                !==============================================

             end do
          end do

          select case(profile_type)

          case('no_path')


             do k = 1,mz
                do j = 1,my


                   Xleft(j, k)  = ax
                   Xright(j, k) = bx

                   !==============================================

                end do
             end do

          case('read_from_memomry_fractal')

             !read fault geometry from memory
             call read_2darray(Fault_Geometry, ny*nz,3,profile_path)
             call reshape_2darray(X1, Fault_Geometry(:,1), ny,nz)

             call reshape_2darray(Y1, Fault_Geometry(:,2), ny,nz)
             call reshape_2darray(Z1, Fault_Geometry(:,3), ny,nz)



             !stop

             do k = 1,mz
                do j = 1,my

                   gg = 0.0_wp



                   if (rc /=0.0_wp) then

                      y =  Yleft(j, k)
                      z =  Zleft(j, k)



                      ! interpolate fault surface to match the free-surface topography
                      !if (minval(Y1(1,:)) .le. minval(Yleft(1, :))) then
                      !call interpolate_fault_surface(X1, Y1, Z1, y, z, ny, nz, gg)

                      !XX1(j, k) = gg
                      XX1(j, k) = X1(j,k)
                      YY1(j, k) = Y1(j,k)
                      ZZ1(j, k) = Z1(j,k)



                      !else

                      !print *, minval(Y1(1,:)),  minval(Yleft(1,:))
                      !STOP 'Fault not compatible with topography'
                      !end if
                   end if

                   if (lc /=0.0_wp) then

                      y =  Yright(j, k)
                      z =  Zright(j, k)

                      ! interpolate fault surface to match the free-surface topography
                      !if (minval(Y1(1,:)) .le. minval(Yright(1, :))) then
                      !   call interpolate_fault_surface(X1, Y1, Z1, y, z, ny, nz, gg)

                      !XX1(j, k) = gg
                      XX1(j, k) = X1(j, k)
                      YY1(j, k) = Y1(j, k)
                      ZZ1(j, k) = Z1(j, k)

                      !print *, gg, y, z


                      !else

                      !   print *, minval(Y1(1,:)),  minval(Yright(1,:))
                      !   STOP 'Fault not compatible with topography'
                      !end if

                   end if

                end do

             end do


             do k = 1,mz
                do j = 1,my

                   Xleft(j, k) = ax
                   Xright(j, k) = bx

                   gg = 0.0_wp

                   if (rc /=0.0_wp) then

                      y =  Yleft(j, k)
                      z =  Zleft(j, k)

                      ! contrain fault against topography
                      !if (minval(Y1(1,:)) .le. minval(Yleft(1, :))) then


                      Xleft(j, k) = ax  + XX1(j, k)
                      Yleft(j, k) = Y1(j, k)
                      Zleft(j, k) = Z1(j, k) 


                      !else

                      !print *, minval(Y1(1,:)),  minval(Yleft(1,:))
                      !STOP 'Fault not compatible with topography'
                      !end if
                   end if

                   if (lc /=0.0_wp) then

                      y =  Yright(j, k)
                      z =  Zright(j, k)

                      ! contrain fault against topography
                      !if (minval(Y1(1,:)) .le. minval(Yright(1, :))) then


                      Xright(j, k) = bx  + XX1(j, k)
                      Yright(j, k) = Y1(j, k)
                      Zright(j, k) = Z1(j, k) 


                      !else

                      !print *, minval(Y1(1,:)),  minval(Yright(1,:))
                      !STOP 'Fault not compatible with topography'
                      !end if

                   end if


                end do


             end do

          case default

             stop 'specify fault profile'

          end select

          ! z-face

          !                                                   front-face: k = 1
          !
          !                                                        Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !y                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom



          ! z-face
          do i = 1,mx

             q = real(i-1, wp)*hx

!!$             if (lc /=0.0_wp) then
!!$                call stretch_hyp(q, delta, stretch)
!!$             end if
!!$             
!!$             if (rc /=0.0_wp) then
!!$                call stretch_hyp(1.0_wp-q, delta, stretch)
!!$                stretch = 1.0_wp-stretch
!!$             end if

             stretch = q

             X0zT(i) = Xleft(my,1)+(Xright(my,1)-Xleft(my,1))*stretch

             Y0zT(i) = Ybottom(i, 1)

             X0zB(i) = Xtop(i, 1)
             Y0zB(i) = Ytop(i, 1) 


          end do

          !stop

          do j = 1,my

             X0zL(j) = Xleft(j,1)
             X0zR(j) = Xright(j,1)  

             Y0zL(j) = Yleft(j,1)   !
             Y0zR(j) = Yright(j,1)  !

          end do


          !                                                   back-face: k = mz
          !
          !                                                        Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !y                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom



          ! z-face
          do i = 1,mx

             q = real(i-1, wp)*hx

!!$             if (lc /=0.0_wp) call stretch_hyp(q, delta, stretch)
!!$             
!!$             if (rc /=0.0_wp) then
!!$                call stretch_hyp(1.0_wp-q, delta, stretch)
!!$                stretch = 1.0_wp-stretch
!!$             end if

             stretch = q

             XnzT(i) = Xleft(my,mz)+ (Xright(my,mz)-Xleft(my,mz))*stretch
             YnzT(i) = Ybottom(i, mz)

             XnzB(i) = Xtop(i, mz)
             YnzB(i) = Ytop(i, mz) 


          end do

          !stop


          do j = 1,my

             XnzL(j) = Xleft(j, mz)
             XnzR(j) = Xright(j, mz)  

             YnzL(j) = Yleft(j, mz)   !
             YnzR(j) = Yright(j, mz)  !

          end do




          call Gen_Surf_Mesh(mx,my,X0zL,X0zR,X0zB,X0zT,Y0zL,Y0zR,Y0zB,Y0zT,Xfront,Yfront)
          call Gen_Surf_Mesh(mx,my,XnzL,XnzR,XnzB,XnzT,YnzL,YnzR,YnzB,YnzT,Xback,Yback)





          do j = 1,my
             do i = 1,mx

                Zfront(i, j) = az
                Zback(i, j) = bz

             end do
          end do


          ! y-face:
          ! boundary curves at the top and bottom surfaces

          !                                                   bottom-face: j = my
          !
          !                                                        Top
          !                                          !-----------------------------!
          !                                          !                             ! 
          !                                     Left !z                            ! Right
          !                                          !                             !
          !                                          !             x               !
          !                                          !-----------------------------!
          !                                                        Bottom


          do i = 1,mx

             XnyT(i) = Xback(i,my)
             XnyB(i) = Xfront(i,my)

             ZnyT(i) = Zback(i,my)
             ZnyB(i) = Zfront(i,my)


          end do

          do k = 1,mz

             XnyL(k) = Xleft(my,k)
             XnyR(k) = Xright(my,k)

             ZnyL(k) = Zleft(my,k)
             ZnyR(k) = Zright(my,k)

          end do



          ! generate remaining boundary curves for the bottom boundary 
          call Gen_Surf_Mesh(mx,mz,XnyL,XnyR,XnyB,XnyT,ZnyL,ZnyR,ZnyB,ZnyT,Xbottom,Zbottom)


       end if



    end if

    ! Generate the x, y and z components of the grid
    call Gen3_Curve_Mesh(G,Xleft,Xright,Xtop,Xbottom,Xfront,Xback,1)
    call Gen3_Curve_Mesh(G,Yleft,Yright,Ytop,Ybottom,Yfront,Yback,2)
    call Gen3_Curve_Mesh(G,Zleft,Zright,Ztop,Zbottom,Zfront,Zback,3)

    select case(trim(profile_type))
    case('analytical_tpv36_a', 'analytical_tpv36_ver_a')
       call TPV36_a(G, lc, rc, az, bz)
    case('analytical_tpv36_b', 'analytical_tpv36_ver_b')
       call TPV36_b(G, lc, rc, az, bz)
    case('analytical_tpv36_c', 'analytical_tpv36_ver_c')
       call TPV36_c(G, lc, rc, ax, bx, az, bz)
    case('analytical_tpv36_d', 'analytical_tpv36_ver_d')
       call TPV36_d(G, lc, rc, ax, bx)
    case('analytical_tpv36_e', 'analytical_tpv36_ver_e')
       call TPV36_e(G, lc, rc, az, bz)
    end select

    ! Generate the x, y and z components of the grid
    !call Gen3_Curve_Mesh(G,Xleft,Xright,Xbottom,Xtop,Xfront,Xback,1)
    !call Gen3_Curve_Mesh(G,Yleft,Yright,Ybottom,Ytop,Yfront,Yback,2)
    !call Gen3_Curve_Mesh(G,Zleft,Zright,Zbottom,Ztop,Zfront,Zback,3)


    ! free memory
    deallocate(Xfront,Xback,Yfront,Yback,Zfront,Zback)
    deallocate(Xtop,Xbottom,Ytop,Ybottom,Ztop,Zbottom)
    deallocate(Xleft,Xright,Yleft,Yright,Zleft,Zright)
    deallocate(Fault_Geometry,Topography, X1,Y1,Z1,XT,YT,ZT)
    deallocate(XX1,YY1,ZZ1)

    deallocate(Y0xL,Y0xR,Z0xL,Z0xR)
    deallocate(Y0xB,Y0xT,Z0xB,Z0xT)

    deallocate(X0yL,X0yR,Z0yL,Z0yR)
    deallocate(X0yB,X0yT,Z0yB,Z0yT)

    deallocate(X0zL,X0zR,Y0zL,Y0zR)
    deallocate(X0zB,X0zT,Y0zB,Y0zT)

    deallocate(YnxL,YnxR,ZnxL,ZnxR)
    deallocate(YnxB,YnxT,ZnxB,ZnxT)

    deallocate(XnyL,XnyR,ZnyL,ZnyR)
    deallocate(XnyB,XnyT,ZnyB,ZnyT)

    deallocate(XnzL,XnzR,YnzL,YnzR)
    deallocate(XnzB,XnzT,YnzB,YnzT)

  end subroutine Curve_Grid_3D3


  subroutine Curve_Grid_3D(G,ax,bx,ay,by,az,bz,lc,rc,use_mms,profile_type,profile_path)
    ! Generate a 3D curvilinear mesh using the transfinite interpolation grid generator
    ! implemented in Gen_Curve_Mesh below. Note that any user defined smooth surfaces can be
    ! used

    type(block_grid_t),intent(inout) :: G
    integer :: mx,my,mz,px,py,pz,nq,nr,ns            ! number of grid points
    real(kind = wp), intent(in) :: ax,bx, ay,by, az,bz          ! dimensions of a 3D computational elastic block
    real(kind = wp), intent(in) :: lc, rc                       ! lc, rc  add boundary curves to the x-boundaries
    logical,intent(in) :: use_mms
    character(64), intent(in) :: profile_type,profile_path

    real(kind = wp), dimension(:,:,:), allocatable:: X, Y, Z      ! x y z components of the grid
    ! Six boundary surfaces for a 3D physical domain
    real(kind = wp),dimension(:,:  ),allocatable :: Xleft,Xright,Xtop,Xbottom,Xfront,Xback
    real(kind = wp),dimension(:,:  ),allocatable :: Yleft,Yright,Ytop,Ybottom,Yfront,Yback
    real(kind = wp),dimension(:,:  ),allocatable :: Zleft,Zright,Ztop,Zbottom,Zfront,Zback

    real(kind = wp),dimension(:,:  ),allocatable :: Fault_Geometry
    real(kind = wp),dimension(:,:  ),allocatable :: X1,Y1,Z1

    real(kind = wp) :: hx, hy, hz
    integer :: i, j, k
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp) :: r1, r2, f
    real(kind = wp) :: r, r_nuke
    real(kind = wp) :: r0y, rny, r0z, rnz

    nq = G%C%nq
    nr = G%C%nr
    ns = G%C%ns

    allocate(X(nq,nr,ns),Y(nq,nr,ns),Z(nq,nr,ns))

    allocate(Xfront(nq,nr),Xback(nq,nr),Yfront(nq,nr),&
         Yback(nq,nr),Zfront(nq,nr),Zback(nq,nr))

    allocate(Xtop(nq,ns),Xbottom(nq,ns),Ytop(nq,ns),&
         Ybottom(nq,ns),Ztop(nq,ns),Zbottom(nq,ns))

    allocate(Xleft(nr,ns),Xright(nr,ns),Yleft(nr,ns),&
         Yright(nr,ns),Zleft(nr,ns),Zright(nr,ns))
    allocate(Fault_Geometry(nr*ns, 3))
    allocate(X1(nr,ns),Y1(nr,ns),Z1(nr,ns))

    hx = G%bhq
    hy = G%bhr
    hz = G%bhs

    r0y = 1.0_wp
    r0z = 1.0_wp
    rny = 1.0_wp
    rnz = 1.0_wp

    if (use_mms .NEQV.  .TRUE.) then

       select case(profile_type)

       case('analytical_tpv28')

          do k = 1,ns
             do j = 1,nr

                Yleft(j, k) = ay+real(j-1, wp)*hy;
                Zleft(j, k) = az+real(k-1, wp)*hz;
                Yright(j, k) = ay+real(j-1, wp)*hy;
                Zright(j, k) = az+real(k-1, wp)*hz;
                Xleft(j, k) = ax
                Xright(j, k) = bx

                ! Geometry for TPV28
                f = 0.0_wp
                r1 = sqrt((Zleft(j, k)+10.5_wp)**2 + (Yleft(j, k)-7.5_wp)**2)
                r2 = sqrt((Zleft(j, k)-10.5_wp)**2 + (Yleft(j, k)-7.5_wp)**2)
                if (r1 .le. 3.0_wp) f = -0.3_wp*(1.0_wp + cos(pi*r1/3.0_wp))
                if (r2 .le. 3.0_wp) f = -0.3_wp*(1.0_wp + cos(pi*r2/3.0_wp))

                Xleft(j, k) =  Xleft(j, k) + rc*f
                Xright(j, k) =  Xright(j, k) + lc*f

                !==============================================

             end do
          end do


       case('analytical_test_problem')

          do k = 1,ns
             do j = 1,nr

                Yleft(j, k) = ay+real(j-1, wp)*hy;
                Zleft(j, k) = az+real(k-1, wp)*hz;
                Yright(j, k) = ay+real(j-1, wp)*hy;
                Zright(j, k) = az+real(k-1, wp)*hz;
                Xleft(j, k) = ax
                Xright(j, k) = bx

                ! Non-planar fault
                Xleft(j, k) = ax&
                     + 0.1_wp*rc*exp(-((1.0_wp*Yleft(j, k)-0.5_wp*(ay+by))**2 + 1.0_wp*(Zleft(j, k)-0.5_wp*(az+bz))**2)/0.025_wp)&
                     + 0.1_wp*rc*sin(2.0_wp/(by-ay)*pi*Yleft(j, k))&
                     *sin(2.0_wp/(bz-az)*pi*Zleft(j, k))

                Xright(j, k) = bx&
                     + 0.1_wp*lc*exp(-(1.0_wp*(Yright(j, k)-0.5_wp*(ay+by))**2 + 1.0_wp*(Zright(j, k)-0.5_wp*(az+bz))**2)/0.025_wp)&
                     + 0.1_wp*lc*sin(2.0_wp/(by-ay)*pi*Yright(j, k))&
                     *sin(2.0_wp/(bz-az)*pi*Zright(j, k))


                !==============================================

             end do
          end do

       case('read_from_memomry_fractal')

          !read fault geometry from memory

          call read_2darray(Fault_Geometry, nr*ns,3,profile_path)

          call reshape_2darray(X1, Fault_Geometry(:,1), nr, ns)
          call reshape_2darray(Y1, Fault_Geometry(:,2), nr, ns)
          call reshape_2darray(Z1, Fault_Geometry(:,3), nr, ns)

          do k = 1,ns
             do j = 1,nr
                !
                Yleft(j, k) = ay+real(j-1, wp)*hy;
                Zleft(j, k) = az+real(k-1, wp)*hz;
                Yright(j, k) = ay+real(j-1, wp)*hy;
                Zright(j, k) = az+real(k-1, wp)*hz;
                Xleft(j, k) = ax
                Xright(j, k) = bx
                !
                ! construct smooth functions
                if (abs(Yleft(j,k)-ay) <= 1d-9) then
                   r0y = 0.0_wp
                else
                   r0y = exp(-2.0_wp/(abs(Yleft(j,k)-ay)))
                end if

                if (abs(Yleft(j,k)-by) <= 1d-9) then
                   rny = 0.0_wp
                else
                   rny = exp(-2.0_wp/(abs(Yleft(j,k)-by)))
                end if


                if (abs(Zleft(j,k)-az) <= 1d-9) then
                   r0z = 0.0_wp
                else
                   r0z = exp(-2.0_wp/(abs(Zleft(j,k)-az)))
                end if

                if (abs(Zleft(j,k)-bz) <= 1d-9) then
                   rnz = 0.0_wp
                else
                   rnz = exp(-2.0_wp/(abs(Zleft(j,k)-bz)))
                end if

                ! Add geometry for fractal faults with smoothed edges

                if (rc /=0.0_wp) then
                   Xleft(j, k) = ax+r0y*rny*r0z*rnz*X1(j,k)
                   Yleft(j, k) = Y1(j,k)
                   Zleft(j, k) = Z1(j,k)
                end if

                if (lc /=0.0_wp) then
                   Xright(j, k) = bx+r0y*rny*r0z*rnz*X1(j,k)
                   Yright(j, k) = Y1(j,k)
                   Zright(j, k) = Z1(j,k)
                end if

                !==============================================
             end do
          end do

       case default

          stop 'specify fault profile'

       end select
       do j = 1, nr
          do i = 1,nq

             Xfront(i, j) = ax+real(i-1, wp)*hx
             Yfront(i, j) = ay+real(j-1, wp)*hy
             Xback(i, j) = ax+real(i-1, wp)*hx
             Yback(i, j) = ay+real(j-1, wp)*hy

             Zfront(i, j) = az&
                  + 0.0_wp*rc*exp(-((1.0_wp*Yfront(i, j)-0.5_wp)**2 + 1.0_wp*(Xfront(i, j)-0.5_wp)**2)/0.025_wp)
             Zback(i, j) = bz&
                  + 0.0_wp*lc*exp(-(1.0_wp*(Yback(i, j)-0.5_wp)**2 + 1.0_wp*(Xback(i, j)+0.5_wp)**2)/0.025_wp)

          end do
       end do

       do k = 1, ns
          do i = 1,nq

             Xtop(i, k) = ax+real(i-1, wp)*hx;
             Ztop(i, k) = az+real(k-1, wp)*hz;
             Xbottom(i, k) = ax+real(i-1, wp)*hx;
             Zbottom(i, k) = az+real(k-1, wp)*hz;

             Ytop(i, k) = ay;
             Ybottom(i, k) = by

          end do
       end do

    end if

    if (use_mms .EQV.  .TRUE.) then
       do k = 1, ns
          do j = 1, nr

             Yleft(j, k) = ay+real(j-1, wp)*hy;
             Zleft(j, k) = az+real(k-1, wp)*hz;
             Yright(j, k) = ay+real(j-1, wp)*hy;
             Zright(j, k) = az+real(k-1, wp)*hz;
             Xleft(j, k) = ax
             Xright(j, k) = bx

             ! Non-planar fault
             Xleft(j, k) = ax&
                  + 0.1_wp*rc*exp(-((1.0_wp*Yleft(j, k)-0.5_wp*(ay+by))**2 + 1.0_wp*(Zleft(j, k)-0.5_wp*(az+bz))**2)/0.025_wp)&
                  + 0.1_wp*rc*sin(2.0_wp/(by-ay)*pi*Yleft(j, k))&
                  *sin(2.0_wp/(bz-az)*pi*Zleft(j, k))

             Xright(j, k) = bx&
                  + 0.1_wp*lc*exp(-(1.0_wp*(Yright(j, k)-0.5_wp*(ay+by))**2 + 1.0_wp*(Zright(j, k)-0.5_wp*(az+bz))**2)/0.025_wp)&
                  + 0.1_wp*lc*sin(2.0_wp/(by-ay)*pi*Yright(j, k))&
                  *sin(2.0_wp/(bz-az)*pi*Zright(j, k))

          end do
       end do

       do j = 1, nr
          do i = 1,nq

             Xfront(i, j) = ax+real(i-1, wp)*hx
             Yfront(i, j) = ay+real(j-1, wp)*hy
             Xback(i, j) = ax+real(i-1, wp)*hx
             Yback(i, j) = ay+real(j-1, wp)*hy

             Zfront(i, j) = az&
                  + 0.0_wp*rc*exp(-((1.0_wp*Yfront(i, j)-0.5_wp)**2 + 1.0_wp*(Xfront(i, j)-0.5_wp)**2)/0.025_wp)
             Zback(i, j) = bz&
                  + 0.0_wp*lc*exp(-(1.0_wp*(Yback(i, j)-0.5_wp)**2 + 1.0_wp*(Xback(i, j)+0.5_wp)**2)/0.025_wp)

          end do
       end do

       do k = 1, ns
          do i = 1, nq

             Xtop(i, k) = ax+real(i-1, wp)*hx;
             Ztop(i, k) = az+real(k-1, wp)*hz;
             Xbottom(i, k) = ax+real(i-1, wp)*hx;
             Zbottom(i, k) = az+real(k-1, wp)*hz;

             Ytop(i, k) = ay;
             Ybottom(i, k) = by

          end do
       end do
    end if

    ! mx = G%C%mq
    ! px = G%C%pq
    ! my = G%C%mr
    ! py = G%C%pr
    ! mz = G%C%ms
    ! pz = G%C%ps

    ! Generate the x, y and z components of the grid
    !call Gen3_Curve_Mesh(G%x(mx:px,my:py,mz:pz,1),Xleft,Xright,Xtop,Xbottom,Xfront,Xback,mx,my,mz,px,py,pz,nq,nr,ns)
    !call Gen3_Curve_Mesh(G%x(mx:px,my:py,mz:pz,2),Yleft,Yright,Ytop,Ybottom,Yfront,Yback,mx,my,mz,px,py,pz,nq,nr,ns)
    !call Gen3_Curve_Mesh(G%x(mx:px,my:py,mz:pz,3),Zleft,Zright,Ztop,Zbottom,Zfront,Zback,mx,my,mz,px,py,pz,nq,nr,ns)

    ! Generate the x, y and z components of the grid

    call Gen_Curve_Mesh(X,Xleft,Xright,Xtop,Xbottom,Xfront,Xback,mx, my, mz)
    call Gen_Curve_Mesh(Y,Yleft,Yright,Ytop,Ybottom,Yfront,Yback,mx, my, mz)
    call Gen_Curve_Mesh(Z,Zleft,Zright,Ztop,Zbottom,Zfront,Zback,mx, my, mz)

    ! Store the grid in a higher rank array
    G%x(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, 1) = X(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps)
    G%x(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, 2) = Y(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps)
    G%x(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, 3) = Z(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps)

    ! Store the grid in a higher rank array
    !G%x(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, 1) = X(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps)
    !G%x(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, 2) = Y(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps)
    !G%x(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps, 3) = Z(G%C%mq:G%C%pq, G%C%mr:G%C%pr, G%C%ms:G%C%ps)

    ! free memory
    !deallocate(X,Y,Z)
    deallocate(Xfront,Xback,Yfront,Yback,Zfront,Zback)
    deallocate(Xtop,Xbottom,Ytop,Ybottom,Ztop,Zbottom)
    deallocate(Xleft,Xright,Yleft,Yright,Zleft,Zright)
    deallocate(Fault_Geometry,X1,Y1,Z1)

  end subroutine Curve_Grid_3D

  subroutine Gen_Curve_Mesh(X,Xleft,Xright,Xtop,Xbottom,Xfront,Xback,mx,my,mz)

    ! Given the boundary surfaces of a physical space in 3D, compute a component
    ! of the grid x(q,r,s), y(q,r,s), or z(q,r,s) using transfinite interpolation

    implicit none

    real(kind = wp), dimension(:,:,:), intent(out) :: X                       ! the grid

    ! Six boundary surfaces for a 3D physical domain
    real(kind = wp), dimension(:,:), intent(in) :: Xleft,Xright,Xtop,Xbottom,Xfront,Xback
    integer, intent(in) :: mx, my, mz          ! number of grid points

    real(kind = wp) :: hq, hr, hs
    integer :: i, j, k
    real(kind = wp) :: U,V,W,UV,UW,VW,UVW1,UVW2

    hq = 1.0_wp/real(mx-1, wp);
    hr = 1.0_wp/real(my-1, wp);
    hs = 1.0_wp/real(mz-1, wp);


    do k = 1,mz
       do j = 1,my
          do i = 1,mx
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             U = (1.0_wp-real(i-1, wp)*hq)*Xleft(j,k) + real(i-1, wp)*hq*Xright(j,k)

             V = (1.0_wp-real(j-1, wp)*hr)*Xtop(i,k) + real(j-1, wp)*hr*Xbottom(i,k)

             W = (1.0_wp-real(k-1, wp)*hs)*Xfront(i,j) + real(k-1, wp)*hs*Xback(i,j)
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             UW = (1.0_wp-real(i-1, wp)*hq)*((1.0_wp-real(j-1, wp)*hr)*Xleft(1,k) &
                  + real(j-1, wp)*hr*Xleft(my,k))&
                  + real(i-1, wp)*hq*((1.0_wp-real(j-1, wp)*hr)*Xright(1, k) &
                  + real(j-1, wp)*hr*Xright(my,k))


             UV = (1.0_wp-real(j-1, wp)*hr)*((1.0_wp-real(k-1, wp)*hs)*Xtop(i,1) &
                  + real(k-1, wp)*hs*Xtop(i,mz))&
                  + real(j-1, wp)*hr*((1.0_wp-real(k-1, wp)*hs)*Xbottom(i,1) &
                  + real(k-1, wp)*hs*Xbottom(i,mz))

             VW = (1.0_wp-real(k-1, wp)*hs)*((1.0_wp-real(i-1, wp)*hq)*Xfront(1,j) &
                  + real(i-1, wp)*hq*Xfront(mx,j))&
                  + real(k-1, wp)*hs*((1.0_wp-real(i-1, wp)*hq)*Xback(1,j) &
                  + real(i-1, wp)*hq*Xback(mx,j))
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             UVW1 = (1.0_wp-real(i-1, wp)*hq)*((1.0_wp-real(j-1, wp)*hr)*((1.0_wp-real(k-1, wp)*hs)*Xtop(1,1) &
                  + real(k-1, wp)*hs*Xtop(1,Mz))&
                  + real(j-1, wp)*hr*((1.0_wp-real(k-1, wp)*hs)*Xbottom(1,1) &
                  + real(k-1, wp)*hs*Xbottom(1,mz)))

             UVW2 = real(i-1, wp)*hq*((1.0_wp-real(j-1, wp)*hr)*((1.0_wp-real(k-1, wp)*hs)*Xtop(mx,1) &
                  + real(k-1, wp)*hs*Xtop(mx,mz))&
                  + real(j-1, wp)*hr*((1.0_wp-real(k-1, wp)*hs)*Xbottom(mx,1) &
                  + real(k-1, wp)*hs*Xbottom(mx,mz)))
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             X(i, j, k) = U + V + W - UV - UW - VW + UVW1 + UVW2

          end do
       end do
    end do

  end subroutine Gen_Curve_Mesh

  pure real(kind=wp) function tpv36_ver_interface(y) result(x)
    ! C2 transition from the 15-degree planar fault to a locked vertical
    ! extension over Y=8--11 km; the padding and PML below remain straight.
    implicit none
    real(kind=wp), intent(in) :: y
    real(kind=wp), parameter :: pi = 3.141592653589793_wp
    real(kind=wp) :: slope, t, integral_s

    slope = 1.0_wp/tan(15.0_wp*pi/180.0_wp)
    if (y <= 8.0_wp) then
       x = y*slope
    else if (y < 11.0_wp) then
       t = (y-8.0_wp)/3.0_wp
       integral_s = t**6 - 3.0_wp*t**5 + 2.5_wp*t**4
       x = 8.0_wp*slope + 3.0_wp*slope*(t-integral_s)
    else
       x = 9.5_wp*slope
    end if
  end function tpv36_ver_interface

  pure real(kind=wp) function tpv36_transition_fraction(u, length, n, ds) result(f)
    ! Monotone C2 redistribution.  It matches ds-sized cells at both ends of
    ! a transition while absorbing unavoidable excess nodes smoothly inside.
    implicit none
    real(kind=wp), intent(in) :: u, length, ds
    integer, intent(in) :: n
    real(kind=wp), parameter :: edge = 0.1_wp
    real(kind=wp) :: endpoint_density, amp, base, tl, tr, hl, hr

    endpoint_density = real(n,wp)*ds/length
    amp = (endpoint_density-1.0_wp)/(1.0_wp-edge)
    base = 1.0_wp-amp*edge
    tl = min(max(u/edge,0.0_wp),1.0_wp)
    hl = tl-2.5_wp*tl**4+3.0_wp*tl**5-tl**6
    hr = 0.0_wp
    if (u > 1.0_wp-edge) then
       tr = min(max((1.0_wp-u)/edge,0.0_wp),1.0_wp)
       hr = 0.5_wp-(tr-2.5_wp*tr**4+3.0_wp*tr**5-tr**6)
    end if
    f = base*u + amp*edge*(hl+hr)
  end function tpv36_transition_fraction

  subroutine TPV36_a(G, lc, rc, az, bz)
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, az, bz
    call apply_tpv36_zoned_x_map(G, lc, rc, az, bz, .true.)
  end subroutine TPV36_a

  subroutine TPV36_b(G, lc, rc, az, bz)
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, az, bz
    call apply_tpv36_zoned_x_map(G, lc, rc, az, bz, .false.)
  end subroutine TPV36_b

  subroutine TPV36_c(G, lc, rc, ax, bx, az, bz)
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, ax, bx, az, bz
    call apply_tpv36_parallel_x_map(G, lc, rc, ax, bx, az, bz)
  end subroutine TPV36_c

  subroutine TPV36_d(G, lc, rc, ax, bx)
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, ax, bx
    call apply_tpv36_planar_x_map(G, lc, rc, ax, bx)
  end subroutine TPV36_d

  subroutine TPV36_e(G, lc, rc, az, bz)
    ! Hybrid fault-normal grid.  Begin with the ver-a zoned mesh, rotate the
    ! transverse direction to the fault normal at the 18 km hypocenter, and
    ! relax it exponentially toward the outer parts of each block.  Rotation
    ! fades independently toward the free surface and lower boundary.
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, az, bz
    call apply_tpv36_zoned_x_map(G, lc, rc, az, bz, .true.)
    call apply_tpv36_fault_normal_buffer(G)
  end subroutine TPV36_e

  pure real(kind=wp) function tpv36_smootherstep(t) result(s)
    implicit none
    real(kind=wp), intent(in) :: t
    real(kind=wp) :: u
    u = min(max(t,0.0_wp),1.0_wp)
    s = u**3*(10.0_wp-15.0_wp*u+6.0_wp*u**2)
  end function tpv36_smootherstep

  subroutine apply_tpv36_fault_normal_buffer(G)
    implicit none
    type(block_grid_t), intent(inout) :: G
    integer :: i, j, k, mx, px, my, py, mz, pz
    real(kind=wp), parameter :: relax_length = 2.0_wp
    real(kind=wp) :: dip, yh, ybottom, y0, xf, d, ad, depth_w, radial_w
    real(kind=wp) :: beta, x0

    mx = G%C%mq; px = G%C%pq
    my = G%C%mr; py = G%C%pr
    mz = G%C%ms; pz = G%C%ps
    dip = 15.0_wp*3.141592653589793_wp/180.0_wp
    yh = 18.0_wp*sin(dip)
    ybottom = 16.0_wp

    do k = mz,pz
       do j = my,py
          y0 = G%x(mx,j,k,2)
          if (y0 <= yh) then
             depth_w = tpv36_smootherstep(y0/yh)
          else
             depth_w = 1.0_wp-tpv36_smootherstep((y0-yh)/(ybottom-yh))
          end if
          xf = tpv36_ver_interface(y0)
          do i = mx,px
             x0 = G%x(i,j,k,1)
             d = x0-xf
             ad = abs(d)
             ! Unit normal with positive X component is
             ! (sin(dip),-cos(dip)).  d*radial_w approaches a bounded 2 km
             ! displacement, avoiding the foldover caused by a compact taper.
             if (ad > tiny(1.0_wp)) then
                radial_w = relax_length/ad*(1.0_wp-exp(-ad/relax_length))
             else
                radial_w = 1.0_wp
             end if
             beta = depth_w*radial_w
             G%x(i,j,k,1) = xf+d*(1.0_wp-beta*(1.0_wp-sin(dip)))
             G%x(i,j,k,2) = y0-d*beta*cos(dip)
          end do
       end do
    end do
  end subroutine apply_tpv36_fault_normal_buffer

  subroutine apply_tpv36_zoned_x_map(G, lc, rc, az, bz, with_buffer)
    ! X-only TPV36 map with straight 2 km padding and 3 km PML zones.
    ! ver-a reserves a 2 km X-parallel buffer beside the fault; ver-b does not.
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, az, bz
    logical, intent(in) :: with_buffer
    integer :: i, j, k, mx, px, my, py, mz, pz, ni
    integer :: nbuf, npad, npml, ntransition, i1, i2, i3
    integer :: nj, nactive_y, ntransition_y, j1, j2
    real(kind=wp) :: ds, u, y, xf, dip

    mx = G%C%mq; px = G%C%pq
    my = G%C%mr; py = G%C%pr
    mz = G%C%ms; pz = G%C%ps
    ni = G%C%nq-1
    ds = abs(bz-az)/real(G%C%ns-1,wp)
    npml = max(1,nint(3.0_wp/ds))
    npad = max(1,nint(2.0_wp/ds))
    nbuf = 0
    if (with_buffer) nbuf = max(1,nint(2.0_wp/ds))
    ntransition = ni-npml-npad-nbuf
    if (ntransition < 1) error stop 'TPV36 zoned X map has too few q intervals'

    ! Redistribute Y so the target fault has <=ds arc-length spacing and the
    ! final 2 km padding and 3 km PML have their stated physical thicknesses.
    dip = 15.0_wp*3.141592653589793_wp/180.0_wp
    nj = G%C%nr-1
    nactive_y = ceiling(8.0_wp/(ds*sin(dip)))
    ntransition_y = nj-nactive_y-npad-npml
    if (ntransition_y < 1) error stop 'TPV36 zoned Y map has too few r intervals'
    j1 = nactive_y
    j2 = j1+ntransition_y

    do k = mz,pz
       do j = my,py
          if (j-1 <= j1) then
             y = 8.0_wp*real(j-1,wp)/real(nactive_y,wp)
          else if (j-1 <= j2) then
             y = 8.0_wp + 3.0_wp*real(j-1-j1,wp)/real(ntransition_y,wp)
          else if (j-1 <= j2+npad) then
             y = 11.0_wp + 2.0_wp*real(j-1-j2,wp)/real(npad,wp)
          else
             y = 13.0_wp + 3.0_wp*real(j-1-j2-npad,wp)/real(npml,wp)
          end if
          G%x(mx:px,j,k,2) = y
          xf = tpv36_ver_interface(y)
          do i = mx,px
             if (rc > 0.0_wp) then
                ! Block 1: outer boundary, PML, padding, transition, buffer, fault.
                i1 = npml; i2 = i1+npad; i3 = i2+ntransition
                if (i-1 <= i1) then
                   u = real(i-1,wp)/real(i1,wp)
                   G%x(i,j,k,1) = -20.0_wp + 3.0_wp*u
                else if (i-1 <= i2) then
                   u = real(i-1-i1,wp)/real(npad,wp)
                   G%x(i,j,k,1) = -17.0_wp + 2.0_wp*u
                else if (i-1 <= i3) then
                   u = real(i-1-i2,wp)/real(ntransition,wp)
                   G%x(i,j,k,1) = -15.0_wp + (xf-real(nbuf,wp)*ds+15.0_wp)* &
                        tpv36_transition_fraction(u, &
                        xf-real(nbuf,wp)*ds+15.0_wp,ntransition,ds)
                else
                   u = real(i-1-i3,wp)/real(max(1,nbuf),wp)
                   G%x(i,j,k,1) = xf-real(nbuf,wp)*ds + real(nbuf,wp)*ds*u
                end if
             else if (lc > 0.0_wp) then
                ! Block 2: fault, buffer, transition, padding, PML, outer boundary.
                i1 = nbuf; i2 = i1+ntransition; i3 = i2+npad
                if (nbuf > 0 .and. i-1 <= i1) then
                   u = real(i-1,wp)/real(nbuf,wp)
                   G%x(i,j,k,1) = xf + real(nbuf,wp)*ds*u
                else if (i-1 <= i2) then
                   u = real(i-1-i1,wp)/real(ntransition,wp)
                   G%x(i,j,k,1) = xf+real(nbuf,wp)*ds + &
                        (53.0_wp-xf-real(nbuf,wp)*ds)* &
                        tpv36_transition_fraction(u, &
                        53.0_wp-xf-real(nbuf,wp)*ds,ntransition,ds)
                else if (i-1 <= i3) then
                   u = real(i-1-i2,wp)/real(npad,wp)
                   G%x(i,j,k,1) = 53.0_wp + 2.0_wp*u
                else
                   u = real(i-1-i3,wp)/real(npml,wp)
                   G%x(i,j,k,1) = 55.0_wp + 3.0_wp*u
                end if
             end if
          end do
       end do
    end do
  end subroutine apply_tpv36_zoned_x_map

  subroutine apply_tpv36_parallel_x_map(G, lc, rc, ax, bx, az, bz)
    ! Entire block follows the fault at a constant X offset.  This gives
    ! uniform transverse spacing, but deliberately warps the side PMLs and
    ! makes the exterior X boundaries parallel to the locked extension.
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, ax, bx, az, bz
    integer :: i, j, k, mx, px, my, py, mz, pz
    integer :: nj, npad, npml, nactive_y, ntransition_y, j1, j2
    real(kind=wp) :: ds, dip, y, xf, q, offset

    mx = G%C%mq; px = G%C%pq
    my = G%C%mr; py = G%C%pr
    mz = G%C%ms; pz = G%C%ps
    ds = abs(bz-az)/real(G%C%ns-1,wp)
    npml = max(1,nint(3.0_wp/ds))
    npad = max(1,nint(2.0_wp/ds))
    dip = 15.0_wp*3.141592653589793_wp/180.0_wp
    nj = G%C%nr-1
    nactive_y = ceiling(8.0_wp/(ds*sin(dip)))
    ntransition_y = nj-nactive_y-npad-npml
    if (ntransition_y < 1) error stop 'TPV36 parallel map has too few r intervals'
    j1 = nactive_y
    j2 = j1+ntransition_y

    do k = mz,pz
       do j = my,py
          if (j-1 <= j1) then
             y = 8.0_wp*real(j-1,wp)/real(nactive_y,wp)
          else if (j-1 <= j2) then
             y = 8.0_wp + 3.0_wp*real(j-1-j1,wp)/real(ntransition_y,wp)
          else if (j-1 <= j2+npad) then
             y = 11.0_wp + 2.0_wp*real(j-1-j2,wp)/real(npad,wp)
          else
             y = 13.0_wp + 3.0_wp*real(j-1-j2-npad,wp)/real(npml,wp)
          end if
          xf = tpv36_ver_interface(y)
          do i = mx,px
             q = real(i-1,wp)/real(G%C%nq-1,wp)
             offset = ax+(bx-ax)*q
             if (rc > 0.0_wp .or. lc > 0.0_wp) then
                G%x(i,j,k,1) = xf+offset
                G%x(i,j,k,2) = y
             end if
          end do
       end do
    end do
  end subroutine apply_tpv36_parallel_x_map

  subroutine apply_tpv36_planar_x_map(G, lc, rc, ax, bx)
    ! Straight 15-degree continuation through the full Y extent.  All X-grid
    ! lines, exterior boundaries, and side PMLs are constant fault offsets.
    implicit none
    type(block_grid_t), intent(inout) :: G
    real(kind=wp), intent(in) :: lc, rc, ax, bx
    integer :: i, j, k, mx, px, my, py, mz, pz
    real(kind=wp) :: q, y, xf, offset, dip

    mx = G%C%mq; px = G%C%pq
    my = G%C%mr; py = G%C%pr
    mz = G%C%ms; pz = G%C%ps
    dip = 15.0_wp*3.141592653589793_wp/180.0_wp
    do k = mz,pz
       do j = my,py
          y = G%x(mx,j,k,2)
          xf = y/tan(dip)
          do i = mx,px
             q = real(i-1,wp)/real(G%C%nq-1,wp)
             offset = ax+(bx-ax)*q
             if (rc > 0.0_wp .or. lc > 0.0_wp) G%x(i,j,k,1) = xf+offset
          end do
       end do
    end do
  end subroutine apply_tpv36_planar_x_map

  subroutine Gen3_Curve_Mesh(G,Xleft,Xright,Xtop,Xbottom,Xfront,Xback,n)

    ! Given the boundary surfaces of a physical space in 3D, compute a component
    ! of the grid x(q,r,s), y(q,r,s), or z(q,r,s) using transfinite interpolation

    implicit none

    type(block_grid_t),intent(inout) :: G

    ! Six boundary surfaces for a 3D physical domain
    real(kind = wp), dimension(:,:), intent(in) :: Xleft,Xright,Xtop,Xbottom,Xfront,Xback
    integer, intent(in) :: n
    integer :: mx,my,mz,px,py,pz,nq,nr,ns          ! number of grid points

    real(kind = wp) :: hq, hr, hs,  q,r,s
    integer :: i, j, k
    real(kind = wp) :: U,V,W,UV,UW,VW,UVW1,UVW2


    nq = G%C%nq
    nr = G%C%nr
    ns = G%C%ns

    hq = 1.0_wp/real(nq-1, wp);
    hr = 1.0_wp/real(nr-1, wp);
    hs = 1.0_wp/real(ns-1, wp);

    mx = G%C%mq
    px = G%C%pq
    my = G%C%mr
    py = G%C%pr
    mz = G%C%ms
    pz = G%C%ps


    do k = mz,pz
       do j = my,py
          do i = mx,px


             q = real(i-1, wp)*hq
             r = real(j-1, wp)*hr
             s = real(k-1, wp)*hs
             
             U = (1.0_wp-q)*Xleft(j,k)  + q*Xright(j,k)

             V = (1.0_wp-r)*Xtop(i,k)   + r*Xbottom(i,k)

             W = (1.0_wp-s)*Xfront(i,j) + s*Xback(i,j)
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             UW = (1.0_wp-q)*((1.0_wp-r)*Xleft(1,k) &
                  + r*Xleft(nr,k))&
                  + q*((1.0_wp-r)*Xright(1, k) &
                  + r*Xright(nr,k))


             UV = (1.0_wp-r)*((1.0_wp-s)*Xtop(i,1) &
                  + s*Xtop(i,ns))&
                  + r*((1.0_wp-s)*Xbottom(i,1) &
                  + s*Xbottom(i,ns))

             VW = (1.0_wp-s)*((1.0_wp-q)*Xfront(1,j) &
                  + q*Xfront(nq,j))&
                  + s*((1.0_wp-q)*Xback(1,j) &
                  + q*Xback(nq,j))
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             UVW1 = (1.0_wp-q)*((1.0_wp-r)*((1.0_wp-s)*Xtop(1,1) &
                  + s*Xtop(1,ns))&
                  + r*((1.0_wp-s)*Xbottom(1,1) &
                  + s*Xbottom(1,ns)))

             UVW2 = q*((1.0_wp-r)*((1.0_wp-s)*Xtop(nq,1) &
                  + s*Xtop(nq,ns))&
                  + r*((1.0_wp-s)*Xbottom(nq,1) &
                  + s*Xbottom(nq,ns)))
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$             U = (1.0_wp-real(i-1, wp)*hq)*Xleft(j,k) + real(i-1, wp)*hq*Xright(j,k)
!!$
!!$             V = (1.0_wp-real(j-1, wp)*hr)*Xtop(i,k) + real(j-1, wp)*hr*Xbottom(i,k)
!!$
!!$             W = (1.0_wp-real(k-1, wp)*hs)*Xfront(i,j) + real(k-1, wp)*hs*Xback(i,j)
!!$             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$
!!$             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$             UW = (1.0_wp-real(i-1, wp)*hq)*((1.0_wp-real(j-1, wp)*hr)*Xleft(1,k) &
!!$                  + real(j-1, wp)*hr*Xleft(nr,k))&
!!$                  + real(i-1, wp)*hq*((1.0_wp-real(j-1, wp)*hr)*Xright(1, k) &
!!$                  + real(j-1, wp)*hr*Xright(nr,k))
!!$
!!$
!!$             UV = (1.0_wp-real(j-1, wp)*hr)*((1.0_wp-real(k-1, wp)*hs)*Xtop(i,1) &
!!$                  + real(k-1, wp)*hs*Xtop(i,ns))&
!!$                  + real(j-1, wp)*hr*((1.0_wp-real(k-1, wp)*hs)*Xbottom(i,1) &
!!$                  + real(k-1, wp)*hs*Xbottom(i,ns))
!!$
!!$             VW = (1.0_wp-real(k-1, wp)*hs)*((1.0_wp-real(i-1, wp)*hq)*Xfront(1,j) &
!!$                  + real(i-1, wp)*hq*Xfront(nq,j))&
!!$                  + real(k-1, wp)*hs*((1.0_wp-real(i-1, wp)*hq)*Xback(1,j) &
!!$                  + real(i-1, wp)*hq*Xback(nq,j))
!!$             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$
!!$             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$             UVW1 = (1.0_wp-real(i-1, wp)*hq)*((1.0_wp-real(j-1, wp)*hr)*((1.0_wp-real(k-1, wp)*hs)*Xtop(1,1) &
!!$                  + real(k-1, wp)*hs*Xtop(1,ns))&
!!$                  + real(j-1, wp)*hr*((1.0_wp-real(k-1, wp)*hs)*Xbottom(1,1) &
!!$                  + real(k-1, wp)*hs*Xbottom(1,ns)))
!!$
!!$             UVW2 = real(i-1, wp)*hq*((1.0_wp-real(j-1, wp)*hr)*((1.0_wp-real(k-1, wp)*hs)*Xtop(nq,1) &
!!$                  + real(k-1, wp)*hs*Xtop(nq,ns))&
!!$                  + real(j-1, wp)*hr*((1.0_wp-real(k-1, wp)*hs)*Xbottom(nq,1) &
!!$                  + real(k-1, wp)*hs*Xbottom(nq,ns)))
!!$             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             G%x(i, j, k,n) = U + V + W - UV - UW - VW + UVW1 + UVW2



          end do
       end do
    end do

  end subroutine Gen3_Curve_Mesh

  subroutine Gen_Surf_Mesh(mx,my,XL,XR,XB,XT,YL,YR,YB,YT,X,Y)

    implicit none

    integer,intent(in) :: mx,my
    real(kind = wp),dimension(:),intent(in) :: XL,XR,YL,YR
    real(kind = wp),dimension(:),intent(in) :: XB,XT,YB,YT
    real(kind = wp),dimension(:,:),intent(out) :: X,Y
    real(kind = wp) :: hq,hr,q,r

    integer :: i,j



    hq = 1.0_wp/real(mx-1, wp);
    hr = 1.0_wp/real(my-1, wp);


    do j = 1,my
       do i = 1,mx

          q = real(i-1, wp)*hq
          r = real(j-1, wp)*hr


          X(i,j) = (1.0_wp-q)*XL(j)+q*XR(j)+(1.0_wp-r)*XB(i)+r*XT(i)-&
               (1.0_wp-q)*(1.0_wp-r)*XL(1)-q*(1.0_wp-r)*XR(1)-r*(1.0_wp-q)*XT(1)-&
               (r*q)*XT(mx)

          Y(i,j) = (1.0_wp-q)*YL(j)+q*YR(j)+(1.0_wp-r)*YB(i)+r*YT(i)-&
               (1.0_wp-q)*(1.0_wp-r)*YL(1)-q*(1.0_wp-r)*YR(1)-r*(1.0_wp-q)*YT(1)-&
               (r*q)*YT(mx);

       end do
    end do


  end subroutine Gen_Surf_Mesh

  subroutine read_2darray(X, m,n,filename)
    ! A subroutine to input a 2d array  from a text file.
    ! The number of rows (m),
    ! the number of columns (n) and
    ! the name of the file are input arguments

    integer, intent(in) :: m,n                  ! number of rows and column
    character(*), intent(in) :: filename        ! name of the file
    real(kind = wp), dimension(m,n),intent(inout) :: X     ! 2d array (out)
    integer :: i


    ! open the file passed in as the string "filename" on unit one
    open(unit=2, file=filename, status="old",action="read")

    do i=1,m                                 ! do for each row

       read(unit=2, fmt=*) X(i,:)                     ! read in row at a time

    end do


    close(unit=2) !***** Close the file

  end subroutine read_2darray


  subroutine reshape_2darray(Y, X, m,n)
    ! A subroutine to reshape a 1d array to a 2d array colunm-wise.
    ! The number of rows (m),
    ! the number of columns (n) and
    ! the name of the file are input arguments
    integer, intent(in) :: m,n                  ! number of rows and column
    real(kind = wp), dimension(m,n),intent(out) :: Y       ! 2d array (out)
    real(kind = wp), dimension(m*n),intent(in) :: X        ! 1d array (in)
    integer :: j


    do j=1,n                                  ! do for each column

       Y(:,j) = X(1+(j-1)*m:j*m)

    end do

  end subroutine reshape_2darray


  subroutine reshape_2darray_topo(Y, X, m,n)
    ! A subroutine to reshape a 1d array to a 2d array colunm-wise.
    ! The number of rows (m),
    ! the number of columns (n) and
    ! the name of the file are input arguments
    integer, intent(in) :: m,n                  ! number of rows and column
    real(kind = wp), dimension(m,n),intent(out) :: Y       ! 2d array (out)
    real(kind = wp), dimension(m*n),intent(in) :: X        ! 1d array (in)
    integer :: j


    do j=1,n                                  ! do for each column

       !Y(j, :) = X(1+(j-1)*m:j*m)
       Y(:,j) = X(1+(j-1)*m:j*m)

    end do

  end subroutine reshape_2darray_topo

  subroutine interpolate_fault_surface(X1, Y1, Z1, y, z, my, mz, gg)

    real(kind = wp), dimension(:,:),intent(in) :: X1, Y1, Z1       ! data
    real(kind = wp),intent(in) :: y, z
    real(kind = wp),intent(inout) :: gg
    integer,intent(in) :: my, mz

    real(kind = wp) :: dz, dy, r, r0

    integer :: jj, kk, nmy, npy, nmz, npz, stenc, ndp
    integer :: j, k

    dy = Y1(2, 1) - Y1(1,1)
    dz = Z1(1, 2) - Z1(1,1)

    r0 = sqrt(dy**2.0_wp + dz**2.0_wp)

    stenc = 1
    ndp = 3


    do kk = 1,mz
       do jj = 1,my

          r = sqrt((y-Y1(jj,kk))**2.0_wp + (z-Z1(jj,kk))**2.0_wp)

          if (r .le. 2.0_wp*r0) then

             if (abs(y-Y1(jj,kk)) .le. dy) then

                nmy = jj-stenc
                npy = jj+stenc


                if (jj == 1) then

                   nmy = jj
                   npy = jj+ndp-1

                end if

                if (jj == 2) then

                   nmy = jj-1
                   npy = jj+ndp-2

                end if

                if (jj == 3) then

                   nmy = jj-2
                   npy = jj+ndp-3

                end if

                if (jj == my-2) then

                   nmy = jj-(ndp-3)
                   npy = jj+2

                end if

                if (jj == my-1) then

                   nmy = jj-(ndp-2)
                   npy = jj+1

                end if

                if (jj == my) then

                   nmy = jj-(ndp-1)
                   npy = jj

                end if


                if  (abs(z-Z1(jj,kk)) .le. dz) then

                   nmz = kk-stenc
                   npz = kk+stenc

                   if (kk == 1) then

                      nmz = kk
                      npz = kk+ndp-1

                   end if

                   if (kk == 2) then

                      nmz = kk-1
                      npz = kk+ndp-2

                   end if

                   if (kk == 3) then

                      nmz = kk-2
                      npz = kk+ndp-3

                   end if

                   if (kk == mz) then

                      nmz = kk-(ndp-1)
                      npz = kk

                   end if

                   if (kk == mz-1) then

                      nmz = kk-(ndp-2)
                      npz = kk+1

                   end if

                   if (kk == mz-2) then

                      nmz = kk-(ndp-3)
                      npz = kk+2

                   end if


                   call interpol2d_dG(nmy, npy, nmz, npz, y, z, Y1(:,kk), Z1(jj, :), X1(:, :), gg)


                end if
                !==============================================
             end if
          end if

          !print *,  jj, kk, Y1(jj,kk), Z1(jj, kk)
       end do
    end do

  end subroutine interpolate_fault_surface



  subroutine interpol2d_dG(nmx, npx, nmy, npy, x, y, xi, yj, f, g)

    !contruct the langrange interpolation g(x,y) of  f(x,y)

    ! g(x,y): lagrange interpolant
    ! nmx, nmy: lower bound
    ! npx, npy: upper bound
    ! xi,yj: data points
    ! x,y : interpolation point
    integer, intent(in) :: nmx, nmy, npx, npy             ! number of rows and column
    real(kind = wp), intent(out) :: g                     ! lagrange interpolant 
    real(kind = wp), intent(in) :: x, y                   ! interpolation point
    real(kind = wp), dimension(:),intent(in) :: xi, yj    ! data points
    real(kind = wp), dimension(:,:),intent(in) :: f       ! data

    integer :: i, j
    real(kind = wp) :: a_i, a_j

    g = 0.0_wp

    do j = nmy, npy

       do i = nmx, npx

          a_i = lagrange(nmx, npx, i, x, xi);
          a_j = lagrange(nmy, npy, j, y, yj);

          g = g + a_i*a_j*f(i,j);

       end do
    end do

  end subroutine interpol2d_dG


  function lagrange(m, p, i, x, xi)  result(h)

    !Function to calculate  Lagrange polynomial for order N and polynomial
    ![nm, np] at location x.

    ! nm: lower bound
    ! np: upper bound
    ! xi: data points
    ! x : interpolation point

    integer, intent(in) :: m, p, i
    real(kind = wp), intent(in) :: x, xi(:)
    integer :: j
    real(kind = wp) :: h, num, den

    h = 1.0_wp
    do j = m,p
       if (j /= i) then

          num = x - xi(j)
          den = xi(i) - xi(j)
          h = h*num/den

       end if
    end do

  end function lagrange


  subroutine stretch_hyp(r, d, z)

    real(kind = wp), intent(in) :: r, d
    real(kind = wp), intent(out) :: z
    !
    z = 1.0_wp + tanh(d*(r-1.0_wp)/2.0_wp)/(tanh(d/2.0_wp))

  end subroutine stretch_hyp

end module grid
