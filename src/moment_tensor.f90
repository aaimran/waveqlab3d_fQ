module moment_tensor

  use common, only : wp
  implicit none

contains

  subroutine init_moment_tensor(G,S,infile)

    use datatypes, only : moment_tensor,block_grid_t
    use seismogram, only : Find_Coordinates_moment
    use mpi

    implicit none

    type(moment_tensor), intent(inout) :: S
    type(block_grid_t), intent(in) :: G
    integer, intent(in) :: infile
    integer :: n_mom,stat,n,p,AllocateStatus
    integer :: mq,mr,ms,pq,pr,ps
    character(256) :: temp


    character(64), dimension(:), allocatable :: source_type
    real(kind = wp), dimension(:), allocatable :: mXX,mXY,mXZ,mYY,mYZ,mZZ
    real(kind = wp), dimension(:), allocatable :: location_x,location_y, location_z
    integer, dimension(:), allocatable :: near_x, near_y, near_z, alpha
    real(kind = wp), dimension(:), allocatable :: near_phys_x, near_phys_y, near_phys_z
    real(kind = wp), dimension(:), allocatable :: duration,t_init 

    integer :: ierr,my_id,csize,status,near_temp,tag,dummy
    integer :: rootX,rootY,rootZ,rootPX,rootPY,rootPZ
    real(kind = wp) :: physX,physY,physZ,physPX,physPY,physPZ


    if (S%block_id == 1) then
       rewind(infile)
       do
          read(infile,'(a)') temp
          if (temp=='!---begin:tensor_listU---') exit
       end do

       ! determine number of moment tensors in list

       n_mom = 0
       do
          read(infile,'(a)') temp
          if (temp=='!---end:tensor_listU---') exit
          n_mom = n_mom+1
       end do
    end if

    if (S%block_id == 2) then
       rewind(infile)
       do
          read(infile,'(a)') temp
          if (temp=='!---begin:tensor_listV---') exit
       end do

       ! determine number of moment tensors in list

       n_mom = 0
       do
          read(infile,'(a)') temp
          if (temp=='!---end:tensor_listV---') exit
          n_mom = n_mom+1
       end do
    end if


    ! Set default values for moment tensor
    allocate(source_type(n_mom))
    allocate(duration(n_mom), t_init(n_mom))
    allocate(mXX(n_mom), mXY(n_mom), mXZ(n_mom),   &
         mYY(n_mom) ,mYZ(n_mom), mZZ(n_mom))
    allocate(location_x(n_mom), location_y(n_mom), location_z(n_mom)) 
    allocate(near_x(n_mom), near_y(n_mom), near_z(n_mom), alpha(n_mom))
    allocate(near_phys_x(n_mom), near_phys_y(n_mom), near_phys_z(n_mom)) 

    source_type = 'gaussian'
    duration = 0.0_wp
    location_x = 0.0_wp
    location_y = 0.0_wp
    location_z = 0.0_wp
    near_x = 0.0_wp
    near_y = 0.0_wp
    near_z = 0.0_wp
    alpha = 0
    near_phys_x = 0.0_wp
    near_phys_y = 0.0_wp
    near_phys_z = 0.0_wp
    t_init = 0.0_wp
    mXX = 0.0_wp
    mXY = 0.0_wp
    mXZ = 0.0_wp
    mYY = 0.0_wp
    mYZ = 0.0_wp
    mZZ = 0.0_wp 


    allocate(S%source_type(n_mom))
    allocate(S%duration(n_mom), S%t_init(n_mom))
    allocate(S%mXX(n_mom), S%mXY(n_mom), S%mXZ(n_mom),   &
         S%mYY(n_mom) ,S%mYZ(n_mom), S%mZZ(n_mom))
    allocate(S%location_x(n_mom), S%location_y(n_mom), S%location_z(n_mom)) 
    allocate(S%near_x(n_mom), S%near_y(n_mom), S%near_z(n_mom), S%alpha(n_mom)) 
    allocate(S%near_phys_x(n_mom), S%near_phys_y(n_mom), S%near_phys_z(n_mom)) 

    S%num_tensor = n_mom

    rewind(infile)

    if(S%block_id == 1) then
       do
          read(infile,'(a)') temp
          if (temp=='!---begin:tensor_listU---') exit
       end do
    end if

    if(S%block_id == 2) then
       do
          read(infile,'(a)') temp
          if (temp=='!---begin:tensor_listV---') exit
       end do
    end if

    if (S%num_tensor > 0) then
       do n = 1,S%num_tensor
          read(infile,*) S%source_type(n),S%duration(n),S%t_init(n),      &
               S%mXX(n),S%mYY(n),S%mZZ(n),                      &
               S%mXY(n),S%mXZ(n),S%mYZ(n),                      &
               S%location_x(n),S%location_y(n),S%location_z(n), &
               S%alpha(n)
       end do
    end if

    mq = G%C%mq 
    mr = G%C%mr 
    ms = G%C%ms 
    pq = G%C%pq 
    pr = G%C%pr 
    ps = G%C%ps 


    ! allocate (S%exact(mq:pq,mr:pr,ms:ps,3), stat = AllocateStatus)
    ! if (AllocateStatus /= 0) stop "*** Not enough memory for exact solution ***" 
    ! if(allocated(S%exact(mq:pq,mr:pr,ms:ps,3))) print*,'ALLOCATED!!!!!!!, block ',S%block_id,&
    !    ' bounds: ',mq,pq,mr,pr,ms,ps

    ! print*,'blockID: ',S%block_id
    ! print*,'init: ',mq,mr,ms,pq,pr,ps


    if (S%num_tensor > 0) then

       call Find_Coordinates_moment(G%X,S%location_x,S%location_y,S%location_z,  &             
            S%near_x,S%near_y,S%near_z,S%num_tensor,G%C%nq,G%C%nr,G%C%ns,mq,mr,ms,pq,pr,ps)          

       call MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD,my_id,ierr) 

       do n = 1,S%num_tensor
          rootX  = -1
          rootY  = -1  
          rootZ  = -1  
          rootPX = -1  
          rootPY = -1  
          rootPZ = -1  
          physX  = -1_wp
          physY  = -1_wp
          physZ  = -1_wp
          physPX = -1_wp
          physPY = -1_wp
          physPZ = -1_wp

          ! Find_Coordinates_moment already resolved whether this rank owns the
          ! nearest physical node.  Rechecking floating-point coordinates at a
          ! partition boundary can reject a source because of roundoff (for
          ! example x=0 versus a computed local minimum slightly above zero).
          if (S%near_x(n) >= mq .and. S%near_x(n) <= pq .and. &
              S%near_y(n) >= mr .and. S%near_y(n) <= pr .and. &
              S%near_z(n) >= ms .and. S%near_z(n) <= ps) then

             rootX = S%near_x(n)
             rootY = S%near_y(n)
             rootZ = S%near_z(n)

             physX = G%X(S%near_x(n),S%near_y(n),S%near_z(n),1)
             physY = G%X(S%near_x(n),S%near_y(n),S%near_z(n),2)
             physZ = G%X(S%near_x(n),S%near_y(n),S%near_z(n),3) 

          end if


          call MPI_Allreduce(rootX,rootPX,1,&
               MPI_INTEGER,MPI_MAX,G%C%comm,ierr)
          call MPI_Allreduce(rootY,rootPY,1,&
               MPI_INTEGER,MPI_MAX,G%C%comm,ierr)
          call MPI_Allreduce(rootZ,rootPZ,1,&
               MPI_INTEGER,MPI_MAX,G%C%comm,ierr)

          call MPI_Allreduce(physX,physPX,1,&
               MPI_DOUBLE,MPI_MAX,G%C%comm,ierr)
          call MPI_Allreduce(physY,physPY,1,&
               MPI_DOUBLE,MPI_MAX,G%C%comm,ierr)
          call MPI_Allreduce(physZ,physPZ,1,&
               MPI_DOUBLE,MPI_MAX,G%C%comm,ierr)

          S%near_x(n) = rootPX
          S%near_y(n) = rootPY
          S%near_z(n) = rootPZ

          S%near_phys_x(n) = physPX
          S%near_phys_y(n) = physPY
          S%near_phys_z(n) = physPZ


          ! print*,'id: ',my_id,'near: ',S%near_x(n),S%near_y(n),S%near_z(n)
          ! print*,'id: ',my_id,'near: ',S%near_phys_x(n),S%near_phys_y(n),S%near_phys_z(n)

       end do

       !    do n = 1,S%num_tensor
       !      if (S%near_x(n) > 0 .and. S%near_y(n) > 0 .and. S%near_z(n) > 0) then
       !        S%near_phys_x(n) = G%X(S%near_x(n),S%near_y(n),S%near_z(n),1)
       !        S%near_phys_y(n) = G%X(S%near_x(n),S%near_y(n),S%near_z(n),2)
       !        S%near_phys_z(n) = G%X(S%near_x(n),S%near_y(n),S%near_z(n),3)
       !      else
       !        S%near_phys_x(n) = -9999_wp
       !      end if
       !    end do

    end if

    !print*, 'block id: ', S%block_id
    !print*, 'source type: ', S%source_type
    !print*, 'location: ', S%location_x,S%location_y,S%location_z
    !print*, 'nearest indices: ',  S%near_x,S%near_y,S%near_z
    !print*, 'nearest physical: ',  S%near_phys_x,S%near_phys_y,S%near_phys_z


  end subroutine init_moment_tensor


  !   subroutine exact_moment_tensor(receiver,source,time)

  !       use datatypes, only : moment_tensor, block_type
  !       use mpi

  !       implicit none

  !       type(block_type), intent(in) :: receiver
  !       type(moment_tensor), intent(inout) :: source
  !       real(kind = wp), intent(in) :: time
  !       integer :: mq,mr,ms,pq,pr,ps,l,m,n,i,AllocateStatus
  !       real(kind = wp) :: m0=1_wp, r, time2 , gamX,gamY,gamZ
  !       real(kind = wp) :: A,B,C,D,E,alpha,beta,rho
  !       real(kind = wp),parameter :: pi = 3.141592653589793_wp
  !       integer :: ierr, rank,num_procs

  !       mq = receiver%G%C%mq 
  !       mr = receiver%G%C%mr 
  !       ms = receiver%G%C%ms 
  !       pq = receiver%G%C%pq 
  !       pr = receiver%G%C%pr 
  !       ps = receiver%G%C%ps 

  !   if(allocated(source%exact(mq:pq,mr:pr,ms:ps,:))) then

  !       alpha = receiver%rho_s_p(3)
  !       beta = receiver%rho_s_p(2)
  !       rho = receiver%rho_s_p(1)

  !      do i = 1,source%num_tensor   
  !        do l = mq,pq
  !          do m = mr,pr
  !            do n = ms,ps

  !           time2 = time - source%t_init(i)

  !          ! Compute source-to-receiver distance

  !            r = sqrt((source%location_x(i) - receiver%G%x(l,m,n,1))**2 + &
  !                     (source%location_y(i) - receiver%G%x(l,m,n,2))**2 + &
  !                     (source%location_z(i) - receiver%G%x(l,m,n,3))**2 )

  !            select case(source%source_type(i))

  !            case default

  !               stop 'invalid moment tensor source type!'

  !            case('gaussian')

  !                call Gaussian(A,time2,r,alpha,source%duration(i))
  !                call Gaussian(B,time2,r,beta,source%duration(i))
  !                call integral_Gaussian(C,time2,r,source%duration(i),alpha,beta)
  !                call d_dt_Gaussian(D,time2,r,alpha,source%duration(i))
  !                call d_dt_Gaussian(E,time2,r,beta,source%duration(i))

  !            end select  


  !     ! Compute unit vectors between source and receiver

  !              if (r /= 0.0_wp) then 
  !                gamX = (receiver%G%x(l,m,n,1) - source%location_x(i))/r
  !                gamY = (receiver%G%x(l,m,n,2) - source%location_y(i))/r
  !                gamZ = (receiver%G%x(l,m,n,3) - source%location_z(i))/r
  !              else
  !                source%exact(l,m,n,:) = 0.0_wp
  !                cycle
  !              end if


  !     ! ----- Displacement component X ----- !

  !          ! Mxx * Gxx,x
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mXX(i)/(4_wp*pi*rho)*    &
  !              
  !              ( ((15_wp*gamX*gamX*gamX - 9_wp*gamX) / r**4_wp) * C                           &

  !              + ((6_wp*gamX*gamX*gamX  - 3_wp*gamX) / (alpha**2_wp * r**2_wp)) * A           &

  !              - ((6_wp*gamX*gamX*gamX  - 4_wp*gamX) / (beta**2_wp  * r**2_wp)) * B           &

  !              + ((  gamX*gamX*gamX             ) / (alpha**3_wp * r   )) * D                    &

  !              - ((( gamX*gamX   - 1_wp ) * gamX) / (beta**3_wp  * r   )) * E                 &

  !              )


  !           ! Myy * Gxy,y
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mYY(i)/(4_wp*pi*rho)*    &
  !              
  !              ( ((15_wp*gamX*gamY*gamY - 3_wp*gamX) / r**4_wp) * C                           &

  !              + ((6_wp*gamX*gamY*gamY  -   gamX) / (alpha**2_wp * r**2_wp)) * A              &

  !              - ((6_wp*gamX*gamY*gamY  -   gamX) / (beta**2_wp  * r**2_wp)) * B              &

  !              + ((  gamX*gamY*gamY          ) / (alpha**3_wp * r   )) * D                    &

  !              - ((  gamX*gamY*gamY          ) / (beta**3_wp  * r   )) * E                    &

  !              )  

  !            ! Mzz * Gxz,z
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mZZ(i)/(4_wp*pi*rho)*    &
  !              
  !              ( ((15_wp*gamX*gamZ*gamZ - 3_wp*gamX) / r**4_wp) * C                           &

  !              + ((6_wp*gamX*gamZ*gamZ  -   gamX) / (alpha**2_wp * r**2_wp)) * A              &

  !              - ((6_wp*gamX*gamZ*gamZ  -   gamX) / (beta**2_wp  * r**2_wp)) * B              &

  !              + ((  gamX*gamZ*gamZ          ) / (alpha**3_wp * r   )) * D                    &

  !              - ((  gamX*gamZ*gamZ          ) / (beta**3_wp  * r   )) * E                    &

  !              ) 

  !             ! Mxy * Gxx,y
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mXY(i)/(4_wp*pi*rho)*    &
  !              
  !              ( ((15_wp*gamX*gamX*gamY - 3_wp*gamY) / r**4_wp) * C                           &

  !              + ((6_wp*gamX*gamX*gamY  -   gamY) / (alpha**2_wp * r**2_wp)) * A              &

  !              - ((6_wp*gamX*gamX*gamY  - 2_wp*gamY) / (beta**2_wp  * r**2_wp)) * B           &

  !              + ((  gamX*gamX*gamY          ) / (alpha**3_wp * r   )) * D                    &

  !              - ((( gamX*gamX   - 1_wp ) * gamY) / (beta**3_wp  * r   )) * E                 &

  !              ) 

  !              ! Mxz * Gxx,z
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mXZ(i)/(4_wp*pi*rho)*    &
  !              
  !              ( ((15_wp*gamX*gamX*gamZ - 3_wp*gamZ) / r**4_wp) * C                           &

  !              + ((6_wp*gamX*gamX*gamZ  -   gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamX*gamX*gamZ  - 2_wp*gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamX*gamX*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamX*gamX   - 1_wp ) * gamZ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Myx * Gxy,x
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mXY(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamX*gamY*gamX - 3_wp*gamY) / r**4_wp) * C                            &

  !              + ((6_wp*gamX*gamY*gamX  -   gamY) / (alpha**2_wp * r**2_wp)) * A               &

  !              - ((6_wp*gamX*gamY*gamX  -   gamY) / (beta**2_wp  * r**2_wp)) * B               &

  !              + ((  gamX*gamY*gamX          ) / (alpha**3_wp * r   )) * D                     &

  !              - ((  gamX*gamY*gamX          ) / (beta**3_wp  * r   )) * E                     &

  !              )

  !              ! Mzx * Gxz,x
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mXZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamX*gamZ*gamX - 3_wp*gamZ) / r**4_wp) * C               &

  !              + ((6_wp*gamX*gamZ*gamX  -   gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamX*gamZ*gamX  -   gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamX*gamZ*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamX*gamZ*gamX          ) / (beta**3_wp  * r   )) * E  &

  !              )

  !              ! Myz * Gxy,z
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mYZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamX*gamY*gamZ         ) / r**4_wp) * C               &

  !              + ((6_wp*gamX*gamY*gamZ          ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamX*gamY*gamZ          ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamX*gamY*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamX*gamY*gamZ          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mzy * Gxz,y
  !          source%exact(l,m,n,1) = source%exact(l,m,n,1) + m0*source%mYZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamX*gamZ*gamY         ) / r**4_wp) * C               &

  !              + ((6_wp*gamX*gamZ*gamY          ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamX*gamZ*gamY          ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamX*gamZ*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamX*gamZ*gamY          ) / (beta**3_wp  * r   )) * E  &

  !              )  



  !     ! ----- Displacement component Y ----- !

  !          ! Mxx * Gyx,x
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mXX(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamX*gamX - 3_wp*gamY) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamX*gamX  -   gamY) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamX*gamX  -   gamY) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamX*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamY*gamX*gamX          ) / (beta**3_wp  * r   )) * E  &

  !              )
  !                 
  !          ! Myy * Gyy,y
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mYY(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamY*gamY - 9_wp*gamY) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamY*gamY  - 3_wp*gamY) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamY*gamY  - 4_wp*gamY) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamY*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamY*gamY   - 1_wp ) * gamY) / (beta**3_wp  * r   )) * E  &

  !              )  

  !            ! Mzz * Gyz,z
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mZZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamZ*gamZ - 3_wp*gamY) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamZ*gamZ  -   gamY) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamZ*gamZ  -   gamY) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamZ*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamY*gamZ*gamZ          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !             ! Mxz * Gyx,z
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mXZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamX*gamZ         ) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamX*gamZ          ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamX*gamZ          ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamX*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamY*gamX*gamZ          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mxy * Gyx,y
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mXY(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamX*gamY - 3_wp*gamX) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamX*gamY  -   gamX) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamX*gamY  -   gamX) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamX*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamY*gamX*gamY         )) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Myx * Gyy,x
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mXY(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamY*gamX - 3_wp*gamX) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamY*gamX  -   gamX) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamY*gamX  - 2_wp*gamX) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamY*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamY*gamY   - 1_wp ) * gamX) / (beta**3_wp  * r   )) * E  &

  !              )

  !              ! Myz * Gyy,z
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mYZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamY*gamZ - 3_wp*gamZ) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamY*gamZ  -   gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamY*gamZ  - 2_wp*gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamY*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamY*gamY   - 1_wp ) * gamZ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mzx * Gyz,x
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mXZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamZ*gamX         ) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamZ*gamX          ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamZ*gamX          ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamZ*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamY*gamZ*gamX          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mzy * Gyz,y
  !          source%exact(l,m,n,2) = source%exact(l,m,n,2) + m0*source%mYZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamY*gamZ*gamY - 3_wp*gamZ) / r**4_wp) * C               &

  !              + ((6_wp*gamY*gamZ*gamY  -   gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamY*gamZ*gamY  -   gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamY*gamZ*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamY*gamZ*gamY          ) / (beta**3_wp  * r   )) * E  &

  !              ) 




  !     ! ----- Displacement component Z ----- !

  !          ! Mzz * Gzz,z
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mZZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamZ*gamZ - 9_wp*gamZ) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamZ*gamZ  - 3_wp*gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamZ*gamZ  - 4_wp*gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamZ*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamZ*gamZ   - 1_wp ) * gamZ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mzy * Gzz,y
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mYZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamZ*gamY - 3_wp*gamY) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamZ*gamY  -   gamY) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamZ*gamY  - 2_wp*gamY) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamZ*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamZ*gamZ   - 1_wp ) * gamY) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mzx * Gzz,x
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mXZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamZ*gamX - 3_wp*gamX) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamZ*gamX  -   gamX) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamZ*gamX  - 2_wp*gamX) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamZ*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamZ*gamZ   - 1_wp ) * gamX) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Myz * Gzy,z
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mYZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamY*gamZ - 3_wp*gamY) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamY*gamZ  -   gamY) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamY*gamZ  -   gamY) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamY*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamZ*gamY*gamZ          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !            ! Myy * Gzy,y
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mYY(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamY*gamY - 3_wp*gamZ) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamY*gamY  -   gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamY*gamY  -   gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamY*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamZ*gamY*gamY          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Myx * Gzy,x
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mXY(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamY*gamX         ) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamY*gamX          ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamY*gamX          ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamY*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamZ*gamY*gamX          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !              ! Mxz * Gzx,z
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mXZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamX*gamZ - 3_wp*gamX) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamX*gamZ  -   gamX) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamX*gamZ  -   gamX) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamX*gamZ          ) / (alpha**3_wp * r   )) * D  &

  !              - ((( gamZ*gamX*gamZ         )) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !             ! Mxy * Gzx,y
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mXZ(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamX*gamY         ) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamX*gamY          ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamX*gamY          ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamX*gamY          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamZ*gamX*gamY          ) / (beta**3_wp  * r   )) * E  &

  !              ) 

  !          ! Mxx * Gzx,x
  !          source%exact(l,m,n,3) = source%exact(l,m,n,3) + m0*source%mXX(i)/(4_wp*pi*rho)*     &
  !              
  !              ( ((15_wp*gamZ*gamX*gamX - 3_wp*gamZ) / r**4_wp) * C               &

  !              + ((6_wp*gamZ*gamX*gamX  -   gamZ) / (alpha**2_wp * r**2_wp)) * A  &

  !              - ((6_wp*gamZ*gamX*gamX  -   gamZ) / (beta**2_wp  * r**2_wp)) * B  &

  !              + ((  gamZ*gamX*gamX          ) / (alpha**3_wp * r   )) * D  &

  !              - ((  gamZ*gamX*gamX          ) / (beta**3_wp  * r   )) * E  &

  !              ) 
  !           ! print*,source%exact(l,m,n,1), source%exact(l,m,n,2), source%exact(l,m,n,3)
  !            end do
  !          end do
  !        end do   
  !      end do
  !    end if


  !   end subroutine exact_moment_tensor


  ! This subroutine will compute a Gaussian of period 'f'

  subroutine Gaussian(A,t,r,c,w)

    implicit none

    real(kind = wp), intent(in) :: t,r,c,w
    real(kind = wp), intent(inout) :: A
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp) :: f

    !  f = 1/w

    !  A = (1.0_wp / (f * sqrt(2.0_wp*pi))) * &

    !           exp(-((t-r/c)**2.0_wp)/(2.0_wp*f*f))

    A = (w/sqrt(2*pi))*exp((-w**2*(t)**2_wp)/2_wp)

  end subroutine Gaussian


  ! This subroutine will compute the derivative of a Gaussian of period 'f'

  subroutine d_dt_Gaussian(A,t,r,c,w)

    implicit none

    real(kind = wp), intent(in) :: t,r,c,w
    real(kind = wp), intent(inout) :: A
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp) :: f

    !f = 1/w

    !A = (1.0_wp / (f * sqrt(2.0_wp*pi))) * &

    !     (-exp(-((t-r/c)**2.0_wp)/(2.0_wp*f*f)) * (t-r/c)) / (2.0_wp*f*f)

    A = -((w**2_wp)*(t)/(sqrt(2_wp*pi)))*exp((-w**2*(t)**2_wp)/2_wp)

  end subroutine d_dt_Gaussian


  ! This subroutine will compute the integral of a Gaussian of
  ! period 'f' over the time interval between r/alpha and r/beta

  subroutine integral_Gaussian(A,t,r,w,alpha,beta)

    implicit none

    real(kind = wp), intent(in) :: t,r,w,alpha,beta
    real(kind = wp), intent(inout) :: A
    real(kind = wp),parameter :: pi = 3.141592653589793_wp
    real(kind = wp) :: f

    f = 1/w

    A = -0.5_wp*t*(erf((t-r/beta)/(sqrt(2.0_wp)*f))        &
         
         - erf((t-r/alpha)/(sqrt(2.0_wp)*f)))               &
         
         - f/sqrt(2.0_wp*pi)*(exp(-(t-r/beta)**2.0_wp/(2.0_wp*f*f))   &  
         
         - exp(-(t-r/alpha)**2.0_wp/(2.0_wp*f*f)))


  end subroutine integral_Gaussian


  subroutine source_time(A,source_type,t,t0,w)

    implicit none

    character(*), intent(in) :: source_type
    real(kind = wp), intent(out) :: A
    real(kind = wp), intent(in) :: t,t0
    real(kind = wp), intent(in), optional :: w
    real(kind = wp),parameter :: pi = 3.141592653589793_wp

    select case (source_type)

    case default

       stop 'invalid source type'

    case('gaussian')

       ! Source time function, derivative of Gaussian

       ! For velocities
       A = (w/(sqrt(2_wp*pi)))*exp((-w**2*(t-t0)**2)/2_wp)

       ! For accelerations
       ! A = ((-w**3*(t-t0))/(sqrt(2_wp*pi)))*exp((-w**2*(t-t0)**2)/2_wp)

    case('LOH_discontinuity')

       A = (t/w**2) * exp(-t/w)

    case('Brune')

       A = (w**2) * t * exp(-w*t)

    case('Cos-Bell')

       if (t < t0) then
           A = (1.0_wp-cos(2.0_wp*pi*t/t0))/t0
       else
           A = 0.0_wp  
       end if 

    case('cosine')

       A = 0.0_wp
       if((t .ge. 0.0_wp) .and. (t .le. w)) A = 1.0_wp-cos(2.0_wp*pi*t/w) 


    case('LOH_moment_time_history')

       A = (1.0_wp-(1.0_wp+t/w)*exp(-t/w))

    case('sine_moment_time_history')

       A = w
       if(t .le. w) A = t - w/(2.0_wp*pi)*sin(2.0_wp*pi*t/w)

    end select

  end subroutine source_time


  subroutine set_moment_tensor(B,t)

    use datatypes, only : block_type
    use seismogram, only : Find_Coordinates
    use mpi

    implicit none

    type(block_type), intent(inout) :: B
    real(kind = wp), intent(in) :: t
    integer :: mq,mr,ms,pq,pr,ps,l,m,n,i,nstations
    real(kind = wp) :: A,hx,hy,hz,ax,ay,az,dx,dy,dz,w,delta
    integer :: rank,num_procs,ierr


    mq = B%G%C%mq 
    mr = B%G%C%mr 
    ms = B%G%C%ms 
    pq = B%G%C%pq 
    pr = B%G%C%pr 
    ps = B%G%C%ps 


    do i = 1,B%MT%num_tensor
       do n = ms,ps
          do m = mr,pr
             do l = mq,pq

                ! Compute local grid spacing
                hx = B%G%hq/B%G%metricx(l,m,n,1)
                hy = B%G%hr/B%G%metricy(l,m,n,2)
                hz = B%G%hs/B%G%metricz(l,m,n,3)

                ! Compute offsets between physical source and nearest grid point, normalized by grid spacing
                ax = (B%MT%location_x(i) - B%MT%near_phys_x(i)) / hx
                ay = (B%MT%location_y(i) - B%MT%near_phys_y(i)) / hy
                az = (B%MT%location_z(i) - B%MT%near_phys_z(i)) / hz

                ! Compute coefficient representing source time function at time t
                w = B%MT%duration(i)
                call source_time(A,B%MT%source_type(i),t,B%MT%t_init(i),w)

                call singular_source(l,B%MT%near_x(i),hx,ax,dx,B%MT%order,B%MT%alpha)
                call singular_source(m,B%MT%near_y(i),hy,ay,dy,B%MT%order,B%MT%alpha)
                call singular_source(n,B%MT%near_z(i),hz,az,dz,B%MT%order,B%MT%alpha)

                delta = dx*dy*dz

                B%F%DF(l,m,n,4) = B%F%DF(l,m,n,4) - A*B%MT%mXX(i)*delta 
                B%F%DF(l,m,n,5) = B%F%DF(l,m,n,5) - A*B%MT%mYY(i)*delta 
                B%F%DF(l,m,n,6) = B%F%DF(l,m,n,6) - A*B%MT%mZZ(i)*delta 
                B%F%DF(l,m,n,7) = B%F%DF(l,m,n,7) - A*B%MT%mXY(i)*delta 
                B%F%DF(l,m,n,8) = B%F%DF(l,m,n,8) - A*B%MT%mXZ(i)*delta 
                B%F%DF(l,m,n,9) = B%F%DF(l,m,n,9) - A*B%MT%mYZ(i)*delta 

             end do
          end do
       end do
    end do

  end subroutine set_moment_tensor


  subroutine singular_source(i,i0,h,a,d,order,alpha)

    implicit none

    integer,intent(in) :: i,i0,order
    integer,intent(in),dimension(:) :: alpha
    real(kind = wp),intent(in) :: h,a
    real(kind = wp),intent(out) :: d
    real(kind = wp) :: beta,gam,gap

    ! d = delta function

    if (order .EQ. 2) then
       ! centered five point discretization, correct for w only with a/=0

       beta = (1_wp-2_wp*a)/4_wp + 0.25_wp;
       gam = (2_wp*a-1_wp)/(16_wp);
       gap = (2_wp*a+1_wp)/(16_wp);

       if     (i==i0-2) then
          d = (-gam)/h

       elseif (i==i0-1) then
          d = (beta + 3_wp*gam - gap)/h

       elseif (i==i0  ) then
          d = (1_wp-a - 2_wp*beta -3_wp*gam + 3_wp*gap)/h

       elseif (i==i0+1) then
          d = (a + beta + gam - 3_wp*gap)/h

       elseif (i==i0+2) then
          d = (gap)/h

       else
          d = 0_wp

       end if
    end if

    if (order .EQ. 4) then
       ! centered 9-point discretization, for alpha = 0

       if (alpha(1) .EQ. 0) then

          if     (i==i0-4) then
             d = -0.01171875_wp  / h

          elseif (i==i0-3) then
             d = -0.03125_wp     / h

          elseif (i==i0-2) then
             d = 0.046875_wp     / h

          elseif (i==i0-1) then
             d = 0.28125_wp      / h

          elseif (i==i0  ) then
             d = 0.4296875_wp    / h

          elseif (i==i0+1) then
             d = 0.28125_wp      / h

          elseif (i==i0+2) then
             d = 0.046875_wp     / h

          elseif (i==i0+3) then
             d = -0.03125_wp     / h

          elseif (i==i0+4) then
             d = -0.01171875_wp  / h

          else
             d = 0_wp 

          end if
       end if

       if (alpha(1) .EQ. 25) then

          if     (i==i0-4) then
             d = -0.00246875_wp  / h

          elseif (i==i0-3) then
             d = -0.02275_wp     / h

          elseif (i==i0-2) then
             d = -0.014125_wp     / h

          elseif (i==i0-1) then
             d = 0.16575_wp      / h

          elseif (i==i0  ) then
             d = 0.4171875_wp    / h

          elseif (i==i0+1) then
             d = 0.38675_wp      / h

          elseif (i==i0+2) then
             d = 0.117875_wp     / h

          elseif (i==i0+3) then
             d = -0.02975_wp     / h

          elseif (i==i0+4) then
             d = -0.01846875_wp  / h

          else
             d = 0_wp 

          end if
       end if

       if (alpha(1) .EQ. 50) then

          if     (i==i0-4) then
             d = -0.00728125_wp  / h

          elseif (i==i0-3) then
             d = -0.028_wp     / h

          elseif (i==i0-2) then
             d = 0.014875_wp     / h

          elseif (i==i0-1) then
             d = 0.224_wp      / h

          elseif (i==i0  ) then
             d = 0.4265625_wp    / h

          elseif (i==i0+1) then
             d = 0.336_wp      / h

          elseif (i==i0+2) then
             d = 0.081375_wp     / h

          elseif (i==i0+3) then
             d = -0.032_wp     / h

          elseif (i==i0+4) then
             d = -0.01553125_wp  / h

          else
             d = 0_wp 

          end if
       end if

       if (alpha(1) .EQ. 100) then


          if     (i==i0-4) then
             d = -0.0095625_wp  / h

          elseif (i==i0-3) then
             d = -0.02990625_wp / h

          elseif (i==i0-2) then
             d = 0.03053125_wp  / h

          elseif (i==i0-1) then
             d = 0.25284375_wp  / h

          elseif (i==i0  ) then
             d = 0.42890625_wp  / h

          elseif (i==i0+1) then
             d = 0.30903125_wp  / h

          elseif (i==i0+2) then
             d = 0.06384375_wp  / h

          elseif (i==i0+3) then
             d = -0.03196875_wp / h

          elseif (i==i0+4) then
             d = -0.01371875_wp / h

          else
             d = 0_wp 

          end if
       end if

    end if

  end subroutine singular_source


  subroutine moment_tensor_body_force(B,t)

    use datatypes, only : block_type
    implicit none

    type(block_type), intent(inout) :: B
    real(kind = wp), intent(in) :: t

    ! moment tensor parameters
    real(kind = wp) :: M0,T0,sigma,x,y,z                                 
    real(kind = wp) :: r2, A, D1, D2, D3, w
    real(kind = wp), dimension(3,3) :: Munit                              ! unit moment tensor

    integer :: i,j,k,m                                           ! grid indices
    real(kind = wp),parameter :: pi = 3.14159265359_wp
    integer :: mq,mr,ms,pq,pr,ps,nx,ny,nz
    real(kind = wp) :: fx, fy, fz, hx,hy,hz

    mq = B%G%C%mq 
    mr = B%G%C%mr 
    ms = B%G%C%ms 
    pq = B%G%C%pq 
    pr = B%G%C%pr 
    ps = B%G%C%ps

    nx = B%G%C%nq
    ny = B%G%C%nr
    nz = B%G%C%ns

    ! unit moment tensor
    Munit(1,1) = 0.0_wp
    Munit(1,2) = 0.0_wp
    Munit(1,3) = 0.0_wp

    Munit(2,1) = 0.0_wp
    Munit(2,2) = 0.0_wp
    Munit(2,3) = 1.0_wp

    Munit(3,1) = 0.0_wp
    Munit(3,2) = 1.0_wp
    Munit(3,3) = 0.0_wp


    !M0 = 1000.0_wp
    !T0 = 0.1_wp
    sigma = 0.2_wp
    !x = 2.0_wp
    !y = 0.0_wp
    !z = 0.0_wp

    ! now evaluate body force at all points
    do m = 1,B%MT%num_tensor
       do k = ms,ps
          do j = mr,pr
             do i = mq,pq


                ! prescribe the moment tensor
                Munit(1,1) = B%MT%mXX(m)
                Munit(1,2) = B%MT%mXY(m)
                Munit(1,3) = B%MT%mXZ(m)

                Munit(2,1) = B%MT%mXY(m)
                Munit(2,2) = B%MT%mYY(m)
                Munit(2,3) = B%MT%mYZ(m)

                Munit(3,1) = B%MT%mXZ(m)
                Munit(3,2) = B%MT%mYZ(m)
                Munit(3,3) = B%MT%mZZ(m)


                ! extract local grid spacing
                hx = B%G%X(i,j,k,1)-B%G%X(i-1,j,k,1)
                hy = B%G%X(i,j,k,2)-B%G%X(i,j-1,k,2)
                hz = B%G%X(i,j,k,3)-B%G%X(i,j,k-1,3)

                if (i == mq) hx = B%G%X(mq+1,j,k,1)-B%G%X(mq,j,k,1)
                if (j == mr) hy = B%G%X(i,mr+1,k,2)-B%G%X(i,mr,k,2)
                if (k == ms) hz = B%G%X(i,j,ms+1,3)-B%G%X(i,ms,1,3)

                ! coordinates of the point source
                x = B%MT%location_x(m)
                y = B%MT%location_y(m)
                z = B%MT%location_z(m)


                ! weights from derivatives in different directions
                D1 = B%G%X(i,j,k,1) - x
                D2 = B%G%X(i,j,k,2) - y
                D3 = B%G%X(i,j,k,3) - z


                ! common amplitude factor
                r2 = D1**2+D2**2+D3**2

                ! Compute coefficient representing source time function at time t
                w = B%MT%duration(m)
                call source_time(A,B%MT%source_type(m),t,B%MT%t_init(m),w)

                ! A = -M0/B%M%M(i,i,k,3)*(1.0_wp-(1.0_wp+t/T0)*exp(-t/T0))* &
                !     exp(-0.5_wp*r2/sigma**2)/(sigma**5*sqrt(8.0_wp*pi**3))

                ! set the standard deviation of the Gaussian
                sigma = 2.0_wp*max(hx, hy, hz)

                A = -1.0_wp/B%M%M(i,i,k,3)*A* &
                     exp(-0.5_wp*r2/sigma**2)/(sigma**5*sqrt(8.0_wp*pi**3))

                ! compute double couple force field
                fx = A*(Munit(1,1)*D1 + Munit(1,2)*D2 + Munit(1,3)*D3)
                fy = A*(Munit(2,1)*D1 + Munit(2,2)*D2 + Munit(2,3)*D3)
                fz = A*(Munit(3,1)*D1 + Munit(3,2)*D2 + Munit(3,3)*D3)

                ! add body forces to rates in momentum balance
                B%F%DF(i,j,k,1) = B%F%DF(i,j,k,1) - fx
                B%F%DF(i,j,k,2) = B%F%DF(i,j,k,2) - fy
                B%F%DF(i,j,k,3) = B%F%DF(i,j,k,3) - fz


             end do
          end do
       end do
    end do

    !stop 
  end subroutine moment_tensor_body_force


  subroutine set_moment_tensor_smooth(B,t)

    use datatypes, only : block_type
    use seismogram, only : Find_Coordinates
    use mpi

    implicit none

    type(block_type), intent(inout) :: B
    real(kind = wp), intent(in) :: t
    integer :: mq,mr,ms,pq,pr,ps,l,m,n,i,nstations
    real(kind = wp) :: A,hx,hy,hz,ax,ay,az,dx,dy,dz,w,delta
    real(kind = wp) :: sigmax,sigmay,sigmaz, x,y,z 
    real(kind = wp),parameter :: pi = 3.14159265359_wp
    integer :: rank,num_procs,ierr


    mq = B%G%C%mq 
    mr = B%G%C%mr 
    ms = B%G%C%ms 
    pq = B%G%C%pq 
    pr = B%G%C%pr 
    ps = B%G%C%ps 


    do i = 1,B%MT%num_tensor
       do n = ms,ps
          do m = mr,pr
             do l = mq,pq

                ! Compute local grid spacing
                ! extract local grid spacing
                hx = B%G%X(l,m,n,1)-B%G%X(l-1,m,n,1)
                hy = B%G%X(l,m,n,2)-B%G%X(l,m-1,n,2)
                hz = B%G%X(l,m,n,3)-B%G%X(l,m,n-1,3)
                
                if (l == mq) hx = B%G%X(mq+1,m,n,1)-B%G%X(l,m,n,1)
                if (m == mr) hy = B%G%X(l,mr+1,n,2)-B%G%X(l,m,n,2)
                if (n == ms) hz = B%G%X(l,m,ms+1,3)-B%G%X(l,m,n,3)

                ! coordinates of the point source
                x = B%MT%location_x(i)
                y = B%MT%location_y(i)
                z = B%MT%location_z(i)

                ! standard deviation of the Gaussian
                sigmax = 2.0_wp*hx
                sigmay = 2.0_wp*hy
                sigmaz = 2.0_wp*hz
                
                ! Gaussian weights in different directions
               
                dx = exp(-0.5_wp*((B%G%X(l,m,n,1) - x)/sigmax)**2.0_wp)/(sigmax*sqrt(2.0_wp*pi))
                dy = exp(-0.5_wp*((B%G%X(l,m,n,2) - y)/sigmay)**2.0_wp)/(sigmay*sqrt(2.0_wp*pi))
                dz = exp(-0.5_wp*((B%G%X(l,m,n,3) - z)/sigmaz)**2.0_wp)/(sigmaz*sqrt(2.0_wp*pi))

                ! Compute coefficient representing source time function at time t
                w = B%MT%duration(i)
                call source_time(A,B%MT%source_type(i),t,B%MT%t_init(i),w)

                
                delta = dx*dy*dz

                B%F%DF(l,m,n,4) = B%F%DF(l,m,n,4) - A*B%MT%mXX(i)*delta 
                B%F%DF(l,m,n,5) = B%F%DF(l,m,n,5) - A*B%MT%mYY(i)*delta 
                B%F%DF(l,m,n,6) = B%F%DF(l,m,n,6) - A*B%MT%mZZ(i)*delta 
                B%F%DF(l,m,n,7) = B%F%DF(l,m,n,7) - A*B%MT%mXY(i)*delta 
                B%F%DF(l,m,n,8) = B%F%DF(l,m,n,8) - A*B%MT%mXZ(i)*delta 
                B%F%DF(l,m,n,9) = B%F%DF(l,m,n,9) - A*B%MT%mYZ(i)*delta 

             end do
          end do
       end do
    end do

  end subroutine set_moment_tensor_smooth


end module moment_tensor
