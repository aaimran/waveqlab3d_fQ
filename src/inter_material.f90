module inter_material

  use common, only : wp
  use datatypes, only : block_type, block_temp_parameters,block_grid_t
  
  implicit none

contains

subroutine interpolatematerials(B,btp)
  ! interpolate elastic material properties
  !
  type(block_type),intent(inout) :: B
  type(block_temp_parameters), intent(in) :: btp
  real(kind = wp), dimension(:), allocatable :: xi, yj, zk
  real(kind = wp):: hx, hy, hz
  integer :: ii, i,j,k, nx, ny, nz
  real(kind = wp), dimension(:,:,:), allocatable :: FF
  integer :: i_m,j_m,k_m,i_p,j_p,k_p
  real(kind = wp) ::  x, y, z
  integer :: i0,j0,k0

  real(kind = wp) :: a_i,a_j,a_k
  real(kind = wp) :: F_ijk
  
  integer :: mx,px, my,py, mz,pz
  integer :: ix, iy, iz, fx, fy, fz
  integer :: ndp, stenc

  
  allocate(FF(B%G%C%mq:B%G%C%pq, B%G%C%mr:B%G%C%pr, B%G%C%ms:B%G%C%ps))
  
  FF(:,:,:) = 0.0_wp

  nx = btp%nqrs(1)
  ny = btp%nqrs(2)
  nz = btp%nqrs(3)
  
  hx = (btp%bqrs(1)-btp%aqrs(1))/real(nx-1, wp)
  hy = (btp%bqrs(2)-btp%aqrs(2))/real(ny-1, wp)
  hz = (btp%bqrs(3)-btp%aqrs(3))/real(nz-1, wp)


  allocate(xi(nx), yj(ny), zk(nz))

  do i = 1, nx
     xi(i) = btp%aqrs(1) + real(i-1, wp)*hx
  end do

   do j = 1, ny
     yj(j) = btp%aqrs(2) + real(j-1, wp)*hy
  end do

   do k = 1, nz
     zk(k) = btp%aqrs(3) + real(k-1, wp)*hz
  end do

  mx = B%G%C%mq
  px = B%G%C%pq

  my = B%G%C%mr
  py = B%G%C%pr

  mz = B%G%C%ms
  pz = B%G%C%ps

  nx = B%G%C%nq
  ny = B%G%C%nr
  nz = B%G%C%ns

  !FF(mx:px,my:py,mz:pz) = f(mx:px,my:py,mz:pz)                                                                                                                                                  

  iz = mz
  fz = pz
  if (mz == 1) iz = mz 
  if (pz == nz) fz = pz 

  iy = my
  fy = py
  if (my == 1) iy= my 
  if (py == ny) fy = py 

  ix = mx
  fx = px
  if (mx == 1 ) ix = mx 
  if (px == nx) fx = px 

  ndp = 3
  stenc = 1
  
  do ii = 1,3
     FF(B%G%C%mq:B%G%C%pq,B%G%C%mr:B%G%C%pr,B%G%C%ms:B%G%C%ps) = &
           B%M%M(B%G%C%mq:B%G%C%pq,B%G%C%mr:B%G%C%pr,B%G%C%ms:B%G%C%ps, ii)

     do i =  ix, fx
        do j = iy, fy
           do k = iz, fz

              
              i_m = i-stenc
              i_p = i+stenc
              
              j_m = j-stenc
              j_p = j+stenc
              
              k_m = k-stenc
              k_p = k+stenc

              ! x-direction
              if (i == mx) then
                 
                 i_m = i
                 i_p = i+ndp-1
                 
              end if
         
              if (i == mx+1) then
                 i_m = i-1
                 i_p = i+ndp-2
              end if
              
              if (i == mx+2) then
                 i_m = i-2
                 i_p = i-3+ndp
              end if
        
              if (i == px) then
                 
                 i_m = i-(ndp-1)
                 i_p = i
                 
              end if
              
              if (i == px-1) then
                 
                 i_m = i+1-(ndp-1)
                 i_p = i+1
                 
              end if
              
              if (i == px-2) then
                 
                 i_m = i+2-(ndp-1)
                 i_p = i+2
                 
              end if

              ! y-direction
              if (j == my) then
                 
                 j_m = j
                 j_p = j+ndp-1
                 
              end if
         
              if (j == my+1) then
                 
                 j_m = j-1
                 j_p = j+ndp-2
                 
              end if
              
              if (j == my+2) then
                 
                 j_m = j-2
                 j_p = j-3+ndp
                 
              end if
        
              if (j == py) then
                 j_m = j-(ndp-1)
                 j_p = j
              end if
              
              if (j == py-1) then
                 
                 j_m = j+1-(ndp-1)
                 j_p = j+1
                 
              end if
              
              if (j == py-2) then
                 
                 j_m = j+2-(ndp-1)
                 j_p = j+2
                 
              end if


               ! z-direction
              if (k == mz) then
                 
                 k_m = k
                 k_p = k+ndp-1
                 
              end if
         
              if (k == mz+1) then
                 
                 k_m = k-1
                 k_p = k+ndp-2
                 
              end if
              
              if (k == mz+2) then
                 
                 k_m = k-2
                 k_p = k-3+ndp
                 
              end if
        
              if (k == pz) then
                 
                 k_m = k-(ndp-1)
                 k_p = k
                 
              end if
              
              if (k == pz-1) then
                 
                 k_m = k+1-(ndp-1)
                 k_p = k+1
                 
              end if
              
              if (k == pz-2) then
                 
                 k_m = k+2-(ndp-1)
                 k_p = k+2
                 
              end if          

           x = B%G%X(i,j,k,1)
           y = B%G%X(i,j,k,2)
           z = B%G%X(i,j,k,3)

           F_ijk = 0.0_wp
           

           ! construct 3rd order lagrange interpolation                                                                                                                                                   
           do i0 =  i_m, i_p
              
              do j0 = j_m, j_p
                 
                 do k0 = k_m, k_p
                                        
                    a_i = lagrange(i_m, i_p, i0, x, xi)
                    a_j = lagrange(j_m, j_p, j0, y, yj)
                    a_k = lagrange(k_m, k_p, k0, z, zk)
                    
                    F_ijk = F_ijk + a_i*a_j*a_k*B%M%M(i0, j0, k0, ii)

                 end do
              end do
           end do
           FF(i,j,k) = F_ijk
        end do
     end do
  end do


  do i =  B%G%C%mq, B%G%C%pq                                                                                                                                                                
     do j = B%G%C%mr, B%G%C%pr
        do k = B%G%C%ms, B%G%C%ps

           if (FF(i,j,k) .le. 0.0 .and.  B%M%M(i,j,k,ii) > 0.0) STOP 'material interpolation is not compatible'
           B%M%M(i,j,k,ii) = FF(i,j,k)       
           
        end do
     end do
  end do
  
end do

  
!free memory
deallocate(FF)

end subroutine interpolatematerials


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


end module inter_material
    
    

      
