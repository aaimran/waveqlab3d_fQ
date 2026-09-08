module parallel_write

  use common, only : wp
  use mpi3dcomm, only : cartesian3d_t
  use grid, only : block_grid_t

  implicit none

contains

  subroutine init_seismogram_parallel(G, point, loc)
    
    type(block_grid_t), intent(in) :: G
    real(kind = wp), dimension(3), intent(in) :: point
    real(kind = wp), dimension(3) :: grid_point
    real(kind = wp) :: max_dist
    integer, dimension(3) :: loc
    integer :: i, j, k


    max_dist = 1000000.0_wp
    do k = G%C%ms, G%C%ps
      do j = G%C%mr, G%C%pr
        do i = G%C%mq, G%C%pq
          grid_point  = G%X(i,j,k,1:3)
          if (dist(grid_point, point) < max_dist) then
            max_dist = dist(grid_point, point)
            loc = [i, j, k]
          end if
        end do
      end do
    end do

  contains
    pure function dist(p1, p2)
      real(kind = wp) :: dist
      real(kind = wp), dimension(3),intent(in) :: p1, p2

      dist = sqrt(dot_product((p1-p2),(p1-p2)))

    end function dist


  end subroutine init_seismogram_parallel

end module parallel_write