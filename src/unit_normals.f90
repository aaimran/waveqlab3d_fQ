module unit_normals
  
  use common, only : wp
  
contains
  
  
  subroutine local_orth_vectors_q(i, G, B)
    
    use mpi3dbasic, only : rank
    use datatypes, only: block_grid_t, block_boundary
    
    type(block_boundary),intent(inout) :: B
    type(block_grid_t),intent(in) :: G
    integer, intent(in) :: i
    integer :: my, mz, py, pz
    real(kind = wp), dimension(1:3) :: l, m, n  !< three local othogonal vectors                                                                     
    integer :: y, z
    real(kind = wp) :: norm, q_x, q_y, q_z
    
    my = G%C%mr
    mz = G%C%ms
    py = G%C%pr
    pz = G%C%ps
    
    l(:) = 0.0_wp
    m(:) = 0.0_wp
    
    do z = mz, pz
       do y = my, py
          
          q_x = G%metricx(i, y, z, 1)
          q_y = G%metricy(i, y, z, 1)
          q_z = G%metricz(i, y, z, 1)
          
          norm = sqrt(q_x**2.0_wp + q_y**2.0_wp + q_z**2.0_wp)
          
          ! outward unit vectors 
          n(1) = q_x/norm
          n(2) = q_y/norm
          n(3) = q_z/norm
          
          ! compute two other orthogonal vectors, l, m 
          call gen_orth_vectors(l, m, n)
          
          B%n_n(y, z, 1) = n(1)
          B%n_n(y, z, 2) = n(2)
          B%n_n(y, z, 3) = n(3)
          
          B%n_m(y, z, 1) = m(1)
          B%n_m(y, z, 2) = m(2)
          B%n_m(y, z, 3) = m(3)

          B%n_l(y, z, 1) = l(1)
          B%n_l(y, z, 2) = l(2)
          B%n_l(y, z, 3) = l(3)
          
       end do
    end do
  end subroutine  Local_orth_vectors_q
  
  subroutine local_orth_vectors_r(i, G, B)
    
    use mpi3dbasic, only : rank
    use datatypes, only: block_grid_t, block_boundary
    
    type(block_boundary),intent(inout) :: B
    type(block_grid_t),intent(in) :: G
    integer, intent(in) :: i
    integer :: mx, mz, px, pz
    real(kind = wp), dimension(1:3) :: l, m, n  !< three local othogonal vectors                                                                       
    integer :: x, z
    real(kind = wp) :: norm, r_x, r_y, r_z

    mx = G%C%mq
    mz = G%C%ms
    px = G%C%pq
    pz = G%C%ps
    
    l = 0.0_wp
    m = 0.0_wp
    
    do z = mz, pz
       do x = mx, px
          
          r_x = G%metricx(x, i, z, 2)
          r_y = G%metricy(x, i, z, 2)
          r_z = G%metricz(x, i, z, 2)
          
          norm = sqrt(r_x**2.0_wp + r_y**2.0_wp + r_z**2.0_wp)
          
          ! outward unit vectors                                                                                                                       
          n(1) = r_x/norm
          n(2) = r_y/norm
          n(3) = r_z/norm
          
          ! compute two other orthogonal vectors, l, m                                                                                                 
          call gen_orth_vectors(l, m, n)
          
          B%n_n(x, z, 1) = n(1)
          B%n_n(x, z, 2) = n(2)
          B%n_n(x, z, 3) = n(3)
          
          B%n_m(x, z, 1) = m(1)
          B%n_m(x, z, 2) = m(2)
          B%n_m(x, z, 3) = m(3)
          
          B%n_l(x, z, 1) = l(1)
          B%n_l(x, z, 2) = l(2)
          B%n_l(x, z, 3) = l(3)
          
       end do
    end do
  end subroutine  Local_orth_vectors_r
  
  subroutine local_orth_vectors_s(i, G, B)

    use mpi3dbasic, only : rank
    use datatypes, only: block_grid_t, block_boundary
    
    type(block_boundary),intent(inout) :: B
    type(block_grid_t),intent(in) :: G
    integer, intent(in) :: i
    integer :: mx, my, px, py
    real(kind = wp), dimension(1:3) :: l, m, n  !< three local othogonal vectors                                                                                                                                                                                                                           
    integer :: x, y
    real(kind = wp) :: norm, s_x, s_y, s_z

    mx = G%C%mq
    my = G%C%mr
    px = G%C%pq
    py = G%C%pr

    l(:) = 0.0_wp
    m(:) = 0.0_wp

    do y = my, py
       do x = mx, px

          s_x = G%metricx(x, y, i, 3)
          s_y = G%metricy(x, y, i, 3)
          s_z = G%metricz(x, y, i, 3)

          norm = sqrt(s_x**2.0_wp + s_y**2.0_wp + s_z**2.0_wp)

          ! outward unit vectors
          n(1) = s_x/norm
          n(2) = s_y/norm
          n(3) = s_z/norm

          ! compute two other orthogonal vectors, l, m
          call gen_orth_vectors(l, m, n)

          B%n_n(x, y, 1) = n(1)
          B%n_n(x, y, 2) = n(2)
          B%n_n(x, y, 3) = n(3)

          B%n_m(x, y, 1) = m(1)
          B%n_m(x, y, 2) = m(2)
          B%n_m(x, y, 3) = m(3)
          
          B%n_l(x, y, 1) = l(1)
          B%n_l(x, y, 2) = l(2)
          B%n_l(x, y, 3) = l(3)

       end do
    end do
  end subroutine  Local_orth_vectors_s

  subroutine local_orth_vectors(i, G, B)
    !subroutine Local_Orth_Vectors(G%metricx, X_y, X_z, X_l, X_m, X_n, &
    !mby, mbz, pby, pbz, my, mz, py, pz)
    
    use mpi3dbasic, only : rank
    use datatypes, only: block_grid_t, block_boundary
    ! Compute local orthogonal unit vectors
    type(block_boundary),intent(inout) :: B
    type(block_grid_t),intent(in) :: G
    integer, intent(in) :: i
    integer :: my, mz, py, pz
    real(kind = wp), dimension(1:3) :: l, m, n  !< three local othogonal vectors

    integer :: y, z
    real(kind = wp) :: norm, q_x, q_y, q_z

    my = G%C%mr
    mz = G%C%ms
    py = G%C%pr
    pz = G%C%ps

    l(:) = 0.0_wp
    m(:) = 0.0_wp

    do z = mz, pz
       do y = my, py

          q_x = G%metricx(i, y, z, 1)
          q_y = G%metricy(i, y, z, 1)
          q_z = G%metricz(i, y, z, 1)

          norm = sqrt(q_x**2.0_wp + q_y**2.0_wp + q_z**2.0_wp)

          ! outward unit vectors
          n(1) = q_x/norm
          n(2) = q_y/norm
          n(3) = q_z/norm

          ! compute two other orthogonal vectors, l, m
          call gen_orth_vectors(l, m, n)

          B%n_n(y, z, 1) = n(1)
          B%n_n(y, z, 2) = n(2)
          B%n_n(y, z, 3) = n(3)

          B%n_m(y, z, 1) = m(1)
          B%n_m(y, z, 2) = m(2)
          B%n_m(y, z, 3) = m(3)

          B%n_l(y, z, 1) = l(1)
          B%n_l(y, z, 2) = l(2)
          B%n_l(y, z, 3) = l(3)

          !==================================================================================
          !==================================================================================

       end do
    end do
    !     stop
  end subroutine  Local_Orth_Vectors


  subroutine  gen_orth_vectors(l, m, n)

    real(kind = wp), dimension(:), intent(in) :: n
    real(kind = wp), dimension(:), intent(inout) :: l, m
    real(kind = wp) :: tol, diff_norm1, diff_norm2

    tol = 1.0e-12_wp

    m(1) = 0.0_wp
    m(2) = 1.0_wp
    m(3) = 0.0_wp

    if (abs(sqrt(n(1)**2.0_wp+n(2)**2+n(3)**2.0_wp)-1.0_wp) .ge. tol) stop 'input vector must be a unit vector'

    diff_norm1 = sqrt((n(1)-m(1))**2.0_wp + (n(2)-m(2))**2.0_wp + (n(3)-m(3))**2.0_wp)

    diff_norm2 = sqrt((n(1)+m(1))**2.0_wp + (n(2)+m(2))**2.0_wp + (n(3)+m(3))**2.0_wp)

    if (diff_norm1 .ge. tol .and. diff_norm2 .ge. tol) then

       call Gram_Schmidt(n, m)

    else
       m(1) = 0.0_wp
       m(2) = 0.0_wp
       m(3) = 1.0_wp
       call Gram_Schmidt(n, m )

    end if

    l(1) = n(2)*m(3)-n(3)*m(2)
    l(2) = n(3)*m(1)-n(1)*m(3)
    l(3) = n(1)*m(2)-n(2)*m(1)

  end subroutine gen_orth_vectors


  subroutine Gram_Schmidt(y, z)

    ! Gram Schmidt orthonormalization

    real(kind = wp), dimension(:), intent(in) :: y
    real(kind = wp), dimension(:), intent(inout) :: z

    real(kind = wp) :: a_yz

    a_yz = y(1)*z(1) + y(2)*z(2) + y(3)*z(3)
    z = z - a_yz*y
    z = 1.0_wp/sqrt(z(1)**2.0_wp + z(2)**2.0_wp + z(3)**2.0_wp)*z

  end subroutine Gram_Schmidt


end module unit_normals
