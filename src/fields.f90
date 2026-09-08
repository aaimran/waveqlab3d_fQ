module fields

  !> fields modules handles fields defined within blocks

  use common, only : wp
  use datatypes, only : block_type, block_fields
   use, intrinsic :: ieee_arithmetic
  implicit none

  !> first few indices are for spatial dimensions, last index is for field component

contains


  subroutine init_fields(F, C, physics)

!!! initialize fields in a block

    use mpi3dcomm

    implicit none

    type(block_fields),intent(out) :: F
    type(cartesian3d_t),intent(in) :: C
    character(*),intent(in) :: physics

!!! fields will be specific to block physics

    select case(physics)

    case default

       stop 'invalid block physics in init_fields'

    case('elastic')

       call allocate_array_body(F%F ,C,9, ghost_nodes=.true., Fval = 0.0_wp)
       call allocate_array_body(F%DF,C,9, ghost_nodes=.true., Fval = 1.0e40_wp)

    case('acoustic')

       call allocate_array_body(F%F , C, 4, ghost_nodes=.true.)
       call allocate_array_body(F%DF, C, 4, ghost_nodes=.true.)

    end select

  end subroutine init_fields


  subroutine scale_rates_interior(F,A)

!!! multiply rates by RK coefficient A

    implicit none

    type(block_type),intent(inout) :: F
    real(kind = wp),intent(in) :: A
    

    F%F%DF = A*F%F%DF

      if (F%M%anelastic) then
         F%M%Deta4 = A*F%M%Deta4
         F%M%Deta5 = A*F%M%Deta5
         F%M%Deta6 = A*F%M%Deta6
         F%M%Deta7 = A*F%M%Deta7
         F%M%Deta8 = A*F%M%Deta8
         F%M%Deta9 = A*F%M%Deta9
      end if
      if (F%M%anelastic_Q) then
         F%M%Deta4Q = A*F%M%Deta4Q
         F%M%Deta5Q = A*F%M%Deta5Q
         F%M%Deta6Q = A*F%M%Deta6Q
         F%M%Deta7Q = A*F%M%Deta7Q
         F%M%Deta8Q = A*F%M%Deta8Q
         F%M%Deta9Q = A*F%M%Deta9Q
      end if
      if (F%M%anelastic_Q8) then
         F%M%Deta4Q8 = A*F%M%Deta4Q8
         F%M%Deta5Q8 = A*F%M%Deta5Q8
         F%M%Deta6Q8 = A*F%M%Deta6Q8
         F%M%Deta7Q8 = A*F%M%Deta7Q8
         F%M%Deta8Q8 = A*F%M%Deta8Q8
         F%M%Deta9Q8 = A*F%M%Deta9Q8
      end if
      if (F%M%anelastic_cQ) then
         F%M%Deta4cQ=A*F%M%Deta4cQ; F%M%Deta5cQ=A*F%M%Deta5cQ
         F%M%Deta6cQ=A*F%M%Deta6cQ; F%M%Deta7cQ=A*F%M%Deta7cQ
         F%M%Deta8cQ=A*F%M%Deta8cQ; F%M%Deta9cQ=A*F%M%Deta9cQ
      end if
      if (allocated(F%M%eta4Qf8)) then
         F%M%Deta4Qf8 = A*F%M%Deta4Qf8
         F%M%Deta5Qf8 = A*F%M%Deta5Qf8
         F%M%Deta6Qf8 = A*F%M%Deta6Qf8
         F%M%Deta7Qf8 = A*F%M%Deta7Qf8
         F%M%Deta8Qf8 = A*F%M%Deta8Qf8
         F%M%Deta9Qf8 = A*F%M%Deta9Qf8
      else if (F%M%anelastic_Qf) then
         F%M%Deta4Qf = A*F%M%Deta4Qf
         F%M%Deta5Qf = A*F%M%Deta5Qf
         F%M%Deta6Qf = A*F%M%Deta6Qf
         F%M%Deta7Qf = A*F%M%Deta7Qf
         F%M%Deta8Qf = A*F%M%Deta8Qf
         F%M%Deta9Qf = A*F%M%Deta9Qf
      end if
      if (allocated(F%M%eta4_8M)) then
         F%M%Deta4_8M = A*F%M%Deta4_8M
         F%M%Deta5_8M = A*F%M%Deta5_8M
         F%M%Deta6_8M = A*F%M%Deta6_8M
         F%M%Deta7_8M = A*F%M%Deta7_8M
         F%M%Deta8_8M = A*F%M%Deta8_8M
         F%M%Deta9_8M = A*F%M%Deta9_8M
      else if (F%M%anelastic_const_Q_4M) then
         F%M%Deta4_4M = A*F%M%Deta4_4M
         F%M%Deta5_4M = A*F%M%Deta5_4M
         F%M%Deta6_4M = A*F%M%Deta6_4M
         F%M%Deta7_4M = A*F%M%Deta7_4M
         F%M%Deta8_4M = A*F%M%Deta8_4M
         F%M%Deta9_4M = A*F%M%Deta9_4M
      end if

    if( F%PMLB(1)%pml .EQV. .TRUE.) then
       if(F%G%C%mq .le.  F%PMLB(1)%N_pml) then
          F%PMLB(1)%DQ = A*F%PMLB(1)%DQ
       end if
    end if

    
    if( F%PMLB(2)%pml .EQV. .TRUE.) then
        if(F%G%C%pq .ge. (F%G%C%nq- F%PMLB(2)%N_pml+1)) then
           F%PMLB(2)%DQ = A*F%PMLB(2)%DQ
        end if
     end if
     
     if( F%PMLB(3)%pml .EQV. .TRUE.) then
        if(F%G%C%mr .le. F%PMLB(3)%N_pml) then
           F%PMLB(3)%DQ = A*F%PMLB(3)%DQ
        end if
     end if
     
     if( F%PMLB(4)%pml .EQV. .TRUE.) then
         if(F%G%C%pr .ge. (F%G%C%nr- F%PMLB(4)%N_pml+1)) then
            F%PMLB(4)%DQ = A*F%PMLB(4)%DQ
         end if
     end if
     
     if( F%PMLB(5)%pml .EQV. .TRUE.) then
        if(F%G%C%ms .le.  F%PMLB(5)%N_pml) then
           F%PMLB(5)%DQ = A*F%PMLB(5)%DQ
        end if
     end if
     
     if( F%PMLB(6)%pml .EQV. .TRUE.) then
         if(F%G%C%ps .ge. (F%G%C%ns- F%PMLB(6)%N_pml+1)) then
            F%PMLB(6)%DQ = A*F%PMLB(6)%DQ
         end if
     end if
     

  end subroutine scale_rates_interior


  subroutine update_fields_interior(F,dt)

!!! update fields using rates

    implicit none

    type(block_type),intent(inout) :: F
    real(kind = wp),intent(in) :: dt

    F%F%F = F%F%F + dt*F%F%DF

      if (F%M%anelastic) then
         F%M%eta4 = F%M%eta4 + dt*F%M%Deta4
         F%M%eta5 = F%M%eta5 + dt*F%M%Deta5
         F%M%eta6 = F%M%eta6 + dt*F%M%Deta6
         F%M%eta7 = F%M%eta7 + dt*F%M%Deta7
         F%M%eta8 = F%M%eta8 + dt*F%M%Deta8
         F%M%eta9 = F%M%eta9 + dt*F%M%Deta9
      end if
      if (F%M%anelastic_Q) then
         F%M%eta4Q = F%M%eta4Q + dt*F%M%Deta4Q
         F%M%eta5Q = F%M%eta5Q + dt*F%M%Deta5Q
         F%M%eta6Q = F%M%eta6Q + dt*F%M%Deta6Q
         F%M%eta7Q = F%M%eta7Q + dt*F%M%Deta7Q
         F%M%eta8Q = F%M%eta8Q + dt*F%M%Deta8Q
         F%M%eta9Q = F%M%eta9Q + dt*F%M%Deta9Q
      end if
      if (F%M%anelastic_Q8) then
         F%M%eta4Q8 = F%M%eta4Q8 + dt*F%M%Deta4Q8
         F%M%eta5Q8 = F%M%eta5Q8 + dt*F%M%Deta5Q8
         F%M%eta6Q8 = F%M%eta6Q8 + dt*F%M%Deta6Q8
         F%M%eta7Q8 = F%M%eta7Q8 + dt*F%M%Deta7Q8
         F%M%eta8Q8 = F%M%eta8Q8 + dt*F%M%Deta8Q8
         F%M%eta9Q8 = F%M%eta9Q8 + dt*F%M%Deta9Q8
      end if
      if (F%M%anelastic_cQ) then
         F%M%eta4cQ=F%M%eta4cQ+dt*F%M%Deta4cQ
         F%M%eta5cQ=F%M%eta5cQ+dt*F%M%Deta5cQ
         F%M%eta6cQ=F%M%eta6cQ+dt*F%M%Deta6cQ
         F%M%eta7cQ=F%M%eta7cQ+dt*F%M%Deta7cQ
         F%M%eta8cQ=F%M%eta8cQ+dt*F%M%Deta8cQ
         F%M%eta9cQ=F%M%eta9cQ+dt*F%M%Deta9cQ
      end if
      if (allocated(F%M%eta4Qf8)) then
         F%M%eta4Qf8 = F%M%eta4Qf8 + dt*F%M%Deta4Qf8
         F%M%eta5Qf8 = F%M%eta5Qf8 + dt*F%M%Deta5Qf8
         F%M%eta6Qf8 = F%M%eta6Qf8 + dt*F%M%Deta6Qf8
         F%M%eta7Qf8 = F%M%eta7Qf8 + dt*F%M%Deta7Qf8
         F%M%eta8Qf8 = F%M%eta8Qf8 + dt*F%M%Deta8Qf8
         F%M%eta9Qf8 = F%M%eta9Qf8 + dt*F%M%Deta9Qf8
      else if (F%M%anelastic_Qf) then
         F%M%eta4Qf = F%M%eta4Qf + dt*F%M%Deta4Qf
         F%M%eta5Qf = F%M%eta5Qf + dt*F%M%Deta5Qf
         F%M%eta6Qf = F%M%eta6Qf + dt*F%M%Deta6Qf
         F%M%eta7Qf = F%M%eta7Qf + dt*F%M%Deta7Qf
         F%M%eta8Qf = F%M%eta8Qf + dt*F%M%Deta8Qf
         F%M%eta9Qf = F%M%eta9Qf + dt*F%M%Deta9Qf
      end if
      if (allocated(F%M%eta4_8M)) then
         F%M%eta4_8M = F%M%eta4_8M + dt*F%M%Deta4_8M
         F%M%eta5_8M = F%M%eta5_8M + dt*F%M%Deta5_8M
         F%M%eta6_8M = F%M%eta6_8M + dt*F%M%Deta6_8M
         F%M%eta7_8M = F%M%eta7_8M + dt*F%M%Deta7_8M
         F%M%eta8_8M = F%M%eta8_8M + dt*F%M%Deta8_8M
         F%M%eta9_8M = F%M%eta9_8M + dt*F%M%Deta9_8M
      else if (F%M%anelastic_const_Q_4M) then
         F%M%eta4_4M = F%M%eta4_4M + dt*F%M%Deta4_4M
         F%M%eta5_4M = F%M%eta5_4M + dt*F%M%Deta5_4M
         F%M%eta6_4M = F%M%eta6_4M + dt*F%M%Deta6_4M
         F%M%eta7_4M = F%M%eta7_4M + dt*F%M%Deta7_4M
         F%M%eta8_4M = F%M%eta8_4M + dt*F%M%Deta8_4M
         F%M%eta9_4M = F%M%eta9_4M + dt*F%M%Deta9_4M
      end if

    if( F%PMLB(1)%pml .EQV. .TRUE.) then
       if(F%G%C%mq .le. F%PMLB(1)%N_pml) then
         F%PMLB(1)%Q = F%PMLB(1)%Q + dt*F%PMLB(1)%DQ
       end if
    end if

    
    if( F%PMLB(2)%pml .EQV. .TRUE.) then
        if(F%G%C%pq .ge. (F%G%C%nq- F%PMLB(2)%N_pml+1)) then
          F%PMLB(2)%Q = F%PMLB(2)%Q + dt*F%PMLB(2)%DQ  
        end if
     end if
     
     if( F%PMLB(3)%pml .EQV. .TRUE.) then
        if(F%G%C%mr .le.  F%PMLB(3)%N_pml) then
           F%PMLB(3)%Q = F%PMLB(3)%Q + dt*F%PMLB(3)%DQ
        end if
     end if
     
     if( F%PMLB(4)%pml .EQV. .TRUE.) then
         if(F%G%C%pr .ge. (F%G%C%nr- F%PMLB(4)%N_pml+1)) then
            F%PMLB(4)%Q = F%PMLB(4)%Q + dt*F%PMLB(4)%DQ 
         end if
     end if
     
     if( F%PMLB(5)%pml .EQV. .TRUE.) then
        if(F%G%C%ms .le.  F%PMLB(5)%N_pml) then
           F%PMLB(5)%Q = F%PMLB(5)%Q + dt*F%PMLB(5)%DQ
        end if
     end if
     
     if( F%PMLB(6)%pml .EQV. .TRUE.) then
         if(F%G%C%ps .ge. (F%G%C%ns- F%PMLB(6)%N_pml+1)) then
            F%PMLB(6)%Q = F%PMLB(6)%Q + dt*F%PMLB(6)%DQ
         end if
     end if

     
         

  end subroutine update_fields_interior


  subroutine check_block_fields_for_invalid_values(F)

    implicit none

    type(block_type), intent(in) :: F

    if (allocated(F%F%F)) call check_real_4d(F%F%F, 'F%F', F%id)
    if (allocated(F%F%DF)) call check_real_4d(F%F%DF, 'F%DF', F%id)

  end subroutine check_block_fields_for_invalid_values


  subroutine check_real_4d(arr, label, block_id)

    implicit none

    real(kind = wp), intent(in) :: arr(:,:,:,:)
    character(*), intent(in) :: label
    integer, intent(in) :: block_id

    integer :: i, j, k, l

    do l = lbound(arr, 4), ubound(arr, 4)
       do k = lbound(arr, 3), ubound(arr, 3)
          do j = lbound(arr, 2), ubound(arr, 2)
             do i = lbound(arr, 1), ubound(arr, 1)
                if (.not. ieee_is_finite(arr(i,j,k,l))) then
                   write(*,'(A,I0,2X,A,2X,A,4(I0,1X),2X,A,ES12.4)') &
                        'ERROR: non-finite value in block', block_id, trim(label), &
                        '(i,j,k,l)=', i, j, k, l, 'value=', arr(i,j,k,l)
                   stop 'check_block_fields_for_invalid_values'
                end if
             end do
          end do
       end do
    end do

  end subroutine check_real_4d


  subroutine exchange_all_fields(F, C)

    use mpi3dcomm
    use mpi3dbasic, only: rank

    implicit none

    type(block_fields), intent(inout) :: F
    type(cartesian3d_t),intent(in) :: C
    integer :: i

    do i = 1, size(F%F, 4)
       call exchange_all_neighbors(C, F%F(:,:,:,i))
    end do

  end subroutine exchange_all_fields
!!! note that there is no set_rates_interior routine in this module;
!!! such a routine, which is quite involved, is stored separately in other modules,
!!! like elastic in elastic.f90, with one for each physics

  subroutine norm_fields(F,mx,my,mz,px, py, pz, n,sum)

    implicit none

    integer, intent(in) :: mx, my, mz, px, py, pz, n
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: F
    real(kind = wp), intent(inout) :: sum

    integer :: i,j,k,l

    sum = 0.0_wp

    do l = 1,n
       do k = mz, pz
          do j = my, py
             do i = mx, px
                sum = sum + F(i,j,k,l)**2
             end do
          end do
       end do
    end do

    sum = sqrt(sum)

  end subroutine norm_fields


end module fields
