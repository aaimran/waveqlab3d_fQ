module slice_output

  use mpi
  use mpi3dbasic, only : pw,ps,MPI_REAL_PW,MPI_REAL_PS,error
  use mpi3dio
  use datatypes, only : slice_type
  implicit none

contains

  subroutine init_slice_output(input,name,slicer, C)

    use mpi3dcomm
    !       use fault_output, only : open_file

    implicit none

    integer, intent(in) :: input
    character(*), intent(in) :: name

    type(cartesian3d_t), intent(in) :: C
    type(slice_type), intent(inout) :: slicer
    character(len=256) :: filename1,filename2
    character(len=16) :: extension1,extension2
    integer :: locationH,locationV,stat
    logical :: horz_slice,vert_slice

    ! Read in slice variables from output_list

    namelist / slice_list / horz_slice,vert_slice,locationH,locationV

    ! Set defaults

    slicer%horz_slice = .false.
    slicer%vert_slice = .false.
    locationH = 1
    locationV = 1

    ! Read from input file

    rewind(input)
    read(input,nml=slice_list,iostat=stat)
    if (stat>0) stop 'error reading namelist slice_list'

    slicer%horz_slice = horz_slice
    slicer%vert_slice = vert_slice
    slicer%locationH = locationH
    slicer%locationV = locationV

    if (slicer%horz_slice .eqv. .true.) then

       extension1 = '.Hslice'
       filename1 = trim(name) // extension1

       call subarray(C%nq,C%ns,C%mq,C%pq,C%ms,C%ps, MPI_REAL_PS, slicer%array_s)

       ! open file
       call open_file_distributed(slicer%Hslice, filename1, &
            'write', C%comm, slicer%array_s, ps)

    end if

    if (slicer%vert_slice .eqv. .true.) then

       extension2 = '.Vslice'
       filename2 = trim(name) // extension2

       call subarray(C%nq,C%nr,C%mq,C%pq,C%mr,C%pr, MPI_REAL_PS, slicer%array_s)

       ! open file
       call open_file_distributed(slicer%Vslice, filename2, &
            'write', C%comm, slicer%array_s, ps)


    end if

  end subroutine init_slice_output

  subroutine write_slice(F, C, slicer)
    use mpi3dcomm

    implicit none

    type(cartesian3d_t), intent(in) :: C
    type(slice_type), intent(in) :: slicer
    real(kind = wp), dimension(:,:,:,:), allocatable, intent(in) :: F

    if (slicer%horz_slice .eqv. .true.) &
         call write_file_distributed(slicer%Hslice, F(C%mq:C%pq,slicer%locationH,C%ms:C%ps,1:9))

    if (slicer%vert_slice .eqv. .true.) &
         call write_file_distributed(slicer%Vslice, F(C%mq:C%pq,C%mr:C%pr, slicer%locationV,1:9))

  end subroutine write_slice

  subroutine end_slice_output(slicer)
        type(slice_type), intent(inout) :: slicer

        if (slicer%horz_slice) call close_file_distributed(slicer%Hslice)
        if (slicer%vert_slice) call close_file_distributed(slicer%Vslice)

  end subroutine end_slice_output
end module slice_output
