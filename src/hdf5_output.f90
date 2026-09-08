!> hdf5 output for waveqlab3D


module hdf5_output

  use mpi
  use hdf5

  implicit none

  ! Variables for writing HDF5
  integer :: comm, info
  character(len=9), parameter :: filename = "output.h5"
  character(len=5), parameter :: groupnameb = "/data"
  character(len=5), parameter :: groupname0 = "var0d", groupname1 = "var1d", &
       groupname2 = "var2d", groupname3 = "var3d", &
       groupname4 = "var4d"
  integer(hid_t) :: file_id, plist_id, groupb_id, group0_id, group1_id, &
       group2_id, group3_id, group4_id, dspace_id, dset_id, loc_id
  integer :: hdferror

contains

  subroutine start_hdf_output()

    ! Open the fortran hdf5 interface
    call h5open_f(hdferror)

    ! Setup file access property list for parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferror)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, hdferror)

    ! Open a new file using default properties
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferror, access_prp = plist_id)
    call h5pclose_f(plist_id, hdferror)

    ! Create the "/data" group in the root
    call h5gcreate_f(file_id, groupnameb, groupb_id, hdferror)
    ! Now create the 0-4D var groups in "/data" group
    call h5gcreate_f(groupb_id, groupname0, group0_id, hdferror)
    call h5gcreate_f(groupb_id, groupname1, group1_id, hdferror)
    call h5gcreate_f(groupb_id, groupname2, group2_id, hdferror)
    call h5gcreate_f(groupb_id, groupname3, group3_id, hdferror)
    call h5gcreate_f(groupb_id, groupname4, group4_id, hdferror)

  end subroutine start_hdf_output

  subroutine finish_hdf_output()

    call h5gclose_f(group4_id, hdferror)
    call h5gclose_f(group3_id, hdferror)
    call h5gclose_f(group2_id, hdferror)
    call h5gclose_f(group1_id, hdferror)
    call h5gclose_f(group0_id, hdferror)
    call h5gclose_f(groupb_id, hdferror)
    call h5fclose_f(file_id, hdferror)
    call h5close_f(hdferror)
  end subroutine finish_hdf_output

  subroutine create_dataset(name, rank, shape)

    character(len=*), intent(in) :: name
    integer, intent(in) :: rank
    integer, dimension(:), intent(in) :: shape

    integer(hsize_t), dimension(rank+1) :: shapep1, maxshape

    if (rank > 0) shapep1(1:rank) = shape(1:rank) !< If there are any cartesian dimensions
    shapep1(rank+1) = 0 !< Initial dimension for time is 0, max will be unlimited.
    maxshape = shapep1 !< Max of cartesian dimensions is same as shape
    maxshape(rank+1) = H5S_UNLIMITED_F

    call h5screate_simple_f(rank+1, shapep1, dspace_id, hdferror, maxshape)

    ! Use chunking for the time dimension
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdferror)
    shapep1(rank+1) = 8
    call h5pset_chunk_f(plist_id, rank+1, shapep1, hdferror)

    ! Decide where to put the data based on rank
    select case(rank)
    case (0)
       loc_id = group1_id
    case (1) 
       loc_id = group2_id
    case (2)
       loc_id = group3_id
    case (3)
       loc_id = group4_id
    end select

    ! Create the dataset
    call h5dcreate_f(loc_id, name, H5T_IEEE_F32LE, dspace_id, dset_id, hdferror, plist_id)

    call h5sclose_f(dspace_id, hdferror)
    call h5pclose_f(plist_id, hdferror)
    call h5dclose_f(dset_id, hdferror)

  end subroutine create_dataset


end module hdf5_output
