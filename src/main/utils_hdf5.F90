!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: utils_hdf5
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: hdf5, iso_c_binding
!+
!--------------------------------------------------------------------------
module utils_hdf5
 use hdf5
 use iso_c_binding, only:c_loc

 implicit none

 public :: write_to_hdf5,    &
           read_from_hdf5,   &
           open_hdf5file,    &
           create_hdf5file,  &
           close_hdf5file,   &
           open_hdf5group,   &
           create_hdf5group, &
           close_hdf5group,  &
           HID_T

 private

 integer, parameter :: compression_level = 9

 interface write_to_hdf5
  module procedure write_real_kind4,              & ! real(4)
                   write_real_kind8,              & ! real(8)
                   write_real_1d_array_kind4,     & ! 1d real(4) arrays
                   write_real_1d_array_kind8,     & ! 1d real(8) arrays
                   write_real_2d_array_kind4,     & ! 2d real(4) arrays
                   write_real_2d_array_kind8,     & ! 2d real(8) arrays
                   write_real_3d_array_kind4,     & ! 3d real(4) arrays
                   write_real_3d_array_kind8,     & ! 3d real(8) arrays
                   write_integer_kind4,           & ! integer(4)
                   write_integer_1d_array_kind1,  & ! 1d integer(1) arrays
                   write_integer_1d_array_kind4,  & ! 1d integer(4) arrays
                   write_string                     ! strings
 end interface write_to_hdf5

 interface read_from_hdf5
  module procedure read_real_kind4,              & ! real(4)
                   read_real_kind8,              & ! real(8)
                   read_real_1d_array_kind4,     & ! 1d real(4) arrays
                   read_real_1d_array_kind8,     & ! 1d real(8) arrays
                   read_real_2d_array_kind4,     & ! 2d real(4) arrays
                   read_real_2d_array_kind8,     & ! 2d real(8) arrays
                   read_real_3d_array_kind4,     & ! 3d real(4) arrays
                   read_real_3d_array_kind8,     & ! 3d real(8) arrays
                   read_integer_kind4,           & ! integer(4)
                   read_integer_1d_array_kind1,  & ! 1d integer(1) arrays
                   read_integer_1d_array_kind4,  & ! 1d integer(4) arrays
                   read_string                     ! strings
 end interface read_from_hdf5

 interface create_hdf5group
  module procedure h5gcreate_f
 end interface create_hdf5group

 interface open_hdf5group
  module procedure h5gopen_f
 end interface open_hdf5group

 interface close_hdf5group
  module procedure h5gclose_f
 end interface close_hdf5group

contains

subroutine create_hdf5file(filename,file_id,error)
 character(len=*), intent(in)  :: filename
 integer(HID_T),   intent(out) :: file_id
 integer,          intent(out) :: error
 integer :: filter_info
 integer :: filter_info_both
 logical :: avail

 ! Initialise HDF5
 call h5open_f(error)
 if (error /= 0) then
    write(*,'("cannot initialise HDF5",/)')
    return
 endif

 ! Check if gzip compression is available.
 call h5zfilter_avail_f(h5z_filter_deflate_f,avail,error)
 if (.not.avail) then
    write(*,'("gzip filter not available.",/)')
    return
 endif
 call h5zget_filter_info_f(h5z_filter_deflate_f,filter_info,error)
 filter_info_both=ior(h5z_filter_encode_enabled_f,h5z_filter_decode_enabled_f)
 if (filter_info /= filter_info_both) then
    write(*,'("gzip filter not available for encoding and decoding.",/)')
    return
 endif

 ! Create file
 call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 file",/)')
    return
 endif

end subroutine create_hdf5file

subroutine open_hdf5file(filename,file_id,error)
 character(len=*), intent(in)  :: filename
 integer(HID_T),   intent(out) :: file_id
 integer,          intent(out) :: error

 ! Initialise HDF5
 call h5open_f(error)
 if (error /= 0) then
    write(*,'("cannot initialise HDF5",/)')
    return
 endif

 ! Open file
 call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
 if (error /= 0) then
    write(*,'("cannot open HDF5 file",/)')
    return
 endif

end subroutine open_hdf5file

subroutine close_hdf5file(file_id,error)
 integer(HID_T), intent(in)  :: file_id
 integer,        intent(out) :: error

 ! Close file
 call h5fclose_f(file_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 file",/)')
    return
 endif

 ! Close HDF5
 call h5close_f(error)
 if (error /= 0) then
    write(*,'("cannot close HDF5",/)')
    return
 endif

end subroutine close_hdf5file

subroutine write_real_kind4(x,name,id,error)
 real(kind=4),   intent(in)  :: x
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T) :: dspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: dtype_id

 dtype_id = H5T_NATIVE_REAL

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_kind4

subroutine write_real_kind8(x,name,id,error)
 real(kind=8),   intent(in)  :: x
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T) :: dspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: dtype_id

 dtype_id = H5T_NATIVE_DOUBLE

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_kind8

subroutine write_real_1d_array_kind4(x,name,id,error)
 real(kind=4),   intent(in)  :: x(:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_REAL

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_1d_array_kind4

subroutine write_real_1d_array_kind8(x,name,id,error)
 real(kind=8),   intent(in)  :: x(:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_DOUBLE

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_1d_array_kind8

subroutine write_real_2d_array_kind4(x,name,id,error)
 real(kind=4),   intent(in)  :: x(:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer, parameter :: ndims = 2
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_REAL

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_2d_array_kind4

subroutine write_real_2d_array_kind8(x,name,id,error)
 real(kind=8),   intent(in)  :: x(:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer, parameter :: ndims = 2
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_DOUBLE

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_2d_array_kind8

subroutine write_real_3d_array_kind4(x,name,id,error)
 real(kind=4),   intent(in)  :: x(:,:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer, parameter :: ndims = 3
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_REAL

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_3d_array_kind4

subroutine write_real_3d_array_kind8(x,name,id,error)
 real(kind=8),   intent(in)  :: x(:,:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 integer,        intent(out) :: error

 integer, parameter :: ndims = 3
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_DOUBLE

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_real_3d_array_kind8

subroutine write_integer_kind4(x,name,id,error)
 integer(kind=4), intent(in)  :: x
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 integer,         intent(out) :: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T) :: dspace_id
 integer(HID_T) :: dset_id
 integer(HID_T) :: dtype_id

 dtype_id = H5T_NATIVE_INTEGER

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_integer_kind4

subroutine write_integer_1d_array_kind4(x,name,id,error)
 integer(kind=4), intent(in)  :: x(:)
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 integer,         intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_NATIVE_INTEGER

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_integer_1d_array_kind4

subroutine write_integer_1d_array_kind1(x,name,id,error)
 integer(kind=1), intent(in)  :: x(:)
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 integer,         intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HSIZE_T)   :: chunk(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(HID_T)     :: prop_id
 integer(HID_T)     :: dtype_id

 xshape = shape(x)
 chunk = shape(x)
 dtype_id = H5T_STD_I8LE

 ! Create dataspace
 call h5screate_simple_f(ndims,xshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset creation property list, add the gzip
 ! compression filter and set the chunk size.
 call h5pcreate_f(H5P_DATASET_CREATE_F,prop_id,error)
 call h5pset_deflate_f(prop_id,compression_level,error)
 call h5pset_chunk_f(prop_id,ndims,chunk,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 property list",/)')
    return
 endif

 ! Create dataset in file
 call h5dcreate_f(id,name,dtype_id,dspace_id,dset_id,error,prop_id)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Write to file
 call h5dwrite_f(dset_id,dtype_id,x,xshape,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close property list
 call h5pclose_f(prop_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 property list",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

end subroutine write_integer_1d_array_kind1

subroutine write_string(str,name,id,error)
 character(*),    intent(in), target :: str
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 integer,         intent(out) :: error

 integer, parameter :: ndims = 0
 integer(HSIZE_T)   :: sshape(ndims)
 integer(HID_T)     :: dspace_id
 integer(HID_T)     :: dset_id
 integer(SIZE_T)    :: slength
 integer(HID_T)     :: filetype
 type(C_PTR)        :: cpointer

 slength = len(str)
 sshape  = shape(str)

 ! Create file datatypes. Save the string as FORTRAN string
 call h5tcopy_f(H5T_FORTRAN_S1,filetype,error)
 call h5tset_size_f(filetype,slength,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 datatype",/)')
    return
 endif

 ! Create dataspace
 call h5screate_simple_f(ndims,sshape,dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataspace",/)')
    return
 endif

 ! Create the dataset in file
 call h5dcreate_f(id,name,filetype,dspace_id,dset_id,error)
 if (error /= 0) then
    write(*,'("cannot create HDF5 dataset",/)')
    return
 endif

 ! Find C pointer
 cpointer = c_loc(str(1:1))

 ! Write to file
 call h5dwrite_f(dset_id,filetype,cpointer,error)
 if (error /= 0) then
    write(*,'("cannot write to HDF5 file",/)')
    return
 endif

 ! Close dataset
 call h5dclose_f(dset_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataset",/)')
    return
 endif

 ! Close dataspace
 call h5sclose_f(dspace_id,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 dataspace",/)')
    return
 endif

 ! Close datatype
 call h5tclose_f(filetype,error)
 if (error /= 0) then
    write(*,'("cannot close HDF5 datatype",/)')
    return
 endif

end subroutine write_string

subroutine read_real_kind4(x,name,id,got,error)
 real(kind=4),   intent(out) :: x
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T) :: dset_id
 integer :: errors(3)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_REAL,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_kind4

subroutine read_real_kind8(x,name,id,got,error)
 real(kind=8),   intent(out) :: x
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T) :: dset_id
 integer :: errors(3)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_kind8

subroutine read_real_1d_array_kind4(x,name,id,got,error)
 real(kind=4),   intent(out) :: x(:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_REAL,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_1d_array_kind4

subroutine read_real_1d_array_kind8(x,name,id,got,error)
 real(kind=8),   intent(out) :: x(:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_1d_array_kind8

subroutine read_real_2d_array_kind4(x,name,id,got,error)
 real(kind=4),   intent(out) :: x(:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer, parameter :: ndims = 2
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_REAL,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_2d_array_kind4

subroutine read_real_2d_array_kind8(x,name,id,got,error)
 real(kind=8),   intent(out) :: x(:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer, parameter :: ndims = 2
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_2d_array_kind8

subroutine read_real_3d_array_kind4(x,name,id,got,error)
 real(kind=4),   intent(out) :: x(:,:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer, parameter :: ndims = 3
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_REAL,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_3d_array_kind4

subroutine read_real_3d_array_kind8(x,name,id,got,error)
 real(kind=8),   intent(out) :: x(:,:,:)
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer, parameter :: ndims = 3
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_real_3d_array_kind8

subroutine read_integer_kind4(x,name,id,got,error)
 integer(kind=4), intent(out) :: x
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 logical,         intent(out) :: got
 integer,         intent(out) :: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T) :: dset_id
 integer :: errors(3)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_integer_kind4

subroutine read_integer_1d_array_kind1(x,name,id,got,error)
 integer(kind=1), intent(out) :: x(:)
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 logical,         intent(out) :: got
 integer,         intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_STD_I8LE,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_integer_1d_array_kind1

subroutine read_integer_1d_array_kind4(x,name,id,got,error)
 integer(kind=4), intent(out) :: x(:)
 character(*),    intent(in)  :: name
 integer(HID_T),  intent(in)  :: id
 logical,         intent(out) :: got
 integer,         intent(out) :: error

 integer, parameter :: ndims = 1
 integer(HSIZE_T)   :: xshape(ndims)
 integer(HID_T)     :: dset_id
 integer :: errors(3)

 xshape = shape(x)

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 ! Open dataset
 call h5dopen_f(id,name,dset_id,errors(1))

 ! Read dataset
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,x,xshape,errors(2))

 ! Close dataset
 call h5dclose_f(dset_id,errors(3))

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_integer_1d_array_kind4

subroutine read_string(str,name,id,got,error)
 character(*),   intent(out) :: str
 character(*),   intent(in)  :: name
 integer(HID_T), intent(in)  :: id
 logical,        intent(out) :: got
 integer,        intent(out) :: error

 integer :: errors(12)

 integer,         parameter :: dim0 = 1
 integer(SIZE_T), parameter :: sdim = 100

 integer(HSIZE_T) :: dims(1) = (/dim0/)
 integer(HSIZE_T) :: maxdims(1)

 integer(HID_T) :: filetype,memtype,space,dset

 character(LEN=sdim), allocatable, target :: rdata(:)
 integer(SIZE_T) :: size
 type(c_ptr) :: f_ptr

 ! Check if dataset exists
 call h5lexists_f(id,name,got,error)
 if (.not.got) return

 call h5dopen_f(id,name,dset,errors(1))

 ! Get the datatype and its size.
 call h5dget_type_f(dset,filetype,errors(2))
 call h5tget_size_f(filetype,size,errors(3))

 ! Make sure the declared length is large enough,
 ! the C string contains the null character.
 if (size > sdim+1) then
    print*,'ERROR: Character LEN is too small'
    stop
 endif

 ! Get dataspace.
 call h5dget_space_f(dset,space,errors(4))
 call h5sget_simple_extent_dims_f(space,dims,maxdims,errors(5))

 allocate(rdata(1:dims(1)))

 ! Create the memory datatype.
 call H5Tcopy_f(H5T_FORTRAN_S1,memtype,errors(6))
 call H5Tset_size_f(memtype,sdim,errors(7))

 ! Read the data.
 f_ptr = C_LOC(rdata(1)(1:1))
 call H5Dread_f(dset,memtype,f_ptr,errors(8),space)

 ! Close and release resources.
 call H5Dclose_f(dset,errors(9))
 call H5Sclose_f(space,errors(10))
 call H5Tclose_f(filetype,errors(11))
 call H5Tclose_f(memtype,errors(12))

 str = rdata(1)

 deallocate(rdata)

 error = maxval(abs(errors))

 if (error /= 0) got = .false.

end subroutine read_string

end module utils_hdf5
