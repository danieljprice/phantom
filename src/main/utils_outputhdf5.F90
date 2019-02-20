module utils_outputhdf5
 use hdf5, only:h5screate_f,h5sclose_f,h5screate_simple_f,h5dcreate_f,h5dclose_f,h5dwrite_f
 use hdf5, only:HID_T,H5F_ACC_TRUNC_F,HSIZE_T,H5S_SCALAR_F,H5T_NATIVE_DOUBLE,H5T_NATIVE_INTEGER
 implicit none
 public :: write_to_hdf5

 private

 interface write_to_hdf5
  module procedure write_scalar, write_array_1d, write_array_2d, write_array_3d, write_array_4d, write_array_6d, &
                   write_scalar_int, write_scalar_intkind8, write_intarray_1d, write_intarray_1dkind8, write_intarray_1dkind1, &
                   write_array_1dkind4, write_array_2dkind4, write_string
 end interface

contains

subroutine write_scalar(x, name, id, error)
 real,           intent(in) :: x
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T)    :: dspace_id
 integer(HID_T)    :: dset_id
 integer :: errors(5)

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_scalar

subroutine write_array_1d(x, name, id, error)
 real,           intent(in) :: x(:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_1d

subroutine write_array_1dkind4(x, name, id, error)
 real(kind=4),   intent(in) :: x(:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_1dkind4

subroutine write_array_2d(x, name, id, error)
 real,           intent(in) :: x(:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 2
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_2d

subroutine write_array_2dkind4(x, name, id, error)
 real(kind=4),   intent(in) :: x(:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 2
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_2dkind4

subroutine write_array_3d(x, name, id, error)
 real,           intent(in) :: x(:,:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 3
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_3d

subroutine write_array_4d(x, name, id, error)
 real,           intent(in) :: x(:,:,:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 4
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_4d

subroutine write_array_6d(x, name, id, error)
 real,           intent(in) :: x(:,:,:,:,:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 6
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_array_6d

subroutine write_scalar_int(x, name, id, error)
 integer,        intent(in) :: x
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T)    :: dspace_id
 integer(HID_T)    :: dset_id
 integer :: errors(5)

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_scalar_int

subroutine write_scalar_intkind8(x, name, id, error)
 integer(kind=8), intent(in) :: x
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id
 integer,         intent(out):: error

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T)    :: dspace_id
 integer(HID_T)    :: dset_id
 integer :: errors(5)

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_scalar_intkind8

subroutine write_intarray_1d(x, name, id, error)
 integer,        intent(in) :: x(:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id
 integer,        intent(out):: error

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_intarray_1d

subroutine write_intarray_1dkind8(x, name, id, error)
 integer(kind=8), intent(in) :: x(:)
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id
 integer,         intent(out):: error

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_intarray_1dkind8

subroutine write_intarray_1dkind1(x, name, id, error)
 integer(kind=1), intent(in) :: x(:)
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id
 integer,         intent(out):: error

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer :: errors(5)

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, errors(1))

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, errors(2))

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, errors(3))

 ! Close dataset
 call h5dclose_f(dset_id, errors(4))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(5))

 error = maxval(abs(errors))

end subroutine write_intarray_1dkind1

subroutine write_string(str, name, id, error)
 use hdf5,          only:SIZE_T,H5T_FORTRAN_S1,C_PTR
 use hdf5,          only:h5tcopy_f,h5tset_size_f,h5tclose_f
 use iso_c_binding, only:c_loc
 character(*),    intent(in), target :: str
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id
 integer,         intent(out):: error

 integer, parameter  :: ndims = 0
 integer(HSIZE_T)    :: sshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer(SIZE_T)     :: slength
 integer(HID_T)      :: filetype
 type(C_PTR)         :: cpointer
 integer :: errors(8)

 slength = len(str)
 sshape  = shape(str)

 ! Create file datatypes. Save the string as FORTRAN string
 call h5tcopy_f(H5T_FORTRAN_S1, filetype, errors(1))
 call h5tset_size_f(filetype, slength, errors(2))

 ! Create dataspace
 call h5screate_simple_f(ndims, sshape, dspace_id, errors(3))

 ! Create the dataset in file
 call h5dcreate_f(id, name, filetype, dspace_id, dset_id, errors(4))

 ! Find C pointer
 cpointer = c_loc(str(1:1))

 ! Write to file
 call h5dwrite_f(dset_id, filetype, cpointer, errors(5))

 ! Close dataset
 call h5dclose_f(dset_id, errors(6))

 ! Closet dataspace
 call h5sclose_f(dspace_id, errors(7))

 ! Close datatype
 call h5tclose_f(filetype, errors(8))

 error = maxval(abs(errors))

end subroutine write_string

end module utils_outputhdf5
