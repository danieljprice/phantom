module utils_outputhdf5
 use hdf5, only:h5screate_f,h5sclose_f,h5screate_simple_f,h5dcreate_f,h5dclose_f,h5dwrite_f
 use hdf5, only:HID_T,H5F_ACC_TRUNC_F,HSIZE_T,H5S_SCALAR_F,H5T_NATIVE_DOUBLE,H5T_NATIVE_INTEGER
implicit none
 public :: write_to_hdf5

 private

 interface write_to_hdf5
    module procedure write_scalar, write_array_1d, write_array_2d, write_array_3d, write_array_4d, write_array_6d, &
                     write_scalar_int, write_scalar_intkind8, write_intarray_1d, write_intarray_1dkind8, write_intarray_1dkind1, &
                     write_array_1dkind4, write_array_2dkind4
 end interface

contains

subroutine write_scalar(x, name, id)
 real,           intent(in) :: x
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T)    :: dspace_id
 integer(HID_T)    :: dset_id
 integer           :: error

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_scalar

subroutine write_array_1d(x, name, id)
 real,           intent(in) :: x(:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_1d

subroutine write_array_1dkind4(x, name, id)
 real(kind=4),   intent(in) :: x(:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_1dkind4

subroutine write_array_2d(x, name, id)
 real,           intent(in) :: x(:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 2
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_2d

subroutine write_array_2dkind4(x, name, id)
 real(kind=4),   intent(in) :: x(:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 2
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_2dkind4

subroutine write_array_3d(x, name, id)
 real,           intent(in) :: x(:,:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 3
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_3d

subroutine write_array_4d(x, name, id)
 real,           intent(in) :: x(:,:,:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 4
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_4d

subroutine write_array_6d(x, name, id)
 real,           intent(in) :: x(:,:,:,:,:,:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 6
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_array_6d

subroutine write_scalar_int(x, name, id)
 integer,        intent(in) :: x
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T)    :: dspace_id
 integer(HID_T)    :: dset_id
 integer           :: error

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_scalar_int

subroutine write_scalar_intkind8(x, name, id)
 integer(kind=8), intent(in) :: x
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id

 integer(HSIZE_T), parameter  :: xshape(0) = 0
 integer(HID_T)    :: dspace_id
 integer(HID_T)    :: dset_id
 integer           :: error

 ! Create dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_scalar_intkind8

subroutine write_intarray_1d(x, name, id)
 integer,        intent(in) :: x(:)
 character(*),   intent(in) :: name
 integer(HID_T), intent(in) :: id

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_intarray_1d

subroutine write_intarray_1dkind8(x, name, id)
 integer(kind=8), intent(in) :: x(:)
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)

end subroutine write_intarray_1dkind8

subroutine write_intarray_1dkind1(x, name, id)
 integer(kind=1), intent(in) :: x(:)
 character(*),    intent(in) :: name
 integer(HID_T),  intent(in) :: id

 integer, parameter  :: ndims = 1
 integer(HSIZE_T)    :: xshape(ndims)
 integer(HID_T)      :: dspace_id
 integer(HID_T)      :: dset_id
 integer             :: error

 xshape = shape(x)

 ! Create dataspace
 call h5screate_simple_f(ndims, xshape, dspace_id, error)

 ! Create dataset in file
 call h5dcreate_f(id, name, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

 ! Write to file
 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, xshape, error)

 ! Close dataset
 call h5dclose_f(dset_id, error)

 ! Closet dataspace
 call h5sclose_f(dspace_id, error)
end subroutine write_intarray_1dkind1

end module utils_outputhdf5
