!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module phantom2hdf5_utils
!
! Utility functions for converting Phantom dumps to HDF5 format
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dump_utils, hdf5
!
 use dump_utils, only:open_dumpfile_r,read_header,read_block_header,ndatatypes, &
                      dump_h,extract,free_header,lentag,lenid,ierr_realsize,&
                      read_global_block_header,i_int,i_int1,i_int2,i_int4,i_int8,&
                      i_real,i_real4,i_real8
 use hdf5
 implicit none

 public :: convert_dump_to_hdf5

 private

contains

!----------------------------------------------------------------------------
!+
!   Convert a Phantom dump file to HDF5 format, while preserving all
!   information in the header and array blocks
!+
!----------------------------------------------------------------------------
subroutine convert_dump_to_hdf5(dumpfile,ierr)
 character(len=*), intent(in)  :: dumpfile
 integer,          intent(out) :: ierr
 integer :: iunit, nblocks, narraylengths
 character(len=lenid) :: fileid
 type(dump_h) :: hdr
 integer(kind=8) :: number8(12)  ! Max 12 array lengths per block
 integer :: nums(ndatatypes,12)  ! ndatatypes data types, max 12 array lengths
 integer :: iblock, iarr, k, i
 character(len=lentag) :: tag
 logical :: iexist,singleprec
 integer(hid_t) :: file_id, header_id, particles_id, sinks_id, block_id
 integer :: hdferr
 integer, allocatable :: int_arr(:)
 integer(kind=1), allocatable :: int1_arr(:)
 integer(kind=2), allocatable :: int2_arr(:)
 integer(kind=4), allocatable :: int4_arr(:)
 integer(kind=8), allocatable :: int8_arr(:)
 real,            allocatable :: real_arr(:)
 real(kind=4),    allocatable :: real4_arr(:)
 real(kind=8),    allocatable :: real8_arr(:)

 ! check if file exists
 inquire(file=trim(dumpfile),exist=iexist)
 if (.not.iexist) then
    print "(1x,a)",'ERROR: dump file '//trim(dumpfile)//' does not exist'
    ierr = 1
    return
 endif

 ! open dump file
 iunit = 12
 call open_dumpfile_r(iunit,dumpfile,fileid,ierr)
 if (ierr == ierr_realsize) then
    close(iunit)
    singleprec = .true.
    call open_dumpfile_r(iunit,dumpfile,fileid,ierr,singleprec=singleprec)
 endif
 if (ierr /= 0) then
    print*,'ERROR: Failed to open dump file'
    ierr = 1
    return
 endif

 print "(1x,a)",'Converting '//trim(dumpfile)//' to HDF5 format...'

 ! initialize HDF5
 call h5open_f(hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to initialize HDF5'
    ierr = 1
    return
 endif

 ! create HDF5 file
 call h5fcreate_f(trim(dumpfile)//'.h5', H5F_ACC_TRUNC_F, file_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 file'
    ierr = 1
    return
 endif

 ! create groups
 call h5gcreate_f(file_id, 'header', header_id, hdferr)
 call h5gcreate_f(file_id, 'particles', particles_id, hdferr)
 call h5gcreate_f(file_id, 'sinks', sinks_id, hdferr)

 ! read and write header
 call read_header(iunit,hdr,ierr)
 if (ierr /= 0) then
    print*,'ERROR: Failed to read header'
    ierr = 1
    return
 endif

 ! write header variables to HDF5
 print "(/,a)",'/header'
 call write_header_to_hdf5(header_id, hdr, hdferr)

 ! read array blocks
 call read_global_block_header(nblocks,narraylengths,hdr,iunit,ierr)
 if (ierr /= 0) return

 ! free header
 call free_header(hdr)

 ! process each block
 do iblock = 1,nblocks
    call read_block_header(narraylengths,number8,nums,iunit,ierr)
    do iarr = 1,narraylengths
       if (iarr == 2) then
          block_id = sinks_id
          print "(/,a)",'/sinks'
       else
          block_id = particles_id
          print "(/,a)",'/particles'
       endif
       do k = 1,ndatatypes  ! loop over data types
          i = 0
          do while (i < nums(k,iarr))
             i = i + 1
             read(iunit,iostat=ierr) tag
             if (.not.is_vector(tag)) print "(a20,i10)",trim(tag),number8(iarr)
             if (ierr /= 0) return

             ! process based on data type
             select case(k)
             case(i_int)  ! default integer
                allocate(int_arr(number8(iarr)))
                read(iunit) int_arr
                call write_array_to_hdf5(block_id, trim(tag), int_arr, hdferr)
                deallocate(int_arr)
             case(i_int1)  ! int*1
                allocate(int1_arr(number8(iarr)))
                read(iunit) int1_arr
                call write_array_to_hdf5(block_id, trim(tag), int1_arr, hdferr)
                deallocate(int1_arr)
             case(i_int2)  ! int*2
                allocate(int2_arr(number8(iarr)))
                read(iunit) int2_arr
                call write_array_to_hdf5(block_id, trim(tag), int2_arr, hdferr)
                deallocate(int2_arr)
             case(i_int4)  ! int*4
                allocate(int4_arr(number8(iarr)))
                read(iunit) int4_arr
                call write_array_to_hdf5(block_id, trim(tag), int4_arr, hdferr)
                deallocate(int4_arr)
             case(i_int8)  ! int*8
                allocate(int8_arr(number8(iarr)))
                read(iunit) int8_arr
                call write_array_to_hdf5(block_id, trim(tag), int8_arr, hdferr)
                deallocate(int8_arr)
             case(i_real)  ! default real
                allocate(real_arr(number8(iarr)))
                if (is_vector(tag)) then
                   call read_vector_components_and_write(iunit,real_arr,block_id,tag,number8(iarr),i,hdferr)
                else
                   read(iunit) real_arr
                   call write_array_to_hdf5(block_id, trim(tag), real_arr, hdferr)
                endif
                deallocate(real_arr)
             case(i_real4)  ! real*4
                allocate(real4_arr(number8(iarr)))
                if (is_vector(tag)) then
                   call read_vector_components_and_write(iunit,real4_arr,block_id,tag,number8(iarr),i,hdferr)
                else
                   read(iunit) real4_arr
                   call write_array_to_hdf5(block_id, trim(tag), real4_arr, hdferr)
                endif
                deallocate(real4_arr)
             case(i_real8)  ! real*8
                allocate(real8_arr(number8(iarr)))
                if (is_vector(tag)) then
                   call read_vector_components_and_write(iunit,real8_arr,block_id,tag,number8(iarr),i,hdferr)
                else
                   read(iunit) real8_arr
                   call write_array_to_hdf5(block_id, trim(tag), real8_arr, hdferr)
                endif
                deallocate(real8_arr)
             end select
          enddo
       enddo
    enddo
 enddo

 ! close HDF5 file and groups
 call h5gclose_f(header_id, hdferr)
 call h5gclose_f(particles_id, hdferr)
 call h5gclose_f(sinks_id, hdferr)
 call h5fclose_f(file_id, hdferr)
 call h5close_f(hdferr)

 close(iunit)

end subroutine convert_dump_to_hdf5

! function to check if a tag represents a vector component
function is_vector(tag) result(isvec)
 character(len=*), intent(in) :: tag
 logical :: isvec

 isvec = (len_trim(tag) >= 1 .and. tag(len_trim(tag):len_trim(tag)) == 'x')

end function is_vector

! read vector components and write as 2D array (polymorphic for real*4 and real*8)
subroutine read_vector_components_and_write(iunit, arr, group_id, tag, n, i, hdferr)
 integer, intent(in) :: iunit
 class(*), intent(inout) :: arr(:)
 integer(hid_t), intent(in) :: group_id
 character(len=*), intent(in) :: tag
 integer(kind=8), intent(in) :: n
 integer, intent(inout) :: i
 integer, intent(out) :: hdferr
 character(len=lentag) :: base_name, tag_tmp
 real(kind=4), allocatable :: vector4_arr(:,:)
 real(kind=8), allocatable :: vector8_arr(:,:)
 integer :: ierr

 ! Determine base name for the vector
 base_name = tag(1:len_trim(tag)-1)//'xyz'

 ! Handle based on type
 select type(arr)
 type is (real(kind=4))
    allocate(vector4_arr(3,n))

    ! Read x component
    read(iunit) vector4_arr(1,:)

    ! Read y component
    read(iunit,iostat=ierr) tag_tmp
    if (ierr /= 0 .or. tag_tmp(len_trim(tag_tmp):len_trim(tag_tmp)) /= 'y') return
    read(iunit) vector4_arr(2,:)

    ! Read z component
    read(iunit,iostat=ierr) tag_tmp
    if (ierr /= 0 .or. tag_tmp(len_trim(tag_tmp):len_trim(tag_tmp)) /= 'z') return
    read(iunit) vector4_arr(3,:)

    i = i + 2

    ! Write combined vector
    call write_2Darray_to_hdf5(group_id, trim(base_name), vector4_arr, hdferr)
    deallocate(vector4_arr)

 type is (real(kind=8))
    allocate(vector8_arr(3,n))

    ! Read x component
    read(iunit) vector8_arr(1,:)

    ! Read y component
    read(iunit,iostat=ierr) tag_tmp
    if (ierr /= 0 .or. tag_tmp(len_trim(tag_tmp):len_trim(tag_tmp)) /= 'y') return
    read(iunit) vector8_arr(2,:)

    ! Read z component
    read(iunit,iostat=ierr) tag_tmp
    if (ierr /= 0 .or. tag_tmp(len_trim(tag_tmp):len_trim(tag_tmp)) /= 'z') return
    read(iunit) vector8_arr(3,:)

    i = i + 2

    ! Write combined vector
    call write_2Darray_to_hdf5(group_id, trim(base_name), vector8_arr, hdferr)
    deallocate(vector8_arr)
 end select

end subroutine read_vector_components_and_write

! function to check for identical tags and return number of matches
function check_for_identical_tags(tags, start_idx, ntags) result(nmatches)
 character(len=lentag), intent(in) :: tags(:)
 integer, intent(in) :: start_idx, ntags
 integer :: nmatches
 character(len=lentag) :: current_tag
 integer :: i

 current_tag = trim(tags(start_idx))
 nmatches = 1
 do i = start_idx + 1, ntags
    if (trim(tags(i)) == trim(current_tag)) then
       nmatches = nmatches + 1
    else
       exit
    endif
 enddo

end function check_for_identical_tags

subroutine write_header_to_hdf5(group_id, hdr, hdferr)
 integer(hid_t), intent(in)  :: group_id
 type(dump_h),   intent(in)  :: hdr
 integer,        intent(out) :: hdferr
 integer :: i, nmatches

 ! write each header variable to HDF5
 do i = 1,hdr%nums(i_int)  ! default integers
    if (i > 1) then
       if (trim(hdr%inttags(i)) == trim(hdr%inttags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%inttags, i, hdr%nums(1))
    call write_to_hdf5(group_id, trim(hdr%inttags(i)), hdr%intvals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_int1)  ! int*1
    if (i > 1) then
       if (trim(hdr%int1tags(i)) == trim(hdr%int1tags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%int1tags, i, hdr%nums(2))
    call write_to_hdf5(group_id, trim(hdr%int1tags(i)), hdr%int1vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_int2)  ! int*2
    if (i > 1) then
       if (trim(hdr%int2tags(i)) == trim(hdr%int2tags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%int2tags, i, hdr%nums(3))
    call write_to_hdf5(group_id, trim(hdr%int2tags(i)), hdr%int2vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_int4)  ! int*4
    if (i > 1) then
       if (trim(hdr%int4tags(i)) == trim(hdr%int4tags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%int4tags, i, hdr%nums(4))
    call write_to_hdf5(group_id, trim(hdr%int4tags(i)), hdr%int4vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_int8)  ! int*8
    if (i > 1) then
       if (trim(hdr%int8tags(i)) == trim(hdr%int8tags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%int8tags, i, hdr%nums(5))
    call write_to_hdf5(group_id, trim(hdr%int8tags(i)), hdr%int8vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_real)  ! default reals
    if (i > 1) then
       if (trim(hdr%realtags(i)) == trim(hdr%realtags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%realtags, i, hdr%nums(6))
    call write_to_hdf5(group_id, trim(hdr%realtags(i)), hdr%realvals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_real4)  ! real*4
    if (i > 1) then
       if (trim(hdr%real4tags(i)) == trim(hdr%real4tags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%real4tags, i, hdr%nums(7))
    call write_to_hdf5(group_id, trim(hdr%real4tags(i)), hdr%real4vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(i_real8)  ! real*8
    if (i > 1) then
       if (trim(hdr%real8tags(i)) == trim(hdr%real8tags(i-1))) cycle
    endif
    nmatches = check_for_identical_tags(hdr%real8tags, i, hdr%nums(8))
    call write_to_hdf5(group_id, trim(hdr%real8tags(i)), hdr%real8vals(i:i+nmatches-1), hdferr)
 enddo

end subroutine write_header_to_hdf5

! generic subroutine to write a scalar to HDF5
subroutine write_scalar_to_hdf5(group_id, name, value, hdferr)
 integer(hid_t), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value
 integer, intent(out) :: hdferr
 integer(hid_t) :: dspace_id, dset_id
 integer(hsize_t) :: dims(1) = (/1/)

 ! create scalar dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
 if (hdferr /= 0) return

 ! create and write dataset based on type
 select type(value)
 type is (integer(kind=1))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I8LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I8LE, value, dims, hdferr)
 type is (integer(kind=2))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I16LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I16LE, value, dims, hdferr)
 type is (integer(kind=4))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I32LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I32LE, value, dims, hdferr)
 type is (integer(kind=8))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I64LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I64LE, value, dims, hdferr)
 type is (real(kind=4))
    call h5dcreate_f(group_id, trim(name), H5T_IEEE_F32LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_IEEE_F32LE, value, dims, hdferr)
 type is (real(kind=8))
    call h5dcreate_f(group_id, trim(name), H5T_IEEE_F64LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_IEEE_F64LE, value, dims, hdferr)
 end select

 ! Close dataset and dataspace
 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(dspace_id, hdferr)

end subroutine write_scalar_to_hdf5

! generic subroutine to write an array to HDF5
subroutine write_array_to_hdf5(group_id, name, value, hdferr)
 integer(hid_t), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value(:)
 integer, intent(out) :: hdferr
 integer(hid_t) :: dspace_id, dset_id
 integer(hsize_t) :: dims(1)

 ! create dataspace
 dims(1) = size(value)
 call h5screate_simple_f(1, dims, dspace_id, hdferr)
 if (hdferr /= 0) return

 ! create and write dataset based on type
 select type(value)
 type is (integer(kind=1))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I8LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I8LE, value, dims, hdferr)
 type is (integer(kind=2))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I16LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I16LE, value, dims, hdferr)
 type is (integer(kind=4))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I32LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I32LE, value, dims, hdferr)
 type is (integer(kind=8))
    call h5dcreate_f(group_id, trim(name), H5T_STD_I64LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_STD_I64LE, value, dims, hdferr)
 type is (real(kind=4))
    call h5dcreate_f(group_id, trim(name), H5T_IEEE_F32LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_IEEE_F32LE, value, dims, hdferr)
 type is (real(kind=8))
    call h5dcreate_f(group_id, trim(name), H5T_IEEE_F64LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_IEEE_F64LE, value, dims, hdferr)
 end select

 ! close dataset and dataspace
 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(dspace_id, hdferr)

end subroutine write_array_to_hdf5

subroutine write_to_hdf5(group_id, name, value, hdferr)
 integer(hid_t), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value(:)
 integer, intent(out) :: hdferr
 character(len=len(name)+4) :: dataset_name

 dataset_name = trim(name)
 ! append _i8 suffix for int*8 arrays so they don't conflict with int arrays
 select type(value)
 type is (integer(kind=8))
    dataset_name = trim(name)//'_i8'
 end select
 print "(a20,i10)",dataset_name,size(value)

 ! If array has dimension 1, treat as scalar
 if (size(value) == 1) then
    call write_scalar_to_hdf5(group_id, dataset_name, value(1), hdferr)
 else
    ! Otherwise treat as array
    call write_array_to_hdf5(group_id, dataset_name, value, hdferr)
 endif

end subroutine write_to_hdf5

! generic subroutine to write a 2D array to HDF5
subroutine write_2Darray_to_hdf5(group_id, name, value, hdferr)
 integer(hid_t), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value(:,:)
 integer, intent(out) :: hdferr
 integer(hid_t) :: dspace_id, dset_id
 integer(hsize_t) :: dims(2)

 print "(a20,i10,i10)",name,size(value,1),size(value,2)

 ! create dataspace
 dims(1) = size(value,1)
 dims(2) = size(value,2)
 call h5screate_simple_f(2, dims, dspace_id, hdferr)
 if (hdferr /= 0) return

 ! create and write dataset based on type
 select type(value)
 type is (real(kind=4))
    call h5dcreate_f(group_id, trim(name), H5T_IEEE_F32LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_IEEE_F32LE, value, dims, hdferr)
 type is (real(kind=8))
    call h5dcreate_f(group_id, trim(name), H5T_IEEE_F64LE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_IEEE_F64LE, value, dims, hdferr)
 end select

 ! close dataset and dataspace
 call h5dclose_f(dset_id, hdferr)
 call h5sclose_f(dspace_id, hdferr)

end subroutine write_2Darray_to_hdf5

end module phantom2hdf5_utils
