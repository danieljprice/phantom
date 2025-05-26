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
 use dump_utils, only:open_dumpfile_r,read_header,read_block_header,ndatatypes, &
                       dump_h,extract,free_header,lentag,lenid,ierr_realsize,&
                       read_global_block_header
 use hdf5
 implicit none

 public :: convert_dump_to_hdf5

 private

contains

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
 integer(HID_T) :: file_id, header_id, particles_id, sinks_id, block_id
 integer :: hdferr
 integer, allocatable :: int_arr(:)
 integer(kind=1), allocatable :: int1_arr(:)
 integer(kind=2), allocatable :: int2_arr(:)
 integer(kind=4), allocatable :: int4_arr(:)
 integer(kind=8), allocatable :: int8_arr(:)
 real,            allocatable :: real_arr(:)
 real(kind=4),    allocatable :: real4_arr(:)
 real(kind=8),    allocatable :: real8_arr(:)

 ! Check if file exists
 inquire(file=trim(dumpfile),exist=iexist)
 if (.not.iexist) then
    print "(1x,a)",'ERROR: dump file '//trim(dumpfile)//' does not exist'
    ierr = 1
    return
 endif

 ! Open dump file
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

 ! Initialize HDF5
 call h5open_f(hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to initialize HDF5'
    ierr = 1
    return
 endif

 ! Create HDF5 file
 call h5fcreate_f(trim(dumpfile)//'.h5', H5F_ACC_TRUNC_F, file_id, hdferr)
 if (hdferr /= 0) then
    print*,'ERROR: Failed to create HDF5 file'
    ierr = 1
    return
 endif

 ! Create groups
 call h5gcreate_f(file_id, 'header', header_id, hdferr)
 call h5gcreate_f(file_id, 'particles', particles_id, hdferr)
 call h5gcreate_f(file_id, 'sinks', sinks_id, hdferr)

 ! Read and write header
 call read_header(iunit,hdr,ierr)
 if (ierr /= 0) then
    print*,'ERROR: Failed to read header'
    ierr = 1
    return
 endif

 ! Write header variables to HDF5
 print "(/,a)",'/header'
 call write_header_to_hdf5(header_id, hdr, hdferr)

 ! Read array blocks
 call read_global_block_header(nblocks,narraylengths,hdr,iunit,ierr)
 if (ierr /= 0) return

 ! Free header
 call free_header(hdr)

 ! Process each block
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
       do k = 1,ndatatypes  ! Loop over data types
          do i = 1,nums(k,iarr)
             read(iunit,iostat=ierr) tag
             print "(a20,i10)",trim(tag),number8(iarr)
             if (ierr /= 0) return

             ! Read and write array based on type
             select case(k)
             case(1)  ! Default integer
                allocate(int_arr(number8(iarr)))
                read(iunit) int_arr
                call write_array_to_hdf5(block_id, trim(tag), int_arr, hdferr)
                deallocate(int_arr)
             case(2)  ! Integer*1
                allocate(int1_arr(number8(iarr)))
                read(iunit) int1_arr
                call write_array_to_hdf5(block_id, trim(tag), int1_arr, hdferr)
                deallocate(int1_arr)
             case(3)  ! Integer*2
                allocate(int2_arr(number8(iarr)))
                read(iunit) int2_arr
                call write_array_to_hdf5(block_id, trim(tag), int2_arr, hdferr)
                deallocate(int2_arr)
             case(4)  ! Integer*4
                allocate(int4_arr(number8(iarr)))
                read(iunit) int4_arr
                call write_array_to_hdf5(block_id, trim(tag), int4_arr, hdferr)
                deallocate(int4_arr)
             case(5)  ! Integer*8
                allocate(int8_arr(number8(iarr)))
                read(iunit) int8_arr
                call write_array_to_hdf5(block_id, trim(tag), int8_arr, hdferr)
                deallocate(int8_arr)
             case(6)  ! Default real
                allocate(real_arr(number8(iarr)))
                read(iunit) real_arr
                call write_array_to_hdf5(block_id, trim(tag), real_arr, hdferr)
                deallocate(real_arr)
             case(7)  ! Real*4
                allocate(real4_arr(number8(iarr)))
                read(iunit) real4_arr
                call write_array_to_hdf5(block_id, trim(tag), real4_arr, hdferr)
                deallocate(real4_arr)
             case(8)  ! Real*8
                allocate(real8_arr(number8(iarr)))
                read(iunit) real8_arr
                call write_array_to_hdf5(block_id, trim(tag), real8_arr, hdferr)
                deallocate(real8_arr)
             end select
          enddo
       enddo
    enddo
 enddo

 ! Close HDF5 file and groups
 call h5gclose_f(header_id, hdferr)
 call h5gclose_f(particles_id, hdferr)
 call h5gclose_f(sinks_id, hdferr)
 call h5fclose_f(file_id, hdferr)
 call h5close_f(hdferr)

 close(iunit)

end subroutine convert_dump_to_hdf5

! Function to check for identical tags and return number of matches
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
 integer(HID_T), intent(in)  :: group_id
 type(dump_h),   intent(in)  :: hdr
 integer,        intent(out) :: hdferr
 integer :: i, nmatches

 ! Write each header variable to HDF5
 do i = 1,hdr%nums(1)  ! Default integers
    if (i > 1 .and. trim(hdr%inttags(i)) == trim(hdr%inttags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%inttags, i, hdr%nums(1))
    call write_to_hdf5(group_id, trim(hdr%inttags(i)), hdr%intvals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(2)  ! Integer*1
    if (i > 1 .and. trim(hdr%int1tags(i)) == trim(hdr%int1tags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%int1tags, i, hdr%nums(2))
    call write_to_hdf5(group_id, trim(hdr%int1tags(i)), hdr%int1vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(3)  ! Integer*2
    if (i > 1 .and. trim(hdr%int2tags(i)) == trim(hdr%int2tags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%int2tags, i, hdr%nums(3))
    call write_to_hdf5(group_id, trim(hdr%int2tags(i)), hdr%int2vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(4)  ! Integer*4
    if (i > 1 .and. trim(hdr%int4tags(i)) == trim(hdr%int4tags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%int4tags, i, hdr%nums(4))
    call write_to_hdf5(group_id, trim(hdr%int4tags(i)), hdr%int4vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(5)  ! Integer*8
    if (i > 1 .and. trim(hdr%int8tags(i)) == trim(hdr%int8tags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%int8tags, i, hdr%nums(5))
    call write_to_hdf5(group_id, trim(hdr%int8tags(i)), hdr%int8vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(6)  ! Default reals
    if (i > 1 .and. trim(hdr%realtags(i)) == trim(hdr%realtags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%realtags, i, hdr%nums(6))
    call write_to_hdf5(group_id, trim(hdr%realtags(i)), hdr%realvals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(7)  ! Real*4
    if (i > 1 .and. trim(hdr%real4tags(i)) == trim(hdr%real4tags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%real4tags, i, hdr%nums(7))
    call write_to_hdf5(group_id, trim(hdr%real4tags(i)), hdr%real4vals(i:i+nmatches-1), hdferr)
 enddo

 do i = 1,hdr%nums(8)  ! Real*8
    if (i > 1 .and. trim(hdr%real8tags(i)) == trim(hdr%real8tags(i-1))) cycle
    nmatches = check_for_identical_tags(hdr%real8tags, i, hdr%nums(8))
    call write_to_hdf5(group_id, trim(hdr%real8tags(i)), hdr%real8vals(i:i+nmatches-1), hdferr)
 enddo

end subroutine write_header_to_hdf5

! Generic subroutine to write a scalar to HDF5
subroutine write_scalar_to_hdf5(group_id, name, value, hdferr)
 integer(HID_T), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value
 integer, intent(out) :: hdferr
 integer(HID_T) :: dspace_id, dset_id
 integer(HSIZE_T) :: dims(1) = (/1/)

 ! Create scalar dataspace
 call h5screate_f(H5S_SCALAR_F, dspace_id, hdferr)
 if (hdferr /= 0) return

 ! Create and write dataset based on type
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

! Generic subroutine to write an array to HDF5
subroutine write_array_to_hdf5(group_id, name, value, hdferr)
 integer(HID_T), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value(:)
 integer, intent(out) :: hdferr
 integer(HID_T) :: dspace_id, dset_id
 integer(HSIZE_T) :: dims(1)

 ! Create dataspace
 dims(1) = size(value)
 call h5screate_simple_f(1, dims, dspace_id, hdferr)
 if (hdferr /= 0) return

 ! Create and write dataset based on type
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

end subroutine write_array_to_hdf5

subroutine write_to_hdf5(group_id, name, value, hdferr)
 integer(HID_T), intent(in) :: group_id
 character(len=*), intent(in) :: name
 class(*), intent(in) :: value(:)
 integer, intent(out) :: hdferr

 print "(a20,i10)",name,size(value)

 ! If array has dimension 1, treat as scalar
 if (size(value) == 1) then
    call write_scalar_to_hdf5(group_id, name, value(1), hdferr)
 else
    ! Otherwise treat as array
    call write_array_to_hdf5(group_id, name, value, hdferr)
 endif

end subroutine write_to_hdf5

end module phantom2hdf5_utils