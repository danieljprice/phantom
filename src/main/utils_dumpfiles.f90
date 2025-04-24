!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dump_utils
!
! Utility routines used when reading and writing the
!   sphNG/Phantom dump file format
!
!  "Every complex file format eventually turns into a
!    badly-designed programming language." - Anon
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: open_dumpfile_w, open_dumpfile_r, get_error_text
 public :: tag,check_tag,match_tag
 public :: skipblock,skip_arrays,skip_headerblock
 public :: get_dumpname,get_dump_size
 public :: add_to_header,add_to_rheader,add_to_iheader
 public :: num_in_header,reset_header,extract
 public :: read_array_from_file
 public :: write_block_header, write_array
 public :: read_block_header, read_array
 public :: print_arrays_in_file
 integer, parameter, public :: lentag = 16    ! tag length
 integer, parameter, public :: lenid  = 100
 integer, parameter, public :: maxphead = 256
 !
 ! Format version. Change this ONLY if the dump format is incompatible
 ! with previous formats
 !
 integer, parameter, public :: iversion = 1

 ! magic numbers
 integer(kind=4), parameter, public :: int1=060769,int2=060878
 integer(kind=4), parameter, public :: int1o=690706,int2o=780806

 ! data types
 integer, parameter, public :: ndatatypes = 8
 integer, parameter, public :: i_int   = 1, &
                               i_int1  = 2, &
                               i_int2  = 3, &
                               i_int4  = 4, &
                               i_int8  = 5, &
                               i_real  = 6, &
                               i_real4 = 7, &
                               i_real8 = 8
 character(len=*), parameter  :: datatype_label(ndatatypes) = &
     (/'integer', &
       'int*1  ', &
       'int*2  ', &
       'int*4  ', &
       'int*8  ', &
       'real   ', &
       'real*4 ', &
       'double '/)

 ! error codes
 integer, parameter, public :: ierr_fileopen  = 1,&
                               ierr_endian    = 2,&
                               ierr_version   = 3,&
                               ierr_realsize  = 4,&
                               ierr_intsize   = 5,&
                               ierr_notags    = 6,&
                               ierr_unknown   = 7,&
                               ierr_notenough = 8,&
                               ierr_arraysize = 9

 type dump_h
    integer :: nums(ndatatypes)
    character(len=lentag), allocatable :: inttags(:)
    character(len=lentag), allocatable :: int1tags(:)
    character(len=lentag), allocatable :: int2tags(:)
    character(len=lentag), allocatable :: int4tags(:)
    character(len=lentag), allocatable :: int8tags(:)
    character(len=lentag), allocatable :: realtags(:)
    character(len=lentag), allocatable :: real4tags(:)
    character(len=lentag), allocatable :: real8tags(:)
    integer,         allocatable :: intvals(:)
    integer(kind=1), allocatable :: int1vals(:)
    integer(kind=2), allocatable :: int2vals(:)
    integer(kind=4), allocatable :: int4vals(:)
    integer(kind=8), allocatable :: int8vals(:)
    real,            allocatable :: realvals(:)
    real(kind=4),    allocatable :: real4vals(:)
    real(kind=8),    allocatable :: real8vals(:)
 end type dump_h

 public :: dump_h

 public :: write_header, read_header
 public :: allocate_header, free_header
 public :: print_header
 public :: get_blocklimits

 ! generic interface to extract quantities from header
 interface extract
  module procedure extract_int4, extract_int8, &
                   extract_real4, extract_real8, &
                   extract_int4arr, extract_int8arr, &
                   extract_real4arr, extract_real8arr
  module procedure extracthdr_int4, extracthdr_int8, &
   extracthdr_real4, extracthdr_real8, &
   extracthdr_int4arr, extracthdr_int8arr, &
   extracthdr_real4arr, extracthdr_real8arr
 end interface extract

 ! generic interface for writing values to header
 interface add_to_header
  module procedure add_to_header_int4, add_to_header_int8, &
    add_to_header_int4arr,  add_to_header_int8arr, &
    add_to_header_real4,    add_to_header_real8, &
    add_to_header_real4arr, add_to_header_real8arr
 end interface add_to_header

 ! add to the default real section of the header
 interface add_to_rheader
  module procedure add_to_rheader, add_to_rheader_arr
 end interface add_to_rheader

 ! add to the default int section of the header
 interface add_to_iheader
  module procedure add_to_iheader, add_to_iheader_arr
 end interface add_to_iheader

 ! generic interface for reset of header
 interface reset_header
  module procedure reset_header_real
 end interface reset_header

 ! generic interface for writing arrays to file
 interface write_array
  module procedure write_array_int1, &
   write_array_int4,  write_array_int8, &
   write_array_real4, write_array_real4arr, &
   write_array_real8, write_array_real8arr, &
   write_array_int4arr
 end interface write_array

 ! generic interface for writing arrays to file
 interface read_array
  module procedure read_array_int1, &
   read_array_int4, read_array_int8, &
   read_array_real4, read_array_real4arr, &
   read_array_real8, read_array_real8arr, &
   read_array_int4arr
 end interface read_array

 ! generic interface for reading arrays from dumpfile
 interface read_array_from_file
  module procedure read_array_from_file_r4, read_array_from_file_r8
 end interface read_array_from_file

 private

contains
!--------------------------------------------------------------------
!+
!  filenaming convention for split MPI dumps
!+
!--------------------------------------------------------------------
function get_dumpname(filename,id)
 character(len=*), intent(in) :: filename
 character(len=len_trim(filename)+8) :: get_dumpname
 integer,          intent(in) :: id

 write(get_dumpname,"(a,a5,i3.3)") trim(filename),'_part',id+1

end function get_dumpname

!--------------------------------------------------------------------
!+
!  extract dump size (full or small) from the fileid string
!+
!--------------------------------------------------------------------
subroutine get_dump_size(fileid,smalldump)
 character(len=lenid), intent(in)  :: fileid
 logical,              intent(out) :: smalldump
 !
 if (fileid(1:1)=='S') then
    smalldump = .true.
 else
    smalldump = .false.
 endif

end subroutine get_dump_size

!--------------------------------------------------------------------
!+
!  small utility to skip an entire block in a file
!+
!-------------------------------------------------------------------
subroutine skipblock(iunit,nums1,nums2,nums3,nums4,tagged,ierr)
 integer, intent(in)  :: iunit
 integer, intent(in)  :: nums1(ndatatypes),nums2(ndatatypes),nums3(ndatatypes),nums4(ndatatypes)
 logical, intent(in)  :: tagged
 integer, intent(out) :: ierr
 integer :: nskip

 ierr = 0
 nskip = sum(nums1)+sum(nums2)+sum(nums3)+sum(nums4)
 call skip_arrays(iunit,nskip,tagged,ierr)

end subroutine skipblock

!--------------------------------------------------------------------

subroutine skip_arrays(iunit,nskip,tagged,ierr,verbose)
 integer, intent(in)  :: iunit, nskip
 logical, intent(in)  :: tagged
 integer, intent(out) :: ierr
 integer, intent(in), optional :: verbose
 character(len=lentag) :: tag
 integer :: i

 ierr = 0
 do i=1,nskip
    if (tagged) read(iunit,iostat=ierr) tag
    if (present(verbose)) then
       if (verbose >= 2) print*,' skipping '//trim(tag)
    endif
    read(iunit,iostat=ierr)
 enddo

end subroutine skip_arrays

!--------------------------------------------------------------------
!+
!  small utility to skip the whole single variable header
!+
!-------------------------------------------------------------------
subroutine skip_headerblock(iunit,ierr)
 integer, intent(in)  :: iunit
 integer, intent(out) :: ierr
 integer :: number

 read(iunit,iostat=ierr) number
 if (ierr /= 0) return
 if (number > 0) then
    read(iunit,iostat=ierr) ! skip tags
    if (ierr /= 0) return
    read(iunit,iostat=ierr) ! skip variables
    if (ierr /= 0) return
 endif

end subroutine skip_headerblock

!---------------------------------------------------------------------
!+
! Construct tag, padded or truncated to 16-characters based on
!  input string
!+
!---------------------------------------------------------------------
elemental function tag(label)
 character(len=lentag) :: tag
 character(len=*), intent(in) :: label

 tag = adjustl(label)

end function tag

!-------------------------------------------
!+
! Check tag against an expected value
!+
!-------------------------------------------
subroutine check_tag(tag,expectedtag)
 character(len=*), intent(in) :: tag, expectedtag

 if (trim(tag)/=trim(expectedtag)) then
    print "(a)",' ERROR reading file: expecting '//trim(expectedtag)//' but got '//trim(tag)
 endif

end subroutine check_tag

!-------------------------------------------
!+
!  check if two tags match
!+
!-------------------------------------------
logical function match_tag(tag,expectedtag)
 character(len=*), intent(in) :: tag, expectedtag

 match_tag = (trim(tag)==trim(expectedtag))

end function match_tag

!------------------------------------------------
!+
!  Read single value from the header structure
!  For default int and default real we first
!  check the default header, but then check
!  the other headers of the same type
!+
!------------------------------------------------
subroutine extracthdr_int4(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 integer(kind=4),  intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer(kind=4),  intent(in), optional :: default
 integer(kind=4) :: def
 integer :: ival,idef

 if (present(default)) then
    def = default
    idef = default
 else
    def = 0
    idef = 0
 endif
 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.,default=idef)
    if (ierr==0) then
       val = ival
    else
       call extract_int4(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr,default=def)
    endif
 else
    call extract_int4(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr,default=def)
 endif

end subroutine extracthdr_int4

!-----------------------------------------------
! extraction of single int*8 value from header
!-----------------------------------------------
subroutine extracthdr_int8(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 integer(kind=8),  intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer(kind=8),  intent(in), optional :: default
 integer(kind=8) :: def
 integer :: ival,idef

 if (present(default)) then
    def = default
    idef = int(default)
 else
    def = 0_8
    idef = 0
 endif
 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.,default=idef)
    if (ierr==0) then
       val = ival
    else
       call extract_int8(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr,default=def)
    endif
 else
    call extract_int8(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr,default=def)
 endif

end subroutine extracthdr_int8

!-----------------------------------------------
! extraction of single real*4 value from header
!-----------------------------------------------
subroutine extracthdr_real4(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 real(kind=4),     intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real(kind=4),     intent(in), optional :: default
 real(kind=4) :: def
 real :: rval,rdef

 if (present(default)) then
    def = default
    rdef = default
 else
    def = 0.d0
    rdef = 0.
 endif
 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,default=rdef,q=.true.)
    if (ierr==0) then
       val = real(rval,kind=4)
    else
       call extract_real4(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr,default=def)
    endif
 else
    call extract_real4(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr,default=def)
 endif
end subroutine extracthdr_real4

!-----------------------------------------------
! extraction of single real*8 value from header
!-----------------------------------------------
subroutine extracthdr_real8(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 real(kind=8),     intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real(kind=8),     intent(in), optional :: default
 real(kind=8) :: def
 real :: rval,rdef

 if (present(default)) then
    def = default
    rdef = real(default)
 else
    def = 0.d0
    rdef = 0.
 endif
 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,q=.true.,default=rdef)
    if (ierr==0) then
       val = rval
    else
       call extract_real8(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr,default=def)
    endif
 else
    call extract_real8(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr,default=def)
 endif

end subroutine extracthdr_real8

!------------------------------------------------
!+
!  Read array from the header structure
!+
!------------------------------------------------
subroutine extracthdr_int4arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 integer(kind=4),  intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer :: ival(size(val))

 ierr = 0
 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.)
    if (ierr==0) then
       val = ival
    else
       call extract_int4arr(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr)
    endif
 else
    call extract_int4arr(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr)
 endif

end subroutine extracthdr_int4arr

!------------------------------------------
! extraction of int*8 arrays from header
!------------------------------------------
subroutine extracthdr_int8arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 integer(kind=8),  intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer :: ival(size(val))

 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.)
    if (ierr==0) then
       val = ival
    else
       call extract_int8arr(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr)
    endif
 else
    call extract_int8arr(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr)
 endif

end subroutine extracthdr_int8arr

!------------------------------------------
! extraction of real*4 arrays from header
!------------------------------------------
subroutine extracthdr_real4arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 real(kind=4),     intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real :: rval(size(val))

 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,q=.true.)
    if (ierr==0) then
       val = real(rval,kind=4)
    else
       call extract_real4arr(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr)
    endif
 else
    call extract_real4arr(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr)
 endif
end subroutine extracthdr_real4arr

!------------------------------------------
! extraction of real*8 arrays from header
!------------------------------------------
subroutine extracthdr_real8arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 real(kind=8),     intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real :: rval(size(val))

 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,q=.true.)
    if (ierr==0 .or. ierr==ierr_notenough) then
       val = rval
    else
       call extract_real8arr(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr)
    endif
 else
    call extract_real8arr(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr)
 endif

end subroutine extracthdr_real8arr

!-------------------------------------------
! Extraction of int*4 variables from header
!-------------------------------------------
subroutine extract_int4(tag,ival,intarr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=4),       intent(out) :: ival
 integer,               intent(in)  :: ntags
 integer(kind=4),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 integer(kind=4),       intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ! default if not found
 if (present(default)) then
    ival = default
 else
    ival = 0
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i) then
          ival = intarr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int4

!-------------------------------------------
! Extraction of int*8 variables from header
!-------------------------------------------
subroutine extract_int8(tag,ival,intarr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=8),       intent(out) :: ival
 integer,               intent(in)  :: ntags
 integer(kind=8),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 integer(kind=8),       intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ! default if not found
 if (present(default)) then
    ival = default
 else
    ival = 0_8
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i) then
          ival = intarr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int8

!---------------------------------------------------
! Extraction of single real*8 variables from header
!---------------------------------------------------
subroutine extract_real8(tag,rval,r8arr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 real(kind=8),          intent(out) :: rval
 real(kind=8),          intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 real(kind=8),          intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 if (present(default)) then
    rval = default
 else
    rval = 0.d0 ! default if not found
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r8arr) >= i) then
          rval = r8arr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real8

!---------------------------------------------------
! Extraction of single real*4 variables from header
!---------------------------------------------------
subroutine extract_real4(tag,rval,r4arr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 real(kind=4),          intent(out) :: rval
 real(kind=4),          intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 real(kind=4),          intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 if (present(default)) then
    rval = default
 else
    rval = 0. ! default if not found
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r4arr) >= i) then
          rval = r4arr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real4

!------------------------------------------
! Extraction of int*4 arrays frmo header
!------------------------------------------
subroutine extract_int4arr(tag,ival,intarr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=4),       intent(out) :: ival(:)
 integer,               intent(in)  :: ntags
 integer(kind=4),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 ival(:) = 0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i .and. size(ival) > nmatched) then
          nmatched = nmatched + 1
          ival(nmatched) = intarr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(ival)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int4arr

!------------------------------------------
! Extraction of int*8 arrays from header
!------------------------------------------
subroutine extract_int8arr(tag,ival,intarr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=8),       intent(out) :: ival(:)
 integer,               intent(in)  :: ntags
 integer(kind=8),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 ival(:) = 0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i .and. size(ival) > nmatched) then
          nmatched = nmatched + 1
          ival(nmatched) = intarr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(ival)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int8arr

!------------------------------------------
! Extraction of real*8 arrays from header
!------------------------------------------
subroutine extract_real8arr(tag,rval,r8arr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 real(kind=8),          intent(out) :: rval(:)
 real(kind=8),          intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 rval = 0.d0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r8arr) >= i .and. size(rval) > nmatched) then
          nmatched = nmatched + 1
          rval(nmatched) = r8arr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(rval)) then
    ierr = 0
 elseif (nmatched > 0) then
    ierr = ierr_notenough
 endif
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real8arr

!------------------------------------------
! extraction of real*4 arrays from header
!------------------------------------------
subroutine extract_real4arr(tag,rval,r4arr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 real(kind=4),          intent(out) :: rval(:)
 real(kind=4),          intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 rval = 0. ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r4arr) >= i .and. size(rval) > nmatched) then
          nmatched = nmatched + 1
          rval(nmatched) = r4arr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(rval)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real4arr

!------------------------------------------
! reset header to all blank entries
!------------------------------------------
subroutine reset_header_real(rheader,rtags)
 real,                  intent(out)   :: rheader(:)
 character(len=lentag), intent(inout) :: rtags(:)
 integer :: i

 do i=1,size(rheader)
    rheader(i) = 0.
    rtags(i)   = ''
 enddo

end subroutine reset_header_real

!------------------------------------------
! add item to int header
!------------------------------------------
subroutine add_to_header_int4(ival,tag,hdr,ierr)
 integer(kind=4),  intent(in)    :: ival
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_int4) + 1
 if (i > size(hdr%int4vals)) then
    ierr = 1
 else
    hdr%nums(i_int4) = i
    hdr%int4vals(i)  = ival
    hdr%int4tags(i)  = tag
 endif

end subroutine add_to_header_int4

!------------------------------------------
! add item to int header
!------------------------------------------
subroutine add_to_header_int8(ival,tag,hdr,ierr)
 integer(kind=8),  intent(in)    :: ival
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_int8) + 1
 if (i > size(hdr%int8vals)) then
    ierr = 1
 else
    hdr%nums(i_int8) = i
    hdr%int8vals(i)  = ival
    hdr%int8tags(i)  = tag
 endif

end subroutine add_to_header_int8

!------------------------------------------
! add array to integer header
!------------------------------------------
subroutine add_to_header_int4arr(ival,tag,hdr,ierr)
 integer(kind=4),  intent(in)    :: ival(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_int4) + 1
 do j=1,size(ival)
    if (i < size(hdr%int4vals)) then
       hdr%nums(i_int4) = i
       hdr%int4vals(i) = ival(j)
       hdr%int4tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_int4arr

!------------------------------------------
! add array to integer*8 header
!------------------------------------------
subroutine add_to_header_int8arr(ival,tag,hdr,ierr)
 integer(kind=8),  intent(in)    :: ival(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_int8) + 1
 do j=1,size(ival)
    if (i < size(hdr%int8vals)) then
       hdr%nums(i_int8) = i
       hdr%int8vals(i) = ival(j)
       hdr%int8tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_int8arr

!------------------------------------------
! add item to real*4 header
!------------------------------------------
subroutine add_to_header_real4(rval,tag,hdr,ierr)
 real(kind=4),     intent(in)    :: rval
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_real4) + 1
 if (i > size(hdr%real4vals)) then
    ierr = 1
 else
    hdr%nums(i_real4) = i
    hdr%real4vals(i)  = rval
    hdr%real4tags(i)  = tag
 endif

end subroutine add_to_header_real4

!------------------------------------------
! add item to real*8 header
!------------------------------------------
subroutine add_to_header_real8(rval,tag,hdr,ierr)
 real(kind=8),     intent(in)    :: rval
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_real8) + 1
 if (i > size(hdr%real8vals)) then
    ierr = 1
 else
    hdr%nums(i_real8) = i
    hdr%real8vals(i)  = rval
    hdr%real8tags(i)  = tag
 endif

end subroutine add_to_header_real8

!------------------------------------------
! add array to real*4 header
!------------------------------------------
subroutine add_to_header_real4arr(rval,tag,hdr,ierr)
 real(kind=4),     intent(in)    :: rval(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_real4) + 1
 do j=1,size(rval)
    if (i < size(hdr%real4vals)) then
       hdr%nums(i_real4) = i
       hdr%real4vals(i) = rval(j)
       hdr%real4tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_real4arr

!------------------------------------------
! add array to real*8 header
!------------------------------------------
subroutine add_to_header_real8arr(rval,tag,hdr,ierr)
 real(kind=8),     intent(in)    :: rval(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_real8) + 1
 do j=1,size(rval)
    if (i < size(hdr%real8vals)) then
       hdr%nums(i_real8) = i
       hdr%real8vals(i) = rval(j)
       hdr%real8tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_real8arr

!------------------------------------------
! add item to default real header
!------------------------------------------
subroutine add_to_rheader(rval,tag,hdr,ierr)
 real,             intent(in)    :: rval
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_real) + 1
 if (i > size(hdr%realvals)) then
    ierr = 1
 else
    hdr%nums(i_real) = i
    hdr%realvals(i)  = rval
    hdr%realtags(i)  = tag
 endif

end subroutine add_to_rheader

!------------------------------------------
! add array to default real header
!------------------------------------------
subroutine add_to_rheader_arr(rval,tag,hdr,ierr)
 real,             intent(in)    :: rval(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_real) + 1
 do j=1,size(rval)
    if (i < size(hdr%realvals)) then
       hdr%nums(i_real) = i
       hdr%realvals(i) = rval(j)
       hdr%realtags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_rheader_arr

!------------------------------------------
! add item to default int header
!------------------------------------------
subroutine add_to_iheader(ival,tag,hdr,ierr)
 integer,          intent(in)    :: ival
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_int) + 1
 if (i > size(hdr%intvals)) then
    ierr = 1
 else
    hdr%nums(i_int) = i
    hdr%intvals(i)  = ival
    hdr%inttags(i)  = tag
 endif

end subroutine add_to_iheader

!------------------------------------------
! add array to default int header
!------------------------------------------
subroutine add_to_iheader_arr(ival,tag,hdr,ierr)
 integer,          intent(in)    :: ival(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_int) + 1
 do j=1,size(ival)
    if (i < size(hdr%intvals)) then
       hdr%nums(i_int) = i
       hdr%intvals(i)  = ival(j)
       hdr%inttags(i)  = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_iheader_arr

!------------------------------------------
! add item to real header
!------------------------------------------
integer function num_in_header(tags)
 character(len=lentag), intent(in) :: tags(:)
 integer :: i

 ! cycle through header backwards until non-blank value is found
 i = size(tags)
 do while(len_trim(tags(i))==0 .and. i > 0)
    i = i - 1
 enddo
 num_in_header = i

end function num_in_header

!----------------------------------------
! open a dump file and write the file id
! and other generic header information
!----------------------------------------
subroutine open_dumpfile_w(iunit,filename,fileid,ierr,singleprec)
 integer,              intent(in)  :: iunit
 character(len=*),     intent(in)  :: filename
 character(len=lenid), intent(in)  :: fileid
 integer,              intent(out) :: ierr
 logical,              intent(in), optional :: singleprec
 integer :: i1
 logical :: r4
 real    :: r1

 open(unit=iunit,file=filename,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) return

 i1 = int1
 r1 = real(int2)

 r4 = .false.
 if (present(singleprec)) r4 = singleprec

 if (r4) then
    write(iunit,iostat=ierr) int1,real(r1,kind=4),int2,iversion,int1o
 else
    write(iunit,iostat=ierr) int1,r1,int2,iversion,int1o
 endif
 if (ierr /= 0) return

 write(iunit,iostat=ierr) fileid

end subroutine open_dumpfile_w

!-----------------------------------------
! open a dump file and read the file id
! and generic header information
!-----------------------------------------
subroutine open_dumpfile_r(iunit,filename,fileid,ierr,singleprec,requiretags,tagged)
 integer,              intent(in)  :: iunit
 character(len=*),     intent(in)  :: filename
 character(len=lenid), intent(out) :: fileid
 integer,              intent(out) :: ierr
 logical,              intent(in),  optional :: singleprec,requiretags
 logical,              intent(out), optional :: tagged
 integer(kind=4) :: int1i,int2i,int3i
 integer         :: iversion_file,ierr1
 logical         :: r4,must_have_tags
 real(kind=4)    :: r1s
 real            :: r1i

 r4 = .false.
 must_have_tags = .false.
 if (present(singleprec)) r4 = singleprec
 if (present(requiretags)) must_have_tags = requiretags
!
!--open dump file
!
 open(unit=iunit,file=filename,status='old',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    ierr = ierr_fileopen
    return
 endif
!
!--read output file
!
 if (r4) then
    read (iunit,iostat=ierr1) int1i,r1s,int2i,iversion_file,int3i
 else
    read (iunit,iostat=ierr1) int1i,r1i,int2i,iversion_file,int3i
 endif
 if (int1i /= int1 .and. int1i /= int1o) then
    ierr = ierr_endian
    return
 endif

! handle version numbers
! version 0 had iversion=690706

 if (iversion_file==int1o) iversion_file = 0
 if (iversion_file > iversion) then
    !write (*,"(a,i2,a,i2)") 'error 4 in readdump: format version is ',iversion, &
    !   ' but this version of Phantom can only read ',maxversion
    ierr = ierr_version
 endif

 read (iunit,iostat=ierr1) fileid

 if (int2i /= int2 .and. int2i /= int2o) then
    ierr = ierr_realsize
    return
 endif
 if (int3i /= int1o) then
    ierr = ierr_intsize
    return
 endif

 ! generic read error, only return this if other errors have been ruled out first
 if (ierr1 /= 0) then
    ierr = ierr_unknown
    return
 endif

 ! return error if not tagged format, if this is required
 if (must_have_tags) then
    if (fileid(2:2) /= 'T' .and. fileid(2:2) /= 't') then
       ierr = ierr_notags
       return
    endif
 endif

 ! return whether or not file is in tagged format
 if (present(tagged)) then
    tagged = (fileid(2:2) == 'T' .or. fileid(2:2) == 't')
 endif

end subroutine open_dumpfile_r

!-------------------------------------------------------
!+
!  error handling routine so that errors can be handled
!  outside of this library
!+
!-------------------------------------------------------
character(len=60) function get_error_text(ierr)
 integer, intent(in) :: ierr

 select case(ierr)
 case(ierr_fileopen)
    get_error_text = 'error opening file'
 case(ierr_endian)
    get_error_text = 'wrong endian?'
 case(ierr_version)
    get_error_text = 'file format version newer than current code can read'
 case(ierr_realsize)
    get_error_text = 'default real size wrong'
 case(ierr_intsize)
    get_error_text = 'default int size wrong'
 case(ierr_notags)
    get_error_text = 'routine requires tagged format but not detected'
 case(ierr_arraysize)
    get_error_text = 'array size too small for requested data'
 case default
    get_error_text = 'unknown error'
 end select

end function get_error_text

!-------------------------------------------------------
!+
!  read the file header into the dump_header structure
!+
!-------------------------------------------------------
subroutine read_header(iunit,hdr,ierr,singleprec,tagged)
 integer,      intent(in) :: iunit
 type(dump_h), intent(out) :: hdr
 integer,      intent(out) :: ierr
 logical,      intent(in), optional :: singleprec
 logical,      intent(in), optional :: tagged
 logical :: convert_prec,tags
 integer :: i,n
 real(kind=4), allocatable :: dumr4(:)

 convert_prec = .false.
 if (present(singleprec)) convert_prec = singleprec

 tags = .true.
 if (present(tagged)) tags = tagged

 do i=1,ndatatypes
    read (iunit,iostat=ierr) n
    if (n < 0) n = 0
    hdr%nums(i) = n
    select case(i)
    case(i_int)
       allocate(hdr%inttags(n),hdr%intvals(n),stat=ierr)
       if (n > 0) then
          hdr%inttags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%inttags
          read(iunit,iostat=ierr) hdr%intvals
       endif
    case(i_int1)
       allocate(hdr%int1tags(n),hdr%int1vals(n),stat=ierr)
       if (n > 0) then
          hdr%int1tags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%int1tags
          read(iunit,iostat=ierr) hdr%int1vals
       endif
    case(i_int2)
       allocate(hdr%int2tags(n),hdr%int2vals(n),stat=ierr)
       if (n > 0) then
          hdr%int2tags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%int2tags
          read(iunit,iostat=ierr) hdr%int2vals
       endif
    case(i_int4)
       allocate(hdr%int4tags(n),hdr%int4vals(n),stat=ierr)
       if (n > 0) then
          hdr%int4tags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%int4tags
          read(iunit,iostat=ierr) hdr%int4vals
       endif
    case(i_int8)
       allocate(hdr%int8tags(n),hdr%int8vals(n),stat=ierr)
       if (n > 0) then
          hdr%int8tags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%int8tags
          read(iunit,iostat=ierr) hdr%int8vals
       endif
    case(i_real)
       allocate(hdr%realtags(n),hdr%realvals(n),stat=ierr)
       if (n > 0) then
          hdr%realtags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%realtags
          if (convert_prec .and. kind(0.) /= 4) then
             allocate(dumr4(n),stat=ierr)
             read(iunit,iostat=ierr) dumr4
             hdr%realvals(1:n) = real(dumr4(1:n))
             deallocate(dumr4)
          else
             read(iunit,iostat=ierr) hdr%realvals
          endif
       endif
    case(i_real4)
       allocate(hdr%real4tags(n),hdr%real4vals(n),stat=ierr)
       if (n > 0) then
          hdr%real4tags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%real4tags
          read(iunit,iostat=ierr) hdr%real4vals
       endif
    case(i_real8)
       allocate(hdr%real8tags(n),hdr%real8vals(n),stat=ierr)
       if (n > 0) then
          hdr%real8tags(:) = ''
          if (tags) read(iunit,iostat=ierr) hdr%real8tags
          read(iunit,iostat=ierr) hdr%real8vals
       endif
    end select
 enddo

end subroutine read_header

!-------------------------------------------------------
!+
!  allocate the dump header structure for writing
!  IN: (all are optional, default size is 256)
!     nint,nint2 : size of buffer for each data type
!               e.g. (/256,0,0,256,256,256,256,256/)
!  OUT:
!     hdr : header structure allocated for each
!           data type
!+
!-------------------------------------------------------
function allocate_header(nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8,err) result(hdr)
 integer, intent(in),  optional :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
 integer, intent(out), optional :: err
 type(dump_h) :: hdr
 integer      :: size_type(ndatatypes)
 integer      :: ierrs(ndatatypes)
 integer      :: ierr

 ! make sure header is deallocated first
 call free_header(hdr,ierr)

 size_type(:) = maxphead
 if (present(nint))   size_type(i_int)  = nint
 if (present(nint1))  size_type(i_int1) = nint1
 if (present(nint2))  size_type(i_int2) = nint2
 if (present(nint4))  size_type(i_int4) = nint4
 if (present(nint8))  size_type(i_int8) = nint8
 if (present(nreal))  size_type(i_real) = nreal
 if (present(nreal4)) size_type(i_real4) = nreal4
 if (present(nreal8)) size_type(i_real8) = nreal8

 if (present(err)) err = 0
 ierrs(:) = 0
 hdr%nums(:) = 0
 if (size_type(i_int) > 0)  then
    allocate(hdr%inttags(size_type(i_int)),hdr%intvals(size_type(i_int)),stat=ierrs(1))
    if (ierrs(1)==0) hdr%inttags(:) = ''
 endif
 if (size_type(i_int1) > 0) then
    allocate(hdr%int1tags(size_type(i_int1)),hdr%int1vals(size_type(i_int1)),stat=ierrs(2))
    if (ierrs(2)==0) hdr%int1tags(:) = ''
 endif
 if (size_type(i_int2) > 0) then
    allocate(hdr%int2tags(size_type(i_int2)),hdr%int2vals(size_type(i_int2)),stat=ierrs(3))
    if (ierrs(3)==0) hdr%int2tags(:) = ''
 endif
 if (size_type(i_int4) > 0) then
    allocate(hdr%int4tags(size_type(i_int4)),hdr%int4vals(size_type(i_int4)),stat=ierrs(4))
    if (ierrs(4)==0) hdr%int4tags(:) = ''
 endif
 if (size_type(i_int8) > 0) then
    allocate(hdr%int8tags(size_type(i_int8)),hdr%int8vals(size_type(i_int8)),stat=ierrs(5))
    if (ierrs(5)==0) hdr%int8tags(:) = ''
 endif
 if (size_type(i_real) > 0)  then
    allocate(hdr%realtags(size_type(i_real)),hdr%realvals(size_type(i_real)),stat=ierrs(6))
    if (ierrs(6)==0) hdr%realtags(:) = ''
 endif
 if (size_type(i_real4) > 0)  then
    allocate(hdr%real4tags(size_type(i_real4)),hdr%real4vals(size_type(i_real4)),stat=ierrs(7))
    if (ierrs(7)==0) hdr%real4tags(:) = ''
 endif
 if (size_type(i_real8) > 0)  then
    allocate(hdr%real8tags(size_type(i_real8)),hdr%real8vals(size_type(i_real8)),stat=ierrs(8))
    if (ierrs(8)==0) hdr%real8tags(:) = ''
 endif

 if (present(err) .and. any(ierrs /= 0)) err = 1

end function allocate_header

!-------------------------------------------------------
!+
!  allocate the dump header structure for writing
!+
!-------------------------------------------------------
subroutine free_header(hdr,ierr)
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out), optional :: ierr

 if (present(ierr)) ierr = 0
 if (allocated(hdr%inttags))   deallocate(hdr%inttags)
 if (allocated(hdr%int1tags))  deallocate(hdr%int1tags)
 if (allocated(hdr%int2tags))  deallocate(hdr%int2tags)
 if (allocated(hdr%int4tags))  deallocate(hdr%int4tags)
 if (allocated(hdr%int8tags))  deallocate(hdr%int8tags)
 if (allocated(hdr%realtags))  deallocate(hdr%realtags)
 if (allocated(hdr%real4tags)) deallocate(hdr%real4tags)
 if (allocated(hdr%real8tags)) deallocate(hdr%real8tags)

 if (allocated(hdr%intvals))   deallocate(hdr%intvals)
 if (allocated(hdr%int1vals))  deallocate(hdr%int1vals)
 if (allocated(hdr%int2vals))  deallocate(hdr%int2vals)
 if (allocated(hdr%int4vals))  deallocate(hdr%int4vals)
 if (allocated(hdr%realvals))  deallocate(hdr%realvals)
 if (allocated(hdr%real4vals)) deallocate(hdr%real4vals)
 if (allocated(hdr%real8vals)) deallocate(hdr%real8vals)

end subroutine free_header

!-------------------------------------------------------
!+
!  print contents of header structure
!+
!-------------------------------------------------------
subroutine print_header(hdr)
 type(dump_h), intent(in) :: hdr
 integer :: i

 if (allocated(hdr%inttags) .and. allocated(hdr%intvals)) then
    do i=1,size(hdr%inttags)
       print*,hdr%inttags(i),hdr%intvals(i)
    enddo
 endif
 if (allocated(hdr%int1tags) .and. allocated(hdr%int1vals)) then
    do i=1,size(hdr%int1tags)
       print*,hdr%int1tags(i),hdr%int1vals(i)
    enddo
 endif
 if (allocated(hdr%int2tags) .and. allocated(hdr%int2vals)) then
    do i=1,size(hdr%int2tags)
       print*,hdr%int2tags(i),hdr%int2vals(i)
    enddo
 endif
 if (allocated(hdr%int4tags) .and. allocated(hdr%int4vals)) then
    do i=1,size(hdr%int4tags)
       print*,hdr%int4tags(i),hdr%int4vals(i)
    enddo
 endif
 if (allocated(hdr%realtags) .and. allocated(hdr%realvals)) then
    do i=1,size(hdr%realtags)
       print*,hdr%realtags(i),hdr%realvals(i)
    enddo
 endif
 if (allocated(hdr%real4tags) .and. allocated(hdr%real4vals)) then
    do i=1,size(hdr%real4tags)
       print*,hdr%real4tags(i),hdr%real4vals(i)
    enddo
 endif
 if (allocated(hdr%real8tags) .and. allocated(hdr%real8vals)) then
    do i=1,size(hdr%real8tags)
       print*,hdr%real8tags(i),hdr%real8vals(i)
    enddo
 endif

end subroutine print_header

!-------------------------------------------------------
!+
!  write the header to file
!+
!-------------------------------------------------------
subroutine write_header(iunit,hdr,ierr,singleprec)
 integer,      intent(in)  :: iunit
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 logical,      intent(in), optional :: singleprec
 integer :: number,i,j,ierrs(17)
 integer,         parameter :: idum  = 0
 integer(kind=1), parameter :: idum1 = 0_1
 integer(kind=2), parameter :: idum2 = 0_2
 integer(kind=4), parameter :: idum4 = 0_4
 integer(kind=8), parameter :: idum8 = 0_8
 real, parameter :: dum = 0.
 real(kind=4), parameter :: dum4 = 0._4
 real(kind=8), parameter :: dum8 = 0._8
 logical :: sing_prec

 ! optional argument to write real header in single precision
 sing_prec = .false.
 if (present(singleprec)) sing_prec = singleprec

 ierrs = 0
 do i=1,ndatatypes
    number = hdr%nums(i)
    write(iunit,iostat=ierrs(1)) number
    if (number > 0) then
       select case(i)
       case(i_int)
          if (allocated(hdr%inttags) .and. allocated(hdr%intvals)) then
             write(iunit,iostat=ierrs(2)) (hdr%inttags(j),j=1,number)
             write(iunit,iostat=ierrs(3)) hdr%intvals(1:number)
          else
             write(iunit,iostat=ierrs(2)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(3)) (idum,j=1,number)
          endif
       case(i_int1)
          if (allocated(hdr%int1tags) .and. allocated(hdr%int1vals)) then
             write(iunit,iostat=ierrs(4)) (hdr%int1tags(j),j=1,number)
             write(iunit,iostat=ierrs(5)) hdr%int1vals(1:number)
          else
             write(iunit,iostat=ierrs(4)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(5)) (idum1,j=1,number)
          endif
       case(i_int2)
          if (allocated(hdr%int2tags) .and. allocated(hdr%int2vals)) then
             write(iunit,iostat=ierrs(6)) (hdr%int2tags(j),j=1,number)
             write(iunit,iostat=ierrs(7)) hdr%int2vals(1:number)
          else
             write(iunit,iostat=ierrs(6)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(7)) (idum2,j=1,number)
          endif
       case(i_int4)
          if (allocated(hdr%int4tags) .and. allocated(hdr%int4vals)) then
             write(iunit,iostat=ierrs(8)) (hdr%int4tags(j),j=1,number)
             write(iunit,iostat=ierrs(9)) hdr%int4vals(1:number)
          else
             write(iunit,iostat=ierrs(8)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(9)) (idum4,j=1,number)
          endif
       case(i_int8)
          if (allocated(hdr%int8tags) .and. allocated(hdr%int8vals)) then
             write(iunit,iostat=ierrs(10)) (hdr%int8tags(j),j=1,number)
             write(iunit,iostat=ierrs(11)) hdr%int8vals(1:number)
          else
             write(iunit,iostat=ierrs(10)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(11)) (idum8,j=1,number)
          endif
       case(i_real)
          if (allocated(hdr%realtags) .and. allocated(hdr%realvals)) then
             write(iunit,iostat=ierrs(12)) (hdr%realtags(j),j=1,number)
             if (sing_prec) then
                write(iunit,iostat=ierrs(13)) real(hdr%realvals(1:number),kind=4)
             else
                write(iunit,iostat=ierrs(13)) hdr%realvals(1:number)
             endif
          else
             write(iunit,iostat=ierrs(12)) (tag('unknown'),j=1,number)
             if (sing_prec) then
                write(iunit,iostat=ierrs(13)) (real(dum,kind=4),j=1,number)
             else
                write(iunit,iostat=ierrs(13)) (dum,j=1,number)
             endif
          endif
       case(i_real4)
          if (allocated(hdr%real4tags) .and. allocated(hdr%real4vals)) then
             write(iunit,iostat=ierrs(14)) (hdr%real4tags(j),j=1,number)
             write(iunit,iostat=ierrs(15)) hdr%real4vals(1:number)
          else
             write(iunit,iostat=ierrs(14)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(15)) (dum4,j=1,number)
          endif
       case(i_real8)
          if (allocated(hdr%real8tags) .and. allocated(hdr%real8vals)) then
             write(iunit,iostat=ierrs(16)) (hdr%real8tags(j),j=1,number)
             write(iunit,iostat=ierrs(17)) hdr%real8vals(1:number)
          else
             write(iunit,iostat=ierrs(16)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(17)) (dum8,j=1,number)
          endif
       end select
    endif
 enddo
 if (any(ierrs /= 0)) ierr = 1

end subroutine write_header

!---------------------------------------------------------------------
!+
!  Write int*1 array to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_int1(ib,iarr,my_tag,len,ikind,ipass,iunit,nums,nerr,func)
 integer(kind=1),  intent(in) :: iarr(:)
 character(len=*), intent(in) :: my_tag
 integer, intent(in)    :: ib,len,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 !procedure(integer(kind=1)), pointer, optional :: func
 interface
  integer(kind=1) pure function func(x)
   integer(kind=1), intent(in) :: x
  end function func
 end interface
 optional :: func
 !integer(kind=1), optional :: func
 integer :: i,ierr

 ierr = 0
 ! check if kind matches
 if (ikind==i_int1) then
    if (ipass==1) then
       nums(i_int1,ib) = nums(i_int1,ib) + 1
    elseif (ipass==2) then
       write(iunit,iostat=ierr) tag(my_tag)
       if (present(func)) then
          write(iunit,iostat=ierr) (func(iarr(i)),i=1,len)
       else
          write(iunit,iostat=ierr) iarr(1:len)
       endif
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_int1

!---------------------------------------------------------------------
!+
!  Write int*4 array to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_int4(ib,iarr,my_tag,len,ikind,ipass,iunit,nums,nerr,func)
 integer(kind=4),  intent(in) :: iarr(:)
 character(len=*), intent(in) :: my_tag
 integer, intent(in)    :: ib,len,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 !procedure(integer(kind=1)), pointer, optional :: func
 interface
  integer(kind=4) pure function func(x)
   integer(kind=4), intent(in) :: x
  end function func
 end interface
 optional :: func
 !integer(kind=1), optional :: func
 integer :: i,ierr

 ierr = 0
 ! check if kind matches
 if (ikind==i_int4) then
    if (ipass==1) then
       nums(i_int4,ib) = nums(i_int4,ib) + 1
    elseif (ipass==2) then
       write(iunit,iostat=ierr) tag(my_tag)
       if (present(func)) then
          write(iunit,iostat=ierr) (func(iarr(i)),i=1,len)
       else
          write(iunit,iostat=ierr) iarr(1:len)
       endif
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_int4

!---------------------------------------------------------------------
!+
!  Write int*4 array to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_int8(ib,iarr,my_tag,len,ikind,ipass,iunit,nums,nerr,func)
 integer(kind=8),  intent(in) :: iarr(:)
 character(len=*), intent(in) :: my_tag
 integer, intent(in)    :: ib,len,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 !procedure(integer(kind=1)), pointer, optional :: func
 interface
  integer(kind=8) pure function func(x)
   integer(kind=8), intent(in) :: x
  end function func
 end interface
 optional :: func
 !integer(kind=1), optional :: func
 integer :: i,ierr

 ierr = 0
 ! check if kind matches
 if (ikind==i_int8) then
    if (ipass==1) then
       nums(i_int8,ib) = nums(i_int8,ib) + 1
    elseif (ipass==2) then
       write(iunit,iostat=ierr) tag(my_tag)
       if (present(func)) then
          write(iunit,iostat=ierr) (func(iarr(i)),i=1,len)
       else
          write(iunit,iostat=ierr) iarr(1:len)
       endif
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_int8

!---------------------------------------------------------------------
!+
!  Write multidimensional real*4 array arr(len1,len2)
!  to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_int4arr(ib,iarr,my_tag,len1,len2,ikind,ipass,iunit,nums,nerr,index)
 integer(kind=4),  intent(in) :: iarr(:,:)
 character(len=*), intent(in) :: my_tag(:)
 integer, intent(in)    :: ib,len1,len2,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 integer, intent(in), optional :: index
 integer :: j,i,istart,iend,ierr

 ierr = 0
 if (present(index)) then
    istart = index
    iend   = index
 else
    istart = 1
    iend   = len1
 endif
 ! check if kind matches
 if (ikind==i_int4) then
    !print*,ipass,' WRITING ',my_tag(istart:iend),' as ',i_int8
    if (ipass==1) then
       nums(i_int4,ib) = nums(i_int4,ib) + (iend - istart) + 1
    elseif (ipass==2) then
       do j=istart,iend
          write(iunit,iostat=ierr) tag(my_tag(j))
          write(iunit,iostat=ierr) (iarr(j,i),i=1,len2)
       enddo
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_int4arr

!---------------------------------------------------------------------
!+
!  Write real*4 array to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_real4(ib,arr,my_tag,len,ikind,ipass,iunit,nums,nerr,func,use_kind,singleprec)
 real(kind=4),     intent(in) :: arr(:)
 character(len=*), intent(in) :: my_tag
 integer, intent(in)    :: ib,len,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 interface
  real(kind=4) pure function func(x)
   real(kind=4), intent(in) :: x
  end function func
 end interface
 optional :: func
 !real(kind=4), optional :: func
 integer, intent(in), optional :: use_kind
 logical, intent(in), optional :: singleprec
 integer :: i,imatch,ierr

 ierr = 0
 ! use default real if it matches, unless kind is specified
 if (kind(0.)==4 .and. .not.present(use_kind)) then
    imatch = i_real
 else
    imatch = i_real4
 endif
 ! check if kind matches
 if (ikind==imatch) then
    if (ipass==1) then
       nums(imatch,ib) = nums(imatch,ib) + 1
    elseif (ipass==2) then
       write(iunit,iostat=ierr) tag(my_tag)
       if (present(func)) then
          write(iunit,iostat=ierr) (func(arr(i)),i=1,len)
       else
          write(iunit,iostat=ierr) arr(1:len)
       endif
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_real4

!---------------------------------------------------------------------
!+
!  Write real*4 array to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_real8(ib,arr,my_tag,len,ikind,ipass,iunit,nums,nerr,func,use_kind,singleprec)
 real(kind=8),     intent(in) :: arr(:)
 character(len=*), intent(in) :: my_tag
 integer, intent(in)    :: ib,len,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 interface
  real(kind=8) pure function func(x)
   real(kind=8), intent(in) :: x
  end function func
 end interface
 optional :: func
 !real(kind=8), optional :: func
 integer, intent(in), optional :: use_kind
 logical, intent(in), optional :: singleprec
 integer :: i,imatch,ierr
 logical :: use_singleprec

 ierr = 0
 use_singleprec = .false.
 if (present(singleprec)) use_singleprec = singleprec
 ! use default real if it matches, unless kind is specified
 if ((kind(0.)==8 .or. use_singleprec).and.(.not.present(use_kind))) then
    imatch = i_real
 elseif (present(use_kind)) then
    if (use_kind==4) then
       imatch = i_real4
    else
       imatch = i_real8
    endif
 else
    imatch = i_real8
 endif
 ! check if kind matches
 if (ikind==imatch) then
    !print*,ipass,' WRITING ',my_tag,' as ',imatch,use_singleprec
    if (ipass==1) then
       nums(imatch,ib) = nums(imatch,ib) + 1
    elseif (ipass==2) then
       write(iunit,iostat=ierr) tag(my_tag)
       if (present(func)) then
          write(iunit,iostat=ierr) (func(arr(i)),i=1,len)
       else
          if (imatch==i_real4 .or. use_singleprec) then
             write(iunit,iostat=ierr) (real(arr(i),kind=4),i=1,len)
          else
             write(iunit,iostat=ierr) arr(1:len)
          endif
       endif
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_real8

!---------------------------------------------------------------------
!+
!  Write multidimensional real*4 array arr(len1,len2)
!  to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_real4arr(ib,arr,my_tag,len1,len2,ikind,ipass,iunit,nums,nerr,use_kind,index,singleprec)
 real(kind=4),     intent(in) :: arr(:,:)
 character(len=*), intent(in) :: my_tag(:)
 integer, intent(in)    :: ib,len1,len2,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 integer, intent(in), optional :: use_kind,index
 logical, intent(in), optional :: singleprec
 integer :: j,i,imatch,istart,iend,ierr

 ierr = 0
 ! use default real if it matches, unless kind is specified
 if (kind(0.)==4 .and. .not.present(use_kind)) then
    imatch = i_real
 else
    imatch = i_real4
 endif
 if (present(index)) then
    istart = index
    iend   = index
 else
    istart = 1
    iend   = len1
 endif
 ! check if kind matches
 if (ikind==imatch) then
    !print*,ipass,' WRITING ',my_tag(istart:iend),' as ',imatch
    if (ipass==1) then
       nums(imatch,ib) = nums(imatch,ib) + (iend - istart) + 1
    elseif (ipass==2) then
       do j=istart,iend
          write(iunit,iostat=ierr) tag(my_tag(j))
          write(iunit,iostat=ierr) (arr(j,i),i=1,len2)
       enddo
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_real4arr

!---------------------------------------------------------------------
!+
!  Write multidimensional real*8 array arr(len1,len2)
!  to block header (ipass=1) or to file (ipass=2)
!+
!---------------------------------------------------------------------
subroutine write_array_real8arr(ib,arr,my_tag,len1,len2,ikind,ipass,iunit,nums,nerr,use_kind,index,singleprec)
 real(kind=8),     intent(in) :: arr(:,:)
 character(len=*), intent(in) :: my_tag(:)
 integer, intent(in)    :: ib,len1,len2,ikind,ipass,iunit
 integer, intent(inout) :: nums(:,:)
 integer, intent(inout) :: nerr
 integer, intent(in), optional :: use_kind,index
 logical, intent(in), optional :: singleprec
 integer :: j,i,imatch,istart,iend,ierr
 logical :: use_singleprec

 ierr = 0
 use_singleprec = .false.
 if (present(singleprec)) use_singleprec = singleprec
 ! use default real if it matches, unless kind is specified
 if ((kind(0.)==8 .or. use_singleprec) .and. .not.present(use_kind)) then
    imatch = i_real
 elseif (present(use_kind)) then
    if (use_kind==4) then
       imatch = i_real4
    else
       imatch = i_real8
    endif
 else
    imatch = i_real8
 endif
 if (present(index)) then
    istart = index
    iend   = index
 else
    istart = 1
    iend   = len1
 endif
 ! check if kind matches
 if (ikind==imatch) then
    !print*,ipass,' WRITING ',my_tag(istart:iend),' as ',imatch,(iend-istart)+1,istart,iend,use_singleprec
    if (ipass==1) then
       nums(imatch,ib) = nums(imatch,ib) + (iend - istart) + 1
    elseif (ipass==2) then
       do j=istart,iend
          write(iunit,iostat=ierr) tag(my_tag(j))
          if (imatch==i_real4 .or. use_singleprec) then
             !print*, "done ", my_tag(j), " | ", tag(my_tag(j))
             write(iunit,iostat=ierr) (real(arr(j,i),kind=4),i=1,len2)
          else
             write(iunit,iostat=ierr) (arr(j,i),i=1,len2)
          endif
       enddo
    endif
 endif
 if (ierr /= 0) nerr = nerr + 1

end subroutine write_array_real8arr

!---------------------------------------------------
!+
!  Write the actual header for each block to disk
!+
!---------------------------------------------------
subroutine write_block_header(nblocks,number,nums,iunit,ierr)
 integer,         intent(in)  :: nblocks
 integer(kind=8), intent(in)  :: number(nblocks)
 integer,         intent(in)  :: nums(ndatatypes,nblocks),iunit
 integer,         intent(out) :: ierr
 integer :: iblock

 do iblock=1,nblocks
    write(iunit,iostat=ierr) number(iblock), nums(1:ndatatypes,iblock)
 enddo

end subroutine write_block_header

!---------------------------------------------------
!+
!  read the header for each block from disk
!+
!---------------------------------------------------
subroutine read_block_header(nblocks,number,nums,iunit,ierr)
 integer,         intent(in)  :: nblocks
 integer(kind=8), intent(out) :: number(nblocks)
 integer,         intent(out) :: nums(ndatatypes,nblocks)
 integer,         intent(in)  :: iunit
 integer,         intent(out) :: ierr
 integer :: iblock

 number(:) = 0
 nums(:,:) = 0
 do iblock=1,nblocks
    read(iunit,iostat=ierr) number(iblock), nums(1:ndatatypes,iblock)
 enddo

end subroutine read_block_header

!--------------------------------------------------------------------
!+
!  utility to determine whether to read a particular block
!  in the dump file, in whole or in part.
!  Allows limited changes to number of threads.
!+
!--------------------------------------------------------------------
subroutine get_blocklimits(npartblock,nblocks,nthreads,id,iblock,noffset,npartread,ierr)
 integer(kind=8), intent(in)  :: npartblock
 integer,         intent(in)  :: nblocks,nthreads,id,iblock
 integer,         intent(out) :: noffset,npartread,ierr
 integer                      :: nblocksperthread,nthreadsperblock
 character(len=15), parameter :: tag = 'get_blocklimits'
!
!--check for errors in input
!
 ierr = 0
 if (npartblock < 0) then
    write(*,*) 'get_blocklimits: block in dump file has npartinblock < 0'
    ierr = 1
 elseif (npartblock > huge(npartread)) then
    write(*,*) 'get_blocklimits: number of particles in block exceeds 32 bit limit'
    ierr = 2
 endif
 if (ierr /= 0) return
!
!--usual situation: nblocks = nprocessors
!  read whole block if id = iblock
!
 if (nblocks==nthreads) then
    if (id==iblock-1) then
       !--read whole block
       npartread = int(npartblock)
       noffset   = 0
    else
       !--do not read block
       npartread = 0
       noffset   = 0
    endif

 elseif (nblocks > nthreads .and. mod(nblocks,nthreads)==0) then
!
!--if more blocks than processes and nblocks exactly divisible by nthreads,
!  then just read more than one block per thread
!
    nblocksperthread = nblocks/nthreads
    if (id==(iblock-1)/nblocksperthread) then
       npartread = int(npartblock)
       noffset   = 0
    else
       npartread = 0
       noffset   = 0
    endif

 elseif (nthreads > nblocks .and. mod(nthreads,nblocks)==0) then
!
!--if more threads than blocks, and exactly divisible, read fractions of blocks only
!
    nthreadsperblock = nthreads/nblocks
    if (id/nthreadsperblock==iblock-1) then
       npartread = int((npartblock-1)/nthreadsperblock) + 1
       noffset   = mod(id,nthreadsperblock)*npartread

       if (mod(id,nthreadsperblock)==nthreadsperblock-1) then
          !--last thread has remainder for non-exactly divisible numbers of particles
          npartread = int(npartblock) - (nthreadsperblock-1)*npartread
          !--die if we would need to load balance between more than the last processor.
          if (npartread < 0) then
             print*,' npart to read from last block =',npartread
             print*,trim(tag)//' error assigning npart to last thread'
             ierr = 3
             return
          endif
       endif
    else
       npartread = 0
       noffset   = 0
    endif
 else
    noffset = 0
    npartread = 0
    ierr = 4
    print*,' ERROR: rearrangement of ',nblocks,' blocks to ',nthreads,' threads not implemented'
    return
 endif

end subroutine get_blocklimits

!--------------------------------------------------------------------
!+
!  Routine for extracting int*1 array from main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_int1(iarr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 integer(kind=1),  intent(inout) :: iarr(:)
 character(len=*), intent(in)    :: arr_tag,tag
 logical,          intent(inout) :: got_arr
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer         :: i
 integer(kind=1) :: dum
 logical         :: match_datatype

 if (matched) return
 match_datatype = (ikind==i_int1)

 if (match_tag(tag,arr_tag) .and. .not.matched) then
    matched    = .true.
    if (match_datatype) then
       got_arr = .true.
       read(iunit,iostat=ierr) (dum,i=1,noffset),iarr(i1:i2)
    else
       print*,'ERROR: wrong datatype for '//trim(tag)//' (is not int1)'
       read(iunit,iostat=ierr)
    endif
 endif

end subroutine read_array_int1

!--------------------------------------------------------------------
!+
!  Routine for extracting int*4 array from main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_int4(iarr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 integer(kind=4),  intent(inout) :: iarr(:)
 character(len=*), intent(in)    :: arr_tag,tag
 logical,          intent(inout) :: got_arr
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer         :: i
 integer(kind=4) :: dum
 logical         :: match_datatype

 if (matched) return
 match_datatype = (ikind==i_int4)

 if (match_tag(tag,arr_tag) .and. .not.matched) then
    matched    = .true.
    if (match_datatype) then
       got_arr = .true.
       read(iunit,iostat=ierr) (dum,i=1,noffset),iarr(i1:i2)
    else
       print*,'ERROR: wrong datatype for '//trim(tag)//' (is not int4)'
       read(iunit,iostat=ierr)
    endif
 endif

end subroutine read_array_int4

!--------------------------------------------------------------------
!+
!  Routine for extracting int*8 array from main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_int8(iarr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 integer(kind=8),  intent(inout) :: iarr(:)
 character(len=*), intent(in)    :: arr_tag,tag
 logical,          intent(inout) :: got_arr
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer         :: i
 integer(kind=8) :: dum
 logical         :: match_datatype

 if (matched) return
 match_datatype = (ikind==i_int8)

 if (match_tag(tag,arr_tag) .and. .not.matched) then
    matched    = .true.
    if (match_datatype) then
       got_arr = .true.
       read(iunit,iostat=ierr) (dum,i=1,noffset),iarr(i1:i2)
    else
       print*,'ERROR: wrong datatype for '//trim(tag)//' (is not int8)'
       read(iunit,iostat=ierr)
    endif
 endif

end subroutine read_array_int8

!--------------------------------------------------------------------
!+
!  Routine for extracting multi-d int*8 array from
!  main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_int4arr(iarr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 integer(kind=4),  intent(inout) :: iarr(:,:)
 character(len=*), intent(in)    :: arr_tag(:),tag
 logical,          intent(inout) :: got_arr(:)
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer         :: i,j,nread
 integer(kind=4) :: dum
 logical         :: match_datatype
 integer(kind=4), allocatable :: dummyi4(:)

 if (matched) return
 match_datatype = (ikind==i_int8)

 do j=1,min(size(iarr(:,1)),size(arr_tag))
    if (match_tag(tag,arr_tag(j)) .and. .not.matched) then
       matched    = .true.
       if (match_datatype) then
          got_arr(j) = .true.
          nread = i2-i1+1
          allocate(dummyi4(nread)) ! to avoid seg fault with ifort
          read(iunit,iostat=ierr) (dum,i=1,noffset),dummyi4(1:nread)
          iarr(j,i1:i2) = dummyi4(:)
          deallocate(dummyi4)
       else
          print*,'ERROR: wrong datatype for '//trim(tag)//' (is not int8)'
          read(iunit,iostat=ierr)
       endif
    endif
 enddo

end subroutine read_array_int4arr


!--------------------------------------------------------------------
!+
!  Routine for extracting real*4 array from main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_real4(arr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 real(kind=4),     intent(inout) :: arr(:)
 character(len=*), intent(in)    :: arr_tag,tag
 logical,          intent(inout) :: got_arr
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer      :: i
 real(kind=4) :: dum
 logical      :: match_datatype

 if (matched .or. ikind < i_real) return
 match_datatype = (ikind==i_real4 .or. (kind(0.)==4 .and. ikind==i_real))

 if (match_tag(tag,arr_tag) .and. .not.matched) then
    matched    = .true.
    if (match_datatype) then
       got_arr = .true.
       if (i2 > size(arr)) then
          print*,'ERROR: array size too small reading array: need ',i2,' got ',size(arr)
          read(iunit,iostat=ierr)
       endif
       read(iunit,iostat=ierr) (dum,i=1,noffset),arr(i1:i2)
    else
       print*,'ERROR: wrong datatype for '//trim(tag)//' (is not real4)'
       read(iunit,iostat=ierr)
    endif
 endif

end subroutine read_array_real4

!--------------------------------------------------------------------
!+
!  Routine for extracting multi-d real*4 array from
!  main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_real4arr(arr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 real(kind=4),     intent(inout) :: arr(:,:)
 character(len=*), intent(in)    :: arr_tag(:),tag
 logical,          intent(inout) :: got_arr(:)
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer      :: i,j,nread
 real(kind=4) :: dum
 real(kind=8) :: dumr8
 real(kind=4), allocatable :: dummy(:)
 real(kind=8), allocatable :: dummyr8(:)
 logical      :: match_datatype

 if (matched .or. ikind < i_real) return
 match_datatype = (ikind==i_real4 .or. (kind(0.)==4 .and. ikind==i_real))

 do j=1,min(size(arr,dim=1),size(arr_tag))
    if (match_tag(tag,arr_tag(j)) .and. .not.matched) then
       matched    = .true.
       if (match_datatype) then
          got_arr(j) = .true.
          if (i2 > size(arr)) then
             print*,'ERROR: array size too small reading '//trim(tag)
             ierr = ierr_arraysize
             return
          endif
          nread = i2-i1+1
          allocate(dummy(nread))
          read(iunit,iostat=ierr) (dum,i=1,noffset),dummy(1:nread)
          arr(j,i1:i2) = dummy
          deallocate(dummy)
          !read(iunit,iostat=ierr) (dum,i=1,noffset),arr(j,i1:i2)
       elseif (ikind==i_real4) then
          got_arr(j) = .true.
          !print*,'WARNING: converting '//trim(tag)//' from real*8->real*4'
          nread = i2-i1+1
          allocate(dummyr8(nread))
          read(iunit,iostat=ierr) (dumr8,i=1,noffset),dummyr8(1:nread)
          if (i2 > size(arr)) then
             print*,'ERROR: array size too small reading '//trim(tag)
             ierr = ierr_arraysize
             return
          endif
          arr(j,i1:i2) = real(dummyr8(:),kind=4)
          deallocate(dummyr8)
       else
          print*,'ERROR: wrong datatype for '//trim(tag)//' (is not real4arr)'
          read(iunit,iostat=ierr)
       endif
    endif
 enddo

end subroutine read_array_real4arr

!--------------------------------------------------------------------
!+
!  Routine for extracting real*8 array from main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_real8(arr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 real(kind=8),     intent(inout) :: arr(:)
 character(len=*), intent(in)    :: arr_tag,tag
 logical,          intent(inout) :: got_arr
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer      :: i,nread
 real(kind=8) :: dum
 real(kind=4) :: dumr4
 real(kind=4), allocatable :: dummyr4(:)
 logical      :: match_datatype

 if (matched .or. ikind < i_real) return
 match_datatype = (ikind==i_real8 .or. (kind(0.)==8 .and. ikind==i_real))

 if (match_tag(tag,arr_tag) .and. .not.matched) then
    matched    = .true.
    if (match_datatype) then
       got_arr = .true.
       if (i2 > size(arr)) then
          print*,'ERROR: array size too small reading '//trim(tag),i2,size(arr)
          ierr = ierr_arraysize
          return
       endif
       read(iunit,iostat=ierr) (dum,i=1,noffset),arr(i1:i2)
    elseif (ikind==i_real4) then
       got_arr = .true.
       !print*,'WARNING: converting '//trim(tag)//' from real*4'
       nread = i2-i1+1
       allocate(dummyr4(nread))
       read(iunit,iostat=ierr) (dumr4,i=1,noffset),dummyr4(1:nread)
       arr(i1:i2) = real(dummyr4(1:nread),kind=8)
       deallocate(dummyr4)
    else
       print*,'ERROR: wrong datatype for '//trim(tag)//' (is not real8)'
       read(iunit,iostat=ierr)
    endif
 endif

end subroutine read_array_real8

!--------------------------------------------------------------------
!+
!  Routine for extracting multi-d real*8 array from
!  main block in dump files
!+
!--------------------------------------------------------------------
subroutine read_array_real8arr(arr,arr_tag,got_arr,ikind,i1,i2,noffset,iunit,tag,matched,ierr)
 real(kind=8),     intent(inout) :: arr(:,:)
 character(len=*), intent(in)    :: arr_tag(:),tag
 logical,          intent(inout) :: got_arr(:)
 integer,          intent(in)    :: ikind,i1,i2,noffset,iunit
 logical,          intent(inout) :: matched
 integer,          intent(out)   :: ierr
 integer      :: i,j,nread
 real(kind=8) :: dum
 real(kind=4) :: dumr4
 real(kind=4), allocatable :: dummyr4(:)
 real(kind=8), allocatable :: dummyr8(:)
 logical      :: match_datatype

 if (matched .or. ikind < i_real) return
 match_datatype = (ikind==i_real8 .or. (kind(0.)==8 .and. ikind==i_real))

 do j=1,min(size(arr,dim=1),size(arr_tag))
    if (match_tag(tag,arr_tag(j)) .and. .not.matched) then
       matched    = .true.
       if (match_datatype) then
          got_arr(j) = .true.
          !print*,'MATCHED '//trim(arr_tag(j))
          nread = i2-i1+1
          allocate(dummyr8(nread)) ! to avoid seg fault with ifort
          read(iunit,iostat=ierr) (dum,i=1,noffset),dummyr8(1:nread)
          arr(j,i1:i2) = dummyr8(:)
          deallocate(dummyr8)
       elseif (ikind==i_real4) then
          got_arr(j) = .true.
          !print*,'WARNING: converting '//trim(tag)//' from real*4'
          nread = i2-i1+1
          allocate(dummyr4(nread))
          read(iunit,iostat=ierr) (dumr4,i=1,noffset),dummyr4(1:nread)
          arr(j,i1:i2) = real(dummyr4(:),kind=8)
          deallocate(dummyr4)
       else
          print*,'ERROR: wrong datatype for '//trim(tag)//' (is not real8arr)'
          read(iunit,iostat=ierr)
       endif
    endif
 enddo

end subroutine read_array_real8arr

!-----------------------------------------------------
!+
!  Internal utility to open file and
!  extract minimum necessary
!  information from the file header in order to
!  read arrays from the main file blocks
!  Use instead of open_dumpfile_r
!+
!-----------------------------------------------------
subroutine open_dumpfile_rh(iunit,filename,nblocks,narraylengths,ierr,singleprec,id)
 integer,          intent(in)  :: iunit
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: nblocks,narraylengths,ierr
 character(len=lenid), intent(out), optional :: id
 logical,          intent(in),      optional :: singleprec
 character(len=lenid)  :: fileid
 character(len=lentag) :: tagarr(maxphead)
 integer :: intarr(maxphead)
 integer :: number,i

 if (present(singleprec)) then
    call open_dumpfile_r(iunit,filename,fileid,ierr,singleprec=singleprec,requiretags=.true.)
 else
    call open_dumpfile_r(iunit,filename,fileid,ierr,requiretags=.true.)
 endif
 if (present(id)) id = fileid
 if (ierr /= 0) return

 ! read nblocks from int header
 read(iunit,iostat=ierr) number
 nblocks = 1
 narraylengths = 1
 if (number >= 5) then
    if (number > maxphead) number = maxphead
    read(iunit,iostat=ierr) tagarr(1:number)
    read(iunit,iostat=ierr) intarr(1:number)
    call extract('nblocks',nblocks,intarr,tagarr,number,ierr)
    if (ierr /= 0) nblocks = 1
 elseif (number > 0) then
    nblocks = 1
    read(iunit,iostat=ierr)
    read(iunit,iostat=ierr)
 endif

 ! no need to read rest of header
 do i=1,ndatatypes-1
    call skip_headerblock(iunit,ierr)
    if (ierr /= 0) print*,' error skipping header block'
    if (ierr /= 0) return
 enddo

 read (iunit,iostat=ierr) number
 if (ierr /= 0) return
 narraylengths = number/nblocks

end subroutine open_dumpfile_rh

!-----------------------------------------------------
!+
!  The following routine can be used to read a single
!  array matching a particular tag from the main blocks
!  in the file
!+
!-----------------------------------------------------
subroutine read_array_from_file_r8(iunit,filename,tag,array,ierr,use_block,iprint_in)
 integer,               intent(in) :: iunit
 character(len=*),      intent(in) :: filename
 character(len=*),      intent(in) :: tag
 real(kind=8),          intent(out) :: array(:)
 integer, intent(out) :: ierr
 integer, intent(in), optional :: use_block
 logical, intent(in), optional :: iprint_in
 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,my_block
 logical :: iprint

 character(len=lentag) :: mytag

 if (present(use_block)) then
    my_block = use_block
 else
    my_block = 1 ! match from block 1 by default
 endif

 ! if printing the tags
 if (present(iprint_in)) then
    iprint = iprint_in
 else
    iprint = .true.
 endif

 array = 0.

 ! open file for read and get minimal information from header
 call open_dumpfile_rh(iunit,filename,nblocks,narraylengths,ierr)
 if (ierr /= 0) return

 do iblock = 1,nblocks
    call read_block_header(narraylengths,number8,nums,iunit,ierr)
    do j=1,narraylengths
       if (j==my_block) then
          do i=1,ndatatypes
             !print*,' data type ',i,' arrays = ',nums(i,j)
             do k=1,nums(i,j)
                if (i==i_real) then
                   read(iunit,iostat=ierr) mytag
                   if (trim(mytag)==trim(tag)) then
                      read(iunit,iostat=ierr) array(1:min(int(number8(j)),size(array)))
                      if (iprint) print*,'->',mytag
                   else
                      if (iprint) print*,'  ',mytag
                      read(iunit,iostat=ierr)
                   endif
                else
                   read(iunit,iostat=ierr) mytag ! tag
                   read(iunit,iostat=ierr) ! array
                endif
             enddo
          enddo
       endif
    enddo
 enddo

 close(iunit)

end subroutine read_array_from_file_r8

!-----------------------------------------------------
!+
!  The following routine can be used to read a single
!  array matching a particular tag from the main blocks
!  in the file
!+
!-----------------------------------------------------
subroutine read_array_from_file_r4(iunit,filename,tag,array,ierr,use_block,iprint_in)
 integer,               intent(in) :: iunit
 character(len=*),      intent(in) :: filename
 character(len=*),      intent(in) :: tag
 real(kind=4), intent(out) :: array(:)
 integer, intent(out) :: ierr
 integer, intent(in), optional :: use_block
 logical, intent(in), optional :: iprint_in
 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,my_block
 character(len=lentag) :: mytag
 logical :: iprint

 if (present(use_block)) then
    my_block = use_block
 else
    my_block = 1 ! match from block 1 by default
 endif

 ! if printing the tags
 if (present(iprint_in)) then
    iprint = iprint_in
 else
    iprint = .true.
 endif

 array = 0.

 ! open file for read
 call open_dumpfile_rh(iunit,filename,nblocks,narraylengths,ierr)
 if (ierr /= 0) return

 do iblock = 1,nblocks
    call read_block_header(narraylengths,number8,nums,iunit,ierr)
    do j=1,narraylengths
       if (j==my_block) then
          do i=1,ndatatypes
             !print*,' data type ',i,' arrays = ',nums(i,j)
             do k=1,nums(i,j)
                if (i==i_real4) then
                   read(iunit,iostat=ierr) mytag
                   if (trim(mytag)==trim(tag)) then
                      read(iunit,iostat=ierr) array(1:min(int(number8(j)),size(array)))
                      if (iprint) print*,'->',mytag
                   else
                      if (iprint) print*,'  ',mytag
                      read(iunit,iostat=ierr)
                   endif
                else
                   read(iunit,iostat=ierr) mytag ! tag
                   read(iunit,iostat=ierr) ! array
                endif
             enddo
          enddo
       endif
    enddo
 enddo

 close(iunit)

end subroutine read_array_from_file_r4

!-------------------------------------------------------
!+
!  print array tags structure
!+
!-------------------------------------------------------
subroutine print_arrays_in_file(iunit,filename)
 integer,          intent(in) :: iunit
 character(len=*), intent(in) :: filename
 integer :: ierr,nblocks,narraylengths
 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths),nread
 character(len=lentag) :: mytag
 character(len=lenid)  :: fileid
 character(len=4) :: str
 integer, parameter :: ndisplay = 5
 real            :: x(ndisplay)
 real(kind=4)    :: x4(ndisplay)
 integer(kind=1) :: i1(ndisplay)
 logical         :: singleprec

 singleprec = .false.

 ! open file for read
 call open_dumpfile_rh(iunit,filename,nblocks,narraylengths,ierr,id=fileid)
 if (ierr == ierr_realsize) then
    close(iunit)
    singleprec = .true.
    call open_dumpfile_rh(iunit,filename,nblocks,narraylengths,ierr,singleprec=singleprec,id=fileid)
 endif
 if (ierr /= 0) return

 print "(a)",trim(fileid)
 print "(a,i3,a,i2)",':: nblocks = ',nblocks,' array lengths per block = ',narraylengths
 do iblock = 1,nblocks
    call read_block_header(narraylengths,number8,nums,iunit,ierr)
    do j=1,narraylengths
       print "(a,i3,a,i3,a,i10)",'Block ',iblock,' array block ',j,': size ',number8(j)
       nread= min(ndisplay,int(number8(j)))
       str = ' ...'
       if (nread >= int(number8(j))) str = ']'
       do i=1,ndatatypes
          do k=1,nums(i,j)
             read(iunit,iostat=ierr) mytag
             select case(i)
             case(i_int1)
                read(iunit,iostat=ierr) i1(1:nread)
                print*,mytag,datatype_label(i),' [',i1(1:nread),str
             case(i_real)
                if (singleprec) then
                   read(iunit,iostat=ierr) x4(1:nread)
                   print*,mytag,datatype_label(i),' [',x4(1:nread),str
                else
                   read(iunit,iostat=ierr) x(1:nread)
                   print*,mytag,datatype_label(i),' [',x(1:nread),str
                endif
             case(i_real4)
                read(iunit,iostat=ierr) x4(1:nread)
                print*,mytag,datatype_label(i),' [',x4(1:nread),str
             case default
                print*,mytag,datatype_label(i)
                read(iunit,iostat=ierr) ! skip actual array
             end select
          enddo
       enddo
    enddo
 enddo
 close(iunit)

end subroutine print_arrays_in_file

end module dump_utils
