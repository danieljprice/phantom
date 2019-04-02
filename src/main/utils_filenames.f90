!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: fileutils
!
!  DESCRIPTION:
!  This module contains useful utilities related to
!  file names and numbering
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------

module fileutils
 implicit none
 public :: getnextfilename,numfromfile,basename,get_ncolumns,skip_header
 public :: strip_extension,is_digit,files_are_sequential
 public :: ucase,lcase,make_tags_unique

 private

contains
!----------------------------------------------------------------
!+
!  this function increments the dumpfile name based on the
!  number at the end of the filename
!+
!----------------------------------------------------------------
function getnextfilename(filename)
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: getnextfilename
 integer :: idot,istartnum,ilen,i,ierr,num
 character(len=10) :: fmtstring
!
!--extract current number from filename
!
 idot = get_idot(filename)
 istartnum = 0
 do i=idot-1,1,-1
    if (istartnum==0) then
       if (.not.is_digit(filename(i:i))) istartnum = i
    endif
 enddo
 if (istartnum /= 0) istartnum = istartnum + 1
 ilen = idot - istartnum
 if (ilen > 0) then
    read(filename(istartnum:istartnum+ilen-1),*,iostat=ierr) num
    if (ierr /= 0) then
       print*,'internal error in getnextfilename'
       getnextfilename = trim(filename)//'001'
       return
    endif
!
!--increment number by one
!
    num = num + 1
    write(fmtstring,"('(i',i1,'.',i1,')')") ilen,ilen
    getnextfilename = trim(filename)
!
!--replace number in new filename
!
    write(getnextfilename(istartnum:istartnum+ilen-1),fmtstring) num
 else
    getnextfilename = trim(filename)//'001'
 endif

end function getnextfilename

!----------------------------------------------------------------
!+
!  this function extracts the number at the end of the filename
!  (cut down version of previous routine)
!+
!----------------------------------------------------------------
integer function numfromfile(filename)
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: string
 integer :: idot,istartnum,ilen,i,ierr
!
!--extract current number from filename
!
 string = basename(filename)
 idot = get_idot(string)
 istartnum = 0
 do i=idot-1,1,-1
    if (istartnum==0) then
       if (.not.is_digit(string(i:i))) istartnum = i
    endif
 enddo
 if (istartnum /= 0) istartnum = istartnum + 1
 ilen = idot - istartnum

 if (ilen > 0) then
    read(string(istartnum:istartnum+ilen-1),*,iostat=ierr) numfromfile
    if (ierr /= 0) then
       !print*,'internal error in numfromfilename'
       numfromfile = -1
    endif
 else
    numfromfile = 0
 endif

 return
end function numfromfile

!----------------------------------------------------------------
!+
!  see if files have been written as a sequence
!  (according to the numbering from getnextfilename routine)
!+
!----------------------------------------------------------------
logical function files_are_sequential(filenames)
 character(len=*) :: filenames(:)
 integer :: i

 if (size(filenames)<=1) then
    files_are_sequential = .false.  ! return false for just one file
 else
    files_are_sequential = .true.
 endif

 do i=2,size(filenames)
    if (.not.(filenames(i)==getnextfilename(filenames(i-1)))) files_are_sequential = .false.
 enddo

end function files_are_sequential

!---------------------------------------------------------------------------
!+
!  whether or not a character is a number
!+
!---------------------------------------------------------------------------
pure logical function is_digit(ch)
 character(len=1), intent(in) :: ch

 is_digit = (iachar(ch) >= iachar('0') .and. iachar(ch) <= iachar('9'))

end function is_digit

!---------------------------------------------------------------------------
!+
!  extract the start of the file extension, if the filename does not
!  end with digits
!+
!---------------------------------------------------------------------------
pure integer function get_idot(string)
 character(len=*), intent(in) :: string
 integer :: ilen

 ilen = len_trim(string)
 get_idot = 0
 !
 ! if file ends in at least two numbers then use the numbers at the end
 ! (two is to avoid problems with .hdf5 etc)
 !
 if (ilen >= 2) then
    if (is_digit(string(ilen:ilen)) .and. is_digit(string(ilen-1:ilen-1))) then
       get_idot = ilen + 1
    endif
 endif
 !
 ! otherwise, look for numbers before the file extension (e.g. _0000.dat)
 !
 if (get_idot==0) then
    get_idot = index(string,'.',back=.true.)
    if (get_idot==0) get_idot = len_trim(string) + 1
 endif

end function get_idot

!---------------------------------------------------------------------------
!
! function stripping the directory off a filename
!
!---------------------------------------------------------------------------
function basename(string)
 character(len=*), intent(in) :: string
 character(len=len(string)) :: basename
 integer :: i,iposmax

 basename = string

 !--find the last forward slash
 iposmax = 0
 i = len_trim(string)
 do while(i >= 2 .and. iposmax==0)
    i = i - 1
    if (string(i:i)=='/') iposmax = i
 enddo
 basename = trim(string(iposmax+1:))

end function basename

!---------------------------------------------------------------------------
!
! routine to strip extension from a filename
!
!---------------------------------------------------------------------------
subroutine strip_extension(string,ext)
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: ext
 integer :: iloc

 ! find extension in the filename
 iloc = index(string,ext)
 if ( iloc > 0 ) then
    string = string(1:iloc-1)
 endif

end subroutine strip_extension

!---------------------------------------------------------------------------
! utility to work out number of columns of real numbers
! in an ascii file
!
! file must already be open and at the start
! slightly ad-hoc but its the best way I could think of!
!---------------------------------------------------------------------------
subroutine get_ncolumns(lunit,ncolumns,nheaderlines)
 integer, intent(in)  :: lunit
 integer, intent(out) :: ncolumns,nheaderlines
 integer :: ierr,ncolprev,ncolsthisline
 character(len=2000) :: line
 logical :: nansinfile,infsinfile

 nheaderlines = 0
 line = ' '
 ierr = 0
 ncolumns = 0
 ncolprev = 666
 ncolsthisline = 0
 nansinfile = .false.
 infsinfile = .false.
!
!--loop until we find two consecutive lines with the same number of columns (but non zero)
!
 do while ((len_trim(line)==0 .or. ncolsthisline /= ncolprev .or. ncolumns <= 0) .and. ierr==0)
    ncolprev = ncolumns
    read(lunit,"(a)",iostat=ierr) line
    if (index(line,'NaN') > 0) nansinfile = .true.
    if (index(line,'Inf') > 0) infsinfile = .true.
    if (ierr==0) ncolsthisline = ncolumnsline(line)
    if (ncolsthisline >= 0) nheaderlines = nheaderlines + 1
    ncolumns = ncolsthisline
 enddo
 !--subtract 2 from the header line count (the last two lines which were the same)
 nheaderlines = max(nheaderlines - 2,0)
 if (ierr  > 0 .or. ncolumns <= 0) then
    ncolumns = 0
 elseif (ierr  <  0) then
    !print*,ncolumns,ncolprev
 endif
 if (nansinfile) print "(a)",' INDIAN BREAD WARNING!! NaNs in file!!'
 if (infsinfile) print "(a)",' WARNING!! Infs in file!!'
 rewind(lunit)

 if (ncolumns==0) then
    print "(a)",' ERROR: no columns of real numbers found'
 else
    print "(a,i3)",' number of data columns = ',ncolumns
 endif

end subroutine get_ncolumns

!---------------------------------------------------------------------------
!
! function returning the number of columns of real numbers from a given line
!
!---------------------------------------------------------------------------
integer function ncolumnsline(line)
 character(len=*), intent(in) :: line
 real :: dummyreal(100)
 integer :: ierr,i

 dummyreal = -666.0

 ierr = 0
 read(line,*,iostat=ierr) (dummyreal(i),i=1,size(dummyreal))
 !if (ierr  >  0) then
 !   ncolumnsline = -1
 !   return
 !endif

 i = 1
 ncolumnsline = 0
 do while(abs(dummyreal(i)+666.) > tiny(0.))
    ncolumnsline = ncolumnsline + 1
    i = i + 1
    if (i > size(dummyreal)) then
       print "(a)",'*** ERROR: too many columns in file'
       ncolumnsline = size(dummyreal)
       return
    endif
 enddo

end function ncolumnsline

!---------------------------------------------------------------------------
!
! skips header lines in file (nheader can be obtained from get_ncolumns)
!
!---------------------------------------------------------------------------
subroutine skip_header(iunit,nheader,ierror)
 integer, intent(in)  :: iunit,nheader
 integer, intent(out), optional :: ierror
 integer :: i,ierr

 do i=1,nheader
    read(iunit,*,iostat=ierr)
    if (present(ierror)) then
       if (ierr /= 0) ierror = ierr
    else
       if (ierr > 0) then
          print*,'ERROR skipping header line ',i
       else
          print*,'ERROR! end-of-file reading header line ',i
       endif
    endif
 enddo

end subroutine skip_header

!---------------------------------------------------------------------------
!
! Converts a string to upper case
!
!---------------------------------------------------------------------------
function ucase(string)
 character(len=*), intent(in) :: string
 character(len=len(string))   :: ucase
 integer :: is,ia
 integer, parameter           :: aoffset = 32

 ucase = string
 do is = 1, len(ucase)
    ia = iachar(ucase(is:is))
    if (ia >= iachar('a').and.ia <= iachar('z')) &
        ucase(is:is) = achar(ia-aoffset)
 enddo

end function ucase

!---------------------------------------------------------------------------
!
! Converts a string to lower case
!
!---------------------------------------------------------------------------
function lcase(string)
 character(len=*), intent(in) :: string
 character(len=len(string))   :: lcase
 integer :: is,ia
 integer, parameter           :: aoffset = 32

 lcase = string
 do is = 1, len(lcase)
    ia = iachar(lcase(is:is))
    if (ia >= iachar('A').and.ia <= iachar('Z')) &
        lcase(is:is) = achar(ia+aoffset)
 enddo

end function lcase

!------------------------------
! Append a number to a string
! e.g. string,2 -> string2
!------------------------------
subroutine append_number(string,j)
 character(len=*), intent(inout) :: string
 integer,          intent(in)    :: j
 character(len=12) :: strj

 write(strj,"(i12)") j
 string = trim(string)//trim(adjustl(strj))

end subroutine append_number

!----------------------------------------------------------------------
! Append numbers to otherwise identical tags to make them unique
! e.g. massoftype1, massoftype2, massoftype3, etc.
!----------------------------------------------------------------------
subroutine make_tags_unique(ntags,tags)
 integer, intent(in) :: ntags
 character(len=*), dimension(ntags), intent(inout) :: tags
 character(len=len(tags)) :: tagprev
 integer :: i,j

 if (ntags < 1) return
 j = 0
 tagprev = tags(1)
 do i=2,ntags
    if (tags(i)==tagprev) then
       j = j + 1
       if (j==1) then
          call append_number(tags(i-1),j)
          j = j + 1
       endif
       call append_number(tags(i),j)
    else
       tagprev = tags(i)
       j = 0
    endif
 enddo

end subroutine make_tags_unique

end module fileutils
