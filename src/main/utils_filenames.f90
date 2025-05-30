!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module fileutils
!
! This module contains useful utilities related to
!  file names and numbering
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io
!

 implicit none
 public :: getnextfilename,numfromfile,basename,get_ncolumns,skip_header,number_of_rows
 public :: read_column_labels,get_column_labels,split
 public :: strip_extension,is_digit,files_are_sequential,load_data_file
 public :: ucase,lcase,make_tags_unique,get_nlines,string_delete,string_replace,nospaces
 integer, parameter :: max_line_length = 10000 ! for finding number of columns

 private

contains
!----------------------------------------------------------------
!+
!  this function increments the dumpfile name based on the
!  number at the end of the filename
!+
!----------------------------------------------------------------
function getnextfilename(filename,ifilename)
 character(len=*), intent(in) :: filename
 character(len=len(filename)) :: getnextfilename
 integer :: idot,istartnum,ilen,i,ierr,num
 integer, optional, intent(out) :: ifilename
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
 if (present(ifilename)) ifilename = num

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
! function to count number of lines in a file. It has optinal arguments for
! number of headerlines and number of columns that can be returned.
!
!---------------------------------------------------------------------------
function get_nlines(string,skip_comments,n_columns,n_headerlines) result(n)
 character(len=*), intent(in) :: string
 integer :: n,iunit,ierr
 logical, optional, intent(in) :: skip_comments
 logical :: do_skip
 integer :: ncolumns,nheaderlines
 integer, optional, intent(out) :: n_columns
 integer, optional, intent(out) :: n_headerlines

 open(newunit=iunit,file=string,status='old',iostat=ierr)
 do_skip = .false.
 if (present(skip_comments)) do_skip = skip_comments

 if (do_skip .or. present(n_columns) .or. present(n_headerlines)) then
    call get_ncolumns(iunit,ncolumns,nheaderlines)
    if (present(n_columns)) n_columns = ncolumns
    if (present(n_headerlines)) n_headerlines = nheaderlines
 endif

 if (do_skip) call skip_header(iunit,nheaderlines,ierr)

 !--reading number of rows in file
 n = 0
 do while (ierr==0)
    n = n + 1
    read(iunit,*,iostat=ierr)
    if (ierr /= 0) n = n - 1
 enddo
 close(iunit)

end function get_nlines

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

!----------------------------------------------------------------------
! Delete a symbol from a string
!----------------------------------------------------------------------
pure subroutine string_delete(string,skey)
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skey
 integer :: ipos,lensub

 ipos = index(string,skey)
 lensub = len(skey)
 do while(ipos > 0)
    string = string(1:ipos-1)//string(ipos+lensub:len_trim(string))
    ipos = index(trim(string),skey)
 enddo
end subroutine string_delete

!---------------------------------------------------------------------------
!
! subroutine to replace a matching section of a string with another
! string, possibly of differing length
!
!---------------------------------------------------------------------------
subroutine string_replace(string,skey,sreplacewith)
 character(len=*), intent(inout) :: string
 character(len=*), intent(in)    :: skey,sreplacewith
 character(len=len(string)) :: remstring
 integer :: ipos,imax,lensub,i

 ipos = index(trim(string),skey)
 lensub = len(skey)
 imax   = len(string)
 i = 0
 do while (ipos > 0 .and. i <= imax)
    i = i + 1  !  only allow as many replacements as characters
    remstring = string(ipos+lensub:len_trim(string))
    string = string(1:ipos-1)//sreplacewith//remstring
    ipos = index(trim(string),skey)
 enddo

end subroutine string_replace

!---------------------------------------------------------------------------
!
! Split a string into substrings based on a delimiter
!
!---------------------------------------------------------------------------
pure subroutine split(string,delim,stringarr,nsplit)
 character(len=*), intent(in)  :: string
 character(len=*), intent(in)  :: delim
 character(len=*), intent(out), dimension(:) :: stringarr
 integer,          intent(out) :: nsplit
 integer :: i,j,imax,iend

 i = 1
 nsplit = 0
 imax = len(string)
 do while(nsplit < size(stringarr) .and. i <= imax)
    ! find next non-blank character
    if (string(i:i)==' ') then
       do while (string(i:i)==' ')
          i = i + 1
          if (i > imax) exit
       enddo
       if (i > imax) exit
    endif

    ! look for next occurrence of delimiter
    j = index(string(i:),delim) - 1
    ! if no delimiter found, use whole rest of string
    if (j < 0) j = imax
    ! set end of substring
    iend = min(i+j-1,imax)
    ! extract the substring
    nsplit = nsplit + 1
    if (nsplit <= size(stringarr)) then
       stringarr(nsplit) = string(i:iend)
    endif
    i = iend + len(delim) + 1
 enddo

end subroutine split

!-----------------------------------------------------------------
!
!  utility to count number of times a character appears in a string
!
!-----------------------------------------------------------------
integer function count_char(string,mychar)
 character(len=*), intent(in) :: string
 character(len=1), intent(in) :: mychar
 integer :: i

 count_char = 0
 do i=1,len(string)
    if (string(i:i)==mychar) count_char = count_char + 1
 enddo

end function count_char

!---------------------------------------------------------------------------
!
! extract a list of labels from the header line of a file
!
!---------------------------------------------------------------------------
subroutine get_column_labels(line,nlabels,labels,method,ndesired,csv)
 character(len=*), intent(in)  :: line
 integer,          intent(out) :: nlabels
 character(len=*), dimension(:), intent(out) :: labels
 integer,          intent(out), optional :: method
 integer,          intent(in),  optional :: ndesired
 logical,          intent(in),  optional :: csv
 integer :: i1,i2,i,nlabelstmp,nlabels_prev,istyle,ntarget
 character(len=1) :: leadingchar
 character(len=4), parameter :: spaces = '    '
 logical :: is_csv

 nlabels = 0
 i1 = 1
 istyle = 0
 ntarget = -1
 is_csv = .false.
 if (present(csv)) is_csv = csv
 if (present(ndesired)) ntarget = ndesired
 !
 ! strip leading comment character ('#')
 !
 leadingchar = trim(adjustl(line))
 if (leadingchar=='#') then
    i1 = index(line,'#') + 1
 endif
 ! strip anything preceding an equals sign
 i1 = max(i1,index(line,'=')+1)
 i2 = i1

 if (index(nospaces(line),'][') > 0 .and. .not.is_csv) then
    !
    ! format style 1: # [ mylabel1 ] [ mylabel2 ] [ mylabel3 ]
    !
    istyle = 1
    i1 = max(index(line,'[')+1,i1)    ! strip leading square bracket
    call split(nospaces(line(i1:)),'][',labels,nlabels)
 elseif (index(line,',') > 1 .or. is_csv) then
    !
    ! format style 2: mylabel1,mylabel2,mylabel3
    !
    istyle = 2
    call split(line(i1:),',',labels,nlabelstmp)
    if (is_csv) then
       nlabels = nlabelstmp  ! allow blank/arbitrary labels in csv format
    else
       nlabels = count_sensible_labels(nlabelstmp,labels)
    endif
 else
    !
    ! format style 3: #     mylabel1     mylabel2     mylabel3
    !
    istyle = 3
    ! try splitting with 4, then 3, then 2 spaces until the number of labels decreases
    nlabels_prev = 0
    over_spaces: do i=4,2,-1
       call split(line(i1:),spaces(1:i),labels,nlabelstmp)
       ! quit if we already have the target number of labels
       if (nlabelstmp == ntarget) exit over_spaces

       ! if the number of labels is > 1 but has decreased, quit, unless nlabels
       ! still exceeds the number of labels we are hoping for (ntarget)
       if ((nlabelstmp < nlabels_prev .or. nlabelstmp >= max(nlabels_prev,2)  &
            .and. i < 4 .and. .not. (ntarget > 0 .and. nlabelstmp > ntarget))) then
          ! take the answer with the previous number of spaces
          call split(line(i1:),spaces(1:i+1),labels,nlabelstmp)
          exit over_spaces
       endif
       nlabels_prev = nlabelstmp
    enddo over_spaces
    !
    ! this style is dangerous, so perform sanity checks
    ! on the labels to ensure they are sensible
    !
    nlabels = count_sensible_labels(nlabelstmp,labels)
    if (nlabels <= 1) then
       !
       ! format style 4: x y z vx vy vz
       ! (this style is also dangerous)
       !
       istyle = 4
       call split(line(i1:),' ',labels,nlabelstmp)
       nlabels = count_sensible_labels(nlabelstmp,labels)
    endif
 endif
 if (present(method)) method = istyle
 !
 ! clean up
 !
 do i=1,nlabels
    ! delete brackets
    if (nlabels <= size(labels)) then
       call string_delete(labels(i),',')
       if (istyle==1 .or. istyle==2) then
          labels(i) = trim(adjustl(labels(i)))
          ! delete leading numbers
          i1 = 1
          do while (isdigit(labels(i)(i1:i1)))
             labels(i)(i1:i1) = ' '
             i1 = i1 + 1
          enddo
       endif
       labels(i) = trim(adjustl(labels(i)))
    endif
 enddo
 ! delete loose trailing square bracket but only if not matching
 if (istyle==1) then
    if (index(labels(nlabels),']') > 0) then
       i1 = count_char(labels(nlabels),'[') ! number of open brackets
       i2 = count_char(labels(nlabels),']') ! number of closed brackets
       if (i2 > i1) then ! if brackets do not match
          ! find last trailing bracket
          i2 = index(labels(nlabels),']',back=.true.)
          ! delete it, but only if followed by spaces
          if (i2==len_trim(labels(nlabels))) then
             labels(nlabels) = labels(nlabels)(1:i2-1)
          endif
       endif
    endif
 endif

end subroutine get_column_labels

!---------------------------------------------------------------------------
!
! interface to the above routine that also searches for the line
! containing the column labels in the list of header lines
!
!---------------------------------------------------------------------------
subroutine read_column_labels(iunit,nheaderlines,ncols,nlabels,labels,csv,debug)
 integer,          intent(in)  :: iunit,nheaderlines,ncols
 integer,          intent(out) :: nlabels
 character(len=*), dimension(:), intent(out) :: labels
 logical, intent(in), optional :: csv,debug
 character(len=len(labels(1))), dimension(size(labels)) :: tmplabel
 character(len=max_line_length) :: line
 logical :: is_csv,verbose,got_labels
 integer :: i,imethod,ierr,nwanted

 is_csv = .false.
 verbose = .false.
 if (present(csv)) is_csv = csv
 if (present(debug)) verbose = debug
 got_labels = .false.
 nlabels = 0
 nwanted = min(ncols,size(labels)) ! can either retrieve all labels or completely fill the labels array
 labels = ''
 rewind(iunit)
 do i=1,nheaderlines
    read(iunit,"(a)",iostat=ierr) line
    !--try to match column labels from this header line, if not already matched (or dubious match)
    call get_column_labels(trim(line),nlabels,tmplabel,method=imethod,ndesired=nwanted,csv=csv)
    !--if we get nlabels > ncolumns, use them, but keep trying for a better match
    if ((got_labels .and. nlabels == nwanted) .or. &
        (.not.got_labels .and. nlabels >= nwanted  & ! only allow single-spaced labels if == ncols
         .and. (.not.(imethod>=4) .or. nlabels==nwanted))) then
       labels(1:nwanted) = tmplabel(1:nwanted)
       got_labels = .true.
    endif
    if (verbose) print "(5(1x,a,i0))",'DEBUG: line ',i,'nlabels = ',nlabels,&
                 'want ',ncols,'method=',imethod,'len_trim(line)=',len_trim(line) !,' LABELS= '//tmplabel(1:ncols)
 enddo

end subroutine read_column_labels

!---------------------------------------------------------------------------
!
! indicate if a character is a digit (number) or not
!
!---------------------------------------------------------------------------
pure elemental logical function isdigit(string)
 character(len=1), intent(in) :: string
 integer :: ia

 isdigit = .false.
 ia = iachar(string)
 if (ia >= iachar('0').and.ia <= iachar('9')) isdigit = .true.

end function isdigit

!---------------------------------------------------------------------------
!
! count the number of sensible labels in a list of possible labels
!
!---------------------------------------------------------------------------
integer function count_sensible_labels(n,labels) result(m)
 integer, intent(in) :: n
 character(len=*), dimension(n), intent(in) :: labels
 integer :: i

 m = 0
 do i=1,n
    if (is_sensible_label(labels(i))) m = m + 1
 enddo

end function count_sensible_labels

!---------------------------------------------------------------------------
!
! determine if a particular string makes sense as a column label or not
!
!---------------------------------------------------------------------------
logical function is_sensible_label(string)
 character(len=*), intent(in) :: string
 real    :: dum
 integer :: ierr
 real, parameter :: dum_prev =  -66666666.

 is_sensible_label = .true.

 ! should not start with a decimal point
 if (string(1:1)=='.') is_sensible_label = .false.

 ! should not contain equals sign
 !if (index(string,'=') > 0) is_sensible_label = .false.

 dum = dum_prev
 ! should not be able to read it as a real number
 read(string,*,iostat=ierr) dum
 if (ierr==0 .and. abs(dum-dum_prev) > tiny(dum)) is_sensible_label = .false.

end function is_sensible_label

!---------------------------------------------------------------------------
!
! function to strip spaces out of a string
!
!---------------------------------------------------------------------------
function nospaces(string)
 character(len=*), intent(in) :: string
 character(len=len(string)) :: nospaces

 nospaces = string
 call string_delete(nospaces,' ')

end function nospaces

subroutine load_data_file(namefile,datafile,nhead)
 use io,      only:fatal
 character(len=*), intent(in) :: namefile
 integer, intent(in), optional :: nhead
 integer           :: nrows,mcolumns,nheadlines,iunit,ierr,i
 character :: c
 real, dimension(:,:), allocatable, intent(inout) :: datafile
 !N.B.: datafile will be deallocated in grids_for_setup.f90:deallocate_sigma()

 ierr=0
 iunit=155
 !--specify number of headlines
 nheadlines=0

 write(*,*) 'Loading data from file:',namefile

 open(unit=iunit,file=namefile,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    if (trim(namefile)=='sigma_grid.dat' .or. trim(namefile)=='ecc_grid.dat') then
       print*,''
       print*,'!!!!!FATAL!!!!!'
       print*,'You chose to initialise sigma or ecc profiles from files, but there are no such files!'
       print*,'Make sure you ran phantomdir/scripts/generate_eccsigma_grid.py before phantomsetup'
    endif
    call fatal('load_from_file','could not open/read '//trim(namefile))
 endif

 !mcolumns=!number_of_columns(iunit,nheadlines)
 call get_ncolumns(iunit,mcolumns,nheadlines)
 nrows=number_of_rows(iunit)

 if (present(nhead)) then
    nheadlines=nhead
 endif
 write(*,*) 'Skipping ',nheadlines,' head lines'


 write(*,*) 'Found nrows, mcolumns:',nrows,mcolumns


 allocate(datafile(nrows-nheadlines,mcolumns))
 do i=1,nheadlines
    read(iunit,*) c
 enddo
 do i=1, nrows-nheadlines
    read(iunit,*) datafile(i,:)
 enddo
 close(iunit)

end subroutine load_data_file

integer function number_of_rows(s) result(nrows)
 !! version: experimental
 !!
 !! determine number or rows
 integer,intent(in)::s

 integer :: ios
 character  :: r

 rewind(s)
 nrows = 0
 do
    read(s, *,iostat=ios) r
    if (ios /= 0) exit
    nrows = nrows + 1
 enddo

 rewind(s)

end function number_of_rows

!subroutine write_in_file_1d(namefile,arraytowrite)
!    character(len=*), intent(in) :: namefile
!    real, intent(in), dimension(:) :: arraytowrite
!    integer :: iunit,i
!
!    iunit=155
!
!    open(unit=iunit,file=namefile,status='replace',action='write')
!    do i=1,size(arraytowrite(:))
!       write(iunit,*) arraytowrite(i)
!    enddo
!    close(unit=iunit)
!
!end subroutine write_in_file_1d
!
!
!subroutine write_in_file_2d(namefile,arraytowrite)
!    character(len=*), intent(in) :: namefile
!    real, intent(in), dimension(:,:) :: arraytowrite
!    integer :: iunit,i
!
!    iunit=1
!
!    open(unit=iunit,file=namefile,status='replace',action='write')
!
!    do i=1,size(arraytowrite(:,1))
!       write(iunit,*) arraytowrite(i,:)
!    enddo
!
!    close(unit=iunit)
!
!end subroutine write_in_file_2d
!
!subroutine write_in_file_1dx2(namefile,arraytowrite1,arraytowrite2)
!    character(len=*), intent(in) :: namefile
!    real, intent(in), dimension(:) :: arraytowrite1,arraytowrite2
!    integer :: iunit,i
!
!    iunit=1
!
!    open(unit=iunit,file=namefile,status='replace',action='write')
!    do i=1,size(arraytowrite1(:))
!       write(iunit,*) arraytowrite1(i),arraytowrite2(i)
!    enddo
!    close(unit=iunit)
!
!end subroutine write_in_file_1dx2
!
!these functions are currently not used, they could be wrapped in an interface
! interface write_in_file
!    module procedure  write_in_file_1d,write_in_file_2d,write_in_file_1dx2
! end interface



end module fileutils
