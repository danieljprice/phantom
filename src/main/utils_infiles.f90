!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module infile_utils
!
! This module contains utility routines for reading
!  and writing of input files
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
 public :: write_inopt, read_inopt
 public :: read_next_inopt, get_inopt
 public :: write_infile_series, check_infile, contains_loop, get_optstring
!
! generic interface write_inopt to write an input option of any type
!
 interface write_inopt
  module procedure write_inopt_int,write_inopt_real4,write_inopt_real8,write_inopt_string,write_inopt_logical
 end interface write_inopt
!
! generic interface read_inopt to read an input option of any type
!
 interface read_inopt
  module procedure read_inopt_int,read_inopt_real,read_inopt_string,read_inopt_logical
 end interface read_inopt

!
! generic interface get_inopt to read a specific input variable
!
 interface get_inopt
  module procedure get_inopt_int,get_inopt_real,get_inopt_string,get_inopt_logical
 end interface get_inopt
!
! maximum length for input strings
! (if you change this, must also change format statements below)
!
 integer, parameter, private :: maxlen = 20 ! max length of string containing variable
 integer, parameter, private :: maxlenval = 100 ! max length of string containing value
 integer, parameter, private :: maxlenstring = 120  ! max length of string variable
 integer, parameter, private :: maxlenline   = 120  ! maximum line length
!
! structure for input options database
!
 type inopts
    character(len=maxlen)    :: tag
    character(len=maxlenval) :: val
    logical :: retrieved
 end type inopts

 public :: inopts,open_db_from_file,close_db,errtext

 integer, parameter, private :: &
    ierr_notfound = -1, &
    ierr_inread = 1, &
    ierr_rangemin = 2, &
    ierr_rangemax = 3

 private

contains

!-----------------------------------------------------------------
!+
!  write integer variable to input file
!+
!-----------------------------------------------------------------
subroutine write_inopt_int(ival,name,descript,iunit,ierr)
 integer,          intent(in)  :: ival
 character(len=*), intent(in)  :: name,descript
 integer,          intent(in)  :: iunit
 integer,          intent(out), optional :: ierr
 integer :: ierror

 write(iunit,"(a20,' = ',1x,i10,4x,'! ',a)",iostat=ierror) name,ival,descript
 if (present(ierr)) ierr = ierror

end subroutine write_inopt_int

!-----------------------------------------------------------------
!+
!  write logical variable to input file
!+
!-----------------------------------------------------------------
subroutine write_inopt_logical(lval,name,descript,iunit,ierr)
 logical,          intent(in)  :: lval
 character(len=*), intent(in)  :: name,descript
 integer,          intent(in)  :: iunit
 integer,          intent(out), optional :: ierr
 integer :: ierror

 write(iunit,"(a20,' = ',1x,l10,4x,'! ',a)",iostat=ierror) name,lval,descript
 if (present(ierr)) ierr = ierror

end subroutine write_inopt_logical

!-----------------------------------------------------------------
!+
!  write real variable to input file
!+
!-----------------------------------------------------------------
subroutine write_inopt_real4(rval,name,descript,iunit,ierr,exp,time)
 real(kind=4),     intent(in)  :: rval
 character(len=*), intent(in)  :: name,descript
 integer,          intent(in)  :: iunit
 logical,          intent(in),  optional :: exp,time
 integer,          intent(out), optional :: ierr
 logical :: doexp,dotime
 integer :: nhr,nmin !,nsec
 character(len=6) :: fmtstring
 character(len=14) :: tmpstring
 real(kind=4) :: trem
 integer :: ierror

 doexp = .false.
 if (present(exp)) then
    if (exp) doexp = .true.
 endif
 dotime = .false.
 if (present(time)) then
    if (time) dotime = .true.
 endif

 if (dotime) then
    trem = rval
    nhr = int(trem/3600.)
    if (nhr > 0) trem = trem - nhr*3600._4
    nmin = int(trem/60.)
    if (nmin > 0) trem = trem - nmin*60._4
    !nsec = int(trem)

    write(iunit,"(a20,' = ',5x,i3.3,':',i2.2,4x,'! ',a)",iostat=ierror) &
          name,nhr,nmin,descript
 else
    if (doexp .or. (abs(rval) < 1.e-4 .and. abs(rval) > tiny(rval)) &
              .or. (abs(rval) >= 1.e4)) then
       fmtstring = 'es10.3'
       write(iunit,"(a20,' = ',1x,"//fmtstring//",4x,'! ',a)",iostat=ierror) name,rval,descript
    else
       if (abs(rval) < 1.e-1) then
          write(tmpstring,"(f10.7)",iostat=ierror) rval
          tmpstring = adjustl(strip_zeros(tmpstring,3))
       elseif (abs(rval) >= 1.e1) then
          write(tmpstring,"(g14.7)",iostat=ierror) rval
          tmpstring = adjustl(strip_zeros(tmpstring,0))
       else
          write(tmpstring,"(g14.7)",iostat=ierror) rval
          tmpstring = adjustl(strip_zeros(tmpstring,3))
       endif
       write(iunit,"(a20,' = ',1x,a10,4x,'! ',a)",iostat=ierror) name,adjustr(trim(tmpstring)),descript
    endif
 endif
 if (present(ierr)) ierr = ierror

end subroutine write_inopt_real4

!-----------------------------------------------------------------
!+
!  write real variable to input file
!+
!-----------------------------------------------------------------
subroutine write_inopt_real8(rval,name,descript,iunit,ierr,exp,time)
 real(kind=8),     intent(in)  :: rval
 character(len=*), intent(in)  :: name,descript
 integer,          intent(in)  :: iunit
 logical,          intent(in),  optional :: exp,time
 integer,          intent(out), optional :: ierr
 logical :: doexp,dotime
 integer :: nhr,nmin !,nsec
 character(len=16) :: tmpstring
 real(kind=8) :: trem
 integer :: ierror

 doexp = .false.
 if (present(exp)) then
    if (exp) doexp = .true.
 endif
 dotime = .false.
 if (present(time)) then
    if (time) dotime = .true.
 endif

 if (dotime) then
    trem = rval
    nhr = int(trem/3600.d0)
    if (nhr > 0) trem = trem - nhr*3600.d0
    nmin = int(trem/60.d0)
    if (nmin > 0) trem = trem - nmin*60.d0
    !nsec = int(trem)

    write(iunit,"(a20,' = ',5x,i3.3,':',i2.2,4x,'! ',a)",iostat=ierror) &
         name,nhr,nmin,descript
 else
    if (doexp .or. (abs(rval) < 1.e-3 .and. abs(rval) > tiny(rval)) &
              .or. (abs(rval) >= 1.e4)) then
       write(iunit,"(a20,' = ',1x,es10.3,4x,'! ',a)",iostat=ierror) &
         name,rval,descript
    else
       if (abs(rval) <= 1.e-1) then
          write(tmpstring,"(f16.13)",iostat=ierror) rval
          tmpstring = adjustl(strip_zeros(tmpstring,3))
       elseif (abs(rval) >= 1.e1) then
          write(tmpstring,"(g16.9)",iostat=ierror) rval
          tmpstring = adjustl(strip_zeros(tmpstring,0))
       else
          write(tmpstring,"(g16.9)",iostat=ierror) rval
          tmpstring = adjustl(strip_zeros(tmpstring,3))
       endif
       if (len_trim(tmpstring) > 10) then
          write(iunit,"(a20,' = ',1x,a,2x,'! ',a)",iostat=ierror) name,adjustr(trim(tmpstring)),descript
       else
          write(iunit,"(a20,' = ',1x,a10,4x,'! ',a)",iostat=ierror) name,adjustr(trim(tmpstring)),descript
       endif
    endif
 endif
 if (present(ierr)) ierr = ierror

end subroutine write_inopt_real8

!-----------------------------------------------------------------
!+
!  internal function to strip zeros from the end of
!  a number after the decimal place
!+
!-----------------------------------------------------------------
function strip_zeros(string,n)
 character(len=*), intent(in) :: string
 integer,          intent(in) :: n ! keep at least this many zeros
 character(len=len(string)) :: strip_zeros
 integer :: idot,i

 ! copy string
 strip_zeros = string

 ! do not strip zeros from exponential notation
 if (index(string,'+') > 0 .or. index(string,'E') > 0) return

 idot = index(string,'.',back=.true.)
 if (idot > 0) then  ! allow max of nstrip zeros after decimal place
    over_string: do i=len_trim(string),idot+n+1,-1
       if (string(i:i)=='0') then
          strip_zeros(i:i) = ' '
       else
          exit over_string
       endif
    enddo over_string
 endif
 !strip_zeros = adjustl(strip_zeros)

end function strip_zeros

!-----------------------------------------------------------------
!+
!  write string to input file
!+
!-----------------------------------------------------------------
subroutine write_inopt_string(sval,name,descript,iunit,ierr)
 character(len=*), intent(in)  :: sval,name,descript
 integer,          intent(in)  :: iunit
 integer,          intent(out), optional :: ierr
 character(len=40) :: fmtstring
 integer :: ierror

 if (len_trim(sval) > 10) then
    fmtstring = '(a20,'' = '',1x,a,3x,''! '',a)'
 else
    fmtstring = '(a20,'' = '',1x,a10,4x,''! '',a)'
 endif

 write(iunit,fmtstring,iostat=ierror) name,trim(sval),trim(descript)
 if (present(ierr)) ierr = ierror

end subroutine write_inopt_string

!-----------------------------------------------------------------
!+
!  open an input file, read options into a database and
!  return the handle
!+
!-----------------------------------------------------------------
subroutine open_db_from_file(db,filename,iunit,ierr)
 character(len=*),          intent(in)  :: filename
 type(inopts), allocatable, intent(out) :: db(:)
 integer,                   intent(out) :: ierr
 integer,                   intent(in), optional :: iunit
 character(len=maxlen) :: name,valstring
 integer :: n,i

 open(unit=iunit,file=trim(filename),status='old',form='formatted',iostat=ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR opening '//trim(filename)
    return
 endif
 !
 ! work out the number of valid entries in the input file
 !
 n = 0
 do while(ierr==0)
    call read_next_inopt(name,valstring,iunit,ierr)
    if (ierr==0) n = n + 1
 enddo

 if (n > 0) then
    !
    ! allocate memory for database
    !
    print "(a,i2,a)",' opening database from '//trim(filename)//' with ',n,' entries'
    if (allocated(db)) deallocate(db)
    allocate(db(n),stat=ierr)
    if (ierr /= 0) then
       close(iunit)
       return
    endif
    !
    ! initialise all entries as blank
    !
    db(:)%tag = ' '
    db(:)%val = ' '
    db(:)%retrieved = .false.
    !
    ! read file
    !
    rewind(iunit)
    ierr = 0
    do i=1,n
       call read_next_inopt(db(i)%tag,db(i)%val,iunit,ierr)
       if (ierr /= 0) then
          print*,'INTERNAL ERROR: open_db_from_file: error: ierr should be zero here'
          if (allocated(db)) deallocate(db)
          close(iunit)
          return
       endif
    enddo
 else
    print "(a)",' ERROR: no valid database entries in '//trim(filename)
    ierr = -1
 endif
 close(unit=iunit)

end subroutine open_db_from_file

!-----------------------------------------------------------------
!+
!  close a database (just deallocates memory)
!+
!-----------------------------------------------------------------
subroutine close_db(db)
 type(inopts), allocatable, intent(inout) :: db(:)

 if (allocated(db)) deallocate(db)

end subroutine close_db

!-----------------------------------------------------
!+
!  translation of ierr value from read_inopt routines
!  into human-readable text
!+
!-----------------------------------------------------
function errtext(ierr)
 character(len=25) :: errtext
 integer, intent(in) :: ierr

 select case(ierr)
 case(ierr_notfound)
    errtext = 'not found'
 case(ierr_rangemin)
    errtext = 'too small'
 case(ierr_rangemax)
    errtext = 'too big'
 case(0)
    errtext = 'read OK'
 case default
    errtext = 'unknown error'
 end select

end function errtext

!---------------------------------------------------------
!+
!  print error message for errors during read_inopt calls
!+
!-=-------------------------------------------------------
subroutine print_error(tag,valstring,chmin,chmax,ierr)
 character(len=*), intent(in) :: tag,valstring,chmin,chmax
 integer,          intent(in) :: ierr

 if (ierr == ierr_rangemin .or. ierr == ierr_rangemax) then
    print "(a,1x,a,' [',a,':',a,']')", &
         ' ERROR: '//trim(adjustl(tag))//' = '//trim(adjustl(valstring)), &
                   trim(errtext(ierr)),trim(adjustl(chmin)),trim(adjustl(chmax))
 else
    print "(a,1x,a)",' ERROR: '//trim(adjustl(tag))//' '//trim(errtext(ierr))
 endif

end subroutine print_error

!-----------------------------------------------------------------
!+
!  read an integer variable from an input options database
!+
!-----------------------------------------------------------------
subroutine read_inopt_int(ival,tag,db,err,errcount,min,max)
 integer,                   intent(out)   :: ival
 character(len=*),          intent(in)    :: tag
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(out),   optional :: err
 integer,                   intent(inout), optional :: errcount
 integer,                   intent(in), optional :: min,max
 character(len=maxlen) :: valstring
 character(len=16)     :: chmin,chmax
 integer :: ioerr,ierr

 chmin = ''
 chmax = ''
 valstring = ''

 ierr = 0
 if (match_inopt_in_db(db,tag,valstring)) then
    read(valstring,*,iostat=ioerr) ival
    if (ioerr /= 0) ierr = ierr_inread
 else
    ierr = ierr_notfound
 endif
 if (ierr==0) then
    if (present(min)) then
       write(chmin,"(g10.0)") min
       if (ival < min) ierr = ierr_rangemin
    endif
    if (present(max)) then
       write(chmax,"(g10.0)") max
       if (ival > max) ierr = ierr_rangemax
    endif
 endif

 ! print error message unless err argument is given
 if (present(err)) then
    err = ierr
 elseif (ierr /= 0) then
    call print_error(tag,valstring,chmin,chmax,ierr)
 endif
 if (present(errcount)) then
    if (ierr /= 0) errcount = errcount + 1
 endif

end subroutine read_inopt_int

!-----------------------------------------------------------------
!+
!  read a real variable from an input options database
!+
!-----------------------------------------------------------------
subroutine read_inopt_real(val,tag,db,err,errcount,min,max)
 real,                      intent(out)   :: val
 character(len=*),          intent(in)    :: tag
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(out),   optional :: err
 integer,                   intent(inout), optional :: errcount
 real,                      intent(in), optional :: min,max
 character(len=maxlen) :: valstring
 character(len=16)     :: chmin,chmax
 integer :: ioerr,ierr

 chmin = ''
 chmax = ''
 valstring = ''

 ierr = 0
 if (match_inopt_in_db(db,tag,valstring)) then
    read(valstring,*,iostat=ioerr) val
    if (ioerr /= 0) ierr = ierr_inread
 else
    ierr = ierr_notfound
 endif
 if (ierr==0) then
    if (present(min)) then
       write(chmin,"(g13.4)") min
       if (val < min) ierr = ierr_rangemin
    endif
    if (present(max)) then
       write(chmax,"(g13.4)") max
       if (val > max) ierr = ierr_rangemax
    endif
 endif
 if (present(err)) then
    err = ierr
 elseif (ierr /= 0) then
    call print_error(tag,valstring,chmin,chmax,ierr)
 endif
 if (present(errcount)) then
    if (ierr /= 0) errcount = errcount + 1
 endif

end subroutine read_inopt_real

!-----------------------------------------------------------------
!+
!  read a string variable from an input options database
!+
!-----------------------------------------------------------------
subroutine read_inopt_string(valstring,tag,db,err,errcount)
 character(len=*),          intent(out)   :: valstring
 character(len=*),          intent(in)    :: tag
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(out),   optional :: err
 integer,                   intent(inout), optional :: errcount
 integer :: ierr

 ierr = 0
 if (.not.match_inopt_in_db(db,tag,valstring)) ierr = -1

 if (present(err)) then
    err = ierr
 elseif (ierr /= 0) then
    call print_error(tag,valstring,'','',ierr)
 endif
 if (present(errcount)) then
    if (ierr /= 0) errcount = errcount + 1
 endif

end subroutine read_inopt_string

!-----------------------------------------------------------------
!+
!  read a logical variable from an input options database
!+
!-----------------------------------------------------------------
subroutine read_inopt_logical(lval,tag,db,err,errcount)
 logical,                   intent(out)   :: lval
 character(len=*),          intent(in)    :: tag
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(out),   optional :: err
 integer,                   intent(inout), optional :: errcount
 character(len=maxlen) :: valstring
 integer :: ierr

 ierr = 0
 if (match_inopt_in_db(db,tag,valstring)) then
    read(valstring,*,iostat=ierr) lval
 else
    ierr = -1
 endif
 if (present(err)) then
    err = ierr
 elseif (ierr /= 0) then
    call print_error(tag,valstring,'','',ierr)
 endif
 if (present(errcount)) then
    if (ierr /= 0) errcount = errcount + 1
 endif

end subroutine read_inopt_logical

!-----------------------------------------------------------------
!+
!  extract a matching entry from an input options database
!  returns the string containing the value for the matching tag
!+
!-----------------------------------------------------------------
logical function match_inopt_in_db(db,tag,valstring)
 type(inopts), allocatable, intent(inout) :: db(:)
 character(len=*),          intent(in)    :: tag
 character(len=*),          intent(out)   :: valstring
 integer :: n

 match_inopt_in_db = .false.
 if (.not.allocated(db)) then
    valstring = ' '
    return
 endif

 do n=1,size(db)
    if (trim(db(n)%tag)==trim(adjustl(tag))) then
       match_inopt_in_db = .true.
       valstring = db(n)%val
       !--mark entry in database as retrieved
       db(n)%retrieved = .true.
    endif
 enddo

end function match_inopt_in_db

!-----------------------------------------------------------------
!+
!  read a line from the input file and extract the tag and the
!  string section containing the data value
!+
!-----------------------------------------------------------------
subroutine read_next_inopt(tag,valstring,iunit,ierr,nlinesread)
 character(len=maxlenline)      :: line
 character(len=*), intent(out) :: tag,valstring
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr ! not optional for reads
 integer,          intent(out), optional :: nlinesread
 integer :: ierrtemp

 tag = ' '
 valstring = ' '
 line = ' '
 ierr = 0
 if (present(nlinesread)) nlinesread = 0
 do while((len_trim(line)==0 .or. line(1:1)=='#') .and. ierr==0)
    if (present(nlinesread)) nlinesread = nlinesread + 1
    read(iunit,"(a)",iostat=ierr,end=88) line
 enddo

 call read_inopt_from_line(line,tag,valstring,ierrtemp)
 if (ierrtemp /= 0) ierr = ierrtemp

 return

88 continue
 if (present(nlinesread)) nlinesread = nlinesread - 1
 return

end subroutine read_next_inopt

!-----------------------------------------------------------------
!+
!  extract the tag and the string section containing the data value
!  from the line read from the file
!+
!-----------------------------------------------------------------
subroutine read_inopt_from_line(line,name,valstring,ierr,comment)
 character(len=*), intent(in)  :: line
 character(len=*), intent(out) :: name
 character(len=*), intent(out) :: valstring
 integer,          intent(out) :: ierr
 character(len=*), intent(out), optional :: comment
 integer :: iequal, icomment,nhr,nmin,nsec

 ierr = 0
 name = ' '
 valstring = ' '
!
!--find position of the equals sign
!
 iequal = index(line,'=')    ! location of first '=' in line
 if (iequal==0) then
    ierr = 1
    return
 endif
!
!--the name is what lies to the left of the equals sign
!
 name = line(1:min(iequal-1,len(name),len(line)))
 name = trim(adjustl(name))
!
!--find starting position of the comment (exclamation mark)
!
 icomment = index(line,'!')  ! location of '!' in line
 if (icomment <= 0) icomment = len_trim(line)+1  ! if no '!', read to end of line
!
!--the value is what lies between the equals sign
!  and either the comment or the end of the line
!
 valstring = trim(adjustl(line(iequal+1:icomment-1)))
! print*,' name= '//trim(name)//' val = '//trim(valstring)
!
!--for time strings, assume they are of the form hh:mm:ss
!  convert to a number of seconds as a real
!
 if (index(valstring,':') /= 0) then
    nsec = 0
    read(valstring,"(i3.3,1x,i2.2,1x,i2.2)",iostat=ierr) nhr,nmin,nsec
    if (ierr/=0) then
       read(valstring,"(i2.2,1x,i2.2,1x,i2.2)",iostat=ierr) nhr,nmin,nsec
    endif
    !print*,'hh,mm,ss=',nhr,nmin,nsec,((nhr*60. + nmin)*60. + nsec),ierr
    write(valstring,"(es11.4)",iostat=ierr) ((nhr*60. + nmin)*60. + nsec)
 endif

 if (present(comment)) then
    comment = ' '
    if (len_trim(line) > icomment+1) then
       comment = line(icomment+1:len_trim(line))
    endif
 endif

 return
end subroutine read_inopt_from_line

!-----------------------------------------------------------------
!+
!  internal subroutine: get the value string matching a
!  supplied input tag: ierr = -1 if string not found
!+
!-----------------------------------------------------------------
subroutine get_inopt_string(valstring,varname,infile,iunit,ierr)
 character(len=*), intent(out) :: valstring
 character(len=*), intent(in)  :: varname,infile
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 character(len=maxlenstring)   :: name
 logical :: imatch

 open(unit=iunit,file=infile,form='formatted',status='old',iostat=ierr)
 imatch = .false.
 do while(.not.imatch .and. ierr==0)
    call read_next_inopt(name,valstring,iunit,ierr)
    if (trim(name)==trim(varname)) imatch = .true.
 enddo
 ierr = 0
 if (.not.imatch) ierr = -1
 close(iunit)

 return
end subroutine get_inopt_string

!-----------------------------------------------------------------
!+
!  read a real variable matching the varname tag from an input file
!+
!-----------------------------------------------------------------
subroutine get_inopt_real(rval,varname,infile,iunit,ierr)
 real,             intent(out) :: rval
 character(len=*), intent(in)  :: varname,infile
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 character(len=maxlenstring)        :: valstring

 call get_inopt_string(valstring,varname,infile,iunit,ierr)
 read(valstring,*,iostat=ierr) rval

 return
end subroutine get_inopt_real

!-----------------------------------------------------------------
!+
!  read an integer variable matching the varname tag from an input file
!+
!-----------------------------------------------------------------
subroutine get_inopt_int(ival,varname,infile,iunit,ierr)
 integer,          intent(out) :: ival
 character(len=*), intent(in)  :: varname,infile
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 character(len=maxlenstring)        :: valstring

 call get_inopt_string(valstring,varname,infile,iunit,ierr)
 read(valstring,*,iostat=ierr) ival

 return
end subroutine get_inopt_int

!-----------------------------------------------------------------
!+
!  read a logical variable matching the varname tag from an input file
!+
!-----------------------------------------------------------------
subroutine get_inopt_logical(lval,varname,infile,iunit,ierr)
 logical,          intent(out) :: lval
 character(len=*), intent(in)  :: varname,infile
 integer,          intent(in)  :: iunit
 integer,          intent(out) :: ierr
 character(len=maxlenstring)        :: valstring

 call get_inopt_string(valstring,varname,infile,iunit,ierr)
 read(valstring,*,iostat=ierr) lval

 return
end subroutine get_inopt_logical

!--------------------------------------------------------------------
!+
!  Function that determines whether or not an input file contains
!  a loop (public)
!+
!--------------------------------------------------------------------
subroutine check_infile(infile,lu_read,containsloop,ierrline)
 character(len=*), intent(in)  :: infile
 integer,          intent(in)  :: lu_read ! unit to read file on
 logical,          intent(out) :: containsloop
 integer,          intent(out) :: ierrline
 character(len=maxlen)       :: name
 character(len=maxlenstring) :: valstring
 integer :: nlinesread,line,ierr

 containsloop = .false.
 ierrline = 0
 line = 0
 open(unit=lu_read,iostat=ierr,file=infile,status='old',form='formatted')
 if (ierr /= 0) ierrline = -1
 do while (ierr == 0)
    call read_next_inopt(name,valstring,lu_read,ierr,nlinesread)
    line = line + nlinesread
    if (ierr > 0) ierrline = line
    if (ierr==0 .and. contains_loop(valstring)) containsloop = .true.
 enddo
 close(unit=lu_read)

end subroutine check_infile

!--------------------------------------------------------------------
!+
!  see if there are any loops implied in the value string
!  (presence of ' to ' or ', step')
!+
!--------------------------------------------------------------------
logical function contains_loop(valstring)
 character(len=*), intent(in) :: valstring

 if (index(valstring,' to ')   > 1 .or. &
     index(valstring,' step ') > 1) then
    contains_loop = .true.
 else
    contains_loop = .false.
 endif

end function contains_loop

!--------------------------------------------------------------------
!+
!  extract information from a loop implied in the value string
!+
!--------------------------------------------------------------------
subroutine get_loopinfo_real(valstring,rvalstart,rvalend,rstep,njobs,log,ierr)
 character(len=*), intent(in)  :: valstring
 real,             intent(out) :: rvalstart,rvalend,rstep
 integer,          intent(out) :: njobs,ierr
 logical,          intent(out) :: log
 integer :: ito,istep,idstep,instep,inlogstep,ilogstep

 ito       = index(valstring,' to ')
 idstep    = index(valstring,' step ')
 instep    = index(valstring,' nstep ')
 inlogstep = index(valstring,' nlogstep ')
 ilogstep  = index(valstring,' logstep ')
 ierr  = 0
 njobs = 0
 log   = .false.

 istep = max(idstep,instep,inlogstep,ilogstep)

 if (ito <= 0 .or. istep  <=  ito) then
    ierr = 1
 else
!--check that only one step variable is in the string
    if (idstep*instep > 0 .or. idstep*inlogstep > 0 .or. instep*inlogstep > 0) then
       ierr = 2
       return
    endif
!--extract the start and end values
    read(valstring(1:ito),*,iostat=ierr) rvalstart
    read(valstring(ito+4:istep-1),*,iostat=ierr) rvalend

    if (instep > 0) then
       !
       !--here the number of steps between the end points has been specified
       !
       read(valstring(istep+7:),*,iostat=ierr) njobs
       if (njobs > 1) then
          rstep = (rvalend - rvalstart)/real(njobs-1)
       else
          rstep = 0.
       endif
    elseif (inlogstep > 0) then
       !
       !--here the number of steps is specified and the loop is in log space
       !
       read(valstring(istep+10:),*,iostat=ierr) njobs
       log = .true.
       if (njobs > 1) then
          rstep = (log10(rvalend/rvalstart))/real(njobs-1)
       else
          rstep = 0.
       endif
    elseif (ilogstep > 0) then
       !
       !--here the step interval has been specified in log space
       !
       read(valstring(istep+9:),*,iostat=ierr) rstep
       log = .true.
       if (rvalend <= rvalstart) then
          njobs = 0
       else
          njobs = int((log10(rvalend) - log10(rvalstart) + tiny(0.))/rstep) + 1
       endif
    else
       !
       !--here the step interval has been specified
       !
       read(valstring(istep+6:),*,iostat=ierr) rstep
       if (rvalend <= rvalstart) then
          njobs = 0
       else
          njobs = int((rvalend-rvalstart + tiny(0.))/rstep) + 1
       endif
    endif
 endif

 return
end subroutine get_loopinfo_real

!--------------------------------------------------------------------
!+
!  extract information from a loop implied in the value string
!+
!--------------------------------------------------------------------
subroutine get_loopinfo_int(valstring,ivalstart,ivalend,istep,ierr)
 character(len=*), intent(in)  :: valstring
 integer,          intent(out) :: ivalstart,ivalend,istep,ierr
 integer :: ito,isteppos

 ito      = index(valstring,' to ')
 isteppos = index(valstring,' step ')
 ierr     = 0

 if (ito <= 0 .or. isteppos  <=  ito) then
    ierr = 1
 else
    read(valstring(1:ito),*,iostat=ierr) ivalstart
    read(valstring(ito+4:isteppos-1),*,iostat=ierr) ivalend
    read(valstring(isteppos+6:),*,iostat=ierr) istep
 endif

 return
end subroutine get_loopinfo_int

!--------------------------------------------------------------------
!+
!  work out whether values supplied for a loop are real variables
!+
!--------------------------------------------------------------------
logical function isrealloop(valstring)
 character(len=*), intent(in) :: valstring
 integer :: ito,istep,instep,inlogstep,ilogstep

 isrealloop = .false.
 ito        = index(valstring,' to ')
 istep      = index(valstring,' step ')
 instep     = index(valstring,' nstep ')
 inlogstep  = index(valstring,' nlogstep ')
 ilogstep   = index(valstring,' logstep ')
 if (ito <= 0 .or. (istep  <=  ito .and. instep <= ito &
     .and. inlogstep <= ito .and. ilogstep <= ito)) then
    return
    !istep = max(istep,instep,inlogstep,ilogstep)
!
!--look for decimal places in the string before and after the ' to '
!  and after the ' step '
!
 elseif (index(valstring(1:ito-1),'.') /= 0 .or. &
         index(valstring(ito+4:istep),'.') /= 0 .or. &
         index(valstring(istep+6:),'.') /= 0) then
    isrealloop = .true.
 endif
end function isrealloop

!--------------------------------------------------------------------
!+
!  work out whether values supplied for a loop are integers
!+
!--------------------------------------------------------------------
logical function isintloop(valstring)
 character(len=*), intent(in) :: valstring
 integer :: ito,istep,inum,i,i1,i2
 logical :: isintvalarr(3)

 isintloop   = .false.
 isintvalarr = .false.
 ito         = index(valstring,' to ')
 istep       = index(valstring,' step ')
 if (ito <= 0 .or. istep  <=  ito) then
    return
!
!--look for numbers
!
 else
    i1 = 1
    i2 = 0
    do inum=1,3
       select case(inum)
       case(1)
          i1 = 1
          i2 = ito-1
       case(2)
          i1 = ito+4
          i2 = istep
       case(3)
          i1 = istep+6
          i2 = len(valstring)
       end select
       do i=i1,i2
          select case(valstring(i:i))
          case('0','1','2','3','4','5','6','7','8','9')
             isintvalarr(inum) = .true.
          end select
       enddo
    enddo
    if (all(isintvalarr)) isintloop = .true.
 endif

end function isintloop

!--------------------------------------------------------------------
!+
!  detects loops in the input file and, if present, proceeds to
!  write a series of input files with the variable changed.
!
!  Loops look something like:
!   var = 1 to 10 step 1
!   var = 1 to 100 logstep 2
!
!  where var is either a real or integer input variable
!+
!--------------------------------------------------------------------
subroutine write_infile_series(lu_read,lu_write,infile,nlines,ierr)
 character(len=*), intent(in) :: infile
 integer,          intent(in) :: lu_read,lu_write,nlines
 character(len=120)            :: infiledata(nlines)
 character(len=40)             :: name,valstring,nameloopval,valstringloop
 character(len=60)             :: commentstring,commentloop
 character(len=len(infile)+30) :: infilenew
 integer, parameter :: maxloop = 1000
 integer :: line,ierr,ierr2,idot
 integer :: nloops,loopline,njobs,i,ivalstart,ivalend,istep
 real    :: rval,rvalstart,rvalend,rstep
 logical :: logloop
!
!--first of all, we read the whole input file as plain text
!
 ierr = 0
 open(unit=lu_read,file=infile,status='old',form='formatted',iostat=ierr)
 line = 0
 nloops = 0
 loopline = 0
 do while (ierr == 0 .and. line < nlines)
    line = line + 1
    !
    !--read each line of input file into character array
    !
    read(lu_read,"(a)",iostat=ierr) infiledata(line)
    !
    !--check each line for loop syntax
    !
    if (ierr == 0) then
       call read_inopt_from_line(infiledata(line),name,valstring,ierr2,comment=commentstring)
       if (contains_loop(valstring)) then
          nloops = nloops + 1
          nameloopval = trim(name)
          valstringloop = trim(valstring)
          commentloop = trim(commentstring)
          loopline = line
          write(*,"(' __',/,'/  \',/,'|',2x,a,i2,/,'\--->')") &
               'LOOP DETECTED IN VARIABLE '//trim(name)//' in input file on line ',loopline
       endif
    endif
 enddo
 close(unit=lu_read)

 if (nloops >= 1) then
    if (nloops > 1) then
       write(*,"(a)") 'ERROR! write_infile_series: more than one loop in input file not yet implemented'
       return
    endif
!
!--extract the base name for the input files from the original input file
!
    idot = index(infile,'.in')
    if (idot <= 1) idot = len_trim(infile)
    infilenew = infile(1:idot-1)
!
!--now extract the loop information
!
    if (isrealloop(valstringloop)) then
!
!--loop in a real variable
!
       call get_loopinfo_real(valstringloop,rvalstart,rvalend,rstep,njobs,logloop,ierr)
       if (ierr /= 0) then
          write(*,"(a)") 'ERROR! write_infile_series: invalid syntax in loop string'
          return
       endif
       if (njobs > maxloop) then ! sanity check
          write(*,"(a)") 'ERROR! write_infile_series: too many input files'
          return
       endif
       write(*,"(a,i4,a)") 'WRITING ',njobs,' INPUT FILES VARYING '//trim(nameloopval)
!
!--loop over all parameters, writing each input file separately
!
       do i=1,njobs
          if (logloop) then
             rval = rvalstart*10**((i-1)*rstep)
          else
             rval = rvalstart + (i-1)*rstep
          endif
          call formatreal(rval,valstring,precision=rstep)
          write(infilenew(idot:),"(a)") '_'//trim(nameloopval)//'_'//trim(valstring)//infile(idot:)
          call write_infile_lines(lu_write,infilenew,infiledata,loopline,nameloopval,commentloop,ierr,rval)
          if (ierr /= 0) return
       enddo
    elseif (isintloop(valstringloop)) then
!
!--loop in an integer variable
!
       call get_loopinfo_int(valstringloop,ivalstart,ivalend,istep,ierr)
       if (ierr /= 0) then
          write(*,"(a)") 'ERROR! write_infile_series: invalid syntax in loop string'
          return
       endif
       if ((ivalend-ivalstart)/istep > maxloop .or. istep <= 0) then
          write(*,"(a)") 'ERROR! write_infile_series: too many input files'
          return
       endif
!
!--loop over all parameters, writing each input file separately
!
       do i=ivalstart,ivalend,istep
          write(valstring,*) i
          write(infilenew(idot:),"(a)") '_'//trim(nameloopval)//'_'//trim(adjustl(valstring))//infile(idot:)
          call write_infile_lines(lu_write,infilenew,infiledata,loopline,nameloopval,commentloop,ierr,ival=i)
          if (ierr /= 0) return
       enddo
    else
       write(*,"(a)") 'ERROR! write_infile_series: cannot determine type of loop variable'
       return
    endif
 else
    write(*,"(a)") 'ERROR! write_infile_series: error extracting loops from input file'
    return
 endif

end subroutine write_infile_series

!----------------------------------------------------------------
!+
! internal copy of formatreal routine to avoid dependencies
! formats a real number into a string as neatly as possible
! this is only used when setting the name of the new input file
!+
!----------------------------------------------------------------
subroutine formatreal(val,string,precision)
 real,             intent(in)  :: val
 character(len=*), intent(out) :: string
 real,             intent(in), optional :: precision
 integer :: ierr,i
 real    :: prec

 if (present(precision)) then
    prec = precision
 else
    prec = 1.e3  ! large number
 endif

 if (abs(val) >= 1.d99) then
    write(string,"(es10.3)",iostat=ierr) val
 elseif (abs(val) < 1.e-3 .or. abs(val) >= 1.e4 .or. abs(prec) < 1.e-3) then
    write(string,"(es9.2)",iostat=ierr) val
 elseif (abs(val) < 0.01 .or. abs(prec) < 0.01) then
    write(string,"(f10.4)",iostat=ierr) val
 elseif (abs(val) < 0.1  .or. abs(prec) < 0.1) then
    write(string,"(f9.3)",iostat=ierr) val
 elseif (abs(val) < 1. .or. abs(prec) < 1.) then
    write(string,"(f8.2)",iostat=ierr) val
 elseif (abs(val) < 10. .or. abs(prec) < 10.) then
    write(string,"(f8.1)",iostat=ierr) val
 elseif (abs(val) < 100. .or. abs(prec) < 100.) then
    write(string,"(f8.0)",iostat=ierr) val
 else
    write(string,"(f8.0)",iostat=ierr) val
 endif
 string = adjustl(trim(string))
 !
 !--strip trailing zeros after the decimal place
 !  (and the decimal place if it is the last character)
 !
 i = len_trim(string)
 if (string(i:i)=='.') string(i:i) = ' '
 string = trim(adjustl(string))
 return
end subroutine formatreal

!--------------------------------------------------------------------
!+
!  replaces the lines containing the looped variable whilst
!  writing the new input file
!+
!--------------------------------------------------------------------
subroutine write_infile_lines(iunit,infilenew,infiledata,linenum,nameval,commentval,ierr,rval,ival)
 character(len=*), intent(in)  :: infilenew
 character(len=*), intent(in)  :: infiledata(:)
 integer,          intent(in)  :: iunit,linenum
 character(len=*), intent(in)  :: nameval,commentval
 integer,          intent(out) :: ierr
 real,             intent(in), optional :: rval
 integer,          intent(in), optional :: ival
 integer, parameter :: stdout = 6
 integer :: i

 ierr = 0
 write(*,"(15('-'),a,15('-'))") trim(infilenew)
 open(unit=iunit,file=trim(infilenew),form='formatted',status='replace',iostat=ierr)
 if (ierr /= 0) then
    write(*,"(/,a)") 'ERROR! write_infile_series: error opening '//trim(infilenew)//' for write'
    return
 endif
 do i=1,size(infiledata)
    if (i==linenum) then
       if (present(rval)) then
          call write_inopt(rval,trim(nameval),trim(commentval),iunit)
          call write_inopt(rval,trim(nameval),trim(commentval),stdout)
       elseif (present(ival)) then
          call write_inopt(ival,trim(nameval),trim(commentval),iunit)
          call write_inopt(ival,trim(nameval),trim(commentval),stdout)
       else
          write(*,"(/,a)") 'ERROR! write_infile_lines: internal error in subroutine call'
          ierr = 1
          return
       endif
    else
       write(iunit,"(a)") trim(infiledata(i))
    endif
 enddo

 close(unit=iunit)

end subroutine write_infile_lines

!---------------------------------------------------------------------------
!
! Creates a string out of a list of options
!
!---------------------------------------------------------------------------
subroutine get_optstring(nopts,optstring,string,maxlen)
 integer,          intent(in)  :: nopts
 character(len=*), intent(in)  :: optstring(nopts)
 character(len=*), intent(out) :: string
 integer,          intent(in), optional :: maxlen
 character(len=len(string)) :: temp
 integer            :: i,maxl,ierr

 if (present(maxlen)) then
    maxl = max(maxlen,1)
 else
    maxl = 4
 endif

 string = ''
 do i=1,nopts
    temp = optstring(i)
    if (i==nopts) then
       write(string(len_trim(string)+1:),"(i2,'=',a)",iostat=ierr) i,trim(temp(1:maxl))
    elseif (i < 10) then
       write(string(len_trim(string)+1:),"(i1,'=',a,',')",iostat=ierr) i,trim(temp(1:maxl))
    else
       write(string(len_trim(string)+1:),"(i2,'=',a,',')",iostat=ierr) i,trim(temp(1:maxl))
    endif
 enddo

end subroutine get_optstring

end module infile_utils
