!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: io
!
!  DESCRIPTION:
!   This module contains utility routines related to input/output
!   of runtime information from the code
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: mpi
!+
!--------------------------------------------------------------------------
module io
 implicit none
 integer, parameter, public :: stdout = 6
 integer, parameter, public :: stderr = 0

 !--MPI stuff but must be set in serial too
 integer, parameter, public :: master = 0
 integer, public :: id,nprocs
 !--the rest
 integer, public :: iprint, ievfile, idump, ireadin, iwritein, idisk1
 integer, public :: imflow, ivmflow, ibinpos, igpos
 integer, public :: ifile,ifdump,ifdumpread,ireadgrid,ireaddrv,ianalysis
 integer, public :: iscfile,ipafile,iskfile,igit,iuniteos

 !--verboseness level is set to zero by default
 !  but can be changed by external routines (e.g. as input file option)
 integer, public :: iverbose = 0

 public :: set_io_unit_numbers,real4,warn,error,fatal,die,comment
 public :: formatreal,formatint,cstring,warning

 !--generic interface to convert either real*8's or integers
 !  to real*4 for output
 interface real4
  module procedure realr4,reald4,reali4
 end interface

 !--generic interface of format int handles both int*4 and int*8
 interface formatint
  module procedure formatint4,formatint8
 end interface
 !
 !--internal storage of warnings database for buffered warnings
 !
 public :: flush_warnings
 logical, private :: buffer_warnings = .true.
 logical, private :: warningdb_initialised = .false.

 integer, parameter, private :: lenmsg = 120
 integer, parameter, private :: lenwhere = 20

 type warningdb_entry
    character(len=lenwhere) :: wherefrom
    character(len=lenmsg)   :: message
    integer(kind=8)         :: ncount
    integer :: level
 end type

 integer, parameter :: maxwarningdb = 20
 type(warningdb_entry) :: warningdb(maxwarningdb)
 integer, parameter :: maxcount = 10

 private

contains
!--------------------------------------------------------------------
!+
!  Sets the logical unit numbers to use for output
!+
!--------------------------------------------------------------------
subroutine set_io_unit_numbers

#ifdef MPI
 iprint = 6     ! only iprint=6 makes sense for MPI runs
#else
 iprint = 6 !8     ! for writing log output
 nprocs = 1
 id     = 0
#endif
 imflow  = 47
 ivmflow = 48
 ibinpos = 49
 igpos   = 46

 ievfile    = 13 ! for writing energies etc
 idump      = 12 ! for writing dump files
 ireadin    = 21 ! for reading input file
 iwritein   = 22 ! for writing input file
 idisk1     =  9 ! for reading dump files
 ifdump     = 51 ! for writing forcing restart file
 ifdumpread = 18 ! for reading forcing restart file
 ireadgrid  = 53 ! for reading gridded density derivative file
 ireaddrv   = 24 ! for reading input file for turbulent driving
 ianalysis  = 25 ! for writing analysis output
 iuniteos   = 28 ! for printing the eos to file
 igit       = 29 ! for reading phantom_version
 iscfile    = 32 ! for writing details of sink creation
 ipafile    = 31 ! for writing details of particle accretion
 iskfile    =407 ! for writing details of the sink particles;
 ! note this actually opens files iskfile+1 to iskfile+nptmass
 iverbose   = 0

end subroutine set_io_unit_numbers


!--------------------------------------------------------------------
!+
!  small function to convert from real8 to real4
!+
!-------------------------------------------------------------------
elemental real(kind=4) function reald4(x)
 real(kind=8), intent(in) :: x

 reald4 = real(x,kind=4)

end function reald4

!--------------------------------------------------------------------
!+
!  small function to convert from real8 to real4
!+
!-------------------------------------------------------------------
elemental real(kind=4) function realr4(x)
 real(kind=4), intent(in) :: x

 realr4 = x

end function realr4

!--------------------------------------------------------------------
!+
!  small function to convert from int to real4
!+
!-------------------------------------------------------------------
elemental real(kind=4) function reali4(ix)
 integer, intent(in) :: ix

 reali4 = ix

end function reali4

!--------------------------------------------------------------------
!+
!  routine that copies a formats a real variable to a string
!  in the neatest possible way
!+
!--------------------------------------------------------------------
subroutine formatreal(val,string,ierror,precision)
 real,             intent(in)  :: val
 character(len=*), intent(out) :: string
 integer,          intent(out), optional :: ierror
 real,             intent(in),  optional :: precision
 integer :: ierr,i !,idot
 !logical :: nonzero
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
 else !if (abs(val) >= 100.) then
    write(string,"(f8.0)",iostat=ierr) val
    !else
    !   write(string,"(f8.2)",iostat=ierr) val
 endif
 string = adjustl(trim(string))

 if (present(ierror)) ierror = ierr

 !
 !--strip trailing zeros after the decimal place
 !  (and the decimal place if it is the last character)
 !
 i = len_trim(string)
 if (string(i:i)=='.') string(i:i) = ' '

! idot = index(string,'.')
! if (idot > 0) then
!    nonzero = .false.
!    do i = len_trim(string),idot,-1
!       if (.not.nonzero .and. string(i:i)=='0') then
!          string(i:i) = ' '
!       elseif (.not.nonzero .and. string(i:i)=='.') then
!          string(i:i) = ' '
!          nonzero = .true.
!       else
!          nonzero = .true.
!       endif
!    enddo
! endif
 string = trim(adjustl(string))

 return
end subroutine formatreal

!--------------------------------------------------------------------
!+
!  routine that copies a formats an integer variable to a string
!  in the neatest possible way
!+
!--------------------------------------------------------------------
subroutine formatint4(ival,string)
 integer(kind=4),  intent(in)  :: ival
 character(len=*), intent(out) :: string

 write(string,"(i12)") ival
 string = trim(adjustl(string))

end subroutine formatint4

subroutine formatint8(ival,string)
 integer(kind=8),  intent(in)  :: ival
 character(len=*), intent(out) :: string

 write(string,"(i12)") ival
 string = trim(adjustl(string))

end subroutine formatint8

!---------------------------------------------------------------------------
!
! function to safely convert a string to c format (ie. with a terminating
! ascii null character)
!
!---------------------------------------------------------------------------
function cstring(string)
 character(len=*), intent(in) :: string
 character(len=len(string)+1) :: cstring

 cstring = trim(string)//achar(0)

end function cstring

!--------------------------------------------------------------------
!+
!  subroutine to initialise warnings database
!+
!--------------------------------------------------------------------
subroutine init_warningdb()
 integer :: j

 do j=1,maxwarningdb
    warningdb(j)%wherefrom = ''
    warningdb(j)%message   = ''
    warningdb(j)%ncount    = 0_8
    warningdb(j)%level     = 0
 enddo
 warningdb_initialised = .true.

end subroutine init_warningdb

!--------------------------------------------------------------------
!+
!  subroutine to buffer warnings in a database
!+
!--------------------------------------------------------------------
subroutine buffer_warning(wherefrom,string,ncount,level)
 character(len=*), intent(in)  :: wherefrom, string
 integer(kind=8),  intent(out) :: ncount
 integer,          intent(in), optional :: level
 integer :: lw,ls,j

 if (.not.warningdb_initialised) call init_warningdb
 ncount = 1
 lw = min(len_trim(wherefrom),lenwhere)
 ls = min(len_trim(string),lenmsg)
 j = 1
 over_db: do while(j <= maxwarningdb)
    if (warningdb(j)%wherefrom(1:lw) == wherefrom(1:lw) &
        .and. warningdb(j)%message(1:ls) == string(1:ls)) then
       !--if warning matches an existing warning in the database
       !  just increase the reference count
!$omp critical
       warningdb(j)%ncount = warningdb(j)%ncount + 1_8
!$omp end critical
       ncount = warningdb(j)%ncount
       exit over_db
    elseif (len_trim(warningdb(j)%message)==0) then
       !--warning not recognised, make new entry in database
       warningdb(j)%wherefrom = wherefrom
       warningdb(j)%message = string
       warningdb(j)%ncount = 1_8
       if (present(level)) then
          warningdb(j)%level = level
       else
          warningdb(j)%level = 2
       endif
       ncount = 1
       exit over_db
    endif
    j = j + 1
 enddo over_db

end subroutine buffer_warning

!--------------------------------------------------------------------
!+
!  flush warnings database, print occurrences of each warning
!+
!--------------------------------------------------------------------
subroutine flush_warnings(all)
 logical, intent(in), optional :: all
 integer :: j
 integer(kind=8) :: ncount
 logical :: print_warnings, flush_all

 if (.not.warningdb_initialised) call init_warningdb

 ! option to flush all warnings, even if count < maxcount
 if (present(all)) then
    flush_all = all
 else
    flush_all = .false.
 endif

 print_warnings = .false.
 do j=1,maxwarningdb
    if (len_trim(warningdb(j)%message) > 0) print_warnings = .true.
 enddo

 if (print_warnings) then
    do j=1,maxwarningdb
       if (len_trim(warningdb(j)%message) > 0) then
          ncount = warningdb(j)%ncount
          if (ncount > maxcount .or. flush_all) then
             if (warningdb(j)%level > 2) then
                write(iprint,"(' ERROR! ',a,': ',a)") &
                      trim(warningdb(j)%wherefrom),trim(warningdb(j)%message)
             else
                write(iprint,"(' WARNING! ',a,': ',a)") &
                      trim(warningdb(j)%wherefrom),trim(warningdb(j)%message)
             endif
             if (ncount < 10000) then
                if (ncount > 1) then
                   write(iprint,"(' (repeated ',i4,' times)')") ncount
                endif
             else
                write(iprint,"(' (repeated ',i12,' times)')") ncount
             endif
          endif
       endif
    enddo
 endif
 !--reset warnings database
 call init_warningdb

end subroutine flush_warnings

!--------------------------------------------------------------------
!+
!  subroutine giving warnings of various severity
!+
!--------------------------------------------------------------------
subroutine warn(wherefrom,string,severity)
 character(len=*), intent(in) :: string,wherefrom
 integer,          intent(in), optional :: severity
 integer :: iseverity
 integer(kind=8) :: ncount

 if (present(severity)) then
    iseverity = severity
 else
    iseverity = 2
 endif

 select case(iseverity)
 case(1)
    if (iverbose >= 1) then
       write(iprint,"(/' COMMENT: ',a,': ',a)") trim(wherefrom),trim(string)
    endif
 case(2)
    if (buffer_warnings) then
       call buffer_warning(trim(wherefrom),trim(string),ncount)
       if (ncount < maxcount) then
          write(iprint,"(' WARNING! ',a,': ',a)") trim(wherefrom),trim(string)
       elseif (ncount==maxcount) then
          write(iprint,"(' (buffering remaining warnings... ',a,') ')") trim(string)
       endif
    else
       write(iprint,"(' WARNING! ',a,': ',a)") trim(wherefrom),trim(string)
    endif
 case(3)
    ncount = 0_8
    if (buffer_warnings) call buffer_warning(trim(wherefrom),trim(string),ncount,level=3)
    if (buffer_warnings .and. ncount==maxcount) then
       write(iprint,"(' (buffering remaining errors... ',a,') ')") trim(string)
    elseif (ncount < maxcount .or. .not.buffer_warnings) then
       write(iprint,"(/' ERROR! ',a,': ',a,/)") trim(wherefrom),trim(string)
       if (iprint /= 6) write(*,"(/' ERROR! ',a,': ',a,/)") trim(wherefrom),trim(string)
    endif
 case(4)
#ifdef MPI
    write(iprint,"(/' FATAL ERROR! on thread ',i4,' in ',a,': ',a)") id,trim(wherefrom),trim(string)
#else
    write(iprint,"(/' FATAL ERROR! ',a,': ',a)") trim(wherefrom),trim(string)
#endif
    if (iprint /= 6) write(*,"(/' FATAL ERROR! ',a,': ',a)") trim(wherefrom),trim(string)
    call die
 case default
    write(iprint,"(/' WARNING(unknown severity)! ',a,': ',a)") trim(wherefrom),trim(string)
 end select

 return
end subroutine warn

!--------------------------------------------------------------------
!+
!  interface to warn for comments
!+
!--------------------------------------------------------------------
subroutine comment(wherefrom,string)
 character(len=*), intent(in) :: string,wherefrom

 call warn(wherefrom,string,1)

end subroutine comment

!--------------------------------------------------------------------
!+
!  interface to warn for warnings
!+
!--------------------------------------------------------------------
subroutine warning(wherefrom,string,i,var,val,ival)
 character(len=*), intent(in) :: string,wherefrom
 integer,          intent(in), optional :: i
 character(len=*), intent(in), optional :: var
 real,             intent(in), optional :: val
 integer,          intent(in), optional :: ival
 character(len=50) :: newstring
 character(len=12) :: stringi

 newstring = ' '
 if (present(i)) then
    write(newstring,"(i10)") i
    newstring = ' on particle '//trim(adjustl(newstring))
 endif
 if (present(var)) then
    write(newstring,"(1x,a,': ',a)") trim(adjustl(newstring)),var
 endif
 if (present(val)) then
    write(newstring,"(1x,a,' = ',es10.3)") trim(adjustl(newstring)),val
 endif
 if (present(ival)) then
    call formatint(ival,stringi)
    write(newstring,"(1x,a,' = ',a)") trim(adjustl(newstring)),trim(stringi)
 endif

 call warn(wherefrom,string//trim(newstring),2)

end subroutine warning

!--------------------------------------------------------------------
!+
!  interface to warn for non-fatal errors
!+
!--------------------------------------------------------------------
subroutine error(wherefrom,string,i,var,val,ival)
 character(len=*), intent(in) :: string,wherefrom
 integer,          intent(in), optional :: i
 character(len=*), intent(in), optional :: var
 real,             intent(in), optional :: val
 integer,          intent(in), optional :: ival
 character(len=50) :: newstring
 character(len=12) :: stringi

 newstring = ' '
 if (present(i)) then
    write(newstring,"(i10)") i
    newstring = ' on particle '//trim(adjustl(newstring))
 endif
 if (present(var)) then
    write(newstring,"(1x,a,': ',a)") trim(adjustl(newstring)),var
 endif
 if (present(val)) then
    write(newstring,"(1x,a,' = ',es10.3)") trim(adjustl(newstring)),val
 endif
 if (present(ival)) then
    call formatint(ival,stringi)
    write(newstring,"(1x,a,' = ',a)") trim(adjustl(newstring)),trim(stringi)
 endif

 call warn(wherefrom,string//trim(newstring),3)

end subroutine error

!--------------------------------------------------------------------
!+
!  interface to fatal errors
!+
!--------------------------------------------------------------------
subroutine fatal(wherefrom,string,i,var,val,ival)
 character(len=*), intent(in) :: string,wherefrom
 integer,          intent(in), optional :: i
 character(len=*), intent(in), optional :: var
 real,             intent(in), optional :: val
 integer,          intent(in), optional :: ival
 character(len=50) :: newstring
 character(len=12) :: stringi

 newstring = ' '
 if (present(i)) then
    write(newstring,"(i10)") i
    newstring = ' on particle '//trim(adjustl(newstring))
 endif
 if (present(var)) then
    write(newstring,"(1x,a,': ',a)") trim(adjustl(newstring)),var
 endif
 if (present(val)) then
    write(newstring,"(1x,a,' = ',es10.3)") trim(adjustl(newstring)),val
 endif
 if (present(ival)) then
    call formatint(ival,stringi)
    write(newstring,"(1x,a,' = ',a)") trim(adjustl(newstring)),trim(stringi)
 endif
 call warn(wherefrom,string//trim(newstring),4)

end subroutine fatal

!--------------------------------------------------------------------
!+
!  MPI shutdown (replaces stop)
!+
!--------------------------------------------------------------------
subroutine die
#ifdef MPI
 use mpi
 integer :: ierr

 call mpi_abort(mpi_comm_world,666,ierr)
 call mpi_finalize(ierr)
#endif

 call exit(1)
end subroutine die

end module io
