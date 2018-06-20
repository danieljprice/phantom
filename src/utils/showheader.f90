!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: showheader
!
!  DESCRIPTION: Utility to print contents of dump file header
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: showheader filename(s)
!
!  DEPENDENCIES: utils_dumpfiles
!+
!--------------------------------------------------------------------------
program showheader
 use dump_utils, only:open_dumpfile_r,read_header,dump_h,lenid,print_header,free_header
 implicit none
 integer, parameter   :: iu = 23
 integer              :: nargs,ierr,i
 character(len=lenid) :: fileid
 character(len=120)   :: dumpfile
 type(dump_h) :: hdr
 !
 ! get filenames from the command line
 !
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",' Usage: showheader file_00000'
    stop
 endif
 do i=1,nargs
    call get_command_argument(i,dumpfile)
    !
    ! open the file
    !
    call open_dumpfile_r(iu,dumpfile,fileid,ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR opening '//trim(dumpfile)
    else
       if (nargs > 1) print "(/,':: ',a,' ::',/)",trim(dumpfile)
       print "(a)",trim(fileid)
       !
       ! read and print the file header
       !
       call read_header(iu,hdr,.true.,ierr)
       call print_header(hdr)
       call free_header(hdr,ierr)
       close(iu)
    endif
 enddo

end program showheader
