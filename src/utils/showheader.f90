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
 use dump_utils, only:open_dumpfile_r,read_header,dump_h,lenid,&
                      print_header,free_header,ierr_realsize
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
       if (ierr == ierr_realsize) then
          ! try to open it as a small dump
          close(iu)
          call open_dumpfile_r(iu,dumpfile,fileid,ierr,singleprec=.true.)
          if (ierr == 0) call read_header(iu,hdr,.true.,ierr,singleprec=.true.)
       else
          print "(a)",' ERROR opening '//trim(dumpfile)
       endif
    else
       !
       ! read and print the file header
       !
       call read_header(iu,hdr,.true.,ierr)
    endif
    if (ierr == 0) then
       if (nargs > 1) print "(/,':: ',a,' ::',/)",trim(dumpfile)
       print "(a)",trim(fileid)
       call print_header(hdr)
       call free_header(hdr,ierr)
       close(iu)
    endif
 enddo

end program showheader
