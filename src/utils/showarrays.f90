!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: showarrays
!
!  DESCRIPTION: Utility to print tags of all array in a dump file
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: showarrays file_00000
!
!  DEPENDENCIES: dump_utils
!+
!--------------------------------------------------------------------------
program showarrays
 use dump_utils, only:print_arrays_in_file
 implicit none
 integer, parameter   :: iu = 23
 integer              :: nargs,ierr,i
 character(len=120)   :: dumpfile
 !
 ! get filenames from the command line
 !
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",' Usage: showarrays file_00000'
    stop
 endif
 do i=1,nargs
    call get_command_argument(i,dumpfile)
    !
    ! open the file
    !
    call print_arrays_in_file(iu,dumpfile)
 enddo

end program showarrays
