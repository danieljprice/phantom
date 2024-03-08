!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program showarrays
!
! Utility to print tags of all array in a dump file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: showarrays file_00000
!
! :Dependencies: dump_utils
!
 use dump_utils, only:print_arrays_in_file
 implicit none
 integer, parameter   :: iu = 23
 integer              :: nargs,i
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
