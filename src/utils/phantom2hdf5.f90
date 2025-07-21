!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantom2hdf5
!
! Utility to convert Phantom dump files to HDF5 format
!
! :References: None
!
! :Owner: David Liptai
!
! :Usage: phantom2hdf5 dumpfile(s)
!
! :Dependencies: phantom2hdf5_utils
!
 use phantom2hdf5_utils, only:convert_dump_to_hdf5
 implicit none
 integer :: nargs,iarg,ierr
 character(len=120) :: dumpfile
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",' Usage: phantom2hdf5 dumpfile(s)'
    stop
 endif

 print "(/,a,/)",' Phantom2hdf5: The best conversion in the west...'

 over_args: do iarg=1,nargs
    call get_command_argument(iarg,dumpfile)
    !
    !--read particle setup from dumpfile
    !
    call convert_dump_to_hdf5(dumpfile,ierr)

    ! If there is still an error, skip to the next file
    if (ierr /= 0) then
       print*,'error reading dumpfile: ',trim(dumpfile)
       print*,'skipping to next one...'
       cycle over_args
    endif

 enddo over_args

 print "(/,a,/)",' Phantom2hdf5: Enjoy a better file format.'

end program phantom2hdf5
