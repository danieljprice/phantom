!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantom2sphNG
!
! This program converts Phantom dumps into sphNG dumps
! (the opposite is not required as Phantom can
!  be started directly from an sphNG dump file).
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantom2sphNG dumpfilein dumpfileout
!
! :Dependencies: dim, io, part, readwrite_dumps
!
 use dim,             only:tagline
 use part,            only:hfact
 use io,              only:set_io_unit_numbers,iprint,idisk1
 use readwrite_dumps, only:read_dump,write_fulldump
 implicit none
 integer :: nargs
 character(len=120) :: dumpfilein,dumpfileout
 real :: time
 integer :: ierr

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs /= 2) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantom2sphNG dumpfilein dumpfileout'
    stop
 endif
 call get_command_argument(1,dumpfilein)
 call get_command_argument(2,dumpfileout)

 print "(/,a,/)",' Phantom2sphNG: The dirtiest conversion in the west...'
!
!--read particle setup from dumpfile
!
 call read_dump(trim(dumpfilein),time,hfact,idisk1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading dumpfile'

 call write_fulldump(time,dumpfileout,sphNG=.true.)

 print "(/,a,/)",' Phantom2sphNG: Good luck with that.'

end program phantom2sphNG

