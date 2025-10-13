!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantom2gadget
!
! This program converts Phantom dumps into GADGET dumps
! (SPLASH can do the opposite)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantom2gadget dumpfilein dumpfileout
!
! :Dependencies: boundary, dim, eos, io, part, readwrite_dumps, units
!
 use dim,             only:tagline
 use part,            only:hfact,massoftype,npart,xyzh,vxyzu,rhoh
 use eos,             only:polyk
 use io,              only:set_io_unit_numbers,iprint,idisk1
 use readwrite_dumps, only:read_dump,write_gadgetdump
 use units,           only:set_units
 use boundary,        only:set_boundary
 implicit none
 integer :: nargs
 character(len=120) :: dumpfilein,dumpfileout
 real :: time
 real, allocatable :: rho(:)
 integer :: ierr,i

 call set_io_unit_numbers
 call set_units
 call set_boundary
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs /= 2) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantom2gadget dumpfilein dumpfileout'
    stop
 endif
 call get_command_argument(1,dumpfilein)
 call get_command_argument(2,dumpfileout)

 print "(/,a,/)",' Phantom2gadget: Why?'
!
!--read particle setup from dumpfile
!
 call read_dump(trim(dumpfilein),time,hfact,idisk1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading dumpfile'

 allocate(rho(npart))
 do i=1,npart
    rho(i) = rhoh(xyzh(4,i),massoftype(1))
 enddo

 call write_gadgetdump(trim(dumpfileout),time,xyzh,massoftype(1),vxyzu,rho,1.5*polyk,npart)
 deallocate(rho)

 print "(/,a,/)",' Phantom2gadget: Good luck with that.'

end program phantom2gadget
