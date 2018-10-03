!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2divv
!
!  DESCRIPTION: This program is a post-processing tool to calculate divv
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom2divv dumpfile(s)
!
!  DEPENDENCIES: deriv, dim, initial, io, kernel, part, readwrite_dumps
!+
!--------------------------------------------------------------------------
program phantom2divv
 use dim,             only:ndivcurlv,maxp,tagline
 use part,            only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
                           hfact,rhoh,dhdrho,igas,isetphase,iphase,maxphase,&
                           dustfrac,ddustevol,temperature,dustprop,ddustprop
 use io,              only:set_io_unit_numbers,iprint,idisk1,idump
 use initial,         only:initialise
 use readwrite_dumps, only:read_dump,write_fulldump
 use deriv,           only:derivs
 use kernel,          only:hfact_default
 implicit none
 integer :: nargs
 character(len=120) :: dumpfile
 real :: time,dtdum
 integer :: ierr,iarg,i

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",trim(tagline)
    print "(a)",' Usage: phantom2divv dumpfile(s)'
    stop
 endif

 print "(/,a,/)",' Phantom2divv: our divergence (and now curl) is your pleasure'

 call initialise()
 if (ndivcurlv < 1) stop 'error: need ndivcurlv=1 for this to do anything'
 if (ndivcurlv < 4) print "(a)",' WARNING: need ndivcurlv=4 in dim file to get curl v as well as div v'

 over_args: do iarg=1,nargs

    call get_command_argument(iarg,dumpfile)
!
!--read particle setup from dumpfile
!
    call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
    if (hfact<epsilon(hfact)) hfact=hfact_default
!
!--calculate derivatives including the divergence of v
!
    if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,&
                temperature,0.,0.,dtdum)
!
!--dump to .divv file
!
    print "(a)",' writing output to file '//trim(dumpfile)//'.divv'
    open(unit=idump,file=trim(dumpfile)//'.divv',form='unformatted',status='replace')
    write(idump) (divcurlv(1,i),i=1,npart)
    if (ndivcurlv >= 4) then
       write(idump) (divcurlv(2,i),i=1,npart)
       write(idump) (divcurlv(3,i),i=1,npart)
       write(idump) (divcurlv(4,i),i=1,npart)
    else
       print*,' skipping curl v (not stored)'
    endif
!    write(idump) (-divv(i)/(dhdrho(xyzh(4,i))*rhoh(xyzh(4,i))),i=1,npart)
    close(idump)

 enddo over_args
 print "(/,a,/)",' Phantom2divv: may your divergence diverge and your curl curl'

end program phantom2divv
