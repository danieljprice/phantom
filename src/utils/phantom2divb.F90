!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2divb
!
!  DESCRIPTION: This program is a post-processing tool to calculate divB
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom2divB dumpfile(s)
!
!  DEPENDENCIES: deriv, dim, initial, io, part, readwrite_dumps
!+
!--------------------------------------------------------------------------
program phantom2divb
 use dim,             only:ndivcurlB,maxp,tagline
 use part,            only:npart,xyzh,vxyzu,fxyzu,Bxyz,fext,divcurlv,divcurlB,Bevol,dBevol, &
                           hfact,rhoh,dhdrho,igas,isetphase,iphase,massoftype,maxphase,&
                           dustfrac,ddustevol,mhd,temperature,dustprop,ddustprop
 use io,              only:set_io_unit_numbers,iprint,idisk1,idump
 use initial,         only:initialise
 use readwrite_dumps, only:read_dump,write_fulldump
 use deriv,           only:derivs
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
    print "(a)",' Usage: phantom2divB dumpfile(s)'
    stop
 endif

 print "(/,a,/)",' Phantom2divb: our divergence (and now curl) is your pleasure'

 call initialise()
 if (ndivcurlB < 1) stop 'error: need ndivcurlB=1 for this to do anything'
 if (ndivcurlB < 4) print "(a)",' WARNING: need ndivcurlB=4 in dim file to get curlB as well as divB'

 over_args: do iarg=1,nargs

    call get_command_argument(iarg,dumpfile)
!
!--read particle setup from dumpfile
!
    call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
!
!--calculate derivatives including the divergence of B
!
    Bevol = 0.
    if (mhd) then
       do i = 1,npart
          Bevol(1:3,i) = Bxyz(1:3,i)/rhoh(xyzh(4,i), massoftype(igas))
       enddo
    endif
    if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,&
                temperature,0.,0.,dtdum)
!
!--dump to .divv file
!
    print "(a)",' writing output to file '//trim(dumpfile)//'.divb'
    open(unit=idump,file=trim(dumpfile)//'.divb',form='unformatted',status='replace')
    write(idump) (divcurlB(1,i),i=1,npart)
    if (ndivcurlB >= 4) then
       write(idump) (divcurlB(2,i),i=1,npart)
       write(idump) (divcurlB(3,i),i=1,npart)
       write(idump) (divcurlB(4,i),i=1,npart)
    else
       print*,' skipping curlB (not stored)'
    endif
!    write(idump) (-divv(i)/(dhdrho(xyzh(4,i))*rhoh(xyzh(4,i))),i=1,npart)
    close(idump)

 enddo over_args
 print "(/,a,/)",' Phantom2divB: may your divergence diverge and your curl curl'

end program phantom2divb
