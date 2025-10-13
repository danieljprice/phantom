!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantom2divv
!
! This program is a post-processing tool to calculate divv
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantom2divv dumpfile(s)
!
! :Dependencies: deriv, dim, initial, io, kernel, part, readwrite_dumps
!
 use dim,             only:maxp,tagline,curlv
 use part,            only:npart,divcurlv,hfact,igas,isetphase,iphase,maxphase
 use io,              only:set_io_unit_numbers,iprint,idisk1,idump
 use initial,         only:initialise
 use readwrite_dumps, only:read_dump,write_fulldump
 use deriv,           only:get_derivs_global
 use kernel,          only:hfact_default
 implicit none
 integer :: nargs
 character(len=120) :: dumpfile
 real :: time
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
 curlv = .true.
 call initialise()

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
    call get_derivs_global()
!
!--dump to .divv file
!
    print "(a)",' writing output to file '//trim(dumpfile)//'.divv'
    open(unit=idump,file=trim(dumpfile)//'.divv',form='unformatted',status='replace')
    write(idump) (divcurlv(1,i),i=1,npart)
    if (curlv) then
       write(idump) (divcurlv(2,i),i=1,npart)
       write(idump) (divcurlv(3,i),i=1,npart)
       write(idump) (divcurlv(4,i),i=1,npart)
    else
       print*,' skipping curl v (not stored)'
    endif
    close(idump)

 enddo over_args
 print "(/,a,/)",' Phantom2divv: may your divergence diverge and your curl curl'

end program phantom2divv
