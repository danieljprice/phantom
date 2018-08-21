!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: combinedustdumps
!
!  DESCRIPTION: This program is a utility for stacking multiple 2-fluid dust
!  dumps onto a single set of gas particles. The gas positions will be taken
!  from the last input dumpfile from the command line before the output dumpfile
!  name. We assume that all dumps have the same number of dust particles.
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  USAGE: combinedustdumps inputdumpfile1 ... inputdumpfileN outputdumpfile
!
!  DEPENDENCIES: dim, io, part, readwrite_dumps, testutils
!+
!--------------------------------------------------------------------------
program combinedustdumps
 use deriv,           only:derivs
 use dim,             only:maxp,tagline
 use initial,         only:initialise
 use io,              only:set_io_unit_numbers,iprint,idisk1
 use part,            only:xyzh,vxyzu,npart,hfact,iphase,npartoftype,&
                           idust,ndusttypes,fxyzu,fext,divcurlv,divcurlB,&
                           Bevol,dBevol,dustfrac,ddustfrac,temperature,dustprop,&
                           ddustprop
 use readwrite_dumps, only:read_dump,write_fulldump
 implicit none
 character(len=120), allocatable :: indumpfiles(:)
 character(len=120) :: outdumpfile
 real, allocatable  :: xyzh_tmp(:,:,:),dustfrac_tmp(:,:)
 integer :: i,j,counter,ierr,nargs
 real    :: time,dtdum

 ndusttypes = 1

 call set_io_unit_numbers
 iprint = 6
 !
 !--get name of run from the command line
 !
 nargs = command_argument_count()
 if (nargs < 3) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: combinedustdumps inputdumpfiles ... outputdumpfile'
    stop
 endif
 allocate(indumpfiles(nargs-1),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store filenames'

 do i=1,nargs-1
    call get_command_argument(i,indumpfiles(i))
 enddo
 call get_command_argument(nargs,outdumpfile)

 print "(/,a,/)",' combinedustdumps: we welcome you'

 !
 !--read particle setup from input dumpfiles
 !
 call read_dump(trim(indumpfiles(1)),time,hfact,idisk1,iprint,0,1,ierr)
 allocate (xyzh_tmp(nargs-1,4,npartoftype(idust)),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store positions'
 allocate (dustfrac_tmp(nargs-1,maxp),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store dustfrac'

 do i=1,nargs-1
    call read_dump(trim(indumpfiles(i)),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
    counter = 0
    do j=1,maxp
       if (iphase(j)==idust) then
          counter = counter + 1
          xyzh_tmp(i,:,counter) = xyzh(:,j)
       endif
    enddo
    if (counter /= npartoftype(idust)) stop 'wrong number of dust particles'
 enddo

 do i=1,nargs-1
    counter = 0
    do j=1,maxp
       if (iphase(j)==idust) then
          counter = counter + 1
          xyzh(:,j) = xyzh_tmp(i,:,counter)
       endif
    enddo
    if (counter /= npartoftype(idust)) stop 'wrong number of dust particles'
    call initialise()
    call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustfrac,temperature,time,0.,dtdum)
    print*, '----------------------------------------'
    print*, dustfrac
    print*, '----------------------------------------'
    dustfrac_tmp(i,:) = dustfrac(1,:)
 enddo

 ndusttypes = nargs-1
 dustfrac = dustfrac_tmp
 call write_fulldump(time,outdumpfile)

 if (allocated(xyzh_tmp)) deallocate(xyzh_tmp)

end program combinedustdumps


