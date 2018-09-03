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
!  from the dumpfile in the first argument on the command line, and the output
!  dumpfile is the last argument. We assume that all dumps have the same number
!  of dust particles.
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  USAGE: combinedustdumps inputdumpfile1 ... inputdumpfileN outputdumpfile
!
!  DEPENDENCIES: deriv, dim, initial, io, part, readwrite_dumps
!+
!--------------------------------------------------------------------------
program combinedustdumps
 use deriv,           only:derivs
 use dim,             only:maxp,tagline
 use initial,         only:initialise
 use io,              only:set_io_unit_numbers,iprint,idisk1
 use part,            only:xyzh,vxyzu,npart,hfact,iphase,npartoftype,massoftype,&
                           igas,idust,&
                           ndusttypes,ndustsmall,ndustlarge,fxyzu,fext,divcurlv,&
                           divcurlB,Bevol,dBevol,dustfrac,ddustevol,temperature,&
                           dustprop,ddustprop,set_particle_type
 use readwrite_dumps, only:read_dump,write_fulldump
 implicit none
 character(len=120), allocatable :: indumpfiles(:)
 character(len=120) :: outdumpfile
 real, allocatable  :: xyzh_tmp(:,:,:),vxyzu_tmp(:,:,:),massoftype_tmp(:),npartoftype_tmp(:)
 integer :: i,j,counter,ipart,itype,ierr,nargs,idust_tmp,ninpdumps
 real    :: time,dtdum

 call set_io_unit_numbers
 iprint = 6

 !
 !--get name of run from the command line
 !
 nargs = command_argument_count()
 ninpdumps = nargs -1
 if (nargs < 3) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: combinedustdumps inputdumpfiles ... outputdumpfile'
    stop
 endif
 allocate(indumpfiles(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store filenames'

 do i=1,ninpdumps
    call get_command_argument(i,indumpfiles(i))
 enddo
 call get_command_argument(nargs,outdumpfile)

 print "(/,a,/)",' combinedustdumps: many grains make light work'

 !
 !--allocate arrays for dust particle positions and dustfrac
 !
 call read_dump(trim(indumpfiles(1)),time,hfact,idisk1,iprint,0,1,ierr)

 !--assume all dumps are from the same phantom version
 if (npartoftype(2) > 0) then
    !--old dumps
    idust_tmp = 2
 else
    !--new dumps
    idust_tmp = 7
 endif

 allocate (xyzh_tmp(ninpdumps,4,npartoftype(idust_tmp)),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store positions'
 allocate (vxyzu_tmp(ninpdumps,4,npartoftype(idust_tmp)),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store velocities'
 allocate (npartoftype_tmp(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store particle mass'
 allocate (massoftype_tmp(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store particle mass'

 !
 !--get dust particle positions and velocities
 !
 do i=1,ninpdumps
    call read_dump(trim(indumpfiles(i)),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
    npartoftype_tmp(i) = npartoftype(idust_tmp)
    massoftype_tmp(i) = massoftype(idust_tmp)
    counter = 0
    do j=1,maxp
       if (iphase(j)==idust_tmp) then
          counter = counter + 1
          xyzh_tmp(i,:,counter) = xyzh(:,j)
          vxyzu_tmp(i,:,counter) = vxyzu(:,j)
       endif
    enddo
    if (counter /= npartoftype(idust_tmp)) stop 'wrong number of dust particles'
 enddo

 !
 !--re-read first dump and set particle properties
 !
 call read_dump(trim(indumpfiles(1)),time,hfact,idisk1,iprint,0,1,ierr)
 do i=1,maxp
    if (iphase(i)==idust_tmp) then
       call set_particle_type(i,idust)
    endif
 enddo
 npartoftype(idust) = npartoftype_tmp(1)
 massoftype(idust) = massoftype_tmp(1)

 !
 !--add dust from other dumps
 !
 do i=2,ninpdumps
    itype = idust + i - 1
    npartoftype(itype) = npartoftype_tmp(i)
    massoftype(itype) = massoftype_tmp(i)
    do j=1,npartoftype_tmp(i)
       ipart = npart + j
       call set_particle_type(ipart,itype)
       xyzh(:,ipart) = xyzh_tmp(i,:,j)
       vxyzu(:,ipart) = vxyzu_tmp(i,:,j)
    enddo
    npart = npart + npartoftype_tmp(i)
 enddo

 ! units?

 !
 !--dust properties
 !
 ! grainsize?
 ! graindens?
 ndusttypes = ninpdumps
 ndustlarge = ndusttypes
 ndustsmall = 0

 !
 !--calculate dustfrac
 !
 call initialise()
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtdum)

 !
 !--write multigrain dump
 !
 call write_fulldump(time,outdumpfile)

 if (allocated(xyzh_tmp)) deallocate(xyzh_tmp)
 if (allocated(vxyzu_tmp)) deallocate(vxyzu_tmp)
 if (allocated(npartoftype_tmp)) deallocate(npartoftype_tmp)
 if (allocated(massoftype_tmp)) deallocate(massoftype_tmp)

end program combinedustdumps


