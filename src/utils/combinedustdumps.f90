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
!  from the dumpfile named in the first argument on the command line, and the
!  output dumpfile name is the last argument. We assume that all dumps have the
!  same number of dust particles.
!
!  REFERENCES: None
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  USAGE: combinedustdumps inputdumpfiles ... outputdumpfile
!
!  DEPENDENCIES: deriv, dim, initial, io, part, readwrite_dumps, units
!+
!--------------------------------------------------------------------------
program combinedustdumps
 use deriv,           only:derivs
 use dim,             only:maxp,tagline
 use initial,         only:initialise
 use io,              only:set_io_unit_numbers,iprint,idisk1
 use part,            only:xyzh,vxyzu,npart,hfact,iphase,npartoftype,massoftype,&
                           igas,idust,ndusttypes,ndustsmall,ndustlarge,fxyzu,fext,&
                           divcurlv,divcurlB,Bevol,dBevol,dustfrac,ddustevol,&
                           temperature,dustprop,ddustprop,set_particle_type,&
                           grainsize,graindens,iamtype,isdead_or_accreted
 use readwrite_dumps, only:read_dump,write_fulldump
 use units,           only:set_units,select_unit,umass,udist,utime
 implicit none

 character(len=120), allocatable :: indumpfiles(:)
 character(len=120) :: outdumpfile
 real, allocatable :: xyzh_tmp(:,:,:),vxyzu_tmp(:,:,:),massofdust_tmp(:)
 real, allocatable :: grainsize_tmp(:),graindens_tmp(:)
 integer, allocatable :: npartofdust_tmp(:)
 integer :: i,j,counter,ipart,itype,ierr,nargs,idust_tmp,ninpdumps
 real    :: time,dtdum
 real(kind=8) :: utime_tmp,udist_tmp,umass_tmp

 call set_io_unit_numbers
 iprint = 6

 !
 !--get dumpfile names from the command line
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
 !--read first dumpfile: check idust, check MAXP
 !  we assume all dumps are from the same phantom version
 !
 call read_dump(trim(indumpfiles(1)),time,hfact,idisk1,iprint,0,1,ierr)
 if (npartoftype(2) > 0) then
    !--old dumps
    idust_tmp = 2
 else
    !--new dumps
    idust_tmp = idust
 endif
 if (npartoftype(1) + ninpdumps*npartoftype(idust_tmp) > maxp) then
    print "(a)",' MAXP too small, set MAXP >= ngas+ndumps*ndust and recompile'
    stop
 endif

 !
 !--allocate temporary arrays
 !
 allocate (xyzh_tmp(ninpdumps,4,npartoftype(idust_tmp)),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store positions'
 allocate (vxyzu_tmp(ninpdumps,4,npartoftype(idust_tmp)),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store velocities'
 allocate (npartofdust_tmp(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store number of dust particles'
 allocate (massofdust_tmp(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store dust particle mass'
 allocate (grainsize_tmp(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store grain size'
 allocate (graindens_tmp(ninpdumps),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store grain dens'

 !
 !--read dumps and get dust particle information
 !
 do i=1,ninpdumps
    call read_dump(trim(indumpfiles(i)),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
    npartofdust_tmp(i) = npartoftype(idust_tmp)
    massofdust_tmp(i) = massoftype(idust_tmp)
    grainsize_tmp(i) = grainsize(1)
    graindens_tmp(i) = graindens(1)
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
 !--re-read first dump and correct dust particle properties
 !
 call read_dump(trim(indumpfiles(1)),time,hfact,idisk1,iprint,0,1,ierr)
 do i=1,maxp
    if (iamtype(iphase(i))==idust_tmp) then
       call set_particle_type(i,idust)
       if (isdead_or_accreted(xyzh(4,i))) then
          iphase(i) = -iphase(i)
       endif
    endif
 enddo
 npartoftype(idust) = npartofdust_tmp(1)
 massoftype(idust) = massofdust_tmp(1)
 npartoftype(2) = 0
 massoftype(2) = 0
 grainsize(1) = grainsize_tmp(1)
 graindens(1) = graindens_tmp(1)

 !
 !--add dust from other dumps
 !
 do i=2,ninpdumps
    itype = idust + i - 1
    npartoftype(itype) = npartofdust_tmp(i)
    massoftype(itype) = massofdust_tmp(i)
    grainsize(i) = grainsize_tmp(i)
    graindens(i) = graindens_tmp(i)
    do j=1,npartofdust_tmp(i)
       ipart = npart + j
       xyzh(:,ipart) = xyzh_tmp(i,:,j)
       vxyzu(:,ipart) = vxyzu_tmp(i,:,j)
       call set_particle_type(ipart,itype)
       if (isdead_or_accreted(xyzh(4,ipart))) then
          iphase(ipart) = -iphase(ipart)
       endif
    enddo
    npart = npart + npartofdust_tmp(i)
 enddo

 !
 !--dust properties
 !
 ndusttypes = ninpdumps
 ndustlarge = ndusttypes
 ndustsmall = 0

 !
 !--store units
 !
 utime_tmp = utime
 umass_tmp = umass
 udist_tmp = udist

 !
 !--calculate dustfrac
 !
 call initialise()
 call set_units(udist_tmp,umass_tmp,utime_tmp)
 call derivs(1,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
             Bevol,dBevol,dustprop,ddustprop,dustfrac,ddustevol,temperature,time,0.,dtdum)

 !
 !--write multigrain dump
 !
 call write_fulldump(time,outdumpfile)

 if (allocated(xyzh_tmp))        deallocate(xyzh_tmp)
 if (allocated(vxyzu_tmp))       deallocate(vxyzu_tmp)
 if (allocated(npartofdust_tmp)) deallocate(npartofdust_tmp)
 if (allocated(massofdust_tmp))  deallocate(massofdust_tmp)
 if (allocated(grainsize_tmp))   deallocate(grainsize_tmp)
 if (allocated(graindens_tmp))   deallocate(graindens_tmp)

end program combinedustdumps
