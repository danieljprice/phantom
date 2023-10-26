!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program combinedustdumps
!
! This program is a utility for stacking multiple 2-fluid dust
!  dumps onto a single set of gas particles. The gas positions will be taken
!  from the dumpfile named in the first argument on the command line, and the
!  output dumpfile name is the last argument. We assume that all dumps have the
!  same number of dust particles.
!
! :References: None
!
! :Owner: Daniel Mentiplay
!
! :Usage: combinedustdumps inputdumpfiles ... outputdumpfile
!
! :Dependencies: checksetup, deriv, dim, initial, io, memory, part,
!   readwrite_dumps, units
!
 use deriv,           only:get_derivs_global
 use dim,             only:maxp,maxvxyzu,tagline
 use initial,         only:initialise
 use io,              only:set_io_unit_numbers,iprint,idisk1,fatal
 use part,            only:xyzh,vxyzu,npart,hfact,iphase,npartoftype,massoftype,&
                           igas,idust,ndusttypes,ndustsmall,ndustlarge,set_particle_type,&
                           grainsize,graindens,iamtype,isdead_or_accreted
 use readwrite_dumps, only:read_dump,write_fulldump,init_readwrite_dumps
 use units,           only:set_units,select_unit,umass,udist,utime
 use memory,          only:allocate_memory
 use checksetup,      only:check_setup
 implicit none
 character(len=120), allocatable :: indumpfiles(:)
 character(len=120) :: outdumpfile
 real, allocatable :: xyzh_tmp(:,:,:),vxyzu_tmp(:,:,:),massofdust_tmp(:)
 real, allocatable :: grainsize_tmp(:),graindens_tmp(:)
 integer, allocatable :: npartofdust_tmp(:)
 integer :: i,j,ipart,itype,ierr,nargs,idust_tmp,ninpdumps
 integer(kind=8) :: counter
 integer :: nwarn,nerror,ndust
 real    :: time
 real(kind=8) :: utime_tmp,udist_tmp,umass_tmp
 logical :: first_gas_only

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
 ! read first dumpfile HEADER ONLY: check idust, check MAXP
 ! we assume all dumps are from the same phantom version
 !

 call init_readwrite_dumps()
 counter = 0
 idust_tmp = idust ! new dumps, location of first dust particle type
 do i=1,ninpdumps
    call read_dump(trim(indumpfiles(i)),time,hfact,idisk1,iprint,0,1,ierr,headeronly=.true.)
    if (ierr /= 0) stop 'error reading dumpfile... aborting'
    if (i==1) then
       counter = counter + npartoftype(1)
       if (npartoftype(2) > 0) idust_tmp = 2  !old dumps
    endif
    counter = counter + npartoftype(idust_tmp)
 enddo
 !
 ! save the number of dust particles for later
 !
 ndust = npartoftype(idust_tmp)

 !
 !--sanity check array sizes
 !
 if (idust+ninpdumps-1 > size(npartoftype)) then
    call fatal('combinedustdumps','not enough particle types: compile with DUST=yes and',&
               var='MAXDUSTLARGE >',ival=ninpdumps-1)
 endif
 !
 ! allocate memory
 !
 call allocate_memory(counter)
 !
 ! read gas particles from first file
 !

 call read_dump(trim(indumpfiles(1)),time,hfact,idisk1,iprint,0,1,ierr)

 !
 ! allocate temporary arrays
 !

 allocate (xyzh_tmp(ninpdumps,4,ndust),stat=ierr)
 print*,shape(xyzh_tmp),npartoftype(idust_tmp),idust_tmp
 if (ierr /= 0) stop 'error allocating memory to store positions'
 allocate (vxyzu_tmp(ninpdumps,maxvxyzu,ndust),stat=ierr)
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
 first_gas_only = .false.
 do i=1,ninpdumps
    call read_dump(trim(indumpfiles(i)),time,hfact,idisk1,iprint,0,1,ierr)
    if (ierr /= 0) stop 'error reading dumpfile'
    npartofdust_tmp(i) = npartoftype(idust_tmp)
    massofdust_tmp(i) = massoftype(idust_tmp)
    grainsize_tmp(i) = grainsize(1)
    graindens_tmp(i) = graindens(1)
    counter = 0
    if (i == 1 .and. npartoftype(idust_tmp) == 0) first_gas_only = .true.
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
    if (first_gas_only) itype = idust + i - 2
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
 if (first_gas_only) ndusttypes = ninpdumps - 1
 ndustlarge = ndusttypes
 ndustsmall = 0

 !
 !--store units
 !
 utime_tmp = utime
 umass_tmp = umass
 udist_tmp = udist

 !
 !--check the setup is OK
 !
 call check_setup(nerror,nwarn,restart=.true.)
 if (nerror > 0) call fatal('combinedustdumps','errors in merged particle setup')

 !
 !--calculate dustfrac
 !
 call initialise()
 call set_units(udist_tmp,umass_tmp,utime_tmp)
 call get_derivs_global()

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
