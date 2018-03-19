!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: get_density_tracer_particles
!
!  DESCRIPTION: This program is a utility for calculating an SPH density
!   on tracer particles from the FLASH code using PHANTOM routines
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: get_density_tracer_particles dumpfilename(s)
!
!  DEPENDENCIES: deriv, dim, eos, flash_hdf5read, initial, io, mpiutils,
!    options, part, readwrite_dumps, setup, timing
!+
!--------------------------------------------------------------------------
program get_density_tracer_particles
 use dim, only:tagline
 use io, only:id,nprocs,iprint,master,set_io_unit_numbers, iverbose
 use initial, only:initialise
 use mpiutils, only:init_mpi
 use options, only:tolh
 use part, only:xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol, &
           massoftype,hfact,npart,npartoftype,rhoh,dustfrac,ddustfrac
 use eos, only:gamma,polyk
 use setup, only:setpart,dumpfile
 use timing, only:getused,printused
 use deriv, only:derivs
 use readwrite_dumps, only:write_smalldump
 use flash_hdf5read, only:write_tracer_particle_density
 use, intrinsic :: iso_c_binding, only:c_double
 implicit none
 integer :: iarg,nargs,i,ierr
 real :: time,dtdum
 real(kind=c_double), allocatable :: rho(:)
 real(kind=4) :: t1
!
!--pre-initialisation
!
 call init_mpi(id,nprocs)
 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: get_density_tracer_particles dumpfilename(s)'
    stop
 endif

 print "(/,a,/)",' get_density_tracer_particles: SPH is great you know.'
 print "(/,a,/)",' get_density_tracer_particles: WARNING! May not work with magnetic fields after 2 Mar 2018!'
!
!--initialise default code options, units etc.
!
 call initialise
 !tolh = huge(1.)
 iverbose = 2
!
!--loop over all filenames on command line
!
 over_args: do iarg=1,nargs

    call get_command_argument(iarg,dumpfile)
!
!--"setup" tracer particles from hdf5 dump file
!
    call setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time)

    print*,'initial density = ',(rhoh(xyzh(4,i),massoftype(i)),i=1,4)
    print*,' tolerance on iterations = ',tolh
!
!--calculate density on the tracer particles
!
    call getused(t1)
    call derivs(1,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
                Bevol,dBevol,dustfrac,ddustfrac,time,0.,dtdum)
    if (id==master) call printused(t1)

!
!--write tracer particle density to file
!
    !call write_smalldump(time,trim(dumpfile)//'.sphNG')
    allocate(rho(npart))
    do i=1,npart
       rho(i) = rhoh(xyzh(4,i),massoftype(1))
    enddo
    print*,'final density  = ',rho(1:4)
    call write_tracer_particle_density(trim(dumpfile)//achar(0),npart,rho(:),ierr)
    if (allocated(rho)) deallocate(rho)

 enddo over_args

 print "(/,a,/)",' get_density_tracer_particles: thank you for using Phantom. Please come again some time.'

end program get_density_tracer_particles

