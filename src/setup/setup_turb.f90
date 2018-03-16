!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  Sets up a calculation of supersonic turbulence in a periodic box.
!  Works for hydro, mhd, and dusty turbulence.
!
!  REFERENCES:
!    Price & Federrath (2010), MNRAS
!    Tricco, Price & Federrath (2016), MNRAS
!    Tricco, Price & Laibe (2017), MNRAS Letters
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, io, mpiutils, options, part, physcon,
!    prompting, setup_params, timestep, unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:use_dust
 use options,      only:use_dustfrac,nfulldump,beta
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:Bxyz,mhd,dustfrac
 use physcon,      only:pi,solarm,pc,km
 use units,        only:set_units
 use prompting,    only:prompt
 use dust,         only:set_dustfrac,grainsizecgs
 use timestep,     only:dtmax,tmax
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=26)                :: filename
 integer :: ipart,i,maxp,maxvxyzu
 logical :: iexist
 real :: totmass,deltax
 real :: Bz_0, dust_to_gas, grainsize

 print *, ''
 print *, 'Setup for turbulence in a periodic box'
 print *, ''
!
!--boundaries
!
 call set_boundary(0.,1.,0.,1.,0.,1.)
!
!--units
!
 ! Molecular cloud conditions:
 !  L = 3 pc, rho = 1e-20 g/cm^3, c_s = 0.2 km/s
 call set_units(dist=3.0*pc, mass=1e-20*(3.0*pc)**3, time=3.0*pc/(0.2*km))
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))
 if (id==master) then
    npartx = 64
    call prompt('Enter number of particles in x ',npartx,16,nint((maxp)**(1/3.)))
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx

 if (id==master) then
    ilattice = 1
    call prompt('Select lattice type (1=cubic, 2=closepacked)',ilattice,1,2)
 endif

 call bcast_mpi(ilattice)
 if (id==master) then
    rhozero = 1.
    call prompt('Enter density (gives particle mass)',rhozero,0.)
 endif
 call bcast_mpi(rhozero)

 if (maxvxyzu < 4) then
    if (id==master) then
       polykset = 1.
       call prompt('Enter sound speed in code units (sets polyk)',polykset,0.)
    endif
    call bcast_mpi(polykset)
    polyk = polykset**2
    print*,' polyk = ',polyk
 else
    polyk = 0.
    polykset = 0.
 endif

 if (use_dust) use_dustfrac = .true.
 if (use_dustfrac) then
    print *, ''
    print *, 'Setting up dusty turbulence:'
    print *, ''

    dust_to_gas = 1.e-2
    if (id==master) call prompt('Enter dust-to-gas ratio ',dust_to_gas)
    call bcast_mpi(dust_to_gas)

    grainsize = 0.1 ! micron
    if (id==master) call prompt('Enter grain size in micron ',grainsize)
    call bcast_mpi(grainsize)
    grainsizecgs = grainsize * 1.0e-4
 endif

 if (mhd) then
    Bz_0 = 1.4142e-5
    print *, ''
    print *, 'Setting up MHD turbulence: (with uniform intial magnetic field in z-direction)'
    print *, ''
    if (id==master) call prompt('Enter initial magnetic field strength ',Bz_0)
 endif


 ! setup preferred values of .in file
 filename= trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax         = 1.00   ! run for 20 turbulent crossing times
    dtmax        = 0.0025
    nfulldump    = 5      ! output 4 full dumps per crossing time
    beta         = 4      ! legacy from Price & Federrath (2010), haven't checked recently if still required
    grainsizecgs = grainsize * 1.0e-4
 endif
 npart = 0
 npart_total = 0


 select case(ilattice)
 case(1)
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 case(2)
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 case default
    print*,' error: chosen lattice not available, using cubic'
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 end select

 npartoftype(:) = 0
 npartoftype(1) = npart
 print *, ' npart = ',ipart,npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 print *, ' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(3,i) = Bz_0
    endif
    if (use_dustfrac) then
       call set_dustfrac(dust_to_gas,dustfrac(:,i))
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup

