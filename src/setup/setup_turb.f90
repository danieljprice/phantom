!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  Sets up a calculation of supersonic turbulence in a periodic box.
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, io, mpiutils, part, physcon,
!    prompting, setup_params, unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: deltax,rhozero,polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:use_dustfrac
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:Bevol,maxvecp,mhd,maxBevol,dustfrac
 use physcon,      only:pi,solarm,pc,km
 use units,        only:set_units, unit_density
 use prompting,    only:prompt
 use dust,         only:set_dustfrac
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real :: totmass,deltax
 integer :: ipart,i,maxp,maxvxyzu
 real :: dust_to_gas
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
    print*,' uniform cubic setup...'
    print*,' enter number of particles in x (max = ',nint((maxp)**(1/3.)),')'
    read*,npartx
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx

 if (id==master) then
    print*,' enter density (gives particle mass)'
    read*,rhozero
 endif
 call bcast_mpi(rhozero)

 if (maxvxyzu < 4) then
    if (id==master) then
       print*,' enter sound speed in code units (sets polyk)'
       read*,polykset
    endif
    call bcast_mpi(polykset)
    polyk = polykset**2
    print*,' polyk = ',polyk
 else
    polyk = 0.
    polykset = 0.
 endif

 if (use_dustfrac) then
    dust_to_gas = 1.e-2
    if (id==master) call prompt(' enter dust-to-gas ratio ',dust_to_gas)
    call bcast_mpi(dust_to_gas)
 endif

 if (id==master) then
    print*,' select lattice type (1=cubic, 2=closepacked)'
    read*,ilattice
 endif
 call bcast_mpi(ilattice)

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
 print*,' npart = ',ipart,npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (mhd) then
       Bevol(:,i) = 0.
       Bevol(3,i) = 1.4142e-5
    endif
    if (use_dustfrac) then
       call set_dustfrac(dust_to_gas,dustfrac(i))
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup

