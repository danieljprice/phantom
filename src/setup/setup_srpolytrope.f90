!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION: Setup a polytrope in flat space, Minkowski metric (special relativity)
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!
!  DEPENDENCIES: infile_utils, io, part, physcon, setbinary, spherical,
!    timestep, units, metric, eos
!+
!--------------------------------------------------------------------------
module setup
 implicit none

 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,        only:igas,set_particle_type,rhoh
 use spherical,   only:set_sphere
 use units,       only:set_units,umass,udist
 use physcon,     only:solarm,solarr
 use io,          only:master,fatal
 use timestep,    only:tmax,dtmax
 use eos,         only:ieos
 use rho_profile, only:rho_polytrope
 use prompting,   only:prompt
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer, parameter :: ntab=5000
 integer :: i,npts,nr
 real    :: psep
 real    :: rtab(ntab),rhotab(ntab)
 real    :: densi,mstar,rstar

!-- general parameters
 time  = 0.
 polyk = 1.e-10
 gamma = 5./3.
 ieos  = 2
 tmax  = 1000.
 dtmax = 10

!-- space available for injected gas particles
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.

!-- set units
 call set_units(mass=1.e6*solarm,c=1.,G=1.)
 mstar = 1.*solarm/umass
 rstar = 1.*solarr/udist

!-- resolution
 nr    = 50
 call prompt('Resolution -- number of radial particles',nr,0)
 psep  = rstar/nr

!-- polytrope
 call rho_polytrope(gamma,polyk,mstar,rtab,rhotab,npts,set_polyk=.true.,Rstar=rstar)
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,xyz_origin=(/0.,0.,0./),rhotab=rhotab(1:npts),rtab=rtab(1:npts))

!-- mass and number of gas particles
 npartoftype(igas) = npart
 massoftype(igas)  = mstar/npart

!-- set thermal energy from density
 do i=1,npart
    call set_particle_type(i,igas)
    densi        = rhoh(xyzh(4,i),massoftype(igas))
    vxyzu(4,i)   = polyk*densi**(gamma-1.) / (gamma-1.)
 enddo

 if (id==master) print "(/,a,i10,/)",' Number of particles setup = ',npart
 if (id==master) print*,' polyk = ',polyk
 if (id==master) print*,' mstar = ',mstar
 if (npart == 0) call fatal('setup','no particles setup')

end subroutine setpart

end module setup
