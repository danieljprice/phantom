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
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: inject, part, physcon, setbinary, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only: xyzmh_ptmass, vxyz_ptmass, nptmass
 use physcon,   only: au, solarm
 use units,     only: udist, umass, utime, set_units
 use inject,    only: wind_init, wind_gamma, mass_of_particles, central_star_mass, companion_star_mass, semi_major_axis
 use inject,    only: icompanion_star, eccentricity, central_star_radius, companion_star_r
 use setbinary, only: set_binary
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 call set_units(dist=au,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 polyk = 0.
 gamma = wind_gamma

 call wind_init(.true.)

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0

 massoftype = mass_of_particles / umass

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 if (icompanion_star > 0) then
    call set_binary(central_star_mass / umass, &
                    companion_star_mass/central_star_mass, &
                    semi_major_axis / udist, &
                    eccentricity, &
                    central_star_radius / udist, &
                    companion_star_r / udist, &
                    xyzmh_ptmass, vxyz_ptmass, nptmass)
 else
    nptmass = 1
    xyzmh_ptmass(4,1) = central_star_mass / umass
    xyzmh_ptmass(5,1) = central_star_radius / udist
 endif

 print *, "udist = ", udist, "; umass = ", umass, "; utime = ", utime

end subroutine setpart

end module setup
