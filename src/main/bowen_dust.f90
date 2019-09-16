!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: bowen_dust
!
!  DESCRIPTION:
!  Acceleration on gas due to dust grains at radiative equilibrium
!  The model is described in the following article:
!  G.H. Bowen - Dynamical modeling of long-period variable star atmospheres (1988)
!
!  REFERENCES: None
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: cooling, eos, physcon, units
!+
!--------------------------------------------------------------------------
module bowen_dust
 implicit none

 public :: setup_bowen, pulsating_wind_profile

 private
 integer, parameter :: wind_emitting_sink = 1
 real :: u_to_temperature_ratio, omega_osc, deltaR_osc, wind_injection_radius, &
      piston_velocity, Rmin, wind_velocity, pulsation_period,wind_temperature
 integer :: nwall_particles

contains


!-----------------------------------------------------------------------
!+
!  Convert parameters into code units.
!+
!-----------------------------------------------------------------------
subroutine setup_bowen(u_to_T_ratio,wind_radius,&
     piston_vamplitude,wind_speed,wind_osc_period,wind_T,nwall)
 use physcon,  only:solarl,c,steboltz,pi,radconst
 use units,    only:udist, umass, utime
 use cooling,  only:init_cooling

 integer, intent(in) :: nwall
 real, intent(in)  :: u_to_T_ratio,wind_radius,wind_speed,&
      piston_vamplitude,wind_osc_period,wind_T
 integer :: ierr

 u_to_temperature_ratio = u_to_T_ratio!/(udist/utime)**2
 pulsation_period = wind_osc_period
 omega_osc = 2.*pi/pulsation_period
 deltaR_osc = piston_vamplitude/omega_osc
 nwall_particles = nwall
 wind_injection_radius = wind_radius
 piston_velocity = piston_vamplitude
 !cls Rmin = Reff - deltaR_osc
 Rmin = wind_injection_radius - deltaR_osc
 wind_velocity = wind_speed
 wind_temperature = wind_T

 call init_cooling(ierr)

end subroutine setup_bowen


!-----------------------------------------------------------------------
!+
!  Oscillating inner boundary : bowen wind
!+
!-----------------------------------------------------------------------
subroutine pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,sphere_number, &
                                  inner_sphere,inner_boundary_sphere,dr3,rho_ini)
 use physcon,     only:pi
 use eos,         only:gamma
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,GM,dr3,rho_ini
 real,    intent(out) :: r, v, u, rho, e

 integer, parameter :: nrho_index = 10
 integer :: k
 real :: surface_radius,r3
 logical :: verbose = .true.

 v = wind_velocity + piston_velocity* cos(omega_osc*time) !same velocity for all wall particles
 surface_radius = wind_injection_radius + deltaR_osc*sin(omega_osc*time)
 !ejected spheres
 if (sphere_number <= inner_sphere) then
    r = surface_radius
    v = max(piston_velocity,wind_velocity)
 else
    !boundary spheres
    r3 = surface_radius**3-dr3
    do k = 2,sphere_number-inner_sphere
       r3 = r3-dr3*(r3/surface_radius**3)**(nrho_index/3.)
    enddo
    r = r3**(1./3)
 endif
 !r = (surface_radius**3-(sphere_number-inner_sphere)*dr3)**(1./3)
 !rho = rho_ini
 u = wind_temperature * u_to_temperature_ratio
 if (gamma > 1.0001) then
    e = .5*v**2 - GM/r + gamma*u
 else
    e = .5*v**2 - GM/r + u
 endif
 rho = rho_ini*(surface_radius/r)**nrho_index
 if (verbose) then
    if (sphere_number > inner_sphere) then
       print '("boundary, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4)', &
             sphere_number,inner_sphere,surface_radius,r,v,time/pulsation_period
    else
       print '("ejected, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4)', &
             sphere_number,inner_sphere,surface_radius,r,v,time/pulsation_period
    endif
 endif

end subroutine pulsating_wind_profile

end module bowen_dust
