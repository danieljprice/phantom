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
!  DEPENDENCIES: dim, eos, kernel, part, physcon, units
!+
!--------------------------------------------------------------------------
module bowen_dust
 implicit none

 public :: radiative_acceleration, setup_bowen, pulsating_wind_profile
 logical, parameter :: use_alpha_wind = .false.

 private
 integer, parameter :: N = 1024
 integer, parameter :: wind_emitting_sink = 1
 logical, parameter :: verbose = .false.
 real, parameter :: alpha_wind = 1.d0
 real :: L, c_light, u_to_temperature_ratio, Teff, Reff, omega_osc, &
      deltaR_osc, usteboltz, wind_injection_radius, piston_velocity, Rmin, &
      wind_velocity, pulsation_period,wind_temperature
 integer :: nwall_particles

contains


!-----------------------------------------------------------------------
!+
!  Computes acceleration due to radiation pressure on dust particles
!  it is assumed that optical depth is only depends on radius
!+
!-----------------------------------------------------------------------
subroutine radiative_acceleration(npart, xyzh, vxyzu, dt, fext, fxyzu_shock, time)
 use part,     only:isdead_or_accreted,rhoh,xyzmh_ptmass,massoftype,igas,divcurlv
 use physcon,  only:pi,mass_proton_cgs
 use dust_formation, only : kappa_dust_bowen
 use eos,      only:gamma
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(in)    :: dt,time
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(inout) :: fext(:,:),fxyzu_shock(:,:)

 real :: r(3, npart), d(npart), z_coord(npart), d2_axis(npart)
 real :: dmin, dmax, dr, part_mass
 integer :: nfound, found(npart)
 real :: rho_grid(N), r_grid(N), Teq(N), a_over_r(N), kappa_dust(N)!, T_grid(N), dlnkap_dlnT(N), tau_prime(N)
 real :: dQ_dt,Mstar,R_star,L_star,Rforce,alpha_surface,dist,frad,frepuls,Teq_part,xr(3),x_star(3)
 real :: force,fac,Trad,kap_dust,alpha_w,rho,T,mass_per_H!,Qcool,dQcool_dlnT
 integer :: i,icooling


 if (npart == nwall_particles) return
 mass_per_H = 1.4 * mass_proton_cgs
 part_mass = massoftype(igas)
 x_star(1:3) = xyzmh_ptmass(1:3,wind_emitting_sink)
 Mstar = xyzmh_ptmass(4,wind_emitting_sink)
 !R_star = Reff + deltaR_osc*sin(omega_osc*time)
 R_star = wind_injection_radius + deltaR_osc*sin(omega_osc*time)
 L_star = 4.*pi*usteboltz*R_star**2*Teff**4
 print *,'time ###  ##',time,R_star,Reff,npart,nwall_particles,L,L_star

 if (use_alpha_wind) then
    force = 1.8
    Rforce = force*R_star
!$omp parallel do default(none) &
!$omp shared(npart,xyzh,vxyzu,fext,divcurlv,x_star,Rforce,Mstar,R_star,gamma,dt,nwall_particles) &
!$omp private(xr,alpha_surface,dist)
    do i=nwall_particles+1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          xr(1:3) = xyzh(1:3,i)-x_star(1:3)
          dist = sqrt(xr(1)**2 + xr(2)**2 + xr(3)**2)
          ! repulsive force, Ayliffe & Bate 2010, MNRAS 408, 876
          if (dist < Rforce) then
             alpha_surface = ((Rforce-dist)/(Rforce-R_star))**4
          else
             alpha_surface = 0.
          endif
          fext(1:3,i) = fext(1:3,i) + Mstar*max(alpha_wind,alpha_surface)/dist**3*xr(1:3)
          !update internal energy: add pdV work
          vxyzu(4,i) = vxyzu(4,i) -(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
       endif
    enddo
!$omp end parallel do
 else
    icooling = 3
    force = 1.2
    Rforce = force*R_star
    fac = L_star/(4.*pi*c_light*Mstar)
!!$omp parallel do default(none) &
!!$omp shared(npart,xyzh,vxyzu,fext,divcurlv,x_star,Rforce,fac,Teff,R_star,Mstar,gamma,dt,nwall_particles) &
!!$omp private(xr,alpha_surface,alpha_w,kap_dust,Trad,dist)
    do i=nwall_particles+1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          xr(1:3) = xyzh(1:3,i)-x_star(1:3)
          dist = sqrt(xr(1)**2 + xr(2)**2 + xr(3)**2)
          ! repulsive force, Ayliffe & Bate 2010, MNRAS 408, 876
          if (dist < Rforce) then
             alpha_surface = ((Rforce-dist)/(Rforce-R_star))**4
          else
             alpha_surface = 0.
          endif
          Trad = Teff*sqrt(R_star/dist)
          !call bowen_dust_opacity(Trad, kap_dust)
          alpha_w = fac*kappa_dust_bowen(Trad)
          !print *,i,alpha_w,alpha_surface,kap_dust,Trad,Teff
          fext(1:3,i) = fext(1:3,i) + Mstar*max(alpha_w,alpha_surface)/dist**3*xr(1:3)
          ! !update internal energy: add pdV work
          ! vxyzu(4,i) = vxyzu(4,i) -(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
          ! !cooling terms
          ! T = vxyzu(4,i)/u_to_temperature_ratio
          ! rho = rhoh(xyzh(4,i), part_mass)
          ! call calc_cooling(icooling,mass_per_H,T,Trad,rho,Qcool,dQcool_dlnT)
          ! !print *,i,vxyzu(4,i),Qcool*dt,Trad,T
          ! if (dt  >=  Cprime/rhoh(xyzh(4,i),part_mass)) then   ! assume thermal equilibrium
          !   vxyzu(4,i) = Trad*u_to_temperature_ratio
          !   fxyzu(4,i) = 0.d0
          ! else
          !   !Qcool = -(vxyzu(4,i) - Trad*u_to_temperature_ratio) &
          !   !     * (rhoh(xyzh(4,i),part_mass)/Cprime)
          !   vxyzu(4,i) = vxyzu(4,i) + Qcool*dt
          !   !fext(4,i) = fext(4,i) + Qcool
          ! endif
       endif
    enddo
!!$omp end parallel do
    !!!!call implicit_wind_cooling(icooling,npart,xyzh,vxyzu,fxyzu_shock,fext,R_star,x_star,dt)
 endif
end subroutine radiative_acceleration

!-----------------------------------------------------------------------
!+
!  Convert parameters into code units.
!+
!-----------------------------------------------------------------------
subroutine setup_bowen(u_to_T_ratio,star_Lum,wind_radius,&
     star_Teff,piston_vamplitude,wind_speed,wind_osc_period,wind_T,nwall)
 use physcon,  only:solarl,c,steboltz,pi,radconst
 use units,    only:udist, umass, utime
 use cooling,  only:init_cooling

 integer, intent(in) :: nwall
 real, intent(in)  :: u_to_T_ratio,star_Lum,wind_radius,wind_speed,&
      star_Teff,piston_vamplitude,wind_osc_period,wind_T
 integer :: ierr

 L = star_Lum /(umass*udist**2/utime**3)
 c_light = c / (udist/utime)
 u_to_temperature_ratio = u_to_T_ratio!/(udist/utime)**2
 Teff  = star_Teff
 Reff = sqrt(star_Lum/(4.*pi*steboltz*star_Teff**4))/udist
 pulsation_period = wind_osc_period
 omega_osc = 2.*pi/pulsation_period
 deltaR_osc = piston_vamplitude/omega_osc
 nwall_particles = nwall
 wind_injection_radius = wind_radius
 piston_velocity = piston_vamplitude
 usteboltz = L/(4.*pi*Reff**2*Teff**4)
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
