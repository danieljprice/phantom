!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: wind_profile
!
!  DESCRIPTION: integrate the 1D wind equation to determine the initial wind profile
!
!  REFERENCES: None
!
!  OWNER: Lionel
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dust_formation, eos, physcon, units
!+
!--------------------------------------------------------------------------
module wind_profile

#ifdef BOWEN
 public :: pulsating_wind_profile
#endif
 public :: evolve_hydro,energy_profile
 public :: stationary_wind_profile
 public :: init_wind_profile

 private

! Wind properties
 real :: Mstar_cgs, Lstar_cgs, Tstar, Rstar_cgs, wind_gamma, wind_mass_rate
 real :: Cprime, u_to_temperature_ratio, expT, alpha, wind_temperature,Rstar
 integer :: wind_type

contains

subroutine init_wind_profile(Mstar_in, Lstar_in, Tstar_in, Rstar_in, Cprime_cgs, &
     expT_in, Mdot_in, CO_ratio, u_to_T, alpha_in, Twind, wind_type_in)
 use units,   only:udist
 use physcon, only:c, solarm, years
 use eos,     only:gamma
 use dust_formation, only: set_abundances,set_cooling
 real, intent(in) :: Mstar_in, Lstar_in, Tstar_in, Rstar_in, Cprime_cgs, expT_in, Mdot_in, CO_ratio, &
          u_to_T,alpha_in, Twind
 integer, intent(in) :: wind_type_in
 logical :: cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan

 Mstar_cgs = Mstar_in*solarm
 Lstar_cgs = Lstar_in
 Tstar = Tstar_in
 wind_type = wind_type_in
 alpha = alpha_in
 Cprime = Cprime_cgs ! keep it in cgs
 if (wind_type == 2) then
    expT = expT_in
    wind_gamma = 1.
 else if (wind_type == 1) then
    wind_gamma = 1.
 else
    wind_gamma = gamma
 endif
 wind_temperature = Twind
 wind_mass_rate = Mdot_in !code units
 Rstar_cgs = Rstar_in
 Rstar = Rstar_in/udist
 u_to_temperature_ratio = u_to_T
 call set_abundances(CO_ratio)

 calc_dust = wind_type > 10

 cool_radiation_H0 = .false.
 if (wind_type == 2) then
    cool_relaxation_Bowen = .false.
 else
    cool_relaxation_Bowen = .true.
 endif
 cool_collisions_dust = .false.
 cool_relaxation_Stefan = .false.
 call set_cooling(cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan)
end subroutine init_wind_profile

! #ifdef NUCLEATION
! !-----------------------------------------------------------------------
! !+
! !  dusty wind model
! !+
! !-----------------------------------------------------------------------
! subroutine dusty_wind_profile(time,local_time,r,v,T,u,rho,e,GM,gamma,Jstar,K,mu,cs)
!  !in/out variables in code units (except Jstar,K,mu)
!  use gailstatwind, only:calc_wind_profile,wind_state
!  use units,        only:udist, utime, unit_velocity, unit_density, unit_pressure
!  real,    intent(in)  :: time,local_time,gamma,GM,T
!  real,    intent(inout) :: r, v
!  real,    intent(out) :: u, rho, e, Jstar, K(:), mu, cs

!  type(wind_state) :: state

!  call calc_wind_profile(r, v, T, local_time*utime, state)
!  r = state%r/udist
!  v = state%v/unit_velocity
!  rho = state%rho/(unit_density)
!  u = state%p/((gamma-1.)*rho)/unit_pressure
!  e = .5*v**2 - GM/r + gamma*u
!  Jstar = state%Jstar
!  K = state%K
!  mu = state%mu
!  cs = state%c/unit_velocity
! end subroutine dusty_wind_profile

! #elif BOWEN
#ifdef BOWEN

!-----------------------------------------------------------------------
!+
!  Oscillating inner boundary : bowen wind
!+
!-----------------------------------------------------------------------
subroutine pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,gamma,sphere_number, &
                                     inner_sphere,inner_boundary_sphere)
 use physcon,     only: pi
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,gamma,GM
 real,    intent(out) :: r, v, u, rho, e

 integer, parameter :: nrho_index = 10
 integer :: k
 real :: surface_radius,r3
 logical :: verbose = .true.

 v = wind_velocity + wind_osc_vamplitude* cos(2.*pi*time/wind_osc_period) !same velocity for all wall particles
 surface_radius = wind_injection_radius + wind_osc_vamplitude*wind_osc_period/(2.*pi)*sin(2.*pi*time/wind_osc_period)
 if (sphere_number <= inner_sphere) then
    r = surface_radius
    v = max(wind_osc_vamplitude,wind_velocity)
 else
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
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)', &
             sphere_number,inner_sphere,surface_radius,r,v,&
             time/wind_osc_period,time_between_spheres/wind_osc_period
    else
       print '("ejected, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)', &
             sphere_number,inner_sphere,surface_radius,r,v,&
             time/wind_osc_period,time_between_spheres/wind_osc_period
    endif
 endif

end subroutine pulsating_wind_profile

#else

!-----------------------------------------------------------------------
!+
!  stationary wind - solution when input velocity is defined
!+
!-----------------------------------------------------------------------
subroutine stationary_wind_profile(local_time, r, v, u, rho, e, GM, gamma, mu)
 use physcon,     only: pi
 use units,       only:udist,unit_velocity
 real, intent(in)  :: local_time, GM, gamma, mu
 real, intent(inout) :: r, v
 real, intent(out) ::  u, rho, e
 real :: dt, r0, v0, T, rvT(3), new_rvT(3), err, Q, dQ_dr, dalpha_dr, numerator, denominator
 integer, parameter :: N = 10000
 integer :: i

 dt = local_time / N
 r0 = r
 v0 = v
 Q = 0.
 dQ_dr = 0.
 rvT(1) = r*udist
 rvT(2) = v*unit_velocity
 rvT(3) = T
 do i=1,N
    call RK4_step_dr(dt, rvT, mu, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
    rvT = new_rvT
 enddo
 r = new_rvt(1)/udist
 v = new_rvt(2)/unit_velocity
 T = new_rvt(3)
 if (gamma > 1.0001) then
    T = wind_temperature * (r0**2 * v0 / (r**2 * v))**(gamma-1.)
    u = T * u_to_temperature_ratio
    e = .5*v**2 - GM/r + gamma*u
 else
    u = T * u_to_temperature_ratio
    e = .5*v**2 - GM/r + u
 endif
 !update radius, velocity and density
 rho = wind_mass_rate / (4.*pi*r**2*v)

end subroutine stationary_wind_profile

#endif

subroutine evolve_hydro(dt, rvT, mu, alpha, dalpha_dr, Q, dQ_dr, spcode, dt_force, dt_max, dt_next)
 use physcon, only:kboltz,atomic_mass_unit,Gg
 logical, intent(in) :: dt_force
 real, intent(in) :: mu, alpha, dalpha_dr, Q, dQ_dr, dt_max
 real, intent(inout) :: dt, rvT(3)
 integer, intent(out) :: spcode
 real, intent(out) :: dt_next

 real :: err, new_rvT(3), numerator, denominator
 real, parameter :: num_tol = 1.e-3
 real, parameter :: denom_tol = 1.e-2

 do
    call RK4_step_dr(dt, rvT, mu, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
    if (dt_force) exit
    if (err > .01) then
       dt = dt * .9
    else
       if (dt*1.05 < dt_max) dt = dt * 1.05
       dt = min(dt,0.03*rvT(1)/(1.d-3+rvT(2)))
       exit
    endif
 enddo
 rvT = new_rvT
 dt_next = dt

 spcode = 0
 if (numerator < -num_tol .and. denominator > -denom_tol) spcode = 1  !no solution for stationary wind
 if (numerator > -num_tol .and. denominator < -denom_tol) spcode = -1 !breeze solution
 if (denominator > denom_tol) spcode = 2

end subroutine evolve_hydro

subroutine RK4_step_dr(dt, rvT, mu, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
 use physcon, only:Gg,kboltz,atomic_mass_unit,pi
 real, intent(in) ::  dt, rvT(3), mu, alpha, dalpha_dr, Q, dQ_dr
 real, intent(out) :: err, new_rvT(3), numerator, denominator

 real :: dv1_dr,dT1_dr,dv2_dr,dT2_dr,dv3_dr,dT3_dr,dv4_dr,dT4_dr,H,r0,v0,T0,r,v,T
 real, parameter :: Rgas = kboltz/atomic_mass_unit

 r0 = rvT(1)
 v0 = rvT(2)
 T0 = rvT(3)
 H = v0*dt
 call calc_dvT_dr(r0, v0, T0, mu, alpha, dalpha_dr, Q, dQ_dr, dv1_dr, dT1_dr, numerator, denominator)
 r = r0+0.5*H
 v = v0+0.5*H*dv1_dr
 T = T0+0.5*H*dT1_dr
 call calc_dvT_dr(r, v, T, mu, alpha, dalpha_dr, Q, dQ_dr, dv2_dr, dT2_dr, numerator, denominator)
 r = r0+0.5*H
 v = v0+0.5*H*dv2_dr
 T = T0+0.5*H*dT2_dr
 call calc_dvT_dr(r, v, T, mu, alpha, dalpha_dr, Q, dQ_dr, dv3_dr, dT3_dr, numerator, denominator)
 r = r0+H
 v = v0+dv3_dr*H
 T = T0+dT3_dr*H
 call calc_dvT_dr(r, v, T, mu, alpha, dalpha_dr, Q, dQ_dr, dv4_dr, dT4_dr, numerator, denominator)
 if (dv2_dr == dv1_dr) then
    err = 0.
 else
    err = 2.*abs((dv3_dr-dv2_dr)/(dv2_dr-dv1_dr))
 endif
 new_rvT(1) = r
 new_rvT(2) = v0 + H*(dv1_dr+2.*(dv2_dr+dv3_dr)+dv4_dr)/6.
 new_rvT(3) = T0 + H*(dT1_dr+2.*(dT2_dr+dT3_dr)+dT4_dr)/6.
 if (wind_type == 2) new_rvT(3) = Tstar*(Rstar_cgs/new_rvT(1))**expT
end subroutine RK4_step_dr

!-------------------------------------------------------------------
!
!  Time derivative of r and v, for Runge-Kutta (stationary solution)
!
!-------------------------------------------------------------------
subroutine calc_dvT_dr(r, v, T, mu, alpha, dalpha_dr, Q, dQ_dr, dv_dr, dT_dr, numerator, denominator)
!all quantities in cgs
 use physcon, only:Gg,kboltz,atomic_mass_unit,pi
 real, intent(in) :: r, v, T, mu, alpha, dalpha_dr, Q, dQ_dr
 real, intent(out) :: dv_dr, dT_dr
 real, intent(out) :: numerator, denominator

 real :: AA, BB, CC, c2, T0
 real, parameter :: Rgas = kboltz/atomic_mass_unit
 real, parameter :: denom_tol = 1.d-2

!Temperature law
 if (wind_type == 2) then
    T0 = Tstar*(Rstar_cgs/r)**expT
    c2 = wind_gamma*Rgas*T0/mu
    denominator = 1.-c2/v**2
    numerator = ((2.+expT)*r*c2 - Gg*Mstar_cgs*(1.-alpha))/(r**2*v)
    if (abs(denominator) < denom_tol) then
       AA = 2.*c2/v**3
       !BB = (2.*r*c2*(1.+expt)-Gg*Mstar_cgs*(1.-alpha))/(r**2*v**2)
       !CC = ((2.+expt)*(1.+expt)*r*c2-Gg*Mstar_cgs*(2.-2.*alpha+r*dalpha_dr))/(r**3*v)
       !dv_dr = solve_q(AA, BB, CC)
       BB = expT*c2/(r*v)
       CC = ((2.+expT)*(1.+expT)*r*v*c2-Gg*Mstar_cgs*v*(2.-2.*alpha+r*dalpha_dr))/(r**3)
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = -expT*T/r
 endif
 !isothermal or adiabatic expansion (no cooling)
 if (wind_type == 1 .or. wind_type == 3) then
    c2 = wind_gamma*Rgas*T/mu
    denominator = 1.-c2/v**2
    numerator = (2.*r*c2 - Gg*Mstar_cgs*(1. - alpha))/(r**2*v)
    if (abs(denominator) < denom_tol) then
       AA = (1.+wind_gamma)*c2/v**3
       BB = ((4.*wind_gamma-2.)*r*c2-Gg*Mstar_cgs*(1.-alpha))/(r**2*v**2)
       CC = ((4.*wind_gamma-2.)*r*c2-Gg*Mstar_cgs*(2.-2.*alpha+r*dalpha_dr))/(r**3*v)
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = (1. - wind_gamma)*T*(2.*v + r*dv_dr)/(r*v)
 endif
!expansion and cooling
 if (wind_type == 4) then
    c2 = wind_gamma*Rgas*T/mu
    denominator = 1.-c2/v**2
    numerator = (2.*r*c2 - Gg*Mstar_cgs*(1.-alpha))/(r**2*v) + Q*(1.-wind_gamma)/v**2
    if (abs(denominator) < denom_tol) then
       AA = (1.+wind_gamma)*c2/v**3
       BB = ((4.*wind_gamma-2.)*r*c2-Gg*Mstar_cgs*(1.-alpha))/(r**2*v**2) &
               + (1.-wind_gamma)*(2.+wind_gamma)*Q/v**3
       CC = ((4.*wind_gamma-2.)*r*c2-Gg*Mstar_cgs*(2.-2.*alpha+r*dalpha_dr))/(r**3*v) &
               - (2.*(wind_gamma-1.)*wind_gamma*Q)/(r*v**2) &
               + (wind_gamma-1.)*dQ_dr/v**2
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = (1. - wind_gamma)*T*(2.*v + r*dv_dr)/(r*v) &
            + (-1.+wind_gamma)*Q*mu/(Rgas*v)
 endif
 numerator = numerator * Rstar_cgs/sqrt(c2)
end subroutine calc_dvT_dr

pure real function solve_q(a, b, c)
 real, intent(in) :: a, b, c
 real :: delta
 if (-4.*a*c/b**2 > 1.d-8) then
    delta = max(b**2-4.*a*c, 0.)
    solve_q = (-b+sqrt(delta))/(2.*a)
 else
    solve_q = -c/b
 endif
end function solve_q

real function energy_profile(xyzh)
 real, intent(in) :: xyzh(4)
 real :: r
 r = sqrt(xyzh(1)**2+xyzh(2)**2+xyzh(3)**2)
 energy_profile = u_to_temperature_ratio*Tstar*(Rstar/r)**expT
end function energy_profile
end module wind_profile
