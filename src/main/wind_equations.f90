!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module wind_equations
!
! integrate the 1D wind equation to determine the initial wind profile
!
! :References: Introduction to stellar winds (Lamers & Cassinelli)
!
! :Owner: Lionel
!
! :Runtime parameters: None
!
! :Dependencies: eos, options, physcon
!

 implicit none

 public :: evolve_hydro,init_wind_equations

 private

! Wind properties
 real :: Mstar_cgs, Tstar, expT, u_to_temperature_ratio

contains

subroutine init_wind_equations (Mstar_in, Tstar_in, u_to_T)
 use physcon, only:solarm
 use eos,     only:qfacdisc
 real, intent(in) :: Mstar_in, Tstar_in, u_to_T
 Mstar_cgs = Mstar_in*solarm
 expT = 2.*qfacdisc
 Tstar = Tstar_in
 u_to_temperature_ratio = u_to_T
end subroutine init_wind_equations

subroutine evolve_hydro(dt, rvT, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, spcode, dt_force, dt_next)
!all quantities in cgs
 logical, intent(in) :: dt_force
 real, intent(in) :: mu, gamma, alpha, dalpha_dr, Q, dQ_dr, Rstar_cgs
 real, intent(inout) :: dt, rvT(3)
 integer, intent(out) :: spcode
 real, intent(out) :: dt_next

 real :: err, new_rvT(3), numerator, denominator,rold
 real, parameter :: num_tol = 1.e-4
 real, parameter :: denom_tol = 1.e-2

 rold = rvT(1)
 do
    call RK4_step_dr(dt, rvT, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
    if (dt_force) exit
    if (err > .01) then
       dt = dt * .9
    else
       !dt = dt * 1.05
       dt = min(dt*1.05,5.*abs(rold-new_rvT(1))/(1.d-3+rvT(2)))
       !dt = min(dt*1.05,0.03*(new_rvT(1))/(1.d-3+rvT(2)))
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

subroutine RK4_step_dr(dt, rvT, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
 use physcon, only:Gg,Rg,pi
 use options, only:ieos
 real, intent(in) ::  dt, rvT(3), Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr
 real, intent(out) :: err, new_rvT(3), numerator, denominator

 real :: dv1_dr,dT1_dr,dv2_dr,dT2_dr,dv3_dr,dT3_dr,dv4_dr,dT4_dr,H,r0,v0,T0,r,v,T

 r0 = rvT(1)
 v0 = rvT(2)
 T0 = rvT(3)
 H = v0*dt
 call calc_dvT_dr(r0, v0, T0, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv1_dr, dT1_dr, numerator, denominator)
 r = r0+0.5*H
 v = v0+0.5*H*dv1_dr
 T = T0+0.5*H*dT1_dr
 call calc_dvT_dr(r, v, T, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv2_dr, dT2_dr, numerator, denominator)
 r = r0+0.5*H
 v = v0+0.5*H*dv2_dr
 T = T0+0.5*H*dT2_dr
 call calc_dvT_dr(r, v, T, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv3_dr, dT3_dr, numerator, denominator)
 r = r0+H
 v = v0+dv3_dr*H
 T = T0+dT3_dr*H
 call calc_dvT_dr(r, v, T, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv4_dr, dT4_dr, numerator, denominator)
 if (dv2_dr == dv1_dr) then
    err = 0.
 else
    err = 2.*abs((dv3_dr-dv2_dr)/(dv2_dr-dv1_dr))
 endif
 new_rvT(1) = r
 new_rvT(2) = v0 + H*(dv1_dr+2.*(dv2_dr+dv3_dr)+dv4_dr)/6.
 new_rvT(3) = T0 + H*(dT1_dr+2.*(dT2_dr+dT3_dr)+dT4_dr)/6.

 ! imposed temperature profile
 if (ieos == 6) new_rvT(3) = Tstar*(Rstar_cgs/new_rvT(1))**expT
end subroutine RK4_step_dr

!--------------------------------------------------------------------------
!
!  Space derivative dv/dr and dT/dr, for Runge-Kutta (stationary solution)
!
!--------------------------------------------------------------------------
subroutine calc_dvT_dr(r, v, T, Rstar_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv_dr, dT_dr, numerator, denominator)
!all quantities in cgs
 use physcon, only:Gg,Rg,pi
 use options, only:icooling,ieos
 real, intent(in) :: r, v, T, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, Rstar_cgs
 real, intent(out) :: dv_dr, dT_dr
 real, intent(out) :: numerator, denominator

 real :: AA, BB, CC, c2, T0
 real, parameter :: denom_tol = 3.d-2 !the solution is very sensitive to this parameter!

!Temperature law
 if (ieos == 6) then
    T0 = Tstar*(Rstar_cgs/r)**expT
    c2 = gamma*Rg*T0/mu
    denominator = 1.-c2/v**2
    numerator = ((2.+expT)*r*c2 - Gg*Mstar_cgs*(1.-alpha))/(r**2*v)
    if (abs(denominator) < denom_tol) then
       AA = 2.*c2/v**3
       BB = expT*c2/(r*v)
       CC = ((2.+expT)*(1.+expT)*r*v*c2-Gg*Mstar_cgs*v*(2.-2.*alpha+r*dalpha_dr))/(r**3)
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = -expT*T0/r
 endif
 if (icooling == 0) then
    !isothermal or adiabatic expansion (no cooling)
    c2 = gamma*Rg*T/mu
    denominator = 1.-c2/v**2
    numerator = (2.*r*c2 - Gg*Mstar_cgs*(1. - alpha))/(r**2*v)
    if (abs(denominator) < denom_tol) then
       AA = ( 1.+gamma)*c2/v**3
       BB = ((4.*gamma-2.)*r*c2-Gg*Mstar_cgs*(1.-alpha))/(r**2*v**2)
       CC = ((4.*gamma-2.)*r*c2-Gg*Mstar_cgs*(2.-2.*alpha+r*dalpha_dr))/(r**3*v)
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = (1. - gamma)*T*(2.*v + r*dv_dr)/(r*v)
 else
!expansion and cooling
    c2 = gamma*Rg*T/mu
    denominator = 1.-c2/v**2
    numerator = (2.*r*c2 - Gg*Mstar_cgs*(1.-alpha))/(r**2*v) + Q*(1.-gamma)/v**2
    if (abs(denominator) < denom_tol) then
       AA = ( 1.+gamma)*c2/v**3
       BB = ((4.*gamma-2.)*r*c2-Gg*Mstar_cgs*(1.-alpha))/(r**2*v**2) &
               + (1.-gamma)*(2.+gamma)*Q/v**3
       CC = ((4.*gamma-2.)*r*c2-Gg*Mstar_cgs*(2.-2.*alpha+r*dalpha_dr))/(r**3*v) &
               - (2.*(gamma-1.)*gamma*Q)/(r*v**2) &
               + (gamma-1.)*dQ_dr/v**2
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = (1.-gamma)*T*(2.*v + r*dv_dr)/(r*v)+(gamma-1.)*Q*mu/(Rg*v)
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

end module wind_equations
