!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module wind_equations
!
! integrate the 1D wind equation to determine the initial wind profile
!
! :References: Introduction to stellar winds (Lamers & Cassinelli)
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: dim, dust_formation, eos, options, physcon
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

subroutine evolve_hydro(dt, rvT, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, spcode, dt_force, dt_next)
!all quantities in cgs
 use physcon, only:au,Rg
 logical, intent(in) :: dt_force
 real, intent(in) :: mu, gamma, alpha, dalpha_dr, Q, dQ_dr, Rstar_cgs, Mdot_cgs
 real, intent(inout) :: dt, rvT(3)
 integer, intent(out) :: spcode
 real, intent(out) :: dt_next

 real :: err, new_rvT(3), numerator, denominator, rold, errmax,cs
 real, parameter :: num_tol = 1.e-4, denom_tol = 1.e-2
 real, parameter :: dt_tol = 1.e-3
 real, parameter :: rvt_tol = 1.e-2, safety = 0.9, pshrnk = -0.25, errcon = 1.89e-4, pgrow = -0.2
 character(len=3), parameter :: RK_solver = 'RK4'

 cs = sqrt(gamma*Rg*rvT(2)/mu)
 rold = rvT(1)
 err = 1.
 do while (err > rvt_tol)
    if (RK_solver == 'RK4') then
       call RK4_step_dr(dt, rvT, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
    else
       call RK6_step_dr(dt, rvT, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
    endif
    if (dt < 1.0) then
       dt_next = 1.
       exit
    endif
    if (dt_force) then
       dt_next = 0.
       exit
    endif

    ! my empirical method requires slightly less iterations
    ! if (err > rvt_tol) then
    !    dt = dt * 0.1
    ! else
    !    if (err < 1.d-3) then
    !       !dt_next = dt * 1.05
    !       dt_next = min(dt*1.25,5.*abs(rold-new_rvT(1))/(1.d-4+rvT(2)))
    !       !dt_next = min(dt*1.05,0.03*(new_rvT(1))/(1.d-3+rvT(2)))
    !    else
    !       dt_next = dt
    !    endif
    !    exit
    ! endif

    !rkqs version
    errmax = err/rvt_tol
    if (errmax > 1.) then
       dt = max(0.1*dt,safety*dt*(errmax**pshrnk))
       !exit
    else
       if (errmax > errcon) then
          dt_next = safety*dt*(errmax**pgrow)
       else
          !limit increase in dt to obtain a smooth solution
          if (abs(log(new_rvT(2)/cs)) > 1.10) then
             dt_next = dt
          else
             dt_next = 5.*dt
          endif
          !print *,new_rvt(1),(new_rvt(2)/cs),abs(log(new_rvt(2)/cs)),dt_next/dt
       endif
       exit
    endif

 enddo

!constrain timestep so the changes in r,v & T do not exceed dt_tol
 dt_next = min(dt_next,1e-2*real(au/new_rvt(2)),&
      dt_tol*dt*abs(rvt(1)/(1e-10+(new_rvt(1)-rvt(1)))),&
      dt_tol*dt*abs(rvt(2)/(1e-10+(new_rvt(2)-rvt(2)))),&
      dt_tol*dt*abs(rvt(3)/(1e-10+(new_rvt(3)-rvt(3)))))
 rvT = new_rvT

 spcode = 0
 if (numerator < -num_tol .and. denominator > -denom_tol) spcode = 1  !no solution for stationary wind
 if (numerator > -num_tol .and. denominator < -denom_tol) spcode = -1 !breeze solution
 if (denominator > denom_tol) spcode = 2                              !supersonic solution

end subroutine evolve_hydro

!--------------------------------------------------------------------------
!
!  Sixth-order Runge-Kutta integrator
!
!--------------------------------------------------------------------------
subroutine RK6_step_dr(dt, rvT, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
 use physcon, only:Gg,Rg,pi
 use options, only:ieos
 real, intent(in) ::  dt, rvT(3), Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr
 real, intent(out) :: err, new_rvT(3), numerator, denominator

 real :: dv1_dr,dv2_dr,dv3_dr,dv4_dr,dv5_dr,dv6_dr
 real :: dT1_dr,dT2_dr,dT3_dr,dT4_dr,dT5_dr,dT6_dr
 real :: H, errv, errT
 real :: r0,v0,T0,r,v,T
 real, parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875
 real, parameter :: B21=.2,B31=3./40.,B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,&
                    B52=2.5,B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,&
                    B63=575./13824.,B64=44275./110592.,B65=253./4096.
 real, parameter :: C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.
 real            :: DC1,DC3,DC4,DC5,DC6

 DC1 = C1 - 2825./27648.
 DC3 = C3 - 18575./48384.
 DC4 = C4 - 13525./55296.
 DC5 = -277./14336.
 DC6 = C6 - 0.25

 r0 = rvT(1)
 v0 = rvT(2)
 T0 = rvT(3)

 call calc_dvT_dr(r0, v0, T0, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv1_dr, dT1_dr, numerator, denominator)
 h = v0*dt
 r = r0+A2*h
 v = v0+B21*h*dv1_dr
 T = T0+B21*h*dT1_dr

 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv2_dr, dT2_dr, numerator, denominator)
 r = r0+A3*h
 v = v0+h*(B31*dv1_dr+B32*dv2_dr)
 T = T0+h*(B31*dT1_dr+B32*dT2_dr)

 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv3_dr, dT3_dr, numerator, denominator)
 r = r0+A4*h
 v = v0+h*(B41*dv1_dr+B42*dv2_dr+B43*dv3_dr)
 T = T0+h*(B41*dT1_dr+B42*dT2_dr+B43*dT3_dr)

 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv4_dr, dT4_dr, numerator, denominator)
 r = r0+A5*h
 v = v0+h*(B51*dv1_dr+B52+dv2_dr+B53*dv3_dr+B54*dv4_dr)
 T = T0+h*(B51*dT1_dr+B52+dT2_dr+B53*dT3_dr+B54*dT4_dr)


 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv5_dr, dT5_dr, numerator, denominator)
 r = r0+A6*h
 v = v0+h*(B61*dv1_dr+B62*dv2_dr+B63*dv3_dr+B64*dv4_dr+B65*dv5_dr)
 T = T0+h*(B61*dT1_dr+B62*dT2_dr+B63*dT3_dr+B64*dT4_dr+B65*dT5_dr)

 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv6_dr, dT6_dr, numerator, denominator)
 new_rvT(1) = r
 new_rvT(2) = v0+h*(C1*dv1_dr+C3*dv3_dr+C4*dv4_dr+C6*dv6_dr)
 new_rvT(3) = T0+h*(C1*dT1_dr+C3*dT3_dr+C4*dT4_dr+C6*dT6_dr)

 errT = abs((h*(DC1*dT1_dr+DC3*dT3_dr+DC4*dT4_dr+DC5*dT5_dr+DC6*dT6_dr))/T0)/10.
 errv = abs((h*(DC1*dv1_dr+DC3*dv3_dr+DC4*dv4_dr+DC5*dv5_dr+DC6*dv6_dr))/v0)
 err  = max(errT,errv)

 ! imposed temperature profile
 if (ieos == 6) new_rvT(3) = Tstar*(Rstar_cgs/new_rvT(1))**expT
end subroutine RK6_step_dr

!--------------------------------------------------------------------------
!
!  Fourth-order Runge-Kutta integrator
!
!--------------------------------------------------------------------------
subroutine RK4_step_dr(dt, rvT, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
 use physcon, only:Gg,Rg,pi
 use options, only:ieos
 real, intent(in) ::  dt, rvT(3), Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr
 real, intent(out) :: err, new_rvT(3), numerator, denominator

 real :: dv1_dr,dT1_dr,dv2_dr,dT2_dr,dv3_dr,dT3_dr,dv4_dr,dT4_dr,H,r0,v0,T0,r,v,T
! default RK4 parameters
 character(len=3), parameter :: method = 'std'
 real, parameter :: A2=.5,A3=.5,A4=1.
 real, parameter :: B21=.5,B31=0.,B32=.5,B41=0.,B42=0.,B43=1.
 real, parameter :: C1=1./6.,C2=1./3.,C3=1./3.,C4=1./6.
 real, parameter :: D1=1./6.,D2=2./3.,D3=0.,D4=1./6.
! version that mimimizes the error
! character(len=3), parameter :: method = 'err'
! real, parameter :: A2=.4,A3=.45573725,A4=1.
! real, parameter :: B21=.4,B31=.29697761,B32=.15875964,B41=.21810040,B42=-3.05096516,B43=3.83286476
! real, parameter :: C1=.17476028,C2=-.55148066,c3=1.20553560,c4=.17118478
! 3/8-rule fourth-order method
 ! character(len=3), parameter :: method = '3/8'
 ! real, parameter :: A2=1./3.,A3=2./3.,A4=1.
 ! real, parameter :: B21=1./3.,B31=-1./3.,B32=1.,B41=1.,B42=-1.,B43=1.
 ! real, parameter :: C1=1./8.,C2=3./8.,C3=3./8.,C4=1./8.
 ! character(len=3), parameter :: method = 'BSm'  !Bogacki–Shampine method
 ! real, parameter :: A2=.5,A3=.75,A4=1.
 ! real, parameter :: B21=.5,B31=0.,B32=.75,B41=2./9.,B42=1./3.,B43=4./9.
 ! real, parameter :: C1=2./9.,C2=1./3.,C3=4./9.,C4=0.
 ! real, parameter :: D1=7./24.,D2=.25,D3=1./3.,D4=1./8.
 real            :: errT,errv,deltas

 !determine RK4 solution
 r0 = rvT(1)
 v0 = rvT(2)
 T0 = rvT(3)
 H = v0*dt
 call calc_dvT_dr(r0, v0, T0, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv1_dr, dT1_dr, numerator, denominator)
 r = r0+A2*H
 v = v0+B21*H*dv1_dr
 T = T0+B21*H*dT1_dr
 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv2_dr, dT2_dr, numerator, denominator)
 r = r0+A3*H
 v = v0+H*(B31*dv1_dr+B32*dv2_dr)
 T = T0+H*(B31*dT1_dr+B32*dT2_dr)
 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv3_dr, dT3_dr, numerator, denominator)
 ! if (method == 'BSm') then
 ! new_rvT(1) = r0+H
 ! new_rvT(2) = v0 + H*(C1*dv1_dr+C2*dv2_dr+C3*dv3_dr)
 ! new_rvT(3) = T0 + H*(C1*dT1_dr+C2*dT2_dr+C3*dT3_dr)
 ! else
 r = r0+A4*H
 v = v0+H*(B41*dv1_dr+B42*dv2_dr+B43*dv3_dr)
 T = T0+H*(B41*dv1_dr+B42*dv2_dr+B43*dT3_dr)
 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv4_dr, dT4_dr, numerator, denominator)
 new_rvT(1) = r
 new_rvT(2) = v0 + H*(C1*dv1_dr+C2*dv2_dr+C3*dv3_dr+C4*dv4_dr)
 new_rvT(3) = T0 + H*(C1*dT1_dr+C2*dT2_dr+C3*dT3_dr+C4*dT4_dr)
 ! endif
 ! imposed temperature profile
 if (ieos == 6) new_rvT(3) = Tstar*(Rstar_cgs/new_rvT(1))**expT

 deltas = min(abs(new_rvt(1)-r0)/(1.d-10+r0),abs(new_rvt(2)-v0)/(1.d-10+v0),abs(new_rvt(3)-T0)/(1.d-10+T0))

 !determine RK3 solution to estimate error
 ! if (method == 'BSm') then
 ! v = new_rvT(2)
 ! T = new_rvT(3)
 ! else
 v = v0+H*(-dv1_dr+2.*dv2_dr)
 T = T0+H*(-dT1_dr+2.*dT2_dr)
 ! endif
 call calc_dvT_dr(r, v, T, Rstar_cgs, Mdot_cgs, mu, gamma, alpha, dalpha_dr, Q, dQ_dr, dv4_dr, dT4_dr, numerator, denominator)
 v = v0 + H*(D1*dv1_dr+D2*dv2_dr+D3*dv3_dr+D4*dv4_dr)
 T = T0 + H*(D1*dT1_dr+D2*dT2_dr+D3*dT3_dr+D4*dT4_dr)
 !estimate errors
 errv = abs((v-new_rvt(2))/new_rvt(2))
 errT = abs((T-new_rvt(3))/new_rvt(3))
 err  = max(deltas,errT,errv)!*100.
 !err = errv!deltas
 !new_rvT(2) = v
 !new_rvT(3) = T
 !print *,r,err,deltas,errv,errt
 !endif
end subroutine RK4_step_dr

!--------------------------------------------------------------------------
!
!  Space derivative dv/dr and dT/dr, for Runge-Kutta (stationary solution)
!
!--------------------------------------------------------------------------
subroutine calc_dvT_dr(r, v, T0, Rstar_cgs, Mdot_cgs, mu0, gamma0, alpha, dalpha_dr, Q, dQ_dr, dv_dr, dT_dr, numerator, denominator)
!all quantities in cgs
 use physcon, only:Gg,Rg,pi
 use dim,     only:update_muGamma
 use options, only:icooling,ieos
 use dust_formation,   only:calc_muGamma,idust_opacity
 real, intent(in) :: r, v, T0, mu0, gamma0, alpha, dalpha_dr, Q, dQ_dr, Rstar_cgs, Mdot_cgs
 real, intent(out) :: dv_dr, dT_dr
 real, intent(out) :: numerator, denominator

 real :: AA, BB, CC, c2, T, mu, gamma, pH, pH_tot, rho_cgs
 real, parameter :: switch_tol = 3.d-2 !the solution is very sensitive to this parameter!

 T = T0
 mu = mu0
 gamma = gamma0
 if (update_muGamma .or. idust_opacity == 2) then
    rho_cgs = Mdot_cgs/(4.*pi*r**2*v)
    call calc_muGamma(rho_cgs, T, mu, gamma, pH, pH_tot)
 endif

!Temperature law
 if (ieos == 6) then
    T = Tstar*(Rstar_cgs/r)**expT
    c2 = gamma*Rg*T/mu
    denominator = 1.-c2/v**2
    numerator = ((2.+expT)*r*c2 - Gg*Mstar_cgs*(1.-alpha))/(r**2*v)
    if (abs(denominator) < switch_tol) then
       AA = 2.*c2/v**3
       BB = (expT*c2+c2*(2.+expT)-Gg*Mstar_cgs*(1.-alpha)/r)/(r*v**2)
       CC = ((2.+expT)*(1.+expT)*c2-Gg*Mstar_cgs*(2.*(1.-alpha)/r+dalpha_dr))/(v*r**2)
       dv_dr = solve_q(AA, BB, CC)
    else
       dv_dr = numerator/denominator
    endif
    dT_dr = -expT*T/r
 endif
 if (icooling == 0) then
    !isothermal or adiabatic expansion (no cooling)
    c2 = gamma*Rg*T/mu
    denominator = 1.-c2/v**2
    numerator = (2.*r*c2 - Gg*Mstar_cgs*(1. - alpha))/(r**2*v)
    if (abs(denominator) < switch_tol) then
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
    if (abs(denominator) < switch_tol) then
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
 numerator = numerator * Rstar_cgs/sqrt(abs(c2))
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
