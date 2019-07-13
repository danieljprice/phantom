!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: dust_formation
!
!  DESCRIPTION: Dust formation and dust cooling terms
!
!  REFERENCES: Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon
!+
!--------------------------------------------------------------------------

module wind_cooling
 implicit none
 character(len=*), parameter :: label = 'wind_cooling'

 public :: init_windcooling,calc_cooling_rate,dust_energy_cooling,write_options_windcooling
 logical, public :: calc_Teq
 real, public::    bowen_Cprime = 3.000d-5

 private
 integer, parameter :: nT = 64
 real, parameter :: Tref = 1.d5, T_floor = 1.
 real :: Tgrid(nT),Cprime
 logical :: cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan


contains

!-----------------------------------------------------------------------
!
!  calculate cooling rates
!
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Teq, gamma, mu, K2, kappa)
! all quantities in cgs
 real, intent(in) :: rho, T
 real, intent(in), optional :: Teq, gamma, mu, K2, kappa
 real, intent(out) :: Q, dlnQ_dlnT
 real :: Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan

 Q_H0 = 0.
 Q_relax_Bowen = 0.
 Q_col_dust = 0.
 Q_relax_Stefan = 0.
 dlnQ_H0 = 0.
 dlnQ_relax_Bowen = 0.
 dlnQ_col_dust = 0.
 dlnQ_relax_Stefan = 0.
 if (cool_radiation_H0) call cooling_neutral_hydrogen(T, rho, Q_H0, dlnQ_H0)
 if (cool_relaxation_Bowen) call cooling_Bowen_relaxation(T, Teq, rho, gamma, mu, Q_relax_Bowen, dlnQ_relax_Bowen)
 if (cool_collisions_dust)  call cooling_collision_dust_grains(T, Teq, rho, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (cool_relaxation_Stefan) call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
 Q = Q_H0 + Q_relax_Bowen+ Q_col_dust+ Q_relax_Stefan
 dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen+ Q_col_dust*dlnQ_col_dust+ Q_relax_Stefan*dlnQ_relax_Stefan)/Q
end subroutine calc_cooling_rate

subroutine init_windcooling(icool)
  use units,        only:umass, utime, udist
  use coolfunc

  integer, intent(in) :: icool
  integer :: iwind
  print *,'init windcool'
  if (icool > 0) then
     iwind = icool
  else
     !dust free wind, only H0 cooling allowed
     if (abs(icool) > 0) iwind = 10
  endif
  Cprime = bowen_Cprime / (umass*utime/udist**3)

  cool_radiation_H0 = .false.
  cool_relaxation_Bowen = .false.
  cool_relaxation_Stefan = .false.
  cool_collisions_dust = .false.
  !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
  if (iwind > 9) cool_radiation_H0 = .true.
  select case (iwind)
  case (2,12)
     cool_relaxation_Stefan = .true.
  case (3,13)
     cool_relaxation_Bowen = .true.
  case (4,14)
     cool_collisions_dust = .true.
  case (6,16)
     cool_collisions_dust = .true.
     cool_relaxation_Stefan = .true.
  case (7,17)
     cool_collisions_dust = .true.
     cool_relaxation_Bowen = .true.
  end select
  !krome calculates its own cooling rate
#ifdef KROME
  cool_radiation_H0 = .false.
  cool_collisions_dust = .false.
#endif
  calc_Teq = cool_relaxation_Bowen .or. cool_relaxation_Stefan .or. cool_collisions_dust

  !initialize grid temperature
  call set_Tcool
  call init_coolfunc(iwind)
end subroutine init_windcooling


!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Teq, rho, wind_gamma, mu, Q, dlnQ_dlnT)
 use physcon, only: Rg
 real, intent(in) :: T, Teq, rho, wind_gamma, mu
 real, intent(out) :: Q,dlnQ_dlnT

 Q = Rg/((wind_gamma-1.)*mu)*rho*(Teq-T)/Cprime
 dlnQ_dlnT = -T/(Teq-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_collision_dust_grains(T, Teq, rho, K2, mu, Q, dlnQ_dlnT)
 use physcon, only: kboltz, mass_proton_cgs, pi
 real, intent(in) :: T, Teq, rho, K2, mu
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter :: f = 0.15, a0 = 1.28e-8
 real :: A

 A = 2. * f * kboltz * a0**2/(mass_proton_cgs**2*mu) &
         * (1.05/1.54) * sqrt(2.*pi*kboltz/mass_proton_cgs) * 2.*K2 * rho
 Q = A * sqrt(T) * (Teq-T)
 if (Q  >  1.d6) then
    print *, f, kboltz, a0, mass_proton_cgs, mu
    print *, mu, K2, rho, T, Teq, A, Q
    stop 'cooling'
 else
    dlnQ_dlnT = 0.5+T/(Teq-T+1.d-10)
 endif
end subroutine cooling_collision_dust_grains

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_radiative_relaxation(T, Teq, kappa, Q, dlnQ_dlnT)
 use physcon, only: steboltz
 real, intent(in) :: T, Teq, kappa
 real, intent(out) :: Q,dlnQ_dlnT

 Q = 4.*steboltz*(Teq**4-T**4)*kappa
 dlnQ_dlnT = -4.*T**4/(Teq**4-T**4+1.d-10)

end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to neutral H (Spitzer)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen(T, rho, Q, dlnQ_dlnT)
 use physcon, only: kboltz, mass_proton_cgs, pi
 real, intent(in) :: T, rho
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter :: f = 0.2! 1.d0 !0.2
 real :: eps_e

 if (T > 3000.) then
    eps_e = calc_eps_e(T)
    !Q = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(mass_per_H)**2
    Q = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(1.4*mass_proton_cgs)**2
    dlnQ_dlnT = 118400.d0/T+log(calc_eps_e(1.001*T)/eps_e)/log(1.001)
 else
    Q = 0.
    dlnQ_dlnT = 0.
 endif
end subroutine cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance (Palla et al 1983)
!+
!-----------------------------------------------------------------------
real function calc_eps_e(T)
 real, intent(in) :: T
 real :: k1, k2, k3, k8, k9, p, q

 if (T > 3000.) then
    k1 = 1.88d-10 / T**6.44e-1
    k2 = 1.83d-18 * T
    k3 = 1.35d-9
    k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
    k9 = 1.7d-4 * k8
    p = .5*k8/k9
    q = k1*(k2+k3)/(k3*k9)
    calc_eps_e = (p + sqrt(q+p**2))/q
 else
    calc_eps_e = 1.d-30
 endif
end function calc_eps_e

subroutine set_Tcool
  integer :: i
  real :: dlnT
  dlnT = log(Tref)/(nT-1)

  do i = 1,nT
     Tgrid(i) = exp((i-1)*dlnT)
  enddo
end subroutine set_Tcool

subroutine dust_energy_cooling1 (u, rho, dt, gam_in, mu_in, Teq, K2, kappa)
  use coolfunc
  real, intent(in) :: rho, dt
  real, intent(in), optional :: Teq, K2, kappa, gam_in, mu_in
  real, intent(inout) :: u
  real :: dudt
  call  energ_coolfunc(u,rho,dt,dudt)
end subroutine

!-----------------------------------------------------------------------
!
!   dust cooling using Townsend method to avoid the timestep constraints
!
!-----------------------------------------------------------------------
subroutine dust_energy_cooling (u, rho, dt, gam_in, mu_in, Teq, K2, kappa)
  use eos,     only:gamma,gmw
  use physcon, only:Rg,kboltz
  use units,   only:unit_density,unit_ergg,utime
  real, intent(in) :: rho, dt
  real, intent(in), optional :: Teq, K2, kappa, gam_in, mu_in
  real, intent(inout) :: u

  real :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,dt_cgs,rho_cgs,gam,mu,du
  integer :: k

  if (.not.present(gam_in)) then
     gam = gamma
     mu = gmw
  else
     gam = gam_in
     mu = mu_in
  endif

  T = (gam-1.)*u/Rg*mu*unit_ergg

  if (T < T_floor) then
     Temp = T_floor
  elseif (T > Tref) then
     call calc_cooling_rate(Q, dlnQ_dlnT, rho_cgs, T, Teq, gam, mu, K2, kappa)
     Temp = T+(gam-1.)*mu*Q*dt_cgs/Rg
  else
     dt_cgs = dt*utime
     rho_cgs = rho*unit_density
     call calc_cooling_rate(Qref,dlnQref_dlnT, rho_cgs, Tref, Teq, gam, mu, K2, kappa)
     Y = 0.
     k = nT
     do while (Tgrid(k) > T)
        k = k-1
        call calc_cooling_rate(Q, dlnQ_dlnT, rho_cgs, Tgrid(k), Teq, gam, mu, K2, kappa)
        if (abs(dlnQ_dlnT-1.) < 1.d-13) then
           y = y - dlnQref_dlnT*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
           !y = y - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
         else
           !y = y - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
           y = y - dlnQref_dlnT*Tgrid(k)/(dlnQ_dlnT*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
        endif
        !print *,'@',y,q*dt_cgs/(u*unit_ergg),tgrid(k),dlnq_dlnt,dlnQref_dlnT
     enddo

     yk = y
     !call calc_cooling_rate(Q, dlnQ_dlnT, rho_cgs, T, Teq, gam, mu, K2, kappa)

     if (abs(dlnQ_dlnT-1.) < 1.d-13) then
        !y = yk + Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
        y = yk + dlnQref_dlnT*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
     else
        !y = yk + Qref*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1))
        y = yk + dlnQref_dlnT*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1))
     endif

     !dy = mu*(gam-1.)/Rg*Qref/Tref*dt_cgs
     dy = mu*(gam-1.)/Rg*dlnQref_dlnT/Tref*dt_cgs
     y = y + dy
!compute Yinv
     if (abs(dlnQ_dlnT-1.) < 1.d-10) then
        Temp = max(Tgrid(k)*exp(-Q*Tref*(y-yk)/(Qref*Tgrid(k))),T_floor)
     else
        Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref/(Qref*Tgrid(k))*(y-yk)
        if (Yinv > 0.) then
           Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
        else
           Temp = T_floor
        endif
     endif
  endif

  du = Rg*(Temp-T)/((gam-1.)*mu*unit_ergg)
  if (abs(du/u) > .1) then
     !print *,u,dlnQ_dlnT,T,temp,rho_cgs
     !stop
  endif
  u = u+du
  !stop

end subroutine

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_windcooling(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

#ifdef NUCLEATION
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
#endif

end subroutine write_options_windcooling

end module wind_cooling
