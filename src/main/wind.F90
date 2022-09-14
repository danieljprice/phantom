!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module wind
!
! driver to integrate the wind equations
!
! :References: Lamers & Cassinelli "Introduction to stellar winds"
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: cooling_solver, dim, dust_formation, eos, io, options,
!   part, physcon, ptmass_radiation, timestep, units, wind_equations
!

!#define CALC_HYDRO_THEN_CHEM

 use part,only: n_nucleation,idJstar,idK0,idK1,idK2,idK3,idmu,idgamma,idsat,idkappa,idalpha

 implicit none
 public :: setup_wind
 public :: wind_state,save_windprofile,interp_wind_profile!,wind_profile

 private
 ! Shared variables
 real, parameter :: Tdust_stop = 1.d-2 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'wind'

 ! input parameters
 real :: Mstar_cgs, Lstar_cgs, wind_gamma, Mdot_cgs, Tstar
 real :: u_to_temperature_ratio
 real, dimension (:,:), allocatable, public :: trvurho_1D, JKmuS_1D

 ! wind properties
 type wind_state
    real :: dt, time, r, r0, Rstar, v, a, time_end, Tg, Tdust
    real :: alpha, rho, p, u, c, dalpha_dr, r_old, Q, dQ_dr
    real :: kappa, mu, gamma, tau_lucy, alpha_Edd, tau
    real :: JKmuS(n_nucleation),K2
    integer :: spcode, nsteps
    logical :: dt_force, error, find_sonic_solution
 end type wind_state
contains

subroutine setup_wind(Mstar_cg, Mdot_code, u_to_T, r0, T0, v0, rsonic, tsonic, stype)
 use ptmass_radiation, only:iget_tdust
 use units,            only:umass,utime
 use physcon,          only:c,years,pi
 use eos,              only:gamma,gmw
 use dust_formation,   only:set_abundances,init_muGamma,idust_opacity

 integer, intent(in) :: stype
 real, intent(in)    :: Mstar_cg, Mdot_code, u_to_T
 real, intent(inout) :: r0,v0,T0
 real, intent(out)   :: rsonic, tsonic
 real :: rho_cgs, wind_mu!, pH, pH_tot

 Mstar_cgs  = Mstar_cg
 wind_gamma = gamma
 Mdot_cgs   = Mdot_code * umass/utime
 u_to_temperature_ratio = u_to_T

 if (idust_opacity == 2) then
    call set_abundances
    rho_cgs    = Mdot_cgs/(4.*pi * r0**2 * v0)
    wind_gamma = gamma
    wind_mu    = gmw
    call init_muGamma(rho_cgs, T0, wind_mu, wind_gamma)
    print *,'reset gamma : ',gamma,wind_gamma
    print *,'reset gmw   : ',gmw,wind_mu
    gamma = wind_gamma
    gmw   = wind_mu
 endif

 if (iget_tdust == -12) then
    !not working
    print *,'get_initial_radius not working'
    call get_initial_radius(r0, T0, v0, rsonic, tsonic, stype)
 else
    call get_initial_wind_speed(r0, T0, v0, rsonic, tsonic, stype)
 endif

end subroutine setup_wind

!-----------------------------------------------------------------------
!
!  Initialize variables for wind integration
!
!-----------------------------------------------------------------------
subroutine init_wind(r0, v0, T0, time_end, state)
! all quantities in cgs
 use physcon,          only:pi,Rg
 use io,               only:fatal
 use eos,              only:gmw
 use ptmass_radiation, only:alpha_rad
 use part,             only:xyzmh_ptmass,iTeff,ilum
 use dust_formation,   only:kappa_gas,init_muGamma,idust_opacity
 use units,            only:umass,unit_energ,utime

 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state

 state%dt = 1000.
 if (time_end > 0.d0) then
    ! integration stops when time = time_end
    state%find_sonic_solution = .false.
    state%time_end = time_end
 else
    ! integration stops once the sonic point is reached
    state%find_sonic_solution = .true.
    state%time_end = -1.d0
 endif
 state%time   = 0.
 state%r_old  = 0.
 state%r0     = r0
 state%r      = r0
 state%Rstar  = r0
 state%v      = v0
 state%a      = 0.
 state%Tg     = T0
 state%Tdust  = T0
 state%tau_lucy = 2./3.
 state%kappa  = kappa_gas
 state%Q      = 0.
 state%dQ_dr  = 0.
 state%rho    = Mdot_cgs/(4.*pi * state%r**2 * state%v)
 state%spcode = 0
 state%nsteps = 1
 state%tau    = 0

 Tstar     = xyzmh_ptmass(iTeff,wind_emitting_sink)
 Lstar_cgs = xyzmh_ptmass(ilum,wind_emitting_sink)*unit_energ/utime
 Mstar_cgs = xyzmh_ptmass(4,wind_emitting_sink)*umass

 state%alpha_Edd = 0.
 state%K2        = 0.
 state%mu        = gmw
 state%gamma     = wind_gamma
 state%JKmuS     = 0.
 if (idust_opacity == 2) call init_muGamma(state%rho, state%Tg, state%mu, state%gamma)
 state%alpha          = state%alpha_Edd + alpha_rad
 state%JKmuS(idalpha) = state%alpha
 state%JKmuS(idmu)    = state%mu
 state%JKmuS(idgamma) = state%gamma
 state%dalpha_dr = 0.
 state%p = state%rho*Rg*state%Tg/state%mu
#ifdef ISOTHERMAL
 state%u = 1.
#else
 state%u = state%p/((state%gamma-1.)*state%rho)
#endif
 state%c = sqrt(state%gamma*Rg*state%Tg/state%mu)
 state%dt_force = .false.
 state%error    = .false.

end subroutine init_wind

#ifdef CALC_HYDRO_THEN_CHEM

!-----------------------------------------------------------------------
!
!  Integrate chemistry, cooling and hydro over one time step
!
!-----------------------------------------------------------------------
subroutine wind_step(state)
! all quantities in cgs

 use wind_equations,   only:evolve_hydro
 use ptmass_radiation, only:alpha_rad,iget_tdust,tdust_exp
 use physcon,          only:pi,Rg
 use dust_formation,   only:evolve_chem,calc_kappa_dust,calc_kappa_bowen,&
      calc_Eddington_factor,idust_opacity
 use part,             only:idK3,idmu,idgamma,idsat,idkappa
 use cooling_solver,   only:calc_cooling_rate
 use options,          only:icooling
 use units,            only:unit_ergg,unit_density

 type(wind_state), intent(inout) :: state
 real :: rvT(3), dt_next, v_old, dlnQ_dlnT, Q_code
 real :: alpha_old, kappa_old, rho_old, Q_old, tau_lucy_bounded, mu_old, dt_old

 rvT(1) = state%r
 rvT(2) = state%v
 rvT(3) = state%Tg
 v_old  = state%v
 dt_old = state%dt
 state%r_old = state%r
 call evolve_hydro(state%dt, rvT, state%Rstar, Mdot_cgs, state%mu, state%gamma, state%alpha, state%dalpha_dr, &
      state%Q, state%dQ_dr, state%spcode, state%dt_force, dt_next)
 state%r   = rvT(1)
 state%v   = rvT(2)
 state%a   = (state%v-v_old)/(state%dt)
 state%Tg  = rvT(3)
 !state%u   = rvT(3)
 rho_old   = state%rho
 state%rho = Mdot_cgs/(4.*pi*state%r**2*state%v)

 kappa_old = state%kappa
 alpha_old = state%alpha
 mu_old    = state%mu
 if (idust_opacity == 2) then
    !state%Tg   = state%u*(state%gamma-1.)/Rg*state%mu
    call evolve_chem(state%dt,state%Tg,state%rho,state%JKmuS)
    state%K2             = state%JKmuS(idK2)
    state%mu             = state%JKmuS(idmu)
    state%gamma          = state%JKmuS(idgamma)
    state%kappa          = calc_kappa_dust(state%JKmuS(idK3), state%Tdust, state%rho)
    state%alpha_Edd      = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, state%kappa)
    state%JKmuS(idalpha) = state%alpha_Edd+alpha_rad
 elseif (idust_opacity == 1) then
    state%kappa     = calc_kappa_bowen(state%Tdust)
    state%alpha_Edd = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, state%kappa)
 endif
 state%alpha     = state%alpha_Edd+alpha_rad
 state%dalpha_dr = (state%alpha_Edd+alpha_rad-alpha_old)/(1.e-10+state%r-state%r_old)

 state%c    = sqrt(state%gamma*Rg*state%Tg/state%mu)
 state%p    = state%rho*Rg*state%Tg/state%mu
#ifndef ISOTHERMAL
 state%u    = state%p/((state%gamma-1.)*state%rho)
#endif
 !state%Tg   = state%p/(state%rho*Rg)*state%mu

 state%tau_lucy = state%tau_lucy &
      - (state%r-state%r_old) * state%r0**2 &
      * (state%kappa*state%rho/state%r**2 + kappa_old*rho_old/state%r_old**2)/2.

 !update dust temperature
 if (iget_tdust == 2) then
    tau_lucy_bounded = max(0., state%tau_lucy)
    state%Tdust = Tstar * (.5*(1.-sqrt(1.-(state%r0/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
 elseif (iget_tdust == 1) then
    state%Tdust = Tstar*(state%r0/state%r)**tdust_exp
 else
    state%Tdust = state%Tg
 endif

 !apply cooling
 if (icooling > 0) then
    Q_old = state%Q
    call calc_cooling_rate(Q_code,dlnQ_dlnT,state%rho/unit_density,state%Tg,state%Tdust,&
         state%mu,state%gamma,state%K2,state%kappa)
    state%Q = Q_code*unit_ergg
    state%dQ_dr = (state%Q-Q_old)/(1.d-10+state%r-state%r_old)
 endif

 state%time = state%time + state%dt
 if (state%time + dt_next > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 else
    state%dt = dt_next
 endif
 state%nsteps = mod(state%nsteps, 65536)+1
 !if  not searching for the sonic point, keep integrating wind equation up to t = time_end
 if (state%time < state%time_end .and. .not.state%find_sonic_solution) state%spcode = 0

end subroutine wind_step

#else

!-----------------------------------------------------------------------
!
!  Integrate chemistry, cooling and hydro over one time step - OLD VERSION
!
!-----------------------------------------------------------------------
subroutine wind_step(state)
! all quantities in cgs

 use wind_equations,   only:evolve_hydro
 use ptmass_radiation, only:alpha_rad,iget_tdust,tdust_exp
 use physcon,          only:pi,Rg
 use dust_formation,   only:evolve_chem,calc_kappa_dust,calc_kappa_bowen,&
      calc_Eddington_factor,idust_opacity
 use part,             only:idK3,idmu,idgamma,idsat,idkappa
 use cooling_solver,   only:calc_cooling_rate
 use options,          only:icooling
 use units,            only:unit_ergg,unit_density
 use dim,              only:itau_alloc

 type(wind_state), intent(inout) :: state
 real :: rvT(3), dt_next, v_old, dlnQ_dlnT, Q_code
 real :: alpha_old, kappa_old, rho_old, Q_old, tau_lucy_bounded

 kappa_old  = state%kappa
 alpha_old  = state%alpha
 if (idust_opacity == 2) then
    call evolve_chem(state%dt,state%Tg,state%rho,state%JKmuS)
    state%K2        = state%JKmuS(idK2)
    state%mu        = state%JKmuS(idmu)
    state%gamma     = state%JKmuS(idgamma)
    state%kappa     = calc_kappa_dust(state%JKmuS(idK3), state%Tdust, state%rho)
    state%alpha_Edd = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, state%kappa)
 elseif (idust_opacity == 1) then
    state%kappa     = calc_kappa_bowen(state%Tdust)
    state%alpha_Edd = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, state%kappa)
 endif
 if (itau_alloc == 1) then
    state%alpha_Edd = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, state%kappa, state%tau)
 else
    state%alpha_Edd = calc_Eddington_factor(Mstar_cgs, Lstar_cgs, state%kappa)
 endif
 state%alpha = state%alpha_Edd+alpha_rad
 if (idust_opacity == 2) state%JKmuS(idalpha) = state%alpha
 if (state%time > 0.)    state%dalpha_dr      = (state%alpha-alpha_old)/(1.e-10+state%r-state%r_old)
 rvT(1) = state%r
 rvT(2) = state%v
 rvT(3) = state%Tg
 v_old  = state%v
 state%r_old = state%r
 call evolve_hydro(state%dt, rvT, state%Rstar, Mdot_cgs, state%mu, state%gamma, state%alpha, state%dalpha_dr, &
      state%Q, state%dQ_dr, state%spcode, state%dt_force, dt_next)
 state%r    = rvT(1)
 state%v    = rvT(2)
 state%a    = (state%v-v_old)/(state%dt)
 state%Tg   = rvT(3)
 state%time = state%time + state%dt
 state%dt   = dt_next
 rho_old    = state%rho

 state%c    = sqrt(state%gamma*Rg*state%Tg/state%mu)
 state%rho  = Mdot_cgs/(4.*pi*state%r**2*state%v)
 state%p    = state%rho*Rg*state%Tg/state%mu
#ifndef ISOTHERMAL
 state%u    = state%p/((state%gamma-1.)*state%rho)
#endif

 if (idust_opacity == 2) then
    state%kappa = calc_kappa_dust(state%JKmuS(idK3), state%Tdust, state%rho)
 elseif (idust_opacity == 1) then
    state%kappa = calc_kappa_bowen(state%Tdust)
 endif
 state%tau      = state%tau + state%kappa*state%rho*(1.e-10+state%r-state%r_old)
 state%tau_lucy = state%tau_lucy &
      - (state%r-state%r_old) * state%r0**2 &
      * (state%kappa*state%rho/state%r**2 + kappa_old*rho_old/state%r_old**2)/2.

 !update dust temperature
 if (iget_tdust == 2) then
    tau_lucy_bounded = max(0., state%tau_lucy)
    state%Tdust = Tstar * (.5*(1.-sqrt(1.-(state%r0/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
 elseif (iget_tdust == 1) then
    state%Tdust = Tstar*(state%r0/state%r)**tdust_exp
 else
    state%Tdust = state%Tg
 endif

 !apply cooling
 if (icooling > 0) then
    Q_old = state%Q
    call calc_cooling_rate(Q_code,dlnQ_dlnT,state%rho/unit_density,state%Tg,state%Tdust,&
         state%mu,state%gamma,state%K2,state%kappa)
    state%Q = Q_code*unit_ergg
    state%dQ_dr = (state%Q-Q_old)/(1.d-10+state%r-state%r_old)
 endif

 if (state%time_end > 0. .and. state%time + state%dt > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 endif
 state%nsteps = mod(state%nsteps, 65536)+1
 !if  not searching for the sonic point, keep integrating wind equation up to t = time_end
 if (state%time < state%time_end .and. .not.state%find_sonic_solution) state%spcode = 0

end subroutine wind_step

#endif

!-----------------------------------------------------------------------
!
!  Integrate the dusty wind equation up to sonic point
!
!-----------------------------------------------------------------------
subroutine calc_wind_profile(r0, v0, T0, time_end, state)
! all quantities in cgs
 use ptmass_radiation, only:iget_tdust
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: tau_lucy_last

 !initialize chemistry and variables
 call init_wind(r0, v0, T0, time_end, state)

 if (state%v > state%c .and. state%find_sonic_solution) then
    print *,'[wind_profile] for trans-sonic solution, the initial velocity cannot exceed the sound speed : v0=',&
         state%v,', cs=',state%c
    return
 endif

! integrate 1D wind solution with dust
 do while(state%dt > dtmin .and. state%Tg > Tdust_stop .and. .not.state%error .and. state%spcode == 0)
    tau_lucy_last = state%tau_lucy

    call wind_step(state)
    if (iget_tdust == 2 .and. (tau_lucy_last-state%tau_lucy)/tau_lucy_last < 1.e-6 .and. state%tau_lucy < .6) exit
    !if (state%r == state%r_old .or. state%tau_lucy < -1.) state%error = .true.
    if (state%r == state%r_old) state%error = .true.
    !print *,state%time,state%r,state%v/state%c,state%dt,dtmin,state%Tg,Tdust_stop,state%error,state%spcode
 enddo

end subroutine calc_wind_profile


#ifdef TO_BE_REMOVED
!-----------------------------------------------------------------------
!+
!  integrate wind equation up to time=local_time
!+
!-----------------------------------------------------------------------
subroutine wind_profile(local_time,r,v,u,rho,e,GM,T0,fdone,JKmuS)
 !in/out variables in code units (except Jstar,K)
 use units,           only:udist, utime, unit_velocity, unit_density, unit_pressure
 use dust_formation,  only:idust_opacity
 real, intent(in)  :: local_time, GM, T0
 real, intent(inout) :: r, v
 real, intent(out) :: u, rho, e, fdone
 real, optional, intent(out) :: JKmuS(:)

 type(wind_state) :: state
 real :: T

 T = T0
 r = r*udist
 v = v*unit_velocity
 if (local_time == 0.) then
    call init_wind(r, v, T, local_time, state)
    fdone = 1.d0
 else
    call calc_wind_profile(r, v, T, local_time*utime, state)
    fdone = state%time/local_time/utime
 endif
 r   = state%r/udist
 v   = state%v/unit_velocity
 rho = state%rho/unit_density
 !u = state%Tg * u_to_temperature_ratio
#ifndef ISOTHERMAL
 u  = state%p/((state%gamma-1.)*rho)/unit_pressure
#endif
 e = .5*v**2 - GM/r + state%gamma*u
 if (idust_opacity == 2) then
    JKmuS(1:n_nucleation-1) = state%JKmuS(1:n_nucleation-1)
    JKmuS(n_nucleation) = state%kappa
 endif

end subroutine wind_profile
#endif

!-----------------------------------------------------------------------
!
!  Determine the initial wind speed for a given stellar radius
!    stype = 0 : do nothing, initial velocity is set
!    stype = 1 : determine the trans-sonic solution
!
!-----------------------------------------------------------------------
subroutine get_initial_wind_speed(r0, T0, v0, rsonic, tsonic, stype)
!all quantities in cgs
 use timestep, only:tmax
 use io,       only:fatal,iverbose
 use units,    only:utime,udist
 use eos,      only:gmw,gamma
 use physcon,  only:Rg,Gg,au,years
 use ptmass_radiation, only:alpha_rad
 integer, intent(in) :: stype
 real, intent(in)    :: r0, T0
 real, intent(inout) :: v0
 real, intent(out)   :: rsonic, tsonic

 type(wind_state) :: state

 real :: v0min, v0max, v0last, vesc, cs, Rs, alpha_max, vin, gmax
 real, parameter :: v_over_cs_min = 1.d-4
 integer, parameter :: ncount_max = 20
 integer :: icount
 character(len=*), parameter :: label = 'get_initial_wind_speed'

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-alpha_rad)/r0)
 cs   = sqrt(gamma*Rg*T0/gmw)
 vin  = cs*(vesc/2./cs)**2*exp(-(vesc/cs)**2/2.+1.5)
 Rs   = Gg*Mstar_cgs*(1.-alpha_rad)/(2.*cs*cs)
 gmax = 1.-2.*cs**2*r0/(Gg*Mstar_cgs)  !Eq 3.38 Lamers & Cassinelli
 if (iverbose>0) then
    if (vesc> 1.d-50) then
       alpha_max = 1.-(2.*cs/vesc)**2
    else
       alpha_max = 0.
    endif
    print *, "[get_initial_wind_speed] searching for initial velocity."
    print *, ' * stype      = ',stype
    print *, ' * unit(au)   = ',udist/au
    print *, ' * Mstar      = ',Mstar_cgs/1.9891d33
    print *, ' * Twind      = ',T0
    print *, ' * Rstar(au)  = ',r0/1.496d13,r0/69600000000.
#ifndef ISOTHERMAL
    print *, ' * gamma      = ',gamma
#endif
    print *, ' * mu         = ',gmw
    print *, ' * cs  (km/s) = ',cs/1e5
    print *, ' * vesc(km/s) = ',vesc/1e5
    if (stype == 1) then
       print *, ' * vin/cs     = ',vin/cs
    else
       print *, ' * v0  (km/s) = ',v0/1e5
    endif
    print *, ' * alpha      = ',alpha_rad
    print *, ' * alpha_max  = ',alpha_max
    print *, ' * tend (s)   = ',tmax*utime,tmax*utime/years
 endif

 !
 ! seach for trans-sonic solution
 !
 if (stype == 1) then
    ! Find lower bound for initial velocity
    v0     = cs*0.991
    v0max  = v0
    v0min  = v_over_cs_min*cs
    icount = 0
    do while (icount < ncount_max .and. v0/cs > v_over_cs_min)
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *,'< v0/cs = ',v0/cs,', spcode = ',state%spcode
       if (state%spcode == -1) then
          v0min = v0
          exit
       else
          v0max = v0
          v0 = v0 / 2.
       endif
       icount = icount+1
    enddo
    if (iverbose>1) print *, 'Lower bound found for v0/cs :',v0min/cs,', icount=',icount
    if (icount == ncount_max .or. v0/cs < v_over_cs_min) &
         call fatal(label,'cannot find v0min, change wind_temperature or wind_injection_radius ?')
    if (v0min/cs > 0.99) call fatal(label,'supersonic wind, set sonic_type = 0 and provide wind_velocity or change alpha_rad')

    ! Find upper bound for initial velocity
    v0 = v0max
    icount = 0
    do while (icount < ncount_max)
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *,'> v0/cs = ',v0/cs,', spcode = ',state%spcode
       if (state%spcode == 1) then
          v0max = v0
          exit
       else
          v0min = max(v0min, v0)
          v0 = v0 * 1.1
          !asymptotically approaching v0max
          !v0 = min(v0, v0max*(v0/(1.+v0)))
       endif
       icount = icount+1
    enddo
    if (iverbose>1) print *, 'Upper bound found for v0/cs :',v0max/cs,', icount=',icount
    if (icount == ncount_max) call fatal(label,'cannot find v0max, change wind_temperature or wind_injection_radius ?')

    ! Find sonic point by dichotomy between v0min and v0max
    do
       v0last = v0
       v0 = (v0min+v0max)/2.
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *, '= v0/cs = ',v0/cs,', spcode = ',state%spcode
       if (state%spcode == -1) then
          v0min = v0
       elseif (state%spcode == 1) then
          v0max = v0
       else
          exit
       endif
       if (abs(v0-v0last)/v0last < 1.e-12) then
          !if (v0 == v0last) then
          exit
       endif
    enddo

    write (*,'("Sonic point properties  cs (km/s) =",f9.3,", Rs/r0 = ",f7.3,&
    &", v0/cs = ",f9.6,", ts =",f8.1,", (Rs-r0)/ts (km/s) = ",f8.4)') &
         state%v/1e5,state%r/r0,v0/state%v,state%time/utime,(state%r-state%r0)/state%time/1.e5

    rsonic = state%r
    tsonic = state%time

 else
    if (v0 >= cs) then
       print *,' supersonic wind : v0/cs = ',v0/cs
    else
       !see discussion p72/73 in Lamers & Cassinelli textbook
       if (r0 > Rs .or. alpha_rad > gmax ) then
          print *,'r0 = ',r0,', Rs = ',rs,', Gamma_max =',gmax,', alpha_rad =',alpha_rad
          print '(/," WARNING : alpha_rad > Gamma_max = ",f7.5," breeze type solution (dv/dr < 0)",/)',gmax
       endif
       print '(" sub-sonic wind : v0/cs = ",f7.5,", Gamma_max = ",f7.5)',v0/cs,gmax
    endif
    !call calc_wind_profile(r0, v0, T0, 0., state)
    rsonic = 0.!state%r
    tsonic = 0.!state%time
 endif
 !
 !store sonic point properties (location, time to reach, ...)
 !
 ! sonic(1) = state%r
 ! sonic(2) = state%v
 ! sonic(3) = state%c
 ! sonic(4) = state%time
 ! sonic(5) = state%Tg
 ! sonic(6) = state%p
 ! sonic(7) = state%alpha
 ! mdot = 4.*pi*rho*v0*ro*ro
end subroutine get_initial_wind_speed

!-----------------------------------------------------------------------
!
!  Determine the initial stellar radius for a trans-sonic solution
!
!-----------------------------------------------------------------------
subroutine get_initial_radius(r0, T0, v0, rsonic, tsonic, stype)
!all quantities in cgs
 use physcon, only:steboltz,pi
 use io,      only:iverbose
 real, parameter :: step_factor = 1.1
 integer, intent(in) :: stype
 real, intent(in) :: T0
 real, intent(inout) :: r0
 real, intent(out) :: v0, rsonic, tsonic

 real, parameter :: time_end = 1.d9
 real :: r0min, r0max, r0best, v0best, tau_lucy_best, initial_guess
 integer :: i, nstepsbest
 type(wind_state) :: state

 ! Find lower bound for initial radius
 !r0 = sqrt(Lstar_cgs/(4.*pi*steboltz*Tstar**4))
 initial_guess = r0
 r0max = r0
 v0 = 0.
 if (iverbose>0) print *, '[get_initial_radius] Searching lower bound for r0'
 do
    call get_initial_wind_speed(r0, T0, v0, rsonic, tsonic, stype)
    call calc_wind_profile(r0, v0, T0, time_end, state)
    if (iverbose>1) print *,' r0 = ', r0, 'tau_lucy = ', state%tau_lucy
    if (state%error) then
       ! something wrong happened!
       r0 = r0 / step_factor
       !elseif (r0 < initial_guess/10.) then
       !   r0min = r0
       !   print *,'radius getting too small!',r0,step_factor
       !   exit
    elseif (state%tau_lucy < 0.) then
       r0min = r0
       exit
    else
       r0max = r0
       r0 = r0 / step_factor
    endif
 enddo
 if (iverbose>1) print *, 'Lower bound found for r0 :', r0min

! Find upper bound for initial radius
 if (iverbose>1) print *, 'Searching upper bound for r0'
 r0 = r0max
 do
    call get_initial_wind_speed(r0, T0, v0, rsonic, tsonic, stype)
    call calc_wind_profile(r0, v0, T0, time_end, state)
    if (iverbose>1) print *, 'r0 = ', r0, 'tau_lucy = ', state%tau_lucy
    if (state%error) then
       ! something wrong happened!
       r0 = r0 * step_factor
    elseif (state%tau_lucy > 0.) then
       r0max = r0
       exit
    else
       r0min = max(r0, r0min)
       r0 = r0 * step_factor
    endif
 enddo
 if (iverbose>1) print *, 'Upper bound found for r0 :', r0max

 ! Find the initial radius by dichotomy between r0min and r0max
 if (iverbose>1) print *, 'Searching r0 by dichotomy'
 r0best = r0
 v0best = v0
 nstepsbest = state%nsteps
 tau_lucy_best = state%tau_lucy
 do i=1,30
    r0 = (r0min+r0max)/2.
    call get_initial_wind_speed(r0, T0, v0, rsonic, tsonic, stype)
    call calc_wind_profile(r0, v0, T0, time_end, state)
    if (iverbose>1) print *, 'r0 = ', r0, 'tau_lucy = ', state%tau_lucy
    if (abs(state%tau_lucy) < abs(tau_lucy_best)) then
       rsonic= state%r
       tsonic= state%time
       r0best = r0
       v0best = v0
       nstepsbest = state%nsteps
       tau_lucy_best = state%tau_lucy
    endif
    if (.not. state%error) then
       if (state%tau_lucy < 0.) then
          r0min = r0
       else
          r0max = r0
       endif
    else
       r0max = r0max + 12345.
    endif
    if (abs(r0min-r0max)/r0max < 1.e-5) exit
 enddo
 r0 = r0best
 v0 = v0best
 if (iverbose>0) then
    print *,'Best initial radius found: r0=',r0,' , initial_guess=',initial_guess
    print *,'with v0=', v0,' , T0=',T0,' , leading to tau_lucy=',tau_lucy_best,' at t=',time_end
 endif
end subroutine get_initial_radius

!-----------------------------------------------------------------------
!
!  Interpolate 1D wind profile
!
!-----------------------------------------------------------------------
subroutine interp_wind_profile(time, local_time, r, v, u, rho, e, GM, fdone, JKmuS)
 !in/out variables in code units (except Jstar,K)
 use units,    only:udist, utime, unit_velocity, unit_density, unit_ergg
 use dust_formation, only:idust_opacity
 use part,     only:idgamma
 use eos,      only:gamma

 real, intent(in)  :: time, local_time, GM
 !real, intent(inout) :: r, v
 real, intent(out) :: u, rho, e, r, v, fdone
 real, intent(out), optional :: JKmuS(:)

 real    :: ltime,gammai
 integer :: indx,j

 ltime = local_time*utime
 call find_nearest_index(trvurho_1D(1,:),ltime,indx)

 r   = interp_1d(ltime,trvurho_1D(1,indx),trvurho_1D(1,indx+1),trvurho_1D(2,indx),trvurho_1D(2,indx+1))/udist
 v   = interp_1d(ltime,trvurho_1D(1,indx),trvurho_1D(1,indx+1),trvurho_1D(3,indx),trvurho_1D(3,indx+1))/unit_velocity
#ifndef ISOTHERMAL
 u   = interp_1d(ltime,trvurho_1D(1,indx),trvurho_1D(1,indx+1),trvurho_1D(4,indx),trvurho_1D(4,indx+1))/unit_ergg
#else
 u   = 0.
#endif
 rho = interp_1d(ltime,trvurho_1D(1,indx),trvurho_1D(1,indx+1),trvurho_1D(5,indx),trvurho_1D(5,indx+1))/unit_density

 if (idust_opacity == 2) then
    do j=1,n_nucleation
       JKmuS(j) = interp_1d(ltime,trvurho_1D(1,indx),trvurho_1D(1,indx+1),JKmuS_1D(j,indx),JKmuS_1D(j,indx+1))
    enddo
    gammai = JKmuS(idgamma)
    gamma  = gammai
 else
    gammai = gamma
 endif

 e = .5*v**2 - GM/r + gammai*u
 if (local_time == 0.) then
    fdone = 1.d0
 else
    fdone = ltime/(local_time*utime)
 endif
end subroutine interp_wind_profile

!-----------------------------------------------------------------------
!+
!  Find index of nearest lower value in array
!+
!-----------------------------------------------------------------------
subroutine find_nearest_index(arr,val,indx)
 real, intent(in)     :: arr(:), val
 integer, intent(out) :: indx
 integer              :: istart,istop,i

 istart = 1
 istop  = size(arr)
 if (val >= arr(istop)) then
    indx = istop-1   ! -1 to avoid array index overflow
 elseif (val <= arr(istart)) then
    indx = istart
 else
    i = istart
    do while (i <= istop)
       if (arr(i)>val) then
          indx = i-1
          exit
       endif
       i = i+1
    enddo
 endif
end subroutine find_nearest_index

!-----------------------------------------------------------------------
!+
!  1D linear interpolation routine
!+
!-----------------------------------------------------------------------
real function interp_1d(x,x1,x2,y1,y2)
 real, intent(in)  :: x, x1, x2, y1, y2

 interp_1d = y1 + (x-x1)*(y2-y1)/(x2-x1)
end function interp_1d

!-----------------------------------------------------------------------
!
!  Integrate the steady wind equation and save wind profile to a file
!
!-----------------------------------------------------------------------
subroutine save_windprofile(r0, v0, T0, rout, tend, tcross, filename)
 use physcon,  only:au
 !use units,    only:utime
 use dust_formation, only:idust_opacity
 real, intent(in) :: r0, v0, T0, tend, rout
 real, intent(out) :: tcross          !time to cross the entire integration domain
 character(*), intent(in) :: filename
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 integer, parameter :: nlmax = 8192   ! maxium number of steps store in the 1D profile
 real :: time_end
 real :: r_incr,v_incr,T_incr,mu_incr,gamma_incr,r_base,v_base,T_base,mu_base,gamma_base,eps
 real, allocatable :: trvurho_temp(:,:)
 real, allocatable :: JKmuS_temp(:,:)
 type(wind_state) :: state
 integer ::iter,itermax,nwrite,writeline

 if (.not. allocated(trvurho_temp)) allocate (trvurho_temp(5,nlmax))
 if (idust_opacity == 2 .and. .not. allocated(JKmuS_temp)) allocate (JKmuS_temp(n_nucleation,nlmax))

 write (*,'("Saving 1D model to ",A)') trim(filename)
 !time_end = tmax*utime
 time_end = tend
 call init_wind(r0, v0, T0, tend, state)

 open(unit=1337,file=filename)
 call filewrite_header(1337,nwrite)
 call filewrite_state(1337,nwrite, state)

 eps       = 0.01
 iter      = 0
 itermax   = int(huge(itermax)/10) !this number is huge but may be needed for RK6 solver
 tcross    = 1.d99
 writeline = 0

 r_base     = state%r
 v_base     = state%v
 T_base     = state%Tg
 mu_base    = state%mu
 gamma_base = state%gamma

 do while(state%time < time_end .and. iter < itermax .and. state%Tg > Tdust_stop .and. writeline < nlmax)
    iter = iter+1
    call wind_step(state)

    r_incr     = state%r
    v_incr     = state%v
    T_incr     = state%Tg
    mu_incr    = state%mu
    gamma_incr = state%gamma

    if (      ( abs((r_incr     -r_base)      /r_base)      > eps ) &
         .or. ( abs((v_incr     -v_base)      /v_base)      > eps ) &
         .or. ( abs((T_incr     -T_base)      /T_base)      > eps ) &
         .or. ( abs((gamma_incr -gamma_base)  /gamma_base)  > eps ) &
         .or. ( abs((mu_incr    -mu_base)     /mu_base)     > eps ) ) then

       writeline = writeline + 1
       call filewrite_state(1337,nwrite, state)

       r_base     = state%r
       v_base     = state%v
       T_base     = state%Tg
       mu_base    = state%mu
       gamma_base = state%gamma
       trvurho_temp(:,writeline) = (/state%time,state%r,state%v,state%u,state%rho/)
       if (idust_opacity == 2) JKmuS_temp(:,writeline) = (/state%JKmuS(1:n_nucleation)/)

    endif
    if (state%r > rout) tcross = min(state%time,tcross)
 enddo
 if (state%time/time_end < .3) then
    write(*,'("[WARNING] wind integration failed : t/tend = ",f7.5,", dt/tend = ",f7.5,&
    &" Tgas = ",f6.0,", r/rout = ",f7.5," iter = ",i7,/)') &
    state%time/time_end,state%dt/time_end,state%Tg,state%r/rout,iter
 else
    print *,'integration successful, #',iter,' iterations required, rout = ',state%r/au
 endif
 close(1337)
 !stop 'save_windprofile'

 if (allocated(trvurho_1D)) deallocate(trvurho_1D)
 allocate (trvurho_1D(5, writeline))
 trvurho_1D(:,:) = trvurho_temp(:,1:writeline)
 deallocate(trvurho_temp)
 if (idust_opacity == 2) then
    if (allocated(JKmuS_1D)) deallocate(JKmuS_1D)
    allocate (JKmuS_1D(n_nucleation, writeline))
    JKmuS_1D(:,:) = JKmuS_temp(:,1:writeline)
    deallocate(JKmuS_temp)
 endif
end subroutine save_windprofile

subroutine filewrite_header(iunit,nwrite)
 use options,        only:icooling
 use dust_formation, only:idust_opacity
 integer, intent(in) :: iunit
 integer, intent(out) :: nwrite
 character (len=20):: fmt

 if (idust_opacity == 2) then
    if (icooling > 0) then
       nwrite = 22
       write(fmt,*) nwrite
       write(iunit,'('// adjustl(fmt) //'(a12))') 't','r','v','T','c','p','u','rho','alpha','a',&
            'mu','S','Jstar','K0','K1','K2','K3','tau_lucy','kappa','Tdust','gamma','Q'
    else
       nwrite = 21
       write(fmt,*) nwrite
       write(iunit,'('// adjustl(fmt) //'(a12))') 't','r','v','T','c','p','u','rho','alpha','a',&
            'mu','S','Jstar','K0','K1','K2','K3','tau_lucy','kappa','Tdust','gamma'
    endif
 else
    if (icooling > 0) then
       nwrite = 14
       write(fmt,*) nwrite
       write(iunit,'('// adjustl(fmt) //'(a12))') 't','r','v','T','c','p','u','rho','alpha','a','mu','kappa','gamma','Q'
    else
       nwrite = 13
       write(fmt,*) nwrite
       write(iunit,'('// adjustl(fmt) //'(a12))') 't','r','v','T','c','p','u','rho','alpha','a','mu','kappa','gamma'
    endif
 endif
end subroutine filewrite_header

subroutine state_to_array(state,nwrite, array)
 use dust_formation, only:idust_opacity
 use options,        only:icooling
 type(wind_state), intent(in) :: state
 integer, intent(in) :: nwrite
 real, intent(out) :: array(:)

 array(1)  = state%time
 array(2)  = state%r
 array(3)  = state%v
 array(4)  = state%Tg
 array(5)  = state%c
 array(6)  = state%p
 array(7)  = state%u
 array(8)  = state%rho
 array(9)  = state%alpha
 array(10) = state%a
 if (idust_opacity == 2) then
    array(11) = state%JKmuS(idmu)
    array(12) = state%JKmuS(idsat)
    array(13) = state%JKmuS(idJstar)
    array(14) = state%JKmuS(idK0)
    array(15) = state%JKmuS(idK1)
    array(16) = state%JKmuS(idK2)
    array(17) = state%JKmuS(idK3)
    array(18) = state%tau_lucy
    array(19) = state%kappa
    array(20) = state%Tdust
    array(21) = state%JKmuS(idgamma)
 else
    array(11) = state%mu
    array(12) = state%kappa
    array(13) = state%gamma
 endif
 if (icooling > 0) array(nwrite) = state%Q
end subroutine state_to_array

subroutine filewrite_state(iunit,nwrite, state)
 integer, intent(in) :: iunit,nwrite
 type(wind_state), intent(in) :: state

 real :: array(nwrite)
 character (len=20):: fmt

 call state_to_array(state,nwrite, array)
 write(fmt,*) nwrite
 write(iunit,'('// adjustl(fmt) //'(1x,es11.3E3:))') array(1:nwrite)
!  write(iunit, '(22(1x,es11.3E3:))') array(1:nwrite)
end subroutine filewrite_state


end module wind
