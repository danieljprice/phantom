!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: wind
!
!  DESCRIPTION: driver to integrate the wind equations
!
!  REFERENCES: Lamers & Cassinelli "Introduction to stellar winds"
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: cooling, dust_formation, eos, io, options, part, physcon,
!    ptmass_radiation, units, wind_equations
!+
!--------------------------------------------------------------------------

module wind
 implicit none
 public :: setup_wind
 public :: wind_state,wind_profile,calc_wind_profile,wind_step,init_wind

 private
 ! Shared variables
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'wind'

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, Lstar_cgs, wind_gamma, Mdot_cgs, wind_temperature, Tstar
 real :: u_to_temperature_ratio

 ! wind properties
 type wind_state
    real :: dt, time, r, v, a, time_end, Tg, Teq, mu
    real :: gamma, alpha, rho, p, c, dalpha_dr, r_old, Q, dQ_dr
    real :: tau_lucy, kappa
#ifdef NUCLEATION
    real :: JKmuS(8)
#endif
    integer :: spcode, nsteps
    logical :: dt_force, error, find_sonic_solution
 end type wind_state
contains

subroutine setup_wind(Mstar_in, Rstar_cg, Mdot_code, u_to_T, Twind)
 use units,          only:umass,utime
 use physcon,        only:c,solarm,years
 use eos,            only:gamma
#ifdef NUCLEATION
 use dust_formation, only:set_abundances
 use part,           only:n_nucleation
 use io,             only:fatal
 type(wind_state) :: state
#endif


 real, intent(in) :: Mstar_in, Rstar_cg, Mdot_code, u_to_T, Twind

 Mstar_cgs = Mstar_in*solarm
 wind_gamma = gamma
 wind_temperature = Twind
 Mdot_cgs = Mdot_code * umass/utime
 Rstar_cgs = Rstar_cg
 u_to_temperature_ratio = u_to_T

#ifdef NUCLEATION
 if (size(state%JKmuS) /= n_nucleation) call fatal(label,'wrong dimension for JKmuS')
 call set_abundances
#endif
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
 use dust_formation,   only:evolve_chem,calc_kappa_dust,kappa_dust_bowen,calc_alpha_dust,idust_opacity
 use units,            only:umass,unit_energ,utime
 use cooling,          only:calc_cooling_rate,calc_Teq
 use options,          only:icooling

 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: tau_lucy_bounded,dlnQ_dlnT

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
 state%time = 0.
 state%r_old = 0.
 state%r = r0
 state%v = v0
 state%a = 0.
 state%Tg = T0
 state%Teq = T0
 state%alpha = alpha_rad
 state%dalpha_dr = 0.
 state%gamma = wind_gamma
#ifdef NUCLEATION
 state%JKmuS = 0.
 state%jKmuS(6) = gmw
#endif
 state%tau_lucy = 2./3.
 state%mu = gmw
 state%kappa = 0.
 state%Q = 0.
 state%dQ_dr = 0.
 state%rho = Mdot_cgs/(4.*pi * state%r**2 * state%v)
 state%spcode = 0
 state%nsteps = 1

 Tstar = xyzmh_ptmass(iTeff,wind_emitting_sink)
 Lstar_cgs = xyzmh_ptmass(ilum,wind_emitting_sink)*unit_energ/utime
 Mstar_cgs = xyzmh_ptmass(4,wind_emitting_sink)*umass

#ifdef NUCLEATION
 call evolve_chem(0., T0, state%rho, state%JKmuS)
 call calc_kappa_dust(state%JKmuS(5), state%Teq, state%rho, state%kappa)
 call calc_alpha_dust(Mstar_cgs, Lstar_cgs, state%kappa, state%alpha)
 state%mu = state%jKmuS(6)
#else
 if (idust_opacity > 0) state%kappa = kappa_dust_bowen(state%Teq)
#endif

 if (icooling >0) then
    if (r0 < Rstar_cgs .and. calc_Teq) then
       call fatal(label,'cannot determine equilibrium temperature because injection_radius < Rstar')
    else
       if (calc_Teq) then
          tau_lucy_bounded = max(0., state%tau_lucy)
          state%Teq = Tstar * (.5*(1.-sqrt(1.-(Rstar_cgs/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
       endif
#ifdef NUCLEATION
       call calc_cooling_rate(state%Q,dlnQ_dlnT,state%rho,state%Tg,state%Teq,state%mu,state%JKmuS(4),state%kappa)
#else
       call calc_cooling_rate(state%Q,dlnQ_dlnT,state%rho,state%Tg,state%Teq,state%mu)
#endif
    endif
 endif
 state%p = state%rho*Rg*state%Tg/state%mu
 state%c = sqrt(wind_gamma*Rg*state%Tg/state%mu)
 state%dt_force = .false.
 state%error = .false.

end subroutine init_wind


!-----------------------------------------------------------------------
!
!  Integrate chemistry, cooling and hydro over one time step
!
!-----------------------------------------------------------------------
subroutine wind_step(state)
! all quantities in cgs

 use wind_equations, only:evolve_hydro
 use physcon,        only:pi,Rg
 use dust_formation, only:evolve_chem,calc_kappa_dust,kappa_dust_bowen,calc_alpha_dust,idust_opacity
 use cooling,        only:calc_cooling_rate,calc_Teq
 use options,        only:icooling

 type(wind_state), intent(inout) :: state
 real :: rvT(3), dt_next, v_old, dlnQ_dlnT
 real :: alpha_old, kappa_old, rho_old, Q_old, tau_lucy_bounded

#ifdef NUCLEATION
 call evolve_chem(state%dt,state%Tg,state%rho,state%JKmuS)
 alpha_old = state%alpha
 state%mu  = state%JKmus(6)
 call calc_kappa_dust(state%JKmuS(5), state%Teq, state%rho, state%kappa)
 call calc_alpha_dust(Mstar_cgs, Lstar_cgs, state%kappa, state%alpha)
 if (state%time > 0.) state%dalpha_dr = (state%alpha-alpha_old)/(1.+state%r-state%r_old)
#endif

 rvT(1) = state%r
 rvT(2) = state%v
 rvT(3) = state%Tg
 v_old = state%v
 state%r_old = state%r
 call evolve_hydro(state%dt, rvT, state%mu, state%gamma, state%alpha, state%dalpha_dr, &
      state%Q, state%dQ_dr, state%spcode, state%dt_force, dt_next)
 state%r    = rvT(1)
 state%v    = rvT(2)
 state%a    = (state%v-v_old)/(state%dt)
 state%Tg   = rvT(3)
 state%time = state%time + state%dt
 state%dt   = dt_next
 rho_old    = state%rho
 state%c    = sqrt(wind_gamma*Rg*state%Tg/state%mu)
 state%rho  = Mdot_cgs/(4.*pi*state%r**2*state%v)
 state%p    = state%rho*Rg*state%Tg/state%mu
 kappa_old = state%kappa

#ifdef NUCLEATION
 call calc_kappa_dust(state%JKmuS(5), state%Teq, state%rho, state%kappa)
#else
 if (idust_opacity > 0) state%kappa = kappa_dust_bowen(state%Teq)
#endif
 state%tau_lucy = state%tau_lucy &
      - (state%r-state%r_old) * Rstar_cgs**2 &
      * (state%kappa*state%rho/state%r**2 + kappa_old*rho_old/state%r_old**2)/2.

 if (icooling > 0) then
    Q_old = state%Q
    if (calc_Teq) then
       tau_lucy_bounded = max(0., state%tau_lucy)
       state%Teq = Tstar * (.5*(1.-sqrt(1.-(Rstar_cgs/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
    endif
#ifdef NUCLEATION
    call calc_cooling_rate(state%Q,dlnQ_dlnT,state%rho,state%Tg,state%Teq,state%JKmuS(6),state%JKmuS(4),state%kappa)
#else
    call calc_cooling_rate(state%Q,dlnQ_dlnT,state%rho,state%Tg,state%Teq,state%mu)
#endif
    if (state%time > 0. .and. state%r /= state%r_old) state%dQ_dr = (state%Q-Q_old)/(1.d-10+state%r-state%r_old)
 else
    !if no cooling or impose temperature profile, assume Tdust = Tgas
    state%Teq = state%Tg
 endif
 if (state%time_end > 0. .and. state%time + state%dt > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 endif
 state%nsteps = state%nsteps + 1
 !if  not searching for the sonic point, keep integrating wind equation up to t = time_end
 if (state%time < state%time_end .and. .not.state%find_sonic_solution) state%spcode = 0

end subroutine wind_step

!-----------------------------------------------------------------------
!
!  Integrate the dusty wind equation up to sonic point
!
!-----------------------------------------------------------------------
subroutine calc_wind_profile(r0, v0, T0, time_end, state)
! all quantities in cgs
 use dust_formation, only:idust_opacity
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: tau_lucy_last
 logical :: tau_test

 !initialize chemistry and variables
 call init_wind(r0, v0, T0, time_end, state)

 if (state%v > state%c .and. state%find_sonic_solution) then
    print *,'[wind_profile] Initial velocity cannot be greater than sound speed'
    return
 endif

 ! integrate 1D wind solution with dust
 do while(state%dt > dtmin .and. state%Tg > Tdust_stop .and. .not.state%error .and. state%spcode == 0)
    tau_lucy_last = state%tau_lucy

    call wind_step(state)

    tau_test = idust_opacity > 0 .and. (tau_lucy_last-state%tau_lucy)/tau_lucy_last < 1.e-6 .and. state%tau_lucy < .6
    if (state%r == state%r_old .or. state%tau_lucy < -1. .or. tau_test) state%error = .true.
 enddo
end subroutine calc_wind_profile

!-----------------------------------------------------------------------
!+
!  integrate wind equation up to time=local_time
!+
!-----------------------------------------------------------------------
subroutine wind_profile(local_time,r,v,u,rho,e,GM,T0,JKmuS)
 !in/out variables in code units (except Jstar,K,mu)
 use units,        only:udist, utime, unit_velocity, unit_density!, unit_pressure
 use eos,          only:gamma
 real, intent(in)  :: local_time, GM, T0
 real, intent(inout) :: r, v
 real, intent(out) :: u, rho, e
 real, optional, intent(out) :: JKmuS(:)

 type(wind_state) :: state
 real :: T

 T = T0
 r = r*udist
 v = v*unit_velocity
 if (local_time == 0.) then
    call init_wind(r, v, T, local_time, state)
 else
    call calc_wind_profile(r, v, T, local_time*utime, state)
 endif
 r = state%r/udist
 v = state%v/unit_velocity
 rho = state%rho/unit_density
 u = state%Tg * u_to_temperature_ratio
 !u = state%p/((gamma-1.)*rho)/unit_pressure
 e = .5*v**2 - GM/r + gamma*u
#ifdef NUCLEATION
 JKmuS = state%JKmuS
#endif
 !cs = state%c/unit_velocity
end subroutine wind_profile



! !-----------------------------------------------------------------------
! !
! !  Determine the initial stellar radius for a trans-sonic solution
! !
! !-----------------------------------------------------------------------
! subroutine profile_findr0(time_end, T0, r0, v0, step_factor, nsteps, sonic)
! !all quantities in cgs
!  use physcon, only:steboltz,pi
!  real, intent(in) :: time_end, T0, step_factor
!  real, intent(inout) :: r0
!  real, intent(out) :: v0, sonic(8)
!  integer, intent(out) :: nsteps
!  logical :: verbose = .false.

!  real :: r0min, r0max, r0best, v0best, tau_lucy_best, initial_guess
!  integer :: i, nstepsbest
!  type(wind_state) :: state

! ! Find lower bound for initial radius
!  !r0 = sqrt(Lstar_cgs/(4.*pi*steboltz*Tstar**4))
!  initial_guess = r0
!  r0max = r0
!  if (verbose) print *, '[profile_findr0] Searching lower bound for r0'
!  do
!     Rstar_cgs = r0
!     call get_initial_wind_speed(r0, T0, v0, sonic, 1)
!     call calc_wind_profile(r0, v0, T0, time_end, state)
!     if (verbose) print *, state%error,' r0 = ', r0, 'tau_lucy = ', state%tau_lucy
!     if (state%error) then
! ! something wrong happened!
!        r0 = r0 / step_factor
!     elseif (r0 < initial_guess/10.) then
!        r0min = r0
!        print *,'radius getting too small!',r0,step_factor
!        exit
!     elseif (state%tau_lucy < 0.) then
!        r0min = r0
!        exit
!     else
!        r0max = r0
!        r0 = r0 / step_factor
!     endif
!  enddo
!  if (verbose) print *, 'Lower bound found for r0 :', r0min

! ! Find upper bound for initial radius
!  r0 = r0max
!  if (verbose) print *, 'Searching upper bound for r0'
!  do
!     Rstar_cgs = r0
!     call get_initial_wind_speed(r0, T0, v0, sonic, 1)
!     call calc_wind_profile(r0, v0, T0, time_end, state)
!     if (verbose) print *, 'r0 = ', r0, 'tau_lucy = ', state%tau_lucy
!     if (state%error) then
!        ! something wrong happened!
!        r0 = r0 * step_factor
!     elseif (state%tau_lucy > 0.) then
!        r0max = r0
!        exit
!     else
!        r0min = max(r0, r0min)
!        r0 = r0 * step_factor
!     endif
!  enddo
!  if (verbose) print *, 'Upper bound found for r0 :', r0max

! ! Find the initial radius by dichotomy between r0min and r0max
!  if (verbose) print *, 'Searching r0 by dichotomy'
!  r0best = r0
!  v0best = v0
!  nstepsbest = state%nsteps
!  tau_lucy_best = state%tau_lucy
!  do i=1,50
!     r0 = (r0min+r0max)/2.
!     Rstar_cgs = r0
!     call get_initial_wind_speed(r0, T0, v0, sonic, 1)
!     call calc_wind_profile(r0, v0, T0, time_end, state)
!     if (verbose) print *, 'r0 = ', r0, 'tau_lucy = ', state%tau_lucy
!     if (abs(state%tau_lucy) < abs(tau_lucy_best)) then
!        r0best = r0
!        v0best = v0
!        nstepsbest = state%nsteps
!        tau_lucy_best = state%tau_lucy
!     endif
!     if (.not. state%error) then
!        if (state%tau_lucy < 0.) then
!           r0min = r0
!        else
!           r0max = r0
!        endif
!     else
!        r0max = r0max + 12345.
!     endif
!     if (abs(r0min-r0max)/r0max < 1.e-5) exit
!  enddo
!  r0 = r0best
!  v0 = v0best
!  nsteps = nstepsbest
!  if (verbose) then
!     print *, 'Best initial radius found: r0 = ', r0,' , initial_guess = ',initial_guess
!     print *, '  with v0 = ', v0,' , T0 = ',T0
!     print *, '  it gives tau_lucy = ', tau_lucy_best
!     print *, '  at t = ', time_end
!  endif
! end subroutine profile_findr0

end module wind
