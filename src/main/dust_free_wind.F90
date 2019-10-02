!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: wind
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: cooling, eos, io, options, physcon, ptmass_radiation,
!    units, wind_equations
!+
!--------------------------------------------------------------------------

module wind
 implicit none
 public :: setup_wind
 public :: wind_state, dust_free_wind_profile,calc_wind_profile,wind_step,init_wind

 private

! Shared variables
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'dust_free_wind'

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, wind_gamma, Mdot_cgs, wind_temperature
 real :: u_to_temperature_ratio

 ! wind properties
 type wind_state
    real :: dt, time, r, v, a, time_end, Tg
    real :: mu, gamma, alpha, rho, p, c, dalpha_dr, r_old, Q, dQ_dr
    integer :: spcode, nsteps
    logical :: dt_force, error, stop_at_tend
 end type wind_state
contains

subroutine setup_wind(Mstar_in, Rstar_cg, Mdot_in, u_to_T, Twind)
 use units,        only:umass,utime
 use physcon,      only:c,solarm,years
 use eos,          only:gamma

 real, intent(in) :: Mstar_in, Rstar_cg, Mdot_in, u_to_T, Twind

 Mstar_cgs = Mstar_in*solarm
 wind_gamma = gamma
 wind_temperature = Twind
 Mdot_cgs = Mdot_in * umass/utime
 Rstar_cgs = Rstar_cg
 u_to_temperature_ratio = u_to_T

end subroutine setup_wind

!-----------------------------------------------------------------------
!
!  Initialize variables for wind integration
!
!-----------------------------------------------------------------------
subroutine init_wind(r0, v0, T0, time_end, state)
! all quantities in cgs
 use physcon,        only:pi,Rg
 use io,             only:fatal
 use eos,            only:gmw
 use ptmass_radiation, only:alpha_rad
#ifndef ISOTHERMAL
 use cooling,        only:calc_cooling_rate
 use options,        only:icooling
#endif
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: dlnQ_dlnT

 state%dt = 1000.
 if (time_end > 0.d0) then
    ! integration stops when time = time_end
    state%stop_at_tend = .true.
    state%time_end = time_end
 else
    ! integration stops once the sonic point is reached
    state%stop_at_tend = .false.
    state%time_end = -1.d0
 endif
 state%time = 0.
 state%r_old = 0.
 state%r = r0
 state%v = v0
 state%a = 0.
 state%Tg = T0
 state%alpha = alpha_rad
 state%dalpha_dr = 0.
 state%gamma = wind_gamma
 state%mu = gmw
 state%Q = 0.
 state%dQ_dr = 0.
 state%rho = Mdot_cgs/(4.*pi * state%r**2 * state%v)

#ifndef ISOTHERMAL
 if (icooling > 0) call calc_cooling_rate(state%Q, dlnQ_dlnT, state%rho, state%Tg)
#endif
 !state%p = state%rho*Rg*state%Tg/(state%mu)
 state%c = sqrt(wind_gamma*Rg*state%Tg/state%mu)
 state%dt_force = .false.
 state%spcode = 0
 state%nsteps = 1
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
#ifndef ISOTHERMAL
 use cooling,        only:calc_cooling_rate
 use options,        only:icooling
#endif

 type(wind_state), intent(inout) :: state
 real :: rvT(3), dt_next, v_old, dlnQ_dlnT

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
 state%c    = sqrt(wind_gamma*Rg*state%Tg/state%mu)
 state%rho  = Mdot_cgs/(4.*pi*state%r**2*state%v)
 !state%p = state%rho*Rg*state%Tg/state%mu

#ifndef ISOTHERMAL
 if (icooling > 0) call calc_cooling_rate(state%Q, dlnQ_dlnT, state%rho, state%Tg)
#endif
 if (state%time_end > 0. .and. state%time + state%dt > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 endif
 state%nsteps = state%nsteps + 1
 if (state%stop_at_tend .and. state%time < state%time_end) state%spcode = 0

end subroutine wind_step

!-----------------------------------------------------------------------
!
!  Integrate the dusty wind equation up to sonic point
!
!-----------------------------------------------------------------------
subroutine calc_wind_profile(r0, v0, T0, time_end, state)
! all quantities in cgs
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state

 !initialize chemistry and variables
 call init_wind(r0, v0, T0, time_end, state)

 if (state%v > state%c) then
    state%spcode = 1
    print *,'[wind_profile] Initial velocity cannot be greater than sound speed'
 endif

 ! integrate 1D wind solution with dust
 do while(state%dt > dtmin .and. state%Tg > Tdust_stop .and. .not.state%error .and. state%spcode == 0)

    call wind_step(state)

    if (state%r == state%r_old) state%error = .true.
 enddo
end subroutine calc_wind_profile

!-----------------------------------------------------------------------
!+
!  integrate wind equation up to time=local_time
!+
!-----------------------------------------------------------------------
subroutine dust_free_wind_profile(local_time, r, v, u, rho, e, GM)
 !in/out variables in code units (except Jstar,K,mu)
 use units,        only:udist,utime,unit_velocity, unit_density
 use eos,          only:gamma
 real, intent(in)    :: local_time, GM
 real, intent(inout) :: r, v
 real, intent(out)   ::  u, rho, e
 real :: T, r0, v0, local_time_cgs
 integer :: iter

 type(wind_state) :: state

 T = wind_temperature
 r0 = r
 v0 = v
 r = r*udist
 v = v*unit_velocity
 local_time_cgs = local_time * utime
 iter = 0
 call init_wind(r, v, T, local_time_cgs, state)
 do while(state%time < local_time_cgs .and. iter < 100000 .and. state%Tg > Tdust_stop)
    iter = iter+1
    call wind_step(state)
 enddo
 r = state%r/udist
 v = state%v/unit_velocity
 rho = state%rho/unit_density
 T = state%Tg
 if (gamma > 1.0001) then
    T = wind_temperature * ((r0**2 * v0) / (r**2 * v))**(gamma-1.)
    u = T * u_to_temperature_ratio
    e = .5*v**2 - GM/r + gamma*u
 else
    u = T * u_to_temperature_ratio
    !u = state%p/((gamma-1.)*rho)/unit_pressure
    e = .5*v**2 - GM/r + u
 endif
 !cs = state%c/unit_velocity
 !rho = Mdot_cgs *utime/(umass*4.*pi*r**2*v)
end subroutine dust_free_wind_profile

end module wind
