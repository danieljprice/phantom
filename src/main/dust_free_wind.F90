!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: dusty_wind
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
!  DEPENDENCIES: dust_formation, eos, io, physcon, timestep, units,
!    wind_profile
!+
!--------------------------------------------------------------------------

module dust_free_wind
 implicit none
 public :: setup_wind
 public :: get_initial_wind_speed!, profile_findr0
 public :: stationary_wind_profile

 private

! Shared variables
 real, parameter :: Tend = 50.d0 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'dust_free_wind'


! Wind properties
 real :: Mstar_cgs, Lstar_cgs, Tstar, Rstar_cgs, wind_gamma, wind_mass_rate, Mdot_cgs
 real :: Cprime, u_to_temperature_ratio, expT, wind_alpha, wind_temperature,Rstar
 integer :: wind_type

! State of the wind
 type wind_state
    real :: dt, t, r, v, a, time_end, Tg, S
    real :: mu, gamma, alpha, rho, p, c, dalpha_dr, r_old, Q, dQ_dr
    integer :: spcode, nsteps
    logical :: dt_force, error, full_integration
 end type wind_state
contains

subroutine setup_wind(Mstar_in, Lstar_in, Tstar_in, Rstar_in, Cprime_cgs, &
     expT_in, Mdot_in, CO_ratio, u_to_T, alpha_in, Twind, wind_type_in)
 use units,   only:udist,umass,utime
 use physcon, only:c, solarm, years
 use eos,     only:gamma
#ifndef ISOTHERMAL
 use dust_formation, only: set_cooling
#endif
 real, intent(in) :: Mstar_in, Lstar_in, Tstar_in, Rstar_in, Cprime_cgs, expT_in, Mdot_in, CO_ratio, &
          u_to_T,alpha_in, Twind
 integer, intent(in) :: wind_type_in
#ifndef ISOTHERMAL
 logical :: cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan
#endif

 Mstar_cgs = Mstar_in*solarm
 Lstar_cgs = Lstar_in
 Tstar = Tstar_in
 wind_type = wind_type_in
 wind_alpha = alpha_in
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
 Mdot_cgs = Mdot_in * umass/utime
 Rstar_cgs = Rstar_in
 Rstar = Rstar_in/udist
 u_to_temperature_ratio = u_to_T

#ifndef ISOTHERMAL
 cool_radiation_H0 = .false.
 if (wind_type == 2) then
    cool_relaxation_Bowen = .false.
 else
    cool_relaxation_Bowen = .true.
 endif
 cool_collisions_dust = .false.
 cool_relaxation_Stefan = .false.
 call set_cooling(cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan)
#endif
end subroutine setup_wind

!-----------------------------------------------------------------------
!
!  Initialize wind variables
!
!-----------------------------------------------------------------------
subroutine init_wind(r0, v0, T0, time_end, state)
! all quantities in cgs
 use physcon, only:pi,kboltz,atomic_mass_unit
 use io,      only:fatal
 !use eos,     only:gmw
 !use dust_formation, only: calc_cooling_rate
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(inout) :: state

 state%dt = 1000.
 if (time_end > 0.d0) then
    state%full_integration = .true.
    state%time_end = time_end
 else
    state%full_integration = .false.
    state%time_end = 1.d99
 endif
 state%t = 0.
 state%r_old = 0.
 state%r = r0
 state%v = v0
 state%a = 0.
 state%Tg = T0
 state%alpha = wind_alpha
 state%dalpha_dr = 0.
 !state%mu = gmw
 !state%gamma = wind_gamma
 state%Q = 0.
 state%dQ_dr = 0.
 state%rho = Mdot_cgs/(4.*pi * state%r**2 * state%v)

 ! if (wind_type == 4) then
 !    !only H0 cooling should be considered
 !    call calc_cooling_rate(state%rho, state%Tg, Teq, wind_gamma, state%mu, Cprime, &
 !            state%K(2)/(state%r**2*state%v), state%kappa_ross, state%Q)
 ! endif
 state%p = state%rho*kboltz*state%Tg/(state%mu*atomic_mass_unit)
 state%c = sqrt(wind_gamma*kboltz*state%Tg/(state%mu*atomic_mass_unit))
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

 use wind_profile,   only:evolve_hydro
 use physcon,        only:pi
 type(wind_state), intent(inout) :: state
 real, parameter :: max_dt = 1.d9
 real :: rvT(3), dt_next, v_old, dt_max

 if (state%time_end > 0.) then
    dt_max = min(state%time_end - state%t, max_dt)
 else
    dt_max = max_dt
 endif
 state%r_old = state%r
 rvT(1) = state%r
 rvT(2) = state%v
 rvT(3) = state%Tg
 v_old = state%v
 call evolve_hydro(state%dt, rvT, state%mu, state%gamma, state%alpha, state%dalpha_dr, &
      state%Q, state%dQ_dr, state%spcode, state%dt_force, dt_max, dt_next)
 state%r = rvT(1)
 state%v = rvT(2)
 state%a = (state%v-v_old)/(state%dt)
 state%Tg = rvT(3)
 state%t = state%t + state%dt
 state%dt = dt_next
 state%rho = Mdot_cgs/(4.*pi*state%r**2*state%v)
 !state%p = state%rho*kboltz*state%Tg/(state%mu*atomic_mass_unit)
 !state%c = sqrt(wind_gamma*kboltz*state%Tg/(state%mu*atomic_mass_unit))
 if (state%time_end > 0. .and. state%t + state%dt > state%time_end) then
    state%dt = state%time_end-state%t
    state%dt_force = .true.
 endif
 state%nsteps = state%nsteps + 1
 if (state%full_integration) state%spcode = 0

end subroutine wind_step

!-----------------------------------------------------------------------
!
!  Integrate the dusty wind equation up to t = t_end
!
!-----------------------------------------------------------------------
subroutine calc_wind_profile(r0, v0, T0, time_end, state)
! all quantities in cgs
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(inout) :: state

!compute chemistry and initialize variables
 call init_wind(r0, v0, T0, time_end, state)

 if (state%v > state%c) then
    state%spcode = 1
    !print *, 'Initial velocity cannot be greater than sound speed'
 endif

!compute 1D wind solution with dust
 do while(state%t < state%time_end .and. state%dt> dtmin .and. state%Tg > Tend .and. .not.state%error .and. state%spcode == 0)

    call wind_step(state)

    if (state%r == state%r_old) state%error = .true.
 enddo
end subroutine calc_wind_profile

!-----------------------------------------------------------------------
!+
!  dusty wind model
!+
!-----------------------------------------------------------------------

!ubroutine dust_free_wind_profile(local_time, r, v, u, rho, e, GM, gamma, mu)
subroutine stationary_wind_profile(local_time, r, v, u, rho, e, GM, gamma, mu)
 !in/out variables in code units (except Jstar,K,mu)
 use units,        only:udist, utime, unit_velocity, unit_density!, unit_pressure
 real, intent(in)  :: local_time, GM, gamma, mu
 real, intent(inout) :: r, v
 real, intent(out) ::  u, rho, e
 real :: dt, T, r0, v0

 type(wind_state) :: state

 T = wind_temperature
 r0 = r
 v0 = v
 r = r*udist
 v = v*unit_velocity
 state%mu = mu
 state%gamma = gamma
 call calc_wind_profile(r, v, T, local_time*utime, state)
 r = state%r/udist
 v = state%v/unit_velocity
 rho = state%rho/(unit_density)
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
end subroutine


!-----------------------------------------------------------------------
!+
!  stationary wind - solution when input velocity is defined
!+
!-----------------------------------------------------------------------
subroutine old_stationary_wind_profile(local_time, r, v, u, rho, e, GM, gamma, mu)
 use physcon,     only: pi
 use units,       only:udist,unit_velocity
 use wind_profile, only:RK4_step_dr
 real, intent(in)  :: local_time, GM, gamma, mu
 real, intent(inout) :: r, v
 real, intent(out) ::  u, rho, e
 real :: dt, T, r0, v0, rvT(3), new_rvT(3), err, Q, dQ_dr, dalpha_dr, numerator, denominator
 integer, parameter :: N = 10000
 integer :: i

 dt = local_time / N
 r0 = r
 v0 = v
 Q = 0.
 dQ_dr = 0.
 rvT(1) = r*udist
 rvT(2) = v*unit_velocity
 rvT(3) = wind_temperature
 do i=1,N
  call RK4_step_dr(dt, rvT, mu, gamma, wind_alpha, dalpha_dr, Q, dQ_dr, err, new_rvT, numerator, denominator)
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

end subroutine

!-----------------------------------------------------------------------
!
!  Integrate the steady wind equation and save variables to file
!
!-----------------------------------------------------------------------
subroutine save_windprofile(r0, T0, v0, time_end, dt_write, filename)
 use units,    only:utime
 real, intent(in) :: r0, v0, T0, time_end, dt_write
 character(*), intent(in) :: filename

 type(wind_state) :: state
 logical :: written
 integer :: i

 call init_wind(r0, v0, T0, time_end, state)
 open(unit=1337,file=filename)
 call filewrite_header(1337)
 call filewrite_state(1337, state)

 i = 0
 do while(state%t < state%time_end .and. state%dt > dtmin .and. state%Tg > Tend)
    call wind_step(state)
    written = .false.
    if (dt_write == 0. .or. mod(i,1000)==0) then
       call filewrite_state(1337, state)
       written = .true.
    else
       if (int((state%t-state%dt)/dt_write) < int(state%t/dt_write)) then
          call filewrite_state(1337, state)
          written = .true.
       endif
    endif
    i = i+1
 enddo
 if (.not. written) then
    call filewrite_state(1337, state) ! write last state
    if (dt_write > 0.) print *, 't/tend = ', state%t/state%time_end, ' (last step) t =',state%time_end/utime
 endif
 close(1337)
end subroutine save_windprofile


subroutine filewrite_header(iunit)
 integer, intent(in) :: iunit

 if (wind_type == 4) then
    write(iunit, '(a)') '# t  r  v  T  c  p  rho  mu  alpha  a  Q'
 else
    write(iunit, '(a)') '# t  r  v  T  c  p  rho  mu  alpha  a'
 endif
end subroutine filewrite_header

subroutine state_to_array(state, array)
 type(wind_state), intent(in) :: state
 real, intent(out) :: array(11)
 real :: f

 f = state%r**2 * state%v
 array(1) = state%t
 array(2) = state%r
 array(3) = state%v
 array(4) = state%Tg
 array(5) = state%c
 array(6) = state%p
 array(7) = state%rho
 array(8) = state%mu
 array(9) = state%alpha
 array(10) = state%a
 if (wind_type == 4) array(11) = state%Q
end subroutine state_to_array

subroutine filewrite_state(iunit, state)
 integer, intent(in) :: iunit
 type(wind_state), intent(in) :: state

 real :: array(11)
 call state_to_array(state, array)
 if (wind_type == 4) then
    write(iunit, '(20E20.10E3)') array(1:11)
 else
    write(iunit, '(19E20.10E3)') array(1:10)
 endif
end subroutine filewrite_state


!-----------------------------------------------------------------------
!
!  Determine the initial wind speed for a trans-sonic solution
!
!-----------------------------------------------------------------------
subroutine get_initial_wind_speed(r0, T0, v0, sonic, verbose)
!all quantities in cgs
 use io,       only:fatal
 use timestep, only:tmax
 use units,    only:utime
 use eos,      only:gmw
 use physcon,  only:kboltz,atomic_mass_unit,Gg
 real, intent(in) :: r0, T0
 real, intent(out) :: v0, sonic(:)
 logical, intent(in) :: verbose


 type(wind_state) :: state

 real :: v0min, v0max, v0last, vesc, cs
 integer, parameter :: ncount_max = 20
 integer :: icount
 character(len=*), parameter :: label = 'get_initial_wind_speed'

 vesc = sqrt(2.*Gg*Mstar_cgs/r0)
 cs = sqrt(wind_gamma*kboltz/atomic_mass_unit*T0/gmw)
 v0 = cs*(vesc/2./cs)**2*exp(-(vesc/cs)**2/2.+1.5)
 if (verbose) then
    print *, "[get_initial_wind_speed] Looking for initial velocity."
#ifndef ISOTHERMAL
    print *, ' * Tstar  = ', Tstar
    print *, ' * Mstar  = ', Mstar_cgs/1.9891d33
    print *, ' * Lstar  = ', Lstar_cgs/3.9d33
    print *, ' * Rstar  = ', Rstar_cgs/1.496d13
    print *, ' * Mdot   = ', Mdot_cgs/6.30303620274d25
    print *, ' * mu     = ', gmw
    print *, ' * T0     = ', T0
    print *, ' * r0(au) = ', r0/1.496d13
    print *, ' * Cprime = ', Cprime
    print *, ' * gamma  = ', wind_gamma
    if (wind_type == 2) print *, ' * expT   = ', expT
#else
    print *, ' * Mstar (Mo) = ', Mstar_cgs/1.9891d33
    print *, ' * Twind      = ', T0
    print *, ' * Rstar (Ro) = ', r0/69600000000.,Rstar_cgs/69600000000.
    print *, ' * mu         = ', gmw
    print *, ' * cs  (km/s) = ', cs/1e5
    print *, ' * vesc(km/s) = ', vesc/1e5
    print *, ' * v0  (km/s) = ', v0/1e5
#endif
 endif

! Find lower bound for initial velocity
 v0 = v0*100.
 v0max = v0
 icount = 0
 state%mu = gmw
 state%gamma = wind_gamma
 do while (icount < ncount_max)
    call calc_wind_profile(r0, v0, T0, 0., state)
    if (verbose) print *,' v0 = ', v0,state%r,state%v,state%c,state%t,icount,state%spcode
    if (state%spcode == -1) then
       v0min = v0
       exit
    else
       v0max = v0
       v0 = v0 / 2.
    endif
    icount = icount+1
 enddo
 if (icount == ncount_max) call fatal(label,'cannot find v0min, change wind_mass_rate, wind_temperature or v0 ?')
 if (verbose) print *, 'Lower bound found for v0 :', v0min

! Find upper bound for initial velocity
 v0 = v0max
 icount = 0
 do while (icount < ncount_max)
    if (verbose) print *, ' v0 = ', v0
    call calc_wind_profile(r0, v0, T0, 0., state)
    if (state%spcode == 1) then
       v0max = v0
       exit
    else
       v0min = max(v0min, v0)
       v0 = v0 * 1.1
    endif
    icount = icount+1
 enddo
 if (icount == ncount_max)  call fatal(label,'cannot find v0max, change wind_mass_rate or wind_temperature ?')
 if (verbose) print *, 'Upper bound found for v0 :', v0max

! Find sonic point by dichotomy between v0min and v0max
 do
    v0last = v0
    v0 = (v0min+v0max)/2.
    if (verbose) print *, 'v0 = ', v0
    call calc_wind_profile(r0, v0, T0, 0., state)
    if (state%spcode == -1) then
       v0min = v0
    elseif (state%spcode == 1) then
       v0max = v0
    else
       exit
    endif
    if (abs(v0-v0last)/v0last < 1.e-5) then
       exit
    endif
 enddo
 sonic(1) = state%r
 sonic(2) = state%v
 sonic(3) = state%c
 sonic(4) = state%t
 sonic(5) = state%Tg
 sonic(6) = state%p
 sonic(7) = state%S
 sonic(8) = state%alpha
!if (verbose) then
 !mdot = 4.*pi*rho*v0*ro*ro
 write (*,'("Initial conditions     v0 (km/s) =",f9.1," r0/R* = ",f7.3)') v0/1e5,r0/Rstar_cgs
 write (*,'("Sonic point properties vs (km/s) =",f9.1," Rs/R* = ",f7.3," Ts =",f7.1," S =",f6.1," alpha =",f5.3,/)')&
            sonic(2)/1e5,sonic(1)/Rstar_cgs,sonic(5),sonic(7),sonic(8)
 !endif
!save 1D initial profile for comparison
 call save_windprofile(R0, T0, v0, tmax*utime, sonic(4)*utime/10., 'gailstatwind1D.dat')

end subroutine get_initial_wind_speed


end module dust_free_wind
