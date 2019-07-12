!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: dust_free_wind
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
!  DEPENDENCIES: dust_physics, eos, io, physcon, timestep, units,
!    wind_profile
!+
!--------------------------------------------------------------------------

module dust_free_wind
 implicit none
 public :: setup_wind
 public :: get_initial_wind_speed!, profile_findr0
 public :: dust_free_wind_profile

 private

! Shared variables
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'dust_free_wind'

! Wind properties
 real :: Mstar_cgs, Rstar_cgs, wind_gamma, Mdot_cgs, wind_temperature
 real :: u_to_temperature_ratio, wind_alpha
 integer :: wind_type

! State of the wind
 type wind_state
    real :: dt, time, r, v, a, time_end, Tg
    real :: mu, gamma, alpha, rho, p, c, dalpha_dr, r_old, Q, dQ_dr
    integer :: spcode, nsteps
    logical :: dt_force, error, stop_at_tmax
 end type wind_state
contains

subroutine setup_wind(Mstar_in, Rstar_cg, Cprime_cgs, Mdot_in, u_to_T, alpha_in, Twind, wind_type_in, wind_cooling)
 use units,        only:umass,utime
 use physcon,      only:c, solarm, years
 use eos,          only:gamma
#ifndef ISOTHERMAL
 use dust_physics, only: init_cooling,set_abundances
#endif

 real, intent(in) :: Mstar_in, Rstar_cg, Cprime_cgs, Mdot_in,u_to_T,alpha_in, Twind
 integer, intent(in) :: wind_type_in, wind_cooling

 Mstar_cgs = Mstar_in*solarm
 wind_type = wind_type_in
 wind_alpha = alpha_in
 if (wind_type == 2) then
    wind_gamma = 1.
 else if (wind_type == 1) then
    wind_gamma = 1.
 else
    wind_gamma = gamma
 endif
 wind_temperature = Twind
 Mdot_cgs = Mdot_in * umass/utime
 Rstar_cgs = Rstar_cg
 u_to_temperature_ratio = u_to_T

#ifndef ISOTHERMAL
 call set_abundances(0.4d0) !needed to initialize mass_per_H
 call init_cooling(wind_cooling,Cprime_cgs)
#endif

end subroutine setup_wind

!-----------------------------------------------------------------------
!
!  Initialize wind variables
!
!-----------------------------------------------------------------------
subroutine init_wind(r0, v0, T0, time_end, state)
! all quantities in cgs
 use physcon,      only:pi,Rg
 use io,           only:fatal
 !use eos,         only:gmw
#ifndef ISOTHERMAL
 use dust_physics, only: calc_cooling_rate
#endif
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(inout) :: state
 real :: dlnq_dlnT

 state%dt = 1000.
 if (time_end > 0.d0) then
    ! compute the full wind profile
    state%stop_at_tmax = .true.
    state%time_end = time_end
 else
    ! integration to find sonic point
    state%stop_at_tmax = .false.
    state%time_end = -1.d0
 endif
 state%time = 0.
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

#ifndef ISOTHERMAL
 if (wind_type == 4) call calc_cooling_rate(state%Q, dlnQ_dlnT, state%rho, state%Tg)
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
 use dust_physics,   only:calc_cooling_rate
#endif

 type(wind_state), intent(inout) :: state
 real :: rvT(3), dt_next, v_old, dlnq_dlnT

 rvT(1) = state%r
 rvT(2) = state%v
 rvT(3) = state%Tg
 v_old = state%v
 state%r_old = state%r
 call evolve_hydro(state%dt, rvT, state%mu, state%gamma, state%alpha, state%dalpha_dr, &
      state%Q, state%dQ_dr, state%spcode, state%dt_force, dt_next)
 state%r = rvT(1)
 state%v = rvT(2)
 state%a = (state%v-v_old)/(state%dt)
 state%Tg = rvT(3)
 state%time = state%time + state%dt
 state%dt = dt_next
 state%rho = Mdot_cgs/(4.*pi*state%r**2*state%v)
 !state%p = state%rho*Rg*state%Tg/(state%mu)
 state%c = sqrt(wind_gamma*state%Tg/state%mu)

#ifndef ISOTHERMAL
 if (wind_type == 4) call calc_cooling_rate(state%Q, dlnq_dlnT, state%rho, state%Tg)
#endif
 if (state%time_end > 0. .and. state%time + state%dt > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 endif
 state%nsteps = state%nsteps + 1
 if (state%stop_at_tmax .and. state%time < state%time_end) state%spcode = 0

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

!integrate 1D wind solution with dust
 do while(state%dt> dtmin .and. state%Tg > Tdust_stop .and. .not.state%error .and. state%spcode == 0)

    call wind_step(state)

    if (state%r == state%r_old) state%error = .true.
 enddo

end subroutine calc_wind_profile

!-----------------------------------------------------------------------
!+
!  dusty wind model
!+
!-----------------------------------------------------------------------
subroutine dust_free_wind_profile(local_time, r, v, u, rho, e, GM, gamma, mu)
 !in/out variables in code units (except Jstar,K,mu)
 use units,        only:udist, utime, unit_velocity, unit_density!, unit_pressure
 real, intent(in)  :: local_time, GM, gamma, mu
 real, intent(inout) :: r, v
 real, intent(out) ::  u, rho, e
 real :: T, r0, v0

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
 !rho = Mdot_cgs *utime/(umass*4.*pi*r**2*v)
end subroutine dust_free_wind_profile

!-----------------------------------------------------------------------
!
!  Integrate the steady wind equation and save variables to file
!
!-----------------------------------------------------------------------
subroutine save_windprofile(r0, T0, v0, tsonic, filename)
 use units,    only:utime
 use timestep, only:tmax
 real, intent(in) :: r0, v0, T0, tsonic
 character(*), intent(in) :: filename

 real :: dt_print,time_end
 type(wind_state) :: state
 logical :: written
 integer :: n

 time_end = tmax*utime*2.
 call init_wind(r0, v0, T0, time_end, state)
 open(unit=1337,file=filename)
 call filewrite_header(1337)
 call filewrite_state(1337, state)

 n = 1
 dt_print = min(tsonic/10.,state%time_end/100.)
 do while(state%time < state%time_end .and. state%dt > dtmin .and. state%Tg > Tdust_stop)
    call wind_step(state)
    written = .false.
    if (state%time > n*dt_print) then
       n = floor(state%time/dt_print)+1
       call filewrite_state(1337, state)
       written = .true.
    endif
 enddo
 if (.not. written) call filewrite_state(1337, state) ! write last state
 write(*,'("t/tend = ",f7.5," (last step) t =",f6.1)') state%time/state%time_end,state%time_end/utime
 close(1337)
end subroutine save_windprofile


subroutine filewrite_header(iunit)
 integer, intent(in) :: iunit

 if (wind_type == 4) then
    write(iunit,'("#",11x,a1,10(a20))') 't','r','v','T','c','p','rho','mu','alpha','a','Q'
 else
    write(iunit,'("#",11x,a1,9(a20))') 't','r','v','T','c','p','rho','mu','alpha','a'
 endif
end subroutine filewrite_header

subroutine state_to_array(state, array)
 type(wind_state), intent(in) :: state
 real, intent(out) :: array(11)

 array(1) = state%time
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
subroutine get_initial_wind_speed(r0, T0, v0, sonic)
!all quantities in cgs
 use timestep, only:tmax
 use io,       only:fatal,iverbose
 use units,    only:utime,udist
 use eos,      only:gmw
 use physcon,  only:Rg,Gg,au,years
 real, intent(in) :: r0, T0
 real, intent(out) :: v0, sonic(:)

 type(wind_state) :: state

 real :: v0min, v0max, v0last, vesc, cs, Rs, alpha_max
 integer, parameter :: ncount_max = 20
 integer :: icount
 character(len=*), parameter :: label = 'get_initial_wind_speed'

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-wind_alpha)/r0)
 cs = sqrt(wind_gamma*Rg*T0/gmw)
 v0 = cs*(vesc/2./cs)**2*exp(-(vesc/cs)**2/2.+1.5)
 Rs = Gg*Mstar_cgs*(1.-wind_alpha)/(2.*cs*cs)
 if (iverbose>0) then
    if (vesc> 1.d-100) then
       alpha_max = 1.-(2.*cs/vesc)**2
    else
       alpha_max = 0.
    endif
    print *, "[get_initial_wind_speed] Looking for initial velocity."
    print *, ' * unit(au)   = ',udist/au
    print *, ' * Mstar      = ',Mstar_cgs/1.9891d33
    print *, ' * Twind      = ',T0
#ifndef ISOTHERMAL
    print *, ' * Rstar(au)  = ',Rstar_cgs/1.496d13
    print *, ' * Mdot       = ',Mdot_cgs/6.30303620274d25
    print *, ' * r0(au)     = ',r0/1.496d13,r0/69600000000.
    print *, ' * gamma      = ',wind_gamma
#else
    print *, ' * Rstar (Ro) = ',r0/69600000000.,Rstar_cgs/69600000000.
#endif
    print *, ' * mu         = ',gmw
    print *, ' * cs  (km/s) = ',cs/1e5
    print *, ' * vesc(km/s) = ',vesc/1e5
    print *, ' * v0  (km/s) = ',v0/1e5
    print *, ' * alpha      = ',wind_alpha
    print *, ' * alpha_max  = ',alpha_max
    print *, ' * tend (s)   = ',tmax*utime,tmax*utime/years
 endif
 write (*,'("Computing 1D model with v0 (km/s) =",f9.3,"  r0/R* = ",f7.3)') cs/1e5,r0/Rstar_cgs

! Find lower bound for initial velocity
 v0 = cs
 v0max = v0
 icount = 0
 state%mu = gmw
 state%gamma = wind_gamma
 do while (icount < ncount_max)
    call calc_wind_profile(r0, v0, T0, 0., state)
    if (iverbose>1) print *,' v0 = ', v0,state%r,state%v,state%c,state%time,icount,state%spcode
    if (state%spcode == -1) then
       v0min = v0
       exit
    else
       v0max = v0
       v0 = v0 / 2.
    endif
    icount = icount+1
 enddo
 if (icount == ncount_max) call fatal(label,'cannot find v0min, change wind_temperature or wind_injection_radius ?')
 if (iverbose>1) print *, 'Lower bound found for v0 :', v0min

! Find upper bound for initial velocity
 v0 = v0max
 icount = 0
 do while (icount < ncount_max)
    if (iverbose>1) print *, ' v0 = ', v0
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
 if (icount == ncount_max)  call fatal(label,'cannot find v0max, change wind_temperature or wind_injection_radius ?')
 if (iverbose>1) print *, 'Upper bound found for v0 :', v0max

! Find sonic point by dichotomy between v0min and v0max
 do
    v0last = v0
    v0 = (v0min+v0max)/2.
    if (iverbose>1) print *, 'v0 = ', v0
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
 !sonic point properties (location, time to reach, ...)
 sonic(1) = state%r
 sonic(2) = state%v
 sonic(3) = state%c
 sonic(4) = state%time
 sonic(5) = state%Tg
 sonic(6) = state%p
 sonic(7) = state%alpha
 !mdot = 4.*pi*rho*v0*ro*ro

 write (*,'("Sonic point properties  vs (km/s) =",f9.3,"  Rs/R* = ",f7.3," theoric = ",f7.3," Ts =",f7.1," alpha =",f5.3,/)')&
            sonic(2)/1e5,sonic(1)/Rstar_cgs,Rs/Rstar_cgs,sonic(5),sonic(7)
 !save 1D initial profile for comparison
 call save_windprofile(R0, T0, v0, sonic(4), 'gailstatwind1D.dat')

end subroutine get_initial_wind_speed


end module dust_free_wind
