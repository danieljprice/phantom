!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
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
! :Dependencies: cooling, dust_formation, eos, io, options, part, physcon,
!   ptmass_radiation, timestep, units, wind_equations
!

 implicit none
 public :: setup_wind
 public :: wind_state,wind_profile,save_windprofile

 private
 ! Shared variables
 real, parameter :: Tdust_stop = 1.d-2 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'wind'

 ! input parameters
 real :: Mstar_cgs, Lstar_cgs, wind_gamma, Mdot_cgs, Tstar
 real :: u_to_temperature_ratio
#ifdef NUCLEATION
 integer, parameter :: nwrite = 19
#else
 integer, parameter :: nwrite = 12
#endif

 ! wind properties
 type wind_state
    real :: dt, time, r, r0, Rstar, v, a, time_end, Tg, Teq, mu
    real :: gamma, alpha, rho, p, c, dalpha_dr, r_old, Q, dQ_dr
    real :: tau_lucy, kappa
#ifdef NUCLEATION
    real :: JKmuS(7)
#endif
    integer :: spcode, nsteps
    logical :: dt_force, error, find_sonic_solution
 end type wind_state
contains

subroutine setup_wind(Mstar_cg, Mdot_code, u_to_T, r0, T0, v0, rsonic, tsonic, stype)
 use ptmass_radiation, only:iget_Tdust
 use units,            only:umass,utime
 use physcon,          only:c,years
 use eos,              only:gamma
#ifdef NUCLEATION
 use dust_formation,   only:set_abundances
 use part,             only:n_nucleation
 use io,               only:fatal
 type(wind_state) :: state
#endif

 integer, intent(in) :: stype
 real, intent(in)    :: Mstar_cg, Mdot_code, u_to_T, T0
 real, intent(inout) :: r0,v0
 real, intent(out)   :: rsonic, tsonic

 Mstar_cgs = Mstar_cg
 wind_gamma = gamma
 Mdot_cgs = Mdot_code * umass/utime
 u_to_temperature_ratio = u_to_T

#ifdef NUCLEATION
 if (size(state%JKmuS) /= n_nucleation-1) call fatal(label,'wrong dimension for JKmuS')
 call set_abundances
#endif

 if (iget_Tdust == 2) then
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
 use dust_formation,   only:kappa_gas,evolve_chem
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
 state%time = 0.
 state%r_old = 0.
 state%r0 = r0
 state%r = r0
 state%Rstar = r0
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
 state%alpha = 0.
#endif
 state%tau_lucy = 2./3.
 state%mu = gmw
 state%kappa = kappa_gas
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
 state%mu = state%jKmuS(6)
#endif
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

 kappa_old = state%kappa
#ifdef NUCLEATION
 call evolve_chem(state%dt,state%Tg,state%rho,state%JKmuS)
 alpha_old = state%alpha
 state%mu  = state%JKmuS(6)
 call calc_kappa_dust(state%JKmuS(5), state%Teq, state%rho, state%kappa)
 call calc_alpha_dust(Mstar_cgs, Lstar_cgs, state%kappa, state%alpha)
 if (state%time > 0.) state%dalpha_dr = (state%alpha-alpha_old)/(1.+state%r-state%r_old)
#endif
 rvT(1) = state%r
 rvT(2) = state%v
 rvT(3) = state%Tg
 v_old = state%v
 state%r_old = state%r
 call evolve_hydro(state%dt, rvT, state%Rstar, state%mu, state%gamma, state%alpha, state%dalpha_dr, &
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

#ifdef NUCLEATION
 call calc_kappa_dust(state%JKmuS(5), state%Teq, state%rho, state%kappa)
#else
 if (idust_opacity > 0) state%kappa = kappa_dust_bowen(state%Teq)
#endif
 state%tau_lucy = state%tau_lucy &
      - (state%r-state%r_old) * state%r0**2 &
      * (state%kappa*state%rho/state%r**2 + kappa_old*rho_old/state%r_old**2)/2.
 if (icooling > 0) then
    Q_old = state%Q
    if (calc_Teq) then
       tau_lucy_bounded = max(0., state%tau_lucy)
       state%Teq = Tstar * (.5*(1.-sqrt(1.-(state%r0/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
    endif
#ifdef NUCLEATION
    call calc_cooling_rate(state%Q,dlnQ_dlnT,state%rho,state%Tg,state%Teq,state%JKmuS(6),state%JKmuS(4),state%kappa)
#else
    call calc_cooling_rate(state%Q,dlnQ_dlnT,state%rho,state%Tg,state%Teq,state%mu)
#endif
    if (state%time > 0. .and. state%r /= state%r_old) state%dQ_dr = (state%Q-Q_old)/(1.d-10+state%r-state%r_old)
 else
    !if cooling disabled or no imposed temperature profile, set Tdust = Tgas
    state%Teq = state%Tg
 endif
 if (state%time_end > 0. .and. state%time + state%dt > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 endif
 state%nsteps = mod(state%nsteps, 65536)+1
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
 use ptmass_radiation, only:iget_Tdust
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: tau_lucy_last

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
    if (iget_tdust == 2 .and. (tau_lucy_last-state%tau_lucy)/tau_lucy_last < 1.e-6 .and. state%tau_lucy < .6) exit
    if (state%r == state%r_old .or. state%tau_lucy < -1.) state%error = .true.
 enddo

end subroutine calc_wind_profile

!-----------------------------------------------------------------------
!+
!  integrate wind equation up to time=local_time
!+
!-----------------------------------------------------------------------
subroutine wind_profile(local_time,r,v,u,rho,e,GM,T0,fdone,JKmuS)
 !in/out variables in code units (except Jstar,K)
 use units,        only:udist, utime, unit_velocity, unit_density!, unit_pressure
 use eos,          only:gamma
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
 r = state%r/udist
 v = state%v/unit_velocity
 rho = state%rho/unit_density
 u = state%Tg * u_to_temperature_ratio
 !u = state%p/((gamma-1.)*rho)/unit_pressure
 e = .5*v**2 - GM/r + gamma*u
#ifdef NUCLEATION
 JKmuS(1:7) = state%JKmuS(1:7)
 JKmuS(8) = state%kappa
#endif

end subroutine wind_profile

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

 real :: v0min, v0max, v0last, vesc, cs, Rs, alpha_max,vin
 integer, parameter :: ncount_max = 10
 integer :: icount
 character(len=*), parameter :: label = 'get_initial_wind_speed'

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-alpha_rad)/r0)
 cs   = sqrt(gamma*Rg*T0/gmw)
 vin  = cs*(vesc/2./cs)**2*exp(-(vesc/cs)**2/2.+1.5)
 Rs   = Gg*Mstar_cgs*(1.-alpha_rad)/(2.*cs*cs)
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
       print *, ' * v0  (km/s) = ',vin/1e5
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
    v0 = cs*0.991
    v0max = v0
    v0min = 0.
    icount = 0
    do while (icount < ncount_max)
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *,' v0/cs = ',v0/cs
       if (state%spcode == -1) then
          v0min = v0
          exit
       else
          v0max = v0
          v0 = v0 / 2.
       endif
       icount = icount+1
    enddo
    if (iverbose>1) print *, 'Lower bound found for v0/cs :',v0min/cs
    if (icount == ncount_max) call fatal(label,'cannot find v0min, change wind_temperature or wind_injection_radius ?')
    if (v0min/cs > 0.99) call fatal(label,'supersonic wind, set sonic_type = 0 and provide wind_velocity or change alpha_rad')

    ! Find upper bound for initial velocity
    v0 = v0max
    icount = 0
    do while (icount < ncount_max)
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *,' v0/cs = ',v0/cs
       if (state%spcode == 1) then
          v0max = v0
          exit
       else
          v0min = max(v0min, v0)
          v0 = v0 * 1.1
       endif
       icount = icount+1
    enddo
    if (icount == ncount_max) call fatal(label,'cannot find v0max, change wind_temperature or wind_injection_radius ?')
    if (iverbose>1) print *, 'Upper bound found for v0/cs :', v0max/cs

    ! Find sonic point by dichotomy between v0min and v0max
    do
       v0last = v0
       v0 = (v0min+v0max)/2.
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *, 'v0/cs = ',v0/cs
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

    write (*,'("Sonic point properties  cs (km/s) =",f9.3,", Rs/r0 = ",f7.3,&
    &", v0/cs = ",f9.6,", ts =",f8.1)') &
         state%v/1e5,state%r/r0,v0/state%v,state%time/utime

 else
    if (v0 >= cs) then
       print *,' supersonic wind : v0/cs = ',v0/cs
    else
       print *,' sub-sonic wind : v0/cs = ',v0/cs
    endif
    call calc_wind_profile(r0, v0, T0, 0., state)
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
 !mdot = 4.*pi*rho*v0*ro*ro
 rsonic = state%r
 tsonic = state%time
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
!  Integrate the steady wind equation and save wind profile to a file
!
!-----------------------------------------------------------------------
subroutine save_windprofile(r0, v0, T0, rout, tsonic, tend, tcross, filename)
 use physcon,  only:au
 use units,    only:utime
 use timestep, only:tmax
 real, intent(in) :: r0, v0, T0, tsonic, tend, rout
 real, intent(out) :: tcross
 character(*), intent(in) :: filename
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 real :: dt_print,time_end
 type(wind_state) :: state
 integer :: n,iter

 write (*,'("Saving 1D model : ")')
 time_end = tmax*utime
 call init_wind(r0, v0, T0, tend, state)

 open(unit=1337,file=filename)
 call filewrite_header(1337)
 call filewrite_state(1337, state)

 n = 1
 iter = 0
 tcross = 1.d99
 dt_print = max(min(tsonic/10.,time_end/256.),time_end/5000.)
 do while(state%time < time_end .and. iter < 10000000 .and. state%Tg > Tdust_stop)
    iter = iter+1
    call wind_step(state)
    if (state%time > n*dt_print) then
       n = floor(state%time/dt_print)+1
       call filewrite_state(1337, state)
    endif
    if (state%r > rout) tcross = min(state%time,tcross)
 enddo
 if (state%time/time_end < .3) then
    write(*,'("[WARNING] wind integration failed : t/tend = ",f7.5,", dt/tend = ",f7.5,&
    &" Tgas = ",f6.0,", r/rout = ",f7.5," iter = ",i7,/)') &
         state%time/time_end,state%dt/time_end,state%Tg,state%r/rout,iter
 endif
 print *,'saveprofile : vinf = ',state%v,state%r/au
 close(1337)

end subroutine save_windprofile

subroutine filewrite_header(iunit)
 use options, only : icooling
 integer, intent(in) :: iunit

#ifdef NUCLEATION
 if (icooling > 0) then
    write(iunit,'("#",11x,a1,19(a20))') 't','r','v','T','c','p','rho','alpha','a',&
         'mu','S','Jstar','K0','K1','K2','K3','tau_lucy','kappa','Q'
 else
    write(iunit,'("#",11x,a1,18(a20))') 't','r','v','T','c','p','rho','alpha','a',&
         'mu','S','Jstar','K0','K1','K2','K3','tau_lucy','kappa'
 endif
#else
 if (icooling > 0) then
    write(iunit,'("#",11x,a1,10(a20))') 't','r','v','T','c','p','rho','alpha','a','mu','kappa','Q'
 else
    write(iunit,'("#",11x,a1,9(a20))') 't','r','v','T','c','p','rho','alpha','a','mu','kappa'
 endif
#endif
end subroutine filewrite_header

subroutine state_to_array(state, array)
 use options, only:icooling
 type(wind_state), intent(in) :: state
 real, intent(out) :: array(:)

 array(1) = state%time
 array(2) = state%r
 array(3) = state%v
 array(4) = state%Tg
 array(5) = state%c
 array(6) = state%p
 array(7) = state%rho
 array(8)  = state%alpha
 array(9) = state%a
#ifdef NUCLEATION
 array(10)  = state%JKmuS(6)
 array(11)  = state%JKmuS(7)
 array(12) = state%JKmuS(1)
 array(13:16) = state%JKmuS(2:5)
 array(17) = state%tau_lucy
 array(18) = state%kappa
#else
 array(10) = state%mu
 array(11) = state%kappa
#endif
 if (icooling > 0) array(nwrite) = state%Q
end subroutine state_to_array

subroutine filewrite_state(iunit, state)
 integer, intent(in) :: iunit
 type(wind_state), intent(in) :: state

 real :: array(nwrite)

 call state_to_array(state, array)
 write(iunit, '(20E20.10E3)') array(1:nwrite)
end subroutine filewrite_state

end module wind
