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
!  DEPENDENCIES: dust_physics, eos, io, part, physcon, timestep, units,
!    wind_profile
!+
!--------------------------------------------------------------------------

module dusty_wind
 implicit none
 public :: setup_dustywind
 public :: get_initial_wind_speed!, profile_findr0
 public :: wind_state, dusty_wind_profile,radiative_acceleration
 public :: evolve_gail

 private

! Shared variables
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 real, parameter :: dtmin = 1.d-3 ! Minimum allowed timsestep (for 1D integration)
 integer, parameter :: wind_emitting_sink = 1
 character(len=*), parameter :: label = 'dusty_wind'

! Wind properties
 real :: Mstar_cgs, Lstar_cgs, Tstar, Rstar_cgs, wind_gamma, expT, Mdot_cgs
 real :: Cprime, u_to_temperature_ratio, wind_alpha
 integer :: wind_type

! State of the wind
 type wind_state
    real :: dt, time, r, v, a, time_end, Tg, S, Jstar, K(0:3)
    real :: mu, gamma, alpha, rho, p, c, dalpha_dr, r_old, Q, dQ_dr
    real :: tau_lucy, kappa_ross, kappa_planck
    integer :: spcode, nsteps
    logical :: dt_force, error, stop_at_tmax
 end type wind_state
contains

subroutine setup_dustywind(Mstar_in, Lstar_in, Tstar_in, Rstar_cg, CO_ratio_in, Cprime_cgs, &
       expT_in, Mdot_in, u_to_T, alpha_in, wind_type_in,wind_cooling)
 use physcon,      only:c, solarm, years
 use eos,          only:gamma
 use units,        only:umass,utime
 use dust_physics, only:set_abundances,init_cooling

 integer, intent(in) :: wind_type_in,wind_cooling
 real, intent(in) :: Mstar_in, Lstar_in, Tstar_in, Rstar_cg, Cprime_cgs, expT_in,&
      Mdot_in, u_to_T, alpha_in, CO_ratio_in

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
 Mdot_cgs = Mdot_in * umass/utime
 Rstar_cgs = Rstar_cg
 u_to_temperature_ratio = u_to_T

 call set_abundances(CO_ratio_in)

 call init_cooling(wind_cooling,Cprime_cgs)

end subroutine setup_dustywind

!-----------------------------------------------------------------------
!
!  set particle dust properties
!
!-----------------------------------------------------------------------
subroutine evolve_gail(dtlast, xyzh, vxyzu, xyzmh_ptmass, vxyz_ptmass, partJstarKmu, npart)
 use units,          only:udist,utime,unit_velocity
 use physcon,        only:pi
!use part,           only:massoftype,igas
 use eos,            only:gmw
 use dust_physics,   only:evolve_chem
 real,    intent(in) :: dtlast
 real,    intent(in) :: xyzh(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 real,    intent(inout) :: vxyzu(:,:), partJstarKmu(:,:)
 integer, intent(in) :: npart

 integer, parameter :: wind_emitting_sink = 1
 integer :: i
 real :: x, y, z, h, vx, vy, vz, rho, dt
 real :: r, v, T, Jstar, K(0:3), mu, S

 if (wind_type <= 1) stop 'problem with wind_type in injject_dust'
 dt = dtlast* utime
 do i=1,npart
    x = xyzh(1,i) - xyzmh_ptmass(1,wind_emitting_sink)
    y = xyzh(2,i) - xyzmh_ptmass(2,wind_emitting_sink)
    z = xyzh(3,i) - xyzmh_ptmass(3,wind_emitting_sink)
    h = xyzh(4,i)
    vx = vxyzu(1,i) - vxyz_ptmass(1,wind_emitting_sink)
    vy = vxyzu(2,i) - vxyz_ptmass(2,wind_emitting_sink)
    vz = vxyzu(3,i) - vxyz_ptmass(3,wind_emitting_sink)
    r = sqrt(x**2+y**2+z**2)*udist
    v = sqrt(vx**2+vy**2+vz**2)*unit_velocity
    !rho = rhoh(xyzh(4,i), part_mass)*(unit_density)
    rho    = Mdot_cgs/(4.*pi*r**2*v)
    Jstar  = partJstarKmu(1,i)
    K(0:3) = partJstarKmu(2:5,i)
    mu = partJstarKmu(6,i)
    if (wind_type == 2) T = Tstar*(Rstar_cgs/r)**expT
    if (wind_type == 3 .or. wind_type == 4) T = mu*vxyzu(4,i)/(u_to_temperature_ratio*gmw)
    call evolve_chem(dt, r, v, T, rho, Jstar, K, mu, S)
    partJstarKmu(1,i)   = Jstar
    partJstarKmu(2:5,i) = K(0:3)
    partJstarKmu(6,i)   = mu
    !even though we consider an isothermal gas, internal energy can change because of mu changes
    !if (wind_type == 1) then
    !   vxyzu(4,i) =  u_to_temperature_ratio*T*gmw/mu
    !endif
 enddo
end subroutine evolve_gail

subroutine radiative_acceleration(npart, xyzh, vxyzu, dt, fext, fxyzu, time)
 use part,         only:rhoh,xyzmh_ptmass,massoftype,igas,partJstarKmu
 !use eos,          only:gmw!,gamma
 use physcon,      only:Rg
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:)
 real, intent(in)    :: dt,time
 real, intent(inout) :: vxyzu(:,:),fxyzu(:,:), fext(:,:)

 integer :: i
 real :: Mstar, Rstar
 real :: part_mass, r, xr(3), x_star(3), K3
 real :: rho, T, mu, Teq, ueq, tau_lucy, Q, alpha

 part_mass = massoftype(igas)
 x_star(1:3) = xyzmh_ptmass(1:3,wind_emitting_sink)
 Mstar = xyzmh_ptmass(4,wind_emitting_sink)
 Rstar = xyzmh_ptmass(5,wind_emitting_sink)
 Q = 0.d0
 !Tstar = stars(attached_to_star)%temperature
 do i=1,npart
    xr(1:3) = xyzh(1:3,i)-x_star(1:3)
    r = sqrt(xr(1)**2 + xr(2)**2 + xr(3)**2)
    K3 = partJstarKmu(5,i)
    call calc_alpha(K3, alpha)
    fext(1:3,i) = fext(1:3,i) + alpha*Mstar/r**3*xr(1:3)
    ! if (wind_type == 2 .or. r < Rstar) then
    !    fxyzu(4,i) = 0.
    ! elseif (wind_type == 4) then
    !    rho = rhoh(xyzh(4,i), part_mass)
    !    mu = partJstarKmu(6,i)
    !    T = mu*vxyzu(4,i)/(u_to_temperature_ratio*gmw)
    !    !T = (gamma-1.)*mu*vxyzu(4,i)/Rg
    !    tau_lucy = 0.
    !    Teq = star_Teff * (.5*(1.-sqrt(1.-(Rstar/r)**2)+3./2.*tau_lucy))**(1./4.)
    !    ueq = u_to_temperature_ratio*gmw*Teq/mu
    !    if (dt >= Cprime/rho) then
    !       vxyzu(4,i) = ueq
    !       if (vxyzu(4,i) < 0.) then
    !          print *, 'aargh! ueq <0'
    !       endif
    !       fxyzu(4,i) = 0.
    !    else
    !       Q = -(vxyzu(4,i) -ueq)*rho/Cprime
    !       fxyzu(4,i) = fxyzu(4,i) + Q
    !       if (vxyzu(4,i)+fxyzu(4,i)*dt < 0.) then
    !          print *, 'oups! u^n+1 < 0 - timestep likely too large',dt,Cprime/rho,rho
    !       endif
    !    endif
    ! endif
 enddo
end subroutine radiative_acceleration

!-----------------------------------------------------------------------
!
!  calculate alpha, reduced gravity factor
!
!-----------------------------------------------------------------------
subroutine calc_alpha(K3, alpha)
!all quantities in cgs
 use physcon,      only:pi,c,Gg
 use dust_physics, only:calc_kappa_dust
 real, intent(in) :: K3
 real, intent(out) :: alpha

 real :: kappa_planck, dummy

 call calc_kappa_dust(K3, Mdot_cgs, kappa_planck, dummy)
 alpha = Lstar_cgs/(4.*pi*c*Gg*Mstar_cgs) * kappa_planck
end subroutine calc_alpha

!-----------------------------------------------------------------------
!
!  Initialize wind variables
!
!-----------------------------------------------------------------------
subroutine init_dustywind(r0, v0, T0, time_end, state)
! all quantities in cgs
 use physcon,      only:pi,Rg
 use io,           only:fatal
 use eos,          only:gmw
 use dust_physics, only:evolve_chem,calc_kappa_dust,calc_cooling_rate,calc_Teq
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: tau_lucy_bounded,Teq, dlnq_dlnT

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
 state%mu = gmw
 state%gamma = wind_gamma
 state%S = 0.
 state%Jstar = 0.
 state%K(:) = 0.
 state%tau_lucy = 2./3.
 state%Q = 0.
 state%dQ_dr = 0.
 state%rho = Mdot_cgs/(4.*pi * state%r**2 * state%v)

 call evolve_chem(0., r0, v0, T0, state%rho, state%Jstar, state%K, state%mu, state%S)
 call calc_alpha(state%K(3), state%alpha)
 call calc_kappa_dust(state%K(3), Mdot_cgs, state%kappa_planck, state%kappa_ross)

 if (wind_type == 4) then
    if (r0 < Rstar_cgs) then
       call fatal(label,'cannot determine equilibrium temperature because injection_radius < Rstar')
    else
       if (calc_Teq) then
          tau_lucy_bounded = max(0., state%tau_lucy)
          Teq = Tstar * (.5*(1.-sqrt(1.-(Rstar_cgs/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
       endif
       call calc_cooling_rate(state%Q, dlnq_dlnT, state%rho, state%Tg, Teq, wind_gamma, state%mu, &
            state%K(2)/(state%r**2*state%v), state%kappa_ross)
    endif
 endif
 state%p = state%rho*Rg*state%Tg/state%mu
 state%c = sqrt(wind_gamma*Rg*state%Tg/state%mu)
 state%dt_force = .false.
 state%spcode = 0
 state%nsteps = 1
 state%error = .false.

end subroutine init_dustywind


!-----------------------------------------------------------------------
!
!  Integrate chemistry, cooling and hydro over one time step
!
!-----------------------------------------------------------------------
subroutine dustywind_step(state)
! all quantities in cgs

 use wind_equations, only:evolve_hydro
 use physcon,        only:pi,Rg
 use dust_physics,   only:evolve_chem,calc_kappa_dust,calc_cooling_rate,calc_Teq
 type(wind_state), intent(inout) :: state

 real :: rvT(3), alpha_old, kappa_ross_old, rho_old, dt_next, Q_old, v_old, tau_lucy_bounded, Teq, dlnQ_dlnT

 call evolve_chem(state%dt,state%r,state%v,state%Tg,state%rho,state%Jstar,state%K,state%mu,state%S)
 alpha_old = state%alpha
 call calc_alpha(state%K(3), state%alpha)
 if (state%time > 0.) state%dalpha_dr = (state%alpha-alpha_old)/(1.+state%r-state%r_old)

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
 rho_old = state%rho
 state%rho = Mdot_cgs/(4.*pi*state%r**2*state%v)
 state%p = state%rho*Rg*state%Tg/state%mu
 state%c = sqrt(wind_gamma*Rg*state%Tg/state%mu)
 kappa_ross_old = state%kappa_ross

 call calc_kappa_dust(state%K(3), Mdot_cgs, state%kappa_planck, state%kappa_ross)
 state%tau_lucy = state%tau_lucy &
      - (state%r-state%r_old) * Rstar_cgs**2 &
      * (state%kappa_ross*state%rho/state%r**2 + kappa_ross_old*rho_old/state%r_old**2)/2.

 if (wind_type == 4) then
    Q_old = state%Q
    if (calc_Teq) then
       tau_lucy_bounded = max(0., state%tau_lucy)
       Teq = Tstar * (.5*(1.-sqrt(1.-(Rstar_cgs/state%r)**2)+3./2.*tau_lucy_bounded))**(1./4.)
    endif
    call calc_cooling_rate(state%Q, dlnQ_dlnT, state%rho, state%Tg, Teq, wind_gamma, state%mu, &
         state%K(2)/(state%r**2*state%v),state%kappa_ross)
    if (state%time > 0. .and. state%r /= state%r_old) state%dQ_dr = (state%Q-Q_old)/(1.d-10+state%r-state%r_old)
 endif
 if (state%time_end > 0. .and. state%time + state%dt > state%time_end) then
    state%dt = state%time_end-state%time
    state%dt_force = .true.
 endif
 state%nsteps = state%nsteps + 1
 if (state%stop_at_tmax .and. state%time < state%time_end) state%spcode = 0

end subroutine dustywind_step

!-----------------------------------------------------------------------
!
!  Integrate the dusty wind equation up to t = t_end
!
!-----------------------------------------------------------------------
subroutine calc_dustywind_profile(r0, v0, T0, time_end, state)
! all quantities in cgs
 real, intent(in) :: r0, v0, T0, time_end
 type(wind_state), intent(out) :: state
 real :: tau_lucy_last
 logical :: tau_test

!compute chemistry and initialize variables
 call init_dustywind(r0, v0, T0, time_end, state)
 !print *,'calc_wind_profile',r0, v0, T0, time_end,state%v > state%c,state%time

 if (state%v > state%c) then
    state%spcode = 1
    !print *, 'Initial velocity cannot be greater than sound speed'
 endif

!integrate 1D wind solution with dust
 do while(state%dt> dtmin .and. state%Tg > Tdust_stop .and. .not.state%error .and. state%spcode == 0)
    tau_lucy_last = state%tau_lucy

    call dustywind_step(state)
    !print *,state%time,state%dt,state%Tg,state%spcode,state%r,state%v,state%tau_lucy
    tau_test = (tau_lucy_last-state%tau_lucy)/tau_lucy_last < 1.e-6 .and. state%tau_lucy < .6
    if (state%r == state%r_old .or. state%tau_lucy < -1. .or. tau_test) state%error = .true.
 enddo

end subroutine calc_dustywind_profile

!-----------------------------------------------------------------------
!+
!  dusty wind model
!+
!-----------------------------------------------------------------------
subroutine dusty_wind_profile(time,local_time,r,v,u,rho,e,GM,gamma,T,Jstar,K,mu,cs)
 !in/out variables in code units (except Jstar,K,mu)
 use units,        only:udist, utime, unit_velocity, unit_density, unit_pressure
 real, intent(in)  :: time,local_time, GM, gamma, T
 real, intent(inout) :: r, v
 real, intent(out) :: u, rho, e, Jstar, K(:), mu, cs

 type(wind_state) :: state

 r = r*udist
 v = v*unit_velocity
 call calc_dustywind_profile(r, v, T, local_time*utime, state)
 r = state%r/udist
 v = state%v/unit_velocity
 rho = state%rho/(unit_density)
 u = state%p/((gamma-1.)*rho)/unit_pressure
 e = .5*v**2 - GM/r + gamma*u
 Jstar = state%Jstar
 K = state%K
 mu = state%mu
 cs = state%c/unit_velocity
end subroutine dusty_wind_profile

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
 call init_dustywind(r0, v0, T0, time_end, state)
 open(unit=1337,file=filename)
 call filewrite_header(1337)
 call filewrite_state(1337, state)

 n = 1
 dt_print = min(tsonic/10.,state%time_end/256.)
 do while(state%time < state%time_end .and. state%dt > dtmin .and. state%Tg > Tdust_stop)
    call dustywind_step(state)
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
    write(iunit,'("#",11x,a1,19(a20))') 't','r','v','T','c','p','rho','mu','S','Jstar',&
         'K0','K1','K2','K3','tau_lucy','kappa_planck','kappa_ross','alpha','a','Q'
 else
    write(iunit,'("#",11x,a1,18(a20))') 't','r','v','T','c','p','rho','mu','S','Jstar',&
         'K0','K1','K2','K3','tau_lucy','kappa_planck','kappa_ross','alpha','a'
 endif
end subroutine filewrite_header

subroutine state_to_array(state, array)
 type(wind_state), intent(in) :: state
 real, intent(out) :: array(20)
 real :: f

 f = state%r**2 * state%v
 array(1) = state%time
 array(2) = state%r
 array(3) = state%v
 array(4) = state%Tg
 array(5) = state%c
 array(6) = state%p
 array(7) = state%rho
 array(8) = state%mu
 array(9) = state%S
 array(10) = state%Jstar/f
 array(11:14) = state%K/f
 array(15) = state%tau_lucy
 array(16) = state%kappa_planck
 array(17) = state%kappa_ross
 array(18) = state%alpha
 array(19) = state%a
 if (wind_type == 4) array(20) = state%Q
end subroutine state_to_array

subroutine filewrite_state(iunit, state)
 integer, intent(in) :: iunit
 type(wind_state), intent(in) :: state

 real :: array(20)
 call state_to_array(state, array)
 if (wind_type == 4) then
    write(iunit, '(20E20.10E3)') array(1:20)
 else
    write(iunit, '(19E20.10E3)') array(1:19)
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
 real, intent(out) :: v0, sonic(8)

 type(wind_state) :: state

 real :: v0min, v0max, v0last, vesc, cs, Rs, alpha_max
 integer, parameter :: ncount_max = 20
 integer :: icount
 character(len=*), parameter :: label = 'get_initial_wind_speed'

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-wind_alpha)/r0)
 cs = sqrt(wind_gamma*Rg*T0/gmw)
 v0 = cs*(vesc/2./cs)**2*exp(-(vesc/cs)**2/2.+1.5)
 Rs = Gg*Mstar_cgs*(1.-wind_alpha)/(2.*cs*cs)
 alpha_max = 1.-(2.*cs/vesc)**2
 if (iverbose>0) then
    print *, "[get_initial_wind_speed] Looking for initial velocity."
    print *, ' * unit(au)   = ',udist/au
    print *, ' * Mstar      = ',Mstar_cgs/1.9891d33
    print *, ' * Twind      = ',T0
    print *, ' * Tstar      = ',Tstar
    print *, ' * Lstar      = ',Lstar_cgs/3.9d33
    print *, ' * Rstar(au)  = ',Rstar_cgs/1.496d13
    print *, ' * Mdot       = ',Mdot_cgs/6.30303620274d25
    print *, ' * r0(au)     = ',r0/1.496d13,r0/69600000000.
    print *, ' * gamma      = ',wind_gamma
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
 do while (icount < ncount_max)
    call calc_dustywind_profile(r0, v0, T0, 0., state)
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
    call calc_dustywind_profile(r0, v0, T0, 0., state)
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
    call calc_dustywind_profile(r0, v0, T0, 0., state)
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
 sonic(8) = state%S
 !mdot = 4.*pi*rho*v0*ro*ro

 write (*,'("Sonic point properties  vs (km/s) =",f9.3,"  Rs/R* = ",f7.3," Ts =",f7.1,"K, alpha =",f6.3," S = ",es9.2)')&
            sonic(2)/1e5,sonic(1)/Rstar_cgs,sonic(5),sonic(7),sonic(8)
 !save 1D initial profile for comparison
 write (*,'("Saving 1D model : ")',advance='no')
 call save_windprofile(R0, T0, v0, sonic(4), 'gailstatwind1D.dat')

end subroutine get_initial_wind_speed

!-----------------------------------------------------------------------
!
!  Determine the initial stellar radius for a trans-sonic solution
!
!-----------------------------------------------------------------------
subroutine profile_findr0(time_end, T0, r0, v0, step_factor, nsteps, sonic)
!all quantities in cgs
 use physcon, only:steboltz,pi
 real, intent(in) :: time_end, T0, step_factor
 real, intent(inout) :: r0
 real, intent(out) :: v0, sonic(8)
 integer, intent(out) :: nsteps
 logical :: verbose = .false.

 real :: r0min, r0max, r0best, v0best, tau_lucy_best, initial_guess
 integer :: i, nstepsbest
 type(wind_state) :: state

! Find lower bound for initial radius
 !r0 = sqrt(Lstar_cgs/(4.*pi*steboltz*Tstar**4))
 initial_guess = r0
 r0max = r0
 if (verbose) print *, '[profile_findr0] Searching lower bound for r0'
 do
    Rstar_cgs = r0
    call get_initial_wind_speed(r0, T0, v0, sonic)
    call calc_dustywind_profile(r0, v0, T0, time_end, state)
    if (verbose) print *, state%error,' r0 = ', r0, 'tau_lucy = ', state%tau_lucy
    if (state%error) then
! something wrong happened!
       r0 = r0 / step_factor
    elseif (r0 < initial_guess/10.) then
       r0min = r0
       print *,'radius getting too small!',r0,step_factor
       exit
    elseif (state%tau_lucy < 0.) then
       r0min = r0
       exit
    else
       r0max = r0
       r0 = r0 / step_factor
    endif
 enddo
 if (verbose) print *, 'Lower bound found for r0 :', r0min

! Find upper bound for initial radius
 r0 = r0max
 if (verbose) print *, 'Searching upper bound for r0'
 do
    Rstar_cgs = r0
    call get_initial_wind_speed(r0, T0, v0, sonic)
    call calc_dustywind_profile(r0, v0, T0, time_end, state)
    if (verbose) print *, 'r0 = ', r0, 'tau_lucy = ', state%tau_lucy
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
 if (verbose) print *, 'Upper bound found for r0 :', r0max

! Find the initial radius by dichotomy between r0min and r0max
 if (verbose) print *, 'Searching r0 by dichotomy'
 r0best = r0
 v0best = v0
 nstepsbest = state%nsteps
 tau_lucy_best = state%tau_lucy
 do i=1,50
    r0 = (r0min+r0max)/2.
    Rstar_cgs = r0
    call get_initial_wind_speed(r0, T0, v0, sonic)
    call calc_dustywind_profile(r0, v0, T0, time_end, state)
    if (verbose) print *, 'r0 = ', r0, 'tau_lucy = ', state%tau_lucy
    if (abs(state%tau_lucy) < abs(tau_lucy_best)) then
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
 nsteps = nstepsbest
 if (verbose) then
    print *, 'Best initial radius found: r0 = ', r0,' , initial_guess = ',initial_guess
    print *, '  with v0 = ', v0,' , T0 = ',T0
    print *, '  it gives tau_lucy = ', tau_lucy_best
    print *, '  at t = ', time_end
 endif
end subroutine profile_findr0

end module dusty_wind
