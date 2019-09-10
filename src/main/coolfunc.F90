!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: coolfunc
!
!  DESCRIPTION: Gas cooling
!
!  REFERENCES:
!   Townsend (2009), ApJS 181, 391-397
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
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

module coolfunc
 implicit none
 character(len=*), parameter :: label = 'cooling_rates'

 public :: init_coolfunc,calc_cooling_rate,energ_coolfunc,read_options_coolfunc,&
      write_options_coolfunc
 logical, public :: calc_Teq
 real, public:: bowen_Cprime = 3.000d-5

 private
 integer, parameter :: nTg = 64
 integer, parameter :: maxt = 1000
 real,    parameter :: Tref = 1.d5, T_floor = 10.
 real :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt)
 real :: habund     = 0.7
 real :: temp_floor = 1.e4
 real :: Tgrid(nTg)
 logical :: cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan, cool_table
 character(len=120) :: cooltable = 'cooltable.dat'


contains

subroutine init_coolfunc(icool)
  !use units,        only:umass, utime, udist

  integer, intent(in) :: icool
  integer :: iwind

  if (icool > 0) then
     iwind = icool
  else
     !dust free wind, only H0 cooling allowed
     if (abs(icool) > 0) iwind = 10
  endif
  !Cprime = bowen_Cprime / (umass*utime/udist**3)

  cool_table = .false.
  cool_radiation_H0 = .false.
  cool_relaxation_Bowen = .false.
  cool_relaxation_Stefan = .false.
  cool_collisions_dust = .false.
  !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
  if (iwind > 9) cool_radiation_H0 = .true.
  select case (iwind)
  case (1)
     cool_table = .true.
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
  if (cool_table) then
     call init_cooltable(ierr)
  endif
  call set_Tgrid

end subroutine init_coolfunc


!-----------------------------------------------------------------------
!+
!  read cooling table from file and initialise arrays
!+
!-----------------------------------------------------------------------
subroutine init_cooltable(ierr)
 use io,        only:fatal
 use datafiles, only:find_phantom_datafile
 integer, intent(out) :: ierr
 integer, parameter :: iu = 127
 integer :: i
 character(len=120) :: filepath

 !
 ! read the cooling table from file
 !
 filepath=find_phantom_datafile(cooltable,'cooling')
 open(unit=iu,file=filepath,status='old',iostat=ierr)
 if (ierr /= 0) call fatal('coolfunc','error opening cooling table')
 i = 0
 do while(ierr==0 .and. i < maxt)
    i = i + 1
    read(iu,*,iostat=ierr) temper(i),lambda(i)
 enddo
 nt = i-1
 if (nt==maxt) call fatal('coolfunc','size of cooling table exceeds array size')
 if (nt < 2) call fatal('coolfunc','size of cooling table is too small',ival=nt,var='nt')
 !
 ! calculate the slope of the cooling function
 !
 do i=1,nt-1
    slope(i) = log(lambda(i+1)/lambda(i))/log(temper(i+1)/temper(i))
 enddo
 slope(nt) = slope(nt-1)

 !
 ! initialise the functions required for Townsend method
 !
 yfunc(nt) = 0.
 do i=nt-1,1,-1
  !Lionel Siess : I think there is an error. in yfunc slope(nt) should be replaced by lambda(nt)
  ! Eq A6
   if (abs(slope(i)-1.) < tiny(0.)) then
       !ori yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
       yfunc(i) = yfunc(i+1) - lambda(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
    else
       !ori yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
       yfunc(i) = yfunc(i+1) - lambda(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
                 *(1.- (temper(i)/temper(i+1))**(slope(i) - 1.))
    endif
 enddo
end subroutine init_cooltable

!-----------------------------------------------------------------------
!
!  calculate cooling rates
!
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Teq, mu, K2, kappa)
 use units,   only:unit_ergg,unit_density
 real, intent(in) :: rho, T !rho in code units
 real, intent(in), optional :: Teq, mu, K2, kappa !cgs
 real, intent(out) :: Q, dlnQ_dlnT !code units
 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, Q_table, rho_cgs
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan

 rho_cgs = rho*unit_density
 Q_H0 = 0.
 Q_relax_Bowen = 0.
 Q_col_dust = 0.
 Q_relax_Stefan = 0.
 dlnQ_H0 = 0.
 dlnQ_relax_Bowen = 0.
 dlnQ_col_dust = 0.
 dlnQ_relax_Stefan = 0.
 if (cool_radiation_H0) call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (cool_relaxation_Bowen) call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, Q_relax_Bowen, dlnQ_relax_Bowen)
 if (cool_collisions_dust)  call cooling_collision_dust_grains(T, Teq, rho_cgs, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (cool_relaxation_Stefan) call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
 Q_cgs = Q_H0 + Q_relax_Bowen+ Q_col_dust+ Q_relax_Stefan
 dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen+ Q_col_dust*dlnQ_col_dust+ Q_relax_Stefan*dlnQ_relax_Stefan)/Q_cgs
 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q = Q_cgs/unit_ergg
end subroutine calc_cooling_rate


!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Teq, rho, mu, Q, dlnQ_dlnT)
! all quantities in cgs
 use eos,     only:gamma
 use physcon, only:Rg
 real, intent(in) :: T, Teq, rho, mu
 real, intent(out) :: Q,dlnQ_dlnT

 Q = Rg/((gamma-1.)*mu)*rho*(Teq-T)/bowen_Cprime
 dlnQ_dlnT = -T/(Teq-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_collision_dust_grains(T, Teq, rho, K2, mu, Q, dlnQ_dlnT)
! all quantities in cgs
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
 use physcon, only: mass_proton_cgs, pi
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

 k1 = 1.88d-10 / T**6.44e-1
 k2 = 1.83d-18 * T
 k3 = 1.35d-9
 k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
 k9 = 1.7d-4 * k8
 p = .5*k8/k9
 q = k1*(k2+k3)/(k3*k9)
 calc_eps_e = (p + sqrt(q+p**2))/q
end function calc_eps_e

subroutine set_Tgrid
  integer :: i
  real :: dlnT
  dlnT = log(Tref)/(nTg-1)

  do i = 1,nTg
     Tgrid(i) = exp((i-1)*dlnT)
  enddo
end subroutine set_Tgrid

!-----------------------------------------------------------------------
!
!   explicit cooling cooling
!
!-----------------------------------------------------------------------
subroutine explicit_cooling (u, dudt, rho, dt, Trad, mu_in, K2, kappa)
  use eos,     only:gamma,gmw
  use physcon, only:Rg
  use units,   only:unit_ergg
  real, intent(in) :: rho, dt, Trad !code units
  real, intent(in), optional :: mu_in, K2, kappa
  real, intent(inout) :: u,dudt !code units

  real :: r,Q,dlnQ_dlnT,T,mu,T_on_u

  if (.not.present(mu_in)) then
     mu = gmw
  else
     mu = mu_in
  endif
  T_on_u = (gamma-1.)*mu*unit_ergg/Rg
  T = T_on_u*u
  call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
  if (-Q*dt  > u) then   ! assume thermal equilibrium
    u = Trad/T_on_u
    dudt = 0.d0
  else
    u = u + Q*dt
    dudt = dudt + Q
  endif
  !print *,T,Teq,T_on_u*u,'dT=',T_on_u*Q*dt,u,Q*dt

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!
!   implicit cooling
!
!-----------------------------------------------------------------------
subroutine implicit_cooling (u, dudt, rho, dt, Trad, mu_in, K2, kappa)
  use eos,     only:gamma,gmw
  use physcon, only:Rg
  use units,   only:unit_ergg
  real, intent(in) :: rho, dt
  real, intent(in), optional :: Trad, mu_in, K2, kappa
  real, intent(inout) :: u,dudt

  real, parameter :: tol = 1.d-4 ! to be adjusted
  integer, parameter :: iter_max = 200
  real :: Q,dlnQ_dlnT,T,mu,T_on_u,delta_u,term1,term2,term3
  integer :: iter

  if (.not.present(mu_in)) then
     mu = gmw
  else
     mu = mu_in
  endif
  T_on_u = (gamma-1.)*mu*unit_ergg/Rg
  delta_u = 1.d-3
  iter = 0
  !The pdv_work also depends on the internal energy and could also be included
  !in this loop provided this contribution was not accounted for in Force.F90
  ! see PP flag : IMPLICIT COOLING - pb: we need div(v) and it is only real*4
  !term2 = 1.-(gamma-1.)*dt*divcurlv !pdv=(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
  term2 = 1.
  term1 = u !updated internal energy without cooling contributions
  do while (abs(delta_u) > tol .and. iter < iter_max)
     T = u*T_on_u
     call calc_cooling_rate(Q,dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
     term3 = u*term2-Q*dt
     delta_u = (term1-term3)/(term2-Q*dlnQ_dlnT*dt/u)
     u = u+delta_u
     iter = iter + 1
  enddo
  ! dudt = dudt+(u-term1)/dt
  if (u < 0. .or. isnan(u)) then
     print *,u
     stop ' u<0'
  endif

end subroutine implicit_cooling

!-----------------------------------------------------------------------
!
!   implement cooling contribution du/dt equation
!
!-----------------------------------------------------------------------
subroutine energ_coolfunc (u, dudt, rho, dt, Trad, mu_in, K2, kappa)
  real, intent(in) :: rho, dt !in code unit
  real, intent(in), optional :: Trad, mu_in, K2, kappa ! in cgs!!!
  real, intent(inout) :: u, dudt

  if (icooling == itable) then
     call exact_cooling_table(u, rho, dt, dudt)
  else
     call exact_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
     !call implicit_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
     !call explicit_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
  endif

end subroutine energ_coolfunc

!-----------------------------------------------------------------------
!
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!
!-----------------------------------------------------------------------
subroutine exact_cooling (u, dudt, rho, dt, Trad, mu_in, K2, kappa)
  use eos,     only:gamma,gmw
  use physcon, only:Rg
  use units,   only:unit_ergg,utime
  use inject,  only:star_Teff,Rstar
  real, intent(in) :: rho, dt, Trad
  real, intent(in), optional :: mu_in, K2, kappa
  real, intent(inout) :: u,dudt

  real, parameter :: tol = 1.d-12
  real :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,mu,du,T_on_u
  integer :: k

  if (.not.present(mu_in)) then
     mu = gmw
  else
     mu = mu_in
  endif
  T_on_u = (gamma-1.)*mu*unit_ergg/Rg
  T = T_on_u*u

  if (T < T_floor) then
     Temp = T_floor
  elseif (T > Tref) then
     call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
     Temp = T+T_on_u*Q*dt
  else
     call calc_cooling_rate(Qref,dlnQref_dlnT, rho, Tref, Trad, mu, K2, kappa)
     Y = 0.
     k = nT
     do while (Tgrid(k) > T)
        k = k-1
        call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Trad, mu, K2, kappa)
        ! eqs A6
        if (abs(dlnQ_dlnT-1.) < tol) then
           y = y - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
         else
           y = y - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
        endif
     enddo
     !eqs A5
     yk = y
     if (abs(dlnQ_dlnT-1.) < tol) then
        y = yk + Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
     else
        y = yk + Qref*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1))
     endif
     !eq 26
     dy = Qref*dt*T_on_u/Tref
     y = y + dy
     !compute Yinv (eqs A7)
     if (abs(dlnQ_dlnT-1.) < tol) then
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

  du = (Temp-T)/T_on_u
  u = u+du
  !note that u = Temp/T_on_u
  dudt = dudt + du/dt

end subroutine exact_cooling

!-----------------------------------------------------------------------
!+
!  cooling using Townsend (2009) method with tabulated rate
!
!   Implements cooling defined using a tabulated cooling table
!   produced e.g. by CLOUDY.
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_table(uu,rho,dt,dudt)
 use eos,     only:gamma,gmw
 use physcon, only:atomic_mass_unit,kboltz,Rg
 use units,   only:unit_density,unit_ergg,utime
 real, intent(in)    :: rho,dt
 real, intent(inout) :: uu,dudt
 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp,durad
 real    :: sloperef,slopek,temp,temp1,tref,yfunx,yinv0
 integer :: k

 gam1 = gamma - 1.
 temp = gam1*uu/Rg*gmw*unit_ergg

 tref     = temper(nt)
 sloperef = slope(nt)

 if (temp < temp_floor) then
    temp1 = temp_floor
 else
    amue = 2.*atomic_mass_unit/(1. + habund)
    amuh = atomic_mass_unit/habund
    density_cgs = rho*unit_density
    dt_cgs      = dt*utime

 !Lionel Siess : I think there is an error. in dtemp sloperef should be replaced by lambda(nt)
    !original dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))* &
    !     sloperef/tref*dt_cgs
    ! Eq 26
    dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))*lambda(nt)/tref*dt_cgs

    k = find_in_table(nt,temper,temp)

    slopek = slope(k)
    ! Eq A5
    if (abs(slopek - 1.) < tiny(0.)) then
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt))*log(temper(k)/temp)
    else
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt)*(1. - slopek)) &
                          *(1. - (temper(k)/temp)**(slopek-1.))
    endif
    yfunx = yfunx + dtemp

    ! Eq A7
    if (abs(slopek - 1.) < tiny(0.)) then
       temp1 = max(temper(k)*exp(-lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))),temp_floor)
    else
       yinv0 = 1. - (1. - slopek)*lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))
       if (yinv0 > 0.) then
          temp1 = max(temper(k)*yinv0**(1./(1. - slopek)),temp_floor)
       else
          temp1 = temp_floor
       endif
    endif
 endif

 durad = (temp1 - temp)*Rg/(gam1*gmw*unit_ergg)
 !dudt = dudt + durad/dt
 uu = uu + durad

end subroutine exact_cooling_table

!-----------------------------------------------------------------------
!+
!  utility to find the index of closest value in a table
!+
!-----------------------------------------------------------------------
pure integer function find_in_table(n,table,val) result(i)
 integer, intent(in) :: n
 real,    intent(in) :: table(n), val
 integer :: i0,i1

 i0 = 0
 i1 = n + 1
 do while (i1 - i0 > 1)
    i = (i0 + i1)/2
    if ((table(n) >= table(1)).eqv.(val >= table(i))) then
       i0 = i
    else
       i1 = i
    endif
 enddo
 if (abs(val-table(1)) < tiny(0.)) then
    i = 1
 elseif (abs(val-table(n)) < tiny(0.)) then
    i = n-1
 else
    i = i0
 endif

end function find_in_table

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_coolfunc(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 if (icooling = itable) then
    call write_inopt(cooltable,'cooltable','data file containing cooling function',iunit)
    call write_inopt(habund,'habund','Hydrogen abundance assumed in cooling function',iunit)
    call write_inopt(temp_floor,'temp_floor','Minimum allowed temperature in K',iunit)
 endif
 if (icooling == idust) then
    call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
 endif

end subroutine write_options_coolfunc

subroutine read_options_coolfunc(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 select case(trim(name))
 case('cooltable')
    read(valstring,*,iostat=ierr) cooltable
    ngot = ngot + 1
 case('habund')
    read(valstring,*,iostat=ierr) habund
    ngot = ngot + 1
 case('temp_floor')
    read(valstring,*,iostat=ierr) temp_floor
    ngot = ngot + 1
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (icooling = itable .and. ngot >= 3) igotall = .true.
 if (icooling = idust  .and. ngot >= 1) igotall = .true.

end subroutine read_options_coolfunc

end module coolfunc
