!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling
!
! Gas cooling
!
! :References:
!   Koyama & Inutsuka (2002), ApJL 564, 97-100
!   Vazquez-Semadeni, et.al (2007), ApJ 657, 870-883
!   Townsend (2009), ApJS 181, 391-397
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - C_cool               : *factor controlling cooling timestep*
!   - beta_cool            : *beta factor in Gammie (2001) cooling*
!   - bowen_Cprime         : *radiative cooling rate (g.s/cm³)*
!   - cooltable            : *data file containing cooling function*
!   - habund               : *Hydrogen abundance assumed in cooling function*
!   - icool_dust_collision : *dust collision on/off*
!   - icool_radiation_H0   : *H0 cooling on/off*
!   - icool_relax_bowen    : *Bowen (diffusive) relaxation on/off*
!   - icool_relax_stefan   : *radiative relaxation on/off*
!   - icooling             : *cooling function (0=off, 1=explicit, 2=Townsend table, 3=Gammie, 5=KI02)*
!   - temp_floor           : *Minimum allowed temperature in K*
!
! :Dependencies: datafiles, eos, h2cooling, infile_utils, io, options,
!   part, physcon, timestep, units
!

 use options,  only:icooling
 use timestep, only:C_cool

 implicit none
 character(len=*), parameter :: label = 'cooling'

 public :: init_cooling,calc_cooling_rate,energ_cooling
 public :: write_options_cooling, read_options_cooling
 public :: find_in_table
 logical, public :: calc_Teq
 logical, public :: cooling_implicit
 logical, public :: cooling_explicit
 real,    public :: bowen_Cprime = 3.000d-5

 private
 integer, parameter :: nTg = 64
 integer, parameter :: maxt = 1000
 real,    parameter :: Tref = 1.d5, T_floor = 10.
 integer :: nt
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt)
 real    :: beta_cool  = 3.
 real    :: habund     = 0.7
 real    :: temp_floor = 1.e4
 real    :: Tgrid(nTg)
 real    :: crate_coef
 integer :: icool_radiation_H0 = 0, icool_relax_Bowen = 0, icool_dust_collision = 0, icool_relax_Stefan = 0
 character(len=120) :: cooltable = 'cooltable.dat'
 !--Minimum temperature (failsafe to prevent u < 0)
 real,    public :: Tfloor = 0. ! [K]; set in .in file.  On if Tfloor > 0.
 real,    public :: ufloor = 0. ! [code units]; set in init_cooling

contains

!-----------------------------------------------------------------------
!+
!  Initialise cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling(id,master,iprint,ierr)
 use dim,       only:maxvxyzu
 use units,     only:utime,umass,udist,unit_ergg
 use physcon,   only:mass_proton_cgs,kboltz
 use io,        only:fatal
 use eos,       only:gamma,gmw
 use part,      only:h2chemistry
 use h2cooling, only:init_h2cooling
 use chem,      only:init_chem
 integer, intent(in)  :: id,master,iprint
 integer, intent(out) :: ierr

 if (h2chemistry) then
    if (id==master) write(iprint,*) 'initialising cooling function...'
    call init_chem()
    call init_h2cooling()
 else
    !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
    if (icool_relax_bowen == 1 .and. icool_relax_stefan == 1) then
       call fatal(label,'you can"t have bowen and stefan cooling at the same time')
    endif

#ifdef KROME
    !krome calculates its own cooling rate
    icool_radiation_H0 = 0
    icool_dust_collision = 0
#else
    !if no cooling flag activated, disable cooling
    if (icooling == 1 .and. (icool_radiation_H0+icool_relax_Bowen+icool_dust_collision+&
          icool_relax_Stefan == 0)) then
       icooling = 0
       calc_Teq = .false.
       return
    endif
#endif
    calc_Teq = (icool_relax_Bowen == 1) .or. (icool_relax_Stefan == 1) .or. (icool_dust_collision == 1)

    !--initialise remaining variables
    if (icooling == 2) then
       call init_cooltable(ierr)
    elseif (icooling == 5) then
       crate_coef = 2.0d-26*umass*utime**3/(mass_proton_cgs**2 * udist**5)
    elseif (icooling > 0) then
       call set_Tgrid
    endif
 endif

 !--Determine if this is implicit or explicit cooling
 cooling_implicit = .false.
 cooling_explicit = .false.
 if (h2chemistry) then
    if (icooling > 0) cooling_implicit = .true.    ! cooling is calculated implicitly in step
 elseif (icooling > 0) then
    if (icooling == 3 .or. icooling == 5) then
       cooling_explicit = .true.                   ! cooling is calculated explicitly in force
    else
       cooling_implicit = .true.                   ! cooling is calculated implicitly in step
    endif
 endif

 !--calculate the energy floor in code units
 if (Tfloor > 0.) then
    if (gamma > 1.) then
       ufloor = kboltz*Tfloor/((gamma-1.)*gmw*mass_proton_cgs)/unit_ergg
    else
       ufloor = 3.0*kboltz*Tfloor/(2.0*gmw*mass_proton_cgs)/unit_ergg
    endif
    if (maxvxyzu < 4) ierr = 1
 else
    ufloor = 0.
 endif

end subroutine init_cooling

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
 if (ierr /= 0) call fatal('cooling','error opening cooling table')
 i = 0
 do while(ierr==0 .and. i < maxt)
    i = i + 1
    read(iu,*,iostat=ierr) temper(i),lambda(i)
 enddo
 nt = i-1
 if (nt==maxt) call fatal('cooling','size of cooling table exceeds array size')
 if (nt < 2) call fatal('cooling','size of cooling table is too small',ival=nt,var='nt')
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
 real, intent(in) :: rho, T, Teq !rho in code units
 real, intent(in), optional :: mu, K2, kappa !cgs
 real, intent(out) :: Q, dlnQ_dlnT !code units
 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, rho_cgs
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
 if (icool_radiation_H0 == 1)   call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (icool_relax_Bowen == 1)    call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, Q_relax_Bowen, dlnQ_relax_Bowen)
 if (icool_dust_collision == 1) call cooling_dust_collision(T, Teq, rho_cgs, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (icool_relax_Stefan == 1)   call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
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
subroutine cooling_dust_collision(T, Teq, rho, K2, mu, Q, dlnQ_dlnT)
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
end subroutine cooling_dust_collision

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
!+
!   Gammie (2001) cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_Gammie(xi,yi,zi,ui,dudti)
 real, intent(in)    :: ui,xi,yi,zi
 real, intent(inout) :: dudti

 real :: omegai,r2,tcool1

 r2     = xi*xi + yi*yi + zi*zi
 Omegai = r2**(-0.75)
 tcool1 = Omegai/beta_cool
 dudti  = dudti - ui*tcool1

end subroutine cooling_Gammie

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska(rhoi,Tgas,dudti)
 real, intent(in)    :: rhoi,Tgas
 real, intent(inout) :: dudti
 real                :: crate

 crate = crate_coef*(1.d7*exp(-118400./(Tgas+1000.))+0.014*sqrt(Tgas)*exp(-92./Tgas))
 dudti = dudti - crate*rhoi

end subroutine cooling_KoyamaInutuska

!-----------------------------------------------------------------------
!
!   explicit cooling
!
!-----------------------------------------------------------------------
subroutine explicit_cooling (ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: ui, rho, dt, Trad !code units
 real, intent(in), optional :: mu_in, K2, kappa
 real, intent(out) :: dudt !code units

 real :: u,Q,dlnQ_dlnT,T,mu,T_on_u

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T = T_on_u*ui
 call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
 if (-Q*dt  > ui) then   ! assume thermal equilibrium
    u = Trad/T_on_u
    dudt = (u-ui)/dt
 else
    dudt = Q
 endif
 !print *,T,Teq,T_on_u*u,'dT=',T_on_u*Q*dt,u,Q*dt

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!
!   implicit cooling
!
!-----------------------------------------------------------------------
subroutine implicit_cooling (ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: ui, rho, dt
 real, intent(in), optional :: Trad, mu_in, K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-4 ! to be adjusted
 integer, parameter :: iter_max = 200
 real :: u,Q,dlnQ_dlnT,T,mu,T_on_u,delta_u,term1,term2,term3
 integer :: iter

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 u = ui
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 delta_u = 1.d-3
 iter = 0
 !The pdv_work also depends on the internal energy and could also be included
 !in this loop provided this contribution was not accounted for in Force.F90
 ! see PP flag : IMPLICIT COOLING - pb: we need div(v) and it is only real*4
 !term2 = 1.-(gamma-1.)*dt*divcurlv !pdv=(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
 term2 = 1.
 term1 = u !initial internal energy without cooling contributions
 do while (abs(delta_u) > tol .and. iter < iter_max)
    T = u*T_on_u
    call calc_cooling_rate(Q,dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
    term3 = u*term2-Q*dt
    delta_u = (term1-term3)/(term2-Q*dlnQ_dlnT*dt/u)
    u = u+delta_u
    iter = iter + 1
 enddo
 dudt =(u-term1)/dt
 if (u < 0. .or. isnan(u)) then
    print *,u
    stop ' u<0'
 endif

end subroutine implicit_cooling

!-----------------------------------------------------------------------
!
!   this routine returns the effective cooling rate du/dt
!
!-----------------------------------------------------------------------
subroutine energ_cooling(xi,yi,zi,ui,dudt,rho,dt,Trad,mu_in,K2,kappa,Tgas)
 use io, only: fatal
 real, intent(in)           :: xi,yi,zi,ui,rho,dt         ! in code units
 real, intent(in), optional :: Tgas,Trad,mu_in,K2,kappa   ! in cgs units
 real, intent(inout)        :: dudt                       ! in code units

 select case (icooling)
 case (3)
    call cooling_Gammie(xi,yi,zi,ui,dudt)
 case (2)
    call exact_cooling_table(ui,rho,dt,dudt)
 case (5)
    if (present(Tgas)) then
       call cooling_KoyamaInutuska(rho,Tgas,dudt)
    else
       call fatal('energ_cooling','Koyama & Inutuska cooling requires gas temperature')
    endif
 case default
    !call exact_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
    !call implicit_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
    if (present(Trad) .and. present(mu_in) .and. present(K2) .and. present(kappa)) then
       call explicit_cooling(ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
    else
       call fatal('energ_cooling','default requires optional arguments; change icooling or ask D Price or L Siess to patch')
    endif
 end select

end subroutine energ_cooling

!-----------------------------------------------------------------------
!
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!
!-----------------------------------------------------------------------
subroutine exact_cooling (u, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: u, rho, dt, Trad
 real, intent(in), optional :: mu_in, K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,mu,T_on_u
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
    k = nTg
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

 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

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
 real, intent(in)  :: uu, rho,dt
 real, intent(out) :: dudt
 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp
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

 dudt = (temp1 - temp)*Rg/(gam1*gmw*unit_ergg)/dt

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
subroutine write_options_cooling(iunit)
 use infile_utils, only:write_inopt
 use h2cooling,    only:write_options_h2cooling
 use part,         only:h2chemistry
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 if (h2chemistry) then
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=on)',iunit)
    if (icooling > 0) then
       call write_options_h2cooling(iunit)
    endif
 else
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=explicit, 2=Townsend table, 3=Gammie, 5=KI02)',iunit)
    select case(icooling)
    case(1)
       call write_inopt(icool_radiation_H0,'icool_radiation_H0','H0 cooling on/off',iunit)
       call write_inopt(icool_relax_bowen,'icool_relax_bowen','Bowen (diffusive) relaxation on/off',iunit)
       call write_inopt(icool_relax_stefan,'icool_relax_stefan','radiative relaxation on/off',iunit)
       call write_inopt(icool_dust_collision,'icool_dust_collision','dust collision on/off',iunit)
       call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
    case(2)
       call write_inopt(cooltable,'cooltable','data file containing cooling function',iunit)
       call write_inopt(habund,'habund','Hydrogen abundance assumed in cooling function',iunit)
       call write_inopt(temp_floor,'temp_floor','Minimum allowed temperature in K',iunit)
    case(3)
       call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)
    end select
 endif
 if (icooling > 0) call write_inopt(Tfloor,'Tfloor','temperature floor (K); on if > 0',iunit)

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads sink particle options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(name,valstring,imatch,igotall,ierr)
 use part,         only:h2chemistry
 use h2cooling,    only:read_options_h2cooling
 use io,           only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 logical :: igotallh2,igotallcf

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 igotallh2 = .true.
 igotallcf = .true.
 select case(trim(name))
 case('icooling')
    read(valstring,*,iostat=ierr) icooling
    ngot = ngot + 1
 case('icool_radiation_H0')
    read(valstring,*,iostat=ierr) icool_radiation_H0
    ngot = ngot + 1
 case('icool_relax_bowen')
    read(valstring,*,iostat=ierr) icool_relax_bowen
    ngot = ngot + 1
 case('icool_relax_stefan')
    read(valstring,*,iostat=ierr) icool_relax_stefan
    ngot = ngot + 1
 case('icool_dust_collision')
    read(valstring,*,iostat=ierr) icool_dust_collision
    ngot = ngot + 1
 case('C_cool')
    read(valstring,*,iostat=ierr) C_cool
    ngot = ngot + 1
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
 case('beta_cool')
    read(valstring,*,iostat=ierr) beta_cool
    ngot = ngot + 1
    if (beta_cool < 1.) call fatal('read_options','beta_cool must be >= 1')
 case('Tfloor')
    ! not compulsory to read in
    read(valstring,*,iostat=ierr) Tfloor
 case default
    imatch = .false.
    if (h2chemistry) then
       call read_options_h2cooling(name,valstring,imatch,igotallh2,ierr)
    endif
 end select
 if (icooling == 3 .and. ngot >= 1) igotall = .true.
 if (icooling == 2 .and. ngot >= 3) igotall = .true.
 if (icooling == 1 .and. ngot >= 5) igotall = .true.
 if (igotallh2 .and. ngot >= 1) igotall = .true.

end subroutine read_options_cooling

end module cooling
