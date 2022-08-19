!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling
!
! Gas cooling
!  Current options:
!     0 = none
!     1 = Wind cooling                    [implicit/exact]
!     2 = Wind cooling                    [explicit]
!     3 = Gammie cooling                  [explicit]
!     4 = Townsend (2009) cooling tables  [implicit]
!     5 = Koyama & Inutuska (2002)        [explicit]
!     6 = Koyama & Inutuska (2002)        [implicit]
!
! :References:
!   Koyama & Inutsuka (2002), ApJL 564, 97-100
!   Vazquez-Semadeni, et.al (2007), ApJ 657, 870-883
!   Townsend (2009), ApJS 181, 391-397
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Lionel Siess, Ward Homan
!
! :Runtime parameters:
!   - C_cool         : *factor controlling cooling timestep*
!   - Tfloor         : *temperature floor (K); on if > 0*
!   - beta_cool      : *beta factor in Gammie (2001) cooling*
!   - bowen_Cprime   : *radiative cooling rate (g.s/cm³)*
!   - cooltable      : *data file containing cooling function*
!   - dust_collision : *dust collision [1=on/0=off]*
!   - excitation_HI  : *cooling via electron excitation of HI [1=on/0=off]*
!   - habund         : *Hydrogen abundance assumed in cooling function*
!   - icooling       : *cooling function (0=off, 1,2=wind cooling, 4=Townsend table, 3=Gammie, 5,6=KI02)*
!   - relax_bowen    : *Bowen (diffusive) relaxation [1=on/0=off]*
!   - relax_stefan   : *radiative relaxation [1=on/0=off]*
!   - temp_floor     : *Minimum allowed temperature in K for Townsend cooling table*
!
! :Dependencies: chem, cooling_molecular, datafiles, dim, eos, h2cooling,
!   infile_utils, io, options, part, physcon, timestep, units
!

 use options,  only:icooling
 use timestep, only:C_cool
 use physcon

 implicit none
 character(len=*), parameter :: label = 'cooling'

 public :: init_cooling,calc_cooling_rate,energ_cooling
 public :: write_options_cooling, read_options_cooling, print_cooling_rates
 public :: find_in_table, implicit_cooling, exact_cooling
 public :: cool_dust_discrete_contact, cool_coulomb, &
           cool_HI, cool_H_ionisation, cool_He_ionisation, &
           cool_H2_rovib, cool_H2_dissociation, cool_CO_rovib, &
           cool_H2O_rovib, cool_OH_rot, heat_dust_friction, &
           heat_dust_photovoltaic_soft, heat_CosmicRays, &
           heat_H2_recombination, calc_Q, calc_dlnQdlnT, &
           cool_dust_full_contact, cool_dust_radiation, &
           heat_dust_photovoltaic_hard
 logical, public :: cooling_in_step  = .true.
 real,    public :: bowen_Cprime     = 3.000d-5
 real,    public :: GammaKI_cgs      = 2.d-26        ! [erg/s] heating rate for Koyama & Inutuska cooling
 real,    public :: lambda_shock_cgs = 1.d0
 real,    public :: T1_factor = 20., T0_value = 0.

 private
 integer, parameter :: nTg  = 64
 integer, parameter :: maxt = 1000
 real,    parameter :: Tref = 1.d5, T_floor = 10.   ! required for exact_cooling
 integer :: nt
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt),rhov4_KI02(2,maxt)
 real    :: beta_cool  = 3.
 real    :: habund     = 0.7
 real    :: temp_floor = 1.e4                       ! required for cooling_Townsend_table
 real    :: Tgrid(nTg)
 real    :: LambdaKI_coef,GammaKI
 real    :: KI02_rho_min_cgs = 1.0d-30              ! minimum density of the KI02 cooling curve
 real    :: KI02_rho_max_cgs = 1.0d-14              ! maximum density of the KI02 cooling curve
 real    :: KI02_rho_min,KI02_rho_max
 real    :: kappa_dust_min = 1e-3                   ! dust opacity value below which dust cooling is not calculated
 integer, public :: excitation_HI = 0, relax_Bowen = 0, dust_collision = 0, relax_Stefan = 0, shock_problem = 0
 integer, public :: icool_method  = 0
 character(len=120) :: cooltable = 'cooltable.dat'
 !--Minimum temperature (failsafe to prevent u < 0); optional for ALL cooling options
 real,    public :: Tfloor = 0.                     ! [K]; set in .in file.  On if Tfloor > 0.
 real,    public :: ufloor = 0.                     ! [code units]; set in init_cooling

contains

!-----------------------------------------------------------------------
!+
!  Initialise cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling(id,master,iprint,ierr)

 use dim,               only:maxvxyzu,h2chemistry
 use units,             only:utime,umass,udist,unit_ergg
 use physcon,           only:mass_proton_cgs,kboltz
 use io,                only:fatal
 use eos,               only:gamma,gmw
 use h2cooling,         only:init_h2cooling
 use chem,              only:init_chem
 use cooling_molecular, only:do_molecular_cooling,init_cooling_molec

 integer, intent(in)  :: id,master,iprint
 integer, intent(out) :: ierr

 if (h2chemistry) then
    if (id==master) write(iprint,*) 'initialising cooling function...'
    call init_chem()
    call init_h2cooling()
 else
    !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
    if (relax_Bowen == 1 .and. relax_Stefan == 1) then
       call fatal(label,'you can"t have bowen and stefan cooling at the same time')
    endif

#ifdef KROME
    !krome calculates its own cooling rate
    excitation_HI  = 0
    dust_collision = 0
#else
    if ( icooling == 1 .or. icooling == 2 ) then
    !if no cooling flag activated, disable cooling
       if ( (excitation_HI+relax_Bowen+dust_collision+relax_Stefan+shock_problem) == 0 .and. &
            .not. do_molecular_cooling) then
          print *,'no cooling prescription activated, reset icooling = 0'
          icooling     = 0
          icool_method = 0
          return
       elseif ( do_molecular_cooling) then
       ! Initialise cooling tables
          call init_cooling_molec
       endif
    endif
#endif

    !--initialise remaining variables
    if (icooling == 4) then
       call init_Townsend_table(ierr)
    elseif (icooling == 5 .or. icooling == 6) then
       LambdaKI_coef = GammaKI_cgs*umass*utime**3/(mass_proton_cgs**2 * udist**5)
       GammaKI       = GammaKI_cgs*utime**3/(mass_proton_cgs*udist**2)
       call init_hv4table(ierr)
       if (ierr > 0) call fatal('init_cooling','Failed to create KI02 cooling table')
    elseif (icool_method == 1) then
       call set_Tgrid
    endif
 endif

 !--Determine if this is implicit or explicit cooling
 if (icooling > 0) then
    if (h2chemistry) then
       cooling_in_step = .true.         ! cooling is calculated implicitly in step
    else
       select case (icooling)
       ! cooling is calculated explicitly in force
       case (2,3,5)
          cooling_in_step = .false.
       ! cooling is calculated implicitly in step
       case default
          cooling_in_step = .true.
       end select
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
subroutine init_Townsend_table(ierr)

 use io,        only:fatal
 use datafiles, only:find_phantom_datafile

 integer, intent(out) :: ierr

 integer, parameter :: iu = 127
 integer            :: i
 character(len=120) :: filepath

 !
 ! read the cooling table from file
 !
 filepath = find_phantom_datafile(cooltable,'cooling')
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
 slope(nt)   = slope(nt-1)

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
end subroutine init_Townsend_table

!-----------------------------------------------------------------------
!+
!  create a h-v4 table based upon the cooling curve of KI02
!+
!-----------------------------------------------------------------------
subroutine init_hv4table(ierr)

 use part,    only:hrho,igas
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_velocity
 use eos,     only:gmw,gamma

 integer, intent(out) :: ierr

 integer              :: i,ctr
 real                 :: nrho0_min,nrho0_max,nrho,dnrho,dGammaKI,Lambda,dLambda
 real                 :: T,Tnew,Trat,fatT,faTdT
 logical              :: iterate
 logical              :: print_cc = .false. ! Print the cooling curve (for testing)

 !--Initialise densities
 KI02_rho_min = KI02_rho_min_cgs/unit_density
 KI02_rho_max = KI02_rho_max_cgs/unit_density
 nrho0_min    = KI02_rho_min_cgs/mass_proton_cgs
 nrho0_max    = KI02_rho_max_cgs/mass_proton_cgs
 dnrho        = (log10(nrho0_max) - log10(nrho0_min))/maxt
 !--Initialise additional variables
 dGammaKI     = 0.0
 ierr         = 0

 if (print_cc) open(unit=1031,file='coolingcurve.dat')

 !--Iterate (in cgs units)!
 T = 20000.
 do i = 1,maxt
    ctr     = 0
    iterate = .true.
    nrho    = 10**(log10(nrho0_min) + (i-1)*dnrho)
    do while ( iterate )
       Lambda  = 1.d7*exp(-1.184d5/(T+1.d3)) + 0.014*sqrt(T)*exp(-92./T) ! This is actually Lamda / Gamma
       dLambda = 0.007*exp(-92./T)*(T+184.)*T**(-1.5) + 1.184d12*exp(-1.184d5/(T+1.d3))*(T+1.d3)**(-2)
       fatT    =  Lambda*GammaKI_cgs*nrho -  GammaKI_cgs
       faTdT   = dLambda*GammaKI_cgs*nrho - dGammaKI
       Tnew    = abs(T - fatT/faTdT)
       Trat    = abs( 1.0 - T/Tnew )
       T       = Tnew
       ctr     = ctr + 1
       !--converged
       if (Trat < 1.0d-6) iterate = .false.
       !--failed to converge
       if (T < 0. .or. ctr > 2000) then
          iterate = .false.
          ierr    = 1
       endif
    enddo
    if (print_cc) write(1031,*) nrho,nrho*mass_proton_cgs,T,T*nrho,Lambda*GammaKI_cgs
    rhov4_KI02(1,i) = nrho
    rhov4_KI02(2,i) = T
 enddo
 if (print_cc) close(1031)

 !--Convert to useful values
 do i = 1,maxt
    rhov4_KI02(1,i) = rhov4_KI02(1,i)*mass_proton_cgs/unit_density                               ! number density (cm^-3) -> mass density (code units)
    rhov4_KI02(2,i) = kboltz*rhov4_KI02(2,i)/(gmw*mass_proton_cgs*(gamma-1.0))/unit_velocity**2  ! T -> internal energy (code units)
 enddo

end subroutine init_hv4table

!-----------------------------------------------------------------------
!
!  calculate cooling rates
!
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Teq, mu, gamma, K2, kappa)
 use units,             only:unit_ergg,unit_density
 !use cooling_molecular, only:do_molecular_cooling,calc_cool_molecular

 real, intent(in)  :: rho, T, Teq     !rho in code units
 real, intent(in)  :: mu, gamma
 real, intent(in)  :: K2, kappa       !cgs
 real, intent(out) :: Q, dlnQ_dlnT    !code units

 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, Q_molec, Q_shock
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan, dlnQ_molec, dlnQ_shock
 real :: rho_cgs, ndens

 rho_cgs           = rho*unit_density
 ndens             = rho_cgs/mass_proton_cgs

 Q_H0              = 0.
 Q_relax_Bowen     = 0.
 Q_col_dust        = 0.
 Q_relax_Stefan    = 0.
 Q_shock           = 0.
 Q_molec           = 0.

 dlnQ_H0           = 0.
 dlnQ_relax_Bowen  = 0.
 dlnQ_col_dust     = 0.
 dlnQ_relax_Stefan = 0.
 dlnQ_shock        = 0.
 dlnQ_molec        = 0.

 if (excitation_HI  == 1) call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (relax_Bowen    == 1) call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, gamma, &
                                                        Q_relax_Bowen, dlnQ_relax_Bowen)
 if (dust_collision == 1) call cooling_dust_collision(T, Teq, rho_cgs, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (relax_Stefan   == 1) call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
 if (shock_problem  == 1) call piecewise_law(T, T0_value, ndens, Q_H0, dlnQ_H0)
 !if (do_molecular_cooling) call calc_cool_molecular(T, r, rho_cgs, Q_molec, dlnQ_molec)

 Q_cgs = Q_H0 + Q_relax_Bowen + Q_col_dust + Q_relax_Stefan + Q_molec + Q_shock
 if (Q_cgs == 0.) then
    dlnQ_dlnT = 0.
 else
    dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen + Q_col_dust*dlnQ_col_dust&
   + Q_relax_Stefan*dlnQ_relax_Stefan + Q_molec*dlnQ_molec + Q_shock*dlnQ_shock)/Q_cgs
 endif
 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q         = Q_cgs/unit_ergg

 !call testfunc()
 !call exit

end subroutine calc_cooling_rate

!-----------------------------------------------------------------------
!+
!  Piecewise cooling law for simple shock problem (Creasey et al. 2011)
!+
!-----------------------------------------------------------------------
subroutine piecewise_law(T, T0, ndens, Q, dlnQ)

 real, intent(in)  :: T, T0, ndens
 real, intent(out) :: Q,dlnQ
 real :: T1,Tmid !,dlnT,fac

 T1 = T1_factor*T0
 Tmid = 0.5*(T0+T1)
 if (T < T0) then
    Q    = 0.
    dlnQ = 0.
 elseif (T >= T0 .and. T <= Tmid) then
    !dlnT = (T-T0)/(T0/100.)
    Q = -lambda_shock_cgs*ndens**2*(T-T0)/T0
    !fac = 2./(1.d0 + exp(dlnT))
    dlnQ = 1./(T-T0+1.d-10)
 elseif (T >= Tmid .and. T <= T1) then
    Q = -lambda_shock_cgs*ndens**2*(T1-T)/T0
    dlnQ = -1./(T1-T+1.d-10)
 else
    Q    = 0.
    dlnQ = 0.
 endif
 !derivatives are discontinuous!

end subroutine piecewise_law

!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling prescription
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Teq, rho, mu, gamma, Q, dlnQ_dlnT)

 use physcon, only:Rg

 real, intent(in)  :: T, Teq, rho, mu, gamma
 real, intent(out) :: Q,dlnQ_dlnT

 Q         = Rg/((gamma-1.)*mu)*rho*(Teq-T)/bowen_Cprime
 dlnQ_dlnT = -T/(Teq-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_dust_collision(T, Teq, rho, K2, mu, Q, dlnQ_dlnT)

 use physcon, only: kboltz, mass_proton_cgs, pi

 real, intent(in)  :: T, Teq, rho, K2, mu
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter   :: f = 0.15, a0 = 1.28e-8
 real              :: A

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

 Q         = 4.*steboltz*(Teq**4-T**4)*kappa
 dlnQ_dlnT = -4.*T**4/(Teq**4-T**4+1.d-10)

end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to electron excitation of neutral H (Spitzer 1978)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen(T, rho_cgs, Q, dlnQ_dlnT)

 use physcon, only: mass_proton_cgs, pi

 real, intent(in)  :: T, rho_cgs
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter   :: f = 1.0d0
 real              :: eps_e

 if (T > 3000.) then
    eps_e = calc_eps_e(T)
    Q = -f*7.3d-19*eps_e*exp(-118400./T)*rho_cgs/(1.4*mass_proton_cgs)**2
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

 real             :: k1, k2, k3, k8, k9, p, q

 k1 = 1.88d-10 / T**6.44e-1
 k2 = 1.83d-18 * T
 k3 = 1.35d-9
 k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
 k9 = 1.7d-4 * k8
 p  = .5*k8/k9
 q  = k1*(k2+k3)/(k3*k9)
 calc_eps_e = (p + sqrt(q+p**2))/q

end function calc_eps_e

!-----------------------------------------------------------------------
!+
!  Set Temperature grid
!+
!-----------------------------------------------------------------------
subroutine set_Tgrid

 integer :: i
 real    :: dlnT

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
!   This is for the explicit calculation
!   In equilibrium, n*LambdaKI = (rho/mp)*LambdaKI = GammaKI
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska_explicit(rhoi,Tgas,dudti)

 real, intent(in)    :: rhoi,Tgas
 real, intent(inout) :: dudti

 real                :: LambdaKI

 ! Derivation to obtain correct units; used Koyama & Inutuska (2002) as the reference
 !LambdaKI = GammaKI_cgs * (1.d7*exp(-118400./(Tgas+1000))+0.014*sqrt(Tgas)*exp(-92./Tgas)) ! The cooling rate in erg cm^3/s = g cm^5/s^3
 !LambdaKI = LambdaKI/mass_proton_cgs**2                                                    ! units are now cm^5/(g s^3) ! since [u] = erg/g = cm^2/s^2
 !LambdaKI = LambdaKI*umass*utime**3/udist**5                                               ! convert to from cm^5/(g s^3) to code units
 !dudti    = dudti - LambdaKI*rhoi*fac                                                      ! multiply by rho (code) to get l^5/(m t^3) * m/l^3 = l^2/s^3 = [u]
 !
 !GammaKI = GammaKI_cgs                                                                     ! The heating rate in erg /s = g cm^2/s^3
 !GammaKI = GammaKI/mass_proton_cgs                                                         ! divide by proton mass.  Units are now g cm^2 / s^3 / g = cm^2/s^3
 !GammaKI = GammaKI*utime**3/udist**2                                                       ! convert from cm^2/s^3 to code units
 !dudti   = dudti + GammaKI                                                                 ! units and dependencies are correct

 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(Tgas+1000.))+0.014*sqrt(Tgas)*exp(-92./Tgas))
 dudti    = dudti - LambdaKI*rhoi + GammaKI

end subroutine cooling_KoyamaInutuska_explicit

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is the implicit method given by (5)-(6) in Vazquez-Semadeni+ (2007)
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska_implicit(eni,rhoi,dt,dudti)

 use eos, only: gamma,temperature_coef,gmw

 real, intent(in)    :: rhoi,eni,dt
 real, intent(out)   :: dudti

 integer             :: i,j,jm1
 real                :: ponrhoi,tempi,eni_equil,eni_final,deni,tau1,LambdaKI

 !--Determine the indicies surrounding the input h
 i = minloc(abs(rhov4_KI02(1,1:maxt)-rhoi), 1)
 if (i==1) then
    !print*, 'min density too large! extrapolating using two smallest densities'
    j = 2
 elseif (i==maxt) then
    !print*, 'max density too small! extrapolating using two largest densities'
    j = maxt
 elseif (rhov4_KI02(1,i-1) <= rhoi .and. rhoi <= rhov4_KI02(1,i  )) then
    j = i
 elseif (rhov4_KI02(1,i  ) <= rhoi .and. rhoi <= rhov4_KI02(1,i+1)) then
    j = i+1
 else
    print*, rhoi,rhov4_KI02(1,i-1:i+1)
    print*, 'this should not happen'
    stop
 endif

 !--Calculate the equilibrium energy by linear interpolation
 jm1       = j - 1
 eni_equil = rhov4_KI02(2,j) + (rhov4_KI02(2,jm1)-rhov4_KI02(2,j))/(rhov4_KI02(1,jm1)-rhov4_KI02(1,j))*(rhoi-rhov4_KI02(1,j))

 !--Determine the inverse time require to radiate/acquire excess/deficit energy & Update energy
 ponrhoi  = (gamma-1.)*eni
 tempi    = temperature_coef*gmw*ponrhoi
 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(tempi+1000.))+0.014*sqrt(tempi)*exp(-92./tempi))
 dudti    = LambdaKI*rhoi - GammaKI
 deni     = eni - eni_equil

 if (abs(deni) > 0.) then
    ! in both limits, this will approach the correct value
    tau1      = abs(dudti/deni)
    eni_final = eni_equil + deni*exp(-dt*tau1)
    dudti     = -(eni - eni_final)/dt
 else
    ! in the unlikly chance deni = 0
    dudti = -dudti
 endif

end subroutine cooling_KoyamaInutuska_implicit

!-----------------------------------------------------------------------
!
!   explicit cooling
!
!-----------------------------------------------------------------------
subroutine explicit_cooling (ui, dudt, rho, dt, mu, gamma, Trad, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, Trad, mu, gamma !code units
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt                                     !code units

 real              :: u,Q,dlnQ_dlnT,T,T_on_u

 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui
 call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, gamma, K2, kappa)
 if (ui - Q*dt < 0.) then   ! assume thermal equilibrium
    u    = Trad/T_on_u      ! set T=Trad
    dudt = (u-ui)/dt
 else
    dudt = Q
 endif

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!
!   implicit cooling
!
!-----------------------------------------------------------------------
subroutine implicit_cooling (ui, dudt, rho, dt, mu, gamma, Trad, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, mu, gamma
 real, intent(in)  :: Trad, K2, kappa
 real, intent(out) :: dudt

 real, parameter    :: tol      = 1.d-4    ! to be adjusted
 integer, parameter :: iter_max = 200
 real               :: u,Q,dlnQ_dlnT,T,T_on_u,delta_u,term1,term2
 integer            :: iter

 u       = ui
 T_on_u  = (gamma-1.)*mu*unit_ergg/Rg
 delta_u = 1.d-3
 iter    = 0
 !term1 = 1.-(gamma-1.)*dt*divcurlv !pdv=(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
 term1 = 1.
 do while (abs(delta_u) > tol .and. iter < iter_max)
    T       = u*T_on_u
    call calc_cooling_rate(Q,dlnQ_dlnT, rho, T, Trad, mu, gamma, K2, kappa)
    term2   = u*term1-Q*dt
    delta_u = (ui-term2)/(term1-Q*dlnQ_dlnT*dt/u)
    u       = u+delta_u
    iter    = iter + 1
 enddo
 dudt =(u-ui)/dt
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
subroutine energ_cooling(xi,yi,zi,ui,dudt,rho,dt,Trad_in,mu_in,gamma_in,K2_in,kappa_in,Tgas_in)
 use io,      only:fatal
 use eos,     only:gmw,gamma
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in)           :: xi,yi,zi,ui,rho,dt                  ! in code units
 real, intent(in), optional :: Tgas_in,Trad_in,mu_in,gamma_in,K2_in,kappa_in   ! in cgs
 real, intent(out)          :: dudt                                ! in code units
 real                       :: mu,polyIndex,T_on_u,Tgas,Trad,K2,kappa

 dudt       = 0.
 mu         = gmw
 polyIndex  = gamma
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 Tgas   = T_on_u*ui
 Trad   = Tgas
 kappa  = 0.
 K2     = 0.
 if (present(gamma_in)) polyIndex = gamma_in
 if (present(mu_in))    mu        = mu_in
 if (present(Trad_in))  Trad      = Trad_in
 if (present(Tgas_in))  Tgas      = Tgas_in
 if (present(K2_in))    K2        = K2_in
 if (present(kappa_in)) kappa     = kappa_in

 select case (icooling)
 case(1,2)
    if (icool_method == 2) then
       call exact_cooling   (ui,dudt,rho,dt,mu,polyIndex,Trad,K2,kappa)
    elseif (icool_method == 0) then
       call implicit_cooling(ui,dudt,rho,dt,mu,polyIndex,Trad,K2,kappa)
    else
       call explicit_cooling(ui,dudt,rho,dt,mu,polyIndex,Trad,K2,kappa)
    endif
 case (3)
    call cooling_Gammie(xi,yi,zi,ui,dudt)
 case (4)
    call cooling_Townsend_table(ui,rho,dt,dudt,mu,polyIndex)
 case (5)
    if (present(Tgas_in)) then
       call cooling_KoyamaInutuska_explicit(rho,Tgas,dudt)
    else
       call fatal('energ_cooling','Koyama & Inutuska cooling requires gas temperature')
    endif
 case (6)
    call cooling_KoyamaInutuska_implicit(ui,rho,dt,dudt)
 case default
    call implicit_cooling(ui,dudt,rho,dt,mu,polyIndex,Trad,K2,kappa)
 end select

end subroutine energ_cooling

!-----------------------------------------------------------------------
!
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!
!-----------------------------------------------------------------------
subroutine exact_cooling(ui, dudt, rho, dt, mu, gamma, Trad, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg

 real, intent(in)  :: ui, rho, dt, Trad, mu, gamma
 real, intent(in)  :: K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u
 integer         :: k

 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui

 if (T < T_floor) then
    Temp = T_floor
 elseif (T > Tref) then
    call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, gamma, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    call calc_cooling_rate(Qref,dlnQref_dlnT, rho, Tref, Trad, mu, gamma, K2, kappa)
    Y         = 0.
    k         = nTg
    Q         = Qref          ! default value if Tgrid < T for all k
    dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
    do while (Tgrid(k) > T)
       k = k-1
       call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Trad, mu, gamma, K2, kappa)
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
    y  = y + dy
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
subroutine cooling_Townsend_table(uu,rho,dt,dudt,mu,gamma)

 use physcon, only:atomic_mass_unit,kboltz,Rg
 use units,   only:unit_density,unit_ergg,utime

 real, intent(in)  :: uu, rho,dt,mu,gamma
 real, intent(out) :: dudt

 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp
 real    :: sloperef,slopek,temp,temp1,tref,yfunx,yinv0
 integer :: k

 gam1 = gamma - 1.
 temp = gam1*uu/Rg*mu*unit_ergg

 tref     = temper(nt)
 sloperef = slope(nt)

 if (temp < temp_floor) then
    temp1 = temp_floor
 else
    amue        = 2.*atomic_mass_unit/(1. + habund)
    amuh        = atomic_mass_unit/habund
    density_cgs = rho*unit_density
    dt_cgs      = dt*utime

    !Lionel Siess : I think there is an error. in dtemp sloperef should be replaced by lambda(nt)
    !original dtemp = gam1*density_cgs*(atomic_mass_unit*mu/(amue*amuh*kboltz))* &
    !     sloperef/tref*dt_cgs
    ! Eq 26
    dtemp = gam1*density_cgs*(atomic_mass_unit*mu/(amue*amuh*kboltz))*lambda(nt)/tref*dt_cgs

    k     = find_in_table(nt,temper,temp)
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

 dudt = (temp1 - temp)*Rg/(gam1*mu*unit_ergg)/dt

end subroutine cooling_Townsend_table

!-----------------------------------------------------------------------
!+
!  utility to find the index of closest value in a table
!+
!-----------------------------------------------------------------------
pure integer function find_in_table(n,table,val) result(i)

 integer, intent(in) :: n
 real,    intent(in) :: table(n), val

 integer             :: i0,i1

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
 use cooling_molecular, only: write_options_molecularcooling
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 if (h2chemistry) then
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=on)',iunit)
    if (icooling > 0) call write_options_h2cooling(iunit)
 else
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=cooling in substep, 2=cooling in force,'// &
                              '3=Gammie, 4=Townsend table, 5,6=KI02)',iunit)
    select case(icooling)
    case(1,2)
       !call write_options_molecularcooling(iunit)
       call write_inopt(icool_method,'icool_method',&
            'integration method (0=implicit, 1=explicit, 2=exact solution)',iunit)
       call write_inopt(excitation_HI,'excitation_HI','cooling via electron excitation of HI (1=on/0=off)',iunit)
       call write_inopt(relax_bowen,'relax_bowen','Bowen (diffusive) relaxation (1=on/0=off)',iunit)
       call write_inopt(relax_stefan,'relax_stefan','radiative relaxation (1=on/0=off)',iunit)
       call write_inopt(dust_collision,'dust_collision','dust collision (1=on/0=off)',iunit)
       call write_inopt(shock_problem,'shock_problem','piecewise formulation for analytic shock solution (1=on/0=off)',iunit)
       if (shock_problem == 1) then
          call write_inopt(lambda_shock_cgs,'lambda_shock','Cooling rate parmaeter for analytic shock solution',iunit)
          call write_inopt(T1_factor,'T1_factor','factor by which T0 is increased (T1= T1_factor*T0)',iunit)
          call write_inopt(T0_value,'T0','temperature to cool towards',iunit)
       endif
       call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
    case(3)
       call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)
    case(4)
       call write_inopt(cooltable,'cooltable','data file containing cooling function',iunit)
       call write_inopt(habund,'habund','Hydrogen abundance assumed in cooling function',iunit)
       call write_inopt(temp_floor,'temp_floor','Minimum allowed temperature in K for Townsend cooling table',iunit)
    end select
 endif
 if (icooling > 0) call write_inopt(Tfloor,'Tfloor','temperature floor (K); on if > 0',iunit)

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(name,valstring,imatch,igotall,ierr)
 use part,              only:h2chemistry
 use h2cooling,         only:read_options_h2cooling
 use io,                only:fatal
 use cooling_molecular, only:read_options_molecular_cooling
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0, nn = 9
 logical :: igotallh2,igotallcf,igotallmol

 imatch     = .true.
 igotall    = .false.  ! cooling options are compulsory
 igotallh2  = .true.
 igotallcf  = .true.
 igotallmol = .true.

 select case(trim(name))
 case('icooling')
    read(valstring,*,iostat=ierr) icooling
    ngot = ngot + 1
 case('icool_method')
    read(valstring,*,iostat=ierr) icool_method
    ngot = ngot + 1
 case('excitation_HI')
    read(valstring,*,iostat=ierr) excitation_HI
    ngot = ngot + 1
 case('relax_bowen')
    read(valstring,*,iostat=ierr) relax_bowen
    ngot = ngot + 1
 case('relax_stefan')
    read(valstring,*,iostat=ierr) relax_stefan
    ngot = ngot + 1
 case('dust_collision')
    read(valstring,*,iostat=ierr) dust_collision
    ngot = ngot + 1
 case('shock_problem')
    read(valstring,*,iostat=ierr) shock_problem
    ngot = ngot + 1
 case('lambda_shock')
    read(valstring,*,iostat=ierr) lambda_shock_cgs
    ngot = ngot + 1
 case('T1_factor')
    read(valstring,*,iostat=ierr) T1_factor
    ngot = ngot + 1
 case('T0')
    read(valstring,*,iostat=ierr) T0_value
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
 ierr = 0
 if (shock_problem == 1) then
    nn = 12
 else
    nn = 9
 endif
 !if (icooling > 0 .and. .not. imatch) call read_options_molecular_cooling(name,valstring,imatch,igotallmol,ierr)
 if (icooling == 0 .and. ngot >= 2) igotall = .true.
 if (icooling == 1 .and. ngot >= nn) igotall = .true.
 if (icooling == 2 .and. ngot >= 8) igotall = .true.
 if (icooling == 3 .and. ngot >= 1) igotall = .true.
 if (icooling == 4 .and. ngot >= 3) igotall = .true.
 if (h2chemistry .and. igotallh2 .and. ngot >= 1) igotall = .true.

end subroutine read_options_cooling








!=======================================================================
!=======================================================================
!=======================================================================
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!
!  Test routine for cooling functions
!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!=======================================================================
!=======================================================================
!=======================================================================

subroutine testfunc()

 use physcon, only: mass_proton_cgs

 real :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL
 real :: n_gas

 ! evaluate parameters using plausible values to test the cooling functions

 T_gas      = 4000.
 rho_gas    = 2.d-16
 mu         = 1.78

 n_gas      = rho_gas/(mu*mass_proton_cgs)

 nH         = 0.5 *n_gas
 nH2        = 0.5 *n_gas
 nHe        = 0.1 *n_gas
 nCO        = 1.d-4*n_gas
 nH2O       = 5.d-5*n_gas
 nOH        = 1.d-7*n_gas
 kappa_gas  = 2.d-4

 T_dust     = 1500.
 v_drift    = 1.d6
 d2g        = 1./200.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_dust = 1.d-4

 JL = 2.5d-12     ! Value taken from Barstow et al. 1997
 call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

end subroutine testfunc

subroutine print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real, intent(in)  :: JL
 real :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Qtot, dlnQ_dlnT, nH_tot

 !nH_tot = nH+2.*nH2
 nH_tot = 1.
 print*, ' '
 print*, ' '
 print*, '-----------------------------------------------------------------------------'
 print*, ' '
 print*, ' '
 Q1  = cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust)
 print*, 'Q1   = ', Q1/nH_tot
 Q2  = cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust)
 print*, 'Q2   = ', Q2/nH_tot
 Q3  = cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)
 print*, 'Q3   = ', Q3/nH_tot
 Q4  = cool_coulomb(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q4   = ', Q4/nH_tot
 Q5  = cool_HI(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q5   = ', Q5/nH_tot
 Q6  = cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q6   = ', Q6/nH_tot
 Q7  = cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q7   = ', Q7/nH_tot
 Q8  = cool_H2_rovib(T_gas, nH, nH2)
 print*, 'Q8   = ', Q8/nH_tot
 Q9  = cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2)
 print*, 'Q9   = ', Q9/nH_tot
 Q10 = cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)
 print*, 'Q10  = ', Q10/nH_tot
 Q11 = cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O)
 print*, 'Q11  = ', Q11/nH_tot
 Q12 = cool_OH_rot(T_gas, rho_gas, mu, nOH)
 print*, 'Q12  = ', Q12/nH_tot
 Q13 = heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust)
 print*, 'Q13  = ', Q13/nH_tot
 Q14 = heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust)
 print*, 'Q14  = ', Q14/nH_tot
 Q15 = heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL)
 print*, 'Q15  = ', Q15/nH_tot
 Q16 = heat_CosmicRays(nH, nH2)
 print*, 'Q16  = ', Q16/nH_tot
 Q17 = heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)
 print*, 'Q17  = ', Q17/nH_tot

 Qtot = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 print*, 'Qtot = ', Qtot/nH_tot

 dlnQ_dlnT = calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 print*, 'dlnQdlnT = ', dlnQ_dlnT

 print*, ' '
 print*, ' '
 print*, '------------------- exit in calc_cooling_rate --------------------------------'
 print*, ' '
 print*, ' '

end subroutine print_cooling_rates






!=======================================================================
!=======================================================================
!=======================================================================
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!
!  Physics required for cooling functions below
!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!=======================================================================
!=======================================================================
!=======================================================================






!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute LTE electron density from SAHA equations
!                      (following D'Angelo & Bodenheimer 2013)
!+
!-----------------------------------------------------------------------
real function n_e(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: kboltz, mass_proton_cgs, mass_electron_cgs, planckhbar, pi

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe

 real, parameter  :: H_ion   = 2.179d-11    ! 13.60 eV in erg
 real, parameter  :: He_ion  = 3.940d-11    ! 24.59 eV in erg
 real, parameter  :: He2_ion = 8.720d-11    ! 54.42 eV in erg
 real             :: n_gas, X, KH, xx, Y, KHe, KHe2, z1, z2, cst

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 X      = nH /n_gas
 Y      = nHe/n_gas
 cst    = mass_proton_cgs/rho_gas * sqrt(mass_electron_cgs*kboltz*T_gas/(2.*pi*planckhbar**2))**3
 if (T_gas > 1.d5) then
   xx = 1.
 else
   KH   = cst/X * exp(-H_ion /(kboltz*T_gas))
   ! solution to quadratic SAHA equations (Eq. 16 in D'Angelo et al 2013)
   xx   = (1./2.) * (-KH    + sqrt(KH**2+4.*KH))
 endif
 if (T_gas > 3.d5) then
   z1 = 1.
   z2 = 1.
 else
   KHe    = 4.*cst * exp(-He_ion/(kboltz*T_gas))
   KHe2   =    cst * exp(-He2_ion/(kboltz*T_gas))

! solution to quadratic SAHA equations (Eq. 17 in D'Angelo et al 2013)
   z1     = (2./Y ) * (-KHe-X + sqrt((KHe+X)**2+KHe*Y))
! solution to quadratic SAHA equations (Eq. 18 in D'Angelo et al 2013)
   z2     = (2./Y ) * (-KHe2-X + sqrt((KHe+X+Y/4.)**2+KHe2*Y))
 endif
 n_e    = xx * nH + z1*(1.+z2) * nHe

end function n_e

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute mean thermal speed of molecules
!+
!-----------------------------------------------------------------------
real function v_th(T_gas,mu)

 use physcon, only: kboltz, mass_proton_cgs

 real, intent(in) :: T_gas, mu

 v_th = sqrt((3.*kboltz*T_gas)/(mu*mass_proton_cgs))

end function v_th


!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute fraction of gas that has speeds lower than v_crit
!                      from the cumulative distribution function of the
!                      Maxwell-Boltzmann distribution
!+
!-----------------------------------------------------------------------
real function MaxBol_cumul(T_gas, mu,  v_crit)

 use physcon, only: kboltz, mass_proton_cgs, pi

 real, intent(in) :: T_gas, mu, v_crit

 real             :: a

 a            = sqrt( kboltz*T_gas/(mu*mass_proton_cgs) )
 MaxBol_cumul = erf(v_crit/(sqrt(2.)*a)) - sqrt(2./pi) * (v_crit*exp(-v_crit**2/(2.*a**2))) / a

end function MaxBol_cumul


!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute dust number density from dust-to-gas mass ratio,
!                      mean grain size a, and specific density of the grain
!+
!-----------------------------------------------------------------------
real function n_dust(rho_gas, d2g, a, rho_grain)

 use physcon, only: pi

 real, intent(in) ::rho_gas,d2g,a,rho_grain

 n_dust = ( rho_gas*d2g ) / ( (4./3.)*pi*a**3.*rho_grain )

end function n_dust








!=======================================================================
!=======================================================================
!=======================================================================
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!
!  Cooling functions
!
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!=======================================================================
!=======================================================================
!=======================================================================








!-----------------------------------------------------------------------
!+
!  DUST:  Full contact cooling (Bowen 1988)
!+
!-----------------------------------------------------------------------
real function cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust)

 use physcon, only:Rg

 real, intent(in)  :: T_gas, rho_gas, mu
 real, intent(in)  :: T_dust, kappa_dust

 if (kappa_dust > kappa_dust_min) then
    cool_dust_full_contact = (3.*Rg)/(2.*mu*bowen_Cprime)*rho_gas*(T_gas-T_dust)
 else
    cool_dust_full_contact = 0.0
 endif
end function cool_dust_full_contact

!-----------------------------------------------------------------------
!+
!  DUST: Discrete contact cooling (Hollenbach & McKee 1979)
!+
!-----------------------------------------------------------------------
real function cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust)

 use physcon, only: kboltz, mass_proton_cgs, pi

 real, intent(in)  :: T_gas, rho_gas, mu
 real, intent(in)  :: T_dust, d2g, a, rho_grain, kappa_dust

 real, parameter   :: alpha = 0.33  ! See Burke & Hollenbach 1983
 real              :: n_gas, sigma_dust

 if (kappa_dust > kappa_dust_min) then
    sigma_dust                 = 2.*pi*a**2
    n_gas                      = rho_gas/(mu*mass_proton_cgs)
    cool_dust_discrete_contact = alpha*n_gas*n_dust(rho_gas,d2g,a,rho_grain)*sigma_dust*v_th(T_gas,mu)*kboltz*(T_gas-T_dust)
 else
    cool_dust_discrete_contact = 0.0
 endif
end function cool_dust_discrete_contact

!-----------------------------------------------------------------------
!+
!  DUST: Radiative cooling (Woitke 2006) - DO NOT USE, PHYSICALLY INCORRECT
!+
!-----------------------------------------------------------------------
real function cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)

 use physcon, only: steboltz

 real, intent(in)  :: T_gas, kappa_gas
 real, intent(in)  :: T_dust, kappa_dust

 if (kappa_dust > kappa_dust_min) then
    cool_dust_radiation = 4.*steboltz*(kappa_gas*T_gas**4-kappa_dust*T_dust**4)
 else
    cool_dust_radiation = 0.0
 endif
end function cool_dust_radiation

!-----------------------------------------------------------------------
!+
!  DUST: Friction heating caused by dust-drift through gas (Golreich & Scoville 1976)
!+
!-----------------------------------------------------------------------
real function heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust)

 real, intent(in)  :: rho_gas
 real, intent(in)  :: v_drift, d2g, a, rho_grain, kappa_dust

 real              :: sigma_dust
 real, parameter   :: alpha = 0.33                            ! see Burke & Hollenbach 1983

 ! Warning, alpha depends on the type of dust
 if (kappa_dust > kappa_dust_min) then
    sigma_dust         = 2.*pi*a**2
    heat_dust_friction = n_dust(rho_gas,d2g,a,rho_grain)*sigma_dust*v_drift*alpha*0.5*rho_gas*v_drift**2
 else
    heat_dust_friction = 0.0
 endif

end function heat_dust_friction

!-----------------------------------------------------------------------
!+
!  DUST: photovoltaic heating by soft UV field (Weingartner & Draine 2001)
!+
!-----------------------------------------------------------------------
real function heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 real, intent(in)  :: kappa_dust

 real              :: x
 real, parameter   :: G=1.68 ! ratio of true background UV field to Habing field
 real, parameter   :: C0=5.45, C1=2.50, C2=0.00945, C3=0.01453, C4=0.147, C5=0.623, C6=0.511 ! see Table 2 in Weingartner & Draine 2001, last line

 if (kappa_dust > kappa_dust_min) then
    x                           = G*sqrt(T_gas)/n_e(T_gas, rho_gas, mu, nH, nHe)
    heat_dust_photovoltaic_soft = 1.d-26*G*nH*(C0+C1*T_gas**C4)/(1.+C2*x**C5*(1.+C3*x**C6))
 else
    heat_dust_photovoltaic_soft = 0.0
 endif

end function heat_dust_photovoltaic_soft

!-----------------------------------------------------------------------
!+
!  DUST: photovoltaic heating by hard UV field (Inoue & Kamaya 2010)
!+
!-----------------------------------------------------------------------
real function heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL)

 real, intent(in)  :: T_gas, nH
 real, intent(in)  :: d2g, kappa_dust
 real, intent(in)  :: JL       ! mean intensity of background UV radiation at hydrogen Lyman limit (91.2 nm)

 if (kappa_dust > kappa_dust_min) then
    heat_dust_photovoltaic_hard = 1.2d-34*(d2g   /1.d-4) &
                                         *(nH   /1.d-5 )**(4. /3.) &
                                         *(T_gas/1.d4  )**(-1./6.) &
                                         *(JL   /1.d-21)**(2. /3.)
 else
    heat_dust_photovoltaic_hard = 0.0
 endif

end function heat_dust_photovoltaic_hard

!-----------------------------------------------------------------------
!+
!  PARTICLE: Coulomb cooling via electron scattering (Weingartner & Draine 2001)
!+
!-----------------------------------------------------------------------
real function cool_coulomb(T_gas, rho_gas, mu, nH, nHe)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe

 real              :: x, ne
 real, parameter   :: G=1.68 ! ratio of true background UV field to Habing field
 real, parameter   :: D0=0.4255, D1=2.457, D2=-6.404, D3=1.513, D4=0.05343 ! see Table 3 in Weingartner & Draine 2001, last line

 if (T_gas > 1000.) then
    ne = n_e(T_gas, rho_gas, mu, nH, nHe)
    x  = log(G*sqrt(T_gas)/ne)
    cool_coulomb = 1.d-28*ne*nH*T_gas**(D0+D1/x)*exp(D2+D3*x-D4*x**2)
 else
    cool_coulomb = 0.0
 endif

end function cool_coulomb

!-----------------------------------------------------------------------
!+
!  PARTICLE: Cosmic ray heating (Jonkheid et al. 2004)
!+
!-----------------------------------------------------------------------
real function heat_CosmicRays(nH, nH2)

 real, intent(in) :: nH, nH2

 real, parameter  :: Rcr = 5.0d-17  !cosmic ray ionisation rate [s^-1]

 heat_CosmicRays = Rcr*(5.5d-12*nH+2.5d-11*nH2)

end function heat_CosmicRays

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to electron excitation of neutral H (Spitzer 1978, Black 1982, Cen 1992)
!+
!-----------------------------------------------------------------------
real function cool_HI(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe

 real              :: n_gas

 ! all hydrogen atomic, so nH = n_gas
 ! Dalgarno & McCray (1972) provide data starting at 3000K
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 3000.) then
    n_gas   = rho_gas/(mu*mass_proton_cgs)
    cool_HI = 7.3d-19*n_e(T_gas, rho_gas, mu, nH, nHe)*n_gas/(1.+sqrt(T_gas/1.d5))*exp(-118400./T_gas)
 else
    cool_HI = 0.0
 endif

end function cool_HI

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to collisional ionisation of neutral H (Black 1982, Cen 1992)
!+
!-----------------------------------------------------------------------
real function cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe

 real              :: n_gas

 ! all hydrogen atomic, so nH = n_gas
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 4000.) then
    n_gas             = rho_gas/(mu*mass_proton_cgs)
    cool_H_ionisation = 1.27d-21*n_e(T_gas, rho_gas, mu, nH, nHe)*n_gas*sqrt(T_gas)/(1.+sqrt(T_gas/1.d5))*exp(-157809./T_gas)
 else
    cool_H_ionisation = 0.0
 endif

end function cool_H_ionisation

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to collisional ionisation of neutral He (Black 1982, Cen 1992)
!+
!-----------------------------------------------------------------------
real function cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe

 real              :: n_gas

 ! all hydrogen atomic, so nH = n_gas
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 4000.) then
    n_gas              = rho_gas/(mu*mass_proton_cgs)
    cool_He_ionisation = 9.38d-22*n_e(T_gas, rho_gas, mu, nH, nHe)*nHe*sqrt(T_gas)*(1+sqrt(T_gas/1.d5))**(-1)*exp(-285335./T_gas)
 else
    cool_He_ionisation = 0.0
 endif

end function cool_He_ionisation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: Cooling due to ro-vibrational excitation of H2 (Lepp & Shull 1983)
!+
!-----------------------------------------------------------------------
real function cool_H2_rovib(T_gas, nH, nH2)

 real, intent(in)  :: T_gas, nH, nH2

 real              :: kH_01, kH2_01
 real              :: Lvh, Lvl, Lrh, Lrl
 real              :: x, Qn

 if (T_gas < 1635.) then
    kH_01 = 1.4d-13*exp((T_gas/125.)-(T_gas/577.)**2)
 else
    kH_01 = 1.0d-12*sqrt(T_gas)*exp(-1000./T_gas)
 endif
 kH2_01 = 1.45d-12*sqrt(T_gas)*exp(-28728./(T_gas+1190.))
 Lvh    = 1.1d-13*exp(-6744./T_gas)
 Lvl    = 8.18d-13*(nH*kH_01+nH2*kH2_01)

 x   = log10(T_gas/1.0d4)
 if (T_gas < 1087.) then
    Lrh = 10.**(-19.24+0.474*x-1.247*x**2)
 else
    Lrh = 3.9d-19*exp(-6118./T_gas)
 endif

 Qn = nH2**0.77+1.2*nH**0.77
 if (T_gas > 4031.) then
    Lrl = 10.**(-22.9-0.553*x-1.148*x**2)*Qn
 else
    Lrl = 1.38d-22*exp(-9243./T_gas)*Qn
 endif

 cool_H2_rovib = nH2*( Lvh/(1.+(Lvh/Lvl)) + Lrh/(1.+(Lrh/Lrl)) )

end function cool_H2_rovib

!-----------------------------------------------------------------------
!+
!  CHEMICAL: H2 dissociation cooling (Shapiro & Kang 1987)
!+
!-----------------------------------------------------------------------
real function cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2)

 use physcon, only: mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2

 real              :: n_gas
 real              :: x, n1, n2, beta
 real              :: kD_H, kD_H2

 n_gas = rho_gas/(mu*mass_proton_cgs)
 x     = log10(T_gas/1.0d4)
 n1    = 10.**(4.0   -0.416*x -0.327*x**2)
 n2    = 10.**(4.845 -1.3*x   +1.62*x**2)
 beta  = 1./(1.+n_gas*(2.*nH2/n_gas*((1./n2)-(1./n1))+1./n1))
 kD_H  = 1.2d-9*exp(-52400/T_gas)*(0.0933*exp(-17950./T_gas))**beta
 kD_H2 = 1.3d-9*exp(-53300/T_gas)*(0.0908*exp(-16200./T_gas))**beta

 cool_H2_dissociation = 7.18d-12*(nH2**2*kD_H2+nH*nH2*kD_H)

end function cool_H2_dissociation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: H2 recombination heating (Hollenbach & Mckee 1979)
!            for an overview, see Valentine Wakelama et al. 2017
!+
!-----------------------------------------------------------------------
real function heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)

 use physcon, only: mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, T_dust

 real              :: n_gas
 real              :: x, n1, n2, beta
 real              :: xi, fa, k_rec

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 x      = log10(T_gas/1.0d4)
 n1     = 10.**(4.0   -0.416*x -0.327*x**2)
 n2     = 10.**(4.845 -1.3*x   +1.62*x**2)
 beta   = 1./(1.+n_gas*(2.*nH2/n_gas*((1./n2)-(1./n1))+1./n1))
 xi     = 7.18d-12*n_gas*nH*(1.-beta)

 fa     = (1.+1.0d4*exp(-600./T_dust))**(-1.)    ! eq 3.4
 k_rec  = 3.0d-1*(sqrt(T_gas)*fa)/(1.+0.04*sqrt(T_gas+T_dust)+2.0d-3*T_gas+8.0d-6*T_gas**2) ! eq 3.8

 heat_H2_recombination = k_rec*xi

end function heat_H2_recombination

!-----------------------------------------------------------------------
!+
!  RADIATIVE: optically thin CO ro-vibrational cooling (Hollenbach & McKee 1979, McKee et al. 1982)
!+
!-----------------------------------------------------------------------
real function cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)

 use physcon, only: kboltz, mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nCO

 real              :: Qrot, QvibH2, QvibH
 real              :: n_gas, n_crit, sigma
 real              :: v_crit, nfCO

! CO bond dissociation energy = 11.11 eV = 1.78e-11 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy CO

 v_crit = sqrt( 2.*1.78d-11/(mu*mass_proton_cgs) )  ! kinetic energy
 nfCO   = MaxBol_cumul(T_gas, mu,  v_crit) * nCO

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 n_crit = 3.3d6*(T_gas/1000.)**0.75                                                                             !McKee et al. 1982 eq. 5.3
 sigma  = 3.0d-16*(T_gas/1000.)**(-1./4.)                                                                       !McKee et al. 1982 eq. 5.4
 Qrot   = n_gas*nfCO*0.5*(kboltz*T_gas*sigma*v_th(T_gas, mu)) / (1. + (n_gas/n_crit) + 1.5*sqrt(n_gas/n_crit))  !McKee et al. 1982 eq. 5.2

 QvibH2 = 1.83d-26*nH2*nfCO*exp(-3080./T_gas)*exp(-68./(T_gas**(1./3.))) !Neufeld & Kaufman 1993
 QvibH  = 1.28d-24*nH *nfCO*exp(-3080./T_gas)*exp(-(2000./T_gas)**3.43)  !Neufeld & Kaufman 1993

 cool_CO_rovib = Qrot+QvibH+QvibH2

end function cool_CO_rovib

!-----------------------------------------------------------------------
!+
!  RADIATIVE: H20 ro-vibrational cooling (Hollenbach & McKee 1989, Neufeld & Kaufman 1993)
!+
!-----------------------------------------------------------------------
real function cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O)

 use physcon, only: mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nH2O

 real              :: Qrot, QvibH2, QvibH
 real              :: alpha, lambdaH2O
 real              :: v_crit, nfH2O

! Binding energy of singular O-H bond = 5.151 eV = 8.25e-12 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy H2O

 v_crit = sqrt( 2.*8.25d-12/(mu*mass_proton_cgs) )  ! kinetic energy
 nfH2O  = MaxBol_cumul(T_gas, mu,  v_crit) * nH2O

 alpha     = 1.35 - 0.3*log10(T_gas/1000.)          ! Neufeld & Kaufmann 1993
 lambdaH2O = 1.32d-23*(T_gas/1000.)**alpha
 Qrot      = (nH2+1.39*nH)*nfH2O*LambdaH2O

 QvibH2    = 1.03d-26*nH2*nfH2O*T_gas*exp(-2352./T_gas)*exp(-47.5/(T_gas**(1./3.)))            !Hollenbach & McKee 1989 eq 2.14b
 QvibH     = 7.40d-27*nH *nfH2O*lambdaH2O*T_gas*exp(-2352./T_gas)*exp(-34.5/(T_gas**(1./3.)))  !Hollenbach & McKee 1989 eq 2.14a

 cool_H2O_rovib = Qrot+QvibH+QvibH2

end function cool_H2O_rovib

!-----------------------------------------------------------------------
!+
!  RADIATIVE: OH rotational cooling (Hollenbach & McKee 1979, McKee et al. 1982)
!+
!-----------------------------------------------------------------------
real function cool_OH_rot(T_gas, rho_gas, mu, nOH)

 use physcon, only: kboltz, mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nOH

 real              :: n_gas
 real              :: sigma, n_crit
 real              :: v_crit, nfOH

! Binding energy of singular O-H bond = 5.151 eV = 8.25e-12 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy OH

 v_crit = sqrt( 2.*8.25d-12/(mu*mass_proton_cgs) )  ! kinetic energy
 nfOH   = MaxBol_cumul(T_gas, mu,  v_crit) * nOH

 n_gas     = rho_gas/(mu*mass_proton_cgs)
 sigma     = 2.0d-16
 n_crit    = 1.33d7*sqrt(T_gas)

 cool_OH_rot = n_gas*nfOH*(kboltz*T_gas*sigma*v_th(T_gas, mu)) / (1 + n_gas/n_crit + 1.5*sqrt(n_gas/n_crit))  !McKee et al. 1982 eq. 5.2

end function cool_OH_rot

!-----------------------------------------------------------------------
!+
!  UTILITY: Total cooling function
!+
!-----------------------------------------------------------------------
 real function calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL

  calc_Q =  cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust) &
!     + cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust) &
!     + cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust) &
    + cool_coulomb(T_gas, rho_gas, mu, nH, nHe) &
    + cool_HI(T_gas, rho_gas, mu, nH, nHe) &
    + cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe) &
    + cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe) &
    + cool_H2_rovib(T_gas, nH, nH2) &
    + cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2) &
    + cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO) &
    + cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O) &
    + cool_OH_rot(T_gas, rho_gas, mu, nOH) &
    - heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust) &
    - heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust) &
!     - heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL) &
    - heat_CosmicRays(nH, nH2)
!     - heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)

#ifdef DEBUG
  call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
#endif

end function calc_Q

!-----------------------------------------------------------------------
!+
!  UTILITY: Calculate dlnQ/dlnT
!+
!-----------------------------------------------------------------------

real function calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 use timestep,          only:bignumber

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real, intent(in)  :: JL

 real, parameter    :: tolQ    = 1.d-4
 real               :: Qtot, dlnQ_dlnT, dT, Q1, Q2, dQdT
 integer, parameter :: itermax = 20
 integer            :: iter

 Qtot      = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                    T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 dlnQ_dlnT = 0.

! centered finite order approximation for the numerical derivative
 dT = T_gas/100.
 Q1 = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
              T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
 Q2 = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
              T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 iter    = 0
 do while ( abs(Q1/Qtot-1.) > tolQ .and. abs(Q2/Qtot-1.) > tolQ .and. iter < itermax )
    dT   = dT/2.
    Q1   = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                  T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    Q2   = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                  T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    iter = iter + 1
 enddo

 dQdT      = (Q1-Q2)/(2.*dT)
 dlnQ_dlnT = (T_gas/Qtot)*dQdT

! gradient can become large at discontinuous physical temperature boundaries (see e.g. atomic and chemical cooling)
 if (dlnQ_dlnT > bignumber) then
    dlnQ_dlnT = 0.
 endif

 calc_dlnQdlnT = dlnQ_dlnT

end function calc_dlnQdlnT

end module cooling
