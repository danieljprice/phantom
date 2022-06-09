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
 public :: write_options_cooling, read_options_cooling
 public :: find_in_table, implicit_cooling, exact_cooling
 logical, public :: cooling_in_step = .true.
 real,    public :: bowen_Cprime    = 3.000d-5
 real,    public :: GammaKI_cgs     = 2.d-26        ! [erg/s] heating rate for Koyama & Inutuska cooling

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
 integer :: excitation_HI = 0, relax_Bowen = 0, dust_collision = 0, relax_Stefan = 0
 integer :: icool_method  = 0
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
       if ( (excitation_HI+relax_Bowen+dust_collision+relax_Stefan) == 0 .and. &
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
subroutine calc_cooling_rate(r, Q, dlnQ_dlnT, rho, T, Teq, mu, gamma, K2, kappa)

 use units,             only:unit_ergg,unit_density
 !use cooling_molecular, only:do_molecular_cooling,calc_cool_molecular

 real, intent(in)           :: rho, T, Teq     !rho in code units
 real, intent(in)           :: r, mu, gamma
 real, intent(in), optional :: K2, kappa       !cgs
 real, intent(out)          :: Q, dlnQ_dlnT    !code units
 
 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, Q_molec, rho_cgs
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan, dlnQ_molec

 rho_cgs           = rho*unit_density
 Q_H0              = 0.
 Q_relax_Bowen     = 0.
 Q_col_dust        = 0.
 Q_relax_Stefan    = 0.
 Q_molec           = 0.

 dlnQ_H0           = 0.
 dlnQ_relax_Bowen  = 0.
 dlnQ_col_dust     = 0.
 dlnQ_relax_Stefan = 0.
 dlnQ_molec        = 0.

 if (excitation_HI  == 1) &
                    call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (relax_Bowen    == 1 ) &
                    call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, gamma, Q_relax_Bowen, dlnQ_relax_Bowen)
 if (dust_collision == 1 .and. present(K2)) &
                    call cooling_dust_collision(T, Teq, rho_cgs, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (relax_Stefan   == 1 .and. present(kappa)) &
                    call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
 !if (do_molecular_cooling) &
 !                   call calc_cool_molecular(T, r, rho_cgs, Q_molec, dlnQ_molec)

 Q_cgs = Q_H0 + Q_relax_Bowen + Q_col_dust + Q_relax_Stefan + Q_molec
 if (Q_cgs == 0.) then
    dlnQ_dlnT = 0.
 else
    dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen + Q_col_dust*dlnQ_col_dust&
   + Q_relax_Stefan*dlnQ_relax_Stefan + Q_molec*dlnQ_molec)/Q_cgs
 endif
 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q         = Q_cgs/unit_ergg
end subroutine calc_cooling_rate


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
subroutine cooling_neutral_hydrogen(T, rho, Q, dlnQ_dlnT)

 use physcon, only: mass_proton_cgs, pi
 
 real, intent(in)  :: T, rho
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter   :: f = 1.0d0
 real              :: eps_e

 if (T > 3000.) then
    eps_e = calc_eps_e(T)
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
subroutine explicit_cooling (xi,yi,zi, ui, dudt, rho, dt, mu, gamma, Trad, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg
 
 real, intent(in)           :: xi, yi, zi, ui, rho, dt, Trad, mu, gamma !code units
 real, intent(in), optional :: K2, kappa
 real, intent(out)          :: dudt                                     !code units
 
 real                       :: u,Q,dlnQ_dlnT,T,T_on_u, r

 r      = sqrt(xi*xi + yi*yi + zi*zi)
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui
 call calc_cooling_rate(r, Q, dlnQ_dlnT, rho, T, Trad, mu, gamma, K2, kappa)
 if (-Q*dt  > ui) then   ! assume thermal equilibrium
    u    = Trad/T_on_u
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
subroutine implicit_cooling (xi,yi,zi, ui, dudt, rho, dt, mu, gamma, Trad, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg
 
 real, intent(in)           :: xi,yi,zi, ui, rho, dt, mu, gamma
 real, intent(in), optional :: Trad, K2, kappa
 real, intent(out)          :: dudt

 real, parameter    :: tol      = 1.d-4    ! to be adjusted
 integer, parameter :: iter_max = 200
 real               :: u,Q,dlnQ_dlnT,T,T_on_u,delta_u,term1,term2
 real               :: r                   ! in au
 integer            :: iter

 u       = ui
 T_on_u  = (gamma-1.)*mu*unit_ergg/Rg
 delta_u = 1.d-3
 iter    = 0
 r       = sqrt(xi**2+yi**2+zi**2)
 !term1 = 1.-(gamma-1.)*dt*divcurlv !pdv=(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
 term1 = 1.
 do while (abs(delta_u) > tol .and. iter < iter_max)
    T       = u*T_on_u
    call calc_cooling_rate(r,Q,dlnQ_dlnT, rho, T, Trad, mu, gamma, K2, kappa)
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
subroutine energ_cooling(xi,yi,zi,ui,dudt,rho,dt,Trad,mu_in,gamma_in,K2,kappa,Tgas)
 use io,  only: fatal
 use eos, only: gmw,gamma
 real, intent(in)           :: xi,yi,zi,ui,rho,dt                  ! in code units
 real, intent(in), optional :: Tgas,Trad,mu_in,gamma_in,K2,kappa   ! in cgs
 real, intent(out)          :: dudt                                ! in code units
 real                       :: mu,polyIndex

 dudt       = 0.
 mu         = gmw
 polyIndex  = gamma
 if (present(gamma_in)) polyIndex = gamma_in
 if (present(mu_in))    mu        = mu_in

 select case (icooling)
 case(1,2)
    if (icool_method == 2) then
       call exact_cooling   (xi,yi,zi, ui, dudt, rho, dt ,mu, polyIndex, Trad, K2, kappa)
    elseif (icool_method == 0) then
       call implicit_cooling(xi,yi,zi, ui, dudt, rho, dt, mu, polyIndex, Trad, K2, kappa)
    else
       call explicit_cooling(xi,yi,zi, ui, dudt, rho, dt, mu, polyIndex, Trad, K2, kappa)
    endif
 case (3)
    call cooling_Gammie(xi,yi,zi,ui,dudt)
 case (4)
    call cooling_Townsend_table(ui,rho,dt,dudt,mu,polyIndex)
 case (5)
    if (present(Tgas)) then
       call cooling_KoyamaInutuska_explicit(rho,Tgas,dudt)
    else
       call fatal('energ_cooling','Koyama & Inutuska cooling requires gas temperature')
    endif
 case (6)
    call cooling_KoyamaInutuska_implicit(ui,rho,dt,dudt)
 case default
    call implicit_cooling(xi,yi,zi,ui, dudt, rho, dt, mu, polyIndex, Trad, K2, kappa)
 end select

end subroutine energ_cooling

!-----------------------------------------------------------------------
!
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!
!-----------------------------------------------------------------------
subroutine exact_cooling    (xi,yi,zi, ui, dudt, rho, dt, mu, gamma, Trad, K2, kappa)

 use physcon, only:Rg
 use units,   only:unit_ergg
 
 real, intent(in)           :: xi,yi,zi, ui, rho, dt, Trad, mu, gamma
 real, intent(in), optional :: K2, kappa
 real, intent(out)          :: dudt

 real, parameter :: tol = 1.d-12
 real            :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,T_on_u
 real            :: r
 integer         :: k

 r      = sqrt(xi*xi + yi*yi + zi*zi)
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T      = T_on_u*ui

 if (T < T_floor) then
    Temp = T_floor
 elseif (T > Tref) then
    call calc_cooling_rate(r,Q, dlnQ_dlnT, rho, T, Trad, mu, gamma, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    call calc_cooling_rate(r,Qref,dlnQref_dlnT, rho, Tref, Trad, mu, gamma, K2, kappa)
    Y         = 0.
    k         = nTg
    Q         = Qref          ! default value if Tgrid < T for all k
    dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
    do while (Tgrid(k) > T)
       k = k-1
       call calc_cooling_rate(r,Q, dlnQ_dlnT, rho, Tgrid(k), Trad, mu, gamma, K2, kappa)
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
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=cooling in substep, 2=cooling in force, &
                               3=Gammie, 4=Townsend table, 5,6=KI02)',iunit)
    select case(icooling)
    case(1,2)
       !call write_options_molecularcooling(iunit)
       call write_inopt(icool_method,'icool_method',&
            'integration method (0=implicit, 1=explicit, 2=exact solution)',iunit)
       call write_inopt(excitation_HI,'excitation_HI','cooling via electron excitation of HI (1=on/0=off)',iunit)
       call write_inopt(relax_bowen,'relax_bowen','Bowen (diffusive) relaxation (1=on/0=off)',iunit)
       call write_inopt(relax_stefan,'relax_stefan','radiative relaxation (1=on/0=off)',iunit)
       call write_inopt(dust_collision,'dust_collision','dust collision (1=on/0=off)',iunit)
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
 integer, save :: ngot = 0
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
 !if (icooling > 0 .and. .not. imatch) call read_options_molecular_cooling(name,valstring,imatch,igotallmol,ierr)
 if (icooling == 0 .and. ngot >= 2) igotall = .true.
 if (icooling == 1 .and. ngot >= 8) igotall = .true.
 if (icooling == 2 .and. ngot >= 7) igotall = .true.
 if (icooling == 3 .and. ngot >= 1) igotall = .true.
 if (icooling == 4 .and. ngot >= 3) igotall = .true.
 if (h2chemistry .and. igotallh2 .and. ngot >= 1) igotall = .true.

end subroutine read_options_cooling








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
!  ADDITIONAL PHYSICS: compute electron equilibrium abundance using LTE network
!                      Warning: Unable to re-derive p and q factors used below,
!                      use with caution!
!+
!-----------------------------------------------------------------------
real function n_e_old(T_gas, rho_gas, mu, nH, nHe)

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe
 
 real             :: n_gas, Te, k1, k2, k3, k8, k9, p, q
 
 ! An up-to-date summary of the relevant electron reations can be found in Vorobyov et al. (2020)
 ! Additionally, a significant overview of primordial reactions, including a large amount of
 ! electron reactions can be found in Glover & Abel (2008)

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 
 k1 = 2.753d-14*(315614/T_gas)**1.5*(1.0+(115188/T_gas))**(-2.242)  ! Ferland et al (1992)
 if (T_gas < 6000.) then                                            ! Wishart (1979)
    k2 = 10**(-17.845+0.762*log10(T_gas)+0.1523*log10(T_gas)**2 &
              -0.03274*log10(T_gas)**3)
 else
    k2 = 10**(-16.4199+0.1998*log10(T_gas)**2 &
              -5.447d-3*log10(T_gas)**4+415d-5*log10(T_gas)**6)
 endif
 k3 = 1.35d-9*(T_gas**9.8493d-2+3.2852d-1*T_gas**5.561d-1&         ! Kreckel et al. (2010)
               +2.771d-7*T_gas**2.1826) / (1.0 &
               +6.191d-3*T_gas**1.0461+8.9712d-11*T_gas**3.0424 &
               +3.2576d-14*T_gas**3.7741)
 Te = T_gas*8.617333262d-5                                         ! gas temperature in eV
 k8 = exp(-3.271396786d1 +1.35365560d+1*log(Te)    &               ! Janev et al. (1987)
                         -5.73932875d+0*log(Te)**2 &
                         +1.56315498d+0*log(Te)**3 &
                         -2.87705600d-1*log(Te)**4 &
                         +3.48255977d-2*log(Te)**5 &
                         -2.63197617d-3*log(Te)**6 &
                         +1.11954395d-4*log(Te)**7 &
                         -2.03914985d-6*log(Te)**8)
 k9 = 1.27d-17*sqrt(T_gas)*exp(-15800/T_gas)                       ! Black (1981)

! original rates from Palla et al. 1983
!  k1 = 1.88d-10 * T_gas**-0.644
!  k2 = 1.83d-18 * T_gas
!  k3 = 1.35d-9
!  k8 = 5.80d-11 * sqrt(T_gas) * exp(-158000/T_gas)
!  k9 = 1.7d-4 * k8
 
 p = .5*k8/k9
 q = k1*(k2+k3)/(k3*k9)
 
 n_e_old = ((p + sqrt(q+p**2))/q)*n_gas
 
end function n_e_old

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute LTE electron density from SAHA equations
!                      (following D'Angelo & Bodenheimer 2013)
!+
!-----------------------------------------------------------------------
real function n_e(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: kboltz, mass_proton_cgs, mass_electron_cgs, planckhbar, pi

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe
 
 real             :: H_ion, He_ion, n_gas
 real             :: X, KH, xx, Y, KHe, z1
 
 H_ion  = 2.179d-11    ! 13.60 eV in erg
 He_ion = 3.940d-11    ! 24.59 eV in erg
 
 n_gas  = rho_gas/(mu*mass_proton_cgs)
 X      = nH /n_gas
 Y      = nHe/n_gas
 KH     = (mass_proton_cgs/(X*rho_gas))*( mass_electron_cgs*kboltz*T_gas /  &
                                         (2.*pi*planckhbar**2.) )**(3./2.)  &
                                       * exp(-H_ion /(kboltz*T_gas))
                                               
 KHe    = (4.*mass_proton_cgs/rho_gas) *( mass_electron_cgs*kboltz*T_gas /  &
                                         (2.*pi*planckhbar**2.) )**(3./2.) &
                                       * exp(-He_ion/(kboltz*T_gas))
                                               
! solution to quadratic SAHA equations (eqns 16 and 17 in D'Angelo et al 2013)
 xx     = (1./2.) * (-KH    + sqrt(KH**2+4.*KH))
 z1     = (2./Y ) * (-KHe-X + sqrt(KHe**2.+2.*KHe*X+X**2.+KHe*Y))
 
 n_e    = xx * nH + z1 * nHe
 
end function n_e

!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute mean thermal speed of molecules
!+
!-----------------------------------------------------------------------
real function vth(T_gas,mu)

 use physcon, only: kboltz, mass_proton_cgs

 real, intent(in) :: T_gas, mu
 
 vth = sqrt((3.*kboltz*T_gas)/(mu*mass_proton_cgs))
 
end function vth


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
 MaxBol_cumul = erf(v_crit/(sqrt(2.)*a)) - sqrt(2./pi) * (v_crit*exp(-v_crit**2./(2.*a**2.))) / a 
 
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
 
 n_dust = ( rho_gas*d2g ) / ( (4./3.)*pi*a**3*rho_grain )
 
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
real function cool_dust_full_contact(T_gas, rho_gas, mu, T_dust)

 use physcon, only:Rg
 
 real, intent(in)  :: T_gas, T_dust, rho_gas, mu
 
 cool_dust_full_contact = (3*Rg*bowen_Cprime)/(2*mu)*rho_gas*(T_dust-T_gas)

end function cool_dust_full_contact

!-----------------------------------------------------------------------
!+
!  DUST: Discrete contact cooling (Hollenbach & McKee 1979)
!+
!-----------------------------------------------------------------------
real function cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust)

 use physcon, only: kboltz, mass_proton_cgs, pi
 
 real, intent(in)  :: T_gas, T_dust, rho_gas, mu
 
 real, parameter   :: alpha = 0.0
 real              :: n_gas, sigma_dust, d2g, a, rho_grain
 
 sigma_dust                 = 2.*pi*a**2.
 n_gas                      = rho_gas/(mu*mass_proton_cgs)
 cool_dust_discrete_contact = alpha*n_gas*n_dust(rho_gas,d2g,a,rho_grain)*sigma_dust*vth(T_gas,mu)*kboltz*(T_dust-T_gas)
 
end function cool_dust_discrete_contact

!-----------------------------------------------------------------------
!+
!  DUST: Radiative cooling (Woitke 2006) - DO NOT USE, PHYSICALLY INCORRECT
!+
!-----------------------------------------------------------------------
real function cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)

 use physcon, only: steboltz
 
 real, intent(in)  :: T_gas, kappa_gas, T_dust, kappa_dust

 cool_dust_radiation = 4.*steboltz*(kappa_dust*T_dust**4-kappa_gas*T_gas**4)

end function cool_dust_radiation

!-----------------------------------------------------------------------
!+
!  DUST: Viscous heating
!+
!-----------------------------------------------------------------------
real function heat_dust_friction(rho_gas)

 real, intent(in)  :: rho_gas
 
 real              :: sigma_dust, v_drift, d2g, a, rho_grain
 real, parameter   :: alpha = 0.33                            ! see Burke & Hollenbach 1983
 
 ! Warning, alpha depends on the type of dust 
 
 sigma_dust         = 2.*pi*a**2.
 heat_dust_friction = n_dust(rho_gas,d2g,a,rho_grain)*sigma_dust*v_drift*alpha*0.5*rho_gas*v_drift**2

end function heat_dust_friction

!-----------------------------------------------------------------------
!+
!  DUST: photovoltaic heating (Weingartner & Draine 2001)
!+
!-----------------------------------------------------------------------
real function heat_dust_photovoltaic(T_gas, rho_gas, mu, nH, nHe)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 
 real              :: n_gas, x
 real, parameter   :: G=1.68 ! ratio of true background UV field to Habing field
 real, parameter   :: C0=5.45, C1=2.50, C2=0.00945, C3=0.01453, C4=0.147, C5=0.623, C6=0.511 ! see Table 2 in Weingartner & Draine 2001, last line
 
 ! only for soft UV background fields
 
 x  = G*sqrt(T_gas)/n_e(T_gas, rho_gas, mu, nH, nHe)
 
 heat_dust_photovoltaic = 1d-26*G*nH*(C0+C1*T_gas**C4)/(1+C2*x**C5*(1+C3*x**C6))
 
 ! for hard UV background fields, such as in the presence of a WD, use

end function heat_dust_photovoltaic

!-----------------------------------------------------------------------
!+
!  PARTICLE: Coulomb cooling
!+
!-----------------------------------------------------------------------
real function cool_coulomb(T_gas, rho_gas, mu, nH, nHe)

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 
 real              :: x
 real, parameter   :: G=1.68 ! ratio of true background UV field to Habing field
 real, parameter   :: D0=0.4255, D1=2.457, D2=-6.404, D3=1.513, D4=0.05343 ! see Table 3 in Weingartner & Draine 2001, last line
 
 if (T_gas > 1000.) then
    x  = G*sqrt(T_gas)/n_e(T_gas, rho_gas, mu, nH, nHe)
    
    cool_coulomb = 1d-28*n_e(T_gas, rho_gas, mu, nH, nHe)*nH*T_gas**(D0+D1/x)*exp(D2+D3*x-D4*x**2)
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
 
 real             :: Rcr
 
 Rcr = 5.0d-17  !cosmic ray ionisation rate [s^-1]
 
 heat_CosmicRays = Rcr*(5.5d-12*nH2+2.5d-11*nH2)

end function heat_CosmicRays

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to electron excitation of neutral H (Spitzer 1978)
!+
!-----------------------------------------------------------------------
real function cool_HI(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: mass_proton_cgs
 
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 
 real              :: n_gas
 
 ! all hydrogen atomic, so nH = n_gas
 ! Dalgarno & McCray (1972) provide data starting at 3000K
 if (T_gas > 3000.) then
    n_gas = rho_gas/(mu*mass_proton_cgs)
    
    cool_HI = -7.3d-19*n_e(T_gas, rho_gas, mu, nH, nHe)*n_gas*exp(-118400./T_gas)
 else
    cool_HI = 0.0
 endif
 
end function cool_HI

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to collisional ionisation of neutral H
!+
!-----------------------------------------------------------------------
real function cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: mass_proton_cgs
 
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 
 real              :: n_gas

 ! all hydrogen atomic, so nH = n_gas
 if (T_gas > 4000.) then    
    n_gas = rho_gas/(mu*mass_proton_cgs) 
    
    cool_H_ionisation = -1.27d-21*n_e(T_gas, rho_gas, mu, nH, nHe)*n_gas*sqrt(T_gas)*exp(-157809./T_gas)
 else
    cool_H_ionisation = 0.0
 endif
 
end function cool_H_ionisation

!-----------------------------------------------------------------------
!+
!  ATOMIC: Cooling due to collisional ionisation of neutral He
!+
!-----------------------------------------------------------------------
real function cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)
 
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 
 real              :: n_gas

 ! all hydrogen atomic, so nH = n_gas
 if (T_gas > 4000.) then
    n_gas = rho_gas/(mu*mass_proton_cgs) 
    
    cool_He_ionisation = -9.38d-22*n_e(T_gas, rho_gas, mu, nH, nHe)*nHe*sqrt(T_gas)*exp(-285335./T_gas)
 else
    cool_He_ionisation = 0.0
 endif
 
end function cool_He_ionisation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: Cooling due to ro-vibrational excitation of H2
!+
!-----------------------------------------------------------------------
real function cool_H2_rovib(T_gas, nH, nH2)
 
 real, intent(in)  :: T_gas, nH, nH2
 
 real              :: kH_01, kH2_01
 real              :: Lvh, Lvl, Lrh, Lrl
 real              :: x, Qn
 
 if (T_gas < 1635.) then
    kH_01 = 1.4d-13*exp((T_gas/125)-(T_gas/577)**2)
 else
    kH_01 = 1.2d-12*sqrt(T_gas)*exp(-1000/T_gas)
 endif
 kH2_01 = 1.45d-12*sqrt(T_gas)*exp(-28728/(T_gas+1190.))
 Lvh    = 1.1d-18*exp(-6744/T_gas)
 Lvl    = 8.18d-13*exp(-6840/T_gas)*(nH*kH_01+nH2*kH2_01)
 
 x   = log10(T_gas/1.0d4)
 if (T_gas < 1087.) then
    Lrh = 10**(-19.24+0.474*x-1.247*x**2)
 else
    Lrh = 3.9d-19*exp(-6118/T_gas)
 endif
 
 Qn = nH2**0.77+1.2*nH**0.77
 if (T_gas > 4031.) then
    Lrl = 10**(-22.9-0.553*x-1.148*x**2)*Qn
 else
    Lrl = 1.38d-22*exp(-9243/T_gas)*Qn
 endif
 
 cool_H2_rovib = nH2*( Lvh/(1+Lvh/Lvl) + Lrh/(1+Lrh/Lrl) )
 
end function cool_H2_rovib

!-----------------------------------------------------------------------
!+
!  CHEMICAL: H2 dissociation cooling
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
 n1    = 10**(4.0   -0.416*x -0.327*x**2)
 n2    = 10**(4.845 -1.3*x   +1.62*x**2)
 beta  = (1+n_gas*(2*nH2/n_gas*((1/n2)-(1/n1))+1/n1))**(-1)
 kD_H  = 1.2d-9*exp(-52400/T_gas)*(0.0933*exp(-17950/T_gas))**beta
 kD_H2 = 1.3d-9*exp(-53300/T_gas)*(0.0908*exp(-16200/T_gas))**beta
 
 cool_H2_dissociation = 7.18d-12*(nH2**2*kD_H2+nH*nH2*kD_H)

end function cool_H2_dissociation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: H2 recombination heating
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
 n1     = 10**(4.0   -0.416*x -0.327*x**2)
 n2     = 10**(4.845 -1.3*x   +1.62*x**2)
 beta   = (1+n_gas*(2*nH2/n_gas*((1/n2)-(1/n1))+1/n1))**(-1)
 xi     = 7.18d-12*n_gas*nH*(1-beta)
 
 fa     = (1+1.0d4*exp(-600/T_dust))**(-1)
 k_rec  = 3.0d-18*(sqrt(T_gas)*fa)/(1+0.04*sqrt(T_gas+T_dust)+2.0d-3*T_gas+8.0d-6*T_gas**2) 
 
 heat_H2_recombination = -k_rec*xi

end function heat_H2_recombination

!-----------------------------------------------------------------------
!+
!  RADIATIVE: CO ro-vibrational cooling
!+
!-----------------------------------------------------------------------
real function cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)
 
 use physcon, only: kboltz, mass_proton_cgs
 
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nCO
 
 real              :: Qrot, QvibH2, QvibH
 real              :: n_gas, n_crit, sigma, v_th
 real              :: v_crit, nfCO
 
! CO bond dissociation energy = 11.11 eV = 1.78e-11 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy CO

 v_crit = sqrt( 2.*1.78d-11/(mu*mass_proton_cgs) )  ! kinetic energy
 nfCO   = MaxBol_cumul(T_gas, mu,  v_crit) * nCO
 
 n_gas  = rho_gas/(mu*mass_proton_cgs)
 n_crit = 3.3d6*(T_gas/1000.)**0.75
 sigma  = 3.0d-16*(T_gas/1000.)**(-1/4)
 Qrot   = n_gas*nfCO*(kboltz*T_gas*sigma*v_th(T_gas, mu)) / (1 + n_gas/n_crit + 1.5*sqrt(n_gas/n_crit))
 
 QvibH2 = 1.83d-26*nH2*nfCO*exp(-3080./T_gas)*exp(-68./T_gas**(1./3.))
 QvibH  = 1.28d-24*nH *nfCO*exp(-3080./T_gas)*exp(-(2000/T_gas)**3.43)
 
 cool_CO_rovib = Qrot+QvibH+QvibH2

end function cool_CO_rovib

!-----------------------------------------------------------------------
!+
!  RADIATIVE: H20 ro-vibrational cooling
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
 
 alpha     = 1.35 - 0.3*log10(T_gas/1000.)
 lambdaH2O = 1.32d-23*(T_gas/1000.)**alpha
 Qrot      = (nH2+1.39*nH)*nfH2O*LambdaH2O
 
 QvibH2    = 1.03d-26*nH2*nfH2O*T_gas*exp(-2352/T_gas)*exp(-47.5/T_gas**(1/3))
 QvibH     = 7.40d-27*nH *nfH2O*lambdaH2O*T_gas*exp(-2352/T_gas)*exp(-34.5/T_gas**(1/3))
 
 cool_H2O_rovib = Qrot+QvibH+QvibH2

end function cool_H2O_rovib

!-----------------------------------------------------------------------
!+
!  RADIATIVE: OH rotational cooling
!+
!-----------------------------------------------------------------------
real function cool_OH_rot(T_gas, rho_gas, mu, nOH)

 use physcon, only: kboltz, mass_proton_cgs
 
 real, intent(in)  :: T_gas, rho_gas, mu, nOH
 
 real              :: n_gas
 real              :: sigma, n_crit, v_th
 real              :: v_crit, nfOH
 
! Binding energy of singular O-H bond = 5.151 eV = 8.25e-12 erg
! use cumulative distribution of Maxwell-Boltzmann
! to account for collisions that destroy OH

 v_crit = sqrt( 2.*8.25d-12/(mu*mass_proton_cgs) )  ! kinetic energy
 nfOH   = MaxBol_cumul(T_gas, mu,  v_crit) * nOH

 n_gas     = rho_gas/(mu*mass_proton_cgs)
 sigma     = 2.0d-16
 n_crit    = 1.33d7*sqrt(T_gas) 
 
 cool_OH_rot = n_gas*nfOH*(kboltz*T_gas*sigma*v_th(T_gas, mu)) / (1 + n_gas/n_crit + 1.5*sqrt(n_gas/n_crit))

end function cool_OH_rot

!-----------------------------------------------------------------------
!+
!  UTILITY: Total cooling function
!+
!-----------------------------------------------------------------------
real function calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)
 
 real, intent(in)  :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)  :: T_dust, kappa_dust
 
  calc_Q =  cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust) &
!     + cool_dust_full_contact(T_gas, rho_gas, mu, T_dust) &
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
    - heat_dust_friction(rho_gas) &
    - heat_dust_photovoltaic(T_gas, rho_gas, mu, nH, nHe) &
    - heat_CosmicRays(nH, nH2) &
    - heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)

end function calc_Q

!-----------------------------------------------------------------------
!+
!  UTILITY: Calculate dlnQ/dlnT
!+
!-----------------------------------------------------------------------

real function calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)

 use timestep,          only:bignumber

 real, intent(in)           :: T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real, intent(in)           :: T_dust
 real, intent(in), optional :: kappa_dust
 real                       :: Qtot, dlnQ_dlnT, dT, dQ1, dQ2, dQdT, tolQ
 integer                    :: iter, itermax

 Qtot      = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)
 dlnQ_dlnT = 0.0
 
! centered finite order approximation for the numerical derivative
 dT  = T_gas/100.
 dQ1 = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)
 dQ2 = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)
 
 iter    = 0
 itermax = 20
 tolQ    = 1.0d-4
 do while ( dQ1/Qtot > tolQ .and. dQ2/Qtot > tolQ .and. iter < itermax )    
    dT   = dT/2.
    dQ1  = calc_Q(T_gas+dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)
    dQ2  = calc_Q(T_gas-dT, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, T_dust, kappa_dust)
    iter = iter + 1    
 enddo

 dQdT      = (dQ1-dQ2)/(2.*dT)
 dlnQ_dlnT = (T_gas/Qtot)*dQdT 
 
! gradient can become large at discontinuous physical temperature boundaries of Q (see e.g. atomic and chemical cooling)
 if (dlnQ_dlnT > bignumber) then
    dlnQ_dlnT = 0.0
 endif   
 
 calc_dlnQdlnT = dlnQ_dlnT
 
end function calc_dlnQdlnT

end module cooling
