!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_functions
!
! A library of cooling functions that can be handled by cooling_solver
!  Contributed by Lionel Siess and Ward Homan
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none

 real, public :: bowen_Cprime     = 3.000d-5
 real, public :: lambda_shock_cgs = 1.d0
 real, public :: T1_factor = 20., T0_value = 0.
 real, public :: kappa_dust_min = 1e-3  ! dust opacity value below which dust cooling is not calculated

 public :: cool_dust_discrete_contact, cool_coulomb, &
           cool_HI, cool_H_ionisation, cool_He_ionisation, &
           cool_H2_rovib, cool_H2_dissociation, cool_CO_rovib, &
           cool_H2O_rovib, cool_OH_rot, heat_dust_friction, &
           heat_dust_photovoltaic_soft, heat_CosmicRays, &
           heat_H2_recombination, &
           cool_dust_full_contact, cool_dust_radiation, &
           cooling_neutral_hydrogen, &
           heat_dust_photovoltaic_hard, &
           piecewise_law, &
           cooling_Bowen_relaxation, &
           cooling_dust_collision, &
           cooling_radiative_relaxation, &
           testing_cooling_functions

 private
 real, parameter  :: xH = 0.7, xHe = 0.28 !assumed H and He mass fractions

contains
!-----------------------------------------------------------------------
!+
!  Piecewise cooling law for simple shock problem (Creasey et al. 2011)
!+
!-----------------------------------------------------------------------
subroutine piecewise_law(T, T0, rho_cgs, ndens, Q, dlnQ)

 real, intent(in)  :: T, T0, rho_cgs, ndens
 real, intent(out) :: Q, dlnQ
 real :: T1,Tmid !,dlnT,fac

 T1 = T1_factor*T0
 Tmid = 0.5*(T0+T1)
 if (T < T0) then
    Q    = 0.
    dlnQ = 0.
 elseif (T >= T0 .and. T <= Tmid) then
    !dlnT = (T-T0)/(T0/100.)
    Q = -lambda_shock_cgs*ndens**2/rho_cgs*(T-T0)/T0
    !fac = 2./(1.d0 + exp(dlnT))
    dlnQ = 1./(T-T0+epsilon(0.))
 elseif (T >= Tmid .and. T <= T1) then
    Q = -lambda_shock_cgs*ndens**2/rho_cgs*(T1-T)/T0
    dlnQ = -1./(T1-T+epsilon(0.))
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
subroutine cooling_Bowen_relaxation(T, Tdust, rho_cgs, mu, gamma, Q_cgs, dlnQ_dlnT)

 use physcon, only:Rg

 real, intent(in)  :: T, Tdust, rho_cgs, mu, gamma
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 Q_cgs     = Rg/((gamma-1.)*mu)*rho_cgs*(Tdust-T)/bowen_Cprime
 dlnQ_dlnT = -T/(Tdust-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_dust_collision(T, Tdust, rho, K2, mu, Q_cgs, dlnQ_dlnT)

 use physcon, only: kboltz, mass_proton_cgs, pi

 real, intent(in)  :: T, Tdust, rho, K2, mu
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 real, parameter   :: f = 0.15, a0 = 1.28e-8
 real              :: A

 A = 2. * f * kboltz * a0**2/(mass_proton_cgs**2*mu) &
         * (1.05/1.54) * sqrt(2.*pi*kboltz/mass_proton_cgs) * 2.*K2 * rho
 Q_cgs = A * sqrt(T) * (Tdust-T)
 if (Q_cgs  >  1.d6) then
    print *, f, kboltz, a0, mass_proton_cgs, mu
    print *, mu, K2, rho, T, Tdust, A, Q_cgs
    stop 'cooling'
 else
    dlnQ_dlnT = 0.5+T/(Tdust-T+1.d-10)
 endif

end subroutine cooling_dust_collision

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_radiative_relaxation(T, Tdust, kappa, Q_cgs, dlnQ_dlnT)

 use physcon, only: steboltz

 real, intent(in) :: T, Tdust, kappa
 real, intent(out) :: Q_cgs, dlnQ_dlnT

 Q_cgs     = 4.*steboltz*(Tdust**4-T**4)*kappa
 dlnQ_dlnT = -4.*T**4/(Tdust**4-T**4+1.d-10)

end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to electron excitation of neutral H (Spitzer 1978)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen(T, rho_cgs, Q_cgs, dlnQ_dlnT)

 use physcon, only: mass_proton_cgs

 real, intent(in)  :: T, rho_cgs
 real, intent(out) :: Q_cgs,dlnQ_dlnT

 real, parameter   :: f = 1.0d0
 real              :: ne,nH

 if (T > 3000.) then
    nH = rho_cgs/(1.4*mass_proton_cgs)
    ne = calc_eps_e(T)*nH
    !the term 1/(1+sqrt(T)) comes from Cen (1992, ApjS, 78, 341)
    Q_cgs  = -f*7.3d-19*ne*nH*exp(-118400./T)/rho_cgs/(1.+sqrt(T/1.d5))
    dlnQ_dlnT = -118400./T+log(nH*calc_eps_e(1.001*T)/ne)/log(1.001) &
         - 0.5*sqrt(T/1.d5)/(1.+sqrt(T/1.d5))
 else
    Q_cgs = 0.
    dlnQ_dlnT = 0.
 endif

end subroutine cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance per nH atom (Palla et al 1983)
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
!  cooling functions with analytical solutions (for analysis)
!+
!-----------------------------------------------------------------------
subroutine testing_cooling_functions(ifunct, T, Q, dlnQ_dlnT)

 integer, intent(in) :: ifunct
 real, intent(in) :: T
 real, intent(out) :: Q,dlnQ_dlnT

 select case(ifunct)
 case (0)
    !test1 : du/dt = cst --> linear decrease in time
    Q = -1e13
    dlnQ_dlnT = 0.
 case(1)
    !test2 : du/dt = -a*u --> exponential decrease in time
    Q = -1e7*T
    dlnQ_dlnT = 1.
 case(3)
    !test3 : du/dt = -a*u**3 --> powerlaw decrease in time**(-1/2)
    Q = -1e-5*T**3
    dlnQ_dlnT = 3.
 case(-3)
    !test4 : du/dt = -a*u**-3 --> powerlaw decrease in time**(1/4)
    Q = -1e20/T**3
    dlnQ_dlnT = -3.
 case default
    Q = 0.
    dlnQ_dlnT = 0.
 end select

end subroutine testing_cooling_functions





!-----------------------------------------------------------------------
!+
!  ADDITIONAL PHYSICS: compute LTE electron density from SAHA equations
!                      (following D'Angelo & Bodenheimer 2013)
!+
!-----------------------------------------------------------------------
real function n_e(T_gas, rho_gas, mu, nH, nHe)

 use physcon, only: kboltz, mass_proton_cgs, mass_electron_cgs, planckhbar, pi

 real, intent(in) :: T_gas, rho_gas, mu, nH, nHe

 real, parameter  :: H2_diss = 7.178d-12    !  4.48 eV in erg
 real, parameter  :: H_ion   = 2.179d-11    ! 13.60 eV in erg
 real, parameter  :: He_ion  = 3.940d-11    ! 24.59 eV in erg
 real, parameter  :: He2_ion = 8.720d-11    ! 54.42 eV in erg
 real             :: KH, KH2, xx, yy, KHe, KHe2, z1, z2, cst

 cst = mass_proton_cgs/rho_gas*sqrt(mass_electron_cgs*kboltz*T_gas/(2.*pi*planckhbar**2))**3
 if (T_gas > 1.d5) then
    xx = 1.
 else
    KH   = cst/xH * exp(-H_ion /(kboltz*T_gas))
    ! solution to quadratic SAHA equations (Eq. 16 in D'Angelo et al 2013)
    xx   = 0.5 * (-KH    + sqrt(KH**2+4.*KH))
 endif
 if (T_gas > 1.d4) then
    yy = 1.
 else
    KH2  = 0.5*sqrt(0.5*mass_proton_cgs/mass_electron_cgs)**3*cst/xH * exp(-H2_diss/(kboltz*T_gas))
    ! solution to quadratic SAHA equations (Eq. 15 in D'Angelo et al 2013)
    yy   = 0.5 * (-KH    + sqrt(KH2**2+4.*KH2))
 endif
 if (T_gas > 3.d5) then
    z1 = 1.
    z2 = 1.
 else
    KHe    = 4.*cst * exp(-He_ion/(kboltz*T_gas))
    KHe2   =    cst * exp(-He2_ion/(kboltz*T_gas))
! solution to quadratic SAHA equations (Eq. 17 in D'Angelo et al 2013)
    z1     = (2./XHe ) * (-KHe-xH + sqrt((KHe+xH)**2+KHe*xHe))
! solution to quadratic SAHA equations (Eq. 18 in D'Angelo et al 2013)
    z2     = (2./xHe ) * (-KHe2-xH + sqrt((KHe+xH+xHe/4.)**2+KHe2*xHe))
 endif
 n_e       = xx * nH + z1*(1.+z2) * nHe
 !mu  = 4./(2.*xH*(1.+xx+2.*xx*yy)+xHe*(1+z1+z1*z2))

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
! doi : 10.4236/ijaa.2020.103010
!-----------------------------------------------------------------------
real function MaxBol_cumul(T_gas, mu,  v_crit)

 use physcon, only: kboltz, mass_proton_cgs, pi

 real, intent(in) :: T_gas, mu, v_crit

 real             :: a

 a            = sqrt(2.*kboltz*T_gas/(mu*mass_proton_cgs))
 MaxBol_cumul = erf(v_crit/a) - 2./sqrt(pi) * v_crit/a *exp(-(v_crit/a)**2)

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
!  Cooling functions    **** ALL IN cgs  ****
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
 use physcon, only:pi
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

 if (T_gas > 1000.) then !. .and. T_gas < 1.e4) then
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
    !nH      = XH*n_gas
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
    n_gas   = rho_gas/(mu*mass_proton_cgs)
    !nH      = XH*n_gas
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

 use physcon, only:mass_proton_cgs

 real, intent(in)  :: T_gas, rho_gas, mu, nH, nHe
 real              :: n_gas

 ! all hydrogen atomic, so nH = n_gas
 ! (1+sqrt(T_gas/1.d5))**(-1) correction factor added by Cen 1992
 if (T_gas > 4000.) then
    n_gas   = rho_gas/(mu*mass_proton_cgs)
    !nH      = XH*n_gas
    cool_He_ionisation = 9.38d-22*n_e(T_gas, rho_gas, mu, nH, nHe)*nHe*sqrt(T_gas)*(1+sqrt(T_gas/1.d5))**(-1)*exp(-285335./T_gas)
 else
    cool_He_ionisation = 0.0
 endif

end function cool_He_ionisation

!-----------------------------------------------------------------------
!+
!  CHEMICAL: Cooling due to ro-vibrational excitation of H2 (Lepp & Shull 1983)
!            (Smith & Rosen, 2003, MNRAS, 339)
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
 Lvh    = 1.1d-18*exp(-6744./T_gas)
 Lvl    = 8.18d-13*(nH*kH_01+nH2*kH2_01)*exp(-6840./T_gas)

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
!  CHEMICAL: H2 dissociation cooling (Shapiro & Kang 1987, Smith & Rosen 2003)
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
!            for an overview, see Wakelam et al. 2017, Smith & Rosen 2003
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

 fa     = 1./(1.+1.d4*exp(-600./T_dust))    ! eq 3.4
 k_rec  = 3.d-18*(sqrt(T_gas)*fa)/(1.+0.04*sqrt(T_gas+T_dust)+2.d-3*T_gas+8.d-6*T_gas**2) ! eq 3.8

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

 if (T_gas > 3000. .or. T_gas < 250.) then
    cool_CO_rovib = 0.
    return
 endif
 v_crit = sqrt( 2.*1.78d-11/(mu*mass_proton_cgs) )  ! kinetic energy
 nfCO   = MaxBol_cumul(T_gas, mu,  v_crit) * nCO

 n_gas  = rho_gas/(mu*mass_proton_cgs)
 n_crit = 3.3d6*(T_gas/1000.)**0.75      !McKee et al. 1982 eq. 5.3
 sigma  = 3.d-16*(T_gas/1000.)**(-0.25)  !McKee et al. 1982 eq. 5.4
 !v_th = sqrt((8.*kboltz*T_gas)/(pi*mH2_cgs)) !3.1
 Qrot   = 0.5*n_gas*nfCO*kboltz*T_gas*sigma*v_th(T_gas, mu) / (1. + (n_gas/n_crit) + 1.5*sqrt(n_gas/n_crit))
!McKee et al. 1982 eq. 5.2

 QvibH2 = 1.83d-26*nH2*nfCO*T_gas*exp(-3080./T_gas)*exp(-68./(T_gas**(1./3.)))  !Smith & Rosen
 QvibH  = 1.28d-24*nH *nfCO*sqrt(T_gas)*exp(-3080./T_gas)*exp(-(2000./T_gas)**3.43) !Smith & Rosen

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
 !n_crit    = 1.33d7*sqrt(T_gas)
 n_crit    = 1.5d10*sqrt(T_gas/1000.) !table 3 Hollenbach & McKee 1989

 cool_OH_rot = n_gas*nfOH*(kboltz*T_gas*sigma*v_th(T_gas, mu)) / (1 + n_gas/n_crit + 1.5*sqrt(n_gas/n_crit))  !McKee et al. 1982 eq. 5.2

end function cool_OH_rot

end module cooling_functions
