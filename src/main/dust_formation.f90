!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: dust_formation
!
!  DESCRIPTION: Dust formation routine using the method of moments
!
!  REFERENCES: Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
!  OWNER: Not Committed Yet
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon
!+
!--------------------------------------------------------------------------

module dust_formation
 implicit none
 character(len=*), parameter :: label = 'dust_formation'

 public :: set_abundances,set_cooling,evolve_chem,calc_kappa_dust,calc_dust_cooling_rate
 logical, public :: calc_Teq

 private
! Indices for elements and molecules:
 integer, parameter :: nElements = 10, nMolecules = 25
 integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
 integer, parameter :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, iNH3=10, iCN=11, &
       iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, iHS=19, iH2S=20, iSiS=21, &
       iSiH=22, iTiO=23, iTiO2=24, iC2 = 25, iTiS=26
 real, parameter :: coefs(5,nMolecules) = reshape([&
       4.25321d+05, -1.07123d+05, 2.69980d+01, 5.48280d-04, -3.81498d-08, & !H2-
       4.15670d+05, -1.05260d+05, 2.54985d+01, 4.78020d-04, -2.82416d-08, & !OH-
       8.66184d+05, -2.27851d+05, 5.61473d+01, 7.62548d-04, -4.95254d-08, & !H2O
       3.30340d+05, -2.59792d+05, 3.23662d+01, 3.33339d-04, -1.69521d-08, & !CO
       5.34072d+05, -3.88536d+05, 6.95497d+01, 0.00000d+00,  0.00000d+00, & !CO2
       1.51591d+06, -4.09711d+05, 1.19952d+02, 5.71281d-04, -4.77554d-08, & !CH4
       6.31402d+05, -2.85395d+05, 5.97427d+01, 1.13748d-04, -2.04586d-08, & !C2H
       8.11613d+05, -3.99041d+05, 9.15250d+01, 2.53889d-05, -2.03949d-08, & !C2H2
       3.80144d+05, -2.28698d+05, 3.09554d+01, 2.71276d-04, -6.21348d-09, & !N2
       1.23626d+06, -2.89588d+05, 8.52686d+01, 5.66638d-04, -3.66862d-08, & !NH3
       4.51527d+05, -1.83350d+05, 2.97771d+01, 1.58267d-04, -1.72510d-08, & !CN
       6.02063d+05, -3.08724d+05, 6.00139d+01, 1.66592d-04, -1.13252d-08, & !HCN
      -1.10204d+05, -7.41179d+04, 2.57470d+01, 1.37542d-04, -2.83506d-09, & !Si2
       6.06066d+04, -1.71746d+05, 5.75688d+01, 1.18547d-04, -6.23584d-08, & !Si3
       1.82780d+05, -1.92839d+05, 3.03804d+01, 4.16079d-04, -1.91035d-08, & !SiO
       2.34924d+05, -2.60546d+05, 6.29047d+01, 2.75378d-04, -3.64630d-08, & !Si2C
       9.83348d+05, -3.16438d+05, 1.13885d+02, 2.90172d-04, -2.37569d-08, & !SiH4
       2.24963d+05, -1.03636d+05, 2.85814d+01, 1.77872d-04, -8.41752d-09, & !S2
       3.65656d+05, -8.77172d-04, 2.41557d+01, 3.40917d-04, -2.09217d-08, & !HS
       7.41911d+05, -1.81063d+05, 5.35114d+01, 4.94594d-04, -3.39231d-08, & !H2S
       1.44507d+05, -1.49916d+05, 2.87381d+01, 3.96186d-04, -2.18002d-08, & !SiS
       2.45408d+05, -7.17752d+04, 2.29449d+01, 4.52999d-04, -2.69915d-08, & !SiH
      -4.38897d+05, -1.58111d+05, 2.49224d+01, 1.08714d-03, -5.62504d-08, & !TiO
      -3.32351d+05, -3.04694d+05, 5.86984d+01, 1.17096d-03, -5.06729d-08, & !TiO2
       2.26786d+05, -1.43775d+05, 2.92429d+01, 1.69434d-04, -1.79867d-08], shape(coefs)) !C2
 real, parameter :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 0., 55.847, 0.] ! Atomic weight for S and Ti missing
 real, parameter :: patm = 1.013250d6 ! Standard atmospheric pressure
 real, parameter :: Scrit = 2. ! Critical saturation ratio

 real :: mass_per_H, eps(nElements)
 logical :: cool_radiation_H0, cool_relaxation_Bowen, cool_collisions_dust, cool_relaxation_Stefan

contains

!-----------------------------------------------------------------------
!
!  calculate cooling rates
!
!-----------------------------------------------------------------------
subroutine calc_dust_cooling_rate(rho, T, Teq, gamma, mu, Cprime, K2, kappa, Q)
! all quantities in cgs
 real, intent(in) :: rho, T, Teq, gamma, mu, Cprime, K2, kappa
 real, intent(out) :: Q
 real :: Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan

 Q_H0 = 0.
 Q_relax_Bowen = 0.
 Q_col_dust = 0.
 Q_relax_Stefan = 0.
 if (cool_radiation_H0)      Q_H0 = cooling_neutral_hydrogen(T, rho)
 if (cool_relaxation_Bowen)  Q_relax_Bowen = cooling_Bowen_relaxation(T, Teq, rho, gamma, mu, Cprime)
 if (cool_collisions_dust)   Q_col_dust = cooling_collision_dust_grains(T, Teq, rho, K2, mu)
 if (cool_relaxation_Stefan) Q_relax_Stefan = cooling_radiative_relaxation(T, Teq, kappa)
 Q = Q_H0 + Q_relax_Bowen+ Q_col_dust+ Q_relax_Stefan
end subroutine calc_dust_cooling_rate

!-----------------------------------------------------------------------
!
!  evolve the chemistry and moments
!
!-----------------------------------------------------------------------
subroutine evolve_chem(dt, r, v, T, rho, Jstar, K, mu, S)
!all quantities in cgs
 use physcon, only:pi,kboltz,atomic_mass_unit
 real, intent(in) :: dt, r, v, T, rho
 real, intent(inout) :: Jstar, K(0:3), S
 real, intent(out) :: mu

 real :: pC, pC2, pC2H, pC2H2, nH_tot, epsC
 real :: JstarS, taustar, taugr
 real :: Jstar_new, K_new(0:3)
 real :: nH, nH2, v1
 real, parameter :: A0 = 20.7d-16
 real, parameter :: alpha2 = 0.34
 real, parameter :: vfactor = sqrt(kboltz/(8.*pi*atomic_mass_unit*12.01))

 nH_tot = rho/mass_per_H
 epsC = eps(iC) - K(3)/(r**2*v*nH_tot)
 if (T > 450.) then
    call chemical_equilibrium_light(rho, T, epsC, pC, pC2, pC2H, pC2H2, mu)
    S = pC/psat_C(T)
    if (S > Scrit) then
       call nucleation(T, pC, 0., 0., 0., pC2H2, S, JstarS, taustar, taugr)
       JstarS = JstarS * r**2*v
       call evol_K(K, Jstar, JstarS, taustar, taugr, dt, Jstar_new, K_new)
    else
       Jstar_new = Jstar
       K_new(:) = K
    endif
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    nH = 0.
    nH2 = nH_tot/2.
    mu = (1.+4.*eps(iHe))*nH_tot/(nH+nH2+eps(iHe)*nH_tot)
    pC2H2 = .5*(epsC-eps(iOx))*nH_tot * kboltz * T
    pC2H = 0.
    S = 0.
    v1 = sqrt(kboltz*T/(8.*pi*atomic_mass_unit*12.01))
!      v1 = vfactor*sqrt(T)
    taugr = kboltz*T/(A0*v1*sqrt(2.)*alpha2*(pC2H+pC2H2))
    call evol_K(K, 0., 0., 1., taugr, dt, Jstar_new, K_new)
 endif
 Jstar = Jstar_new
 K(:) = K_new(:)
end subroutine evolve_chem

!-----------------------------------------------------------------------
!
!  calculate dust opacity
!
!-----------------------------------------------------------------------
subroutine calc_kappa_dust(K3, Mdot_cgs, kappa_planck, kappa_rosseland)
!all quantities in cgs
 use physcon, only:pi
 real, intent(in) :: K3, Mdot_cgs
 real, intent(out) :: kappa_planck, kappa_rosseland

 real :: fC
 real, parameter :: rho_Cdust = 2.62
 real, parameter :: Qplanck_abs = 1.6846124267740528e+04
 real, parameter :: Qross_ext = 9473.2722815583311

 !carbon fraction
 fC = (4.*pi*K3*mass_per_H)/(Mdot_cgs*eps(iC))
 fC = (3./4.)*eps(iC)*12./(1.4*rho_Cdust) * max(fC, 1.d-15)
 kappa_planck    = Qplanck_abs * fC
 kappa_rosseland = Qross_ext * fC
 !should add gas contribution
 !!!kappa_rosseland = kappa_rosseland + 2.d-4
end subroutine calc_kappa_dust

!----------------------------
!
!  Compute nucleation rate
!
!----------------------------
subroutine nucleation(T, pC, pC2, pC3, pC2H, pC2H2, S, Jstar, taustar, taugr)
! all quantities are in cgs
 use physcon, only:kboltz,pi,atomic_mass_unit
 real, intent(in) :: T, pC, pC2, pC3, pC2H, pC2H2, S
 real, intent(out) :: Jstar, taustar, taugr
 real, parameter :: A0 = 20.7d-16
 real, parameter :: sigma = 1400.
 real, parameter :: theta_inf = A0*sigma/kboltz
 real, parameter :: alpha1 = 0.37
 real, parameter :: alpha2 = 0.34
 real, parameter :: alpha3 = 0.08
 real, parameter :: Nl = 5.
 real, parameter :: mproton = 1.6726485d-24

 real :: ln_S_g, Nstar_inf_13, Nstar, Nstar_m1_13, Nstar_m1_23, theta_Nstar
 real :: dtheta_dNstar, d2lnc_dN2star, Z, A_Nstar, v1, beta, c_star, expon

 ln_S_g = log(S)
 Nstar_inf_13 = 2.*theta_inf/(3.*T*ln_S_g)
 Nstar = 1. + (Nstar_inf_13**3)/8. * (1. + sqrt(1. + 2.*Nl**(1./3.)/Nstar_inf_13) - 2.*Nl**(1./3.)/Nstar_inf_13)**3
 Nstar_m1_13 = (Nstar - 1.)**(1./3.)
 Nstar_m1_23 = Nstar_m1_13**2
 theta_Nstar = theta_inf/(1. + Nl**(1./3.)/Nstar_m1_13)
 dtheta_dNstar = Nl**(1./3.)/3. * theta_Nstar**2/(theta_inf * Nstar_m1_23**2)
 d2lnc_dN2star = -2./T * Nstar_m1_23 * (dtheta_dNstar**2/theta_Nstar - theta_Nstar/(9.*(Nstar-1.)**2))
 Z = sqrt(d2lnc_dN2star/(2.*pi))
 A_Nstar = A0 * Nstar**(2./3.)
 v1 = sqrt(kboltz*T/(8.*pi*12.01*mproton))
 beta = v1/(kboltz*T) * (pC*alpha1 + 4.*alpha2/sqrt(2.)*(pC2 + pC2H + pC2H2) + 9.*alpha3/sqrt(3.)*pC3)
 expon = (Nstar-1.)*ln_S_g - theta_Nstar*Nstar_m1_23/T
 if (expon < -100.) then
    c_star = 1.d-99
 else
    c_star = pC/(kboltz*T) * exp(expon)
 endif
 Jstar = beta * A_Nstar * Z * c_star
 taustar = 1./(d2lnc_dN2star*beta*A_Nstar)
 v1 = sqrt(kboltz*T/(8.*pi*12.01*atomic_mass_unit))
 taugr = kboltz*T/(A0*v1*(alpha1*pC*(1.-1./S) + 2.*alpha2/sqrt(2.)*(pC2+pC2H+pC2H2)*(1.-1./S**2)))
end subroutine nucleation

!------------------------------------
!
!  Compute evolution of the moments
!
!------------------------------------
subroutine evol_K(K, Jstar, JstarS, taustar, taugr, dt, Jstar_new, K_new)
! all quantities are in cgs
 real, intent(in) :: K(0:3), Jstar, JstarS, taustar, taugr, dt
 real, intent(out) :: Jstar_new, K_new(0:3)

 real :: d, i0, i1, i2, i3, i4, i5, dK3

 d = dt/taustar
 if (d > 500.) then
    i0 = 0.
 else
    i0 = exp(-d)
 endif
 i1 = 1. - i0
 i2 = d - i1
 i3 = d**2/2. - i2
 i4 = d**3/6. - i3
 i5 = d**4/24. - i4
 Jstar_new = Jstar*i0 + JstarS*i1
 K_new(0) = K(0) + taustar*(Jstar*i1 + JstarS*i2)
 K_new(1) = K(1) + K(0)*dt/(3.*taugr) + taustar**2/(3.*taugr)*(Jstar*i2 + JstarS*i3)
 K_new(2) = K(2) + 2.*K(1)*dt/(3.*taugr) + (dt/(3.*taugr))**2*K(0) + 2.*taustar**3/(3.*taugr)**2 * (Jstar*i3 + JstarS*i4)
 dK3 = 3.*dt/(3.*taugr)*K(2) + 3.*(dt/(3.*taugr))**2*K(1) &
         + (dt/(3.*taugr))**3*K(0) + (6.*taustar**4)/(3.*taugr)**3*(Jstar*i4+JstarS*i5)
 K_new(3) = K(3) + dK3
end subroutine evol_K

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance in the gas phase
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light(rho, T, epsC, pC, pC2, pC2H, pC2H2, mu)
! all quantities are in cgs
 use physcon,only:kboltz
 real, intent(in) :: rho, T, epsC
 real, intent(out) :: pC, pC2, pC2H, pC2H2, mu
 real :: pH_tot, Kd(nMolecules+1), err, a, b, c, d
 real :: pH, pCO, pO, pSi, pS, pTi, pN
 real :: pC_old, pO_old, pSi_old, pS_old, pTi_old
 integer :: i, nit

! Total partial pressure of H (in atm)
 pH_tot = rho*T*kboltz/(patm*mass_per_H)
! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)

! H
 pH = solve_q(2.*Kd(iH2), &
         1., &
         -pH_tot)
 pCO = epsC*pH_tot

 pC = 0.
 pC_old = 0.
 pO = 0.
 pO_old = 0.
 pSi = 0.
 pSi_old = 0.
 pS = 0.
 pS_old = 0.
 pTi = 0.
 pTi_old = 0.
 err = 1.
 nit = 0

 do while (err > 1.e-5)
! N
    pN = solve_q(2.*Kd(iN2), &
            1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+Kd(iHCN)*pH), &
            -eps(iN)*pH_tot)

! C
    pC = solve_q(2.*(Kd(iC2H)*pH + Kd(iC2H2)*pH**2 + Kd(iC2)), &
            1.+.02*pO*Kd(iCO)+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C), &
            .98*pCO-epsC*pH_tot)

! O
    pO = eps(iOx)*pH_tot/(1.+pC*Kd(iCO)+pH*Kd(iOH)+pH**2*Kd(iH2O)+2.*pO*pC*Kd(iCO2)+pSi*Kd(iSiO))
    pCO = Kd(iCO)*pC*pO

! Si
    a = 3.*Kd(iSi3)
    b = 2.*Kd(iSi2)+2.*Kd(iSi2C)*pC
    c = 1.+Kd(iSiH)*pH+Kd(iSiH4)*pH**4+Kd(iSiO)*pO+Kd(iSiS)*pS
    d = -eps(iSi)*pH_tot
    pSi = -d/(a*pSi**2+b*pSi+c)

! S
    pS = solve_q(2.*Kd(iS2), &
            1.+Kd(iSiS)*pSi+Kd(iHS)*pH+Kd(iH2S)*pH**2, &
            -eps(iS)*pH_tot)

! Ti
    pTi = eps(iTi)*pH_tot/(1.+Kd(iTiO)*pO+Kd(iTiO2)*pO**2+Kd(iTiS)*pS)

    err = 0.
    if (pC > 1.e-50)  err = err + abs((pC-pC_old)/pC)
    if (pO > 1.e-50)  err = err + abs((pO-pO_old)/pO)
    if (pSi > 1.e-50) err = err + abs((pSi-pSi_old)/pSi)
    if (pS > 1.e-50)  err = err + abs((pS-pS_old)/pS)
    if (pTi > 1.e-50) err = err + abs((pTi-pTi_old)/pTi)

    nit = nit + 1
    if (nit == 2000) exit

    pC_old = pC
    pO_old = pO
    pSi_old = pSi
    pS_old = pS
    pTi_old = pTi
 enddo

 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2
 mu = (1.+4.*eps(iHe))*pH_tot/(pH+Kd(iH2)*pH**2+eps(iHe)*pH_tot)

! Convert partial pressures from atm to cgs
 pC    = pC*patm
 pC2   = pC2*patm
 pC2H  = pC2H*patm
 pC2H2 = pC2H2*patm
end subroutine chemical_equilibrium_light

subroutine set_abundances(CO_ratio)
! all quantities in cgs
 use physcon, only:atomic_mass_unit
 real, intent(in) :: CO_ratio
 eps(iH)  = 1.d0
 eps(iHe) = 1.04d-1
 eps(iOx) = 6.d-4
 eps(iN)  = 2.52d-4
 eps(iNe) = 1.17d-4
 eps(iSi) = 3.58d-5
 eps(iS)  = 1.85d-5
 eps(iFe) = 3.24d-5
 eps(iTi) = 8.6d-8
 eps(iC)  = eps(iOx) * CO_ratio
 mass_per_H = atomic_mass_unit*dot_product(Aw,eps)
end subroutine set_abundances

!-----------------------------
!
!  solve 2nd order polynomical
!
!-----------------------------
pure real function solve_q(a, b, c)
 real, intent(in) :: a, b, c
 real :: delta
 if (-4.*a*c/b**2 > 1.d-8) then
    delta = max(b**2-4.*a*c, 0.)
    solve_q = (-b+sqrt(delta))/(2.*a)
 else
    solve_q = -c/b
 endif
end function solve_q

!-------------------------------------------------------------------------
!
!  Compute saturation pressure of carbon clusters C_1, ..., C_5 over graphite
!
!-------------------------------------------------------------------------
pure real function psat_C(T)
! all quantities are in cgs
 real, intent(in) :: T

 real, parameter :: f = 13.8287
 real :: pC1,pC2,pC3,pC4,pC5,T2,T3
 T2 = T*T
 T3 = T*T2
 pC1 = exp(-8.61240d4/T + 1.85106d1 + 5.87980d-4*T - 2.51549d-7*T2 + 3.24892d-11*T3 + f)
 pC2 = exp(-1.01124d5/T + 2.35611d1 + 3.37807d-4*T - 2.94959d-7*T2 + 4.41801D-11*T3 + f)
 pC3 = exp(-9.90261d4/T + 2.81161d1 - 1.55820d-3*T + 1.60002d-7*T2 - 4.47171D-12*T3 + f)
 pC4 = exp(-1.17037d5/T + 2.55579d1 - 5.63869d-6*T - 2.13596d-7*T2 + 3.39660D-11*T3 + f)
 pC5 = exp(-1.18080d5/T + 2.65798d1 + 1.20285d-4*T - 2.68583d-7*T2 + 4.12365D-11*T3 + f)
 psat_C = pC1
end function psat_C

!------------------------------------
!
!  Compute dissociation coefficients
!
!------------------------------------
pure real function calc_Kd(coefs, T)
! all quantities are in cgs
 real, intent(in) :: coefs(5), T
 real, parameter :: R = 1.987165
 real :: G, d
 G = coefs(1)/T + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T)*T)*T
 d = -G/(R*T)
 calc_Kd = exp(d)
end function calc_Kd

pure real function calc_Kd_TiS(T)
! all quantities are in cgs
 real, intent(in) :: T
 real, parameter :: a = 1.3316d1, b = -6.2216, c = 4.5829d-1, d = -6.4903d-2, e = 3.2788d-3
 real :: theta, logKd
 theta = 5040./T
 logKd = a+(b+(c+(d+e*theta)*theta)*theta)*theta
 calc_Kd_TiS = 10.**(-logKd)*patm
end function calc_Kd_TiS

subroutine set_cooling(cool_radiation_H0_in, cool_relaxation_Bowen_in, cool_collisions_dust_in, cool_relaxation_Stefan_in)
 logical, intent(in) :: cool_radiation_H0_in, cool_relaxation_Bowen_in, cool_collisions_dust_in, cool_relaxation_Stefan_in
 cool_radiation_H0 = cool_radiation_H0_in
 cool_relaxation_Bowen = cool_relaxation_Bowen_in
 cool_collisions_dust = cool_collisions_dust_in
 cool_relaxation_Stefan = cool_relaxation_Stefan_in
 calc_Teq = cool_relaxation_Bowen .or. cool_relaxation_Stefan .or. cool_collisions_dust
end subroutine set_cooling


!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling term
!+
!-----------------------------------------------------------------------
real function cooling_Bowen_relaxation(T, Teq, rho, wind_gamma, mu, Cprime)
 use physcon, only: atomic_mass_unit, kboltz
 real, parameter :: Rgas = kboltz/atomic_mass_unit
 real, intent(in) :: T, Teq, rho, wind_gamma, mu, Cprime

 cooling_Bowen_relaxation = Rgas/((wind_gamma-1.)*mu)*(Teq-T)*rho/Cprime
end function cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
real function cooling_collision_dust_grains(T, Teq, rho, K2, mu)
 use physcon, only: kboltz, mass_proton_cgs, pi
 real, intent(in) :: T, Teq, rho, K2, mu
 real, parameter :: f = 0.15, a0 = 1.28e-8
 real :: A

 A = 2. * f * kboltz * a0**2/(mass_proton_cgs**2*mu) &
         * (1.05/1.54) * sqrt(2.*pi*kboltz/mass_proton_cgs) &
         * 2.*K2 * rho
 cooling_collision_dust_grains = A * sqrt(T) * (Teq-T)
 if (cooling_collision_dust_grains  >  1000000.) then
    print *, f, kboltz, a0, mass_proton_cgs, mu
    print *, mu, K2, rho, T, Teq
    print *, A, cooling_collision_dust_grains
    stop 'cooling'
 endif
end function cooling_collision_dust_grains

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
real function cooling_radiative_relaxation(T, Teq, kappa)
 use physcon, only: steboltz
 real, intent(in) :: T, Teq, kappa

 cooling_radiative_relaxation = 4.*steboltz*(Teq**4-T**4)*kappa
end function cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to neutral H (Spitzer)
!+
!-----------------------------------------------------------------------
real function cooling_neutral_hydrogen(T, rho)
 use physcon, only: mass_proton_cgs
 real, intent(in) :: T, rho
 real, parameter :: f = 0.2
 real :: eps_e

 if (T > 3000.) then
    eps_e = calc_eps_e(T)
    !cooling_neutral_hydrogen = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(1.4*mass_proton_cgs)**2
    cooling_neutral_hydrogen = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(mass_per_H)**2
 else
    cooling_neutral_hydrogen = 0.
 endif
end function cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance (Palla et al 1983)
!+
!-----------------------------------------------------------------------
real function calc_eps_e(T)
 real, intent(in) :: T
 real :: k1, k2, k3, k8, k9, p, q

 if (T > 3000.) then
    k1 = 1.88e-10 / T**6.44e-1
    k2 = 1.83e-18 * T
    k3 = 1.35e-9
    k8 = 5.80e-11 * sqrt(T) * exp(-1.58e5/T)
    k9 = 1.7e-4 * k8
    p = .5*k8/k9
    q = k1*(k2+k3)/(k3*k9)
    calc_eps_e = (p + sqrt(q+p**2))/q
 else
    calc_eps_e = 0.
 endif
end function calc_eps_e

end module dust_formation
