!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dust_formation
!
! Dust formation routine : theory of moments
!
! :References: Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - bowen_Tcond   : *dust condensation temperature (K)*
!   - bowen_delta   : *condensation temperature range (K)*
!   - bowen_kmax    : *maximum dust opacity (cm²/g)*
!   - idust_opacity : *compute dust opacity (0=off, 1=bowen)*
!   - kappa_gas     : *constant gas opacity (cm²/g)*
!   - wind_CO_ratio : *wind initial C/O ratio (> 1)*
!
! :Dependencies: dim, dump_utils, eos, infile_utils, io, part, physcon,
!   units
!

 use part,    only:idJstar,idK0,idK1,idK2,idK3,idmu,idgamma,idsat,idkappa
 use physcon, only:kboltz,pi,atomic_mass_unit,mass_proton_cgs,patm
 use dim,     only:nElements

 implicit none
 integer, public :: idust_opacity = 0

 public :: set_abundances,evolve_dust,evolve_chem,calc_kappa_dust,&
      calc_kappa_bowen,chemical_equilibrium_light,psat_C,calc_nucleation,&
      read_options_dust_formation,write_options_dust_formation,&
      calc_Eddington_factor,calc_muGamma,init_muGamma,init_nucleation,&
      write_headeropts_dust_formation,read_headeropts_dust_formation
!
!--runtime settings for this module
!

 real, public :: kappa_gas   = 2.d-4
 real, public, parameter :: Scrit = 2. ! Critical saturation ratio
 real, public :: mass_per_H, eps(nElements)
 real, public :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]

 private

 character(len=*), parameter :: label = 'dust_formation'
 real :: wind_CO_ratio = 2.
 real :: bowen_kmax  = 2.7991
 real :: bowen_Tcond = 1500.
 real :: bowen_delta = 60.

! Indices for elements and molecules:
 integer, parameter :: nMolecules = 25
 integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
 integer, parameter :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, &
      iNH3=10, iCN=11, iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, &
      iHS=19, iH2S=20, iSiS=21, iSiH=22, iTiO=23, iTiO2=24,iC2 = 25, iTiS=26
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
 real, parameter :: vfactor = sqrt(kboltz/(2.*pi*atomic_mass_unit*12.01))
 !real, parameter :: vfactor = sqrt(kboltz/(8.*pi*atomic_mass_unit*12.01))

contains

subroutine init_nucleation
 use part,  only:npart,nucleation,n_nucleation
 use eos,   only:gamma,gmw
 integer :: i
 real :: JKmuS(n_nucleation)

 call set_abundances

 !initialize nucleation array
 gamma = 5./3.
 JKmuS = 0.
 JKmuS(idmu)    = gmw
 JKmuS(idgamma) = gamma
 do i=1,npart
    nucleation(:,i) = JKmuS(:)
 enddo

end subroutine init_nucleation

subroutine set_abundances
! all quantities in cgs
 eps(iH)  = 1.0
 eps(iHe) = 1.04d-1
 eps(iOx) = 6.d-4
 eps(iN)  = 2.52d-4
 eps(iNe) = 1.17d-4
 eps(iSi) = 3.58d-5
 eps(iS)  = 1.85d-5
 eps(iFe) = 3.24d-5
 eps(iTi) = 8.6d-8
 eps(iC)  = eps(iOx) * wind_CO_ratio
 mass_per_H = atomic_mass_unit*dot_product(Aw,eps)
 !XH  = atomic_mass_unit*eps(iH)/mass_per_H  ! H mass fraction
 !XHe = atomic_mass_unit*eps(iHe)/mass_per_H ! He mass fraction
end subroutine set_abundances

!-----------------------------------------------------------------------
!
!  set particle dust properties (particle's call)
!
!-----------------------------------------------------------------------
subroutine evolve_dust(dtsph, xyzh, u, JKmuS, Tdust, rho)
 use units,          only:utime,unit_density
 use eos,            only:ieos,get_temperature

 real,    intent(in) :: dtsph,Tdust,rho,u,xyzh(4)
 real,    intent(inout) :: JKmuS(:)

 integer, parameter :: wind_emitting_sink = 1
 real :: dt_cgs, T, rho_cgs, vxyzui(4)

 dt_cgs    = dtsph* utime
 rho_cgs   = rho*unit_density
 vxyzui(4) = u
 T         = get_temperature(ieos,xyzh,rho,vxyzui,gammai=JKmuS(idgamma),mui=JKmuS(idmu))
 call evolve_chem(dt_cgs, T, rho_cgs, JKmuS)
 JKmuS(idkappa)     = calc_kappa_dust(JKmuS(idK3), Tdust, rho_cgs)

end subroutine evolve_dust

!-----------------------------------------------------------------------
!
!  evolve the chemistry and moments
!
!-----------------------------------------------------------------------
subroutine evolve_chem(dt, T, rho_cgs, JKmuS)
!all quantities in cgs
 real, intent(in)    :: dt, rho_cgs
 real, intent(inout) :: T, JKmuS(:)

 real :: pC, pC2, pC2H, pC2H2, nH_tot, epsC
 real :: JstarS, taustar, taugr, S
 real :: Jstar_new, K_new(0:3)
 real :: nH, nH2, v1
 real, parameter :: A0 = 20.7d-16
 real, parameter :: alpha2 = 0.34

 nH_tot = rho_cgs/mass_per_H
 epsC   = eps(iC) - JKmuS(idK3)
 if (epsC < 0.) then
    print *,'eps(C) =',eps(iC),', K3=',JKmuS(idK3),', epsC=',epsC,', T=',T,' rho=',rho_cgs
    print *,'JKmuS=',JKmuS
    stop '[S-dust_formation] epsC < 0!'
 endif
 if (T > 450.) then
    call chemical_equilibrium_light(rho_cgs, T, epsC, pC, pC2, pC2H, pC2H2, JKmuS(idmu), JKmuS(idgamma))
    S = pC/psat_C(T)
    if (S > Scrit) then
       !call nucleation(T, pC, pC2, 0., pC2H, pC2H2, S, JstarS, taustar, taugr)
       call calc_nucleation(T, pC, 0., 0., 0., pC2H2, S, JstarS, taustar, taugr)
       JstarS = JstarS/ nH_tot
       call evol_K(JKmuS(idJstar), JKmuS(idK0:idK3), JstarS, taustar, taugr, dt, Jstar_new, K_new)
    else
       Jstar_new  = JKmuS(idJstar)
       K_new(0:3) = JKmuS(idK0:idK3)
    endif
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules, all O in CO
    nH  = 0.
    nH2 = nH_tot/2.
    JKmuS(idmu)    = (1.+4.*eps(iHe))*nH_tot/(nH+nH2+eps(iHe)*nH_tot)
    JKmuS(idgamma) = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
    pC2H2 = .5*(epsC-eps(iOx))*nH_tot * kboltz * T
    pC2H  = 0.
    S     = 1.d-3
    v1    = vfactor*sqrt(T)
    taugr = kboltz*T/(A0*v1*sqrt(2.)*alpha2*(pC2H+pC2H2))
    call evol_K(0., JKmuS(idK0:idK3), 0., 1., taugr, dt, Jstar_new, K_new)
 endif
 JKmuS(idJstar)   = Jstar_new
 JKmuS(idK0:idK3) = K_new(0:3)
 JKmuS(idsat)     = S

end subroutine evolve_chem

!------------------------------------
!
!  Bowen dust opacity formula
!
!------------------------------------
pure elemental real function calc_kappa_bowen(Teq)
!all quantities in cgs
 real,    intent(in)  :: Teq
 real :: dlnT

 dlnT = (Teq-bowen_Tcond)/bowen_delta
 if (dlnT > 50.) then
    calc_kappa_bowen = 0.
 else
    calc_kappa_bowen = bowen_kmax/(1.0 + exp(dlnT)) + kappa_gas
 endif

end function calc_kappa_bowen

!-----------------------------------------------------------------------
!
!  calculate dust opacity
!
!-----------------------------------------------------------------------
pure real function calc_kappa_dust(K3, Tdust, rho_cgs)
!all quantities in cgs
 real, intent(in) :: K3, Tdust, rho_cgs

 real :: kappa_cgs, fac
 real, parameter :: rho_Cdust = 2.62, mc = 12.*atomic_mass_unit
 real, parameter :: Qplanck_abs = 1.6846124267740528e+04
 real, parameter :: Qross_ext = 9473.2722815583311

 fac = max(0.75*K3*mc/(mass_per_H*rho_Cdust),1.e-15)
 !kappa_cgs = Qplanck_abs *fac ! planck
 !kappa_cgs = Qross_ext * fac  ! Rosseland

 ! Gail & Sedlmayr, 1985, A&A, 148, 183, eqs 23,24
 !kappa_cgs = 6.7d0 * fac * Tdust  ! planck
 kappa_cgs = 5.9d0 * fac * Tdust  ! Rosseland

 calc_kappa_dust = kappa_cgs + kappa_gas
end function calc_kappa_dust

!-----------------------------------------------------------------------
!
!  calculate alpha, reduced gravity factor
!
!-----------------------------------------------------------------------
pure real function calc_Eddington_factor(Mstar_cgs, Lstar_cgs, kappa_cgs, tau)
!all quantities in cgs
 use physcon, only:c,Gg
 real, intent(in) :: Mstar_cgs,Lstar_cgs,kappa_cgs
 real, intent(in), optional :: tau

 if (present(tau)) then
    calc_Eddington_factor = Lstar_cgs*exp(-tau)/(4.*pi*c*Gg*Mstar_cgs) * kappa_cgs
 else
    calc_Eddington_factor = Lstar_cgs/(4.*pi*c*Gg*Mstar_cgs) * kappa_cgs
 endif
end function calc_Eddington_factor

!----------------------------
!
!  Compute nucleation rate
!
!----------------------------
subroutine calc_nucleation(T, pC, pC2, pC3, pC2H, pC2H2, S, JstarS, taustar, taugr)
! all quantities are in cgs
 real, intent(in)  :: T, pC, pC2, pC3, pC2H, pC2H2, S
 real, intent(out) :: JstarS, taustar, taugr
 real, parameter   :: A0 = 20.7d-16
 real, parameter   :: sigma = 1400.
 real, parameter   :: theta_inf = A0*sigma/kboltz
 real, parameter   :: alpha1 = 0.37 !sticking coef for C
 real, parameter   :: alpha2 = 0.34 !sticking coef for C2,C2H,C2H2
 real, parameter   :: alpha3 = 0.08 !sticking coef for C3
 real, parameter   :: Nl_13 = 5.**(1./3.)
 real, parameter   :: mproton = 1.6726485d-24

 real :: ln_S_g, Nstar_inf_13, Nstar_m1_13, Nstar_m1_23, theta_Nstar,Nstar_m1,Nstar_tmp
 real :: dtheta_dNstar, d2lnc_dN2star, Z, A_Nstar, v1, beta, c_star, expon

 ln_S_g = log(S)
 v1     = vfactor*sqrt(T)
 Nstar_inf_13  = 2.*theta_inf/(3.*T*ln_S_g)
 Nstar_tmp     = 1. + sqrt(1. + 2.*Nl_13/Nstar_inf_13) - 2.*Nl_13/Nstar_inf_13
 if (Nstar_tmp > 0.) then
    Nstar_m1      = (Nstar_inf_13**3)/8. * Nstar_tmp**3
    Nstar_m1_13   = 0.5*Nstar_inf_13*Nstar_tmp
    Nstar_m1_23   = Nstar_m1_13**2
    theta_Nstar   = theta_inf/(1. + Nl_13/Nstar_m1_13)
    dtheta_dNstar = Nl_13 * theta_Nstar**2/(3.*theta_inf * Nstar_m1_23**2)
    d2lnc_dN2star = 2.*theta_Nstar/(9.*T*Nstar_m1)*(Nstar_m1_13+2.*Nl_13)/(Nstar_m1_13+Nl_13)**2
    Z             = sqrt(d2lnc_dN2star/(2.*pi))
    A_Nstar       = A0 * (1.+Nstar_m1)**(2./3.)
    beta          = v1/(kboltz*T) * (pC*alpha1 + 4.*alpha2/sqrt(2.)*(pC2 + pC2H + pC2H2) + 9.*alpha3/sqrt(3.)*pC3)
    expon         = Nstar_m1*ln_S_g - theta_Nstar*Nstar_m1_23/T
    if (expon < -100.) then
       c_star = 1.d-70
    else
       c_star = pC/(kboltz*T) * exp(expon)
    endif
    JstarS  = beta * A_Nstar * Z * c_star
    taustar = 1./(d2lnc_dN2star*beta*A_Nstar)
    ! if (isnan(JstarS)) then
    !   print*,i,'(N-1)^1/3=',Nstar_m1_13,'exp=',expon,'T=',T,'theta_N=',theta_Nstar,'d2lnc/dN2=',d2lnc_dN2star,ddd,&
    !        'beta=',beta,'Z=',Z,'c_star=',c_star,'JstarS=',JstarS,'tau*=',taustar
    !   if (isnan(JstarS)) stop
    ! endif
 else
    JstarS  = 0.d0
    taustar = 1.d-30
 endif
 taugr = kboltz*T/(A0*v1*(alpha1*pC*(1.-1./S) + 2.*alpha2/sqrt(2.)*(pC2+pC2H+pC2H2)*(1.-1./S**2)))
end subroutine calc_nucleation

!------------------------------------
!
!  Compute evolution of the moments
!
!------------------------------------
subroutine evol_K(Jstar, K, JstarS, taustar, taugr, dt, Jstar_new, K_new)
! all quantities are in cgs, K and K_new are the *normalized* moments (K/n<H>)
 real, intent(in) :: Jstar, K(0:3), JstarS, taustar, taugr, dt
 real, intent(out) :: Jstar_new, K_new(0:3)

 real, parameter :: Nl_13 = 10. !(lower grain size limit)**1/3
 real :: d, i0, i1, i2, i3, i4, i5, dK0, dK1, dK2, DK3

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
 dK0 = taustar*(Jstar*i1 + JstarS*i2)
 K_new(0) = K(0) + dK0
 dK1 = taustar**2/(3.*taugr)*(Jstar*i2 + JstarS*i3)
 K_new(1) = K(1) + K(0)*dt/(3.*taugr) + dK1 + Nl_13*dK0
 dK2 = 2.*taustar**3/(3.*taugr)**2 * (Jstar*i3 + JstarS*i4)
 K_new(2) = K(2) + 2.*K(1)*dt/(3.*taugr) + (dt/(3.*taugr))**2*K(0) + dK2 + 2.*Nl_13*dK1 + Nl_13**2*dK0
 dK3 = 3.*dt/(3.*taugr)*K(2) + 3.*(dt/(3.*taugr))**2*K(1) + (dt/(3.*taugr))**3*K(0)  &
     + (6.*taustar**4)/(3.*taugr)**3*(Jstar*i4+JstarS*i5)
 K_new(3) = K(3) + dK3 + Nl_13**3*dK0 + 3.*Nl_13**2*dK1 + 3.*Nl_13*dK2
 !if (any(isnan(K_new))) then
 !  print*,'NaNs in K_new for particle #',i
 !  print *,'dt=',dt,'tau*=',taustar,'taug=',taugr,'d=',d,'i0=',i0,'i1=',i1,'Jstar=',Jstar,'JstarS=',JstarS,&
 !      'k1=',k(1),'dk1=',dk1,'Kn1=',k_new(1),'k2=',k(2),'dk2=',dk2,'Kn2=',k_new(2),'k3=',k(3),'dk3=',dk3,'Kn3=',k_new(3)
 !  stop
 !endif

end subroutine evol_K

!----------------------------------------
!
!  Calculate mean molecular weight, gamma
!
!----------------------------------------
subroutine calc_muGamma(rho_cgs, T, mu, gamma, pH, pH_tot)
! all quantities are in cgs
 use io,  only:fatal
 use eos, only:ieos

 real, intent(in)    :: rho_cgs
 real, intent(inout) :: T, mu, gamma
 real, intent(out)   :: pH, pH_tot
 real :: KH2, pH2, x
 real :: T_old, mu_old, gamma_old, tol
 logical :: converged
 integer :: i,isolve
 integer, parameter :: itermax = 100
 character(len=30), parameter :: label = 'calc_muGamma'

 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
 T_old = T
 if (T > 1.d4) then
    mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
    pH     = pH_tot
    if (ieos /= 17) gamma  = 5./3.
 elseif (T > 450.) then
! iterate to get consistently pH, T, mu and gamma
    tol       = 1.d-3
    converged = .false.
    isolve    = 0
    pH        = pH_tot ! initial value, overwritten below, to avoid compiler warning
    i = 0
    do while (.not. converged .and. i < itermax)
       i = i+1
       pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H)
       KH2       = calc_Kd(coefs(:,iH2), T)
       pH        = solve_q(2.*KH2, 1., -pH_tot)
       pH2       = KH2*pH**2
       mu        = (1.+4.*eps(iHe))/(.5+eps(iHe)+0.5*pH/pH_tot)
       if (ieos == 17) exit !only update mu, keep gamma constant
       x         = 2.*(1.+4.*eps(iHe))/mu
       gamma     = (3.*x+4.+4.*eps(iHe))/(x+4.+4.*eps(iHe))
       converged = (abs(T-T_old)/T_old) < tol
       if (i == 1) then
          mu_old = mu
          gamma_old = gamma
       else
          T = T_old*mu/mu_old/(gamma_old-1.)*2.*x/(x+4.+4.*eps(iHe))
          if (i>=itermax .and. .not.converged) then
             if (isolve==0) then
                isolve = isolve+1
                i      = 0
                tol    = 1.d-2
                print *,'[dust_formation] cannot converge on T(mu,gamma). Trying with lower tolerance'
             else
                print *,'Told=',T_old,',T=',T,',gamma_old=',gamma_old,',gamma=',gamma,',mu_old=',&
                  mu_old,',mu=',mu,',dT/T=',abs(T-T_old)/T_old,', rho=',rho_cgs
                call fatal(label,'cannot converge on T(mu,gamma)')
             endif
          endif
       endif
    enddo
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2    = pH_tot/2.
    pH     = 0.
    mu     = (1.+4.*eps(iHe))/(0.5+eps(iHe))
    if (ieos /= 17)  gamma  = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
 endif

end subroutine calc_muGamma

!--------------------------------------------
!
!  Initialise mean molecular weight and gamma
!
!--------------------------------------------
subroutine init_muGamma(rho_cgs, T, mu, gamma, ppH, ppH2)
! all quantities are in cgs
 real, intent(in)              :: rho_cgs
 real, intent(inout)           :: T
 real, intent(out)             :: mu, gamma
 real, intent(out), optional   :: ppH, ppH2
 real :: KH2, pH_tot, pH, pH2

 pH_tot = rho_cgs*kboltz*T/(patm*mass_per_H)
 if (T > 1.d5) then
    pH2 = 0.
    pH  = pH_tot
 elseif (T > 450.) then
    KH2 = calc_Kd(coefs(:,iH2), T)
    pH  = solve_q(2.*KH2, 1., -pH_tot)
    pH2 = KH2*pH**2
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2 = pH_tot/2.
    pH  = 0.
 endif
 mu    = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
 gamma = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
 call calc_muGamma(rho_cgs, T, mu, gamma, pH, pH_tot)
 if (present(ppH))  ppH = pH
 if (present(ppH2)) ppH2 = pH2

end subroutine init_muGamma

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance in the gas phase
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light(rho_cgs, T, epsC, pC, pC2, pC2H, pC2H2, mu, gamma,&
     nH, nH2, nHe, nCO, nH2O, nOH)
! all quantities are in cgs
 real, intent(in)    :: rho_cgs, epsC
 real, intent(inout) :: T, mu, gamma
 real, intent(out)   :: pC, pC2, pC2H, pC2H2
 real, intent(out), optional :: nH, nH2, nHe, nCO, nH2O, nOH
 real    :: pH_tot, Kd(nMolecules+1), err, a, b, c, d
 real    :: pH, pCO, pO, pSi, pS, pTi, pN
 real    :: pC_old, pO_old, pSi_old, pS_old, pTi_old, cst
 integer :: i, nit

 call calc_muGamma(rho_cgs, T, mu, gamma, pH, pH_tot)
 if (T > 1.d4) then
    pC    = eps(iC)*pH_tot
    pC2   = 0.
    pC2H  = 0.
    pC2H2 = 0.
    if (present(nH)) then
       nH   = pH  *(patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
       nH2  = 1.d-50
       nHe  = eps(ihe)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
       nCO  = 1.d-50
       nH2O = 1.d-50
       nOH  = 1.d-50
    endif
    return
 endif

! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)
 pCO      = epsC*pH_tot

 pC       = 0.
 pC_old   = 0.
 pO       = 0.
 pO_old   = 0.
 pSi      = 0.
 pSi_old  = 0.
 pS       = 0.
 pS_old   = 0.
 pTi      = 0.
 pTi_old  = 0.
 err      = 1.
 nit      = 0

 do while (err > 1.e-5)
! N
    pN  = solve_q(2.*Kd(iN2), &
            1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+Kd(iHCN)*pH), &
            -eps(iN)*pH_tot)
! C
    pC  = solve_q(2.*(Kd(iC2H)*pH + Kd(iC2H2)*pH**2 + Kd(iC2)), &
            1.+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C), &
            pCO-epsC*pH_tot)
! O
    pO  = eps(iOx)*pH_tot/(1.+pC*Kd(iCO)+pH*Kd(iOH)+pH**2*Kd(iH2O)+2.*pO*pC*Kd(iCO2)+pSi*Kd(iSiO))
    pCO = Kd(iCO)*pC*pO
! Si
    a   = 3.*Kd(iSi3)
    b   = 2.*Kd(iSi2)+2.*Kd(iSi2C)*pC
    c   = 1.+Kd(iSiH)*pH+Kd(iSiH4)*pH**4+Kd(iSiO)*pO+Kd(iSiS)*pS
    d   = -eps(iSi)*pH_tot
    pSi = -d/(a*pSi**2+b*pSi+c)
! S
    pS  = solve_q(2.*Kd(iS2), &
            1.+Kd(iSiS)*pSi+Kd(iHS)*pH+Kd(iH2S)*pH**2, &
            -eps(iS)*pH_tot)
! Ti
    pTi = eps(iTi)*pH_tot/(1.+Kd(iTiO)*pO+Kd(iTiO2)*pO**2+Kd(iTiS)*pS)

    err = 0.
    if (pC  > 1.e-50) err = err + abs((pC-pC_old)/pC)
    if (pO  > 1.e-50) err = err + abs((pO-pO_old)/pO)
    if (pSi > 1.e-50) err = err + abs((pSi-pSi_old)/pSi)
    if (pS  > 1.e-50) err = err + abs((pS-pS_old)/pS)
    if (pTi > 1.e-50) err = err + abs((pTi-pTi_old)/pTi)

    nit = nit + 1
    if (nit == 200) exit

    pC_old  = pC
    pO_old  = pO
    pSi_old = pSi
    pS_old  = pS
    pTi_old = pTi
 enddo

 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2

! Convert partial pressures from atm to cgs
 pC    = pC*patm
 pC2   = pC2*patm
 pC2H  = pC2H*patm
 pC2H2 = pC2H2*patm
 if (present(nH)) then
    cst  = mass_per_H/(mu*mass_proton_cgs*kboltz*T)
    if (T < 450.) then
       nH2 = pH_tot/2.      *patm*cst
       nH  = 1.d-99
    else
       nH   = pH            *patm*cst
       nH2  = Kd(iH2)*pH**2 *patm*cst
    endif
    nHe  = eps(ihe)*pH_tot  *patm*cst
    nCO  = Kd(iCO) *pC*pO   *patm*cst
    nH2O = Kd(iH2O)*pH**2*pO*patm*cst
    nOH  = Kd(iOH) *pH*pO   *patm*cst
 endif
end subroutine chemical_equilibrium_light

!-----------------------------
!
!  solve 2nd order polynomial
!
!-----------------------------
pure real function solve_q(a, b, c)
 real, intent(in) :: a, b, c
 real :: delta
 if (-4.*a*c/b**2 > epsilon(0.)) then
    delta = max(b**2-4.*a*c, 0.)
    solve_q = (-b+sqrt(delta))/(2.*a)
 else
    solve_q = -c/b
 endif
 solve_q = max(solve_q,1e-50)
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
 real :: T2,T3,pC1!,pC2,pC3,pC4,pC5

 if (T > 1.d4) then
    Psat_C = huge(Psat_C)
 else
    T2 = T*T
    T3 = T*T*T
    pC1 = exp(-8.61240d4/T + 1.85106d1 + 5.87980d-4*T - 2.51549d-7*T2 + 3.24892d-11*T3 + f)
    !pC2 = exp(-1.01124d5/T + 2.35611d1 + 3.37807d-4*T - 2.94959d-7*T2 + 4.41801D-11*T3 + f)
    !pC3 = exp(-9.90261d4/T + 2.81161d1 - 1.55820d-3*T + 1.60002d-7*T2 - 4.47171D-12*T3 + f)
    !pC4 = exp(-1.17037d5/T + 2.55579d1 - 5.63869d-6*T - 2.13596d-7*T2 + 3.39660D-11*T3 + f)
    !pC5 = exp(-1.18080d5/T + 2.65798d1 + 1.20285d-4*T - 2.68583d-7*T2 + 4.12365D-11*T3 + f)
    psat_C = pC1
 endif
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
 d = min(-G/(R*T),222.)
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

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_dust_formation(hdr,ierr)
 use dump_utils,        only:dump_h,add_to_rheader
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out)   :: ierr

! initial gas composition for dust formation
 call set_abundances
 call add_to_rheader(eps,'epsilon',hdr,ierr) ! array
 call add_to_rheader(Aw,'Amean',hdr,ierr)    ! array
 call add_to_rheader(mass_per_H,'mass_per_H',hdr,ierr) ! real

end subroutine write_headeropts_dust_formation

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_dust_formation(hdr,ierr)
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 real :: dum(nElements)


 ierr = 0
 call extract('mass_per_H',mass_per_H,hdr,ierr) ! real
 ! it is likely that your dump was generated with an old version of phantom
 ! and the chemical properties not stored. restore and save the default values
 if (mass_per_H < tiny(0.)) then
    print *,'reset dust chemical network properties'
    call set_abundances
    call extract('epsilon',dum(1:nElements),hdr,ierr) ! array
    call extract('Amean',dum(1:nElements),hdr,ierr) ! array
 else
    call extract('epsilon',eps(1:nElements),hdr,ierr) ! array
    call extract('Amean',Aw(1:nElements),hdr,ierr) ! array
 endif


end subroutine read_headeropts_dust_formation

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_dust_formation(iunit)
 use dim,          only:nucleation
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling dust'
 if (nucleation) then
    call write_inopt(idust_opacity,'idust_opacity','compute dust opacity (0=off, 1=bowen, 2=nucleation)',iunit)
 else
    call write_inopt(idust_opacity,'idust_opacity','compute dust opacity (0=off, 1=bowen)',iunit)
 endif
 if (idust_opacity == 1) then
    call write_inopt(kappa_gas,'kappa_gas','constant gas opacity (cm²/g)',iunit)
    call write_inopt(bowen_kmax,'bowen_kmax','maximum dust opacity (cm²/g)',iunit)
    call write_inopt(bowen_Tcond,'bowen_Tcond','dust condensation temperature (K)',iunit)
    call write_inopt(bowen_delta,'bowen_delta','condensation temperature range (K)',iunit)
 endif
 if (nucleation .and. idust_opacity == 2) then
    call write_inopt(kappa_gas,'kappa_gas','constant gas opacity (cm²/g)',iunit)
    call write_inopt(wind_CO_ratio ,'wind_CO_ratio','wind initial C/O ratio (> 1)',iunit)
 endif

end subroutine write_options_dust_formation

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_dust_formation(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use dim,     only:do_nucleation,inucleation,store_dust_temperature
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_nucleation'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('idust_opacity')
    read(valstring,*,iostat=ierr) idust_opacity
    ngot = ngot + 1
    if (idust_opacity == 2) then
       do_nucleation = .true.
       inucleation = 1
    endif
 case('wind_CO_ratio')
    read(valstring,*,iostat=ierr) wind_CO_ratio
    ngot = ngot + 1
    if (wind_CO_ratio < 0.) call fatal(label,'invalid setting for wind_CO_ratio (must be > 0)')
 case('kappa_gas')
    read(valstring,*,iostat=ierr) kappa_gas
    ngot = ngot + 1
    if (kappa_gas < 0.)    call fatal(label,'invalid setting for kappa_gas (<0)')
    !kgas = kappa_gas / (udist**2/umass)
 case('bowen_kmax')
    read(valstring,*,iostat=ierr) bowen_kmax
    ngot = ngot + 1
    if (bowen_kmax < 0.)    call fatal(label,'invalid setting for bowen_kmax (<0)')
 case('bowen_Tcond')
    read(valstring,*,iostat=ierr) bowen_Tcond
    ngot = ngot + 1
    if (bowen_Tcond < 0.) call fatal(label,'invalid setting for bowen_Tcond (<0)')
 case('bowen_delta')
    read(valstring,*,iostat=ierr) bowen_delta
    ngot = ngot + 1
    if (bowen_delta < 0.) call fatal(label,'invalid setting for bowen_delta (<0)')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)
 if (idust_opacity == 1) igotall = (ngot >= 5)
 if (idust_opacity == 2) igotall = (ngot >= 3)
 if (idust_opacity > 0) store_dust_temperature = .true.

end subroutine read_options_dust_formation

end module dust_formation
