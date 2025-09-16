!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
 use dim,     only:nElements,nabn_AGB

 implicit none
 integer, public :: idust_opacity = 0

 public :: set_abundances,evolve_dust,evolve_chem,calc_kappa_dust,&
      calc_kappa_bowen,chemical_equilibrium_light,psat_C,calc_nucleation,&
      read_options_dust_formation,write_options_dust_formation,&
      calc_Eddington_factor,calc_muGamma,init_muGamma,init_nucleation,&
      write_headeropts_dust_formation,read_headeropts_dust_formation,&
      chemical_equilibrium_light_fixed_mu_gamma,&
      chemical_equilibrium_light_fixed_mu_gamma_broyden, &
      chemical_equilibrium_Cspecies
!
!--runtime settings for this module
!

 real, public :: kappa_gas   = 2.d-4
 real, public :: wind_CO_ratio = 1.7 
 real, public :: Tmol = 450.
 real, public, parameter :: Scrit = 1. ! Critical saturation ratio
 real, public :: mass_per_H, eps(nElements)
 real, public :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
 integer, public :: icoolH=1, icoolC=2, icoolO=3, icoolSi=4, icoolH2=5, icoolCO=6, &
                       icoolH2O=7, icoolOH=8, icoolC2=9, icoolC2H=10, icoolC2H2=11, &
                       icoolHe=12, icoolSiO=13, icoolCH4=14, icoolS=15, icoolTi=16

 private

 character(len=*), parameter :: label = 'dust_formation'
 real :: bowen_kmax  = 2.7991
 real :: bowen_Tcond = 1500.
 real :: bowen_delta = 60.

! Indices for elements and molecules:
 integer, parameter :: nMolecules = 25
 integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
 integer, parameter :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, &
      iNH3=10, iCN=11, iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, &
      iHS=19, iH2S=20, iSiS=21, iSiH=22, iTiO=23, iTiO2=24,iC2 = 25, iTiS=26
! Indices for cooling species:
 
 real(kind=16), parameter :: coefs(5,nMolecules) = reshape([&
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
       3.65656d+05, -8.77172d+04, 2.41557d+01, 3.40917d-04, -2.09217d-08, & !HS
       7.41911d+05, -1.81063d+05, 5.35114d+01, 4.94594d-04, -3.39231d-08, & !H2S
       1.44507d+05, -1.49916d+05, 2.87381d+01, 3.96186d-04, -2.18002d-08, & !SiS
       2.45408d+05, -7.17752d+04, 2.29449d+01, 4.52999d-04, -2.69915d-08, & !SiH
      -4.38897d+05, -1.58111d+05, 2.49224d+01, 1.08714d-03, -5.62504d-08, & !TiO
      -3.32351d+05, -3.04694d+05, 5.86984d+01, 1.17096d-03, -5.06729d-08, & !TiO2
       2.26786d+05, -1.43775d+05, 2.92429d+01, 1.69434d-04, -1.79867d-08], shape(coefs)) !C2
!  real, parameter :: vfactor = sqrt(kboltz/(2.*pi*atomic_mass_unit*12.01))
 real, parameter :: vfactor = sqrt(kboltz/(8.*pi*atomic_mass_unit*12.01))

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
 real :: cst
 real :: abundi(nabn_AGB)
 ! Needed to write output abundances to external file:
 integer :: iunit_out, j
 character(len=100) :: filename
 logical :: file_exists

 nH_tot = rho_cgs/mass_per_H
 epsC   = eps(iC) - JKmuS(idK3)
 if (epsC < 0.) then
    print *,'eps(C) =',eps(iC),', K3=',JKmuS(idK3),', epsC=',epsC,', T=',T,' rho=',rho_cgs
    print *,'JKmuS=',JKmuS
    stop '[S-dust_formation] epsC < 0!'
 endif
!  if (T> 3000.) then
!     K_new = 0.
!     Jstar_new = 0.
!     S = 0.
!     JKmuS(idkappa) = 0.
 if (T > 450.) then
    call chemical_equilibrium_light(rho_cgs, T, epsC, JKmuS(idmu), JKmuS(idgamma), abundi)
    cst = mass_per_H/(JKmuS(idmu)*mass_proton_cgs*kboltz*T)
    pC = abundi(icoolC)/cst 
    pC2 = abundi(icoolC2)/cst 
    pC2H  = abundi(icoolC2H)/cst 
    pC2H2  = abundi(icoolC2H2)/cst 
    S = pC/psat_C(T)

   ! ! Write output abundances to external file
   ! filename = 'abundances_output.txt'
   ! inquire(file=filename, exist=file_exists)
   ! if (.not. file_exists) then
   !    open(newunit=iunit_out, file=filename, status='unknown', action='write', position='append')
   ! else
   !    open(newunit=iunit_out, file=filename, status='old', action='write', position='append')
   ! endif
   ! write(iunit_out, '(100(1p,es14.6e3,1x))') (abundi(j)/nH_tot, j=1,nabn_AGB), psat_C(T)*cst/nH_tot
   ! close(iunit_out)


    if (S > Scrit) then
       call calc_nucleation(T, pC, pC2, 0.0, pC2H, pC2H2, S, JstarS, taustar, taugr)
       JstarS = JstarS/ nH_tot
       call evol_K(JKmuS(idJstar), JKmuS(idK0:idK3), JstarS, taustar, taugr, dt, Jstar_new, K_new)
    else
       if (any(JKmuS(idK0:idK3) > 0.0)) then
          call calc_nucleation(T, pC, pC2, 0.0, pC2H, pC2H2, S, JstarS, taustar, taugr)
          JstarS = JstarS / nH_tot
          call evol_K_ev(JKmuS(idJstar), JKmuS(idK0:idK3), taustar, taugr, dt, Jstar_new, K_new)
       else
          Jstar_new = 0.0
          K_new(0:3) = JKmuS(idK0:idK3)
       endif
   !  else
   !     Jstar_new  = JKmuS(idJstar)
   !     K_new(0:3) = JKmuS(idK0:idK3)
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
    taugr = kboltz*T/(A0*v1*sqrt(2.)*alpha2*(pC2+pC2H+pC2H2))
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

 v1     = vfactor*sqrt(T)
 ln_S_g = abs(log(S))

!  if (ln_S_g <= 0.d0) then
!     JstarS  = 0.d0
!     taustar = 1.d-30
!     taugr = kboltz*T / (A0*v1 * ( alpha1*pC*(1 - 1/S) + &
!              2*alpha2/sqrt(2.)*(pC2 + pC2H + pC2H2)*(1 - 1/S**2) ))
!     return
!  endif
 
 Nstar_inf_13  = 2.*theta_inf/(3.*T*ln_S_g)
 Nstar_tmp     = 1. + sqrt(1. + 2.*Nl_13/Nstar_inf_13) - 2.*Nl_13/Nstar_inf_13
 if (Nstar_tmp > 0.) then
    Nstar_m1      = (Nstar_inf_13**3)/8. * Nstar_tmp**3
    Nstar_m1_13   = 0.5*Nstar_inf_13*Nstar_tmp
    Nstar_m1_23   = Nstar_m1_13**2
    theta_Nstar   = theta_inf/(1. + Nl_13/Nstar_m1_13)
    dtheta_dNstar = Nl_13 * theta_Nstar**2/(3.*theta_inf * Nstar_m1_23**2)
    d2lnc_dN2star = 2.*theta_Nstar/(9.*T*Nstar_m1)*(Nstar_m1_13+2.*Nl_13)/(Nstar_m1_13+Nl_13)**2
   !  d2lnc_dN2star = -2./T * Nstar_m1_23 * (dtheta_dNstar**2/theta_Nstar - theta_Nstar/(9.*Nstar_m1**2))
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
    taustar = 1.d-50
 endif
 taugr = kboltz*T/(A0*v1*(alpha1*pC*(1.-1./S) + 2.*alpha2/sqrt(2.)*(pC2+pC2H+pC2H2)*(1.-1./S**2)))
end subroutine calc_nucleation

!------------------------------------
!
!  Compute evolution of the moments
!
!------------------------------------
subroutine evol_K(Jstar_in, K, JstarS_in, taustar, taugr, dt, Jstar_new, K_new)
! all quantities are in cgs, K and K_new are the *normalized* moments (K/n<H>)
 real, intent(in) :: Jstar_in, K(0:3), JstarS_in, taustar, taugr, dt
 real, intent(out) :: Jstar_new, K_new(0:3)

 real, parameter :: Nl_13 = 10. !(lower grain size limit)**1/3
 real :: d, i0, i1, i2, i3, i4, i5, dK0, dK1, dK2, dK3
 real :: Jstar, JstarS
!  real(kind=16)  :: dK3

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
 Jstar_new = Jstar_in*i0 + JstarS_in*i1
 ! When Jstar or JstarS are < 1.d-50 (around 1d-70), they make the moments grow in an unrealistic way,
 ! for example we are unable to reproduce the average radius.
 ! The correct values to be considered (to reproduce Gauger 1990) are only larger than 1d-50
 if (Jstar_in < 1.d-50) then
    Jstar = 0.0
 else
    Jstar = Jstar_in
 endif
 if (JstarS_in < 1.d-50) then
    JstarS = 0.0
 else
    JstarS = JstarS_in
 endif
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
K_new = max(K_new, 0.0)

end subroutine evol_K

!------------------------------------------------------------
!
!  Compute evolution of the moments during dust evaporation
!
!------------------------------------------------------------
subroutine evol_K_ev(Jstar, K, taustar, taugr, dt, Jstar_new, K_new)
! all quantities are in cgs, K and K_new are the *normalized* moments (K/n<H>)
 implicit none
 real, intent(in) :: K(0:3), Jstar, taustar, taugr, dt
 real, intent(out) :: Jstar_new, K_new(0:3)

 real, parameter :: Nl_13 = 10.0
 real :: J_ev_new, J_ev_exp, dK0, dK1, dK2, dK3

 J_ev_exp = dt/taustar
 J_ev_exp = min(J_ev_exp, 50.0)
 Jstar_new = Jstar * exp(J_ev_exp)
 J_ev_new = -Jstar_new
 
 dK0 = J_ev_new * dt
 K_new(0) = K(0) + dK0
 dK1 = J_ev_new * dt**2 / (6.0*taugr)
 K_new(1) = K(1) + K(0)*dt/(3.0*taugr) + dK1 + Nl_13*dK0
 dK2 = J_ev_new * dt**3 / (27.0*taugr**2)
 K_new(2) = K(2) + 2.0*K(1)*dt/(3.0*taugr) + (dt/(3.0*taugr))**2*K(0) + dK2 + 2.0*Nl_13*dK1 + Nl_13**2*dK0
 dK3 = J_ev_new * dt**4 / (108.0*taugr**3)
 K_new(3) = K(3) + dt/taugr*K(2) + (dt/taugr)**2*K(1)/3.0 + (dt/(3.0*taugr))**3*K(0) + dK3 + &
            Nl_13**3*dK0 + 3.0*Nl_13**2*dK1 + 3.0*Nl_13*dK2

 K_new = max(K_new, 0.0)

end subroutine evol_K_ev

!----------------------------------------
!
!  Calculate mean molecular weight, gamma
!
!----------------------------------------
subroutine calc_muGamma(rho_cgs, T, mu, gamma, pH_out, pH_tot_out, ppH2)
! all quantities are in cgs
 use io,  only:fatal
 use eos, only:ieos

 real, intent(in)    :: rho_cgs
 real, intent(inout) :: T, mu, gamma
 real, intent(out)   :: pH_out, pH_tot_out
 real, intent(out), optional :: ppH2
 real(kind=16) :: KH2, pH2, x
 real :: T_old, mu_old, gamma_old, tol, mu_original
 logical :: converged
 integer :: i,isolve
 integer, parameter :: itermax = 100
 character(len=30), parameter :: label = 'calc_muGamma'
 real(kind=16) :: pH, pH_tot


 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
 T_old = T
 mu_original = mu
 if (T > 1.d4) then
    mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
    pH     = pH_tot
   !  pH2    = 1.d-50
    if (ieos /= 17) gamma  = 5./3.
 elseif (T > Tmol) then
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
       pH        = solve_q(2.*KH2, 1._16, -pH_tot)
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

!  pH        = solve_q(2.*KH2, 1., -pH_tot)
 if (present(ppH2)) ppH2 = pH2

 pH_out = real(pH, kind=8)
 pH_tot_out = real(pH_tot, kind=8)

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
 real(kind=16) :: KH2, pH_tot, pH, pH2, mu_old
 real          :: pH_double, pH_tot_double


 pH_tot = rho_cgs*kboltz*T/(patm*mass_per_H)
 if (T > 1.d5) then
    pH2 = 0.
    pH  = pH_tot
 elseif (T > Tmol) then
    KH2 = calc_Kd(coefs(:,iH2), T)
    pH  = solve_q(2.*KH2, 1._16, -pH_tot)
    pH2 = KH2*pH**2
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2 = pH_tot/2.
    pH  = 0.
 endif
 mu    = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
 gamma = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
 call calc_muGamma(rho_cgs, T, mu, gamma, pH_double, pH_tot_double)
 if (present(ppH))  ppH = pH_double
 if (present(ppH2)) ppH2 = pH2


end subroutine init_muGamma

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance in the gas phase
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light(rho_cgs, T_in, epsC, mu, gamma, abundi)
! all quantities are in cgs
 real, intent(in)    :: rho_cgs, T_in, epsC
 real, intent(inout) :: mu, gamma
 real, intent(out)   :: abundi(nabn_AGB)
 real(kind=16)   :: pC, pC2, pC2H, pC2H2
 real(kind=16)    :: pH_tot, err, a, b, c, d
 real(kind=16) :: Kd(nMolecules+1)
 real(kind=16)    :: pH, pCO, pO, pSi, pS, pTi, pN, pH2, pSiO, pCH4
 real(kind=16)    :: pN_old, pC_old, pO_old, pSi_old, pS_old, pTi_old
 real(kind=16)    :: cst
 integer :: i, nit
 real(kind=16)    :: X, AA, BB
 real             :: T
 real             :: pH_mugamma, pH_tot_mugamma, pH2_mugamma

 T = max(T_in, 20.d0)

 pH_mugamma = real(pH, kind=8)
 pH_tot_mugamma = real(pH_tot, kind=8)
 pH2_mugamma = real(pH2, kind=8)
 call calc_muGamma(rho_cgs, T, mu, gamma, pH_mugamma, pH_tot_mugamma, pH2_mugamma)
 pH = real(pH_mugamma, kind=16)
 pH_tot = real(pH_tot_mugamma, kind=16)
 pH2 = real(pH2_mugamma, kind=16)
 if (T > 1.d4) then
    abundi(icoolC)    = eps(iC)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolC2)   = 0.
    abundi(icoolC2H)  = 0.
    abundi(icoolC2H2) = 0.
    abundi(icoolH)   = pH  *(patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolH2)  = 1.d-50
    abundi(icoolHe)  = eps(ihe)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolCO)  = 1.d-50
    abundi(icoolH2O) = 1.d-50
    abundi(icoolOH)  = 1.d-50
    abundi(icoolO)   = eps(iOx)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolSi)  = eps(iSi)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolSiO) = 1.d-50
    abundi(icoolCH4) = 1.d-50
    abundi(icoolS)   = eps(iS)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolTi)  = eps(iTi)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    return
 endif
 
 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)

! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)
 pCO      = epsC*pH_tot
 pH       = solve_q(2.*Kd(iH2), 1._16, -pH_tot)

 pN       = 0.
 pN_old   = 0.
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

 do while (err > 1.d-6)
! N
    X = 1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+pH*Kd(iHCN))
    AA = .25*X/Kd(iN2)
    BB = 8.*eps(iN)*pH_tot*Kd(iN2)/X**2
    if (BB > 1.d-8) then
       pN = AA*(sqrt(1.+BB)-1.)  
    else 
       pN = 0.5*AA*BB
    endif
! C
   !  pC  = solve_q(2.*(Kd(iC2H)*pH + Kd(iC2H2)*pH**2 + Kd(iC2)), &
   !          1.+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C)+ &
   !          pO*Kd(iCO)+pN*Kd(iCN)+pH*pN*Kd(iHCN), &
   !          -epsC*pH_tot)
    X = 1.+pO*Kd(iCO)+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C) + &
         pN*Kd(iCN)+pH*pN*Kd(iHCN)
    AA = .25*X/(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))
    BB = 8.*(epsC*pH_tot)*(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))/(X*X)
    if (BB > 1.d-8) then
       pC = AA*(sqrt(1.+BB)-1.)
    else 
       pC = 0.5*AA*BB
    endif
! O
    pO  = eps(iOx)*pH_tot/(1.+pC*Kd(iCO)+pH*Kd(iOH)+pH**2*Kd(iH2O)+ &
          2.*pO*pC*Kd(iCO2)+pSi*Kd(iSiO) &
          +pTi*Kd(iTiO)+2.*pO*pTi*Kd(iTiO2))
   !  pO  = solve_q(2.*(pC*Kd(iCO2)+pTi*Kd(iTiO2)), & 
   !          1.+pH*Kd(iOH)+pH**2*Kd(iH2O)+pC*Kd(iCO)+pSi*Kd(iSiO)+ &
   !          pTi*Kd(iTiO), &
   !          -eps(iOx)*pH_tot)
    pCO = Kd(iCO)*pC*pO
! Si
    a   = 3.*Kd(iSi3)
    b   = 2.*Kd(iSi2)+2.*Kd(iSi2C)*pC
    c   = 1.+Kd(iSiH)*pH+Kd(iSiH4)*pH**4+Kd(iSiO)*pO+Kd(iSiS)*pS
    d   = -eps(iSi)*pH_tot
    pSi = -d/(a*pSi**2+b*pSi+c)
! S
   !  pS  = solve_q(2.*Kd(iS2), &
   !          1.+Kd(iSiS)*pSi+Kd(iHS)*pH+Kd(iH2S)*pH**2+Kd(iTiS)*pTi, &
   !          -eps(iS)*pH_tot)
    X = 1.+pSi*Kd(iSiS)+pH*Kd(iHS)+pH**2*Kd(iH2S)+pTi*Kd(iTiS)
    AA = .25*X/Kd(iS2)
    BB = 8.*eps(iS)*pH_tot*Kd(iS2)/X**2
    if (BB > 1.d-8) then
       pS = AA*(sqrt(1.+BB)-1.)
    else
       pS = 0.5*AA*BB
    endif
! Ti
    pTi = eps(iTi)*pH_tot/(1.+Kd(iTiO)*pO+Kd(iTiO2)*pO**2+Kd(iTiS)*pS)
   !  pTi = solve_q(0._16, &
   !          1.+Kd(iTiO)*pO+Kd(iTiS)*pS+Kd(iTiO2)*pO**2, &
   !          -eps(iTi)*pH_tot)

    err = 0.
    if (pN  > 1.d-50) err = err + abs((pN-pN_old)/pN)
    if (pC  > 1.d-50) err = err + abs((pC-pC_old)/pC)
    if (pO  > 1.d-50) err = err + abs((pO-pO_old)/pO)
    if (pSi > 1.d-50) err = err + abs((pSi-pSi_old)/pSi)
    if (pS  > 1.d-50) err = err + abs((pS-pS_old)/pS)
    if (pTi > 1.d-50) err = err + abs((pTi-pTi_old)/pTi)

    nit = nit + 1
    if (nit == 2000) exit

    pN_old  = pN
    pC_old  = pC
    pO_old  = pO
    pSi_old = pSi
    pS_old  = pS
    pTi_old = pTi
 enddo
 
 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2
 pCH4  = Kd(iCH4)*pC*pH**4
 pSiO  = Kd(iSiO)*pO*pSi

 cst = mass_per_H/(mu*mass_proton_cgs*kboltz*T)
 abundi(icoolH)   = pH               *patm*cst
 abundi(icoolH2)  = Kd(iH2)*pH**2    *patm*cst
 abundi(icoolHe)  = eps(ihe)*pH_tot  *patm*cst  ! pH_tot is not changing, but helium probably changes following change in mu
 abundi(icoolCO)  = Kd(iCO)*pC*pO    *patm*cst
 abundi(icoolH2O) = Kd(iH2O)*pH**2*pO*patm*cst
 abundi(icoolOH)  = Kd(iOH) *pH*pO   *patm*cst
 abundi(icoolO)   = pO               *patm*cst
 abundi(icoolSi)  = pSi              *patm*cst
 abundi(icoolC2)  = pC2              *patm*cst
 abundi(icoolC)   = pC               *patm*cst
 abundi(icoolC2H2) = pC2H2           *patm*cst
 abundi(icoolC2H) = pC2H             *patm*cst
 abundi(icoolSiO) = pSiO             *patm*cst
 abundi(icoolCH4) = pCH4             *patm*cst

 ! These abundances are number densities, so in units of cm^{-3}

!  endif


end subroutine chemical_equilibrium_light

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance with fixed mu and gamma
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light_fixed_mu_gamma(rho_cgs, T_in, epsC, mu, gamma, abundi)
! all quantities are in cgs
 real, intent(in)    :: rho_cgs, T_in, epsC
 real, intent(in) :: mu, gamma
 real, intent(out)   :: abundi(nabn_AGB)
 real(kind=16)   :: pC, pC2, pC2H, pC2H2
 real(kind=16)    :: pH_tot, err, a, b, c, d
 real(kind=16) :: Kd(nMolecules+1)
 real(kind=16)    :: pH, pCO, pO, pSi, pS, pTi, pN, pH2, pSiO, pCH4
 real(kind=16)    :: pN_old, pC_old, pO_old, pSi_old, pS_old, pTi_old
 real(kind=16)    :: cst
 integer :: i, nit
 real(kind=16)    :: X, AA, BB
 real             :: T
 real             :: pH_mugamma, pH_tot_mugamma, pH2_mugamma

 T = max(T_in, 20.d0)

 cst = mass_per_H/(mu*mass_proton_cgs*kboltz*T)
 ! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
 pH     = solve_q(2.*Kd(iH2), 1._16, -pH_tot)

!  pH_mugamma = real(pH, kind=8)
!  pH_tot_mugamma = real(pH_tot, kind=8)
!  pH2_mugamma = real(pH2, kind=8)
!  call calc_muGamma(rho_cgs, T, mu, gamma, pH_mugamma, pH_tot_mugamma, pH2_mugamma)
!  pH = real(pH_mugamma, kind=16)
!  pH_tot = real(pH_tot_mugamma, kind=16)
!  pH2 = real(pH2_mugamma, kind=16)
 if (T > 1.d4) then
    abundi(icoolC)    = eps(iC)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolC2)   = 0.
    abundi(icoolC2H)  = 0.
    abundi(icoolC2H2) = 0.

    abundi(icoolH)   = pH  *(patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolH2)  = 1.d-50
    abundi(icoolHe)  = eps(ihe)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolCO)  = 1.d-50
    abundi(icoolH2O) = 1.d-50
    abundi(icoolOH)  = 1.d-50
    abundi(icoolO)   = eps(iOx)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolSi)  = eps(iSi)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolSiO) = 1.d-50
    abundi(icoolCH4) = 1.d-50
    abundi(icoolS)   = eps(iS)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolTi)  = eps(iTi)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    return
 endif


 Kd(iTiS) = calc_Kd_TiS(T)
 pCO      = epsC*pH_tot

!  pC       = abundi(icoolC) / (patm*cst)
!  pC_old   = abundi(icoolC) / (patm*cst)
!  pO       = abundi(icoolO) / (patm*cst)
!  pO_old   = abundi(icoolO) / (patm*cst)
!  pSi      = abundi(icoolSi) / (patm*cst)
!  pSi_old  = abundi(icoolSi) / (patm*cst)
!  pS       = abundi(icoolS) / (patm*cst)
!  pS_old   = abundi(icoolS) / (patm*cst)
!  pTi      = abundi(icoolTi) / (patm*cst)
!  pTi_old  = abundi(icoolTi) / (patm*cst)
 pN       = 0.
 pN_old   = 0.
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
 err      = 1.d0
 nit      = 0

 do while (err > 1.d-6)
! N
!  pN  = solve_q(2.*Kd(iN2), &
!             1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+Kd(iHCN)*pH), &
!             -eps(iN)*pH_tot)
    X = 1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+pH*Kd(iHCN))
    AA = .25*X/Kd(iN2)
    BB = 8.*eps(iN)*pH_tot*Kd(iN2)/X**2
    if (BB > 1.d-8) then
       pN = AA*(sqrt(1.+BB)-1.)  
    else 
       pN = 0.5*AA*BB
    endif
! C
   !  pC  = solve_q(2.*(Kd(iC2H)*pH + Kd(iC2H2)*pH**2 + Kd(iC2)), &
   !          1.+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C)+ &
   !          pO*Kd(iCO)+pN*Kd(iCN)+pH*pN*Kd(iHCN), &
   !          -epsC*pH_tot)
    X = 1.+pO*Kd(iCO)+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C) &
          +pN*Kd(iCN)+pH*pN*Kd(iHCN)
    AA = .25*X/(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))
    BB = 8.*(epsC*pH_tot)*(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))/(X*X)
    if (BB > 1.d-8) then
       pC = AA*(sqrt(1.+BB)-1.)
    else 
       pC = 0.5*AA*BB
    endif
! O
    pO  = eps(iOx)*pH_tot/(1.+pC*Kd(iCO)+pH*Kd(iOH)+pH**2*Kd(iH2O)+ &
          pO*pC*Kd(iCO2)+pSi*Kd(iSiO) &
           +pTi*Kd(iTiO)+ pO*pTi*Kd(iTiO2))
   !  pO  = solve_q(2.*(pC*Kd(iCO2)+pTi*Kd(iTiO2)), & 
   !          1.+pH*Kd(iOH)+pH**2*Kd(iH2O)+pC*Kd(iCO)+pSi*Kd(iSiO)+ &
   !          pTi*Kd(iTiO), &
   !          -eps(iOx)*pH_tot)
    pCO = Kd(iCO)*pC*pO
! Si
    a   = 3.*Kd(iSi3)
    b   = 2.*Kd(iSi2)+2.*Kd(iSi2C)*pC
    c   = 1.+Kd(iSiH)*pH+Kd(iSiH4)*pH**4+Kd(iSiO)*pO+Kd(iSiS)*pS
    d   = -eps(iSi)*pH_tot
    pSi = -d/(a*pSi**2+b*pSi+c)
! S
   !  pS  = solve_q(2.*Kd(iS2), &
   !          1.+Kd(iSiS)*pSi+Kd(iHS)*pH+Kd(iH2S)*pH**2+Kd(iTiS)*pTi, &
   !          -eps(iS)*pH_tot)
    X = 1.+pSi*Kd(iSiS)+pH*Kd(iHS)+pH**2*Kd(iH2S)+pTi*Kd(iTiS)
    AA = .25*X/Kd(iS2)
    BB = 8.*eps(iS)*pH_tot*Kd(iS2)/X**2
    if (BB > 1.d-8) then
       pS = AA*(sqrt(1.+BB)-1.)
    else
       pS = 0.5*AA*BB
    endif
! Ti
    pTi = eps(iTi)*pH_tot/(1.+Kd(iTiO)*pO+Kd(iTiO2)*pO**2+Kd(iTiS)*pS)
   !  pTi = solve_q(0._16, &
   !          1.+Kd(iTiO)*pO+Kd(iTiS)*pS+Kd(iTiO2)*pO**2, &
   !          -eps(iTi)*pH_tot)

    err = 0.
    if (pN  > 1.d-50) err = err + abs((pN-pN_old)/pN)
    if (pC  > 1.d-50) err = err + abs((pC-pC_old)/pC)
    if (pO  > 1.d-50) err = err + abs((pO-pO_old)/pO)
    if (pSi > 1.d-50) err = err + abs((pSi-pSi_old)/pSi)
    if (pS  > 1.d-50) err = err + abs((pS-pS_old)/pS)
    if (pTi > 1.d-50) err = err + abs((pTi-pTi_old)/pTi)

    nit = nit + 1
    if (nit == 2000) exit

    pN_old  = pN
    pC_old  = pC
    pO_old  = pO
    pSi_old = pSi
    pS_old  = pS
    pTi_old = pTi
 enddo
 
 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2
 pCH4  = Kd(iCH4)*pC*pH**4
 pSiO  = Kd(iSiO)*pO*pSi

 abundi(icoolH)   = pH               *patm*cst
 abundi(icoolH2)  = Kd(iH2)*pH**2    *patm*cst
 abundi(icoolHe)  = eps(ihe)*pH_tot  *patm*cst  ! pH_tot is not changing, but helium probably changes following change in mu
 abundi(icoolCO)  = Kd(iCO)*pC*pO    *patm*cst
 abundi(icoolH2O) = Kd(iH2O)*pH**2*pO*patm*cst
 abundi(icoolOH)  = Kd(iOH) *pH*pO   *patm*cst
 abundi(icoolO)   = pO               *patm*cst
 abundi(icoolSi)  = pSi              *patm*cst
 abundi(icoolC2)  = pC2              *patm*cst
 abundi(icoolC)   = pC               *patm*cst
 abundi(icoolC2H2) = pC2H2           *patm*cst
 abundi(icoolC2H) = pC2H             *patm*cst
 abundi(icoolSiO) = pSiO             *patm*cst
 abundi(icoolCH4) = pCH4             *patm*cst
 abundi(icoolS)   = pS               *patm*cst
 abundi(icoolTi)  = pTi              *patm*cst


end subroutine chemical_equilibrium_light_fixed_mu_gamma

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance with fixed mu and gamma
!  Including Broyden method
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light_fixed_mu_gamma_broyden(rho_cgs, T_in, epsC, mu, gamma, abundi)
use linalg,         only:inverse
! all quantities are in cgs
 real, intent(in)    :: rho_cgs, T_in, epsC
 real, intent(in) :: mu, gamma
 real, intent(out)   :: abundi(nabn_AGB)
 real(kind=16)   :: pC, pC2, pC2H, pC2H2
 real(kind=16)    :: pH_tot, err, a, b, c, d
 real(kind=16) :: Kd(nMolecules+1)
 real(kind=16)    :: pH, pCO, pO, pSi, pS, pTi, pN, pH2, pSiO, pCH4
 real(kind=16)    :: pN_old, pC_old, pO_old, pSi_old, pS_old, pTi_old
 real(kind=16)    :: cst
 integer :: i, nit
 real(kind=16)    :: X, AA, BB
 real             :: T
 real             :: pH_mugamma, pH_tot_mugamma, pH2_mugamma
 real(kind=16) :: denom
 logical :: first_iter
 integer :: info

 real(kind=16), dimension(6) :: x_new, x_old, f, f_old, delta, s, y, residual
 real(kind=16), dimension(6,6) :: A_b, A_inv, update, sTA
 real(kind=16) :: unity(6,6), frob_norm
 integer :: i_f, j_f

 T = max(T_in, 20.d0)

 cst = mass_per_H/(mu*mass_proton_cgs*kboltz*T)
 ! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
 pH     = solve_q(2.*Kd(iH2), 1._16, -pH_tot)

!  pH_mugamma = real(pH, kind=8)
!  pH_tot_mugamma = real(pH_tot, kind=8)
!  pH2_mugamma = real(pH2, kind=8)
!  call calc_muGamma(rho_cgs, T, mu, gamma, pH_mugamma, pH_tot_mugamma, pH2_mugamma)
!  pH = real(pH_mugamma, kind=16)
!  pH_tot = real(pH_tot_mugamma, kind=16)
!  pH2 = real(pH2_mugamma, kind=16)
 if (T > 1.d4) then
    abundi(icoolC)    = eps(iC)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolC2)   = 0.
    abundi(icoolC2H)  = 0.
    abundi(icoolC2H2) = 0.

    abundi(icoolH)   = pH  *(patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolH2)  = 1.d-50
    abundi(icoolHe)  = eps(ihe)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolCO)  = 1.d-50
    abundi(icoolH2O) = 1.d-50
    abundi(icoolOH)  = 1.d-50
    abundi(icoolO)   = eps(iOx)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolSi)  = eps(iSi)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolSiO) = 1.d-50
    abundi(icoolCH4) = 1.d-50
    abundi(icoolS)   = eps(iS)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolTi)  = eps(iTi)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    return
 endif


 Kd(iTiS) = calc_Kd_TiS(T)
 pCO      = epsC*pH_tot
 

 pN       = 0.
 pN_old   = 0.
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

 do while ((err > 1.d-6) .and. (nit<10))
! N
    X = 1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+pH*Kd(iHCN))
    AA = .25*X/Kd(iN2)
    BB = 8.*eps(iN)*pH_tot*Kd(iN2)/X**2
    if (BB > 1.d-8) then
       pN = AA*(sqrt(1.+BB)-1.)  
    else 
       pN = 0.5*AA*BB
    endif
! C
    X = 1.+pO*Kd(iCO)+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C) &
          +pN*Kd(iCN)+pH*pN*Kd(iHCN)
    AA = .25*X/(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))
    BB = 8.*(epsC*pH_tot)*(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))/(X*X)
    if (BB > 1.d-8) then
       pC = AA*(sqrt(1.+BB)-1.)
    else 
       pC = 0.5*AA*BB
    endif
! O
    pO  = eps(iOx)*pH_tot/(1.+pC*Kd(iCO)+pH*Kd(iOH)+pH**2*Kd(iH2O)+ &
          pO*pC*Kd(iCO2)+pSi*Kd(iSiO) &
           +pTi*Kd(iTiO)+ pO*pTi*Kd(iTiO2))
   !  pCO = Kd(iCO)*pC*pO
! Si
    a   = 3.*Kd(iSi3)
    b   = 2.*Kd(iSi2)+2.*Kd(iSi2C)*pC
    c   = 1.+Kd(iSiH)*pH+Kd(iSiH4)*pH**4+Kd(iSiO)*pO+Kd(iSiS)*pS
    d   = -eps(iSi)*pH_tot
    pSi = -d/(a*pSi**2+b*pSi+c)
! S
    X = 1.+pSi*Kd(iSiS)+pH*Kd(iHS)+pH**2*Kd(iH2S)+pTi*Kd(iTiS)
    AA = .25*X/Kd(iS2)
    BB = 8.*eps(iS)*pH_tot*Kd(iS2)/X**2
    if (BB > 1.d-8) then
       pS = AA*(sqrt(1.+BB)-1.)
    else
       pS = 0.5*AA*BB
    endif
! Ti
    pTi = eps(iTi)*pH_tot/(1.+Kd(iTiO)*pO+Kd(iTiO2)*pO**2+Kd(iTiS)*pS)

    err = 0.
    if (pN  > 1.d-50) err = err + abs((pN-pN_old)/pN)
    if (pC  > 1.d-50) err = err + abs((pC-pC_old)/pC)
    if (pO  > 1.d-50) err = err + abs((pO-pO_old)/pO)
    if (pSi > 1.d-50) err = err + abs((pSi-pSi_old)/pSi)
    if (pS  > 1.d-50) err = err + abs((pS-pS_old)/pS)
    if (pTi > 1.d-50) err = err + abs((pTi-pTi_old)/pTi)

    nit = nit + 1
    if (nit == 2000) exit

    pN_old  = pN
    pC_old  = pC
    pO_old  = pO
    pSi_old = pSi
    pS_old  = pS
    pTi_old = pTi
 enddo
 ! Initial guess
 x_new = (/ pN, pC, pO, pSi, pS, pTi /)
 x_old = (/ pN_old, pC_old, pO_old, pSi_old, pS_old, pTi_old /)
!  nit = 0
 first_iter = .true.

 do while (err > 1.d-6)
    if (first_iter) then
       call functions(x_new, pH, Kd, pH_tot, f)
       call jacobian(x_new, pH, Kd, A_b)
      !  A_inv = real(inverse(real(A_b, kind=8), 6), kind=16)
       call inverse_m(A_b, A_inv, 6)
      !  unity = matmul(A_b, A_inv)

       delta = matmul(-A_inv, f)
       first_iter = .false.

    else
       f_old = f
       call functions(x_new, pH, Kd, pH_tot, f)
       y = f - f_old
       s = x_new - x_old
       x_old = x_new
       

       ! residual = s - A_inv @ y
       residual = s - matmul(A_inv, y)

       ! sTA = s^T @ A_inv
       do i = 1,6
         sTA(1,i) = sum(s(:) * A_inv(:,i))
       end do

       denom = dot_product(s, matmul(A_inv, y))

       do i = 1,6
          update(i,:) = residual(i) * sTA(1,:)
       end do

       A_inv = A_inv + (1.0d0 / denom) * update

       delta = -matmul(A_inv, f)
    end if

    x_new = x_new + delta

    err = 0.0d0
    do i = 1,6
      if (x_new(i) > 1.0d-50) then
        err = err + abs((x_new(i) - x_old(i)) / x_new(i))
      end if
    end do

    nit = nit + 1
    pN  = x_new(1)
    pC  = x_new(2)
    pO  = x_new(3)
    pSi = x_new(4)
    pS  = x_new(5)
    pTi = x_new(6)
    if (nit == 2000) then
      print *, 'Not converged! err = ', err, ', T = ', T
      print *, abs((pC - pC_old)/pC), abs((pO - pO_old)/pO), &
               abs((pSi - pSi_old)/pSi), abs((pS - pS_old)/pS), &
               abs((pTi - pTi_old)/pTi)
      print *, pC, pO, pSi, pS, pTi
      exit
    end if

 enddo

 pN = x_new(1)
 pC = x_new(2)
 pO = x_new(3)
 pSi = x_new(4)
 pS = x_new(5)
 pTi = x_new(6)

 
 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2
 pCH4  = Kd(iCH4)*pC*pH**4
 pSiO  = Kd(iSiO)*pO*pSi

 abundi(icoolH)   = pH               *patm*cst
 abundi(icoolH2)  = Kd(iH2)*pH**2    *patm*cst
 abundi(icoolHe)  = eps(ihe)*pH_tot  *patm*cst  ! pH_tot is not changing, but helium probably changes following change in mu
 abundi(icoolCO)  = Kd(iCO)*pC*pO    *patm*cst
 abundi(icoolH2O) = Kd(iH2O)*pH**2*pO*patm*cst
 abundi(icoolOH)  = Kd(iOH) *pH*pO   *patm*cst
 abundi(icoolO)   = pO               *patm*cst
 abundi(icoolSi)  = pSi              *patm*cst
 abundi(icoolC2)  = pC2              *patm*cst
 abundi(icoolC)   = pC               *patm*cst
 abundi(icoolC2H2) = pC2H2           *patm*cst
 abundi(icoolC2H) = pC2H             *patm*cst
 abundi(icoolSiO) = pSiO             *patm*cst
 abundi(icoolCH4) = pCH4             *patm*cst
 abundi(icoolS)   = pS               *patm*cst
 abundi(icoolTi)  = pTi              *patm*cst
!  endif


end subroutine chemical_equilibrium_light_fixed_mu_gamma_broyden

!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance in the gas phase
!  Only including Carbon species
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_Cspecies(rho_cgs, T_in, epsC, mu, gamma, abundi)
! all quantities are in cgs
 real, intent(in)    :: rho_cgs, T_in, epsC
 real, intent(inout) :: mu, gamma
 real, intent(out)   :: abundi(nabn_AGB)
 real(kind=16)   :: pC, pC2, pC2H, pC2H2
 real(kind=16)    :: pH_tot, err, a, b, c, d
 real(kind=16) :: Kd(nMolecules+1)
 real(kind=16)    :: pH, pH2, pCH4
 real(kind=16)    :: pC_old
 real(kind=16)    :: cst
 integer :: i, nit
 real(kind=16)    :: X, AA, BB
 real             :: T
 real             :: pH_mugamma, pH_tot_mugamma, pH2_mugamma

 T = max(T_in, 20.d0)

 pH_mugamma = real(pH, kind=8)
 pH_tot_mugamma = real(pH_tot, kind=8)
 pH2_mugamma = real(pH2, kind=8)
 call calc_muGamma(rho_cgs, T, mu, gamma, pH_mugamma, pH_tot_mugamma, pH2_mugamma)
 pH = real(pH_mugamma, kind=16)
 pH_tot = real(pH_tot_mugamma, kind=16)
 pH2 = real(pH2_mugamma, kind=16)
 if (T > 1.d4) then
    abundi(icoolC)    = eps(iC)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    abundi(icoolC2)   = 0.
    abundi(icoolC2H)  = 0.
    abundi(icoolC2H2) = 0.
    abundi(icoolH)   = pH  *(patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
    return
 endif
 
 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)

! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)
 pH       = solve_q(2.*Kd(iH2), 1._16, -pH_tot)

 pC       = 0.
 pC_old   = 0.

 err      = 1.
 nit      = 0

 do while (err > 1.d-6)
! C
    X = 1.+pH**4*Kd(iCH4)
    AA = .25*X/(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))
    BB = 8.*(epsC*pH_tot)*(Kd(iC2)+pH*Kd(iC2H)+pH**2*Kd(iC2H2))/(X*X)
    if (BB > 1.d-8) then
       pC = AA*(sqrt(1.+BB)-1.)
    else 
       pC = 0.5*AA*BB
    endif

    err = 0.
    if (pC  > 1.d-50) err = err + abs((pC-pC_old)/pC)

    nit = nit + 1
    if (nit == 2000) exit

    pC_old  = pC
 enddo
 
 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2
 pCH4  = Kd(iCH4)*pC*pH**4

 cst = mass_per_H/(mu*mass_proton_cgs*kboltz*T)
 abundi(icoolH)   = pH               *patm*cst
 abundi(icoolH2)  = Kd(iH2)*pH**2    *patm*cst
 abundi(icoolC2)  = pC2              *patm*cst
 abundi(icoolC)   = pC               *patm*cst
 abundi(icoolC2H2) = pC2H2           *patm*cst
 abundi(icoolC2H) = pC2H             *patm*cst
 abundi(icoolCH4) = pCH4             *patm*cst

 ! These abundances are number densities, so in units of cm^{-3}

!  endif


end subroutine chemical_equilibrium_Cspecies

subroutine inverse_m(a_in,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none
integer, intent(in) :: n
real(kind=16), intent(in) :: a_in(n,n)
real(kind=16), intent(out) :: c(n,n)
real(kind=16) :: L(n,n), U(n,n), b(n), d(n), x(n)
real(kind=16) :: coeff
real(kind=16), dimension(n,n) :: a
integer :: i, j, k

a = a_in

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse_m

!-------------------------------------------------------
!
!  Subroutine to compute equilibrium chemistry equations
!
!-------------------------------------------------------
subroutine functions(p, pH, K, pH_tot, f)
  implicit none
  real(kind=16), intent(in)  :: p(6), pH, K(:), pH_tot
  real(kind=16), intent(out) :: f(6)

  real(kind=16) :: pN, pC, pO, pSi, pS, pTi
  pN = p(1); pC = p(2); pO = p(3)
  pSi = p(4); pS = p(5); pTi = p(6)

  f(1) = pN + pN * pH**3 * K(iNH3) + pC * pN * K(iCN) + pH * pC * pN * K(iHCN) + 2 * pN**2 * K(iN2) - eps(iN) * pH_tot
  f(2) = pC + 2 * pC**2 * K(iC2) + 2 * pC**2 * pH * K(iC2H) + 2 * pC**2 * pH**2 * K(iC2H2) + pC * pO * K(iCO) + &
         pC * pO**2 * K(iCO2) + pC * pH**4 * K(iCH4) + pSi**2 * pC * K(iSi2C) + pC * pN * K(iCN) + &
         pH * pC * pN * K(iHCN) - eps(iC) * pH_tot
  f(3) = pO + pC * pO * K(iCO) + pO * pH * K(iOH) + pH**2 * pO * K(iH2O) + pSi * pO * K(iSiO) + &
         2 * pC * pO**2 * K(iCO2) + pTi * pO * K(iTiO) + 2 * pTi * pO**2 * K(iTiO2) - eps(iOx) * pH_tot
  f(4) = pSi + 2 * pSi**2 * K(iSi2) + 3 * pSi**3 * K(iSi3) + pSi * pH * K(iSiH) + 2 * pSi**2 * pC * K(iSi2C) + &
         pSi * pH**4 * K(iSiH4) + pSi * pO * K(iSiO) + pSi * pS * K(iSiS) - eps(iSi) * pH_tot
  f(5) = pS + 2 * pS**2 * K(iS2) + pH * pS * K(iHS) + pH**2 * pS * K(iH2S) + pSi * pS * K(iSiS) + &
         pTi * pS * K(iTiS) - eps(iS) * pH_tot
  f(6) = pTi + pTi * pO * K(iTiO) + pTi * pO**2 * K(iTiO2) + pTi * pS * K(iTiS) - eps(iTi) * pH_tot
end subroutine functions

!---------------------------------------------------------------
!
!  Subroutine to compute the jacobian of the chemistry equations
!
!---------------------------------------------------------------
subroutine jacobian(p, pH, K, J)
  implicit none
  real(kind=16), intent(in) :: p(6), pH, K(:)
  real(kind=16), intent(out) :: J(6,6)

  real(kind=16) :: pN, pC, pO, pSi, pS, pTi
  pN = p(1); pC = p(2); pO = p(3)
  pSi = p(4); pS = p(5); pTi = p(6)

  ! Initialize matrix
  J = 0.0d0

  J(1,1) = 1 + pH**3 * K(iNH3) + pC * K(iCN) + pH * pC * K(iHCN) + 4 * pN * K(iN2)
  J(1,2) = pN * K(iCN) + pH * pN * K(iHCN)

  J(2,1) = pC * K(iCN) + pH * pC * K(iHCN)
  J(2,2) = 1 + 4 * pC * K(iC2) + 4 * pC * pH * K(iC2H) + 4 * pC * pH**2 * K(iC2H2) + pO * K(iCO) + &
               pO**2 * K(iCO2) + pH**4 * K(iCH4) + pSi**2 * K(iSi2C) + pN * K(iCN) + pH * pN * K(iHCN)
  J(2,3) = pC * K(iCO) + 2 * pC * pO * K(iCO2)
  J(2,4) = 2 * pSi * pC * K(iSi2C)

  J(3,2) = pO * K(iCO) + 2 * pO**2 * K(iCO2)
  J(3,3) = 1 + pC * K(iCO) + pH * K(iOH) + pH**2 * K(iH2O) + pSi * K(iSiO) + 4 * pC * pO * K(iCO2) + &
               pTi * K(iTiO) + 4 * pTi * pO * K(iTiO2)
  J(3,4) = pO * K(iSiO)
  J(3,6) = pO * K(iTiO) + 2 * pO**2 * K(iTiO2)

  J(4,2) = 2 * pSi**2 * K(iSi2C)
  J(4,3) = pSi * K(iSiO)
  J(4,4) = 1 + 4 * pSi * K(iSi2) + 9 * pSi**2 * K(iSi3) + pH * K(iSiH) + 4 * pSi * pC * K(iSi2C) + &
               pH**4 * K(iSiH4) + pO * K(iSiO) + pS * K(iSiS)
  J(4,5) = pSi * K(iSiS)

  J(5,4) = pS * K(iSiS)
  J(5,5) = 1 + 4 * pS * K(iS2) + pH * K(iHS) + pH**2 * K(iH2S) + pSi * K(iSiS) + pTi * K(iTiS)
  J(5,6) = pS * K(iTiS)

  J(6,3) = 2 * pTi * pO * K(iTiO2) + pTi * K(iTiO)
  J(6,5) = pTi * K(iTiS)
  J(6,6) = 1 + pO * K(iTiO) + pO**2 * K(iTiO2) + pS * K(iTiS)
end subroutine jacobian



!-----------------------------
!
!  solve 2nd order polynomial
!
!-----------------------------
pure real(kind=16) function solve_q(a, b, c)
  ! Inputs in quad precision
  real(kind=16), intent(in) :: a, b, c
  real(kind=16) :: delta

  if (-4.0_16 * a * c / (b**2) > epsilon(0.0_16)) then
   !   delta = max(b**2 - 4.0_16 * a * c, 0.0_16)
     delta = b**2 - 4.0_16 * a * c
     solve_q = (-b + sqrt(delta)) / (2.0_16 * a)
  else
     solve_q = -c / b
  end if

end function solve_q

!-------------------------------------------------------------------------
!
!  Compute saturation pressure of carbon clusters C_1, ..., C_5 over graphite
!
!-------------------------------------------------------------------------
pure real function psat_C(T)
! all quantities are in cgs
 real, intent(in) :: T

 real, parameter :: f = 13.8287  ! Conversion factor from atm to cgs
 real :: T2,T3,pC1!,pC2,pC3,pC4,pC5

 if (T > 1.d4) then
    Psat_C = 1.d50
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
pure real(16) function calc_Kd(coefs, T)
  ! all quantities are in cgs, computed in quad precision
  implicit none
  real(16), intent(in) :: coefs(5)
  real, intent(in) :: T
  real(16), parameter :: R = 1.987165_16
  real(16) :: G, d

  G = coefs(1)/T + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T)*T)*T
  d = -G/(R*T)
  calc_Kd = exp(d)
end function calc_Kd

pure real(16) function calc_Kd_TiS(T)
! all quantities are in cgs
 real, intent(in) :: T
 real(16), parameter :: a = 1.3316d1, b = -6.2216, c = 4.5829d-1, d = -6.4903d-2, e = 3.2788d-3
 real(16) :: theta, logKd
 theta = 5040./T
 logKd = a+(b+(c+(d+e*theta)*theta)*theta)*theta
!  calc_Kd_TiS = 10.**(-min(222.,max(logKd,-222.)))*patm
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
