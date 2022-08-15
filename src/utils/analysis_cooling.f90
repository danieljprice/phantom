!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis

   use physcon,        only: mass_proton_cgs, kboltz, atomic_mass_unit
   use cooling ,       only: cool_dust_discrete_contact, &
                             cool_dust_full_contact, &
                             cool_dust_radiation, &
                             cool_coulomb, &
                             cool_HI, &
                             cool_H_ionisation, &
                             cool_He_ionisation, &
                             cool_H2_rovib, &
                             cool_H2_dissociation, &
                             cool_CO_rovib, &
                             cool_H2O_rovib, &
                             cool_OH_rot, &
                             heat_dust_friction, &
                             heat_dust_photovoltaic_soft, &
                             heat_dust_photovoltaic_hard, &
                             heat_CosmicRays, &
                             heat_H2_recombination, &
                             calc_Q, &
                             calc_dlnQdlnT, &
                             print_cooling_rates
   use dust_formation, only: init_muGamma, set_abundances, kappa_gas, &
                             calc_kappa_bowen, chemical_equilibrium_light, mass_per_H
   use dim,            only:nElements

   implicit none

   character(len=20), parameter, public :: analysistype = 'cooling'
   public :: do_analysis
   private
   real, parameter :: patm = 1.013250d6
   real    :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
   real    :: eps(nElements) = [1.d0, 1.04d-1, 0.0,  6.d-4, 2.52d-4, 1.17d-4, 3.58d-5, 1.85d-5, 3.24d-5, 8.6d-8]


contains

real function MPH(eps, Aw)

  real, dimension(nElements), intent(inout) :: eps, Aw
  real :: wind_CO_ratio

  wind_CO_ratio = 2.0
  eps(3)        = eps(4) * wind_CO_ratio
  MPH           = atomic_mass_unit*dot_product(Aw,eps)

end function MPH

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:)
 real(kind=8),     intent(in) :: particlemass,time

 real :: T_gas, rho_gas, mu, gamma, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real :: pH, pH2
 real :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real :: JL
 real :: n_gas

 T_gas      = 1500.
 rho_gas    = 1.d-15

 call set_abundances
 call init_muGamma(rho_gas, T_gas, mu, gamma, pH, pH2)
 nH         = pH  *(patm*MPH(eps, Aw))/(mu*mass_proton_cgs*kboltz*T_gas)
 nH2        = pH2 *(patm*MPH(eps, Aw))/(mu*mass_proton_cgs*kboltz*T_gas)

 n_gas      = rho_gas/(mu*mass_proton_cgs)
 nHe        = 1.d-1*n_gas
 nCO        = 1.d-4*n_gas
 nH2O       = 5.d-5*n_gas
 nOH        = 1.d-7*n_gas

 kappa_gas  = 2.d-4

 T_dust     = 1000.
 v_drift    = 1.d6
 d2g        = 1./200.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_dust = calc_kappa_bowen(T_dust)

 JL         = 2.5d-12     ! Value taken from Barstow et al. 1997

 call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

 !c build rho-T grid of cooling rates
 call generate_grid()


end subroutine do_analysis

subroutine generate_grid

 real :: logtmin,logtmax,logT,dlogt,T,crate,nH_tot,rho_cgs
 real :: pC, pC2, pC2H, pC2H2, mu, gamma, T_dust, d2g, v_drift
 real :: nH, nH2, nHe, nCO, nH2O, nOH, a, rho_grain, kappa_g, n_gas, rho_gas, kappa_dust, JL
 integer :: i,iunit
 integer, parameter :: nt = 400, iC=3

 logtmax = log10(1.e8)
 logtmin = log10(1.e0)

 call set_abundances()

 open(newunit=iunit,file='new_cooltable.txt',status='replace')
 write(iunit,"(a)") '#   T   \Lambda_E(T) erg s^{-1} cm^3   \Lambda erg s^{-1} cm^{-3}'
 dlogt = (logtmax - logtmin)/real(nt)
 d2g        = 1.d-2
 v_drift    = 0.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_g    = 2.d-4
 JL         = 2.5d-12     ! Value taken from Barstow et al. 1997
 rho_cgs    = 1.d-15

 do i=1,nt
    logT   = logtmin + (i-1)*dlogt
    T      = 10**logT
    T_dust = T
    kappa_dust = calc_kappa_bowen(T)
    nH_tot     = rho_cgs/mass_per_H
    n_gas      = rho_cgs/(mu*mass_proton_cgs)
    call chemical_equilibrium_light(rho_cgs, T, eps(iC), pC, pC2, pC2H, pC2H2, mu, gamma, nH, nH2, nHe, nCO, nH2O, nOH)
    crate = calc_Q(T, rho_cgs, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_g, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    !ndens = (rho_cgs/mass_proton_cgs)*5.d0/7.d0
    !print *,rho_cgs, T, mu, gamma, nH, nH2, nHe, nCO, nH2O, nOH, crate
    write(iunit,*) t,crate/nH_tot**2,crate
 enddo
 close(iunit)
end subroutine generate_grid


end module analysis
