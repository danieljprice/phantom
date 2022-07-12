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
                             calc_dlnQdlnT
   use dust_formation, only: init_muGamma, set_abundances, kappa_gas, &
                             calc_kappa_bowen
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
 real :: Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17, Qtot, dlnQ_dlnT
 
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
 
 print*, ' '
 print*, ' '
 print*, '-----------------------------------------------------------------------------'
 print*, ' '
 print*, ' '
 Q1  = cool_dust_discrete_contact(T_gas, rho_gas, mu, T_dust, d2g, a, rho_grain, kappa_dust)
 print*, 'Q1   = ', Q1
!  Q2  = cool_dust_full_contact(T_gas, rho_gas, mu, T_dust, kappa_dust)
!  print*, 'Q2   = ', Q2
!  Q3  = cool_dust_radiation(T_gas, kappa_gas, T_dust, kappa_dust)
!  print*, 'Q3   = ', Q3
 Q4  = cool_coulomb(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q4   = ', Q4
 Q5  = cool_HI(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q5   = ', Q5
 Q6  = cool_H_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q6   = ', Q6
 Q7  = cool_He_ionisation(T_gas, rho_gas, mu, nH, nHe)
 print*, 'Q7   = ', Q7
 Q8  = cool_H2_rovib(T_gas, nH, nH2)
 print*, 'Q8   = ', Q8
 Q9  = cool_H2_dissociation(T_gas, rho_gas, mu, nH, nH2)
 print*, 'Q9   = ', Q9
 Q10 = cool_CO_rovib(T_gas, rho_gas, mu, nH, nH2, nCO)
 print*, 'Q10  = ', Q10
 Q11 = cool_H2O_rovib(T_gas, rho_gas, mu, nH, nH2, nH2O)
 print*, 'Q11  = ', Q11
 Q12 = cool_OH_rot(T_gas, rho_gas, mu, nOH)
 print*, 'Q12  = ', Q12
 Q13 = heat_dust_friction(rho_gas, v_drift, d2g, a, rho_grain, kappa_dust)
 print*, 'Q13  = ', Q13
 Q14 = heat_dust_photovoltaic_soft(T_gas, rho_gas, mu, nH, nHe, kappa_dust)
 print*, 'Q14  = ', Q14
!  Q15 = heat_dust_photovoltaic_hard(T_gas, nH, d2g, kappa_dust, JL)
!  print*, 'Q15  = ', Q15
 Q16 = heat_CosmicRays(nH, nH2)
 print*, 'Q16  = ', Q16
!  Q17 = heat_H2_recombination(T_gas, rho_gas, mu, nH, nH2, T_dust)
!  print*, 'Q17  = ', Q17
 
 Qtot = calc_Q(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, & 
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, &
                     JL)
 print*, 'Qtot = ', Qtot
                     
 dlnQ_dlnT = calc_dlnQdlnT(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, & 
                            T_dust, v_drift, d2g, a, rho_grain, kappa_dust, &
                            JL)
 print*, 'dlnQdlnT = ', dlnQ_dlnT 
 
 print*, ' '
 print*, ' '
 print*, '------------------- exit --------------------------------'
 print*, ' '
 print*, ' '

end subroutine do_analysis

end module analysis
