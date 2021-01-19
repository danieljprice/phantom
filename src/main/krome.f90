!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module krome_interface
!
! This module contains all the necessary subroutines to establish
!   the coupling between phantom and KROME
!
! :References: None
!
! :Owner: Lionel
!
! :Runtime parameters: None
!
! :Dependencies: eos, krome_main, krome_user, part, units
!

 implicit none

 public :: initialise_krome,update_krome,write_KromeSetupFile

 private
 real  :: cosmic_ray_rate
 real  :: H_init, He_init, C_init, N_init, O_init
 real  :: S_init, Fe_init, Si_init, Mg_init
 real  :: Na_init, P_init, F_init


contains
!----------------------------------------------------------------
!+
!  short initialisation routine that initialises krome,
!  the arrays of chemical species labels, and other
!  necessary additional variables.
!+
!----------------------------------------------------------------
subroutine initialise_krome()

 use krome_main, only:krome_init
 use krome_user, only:krome_idx_He,krome_idx_C,krome_idx_N,krome_idx_O,krome_idx_H,&
       krome_set_user_crflux,krome_get_names,krome_get_mu_x,krome_get_gamma_x,&
       krome_idx_S,krome_idx_Fe,krome_idx_Si,krome_idx_Mg,krome_idx_Na,&
       krome_idx_P,krome_idx_F
 use part,       only:abundance,abundance_label,mu_chem,gamma_chem,T_chem
 real :: wind_temperature

 print *, ""
 print *, "==================================================="
 print *, "=                                                 ="
 print *, "=             INITIALISING KROME                  ="
 print *, "=                                                 ="
 print *, "==================================================="
 call krome_init()
 print *, ""
 print *, "========================================================="
 print *, "=            KROME INITIALISATION SUCCESSFUL            ="
 print *, "========================================================="
 print *, ""

 cosmic_ray_rate = 1.36e-17 ! in s^-1
 call krome_set_user_crflux(cosmic_ray_rate)

 abundance_label(:) = krome_get_names()

 ! Initial chemical abundance value for AGB surface
 He_init = 3.11e-1 ! mass fraction
 C_init  = 2.63e-3 ! mass fraction
 N_init  = 1.52e-3 ! mass fraction
 O_init  = 9.60e-3 ! mass fraction

 S_init  = 3.97e-4 ! mass fraction
 Fe_init = 1.17e-3 ! mass fraction
 Si_init = 6.54e-4 ! mass fraction
 Mg_init = 5.16e-4

 Na_init = 3.38e-5
 P_init  = 8.17e-6
 F_init  = 4.06e-7

 H_init = 1.0 - He_init - C_init - N_init - O_init - S_init - Fe_init - &
          Si_init - Mg_init - Na_init - P_init - F_init

 abundance(krome_idx_He,:) = He_init
 abundance(krome_idx_C,:)  = C_init
 abundance(krome_idx_N,:)  = N_init
 abundance(krome_idx_O,:)  = O_init
 abundance(krome_idx_S,:)  = S_init
 abundance(krome_idx_Fe,:) = Fe_init
 abundance(krome_idx_Si,:) = Si_init
 abundance(krome_idx_Mg,:) = Mg_init
 abundance(krome_idx_Na,:) = Na_init
 abundance(krome_idx_P,:)  = P_init
 abundance(krome_idx_F,:)  = F_init
 abundance(krome_idx_H,:)  = H_init

 !set initial wind temperature to star's effective temperature
 mu_chem(:)    = krome_get_mu_x(abundance(:,1))
 gamma_chem(:) = krome_get_gamma_x(abundance(:,1),wind_temperature)
 T_chem(:)     = wind_temperature

end subroutine initialise_krome

subroutine update_krome(dt,xyzh,u,rho,xchem,gamma_chem,mu_chem,T_chem)

 use krome_main, only: krome
 use krome_user,    only:krome_consistent_x,krome_get_mu_x,krome_get_gamma_x
 use units,         only:unit_density,utime
 use eos,           only:ieos,get_local_temperature,get_local_u_internal!equationofstate

 real, intent(in)    :: dt,xyzh(4),rho
 real, intent(inout) :: u,gamma_chem,mu_chem,xchem(:)
 real, intent(out)   :: T_chem
 real :: T_local, dt_cgs, rho_cgs

 dt_cgs = dt*utime
 rho_cgs = rho*unit_density
 call get_local_temperature(ieos,xyzh(1),xyzh(2),xyzh(3),rho,mu_chem,u,gamma_chem,T_local)
 T_local=max(T_local,20.0d0)
! normalise abudances and balance charge conservation with e-
 call krome_consistent_x(xchem)
! evolve the chemistry and update the abundances
 call krome(xchem,rho_cgs,T_local,dt_cgs)
! update the particle's mean molecular weight
 mu_chem =  krome_get_mu_x(xchem)
! update the particle's adiabatic index
 gamma_chem = krome_get_gamma_x(xchem,T_local)
! update the particle's temperature
 T_chem = T_local
! get the new internal energy
 u = get_local_u_internal(gamma_chem,mu_chem,T_local)

end subroutine update_krome

!----------------------------------------------------------------
!+
!  write Krome parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_KromeSetupFile

 integer, parameter           :: iunit = 21

 print "(a)",' writing krome setup options in krome.setup'
 open(unit=iunit,file='krome.setup',status='replace',form='formatted')
 write (iunit,'("-n=networks/react_AGB_full_noNucl")')
 write (iunit,'("#-compact")')
 write (iunit,'("-cooling=ATOMIC,CHEM,H2,CIE,Z,CI,CII,OI,OII,CO,OH,H2O,HCN")')
 write (iunit,'("-heating=CHEM,CR")')
 write (iunit,'("-H2opacity=RIPAMONTI")')
 write (iunit,'("-gamma=EXACT")')
 write (iunit,'("-noSinkCheck")')
 write (iunit,'("-noRecCheck")')
 write (iunit,'("-noTlimits")')
 write (iunit,'("-useX")')
 write (iunit,'("-conserveLin")')
 write (iunit,'("-useTabs")')
 close(iunit)

end subroutine write_KromeSetupFile

! subroutine get_local_temperature(eos_type,xi,yi,zi,rhoi,gmwi,intenerg,gammai,local_temperature)
!  use dim, only:maxvxyzu
!  integer,      intent(in)    :: eos_type
!  real,         intent(in)    :: xi,yi,zi,rhoi,gmwi,gammai
!  real,         intent(inout) :: intenerg
!  real,         intent(out)   :: local_temperature
!  real :: spsoundi,ponrhoi

!  if (maxvxyzu==4) then
!     call equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,eni=intenerg,gamma_local=gammai)
!  else
!     print *, "CHEMISTRY PROBLEM: ISOTHERMAL SETUP USED, INTERNAL ENERGY NOT STORED"
!  endif
!  local_temperature = temperature_coef*gmwi*ponrhoi

! end subroutine get_local_temperature
end module krome_interface
