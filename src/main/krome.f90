!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: krome_phantom_coupling
!
!  DESCRIPTION:
!   This module contains all the necessary subroutines to establish
!   the coupling between phantom and KROME
!
!  REFERENCES: None
!
!  OWNER: Lionel
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: krome_main, krome_user, part
!+
!--------------------------------------------------------------------------
module krome_interface

 use krome_user
 use krome_main
 use part, only: species_abund_label,mu_chem,gamma_chem,krometemperature

 implicit none

 public :: initialise_krome,update_krome

 real  :: cosmic_ray_rate
 real  :: H_init, He_init, C_init, N_init, O_init

contains
!----------------------------------------------------------------
!+
!  short initialisation routine that initialises krome,
!  the arrays of chemical species labels, and other
!  necessary additional variables.
!+
!----------------------------------------------------------------
subroutine initialise_krome()

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

 species_abund_label(:) = krome_get_names()
 mu_chem(:)            = 2.12444   ! for composition below
 gamma_chem(:)         = 1.66667
 krometemperature(:)   = 0

 ! Initial chemical abundance value for AGB surface
 He_init = 3.11e-1 ! mass fraction
 C_init  = 2.63e-3 ! mass fraction
 N_init  = 1.52e-3 ! mass fraction
 O_init  = 9.60e-3 ! mass fraction

 H_init = 1.0 - He_init - C_init - N_init - O_init

end subroutine initialise_krome

subroutine update_krome(dt,npart,xyzh,vxyzu)

 use part,      only:rhoh,massoftype,isdead_or_accreted,igas,&
                     species_abund,mu_chem,gamma_chem,krometemperature,kromecool
 use units,     only:unit_density, utime
 use eos,       only:ieos,get_local_temperature!equationofstate

 integer, intent(in) :: npart
 real, intent(in) :: dt,xyzh(:,:)
 real, intent(inout) :: vxyzu(:,:)
 integer :: i
 real :: T_local, dt_cgs, hi, rhoi, rho_cgs

 dt_cgs = dt*utime
!$omp parallel do schedule(runtime) &
!$omp default(none) &
!$omp shared(ieos,dt_cgs, npart, unit_density) &
!$omp shared(species_abund) &
!$omp shared(krometemperature) &
!$omp shared(kromecool) &
!$omp shared(mu_chem) &
!$omp shared(gamma_chem) &
!$omp shared(vxyzu,xyzh, massoftype) &
!$omp private(i,T_local, hi, rhoi, rho_cgs)
 
 do i = 1,npart
   hi = xyzh(4,i)
   if (.not.isdead_or_accreted(hi))  then
     rhoi = rhoh(hi,massoftype(igas))
     rho_cgs = rhoi**unit_density
     call get_local_temperature(ieos,xyzh(1,i),xyzh(2,i),xyzh(3,i),rhoi,mu_chem(i),vxyzu(4,i),gamma_chem(i),T_local)
     ! evolve the chemistry and update the abundances
     call krome(species_abund(:,i),rho_cgs,T_local,dt_cgs)
     ! calculate the cooling contribution, needed for force.F90
     kromecool(i) = krome_get_cooling(krome_x2n(species_abund(:,i),rho_cgs),T_local)
     ! update the gas temperature array for the dumpfiles
     krometemperature(i) = T_local
     ! update the particle's mean molecular weight
     mu_chem(i) =  krome_get_mu(krome_x2n(species_abund(:,i),rho_cgs))
     ! update the particle's adiabatic index
     gamma_chem(i) = krome_get_gamma_x(species_abund(:,i),T_local)
   endif
 enddo
!$omp end parallel do

end subroutine update_krome

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
