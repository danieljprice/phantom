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
 use krome_getphys
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
 use eos,       only:ieos,get_local_temperature,get_local_u_internal!equationofstate

 integer, intent(in) :: npart
 real, intent(in) :: dt,xyzh(:,:)
 real, intent(inout) :: vxyzu(:,:)
 integer :: i
 real :: T_local, dt_cgs, hi, rhoi, rho_cgs

 dt_cgs = dt*utime
!$omp parallel do schedule(runtime) &
!$omp default(none) &
!$omp shared(ieos,dt,dt_cgs, npart, unit_density) &
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
     rho_cgs = rhoi*unit_density
     call get_local_temperature(ieos,xyzh(1,i),xyzh(2,i),xyzh(3,i),rhoi,mu_chem(i),vxyzu(4,i),gamma_chem(i),T_local)
     ! normalise abudances and balance charge conservation with e-
     call krome_consistent_x(species_abund(:,i))
     ! evolve the chemistry and update the abundances
     call evolve_chemistry(species_abund(:,i),rho_cgs,T_local,dt_cgs)
     ! update the gas temperature array for the dumpfiles
     krometemperature(i) = T_local
     ! update the particle's mean molecular weight
     mu_chem(i) =  krome_get_mu_x(species_abund(:,i))
     ! update the particle's adiabatic index
     gamma_chem(i) = krome_get_gamma_x(species_abund(:,i),T_local)
     ! calculate the cooling contribution, needed for force.F90
     kromecool(i) = (get_local_u_internal(gamma_chem(i), mu_chem(i), T_local)-vxyzu(4,i))/dt
   endif
 enddo
!$omp end parallel do

end subroutine update_krome

subroutine evolve_chemistry(species, dens, temp, time)

 real, intent(inout) :: species(:), temp
 real, intent(in)    :: dens, time
 real, allocatable   :: test_species1(:), test_species2(:)
 real                :: test_temp1, test_temp2, test_dens, test_time
 real                :: dudt, dt_cool
 integer             :: i, N
 
 test_species1 = species
 test_species2 = species
 test_dens     = dens
 test_temp1    = temp
 test_temp2    = temp
 test_time     = time
 
 call krome(test_species2,test_dens,test_temp2,test_time)
 
 ! Calculate cooling timescale
 dudt = abs(test_temp2 - test_temp1)/test_time
 dt_cool = abs(test_temp1/dudt)
 
 ! Substepping if dt_cool < input timestep
 !!!! explicit addition of krome cooling in force.F90 to determine hydro timestep not needed anymore
 !!!! skipping the contribution to fxyz4 and updating the final particle energy still to be implemented
 if (dt_cool < test_time) then 
    N = ceiling(test_time/dt_cool)
    do i = 1,N       
       call krome(test_species1,test_dens,test_temp1,test_time/N)
    enddo
    species = test_species1
    temp    = test_temp1
 else
    species = test_species2
    temp    = test_temp2
 endif

end subroutine evolve_chemistry

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
