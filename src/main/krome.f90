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
                     species_abund,mu_chem,gamma_chem,krometemperature
 use units,     only:unit_density,utime,udist
 use eos,       only:ieos,get_local_temperature,get_local_u_internal!equationofstate

 integer, intent(in) :: npart
 real, intent(in) :: dt
 real, intent(inout) :: vxyzu(:,:),xyzh(:,:)
 integer :: i
 real :: T_local, dt_cgs, hi, rhoi, rho_cgs
 real :: particle_position_sq

 dt_cgs = dt*utime
!$omp parallel do schedule(runtime) &
!$omp default(none) &
!$omp shared(ieos,dt,dt_cgs, npart, unit_density, udist) &
!$omp shared(species_abund) &
!$omp shared(krometemperature) &
!$omp shared(mu_chem) &
!$omp shared(gamma_chem) &
!$omp shared(vxyzu,xyzh, massoftype) &
!$omp private(i,T_local, hi, rhoi, rho_cgs, particle_position_sq)
 
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
     ! when modelling stellar wind: remove partiles beyond r_max=100*AU => resolution too low
     particle_position_sq = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
     if (particle_position_sq > (1.496e15/udist)**2) xyzh(4,i) = -abs(xyzh(4,i))
   endif
 enddo
!$omp end parallel do

end subroutine update_krome

subroutine evolve_chemistry(species, dens, temp, time)

 real, intent(inout) :: species(:), temp
 real, intent(in)    :: dens, time
 real, allocatable   :: dupl_species1(:), dupl_species2(:)
 real                :: dupl_temp1, dupl_temp2, dupl_dens, dupl_time
 real                :: dudt, dt_cool
 integer             :: i, N
 
 ! Duplicate input arrays 
 dupl_species1 = species
 dupl_species2 = species
 dupl_dens     = dens
 dupl_temp1    = temp
 dupl_temp2    = temp
 dupl_time     = time
 
 call krome(dupl_species2,dupl_dens,dupl_temp2,dupl_time)
 
 ! Calculate cooling timescale
 dudt = abs(dupl_temp2 - dupl_temp1)/dupl_time
 dt_cool = abs(dupl_temp1/dudt)
 
 ! Substepping if dt_cool < input timestep
 !!!! explicit addition of krome cooling in force.F90 to determine hydro timestep not needed anymore
 !!!! skipping the contribution to fxyz4 and updating the final particle energy still to be implemented
 if (dt_cool < dupl_time) then 
    N = ceiling(dupl_time/dt_cool)
    do i = 1,N       
       call krome(dupl_species1,dupl_dens,dupl_temp1,dupl_time/N)
    enddo
    species = dupl_species1
    temp    = dupl_temp1
 else
    species = dupl_species2
    temp    = dupl_temp2
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
