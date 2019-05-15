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
module krome_phantom_coupling

 use krome_user
 use krome_main
 use part, only: species_abund_label,mu_chem,gamma_chem,krometemperature

 implicit none

 public :: initialise_krome

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
 mu_chem(:)           = 2.12444   ! for composition below
 gamma_chem(:)        = 1.66667
 krometemperature(:)  = 3000

 ! Initial chemical abundance value for AGB surface
 He_init = 3.11e-1 ! mass fraction
 C_init = 2.63e-3  ! mass fraction
 N_init = 1.52e-3  ! mass fraction
 O_init = 9.60e-3  ! mass fraction

 H_init = 1.0 - He_init - C_init - N_init - O_init

end subroutine initialise_krome

end module krome_phantom_coupling
