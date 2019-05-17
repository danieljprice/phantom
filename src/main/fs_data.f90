!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: fs_data
!
!  DESCRIPTION:
!  Atomic data (transition probabilities, energies) required
!  for computing the fine structure cooling rates.
!
!  Written by S. Glover, AMNH, 2004-2005
!  Adapted for use in Phantom by D. Price, MoCA, 2011
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon
!+
!--------------------------------------------------------------------------
module fs_data
 use physcon, only:kb=>kboltz
 implicit none

 real(kind=8), parameter :: oxa10 = 8.865d-5
 real(kind=8), parameter :: oxa20 = 1.275d-10
 real(kind=8), parameter :: oxa21 = 1.772d-5
 real(kind=8), parameter :: oxe10 = 2.2771d2 * kb
 real(kind=8), parameter :: oxe21 = 9.886d1  * kb
 real(kind=8), parameter :: oxe20 = oxe10 + oxe21

 real(kind=8), parameter :: cIa10 = 7.932d-8
 real(kind=8), parameter :: cIa20 = 2.054d-14
 real(kind=8), parameter :: cIa21 = 2.654d-7
 real(kind=8), parameter :: cIe10 = 2.360d1 * kb
 real(kind=8), parameter :: cIe21 = 3.884d1 * kb
 real(kind=8), parameter :: cIe20 = cIe10 + cIe21

 real(kind=8), parameter :: cIIa10 = 2.291d-6
 real(kind=8), parameter :: cIIe10 = 9.125d1 * kb

 real(kind=8), parameter :: siIa10 = 8.4d-6
 real(kind=8), parameter :: siIa20 = 2.4d-10
 real(kind=8), parameter :: siIa21 = 4.2d-5
 real(kind=8), parameter :: siIe10 = 1.1d2 * kb
 real(kind=8), parameter :: siIe21 = 2.1d2 * kb
 real(kind=8), parameter :: siIe20 = siIe10 + siIe21

 real(kind=8), parameter :: siIIa10 = 2.17d-4
 real(kind=8), parameter :: siIIe10 = 4.1224d2 * kb

end module fs_data
