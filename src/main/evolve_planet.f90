!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module evolveplanet
!
! This module evolve the embedded planet
!   calculating the accretion and wind rate according
!   to the pressure at the Bondi radius of the planet
!
! :References: Rogers et al. 2024 for the boil-off model
!
! :Owner: Cristiano Longarini
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units
!
 implicit none
 public :: evolve_planet

 private

contains

subroutine evolve_planet(pbondi,rbondi,mdotacc,mdotwind)
! this routine should call James' code to calculate mdotacc and mdotwind for boil-off
 use units, only:umass,utime
 use physcon, only:jupiterm,years
 real, intent(in) :: pbondi,rbondi
 real, intent(out) :: mdotacc,mdotwind

 mdotacc = 0.
 mdotwind = 1.e-2*jupiterm/umass/(years/utime)

end subroutine evolve_planet

end module evolveplanet
