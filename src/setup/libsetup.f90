!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module libsetup
!
! This module contains the public API for setup library routines
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: stretchmap, unifdis
!
 use stretchmap, only:set_density_profile
 use unifdis,    only:set_unifdis, get_ny_nz_closepacked
 public

end module libsetup
