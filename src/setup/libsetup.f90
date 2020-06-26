!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: libsetup
!
!  DESCRIPTION:
!   This module contains the public API for setup library routines
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  DEPENDENCIES: stretchmap
!+
!--------------------------------------------------------------------------
module libsetup
 use stretchmap, only:set_density_profile
 use unifdis,    only:set_unifdis, get_ny_nz_closepacked
 public

end module libsetup
