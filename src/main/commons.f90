!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup_params
!
!  DESCRIPTION:
! This file contains all of the main "common blocks"
!  for the code
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim
!+
!--------------------------------------------------------------------------

!----------------------------------------------------------------
!+
!  Parameters related the initial conditions
!  These will be calculated upon first initialisation, and will be
!  used for subsequent restarts
!+
!---------------------------------------------------------------
module initial_params
 use dim, only:maxdusttypes
 implicit none
 real,    public :: get_conserv = 1.0 ! to track when we have initial values for conservation laws
 real,    public :: etot_in,angtot_in,totmom_in,mdust_in(maxdusttypes)

end module initial_params

!----------------------------------------------------------------
!+
!  Parameters related to initial setup
!+
!---------------------------------------------------------------
module setup_params
 implicit none
 real            :: rmax,rhozero
 logical         :: ihavesetupB = .false.
 integer(kind=8) :: npart_total

end module setup_params
