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
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------

!----------------------------------------------------------------
!+
!  Parameters related to timestepping
!+
!---------------------------------------------------------------
module timestep
 implicit none
 real    :: tmax,dtmax
 real    :: C_cour,C_force,C_cool
 integer :: nmax,nout
 integer :: nsteps
 real, parameter :: bignumber = 1.e29

 real    :: dt, dtcourant, dtforce, dtextforce, dterr, dtdiff, time
 logical :: restartonshortest

 ! When rho_max > rho_dtthresh, then decrease dtmax by dtmax_rat;
 ! this can currently be done only once per simulation.
 ! The logicals will be automatically adjusted, so do not modify.
 real,    public :: rho_dtthresh
 real,    public :: rho_dtthresh_cgs  = 0.0
 integer, public :: dtmax_rat0        = 1
 integer, public :: dtmax_rat
 logical, public :: mod_dtmax         = .true.
 logical, public :: mod_dtmax_now     = .false.
 logical, public :: mod_dtmax_in_step = .false.

end module timestep

!----------------------------------------------------------------
!+
!  Parameters related the initial conditions
!  These will be calculated upon first initialisation, and will be
!  used for subsequent restarts
!+
!---------------------------------------------------------------
module initial_params
use dim, only:ndusttypes
 implicit none
 real,    public :: get_conserv = 1.0 ! to track when we have initial values for conservation laws
 real,    public :: etot_in,angtot_in,totmom_in,mdust_in(ndusttypes)

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
