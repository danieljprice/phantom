!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  Reads in sphNG data without B-fields, spits out phantom data with B-fields
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, setup_params
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setup_params, only: ihavesetupB
 use part,         only: hfact,mhd
! use timestep,      only: dtmax
 !traced dtmax to "timestep" but can't find "timestep"?
 !also can't find "dim"?

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)


 ! Ensure that sphNG data correctly conforms to phantom standards
 ! Check variables are properly defined:
 ! ntypes
 ! npartoftype
 ! isink (?)
 ! nparttotal
 ! = number of particles = npart = sum(npartoftype)?
 ! system time (check dtmax as well)
 ! - from 'set_default_options' - options which comes from 'timestep' ?

 ! Define hfact
 hfact = 1.2


 print*,'sphNG data reformatted for phantom write'


 ! Now add B field
 if (mhd) then
    ihavesetupB = .false.
 endif

end subroutine modify_dump

end module moddump

