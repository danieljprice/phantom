!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Reads in sphNG data without B-fields, spits out phantom data with B-fields
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: kernel, part, setup_params
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use setup_params, only:ihavesetupB
 use part,         only:hfact,mhd
 use kernel,       only:hfact_default
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
 hfact = hfact_default

 print*,'sphNG data reformatted for phantom write'


 ! Now add B field
 if (mhd) then
    ihavesetupB = .false.
 endif

end subroutine modify_dump

end module moddump
