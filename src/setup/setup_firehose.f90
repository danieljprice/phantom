!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! setup for GR "firehose" problem
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  empty setup for driven simulation
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use units,   only:set_units
 use physcon, only:au,solarm
 use part,    only:igas
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)

 call set_units(mass=1e6*solarm,G=1.d0,c=1.d0)
 time = 0.
 polyk = 0.
 gamma = 5./3.

 npart = 0
 npartoftype(:) = 0
 massoftype = 0.
 massoftype(igas) = 1.e-6

end subroutine setpart

end module setup
