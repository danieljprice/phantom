!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: kernel, part, physcon, setbinary, units
!
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use setbinary, only:set_binary
 use units,     only:umass,udist,set_units
 use physcon,   only:au,solarm,solarr
 use kernel,    only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: m1,m2,a,hacc1,hacc2,ecc
 integer :: ierr

 call set_units(dist=au,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 hfact = hfact_default
 polyk = 0.
 gamma = 1.
!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 massoftype = 1d-9

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 m1    = 6./5.*solarm/umass
 m2    = 1./5.*m1
 a     = 180*solarr/udist
 !a     = 2.9 ! 70 AU
 ecc   = 0.
 hacc1  = 0.017 ! 1.7 AU
 hacc2  = 0.017
 call set_binary(m1,m2,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)

end subroutine setpart

end module setup
