!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!    Setup of "Chinese coin" orbital dynamics problem
!    from Chin & Chen (2005)
!
!  REFERENCES: Chin & Chen (2005), Celest. Mech. Dyn. Astron. 91, 301
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: extern_binary, externalforces, io, options, part, physcon,
!    timestep, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use units,          only:set_units
 use physcon,        only:solarm,au,pi
 use options,        only:iexternalforce
 use externalforces, only:iext_binary
 use extern_binary,  only:binarymassr
 use io,             only:master
 use timestep,       only:dtmax
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer            :: ierr
 logical            :: iexist
 real               :: m1
!
!--units
!
 call set_units(mass=solarm,dist=au,G=1.d0)
!
!--general parameters
!
 time = 0.
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

 m1    = 1.
 nptmass = 1
 xyzmh_ptmass = 0.
 xyzmh_ptmass(2,1) = 0.0580752367
 xyzmh_ptmass(4,1) = m1
 xyzmh_ptmass(ihacc,1) = 0.1

 vxyz_ptmass = 0.
 vxyz_ptmass(1,1) = 0.489765446

 iexternalforce = iext_binary
 binarymassr = 0.5
 dtmax = 0.1*(9.*pi)

end subroutine setpart

end module setup
