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
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: inject, part, physcon, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for 3D Bondi Hoyle Lyttleton simulation
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,       only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use physcon,    only:pi,au,solarm
 use units,      only:udist,umass,utime,set_units
 use inject,     only:BHL_init,BHL_r_star,BHL_m_star,BHL_pmass
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real            :: m,hacc
 real, parameter :: irrational_number_close_to_one = 3./pi

 call set_units(dist=1.,time=1.,G=1.)

!
!--general parameters
!
 time  = 0.
 polyk = 0.
 gamma = 5./3.

 call BHL_init(.true.)
 m     = BHL_m_star / umass ! Depends on parameters in the input file
 hacc  = BHL_r_star / udist ! Depends on parameters in the input file

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 massoftype = BHL_pmass/umass

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 nptmass = 1
 xyzmh_ptmass(4,nptmass) = m
 xyzmh_ptmass(ihacc,nptmass) = hacc
 xyzmh_ptmass(ihsoft,nptmass) = 0.0

 print *, "udist = ", udist, "; umass = ", umass, "; utime = ", utime
end subroutine setpart

end module setup
