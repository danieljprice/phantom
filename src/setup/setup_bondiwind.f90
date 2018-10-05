!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does bondi wind (accretion in reverse)
!
!  REFERENCES: None
!
!  OWNER: David Liptai
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
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,           only: xyzmh_ptmass, vxyz_ptmass, nptmass, gr
 use units,          only: set_units
 use options,        only: iexternalforce
 use timestep,       only: tmax,nmax
 use io,             only: fatal
 use prompting,      only: prompt
 use metric,         only: mass1, imetric
 use metric_tools,   only: imet_schwarzschild
 use externalforces, only: accradius1,accradius1_hard
 use inject,         only: wind_init,inject_particles,gammawind,masspart
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 if (.not.gr) call fatal('setup_bondiwind','This setup only works with GR on')
 if (imetric/=imet_schwarzschild) call fatal('setup_bondiwind','This setup is meant for use with the Schwarzschild metric')
 call set_units(G=1.,c=1.)

 ! General parameters
 time  = 0.
 tmax  = 1000.
 polyk = 0.
 gamma = gammawind
 iexternalforce  = 1
 mass1           = 1.
 accradius1      = 2.*mass1
 accradius1_hard = accradius1

 call wind_init(.true.)

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0

 massoftype = masspart

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.


 ! --- Put sink with no mass (needed to get around bug where phantom doesn't want to read an empty dumpfile)
 nptmass = 1
 xyzmh_ptmass(1,1) = 0.
 xyzmh_ptmass(2,1) = 0.
 xyzmh_ptmass(3,1) = 0.
 xyzmh_ptmass(4,1) = 0.
 xyzmh_ptmass(5,1) = 3.

 ! call inject_particles(time,0.,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 ! nmax = 0

end subroutine setpart

end module setup
