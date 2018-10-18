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
 use part,           only:gr,igas
 use units,          only:set_units
 use options,        only:iexternalforce,alpha,alphau
 use timestep,       only:tmax
 use io,             only:fatal
 use prompting,      only:prompt
 use metric,         only:imetric
 use metric_tools,   only:imet_schwarzschild
 use externalforces, only:accradius1,accradius1_hard
 use inject,         only:inject_init,inject_particles,gammawind,masspart,choose_inject
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=40) :: infile
 integer :: isol
 logical :: iexist

 infile = trim(fileprefix)//'.in'

 if (.not.gr) call fatal('setup_bondiwind','This setup only works with GR on')
 if (imetric/=imet_schwarzschild) call fatal('setup_bondiwind','This setup is meant for use with the Schwarzschild metric')
 call set_units(G=1.,c=1.)

 ! General parameters
 time  = 0.
 tmax  = 360.
 polyk = 0.
 iexternalforce  = 1

 iexist = .false.
 inquire(file=trim(infile),exist=iexist)
 if (.not.iexist) then
    call choose_inject
    accradius1 = 2.5
    call prompt('Enter accretion radius of black hole ',accradius1,0.)
    accradius1_hard = accradius1
 endif

 call inject_init(setup=.true.,sol=isol)
 !-- Geodesic flow
 if (isol == 1) then
    alpha         = 0.
    alphau        = 0.
    tmax          = 90.
 endif

 massoftype(igas) = masspart
 gamma            = gammawind
 npart            = 0
 npartoftype(:)   = 0

 print "(50('*'))"
 print*,'Edit options in '//trim(infile)//' and re-run phantomsetup'
 print "(50('*'),/)"

end subroutine setpart

end module setup
