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
!  this module does general accretion disc setups
!  Modified from an original routine by Giuseppe Lodato
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, extern_lensethirring, externalforces, io, options,
!    part, physcon, setdisc, setup_params, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!
! This subroutine is a utility for setting up discs
!
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setdisc,        only:set_disc
 use part,           only:igas
 use io,             only:master
 use externalforces, only:accradius1,accradius1_hard
 use options,        only:iexternalforce,alpha,alphau,iexternalforce
 use units,          only:set_units,umass
 use physcon,        only:solarm,pi
#ifdef GR
 use metric,         only:mass1,a
#else
 use externalforces, only:mass1,iext_einsteinprec
 use extern_lensethirring, only:blackhole_spin
#endif
 use prompting,      only:prompt
 use timestep,       only:tmax,dtmax
 use eos,            only:ieos
 use kernel,         only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: r_in,r_out,honr,theta,mhole,mdisc,spin

 time            = 0.
 alphau          = 0.0
 npartoftype(:)  = 0
 iexternalforce  = 1
 hfact           = hfact_default
#ifdef GR
 ieos  = 4
 gamma = 5./3.
#else
 ieos  = 2
 gamma = 1.
 iexternalforce = iext_einsteinprec
#endif

 tmax  = 1000.
 dtmax = 10.

 mhole = 1.e6*solarm
 call set_units(G=1.,c=1.,mass=mhole)

 mdisc           = 10.*solarm/umass
 accradius1      = 4.
 accradius1_hard = accradius1

!
! Set problem parameters (defaults)
!
 npart   = 1e5
 r_in    = 40.
 r_out   = 160.
 spin    = 0.
 theta   = 0.      ! inclination angle (degrees)
 honr    = 0.02

 call prompt('Enter number of particles ',npart)
 call prompt('Enter H on R ',honr)
 call prompt('Enter spin of black hole ',spin,-1.,1.)
 call prompt('Enter inner radius of disc ',r_in,2.)
 call prompt('Enter outer radius of disc ',r_out,r_in)
 call prompt('Enter inclination angle (degrees) ',theta,0.,90.)

 npartoftype(1)  = npart

#ifdef GR
 a = spin
#else
 blackhole_spin = spin
#endif

 theta = theta/180. * pi ! convert to radians

 call set_disc(id,master,&
               npart         = npart,                &
               rmin          = r_in,                 &
               rmax          = r_out,                &
               p_index       = -1.0,                 &
               q_index       = 0.0,                  &
               HoverR        = honr,                 &
               gamma         = gamma,                &
               hfact         = hfact,                &
               xyzh          = xyzh,                 &
               vxyzu         = vxyzu,                &
               polyk         = polyk,                &
               particle_mass = massoftype(igas),     &
               ! star_mass     = 1.0,                &
               disc_mass     = mdisc,                &
               inclination   = theta,                &
               ! bh_spin       = spin,               &
               ! alpha         = alpha,              &
               prefix        = fileprefix)
#ifdef GR
 polyk = vxyzu(4,1)
#else
 polyk = (5./3. -1.)*7.2e-6
#endif

 return
end subroutine setpart

end module setup
