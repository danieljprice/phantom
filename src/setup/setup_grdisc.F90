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
 use options,        only:iexternalforce,alpha,alphau,icooling,iexternalforce
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
 real    :: R_in,R_out,HonR,theta,mhole,mdisc,spin

#ifdef GR
 ieos  = 4
 gamma = 5./3.
#else
 ieos  = 2
 gamma = 1.
 iexternalforce = iext_einsteinprec
#endif

 mhole = 200.*solarm!1.e6*solarm
 call set_units(G=1.,c=1.,mass=mhole)
 mdisc = 1.e-6 !10.*solarm/umass
 hfact = hfact_default

 tmax = 1000.
 dtmax = 10.

 !
 !  Set problem parameters
 !
 !--disc inner and outer radius

 spin    = 0.
 R_in    = 40.
 R_out   = 160.
 theta   = 0.          ! inclination angle (degrees)
 HonR    = 0.02
 npart   = 1e5

 call prompt('Enter number of particles ',npart)
 call prompt('Enter H on R ',HonR)
 call prompt('Enter spin of black hole ',spin,-1.,1.)
 call prompt('Enter inner radius of disc ',r_in,2.*mass1)
 call prompt('Enter outer radius of disc ',r_out,r_in)
 call prompt('Enter inclination angle (degrees) ',theta,0.,90.)
 call prompt('Cooling ',icooling)

#ifdef GR
 a = spin
#else
 blackhole_spin = spin
#endif

 theta = theta/180. * pi ! convert to radians

 npartoftype(:) = 0
 npartoftype(1) = npart
 time    = 0.

 alphau  = 0.0

 iexternalforce  = 1
 accradius1      = 2.5 !4.*mass1
 accradius1_hard = accradius1 !- (0.5*(accradius1-2.*mass1))

 call set_disc(id,master,&
               npart         = npart,                &
               rmin          = R_in,                 &
               rmax          = R_out,                &
               p_index       = -1.0,                 &
               q_index       = 0.0,                  &
               HoverR        = HonR,                 &
               gamma         = gamma,                &
               hfact         = hfact_default,        &
               xyzh          = xyzh,                 &
               vxyzu         = vxyzu,                &
               polyk         = polyk,                &
               particle_mass = massoftype(igas),     &
               ! star_mass     = 1.0,                 &
               disc_mass     = mdisc,                &
               inclination   = theta,                &
               ! bh_spin       = spin,                 &
               ! alpha         = alpha,               &
               prefix        = fileprefix)
#ifdef GR
 polyk = vxyzu(4,1)
#endif

 return
end subroutine setpart

end module setup
