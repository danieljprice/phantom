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
 use options,        only:iexternalforce,alpha,alphau,icooling
 use units,          only:set_units
 use physcon,        only:solarm,pi
 use metric,         only:mass1
 use metric,         only:a
 use prompting,      only:prompt
 use timestep,       only:tmax,dtmax
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: R_in,R_out,HonR,theta

 call set_units(mass=1.*solarm,c=1.)

 tmax = 1000.
 dtmax = 10.

 !
 !  Set problem parameters
 !
 !--disc inner and outer radius

 a       = 0.
 R_in    = 4.233
 R_out   = 100.
 theta   = 0.          ! inclination angle (degrees)

 call prompt('Enter spin of black hole ',a,-1.,1.)
 call prompt('Enter inner radius of disc ',r_in,2.*mass1)
 call prompt('Enter outer radius of disc ',r_out,r_in)
 call prompt('Enter inclination angle (degrees) ',theta,0.,90.)
 call prompt('Cooling ',icooling)

 theta = theta/180. * pi ! convert to radians

! npart   = size(xyzh(1,:))
 npart = 1e5
 npartoftype(:) = 0
 npartoftype(1) = npart
 gamma   = 5./3.          !Unless otherwise specified, 1.0
 time    = 0.

 ! alpha   = 0.!1.0
 alphau  = 0.1
 HonR = 0.01

 iexternalforce = 1
 accradius1 = 5.*mass1 !r_in!+2.!4.
 accradius1_hard = accradius1 - (0.5*(accradius1-2.*mass1))

 call set_disc(id,master=master,&
                nparttot  = npart,                &
                npart     = npart,                &
                rmin      = R_in,                 &
                rmax      = R_out,                &
                p_index   = -1.0,                 &
                q_index   = 0.5,                  &
                HoverR    = HonR,                 &
                disc_Q    = 168.,                 &
                star_mass = 1.0,                  &
                gamma     = gamma,                &
                particle_mass = massoftype(igas), &
                hfact     = 1.0,                  &
                xyzh      = xyzh,                 &
                vxyzu     = vxyzu,                &
                polyk     = polyk,                &
                twist     = .false.,              &
                ismooth   = .true.,               &
                inclination = theta,              &
                warp_smoothl = 0.,                &
                bh_spin = 0.,                     &
                prefix = fileprefix)

 return
end subroutine setpart

end module setup
