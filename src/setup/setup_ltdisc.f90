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
!  DEPENDENCIES: extern_lensethirring, externalforces, io, options,
!    physcon, setdisc, units
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
 use setdisc,   only:set_disc
 use physcon,   only:pi,solarm
 use io,             only:master
 use externalforces, only:accradius1
 use options,        only:iexternalforce,alpha
 use extern_lensethirring, only:blackhole_spin
 use units,          only:set_units
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real    :: R_in,R_out,xinc

 call set_units(mass=1.d9*solarm,c=1.)
 !
 !  Set problem parameters
 !
 !--disc inner and outer radius
 R_in    = 50.0
 R_out   = 100.0
 npart   = size(xyzh(1,:))
 npartoftype(:) = 0
 npartoftype(1) = npart
 gamma   = 1.0
 time    = 0.
 hfact   = 1.2
 xinc    = 0.0*(pi/180.0) ! Must be in radians
 alpha   = 0.1
 blackhole_spin=1.0

 call set_disc(id,master=master,&
                npart   = npartoftype(1),&
                rmin    = R_in, &
                rmax    = R_out,&
                p_index = 1.5,    &
                q_index = 0.75,   &
                HoverR  = 0.01,   &
                disc_mass = 0.001,   &
                star_mass = 1.0,    &
                gamma   = gamma,  &
                particle_mass = massoftype(1), &
                hfact = hfact, &
                xyzh = xyzh, &
                vxyzu = vxyzu, &
                polyk = polyk, &
                inclination = xinc, &
                twist = .false., &
                alpha = alpha, &
                bh_spin = blackhole_spin, &
                prefix = fileprefix)
!
!--set default options for the input file
!
 iexternalforce = 9
 accradius1 = 1.0
!
 return
end subroutine setpart

end module setup
