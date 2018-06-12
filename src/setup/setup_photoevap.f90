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
!  this module sets up photoevaporation problem
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: externalforces, io, options, physcon, setdisc, units
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
 use units,          only:umass,set_units
 use setdisc,        only:set_disc
 use physcon,        only:pi,au,solarm
 use io,             only:master
 use externalforces, only:accradius1,iext_star
 use options,        only:iexternalforce
 integer, intent(in)              :: id
 integer, intent(out)             :: npart
 integer, intent(out)             :: npartoftype(:)
 real,    intent(out)             :: xyzh(:,:)
 real,    intent(out)             :: polyk,gamma,hfact
 real,    intent(out)             :: vxyzu(:,:)
 real,    intent(out)             :: massoftype(:)
 real,    intent(inout)           :: time
 character (len=20) , intent (in) :: fileprefix
 real :: R_in,R_out
 real :: xinc
 real :: discmass

 !
 !  Set problem parameters
 !

 !--disc inner and outer radius
 R_in     = 2.
 R_out    = 10.

 npart    = size(xyzh(1,:))

 call set_units(dist=au,mass=solarm,G=1.)

 npartoftype(:) = 0
 npartoftype(1) = npart
 gamma    = 5./3.
 time     = 0.
 xinc     = 0.0*(pi/180.0) ! Must be in radians
 discmass = 5.d22/umass

 call set_disc(id,master        = master,            &
                  npart         = npartoftype(1),    &
                  rmin          = R_in,              &
                  rmax          = R_out,             &
                  p_index       = 1.5,               &
                  q_index       = 0.75,              &
                  HoverR        = 0.05,              &
                  disc_mass     = discmass,          &
                  star_mass     = 1.0,               &
                  gamma         = gamma,             &
                  particle_mass = massoftype(1),     &
                  hfact         = hfact,             &
                  xyzh          = xyzh,              &
                  vxyzu         = vxyzu,             &
                  polyk         = polyk,             &
                  inclination   = xinc,              &
                  prefix        = fileprefix)

 !
 !--set default options for the input file
 !
 iexternalforce = iext_star
 accradius1 = R_in

 return
end subroutine setpart

end module setup
