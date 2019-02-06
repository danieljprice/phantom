!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!  this module sets up a neutron star accretion disc
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
 use setdisc, only:set_disc
 use units,   only:udist, umass, utime, set_units
 use physcon, only:solarm,km
 use io,      only:master
 use externalforces, only:accradius1,mass1,iext_prdrag
 use options,        only:iexternalforce, ieos, alpha
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real(kind=8) :: udist_km
 real :: R_in, R_out
 real :: Mdisc,Mstar

 call set_units(dist=10.*km,c=1.)

 udist_km = 1.e5/udist   ! code units are cm
 gamma = 1.0
 npart = size(xyzh(1,:))
 npartoftype(1) = npart
 hfact = 1.2
 time  = 0.
 R_in  = 10.*udist_km
 R_out  = 8000.*udist_km
 alpha = 0.2

 Mstar = 1.4d0*solarm/umass !0.207
 Mdisc = Mstar*5.e-16

 print*,' 1km is ',udist_km,' in code units'
 print*,' Mstar is ', Mstar, ' in code units, or ', Mstar*umass/solarm, 'solar masses'
 print*,' Mdisc is ', Mdisc

 call set_disc(id,master=master,&
               npart   = npartoftype(1),&
               rmin    = R_in, &
               rmax    = R_out,&
               p_index = 1.5,    &
               q_index = 3./14.,   &
               HoverR  = 0.02,   &
               disc_mass = Mdisc,   &
               star_mass = Mstar,   &
               gamma   = gamma,  &
               particle_mass = massoftype(1), &
               ismooth = .true., &
               hfact = hfact, &
               xyzh = xyzh, &
               vxyzu = vxyzu, &
               polyk = polyk, &
               alpha = alpha, &
               prefix = fileprefix)
 !
 !--set default options for the input file
 !
 iexternalforce = iext_prdrag
 ieos = 3
 accradius1 = R_in
 mass1 = Mstar

 return
end subroutine setpart

end module setup
