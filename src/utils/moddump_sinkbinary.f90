!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  default moddump routine: does not make any modifications
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,      only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas
 use setbinary, only:set_binary
 use units,     only:umass,udist
 use physcon,   only:solarm,solarr,pi
 use prompting, only:prompt
 use centreofmass, only:reset_centreofmass
 use timestep, only:dtmax,tmax
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: m1,m2,a,e,racc,period
 integer :: i

 ! find current mass from existing particles
 m2 = npart*massoftype(igas)
 print*,' Mass of star from existing file in Msun = ',m2*umass/solarm
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! set up a sink particle binary
 !
 m1 = 1.4
 call prompt('enter mass of point mass in Msun ',m1,0.)
 print*,' code unit of distance in Rsun = ',udist/solarr
 
 a = 2500.
 call prompt('enter semi-major axis in code units',a,0.)
 
 racc = a/200.
 call prompt('enter accretion radius of point mass',racc,0.01*a,0.5*a)
 print*,' accretion radius in solar radii = ',racc*udist/solarr

 e = 0.
 call prompt('enter eccentricity ',e,0.)

 !
 ! add sink particle binary
 !
 nptmass = 0
 call set_binary(m2,m1/m2,a,e,racc,racc,xyzmh_ptmass,vxyz_ptmass,nptmass)
 !
 ! delete one of the sink particles and replace it with our polytrope
 !
 nptmass = nptmass - 1
 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,2)
    vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,2)
 enddo
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 period = 2.*pi*sqrt(a**3/(m1 + m2))
 print*,' orbital period = ',period
 tmax = 1000.*period
 dtmax = 0.1*period

 return
end subroutine modify_dump

end module moddump

