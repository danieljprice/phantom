!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump to add a sink particle
!
! :References: None
!
! :Owner: Taj Jankovič & Aleksej Jurca
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, physcon, prompting, setflyby, units,
!   vectorutils
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:nptmass,xyzmh_ptmass,vxyz_ptmass!,igas,ihacc,ihsoft
 use prompting,         only:prompt
 use units,             only:umass,utime,udist,print_units
 use physcon,           only:au,solarm,pi,years
 use centreofmass,      only:reset_centreofmass,get_centreofmass
 use vectorutils,       only:rotatevec
 use setflyby,          only:get_T_flyby
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: msink, rsink, z0, incl, vz0
 real                   :: xp(3), vp(3), rot_axis(3)
 real, parameter        :: disk_H = 3.0 

 ! sink particle properties
 msink = solarm/umass ! mass of the perturber - 1Msun
 !msink = 0.925 ! for Sun-like polytrope with Rsink = 0.8Rstar
 z0 = 5 ! initial distance of the perturber
 incl = 90./180*3.14192 ! inclination
 vz0   = 1 ! 0 for core
 rsink = 1 !0.8 for core

 print *,'-----Sink particle mass is:', msink
 print *,'-----Sink particle radius is:', rsink
 print *,'-----Sink particle distance alog z-direction is:', z0
 print *,'-----Sink particle velocity is:', vz0

 nptmass = nptmass + 1

 ! --- Compute initial position ---
 if (incl < 0.5*pi) then
     xp = (/ -z0 * (cos(incl)/sin(incl)), 0.0, z0 /) 
 else
     xp = (/ 0.0, 0.0, z0 /)
 endif

 print *,'-----Sink particle position is:', xp

 ! perturber initial velocity
 vp  = (/vz0*cos(incl),0.0,-vz0*sin(incl)/)
 print *,'-----Sink particle velocity is:', vp

 ! assigning position, velocity and accretion radius to the sink particle
 xyzmh_ptmass(1:3,nptmass)   = xp
 xyzmh_ptmass(4,nptmass)     = msink
 xyzmh_ptmass(5,nptmass)     = rsink
 xyzmh_ptmass(6,nptmass)     = rsink ! softening length
 vxyz_ptmass(:,nptmass)      = vp

 ! test - 1 particle
 !xyzh(1:3,1) = (/0.0,0.0,7./)

 
! call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return

end subroutine modify_dump

end module moddump

