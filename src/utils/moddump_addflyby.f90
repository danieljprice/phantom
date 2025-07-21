!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump to add a flyby (sink particle on a parabolic orbit)
!
! :References: None
!
! :Owner: Cristiano Longarini
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, part, physcon, prompting, setflyby, units,
!   vectorutils
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

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
 real                   :: star_m, mperturber, accrperturber, dma, n0, big_omega, incl
 real                   :: mtot, eccentricity, reducedmass, pf
 real                   :: x0, y0, z0, vx0, vy0, vz0, r0, m0, period, dtmax
 real                   :: xp(3), vp(3), rot_axis(3)

 star_m = xyzmh_ptmass(4,1) ! star mass
 mperturber = 1. ! mass of the perturber
 accrperturber = 1. ! accretion radius of the perturber
 dma = 100. ! distance of minimum approach
 n0 = 1000. ! initial distance of the perturber
 big_omega = 0. ! position angle of the ascending node
 incl = 0. ! inclination
 dtmax = 0.05

! user's questions
 call prompt('Enter the mass of the perturber: ', mperturber, 0.)
 call prompt('Enter accretion radius of the perturber: ', accrperturber, 0.)
 call prompt('Enter distance of minimum approach of the perturber: ', dma, 0.)
 call prompt('Enter initial distance of the perturber in units of minimum approach: ', n0, 0.)
 call prompt('Enter position angle of the ascending node of the perturber: ', big_omega, 0.)
 call prompt('Enter inclination of the perturber: ', incl, 0.)
 call prompt('Enter time between dumps as fraction of flyby time:', dtmax, 0.)


 nptmass = nptmass + 1
 mtot = star_m + mperturber ! total mass
 eccentricity = 1. ! eccentricity = 1 -> parabolic orbit
 reducedmass = star_m*mperturber/mtot ! reduced mass
 pf = 2*dma ! focal parameter -- dma=pf/2

 ! define m0 = -x0/dma such that r0 = n0*dma
 !  companion starts at negative x and y
 !  positive root of 1/8*m**4 + m**2 + 2(1-n0**2) = 0
 !  for n0 > 1
 ! copied from set_flyby.f90
 m0 = 2*sqrt(n0-1.0)

 ! perturber initial position
 x0 = -m0*dma
 y0 = dma*(1.0-(x0/pf)**2)
 z0 = 0.0
 xp = (/x0,y0,z0/)

 ! perturber initial velocity
 r0  = sqrt(x0**2+y0**2+z0**2)
 vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
 vy0 = -(x0/r0)*sqrt(mtot/pf)
 vz0 = 0.0
 vp  = (/vx0,vy0,vz0/)

 ! assigning position, velocity and accretion radius to the sink particle
 xyzmh_ptmass(1:3,nptmass)   = xp
 xyzmh_ptmass(4,nptmass)     = mperturber
 xyzmh_ptmass(5,nptmass)     = accrperturber
 xyzmh_ptmass(6,nptmass)     = 0.
 vxyz_ptmass(:,nptmass)      = vp

 incl = pi-incl*pi/180.
 big_omega = big_omega*pi/180.
 rot_axis = (/sin(big_omega),-cos(big_omega),0./)
 call rotatevec(xyzmh_ptmass(1:3,nptmass),rot_axis,incl)
 call rotatevec(vxyz_ptmass(1:3,nptmass), rot_axis,incl)

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 period = get_T_flyby(star_m,mperturber,dma,n0) * dtmax ! computing time between dumps in code units

 write(*,*) 'Time between dumps in code units to put in *.in file', period

end subroutine modify_dump

end module moddump

