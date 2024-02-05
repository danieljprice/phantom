!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setflyby
!
! This module is contains utilities for setting up flyby.
! Our conventions for angles are the same as in Xiang-Gruess (2016).
! Eccentricity is set to unity, i.e. for a parabolic orbit:
!
!  - minimum_approach : *distance of minimum approach (pericentre)*
!  - initial_dist     : *the initial separation distance (in units of minimum_approach)*
!  - posang_ascnode   : *angle counter-clockwise (East) from y-axis (North)*
!  - inclination      : *angle of rotation of orbital plane around axis defined by the position angle (for posang_ascnode=0 this is a roll angle)*
!
! :References:
!   - Xiang-Gruess (2016), MNRAS 455, 3086-3100
!   - Cuello et al. (2019), MNRAS 483, 4114-4139
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters: None
!
! :Dependencies: physcon, vectorutils
!
 use physcon, only:pi
 implicit none
 public :: set_flyby,get_T_flyby

 integer, parameter :: &
   ierr_m1   = 1, &
   ierr_m2   = 2, &
   ierr_ecc  = 3, &
   ierr_peri = 4

 private

contains

!----------------------------------------------------------------
!+
!  setup for a flyby
!+
!----------------------------------------------------------------
subroutine set_flyby(m1,m2,minimum_approach,initial_dist, &
                     accretion_radius1,accretion_radius2,xyzmh_ptmass, &
                     vxyz_ptmass,nptmass,ierr,posang_ascnode,inclination,verbose)
 use vectorutils, only:rotatevec
 real,    intent(in)    :: m1,m2
 real,    intent(in)    :: minimum_approach,initial_dist
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 integer, intent(out)   :: ierr
 real,    intent(in), optional :: posang_ascnode,inclination
 logical, intent(in), optional :: verbose

 integer :: i1,i2,i
 real    :: mtot,n0
 real    :: xp(3),vp(3),reducedmass
 real    :: dma,ecc,eccentricity,m0,pf,r0,x0,y0,z0,vx0,vy0,vz0
 real    :: big_omega,incl,rot_axis(3)
 logical :: do_verbose

 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose

 !--inclination about zero position angle corresponds to a roll angle
 big_omega = 0.
 incl = 0.
 if (present(posang_ascnode)) big_omega = posang_ascnode
 if (present(inclination)) incl = inclination

 i1 = nptmass + 1
 i2 = nptmass + 2
 nptmass = nptmass + 2

 !--total mass
 mtot = m1 + m2

 !--eccentricity is set to 1, i.e. parabolic orbit
 eccentricity = 1.0

 reducedmass = m1*m2/mtot

 if (do_verbose) then
    print "(/,2x,a)",'---------- flyby parameters ----------- '
    print "(8(2x,a,g12.3,/),2x,a,g12.3)", &
        'primary mass            :',m1, &
        'secondary mass          :',m2, &
        'mass ratio              :',m2/m1, &
        'reduced mass            :',reducedmass, &
        'dist of min. app.       :',minimum_approach, &
        'pos. angle ascen. node  :',big_omega, &
        'inclination             :',incl, &
        'eccentricity            :',eccentricity
 endif

 !--check for bad parameter choices
 ierr = 0
 if (m1 <= 0.) then
    print "(1x,a)",'ERROR: set_flyby: primary mass <= 0'
    ierr = ierr_m1
 endif
 if (m2 < 0.) then
    print "(1x,a)",'ERROR: set_flyby: secondary mass <= 0'
    ierr = ierr_m2
 endif
 if (eccentricity > 1. .or. eccentricity < 0.) then
    print "(1x,a)",'ERROR: set_flyby: eccentricity must be between 0 and 1'
    ierr = ierr_ecc
 endif
 if (minimum_approach <= 0.) then
    print "(1x,a)",'ERROR: set_flyby: distance of min approach <= 0'
    ierr = ierr_peri
 endif
 if (ierr /= 0) return

 dma = minimum_approach
 ecc = eccentricity
 n0  = initial_dist

 !--focal parameter dma = pf/2
 pf = 2*dma

 !--define m0 = -x0/dma such that r0 = n0*dma
 !  companion starts at negative x and y
 !  positive root of 1/8*m**4 + m**2 + 2(1-n0**2) = 0
 !  for n0 > 1
 m0 = 2*sqrt(n0-1.0)

 !--perturber initial position
 x0 = -m0*dma
 y0 = dma*(1.0-(x0/pf)**2)
 z0 = 0.0
 xp = (/x0,y0,z0/)

 !--perturber initial velocity
 r0  = sqrt(x0**2+y0**2+z0**2)
 vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
 vy0 = -(x0/r0)*sqrt(mtot/pf)
 vz0 = 0.0
 vp  = (/vx0,vy0,vz0/)

 !--positions and accretion radii
 xyzmh_ptmass(:,i1:i2) = 0.
 xyzmh_ptmass(1:3,i1) = 0.
 xyzmh_ptmass(1:3,i2) = xp
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(5,i1) = accretion_radius1
 xyzmh_ptmass(5,i2) = accretion_radius2
 xyzmh_ptmass(6,i1) = 0.
 xyzmh_ptmass(6,i2) = 0.

 !--velocities
 vxyz_ptmass(:,i1) = 0.
 vxyz_ptmass(:,i2) = vp

 !--incline orbit about ascending node
 ! if incl 0 = prograde orbit
 ! if incl 180 = retrograde orbit
 ! Convention: clock-wise rotation in the zx-plane
 incl = pi-incl*pi/180.
 big_omega = big_omega*pi/180.
 do i=i1,i2
    rot_axis = (/sin(big_omega),-cos(big_omega),0./)
    call rotatevec(xyzmh_ptmass(1:3,i),rot_axis,incl)
    call rotatevec(vxyz_ptmass(1:3,i), rot_axis,incl)
 enddo

end subroutine set_flyby

!-------------------------------------------------------------
!
! Function to determine the time from initial separation,
! n0 [dma], to the equivalent separation past periastron given
! the distance of minimum approach (dma) using Barker's equation
!
!-------------------------------------------------------------
function get_T_flyby(m1,m2,dma,n0) result(T)
 real, intent(in) :: m1,m2,dma,n0
 real :: T,nu,xi,yi,Di,Df,p,mu,G

 !--semi-latus rectum
 p = 2*dma

 !--initial position
 xi = -2*sqrt(n0-1.0)*dma
 yi = dma*(1.0-(xi/p)**2)

 !--gravitational parameter
 G = 1.0 ! we assume code units where G=1
 mu = G*(m1+m2)

 !--true anomaly
 nu = pi - atan(abs(xi/yi))

 !--Barker's equation
 Di = tan(-nu/2)
 Df = tan(nu/2)
 T = 1./2.*sqrt(p**3/mu) * (Df + 1./3.*Df**3 - Di - 1./3.*Di**3)

end function get_T_flyby

end module setflyby
