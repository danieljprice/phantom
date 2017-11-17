!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setflyby
!
!  DESCRIPTION:
!   This module is contains utilities for setting up flyby
!   Our conventions for angles are the same as in Xiang-Gruess 2016
!   Eccentricity is set = 1, i.e. parabolic orbits
!   dma = distance of minimum approach
!   r0 = sqrt(y0*y0+x0+x0) = n0 * dma
!   inclination of the orbit --> roll angle
!
!  REFERENCES:
!   Eggleton (1983) ApJ 268, 368-369 (ref:eggleton83)
!   Xiang-Gruess (2016) MNRAS, Volume 455, Issue 3, p.3086-3100
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon, setbinary
!+
!--------------------------------------------------------------------------
module setflyby
 use physcon, only:pi
 implicit none
 public :: set_flyby,get_T_flyby

 private

contains

!----------------------------------------------------------------
!+
!  setup for a flyby
!+
!----------------------------------------------------------------
subroutine set_flyby(mprimary,massratio,dma,n0, &
                     accretion_radius1,accretion_radius2, &
                     xyzmh_ptmass,vxyz_ptmass,nptmass,verbose, &
                     roll)

 use part,      only:ihacc,ihsoft
 use setbinary, only:Rochelobe_estimate
 real,    intent(in)    :: mprimary,massratio
 real,    intent(in)    :: dma,n0
 real,    intent(in)    :: accretion_radius1,accretion_radius2
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass
 real,    intent(in)    :: roll
 logical, intent(in),  optional :: verbose

 integer :: i1,i2,i
 real    :: m1,m2,mtot,Rochelobe,Rochelobe2
 real    :: xp(3),vp(3),cosi,sini,xangle,reducedmass
 real    :: a,ecc,eccentricity,m0,pf,c,r0,vx0,vy0,vz0,x0,y0,z0
 logical :: do_verbose

 do_verbose = .true.
 if (present(verbose)) do_verbose = verbose

 i1 = nptmass + 1
 i2 = nptmass + 2
 nptmass = nptmass + 2

 ! mass of the central star (m1) and the perturber (m2)
 m1 = mprimary
 m2 = mprimary*massratio
 mtot = m1 + m2

 ! eccentricity is set to 1, i.e. parabolic orbit
 eccentricity = 1.0

 Rochelobe = Rochelobe_estimate(m1,m2,dma)
 Rochelobe2 = Rochelobe_estimate(m2,m1,dma)
 reducedmass = m1*m2/mtot

 if (do_verbose) then
    print "(/,2x,a)",'---------- flyby parameters ----------- '
    print "(8(2x,a,g12.3,/),2x,a,g12.3)", &
        'primary mass     :',mprimary, &
        'secondary mass   :',massratio*mprimary, &
        'mass ratio       :',massratio, &
        'reduced mass     :',reducedmass, &
        'dist of min app  :',dma, &
        'roll angle       :',roll, &
        'eccentricity     :',eccentricity
 endif

 if (accretion_radius1 >  Rochelobe2) then
    print "(/,a,/)",'*** WARNING: accretion radius of primary > Roche lobe'
 endif
 if (accretion_radius2 >  Rochelobe) then
    print "(/,a,/)",'*** WARNING: accretion radius of primary > Roche lobe'
 endif
 !
 !--check for stupid parameter choices
 !
 if (mprimary <= 0.) stop 'ERROR: primary mass <= 0'
 if (massratio < 0.) stop 'ERROR: mass ratio < 0'
 if (dma <= 0.)      stop 'ERROR: distance of min approach <= 0'

 a = dma
 ecc = eccentricity
 c = sqrt(2.*mtot*a)
 pf = (c**2)/mtot !focal parameter a = pf/2

 ! We look the value of m that verifies x0 = - m0*a and r0= n0*a
 ! The companion should start at negative x and y coordinates (z=0, roll rotation at the end)
 ! Positive root of: 1/8*m**4 + m**2 + 2(1-n**2) = 0
 m0 = 2*sqrt(n0-1.0) ! this number should be positive, so n0 > 1 (warning in the setup?)

 !--perturber initial position
 x0=-m0*a
 y0=a*(1.0-(x0/pf)**2)
 z0 = 0.0
 r0=sqrt(x0**2+y0**2+z0**2)
 xp(:) = (/x0,y0,z0/)

 !--perturber initial velocity
 vx0 = (1.+ (y0/r0))*(mtot/c)
 vy0 = -(x0/r0)*(mtot/c)
 vz0 = 0.0
 vp(:) = (/vx0,vy0,vz0/)

 !
 !--positions and accretion radii
 !
 xyzmh_ptmass(:,i1:i2) = 0.
 xyzmh_ptmass(1:3,i1) = 0.
 xyzmh_ptmass(1:3,i2) = xp
 xyzmh_ptmass(4,i1) = m1
 xyzmh_ptmass(4,i2) = m2
 xyzmh_ptmass(ihacc,i1) = accretion_radius1
 xyzmh_ptmass(ihacc,i2) = accretion_radius2
 xyzmh_ptmass(ihsoft,i1) = 0.
 xyzmh_ptmass(ihsoft,i2) = 0.
 !
 !--velocities
 !
 vxyz_ptmass(:,i1) = 0.
 vxyz_ptmass(:,i2) = vp
 !
 !--roll angle rotation
 !
 if (roll /= 0.) then
    xangle = roll*pi/180.
    cosi = cos(xangle)
    sini = sin(xangle)
    do i=i1,i2
       call rotate(xyzmh_ptmass(1:3,i),cosi,sini)
       call rotate(vxyz_ptmass(1:3,i),cosi,sini)
    enddo
 endif

end subroutine set_flyby

!--------------------------------------------
! Rotation around the y-axis (roll angle)
!--------------------------------------------
pure subroutine rotate(xyz,cosi,sini)
 real, intent(inout) :: xyz(3)
 real, intent(in)    :: cosi,sini
 real :: xi,yi,zi

 xi = xyz(1)
 yi = xyz(2)
 zi = xyz(3)
 xyz(1) =  xi*cosi + zi*sini  !xi
 xyz(2) =  yi                 !yi*cosi - zi*sini
 xyz(3) =  -xi*sini + zi*cosi !yi*sini + zi*cosi

end subroutine rotate

!-------------------------------------------------------------
! Function to determine the semi-major axis given the period
!-------------------------------------------------------------
function get_T_flyby(m1,m2,a,n0) result(T)
 real, intent(in) :: m1,m2,a,n0
 real :: T,m

 m = 2*sqrt(n0-1.0)
 T = ( 4*a**2 * (8./3. - m + m**3/12) ) / sqrt(2*(m1+m2)*a)

end function get_T_flyby

end module setflyby
