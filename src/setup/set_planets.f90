!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setplanets
!
!  DESCRIPTION:
! this module is a setup utility for adding planets
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part, physcon
!+
!--------------------------------------------------------------------------
module setplanets
 implicit none
 public :: set_planets

 private

contains

!----------------------------------------------------------------
!+
!  setup for planets
!+
!----------------------------------------------------------------
subroutine set_planets(xyzmh_ptmass,vxyz_ptmass,nptmass)
 use part, only:ihacc,ihsoft
 use physcon,only:pi
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: nptmass

 integer,parameter :: nplanets = 1
 real,parameter :: G=1.0
 integer,dimension(nplanets) :: isink
 real,dimension(nplanets) :: r,phi,inc,mp,phase

 integer :: nptmass_old,i,j,mysink
 real :: dr,dphi
 real :: msink,rp,phip,vphip,vrp,tiltp,twistp,rin,rout
 real,dimension(3) :: xsink,vsink,xp,vp

 nptmass_old = nptmass

 isink(:) = 1
 rin = 0.15
 rout = 0.3
! dr = (rout-rin)/real(nplanets-1)
! dphi = 2*pi/real(nplanets-1)
 dr=0.0
 dphi=0.0
 do i=1,nplanets
    mp(i) = 0.0001
    r(i) = rin + real(i-1)*dr
    phi(i) = real(i-1)*dphi
    inc(i) = 60.0*pi/180.0
    phase(i) = 0.0*pi/180 ! phase of inclination
 enddo

 do i=1,nplanets
    j = nptmass_old+i

! set central sink information
    mysink = isink(i)
    if (mysink > 0) then
       msink = xyzmh_ptmass(4,mysink)
       xsink(1:3) = xyzmh_ptmass(1:3,mysink)
       vsink(1:3) = vxyz_ptmass(1:3,mysink)
    else
       msink = xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2)
       xsink(1:3) = (xyzmh_ptmass(4,1)*xyzmh_ptmass(1:3,1) + xyzmh_ptmass(4,2)*xyzmh_ptmass(1:3,2))/msink
       vsink(1:3) = (xyzmh_ptmass(4,1)*vxyz_ptmass(1:3,1) + xyzmh_ptmass(4,2)*vxyz_ptmass(1:3,2))/msink
    endif

! set masses
    xyzmh_ptmass(4,j) = mp(i)

! Set positions
    rp = r(i)
    phip = phi(i)
    xp(1) = rp*cos(phip)
    xp(2) = rp*sin(phip)
    xp(3) = 0.0

! Set velocities
    vphip = sqrt(G*msink/rp)
    vrp   = 0.0
    vp(1)  = -vphip*sin(phip) + vrp*cos(phip)
    vp(2)  =  vphip*cos(phip) + vrp*sin(phip)
    vp(3)  =  0.0

! rotate to get inclination
    tiltp = inc(i)
    twistp = phase(i)

    call rotate_y(xp,tiltp)
    call rotate_y(vp,tiltp)

    call rotate_z(xp,twistp)
    call rotate_z(vp,twistp)

! Move to relevant sink particle
    xyzmh_ptmass(1:3,j) = xp(1:3) + xsink(1:3)
    vxyz_ptmass(1:3,j) = vp(1:3) + vsink(1:3)

! Set planet accretion radius and (used) softening radius
    xyzmh_ptmass(ihacc,j) = 0.005
    xyzmh_ptmass(ihsoft,j) = 0.001

! increment nptmass
    nptmass = nptmass+1
 enddo

end subroutine set_planets
!--------------------------------------------------------------------------
! Rotate a 3D vector around y-axis, counter-clockwise if axis points at you
!--------------------------------------------------------------------------
subroutine rotate_y(vec,angle)
 implicit none

 real,intent(in) :: angle
 real,dimension(3), intent(inout) :: vec

 real :: c,s
 real,dimension(3) :: dummy

 dummy(:) = vec(:)

 c = cos(angle)
 s = sin(angle)

 vec(1) = c*dummy(1) + s*dummy(3)
 vec(2) = dummy(2)
 vec(3) = -s*dummy(1) + c*dummy(3)

end subroutine rotate_y
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Rotate a 3D vector around z-axis, counter-clockwise if axis points at you
!--------------------------------------------------------------------------
subroutine rotate_z(vec,angle)
 implicit none

 real,intent(in) :: angle
 real,dimension(3), intent(inout) :: vec

 real :: c,s
 real,dimension(3) :: dummy

 dummy(:) = vec(:)

 c = cos(angle)
 s = sin(angle)

 vec(1) = c*dummy(1) - s*dummy(2)
 vec(2) = s*dummy(1) + c*dummy(2)
 vec(3) = dummy(3)

end subroutine rotate_z
!--------------------------------------------------------------------------
end module setplanets
