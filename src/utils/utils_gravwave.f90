!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module gravwaveutils
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units
!
 implicit none

 public :: calculate_strain

 private

contains

! This subroutine computes the gravitational wave strain at a distance of 1Mpc
subroutine calculate_strain(hx,hp,hxx,hpp,xyzh,vxyz,axyz,pmass,npart)
 use units,   only:umass,udist,utime
 use physcon, only:gg,c
 real, intent(out)   :: hx,hp,hxx,hpp
 real, intent(in)    :: xyzh(:,:), vxyz(:,:), axyz(:,:), pmass
 integer, intent(in) :: npart
 real       :: q(6),ddq(6),x,y,z,vx,vy,vz,ax,ay,az,fac,r2,gc4,d
 integer    :: i

 ! initialise quadrupole and second deriv to zero
 q(:)   = 0.
 ddq(:) = 0.

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyz,axyz,pmass) &
 !$omp private(i,x,y,z,vx,vy,vz,ax,ay,az,r2) &
 !$omp reduction(+:q,ddq)
 do i=1,npart
    if (xyzh(4,i)>tiny(xyzh)) then  !if not accreted
       x  = xyzh(1,i)
       y  = xyzh(2,i)
       z  = xyzh(3,i)
       vx = vxyz(1,i)
       vy = vxyz(2,i)
       vz = vxyz(3,i)
       ax = axyz(1,i)
       ay = axyz(2,i)
       az = axyz(3,i)

       r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))

       ! calculate the components of the traceless quadrupole--not necessary but maybe useful
       q(1) = q(1) + pmass*(x*x-0.3*r2) !qxx
       q(2) = q(2) + pmass*(x*y-0.3*r2) !qxy
       q(3) = q(3) + pmass*(x*z-0.3*r2) !qxz
       q(4) = q(4) + pmass*(y*y-0.3*r2) !qyy
       q(5) = q(5) + pmass*(y*z-0.3*r2) !qyz
       q(6) = q(6) + pmass*(z*z-0.3*r2) !qzz

       ! calculate the second time derivative of the traceless quadrupole
       ddq(1) = ddq(1) + pmass*(2.*vx*vx+x*ax+x*ax) !ddqxx
       ddq(2) = ddq(2) + pmass*(2.*vx*vy+x*ay+y*ax) !ddqxy
       ddq(3) = ddq(3) + pmass*(2.*vx*vz+x*az+z*ax) !ddqxz
       ddq(4) = ddq(4) + pmass*(2.*vy*vy+y*ay+y*ay) !ddqyy
       ddq(5) = ddq(5) + pmass*(2.*vy*vz+y*az+z*ay) !ddqyz
       ddq(6) = ddq(6) + pmass*(2.*vz*vz+z*az+z*az) !ddqzz
    endif
 enddo
 !omp end parallel do

 d   = 3.10e24/udist                        ! 1Mpc in code units
 gc4 = (gg/c**4) / (utime**2/(umass*udist)) ! G/c^4 in code units
 fac = gc4/d

 ! gw strain in the direction perpendicular to the orbit, theta=0
 hp = fac*(ddq(1)-ddq(4))
 hx = 2.*fac*ddq(2)

 ! gw strain in the plane of the orbit, theta=pi/2
 hpp = fac*(ddq(1)-ddq(6))
 hxx = -2.*fac*ddq(3)

end subroutine calculate_strain

end module gravwaveutils
