module gravwaveutils
 implicit none

 public :: calculate_strain

 private

contains

! This subroutine computes the gravitational wave strain at a distance of 1Mpc
subroutine calculate_strain(hx,hp,pmass,x0,v0,a0,npart,xyzh,vxyz,axyz,axyz1)
 use units,                    only:umass,udist,utime
 use physcon,                  only:gg,c,Mpc,pi
 real, intent(out)             :: hx(4),hp(4)
 real, intent(in)              :: xyzh(:,:), vxyz(:,:), axyz(:,:), pmass,x0(3),v0(3),a0(3)
 real, intent(inout), optional :: axyz1(:,:) !optional, only if there are external forces
 integer, intent(in)           :: npart
 real                          :: q(6),ddq(6),x,y,z,vx,vy,vz,ax,ay,az,fac,r2,gc4,d
 real, parameter               :: onethird = 1./3.
 integer                       :: i
 real                          :: theta,phi,sinphi,cosphi,sintheta,costheta,cosphi2,sinphi2,&
                                  cos2phi,sin2phi,cos2theta,sin2theta,sintheta2,costheta2

!change this line if you want a different value of phi
 phi=0.
 theta=0.

!initialise quadrupole and second deriv to zero
 q(:)   = 0.
 ddq(:) = 0.

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyz,axyz,axyz1,pmass,x0,v0,a0) &
 !$omp private(i,x,y,z,vx,vy,vz,ax,ay,az,r2) &
 !$omp reduction(+:q,ddq)
 do i=1,npart
    if (xyzh(4,i) > tiny(xyzh)) then  !if not accreted
      x  = xyzh(1,i) - x0(1)
      y  = xyzh(2,i) - x0(2)
      z  = xyzh(3,i) - x0(3)
      vx = vxyz(1,i) - v0(1)
      vy = vxyz(2,i) - v0(2)
      vz = vxyz(3,i) - v0(3)
      ax = axyz(1,i) - a0(1)
      ay = axyz(2,i) - a0(2)
      az = axyz(3,i) - a0(3)
      if (present(axyz1)) then ! to avoid memory allocation in phantom
         ax = ax + axyz1(1,i)
         ay = ay + axyz1(2,i)
         az = az + axyz1(3,i)
      endif
      r2 = dot_product(xyzh(1:3,i),xyzh(1:3,i))

      ! calculate the components of the traceless quadrupole--not necessary but maybe useful
       q(1) = q(1) + pmass*(x*x-onethird*r2) !qxx
       q(2) = q(2) + pmass*(x*y-onethird*r2) !qxy
       q(3) = q(3) + pmass*(x*z-onethird*r2) !qxz
       q(4) = q(4) + pmass*(y*y-onethird*r2) !qyy
       q(5) = q(5) + pmass*(y*z-onethird*r2) !qyz
       q(6) = q(6) + pmass*(z*z-onethird*r2) !qzz

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

!define some parameters
 d   = Mpc/udist                        ! 1Mpc in code units
 gc4 = (gg/c**4) / (utime**2/(umass*udist)) ! G/c^4 in code units
 fac = gc4/d

!derive the GW waveforms
 do i=1,4
 !define new functions instead of sin and cos
 sinphi=sin(phi)
 cosphi=cos(phi)
 sinphi2=sinphi*sinphi
 cosphi2=cosphi*cosphi
 sin2phi=sin(2*phi)
 cos2phi=cos(2*phi)
 sintheta=sin(theta)
 costheta=cos(theta)
 sintheta2=sintheta*sintheta
 costheta2=costheta*costheta
 sin2theta=sin(2*theta)
 cos2theta=cos(2*theta)

  hp(i) =fac*(ddq(1)*(cosphi2 - sinphi2*costheta2)+ddq(4)*(sinphi2 &
  -cosphi2*costheta2)-ddq(6)*sintheta2-ddq(2)*(sin2phi)*(1+costheta2) &
  +ddq(3)*(sinphi*sin2theta)+ddq(5)*cosphi*sin2theta)
  hx(i)=2.*fac*(0.5*(ddq(1)-ddq(4))*sin2phi*costheta+ddq(2)*(cos2phi*costheta)&
  -ddq(3)*(cosphi*sintheta)+ddq(5)*sinphi*sintheta)
  theta=theta+pi/6.
 enddo
print*, 'collo'

end subroutine calculate_strain

end module gravwaveutils
