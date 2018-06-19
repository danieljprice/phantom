!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:                                                            !
!  Tilts and puts a (rotating) star in a parabolic orbit around a BH       !
!                                                                          !
!  adapted by Andrea Sacchi                                                !
!                                                                          !
!  REFERENCES: None                                                        !
!                                                                          !
!  OWNER: Daniel Price                                                     !
!                                                                          !
!  $Id$                         !
!                                                                          !
!  RUNTIME PARAMETERS: None                                                !
!                                                                          !
!  DEPENDENCIES: centreofmass, externalforces, options, prompting, physcon !
!
!  REFERENCES: None                                                        !
!                                                                          !
!  OWNER: Daniel Price                                                     !
!                                                                          !
!  $Id$                         !
!                                                                          !
!  RUNTIME PARAMETERS: None                                                !
!                                                                          !
!  DEPENDENCIES: centreofmass, externalforces, options, prompting, physcon !
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: centreofmass, externalforces, options, physcon, prompting
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use centreofmass
 use externalforces, only:mass1
 use externalforces, only:accradius1
 use options,        only:iexternalforce
 use prompting,      only:prompt
 use physcon,        only:pi
 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                 :: i
 real                    :: beta, rt, rp, rs, Ms, Mh
 real                    :: Lx,Ly,Lz,L,Lp
 real                    :: phi,theta
 real                    :: x,y,z,vx,vy,vz
 real                    :: x0,y0,vx0,vy0,alpha,r0

 !--Reset center of mass
 call reset_centreofmass(npart,xyzh,vxyzu)

 phi   = 0.
 theta = 0.

 !--Calculating anuglar momentum (blame on Andrea Sacchi)

 Lx=0.0
 Ly=0.0
 Lz=0.0

 do i = 1, npart
    !
    Lx = Lx+xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i)
    Ly = Ly+xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i)
    Lz = Lz+xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i)
    !
 enddo

 L=sqrt(Lx**2.0+Ly**2.0+Lz**2.0)

 print*,'Checking angular momentum orientation and magnitude...'
 print*,'Angular momentum is L=(',Lx,Ly,Lz,')'
 print*,'Angular momentum modulus is L=',L

 Lp=sqrt(Lx**2.0+Lz**2.0)

 if (Lx > 0.) then
    phi=acos(Lz/Lp)
 elseif (Lx < 0.) then
    phi=-acos(Lz/Lp)
 endif

 print*,'tilting along y axis: ',(phi*180/pi),'degrees'

!
!--Rotate the star so the momentum lies in the yz plan
!

 do i=1,npart
    !
    x=xyzh(1,i)
    z=xyzh(3,i)
    xyzh(1,i)=x*cos(-phi)+z*sin(-phi)
    xyzh(3,i)=-x*sin(-phi)+z*cos(-phi)
    vx=vxyzu(1,i)
    vz=vxyzu(3,i)
    vxyzu(1,i)=vx*cos(-phi)+vz*sin(-phi)
    vxyzu(3,i)=-vx*sin(-phi)+vz*cos(-phi)
    !
 enddo

!
!--Recheck the stellar angular momentum
!

 Lx=0.0
 Ly=0.0
 Lz=0.0

 do i = 1, npart
    !
    Lx = Lx+xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i)
    Ly = Ly+xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i)
    Lz = Lz+xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i)
    !
 enddo

 L=sqrt(Lx**2.0+Ly**2.0+Lz**2.0)

! print*,'Angular momentum is L=(',Lx,Ly,Lz,')'
! print*,'Angular momentum modulus is L=',L


 if (Ly < 0.) then
    theta=acos(Lz/L)
 elseif (Ly > 0.) then
    theta=-acos(Lz/L)
 endif

 print*, 'tilting along x axis: ',(theta*180/pi),'degrees'

!
!--Rotate the star so the momentum lies along the z axis
!

 do i=1,npart
    y=xyzh(2,i)
    z=xyzh(3,i)
    xyzh(2,i)=y*cos(-theta)-z*sin(-theta)
    xyzh(3,i)=y*sin(-theta)+z*cos(-theta)
    vy=vxyzu(2,i)
    vz=vxyzu(3,i)
    vxyzu(2,i)=vy*cos(-theta)-vz*sin(-theta)
    vxyzu(3,i)=vy*sin(-theta)+vz*cos(-theta)
 enddo

!
!--Recheck the stellar angular momentum
!

 Lx=0.0
 Ly=0.0
 Lz=0.0

 do i = 1, npart
    !
    Lx = Lx+xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i)
    Ly = Ly+xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i)
    Lz = Lz+xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i)
    !
 enddo

 L=sqrt(Lx**2.0+Ly**2.0+Lz**2.0)

 print*,'Stellar spin is now along z axis.'
 print*,'Angular momentum is L=(',Lx,Ly,Lz,')'
 print*,'Angular momentum modulus is L=',L

 !--Defaults
 beta = 1.0                   ! penetration factor
 Mh = 1.e6                    ! BH mass
 Ms = 1.0                     ! stellar mass
 rs = 1.0                     ! stellar radius
 theta = 0.0                  ! stellar tilting along x
 phi = 0.0                    ! stellar tilting along y

 !--User enter values
 call prompt(' Enter a value for the penetration factor (beta): ',beta,0.)
 call prompt(' Enter a value for blackhole mass (in code units): ',Mh,0.)
 call prompt(' Enter a value for the stellar mass (in code units): ',Ms,0.)
 call prompt(' Enter a value for the stellar radius (in code units): ',rs,0.)
 call prompt(' Enter a value for the stellar rotation with respect to x-axis (in degrees): ',theta,0.)
 call prompt(' Enter a value for the stellar rotation with respect to y-axis (in degrees): ',phi,0.)


 rt = (Mh/Ms)**(1./3.) * rs   ! tidal radius
 rp = rt/beta                 ! pericenter distance

 alpha = 3./4.*pi  ! Starting angle from x-axis (anti-clockwise)
 ! alpha = 2.2      ! Most time efficient (Elen)
 r0    = 2.*rt/beta*1./(1.+cos(alpha)) ! Starting radius

 x0    = r0*cos(alpha)
 y0    = r0*sin(alpha)
 vx0   = sqrt(mh*beta/(2.*rt)) * sin(alpha)
 vy0   = -sqrt(mh*beta/(2.*rt)) * (cos(alpha)+1.)

 !--Set input file parameters
 mass1 = Mh
 iexternalforce = 1
 accradius1 = (2*Mh*rs)/((6.8565e2)**2) ! R_sch = 2*G*Mh*rs/c**2

 !--Tilting the star
 theta=theta*pi/180.0
 phi=phi*pi/180.0

 if (theta  /=  0.0) then
    do i=1,npart
       y=xyzh(2,i)
       z=xyzh(3,i)
       xyzh(2,i)=y*cos(theta)-z*sin(theta)
       xyzh(3,i)=y*sin(theta)+z*cos(theta)
       vy=vxyzu(2,i)
       vz=vxyzu(3,i)
       vxyzu(2,i)=vy*cos(theta)-vz*sin(theta)
       vxyzu(3,i)=vy*sin(theta)+vz*cos(theta)
    enddo
 endif
 if (phi  /=  0.0) then
    do i=1,npart
       x=xyzh(1,i)
       z=xyzh(3,i)
       xyzh(1,i)=x*cos(phi)+z*sin(phi)
       xyzh(3,i)=-x*sin(phi)+z*cos(phi)
       vx=vxyzu(1,i)
       vz=vxyzu(3,i)
       vxyzu(1,i)=vx*cos(phi)+vz*sin(phi)
       vxyzu(3,i)=-vx*sin(phi)+vz*cos(phi)
    enddo
 endif

 !--Putting star into orbit
 do i = 1, npart
    xyzh(1,i)  = xyzh(1,i)  + x0
    xyzh(2,i)  = xyzh(2,i)  + y0
    vxyzu(1,i) = vxyzu(1,i) + vx0
    vxyzu(2,i) = vxyzu(2,i) + vy0
 enddo

 theta=theta*pi/180.0
 phi=phi*pi/180.0

 write(*,'(a)') "======================================================================"
 write(*,'(a,Es12.5,a)') ' Pericenter distance = ',rp,' R_sun'
 write(*,'(a,Es12.5,a)') ' Tidal radius        = ',rt,' R_sun'
 write(*,'(a,Es12.5,a)') ' Radius of star      = ',rs,' R_sun'
 write(*,'(a,Es12.5,a)') ' Starting distance   = ',r0,' R_sun'
 write(*,'(a,Es12.5,a)') ' Stellar mass        = ',Ms,' M_sun'
 write(*,'(a,Es12.5,a)') ' Tilting along x     = ',theta,' degrees'
 write(*,'(a,Es12.5,a)') ' Tilting along y     = ',phi,' degrees'

 write(*,'(a)') "======================================================================"

 return
end subroutine modify_dump

end module moddump
