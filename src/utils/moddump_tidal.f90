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
 use options,        only:iexternalforce,damp
 use prompting,      only:prompt
 use physcon,        only:pi
 integer,  intent(inout) :: npart
 integer,  intent(inout) :: npartoftype(:)
 real,     intent(inout) :: massoftype(:)
 real,     intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                 :: i
 real                    :: beta, rt, rp, rs, Ms, Mh
 real                    :: Lx,Ly,Lz,L,Lp,Ltot(3)
 real                    :: phi,theta
 real                    :: x,y,z,vx,vy,vz
 real                    :: x0,y0,vx0,vy0,alpha,r0,ecc

 !--Reset center of mass
 call reset_centreofmass(npart,xyzh,vxyzu)

 phi   = 0.
 theta = 0.

 call get_angmom(ltot,npart,xyzh,vxyzu)
 Lx = ltot(1)
 Ly = ltot(2)
 Lz = ltot(3)
 Lp = sqrt(Lx**2.0+Lz**2.0)

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
    x=xyzh(1,i)
    z=xyzh(3,i)
    xyzh(1,i)=x*cos(-phi)+z*sin(-phi)
    xyzh(3,i)=-x*sin(-phi)+z*cos(-phi)
    vx=vxyzu(1,i)
    vz=vxyzu(3,i)
    vxyzu(1,i)=vx*cos(-phi)+vz*sin(-phi)
    vxyzu(3,i)=-vx*sin(-phi)+vz*cos(-phi)
 enddo

!
!--Recheck the stellar angular momentum
!

 call get_angmom(ltot,npart,xyzh,vxyzu)
 lx = ltot(1)
 ly = ltot(2)
 lz = ltot(3)
 L  = sqrt(Lx**2.0+Ly**2.0+Lz**2.0)
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

 call get_angmom(ltot,npart,xyzh,vxyzu)
 print*,'Stellar spin should now be along the z axis.'




 !--Defaults
 beta = 1.0                   ! penetration factor
 Mh = 1.e6                    ! BH mass
 Ms = 1.0                     ! stellar mass
 rs = 1.0                     ! stellar radius
 theta = 0.0                  ! stellar tilting along x
 phi = 0.0                    ! stellar tilting along y
 ecc = 1.                     ! eccentricity

 !--User enter values
 call prompt(' Enter a value for the penetration factor (beta): ',beta,0.)
 call prompt(' Enter a value for blackhole mass (in code units): ',Mh,0.)
 call prompt(' Enter a value for the stellar mass (in code units): ',Ms,0.)
 call prompt(' Enter a value for the stellar radius (in code units): ',rs,0.)
 call prompt(' Enter a value for the eccentricity: ',ecc,0.,1.)
 rt = (Mh/Ms)**(1./3.) * rs         ! tidal radius
 rp = rt/beta                       ! pericenter distance
 r0 = 4.9*rt                        ! starting radius
 call prompt(' Enter a value for the stellar rotation with respect to x-axis (in degrees): ',theta,0.)
 call prompt(' Enter a value for the stellar rotation with respect to y-axis (in degrees): ',phi,0.)
 call prompt(' Enter a value for the starting distance (in code units): ',r0,0.)

 alpha = acos((rt*(1.+ecc)/(r0*beta)-1.)/ecc)         ! starting angle anti-clockwise from positive x-axis
 x0    = r0*cos(alpha)
 y0    = r0*sin(alpha)
 vx0   = sqrt(mh*beta/((1.+ecc)*rt)) * sin(alpha)
 vy0   = -sqrt(mh*beta/((1.+ecc)*rt)) * (cos(alpha)+ecc)

 !--Set input file parameters
 mass1          = Mh
 iexternalforce = 1
 damp           = 0.
 accradius1     = (2*Mh*rs)/((6.8565e2)**2) ! R_sch = 2*G*Mh*rs/c**2

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


subroutine get_angmom(ltot,npart,xyzh,vxyzu)
 real, intent(out)   :: ltot(3)
 integer, intent(in) :: npart
 real, intent(in)    :: xyzh(:,:), vxyzu(:,:)
 integer :: i
 real    :: L

 !--Calculating anuglar momentum (blame on Andrea Sacchi [and David Liptai :P])

 ltot = 0.
 do i=1,npart
     ltot(1) = ltot(1)+xyzh(2,i)*vxyzu(3,i)-xyzh(3,i)*vxyzu(2,i)
     ltot(2) = ltot(2)+xyzh(3,i)*vxyzu(1,i)-xyzh(1,i)*vxyzu(3,i)
     ltot(3) = ltot(3)+xyzh(1,i)*vxyzu(2,i)-xyzh(2,i)*vxyzu(1,i)
 enddo

 L = sqrt(dot_product(ltot,ltot))

 print*,''
 print*,'Checking angular momentum orientation and magnitude...'
 print*,'Angular momentum is L = (',ltot(1),ltot(2),ltot(3),')'
 print*,'Angular momentum modulus is |L| = ',L
 print*,''

end subroutine get_angmom

end module moddump
