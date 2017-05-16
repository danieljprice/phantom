!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
! this module does setup
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: externalforces, io, part, physcon, setup_params,
!    spherical, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,maxvxyzu
 use spherical, only:set_sphere
 use io,        only:master
 use setup_params, only:rhozero,npart_total
 use externalforces, only:mass1
 use units,          only:set_units
 use physcon,        only:solarm,pc
 real, parameter :: pi = 3.1415926536
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: massr,m1,m2,a,hacc1,hacc2,ecc
 real :: P1,P2,P3,Q1,Q2,Q3,big_omega,omega,inc,E,v1,v2,E_dot
 real :: rmin,rmax,psep,totvol,dr,t,dt,totmass,r
 integer :: i,nx,np
 integer, parameter :: ng  = 1024
! real, dimension(ng) :: r,den

 call set_units(dist=1.d-3*pc,mass=solarm,G=1.)
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 polyk = 412.4951
 gamma = 5./3.
 rmin = 0.3
 rmax = 81.0
!
!--space available for injected gas particles
!
 m1    = 4.31E6
 massr = 2.09d-12
 a     = 42.227
 ecc   = 0.9762
 hacc1  = rmin ! 1.7 AU
 hacc2  = 0.017
 m2   = m1*massr

!--set input file options
 mass1 = m1

!*** SET UP MEDIUM***!

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.

 np = 10000 !size(xyzh(1,:))
 totvol = 4./3.*pi*rmax**3
 nx = int(np**(1./3.))
 psep = totvol**(1./3.)/real(nx)
 print*,' Setup for medium around black hole '
 print*,' total volume = ',totvol,' particle separation = ',psep
 print*,totvol/psep**3

 npart = 0
 npart_total = 0

 call set_sphere('closepacked',id,master,rmin,rmax,psep,hfact,npart,xyzh,rhofunc=rhofunc,nptot=npart_total)

!
!--set particle properties
!
 totmass = 0.74985E-3!1.87d-3
 rhozero = totmass/totvol
 polyk   = 41249.51
 print*,' total mass = ',totmass,' mean density = ',rhozero,' polyk = ',polyk
 print*,' free fall time = ',sqrt(3.*pi/(32.*rhozero))

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 massoftype(1) = totmass/npart_total
 print*,' particle mass = ',massoftype(1)

 do i=1,npart
    vxyzu(1:3,i) = 0.
    !
    !--Note: Running the polytrope with u stored is not quite
    !  the same as using P = K rho^gamma because we really
    !  should use the actual rho, not rhozero.
    !
    r = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))
    if (maxvxyzu >= 4) then
       vxyzu(4,i) = m1/(2*(gamma-1)*r)+ 3.0E4
    endif
 enddo


!*** SET UP THE POINT MASSES***!

! Set up the particles
 nptmass = 0

 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 ! Parameters for tilted orbit
 big_omega = 1.4294   ! Position angle of ascending node
 omega = 1.69646      ! Angle to periapse
 inc = 2.0612         ! Inclination of orbit
 E = 5.2457           ! Found externally for 2002 when periapse at 2014.25

 ! Positions in plane
 P1 = cos(omega)*cos(big_omega) - sin(omega)*cos(inc)*sin(big_omega)
 P2 = cos(omega)*sin(big_omega) + sin(omega)*cos(inc)*cos(big_omega)
 P3 = sin(omega)*sin(inc)
 Q1 = -sin(omega)*cos(big_omega) - cos(omega)*cos(inc)*sin(big_omega)
 Q2 = -sin(omega)*sin(big_omega) + cos(omega)*cos(inc)*cos(big_omega)
 Q3 = sin(inc)*cos(omega)

 v1 = cos(E)-ecc
 v2 = sqrt(1-(ecc*ecc))*sin(E)
 E_dot = sqrt((m1 + m2)/(a**3))/(1-ecc*cos(E))

 xyzmh_ptmass(4,1) = m1
! xyzmh_ptmass(4,2:nptmass) = m2

! Create a ring in the x-y and x-z plane
! dt = pi/4  ! Angle in radians
! t = 0                  ! To start
!  do i=2,9
!   t = t + ((i-2)*dt)
!   xyzmh_ptmass(1,i) = rmax*cos(t)
!   xyzmh_ptmass(2,i) = rmax*sin(t)
!  enddo

! Horrible hack
! t = 0                  ! To start
!  do i=10,nptmass
!   t = t + ((i-1)*dt)
!   xyzmh_ptmass(1,i) = rmax*cos(t)
!   xyzmh_ptmass(3,i) = rmax*sin(t)
!  enddo

 ! Rotating everything

!do i=1,npart
! xyzh(1,i) = a*(v1*P1 + v2*Q1) + xyzh(1,i)
! xyzh(2,i) = a*(v1*P2 + v2*Q2) + xyzh(2,i)
! xyzh(3,i) = a*(v1*P3 + v2*Q3) + xyzh(3,i)

 ! Set the velocities
! vxyzu(1,i) = -a*sin(E)*E_dot*P1 + a*sqrt(1-(ecc*ecc))*cos(E)*E_dot*Q1
! vxyzu(2,i) = -a*sin(E)*E_dot*P2 + a*sqrt(1-(ecc*ecc))*cos(E)*E_dot*Q2
! vxyzu(3,i) = -a*sin(E)*E_dot*P3 + a*sqrt(1-(ecc*ecc))*cos(E)*E_dot*Q3
!enddo

!do i=2,nptmass
! xyzmh_ptmass(1,i) = a*(v1*P1 + v2*Q1) + xyzmh_ptmass(1,i)
! xyzmh_ptmass(2,i) = a*(v1*P2 + v2*Q2) + xyzmh_ptmass(2,i)
! xyzmh_ptmass(3,i) = a*(v1*P3 + v2*Q3) + xyzmh_ptmass(3,i)

 ! Set the velocities
! vxyz_ptmass(1,i) = -a*sin(E)*E_dot*P1 + a*sqrt(1-(ecc*ecc))*cos(E)*E_dot*Q1
! vxyz_ptmass(2,i) = -a*sin(E)*E_dot*P2 + a*sqrt(1-(ecc*ecc))*cos(E)*E_dot*Q2
! vxyz_ptmass(3,i) = -a*sin(E)*E_dot*P3 + a*sqrt(1-(ecc*ecc))*cos(E)*E_dot*Q3
!enddo

 !As in set_binary
 xyzmh_ptmass(ihacc,1) = hacc1
! xyzmh_ptmass(ihacc,2:nptmass) = hacc2
 xyzmh_ptmass(ihsoft,1) = 0.0
! xyzmh_ptmass(ihsoft,2:nptmass) = 0.0

! print*,'Positions',xyzmh_ptmass
! print*,'Velocities',vxyz_ptmass
end subroutine setpart

real function rhofunc(r)
 real, intent(in) :: r

 rhofunc = 1./r

end function rhofunc

end module setup
