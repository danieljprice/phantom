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
!  DEPENDENCIES: boundary
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 private

contains

subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,fileprefix)
 use boundary, only:xmin,ymin,zmin,dxbound,dybound,dzbound
 integer,           intent(in)  :: id
 integer,           intent(out) :: npart
 integer,           intent(out) :: npartoftype(:)
 real,              intent(out) :: xyzh(:,:)
 real,              intent(out) :: polyk,gamma,hfact
 real,              intent(out) :: vxyzu(:,:)
 real,              intent(out) :: massoftype(:)
 character(len=20), intent(in)  :: fileprefix
 integer, parameter :: mx=100,my=100,mz=100
 real :: rhogrid(mx,my,mz)
 real :: xi,yi,zi,dx,dy,dz
 integer :: i,j,k,maxp

 !--set a rho distribution on a grid, as a test case
 maxp = size(xyzh(1,:))
 dz = dzbound/mz
 dy = dybound/my
 dx = dxbound/mx
 do k=1,mz
    zi = zmin + (k-0.5)*dz
    do j=1,my
       yi = ymin + (j-0.5)*dy
       do i=1,mx
          xi = xmin + (i-0.5)*dx
          rhogrid(i,j,k) = rhofunc(xi,yi,zi)
       enddo
    enddo
 enddo
 call positions(mx,my,mz,rhogrid,maxp,npart,xyzh(1,:),xyzh(2,:),xyzh(3,:),xyzh(4,:),massoftype(1))

 !!$omp parallel do private(i)
 do i=1,npart  ! normalize to unit interval
    xyzh(1,i) = xyzh(1,i)*(1./mx)
    xyzh(2,i) = xyzh(2,i)*(1./my)
    xyzh(3,i) = xyzh(3,i)*(1./mz)
 enddo
 print*,' setting particle mass = ',massoftype(1)
 npartoftype(:) = 0
 npartoftype(1) = npart

! utherm = 1.
! print*,'setting utherm = ',utherm
 polyk = 1.  !!2./3.*utherm
 print*,'polyk = ',polyk
! polyk = 0.
 hfact = 1.2
!--set h based on rho
 do i=1,npart
    xyzh(4,i) = hfact*(massoftype(1)/xyzh(4,i))**(1./3)
 enddo
 vxyzu(:,i) = 0.
 gamma = 1.

 return
end subroutine setpart

real function rhofunc(xi,yi,zi)
 real, parameter :: pi = 3.1415926536
 real, intent(in) :: xi,yi,zi

 rhofunc = 2. + sin(2.*pi*xi)
 !rhofunc = max(2. - dot_product(xi,xi),1.)
 !rhofunc = 1./exp(5.*((xi(1)-0.5)**2 + (xi(2)-0.5)**2))

end function rhofunc

end module setup

