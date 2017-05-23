!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2015 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: utils_gr
!
!  DESCRIPTION:
!   Contains utility routines for general relativity
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module utils_gr
implicit none

public :: dot_product_gr, get_metric3plus1, get_u0, get_rderivs, get_ev, rho_to_dens, h2dens_all

private

interface get_metric3plus1
module procedure get_metric3plus1_only, get_metric3plus1_both
end interface get_metric3plus1

interface rho_to_dens
module procedure rho2dens, h2dens
end interface rho_to_dens

contains

!----------------------------------------------------------------
!+
!  Function to perform a dot product in general relativity
!  i.e. g_\mu\nu v^\mu \v^nu
!+
!----------------------------------------------------------------
pure real function dot_product_gr(vec1,vec2,gcov)
   real, intent(in) :: vec1(:)
   real, intent(in) :: vec2(size(vec1))
   real, intent(in) :: gcov(size(vec1),size(vec2))
   real :: vec1i
   integer :: i,j

   dot_product_gr = 0.
   do i=1,size(vec1)
      vec1i = vec1(i)
      do j=1,size(vec2)
         dot_product_gr = dot_product_gr + gcov(j,i)*vec1i*vec2(j)
      enddo
   enddo

   return
end function dot_product_gr

subroutine get_metric3plus1_only(x,alpha,beta,gammaijdown,gammaijUP)
   real, intent(in)  :: x(1:3)
   real, intent(out) :: alpha, beta(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
   real              :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg

   call metric3p1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)

end subroutine get_metric3plus1_only

subroutine get_metric3plus1_both(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
   real, intent(in)  :: x(1:3)
   real, intent(out) :: alpha,beta(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
   real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg

   call metric3p1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)

end subroutine get_metric3plus1_both

subroutine metric3p1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
   use metric_tools, only: get_metric
   real, intent(in)  :: x(1:3)
   real, intent(out) :: alpha,beta(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
   real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
   real :: betaUP(1:3)
   integer :: i,j

   call get_metric(x,gcov,gcon,sqrtg)
   beta  = gcov(0,1:3)
   gammaijdown   = gcov(1:3,1:3)
   alpha = sqrt(-1./gcon(0,0))
   betaUP = gcon(0,1:3)*alpha**2
   gammaijUP = 0.
   do i=1,3
      do j=1,3
         gammaijUP(i,j) = gcon(i,j) + betaUP(i)*betaUP(j)/alpha**2
      enddo
   enddo
end subroutine metric3p1

subroutine get_u0(x,v,U0)
   use metric_tools, only: get_metric
   real, intent(in) :: x(1:3),v(1:3)
   real, intent(out) :: U0
   real :: v4(0:3), gcov(0:3,0:3), gcon(0:3,0:3), sqrtg

   call get_metric(x,gcov,gcon,sqrtg)
   v4(0) = 1.
   v4(1:3) = v(1:3)
   U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))
end subroutine get_u0

subroutine get_rderivs(position,dr)
   use metric, only: a
   real, dimension(3), intent(in) :: position
   real, dimension(3), intent(out) :: dr
   real :: x,y,z,dxdr,dydr,dzdr,r,r2sphere

   x=position(1)
   y=position(2)
   z=position(3)

   r2sphere = x**2+y**2+z**2

   r = sqrt(0.5*(r2sphere-a**2+sqrt((r2sphere-a**2))**2 + 4.*a**2*z**2))

   dxdr=x*r/(r**2+a**2)
   dydr=y*r/(r**2+a**2)
   dzdr=z/r

   dr=(/dxdr,dydr,dzdr/)

end subroutine get_rderivs

subroutine get_ev(x,v,energy,angmom)
   use metric, only: metric_type, rs
   use metric_tools, only: coordinate_sys
   real, intent(in), dimension(3) :: x,v
   real, intent(out) :: energy, angmom
   real :: r, U0
   integer, save :: i = 0
   character(len=*), parameter :: force_type = 'GR'

   ! For Schwarzschild only
   if (metric_type=='Schwarzschild') then
      call get_u0(x,v,U0)
      if (coordinate_sys=='Cartesian') then
         r      = sqrt(dot_product(x,x))
         energy = (1. - rs/r)*U0
         angmom = (x(1)*v(2)-x(2)*v(1))*U0
      else if (coordinate_sys=='Spherical') then
         energy = (1. - rs/x(1))*U0
         angmom = x(1)**2*v(3)*U0
      endif
   else if (metric_type == 'Minkowski' .and. force_type == 'Newtonian') then
      energy = 0.5*dot_product(v,v)
      angmom = x(1)*v(2)-x(2)*v(1)
   else
      if (i==0) then
         i = i+1
         print*,'WARNING: Energy and angular momentum are not being calculated for this metric. They will just be set to zero.'
         print*,'Continue?'
         read*
         energy = 0.
         angmom = 0.
      endif
   endif

end subroutine get_ev

subroutine rho2dens(dens,rho,position,v)
   use metric_tools, only: get_metric
   real, intent(in) :: rho,position(1:3),v(1:3)
   real, intent(out):: dens
   real :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg, U0

   call get_metric(position,gcov,gcon,sqrtg)
   call get_u0(position,v,U0)

   dens = rho/(sqrtg*U0)

end subroutine rho2dens

subroutine h2dens(dens,xyzh,v)
   use metric_tools, only: get_metric
   use part,         only: rhoh,massoftype,igas
   real, intent(in) :: xyzh(1:4),v(1:3)
   real, intent(out):: dens
   real :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg, U0, rho

   rho =  rhoh(xyzh(4),massoftype(igas))
   call get_metric(xyzh(1:3),gcov,gcon,sqrtg)
   call get_u0(xyzh(1:3),v,U0)

   dens = rho/(sqrtg*U0)

end subroutine h2dens

subroutine h2dens_all(npart,dens,xyzh,vxyzu)
   use part, only:isdead_or_accreted
   real, intent(in) :: xyzh(:,:),vxyzu(:,:)
   real, intent(out):: dens(:)
   integer, intent(in) :: npart
   integer :: i

!$omp parallel do default (none) &
!$omp shared(xyzh,vxyzu,dens,npart) &
!$omp private(i)
   do i=1,npart
      if (.not.isdead_or_accreted(xyzh(4,i))) call h2dens(dens(i),xyzh(:,i),vxyzu(1:3,i))
   enddo
!$omp end parallel do
end subroutine h2dens_all

subroutine dens_to_rho(rho,dens,position,v)
   use metric_tools, only: get_metric
   real, intent(in) :: dens,position(1:3),v(1:3)
   real, intent(out):: rho
   real :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg, U0

   call get_metric(position,gcov,gcon,sqrtg)
   call get_u0(position,v,U0)

   rho = sqrtg*U0*dens

end subroutine dens_to_rho

end module utils_gr
