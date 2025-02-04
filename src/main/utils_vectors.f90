!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module vectorutils
!
! This module contains utilities for manipulating vectors
!  and 3x3 matrices
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: minmaxave,cross_product3D,curl3D_epsijk,det
 public :: matrixinvert3D,rotatevec,unitvec,mag

 private

contains
!-------------------------------------------------------------------
!+
!  find min, max and average of an array
!+
!-------------------------------------------------------------------
subroutine minmaxave(x,xmin,xmax,xav,npts)
 integer :: i
 integer, intent(in)  :: npts
 real,    intent(in)  :: x(npts)
 real,    intent(out) :: xmin,xmax,xav

 xav = 0.
 xmin = huge(xmin)
 xmax = -xmin
 do i=1,npts
    xav = xav + x(i)
    xmin = min(xmin,x(i))
    xmax = max(xmax,x(i))
 enddo
 xav = xav/real(npts)

end subroutine minmaxave

!-------------------------------------------------------------------
!+
!  vector cross product
!+
!-------------------------------------------------------------------
pure subroutine cross_product3D(veca,vecb,vecc)
 real, intent(in)  :: veca(3),vecb(3)
 real, intent(out) :: vecc(3)

 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end subroutine cross_product3D

!-------------------------------------------------------------------
!+
!  curl from the 3 x 3 gradient matrix
!+
!-------------------------------------------------------------------
pure subroutine curl3D_epsijk(gradAvec,curlA)
 real, intent(in)  :: gradAvec(3,3)
 real, intent(out) :: curlA(3)

 curlA(1) = gradAvec(2,3) - gradAvec(3,2)
 curlA(2) = gradAvec(3,1) - gradAvec(1,3)
 curlA(3) = gradAvec(1,2) - gradAvec(2,1)

end subroutine curl3D_epsijk

!----------------------------------------------------------------
!+
!  Inverts a 3x3 matrix
!+
!----------------------------------------------------------------
subroutine matrixinvert3D(A,Ainv,ierr)
 real,    intent(in)  :: A(3,3)
 real,    intent(out) :: Ainv(3,3)
 integer, intent(out) :: ierr
 real :: x0(3),x1(3),x2(3),result(3)
 real    :: det, ddet

 ierr = 0

 x0 = A(1,:)
 x1 = A(2,:)
 x2 = A(3,:)

 call cross_product3D(x1,x2,result)
 det = dot_product(x0,result)

 if (abs(det) > tiny(det)) then
    ddet = 1./det
 else
    ddet = 0.
    Ainv = 0.
    ierr = 1
    return
 endif

 Ainv(:,1) = result(:)*ddet
 call cross_product3D(x2,x0,result)
 Ainv(:,2) = result(:)*ddet
 call cross_product3D(x0,x1,result)
 Ainv(:,3) = result(:)*ddet

end subroutine matrixinvert3D

!----------------------------------------------------------------
!+
!  Determinant of a 3x3 matrix
!+
!----------------------------------------------------------------
real function det(A)
 real, intent(in) :: A(3,3)
 real :: x0(3),x1(3),x2(3),result(3)

 x0 = A(1,:)
 x1 = A(2,:)
 x2 = A(3,:)

 call cross_product3D(x1,x2,result)
 det = dot_product(x0,result)

end function det

!------------------------------------------------------------------------
!+
!  rotate a vector (u) around an axis defined by another vector (v)
!  by an angle (theta) using the Rodrigues rotation formula
!+
!------------------------------------------------------------------------
pure subroutine rotatevec(u,v,theta)
 real, dimension(3), intent(inout) :: u
 real, dimension(3), intent(in)    :: v
 real, intent(in)   :: theta
 real, dimension(3) :: k,w

 !--normalise v
 k = v/sqrt(dot_product(v,v))
 !--Rodrigues rotation formula
 call cross_product3D(k,u,w)
 u = u*cos(theta) + w*sin(theta) + k*dot_product(k,u)*(1-cos(theta))
end subroutine rotatevec

!------------------------------------------------------------------------
!+
!  return unit vector given a vector
!+
!------------------------------------------------------------------------
pure function unitvec(u) result(uhat)
 real, intent(in) :: u(3)
 real :: uhat(3),u2

 u2 = dot_product(u,u)
 if (u2 > tiny(0.)) then
    uhat = u/sqrt(u2)
 else
    uhat = (/0.,0.,1./)  ! arbitrary if vector is zero
 endif

end function unitvec

!------------------------------------------------------------------------
!+
!  magnitude of a vector
!+
!------------------------------------------------------------------------
pure function mag(u) result(umag)
 real, intent(in) :: u(3)
 real :: umag

 umag = sqrt(dot_product(u,u))

end function mag

end module vectorutils
