!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inverse4x4
!
! Invert a 4 x 4 matrix. Returns zero if determinant is zero
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
contains

pure subroutine inv4x4(A,B,det)
 real, intent(in), dimension(4,4) :: A
 real, intent(out), dimension(4,4) :: B
 real, intent(out) :: det
 real :: a11,a12,a13,a14
 real :: a21,a22,a23,a24
 real :: a31,a32,a33,a34
 real :: a41,a42,a43,a44
 real :: det1

 a11 = A(1,1)
 a21 = A(2,1)
 a31 = A(3,1)
 a41 = A(4,1)
 a12 = A(1,2)
 a22 = A(2,2)
 a32 = A(3,2)
 a42 = A(4,2)
 a13 = A(1,3)
 a23 = A(2,3)
 a33 = A(3,3)
 a43 = A(4,3)
 a14 = A(1,4)
 a24 = A(2,4)
 a34 = A(3,4)
 a44 = A(4,4)

 det = a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 + &
       a13*a22*a34*a41 - a12*a23*a34*a41 - a14*a23*a31*a42 + a13*a24*a31*a42 + &
       a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 + &
       a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 + &
       a12*a21*a34*a43 - a11*a22*a34*a43 - a13*a22*a31*a44 + a12*a23*a31*a44 + &
       a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44

 if (abs(det) < tiny(0.)) then
    B = 0.
    return
 endif

 det1 = 1./det

 B(1,1) = -(a24*a33*a42) + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 + a22*a33*a44
 B(2,1) =    a14*a33*a42 - a13*a34*a42 - a14*a32*a43 + a12*a34*a43 + a13*a32*a44 - a12*a33*a44
 B(3,1) = -(a14*a23*a42) + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 + a12*a23*a44
 B(4,1) =    a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*a34 - a12*a23*a34
 B(1,2) =    a24*a33*a41 - a23*a34*a41 - a24*a31*a43 + a21*a34*a43 + a23*a31*a44 - a21*a33*a44
 B(2,2) = -(a14*a33*a41) + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 + a11*a33*a44
 B(3,2) =    a14*a23*a41 - a13*a24*a41 - a14*a21*a43 + a11*a24*a43 + a13*a21*a44 - a11*a23*a44
 B(4,2) = -(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 + a11*a23*a34
 B(1,3) = -(a24*a32*a41) + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 + a21*a32*a44
 B(2,3) =    a14*a32*a41 - a12*a34*a41 - a14*a31*a42 + a11*a34*a42 + a12*a31*a44 - a11*a32*a44
 B(3,3) = -(a14*a22*a41) + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 + a11*a22*a44
 B(4,3) =    a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*a34 - a11*a22*a34
 B(1,4) =    a23*a32*a41 - a22*a33*a41 - a23*a31*a42 + a21*a33*a42 + a22*a31*a43 - a21*a32*a43
 B(2,4) = -(a13*a32*a41) + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 + a11*a32*a43
 B(3,4) =    a13*a22*a41 - a12*a23*a41 - a13*a21*a42 + a11*a23*a42 + a12*a21*a43 - a11*a22*a43
 B(4,4) = -(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33

 B = B*det1

end subroutine inv4x4


end module inverse4x4
