!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: fastmath
!
!  DESCRIPTION:
!   This module computes a fast inverse square root
!   Algorithm originally in C by Chris Lomont,
!   implemented here in Fortran by Daniel Price
!
!   Provides generic interface finvsqrt referring to
!   both the single and double precision versions
!
!   This is a standalone module which has no dependencies
!
!   For fast sqrt, just use sqrt(x) = x*finvsqrt(x)
!
!  REFERENCES:
!    Chris Lomont, "Fast Inverse Square Root"
!    www.math.purdue.edu/~clomont
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
module fastmath
 implicit none

 private
 public :: finvsqrt,testsqrt,checksqrt

 interface finvsqrt
  module procedure fastinvsqrt,dfastinvsqrt
 end interface

contains

!----------------------------------------------------------------
!+
!  For fast sqrt, just use sqrt(x) = x*(1/sqrt(x))
!+
!----------------------------------------------------------------
real(kind=4) pure elemental function fastinvsqrt(y)
 real(kind=4), intent(in) :: y
 real(kind=4) :: x,xhalf
 integer(kind=4) :: ix
 !equivalence(ix,x)

 x = y
 xhalf = 0.5_4*x
 ix = transfer(x,ix)
 ix = 1597463174 - ISHFT(ix,-1)
 x = transfer(ix,x)
 x = x*(1.5_4 - xhalf*x*x)
 x = x*(1.5_4 - xhalf*x*x)
 x = x*(1.5_4 - xhalf*x*x)
 fastinvsqrt = x

 return
end function

!-----------------------------------------
!+
!  fast inverse sqrt in double precision
!+
!-----------------------------------------
real(kind=8) pure elemental function dfastinvsqrt(y)
 real(kind=8), intent(in) :: y
 real(kind=8) :: x,xhalf
 integer(kind=8) :: ix

 x = y
 xhalf = 0.5d0*x
 ix = transfer(x,ix)
 ix = 6910470738111508698_8 - ISHFT(ix,-1)
 x = transfer(ix,x)
 x = x*(1.5d0 - xhalf*x*x)
 x = x*(1.5d0 - xhalf*x*x)
 x = x*(1.5d0 - xhalf*x*x)
 x = x*(1.5d0 - xhalf*x*x)
 dfastinvsqrt = x

 return
end function

!----------------------------------------------------------------
!+
!  These routines test the sqrt functions *just in case*
!+
!----------------------------------------------------------------
subroutine testsqrt(ierr,output)
 integer, intent(out) :: ierr
 logical, intent(in)  :: output
 real, parameter :: tol = 5.e-7

 ierr = 0
 call checksqrt(4.,tol,ierr,output)
 call checksqrt(2.,tol,ierr,output)
 call checksqrt(16.,tol,ierr,output)
 call checksqrt(3.141592653589,tol,ierr,output)
 call checksqrt(6878.,tol,ierr,output)

 return
end subroutine testsqrt

subroutine checksqrt(x,tol,ierr,output)
 real,    intent(in)    :: x,tol
 integer, intent(inout) :: ierr
 logical, intent(in)    :: output
 real :: fx,fexact,err

 fx = finvsqrt(x)
 fexact = 1./sqrt(x)
 err = abs(fx - fexact)
 if (output) print*,' 1./sqrt(',x,') = ',fx,', error = ',err

 if (err > tol) then
    if (output) print*,'error with fast sqrt on this architecture'
    ierr = ierr + 1
 endif

 return
end subroutine checksqrt

end module fastmath
