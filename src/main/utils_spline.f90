!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: splineutils
!
!  DESCRIPTION:
!  This module contains utilities for fitting/ evaluating
!  cubic spline fits to data
!
!  Routines are originally by Simon Glover,
!  Translated to Fortran 90 and adapted
!  for use in Phantom by Daniel Price (2011)
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
module splineutils
 implicit none

 public :: spline_eval

 private

contains

!----------------------------------------------------------------
!+
!  Written by: Simon Glover, January 2004
!
!  PURPOSE: Compute a spline fit to a set of nval data points.
!  Evaluate this fit at nnew specified locations and return the
!  nnew computed values.
!+
!---------------------------------------------------------------
subroutine spline_eval(nval, positions, values, nnew, new_positions, new_values)
 !use cooling  ! surely not necessary???
 integer, intent(in)  :: nval, nnew
 real,    intent(in)  :: positions(nval),values(nval)
 real,    intent(in)  :: new_positions(nnew)
 real,    intent(out) :: new_values(nnew)
 real    :: coefficients(4,nval-1)
 real    :: dpt, pt, a, b, c, d
 integer :: i, j, index

 call spline_coefficients(nval, values, coefficients)

 do i = 1, nnew
    pt = new_positions(i)

! We're using -1 to indicate that we've been asked for a value outside of
! the range of our data, and setting the value for that point to zero;
! anything fancier should probably be handled in the caller.
    if (pt  <  positions(1) .or. pt  >  positions(nval)) then
       index         = -1
       new_values(i) = 0d0
    elseif (pt == positions(nval)) then
! Another special case, but this one is easy to deal with
       index         = -1
       new_values(i) = values(nval)
    else
       index = -1
       dpt   = 0.
       do j = 1, nval-1
          if (pt  >=  positions(j) .and. pt  <  positions(j+1)) then
             index = j
             dpt = (pt - positions(index)) / (positions(index+1) - positions(index))
          endif
       enddo
    endif

    if (index  /=  -1) then
       a = coefficients(1,index)
       b = coefficients(2,index)
       c = coefficients(3,index)
       d = coefficients(4,index)
       new_values(i) = a + b * dpt + c * dpt**2 + d * dpt**3
    endif

 enddo
end subroutine spline_eval

!---------------------------------------------------------------
!+
!  Written by: Simon Glover, January 2004
!
!  PURPOSE: Compute coefficients for a cubic spline fit to the set
!  of nval datapoints specified in values.
!
!  The solution follows the procedure of Bartels et al., 1998, "An
!  Introduction to Splines for Use in Computer Graphics and Geometric
!  Modelling", Morgan Kaufmann, ch. 3, pp 9-17, as outlined at
!  http://mathworld.wolfram.com/CubicSpline.html
!
!  NB This routine does little or no error-checking, since it is only
!  intended to be used with a few small sets of known good data.
!+
!---------------------------------------------------------------
subroutine spline_coefficients(nval, values, coefficients)
 integer, intent(in)  :: nval
 real,    intent(in)  :: values(nval)
 real,    intent(out) :: coefficients(4,nval-1)
 real :: derivatives(nval)
 integer :: i

! We first calculate the first derivative of the curve at each of
! our data points
 call spline_derivatives(nval, values, derivatives)

! The spline coefficients can then be expressed in terms of the derivatives
! and the values at the data points
 do i = 1, nval-1
    coefficients(1,i) = values(i)
    coefficients(2,i) = derivatives(i)
    coefficients(3,i) = 3d0 * (values(i+1) - values(i)) - &
                        2d0 * derivatives(i) - derivatives(i+1)
    coefficients(4,i) = 2d0 * (values(i) - values(i+1)) +  &
                        derivatives(i) + derivatives(i+1)
 enddo

end subroutine spline_coefficients

!---------------------------------------------------------------
!+
!  Written by: Simon Glover, January 2004
!
!  PURPOSE: Compute the derivates of the spline curve required by
!  the spline_coefficient subroutine. For more details, see the
!  preamble to that subroutine.
!+
!---------------------------------------------------------------
!
subroutine spline_derivatives(nval, values, derivatives)
 integer, intent(in)  :: nval
 real,    intent(in)  :: values(nval)
 real,    intent(out) :: derivatives(nval)
 real :: Cc(nval,nval+1)
 real    :: f,sum
 integer :: i, j

! The derivatives we require are given by the solution of the matrix equation
! Ax = B, where:
!
!      |2 1        |       |D_1  |           |3*(y_1 - y_0)  |
!      |1 4 1      |       |D_2  |           |3*(y_2 - y_0)  |
!  A = |  .....    |,  x = |...  |  and  B = |3*(y_3 - y_1)  |
!      |    .....  |       |...  |           | ............. |
!      |      1 4 1|       |D_n-1|           |3*(y_n - y_n-2)|
!      |        1 2|       |D_n  |           |3*(y_n - y_n-1)|
!
! To solve this equation, we first construct the 'augmented', n by n+1 matrix
! Cc = AB:

 do i = 1, nval
    do j = 1, nval+1
       Cc(I,J) = 0d0
    enddo
    if (i == 1) then
       Cc(i,i)      = 2d0
       Cc(i,i+1)    = 1d0
       Cc(i,nval+1) = 3d0 * (values(2) - values(1))
    elseif (i == nval) then
       Cc(i,i-1)    = 1d0
       Cc(i,i)      = 2d0
       Cc(i,nval+1) = 3d0 * (values(nval) - values(nval-1))
    else
       Cc(i,i-1)    = 1d0
       Cc(i,i)      = 4d0
       Cc(i,i+1)    = 1d0
       Cc(i,nval+1) = 3d0 * (values(i+1) - values(i-1))
    endif
 enddo

! We then reduce this matrix to upper triangular form by Gaussian
! elimination...
 do i = 2, nval
    f = Cc(i,i-1) / Cc(i-1,i-1)
    do j = 1, nval+1
       Cc(i,j) = Cc(i,j) - f * Cc(i-1,j)
    enddo
 enddo

! ...and finally solve for the derivatives by back substitution
 derivatives(nval) = Cc(nval,nval+1) / Cc(nval,nval)

 do i = nval-1, 1, -1
    sum = 0d0
    do j = i+1, nval
       sum = sum + Cc(i,j) * derivatives(j)
    enddo
    derivatives(i) = (Cc(i,nval+1) - sum) / Cc(i,i)
 enddo

end subroutine spline_derivatives

end module splineutils
