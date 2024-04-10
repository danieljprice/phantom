!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module table_utils
!
! This module contains low-level utilities related to tabulated data
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

 public :: yinterp, linspace, logspace, diff, flip_array, interpolator
 public :: find_nearest_index, interp_1d

 private

contains

!----------------------------------------------------
!+
!  Function to linearly interpolate values from
!  from an array. Given x it finds the position
!  in the table of x values (xtab) and interpolates
!  the value of y from the table of y values (ytab)
!  to return the y value corresponding to the
!  position x
!+
!----------------------------------------------------
pure real function yinterp(ytab,xtab,x)
 real, intent(in) :: ytab(:),xtab(:)
 real, intent(in) :: x
 integer :: i,n
 real    :: xprev,xi,dx,dydx

 n = size(xtab)

 yinterp = 0.
 if (n <= 0) return

 if (x <= xtab(1)) then
    yinterp = ytab(1)
    return
 elseif (x >= xtab(n)) then
    yinterp = ytab(n)
    return
 endif

 xprev = xtab(1)
 over_tab: do i=2,n
    xi = xtab(i)
    if (x >= xprev .and. x < xi) then
       dx      = x - xprev
       dydx    = (ytab(i) - ytab(i-1))/(xi - xprev)
       yinterp = ytab(i-1) + dx*dydx
       exit over_tab
    endif
    xprev = xi
 enddo over_tab

end function yinterp

!--------------------------------------------------------
!+
!  Function to fill an array with equally spaced points
!  e.g. call linspace(rgrid,rmin,rmax)
!
!  the argument dx is optional giving the grid spacing
!+
!--------------------------------------------------------
pure subroutine linspace(x,xmin,xmax,dx)
 real, intent(out) :: x(:)
 real, intent(in)  :: xmin,xmax
 real, intent(out), optional :: dx
 integer :: i, n
 real    :: dxi

 n  = size(x)
 dxi = (xmax - xmin)/real(n-1)
 do i=1,n
    x(i) = xmin + (i-1)*dxi
 enddo
 if (present(dx)) dx = dxi

end subroutine linspace

!--------------------------------------------------------
!+
!  Function to fill an array with equally log-spaced points
!+
!--------------------------------------------------------
pure subroutine logspace(x,xmin,xmax)
 real, intent(out) :: x(:)
 real, intent(in)  :: xmin,xmax
 integer :: i, n
 real    :: dx

 n  = size(x)

 dx = log10(xmax/xmin)/real(n-1)
 do i=1,n
    x(i) = log10(xmin) + (i-1)*dx
 enddo

 x = 10.**x

end subroutine logspace

!----------------------------------------------------------------
!+
!  Finds index of the array value closest to a given value for an
!  ordered array
!+
!----------------------------------------------------------------
subroutine interpolator(array,value,valueidx)
 real, intent(in)     :: array(:)
 real, intent(in)     :: value
 integer, intent(out) :: valueidx

 valueidx = minloc(abs(array - value), dim = 1)

end subroutine interpolator

!----------------------------------------------------------------
!+
!  Reverses the elements of a 1-d array
!+
!----------------------------------------------------------------
subroutine flip_array(array)
 real, intent(inout) :: array(:)
 real, allocatable   :: flipped_array(:)
 integer             :: i

 allocate(flipped_array(size(array)))
 do i = 1, size(array)
    flipped_array(i) = array(size(array) - i + 1)
 enddo
 array = flipped_array

end subroutine flip_array

!----------------------------------------------------------------
!+
!  Takes a n-dim array and produces a (n-1)-dim array with the
!  ith element being the (i+1)th element minus the ith element of
!  the original array
!+
!----------------------------------------------------------------
subroutine diff(array, darray)
 real, intent(in)               :: array(:)
 real, allocatable, intent(out) :: darray(:)
 integer                        :: i

 allocate(darray(size(array)-1))
 do i = 1, size(array)-1
    darray(i) = array(i+1) - array(i)
 enddo

end subroutine diff

!-----------------------------------------------------------------------
!+
!  Find index of nearest lower value in array
!+
!-----------------------------------------------------------------------
subroutine find_nearest_index(arr,val,indx)
 real, intent(in)     :: arr(:), val
 integer, intent(out) :: indx
 integer              :: istart,istop,i

 istart = 1
 istop  = size(arr)
 indx = istart
 if (val >= arr(istop)) then
    indx = istop-1   ! -1 to avoid array index overflow
 elseif (val <= arr(istart)) then
    indx = istart
 else
    i = istart
    do while (i <= istop)
       if (arr(i)>val) then
          indx = i-1
          exit
       endif
       i = i+1
    enddo
 endif

end subroutine find_nearest_index

!-----------------------------------------------------------------------
!+
!  1D linear interpolation routine
!+
!-----------------------------------------------------------------------
real function interp_1d(x,x1,x2,y1,y2)
 real, intent(in)  :: x, x1, x2, y1, y2

 interp_1d = y1 + (x-x1)*(y2-y1)/(x2-x1)

end function interp_1d

end module table_utils
