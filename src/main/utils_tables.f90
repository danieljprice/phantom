!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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

 private

contains

!----------------------------------------------------
!+
!  Function to linearly interpolate values from
!  from an array
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
!+
!--------------------------------------------------------
pure subroutine linspace(x,xmin,xmax)
 real, intent(out) :: x(:)
 real, intent(in)  :: xmin,xmax
 integer :: i, n
 real    :: dx

 n  = size(x)
 dx = (xmax - xmin)/real(n-1)
 do i=1,n
    x(i) = xmin + (i-1)*dx
 enddo

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
subroutine interpolator(array, value, valueidx)
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

end module table_utils
