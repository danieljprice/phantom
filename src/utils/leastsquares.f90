!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: leastsquares
!
!  DESCRIPTION:
!  module containing routine(s)
!  to compute least squares fit of a linear function
!  for arbitrary data, optionally in log-space
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
module leastsquares
 implicit none

contains

subroutine fit_slope(nvals,xval,yval,slope,yint,err,errslope,erryint,xmin,xmax,logplot,fixslope)
 integer,      intent(in)  :: nvals
 real,         intent(in)  :: xval(nvals)
 real(kind=8), intent(in)  :: yval(nvals)
 real,         intent(out) :: slope,yint,err,errslope,erryint
 real,         intent(in), optional :: xmin,xmax
 logical,      intent(in), optional :: logplot,fixslope
 real(kind=8) :: sumx,sumy,sumxx,sumxy,sumyy
 real(kind=8) :: xvali,xi,yi,xmean,ymean
 real(kind=8) :: xmini,xmaxi,erri,sfac
 logical :: islog,fix_slope
 integer :: i,npts
!
!--set x limits if present
!
 if (present(xmin)) then
    xmini = xmin
 else
    xmini = -huge(xmini)
 endif
 if (present(xmax)) then
    xmaxi = xmax
 else
    xmaxi = huge(xmaxi)
 endif
 islog = .false.
 if (present(logplot)) islog = logplot
 fix_slope = .false.
 if (present(fixslope)) fix_slope = fixslope
!
!--least squares fit of a straight line
!
 npts = 0
 sumx = 0.
 sumy = 0.
 sumxy = 0.
 sumxx = 0.
 sumyy = 0.
 do i=1,nvals
    xvali = xval(i)
    if (xvali > xmini .and. xvali < xmaxi) then
       npts = npts + 1
       if (islog) then
          xi = log10(xvali)
          yi = log10(yval(i))
       else
          xi = xvali
          yi = yval(i)
       endif
       sumx = sumx + xi
       sumy = sumy + yi
       sumxy = sumxy + xi*yi
       sumxx = sumxx + xi*xi
       sumyy = sumyy + yi*yi
    endif
 enddo

 if (npts > 0 .and. sumxx > 0. .and. sumyy > 0.) then
    if (.not.fix_slope) slope = (npts*sumxy - sumy*sumx)/(npts*sumxx - sumx*sumx)
    yint  = (sumy - slope*sumx)/dble(npts)
!
!--calculate errors
!  see http://mathworld.wolfram.com/LeastSquaresFitting.html
!
!  err      is the correlation coefficient
!  errslope is the standard error in the slope
!  erryint  is the standard error in the y intercept
!
    xmean = sumx/dble(npts)
    ymean = sumy/dble(npts)
    sumxx = sumxx - npts*xmean*xmean
    sumxy = sumxy - npts*xmean*ymean
    sumyy = sumyy - npts*ymean*ymean

    err   = (sumxy*sumxy)/(sumxx*sumyy)
    sfac  = sqrt((sumyy - sumxy*sumxy/sumxx)/dble(npts - 2))
    if (.not.fix_slope) errslope = sfac/sqrt(sumxx)
    erryint  = sfac*sqrt(1./dble(npts) + xmean**2/sumxx)

    if (present(xmin)) then
       xmini = xmin
       do i=1,nvals
          xvali = xval(i)
          if (xvali < xmin) then
             if (islog) then
                xi = log10(xvali)
                yi = log10(yval(i))
             else
                xi = xvali
                yi = yval(i)
             endif
             erri = sqrt((yi - (slope*xi + yint))**2)
             if (erri  <  0.1) then
                !print*,' xi = ',xi,' err = ',erri
                xmini = xvali
             endif
          endif
       enddo
       if (xmini < xmin .and. .not.fix_slope) then
          print*,' slope within error down to x = ',xmini
       endif
    endif

 else
    slope = 0.
    yint  = 0.
    err   = -1.
    errslope = -1.
    erryint  = -1.
 endif

end subroutine fit_slope

end module leastsquares
