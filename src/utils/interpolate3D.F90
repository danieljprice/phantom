!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: interpolations3D
!
!  DESCRIPTION:
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

module interpolations3D
 implicit none
 real, parameter, private :: dpi = 1./3.1415926536d0
 public :: interpolate3D

contains
!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is interpolated according to the formula
!
!     datsmooth(pixel) = sum_b weight_b dat_b W(r-r_b, h_b)
!
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline.
!
!     For a standard SPH smoothing the weight function for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     this version is written for slices through a rectangular volume, ie.
!     assumes a uniform pixel size in x,y, whilst the number of pixels
!     in the z direction can be set to the number of cross-section slices.
!
!     Input: particle coordinates  : x,y,z (npart)
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Institute of Astronomy, Cambridge 16/7/03
!--------------------------------------------------------------------------

subroutine interpolate3D(xyzh,weight,pmass,vxyzu,npart, &
     xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidth,normalise,ilendat4,dat4)
 integer,      intent(in)  :: npart,npixx,npixy,npixz,ilendat4
 real,         intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real,         intent(in)  :: weight,pmass,xmin,ymin,zmin,pixwidth
 real,         intent(out) :: datsmooth(4+ilendat4,npixx,npixy,npixz) !rho, vx, vy, vz
 real(kind=4), intent(in),              optional :: dat4(:,:)
 logical,      intent(in)  :: normalise
 real :: datnorm(npixx,npixy,npixz)

 integer :: i,ipix,jpix,kpix
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nsubgrid,nfull
 real :: xminpix,yminpix,zminpix,hmin,hminall !,dhmin3
 real :: xpix(npixx),dx2i(npixx)
 real :: xi,yi,zi,hi,hi1,hi21,radkern,qq,wab,q2,const,dyz2,dz2
 real :: t_start,t_end
 real :: term1,termnorm,dy,dz,ypix,zpix,xpixi
 real :: termdat(4+ilendat4)
 logical :: iprintprogress
!$ integer :: omp_get_num_threads,j
#ifndef _OPENMP
 integer(kind=8) :: iprogress
#endif

 datsmooth = 0.
 datnorm = 0.
 if (normalise) then
    print "(1x,a,i4,a,i4,a,i4,a)",'interpolating from particles to ',npixx,' x',npixy,' x',npixz,' 3D grid (normalised) ...'
 else
    print "(1x,a,i4,a,i4,a,i4,a)",'interpolating from particles to ',npixx,' x',npixy,' x',npixz,' 3D grid (non-normalised) ...'
 endif
 if (pixwidth <= 0.) then
    print "(1x,a)",'interpolate3D: error: pixel width <= 0'
    return
 endif
 if (present(dat4)) then
    if (ilendat4 /= size(dat4(:,1))) then
       print "(1x,a)",'interpolate3D: error in interface: size(dat4) /= ilendat4'
       return
    endif
 elseif (ilendat4 /= 0) then
    print "(1x,a)",'interpolate3D: error in interface: dat4 has non-zero length but is not present'
    return
 endif
 !if (any(hh(1:npart) <= tiny(hh))) then
 !   print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
 !endif

 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000)
 !
 !--loop over particles
 !
 iprintinterval = 25
 if (npart >= 1e6) iprintinterval = 10
 iprintnext = iprintinterval
 !
 !--get starting CPU time
 !
 call cpu_time(t_start)

 xminpix = xmin - 0.5*pixwidth
 yminpix = ymin - 0.5*pixwidth
 zminpix = zmin - 0.5*pixwidth
 !zminpix = zmin - 0.5*zpixwidth
!  xmax = xmin + npixx*pixwidth
!  ymax = ymin + npixy*pixwidth
 !
 !--use a minimum smoothing length on the grid to make
 !  sure that particles contribute to at least one pixel
 !
 hmin     = 0.5*pixwidth
 !dhmin3   = 1./(hmin*hmin*hmin)
 hminall  = huge(hminall)
 nsubgrid = 0
 !
 !--store x value for each pixel (for optimisation)
 !
 do ipix=1,npixx
    xpix(ipix) = xminpix + ipix*pixwidth
 enddo
 const = dpi  ! normalisation constant (3D)
 !
 !--loop over particles
 !
 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,vxyzu,dat4,datsmooth,datnorm)            &
 !$omp firstprivate(xpix,hmin,const,ilendat4,weight)              &
 !$omp firstprivate(xmin,ymin,zmin,pmass,xminpix,yminpix,zminpix) &
 !$omp firstprivate(npixx,npixy,npixz,pixwidth,normalise)         &
 !$omp private(i,j,hi,hi1,hi21,radkern,termnorm,termdat,term1)    &
 !$omp private(xpixi,dx2i,ypix,zpix,dy,dz,dz2,dyz2,qq,q2,wab)     &
 !$omp private(xi,yi,zi,ipix,jpix,kpix,ipixi,jpixi,kpixi,nxpix)   &
 !$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax)   &
 !$omp reduction(+:nsubgrid) reduction(min:hminall)
 !$omp master
 !$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
 !$omp end master
 !$omp do schedule(guided,10)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
#ifndef _OPENMP
    if (iprintprogress) then
       iprogress = nint(100.*i/npart)
       if (iprogress >= iprintnext) then
          write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
          iprintnext = iprintnext + iprintinterval
       endif
    endif
#endif
    !
    !--skip particles with itype < 0
    !
    !if (itype(i) < 0) cycle over_parts
    !hi = hh(i)
    !if (hi <= 0.) cycle over_parts

    !
    !--set kernel related quantities
    !
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (hi < hmin) then
       !--best results are found if we set h= hmin
       !  but do *not* adjust the weight
       hminall  = min(hi,hminall)
       term1    = const*pmass/(hi*hi*hi)
       termnorm = const*weight !*(hi*hi*hi)*dhmin3
       hi = hmin
       hi1 = 1./hi
       hi21 = hi1*hi1
       nsubgrid = nsubgrid + 1
    else
       hi1 = 1./hi
       hi21 = hi1*hi1
       term1    = const*pmass*hi21*hi1
       termnorm = const*weight
    endif

    radkern = 2.*hi   ! radius of the smoothing kernel
    termdat(1) = term1  ! weight for density calculation
    termdat(2:4) = termnorm*vxyzu(1:3,i)
    if (present(dat4) .and. ilendat4 > 0) then
       termdat(5:5+ilendat4-1) = termnorm*dat4(1:ilendat4,i)
    endif
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((xi - radkern - xmin)/pixwidth)
    jpixmin = int((yi - radkern - ymin)/pixwidth)
    kpixmin = int((zi - radkern - zmin)/pixwidth)
    ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
    jpixmax = int((yi + radkern - ymin)/pixwidth) + 1
    kpixmax = int((zi + radkern - zmin)/pixwidth) + 1

#ifndef PERIODIC
    if (ipixmin < 1) ipixmin = 1  ! make sure they only contribute
    if (jpixmin < 1) jpixmin = 1  ! to pixels in the image
    if (kpixmin < 1) kpixmin = 1
    if (ipixmax > npixx) ipixmax = npixx
    if (jpixmax > npixy) jpixmax = npixy
    if (kpixmax > npixz) kpixmax = npixz
#endif
    !
    !--precalculate an array of dx2 for this particle (optimisation)
    !
    nxpix = 0
    do ipix=ipixmin,ipixmax
       nxpix = nxpix + 1
       ipixi = ipix
#ifdef PERIODIC
       if (ipixi < 1) ipixi = ipixi + npixx
       if (ipixi > npixx) ipixi = ipixi - npixx
#endif
       xpixi = xminpix + ipix*pixwidth

       dx2i(nxpix) = ((xpixi - xi)**2)*hi21
    enddo
    !
    !--loop over pixels, adding the contribution from this particle
    !
    do kpix = kpixmin,kpixmax
       kpixi = kpix
#ifdef PERIODIC
       if (kpixi < 1) kpixi = kpixi + npixz
       if (kpixi > npixz) kpixi = kpixi - npixz
#endif
       zpix = zminpix + kpix*pixwidth
       dz = zpix - zi
       dz2 = dz*dz*hi21

       do jpix = jpixmin,jpixmax
          jpixi = jpix
#ifdef PERIODIC
          if (jpixi < 1) jpixi = jpixi + npixy
          if (jpixi > npixy) jpixi = jpixi - npixy
#endif
          ypix = yminpix + jpix*pixwidth
          dy = ypix - yi
          dyz2 = dy*dy*hi21 + dz2

          nxpix = 0
          do ipix = ipixmin,ipixmax
             ipixi = ipix
#ifdef PERIODIC
             if (ipixi < 1) ipixi = ipixi + npixx
             if (ipixi > npixx) ipixi = ipixi - npixx
#endif
             nxpix = nxpix + 1
             q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
             !
             !--SPH kernel - standard cubic spline
             !
             if (q2 < 4.0) then
                if (q2 < 1.0) then
                   qq = sqrt(q2)
                   wab = 1.-1.5*q2 + 0.75*q2*qq
!
!-- the following lines use a fast inverse sqrt function
!
!                    if (q2 > epsilon(q2)) then
!                       qq = q2*finvsqrt(q2)
!                       wab = 1.-1.5*q2 + 0.75*q2*qq
!                    else
!                       wab = 1.
!                    endif
                else
                   qq = sqrt(q2)
!                    qq = q2*finvsqrt(q2)
                   wab = 0.25*(2.-qq)**3
                endif
                !
                !--calculate data value at this pixel using the summation interpolant
                !
#ifdef _OPENMP
                do j=1,4+ilendat4
                   !$omp atomic
                   datsmooth(j,ipixi,jpixi,kpixi) = datsmooth(j,ipixi,jpixi,kpixi) + termdat(j)*wab
                enddo
#else
                datsmooth(:,ipixi,jpixi,kpixi) = datsmooth(:,ipixi,jpixi,kpixi) + termdat(:)*wab
#endif
                if (normalise) then
                   !$omp atomic
                   datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                endif
             endif
          enddo
       enddo
    enddo
 enddo over_parts
 !$omp enddo
 !$omp end parallel
 !
 !--normalise dat array
 !
 if (normalise) then
    do i=2,4+ilendat4  !--IN GENERAL BETTER TO NOT NORMALISE THE DENSITY (start from 2)
       where (datnorm > tiny(datnorm))
          datsmooth(i,:,:,:) = datsmooth(i,:,:,:)/datnorm(:,:,:)
       end where
    enddo
 endif
 !
 !--warn about subgrid interpolation
 !
 if (nsubgrid > 1) then
    nfull = int((npixx*pixwidth)/hminall) + 1
    print "(a,i9,a,/,a,i8,a)",' Warning: pixel size > 2h for ',nsubgrid,' particles', &
                               '          need',nfull,' pixels for full resolution'
 endif
 !
 !--get ending CPU time
 !
 call cpu_time(t_end)
 print*,' completed in ',t_end-t_start,'s'

 return

end subroutine interpolate3D

end module interpolations3D
