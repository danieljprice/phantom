!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: libphantomsplash
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
!  DEPENDENCIES: kernel
!+
!--------------------------------------------------------------------------
module libphantomsplash
 implicit none
contains

subroutine interpolate3D_fastxsec(xyzh,weight,dat,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise)
 use kernel, only:cnormk,radkern,radkern2,wkern
 integer, intent(in)  :: npart,npixx,npixy
 real,    intent(in)  :: xyzh(4,npart),weight,dat(npart)
 real,    intent(in)  :: xmin,ymin,pixwidthx,pixwidthy,zslice
 real,    intent(out) :: datsmooth(npixx,npixy)
 logical, intent(in)  :: normalise
 real :: datnorm(npixx,npixy)

 integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
 real :: hi,hi1,radkern_local,q,q2,wab,const,xi,yi,hi21
 real :: termnorm,term,dy,dy2,dz,dz2,ypix,rescalefac
 real :: dx2i(npixx)

 datsmooth = 0.
 datnorm = 0.
 if (pixwidthx <= 0. .or. pixwidthy <= 0.) then
    print*,'interpolate3D_xsec: error: pixel width <= 0'
    return
 elseif (npart <= 0) then
    print*,'interpolate3D_xsec: error: npart = 0'
    return
 endif
 const = cnormk
 !
 !--renormalise dat array by first element to speed things up
 !
 if (dat(1) > tiny(dat)) then
    rescalefac = dat(1)
 else
    rescalefac = 1.0
 endif

 !
 !--loop over particles
 !
 over_parts: do i=1,npart
    !
    !--set kernel related quantities
    !
    hi = xyzh(4,i)
    if (hi <= 0.) cycle over_parts
    hi1 = 1./hi
    hi21 = hi1*hi1
    radkern_local = radkern*hi    ! radius of the smoothing kernel
    !
    !--for each particle, work out distance from the cross section slice.
    !
    dz = zslice - xyzh(3,i)
    dz2 = dz**2*hi21
    !
    !--if this is < 2h then add the particle's contribution to the pixels
    !  otherwise skip all this and start on the next particle
    !
    if (dz2  <  radkern2) then

       xi = xyzh(1,i)
       yi = xyzh(2,i)
       termnorm = const*weight
       term = termnorm*dat(i)/rescalefac
       !
       !--for each particle work out which pixels it contributes to
       !
       ipixmin = int((xi - radkern_local - xmin)/pixwidthx)
       jpixmin = int((yi - radkern_local - ymin)/pixwidthy)
       ipixmax = int((xi + radkern_local - xmin)/pixwidthx) + 1
       jpixmax = int((yi + radkern_local - ymin)/pixwidthy) + 1

       if (ipixmin < 1) ipixmin = 1 ! make sure they only contribute
       if (jpixmin < 1) jpixmin = 1 ! to pixels in the image
       if (ipixmax > npixx) ipixmax = npixx
       if (jpixmax > npixy) jpixmax = npixy
       !
       !--precalculate an array of dx2 for this particle (optimisation)
       !
       do ipix=ipixmin,ipixmax
          dx2i(ipix) = ((xmin + (ipix-0.5)*pixwidthx - xi)**2)*hi21 + dz2
       enddo
       !
       !--loop over pixels, adding the contribution from this particle
       !
       do jpix = jpixmin,jpixmax
          ypix = ymin + (jpix-0.5)*pixwidthy
          dy = ypix - yi
          dy2 = dy*dy*hi21
          do ipix = ipixmin,ipixmax
             q2 = dx2i(ipix) + dy2
             q = sqrt(q2)
             !
             !--SPH kernel
             !
             if (q2 < radkern2) then
                wab = wkern(q2,q)
                !
                !--calculate data value at this pixel using the summation interpolant
                !
                datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab

             endif

          enddo
       enddo

    endif                  ! if particle within 2h of slice
 enddo over_parts                    ! over particles
 !
 !--normalise dat array
 !
 if (normalise) then
    !--normalise everywhere (required if not using SPH weighting)
    where (datnorm > tiny(datnorm))
       datsmooth = datsmooth/datnorm
    end where
 endif
 datsmooth = datsmooth*rescalefac

 return

end subroutine interpolate3D_fastxsec




end module
