!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: interpolations3D_amr
!
!  DESCRIPTION:
!  Module containing routine for interpolation from PHANTOM data
!  to 3D adaptive mesh
!
!  Requires adaptivemesh.f90 module
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: adaptivemesh
!+
!--------------------------------------------------------------------------

module interpolations3D_amr
 implicit none
 real, parameter, private :: dpi = 1./3.1415926536d0
 public :: interpolate3D_amr
 !$ integer(kind=8), dimension(:), private, allocatable :: ilock

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
!     Input: particle coordinates and h : xyzh(4,npart)
!            weight for each particle   : weight [ same on all parts in PHANTOM ]
!            scalar data to smooth      : dat   (npart)
!
!     Output: smoothed data             : datsmooth (npixx,npixy,npixz)
!
!     Daniel Price, Monash University 2010
!     daniel.price@monash.edu
!--------------------------------------------------------------------------

subroutine interpolate3D_amr(xyzh,weight,pmass,vxyzu,npart, &
                             xmin,datsmooth,nnodes,dxmax,normalise,ilendat,dat)
 use adaptivemesh, only:ifirstlevel,nsub,ndim,gridnodes
 integer,      intent(in)  :: npart,nnodes,ilendat
 real,         intent(in)  :: xyzh(:,:),vxyzu(:,:)
 real,         intent(in)  :: weight,pmass
 real,         intent(in)  :: xmin(3),dxmax(3)
 real,         intent(out) :: datsmooth(4+ilendat,nsub**ndim,nnodes)
 logical,      intent(in)  :: normalise
 real,         intent(in), optional :: dat(:,:)
 real, allocatable :: datnorm(:,:)
!  real, dimension(nsub**ndim,nnodes) :: datnorm

 integer :: i,ipix,jpix,kpix,isubmesh,imesh,level,icell
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,npixx,npixy,npixz
 real :: xi,yi,zi,hi,hi1,hi21,radkern,qq,wab,q2,const,dyz2,dz2
 real :: xorigi,yorigi,zorigi,xpix,ypix,zpix,dx,dy,dz
 real :: dxcell(ndim),xminnew(ndim)
 real :: t_start,t_end
 real :: termnorm
 real :: termdat(4+ilendat)
 logical :: iprintprogress
!$ integer :: omp_get_num_threads,j
#ifndef _OPENMP
 integer(kind=8) :: iprogress
#endif

 datsmooth = 0.
 !datnorm = 0.
 if (normalise) then
    print "(1x,a)",'interpolating from particles to 3D AMR grid (normalised) ...'
 else
    print "(1x,a)",'interpolating from particles to 3D AMR grid (non-normalised) ...'
 endif
 if (any(dxmax(:) <= 0.)) then
    print "(1x,a)",'interpolate3D_amr: error: grid size <= 0'
    return
 endif
 if (present(dat)) then
    if (ilendat /= size(dat(:,1))) then
       print "(1x,a,i2,i2)",'interpolate3D_amr: error in interface: size(dat) /= ilendat ',size(dat(:,1)),ilendat
       return
    endif
 elseif (ilendat /= 0) then
    print "(1x,a)",'interpolate3D_amr: error in interface: dat has non-zero length but is not present'
    return
 endif
 if (normalise) then
    allocate(datnorm(nsub**ndim,nnodes))
    datnorm = 0.
 endif

!$  allocate(ilock(0:nnodes))
!$  do i=0,nnodes
!$     call omp_init_lock(ilock(i))
!$  enddo

 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = (npart  >=  100000) .or. (nnodes > 10000)
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

 imesh      = 1
 level      = ifirstlevel
 dxcell(:)  = dxmax(:)/real(nsub**level)
!  xminpix(:) = xmin(:) - 0.5*dxcell(:)
 npixx = nsub**level
 npixy = nsub**level
 npixz = nsub**level
 print "(3(a,i4))",' root grid: ',npixx,' x ',npixy,' x ',npixz

 const = dpi  ! kernel normalisation constant (3D)
 !
 !--loop over particles
 !
 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,vxyzu,dat,datsmooth,datnorm,gridnodes) &
 !$omp firstprivate(const,ilendat,weight)                       &
 !$omp firstprivate(xmin,pmass,imesh,nnodes,level)               &
 !$omp firstprivate(npixx,npixy,npixz,dxmax,dxcell,normalise)    &
 !$omp private(i,j,hi,hi1,hi21,radkern,termnorm,termdat)         &
 !$omp private(xpix,ypix,zpix,dx,dy,dz,dz2,dyz2,qq,q2,wab)       &
 !$omp private(xi,yi,zi,xorigi,yorigi,zorigi,xminnew)            &
 !$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi,icell,isubmesh)  &
 !$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax)
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
    !--set kernel related quantities
    !
    xi = xyzh(1,i); xorigi = xi
    yi = xyzh(2,i); yorigi = yi
    zi = xyzh(3,i); zorigi = zi
    hi = xyzh(4,i)
    if (hi <= 0.) cycle over_parts
    hi1 = 1./hi; hi21 = hi1*hi1
    termnorm = const*weight

    radkern      = 2.*hi   ! radius of the smoothing kernel
    termdat(1)   = const*pmass*hi21*hi1  ! weight for density calculation
    termdat(2:4) = termnorm*vxyzu(1:3,i)
    if (present(dat) .and. ilendat > 0) then
       termdat(5:5+ilendat-1) = termnorm*dat(1:ilendat,i)
    endif
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((xi - radkern - xmin(1))/dxcell(1))
    jpixmin = int((yi - radkern - xmin(2))/dxcell(2))
    kpixmin = int((zi - radkern - xmin(3))/dxcell(3))

    ipixmax = int((xi + radkern - xmin(1))/dxcell(1)) + 1
    jpixmax = int((yi + radkern - xmin(2))/dxcell(2)) + 1
    kpixmax = int((zi + radkern - xmin(3))/dxcell(3)) + 1

#ifndef PERIODIC
    if (ipixmin < 1)     ipixmin = 1  ! make sure they only contribute
    if (jpixmin < 1)     jpixmin = 1  ! to pixels in the image
    if (kpixmin < 1)     kpixmin = 1
    if (ipixmax > npixx) ipixmax = npixx
    if (jpixmax > npixy) jpixmax = npixy
    if (kpixmax > npixz) kpixmax = npixz
#endif
    !print*,' part ',i,' lims = ',ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
    !
    !--loop over pixels, adding the contribution from this particle
    !  (note that we handle the periodic boundary conditions
    !   entirely on the root grid)
    !
    do kpix = kpixmin,kpixmax
       kpixi = kpix
#ifdef PERIODIC
       if (kpixi < 1) then
          kpixi = kpixi  + npixz
          zi    = zorigi + dxmax(3)
       elseif (kpixi > npixz) then
          kpixi = kpixi  - npixz
          zi    = zorigi - dxmax(3)
       else
          zi    = zorigi
       endif
#endif
       zpix = xmin(3) + (kpixi-0.5)*dxcell(3)
       dz   = zpix - zi
       dz2  = dz*dz*hi21

       do jpix = jpixmin,jpixmax
          jpixi = jpix
#ifdef PERIODIC
          if (jpixi < 1) then
             jpixi = jpixi  + npixy
             yi    = yorigi + dxmax(2)
          elseif (jpixi > npixy) then
             jpixi = jpixi  - npixy
             yi    = yorigi - dxmax(2)
          else
             yi    = yorigi
          endif
#endif
          ypix = xmin(2) + (jpixi-0.5)*dxcell(2)
          dy   = ypix - yi
          dyz2 = dy*dy*hi21 + dz2

          do ipix = ipixmin,ipixmax
             ipixi = ipix
#ifdef PERIODIC
             if (ipixi < 1) then
                ipixi = ipixi  + npixx
                xi    = xorigi + dxmax(1)
             elseif (ipixi > npixx) then
                ipixi = ipixi  - npixx
                xi    = xorigi - dxmax(1)
             else
                xi    = xorigi
             endif
#endif
             icell    = ((kpixi-1)*nsub + (jpixi-1))*nsub + ipixi
             isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)

             if (isubmesh > 0) then
                xminnew(1) = xmin(1) + (ipixi-1)*dxcell(1)
                xminnew(2) = xmin(2) + (jpixi-1)*dxcell(2)
                xminnew(3) = xmin(3) + (kpixi-1)*dxcell(3)
                call interpolate_submesh(xi,yi,zi,radkern,hi21,termnorm,4+ilendat,termdat, &
                         datsmooth,nnodes,datnorm,normalise,isubmesh,level+1,xminnew,dxmax)
             else
                !
                !--particle interpolates directly onto the root grid
                !
                !print*,'onto root grid ',ipixi,jpixi,kpixi
                xpix = xmin(1) + (ipixi-0.5)*dxcell(1)
                dx   = xpix - xi
                q2   = dx*dx*hi21 + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
                !
                !--SPH kernel - standard cubic spline
                !
                if (q2 < 4.0) then
                   if (q2 < 1.0) then
                      qq = sqrt(q2)
                      wab = 1.-1.5*q2 + 0.75*q2*qq
                   else
                      qq = sqrt(q2)
                      wab = 0.25*(2.-qq)**3
                   endif
                   !
                   !--calculate data value at this pixel using the summation interpolant
                   !
                   !$call omp_set_lock(ilock(imesh))
                   datsmooth(:,icell,imesh) = datsmooth(:,icell,imesh) + termdat(:)*wab

                   if (normalise) then
                      datnorm(icell,imesh) = datnorm(icell,imesh) + termnorm*wab
                   endif
                   !$call omp_unset_lock(ilock(imesh))
                endif
             endif
          enddo
       enddo
    enddo
 enddo over_parts
 !$omp enddo
 !$omp end parallel

!$   do i=0,nnodes
!$      call omp_destroy_lock(ilock(i))
!$   enddo
!$   if (allocated(ilock)) deallocate(ilock)

 !
 !--normalise dat array
 !
 if (normalise) then
    do i=2,4+ilendat  !--IN GENERAL BETTER TO NOT NORMALISE THE DENSITY (start from 2)
       where (datnorm > tiny(datnorm))
          datsmooth(i,:,:) = datsmooth(i,:,:)/datnorm(:,:)
       end where
    enddo
 endif
 if (allocated(datnorm)) deallocate(datnorm)
 !
 !--get ending CPU time
 !
 call cpu_time(t_end)
 print*,'completed in ',t_end-t_start,'s'

 return

end subroutine interpolate3D_amr

recursive subroutine interpolate_submesh(xi,yi,zi,radkern,hi21,termnorm,ilendat,termdat,&
                     datsmooth,nnodes,datnorm,normalise,imesh,level,xminl,dxmax)
 use adaptivemesh, only:ndim,nsub,gridnodes
 real,    intent(in)    :: xi,yi,zi,radkern,hi21,termnorm
 integer, intent(in)    :: ilendat,nnodes
 real,    intent(in)    :: termdat(ilendat)
 real,    intent(inout) :: datsmooth(ilendat,nsub**ndim,nnodes)
 real,    intent(inout) :: datnorm(nsub**ndim,nnodes)
 logical, intent(in)    :: normalise
 integer, intent(in)    :: imesh,level
 real,    intent(in)    :: xminl(ndim),dxmax(ndim)
 real :: xminnew(ndim),dxcell(ndim)
 real :: dxpix2(4,nsub)
 real                    :: q2,qq,wab,dxpix,dypix,dzpix,dlevel,qterm
 integer                 :: icell,icellz,icellyz,ipix,jpix,kpix,isubmesh
 integer                 :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax

 dlevel = 1./real(nsub**level)
 dxcell(:) = dxmax(:)*dlevel
 !print*,' submesh, level ',level,' dx = ',dxcell(:)

 ipixmin = int((xi - radkern - xminl(1))/dxcell(1))
 jpixmin = int((yi - radkern - xminl(2))/dxcell(2))
 kpixmin = int((zi - radkern - xminl(3))/dxcell(3))

 ipixmax = int((xi + radkern - xminl(1))/dxcell(1)) + 1
 jpixmax = int((yi + radkern - xminl(2))/dxcell(2)) + 1
 kpixmax = int((zi + radkern - xminl(3))/dxcell(3)) + 1

 if (ipixmin < 1)    ipixmin = 1
 if (ipixmax > nsub) ipixmax = nsub
 if (jpixmin < 1)    jpixmin = 1
 if (jpixmax > nsub) jpixmax = nsub
 if (kpixmin < 1)    kpixmin = 1
 if (kpixmax > nsub) kpixmax = nsub
 !
 !--this is an optimisation: we pre-calculate
 !  the distances in each direction between
 !  the particle and the grid cells
 !  (note that subgrids are very small, so
 !   total array size here should be small)
 !
 do ipix=1,nsub
    dxpix = (xminl(1) + (ipix-0.5)*dxcell(1)) - xi
    dypix = (xminl(2) + (ipix-0.5)*dxcell(2)) - yi
    dzpix = (xminl(3) + (ipix-0.5)*dxcell(3)) - zi
    dxpix2(1,ipix) = dxpix*dxpix*hi21
    dxpix2(2,ipix) = dypix*dypix*hi21
    dxpix2(3,ipix) = dzpix*dzpix*hi21
 enddo

 do kpix=kpixmin,kpixmax
    icellz = (kpix-1)*nsub
    !zpix   = xminl(3) + (kpix-0.5)*dxcell(3)
    !dz     = zpix - zi
    !dz2    = dz*dz*hi21

    do jpix=jpixmin,jpixmax
       !ypix    = xminl(2) + (jpix-0.5)*dxcell(2)
       !dy      = ypix - yi
       !dyz2    = dy*dy*hi21 + dz2
       icellyz = (icellz + (jpix-1))*nsub

       over_x: do ipix=ipixmin,ipixmax
          icell = icellyz + ipix
          !print*,'level ',level,' icell = ',icell

          isubmesh = gridnodes(icell,imesh) !grid(imesh)%daughter(icell)
          if (isubmesh > 0) then
             !print*,' splitting...'
             xminnew(1) = xminl(1) + (ipix-1)*dxcell(1)
             xminnew(2) = xminl(2) + (jpix-1)*dxcell(2)
             xminnew(3) = xminl(3) + (kpix-1)*dxcell(3)
             call interpolate_submesh(xi,yi,zi,radkern,hi21,termnorm,ilendat,termdat, &
                                      datsmooth,nnodes,datnorm,normalise,&
                                      isubmesh,level+1,xminnew,dxmax)
          else
             !print*,' interpolating...'
             !if (dyz2 > 4.) cycle over_x

             !xpix = xminl(1) + (ipix-0.5)*dxcell(1)
             !dx   = xpix - xi
             !dx2  = dx*dx*hi21
             !q2   = dx2 + dyz2

             !--optimisation: use pre-calculated distances
             q2 = dxpix2(1,ipix) + dxpix2(2,jpix) + dxpix2(3,kpix)

             if (q2 < 4.) then
                !
                !--SPH kernel - standard cubic spline
                !
                if (q2 < 4.0) then
                   if (q2 < 1.0) then
                      qq = sqrt(q2)
                      wab = 1.-1.5*q2 + 0.75*q2*qq
                   else
                      qq = sqrt(q2)
                      qterm = 2.-qq
                      wab = 0.25*qterm*qterm*qterm
                   endif
                endif
                !
                !--calculate data value at this pixel using the summation interpolant
                !
                !$ call omp_set_lock(ilock(imesh))
                datsmooth(:,icell,imesh) = datsmooth(:,icell,imesh) + termdat(:)*wab

                if (normalise) then
                   datnorm(icell,imesh) = datnorm(icell,imesh) + termnorm*wab
                endif
                !$ call omp_unset_lock(ilock(imesh))
             endif
          endif
       enddo over_x
    enddo
 enddo

 return
end subroutine interpolate_submesh

end module interpolations3D_amr
