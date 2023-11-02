!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module interpolations3D
!
! Module containing routine for interpolation from PHANTOM data
!  to 3D adaptive mesh
!
!  Requires adaptivemesh.f90 module
!
! :References: None
!
! :Owner: Spencer Magnall
!
! :Runtime parameters: None
!
! :Dependencies: kernel
!

 implicit none
 real, parameter, private :: dpi = 1./3.1415926536d0
 public :: interpolate3D
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

subroutine interpolate3D(xyzh,weight,npart, &
                             xmin,datsmooth,nnodes,dxgrid,normalise,dat,ngrid,vertexcen)
 use kernel,        only:wkern, radkern, radkern2, cnormk
 !use adaptivemesh, only:ifirstlevel,nsub,ndim,gridnodes
 integer,      intent(in)  :: npart,nnodes,ngrid(3)
 real,         intent(in)  :: xyzh(:,:)! ,vxyzu(:,:)
 real,         intent(in)  :: weight !,pmass
 real,         intent(in)  :: xmin(3),dxgrid(3)
 real,         intent(out) :: datsmooth(:,:,:)
 logical,      intent(in)  :: normalise, vertexcen
 real,         intent(in), optional :: dat(:)
 real, allocatable :: datnorm(:,:,:)
!  real, dimension(nsub**ndim,nnodes) :: datnorm
 integer, parameter :: ndim = 3, nsub=1
 integer :: i,ipix,jpix,kpix,isubmesh,imesh,level,icell
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,npixx,npixy,npixz
 real :: xi,yi,zi,hi,hi1,hi21,radkernh,qq,wab,q2,const,dyz2,dz2
 real :: xorigi,yorigi,zorigi,xpix,ypix,zpix,dx,dy,dz
 real :: dxcell(ndim),xminnew(ndim), dxmax(ndim)
 real :: t_start,t_end
 real :: termnorm
 real :: term
 logical :: iprintprogress
!$ integer :: omp_get_num_threads,j
#ifndef _OPENMP
 integer(kind=8) :: iprogress
#endif

 print*, "size: ", size(datsmooth)
 print*, "datsmooth out of bounds: ", datsmooth(35,1,1)
 datsmooth = 0.
 dxmax(:) = dxgrid(:)
 !datnorm = 0.
 if (normalise) then
    print "(1x,a)",'interpolating from particles to Einstein toolkit grid (normalised) ...'
 else
    print "(1x,a)",'interpolating from particles to Einstein toolkit grid (non-normalised) ...'
 endif
!  if (any(dxmax(:) <= 0.)) then
!     print "(1x,a)",'interpolate3D: error: grid size <= 0'
!     return
!  endif
!  if (ilendat /= 0) then
!     print "(1x,a)",'interpolate3D: error in interface: dat has non-zero length but is not present'
!     return
!  endif
 if (normalise) then
    allocate(datnorm(ngrid(1),ngrid(2),ngrid(3)))
    datnorm = 0.
 endif

!$ allocate(ilock(0:nnodes))
!$ do i=0,nnodes
!$  call omp_init_lock(ilock(i))
!$ enddo

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
 level      = 1
 dxcell(:)  = dxgrid(:)/real(nsub**level)
!  xminpix(:) = xmin(:) - 0.5*dxcell(:)
 npixx = ngrid(1)
 npixy = ngrid(2)
 npixz = ngrid(3)
 print "(3(a,i4))",' root grid: ',npixx,' x ',npixy,' x ',npixz
 print*, "position of i cell  is: ", 1*dxcell(1) + xmin(1)
 print*, "npart: ", npart

 const = cnormk  ! kernel normalisation constant (3D)
 print*,"const: ", const
 !stop

 !
 !--loop over particles
 !
 !$omp parallel default(none) &
 !$omp shared(npart,xyzh,dat,datsmooth,datnorm,vertexcen,const,weight)  &
 !$omp shared(xmin,imesh,nnodes,level)               &
 !$omp shared(npixx,npixy,npixz,dxmax,dxcell,normalise)    &
 !$omp private(i,j,hi,hi1,hi21,termnorm,term)    &
 !$omp private(xpix,ypix,zpix,dx,dy,dz,dz2,dyz2,qq,q2,wab,radkernh)       &
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
    !print*, i
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
    radkernh = radkern*hi
    !print*, "hi: ", hi
    if (hi <= 0.) cycle over_parts
    hi1 = 1./hi; hi21 = hi1*hi1
    termnorm = const*weight
    !  print*, "const: ", const
    !  print*, "weight: ", weight
    !  print*, "termnorm: ", termnorm

    !radkern      = 2.*hi   ! radius of the smoothing kernel
    !print*, "radkern: ", radkern
    !print*, "part pos: ", xi,yi,zi
    term   = termnorm*dat(i)  ! weight for density calculation
    ! I don't understand why this doesnt involve any actual smoothing?
    !dfac = hi**3/(dxcell(1)*dxcell(2)*dxcell(3)*const)
    !
    !--for each particle work out which pixels it contributes to
    !
    !print*, "radkern: ", radkern
    ipixmin = int((xi - radkernh - xmin(1))/dxcell(1))
    jpixmin = int((yi - radkernh - xmin(2))/dxcell(2))
    kpixmin = int((zi - radkernh - xmin(3))/dxcell(3))

    ipixmax = int((xi + radkernh - xmin(1))/dxcell(1)) + 1
    jpixmax = int((yi + radkernh - xmin(2))/dxcell(2)) + 1
    kpixmax = nint((zi + radkernh - xmin(3))/dxcell(3)) + 1

    !if (ipixmax == 33) stop


    !if (ipixmin == 4 .and. jpixmin == 30 .and. kpixmin == 33) print*, "particle (min): ", i
    !if (ipixmax == 4 .and. jpixmax == 30 .and. kpixmax == 33)  print*, "particle (max): ", i
#ifndef PERIODIC
    if (ipixmin < 1)     ipixmin = 1  ! make sure they only contribute
    if (jpixmin < 1)     jpixmin = 1  ! to pixels in the image
    if (kpixmin < 1)     kpixmin = 1
    if (ipixmax > npixx) ipixmax = npixx
    if (jpixmax > npixy) jpixmax = npixy
    if (kpixmax > npixz) kpixmax = npixz
    !print*, "ipixmin: ", ipixmin
    !print*, "ipixmax: ", ipixmax
    !print*, "jpixmin: ", jpixmin
    !print*, "jpixmax: ", jpixmax
    !print*, "kpixmin: ", kpixmin
    !print*, "kpixmax: ", kpixmax
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
          zi    = zorigi !+ dxmax(3)
       elseif (kpixi > npixz) then
          kpixi = kpixi  - npixz
          zi    = zorigi !- dxmax(3)
       else
          zi    = zorigi
       endif
#endif
       if (vertexcen) then
          zpix = xmin(3) + (kpixi-1)*dxcell(3)
       else
          zpix = xmin(3) + (kpixi-0.5)*dxcell(3)
       endif
       dz   = zpix - zi
       dz2  = dz*dz*hi21

       do jpix = jpixmin,jpixmax
          jpixi = jpix
#ifdef PERIODIC
          if (jpixi < 1) then
             jpixi = jpixi  + npixy
             yi    = yorigi !+ dxmax(2)
          elseif (jpixi > npixy) then
             jpixi = jpixi  - npixy
             yi    = yorigi !- dxmax(2)
          else
             yi    = yorigi
          endif
#endif
          if (vertexcen) then
             ypix = xmin(2) + (jpixi-1)*dxcell(2)
          else
             ypix = xmin(2) + (jpixi-0.5)*dxcell(2)
          endif
          dy   = ypix - yi
          dyz2 = dy*dy*hi21 + dz2

          do ipix = ipixmin,ipixmax
             ipixi = ipix
#ifdef PERIODIC
             if (ipixi < 1) then
                ipixi = ipixi  + npixx
                xi    = xorigi !+ dxmax(1)
             elseif (ipixi > npixx) then
                if (ipixi == 33) then
                   print*,"xi old: ", xorigi
                   print*, "xi new: ", xorigi-dxmax(1)
                   print*, "ipixi new: ", ipixi - npixx
                endif
                ipixi = ipixi  - npixx
                xi    = xorigi !- dxmax(1)
             else
                xi    = xorigi
             endif
#endif
             icell    = ((kpixi-1)*nsub + (jpixi-1))*nsub + ipixi
             !
             !--particle interpolates directly onto the root grid
             !
             !print*,'onto root grid ',ipixi,jpixi,kpixi
             if (vertexcen) then
                xpix = xmin(1) + (ipixi-1)*dxcell(1)
             else
                xpix = xmin(1) + (ipixi-0.5)*dxcell(1)
             endif
             !print*, "xpix: ", xpix
             !xpix = xmin(1) + (ipixi-1)*dxcell(1) ! Since we are vertex centered from Et
             dx   = xpix - xi
             q2   = dx*dx*hi21 + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
             !
             !--SPH kernel - standard cubic spline
             !
             if (q2 < radkern2) then
                !  if (q2 < 1.0) then
                !     qq = sqrt(q2)
                !     wab = 1.-1.5*q2 + 0.75*q2*qq
                !  else
                !     qq = sqrt(q2)
                !     wab = 0.25*(2.-qq)**3
                !  endif
                ! Call the kernel routine
                qq  = sqrt(q2)
                wab = wkern(q2,qq)
                !
                !--calculate data value at this pixel using the summation interpolant
                !
                ! Change this to the access the pixel coords x,y,z
                !$omp critical
                datsmooth(ipixi,jpixi,kpixi) = datsmooth(ipixi,jpixi,kpixi) + term*wab

                !if (ipixi==1 .and. jpixi==1 .and. kpixi==1) print*, "x position of 1,1,1", xi,yi,zi
                if (normalise) then
                   datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                endif
                !$omp end critical
             endif
          enddo
       enddo
    enddo
 enddo over_parts
 !$omp enddo
 !$omp end parallel

!$ do i=0,nnodes
!$  call omp_destroy_lock(ilock(i))
!$ enddo
!$ if (allocated(ilock)) deallocate(ilock)

 !
 !--normalise dat array
 !
 if (normalise) then
    where (datnorm > tiny(datnorm))
       datsmooth = datsmooth/datnorm
    end where
 endif
 if (allocated(datnorm)) deallocate(datnorm)
 !
 !--get ending CPU time
 !
 call cpu_time(t_end)
 print*,'completed in ',t_end-t_start,'s'

 return

end subroutine interpolate3D

end module interpolations3D
