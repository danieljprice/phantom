!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module interpolations3D
!
! interpolations3D
!
! :References: None
!
! :Owner: Spencer Magnall
!
! :Runtime parameters: None
!
! :Dependencies: einsteintk_utils, kernel
!

!----------------------------------------------------------------------
!
!  Module containing all of the routines required for interpolation
!  from 3D data to a 3D grid (SLOW!)
!
!----------------------------------------------------------------------

 use einsteintk_utils,  only:exact_rendering
 use kernel,       only:radkern2,radkern,cnormk,wkern!,wallint  ! Moved to this module
 !use interpolation, only:iroll ! Moved to this module

 !use timing,        only:wall_time,print_time ! Using cpu_time for now
 implicit none
 integer, parameter :: doub_prec = kind(0.d0)
 real :: cnormk3D = cnormk
 public :: interpolate3D,interpolate3D_vecexact

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
 !     Revised for "splash to grid", Monash University 02/11/09
 !     Maya Petkova contributed exact subgrid interpolation, April 2019
 !--------------------------------------------------------------------------
subroutine interpolate3D(xyzh,weight,dat,itype,npart,&
   xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidthx,pixwidthy,pixwidthz,&
   normalise,periodicx,periodicy,periodicz)

integer, intent(in) :: npart,npixx,npixy,npixz
real, intent(in)    :: xyzh(4,npart)
!real, intent(in), dimension(npart) :: x,y,z,hh ! change to xyzh()
real, intent(in), dimension(npart) :: weight,dat
integer, intent(in), dimension(npart) :: itype
real, intent(in) :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
real, intent(out), dimension(npixx,npixy,npixz) :: datsmooth
logical, intent(in) :: normalise,periodicx,periodicy,periodicz
!logical, intent(in), exact_rendering
real, allocatable :: datnorm(:,:,:)

integer :: i,ipix,jpix,kpix
integer :: iprintinterval,iprintnext
integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
integer :: ipixi,jpixi,kpixi,nxpix,nwarn,threadid
real :: xminpix,yminpix,zminpix,hmin !,dhmin3
real, dimension(npixx) :: dx2i
real :: xi,yi,zi,hi,hi1,hi21,wab,q2,const,dyz2,dz2
real :: term,termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax,dfac
real :: t_start,t_end,t_used
logical :: iprintprogress
real, dimension(npart) :: x,y,z,hh
real :: radkernel, radkernel2, radkernh

! Exact rendering
real :: pixint, wint
!logical, parameter :: exact_rendering = .true.   ! use exact rendering y/n
integer :: usedpart, negflag


!$ integer :: omp_get_num_threads,omp_get_thread_num
integer(kind=selected_int_kind(10)) :: iprogress,j  ! up to 10 digits

! Fill the particle data with xyzh
x(:) = xyzh(1,:)
y(:) = xyzh(2,:)
z(:) = xyzh(3,:)
hh(:) = xyzh(4,:)
!print*, "smoothing length: ", hh(1:10)
! cnormk3D set the value from the kernel routine
cnormk3D = cnormk
radkernel = radkern
radkernel2 = radkern2
! print*, "radkern: ", radkern
! print*, "radkernel: ",radkernel
! print*, "radkern2: ", radkern2

! print*, "npix: ", npixx, npixy,npixz

if (exact_rendering) then
print "(1x,a)",'interpolating to 3D grid (exact/Petkova+2018 on subgrid) ...'
elseif (normalise) then
print "(1x,a)",'interpolating to 3D grid (normalised) ...'
else
print "(1x,a)",'interpolating to 3D grid (non-normalised) ...'
endif
if (pixwidthx <= 0. .or. pixwidthy <= 0 .or. pixwidthz <= 0) then
print "(1x,a)",'interpolate3D: error: pixel width <= 0'
return
endif
if (any(hh(1:npart) <= tiny(hh))) then
print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
endif

!call wall_time(t_start)

datsmooth = 0.
if (normalise) then
allocate(datnorm(npixx,npixy,npixz))
datnorm = 0.
endif
!
!--print a progress report if it is going to take a long time
!  (a "long time" is, however, somewhat system dependent)
!
iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000) !.or. exact_rendering
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

usedpart = 0

xminpix = xmin !- 0.5*pixwidthx
yminpix = ymin !- 0.5*pixwidthy
zminpix = zmin !- 0.5*pixwidthz
! print*, "xminpix: ", xminpix
! print*, "yminpix: ", yminpix
! print*, "zminpix: ", zminpix
! print*, "dat: ", dat(1:10)
! print*, "weights: ", weight(1:10)
pixwidthmax = max(pixwidthx,pixwidthy,pixwidthz)
!
!--use a minimum smoothing length on the grid to make
!  sure that particles contribute to at least one pixel
!
hmin = 0.5*pixwidthmax
!dhmin3 = 1./(hmin*hmin*hmin)

const = cnormk3D  ! normalisation constant (3D)
!print*, "const: ", const
nwarn = 0
j = 0_8
threadid = 1
!
!--loop over particles
!
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
!$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
!$omp shared(xminpix,yminpix,zminpix,pixwidthx,pixwidthy,pixwidthz) &
!$omp shared(npixx,npixy,npixz,const) &
!$omp shared(datnorm,normalise,periodicx,periodicy,periodicz,exact_rendering) &
!$omp shared(hmin,pixwidthmax) &
!$omp shared(iprintprogress,iprintinterval,j) &
!$omp private(hi,xi,yi,zi,radkernh,hi1,hi21) &
!$omp private(term,termnorm,xpixi,iprogress) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
!$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
!$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
!$omp private(pixint,wint,negflag,dfac,threadid) &
!$omp firstprivate(iprintnext) &
!$omp reduction(+:nwarn,usedpart)
!$omp master
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end master

!$omp do schedule (guided, 2)
over_parts: do i=1,npart
!
!--report on progress
!
if (iprintprogress) then
  !$omp atomic
  j=j+1_8
!$     threadid = omp_get_thread_num()
  iprogress = 100*j/npart
  if (iprogress >= iprintnext .and. threadid==1) then
     write(*,"(i3,'%.')",advance='no') iprogress
     iprintnext = iprintnext + iprintinterval
  endif
endif
!
!--skip particles with itype < 0
!
if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle over_parts

hi = hh(i)
if (hi <= 0.) then
  cycle over_parts
elseif (hi < hmin) then
  !
  !--use minimum h to capture subgrid particles
  !  (get better results *without* adjusting weights)
  !
  termnorm = const*weight(i) !*(hi*hi*hi)*dhmin3
  if (.not.exact_rendering) hi = hmin
else
  termnorm = const*weight(i)
endif

!
!--set kernel related quantities
!
xi = x(i)
yi = y(i)
zi = z(i)

hi1 = 1./hi
hi21 = hi1*hi1
radkernh = radkernel*hi   ! radius of the smoothing kernel
!termnorm = const*weight(i)
term = termnorm*dat(i)
dfac = hi**3/(pixwidthx*pixwidthy*pixwidthz*const)
!dfac = hi**3/(pixwidthx*pixwidthy*const)
!
!--for each particle work out which pixels it contributes to
!
ipixmin = int((xi - radkernh - xmin)/pixwidthx)
jpixmin = int((yi - radkernh - ymin)/pixwidthy)
kpixmin = int((zi - radkernh - zmin)/pixwidthz)
ipixmax = int((xi + radkernh - xmin)/pixwidthx) + 1
jpixmax = int((yi + radkernh - ymin)/pixwidthy) + 1
kpixmax = int((zi + radkernh - zmin)/pixwidthz) + 1

if (.not.periodicx) then
  if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
  if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
endif
if (.not.periodicy) then
  if (jpixmin < 1)     jpixmin = 1
  if (jpixmax > npixy) jpixmax = npixy
endif
if (.not.periodicz) then
  if (kpixmin < 1)     kpixmin = 1
  if (kpixmax > npixz) kpixmax = npixz
endif

negflag = 0

!
!--precalculate an array of dx2 for this particle (optimisation)
!
! Check the x position of the grid cells
!open(unit=677,file="posxgrid.txt",action='write',position='append')
nxpix = 0
do ipix=ipixmin,ipixmax
  nxpix = nxpix + 1
  ipixi = ipix
  if (periodicx) ipixi = iroll(ipix,npixx)
  xpixi = xminpix + ipix*pixwidthx
  !write(677,*) ipix, xpixi
  !--watch out for errors with periodic wrapping...
  if (nxpix <= size(dx2i)) then
     dx2i(nxpix) = ((xpixi - xi)**2)*hi21
  endif
enddo

!--if particle contributes to more than npixx pixels
!  (i.e. periodic boundaries wrap more than once)
!  truncate the contribution and give warning
if (nxpix > npixx) then
  nwarn = nwarn + 1
  ipixmax = ipixmin + npixx - 1
endif
!
!--loop over pixels, adding the contribution from this particle
!
do kpix = kpixmin,kpixmax
  kpixi = kpix
  if (periodicz) kpixi = iroll(kpix,npixz)

  zpix = zminpix + kpix*pixwidthz
  dz = zpix - zi
  dz2 = dz*dz*hi21

  do jpix = jpixmin,jpixmax
     jpixi = jpix
     if (periodicy) jpixi = iroll(jpix,npixy)

     ypix = yminpix + jpix*pixwidthy
     dy = ypix - yi
     dyz2 = dy*dy*hi21 + dz2

     nxpix = 0
     do ipix = ipixmin,ipixmax
        if ((kpix==kpixmin).and.(jpix==jpixmin).and.(ipix==ipixmin)) then
           usedpart = usedpart + 1
        endif

        nxpix = nxpix + 1
        ipixi = ipix
        if (periodicx) ipixi = iroll(ipix,npixx)

        q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21

        if (exact_rendering .and. ipixmax-ipixmin <= 4) then
           if (q2 < radkernel2 + 3.*pixwidthmax**2*hi21) then
              xpixi = xminpix + ipix*pixwidthx

              ! Contribution of the cell walls in the xy-plane
              pixint = 0.0
              wint = wallint(zpix-zi+0.5*pixwidthz,xi,yi,xpixi,ypix,pixwidthx,pixwidthy,hi)
              pixint = pixint + wint

              wint = wallint(zi-zpix+0.5*pixwidthz,xi,yi,xpixi,ypix,pixwidthx,pixwidthy,hi)
              pixint = pixint + wint

              ! Contribution of the cell walls in the xz-plane
              wint = wallint(ypix-yi+0.5*pixwidthy,xi,zi,xpixi,zpix,pixwidthx,pixwidthz,hi)
              pixint = pixint + wint

              wint = wallint(yi-ypix+0.5*pixwidthy,xi,zi,xpixi,zpix,pixwidthx,pixwidthz,hi)
              pixint = pixint + wint

              ! Contribution of the cell walls in the yz-plane
              wint = wallint(xpixi-xi+0.5*pixwidthx,zi,yi,zpix,ypix,pixwidthz,pixwidthy,hi)
              pixint = pixint + wint

              wint = wallint(xi-xpixi+0.5*pixwidthx,zi,yi,zpix,ypix,pixwidthz,pixwidthy,hi)
              pixint = pixint + wint

              wab = pixint*dfac ! /(pixwidthx*pixwidthy*pixwidthz*const)*hi**3

              if (pixint < -0.01d0) then
                 print*, "Error: (",ipixi,jpixi,kpixi,") -> ", pixint, term*wab
              endif

              !
              !--calculate data value at this pixel using the summation interpolant
              !
              !$omp atomic
              datsmooth(ipixi,jpixi,kpixi) = datsmooth(ipixi,jpixi,kpixi) + term*wab
              if (normalise) then
                 !$omp atomic
                 datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
              endif
           endif
        else
           if (q2 < radkernel2) then

              !
              !--SPH kernel - standard cubic spline
              !
              wab = wkernel(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !
              !$omp atomic
              datsmooth(ipixi,jpixi,kpixi) = datsmooth(ipixi,jpixi,kpixi) + term*wab
              if (normalise) then
                 !$omp atomic
                 datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
              endif
           endif
        endif
     enddo
  enddo
enddo
enddo over_parts
!$omp enddo
!$omp end parallel

if (nwarn > 0) then
print "(a,i11,a,/,a)",' interpolate3D: WARNING: contributions truncated from ',nwarn,' particles',&
                         '                that wrap periodic boundaries more than once'
endif
!
!--normalise dat array
!
if (normalise) then
where (datnorm > tiny(datnorm))
  datsmooth = datsmooth/datnorm
end where
endif
if (allocated(datnorm)) deallocate(datnorm)

!call wall_time(t_end)
call cpu_time(t_end)
t_used = t_end - t_start
print*, 'Interpolate3D completed in ',t_end-t_start,'s'
!if (t_used > 10.) call print_time(t_used)

!print*, 'Number of particles in the volume: ', usedpart
!  datsmooth(1,1,1) = 3.14159
!  datsmooth(32,32,32) = 3.145159
!  datsmooth(11,11,11) = 3.14159
!  datsmooth(10,10,10) = 3.145159

end subroutine interpolate3D

subroutine interpolate3D_vecexact(xyzh,weight,dat,ilendat,itype,npart,&
        xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidthx,pixwidthy,pixwidthz,&
        normalise,periodicx,periodicy,periodicz)

 integer, intent(in) :: npart,npixx,npixy,npixz,ilendat
 real, intent(in)    :: xyzh(4,npart)
 !real, intent(in), dimension(npart) :: x,y,z,hh ! change to xyzh()
 real, intent(in), dimension(npart) :: weight
 real, intent(in),dimension(npart,ilendat) :: dat
 integer, intent(in), dimension(npart) :: itype
 real, intent(in) :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
 real, intent(out), dimension(ilendat,npixx,npixy,npixz) :: datsmooth
 logical, intent(in) :: normalise,periodicx,periodicy,periodicz
 !logical, intent(in), exact_rendering
 real, allocatable :: datnorm(:,:,:)

 integer :: i,ipix,jpix,kpix,lockindex,smoothindex
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nwarn,threadid
 real :: xminpix,yminpix,zminpix,hmin !,dhmin3
 real, dimension(npixx) :: dx2i
 real :: xi,yi,zi,hi,hi1,hi21,wab,q2,const,dyz2,dz2
 real :: term(ilendat),termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax,dfac
 real :: t_start,t_end,t_used
 logical :: iprintprogress
 real, dimension(npart) :: x,y,z,hh
 real :: radkernel, radkernel2, radkernh

 ! Exact rendering
 real :: pixint, wint
 !logical, parameter :: exact_rendering = .true.   ! use exact rendering y/n
 integer :: usedpart, negflag


!$ integer :: omp_get_num_threads,omp_get_thread_num
 integer(kind=selected_int_kind(10)) :: iprogress,j  ! up to 10 digits

 ! Fill the particle data with xyzh
 x(:) = xyzh(1,:)
 y(:) = xyzh(2,:)
 z(:) = xyzh(3,:)
 hh(:) = xyzh(4,:)
 !print*, "smoothing length: ", hh(1:10)
 ! cnormk3D set the value from the kernel routine
 cnormk3D = cnormk
 radkernel = radkern
 radkernel2 = radkern2
!  print*, "radkern: ", radkern
!  print*, "radkernel: ",radkernel
!  print*, "radkern2: ", radkern2

 !print*, "npix: ", npixx, npixy,npixz

 if (exact_rendering) then
    print "(1x,a)",'interpolating to 3D grid (exact/Petkova+2018 on subgrid) ...'
 elseif (normalise) then
    print "(1x,a)",'interpolating to 3D grid (normalised) ...'
 else
    print "(1x,a)",'interpolating to 3D grid (non-normalised) ...'
 endif
 if (pixwidthx <= 0. .or. pixwidthy <= 0 .or. pixwidthz <= 0) then
    print "(1x,a)",'interpolate3D: error: pixel width <= 0'
    return
 endif
 if (any(hh(1:npart) <= tiny(hh))) then
    print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
 endif

 !call wall_time(t_start)

!! $ allocate(ilock(npixx*npixy*npixz))
!! $ do i=1,npixx*npixy*npixz
!! $  call omp_init_lock(ilock(i))
!! $ enddo

 datsmooth = 0.
 if (normalise) then
    allocate(datnorm(npixx,npixy,npixz))
    datnorm = 0.
 endif
 !
 !--print a progress report if it is going to take a long time
 !  (a "long time" is, however, somewhat system dependent)
 !
 iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000) !.or. exact_rendering
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

 usedpart = 0

 xminpix = xmin !- 0.5*pixwidthx
 yminpix = ymin !- 0.5*pixwidthy
 zminpix = zmin !- 0.5*pixwidthz
!  print*, "xminpix: ", xminpix
!  print*, "yminpix: ", yminpix
!  print*, "zminpix: ", zminpix
!  print*, "dat: ", dat(1:10)
!  print*, "weights: ", weight(1:10)
 pixwidthmax = max(pixwidthx,pixwidthy,pixwidthz)
 !
 !--use a minimum smoothing length on the grid to make
 !  sure that particles contribute to at least one pixel
 !
 hmin = 0.5*pixwidthmax
 !dhmin3 = 1./(hmin*hmin*hmin)

 const = cnormk3D  ! normalisation constant (3D)
 !print*, "const: ", const
 nwarn = 0
 j = 0_8
 threadid = 1
 !
 !--loop over particles
 !
 !$omp parallel default(none) &
 !$omp shared(hh,z,x,y,weight,dat,itype,datsmooth,npart) &
 !$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
 !$omp shared(xminpix,yminpix,zminpix,pixwidthx,pixwidthy,pixwidthz) &
 !$omp shared(npixx,npixy,npixz,const,ilendat) &
 !$omp shared(datnorm,normalise,periodicx,periodicy,periodicz,exact_rendering) &
 !$omp shared(hmin,pixwidthmax) &
 !$omp shared(iprintprogress,iprintinterval,j) &
 !$omp private(hi,xi,yi,zi,radkernh,hi1,hi21) &
 !$omp private(term,termnorm,xpixi,iprogress) &
 !$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
 !$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
 !$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
 !$omp private(pixint,wint,negflag,dfac,threadid,lockindex,smoothindex) &
 !$omp firstprivate(iprintnext) &
 !$omp reduction(+:nwarn,usedpart)
 !$omp master
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
 !$omp end master

 !$omp do schedule (guided, 2)
 over_parts: do i=1,npart
    !
    !--report on progress
    !
    if (iprintprogress) then
       !$omp atomic
       j=j+1_8
!$     threadid = omp_get_thread_num()
       iprogress = 100*j/npart
       if (iprogress >= iprintnext .and. threadid==1) then
          write(*,"(i3,'%.')",advance='no') iprogress
          iprintnext = iprintnext + iprintinterval
       endif
    endif
    !
    !--skip particles with itype < 0
    !
    if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle over_parts

    hi = hh(i)
    if (hi <= 0.) then
       cycle over_parts
    elseif (hi < hmin) then
       !
       !--use minimum h to capture subgrid particles
       !  (get better results *without* adjusting weights)
       !
       termnorm = const*weight(i) !*(hi*hi*hi)*dhmin3
       if (.not.exact_rendering) hi = hmin
    else
       termnorm = const*weight(i)
    endif

    !
    !--set kernel related quantities
    !
    xi = x(i)
    yi = y(i)
    zi = z(i)

    hi1 = 1./hi
    hi21 = hi1*hi1
    radkernh = radkernel*hi   ! radius of the smoothing kernel
    !termnorm = const*weight(i)
    term(:) = termnorm*dat(i,:)
    dfac = hi**3/(pixwidthx*pixwidthy*pixwidthz*const)
    !dfac = hi**3/(pixwidthx*pixwidthy*const)
    !
    !--for each particle work out which pixels it contributes to
    !
    ipixmin = int((xi - radkernh - xmin)/pixwidthx)
    jpixmin = int((yi - radkernh - ymin)/pixwidthy)
    kpixmin = int((zi - radkernh - zmin)/pixwidthz)
    ipixmax = int((xi + radkernh - xmin)/pixwidthx) + 1
    jpixmax = int((yi + radkernh - ymin)/pixwidthy) + 1
    kpixmax = int((zi + radkernh - zmin)/pixwidthz) + 1

    if (.not.periodicx) then
       if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
       if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
    endif
    if (.not.periodicy) then
       if (jpixmin < 1)     jpixmin = 1
       if (jpixmax > npixy) jpixmax = npixy
    endif
    if (.not.periodicz) then
       if (kpixmin < 1)     kpixmin = 1
       if (kpixmax > npixz) kpixmax = npixz
    endif

    negflag = 0

    !
    !--precalculate an array of dx2 for this particle (optimisation)
    !
    ! Check the x position of the grid cells
    !open(unit=677,file="posxgrid.txt",action='write',position='append')
    nxpix = 0
    do ipix=ipixmin,ipixmax
       nxpix = nxpix + 1
       ipixi = ipix
       if (periodicx) ipixi = iroll(ipix,npixx)
       xpixi = xminpix + ipix*pixwidthx
       !write(677,*) ipix, xpixi
       !--watch out for errors with periodic wrapping...
       if (nxpix <= size(dx2i)) then
          dx2i(nxpix) = ((xpixi - xi)**2)*hi21
       endif
    enddo

    !--if particle contributes to more than npixx pixels
    !  (i.e. periodic boundaries wrap more than once)
    !  truncate the contribution and give warning
    if (nxpix > npixx) then
       nwarn = nwarn + 1
       ipixmax = ipixmin + npixx - 1
    endif
    !
    !--loop over pixels, adding the contribution from this particle
    !
    do kpix = kpixmin,kpixmax
       kpixi = kpix
       if (periodicz) kpixi = iroll(kpix,npixz)

       zpix = zminpix + kpix*pixwidthz
       dz = zpix - zi
       dz2 = dz*dz*hi21

       do jpix = jpixmin,jpixmax
          jpixi = jpix
          if (periodicy) jpixi = iroll(jpix,npixy)

          ypix = yminpix + jpix*pixwidthy
          dy = ypix - yi
          dyz2 = dy*dy*hi21 + dz2

          nxpix = 0
          do ipix = ipixmin,ipixmax
             if ((kpix==kpixmin).and.(jpix==jpixmin).and.(ipix==ipixmin)) then
                usedpart = usedpart + 1
             endif

             nxpix = nxpix + 1
             ipixi = ipix
             if (periodicx) ipixi = iroll(ipix,npixx)

             q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21

             if (exact_rendering .and. ipixmax-ipixmin <= 4) then
                if (q2 < radkernel2 + 3.*pixwidthmax**2*hi21) then
                   xpixi = xminpix + ipix*pixwidthx

                   ! Contribution of the cell walls in the xy-plane
                   pixint = 0.0
                   wint = wallint(zpix-zi+0.5*pixwidthz,xi,yi,xpixi,ypix,pixwidthx,pixwidthy,hi)
                   pixint = pixint + wint

                   wint = wallint(zi-zpix+0.5*pixwidthz,xi,yi,xpixi,ypix,pixwidthx,pixwidthy,hi)
                   pixint = pixint + wint

                   ! Contribution of the cell walls in the xz-plane
                   wint = wallint(ypix-yi+0.5*pixwidthy,xi,zi,xpixi,zpix,pixwidthx,pixwidthz,hi)
                   pixint = pixint + wint

                   wint = wallint(yi-ypix+0.5*pixwidthy,xi,zi,xpixi,zpix,pixwidthx,pixwidthz,hi)
                   pixint = pixint + wint

                   ! Contribution of the cell walls in the yz-plane
                   wint = wallint(xpixi-xi+0.5*pixwidthx,zi,yi,zpix,ypix,pixwidthz,pixwidthy,hi)
                   pixint = pixint + wint

                   wint = wallint(xi-xpixi+0.5*pixwidthx,zi,yi,zpix,ypix,pixwidthz,pixwidthy,hi)
                   pixint = pixint + wint

                   wab = pixint*dfac ! /(pixwidthx*pixwidthy*pixwidthz*const)*hi**3

                   if (pixint < -0.01d0) then
                      print*, "Error: (",ipixi,jpixi,kpixi,") -> ", pixint, term*wab
                   endif

                   !
                   !--calculate data value at this pixel using the summation interpolant
                   !
                   ! Find out where this pixel sits in the lock array 
                   ! lockindex = (k-1)*nx*ny + (j-1)*nx + i 
                   !lockindex = (kpixi-1)*npixx*npixy + (jpixi-1)*npixx + ipixi
                   !!$call omp_set_lock(ilock(lockindex))
                   !!$omp critical (datsmooth)
                   do smoothindex=1, ilendat
                     !$omp atomic 
                     datsmooth(smoothindex,ipixi,jpixi,kpixi) = datsmooth(smoothindex,ipixi,jpixi,kpixi) + term(smoothindex)*wab
                   enddo 
                   !!$omp end critical (datsmooth)
                   if (normalise) then
                      !$omp atomic
                      !!$omp critical (datnorm)
                      datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                      !!$omp end critical (datnorm)
                   endif
                   
                   !!$call omp_unset_lock(ilock(lockindex))
                endif
             else
                if (q2 < radkernel2) then

                   !
                   !--SPH kernel - standard cubic spline
                   !
                   wab = wkernel(q2)
                   !
                   !--calculate data value at this pixel using the summation interpolant
                   !
                   !!$omp atomic ! Atomic statmements only work with scalars 
                   !!$omp set lock ! Does this work with an array?
                   ! Find out where this pixel sits in the lock array 
                   ! lockindex = (k-1)*nx*ny + (j-1)*nx + i 
                   !lockindex = (kpixi-1)*npixx*npixy + (jpixi-1)*npixx + ipixi
                   !!$call omp_set_lock(ilock(lockindex)) 
                   !!$omp critical (datsmooth)
                   do smoothindex=1,ilendat
                     !$omp atomic 
                     datsmooth(smoothindex,ipixi,jpixi,kpixi) = datsmooth(smoothindex,ipixi,jpixi,kpixi) + term(smoothindex)*wab
                   enddo 
                   !!$omp end critical (datsmooth)
                   if (normalise) then
                      !$omp atomic
                      !!$omp critical (datnorm)
                      datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
                      !!$omp end critical (datnorm)
                   endif
                  !!$call omp_unset_lock(ilock(lockindex)) 
                  
                endif
             endif
          enddo
       enddo
    enddo
 enddo over_parts
 !$omp enddo
 !$omp end parallel

!!$ do i=1,npixx*npixy*npixz
!!$  call omp_destroy_lock(ilock(i))
!!$ enddo
!!$ if (allocated(ilock)) deallocate(ilock)

 if (nwarn > 0) then
    print "(a,i11,a,/,a)",' interpolate3D: WARNING: contributions truncated from ',nwarn,' particles',&
                              '                that wrap periodic boundaries more than once'
 endif
 !
 !--normalise dat array
 !
 if (normalise) then
   do i=1, ilendat
    where (datnorm > tiny(datnorm))
      
       datsmooth(i,:,:,:) = datsmooth(i,:,:,:)/datnorm(:,:,:) 
    end where
   enddo 
 endif
 if (allocated(datnorm)) deallocate(datnorm)

 !call wall_time(t_end)
 call cpu_time(t_end)
 t_used = t_end - t_start
 print*, 'Interpolate3DVec completed in ',t_end-t_start,'s'
 !if (t_used > 10.) call print_time(t_used)

 !print*, 'Number of particles in the volume: ', usedpart
 !  datsmooth(1,1,1) = 3.14159
 !  datsmooth(32,32,32) = 3.145159
 !  datsmooth(11,11,11) = 3.14159
 !  datsmooth(10,10,10) = 3.145159

end subroutine interpolate3D_vecexact

 ! subroutine interpolate3D_vec(x,y,z,hh,weight,datvec,itype,npart,&
 !      xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidthx,pixwidthy,pixwidthz,&
 !      normalise,periodicx,periodicy,periodicz)

 !  integer, intent(in) :: npart,npixx,npixy,npixz
 !  real, intent(in), dimension(npart)    :: x,y,z,hh,weight
 !  real, intent(in), dimension(npart,3)  :: datvec
 !  integer, intent(in), dimension(npart) :: itype
 !  real, intent(in) :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
 !  real(doub_prec), intent(out), dimension(3,npixx,npixy,npixz) :: datsmooth
 !  logical, intent(in) :: normalise,periodicx,periodicy,periodicz
 !  real(doub_prec), dimension(npixx,npixy,npixz) :: datnorm

 !  integer :: i,ipix,jpix,kpix
 !  integer :: iprintinterval,iprintnext
 !  integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 !  integer :: ipixi,jpixi,kpixi,nxpix,nwarn
 !  real :: xminpix,yminpix,zminpix
 !  real, dimension(npixx) :: dx2i
 !  real :: xi,yi,zi,hi,hi1,hi21,radkern,wab,q2,const,dyz2,dz2
 !  real :: termnorm,dy,dz,ypix,zpix,xpixi,ddatnorm
 !  real, dimension(3) :: term
 !  !real :: t_start,t_end
 !  logical :: iprintprogress
 !  !$ integer :: omp_get_num_threads
 !  integer(kind=selected_int_kind(10)) :: iprogress  ! up to 10 digits

 !  datsmooth = 0.
 !  datnorm = 0.
 !  if (normalise) then
 !     print "(1x,a)",'interpolating to 3D grid (normalised) ...'
 !  else
 !     print "(1x,a)",'interpolating to 3D grid (non-normalised) ...'
 !  endif
 !  if (pixwidthx <= 0. .or. pixwidthy <= 0. .or. pixwidthz <= 0.) then
 !     print "(1x,a)",'interpolate3D: error: pixel width <= 0'
 !     return
 !  endif
 !  if (any(hh(1:npart) <= tiny(hh))) then
 !     print*,'interpolate3D: WARNING: ignoring some or all particles with h < 0'
 !  endif

 !  !
 !  !--print a progress report if it is going to take a long time
 !  !  (a "long time" is, however, somewhat system dependent)
 !  !
 !  iprintprogress = (npart  >=  100000) .or. (npixx*npixy  > 100000)
 !  !$ iprintprogress = .false.
 !  !
 !  !--loop over particles
 !  !
 !  iprintinterval = 25
 !  if (npart >= 1e6) iprintinterval = 10
 !  iprintnext = iprintinterval
 !  !
 !  !--get starting CPU time
 !  !
 !  !call cpu_time(t_start)

 !  xminpix = xmin - 0.5*pixwidthx
 !  yminpix = ymin - 0.5*pixwidthy
 !  zminpix = zmin - 0.5*pixwidthz

 !  const = cnormk3D  ! normalisation constant (3D)
 !  nwarn = 0

 ! !$omp parallel default(none) &
 ! !$omp shared(hh,z,x,y,weight,datvec,itype,datsmooth,npart) &
 ! !$omp shared(xmin,ymin,zmin,radkernel,radkernel2) &
 ! !$omp shared(xminpix,yminpix,zminpix,pixwidthx,pixwidthy,pixwidthz) &
 ! !$omp shared(npixx,npixy,npixz,const) &
 ! !$omp shared(iprintprogress,iprintinterval) &
 ! !$omp shared(datnorm,normalise,periodicx,periodicy,periodicz) &
 ! !$omp private(hi,xi,yi,zi,radkern,hi1,hi21) &
 ! !$omp private(term,termnorm,xpixi) &
 ! !$omp private(iprogress,iprintnext) &
 ! !$omp private(ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax) &
 ! !$omp private(ipix,jpix,kpix,ipixi,jpixi,kpixi) &
 ! !$omp private(dx2i,nxpix,zpix,dz,dz2,dyz2,dy,ypix,q2,wab) &
 ! !$omp reduction(+:nwarn)
 ! !$omp master
 ! !$  print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
 ! !$omp end master
 !  !
 !  !--loop over particles
 !  !
 ! !$omp do schedule (guided, 2)
 !  over_parts: do i=1,npart
 !     !
 !     !--report on progress
 !     !
 !     if (iprintprogress) then
 !        iprogress = 100*i/npart
 !        if (iprogress >= iprintnext) then
 !           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
 !           iprintnext = iprintnext + iprintinterval
 !        endif
 !     endif
 !     !
 !     !--skip particles with itype < 0
 !     !
 !     if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle over_parts

 !     hi = hh(i)
 !     if (hi <= 0.) cycle over_parts

 !     !
 !     !--set kernel related quantities
 !     !
 !     xi = x(i)
 !     yi = y(i)
 !     zi = z(i)

 !     hi1 = 1./hi
 !     hi21 = hi1*hi1
 !     radkern = radkernel*hi   ! radius of the smoothing kernel
 !     termnorm = const*weight(i)
 !     term(:) = termnorm*datvec(i,:)
 !     !
 !     !--for each particle work out which pixels it contributes to
 !     !
 !     ipixmin = int((xi - radkern - xmin)/pixwidthx)
 !     jpixmin = int((yi - radkern - ymin)/pixwidthy)
 !     kpixmin = int((zi - radkern - zmin)/pixwidthz)
 !     ipixmax = int((xi + radkern - xmin)/pixwidthx) + 1
 !     jpixmax = int((yi + radkern - ymin)/pixwidthy) + 1
 !     kpixmax = int((zi + radkern - zmin)/pixwidthz) + 1

 !     if (.not.periodicx) then
 !        if (ipixmin < 1)     ipixmin = 1      ! make sure they only contribute
 !        if (ipixmax > npixx) ipixmax = npixx  ! to pixels in the image
 !     endif
 !     if (.not.periodicy) then
 !        if (jpixmin < 1)     jpixmin = 1
 !        if (jpixmax > npixy) jpixmax = npixy
 !     endif
 !     if (.not.periodicz) then
 !        if (kpixmin < 1)     kpixmin = 1
 !        if (kpixmax > npixz) kpixmax = npixz
 !     endif
 !     !
 !     !--precalculate an array of dx2 for this particle (optimisation)
 !     !
 !     nxpix = 0
 !     do ipix=ipixmin,ipixmax
 !        nxpix = nxpix + 1
 !        ipixi = ipix
 !        if (periodicx) ipixi = iroll(ipix,npixx)
 !        xpixi = xminpix + ipix*pixwidthx
 !        !--watch out for errors with perioic wrapping...
 !        if (nxpix <= size(dx2i)) then
 !           dx2i(nxpix) = ((xpixi - xi)**2)*hi21
 !        endif
 !     enddo

 !     !--if particle contributes to more than npixx pixels
 !     !  (i.e. periodic boundaries wrap more than once)
 !     !  truncate the contribution and give warning
 !     if (nxpix > npixx) then
 !        nwarn = nwarn + 1
 !        ipixmax = ipixmin + npixx - 1
 !     endif
 !     !
 !     !--loop over pixels, adding the contribution from this particle
 !     !
 !     do kpix = kpixmin,kpixmax
 !        kpixi = kpix
 !        if (periodicz) kpixi = iroll(kpix,npixz)
 !        zpix = zminpix + kpix*pixwidthz
 !        dz = zpix - zi
 !        dz2 = dz*dz*hi21

 !        do jpix = jpixmin,jpixmax
 !           jpixi = jpix
 !           if (periodicy) jpixi = iroll(jpix,npixy)
 !           ypix = yminpix + jpix*pixwidthy
 !           dy = ypix - yi
 !           dyz2 = dy*dy*hi21 + dz2

 !           nxpix = 0
 !           do ipix = ipixmin,ipixmax
 !              ipixi = ipix
 !              if (periodicx) ipixi = iroll(ipix,npixx)
 !              nxpix = nxpix + 1
 !              q2 = dx2i(nxpix) + dyz2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
 !              !
 !              !--SPH kernel - standard cubic spline
 !              !
 !              if (q2 < radkernel2) then
 !                 wab = wkernel(q2)
 !                 !
 !                 !--calculate data value at this pixel using the summation interpolant
 !                 !
 !                 !$omp atomic
 !                 datsmooth(1,ipixi,jpixi,kpixi) = datsmooth(1,ipixi,jpixi,kpixi) + term(1)*wab
 !                 !$omp atomic
 !                 datsmooth(2,ipixi,jpixi,kpixi) = datsmooth(2,ipixi,jpixi,kpixi) + term(2)*wab
 !                 !$omp atomic
 !                 datsmooth(3,ipixi,jpixi,kpixi) = datsmooth(3,ipixi,jpixi,kpixi) + term(3)*wab
 !                 if (normalise) then
 !                    !$omp atomic
 !                    datnorm(ipixi,jpixi,kpixi) = datnorm(ipixi,jpixi,kpixi) + termnorm*wab
 !                 endif
 !              endif
 !           enddo
 !        enddo
 !     enddo
 !  enddo over_parts
 ! !$omp enddo
 ! !$omp end parallel

 !  if (nwarn > 0) then
 !     print "(a,i11,a,/,a)",' interpolate3D: WARNING: contributions truncated from ',nwarn,' particles',&
 !                            '                that wrap periodic boundaries more than once'
 !  endif
 !  !
 !  !--normalise dat array
 !  !
 !  if (normalise) then
 !     !$omp parallel do default(none) schedule(static) &
 !     !$omp shared(datsmooth,datnorm,npixz,npixy,npixx) &
 !     !$omp private(kpix,jpix,ipix,ddatnorm)
 !     do kpix=1,npixz
 !        do jpix=1,npixy
 !           do ipix=1,npixx
 !              if (datnorm(ipix,jpix,kpix) > tiny(datnorm)) then
 !                 ddatnorm = 1./datnorm(ipix,jpix,kpix)
 !                 datsmooth(1,ipix,jpix,kpix) = datsmooth(1,ipix,jpix,kpix)*ddatnorm
 !                 datsmooth(2,ipix,jpix,kpix) = datsmooth(2,ipix,jpix,kpix)*ddatnorm
 !                 datsmooth(3,ipix,jpix,kpix) = datsmooth(3,ipix,jpix,kpix)*ddatnorm
 !              endif
 !           enddo
 !        enddo
 !     enddo
 !     !$omp end parallel do
 !  endif

 !  return

 ! end subroutine interpolate3D_vec

 !------------------------------------------------------------
 ! interface to kernel routine to avoid problems with openMP
 !-----------------------------------------------------------
real function wkernel(q2)
 use kernel, only:wkern
 real, intent(in) :: q2
 real :: q
 q = sqrt(q2)
 wkernel = wkern(q2,q)

end function wkernel

 !------------------------------------------------------------
 ! 3D functions to evaluate exact overlap of kernel with wall boundaries
 ! see Petkova, Laibe & Bonnell (2018), J. Comp. Phys
 !------------------------------------------------------------
real function wallint(r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi)
 real, intent(in) :: r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi
 real(doub_prec) :: R_0, d1, d2, dx, dy, h

 wallint = 0.0
 dx = xc - xp
 dy = yc - yp
 h = hi

 !
 ! Contributions from each of the 4 sides of a cell wall
 !
 R_0 = 0.5*pixwidthy + dy
 d1 = 0.5*pixwidthx - dx
 d2 = 0.5*pixwidthx + dx
 wallint = wallint + pint3D(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthy - dy
 d1 = 0.5*pixwidthx + dx
 d2 = 0.5*pixwidthx - dx
 wallint = wallint + pint3D(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthx + dx
 d1 = 0.5*pixwidthy + dy
 d2 = 0.5*pixwidthy - dy
 wallint = wallint + pint3D(r0, R_0, d1, d2, h)

 R_0 = 0.5*pixwidthx - dx
 d1 = 0.5*pixwidthy - dy
 d2 = 0.5*pixwidthy + dy
 wallint = wallint + pint3D(r0, R_0, d1, d2, h)

end function wallint


real function pint3D(r0, R_0, d1, d2, hi)

 real(doub_prec), intent(in) :: R_0, d1, d2, hi
 real, intent(in) :: r0
 real(doub_prec) :: ar0, aR_0
 real(doub_prec) :: int1, int2
 !integer :: fflag = 0

 if (abs(r0) < tiny(0.)) then
    pint3D = 0.d0
    return
 endif

 if (r0  >  0.d0) then
    pint3D = 1.d0
    ar0 = r0
 else
    pint3D = -1.d0
    ar0 = -r0
 endif

 if (R_0  >  0.d0) then
    aR_0 = R_0
 else
    pint3D = -pint3D
    aR_0 = -R_0
 endif

 int1 = full_integral_3D(d1, ar0, aR_0, hi)
 int2 = full_integral_3D(d2, ar0, aR_0, hi)

 if (int1 < 0.d0) int1 = 0.d0
 if (int2 < 0.d0) int2 = 0.d0

 if (d1*d2  >=  0) then
    pint3D = pint3D*(int1 + int2)
    if (int1 + int2 < 0.d0) print*, 'Error: int1 + int2 < 0'
 elseif (abs(d1)  <  abs(d2)) then
    pint3D = pint3D*(int2 - int1)
    if (int2 - int1 < 0.d0) print*, 'Error: int2 - int1 < 0: ', int1, int2, '(', d1, d2,')'
 else
    pint3D = pint3D*(int1 - int2)
    if (int1 - int2 < 0.d0) print*, 'Error: int1 - int2 < 0: ', int1, int2, '(', d1, d2,')'
 endif

end function pint3D

real(doub_prec) function full_integral_3D(d, r0, R_0, h)

 real(doub_prec), intent(in) :: d, r0, R_0, h
 real(doub_prec) :: B1, B2, B3, a, h2
 real(doub_prec), parameter :: pi = 4.*atan(1.)
 real(doub_prec) :: tanphi, phi, a2, cosp, r0h, r03, r0h2, r0h3, r0h_2, r0h_3
 real(doub_prec) :: r2, R_, linedist2, cosphi
 real(doub_prec) :: I0, I1, I_2, I_3, I_4, I_5
 real(doub_prec) :: D2, D3

 r0h = r0/h
 tanphi = abs(d)/R_0
 phi = atan(tanphi)

 if (abs(r0h) < tiny(0.) .or. abs(R_0/h) < tiny(0.) .or. abs(phi) < tiny(0.)) then
    full_integral_3D = 0.0
    return
 endif

 h2 = h*h
 r03 = r0*r0*r0
 r0h2 = r0h*r0h
 r0h3 = r0h2*r0h
 r0h_2 = 1./r0h2
 r0h_3 = 1./r0h3

 ! Avoid Compiler warnings
 B1 = 0.
 B2 = 0.

 if (r0 >= 2.0*h) then
    B3 = 0.25*h2*h
 elseif (r0 > h) then
    B3 = 0.25*r03 *(-4./3. + (r0h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3+ 8./5.*r0h_2)
    B2 = 0.25*r03 *(-4./3. + (r0h) - 0.3*r0h2 + 1./30.*r0h3 - 1./15. *r0h_3)
 else
    B3 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3 + 7./5.*r0h_2)
    B2 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3 - 1./5.*r0h_2)
    B1 = 0.25*r03 *(-2./3. + 0.3*r0h2 - 0.1*r0h3)
 endif

 a = R_0/r0
 a2 = a*a

 linedist2 = (r0*r0 + R_0*R_0)
 cosphi = cos(phi)
 R_ = R_0/cosphi
 r2 = (r0*r0 + R_*R_)

 D2 = 0.0
 D3 = 0.0

 if (linedist2 < h2) then
    !////// phi1 business /////
    cosp = R_0/sqrt(h2-r0*r0)
    call get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)

    D2 = -1./6.*I_2 + 0.25*(r0h) *I_3 - 0.15*r0h2 *I_4 + 1./30.*r0h3 *I_5 - 1./60. *r0h_3 *I1 + (B1-B2)/r03 *I0
 endif
 if (linedist2 < 4.*h2) then
    !////// phi2 business /////
    cosp = R_0/sqrt(4.0*h2-r0*r0)
    call get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5)

    D3 = 1./3.*I_2 - 0.25*(r0h) *I_3 + 3./40.*r0h2 *I_4 - 1./120.*r0h3 *I_5 + 4./15. *r0h_3 *I1 + (B2-B3)/r03 *I0 + D2
 endif

 !//////////////////////////////
 call get_I_terms(cosphi,a2,a,I0,I1,I_2,I_3,I_4,I_5,phi=phi,tanphi=tanphi)

 if (r2 < h2) then
    full_integral_3D = r0h3/pi  * (1./6. *I_2 - 3./40.*r0h2 *I_4 + 1./40.*r0h3 *I_5 + B1/r03 *I0)
 elseif (r2 < 4.*h2) then
    full_integral_3D=  r0h3/pi  * (0.25 * (4./3. *I_2 - (r0/h) *I_3 + 0.3*r0h2 *I_4 - &
    &   1./30.*r0h3 *I_5 + 1./15. *r0h_3 *I1) + B2/r03 *I0 + D2)
 else
    full_integral_3D = r0h3/pi  * (-0.25*r0h_3 *I1 + B3/r03 *I0 + D3)
 endif

end function full_integral_3D

subroutine get_I_terms(cosp,a2,a,I0,I1,I_2,I_3,I_4,I_5,phi,tanphi)
 real(doub_prec), intent(in) :: cosp,a2,a
 real(doub_prec), intent(out) :: I0,I1,I_2,I_3,I_4,I_5
 real(doub_prec), intent(in), optional :: phi,tanphi
 real(doub_prec) :: cosp2,p,tanp,u2,u,logs,I_1,mu2_1,fac

 cosp2 = cosp*cosp
 if (present(phi)) then
    p = phi
    tanp = tanphi
 else
    p = acos(cosp)
    tanp = sqrt(1.-cosp2)/cosp ! tan(p)
 endif

 mu2_1 = 1. / (1. + cosp2/a2)
 I0  = p
 I_2 = p +    a2 * tanp
 I_4 = p + 2.*a2 * tanp + 1./3.*a2*a2 * tanp*(2. + 1./cosp2)

 u2 = (1.-cosp2)*mu2_1
 u = sqrt(u2)
 logs = log((1.+u)/(1.-u))
 I1 = atan2(u,a)

 fac = 1./(1.-u2)
 I_1 = 0.5*a*logs + I1
 I_3 = I_1 + a*0.25*(1.+a2)*(2.*u*fac + logs)
 I_5 = I_3 + a*(1.+a2)*(1.+a2)/16. *( (10.*u - 6.*u*u2)*fac*fac + 3.*logs)

end subroutine get_I_terms

 !------------------------------------------------------------
 ! function to return a soft maximum for 1/x with no bias
 ! for x >> eps using the cubic spline kernel softening
 ! i.e. something equivalent to 1/sqrt(x**2 + eps**2) but
 ! with compact support, i.e. f=1/x when x > 2*eps
 !------------------------------------------------------------
pure elemental real function soft_func(x,eps) result(f)
 real, intent(in)  :: x,eps
 real :: q,q2,q4

 q = x/eps
 q2 = q*q
 if (q < 1.) then
    q4 = q2*q2
    f = (1./eps)*(q4*q/10. - 3.*q4/10. + 2.*q2/3. - 7./5.)
 elseif (q < 2.) then
    q4 = q2*q2
    f = (1./eps)*(q*(-q4*q + 9.*q4 - 30.*q2*q + 40.*q2 - 48.) + 2.)/(30.*q)
 else
    f = -1./x
 endif
 f = -f

end function soft_func

 !--------------------------------------------------------------------------
 !
 !  utility to wrap pixel index around periodic domain
 !  indices that roll beyond the last position are re-introduced at the first
 !
 !--------------------------------------------------------------------------
pure integer function iroll(i,n)
 integer, intent(in) :: i,n

 if (i > n) then
    iroll = mod(i-1,n) + 1
 elseif (i < 1) then
    iroll = n + mod(i,n) ! mod is negative
 else
    iroll = i
 endif

end function iroll
end module interpolations3D

