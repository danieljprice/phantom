!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
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
! :Dependencies: einsteintk_utils, kernel, kernel_integrals
!
 use einsteintk_utils,  only:exact_rendering
 use kernel,            only:radkern2,radkern,cnormk,wkern
 use kernel_integrals,  only:maxpieces,maxdeg,maxpow
 implicit none
 integer, parameter :: doub_prec = kind(0.d0)
 real :: cnormk3D = cnormk
 ! coefficient tables for the exact integrals over cells (set_exact_tables)
 logical :: tables_set = .false.
 integer :: npieces_int = 0, nmax_int = 0
 real(doub_prec) :: qb_int(maxpieces) = 0., a3_int(-1:maxpow,maxpieces+1) = 0.
 public :: interpolate3D,interpolate3D_vecexact,interpolate3D_nodes

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

 integer, intent(in)  :: npart,npixx,npixy,npixz
 real,    intent(in)  :: xyzh(4,npart)
 real,    intent(in)  :: weight(npart),dat(npart)
 integer, intent(in)  :: itype(npart)
 real,    intent(in)  :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
 real,    intent(out) :: datsmooth(npixx,npixy,npixz)
 logical, intent(in)  :: normalise,periodicx,periodicy,periodicz
 real, allocatable :: datnorm(:,:,:)

 integer :: i,ipix,jpix,kpix
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nwarn,threadid
 real :: xminpix,yminpix,zminpix,hmin !,dhmin3
 real :: dx2i(npixx)
 real :: xi,yi,zi,hi,hi1,hi21,wab,q2,const,dyz2,dz2
 real :: term,termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax,dfac
 real :: t_start,t_end,t_used
 logical :: iprintprogress
 real :: x(npart),y(npart),z(npart),hh(npart)
 real :: radkernel, radkernel2, radkernh

! Exact rendering
 real :: pixint, wint
 integer :: usedpart, negflag

!$ integer, external :: omp_get_num_threads,omp_get_thread_num
 integer(kind=selected_int_kind(10)) :: iprogress,j  ! up to 10 digits

! Fill the particle data with xyzh
 x(:) = xyzh(1,:)
 y(:) = xyzh(2,:)
 z(:) = xyzh(3,:)
 hh(:) = xyzh(4,:)
 cnormk3D = cnormk
 radkernel = radkern
 radkernel2 = radkern2

 if (exact_rendering) call set_exact_tables()
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
!$omp single
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
!$omp end single nowait

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
    term = termnorm*dat(i)
    dfac = hi**3/(pixwidthx*pixwidthy*pixwidthz*const)
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
    nxpix = 0
    do ipix=ipixmin,ipixmax
       nxpix = nxpix + 1
       ipixi = ipix
       if (periodicx) ipixi = iroll(ipix,npixx)
       xpixi = xminpix + ipix*pixwidthx
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

end subroutine interpolate3D

subroutine interpolate3D_vecexact(xyzh,weight,dat,ilendat,itype,npart,&
        xmin,ymin,zmin,datsmooth,npixx,npixy,npixz,pixwidthx,pixwidthy,pixwidthz,&
        normalise,periodicx,periodicy,periodicz)

 integer, intent(in)  :: npart,npixx,npixy,npixz,ilendat
 real,    intent(in)  :: xyzh(4,npart)
 real,    intent(in)  :: weight(npart)
 real,    intent(in)  :: dat(npart,ilendat)
 integer, intent(in)  :: itype(npart)
 real,    intent(in)  :: xmin,ymin,zmin,pixwidthx,pixwidthy,pixwidthz
 real,    intent(out) :: datsmooth(ilendat,npixx,npixy,npixz)
 logical, intent(in)  :: normalise,periodicx,periodicy,periodicz
 !logical, intent(in), exact_rendering
 real, allocatable :: datnorm(:,:,:)

 integer :: i,ipix,jpix,kpix,smoothindex
 integer :: iprintinterval,iprintnext
 integer :: ipixmin,ipixmax,jpixmin,jpixmax,kpixmin,kpixmax
 integer :: ipixi,jpixi,kpixi,nxpix,nwarn,threadid
 real :: xminpix,yminpix,zminpix,hmin !,dhmin3
 real :: dx2i(npixx)
 real :: xi,yi,zi,hi,hi1,hi21,wab,q2,const,dyz2,dz2
 real :: term(ilendat),termnorm,dy,dz,ypix,zpix,xpixi,pixwidthmax,dfac
 real :: t_start,t_end,t_used
 logical :: iprintprogress
 real :: x(npart),y(npart),z(npart),hh(npart)
 real :: radkernel, radkernel2, radkernh

 ! Exact rendering
 real :: pixint, wint
 !logical, parameter :: exact_rendering = .true.   ! use exact rendering y/n
 integer :: usedpart, negflag

!$ integer, external :: omp_get_num_threads,omp_get_thread_num
 integer(kind=selected_int_kind(10)) :: iprogress,j  ! up to 10 digits

 ! Fill the particle data with xyzh
 x(:) = xyzh(1,:)
 y(:) = xyzh(2,:)
 z(:) = xyzh(3,:)
 hh(:) = xyzh(4,:)
 cnormk3D = cnormk
 radkernel = radkern
 radkernel2 = radkern2

 if (exact_rendering) call set_exact_tables()
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
 !$omp private(pixint,wint,negflag,dfac,threadid,smoothindex) &
 !$omp firstprivate(iprintnext) &
 !$omp reduction(+:nwarn,usedpart)
 !$omp single
!$ print "(1x,a,i3,a)",'Using ',omp_get_num_threads(),' cpus'
 !$omp end single nowait

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
    nxpix = 0
    do ipix=ipixmin,ipixmax
       nxpix = nxpix + 1
       ipixi = ipix
       if (periodicx) ipixi = iroll(ipix,npixx)
       xpixi = xminpix + ipix*pixwidthx
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
                   do smoothindex=1, ilendat
                      !$omp atomic
                      datsmooth(smoothindex,ipixi,jpixi,kpixi) = datsmooth(smoothindex,ipixi,jpixi,kpixi) + term(smoothindex)*wab
                   enddo
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
                   do smoothindex=1,ilendat
                      !$omp atomic
                      datsmooth(smoothindex,ipixi,jpixi,kpixi) = datsmooth(smoothindex,ipixi,jpixi,kpixi) + term(smoothindex)*wab
                   enddo
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
    do i=1, ilendat
       where (datnorm > tiny(datnorm))

          datsmooth(i,:,:,:) = datsmooth(i,:,:,:)/datnorm(:,:,:)
       end where
    enddo
 endif
 if (allocated(datnorm)) deallocate(datnorm)

 call cpu_time(t_end)
 t_used = t_end - t_start
 print*, 'Interpolate3DVec completed in ',t_end-t_start,'s'

end subroutine interpolate3D_vecexact

!--------------------------------------------------------------------------
!+
!  interpolate vector data from particles to the nodes of a uniform
!  periodic grid (node (i,j,k) at xmin + (i-1,j-1,k-1)*dx), by sampling
!  the kernel at the nodes,
!
!     datsmooth(:,node) = sum_b weight_b dat_b W(r_node - r_b, h_b)
!
!  with weight = m/(rho h^3). If zero_moment is set, the kernel is
!  B W(q) (1 - A q^2) with A, B chosen so that it is normalised and its
!  second moment vanishes: the smoothing then changes a Fourier mode by
!  O((kh)^4) rather than O((kh)^2).
!
!  The smoothing length used is max(h, hmin). A kernel sampled at the
!  nodes only conserves each particle's contribution if it spans a few
!  nodes (the error in the sum of the sampled kernel is < 1e-2 for
!  h >= dx, but up to O(1) for h = 0.6 dx), so if exact is set, particles
!  whose smoothing length is still below the grid spacing contribute
!  instead the exact integral of the kernel over each cell centred on a
!  node, divided by the cell volume
!+
!--------------------------------------------------------------------------
subroutine interpolate3D_nodes(xyzh,weight,dat,ilendat,itype,npart,xmin,dx,n,zero_moment,hmin,exact,datsmooth)
 use kernel,           only:kernelname
 use kernel_integrals, only:maxpieces,maxdeg,maxpow,get_kernel_pieces,zero_moment_pieces,&
                            set_kernel_integral_coeffs,cell_integral_3d
 integer, intent(in)  :: ilendat,npart,n(3),itype(npart)
 real,    intent(in)  :: xyzh(:,:),weight(npart),dat(npart,ilendat),xmin(3),dx(3),hmin
 logical, intent(in)  :: zero_moment,exact
 real,    intent(out) :: datsmooth(ilendat,n(1),n(2),n(3))
 real    :: a,b,hi,hi1,q2,wi,dr(3),xnode(3),wfac,radi
 real(kind=doub_prec) :: qb(maxpieces),wpoly(0:maxdeg,maxpieces),g2(-1:maxpow,maxpieces+1)
 real(kind=doub_prec) :: a3(-1:maxpow,maxpieces+1),cnorm2,cnorm3,aint,vint
 integer :: i,ip,jp,kp,ii,jj,kk,lo(3),hi_idx(3),c,npieces,nmax,ierr
 logical :: use_exact,exacti

 a = 0.
 b = 1.
 if (zero_moment) call zero_moment_constants(a,b)
 ! coefficient tables for the exact integrals of the same kernel over cells
 use_exact = exact
 nmax = 0
 if (use_exact) then
    call get_kernel_pieces(kernelname,npieces,qb,wpoly,ierr)
    if (ierr /= 0) then
       use_exact = .false.
    else
       if (zero_moment) call zero_moment_pieces(npieces,qb,wpoly,aint)
       call set_kernel_integral_coeffs(npieces,qb,wpoly,nmax,g2,a3,cnorm2,cnorm3)
    endif
 endif
 datsmooth = 0.
 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,weight,dat,ilendat,itype,xmin,dx,n,a,b,hmin,datsmooth) &
 !$omp shared(use_exact,npieces,nmax,qb,a3) &
 !$omp private(i,ip,jp,kp,ii,jj,kk,lo,hi_idx,hi,hi1,q2,wi,dr,xnode,c,wfac,radi,exacti,vint) &
 !$omp schedule(guided)
 do i=1,npart
    if (itype(i) < 0 .or. weight(i) < tiny(0.)) cycle
    hi = max(xyzh(4,i),hmin)
    exacti = use_exact .and. hi < minval(dx)
    if (exacti) then
       ! cell integrals: weight = m/(rho h^3), so m/rho over the cell volume
       wfac = weight(i)*xyzh(4,i)**3/product(dx)
       radi = radkern*hi + 0.5*maxval(dx)
    else
       ! the weight holds 1/h^3 for the particle's own h
       wfac = weight(i)*(xyzh(4,i)/hi)**3*cnormk*b
       radi = radkern*hi
    endif
    hi1 = 1./hi
    ! range of nodes whose cells or points the kernel reaches (may extend
    ! over the periodic boundary)
    lo     = ceiling((xyzh(1:3,i) - radi - xmin)/dx) + 1
    hi_idx = floor((xyzh(1:3,i) + radi - xmin)/dx) + 1
    do kp=lo(3),hi_idx(3)
       kk = modulo(kp-1,n(3)) + 1
       xnode(3) = xmin(3) + (kp-1)*dx(3)
       do jp=lo(2),hi_idx(2)
          jj = modulo(jp-1,n(2)) + 1
          xnode(2) = xmin(2) + (jp-1)*dx(2)
          do ip=lo(1),hi_idx(1)
             ii = modulo(ip-1,n(1)) + 1
             xnode(1) = xmin(1) + (ip-1)*dx(1)
             if (exacti) then
                vint = cell_integral_3d(real(xyzh(1:3,i),doub_prec),real(xnode,doub_prec),&
                                        real(dx,doub_prec),real(hi,doub_prec),npieces,nmax,qb,a3)
                wi = wfac*real(vint)
             else
                dr = xnode - xyzh(1:3,i)
                q2 = dot_product(dr,dr)*hi1*hi1
                if (q2 >= radkern2) cycle
                wi = wfac*wkern(q2,sqrt(q2))*(1. - a*q2)
             endif
             do c=1,ilendat
                !$omp atomic
                datsmooth(c,ii,jj,kk) = datsmooth(c,ii,jj,kk) + wi*dat(i,c)
             enddo
          enddo
       enddo
    enddo
 enddo
 !$omp end parallel do

end subroutine interpolate3D_nodes

!--------------------------------------------------------------------------
!+
!  constants A, B of the zero-second-moment kernel B W(q)(1 - A q^2) for
!  the kernel in use: A = <q^2>/<q^4> and B = 1/(1 - A <q^2>), with the
!  moments <q^p> = 4 pi cnormk int q^(p+2) W(q) dq (midpoint rule)
!+
!--------------------------------------------------------------------------
subroutine zero_moment_constants(a,b)
 real, intent(out) :: a,b
 integer, parameter :: nq = 20000
 real    :: dq,q,w,m2,m4
 integer :: i

 dq = radkern/nq
 m2 = 0.
 m4 = 0.
 do i=1,nq
    q = (i - 0.5)*dq
    w = 4.*acos(-1.)*cnormk*wkern(q*q,q)*q**2*dq
    m2 = m2 + w*q**2
    m4 = m4 + w*q**4
 enddo
 a = m2/m4
 b = 1./(1. - a*m2)

end subroutine zero_moment_constants

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
!+
!  coefficient tables for the exact integrals of the kernel in use
!  over cell walls (set once, before any parallel region)
!+
!------------------------------------------------------------
subroutine set_exact_tables()
 use kernel,           only:kernelname
 use kernel_integrals, only:get_kernel_pieces,set_kernel_integral_coeffs
 real(doub_prec) :: wpoly(0:maxdeg,maxpieces),g2(-1:maxpow,maxpieces+1),cnorm2,cnorm3
 integer :: ierr

 if (tables_set) return
 call get_kernel_pieces(kernelname,npieces_int,qb_int,wpoly,ierr)
 if (ierr /= 0) then
    print "(a)",' ERROR: exact interpolation not available for kernel '//trim(kernelname)
    npieces_int = 0
 else
    call set_kernel_integral_coeffs(npieces_int,qb_int,wpoly,nmax_int,g2,a3_int,cnorm2,cnorm3)
 endif
 tables_set = .true.

end subroutine set_exact_tables

!------------------------------------------------------------
!+
!  exact integral of the kernel over a cell wall, the sum of the
!  contributions from its four sides (Petkova, Laibe & Bonnell 2018,
!  generalised to any kernel as in splash): r0 is the signed distance
!  from the particle to the wall plane, (xp,yp) the particle and (xc,yc)
!  the wall centre in the plane of the wall
!+
!------------------------------------------------------------
real function wallint(r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi)
 real, intent(in) :: r0, xp, yp, xc, yc, pixwidthx, pixwidthy, hi
 real :: dx, dy

 dx = xc - xp
 dy = yc - yp
 wallint = pint3D(r0, 0.5*pixwidthy + dy, 0.5*pixwidthx - dx, 0.5*pixwidthx + dx, hi) &
         + pint3D(r0, 0.5*pixwidthy - dy, 0.5*pixwidthx + dx, 0.5*pixwidthx - dx, hi) &
         + pint3D(r0, 0.5*pixwidthx + dx, 0.5*pixwidthy + dy, 0.5*pixwidthy - dy, hi) &
         + pint3D(r0, 0.5*pixwidthx - dx, 0.5*pixwidthy - dy, 0.5*pixwidthy + dy, hi)

end function wallint

!------------------------------------------------------------
!+
!  contribution from one side of a cell wall: r0 is the distance to the
!  wall plane, R_0 the distance in the plane to the side and d1, d2 the
!  distances along the side to each vertex
!+
!------------------------------------------------------------
real function pint3D(r0, R_0, d1, d2, h)
 use kernel_integrals, only:edge_int_3d
 real, intent(in) :: r0, R_0, d1, d2, h
 real(doub_prec) :: hi1, int1, int2

 if (abs(r0) < tiny(0.) .or. npieces_int == 0) then
    pint3D = 0.
    return
 endif
 pint3D = sign(1.,r0)*sign(1.,R_0)
 hi1 = 1.d0/h
 call edge_int_3d(abs(d1)*hi1,abs(d2)*hi1,abs(R_0)*hi1,abs(r0)*hi1,npieces_int,nmax_int,qb_int,a3_int,int1,int2)
 int1 = max(int1,0.d0)
 int2 = max(int2,0.d0)
 ! combine the two vertex integrals depending on where the foot of the side lies
 if (d1*d2 >= 0.) then
    pint3D = pint3D*real(int1 + int2)
 elseif (abs(d1) < abs(d2)) then
    pint3D = pint3D*real(max(int2 - int1,0.d0))
 else
    pint3D = pint3D*real(max(int1 - int2,0.d0))
 endif

end function pint3D

 !------------------------------------------------------------
 ! function to return a soft maximum for 1/x with no bias
 ! for x >> eps using the cubic spline kernel softening
 ! i.e. something equivalent to 1/sqrt(x**2 + eps**2) but
 ! with compact support, i.e. f=1/x when x > 2*eps
 !------------------------------------------------------------
pure elemental real function soft_func(x,eps) result(f)
 real, intent(in) :: x,eps
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
