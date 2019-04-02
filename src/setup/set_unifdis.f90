!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: unifdis
!
!  DESCRIPTION:
!   Setup of uniform particle distributions on various lattices
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: domain, part, random, stretchmap
!+
!--------------------------------------------------------------------------
module unifdis
 implicit none
 public :: set_unifdis, get_ny_nz_closepacked

 private

contains

!-------------------------------------------------------------
!+
!  This subroutine positions particles on a uniform lattice
!  in three dimensions, either cubic or close packed.
!+
!-------------------------------------------------------------
subroutine set_unifdis(lattice,id,master,xmin,xmax,ymin,ymax, &
                       zmin,zmax,delta,hfact,np,xyzh,rmin,rmax,rcylmin,rcylmax,&
                       nptot,npy,npz,rhofunc,inputiseed,verbose,dir,geom)
 use random,     only:ran2
 use part,       only:periodic
 use stretchmap, only:set_density_profile
 use domain,     only:i_belong
 character(len=*), intent(in)    :: lattice
 integer,          intent(in)    :: id,master
 integer,          intent(inout) :: np
 real,             intent(in)    :: xmin,xmax,ymin,ymax,zmin,zmax,delta,hfact
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(in),    optional :: rmin,rmax
 real,             intent(in),    optional :: rcylmin,rcylmax
 integer(kind=8),  intent(inout), optional :: nptot
 integer,          intent(in),    optional :: npy,npz,dir,geom
 real, external,                  optional :: rhofunc
 integer,          intent(in),    optional :: inputiseed
 logical,          intent(in),    optional :: verbose

 integer            :: i,j,k,l,m,nx,ny,nz,npnew,npin
 integer            :: jy,jz,ipart,maxp,iseed,icoord,igeom
 integer(kind=8)    :: iparttot
 real               :: delx,dely
 real               :: deltax,deltay,deltaz,dxbound,dybound,dzbound
 real               :: xstart,ystart,zstart,xi,yi,zi,rcyl2,rr2
 !real               :: xcentre,ycentre,zcentre
 real               :: rmin2,rmax2,rcylmin2,rcylmax2
 real               :: xpartmin,ypartmin,zpartmin
 real               :: xpartmax,ypartmax,zpartmax,xmins,xmaxs
 logical            :: is_verbose
 character(len=*), parameter :: fmt1 = "(/,1x,16('-'),' particles set on ',i3,2(' x ',i3),"// &
                                       "' uniform ',a,' lattice ',14('-'))"
 character(len=*), parameter :: fmt2 = "(/,1x,13('-'),' particles set on',i6,2(' x',i6),"// &
                                       "' uniform ',a,' lattice ',14('-'))"
 character(len=*), parameter :: fmtnx = "(1x,'enforcing periodicity: nx, ny, nz = ',i6,2(' x',i6),',  np = ',i8)"
 character(len=*), parameter :: fmtxx = "(3(2x,a,':',1pg12.3,'->',1pg11.3,' '))"
 character(len=*), parameter :: fmtdx = "(3(1x,a,':',1pg12.3,14x))"
 character(len=*), parameter :: fmtdy = "(28x,a,f8.4,14x,a,f8.4)"
 character(len=*), parameter :: sep   = "(/,1x,89('-'))"

 xpartmin = huge(0.)
 ypartmin = huge(0.)
 zpartmin = huge(0.)
 xpartmax = -huge(0.)
 ypartmax = -huge(0.)
 zpartmax = -huge(0.)

 dxbound = xmax-xmin
 dybound = ymax-ymin
 dzbound = zmax-zmin
 maxp = size(xyzh(1,:))
 npin = np

 if (present(rmin)) then
    rmin2 = rmin*rmin
 else
    rmin2 = 0.
 endif
 if (present(rmax)) then
    rmax2 = rmax*rmax
 else
    rmax2 = huge(0.)
 endif
 if (present(rcylmin)) then
    rcylmin2 = rcylmin*rcylmin
 else
    rcylmin2 = 0.
 endif
 if (present(rcylmax)) then
    rcylmax2 = rcylmax*rcylmax
 else
    rcylmax2 = huge(0.)
 endif
 if (present(nptot)) then
    iparttot = nptot
 else
    iparttot = 0
 endif

 ! Suppress output to the terminal if wished - handy for setups which call this subroutine frequently

 is_verbose = .true.
 if (present(verbose)) then
    is_verbose = verbose
 endif

 select case(trim(lattice))
 case('cubic')
    nx = nint(dxbound/delta)
    ny = nint(dybound/delta)
    nz = nint(dzbound/delta)
    deltaz = dzbound/nz
    deltay = dybound/ny
    deltax = dxbound/nx
    if (id==master .and. is_verbose) then
       if (max(nx,ny,nz) < 1e3) then
          print fmt1,nx,ny,nz,trim(lattice)
       else
          print fmt2,nx,ny,nz,trim(lattice)
       endif
       print fmtxx,'x',xmin,xmax,'y',ymin,ymax,'z',zmin,zmax
       print fmtdx,'dx',deltax,'dy',deltay,'dz',deltaz
    endif
    npnew=nx*ny*nz

    xstart = 0.
    ystart = 0.
    zstart = 0.

    ipart = np
    do k=1,nz
       zi = zmin + (k-0.5)*deltaz + zstart
       !print*,' z = ',zi
       do j=1,ny
          yi = ymin + (j-0.5)*deltay + ystart
          do i=1,nx
             xi = xmin + (i-0.5)*deltax + xstart

             rcyl2 = xi*xi + yi*yi
             rr2   = rcyl2 + zi*zi
             if (in_range(rr2,rmin2,rmax2) .and. in_range(rcyl2,rcylmin2,rcylmax2)) then
                iparttot = iparttot + 1
                if (i_belong(iparttot)) then
                   ipart = ipart + 1
                   if (ipart > maxp) stop 'ipart > maxp: re-compile with MAXP=bigger number'
                   xyzh(1,ipart) = xi
                   xyzh(2,ipart) = yi
                   xyzh(3,ipart) = zi
                   xyzh(4,ipart) = hfact*deltax
                endif
             endif
          enddo
       enddo
    enddo
    np = ipart
    if (present(nptot)) then
       nptot = iparttot
    endif

 case('hcp','hexagonal')
!
!--set uniform particle distribution, centred at the origin
!
    deltax = delta
    deltay = delta*sqrt(3./4.)
    deltaz = delta*sqrt(6.)/3.

    delx = 0.5*delta
    dely = 1./3.*deltay  !psep*sqrt(3.)/6.

    nx = int((1.-epsilon(0.))*dxbound/deltax) + 1
    ny = int((1.-epsilon(0.))*dybound/deltay) + 1
    nz = int((1.-epsilon(0.))*dzbound/deltaz) + 1
    npnew = nx*ny*nz

    if (id==master .and. is_verbose) then
       if (max(nx,ny,nz) < 1e3) then
          print fmt1,nx,ny,nz,trim(lattice)
       else
          print fmt2,nx,ny,nz,trim(lattice)
       endif
       print fmtxx,'x',xmin,xmax,'y',ymin,ymax,'z',zmin,zmax
       print fmtdx,'dx',deltax,'dy',deltay,'dz',deltaz
       print fmtdy,'dy/dx: ',deltay/deltax,' dz/dx: ',deltaz/deltax
    endif
    if (periodic) then
       ny = 2*int(ny/2)
       nz = 2*int(nz/2)
       deltax = dxbound/nx
       deltay = dybound/ny
       deltaz = dzbound/nz
       npnew = nx*ny*nz
       if (id==master .and. is_verbose) then
          print fmtnx,  nx,ny,nz,npnew
          print fmtdx,'dx',deltax,'dy',deltay,'dz',deltaz
          write(*,*) 'boundary check: (1. if lattice is periodic across boundaries)',deltay/(deltax*sqrt(3./4.)), &
                     deltaz/(deltax*sqrt(6.)/3.)
       endif
    endif

    k = 0
    l = 1
    m = 1
    ipart = np
    do i = 1, npnew
       k = k + 1
       if (k > nx) then
          k = 1
          l = l + 1
          if (l > ny) then
             l = 1
             m = m + 1
             if (m > nz) then
                k = 1
                l = 1
                nz = nz + 1
             endif
          endif
       endif

       xstart = xmin + 0.5*delx
       ystart = ymin + 0.5*dely
       zstart = zmin + 0.5*deltaz

       jy = mod(l, 2)
       jz = mod(m, 2)

       if (jz==0) then  ! 2nd layer
          ystart = ystart + dely
          if (jy==1) xstart = xstart + delx
       elseif (jy==0) then  ! first layer, jz=1
          xstart = xstart + delx
       endif

       xi = xstart + float(k - 1)*deltax
       yi = ystart + float(l - 1)*deltay
       zi = zstart + float(m - 1)*deltaz

       xpartmin = min(xpartmin,xi)
       ypartmin = min(ypartmin,yi)
       zpartmin = min(zpartmin,zi)
       xpartmax = max(xpartmax,xi)
       ypartmax = max(ypartmax,yi)
       zpartmax = max(zpartmax,zi)

       !
       !--trim to fit radius. do not allow particles to have *exactly* rmax
       !  (this stops round-off error from giving non-zero centre of mass)
       !
       rcyl2 = xi*xi + yi*yi
       rr2   = rcyl2 + zi*zi
       if (in_range(rr2,rmin2,rmax2) .and. in_range(rcyl2,rcylmin2,rcylmax2)) then
          iparttot = iparttot + 1
          if (i_belong(iparttot)) then
             ipart = ipart + 1
             if (ipart > maxp) stop 'ipart > maxp: re-compile with MAXP=bigger number'
             xyzh(1,ipart) = xi
             xyzh(2,ipart) = yi
             xyzh(3,ipart) = zi
             xyzh(4,ipart) = hfact*deltax
          endif
       endif

    enddo

    if (id==master .and. is_verbose) then
       print*,'part boundaries',xpartmin,xpartmax,ypartmin,ypartmax,zpartmin,zpartmax

       print*,'part spacing with the edges of the box ','x',(xpartmin-xmin)/deltax,(xmax-xpartmax)/deltax, &
           'y',(ypartmin-ymin)/deltay,(ymax-ypartmax)/deltay, &
           'z',(zpartmin-zmin)/deltaz,(zmax-zpartmax)/deltaz
    endif
    np = ipart
    if (present(nptot)) then
       nptot = iparttot
    endif

 case('closepacked')
!
!--set uniform particle distribution, centred at the origin
!
    !xcentre = 0.5*(xmin + xmax)
    !ycentre = 0.5*(ymin + ymax)
    !zcentre = 0.5*(zmin + zmax)
    ! xmin = -rmax
    ! ymin = -rmax
    ! zmin = -rmax
    ! xmax = rmax
    ! ymax = rmax
    ! zmax = rmax

    deltax = delta
    deltay = delta*sqrt(3./4.)
    deltaz = delta*sqrt(6.)/3.

    delx = 0.5*delta
    dely = 1./3.*deltay  !psep*sqrt(3.)/6.

    nx = int((1.-epsilon(0.))*dxbound/deltax) + 1
    ny = int((1.-epsilon(0.))*dybound/deltay) + 1
    nz = int((1.-epsilon(0.))*dzbound/deltaz) + 1
    if (present(npy)) then
       if (npy > 0) ny = npy
    endif
    if (present(npz)) then
       if (npz > 0) nz = npz
    endif

    npnew = nx*ny*nz

    if (id==master .and. is_verbose) then
       if (max(nx,ny,nz) < 1e3) then
          print fmt1,nx,ny,nz,trim(lattice)
       else
          print fmt2,nx,ny,nz,trim(lattice)
       endif
       print fmtxx,'x',xmin,xmax,'y',ymin,ymax,'z',zmin,zmax
       print fmtdx,'dx',deltax,'dy',deltay,'dz',deltaz
       print fmtdy,'dy/dx: ',deltay/deltax,' dz/dx: ',deltaz/deltax
    endif
    if (periodic) then
       ny = 2*int(ny/2)
       nz = 3*int(nz/3)
       if (.not.present(npy)) then  ! if fixing ny and nz keep deltax = input value
          deltax = dxbound/real(nx)
       endif
       deltay = dybound/real(ny)
       deltaz = dzbound/real(nz)
       dely = 1./3.*deltay
       npnew = nx*ny*nz

       if (id==master .and. is_verbose) then
          print fmtnx,nx,ny,nz,npnew
          print fmtdx,'dx',deltax,'dy',deltay,'dz',deltaz
          print fmtdy,'dy/dx: ',deltay/deltax,' dz/dx: ',deltaz/deltax
          !write(*,*) 'boundary check: (1. if lattice is periodic across boundaries)',deltay/(deltax*sqrt(3./4.)), &
          !        deltaz/(deltax*sqrt(6.)/3.)
       endif
    endif

    !
    !--set the limits so that the particles are
    !  exactly centred on the origin
    !
    ! xmin = xcentre - (nx-1)/2*psepx
    ! xmax = xcentre + (nx-1)/2*psepx
    ! ymin = ycentre - (ny-1)/2*psepy
    ! ymax = ycentre + (ny-1)/2*psepy
    ! zmin = zcentre - (nz-1)/2*psepz
    ! zmax = zcentre + (nz-1)/2*psepz

    k = 0
    l = 1
    m = 1
    ipart = np
    do i = 1, npnew
       k = k + 1
       if (k > nx) then
          k = 1
          l = l + 1
          if (l > ny) then
             l = 1
             m = m + 1
             if (m > nz) then
                k = 1
                l = 1
                nz = nz + 1
             endif
          endif
       endif

       !xstart = xcentre - 0.5*nx*deltax
       !ystart = ycentre - 0.5*ny*deltay
       !zstart = zcentre - 0.5*nz*deltaz

       !   xstart = xmin
       !   ystart = ymin + 0.5*dely
       !   zstart = zmin + 0.5*deltaz

       xstart = xmin + 0.5*delx
       ystart = ymin + 0.5*dely
       zstart = zmin + 0.5*deltaz

       jy = mod(l, 2)
       jz = mod(m, 3)
       !if (mod(m,2)==0) then
       !   if (jy==1) xstart = xstart + delx
       !else
       !   if (jy==0) xstart = xstart + delx
       !endif

       if (jz==0) then  ! 3rd layer
          ystart = ystart + 2.*dely
          if (jy==0) xstart = xstart + delx
       elseif (jz==2) then  ! 2nd layer
          ystart = ystart + dely
          if (jy==1) xstart = xstart + delx
       elseif (jy==0) then  ! first layer, jz=1
          xstart = xstart + delx
       endif

       xi = xstart + float(k - 1)*deltax
       yi = ystart + float(l - 1)*deltay
       zi = zstart + float(m - 1)*deltaz

       xpartmin = min(xpartmin,xi)
       ypartmin = min(ypartmin,yi)
       zpartmin = min(zpartmin,zi)
       xpartmax = max(xpartmax,xi)
       ypartmax = max(ypartmax,yi)
       zpartmax = max(zpartmax,zi)

       !
       !--trim to fit radius. do not allow particles to have *exactly* rmax
       !  (this stops round-off error from giving non-zero centre of mass)
       !
       rcyl2 = xi*xi + yi*yi
       rr2   = rcyl2 + zi*zi
       if (in_range(rr2,rmin2,rmax2) .and. in_range(rcyl2,rcylmin2,rcylmax2)) then
          iparttot = iparttot + 1
          if (i_belong(iparttot)) then
             ipart = ipart + 1
             if (ipart > maxp) stop 'ipart > maxp: re-compile with MAXP=bigger number'
             xyzh(1,ipart) = xi
             xyzh(2,ipart) = yi
             xyzh(3,ipart) = zi
             xyzh(4,ipart) = hfact*deltax
          endif
       endif

    enddo

    !if (id==master .and. periodic) then
    !print*,'part boundaries',xpartmin,xpartmax,ypartmin,ypartmax,zpartmin,zpartmax

    !print*,'part spacing with the edges of the box ','x',(xpartmin-xmin)/deltax,(xmax-xpartmax)/deltax, &
    !    'y',(ypartmin-ymin)/deltay,(ymax-ypartmax)/deltay, &
    !    'z',(zpartmin-zmin)/deltaz,(zmax-zpartmax)/deltaz
    !endif
    np = ipart
    if (present(nptot)) then
       nptot = iparttot
    endif

 case('random')
!
!--initialise random number generator
!
    if(present(inputiseed))then
       iseed = inputiseed
    else
       iseed = -43587
    endif

    if (id==master .and. is_verbose) then
       print sep
       write(*,*) 'random seed = ',iseed
    endif

    ipart = np
    nx = nint(dxbound/delta)
    ny = nint(dybound/delta)
    nz = nint(dzbound/delta)
    npnew = nx*ny*nz

    do i=1,npnew
       xi = xmin + ran2(iseed)*dxbound
       yi = ymin + ran2(iseed)*dybound
       zi = zmin + ran2(iseed)*dzbound
!--do not use if not within radial cuts
       rcyl2 = xi*xi + yi*yi
       rr2   = rcyl2 + zi*zi
       if (in_range(rr2,rmin2,rmax2) .and. in_range(rcyl2,rcylmin2,rcylmax2)) then
          iparttot = iparttot + 1
          if (i_belong(iparttot)) then
             ipart = ipart + 1
             if (ipart > maxp) stop 'ipart > maxp: re-compile with MAXP=bigger number'
             xyzh(1,ipart) = xi
             xyzh(2,ipart) = yi
             xyzh(3,ipart) = zi
             xyzh(4,ipart) = hfact*delta
          endif
       endif
    enddo
    np = ipart
    if (id==master .and. is_verbose) then
       print*,np, 'particles set in uniform random distribution'
    endif
    if (present(nptot)) then
       nptot = iparttot
    endif

 case default
    if (id==master) print*,'unknown lattice choice '''//trim(lattice)//''' in set_unifdis'
    stop
 end select
 if (id==master .and. is_verbose) print "(1x,89('-'),/)"

 !
 ! allow stretch mapping to give arbitrary non-uniform distributions along one coordinate direction
 !
 if (present(rhofunc)) then
    icoord = 1 ! default is x direction
    if (present(dir)) then
       if (dir >= 1 .and. dir <= 3) icoord = dir
    endif
    if (present(geom)) then
       igeom = geom
    else
       igeom = 1
    endif
    xmins = xmin
    xmaxs = xmax
    if (igeom==1) then
       if (icoord==3) then
          xmins = zmin
          xmaxs = zmax
       elseif (icoord==2) then
          xmins = ymin
          xmaxs = ymax
       endif
    endif
    call set_density_profile(np,xyzh,min=xmins,max=xmaxs,rhofunc=rhofunc,start=npin,geom=igeom,coord=icoord)
 endif

 return
end subroutine set_unifdis
!
!-------------------------------------------------------------
!+
!  Supplementary function & subroutine
!+
!-------------------------------------------------------------
pure logical function in_range(x,xmin,xmax)
 real, intent(in) :: x,xmin,xmax

 in_range = (x >= xmin .and. x <= xmax)

end function in_range

subroutine get_ny_nz_closepacked(delta,ymin,ymax,zmin,zmax,ny,nz)
 real,     intent(in) :: delta,ymin,ymax,zmin,zmax
 integer, intent(out) :: ny,nz
 real :: deltay,deltaz

 deltay = delta*sqrt(3./4.)
 deltaz = delta*sqrt(6.)/3.

 ny = int((1.-epsilon(0.))*(ymax-ymin)/deltay) + 1
 nz = int((1.-epsilon(0.))*(zmax-zmin)/deltaz) + 1

 ny = 2*int(ny/2)
 nz = 3*int(nz/3)

end subroutine get_ny_nz_closepacked

!-------------------------------------------------------------
end module unifdis
