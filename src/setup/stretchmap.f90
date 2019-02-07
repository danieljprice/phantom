!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: stretchmap
!
!  DESCRIPTION:
!   This module implements stretch mapping to create one dimensional
!   density profiles from uniform arrangements of SPH particles
!
!   The implementation is quite general, allowing the density function
!   to be specified in any of the coordinate dimensions in either
!   cartesian, cylindrical, spherical or toroidal coordinates. It is a
!   generalisation of the spherical stretch map described in Herant (1994)
!   and of the cartesian version used in Price (2004)
!
!  REFERENCES:
!    Herant, M. (1994) "Dirty Tricks for SPH", MmSAI 65, 1013
!    Price (2004), "Magnetic fields in Astrophysics", PhD Thesis, University of Cambridge
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: geometry, io, physcon, table_utils
!+
!--------------------------------------------------------------------------
module stretchmap
 use physcon, only:pi
 implicit none
 public :: set_density_profile
 public :: get_mass_r

 integer, private :: ngrid = 1024 ! number of points used when integrating rho to get mass
 integer, parameter, private :: maxits = 100  ! max number of iterations
 integer, parameter, private :: maxits_nr = 30  ! max iterations with Newton-Raphson
 real,    parameter, private :: tol = 1.e-9  ! tolerance on iterations
 private

contains

!--------------------------------------------------------------------------
!+
!  Subroutine to implement the stretch mapping procedure
!
!  IN/OUT:
!    xyzh : particle coordinates and smoothing length
!
!  IN:
!    np : number of particles
!    cmin, cmax         : range in the coordinate to apply transformation
!    start (optional)   : only consider particles between start and np
!    geom  (optional)   : geometry in which stretch mapping is to be performed
!                          1 - cartesian
!                          2 - cylindrical
!                          3 - spherical
!                          4 - toroidal
!                         (if not specified, assumed to be cartesian)
!    coord (optional)   : coordinate direction in which stretch mapping is to be performed
!                         (if not specified, assumed to be the first coordinate)
!    rhofunc (optional) : function containing the desired density function rho(r) or rho(x)
!    rhotab  (optional) : array of tabulated density profile
!    ctab    (optional) : array of tabulated coordinate positions for density bins in table
!
!  The function rhofunc is assumed to be a real function with a single argument:
!
!  real function rhofunc(r)
!   real, intent(in) :: r
!
!   rhofunc = 1./r**2
!
!  end function rhofunc
!
!  If the ctab array is not present, the table rhotab(:) is assumed to
!  contain density values on equally spaced bins between cmin and cmax
!+
!----------------------------------------------------------------
subroutine set_density_profile(np,xyzh,min,max,rhofunc,rhotab,xtab,start,geom,coord)
 use io,          only:error,fatal,iprint,id,master
 use geometry,    only:coord_transform,maxcoordsys,labelcoord,igeom_cartesian!,labelcoordsys
 use table_utils, only:yinterp,linspace
 integer, intent(in)    :: np
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: min,max
 real,    external,   optional :: rhofunc
 real,    intent(in), optional :: rhotab(:),xtab(:)
 integer, intent(in), optional :: start, geom, coord
 real :: totmass,rhozero,hi,fracmassold
 real :: x(3),xt(3),xmin,xmax,xold,xi,xminbisect,xmaxbisect
 real :: xprev,func,dfunc,rhoi,rho_at_min
 real, allocatable  :: xtable(:),masstab(:)
 integer            :: i,its,igeom,icoord,istart,ierr,nt
 logical            :: is_r, is_rcyl, bisect

 if (present(rhofunc) .or. present(rhotab)) then
    if (id==master) print "(a)",' >>>>>>  s  t  r  e   t    c     h       m     a    p   p  i  n  g  <<<<<<'
    !
    ! defaults for optional arguments
    !
    icoord = 1
    if (present(coord)) then
       if (coord >= 1 .and. coord <= 3) icoord = coord
    endif
    igeom = 1     ! default geometry is cartesian
    if (present(geom)) then
       if (geom >= 1 .and. geom <= maxcoordsys) igeom = geom
    endif
    !print "(a)",' stretching in '//trim(labelcoordsys(igeom))//' coordinate system'
    is_r    = is_rspherical(igeom,icoord)
    is_rcyl = is_rcylindrical(igeom,icoord)
    xmin    = min
    xmax    = max
    !
    !--get total mass integrated along coordinate direction to normalise the density profile
    !  (we assume total mass is 1 for both particles and the desired profile,
    !   i.e. that the same total mass is in both)
    !
    if (present(rhotab)) then
       if (id==master) write(iprint,*) 'stretching to match tabulated density profile in '//&
                                       trim(labelcoord(icoord,igeom))//' direction'
       nt = size(rhotab)
       if (nt <= 0) call fatal('set_density_profile','size of density table <= 0')
       allocate(xtable(nt),masstab(nt),stat=ierr)
       if (ierr /= 0) call fatal('set_density_profile','cannot allocate memory for mass/coordinate table')
       !
       ! if no coordinate table is passed, create a table
       ! of equally spaced points between xmin and xmax
       !
       if (present(xtab)) then
          if (size(xtab) < nt) call fatal('set_density_profile','coordinate table different size to density table')
          xtable(1:nt) = xtab(1:nt)
       else
          call linspace(xtable(1:nt),xmin,xmax)
       endif
       if (is_r) then
          call get_mass_tab_r(masstab,rhotab,xtable)
       elseif (is_rcyl) then
          call get_mass_tab_rcyl(masstab,rhotab,xtable)
       else
          call get_mass_tab(masstab,rhotab,xtable)
       endif
       totmass    = yinterp(masstab,xtable(1:nt),xmax)
       rho_at_min = yinterp(rhotab,xtable(1:nt),xmin)
    else
       if (id==master) write(iprint,*) 'stretching to match density profile in '&
                       //trim(labelcoord(icoord,igeom))//' direction'
       if (is_r) then
          totmass = get_mass_r(rhofunc,xmax,xmin)
       elseif (is_rcyl) then
          totmass = get_mass_rcyl(rhofunc,xmax,xmin)
       else
          totmass = get_mass(rhofunc,xmax,xmin)
       endif
       rho_at_min = rhofunc(xmin)
       nt = 1 ! not used, to avoid compiler warnings
    endif
    if (is_r) then
       rhozero = totmass/(4./3.*pi*(xmax**3 - xmin**3))
    elseif (is_rcyl) then
       rhozero = totmass/(pi*(xmax**2 - xmin**2))
    else
       rhozero = totmass/(xmax - xmin)
    endif

    if (id==master) then
       write(iprint,*) 'density at '//trim(labelcoord(icoord,igeom))//' = ',xmin,' is ',rho_at_min
       write(iprint,*) 'total mass      = ',totmass
    endif

    if (present(start)) then
       istart = start
    else
       istart = 0
    endif

    !$omp parallel do default(none) &
    !$omp shared(np,xyzh,rhozero,igeom,rhotab,xtable,masstab,nt) &
    !$omp shared(xmin,xmax,totmass,icoord,is_r,is_rcyl,istart) &
    !$omp private(x,xold,xt,fracmassold,its,xprev,xi,hi,rhoi) &
    !$omp private(func,dfunc,xminbisect,xmaxbisect,bisect)
    do i=istart+1,np
       x(1) = xyzh(1,i)
       x(2) = xyzh(2,i)
       x(3) = xyzh(3,i)
       hi = xyzh(4,i)
       if (igeom > 1) then
          call coord_transform(x,3,igeom_cartesian,xt,3,igeom)
          xold = xt(icoord)
       else
          xt = x
          xold = x(icoord)
       endif
       if (is_r) then
          fracmassold = 4./3.*pi*rhozero*(xold**3 - xmin**3)
       elseif (is_rcyl) then
          fracmassold = pi*rhozero*(xold**2 - xmin**2)
       else
          fracmassold = rhozero*(xold - xmin)
       endif
       if (xold > xmin) then  ! if x=0 do nothing
          its   = 0
          xprev = 0.
          xi    = xold   ! starting guess
          ! calc func to determine if tol is met
          if (present(rhotab)) then
             func  = yinterp(masstab,xtable(1:nt),xi)
          else
             if (is_r) then
                func  = get_mass_r(rhofunc,xi,xmin)
             elseif (is_rcyl) then
                func  = get_mass_rcyl(rhofunc,xi,xmin)
             else
                func  = get_mass(rhofunc,xi,xmin)
             endif
          endif
          func = func - fracmassold
          xminbisect = xmin
          xmaxbisect = xmax
          bisect = .false.  ! use Newton-Raphson by default
          do while ((abs(func/totmass) > tol .and. its < maxits))
             xprev = xi
             its   = its + 1
             if (present(rhotab)) then
                func  = yinterp(masstab,xtable(1:nt),xi) - fracmassold
                if (is_r) then
                   dfunc = 4.*pi*xi**2*yinterp(rhotab,xtable(1:nt),xi)
                elseif (is_rcyl) then
                   dfunc = 2.*pi*xi*yinterp(rhotab,xtable(1:nt),xi)
                else
                   dfunc = yinterp(rhotab,xtable(1:nt),xi)
                endif
             else
                if (is_r) then
                   func  = get_mass_r(rhofunc,xi,xmin) - fracmassold
                   dfunc = 4.*pi*xi**2*rhofunc(xi)
                elseif (is_rcyl) then
                   func  = get_mass_rcyl(rhofunc,xi,xmin) - fracmassold
                   dfunc = 2.*pi*xi*rhofunc(xi)
                else
                   func  = get_mass(rhofunc,xi,xmin) - fracmassold
                   dfunc = rhofunc(xi)
                endif
             endif

             if (bisect) then
                if (func > 0.) then
                   xmaxbisect = xi
                else
                   xminbisect = xi
                endif
                xi = 0.5*(xminbisect + xmaxbisect)
                !print*,i,its,' bisect ',xi,xprev,xold,func
             else
                ! Newton-Raphson
                xi = xprev - func/dfunc
                ! do not allow N-R to take big jumps
                if (abs(xi) < 0.8*abs(xprev)) then
                   xi = 0.8*xprev
                elseif (abs(xi) > 1.2*abs(xprev)) then
                   xi = 1.2*xprev
                endif
                !if (its > maxits_nr) print*,i,its,'n-r',xi,xprev,xold,func

                ! Use bisection if Newton-Raphson is failing
                if (xi > xmax .or. xi < xmin .or. its > maxits_nr) then
                   bisect = .true.
                   xi = 0.5*(xminbisect + xmaxbisect)
                endif
             endif
          enddo
          if (present(rhotab)) then
             rhoi = yinterp(rhotab,xtable(1:nt),xi)
          else
             rhoi = rhofunc(xi)
          endif
          xt(icoord) = xi
          call coord_transform(xt,3,igeom,x,3,igeom_cartesian)
          xyzh(1,i) = x(1)
          xyzh(2,i) = x(2)
          xyzh(3,i) = x(3)
          xyzh(4,i) = hi*(rhozero/rhoi)**(1./3.)
          if (its >= maxits) call error('set_density_profile','Stretch mapping not converged')
       endif
    enddo
    !$omp end parallel do

    if (allocated(xtable))  deallocate(xtable)
    if (allocated(masstab)) deallocate(masstab)
    if (id==master) print "(a,/)",' >>>>>> done'
 endif

end subroutine set_density_profile

!--------------------------------------------------------------
!+
!  query function for whether we have spherical r direction
!+
!--------------------------------------------------------------
logical function is_rspherical(igeom,icoord)
 use geometry, only:igeom_spherical
 integer, intent(in) :: igeom,icoord

 is_rspherical = .false.
 if (igeom==igeom_spherical .and. icoord==1) is_rspherical = .true.

end function is_rspherical

!--------------------------------------------------------------
!+
!  query function for whether we have cylindrical r direction
!+
!--------------------------------------------------------------
logical function is_rcylindrical(igeom,icoord)
 use geometry, only:igeom_cylindrical,igeom_toroidal
 integer, intent(in) :: igeom,icoord

 is_rcylindrical = .false.
 if (igeom==igeom_cylindrical .or. igeom==igeom_toroidal .and. icoord==1) is_rcylindrical = .true.

end function is_rcylindrical

!--------------------------------------------------------------
!+
!  Integrate to get total mass along the coordinate direction
!+
!--------------------------------------------------------------
!
! mass integrated along spherical radius
!
real function get_mass_r(rhofunc,r,rmin)
 real, intent(in) :: r,rmin
 real, external   :: rhofunc
 real :: dr,ri,dmi,dmprev
 integer :: i

 dr = (r - rmin)/real(ngrid)
 dmprev     = 0.
 get_mass_r = 0.
 do i=1,ngrid
    ri         = rmin + i*dr
    dmi        = ri*ri*rhofunc(ri)*dr
    get_mass_r = get_mass_r + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev     = dmi
 enddo
 get_mass_r = 4.*pi*get_mass_r

end function get_mass_r
!
! mass integrated along cylindrical radius
!
real function get_mass_rcyl(rhofunc,rcyl,rmin)
 real, intent(in) :: rcyl,rmin
 real, external   :: rhofunc
 real :: dr,ri,dmi,dmprev
 integer :: i

 dr = (rcyl - rmin)/real(ngrid)
 dmprev   = 0.
 get_mass_rcyl = 0.
 do i=1,ngrid
    ri            = rmin + i*dr
    dmi           = ri*rhofunc(ri)*dr
    get_mass_rcyl = get_mass_rcyl + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev        = dmi
 enddo
 get_mass_rcyl = 2.*pi*get_mass_rcyl

end function get_mass_rcyl
!
! mass integrated along cartesian direction
!
real function get_mass(rhofunc,x,xmin)
 real, intent(in) :: x,xmin
 real, external   :: rhofunc
 real :: dx,xi,dmi,dmprev
 integer :: i

 dx = (x - xmin)/real(ngrid)
 dmprev   = 0.
 get_mass = 0.
 do i=1,ngrid
    xi       = xmin + i*dx
    dmi      = rhofunc(xi)*dx
    get_mass = get_mass + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev   = dmi
 enddo

end function get_mass

!------------------------------------
!+
!  Same as above, but fills a table
!+
!------------------------------------
!
! version that integrates along spherical radius
!
subroutine get_mass_tab_r(masstab,rhotab,rtab)
 real, intent(in)  :: rhotab(:),rtab(:)
 real, intent(out) :: masstab(size(rhotab))
 real :: dr,ri,dmi,dmprev,rprev
 integer :: i

 masstab(1) = 0.
 dmprev     = 0.
 rprev      = rtab(1)
 do i=2,size(rhotab)
    ri     = rtab(i)
    dr     = ri - rprev
    dmi    = ri*ri*rhotab(i)*dr
    masstab(i) = masstab(i-1) + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev = dmi
    rprev  = ri
 enddo
 masstab(:) = 4.*pi*masstab(:)

end subroutine get_mass_tab_r
!
! version that integrates along cylindrical radius
!
subroutine get_mass_tab_rcyl(masstab,rhotab,rtab)
 real, intent(in)  :: rhotab(:),rtab(:)
 real, intent(out) :: masstab(size(rhotab))
 real :: dr,ri,dmi,dmprev,rprev
 integer :: i

 masstab(1) = 0.
 dmprev     = 0.
 rprev      = rtab(1)
 do i=2,size(rhotab)
    ri     = rtab(i)
    dr     = ri - rprev
    dmi    = ri*rhotab(i)*dr
    masstab(i) = masstab(i-1) + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev = dmi
    rprev  = ri
 enddo
 masstab(:) = 2.*pi*masstab(:)

end subroutine get_mass_tab_rcyl
!
! version that integrates along a cartesian direction
!
subroutine get_mass_tab(masstab,rhotab,xtab)
 real, intent(in)  :: rhotab(:),xtab(:)
 real, intent(out) :: masstab(size(rhotab))
 real :: dx,xi,dmi,dmprev,xprev
 integer :: i

 masstab(1) = 0.
 dmprev     = 0.
 xprev      = xtab(1)
 do i=2,size(rhotab)
    xi     = xtab(i)
    dx     = xi - xprev
    dmi    = rhotab(i)*dx
    masstab(i) = masstab(i-1) + 0.5*(dmi + dmprev) ! trapezoidal rule
    dmprev = dmi
    xprev  = xi
 enddo

end subroutine get_mass_tab

end module stretchmap
