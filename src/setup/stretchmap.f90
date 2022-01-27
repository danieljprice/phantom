!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module stretchmap
!
! This module implements stretch mapping to create one dimensional
! density profiles from uniform arrangements of SPH particles
!
! The implementation is quite general, allowing the density function
! to be specified in any of the coordinate dimensions in either
! cartesian, cylindrical, spherical or toroidal coordinates. It is a
! generalisation of the spherical stretch map described in Herant (1994)
! and of the cartesian version used in Price (2004)
!
! :References:
!    - Herant, M. (1994) "Dirty Tricks for SPH", MmSAI 65, 1013
!    - Price (2004), "Magnetic fields in Astrophysics", PhD Thesis, University of Cambridge
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: geometry, table_utils
!
 implicit none
 real, parameter :: pi = 4.*atan(1.) ! the circle of life
 public :: set_density_profile
 public :: get_mass_r
 public :: rho_func

 integer, private :: ngrid = 1024 ! number of points used when integrating rho to get mass
 integer, parameter, private :: maxits = 100  ! max number of iterations
 integer, parameter, private :: maxits_nr = 30  ! max iterations with Newton-Raphson
 real,    parameter, private :: tol = 1.e-9  ! tolerance on iterations
 integer, parameter, public :: ierr_zero_size_density_table = 1, & ! error code
                               ierr_memory_allocation = 2, & ! error code
                               ierr_table_size_differs = 3, & ! error code
                               ierr_not_converged = -1 ! error code
 abstract interface
  real function rho_func(x)
   real, intent(in) :: x
  end function rho_func
 end interface

 private

contains

subroutine set_density_profile(np,xyzh,min,max,rhofunc,rhotab,xtab,start,geom,coord,verbose,err)
!
!  Subroutine to implement the stretch mapping procedure
!
!  The function rhofunc is assumed to be a real function with a single argument::
!
!     real function rho(r)
!      real, intent(in) :: r
!
!       rho = 1./r**2
!
!     end function rho
!
!  If the ctab array is not present, the table rhotab(:) is assumed to
!  contain density values on equally spaced bins between cmin and cmax
!
!    xyzh    : particle coordinates and smoothing length
!    np      : number of particles
!    min_bn  : min range in the coordinate to apply transformation
!    max_bn  : max range in the coordinate to apply transformation
!    start   : only consider particles between start and np (optional)
!    geom    : geometry in which stretch mapping is to be performed (optional)
!                      1 - cartesian
!                      2 - cylindrical
!                      3 - spherical
!                      4 - toroidal
!                     (if not specified, assumed to be cartesian)
!    coord   : coordinate direction in which stretch mapping is to be performed (optional)
!              (if not specified, assumed to be the first coordinate)
!    rhofunc : function containing the desired density function rho(r) or rho(x) (optional)
!    xtab    : tabulated coordinate values (optional)
!    rhotab  : tabulated density profile (optional)
!    ctab    : tabulated coordinate positions for density bins in table (optional)
!    verbose : turn on/off verbose output (optional)
!    err     : error code (0 on successful run)
!
 use geometry,    only:coord_transform,maxcoordsys,labelcoord,igeom_cartesian!,labelcoordsys
 use table_utils, only:yinterp,linspace
 integer, intent(in)    :: np
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: min,max
 procedure(rho_func), pointer, optional :: rhofunc
 real,    intent(in), optional :: rhotab(:),xtab(:)
 integer, intent(in), optional :: start, geom, coord
 logical, intent(in), optional :: verbose
 integer, intent(out),optional :: err
 real :: totmass,rhozero,hi,fracmassold
 real :: x(3),xt(3),xmin,xmax,xold,xi,xminbisect,xmaxbisect
 real :: xprev,func,dfunc,rhoi,rho_at_min
 real, allocatable  :: xtable(:),masstab(:)
 integer            :: i,its,igeom,icoord,istart,nt,nerr,ierr
 logical            :: is_r, is_rcyl, bisect, isverbose
 logical            :: use_rhotab

 isverbose = .true.
 use_rhotab = .false.
 if (present(verbose)) isverbose = verbose
 if (present(rhotab)) use_rhotab = .true.

 if (present(rhofunc) .or. present(rhotab)) then
    if (isverbose) print "(a)",' >>>>>>  s  t  r  e   t    c     h       m     a    p   p  i  n  g  <<<<<<'
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
       if (isverbose) write(*,*) 'stretching to match tabulated density profile in '//&
                                 trim(labelcoord(icoord,igeom))//' direction'
       nt = size(rhotab)
       if (nt <= 0) then
          if (isverbose) write(*,*) 'ERROR: zero size density table'
          if (present(err)) err = ierr_zero_size_density_table
          return
       endif
       allocate(xtable(nt),masstab(nt),stat=ierr)
       if (ierr /= 0) then
          if (isverbose) write(*,*) 'ERROR: cannot allocate memory for mass/coordinate table'
          if (present(err)) err = ierr_memory_allocation
          return
       endif
       !
       ! if no coordinate table is passed, create a table
       ! of equally spaced points between xmin and xmax
       !
       if (present(xtab)) then
          if (size(xtab) < nt) then
             if (isverbose) write(*,*) 'ERROR: coordinate table different size to density table'
             if (present(err)) err = ierr_table_size_differs
             return
          endif
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
       if (isverbose) write(*,*) 'stretching to match density profile in '&
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

    if (isverbose) then
       write(*,*) 'density at '//trim(labelcoord(icoord,igeom))//' = ',xmin,' is ',rho_at_min
       write(*,*) 'total mass      = ',totmass
    endif

    if (present(start)) then
       istart = start
    else
       istart = 0
    endif

    nerr = 0
    !$omp parallel do default(none) &
    !$omp shared(np,xyzh,rhozero,igeom,use_rhotab,rhotab,xtable,masstab,nt) &
    !$omp shared(xmin,xmax,totmass,icoord,is_r,is_rcyl,istart,rhofunc) &
    !$omp private(x,xold,xt,fracmassold,its,xprev,xi,hi,rhoi) &
    !$omp private(func,dfunc,xminbisect,xmaxbisect,bisect) &
    !$omp reduction(+:nerr)
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
          if (use_rhotab) then
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
             if (use_rhotab) then
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
          if (use_rhotab) then
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
          if (its >= maxits) nerr = nerr + 1
       endif
    enddo
    !$omp end parallel do

    if (isverbose .and. nerr > 0) then
       if (present(err)) err = ierr_not_converged
       if (isverbose) write(*,*) 'ERROR: stretch mapping not converged on ',nerr,' particles'
    endif

    if (allocated(xtable))  deallocate(xtable)
    if (allocated(masstab)) deallocate(masstab)
    if (isverbose) print "(a,/)",' >>>>>> done'
 endif

end subroutine set_density_profile

logical function is_rspherical(igeom,icoord)
!
! query function for whether we have spherical r direction
!
 use geometry, only:igeom_spherical
 integer, intent(in) :: igeom,icoord

 is_rspherical = .false.
 if (igeom==igeom_spherical .and. icoord==1) is_rspherical = .true.

end function is_rspherical

logical function is_rcylindrical(igeom,icoord)
!
! query function for whether we have cylindrical r direction
!
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
real function get_mass_r(rhofunc,r,rmin)
!
! mass integrated along spherical radius
!
 real, intent(in) :: r,rmin
 procedure(rho_func), pointer :: rhofunc
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

real function get_mass_rcyl(rhofunc,rcyl,rmin)
!
! mass integrated along cylindrical radius
!
 real, intent(in) :: rcyl,rmin
 procedure(rho_func), pointer :: rhofunc
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

real function get_mass(rhofunc,x,xmin)
!
! mass integrated along cartesian direction
!
 real, intent(in) :: x,xmin
 procedure(rho_func), pointer :: rhofunc
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
subroutine get_mass_tab_r(masstab,rhotab,rtab)
!
! version that integrates along spherical radius
!
 real, intent(in)  :: rhotab(:),rtab(:)
 real, intent(out) :: masstab(size(rhotab))
 real :: dmi
 integer :: i

 masstab(1) = 0.
 do i=2,size(rhotab)
    dmi = 4./3. * pi * (rtab(i)**3 - rtab(i-1)**3) * rhotab(i)
    masstab(i) = masstab(i-1) + dmi
 enddo

end subroutine get_mass_tab_r

subroutine get_mass_tab_rcyl(masstab,rhotab,rtab)
!
! version that integrates along cylindrical radius
!
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

subroutine get_mass_tab(masstab,rhotab,xtab)
!
! version that integrates along a cartesian direction
!
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
