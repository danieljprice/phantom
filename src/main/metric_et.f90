!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric
!
! Generic module for a tabulated metric, e.g. from the Einstein Toolkit
!
! :References: None
!
! :Owner: Spencer Magnall
!
! :Runtime parameters: None
!
! :Dependencies: metric_et_utils, table_utils
!
 implicit none
 character(len=*), parameter :: metric_type = 'et'
 integer,          parameter :: imetric     = 6
 ! This are dummy parameters to stop the compiler complaing
 ! Not used anywhere in the code - Needs a fix!
 real, public  :: mass1 = 1.       ! mass of central object
 real, public  :: a     = 0.0       ! spin of central object
contains

!----------------------------------------------------------------
!+
!  Compute the metric tensor in both covariant (gcov) and
!  contravariant (gcon) form. Here we merely interpolate
!  these values from the global grid.
!+
!----------------------------------------------------------------
pure subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
 use metric_et_utils, only:gridinit
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 integer :: ierr

 ! The subroutine that computes the metric tensor for a given position
 ! In this case it is interpolated from the global grid values
 ! Perform trilinear interpolation
 if ( .not. gridinit) then
    ierr = 1
    ! This is required for phantomsetup
    ! As no grid information has been passed to phantom from ET
    ! So interpolation cannot be performed
    if (ierr /= 0) then
       gcov = 0.
       gcov(0,0) = -1.
       gcov(1,1) = 1.
       gcov(2,2) = 1.
       gcov(3,3) = 1.
       if (present(gcon)) then
          gcon      = 0.
          gcon(0,0) = -1.
          gcon(1,1) = 1.
          gcon(2,2) = 1.
          gcon(3,3) = 1.
       endif
       if (present(sqrtg)) sqrtg   = -1.
    endif
 elseif (present(gcon) .and. present(sqrtg)) then
    call interpolate_metric(position,gcov,gcon,sqrtg)
 else
    call interpolate_metric(position,gcov)
 endif

end subroutine get_metric_cartesian

!-----------------------------------------------------------------------
!+
!  dummy routine to get the metric in spherical coordinates (not used)
!+
!-----------------------------------------------------------------------
pure subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional :: gcon(0:3,0:3)
 real, intent(out), optional :: sqrtg
 real :: r2,sintheta

 gcov = 0.

 r2       = position(1)**2
 sintheta = sin(position(2))

 gcov(0,0) = -1.
 gcov(1,1) = 1.
 gcov(2,2) = r2
 gcov(3,3) = r2*sintheta**2

 if (present(gcon)) then
    gcon      = 0.
    gcon(0,0) = -1.
    gcon(1,1) = 1.
    gcon(2,2) = 1./r2
    gcov(3,3) = 1./gcov(3,3)
 endif

 if (present(sqrtg)) sqrtg = r2*sintheta

end subroutine get_metric_spherical

!-----------------------------------------------------------------------
!+
!  cartesian metric derivatives, interpolates the derivatives from
!  the grid
!+
!-----------------------------------------------------------------------
pure subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 use metric_et_utils, only:gridinit
 real,    intent(in)  :: position(3)
 real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 integer :: ierr

 if (.not. gridinit) then
    ierr  = 1
    if (ierr /= 0) then
       dgcovdx = 0.
       dgcovdy = 0.
       dgcovdz = 0.
    else
       !    gridinit = .true.
       call interpolate_metric_derivs(position,dgcovdx,dgcovdy,dgcovdz)
    endif
 else
    call interpolate_metric_derivs(position,dgcovdx,dgcovdy,dgcovdz)
 endif

end subroutine metric_cartesian_derivatives

!-----------------------------------------------------------------------
!+
!  dummy routine for spherical metric derivatives, not used
!+
!-----------------------------------------------------------------------
pure subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi
 real :: r, theta

 r     = position(1)
 theta = position(2)

 dgcovdr     = 0.
 dgcovdtheta = 0.
 dgcovdphi   = 0.

 dgcovdr(2,2) = 2*r
 dgcovdr(3,3) = 2*r*sin(theta)**2

 dgcovdtheta(3,3) = 2*r**2*cos(theta)*sin(theta)

end subroutine metric_spherical_derivatives

!-----------------------------------------------------------------------
!+
!  dummy routine to convert cartesian to spherical coordinates
!+
!-----------------------------------------------------------------------
pure subroutine cartesian2spherical(xcart,xspher)
 real, intent(in)  :: xcart(3)
 real, intent(out) :: xspher(3)
 real :: x,y,z
 real :: r,theta,phi

 x  = xcart(1)
 y  = xcart(2)
 z  = xcart(3)

 r     = sqrt(x**2+y**2+z**2)
 theta = acos(z/r)
 phi   = atan2(y,x)

 xspher   = (/r,theta,phi/)

end subroutine cartesian2spherical

!-----------------------------------------------------------------------
!+
!  writes metric options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_metric(iunit)
 !use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 !write(iunit,"(/,a)") '# There are no options relating to the '//trim(metric_type)//' metric'
 !call write_inopt(metric_file,'metric_file','file from which to read tabulated metric (blank if used with einsteintk)',iunit)

end subroutine write_options_metric

!-----------------------------------------------------------------------
!+
!  reads metric options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_metric(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 !integer, save :: ngot = 0

 select case(trim(name))
    !case('metric_file')
    !   read(valstring,*,iostat=ierr) metric_file
    !   ngot = ngot + 1
 case default
    imatch = .false.
 end select
 !igotall = (ngot >= 1)
 igotall = .true.

end subroutine read_options_metric

!-----------------------------------------------------------------------
!+
! Interpolates value from grid to position
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_metric(position,gcov,gcon,sqrtg)
 ! linear and cubic interpolators should be moved to their own subroutine
 ! away from eos_shen
 use table_utils,     only:linear_interpolator_one_d
 use metric_et_utils, only:gcovgrid,gcongrid,sqrtggrid,dxgrid,gridorigin!,gridsize
 real, intent(in)  :: position(3)
 real, intent(out) :: gcov(0:3,0:3)
 real, intent(out), optional ::  gcon(0:3,0:3), sqrtg
 integer :: xlower,ylower,zlower!,xupper,yupper,zupper
 real    :: xlowerpos,ylowerpos,zlowerpos
 real :: xd,yd,zd
 real :: interptmp(7)
 integer :: i,j

 ! If the issue is that the metric vals are undefined on
 ! Setup since we have not recieved anything about the metric
 ! from ET during phantomsetup
 ! Then simply set gcov and gcon to 0
 ! as these values will be overwritten during the run anyway
 ! Get neighbours
 call get_grid_neighbours(position, dxgrid, xlower, ylower, zlower)

 xlowerpos = gridorigin(1) + (xlower-1)*dxgrid(1)
 ylowerpos = gridorigin(2) + (ylower-1)*dxgrid(2)
 zlowerpos = gridorigin(3) + (zlower-1)*dxgrid(3)

 xd = (position(1) - xlowerpos)/(dxgrid(1))
 yd = (position(2) - ylowerpos)/(dxgrid(2))
 zd = (position(3) - zlowerpos)/(dxgrid(3))

 interptmp = 0.
 ! All the interpolation should go into an interface, then you should just call trilinear_interp
 ! interpolate for gcov
 do i=0, 3
    do j=0, 3
       ! Interpolate along x
       call linear_interpolator_one_d(gcovgrid(i,j,xlower,ylower,zlower), &
                gcovgrid(i,j,xlower+1,ylower,zlower),xd,interptmp(1))
       call linear_interpolator_one_d(gcovgrid(i,j,xlower,ylower,zlower+1), &
                gcovgrid(i,j,xlower+1,ylower,zlower+1),xd,interptmp(2))
       call linear_interpolator_one_d(gcovgrid(i,j,xlower,ylower+1,zlower), &
                gcovgrid(i,j,xlower+1,ylower+1,zlower),xd,interptmp(3))
       call linear_interpolator_one_d(gcovgrid(i,j,xlower,ylower+1,zlower+1), &
                gcovgrid(i,j,xlower+1,ylower+1,zlower+1),xd,interptmp(4))
       ! Interpolate along y
       call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
       call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
       ! Interpolate along z
       call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

       gcov(i,j) = interptmp(7)
    enddo
 enddo

 if (present(gcon)) then
    ! interpolate for gcon
    do i=0, 3
       do j=0, 3
          ! Interpolate along x
          call linear_interpolator_one_d(gcongrid(i,j,xlower,ylower,zlower), &
                gcongrid(i,j,xlower+1,ylower,zlower),xd,interptmp(1))
          call linear_interpolator_one_d(gcongrid(i,j,xlower,ylower,zlower+1), &
                gcongrid(i,j,xlower+1,ylower,zlower+1),xd,interptmp(2))
          call linear_interpolator_one_d(gcongrid(i,j,xlower,ylower+1,zlower), &
                gcongrid(i,j,xlower+1,ylower+1,zlower),xd,interptmp(3))
          call linear_interpolator_one_d(gcongrid(i,j,xlower,ylower+1,zlower+1), &
                gcongrid(i,j,xlower+1,ylower+1,zlower+1),xd,interptmp(4))
          ! Interpolate along y
          call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
          call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
          ! Interpolate along z
          call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

          gcon(i,j) = interptmp(7)
       enddo
    enddo
 endif

 if (present(sqrtg)) then
    ! Interpolate for sqrtg
    ! Interpolate along x
    call linear_interpolator_one_d(sqrtggrid(xlower,ylower,zlower), &
                sqrtggrid(xlower+1,ylower,zlower),xd,interptmp(1))
    call linear_interpolator_one_d(sqrtggrid(xlower,ylower,zlower+1), &
                sqrtggrid(xlower+1,ylower,zlower+1),xd,interptmp(2))
    call linear_interpolator_one_d(sqrtggrid(xlower,ylower+1,zlower), &
                sqrtggrid(xlower+1,ylower+1,zlower),xd,interptmp(3))
    call linear_interpolator_one_d(sqrtggrid(xlower,ylower+1,zlower+1), &
                sqrtggrid(xlower+1,ylower+1,zlower+1),xd,interptmp(4))
    ! Interpolate along y
    call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
    call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
    ! Interpolate along z
    call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

    sqrtg = interptmp(7)
 endif

end subroutine interpolate_metric

!-----------------------------------------------------------------------
!+
! Interpolates derivatives of the metric from the grid to the position
!+
!-----------------------------------------------------------------------
pure subroutine interpolate_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz)
 use table_utils,     only:linear_interpolator_one_d
 use metric_et_utils, only:metricderivsgrid, dxgrid,gridorigin
 real, intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3),dgcovdz(0:3,0:3)
 real, intent(in)  :: position(3)
 integer :: xlower,ylower,zlower!,xupper,yupper,zupper
 real :: xd,yd,zd,xlowerpos, ylowerpos,zlowerpos
 real :: interptmp(7)
 integer :: i,j

 call get_grid_neighbours(position, dxgrid, xlower, ylower, zlower)

 xlowerpos = gridorigin(1) + (xlower-1)*dxgrid(1)
 ylowerpos = gridorigin(2) + (ylower-1)*dxgrid(2)
 zlowerpos = gridorigin(3) + (zlower-1)*dxgrid(3)

 xd = (position(1) - xlowerpos)/(dxgrid(1))
 yd = (position(2) - ylowerpos)/(dxgrid(2))
 zd = (position(3) - zlowerpos)/(dxgrid(3))

 interptmp = 0.

 ! Interpolate for dx
 do i=0, 3
    do j=0, 3
       ! Interpolate along x
       call linear_interpolator_one_d(metricderivsgrid(i,j,1,xlower,ylower,zlower), &
                metricderivsgrid(i,j,1,xlower+1,ylower,zlower),xd,interptmp(1))
       call linear_interpolator_one_d(metricderivsgrid(i,j,1,xlower,ylower,zlower+1), &
                metricderivsgrid(i,j,1,xlower+1,ylower,zlower+1),xd,interptmp(2))
       call linear_interpolator_one_d(metricderivsgrid(i,j,1,xlower,ylower+1,zlower), &
                metricderivsgrid(i,j,1,xlower+1,ylower+1,zlower),xd,interptmp(3))
       call linear_interpolator_one_d(metricderivsgrid(i,j,1,xlower,ylower+1,zlower+1), &
                metricderivsgrid(i,j,1,xlower+1,ylower+1,zlower+1),xd,interptmp(4))
       ! Interpolate along y
       call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
       call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
       ! Interpolate along z
       call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

       dgcovdx(i,j) = interptmp(7)
    enddo
 enddo
 ! Interpolate for dy
 do i=0, 3
    do j=0, 3
       ! Interpolate along x
       call linear_interpolator_one_d(metricderivsgrid(i,j,2,xlower,ylower,zlower), &
                metricderivsgrid(i,j,2,xlower+1,ylower,zlower),xd,interptmp(1))
       call linear_interpolator_one_d(metricderivsgrid(i,j,2,xlower,ylower,zlower+1), &
                metricderivsgrid(i,j,2,xlower+1,ylower,zlower+1),xd,interptmp(2))
       call linear_interpolator_one_d(metricderivsgrid(i,j,2,xlower,ylower+1,zlower), &
                metricderivsgrid(i,j,2,xlower+1,ylower+1,zlower),xd,interptmp(3))
       call linear_interpolator_one_d(metricderivsgrid(i,j,2,xlower,ylower+1,zlower+1), &
                metricderivsgrid(i,j,2,xlower+1,ylower+1,zlower+1),xd,interptmp(4))
       ! Interpolate along y
       call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
       call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
       ! Interpolate along z
       call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

       dgcovdy(i,j) = interptmp(7)
    enddo
 enddo

 ! Interpolate for dz
 do i=0, 3
    do j=0, 3
       ! Interpolate along x
       call linear_interpolator_one_d(metricderivsgrid(i,j,3,xlower,ylower,zlower), &
                metricderivsgrid(i,j,3,xlower+1,ylower,zlower),xd,interptmp(1))
       call linear_interpolator_one_d(metricderivsgrid(i,j,3,xlower,ylower,zlower+1), &
                metricderivsgrid(i,j,3,xlower+1,ylower,zlower+1),xd,interptmp(2))
       call linear_interpolator_one_d(metricderivsgrid(i,j,3,xlower,ylower+1,zlower), &
                metricderivsgrid(i,j,3,xlower+1,ylower+1,zlower),xd,interptmp(3))
       call linear_interpolator_one_d(metricderivsgrid(i,j,3,xlower,ylower+1,zlower+1), &
                metricderivsgrid(i,j,3,xlower+1,ylower+1,zlower+1),xd,interptmp(4))
       ! Interpolate along y
       call linear_interpolator_one_d(interptmp(1),interptmp(3),yd,interptmp(5))
       call linear_interpolator_one_d(interptmp(2),interptmp(4),yd,interptmp(6))
       ! Interpolate along z
       call linear_interpolator_one_d(interptmp(5),interptmp(6),zd,interptmp(7))

       dgcovdz(i,j) = interptmp(7)
    enddo
 enddo

end subroutine interpolate_metric_derivs

!-----------------------------------------------------------------------
!+
!  Utility routine to get the lower grid neighbours of a position
!+
!-----------------------------------------------------------------------
pure subroutine get_grid_neighbours(position,dx,xlower,ylower,zlower)
 use metric_et_utils, only:gridorigin
 real, intent(in) :: position(3)
 real, intent(in) :: dx(3)
 integer, intent(out) :: xlower,ylower,zlower

 ! Get the lower grid neighbours of the position
 ! If this is broken change from floor to int
 ! How are we handling the edge case of a particle being
 ! in exactly the same position as the grid?
 ! Hopefully having different grid sizes in each direction
 ! Doesn't break the lininterp
 xlower = floor((position(1)-gridorigin(1))/dx(1))
 ylower = floor((position(2)-gridorigin(2))/dx(2))
 zlower = floor((position(3)-gridorigin(3))/dx(3))

 ! +1 because fortran
 xlower = xlower + 1
 ylower = ylower + 1
 zlower = zlower + 1

end subroutine get_grid_neighbours

end module metric

