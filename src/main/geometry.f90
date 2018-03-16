!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: geometry
!
!  DESCRIPTION: See below
!
!  REFERENCES: Appendix B of SPLASH user guide
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

!-----------------------------------------------------------------
! Standalone module containing subroutines to transform between
! different co-ordinate systems, for co-ordinates and vectors
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! itype must be one of the following:
!  itype = 1    : cartesian (default)
!  itype = 2    : cylindrical
!  itype = 3    : spherical
!  itype = 4    : toroidal
!  itype = 5    : rotated cartesian
!
! Currently handles:
!
!  cartesian -> cylindrical, spherical polar
!  cylindrical <-> cartesian
!  spherical polar <-> cartesian
!  toroidal r,theta,phi <-> cartesian
!  rotated cartesian <-> cartesian
!
! written by Daniel Price 2004-2011
! as part of the SPLASH SPH visualisation package
!-----------------------------------------------------------------
module geometry
 implicit none
 integer, parameter, public :: maxcoordsys = 5
 integer, parameter, public :: igeom_cartesian   = 1
 integer, parameter, public :: igeom_cylindrical = 2
 integer, parameter, public :: igeom_spherical   = 3
 integer, parameter, public :: igeom_toroidal    = 4
 integer, parameter, public :: igeom_rotated     = 5

 character(len=24), dimension(maxcoordsys), parameter, public :: labelcoordsys = &
    (/'cartesian   x,y,z      ', &
      'cylindrical r,phi,z    ', &
      'spherical   r,phi,theta', &
      'toroidal    r,theta,phi', &
      'rotated     x_1,x_2,x_3'/)
 character(len=6), dimension(3,maxcoordsys), parameter, public :: labelcoord = &
    reshape((/'x    ','y    ','z    ', &
              'r    ','phi  ','z    ', &
              'r    ','phi  ','theta', &
              'r_t  ','theta','phi  ', &
              'x_1  ','x_2  ','x_3  '/),shape=(/3,maxcoordsys/))

 public :: coord_transform, vector_transform, coord_transform_limits
 public :: coord_is_length, print_error
 public :: set_rotation_angles

 real, parameter, private :: pi = 3.1415926536
 real, parameter, private :: Rtorus = 1.0
 real, private :: sina = 2./3.
 real, private :: cosa = sqrt(5.)/3.
 real, private :: sinb = 2./sqrt(5.)
 real, private :: cosb = 1./sqrt(5.)

 integer, parameter, public :: ierr_invalid_dimsin  = 1
 integer, parameter, public :: ierr_invalid_dimsout = 2
 integer, parameter, public :: ierr_invalid_dims    = 3
 integer, parameter, public :: ierr_r_is_zero       = 4
 integer, parameter, public :: ierr_warning_assuming_cartesian = -1

 private

contains
!-----------------------------------------------------------------
! utility that returns whether or not a particular coordinate
! in a given coordinate system has dimensions of length or not
!-----------------------------------------------------------------
pure logical function coord_is_length(ix,igeom)
 integer, intent(in) :: ix,igeom

 coord_is_length = .false.
 select case(igeom)
 case(igeom_toroidal, igeom_spherical)
    if (ix==1) coord_is_length = .true.
 case(igeom_cylindrical)
    if (ix==1 .or. ix==3) coord_is_length = .true.
 case default
    coord_is_length = .true.
 end select

end function coord_is_length

!--------------------------------------------------------
! utility to handle error printing so transform routines
! do not generate verbose output
!--------------------------------------------------------

subroutine print_error(ierr)
 integer, intent(in) :: ierr

 select case(ierr)
 case(ierr_invalid_dimsin)
    print*,'Error: coord transform: invalid number of dimensions on input'
 case(ierr_invalid_dimsout)
    print*,'Error: coord transform: invalid number of dimensions on output'
 case(ierr_invalid_dims)
    print*,'Error: coord transform: ndimout must be <= ndimin'
 case(ierr_r_is_zero)
    print*,'Warning: coord transform: r=0 on input: cannot return angle'
 case(ierr_warning_assuming_cartesian)
    print*,'warning: using default cartesian output'
 case default
    print*,' unknown error'
 end select

end subroutine print_error

!-------------------------------------------------------------
! utility to set rotation angles for rotated cartesian system
!-------------------------------------------------------------
subroutine set_rotation_angles(a,b,sin_a,sin_b,cos_a,cos_b)
 real, intent(in), optional :: a,b,sin_a,sin_b,cos_a,cos_b

 if (present(a)) then
    sina = sin(a)
    cosa = cos(a)
 endif
 if (present(b)) then
    sinb = sin(b)
    cosb = cos(b)
 endif
 if (present(sin_a)) then
    sina = sin_a
    cosa = sqrt(1. - sina**2)
 endif
 if (present(sin_b)) then
    sinb = sin_b
    cosb = sqrt(1. - sinb**2)
 endif
 if (present(cos_a)) then
    cosa = cos_a
    if (.not.present(sin_a)) sina = sqrt(1. - cosa**2)
 endif
 if (present(cos_b)) then
    cosb = cos_b
    if (.not.present(cos_b)) sinb = sqrt(1. - cosb**2)
 endif

end subroutine set_rotation_angles

!-----------------------------------------------------------------
! Subroutine to transform between different co-ordinate systems
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! xin(ndimin)   : input co-ordinates, in ndimin dimensions
! itypein       : input co-ordinate type
!
! xout(ndimout) : output co-ordinates, in ndimout dimensions
! itypeout      : output co-ordinate type
!

!
!-----------------------------------------------------------------
pure subroutine coord_transform(xin,ndimin,itypein,xout,ndimout,itypeout,err)
 integer, intent(in)  :: ndimin,ndimout,itypein,itypeout
 real,    intent(in)  :: xin(ndimin)
 real,    intent(out) :: xout(ndimout)
 integer, intent(out), optional :: err
 real    :: rcyl,xi(3),xouti(3)
 integer :: ierr
!
!--check for errors in input
!
 ierr = 0
 if (itypeout==itypein) then
    xout(1:ndimout) = xin(1:ndimout)
    return
 elseif (ndimin < 1.or.ndimin > 3) then
    ierr = ierr_invalid_dimsin
    if (present(err)) err = ierr
    return
 elseif (ndimout < 1.or.ndimout > 3) then
    ierr = ierr_invalid_dimsout
    if (present(err)) err = ierr
    return
 elseif (ndimout > ndimin) then
    ierr = ierr_invalid_dims
    if (present(err)) err = ierr
    return
 elseif (abs(xin(1)) < 1e-8 .and. ndimout >= 2 .and. &
       (itypein==2 .or. itypein==3)) then
    ierr = ierr_r_is_zero
    if (present(err)) err = ierr
    xout(1:ndimout) = xin(1:ndimout)
    return
 endif
!
!--now do transformation
!
 select case(itypein)
!
!--input is cylindrical polars
!
 case(2)
    !
    ! output is cartesian (default)
    !
    if (itypeout /= 1) ierr = ierr_warning_assuming_cartesian
    if (ndimout==1) then
       xout(1) = xin(1)
    else  ! r,phi,z -> x,y,z
       xout(1) = xin(1)*cos(xin(2))
       xout(2) = xin(1)*sin(xin(2))
       if (ndimout > 2) xout(3) = xin(3)
    endif
!
!--input is spherical polars
!
 case(3)
    !
    ! output is cartesian (default)
    !
    if (itypeout /= 1) ierr = ierr_warning_assuming_cartesian
    select case(ndimout)
    case(1) ! r -> x
       xout(1) = xin(1)
    case(2) ! r,phi -> x,y
       xout(1) = xin(1)*cos(xin(2))
       xout(2) = xin(1)*sin(xin(2))
    case(3) ! r,phi,theta -> x,y,z
       xout(1) = xin(1)*cos(xin(2))*sin(xin(3))
       xout(2) = xin(1)*sin(xin(2))*sin(xin(3))
       xout(3) = xin(1)*cos(xin(3))
    end select
!
!--input is torus co-ordinates
!
 case(4)
    !
    ! output is cartesian (default)
    !
    if (itypeout /= 1) ierr = ierr_warning_assuming_cartesian
    if (ndimin /= 3) then
       xout(1:ndimout) = xin(1:ndimout)
    else
       rcyl = xin(1)*cos(xin(2)) + Rtorus
       xout(1) = rcyl*cos(xin(3))
       if (ndimout >= 2) xout(2) = rcyl*sin(xin(3))
       if (ndimout >= 3) xout(3) = xin(1)*sin(xin(2))
    endif

 case(5)
!
!--input is rotated cartesian coordinates
!
    !
    ! output is cartesian
    !
    if (itypeout /= 1) ierr = ierr_warning_assuming_cartesian
    xi = 0.
    xi(1:ndimin) = xin
    xouti(1) = xi(1)*cosa*cosb - xi(2)*sinb - xi(3)*sina*cosb
    xouti(2) = xi(1)*cosa*sinb + xi(2)*cosb - xi(3)*sina*sinb
    xouti(3) = xi(1)*sina + xi(3)*cosa
    xout(1:ndimout) = xouti
!
!--input is cartesian co-ordinates
!
 case default
    select case(itypeout)
    case(2)
       !
       !--output is cylindrical
       !
       if (ndimin==1) then
          xout(1) = abs(xin(1))   ! cylindrical r
       else
          xout(1) = sqrt(dot_product(xin(1:2),xin(1:2)))
          if (ndimout >= 2) xout(2) = atan2(xin(2),xin(1)) ! phi
          if (ndimout==3) xout(3) = xin(3) ! z
       endif
    case(3)
       !
       !--output is spherical
       !
       xout(1) = sqrt(dot_product(xin,xin))! r
       if (ndimout >= 2) xout(2) = atan2(xin(2),xin(1)) ! phi
       if (ndimout >= 3) then
          ! theta = acos(z/r)
          xout(3) = acos(xin(3)/xout(1))
       endif
    case(4)
       !
       !--output is torus r,theta,phi co-ordinates
       !
       if (ndimin /= 3) then
          ! not applicable if ndim < 3
          xout(1:ndimout) = xin(1:ndimout)
       else
          rcyl = sqrt(xin(1)**2 + xin(2)**2)
          xout(1) = sqrt(xin(3)**2 + (rcyl - Rtorus)**2)
          if (ndimout >= 2) xout(2) = atan2(xin(3),rcyl-Rtorus) ! asin(xin(3)/xout(1))
          if (ndimout >= 3) xout(3) = atan2(xin(2),xin(1))
       endif
    case(5)
       !
       !--output is rotated cartesian x1, x2, x3
       !
       xi = 0.
       xi(1:ndimin) = xin
       xouti(1) = cosa*(xi(1)*cosb + xi(2)*sinb) + xi(3)*sina
       xouti(2) = xi(2)*cosb - xi(1)*sinb
       xouti(3) = -sina*(xi(1)*cosb + xi(2)*sinb) + xi(3)*cosa
       xout(1:ndimout) = xouti(1:ndimout)
    case default
       !
       ! just copy
       !
       xout(1:ndimout) = xin(1:ndimout)
    end select
 end select

 if (present(err)) err = ierr
 return
end subroutine coord_transform

!-----------------------------------------------------------------
! Subroutine to transform vector components
!  between different co-ordinate systems
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! Arguments:
!   xin(ndimin)   : input co-ordinates, in ndimin dimensions
!   vecin(ndimin) : components of vector in input co-ordinate basis
!   itypein      : input co-ordinate type
!
!   vecout(ndimout) : components of vector in output co-ordinate basis
!   itypeout      : output co-ordinate type
!
! coords must be one of the following:
!   'cartesian' (default)
!   'cylindrical'
!   'spherical'
!
! Currently handles:
!
!  cartesian -> cylindrical, spherical polar
!  cylindrical -> cartesian
!  spherical polar -> cartesian
!
!-----------------------------------------------------------------
pure subroutine vector_transform(xin,vecin,ndimin,itypein,vecout,ndimout,itypeout,err)
 integer, intent(in)  :: ndimin,ndimout,itypein,itypeout
 real,    intent(in)  :: xin(ndimin),vecin(ndimin)
 real,    intent(out) :: vecout(ndimout)
 integer, intent(out), optional :: err
 integer :: i,ierr
 real :: dxdx(3,3)
 real :: sinphi, cosphi
 real :: rr,rr1,rcyl,rcyl2,rcyl1

 ierr = 0
!
!--check for errors in input
!
 if (ndimout > ndimin) then
    ierr = ierr_invalid_dims
    if (present(err)) err = ierr
    return
 elseif (itypein==itypeout) then
    vecout(1:ndimout) = vecin(1:ndimout)
    return
 elseif (ndimin < 1.or.ndimin > 3) then
    ierr = ierr_invalid_dimsin
    if (present(err)) err = ierr
    return
 elseif (ndimout < 1.or.ndimout > 3) then
    ierr = ierr_invalid_dimsout
    if (present(err)) err = ierr
    return
 elseif (abs(xin(1)) < 1e-8 .and. &
       (itypein==2 .or. itypein==3)) then
    ierr = ierr_r_is_zero
    if (present(err)) err = ierr
    vecout = 0.
    return
 endif
!
!--set Jacobian matrix to zero
!
 dxdx = 0.
!
!--calculate non-zero components of Jacobian matrix for the transformation
!
 select case(itypein)
!
!--rotated cartesian
!
 case(5)
    call coord_transform(vecin,ndimin,itypein,vecout,ndimout,itypeout,err=ierr)
    if (present(err)) err = ierr
    return
!
!--input is toroidal
!
 case(4)
    select case(itypeout)
    case default
       dxdx(1,1) = cos(xin(2))*cos(xin(3))         ! dx/dr
       dxdx(1,2) = -sin(xin(2))*cos(xin(3))        ! 1/r dx/dtheta
       dxdx(1,3) = sin(xin(3))                     ! 1/rcyl dx/dphi
       dxdx(2,1) = cos(xin(2))*sin(xin(3))         ! dy/dr
       dxdx(2,2) = -sin(xin(2))*sin(xin(3))        ! 1/r dy/dtheta
       dxdx(2,3) = cos(xin(3))                     ! 1/rcyl dy/dphi
       dxdx(3,1) = sin(xin(3))                     ! dz/dr
       dxdx(3,2) = cos(xin(3))                     ! 1/r dz/dtheta
!        dxdx(3,3) = 0.                             ! dz/dphi
    end select
!
!--input is spherical polars
!
 case(3)
    select case(itypeout)
    case default
       !
       ! output is cartesian (default)
       !
       dxdx(1,1) = cos(xin(2))*sin(xin(3))         ! dx/dr
       dxdx(1,2) = -sin(xin(2))                    ! 1/rcyl dx/dphi
       dxdx(1,3) = cos(xin(2))*cos(xin(3))         ! 1/r dx/dtheta
       dxdx(2,1) = sin(xin(2))*sin(xin(3))         ! dy/dr
       dxdx(2,2) = cos(xin(2))                     ! 1/rcyl dy/dphi
       dxdx(2,3) = sin(xin(2))*cos(xin(3))         ! 1/r dy/dtheta
       dxdx(3,1) = cos(xin(3))                     ! dz/dr
       dxdx(3,3) = -sin(xin(3))                    ! 1/r dz/dtheta
    end select
!
!--input is cylindrical polars
!
 case(2)
    select case(itypeout)
    case default
       !
       ! output is cartesian (default)
       !
       sinphi = sin(xin(2))
       cosphi = cos(xin(2))
       dxdx(1,1) = cosphi            ! dx/dr
       dxdx(1,2) = -sinphi           ! 1/r*dx/dphi
       dxdx(2,1) = sinphi            ! dy/dr
       dxdx(2,2) = cosphi            ! 1/r*dy/dphi
       dxdx(3,3) = 1.                ! dz/dz
    end select
!
!--input is cartesian co-ordinates (default)
!
 case default
    select case(itypeout)
    case(5)
       call coord_transform(vecin,ndimin,itypein,vecout,ndimout,itypeout,err=ierr)
       return
    case(4)
       !
       ! output is toroidal
       !
       rcyl = sqrt(xin(1)**2 + xin(2)**2)
       if (rcyl > tiny(rcyl)) then
          rcyl1 = 1./rcyl
       else
          rcyl1 = 0.
       endif
       rr = sqrt((rcyl - Rtorus)**2 + xin(3)**2)
       if (rr > tiny(rr)) then
          rr1 = 1./rr
       else
          rr1 = 0.
       endif
       dxdx(1,1) = (rcyl - Rtorus)*xin(1)*rr1*rcyl1 ! dr/dx
       dxdx(1,2) = (rcyl - Rtorus)*xin(2)*rr1*rcyl1 ! dr/dy
       dxdx(1,3) = xin(3)*rr1                       ! dr/dz
       dxdx(2,1) = -xin(3)*xin(1)*rr1*rcyl1         ! dtheta/dx
       dxdx(2,2) = -xin(3)*xin(2)*rr1*rcyl1         ! dtheta/dy
       dxdx(2,3) = (rcyl - Rtorus)*rr1              ! dtheta/dz
       dxdx(3,1) = -xin(2)*rcyl1                    ! dphi/dx
       dxdx(3,2) = xin(1)*rcyl1                     ! dphi/dy
!        dxdx(3,3) = 0.                              ! dphi/dz
    case(3)
       !
       ! output is spherical
       !
       rr = sqrt(dot_product(xin,xin))
       if (rr > tiny(rr)) then
          rr1 = 1./rr
       else
          rr1 = 0.
       endif
       dxdx(1,1) = xin(1)*rr1  ! dr/dx
       if (ndimin >= 2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
       if (ndimin==3) dxdx(1,3) = xin(3)*rr1  ! dr/dz
       if (ndimin >= 2) then
          rcyl2 = dot_product(xin(1:2),xin(1:2))
          rcyl = sqrt(rcyl2)
          if (rcyl > tiny(rcyl)) then
             rcyl1 = 1./rcyl
          else
             rcyl1 = 0.
          endif
          dxdx(2,1) = -xin(2)*rcyl1 ! rcyl dphi/dx
          dxdx(2,2) = xin(1)*rcyl1  ! rcyl dphi/dy
          dxdx(2,3) = 0.
          if (ndimin >= 3) then
             dxdx(3,1) = xin(1)*xin(3)*rr1*rcyl1 ! r dtheta/dx
             dxdx(3,2) = xin(2)*xin(3)*rr1*rcyl1 ! r dtheta/dy
             dxdx(3,3) = -rcyl2*rr1*rcyl1 ! r dtheta/dz
          endif
       endif
    case(2)
       !
       !--output is cylindrical
       !
       rr = sqrt(dot_product(xin(1:min(ndimin,2)),xin(1:min(ndimin,2))))
       if (rr > tiny(rr)) then
          rr1 = 1./rr
       else
          rr1 = 0.
       endif
       dxdx(1,1) = xin(1)*rr1  ! dr/dx
       if (ndimin >= 2) dxdx(1,2) = xin(2)*rr1  ! dr/dy
       if (ndimout >= 2) then
          dxdx(2,1) = -xin(2)*rr1 ! r*dphi/dx
          dxdx(2,2) = xin(1)*rr1  ! r*dphi/dy
          if (ndimout==3) dxdx(3,3) = 1.  ! dz/dz
       endif
    case default
       ierr = ierr_warning_assuming_cartesian
       vecout(1:ndimout) = vecin(1:ndimout)
       return
    end select
 end select
!
!--now perform transformation using Jacobian matrix
!
 do i=1,ndimout
    vecout(i) = dot_product(dxdx(i,1:ndimin),vecin(1:ndimin))
 enddo

 if (present(err)) err = ierr
 return
end subroutine vector_transform

!------------------------------------------------------------------
! this subroutine attempts to switch plot limits / boundaries
! between various co-ordinate systems.
!------------------------------------------------------------------
subroutine coord_transform_limits(xmin,xmax,itypein,itypeout,ndim)
 integer, intent(in) :: itypein,itypeout,ndim
 real, dimension(ndim), intent(inout) :: xmin,xmax
 real, dimension(ndim) :: xmaxtemp,xmintemp
!
!--check for errors in input
!
 if (ndim < 1 .or. ndim > 3) then
    print*,'Error: limits coord transform: ndim invalid on input'
    return
 endif
 !print*,'modifying plot limits for new coordinate system'
!
!--by default do nothing
!
 xmintemp(1:ndim) = xmin(1:ndim)
 xmaxtemp(1:ndim) = xmax(1:ndim)

 select case(itypein)
 case(5)
!
!--rotated cartesian
!
    call coord_transform(xmin,ndim,itypein,xmintemp,ndim,itypeout)
    call coord_transform(xmax,ndim,itypein,xmaxtemp,ndim,itypeout)
!
!--input is toroidal
!
 case(4)
    select case(itypeout)
    case default
       !
       !--cartesian output
       !
       xmintemp(1:min(ndim,2)) = -Rtorus - xmax(1)
       xmaxtemp(1:min(ndim,2)) = Rtorus + xmax(1)
       if (ndim==3) then
          xmintemp(3) = -xmax(1)
          xmaxtemp(3) = xmax(1)
       endif

    end select
!
!--input is spherical
!
 case(3)
    select case(itypeout)
    case default
       !
       !--cartesian output
       !
       xmintemp(1:ndim) = -xmax(1)
       xmaxtemp(1:ndim) = xmax(1)

    end select
!
!--input is cylindrical
!
 case(2)
    select case(itypeout)
    case default
       !
       !--cartesian output
       !
       xmintemp(1:max(ndim,2)) = -xmax(1)
       xmaxtemp(1:max(ndim,2)) = xmax(1)

    end select
!
!--input is cartesian
!
 case default
    select case(itypeout)
    case(5)
       !
       !--rotated cartesian
       !
       call coord_transform(xmin,ndim,itypein,xmintemp,ndim,itypeout)
       call coord_transform(xmax,ndim,itypein,xmaxtemp,ndim,itypeout)

    case(4)
       !
       !--output is toroidal
       !
       xmintemp(1) = 0.
       xmaxtemp(1) = max(maxval(abs(xmax(1:min(ndim,2))))-Rtorus, &
                      maxval(abs(xmin(1:min(ndim,2))))-Rtorus)
       if (ndim >= 2) then
          xmintemp(2) = -0.5*pi
          xmaxtemp(2) = 0.5*pi
          if (ndim >= 3) then
             xmintemp(3) = -pi
             xmaxtemp(3) = pi
          endif
       endif
       !
       !--output is spherical
       !
    case(3)
       !--rmin, rmax
       xmintemp(1) = 0.
       xmaxtemp(1) = max(maxval(abs(xmin(1:ndim))), &
                         maxval(abs(xmax(1:ndim))))
       if (ndim >= 2) then
          xmintemp(2) = -pi
          xmaxtemp(2) = pi
          if (ndim >= 3) then
             xmintemp(3) = 0.
             xmaxtemp(3) = pi
          endif
       endif
       !
       !--output is cylindrical
       !
    case(2)
       !--rmin, rmax
       xmintemp(1) = 0.
       if (ndim >= 2) then
          xmaxtemp(1) = sqrt(max((xmin(1)**2 + xmin(2)**2), &
                                 (xmax(1)**2 + xmax(2)**2)))
          xmintemp(2) = -pi
          xmaxtemp(2) = pi
       else
          xmaxtemp(1) = max(abs(xmin(1)),abs(xmax(1)))
       endif
    end select
 end select

 xmin(:) = min(xmintemp(:),xmaxtemp(:))
 xmax(:) = max(xmintemp(:),xmaxtemp(:))

 return
end subroutine coord_transform_limits

end module geometry
