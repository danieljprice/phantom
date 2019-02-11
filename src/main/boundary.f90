!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: boundary
!
!  DESCRIPTION:
!  This module contains variables and subroutines
!  relating to boundaries
!
!  REFERENCES: None
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

module boundary
 implicit none
 real, public :: xmin,xmax,ymin,ymax,zmin,zmax
 real, public :: dxbound, dybound, dzbound
 real, public :: totvol

 public :: set_boundary
 public :: cross_boundary

 private

contains

!---------------------------------------------------------------
!+
!  Routine to set the (cartesian) boundaries for the code
!  when using periodic boundary conditions
!
!  This sets the module variables xmin,xmax,ymin,ymax,zmin,zmax
!  as well as the subsidiary variables dxbound,dybound,dzbound
!  and totvol.
!
!  Can be called with no arguments, which gives the defaults:
!
!   call set_boundary()
!
!  Or with each boundary set individually:
!
!   call set_boundary(0.,1.,0.,1.,0.,1.)
!
!  Or with an array of values array=(/xmin,xmax,ymin,ymax,zmin,zmax/)
!
!   call set_boundary(pos=array)
!+
!---------------------------------------------------------------
subroutine set_boundary(x_min,x_max,y_min,y_max,z_min,z_max,pos)
 real, intent(in), optional :: x_min, x_max, y_min, y_max, z_min, z_max
 real, intent(in), optional :: pos(6)
 !
 ! give default values if no settings given
 !
 xmin = -0.5
 xmax = 0.5
 ymin = -0.5
 ymax = 0.5
 zmin = -0.5
 zmax = 0.5
 !
 ! set the min and max position in each coordinate
 !
 if (present(pos)) then
    xmin = pos(1)
    xmax = pos(2)
    ymin = pos(3)
    ymax = pos(4)
    zmin = pos(5)
    zmax = pos(6)
 else
    if (present(x_min)) xmin = x_min
    if (present(x_max)) xmax = x_max
    if (present(y_min)) ymin = y_min
    if (present(y_max)) ymax = y_max
    if (present(z_min)) zmin = z_min
    if (present(z_max)) zmax = z_max
 endif
 !
 ! set subsidiary quantities that depend on above settings
 ! these are to save computation in periodicity calculations
 !
 dxbound = xmax - xmin
 dybound = ymax - ymin
 dzbound = zmax - zmin

 totvol = dxbound*dybound*dzbound

end subroutine set_boundary

!----------------------------------------------------------------
!+
!  This subroutine determines whether particles should cross
!  the boundary
!+
!---------------------------------------------------------------
subroutine cross_boundary(isperiodic,xyz,ncross)
 logical, intent(in)    :: isperiodic(3)
 real,    intent(inout) :: xyz(:)
 integer, intent(inout) :: ncross
!
!--check for crossings of periodic boundaries
!
 if (isperiodic(1)) then
    if (xyz(1) < xmin) then
       xyz(1) = xyz(1) + dxbound
       ncross = ncross + 1
    elseif (xyz(1) > xmax) then
       xyz(1) = xyz(1) - dxbound
       ncross = ncross + 1
    endif
 endif

 if (isperiodic(2)) then
    if (xyz(2) < ymin) then
       xyz(2) = xyz(2) + dybound
       ncross = ncross + 1
    elseif (xyz(2) > ymax) then
       xyz(2) = xyz(2) - dybound
       ncross = ncross + 1
    endif
 endif

 if (isperiodic(3)) then
    if (xyz(3) < zmin) then
       xyz(3) = xyz(3) + dzbound
       ncross = ncross + 1
    elseif (xyz(3) > zmax) then
       xyz(3) = xyz(3) - dzbound
       ncross = ncross + 1
    endif
 endif

 return
end subroutine cross_boundary

end module boundary
