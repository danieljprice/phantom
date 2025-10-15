!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module metric_interp
!
! Interpolate a tabulated metric onto the particle positions
!
! :References:
!  Magnall, Price, Lasky & Macpherson (2023), Phys. Rev D. 108, 103534
!
! :Owner: Spencer Magnall
!
! :Runtime parameters: None
!
! :Dependencies: einsteintk_utils
!
 implicit none

 interface trilinear_interp
  module procedure interp_g, interp_sqrtg, interp_gderiv
 end interface trilinear_interp

contains

subroutine interp_g()

end subroutine interp_g

subroutine interp_sqrtg()

end subroutine interp_sqrtg

subroutine interp_gderiv()

end subroutine interp_gderiv

pure subroutine get_grid_neighbours(position,dx,xlower,ylower,zlower)
 use einsteintk_utils, only:gridorigin
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

end module metric_interp
