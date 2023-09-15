!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr_region
  !
  ! Contains everything for live adaptive particle refinement
  !
  ! :References: None
  !
  ! :Owner: Rebecca Nealon
  !
  ! :Runtime parameters:
  !   - apr_max_in        : number of refinement levels (3 -> 2x resolution)
  !   - ref_dir           : increase (1) or decrease (-1) resolution from the base resolution
  !   - [x,y,z]_centre    : centre coordinates of the region to be more highly resolved
  !   - apr_rad           : radius of the region to be more highly resolved
  !
  ! :Dependencies: None
  !
  implicit none

  logical, public :: dynamic_apr = .false.
  public :: set_apr_region

  private

contains

  !-----------------------------------------------------------------------
  !+
  !  Initialising all the apr arrays and properties
  !+
  !-----------------------------------------------------------------------

subroutine set_apr_region(apr_type,apr_centre,apr_rad)
  integer, intent(in) :: apr_type
  real,    intent(out) :: apr_centre(3),apr_rad

  select case (apr_type)

  case(1) ! a static circle
    dynamic_apr = .false.
    apr_centre(1) = 0.0
    apr_centre(2) = 0.0
    apr_centre(3) = 0.0
    apr_rad = 0.2

  case default
    dynamic_apr = .false.
    apr_centre(:) = 0.
    apr_rad = 0.2

  end select

end subroutine set_apr_region

end module apr_region
