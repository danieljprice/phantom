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

subroutine set_apr_region(apr_type,apr_centre,apr_rad,apr_blend)
  use part, only: xyzmh_ptmass
  integer, intent(in) :: apr_type
  real,    intent(out) :: apr_centre(3),apr_rad,apr_blend

  select case (apr_type)

  case(1) ! a static circle
    dynamic_apr = .false.
    apr_centre(1) = 0.0
    apr_centre(2) = 0.0
    apr_centre(3) = 0.0
    apr_rad = 0.2
    apr_blend = 0.1

  case(2) ! around sink particle 2 - e.g. a planet
    dynamic_apr = .true.
    apr_centre(1) = xyzmh_ptmass(1,2)
    apr_centre(2) = xyzmh_ptmass(2,2)
    apr_centre(3) = xyzmh_ptmass(3,2)
    apr_rad = 2.0
    apr_blend = 0.1

  case default
    dynamic_apr = .false.
    apr_centre(:) = 0.
    apr_rad = 0.2
    apr_blend = 0.1

  end select

end subroutine set_apr_region

end module apr_region
