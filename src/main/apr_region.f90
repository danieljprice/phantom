!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr_region
  !
  ! Contains everything for setting the adaptive particle refinement regions
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
  public :: set_apr_centre, set_apr_regions

  private

contains

  !-----------------------------------------------------------------------
  !+
  !  Setting/updating the centre of the apr region (as it may move)
  !+
  !-----------------------------------------------------------------------

subroutine set_apr_centre(apr_type,apr_centre,apr_blend)
  use part, only: xyzmh_ptmass
  integer, intent(in)  :: apr_type
  real,    intent(out) :: apr_centre(3),apr_blend

  select case (apr_type)

  case(1) ! a static circle
    dynamic_apr = .false.
    apr_centre(1) = 0.0
    apr_centre(2) = 0.0
    apr_centre(3) = 0.0
    apr_blend = 0.1

  case(2) ! around sink particle 2 - e.g. a planet
    dynamic_apr = .true.
    apr_centre(1) = xyzmh_ptmass(1,2)
    apr_centre(2) = xyzmh_ptmass(2,2)
    apr_centre(3) = xyzmh_ptmass(3,2)
    apr_blend = 0.1

  case default
    dynamic_apr = .false.
    apr_centre(:) = 0.
    apr_blend = 0.1

  end select

end subroutine set_apr_centre

!-----------------------------------------------------------------------
!+
!  Initialising all the apr region arrays that decide
!  the spatial arrangement of the regions
!+
!-----------------------------------------------------------------------

subroutine set_apr_regions(ref_dir,apr_max,apr_regions,apr_rad,apr_drad)
  integer, intent(in) :: ref_dir,apr_max
  real, intent(in)    :: apr_rad,apr_drad
  real, intent(inout) :: apr_regions(apr_max)
  integer :: ii,kk

  if (ref_dir == 1) then
    apr_regions(1) = 1000. ! this needs to be a number that encompasses the whole domain
    do ii = 2,apr_max
      kk = apr_max - ii + 2
      apr_regions(kk) = apr_rad + (ii-1)*apr_drad
    enddo
  else
    apr_regions(apr_max) = 1000. ! again this just needs to encompass the whole domain
    do ii = 1,apr_max-1
      apr_regions(ii) = apr_rad + (ii-1)*apr_drad
    enddo
  endif

end subroutine set_apr_regions

end module apr_region
