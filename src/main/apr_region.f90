!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module apr_region
!
! Everything for setting the adaptive particle refinement regions
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: part
!
 implicit none

 logical, public :: dynamic_apr = .false., apr_region_is_circle = .false.
 public :: set_apr_centre, set_apr_regions

 private

contains

!-----------------------------------------------------------------------
!+
!  Setting/updating the centre of the apr region (as it may move)
!+
!-----------------------------------------------------------------------
subroutine set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 use part, only: xyzmh_ptmass,xyzh
 integer, intent(in)  :: apr_type
 real,    intent(out) :: apr_centre(3)
 integer, optional, intent(in) :: ntrack,track_part

 select case (apr_type)

 case(1) ! a static circle
    ! do nothing here

 case(2) ! around sink particle named track_part
    dynamic_apr = .true.
    apr_centre(1) = xyzmh_ptmass(1,track_part)
    apr_centre(2) = xyzmh_ptmass(2,track_part)
    apr_centre(3) = xyzmh_ptmass(3,track_part)

 case(3) ! to derefine a clump - only activated when the centre of the clump
    ! has been found
    dynamic_apr = .true.
    if (present(ntrack)) then
       apr_centre(1) = xyzh(1,track_part)
       apr_centre(2) = xyzh(2,track_part)
       apr_centre(3) = xyzh(3,track_part)
    else
       apr_centre = tiny(apr_centre) ! this *might* be safe? Just want it to be irrelevant
    endif

 case default ! used for the test suite
    apr_centre(:) = 0.

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
    apr_regions(1) = huge(apr_regions(1)) ! this needs to be a number that encompasses the whole domain
    do ii = 2,apr_max
       kk = apr_max - ii + 2
       apr_regions(kk) = apr_rad + (ii-1)*apr_drad
    enddo
 else
    apr_regions(apr_max) = huge(apr_regions(apr_max)) ! again this just needs to encompass the whole domain
    do ii = 1,apr_max-1
       apr_regions(ii) = apr_rad + (ii-1)*apr_drad
    enddo
 endif

end subroutine set_apr_regions


end module apr_region
