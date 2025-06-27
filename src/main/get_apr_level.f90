!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module get_apr_level
!
! Module that holds the routines to return get_apr. This is where you set
! the shape of your APR region.
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: apr_region, dim, io, utils_apr
!
 use dim, only:use_apr
 use apr_region
 use utils_apr

 implicit none

 public :: get_apr
 public :: create_or_update_apr_clump

 procedure(get_apr_sphere), pointer :: get_apr => get_apr_sphere

contains

!-----------------------------------------------------------------------
!+
!  routine to set the get_apr function correctly
!+
!-----------------------------------------------------------------------
subroutine set_get_apr()

 ! For the apr type, chose the region shape
 if (apr_type == 6) then
    ref_dir = -1 ! need to enforce this for this one
 else
    get_apr => get_apr_sphere
 endif

 ! Here set the requirements for the apr_type to read in the right values
! for apr_types that read in a particle number
 if (apr_type == 2) then
    track_part(1) = track_part_in
 endif

 ! for apr_types that read in the centre values from the *.in file
 if (apr_type == 1) apr_centre(1:3,1) = apr_centre_in(1:3) ! from the .in file

end subroutine set_get_apr


!-----------------------------------------------------------------------
!+
!  routine to return the adaptive particle refinement level based on position
!  and the boundaries set by the apr_* arrays for a spherical region
!+
!-----------------------------------------------------------------------
subroutine get_apr_sphere(pos,icentre,apri)
 use io, only:fatal
 use apr_region, only:apr_region_is_circle
 real, intent(in)     :: pos(3)
 integer, intent(in)  :: icentre
 integer, intent(out) :: apri
 integer :: jj, kk
 real :: dx,dy,dz,r

 apri = -1 ! to prevent compiler warnings

 do jj = 1,apr_max
    if (ref_dir == 1) then
       kk = apr_max - jj + 1       ! going from apr_max -> 1
    else
       kk = jj                    ! going from 1 -> apr_max
    endif
    dx = pos(1) - apr_centre(1,icentre)
    dy = pos(2) - apr_centre(2,icentre)
    dz = pos(3) - apr_centre(3,icentre)

    if (apr_region_is_circle) dz = 0.

    r = sqrt(dx**2 + dy**2 + dz**2)
    if (r < apr_regions(kk)) then
       apri = kk
       return
    endif
 enddo

 if (apri == -1) call fatal('apr_region, get_apr','could not find apr level')

end subroutine get_apr_sphere

end module get_apr_level
