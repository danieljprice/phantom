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
! Current options include:
!
! 1: A static sphere; absolute position set in *.in file
!
! 2: A sphere around a sink particle; sink is identified as track_part in *.in
!
! 3: A sphere around fragments; designed for use in an SG disc
!
! 4: A sphere centred on two sinks; this is centred between two sequential sinks
!    (i.e. modelling the region around a binary star)
!
! 5: A sphere centered on the CoM of the system
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, io, part, ptmass, units, utils_apr
!

 use utils_apr

 implicit none

contains

!-----------------------------------------------------------------------
!+
!  Setting/updating the centre of the apr region (as it may move)
!+
!-----------------------------------------------------------------------
subroutine set_apr_centre(apr_type,apr_centre,ntrack,track_part)
 use part, only: xyzmh_ptmass,xyzh,npart,vxyzu,nptmass,vxyz_ptmass,apr_level
 use part, only: aprmassoftype, poten
 use centreofmass, only:get_centreofmass
 integer, intent(in)  :: apr_type
 real,    intent(out) :: apr_centre(3,ntrack_max)
 integer, optional, intent(in) :: ntrack,track_part(:)
 real :: xcom(3), vcom(3)
 integer :: ii, ntrack_temp, track_part_temp(ntrack_max)
 integer, save :: count = 0

 select case (apr_type)

 case(1) ! a static circle
    ! do nothing here

 case(2) ! around sink particle named track_part
    apr_centre(1,1) = xyzmh_ptmass(1,track_part(1))
    apr_centre(2,1) = xyzmh_ptmass(2,track_part(1))
    apr_centre(3,1) = xyzmh_ptmass(3,track_part(1))

 case(3) ! to derefine a clump

    ! to speed things up, we only call this every 10 steps
    count = count + 1
    if (count == 10) then
       ! find potential clumps
       call identify_clumps(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,aprmassoftype, &
                              ntrack_temp,track_part_temp)

       ! update or create from existing clumps
       call create_or_update_apr_clump(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,&
                                      aprmassoftype,ntrack_temp,track_part_temp)

       ! now update the clump locations
       do ii = 1,ntrack
          apr_centre(1:3,ii) = xyzh(1:3,track_part(ii))
       enddo

       ! for next time
       count = 0
    endif

 case(4) ! averaging two sequential sinks
    apr_centre(1,1) = 0.5*(xyzmh_ptmass(1,track_part(1)) + xyzmh_ptmass(1,track_part(1) + 1))
    apr_centre(2,1) = 0.5*(xyzmh_ptmass(2,track_part(1)) + xyzmh_ptmass(2,track_part(1) + 1))
    apr_centre(3,1) = 0.5*(xyzmh_ptmass(3,track_part(1)) + xyzmh_ptmass(3,track_part(1) + 1))

 case(5) ! centre of mass
    call get_centreofmass(xcom,vcom,npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
    apr_centre(:,1) = xcom

 case(6) ! vertically uniformly resolved disc
    !call run_apr_disc_analysis(100,xyzh,vxyzu,apr_H)

 case default ! used for the test suite
    ! do nothing

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
       apr_regions(kk) = apr_rad + (ii-2)*apr_drad
    enddo
 else
    apr_regions(apr_max) = huge(apr_regions(apr_max)) ! again this just needs to encompass the whole domain
    do ii = 1,apr_max-1
       apr_regions(ii) = apr_rad + (ii-1)*apr_drad
    enddo
 endif

end subroutine set_apr_regions

!-----------------------------------------------------------------------
!+
!  Analysis to see where the s-g clumps are - this function identifies
!  *potential* track_part, these are refined and created in
!  create_or_update_clumps
!+
!-----------------------------------------------------------------------
subroutine identify_clumps(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,aprmassoftype, &
                           ntrack_temp,track_part_temp)
 use utils_apr, only: find_inner_and_outer_radius
 use part,      only:igas,rhoh,isdead_or_accreted
 use ptmass,    only:rho_crit_cgs
 use units,     only:unit_density
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: apr_level(:)
 real, intent(in) :: xyzh(:,:), vxyzu(:,:), aprmassoftype(:,:),xyzmh_ptmass(:,:)
 real(kind=4), intent(in) :: poten(:)
 integer, intent(out) :: ntrack_temp, track_part_temp(:)
 integer :: ii, kk
 real :: pmassi, rhoi, r2test, rminlimit2, rin, rout

 ! find the inner and outer radius each time
 call find_inner_and_outer_radius(npart,xyzh,rin,rout)

 ! initialise
 ntrack_temp = 0
 track_part_temp(:) = 0

 ! we shouldn't have clumps too close to the inner edge of the disc, so set a limit
 rminlimit2 = (1.2*rin)**2

 !iterate over all particles and find the ones that are above a certain density threshold

 !$omp parallel do schedule(guided) default(none) &
 !$omp shared(npart,xyzh,apr_level,aprmassoftype,unit_density) &
 !$omp shared(rho_crit_cgs,ntrack_temp,track_part_temp) &
 !$omp shared(rminlimit2,ntrack,track_part) &
 !$omp private(pmassi,rhoi,ii,kk,r2test)
 over_dens: do ii = 1, npart

    ! check the particle isn't dead or accreted
    if (isdead_or_accreted(xyzh(4,ii))) cycle over_dens

    ! Obtain particle density
    pmassi = aprmassoftype(igas,apr_level(ii))
    rhoi = rhoh(xyzh(4,ii),pmassi)

    ! does it have the density we need?
    if (rhoi*unit_density < rho_crit_cgs) cycle over_dens

    ! check we haven't seen this particle already and that we're not already tracking it
    do kk = 1,ntrack_temp
       if (track_part_temp(kk) == ii) cycle over_dens
    enddo
    do kk = 1,ntrack
       if (track_part(kk) == ii) cycle over_dens
    enddo

    ! check that it's not too close to the inner edge of the disc
    r2test = dot_product(xyzh(1:3,ii),xyzh(1:3,ii))
    if (r2test < rminlimit2) cycle over_dens

    ! if we've met the above criteria, add it to the list
    ntrack_temp = ntrack_temp + 1
    track_part_temp(ntrack_temp) = ii

 enddo over_dens
 !$omp end parallel do

end subroutine identify_clumps

 !-----------------------------------------------------------------------
!+
!  Given a list track_part of ntrack particles to track, create or update
!  clumps as appropriate
!+
!-----------------------------------------------------------------------
subroutine create_or_update_apr_clump(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,&
                                      aprmassoftype,ntrack_temp,track_part_temp)
 use utils_apr, only:find_closest_region
 use part,      only:igas,rhoh
 use io,        only:fatal
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: apr_level(:)
 real, intent(in) :: xyzh(:,:), vxyzu(:,:), aprmassoftype(:,:),xyzmh_ptmass(:,:)
 real(kind=4), intent(in) :: poten(:)
 integer, intent(in) :: ntrack_temp, track_part_temp(:)
 integer :: ii, ll, jj, kk
 real :: pmassi, rhoitest, rhoiexisting, rtest, xi(3)

 over_mins: do jj = 1,ntrack_temp
    ii = track_part_temp(jj) ! this is the potential particle to track

    ! density at this location
    pmassi = aprmassoftype(igas,apr_level(ii))
    rhoitest = rhoh(xyzh(4,ii),pmassi)

    ! check if its inside an existing region
    call find_closest_region(xyzh(1:3,ii),ll)
    if (ll > 0) then
       xi = xyzh(1:3,ii) - apr_centre(1:3,ll)
       kk = track_part(ll) ! this is the particle at the centre of the closest region
       rtest = sqrt(dot_product(xi(1:3),xi(1:3)))
    else
       rtest = 0. ! to prevent warnings, but if this occurs then we're not doing the first option next
    endif

    if ((rtest < apr_rad) .and. (ll > 0)) then ! it's already part of the region

       ! is it's density higher than the existing centre?
       pmassi = aprmassoftype(igas,apr_level(kk))
       rhoiexisting = rhoh(xyzh(4,kk),pmassi)

       if (rhoitest > rhoiexisting) then
          ! we replace the existing centre but keep ntrack the same
          !print*, 'UPDATING CLUMP:', ll, ' from', apr_centre(1:2,ll), 'to ', xyzh(1:2,ii)
          track_part(ll) = ii
          apr_centre(1:3,ll) = xyzh(1:3,ii)
       endif

       ! if it's not, we can safely discard it anyway

    else
       ! this is a new region
       ntrack = ntrack + 1
       track_part(ntrack) = ii
       apr_centre(1:3,ntrack) = xyzh(1:3,ii)
       !print*, 'NEW CLUMP:', ntrack, 'with density ', rhoitest*unit_density, 'assigning clump centre at', xyzh(1:2,ii)

    endif

    if (ntrack > ntrack_max) call fatal('create_or_update_clumps',&
        'too many clumps found, increase ntrack_max',var='ntrack',ival=ntrack)

 enddo over_mins

end subroutine create_or_update_apr_clump

end module apr_region
