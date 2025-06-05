!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module get_apr_level
!
! Module that holds the routines to return get_apr
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: 
!
!
!
!
 use dim, only:use_apr
 use apr_region
 implicit none

 public :: get_apr, get_apr_multiple
 public :: create_or_update_apr_clump

 procedure(get_apr_sphere), pointer :: get_apr => get_apr_sphere

 contains

!-----------------------------------------------------------------------
!+
!  routine to set the get_apr function correctly
!+
!-----------------------------------------------------------------------
subroutine set_get_apr()

 if (apr_type == 6) then
    print*,'not using that now'
    ref_dir = -1 ! need to enforce this for this one
 elseif (apr_type == 3) then
    get_apr => get_apr_multiple
 else
    get_apr => get_apr_sphere
 endif

end subroutine set_get_apr


!-----------------------------------------------------------------------
!+
!  routine to return the adaptive particle refinement level based on position
!  and the boundaries set by the apr_* arrays for a spherical region
!+
!-----------------------------------------------------------------------
subroutine get_apr_sphere(pos,apri)
 use io, only:fatal
 use apr_region, only:apr_region_is_circle
 real, intent(in)     :: pos(3)
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
    dx = pos(1) - apr_centre(1,1)
    dy = pos(2) - apr_centre(2,1)
    dz = pos(3) - apr_centre(3,1)

    if (apr_region_is_circle) dz = 0.
    
    r = sqrt(dx**2 + dy**2 + dz**2)
    
    if (r < apr_regions(kk)) then
       apri = kk
       return
    endif
 enddo

 if (apri == -1) call fatal('apr_region, get_apr','could not find apr level')

end subroutine get_apr_sphere

!-----------------------------------------------------------------------
!+
!  same as get_apr_sphere but for multiple regions
!+
!-----------------------------------------------------------------------
subroutine get_apr_multiple(pos,apri)
 use io, only:fatal
 use apr_region, only:apr_region_is_circle, ntrack
 real, intent(in)     :: pos(3)
 integer, intent(out) :: apri
 integer :: jj, kk, ii, apri_test(ntrack), extra_apri_test
 real :: dx,dy,dz,r

 apri = -1 ! to prevent compiler warnings
 
 if (ntrack == 0) then
   if (ref_dir == 1) apri = 1 ! base level
   if (ref_dir == -1) apri = apr_max
   return
 endif

 over_napr: do ii = 1,ntrack
    do jj = 1,apr_max
     if (ref_dir == 1) then
        kk = apr_max - jj + 1       ! going from apr_max -> 1
     else
        kk = jj                    ! going from 1 -> apr_max
     endif
     dx = pos(1) - apr_centre(1,ii)
     dy = pos(2) - apr_centre(2,ii)
     dz = pos(3) - apr_centre(3,ii)

     if (apr_region_is_circle) dz = 0.
    
     r = sqrt(dx**2 + dy**2 + dz**2)
    
     if (r < apr_regions(kk)) then
        apri_test(ii) = kk
        cycle over_napr
     endif
    enddo
 enddo over_napr
 
 if (ref_dir == 1) then
   extra_apri_test = -1
   do ii = 1,ntrack
      if (apri_test(ii) > extra_apri_test) extra_apri_test = apri_test(ii) ! take the biggest
   enddo
 else
   extra_apri_test = 100
   do ii = 1,ntrack
      if (apri_test(ii) < extra_apri_test) extra_apri_test = apri_test(ii) ! take the lowest
   enddo   
 endif

 apri = extra_apri_test

   !print*,ntrack
   !print*,r,apri_test
   !print*,extra_apri_test
   !read*

 if (apri == -1) call fatal('apr_region, get_apr','could not find apr level')

end subroutine get_apr_multiple

 !-----------------------------------------------------------------------
!+
!  Given a list track_part of ntrack particles to track, create or update
!  clumps as appropriate - this function has to be split with 
!  identify_clumps because it uses get_apr
!+
!-----------------------------------------------------------------------
subroutine create_or_update_apr_clump(npart,xyzh,vxyzu,poten,apr_level,xyzmh_ptmass,&
                                      aprmassoftype,ntrack_temp,track_part_temp)
 use part, only:igas,rhoh
 use apr_region, only : ntrack, track_part, apr_rad
 use io, only: fatal
 integer, intent(in) :: npart
 integer(kind=1), intent(in) :: apr_level(:)
 real, intent(in) :: xyzh(:,:), vxyzu(:,:), aprmassoftype(:,:),xyzmh_ptmass(:,:)
 real(kind=4), intent(in) :: poten(:)
 integer, intent(in) :: ntrack_temp, track_part_temp(:)
 integer :: ii, jj, ll, mm, ll_found
 real :: r2test, r2, xi, yi, zi, pmassmm, pmassii


!print*,'ntrack found',ntrack_temp

 over_mins: do jj = 1,ntrack_temp
   ii = track_part_temp(jj)

   ! get the refinement level of the potential particle to track
   !call get_apr(xyzh(1:3,ii),apri)

   ! check if this particle is inside an existing region
   !if (apri > 1) then

      ! work out which region it currently belongs to
      r2 = huge(r2)
      ll_found = -1
      do ll = 1, ntrack ! these are the existing regions
         xi = xyzh(1,ii) - apr_centre(1,ll)
         yi = xyzh(2,ii) - apr_centre(2,ll)
         zi = xyzh(3,ii) - apr_centre(3,ll)
         r2test = xi**2 + yi**2 + zi**2
         if (r2test < r2 .and. (sqrt(r2test) < apr_rad)) then
            ll_found = ll
            r2 = r2test
         endif
      enddo

      !if (ll_found < 0) then
      !   print*,'stop, cannot find the region this particle belongs to'
      !endif
      
   if (ll_found > 0) then
      ! the particle at the centre of that is
      mm = track_part(ll_found)

      ! For that clump, establish if it has a lower potential energy than
      ! the current centre of the clump - careful to take into account mass
      pmassmm = aprmassoftype(igas,apr_level(mm))
      pmassii = aprmassoftype(igas,apr_level(ii))
      ! If it does, reset the particle for that clump to be centred on
      ! this current particle
      !if ((poten(mm)/pmassmm) < (poten(ii)/pmassii)) then
      !   track_part(ll_found) = ii
      !endif
   else

      ! create a new particle to track
      ntrack = ntrack + 1
      track_part(ntrack) = ii
   endif
 enddo over_mins

 if (ntrack > ntrack_max) then
   print*,ntrack,ntrack_max
   call fatal('create_or_update_clumps','too many clumps')
 endif

 ! hack for testing
 ntrack = 1
 apr_centre(:,:) = 0.
 apr_centre(1,1) = 10.0
 !apr_centre(2,2) = 10.0

 ! now, reset the regions
 call set_apr_centre(apr_type,apr_centre,ntrack,track_part)

 !print*,'ntrack: ',ntrack

 end subroutine create_or_update_apr_clump

 end module get_apr_level