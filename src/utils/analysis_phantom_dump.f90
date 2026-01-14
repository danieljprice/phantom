!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! print stellar radius and mass computed with respect to the
! location of the particle with maximum density
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: readwrite_dumps, sortutils, units
!
 implicit none
 character(len=*), parameter, public :: analysistype = 'rstar_and_mstar'
 public :: do_analysis

 private

contains
!----------------------------------------------------------------
!+
!  print stellar radius and mass computed with respect to the
!  location of the particle with maximum density
!+
!----------------------------------------------------------------
subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
 use units,           only:udist,umass
 use readwrite_dumps, only:opened_full_dump
 use sortutils,       only:sort_by_radius
 real,     intent(in) :: pmass,time
 real,     intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer,  intent(in) :: numfile,npart,iunit
 character(len=*), intent(in) :: dumpfile
 integer :: i,j,location
 real    :: xpos(3),pos(3)
 integer, allocatable :: iorder(:)
 real,    allocatable :: radius(:)

 ! if dumpfile is not a full dump, skip it
 if (.not.opened_full_dump) then
    write(*,'("SKIPPING FILE -- (Not a full dump)")')
    return
 endif

 ! print the analysis being done

 ! sort particles by radius
 allocate(iorder(npart),radius(npart))
 call sort_by_radius(npart,xyzh,iorder)

 ! find the particle with the highest density
 location = minloc(xyzh(4,:),dim=1)
 xpos(:) = xyzh(1:3,location)

 do j = 1, npart
    i  = iorder(j) ! access the rank of each particle in radius.

    ! the position of the particle is calculated with respect to the particle of highest density.
    ! xyzh is position wrt the black hole assumed to be at the origin
    pos(:) = xyzh(1:3,i) - xpos(:)

    ! calculate the radial position of the particle
    radius(j) = sqrt(dot_product(pos(:),pos(:)))
 enddo
 print*,'max radius =',maxval(radius)*udist,'cm; total mass =', npart*pmass*umass,'g'

end subroutine do_analysis

end module analysis
