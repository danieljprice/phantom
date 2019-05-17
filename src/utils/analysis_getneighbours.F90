!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine which computes neighbour lists for all particles
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, kernel, linklist, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'getneighbours'

 integer, parameter :: neighmax = 350
 integer, parameter :: maxcellcache = 50000

 integer, allocatable, dimension(:) :: neighcount
 integer, allocatable, dimension(:,:) :: neighb

 real :: meanneigh, sdneigh, neighcrit

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use dim, only: maxp, maxneigh
 use kernel, only: radkern2
 use linklist, only: ncells, ifirstincell, set_linklist, get_neighbour_list
 use part, only: get_partinfo, igas, iboundary,maxphase, ll, iphase
 implicit none

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 real,allocatable,dimension(:,:) :: dumxyzh

 integer :: i, j, iamtypei, icell, ineigh, nneigh
 real :: dx,dy,dz, rij2
 real :: hi1,hj1,hi21,hj21, q2i, q2j

 integer,save :: listneigh(maxneigh)
 real, save:: xyzcache(4,maxcellcache)
 !$omp threadprivate(xyzcache,listneigh)
 real :: fgrav(20)

 logical :: iactivei, iamdusti,iamgasi, ifilledcellcache

 character(100) :: neighbourfile

 !****************************************
 ! 1. Build kdtree and linklist
 ! --> global (shared) neighbour lists for all particles in tree cell
 !****************************************

 print*, 'Building kdtree and linklist: '

 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 print*, '- Done'

 print*, 'Allocating arrays for neighbour storage : '

 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))

 neighcount(:) = 0
 neighb(:,:) = 0

 print*, '- Done'
 print "(A,I3)", 'Maximum neighbour number allocated:  ', neighmax

 !***************************************
 ! 2. Assign neighbour lists to particles by searching shared list of host cell
 !***************************************

 print*, 'Creating neighbour lists for particles'

! Loop over cells

 !$omp parallel default(none) &
 !$omp shared(ncells,ll,ifirstincell,npart) &
 !$omp shared(xyzh,vxyzu,iphase) &
 !$omp shared(neighcount,neighb) &
 !$omp private(icell,i, j)&
 !$omp private(iamtypei,iamgasi,iamdusti,iactivei) &
 !$omp private(ifilledcellcache,nneigh) &
 !$omp private(hi1,hi21,hj1,hj21,rij2,q2i,q2j) &
 !$omp private(fgrav, dx,dy,dz)
 !$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)

    i = ifirstincell(icell)

    ! Skip empty/inactive cells
    if(i<=0) cycle over_cells

    ! Get neighbour list for the cell

#ifdef GRAVITY
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.,f=fgrav)
#else
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)
#endif
    ifilledcellcache = .true.

    ! Loop over particles in the cell

    over_parts: do while(i /=0)
       !print*, i, icell, ncells
       if(i<0) then ! i<0 indicates inactive particles
          i = ll(abs(i))
          cycle over_parts
       endif

       if(maxphase==maxp) then
          call get_partinfo(iphase(i), iactivei,iamdusti,iamtypei)
          iamgasi = (iamtypei ==igas)
       else
          iactivei = .true.
          iamtypei = igas
          iamdusti = .false.
          iamgasi = .true.
       endif

       ! Catches case where first particle is inactive
       if (.not.iactivei) then
          i = ll(i)
          cycle over_parts
       endif

       ! do not compute neighbours for boundary particles
       if(iamtypei ==iboundary) cycle over_parts


       ! Fill neighbour list for this particle

       neighcount(i) = 0

       over_neighbours: do ineigh = 1,nneigh
          !print*, i,ineigh, listneigh(ineigh)
          j = abs(listneigh(ineigh))

          ! Skip self-references
          if(i==j) cycle over_neighbours

          dx = xyzh(1,i) - xyzh(1,j)
          dy = xyzh(2,i) - xyzh(2,j)
          dz = xyzh(3,i) - xyzh(3,j)

          rij2 = dx*dx + dy*dy +dz*dz

          hi1 = 1.0/xyzh(4,i)
          hi21 = hi1*hi1

          q2i = rij2*hi21

          hj1 = 1.0/xyzh(4,j)
          hj21 = hj1*hj1
          q2j = rij2*hj21

          is_sph_neighbour: if(q2i < radkern2 .or. q2j < radkern2) then
             !$omp critical
             neighcount(i) = neighcount(i) + 1
             if(neighcount(i) <=neighmax) neighb(i,neighcount(i)) = j
             !$omp end critical
          endif is_sph_neighbour

       enddo over_neighbours
       ! End loop over neighbours

       i = ll(i)
    enddo over_parts
    ! End loop over particles in the cell

 enddo over_cells
 !$omp enddo
 !$omp end parallel

 ! End loop over cells in the kd-tree

 ! Do some simple stats on neighbour numbers

 meanneigh = 0.0
 sdneigh = 0.0
 neighcrit = 0.0

 call neighbours_stats(npart)

 !**************************************
 ! 3. Output neighbour lists to file
 !**************************************

 neighbourfile = 'neigh_'//TRIM(dumpfile)
 call write_neighbours(neighbourfile, npart)

 print*, 'Neighbour finding complete for file ', TRIM(dumpfile)

 ! Deallocate memory for next dump
 deallocate(dumxyzh, neighcount,neighb)

end subroutine do_analysis

!--------------------------------------------------------------------
!+
! Calculates the mean and standard deviation of the neighbour number
! Also calculates a 5 sigma deviation from meanneigh
! (This is principally used as a diagnostic aid for structure finding
! algorithms that rely on the nearest neighbours, like CLUMPFIND)
!+
!--------------------------------------------------------------------
subroutine neighbours_stats(npart)

 implicit none
 integer, intent(in) :: npart
 integer :: ipart

 real :: minimum, maximum

 ! Calculate mean and standard deviation of neighbour counts

 maximum = maxval(neighcount)
 minimum = minval(neighcount)
 print*, 'The maximum neighbour count is ', maximum
 print*, 'The minimum neighbour count is ', minimum

 if(maximum > neighmax) then
    print*, 'WARNING! Neighbour count too large for allocated arrays'
 endif

 meanneigh = sum(neighcount)/REAL(npart)
 sdneigh = 0.0

!$omp parallel default(none) &
!$omp shared(neighcount,meanneigh,npart)&
!$omp private(ipart) &
!$omp reduction(+:sdneigh)
!$omp do schedule(runtime)
 do ipart=1,npart
    sdneigh = sdneigh+(neighcount(ipart)-meanneigh)**2
 enddo
 !$omp enddo
 !$omp end parallel

 sdneigh = sqrt(sdneigh/REAL(npart))

 print*, 'Mean neighbour number is ', meanneigh
 print*, 'Standard Deviation: ', sdneigh
 neighcrit = meanneigh-5.0*sdneigh

 print*, 'Clumps created only if neighbour number greater than ', neighcrit

 return
end subroutine neighbours_stats


!--------------------------------------------------------------------
!+
! Writes neighbour data to binary file
!+
!--------------------------------------------------------------------
subroutine write_neighbours(neighbourfile,npart)

 implicit none

 integer, intent(in) :: npart

 integer :: i,j
 character(100)::neighbourfile
 ! This is a dummy parameter, used to keep file format similar to other codes
 ! (Will probably delete this later)

 real, parameter :: tolerance = 2.0e0
 logical :: neigh_overload

 neigh_overload = .false.

 neighbourfile = TRIM(neighbourfile)

 print*, 'Writing neighbours to file ', neighbourfile

 OPEN (2, file=neighbourfile, form='unformatted')

 write(2)  neighmax, tolerance, meanneigh,sdneigh,neighcrit
 write(2) (neighcount(i), i=1,npart)
 do i=1,npart
    if(neighcount(i) > neighmax) then
       neigh_overload = .true.
       write(2) (neighb(i,j), j=1,neighmax)
    else
       write(2) (neighb(i,j), j=1,neighcount(i))
    endif
 enddo

 close(2)


 if(neigh_overload) then
    print*, 'WARNING! File write incomplete: neighbour count exceeds array size'
 else
    print*, 'File Write Complete'
 endif

 return
end subroutine write_neighbours

end module
