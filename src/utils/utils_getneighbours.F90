!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module getneighbours
!
! A set of routines generate neighbour lists for all particles, and
!  read/write the list to file
!
! :References:
!
! :Owner: James Wurster
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, kdtree, kernel, linklist, part
!
 implicit none

 public :: generate_neighbour_lists, neighbours_stats, read_neighbours, write_neighbours
 integer, public, allocatable, dimension(:)   :: neighcount
 integer, public, allocatable, dimension(:,:) :: neighb
 real,    public            :: meanneigh, sdneigh, neighcrit
 logical                    :: neigh_overload
 integer,         parameter :: maxcellcache =  50000
 integer,         parameter :: neighall     = 100000 ! maximum number of neighbours to test
 integer, public, parameter :: neighmax     =   2000 ! maximum number of neighbours to store

 private

contains

!-----------------------------------------------------------------------
!+
! Generate neighbour lists for all particles
!+
!-----------------------------------------------------------------------
subroutine generate_neighbour_lists(xyzh,vxyzu,npart,dumpfile,write_neighbour_list)
 use dim,      only:maxp
 use kernel,   only:radkern2
 use linklist, only:ncells, ifirstincell, set_linklist, get_neighbour_list
 use part,     only:get_partinfo, igas, maxphase, iphase, iamboundary, iamtype
 use kdtree,   only:inodeparts,inoderange
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 real,             intent(in)     :: xyzh(:,:),vxyzu(:,:)
 integer,          intent(in)     :: npart
 character(len=*), intent(in)     :: dumpfile
 logical,          intent(in)     :: write_neighbour_list
 real,allocatable, dimension(:,:) :: dumxyzh

 integer      :: i,j,k,p,ip,icell,ineigh,nneigh,dummynpart
 integer      :: ineigh_all(neighall)
 real         :: dx,dy,dz,rij2
 real         :: hi1,hj1,hi21,hj21,q2i,q2j
 integer, allocatable, save :: listneigh(:)
 real,    allocatable, save :: xyzcache(:,:)
 real         :: rneigh_all(neighall)
 !$omp threadprivate(xyzcache,listneigh)
 character(len=100) :: neighbourfile

 if (.not.allocated(listneigh)) allocate(listneigh(maxp))
 if (.not.allocated(xyzcache)) allocate(xyzcache(maxcellcache,4))

 !****************************************
 ! 1. Build kdtree and linklist
 ! --> global (shared) neighbour lists for all particles in tree cell
 !****************************************

 print*, 'Building kdtree and linklist: '
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 dummynpart = npart
 call set_linklist(dummynpart,npart,dumxyzh,vxyzu(:,1:npart))

 print*, 'Allocating arrays for neighbour storage : '
 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))

 neighcount(:) = 0
 neighb(:,:)   = 0

 print "(A,I5)", 'Maximum neighbour number allocated:  ', neighmax

 !***************************************
 ! 2. Assign neighbour lists to particles by searching shared list of host cell
 !***************************************

 print*, 'Creating neighbour lists for particles'

 !$omp parallel default(none) &
 !$omp shared(ncells,ifirstincell,npart,maxphase,maxp,inodeparts,inoderange) &
 !$omp shared(xyzh,vxyzu,iphase,neighcount,neighb) &
#ifdef PERIODIC
 !$omp shared(dxbound,dybound,dzbound) &
#endif
 !$omp private(icell,i,j,k,p,ip)&
 !$omp private(nneigh,ineigh_all,rneigh_all) &
 !$omp private(hi1,hi21,hj1,hj21,rij2,q2i,q2j) &
 !$omp private(dx,dy,dz)
 !$omp do schedule(runtime)
 over_cells: do icell=1,int(ncells)

    k = ifirstincell(icell)

    ! Skip empty/inactive cells
    if (k <= 0) cycle over_cells

    ! Get neighbour list for the cell
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,maxcellcache,getj=.true.)

    ! Loop over particles incellsn the cell
    over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
       i = inodeparts(ip)
       if (maxphase==maxp) then
          if (iamboundary( iamtype(iphase(i)) )) cycle over_parts
       endif

       ! Fill neighbour list for this particle
       neighcount(i) = 0
       ineigh_all    = 0
       rneigh_all    = 0.
       hi1  = 1.0/xyzh(4,i)
       hi21 = hi1*hi1

       over_neighbours: do ineigh = 1,nneigh
          j = abs(listneigh(ineigh))

          ! Skip self
          if (i==j) cycle over_neighbours

          dx = xyzh(1,i) - xyzh(1,j)
          dy = xyzh(2,i) - xyzh(2,j)
          dz = xyzh(3,i) - xyzh(3,j)
#ifdef PERIODIC
          if (abs(dx) > 0.5*dxbound) dx = dx - dxbound*SIGN(1.0,dx)
          if (abs(dy) > 0.5*dybound) dy = dy - dybound*SIGN(1.0,dy)
          if (abs(dz) > 0.5*dzbound) dz = dz - dzbound*SIGN(1.0,dz)
#endif
          rij2 = dx*dx + dy*dy + dz*dz
          q2i  = rij2*hi21

          hj1  = 1.0/xyzh(4,j)
          hj21 = hj1*hj1
          q2j  = rij2*hj21

          is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then
             neighcount(i) = neighcount(i) + 1
             if (neighcount(i) <= neighall) then
                ineigh_all(neighcount(i)) = j
                rneigh_all(neighcount(i)) = rij2
             else
                print*, 'neighbour finding.  neighcount > neighall.  aborting'
                stop
             endif
          endif is_sph_neighbour
       enddo over_neighbours ! End loop over neighbours
       ! Failsafe if too many neighbours
       if (neighcount(i) <= neighmax) then
          neighb(i,1:neighcount(i)) = ineigh_all(1:neighcount(i))
       else
          print*, 'Neighbour finding: There are ',neighcount(i),' neighbours for i = ',i
          print*, 'Neighbour finding: Keeping the ',neighmax,' closest neighbours.'
          neighcount(i) = neighmax
          do p = 1,neighmax
             j = minloc(rneigh_all(1:neighcount(i)),1)
             neighb(i,p) = ineigh_all(j)
             rneigh_all(j) = huge(rij2)
          enddo
       endif
    enddo over_parts         ! End loop over particles in the cell
 enddo over_cells            ! End loop over cells in the kd-tree
 !$omp enddo
 !$omp end parallel

 ! Do some simple stats on neighbour numbers
 meanneigh = 0.0
 sdneigh   = 0.0
 neighcrit = 0.0
 call neighbours_stats(npart)

 !**************************************
 ! 3. Output neighbour lists to file (if requested; these files can become very big)
 !**************************************
 if (write_neighbour_list) then
    neighbourfile = 'neigh_'//trim(dumpfile)
    call write_neighbours(neighbourfile, npart)
    print*, 'Neighbour finding complete for file ', trim(dumpfile)
 endif

 deallocate(dumxyzh)
end subroutine generate_neighbour_lists
!-----------------------------------------------------------------------
!+
! Calculates the mean and standard deviation of the neighbour number
! Also calculates a 5 sigma deviation from meanneigh
! (This is principally used as a diagnostic aid for structure finding
! algorithms that rely on the nearest neighbours, like CLUMPFIND)
!+
!-----------------------------------------------------------------------
subroutine neighbours_stats(npart)
 integer, intent(in) :: npart
 integer             :: ipart, minimum, maximum

 ! Calculate mean and standard deviation of neighbour counts

 maximum = maxval(neighcount)
 minimum = minval(neighcount)
 print*, 'The maximum neighbour count is ', maximum
 print*, 'The minimum neighbour count is ', minimum

 if (maximum > neighmax) then
    print*, 'WARNING! Neighbour count too large for allocated arrays.  aborting'
    stop
 endif

 meanneigh = sum(neighcount)/real(npart)
 sdneigh   = 0.0

!$omp parallel default(none) &
!$omp shared(neighcount,meanneigh,npart) &
!$omp private(ipart) &
!$omp reduction(+:sdneigh)
!$omp do schedule(runtime)
 do ipart=1,npart
    sdneigh = sdneigh + (neighcount(ipart)-meanneigh)**2
 enddo
 !$omp enddo
 !$omp end parallel

 sdneigh = sqrt(sdneigh/real(npart))

 print*, 'Mean neighbour number is ', meanneigh
 print*, 'Standard Deviation: ', sdneigh

end subroutine neighbours_stats
!-----------------------------------------------------------------------
!+
! Reads in a pre-written neighbours file
!+
!-----------------------------------------------------------------------
subroutine read_neighbours(neighbourfile,npart)
 integer,            intent(in) :: npart
 character(len=100), intent(in) :: neighbourfile
 integer                        :: i,j,neighcheck,tolcheck

 neigh_overload = .false.
 allocate(neighcount(npart))
 allocate(neighb(npart,neighmax))
 neighcount(:) = 0
 neighb(:,:)   = 0

 print*, 'Reading neighbour file ', trim(neighbourfile)
 open(2,file= neighbourfile,  form = 'UNFORMATTED')
 read(2)  neighcheck, tolcheck, meanneigh,sdneigh,neighcrit
 if (neighcheck/=neighmax) print*, 'WARNING: mismatch in neighmax: ', neighmax, neighcheck
 read(2) (neighcount(i), i=1,npart)
 do i=1,npart
    if (neighcount(i) > neighmax) then
       neigh_overload = .true.
       read(2) (neighb(i,j), j=1,neighmax)
    else
       read(2) (neighb(i,j), j=1,neighcount(i))
    endif
 enddo
 close(2)

 call neighbours_stats(npart)

 if (neigh_overload) then
    print*, 'WARNING! File Read incomplete: neighbour count exceeds array size.  aborting.'
    stop
 else
    print*, 'File Read Complete'
 endif

end subroutine read_neighbours
!-----------------------------------------------------------------------
!+
! Writes neighbour data to binary file
!+
!-----------------------------------------------------------------------
subroutine write_neighbours(neighbourfile,npart)
 integer, intent(in) :: npart
 integer             :: i,j
 character(len=100)  :: neighbourfile
 real, parameter     :: tolerance = 2.0e0  ! A dummy parameter used to keep file format similar to other codes (Probably delete later)

 neigh_overload = .false.
 neighbourfile  = trim(neighbourfile)
 print*, 'Writing neighbours to file ', neighbourfile

 open(2,file=neighbourfile,form='unformatted')
 write(2)  neighmax, tolerance, meanneigh,sdneigh,neighcrit
 write(2) (neighcount(i), i=1,npart)
 do i=1,npart
    if (neighcount(i) > neighmax) then
       neigh_overload = .true.
       print*, 'neighbour overload: ', neighcount(i), neighmax
       write(2) (neighb(i,j), j=1,neighmax)
    else
       write(2) (neighb(i,j), j=1,neighcount(i))
    endif
 enddo
 close(2)

 if (neigh_overload) then
    print*, 'WARNING! File write incomplete: neighbour count exceeds array size'
 else
    print*, 'File Write Complete'
 endif

end subroutine write_neighbours
!-----------------------------------------------------------------------
end module getneighbours
