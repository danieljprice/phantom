!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module neighcompact
!
! Compacted (CSR) neighbour lists shared by matrix-free particle solvers.
! The list is symmetric: a pair is kept when either particle lies inside
! the kernel support of the other.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, io, kdtree, kernel, neighkdtree, omputils,
!   part, physcon
!
 use io, only:fatal
 implicit none
 public :: compact_neigh_counts,get_compacted_neighbour_list,is_solver_particle

 private

contains
!-----------------------------------------------------------------------
!+
!  Size of the compacted neighbour storage.
!  icompactmax follows the same factor used by the implicit radiation
!  solver: a generous multiple of the mean neighbour count.
!+
!-----------------------------------------------------------------------
subroutine compact_neigh_counts(npart,radkern,hfact,nneigh_average,icompactmax)
 use physcon, only:pi
 integer, intent(in)  :: npart
 real,    intent(in)  :: radkern,hfact
 integer, intent(out) :: nneigh_average,icompactmax

 nneigh_average = int(4./3.*pi*(radkern*hfact)**3) + 1
 icompactmax = int(1.2*10.*nneigh_average*npart)

end subroutine compact_neigh_counts

!-----------------------------------------------------------------------
!+
!  Build a compacted neighbour list for active gas and boundary particles.
!  ivar(1,n) = neighbour count
!  ivar(2,n) = offset into ijvar (the first neighbour is at offset+1)
!  ivar(3,n) = particle index
!  ijvar     = neighbour indices
!
!  If nedges is present, running out of space in ijvar is not fatal:
!  the list is left incomplete and nedges returns the size required,
!  so the caller can reallocate and call again.
!
!  Optionally, neighbours not in xyzh (e.g. particles on other MPI tasks)
!  are added from the candidates of each leaf, icand(icand_start(icell):
!  icand_start(icell+1)-1), with positions xyzh_extra; they are stored in
!  ijvar as joffset + their index in xyzh_extra.
!+
!-----------------------------------------------------------------------
subroutine get_compacted_neighbour_list(xyzh,ivar,ijvar,ncompact,ncompactlocal,nedges, &
                                        icand_start,icand,xyzh_extra,joffset)
 use dim,         only:periodic,maxpsph
 use neighkdtree, only:ncells,get_neighbour_list,listneigh,leaf_is_active
 use kdtree,      only:inodeparts,inoderange
 use boundary,    only:mod_boundary
 use kernel,      only:radkern2
 use omputils,    only:omp_chunk
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(out) :: ivar(:,:),ijvar(:)
 integer, intent(out) :: ncompact,ncompactlocal
 integer, intent(out), optional :: nedges
 integer, intent(in),  optional :: icand_start(:),icand(:),joffset
 real,    intent(in),  optional :: xyzh_extra(:,:)
 integer                    :: icell,nneigh,i,j,k,n,ip,ncompact_private
 integer                    :: icompact_private,icompact,icompactmax,nneigh_trial,ncand,nchunk
 integer, parameter         :: maxcellcache = 10000
 integer, save, allocatable :: neighlist(:)
 real                       :: dx,dy,dz,xi,yi,zi,hi,hi21,hj1,rij2,q2i,q2j
 real, save, allocatable    :: xyzcache(:,:)
 !$omp threadprivate(xyzcache,neighlist)
 logical                    :: growable,extra

 if (.not. allocated(neighlist)) then
    !$omp parallel default(none) shared(xyzh)
    allocate(neighlist(size(xyzh,2)),xyzcache(4,maxcellcache))
    !$omp end parallel
 endif

 ncompact = 0
 ncompactlocal = 0
 icompact = 0
 icompactmax = size(ijvar)
 growable = present(nedges)
 extra = present(icand_start) .and. present(icand) .and. present(xyzh_extra) .and. present(joffset)
 ! chunks of up to 64 cells: one at a time contends for iterations, larger chunks load
 ! balance poorly, and small trees still need several chunks per thread
 nchunk = omp_chunk(int(ncells),64)
 !$omp parallel do default(none) schedule(dynamic,nchunk)&
 !$omp shared(nchunk,ncells,xyzh,inodeparts,inoderange,leaf_is_active)&
 !$omp shared(ivar,ijvar,ncompact,icompact,icompactmax,maxpsph,growable)&
 !$omp shared(extra,icand_start,icand,xyzh_extra,joffset)&
 !$omp private(icell,i,j,k,n,ip,dx,dy,dz,rij2,q2i,q2j)&
 !$omp private(hi,xi,yi,zi,hi21,hj1,ncompact_private,icompact_private,nneigh_trial,nneigh,ncand)

 over_cells: do icell=1,int(ncells)

    !--skip empty cells and inactive cells
    if (leaf_is_active(icell) <= 0) cycle over_cells

    !
    !--get the neighbour list and fill the cell cache
    !
    call get_neighbour_list(icell,listneigh,nneigh_trial,xyzh,xyzcache,maxcellcache,getj=.true.)
    !
    !--a particle can have more neighbours, local plus on other tasks, than
    !  there are particles on this task, so grow this thread's list if needed
    !
    ncand = 0
    if (extra) ncand = icand_start(icell+1) - icand_start(icell)
    if (nneigh_trial + ncand > size(neighlist)) then
       deallocate(neighlist)
       allocate(neighlist(2*(nneigh_trial + ncand)))
    endif

    over_parts: do ip = inoderange(1,icell),inoderange(2,icell)
       i = abs(inodeparts(ip))  ! the tree stores inactive particles as -i
       if (.not.is_solver_particle(i,xyzh)) cycle over_parts

       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       hi = xyzh(4,i)
       nneigh = 0
       hi21 = 1./hi**2

       loop_over_neigh: do n = 1,nneigh_trial

          j = listneigh(n)
          !--skip self; the divergence and gradient are pairwise
          if (j==i .or. j>maxpsph) cycle loop_over_neigh

          if (n <= maxcellcache) then
             ! positions from cache are already mod boundary
             dx = xi - xyzcache(1,n)
             dy = yi - xyzcache(2,n)
             dz = zi - xyzcache(3,n)
             hj1 = xyzcache(4,n)
          else
             dx = xi - xyzh(1,j)
             dy = yi - xyzh(2,j)
             dz = zi - xyzh(3,j)
             hj1 = 1./xyzh(4,j)
          endif
          if (periodic) call mod_boundary(dx,dy,dz)
          rij2 = dx*dx + dy*dy + dz*dz
          q2i = rij2*hi21
          q2j = rij2*hj1*hj1
          !
          !--keep the pair if either particle lies inside the other's kernel
          !
          is_sph_neighbour: if (q2i < radkern2 .or. q2j < radkern2) then
             nneigh = nneigh + 1
             neighlist(nneigh) = j
          endif is_sph_neighbour

       enddo loop_over_neigh
       !
       !--candidates not in xyzh (e.g. on other MPI tasks)
       !
       if (extra) then
          do k = icand_start(icell),icand_start(icell+1)-1
             n = icand(k)
             dx = xi - xyzh_extra(1,n)
             dy = yi - xyzh_extra(2,n)
             dz = zi - xyzh_extra(3,n)
             if (periodic) call mod_boundary(dx,dy,dz)
             rij2 = dx*dx + dy*dy + dz*dz
             q2i = rij2*hi21
             q2j = rij2/xyzh_extra(4,n)**2
             if (q2i < radkern2 .or. q2j < radkern2) then
                nneigh = nneigh + 1
                neighlist(nneigh) = joffset + n
             endif
          enddo
       endif

       !$omp critical(listcompact)
       ncompact = ncompact + 1
       ncompact_private = ncompact
       icompact_private = icompact
       icompact = icompact + nneigh
       !$omp end critical (listcompact)
       ivar(1,ncompact_private) = nneigh
       ivar(2,ncompact_private) = icompact_private
       ivar(3,ncompact_private) = i
       if (icompact_private+nneigh > icompactmax) then
          ! out of space: keep counting so the caller can resize
          if (growable) cycle over_parts
          call fatal('neighcompact','not enough memory allocated for neighbour list', &
                     var='icompactmax',ival=icompactmax)
       endif

       do k = 1, nneigh
          j = neighlist(k)
          ijvar(icompact_private + k) = j
       enddo
    enddo over_parts
 enddo over_cells
 !$omp end parallel do

 ncompactlocal = ncompact
 if (growable) nedges = icompact

end subroutine get_compacted_neighbour_list

!-----------------------------------------------------------------------
!+
!  whether particle i is solved for: active gas particles and
!  boundary particles that are not dead or accreted
!+
!-----------------------------------------------------------------------
logical function is_solver_particle(i,xyzh)
 use dim,  only:maxphase,maxp,maxpsph
 use part, only:iphase,igas,iboundary,get_partinfo,isdead_or_accreted
 integer, intent(in) :: i
 real,    intent(in) :: xyzh(:,:)
 logical :: iactivei,iamgasi,iamdusti
 integer :: iamtypei

 is_solver_particle = .false.
 !--sinks are stored beyond maxpsph
 if (i <= 0 .or. i > maxpsph) return
 if (maxphase==maxp) then
    call get_partinfo(iphase(i),iactivei,iamgasi,iamdusti,iamtypei)
 else
    iactivei = .true.
    iamtypei = igas
    iamgasi  = .true.
 endif
 !--skip inactive and non-gas particles, but keep boundary particles
 if ((.not.iactivei .or. .not.iamgasi) .and. .not.iamtypei==iboundary) return
 is_solver_particle = .not.isdead_or_accreted(xyzh(4,i))

end function is_solver_particle

end module neighcompact
