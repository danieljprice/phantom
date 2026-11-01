!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module hilbertdomain
!
! Domain decomposition along a Peano-Hilbert curve: the particles are
! ordered by their Hilbert key in a cube enclosing all of them, the cube
! is divided into cells (sections) that are key ranges, refined where the
! work is concentrated, and each task gets a contiguous run of sections
! with about the same work.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: hilbert, io, mpiutils, sortutils
!
 use hilbert, only:nbits_key
 implicit none
 integer, parameter :: nbits_total = 3*int(nbits_key)   ! bits in a key (63)

 public :: get_key_cube,choose_sections,print_hilbert_diagnostic

 private

contains

!--------------------------------------------------------------------------
!+
!  the cube enclosing the particles on all tasks (those with h > 0),
!  slightly enlarged so that no particle lies on its upper faces
!+
!--------------------------------------------------------------------------
subroutine get_key_cube(np,xyzh,xmin,dxcube)
 use mpiutils, only:reduceall_mpi
 integer, intent(in)  :: np
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: xmin(3),dxcube
 real    :: xmax(3)
 integer :: i

 xmin =  huge(xmin)
 xmax = -huge(xmax)
 !$omp parallel do default(none) schedule(static) shared(np,xyzh) private(i) &
 !$omp reduction(min:xmin) reduction(max:xmax)
 do i = 1,np
    if (xyzh(4,i) > 0.) then
       xmin = min(xmin,xyzh(1:3,i))
       xmax = max(xmax,xyzh(1:3,i))
    endif
 enddo
 !$omp end parallel do
 xmin = reduceall_mpi('min',xmin)
 xmax = reduceall_mpi('max',xmax)
 dxcube = maxval(xmax - xmin)
 dxcube = dxcube*(1. + 1.e-6) + tiny(dxcube)

end subroutine get_key_cube

!--------------------------------------------------------------------------
!+
!  sections and their owners. The local keys are sorted (keysort) with
!  the running sums of their weights (wcum(0:n), wcum(0)=0). Starting from
!  the octree cells at level lev0, every cell whose work (summed over all
!  tasks) exceeds frac times a task's share is split into its 8 sub-cubes,
!  until none does or the keys run out of bits. The sections, in key
!  order, are then cut into nprocs runs of about equal work.
!
!  Each section is the key range [ilo,ihi]; secw and secn are its work
!  and number of particles summed over all tasks, and owner its task.
!+
!--------------------------------------------------------------------------
subroutine choose_sections(n,keysort,wcum,nprocs,frac,ilo,ihi,secw,secn,owner,nsec,nround)
 use mpiutils, only:reduce_in_place_mpi
 integer,                 intent(in)  :: n,nprocs
 integer(8),              intent(in)  :: keysort(:)
 real(8),                 intent(in)  :: wcum(0:)
 real,                    intent(in)  :: frac
 integer(8), allocatable, intent(out) :: ilo(:),ihi(:)
 real(8),    allocatable, intent(out) :: secw(:),secn(:)
 integer,    allocatable, intent(out) :: owner(:)
 integer,                 intent(out) :: nsec,nround
 integer(8), allocatable :: newlo(:),newhi(:)
 real(8),    allocatable :: sums(:,:)
 integer(8) :: width,step
 real(8)    :: wtot,wmax,wshare,wbefore
 integer    :: lev0,k,j,m,nsplit,i1,i2,nnew

 ! starting level: at least 8 cells per task
 lev0 = 1
 do while (8**lev0 < 8*nprocs .and. 3*lev0 < nbits_total - 3)
    lev0 = lev0 + 1
 enddo
 nsec = 8**lev0
 width = shiftl(1_8,nbits_total - 3*lev0)
 allocate(ilo(nsec),ihi(nsec))
 do k = 1,nsec
    ilo(k) = (k-1)*width
    ihi(k) = ilo(k) + width - 1
 enddo

 nround = 0
 refine: do
    nround = nround + 1
    ! the work and number of particles in each section, summed over all tasks
    if (allocated(sums)) deallocate(sums)
    allocate(sums(2,nsec))
    do k = 1,nsec
       i1 = nkeys_below(n,keysort,ilo(k))
       if (ihi(k) == huge(ihi(k))) then
          i2 = n
       else
          i2 = nkeys_below(n,keysort,ihi(k)+1)
       endif
       sums(1,k) = wcum(i2) - wcum(i1)
       sums(2,k) = real(i2 - i1,kind=8)
    enddo
    call reduce_in_place_mpi('+',sums)
    wtot = sum(sums(1,:))
    wmax = frac*wtot/nprocs
    ! split the sections with too much work into their 8 sub-cubes
    nsplit = 0
    do k = 1,nsec
       if (sums(1,k) > wmax .and. ihi(k) > ilo(k)) nsplit = nsplit + 1
    enddo
    if (nsplit == 0) exit refine
    nnew = nsec + 7*nsplit
    allocate(newlo(nnew),newhi(nnew))
    m = 0
    do k = 1,nsec
       if (sums(1,k) > wmax .and. ihi(k) > ilo(k)) then
          step = (ihi(k) - ilo(k) + 1)/8
          do j = 0,7
             m = m + 1
             newlo(m) = ilo(k) + j*step
             newhi(m) = newlo(m) + step - 1
          enddo
       else
          m = m + 1
          newlo(m) = ilo(k)
          newhi(m) = ihi(k)
       endif
    enddo
    call move_alloc(newlo,ilo)
    call move_alloc(newhi,ihi)
    nsec = nnew
 enddo refine

 ! cut the sections, in key order, into nprocs runs of about equal work:
 ! each section goes to the task whose share contains its middle
 allocate(secw(nsec),secn(nsec),owner(nsec))
 secw = sums(1,1:nsec)
 secn = sums(2,1:nsec)
 wshare  = max(wtot/nprocs,tiny(wtot))
 wbefore = 0.
 do k = 1,nsec
    owner(k) = min(nprocs-1,int((wbefore + 0.5*secw(k))/wshare))
    wbefore  = wbefore + secw(k)
 enddo

end subroutine choose_sections

!--------------------------------------------------------------------------
!+
!  the number of sorted keys below kval (binary search)
!+
!--------------------------------------------------------------------------
pure integer function nkeys_below(n,keysort,kval)
 integer,    intent(in) :: n
 integer(8), intent(in) :: keysort(:),kval
 integer :: lo,hi,mid

 lo = 0
 hi = n
 do while (lo < hi)
    mid = (lo + hi)/2
    if (keysort(mid+1) < kval) then
       lo = mid + 1
    else
       hi = mid
    endif
 enddo
 nkeys_below = lo

end function nkeys_below

!--------------------------------------------------------------------------
!+
!  diagnostic: the decomposition along the Hilbert curve that the work
!  weights wt (per particle, h > 0) would give, against the current one
!+
!--------------------------------------------------------------------------
subroutine print_hilbert_diagnostic(np,xyzh,wt,iprint)
 use io,       only:id,master,nprocs
 use mpiutils, only:reduce_in_place_mpi
 use hilbert,  only:get_hilbert_keys
 use sortutils,only:indexx
 integer, intent(in) :: np,iprint
 real,    intent(in) :: xyzh(:,:),wt(:)
 integer(8), allocatable :: key(:),keysort(:),ilo(:),ihi(:)
 integer,    allocatable :: indx(:),owner(:)
 real(8),    allocatable :: wcum(:),secw(:),secn(:),pred(:,:),now(:,:)
 real    :: xmin(3),dxcube
 integer :: i,n,k,nsec,nround

 call get_key_cube(np,xyzh,xmin,dxcube)
 allocate(key(np),indx(np))
 call get_hilbert_keys(np,xyzh,xmin,dxcube,key)
 call indexx(np,key,indx)
 ! sorted keys of the particles with h > 0, and the running sum of their weights
 allocate(keysort(np),wcum(0:np))
 n = 0
 wcum(0) = 0.
 do k = 1,np
    i = indx(k)
    if (xyzh(4,i) <= 0.) cycle
    n = n + 1
    keysort(n) = key(i)
    wcum(n) = wcum(n-1) + wt(i)
 enddo

 call choose_sections(n,keysort,wcum,nprocs,0.25,ilo,ihi,secw,secn,owner,nsec,nround)

 ! work and particles per task: as cut, and as now
 allocate(pred(2,nprocs),now(2,nprocs))
 pred = 0.
 do k = 1,nsec
    pred(1,owner(k)+1) = pred(1,owner(k)+1) + secw(k)
    pred(2,owner(k)+1) = pred(2,owner(k)+1) + secn(k)
 enddo
 now = 0.
 now(1,id+1) = wcum(n)
 now(2,id+1) = n
 call reduce_in_place_mpi('+',now)

 if (id == master) then
    write(iprint,"(a,i8,a,i3,a,i3)") ' Hilbert decomposition: ',nsec,' sections in ',nround, &
         ' rounds, finest octree level ',maxval((nbits_total - trailz(ihi-ilo+1))/3)
    write(iprint,"(a,2(2x,a,f6.2,a,f6.2))") ' Hilbert decomposition max/mean:', &
         'now work',maxval(now(1,:))/(sum(now(1,:))/nprocs),' particles',maxval(now(2,:))/(sum(now(2,:))/nprocs), &
         'cut work',maxval(pred(1,:))/(sum(pred(1,:))/nprocs),' particles',maxval(pred(2,:))/(sum(pred(2,:))/nprocs)
    write(iprint,"(a)") ' task   particles now  particles cut   work now/mean   work cut/mean'
    do k = 1,nprocs
       write(iprint,"(i5,2(i15),2(f16.3))") k-1,nint(now(2,k)),nint(pred(2,k)), &
            now(1,k)/(sum(now(1,:))/nprocs),pred(1,k)/(sum(pred(1,:))/nprocs)
    enddo
 endif

end subroutine print_hilbert_diagnostic

end module hilbertdomain
