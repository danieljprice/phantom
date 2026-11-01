!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpisections
!
! Domain decomposition by sections of a geometric tree: a box enclosing
! all particles is split at the midpoint of its longest axis, and so on
! down, as in the tree used for gravity. The path of a particle down this
! tree is its key, so every cell of the tree is a range of keys and the
! depth first order of the cells is the order of their keys. The cells
! where the work is concentrated are refined (sections), and each task
! gets a run of sections, in key order, with about the same work.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, mpiutils
!
 implicit none
 integer, parameter :: nbits_total = 63    ! bits in a key: levels of the tree
 real,    parameter :: frac_section = 0.25 ! largest work in a section, as a fraction of a task's share
 real,    parameter :: tol_balance  = 1.1  ! new sections only if a task has more work than this times the mean
 real,    parameter :: frac_margin  = 0.1  ! margin of the box of the sections around the particles, each axis

 public :: get_key_box,get_geometric_keys,choose_sections
 public :: assign_section_owners,local_section_index,owned_section_keys,print_section_diagnostic

 private

 ! the current decomposition: the box and the split axis at each level,
 ! the section boundaries, and the run of sections owned by this task
 ! (kfirst_own to kfirst_own+nown-1)
 integer(8), allocatable :: ilo_cur(:)
 integer, allocatable :: owner_cur(:)
 real    :: boxmin_cur(3),boxsize_cur(3)
 integer :: axis_cur(nbits_total)
 integer :: nsec_cur = 0,kfirst_own = 1,nown = 0
 integer :: nlev_cur = nbits_total  ! depth of the deepest section: keys are only needed to there

contains

!--------------------------------------------------------------------------
!+
!  the box enclosing the particles on all tasks (those with h > 0),
!  slightly enlarged so that no particle lies on its upper faces, and the
!  axis along which the cells of each level are split: their longest
!  (all the cells of a level have the same shape), the lowest on ties
!+
!--------------------------------------------------------------------------
subroutine get_key_box(np,xyzh,boxmin,boxsize,axis)
 use mpiutils, only:reduceall_mpi
 integer, intent(in)  :: np
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: boxmin(3),boxsize(3)
 integer, intent(out) :: axis(nbits_total)
 real    :: xmax(3),dims(3)
 integer :: i,d

 boxmin =  huge(boxmin)
 xmax   = -huge(xmax)
 !$omp parallel do default(none) schedule(static) shared(np,xyzh) private(i) &
 !$omp reduction(min:boxmin) reduction(max:xmax)
 do i = 1,np
    if (xyzh(4,i) > 0.) then
       boxmin = min(boxmin,xyzh(1:3,i))
       xmax   = max(xmax,xyzh(1:3,i))
    endif
 enddo
 !$omp end parallel do
 boxmin = reduceall_mpi('min',boxmin)
 xmax   = reduceall_mpi('max',xmax)
 boxsize = (xmax - boxmin)*(1. + 1.e-6) + tiny(boxsize)

 dims = boxsize
 do d = 1,nbits_total
    axis(d) = maxloc(dims,1)
    dims(axis(d)) = 0.5*dims(axis(d))
 enddo

end subroutine get_key_box

!--------------------------------------------------------------------------
!+
!  the key of each position xyz(1:3,i): its path down the tree of
!  midpoint splits of the box, one bit per level (1 for the upper half
!  along axis(d)). Positions outside the box are moved onto its surface
!+
!--------------------------------------------------------------------------
subroutine get_geometric_keys(np,xyz,boxmin,boxsize,axis,key,nlevels)
 integer,    intent(in)  :: np
 real,       intent(in)  :: xyz(:,:),boxmin(3),boxsize(3)
 integer,    intent(in)  :: axis(nbits_total)
 integer(8), intent(out) :: key(:)
 integer,    intent(in), optional :: nlevels  ! only the first nlevels bits, the rest zero
 real    :: x(3),lo(3),hi(3),xmid
 integer :: i,d,ia,nlev

 nlev = nbits_total
 if (present(nlevels)) nlev = min(nlevels,nbits_total)
 !$omp parallel do default(none) schedule(static) shared(np,xyz,boxmin,boxsize,axis,key,nlev) &
 !$omp private(i,d,ia,x,lo,hi,xmid)
 do i = 1,np
    lo = boxmin
    hi = boxmin + boxsize
    x  = min(max(xyz(1:3,i),lo),hi)
    key(i) = 0
    do d = 1,nlev
       ia = axis(d)
       xmid = 0.5*(lo(ia) + hi(ia))
       if (x(ia) >= xmid) then
          key(i) = ior(shiftl(key(i),1),1_8)
          lo(ia) = xmid
       else
          key(i) = shiftl(key(i),1)
          hi(ia) = xmid
       endif
    enddo
    key(i) = shiftl(key(i),nbits_total-nlev)
 enddo
 !$omp end parallel do

end subroutine get_geometric_keys

!--------------------------------------------------------------------------
!+
!  sections and their owners, from the keys and work weights wt of the n
!  particles on this task (in any order). Starting from the cells at
!  level 3*lev0, every cell whose work (summed over all tasks) exceeds
!  frac times a task's share is split into the 8 cells three levels
!  below it, until none does or the keys run out of bits. The sections,
!  in key order, are then cut into nprocs runs of about equal work.
!
!  Each section is the key range [ilo,ihi]; secw and secn are its work
!  and number of particles summed over all tasks, and owner its task.
!+
!--------------------------------------------------------------------------
subroutine choose_sections(n,key,wt,nprocs,frac,ilo,ihi,secw,secn,owner,nsec,nround)
 use mpiutils, only:reduce_in_place_mpi
 integer,                 intent(in)  :: n,nprocs
 integer(8),              intent(in)  :: key(:)
 real,                    intent(in)  :: wt(:)
 real,                    intent(in)  :: frac
 integer(8), allocatable, intent(out) :: ilo(:),ihi(:)
 real(8),    allocatable, intent(out) :: secw(:),secn(:)
 integer,    allocatable, intent(out) :: owner(:)
 integer,                 intent(out) :: nsec,nround
 integer(8), allocatable :: newlo(:),newhi(:)
 real(8),    allocatable :: sums(:,:)
 integer(8) :: width,step
 real(8)    :: wtot,wmax,wshare,wbefore
 integer    :: lev0,k,j,m,nsplit,nnew,i

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
    sums = 0.
    !$omp parallel do default(none) schedule(static) shared(n,key,wt,nsec,ilo) private(i,k) &
    !$omp reduction(+:sums)
    do i = 1,n
       k = find_section(nsec,ilo,key(i))
       sums(1,k) = sums(1,k) + wt(i)
       sums(2,k) = sums(2,k) + 1.
    enddo
    !$omp end parallel do
    call reduce_in_place_mpi('+',sums)
    wtot = sum(sums(1,:))
    wmax = frac*wtot/nprocs
    ! split the sections with too much work into the 8 cells three levels below
    nsplit = 0
    do k = 1,nsec
       if (sums(1,k) > wmax .and. ihi(k) - ilo(k) >= 7) nsplit = nsplit + 1
    enddo
    if (nsplit == 0) exit refine
    nnew = nsec + 7*nsplit
    allocate(newlo(nnew),newhi(nnew))
    m = 0
    do k = 1,nsec
       if (sums(1,k) > wmax .and. ihi(k) - ilo(k) >= 7) then
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
!  a new decomposition in the box with split axes axis, for the work
!  weights wt (per particle, those with h > 0 only): the keys of the
!  particles, and the sections and their owners
!+
!--------------------------------------------------------------------------
subroutine new_sections(np,xyzh,wt,boxmin,boxsize,axis,key,ilo,ihi,secw,secn,owner,nsec,nround)
 use io, only:nprocs
 integer,                 intent(in)  :: np
 real,                    intent(in)  :: xyzh(:,:),wt(:)
 real,                    intent(in)  :: boxmin(3),boxsize(3)
 integer,                 intent(in)  :: axis(nbits_total)
 integer(8),              intent(out) :: key(:)
 integer(8), allocatable, intent(out) :: ilo(:),ihi(:)
 real(8),    allocatable, intent(out) :: secw(:),secn(:)
 integer,    allocatable, intent(out) :: owner(:)
 integer,                 intent(out) :: nsec,nround
 integer(8), allocatable :: keylive(:)
 real,       allocatable :: wtlive(:)
 integer :: i,n

 call get_geometric_keys(np,xyzh,boxmin,boxsize,axis,key)
 ! keys and weights of the particles with h > 0
 allocate(keylive(np),wtlive(np))
 n = 0
 do i = 1,np
    if (xyzh(4,i) <= 0.) cycle
    n = n + 1
    keylive(n) = key(i)
    wtlive(n)  = wt(i)
 enddo

 call choose_sections(n,keylive,wtlive,nprocs,frac_section,ilo,ihi,secw,secn,owner,nsec,nround)

end subroutine new_sections

!--------------------------------------------------------------------------
!+
!  the task each particle belongs to (ibelong), from the section its key
!  falls in. The sections are chosen for the work weights wt on the
!  first call, when the particles leave the box of the sections or fill
!  less than half of it along an axis, and when rebalance is set and the
!  work of the tasks under the current sections is out of balance;
!  otherwise they are kept. Dead or accreted particles stay on this task.
!  Optionally the sink particles (xyzmh_ptmass, held by every task) are
!  given the owner of their section, or -1 if they have been accreted,
!  the same on every task
!+
!--------------------------------------------------------------------------
subroutine assign_section_owners(np,xyzh,wt,rebalance,ibelong,nptmass,xyzmh_ptmass,ibelong_sink)
 use io, only:id
 integer,           intent(in)    :: np
 real,              intent(in)    :: xyzh(:,:),wt(:)
 logical,           intent(in)    :: rebalance
 integer,           intent(inout) :: ibelong(:)
 integer, optional, intent(in)    :: nptmass
 real,    optional, intent(in)    :: xyzmh_ptmass(:,:)
 integer, optional, intent(out)   :: ibelong_sink(:)
 integer(8), allocatable :: key(:),ihi(:),keysink(:)
 real(8),    allocatable :: secw(:),secn(:)
 real    :: boxmin(3),boxsize(3)
 integer :: axis(nbits_total),i,nround
 logical :: new

 allocate(key(np))
 ! new sections if the particles have left the box, or fill much less of it
 call get_key_box(np,xyzh,boxmin,boxsize,axis)
 new = (nsec_cur == 0)
 if (.not.new) new = any(boxmin < boxmin_cur) .or. any(boxmin + boxsize > boxmin_cur + boxsize_cur) &
                     .or. any(boxsize < 0.5*boxsize_cur)
 if (.not.new) then
    call get_geometric_keys(np,xyzh,boxmin_cur,boxsize_cur,axis_cur,key,nlev_cur)
    call set_owners(np,xyzh,key,ibelong)
    ! when asked to rebalance, keep the sections unless the work is out of balance
    if (rebalance) new = (work_imbalance(np,xyzh,wt,ibelong) > tol_balance)
 endif
 if (new) then
    ! a margin around the particles, so that the sections last while they spread out
    boxmin_cur  = boxmin - 0.5*frac_margin*boxsize
    boxsize_cur = boxsize*(1. + frac_margin)
    axis_cur    = axis
    call new_sections(np,xyzh,wt,boxmin_cur,boxsize_cur,axis_cur,key,ilo_cur,ihi,secw,secn,&
                      owner_cur,nsec_cur,nround)
    nown       = count(owner_cur == id)
    kfirst_own = max(findloc(owner_cur,id,1),1)
    nlev_cur   = maxval(nbits_total - trailz(ihi - ilo_cur + 1))
    call set_owners(np,xyzh,key,ibelong)
 endif

 if (present(nptmass) .and. present(xyzmh_ptmass) .and. present(ibelong_sink)) then
    if (nptmass > 0) then
       allocate(keysink(nptmass))
       call get_geometric_keys(nptmass,xyzmh_ptmass,boxmin_cur,boxsize_cur,axis_cur,keysink,nlev_cur)
       do i = 1,nptmass
          if (xyzmh_ptmass(4,i) > 0.) then
             ibelong_sink(i) = owner_cur(find_section(nsec_cur,ilo_cur,keysink(i)))
          else
             ibelong_sink(i) = -1
          endif
       enddo
    endif
 endif

end subroutine assign_section_owners

!--------------------------------------------------------------------------
!+
!  the task each particle belongs to, from the current sections and the
!  keys of the particles. Dead or accreted particles stay on this task
!+
!--------------------------------------------------------------------------
subroutine set_owners(np,xyzh,key,ibelong)
 use io, only:id
 integer,    intent(in)    :: np
 real,       intent(in)    :: xyzh(:,:)
 integer(8), intent(in)    :: key(:)
 integer,    intent(inout) :: ibelong(:)
 integer :: i

 !$omp parallel do default(none) schedule(static) shared(np,xyzh,key,ilo_cur,owner_cur,nsec_cur,ibelong,id) &
 !$omp private(i)
 do i = 1,np
    if (xyzh(4,i) > 0.) then
       ibelong(i) = owner_cur(find_section(nsec_cur,ilo_cur,key(i)))
    else
       ibelong(i) = id
    endif
 enddo
 !$omp end parallel do

end subroutine set_owners

!--------------------------------------------------------------------------
!+
!  the largest work of a task, over the mean, if the particles moved to
!  the tasks given by ibelong
!+
!--------------------------------------------------------------------------
real function work_imbalance(np,xyzh,wt,ibelong)
 use io,       only:nprocs
 use mpiutils, only:reduce_in_place_mpi
 integer, intent(in) :: np
 real,    intent(in) :: xyzh(:,:),wt(:)
 integer, intent(in) :: ibelong(:)
 real(8) :: work(nprocs)
 integer :: i

 work = 0.
 do i = 1,np
    if (xyzh(4,i) > 0.) work(ibelong(i)+1) = work(ibelong(i)+1) + wt(i)
 enddo
 call reduce_in_place_mpi('+',work,nprocs)
 work_imbalance = real(maxval(work)/max(sum(work)/nprocs,tiny(1.d0)))

end function work_imbalance

!--------------------------------------------------------------------------
!+
!  the index, among the sections this task owns (1 to nsecown, in key
!  order), of the section of each position xyz(1:3,i), from the last call
!  to assign_section_owners. Positions outside them (none, if called after
!  the particles moved to their owners) are put in the nearest one
!+
!--------------------------------------------------------------------------
subroutine local_section_index(n,xyz,isec,nsecown)
 integer, intent(in)  :: n
 real,    intent(in)  :: xyz(:,:)
 integer, intent(out) :: isec(:),nsecown
 integer(8), allocatable :: key(:)
 integer :: i

 nsecown = nown
 if (nown <= 1) then
    isec(1:n) = 1
    return
 endif
 allocate(key(n))
 call get_geometric_keys(n,xyz,boxmin_cur,boxsize_cur,axis_cur,key,nlev_cur)
 !$omp parallel do default(none) schedule(static) shared(n,key,isec,nsec_cur,ilo_cur,kfirst_own,nown) private(i)
 do i = 1,n
    isec(i) = min(max(find_section(nsec_cur,ilo_cur,key(i)) - kfirst_own + 1,1),nown)
 enddo
 !$omp end parallel do

end subroutine local_section_index

!--------------------------------------------------------------------------
!+
!  the first key of each section this task owns, in key order (from the
!  last call to assign_section_owners)
!+
!--------------------------------------------------------------------------
subroutine owned_section_keys(keyfirst)
 integer(8), allocatable, intent(out) :: keyfirst(:)

 allocate(keyfirst(max(nown,1)))
 keyfirst = 0
 if (nown > 0) keyfirst(1:nown) = ilo_cur(kfirst_own:kfirst_own+nown-1)

end subroutine owned_section_keys

!--------------------------------------------------------------------------
!+
!  the section containing key kval: the last one starting at or below it
!  (the sections cover all keys, in increasing order)
!+
!--------------------------------------------------------------------------
pure integer function find_section(nsec,ilo,kval)
 integer,    intent(in) :: nsec
 integer(8), intent(in) :: ilo(:),kval
 integer :: lo,hi,mid

 lo = 1
 hi = nsec
 do while (lo < hi)
    mid = (lo + hi + 1)/2
    if (ilo(mid) <= kval) then
       lo = mid
    else
       hi = mid - 1
    endif
 enddo
 find_section = lo

end function find_section

!--------------------------------------------------------------------------
!+
!  diagnostic: the work and particles of each task now, against those of
!  a new decomposition for the work weights wt (per particle, h > 0)
!+
!--------------------------------------------------------------------------
subroutine print_section_diagnostic(np,xyzh,wt,iprint)
 use io,       only:id,master,nprocs
 use mpiutils, only:reduce_in_place_mpi
 integer, intent(in) :: np,iprint
 real,    intent(in) :: xyzh(:,:),wt(:)
 integer(8), allocatable :: key(:),ilo(:),ihi(:)
 integer,    allocatable :: owner(:)
 real(8),    allocatable :: secw(:),secn(:),pred(:,:),now(:,:)
 real    :: boxmin(3),boxsize(3)
 integer :: axis(nbits_total),k,nsec,nround

 allocate(key(np))
 call get_key_box(np,xyzh,boxmin,boxsize,axis)
 call new_sections(np,xyzh,wt,boxmin,boxsize,axis,key,ilo,ihi,secw,secn,owner,nsec,nround)

 ! work and particles per task: as cut, and as now
 allocate(pred(2,nprocs),now(2,nprocs))
 pred = 0.
 do k = 1,nsec
    pred(1,owner(k)+1) = pred(1,owner(k)+1) + secw(k)
    pred(2,owner(k)+1) = pred(2,owner(k)+1) + secn(k)
 enddo
 now = 0.
 now(1,id+1) = sum(wt(1:np),mask=(xyzh(4,1:np) > 0.))
 now(2,id+1) = count(xyzh(4,1:np) > 0.)
 call reduce_in_place_mpi('+',now)

 if (id == master) then
    write(iprint,"(a,i8,a,i3,a,i3)") ' MPI sections: ',nsec,' in ',nround, &
         ' rounds, deepest level ',maxval(nbits_total - trailz(ihi-ilo+1))
    write(iprint,"(a,2(2x,a,f6.2,a,f6.2))") ' MPI sections max/mean:', &
         'now work',maxval(now(1,:))/(sum(now(1,:))/nprocs),' particles',maxval(now(2,:))/(sum(now(2,:))/nprocs), &
         'cut work',maxval(pred(1,:))/(sum(pred(1,:))/nprocs),' particles',maxval(pred(2,:))/(sum(pred(2,:))/nprocs)
    write(iprint,"(a)") ' task   particles now  particles cut   work now/mean   work cut/mean'
    do k = 1,nprocs
       write(iprint,"(i5,2(i15),2(f16.3))") k-1,nint(now(2,k)),nint(pred(2,k)), &
            now(1,k)/(sum(now(1,:))/nprocs),pred(1,k)/(sum(pred(1,:))/nprocs)
    enddo
 endif

end subroutine print_section_diagnostic

end module mpisections
