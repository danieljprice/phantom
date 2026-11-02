!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpisections
!
! Domain decomposition by sections of a geometric tree, built as the tree
! used for gravity: each cell is split at the midpoint of the bounds of
! its particles (on all tasks) along their longest axis, where the work is
! concentrated. The leaves of this tree are the sections. A cell's key is
! its path from the root (one bit per level, left aligned), so the depth
! first order of the sections is the order of their keys, and each task
! gets a run of sections, in that order, with about the same work.
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
 real,    parameter :: frac_maxp    = 0.75 ! most particles a cut puts on a task, as a fraction of the
                                           ! particle arrays: room for migration, drift and imports

 ! a tree of sections: for each cell, its split axis and position, its
 ! left child (the right one is next to it; 0 for a section), its key and
 ! depth, and for a section its index; then the sections in key order,
 ! with their cells, work, number of particles and owners
 type sectiontree
    integer :: ncell = 0, nsec = 0
    integer,    allocatable :: axis(:),left(:),depth(:),isec(:)
    real,       allocatable :: xsplit(:)
    integer(8), allocatable :: key(:)
    integer,    allocatable :: seccell(:),owner(:)
    real(8),    allocatable :: secw(:),secn(:)
 end type sectiontree

 public :: sectiontree,build_sections,locate_section,cut_sections
 public :: assign_section_owners,local_section_index,owned_section_keys,print_section_diagnostic
 public :: get_section_tree,get_work_weights

 private

 ! the current decomposition, and the run of sections owned by this task
 ! (kfirst_own to kfirst_own+nown-1)
 type(sectiontree) :: cur
 integer :: kfirst_own = 1,nown = 0

 ! each task's share of the work in a cut (mean 1), corrected from the
 ! measured busy time of the tasks; their busy time and the checks since
 ! the last cut
 real, allocatable :: share(:),busy_cut(:)
 integer :: nchecks = 0

contains

!--------------------------------------------------------------------------
!+
!  build the tree of sections for the work weights wt of the particles
!  (those with h > 0): starting from the root, every cell whose work,
!  summed over all tasks, exceeds frac_section times a task's share is
!  split, one level per round, at the midpoint of the bounds of its
!  particles along their longest axis (as is every cell with more than
!  frac_section of a task's share of the particles). The sections, in key
!  order, are then cut into nprocs runs (cut_sections), with at most
!  frac_maxp of the particle arrays (the size of xyzh) on each task
!+
!--------------------------------------------------------------------------
subroutine build_sections(np,xyzh,wt,st,nround)
 use io,       only:nprocs
 use mpiutils, only:reduce_in_place_mpi,reduceall_mpi
 integer,           intent(in)  :: np
 real,              intent(in)  :: xyzh(:,:),wt(:)
 type(sectiontree), intent(out) :: st
 integer,           intent(out) :: nround
 integer, allocatable :: icell(:),cand(:),pos(:),newcand(:)
 real(8), allocatable :: sums(:,:),bnd(:,:),cellw(:),celln(:)
 real(8) :: wtot,wmax,ntot,nmax,npmax
 real    :: ext(3)
 integer :: i,j,k,c,nc,nnew,ia

 call grow_cells(st,64)
 allocate(cellw(size(st%axis)),celln(size(st%axis)))
 st%ncell = 1
 st%left(1)  = 0
 st%depth(1) = 0
 st%key(1)   = 0
 allocate(icell(np))
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(np,xyzh,icell) &
 !$omp private(i)
 do i = 1,np
    icell(i) = merge(1,0,xyzh(4,i) > 0.)
 enddo
 !$omp end parallel do

 ! the cells created in the last round, whose work and bounds are needed
 allocate(cand(1),pos(1))
 cand(1) = 1
 nc = 1
 wtot = 0.
 ntot = 0.
 nround = 0
 refine: do
    nround = nround + 1
    if (size(pos) < st%ncell) then
       deallocate(pos)
       allocate(pos(size(st%axis)))
    endif
    pos(1:st%ncell) = 0
    do j = 1,nc
       pos(cand(j)) = j
    enddo
    ! work, number of particles and bounds (min of x and of -x) of each candidate
    allocate(sums(2,nc),bnd(6,nc))
    sums = 0.
    bnd  = huge(1.d0)
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(np,icell,pos,xyzh,wt) &
    !$omp private(i,j) &
    !$omp reduction(+:sums) &
    !$omp reduction(min:bnd)
    do i = 1,np
       if (icell(i) <= 0) cycle
       j = pos(icell(i))
       if (j <= 0) cycle
       sums(1,j) = sums(1,j) + wt(i)
       sums(2,j) = sums(2,j) + 1.
       bnd(1:3,j) = min(bnd(1:3,j),real(xyzh(1:3,i),kind=8))
       bnd(4:6,j) = min(bnd(4:6,j),-real(xyzh(1:3,i),kind=8))
    enddo
    !$omp end parallel do
    call reduce_in_place_mpi('+',sums)
    bnd = reduceall_mpi('min',bnd)
    if (nround == 1) then
       wtot = sums(1,1)
       ntot = sums(2,1)
    endif
    wmax = frac_section*wtot/nprocs
    nmax = frac_section*ntot/nprocs

    ! split the cells with too much work in two
    if (size(cellw) < st%ncell + 2*nc) call grow_real8(cellw,celln,st%ncell + 2*nc)
    allocate(newcand(2*nc))
    nnew = 0
    do j = 1,nc
       c = cand(j)
       cellw(c) = sums(1,j)
       celln(c) = sums(2,j)
       st%left(c) = 0
       if ((sums(1,j) <= wmax .and. sums(2,j) <= nmax) .or. sums(2,j) < 2. .or. st%depth(c) >= nbits_total) cycle
       ext = real(-bnd(4:6,j) - bnd(1:3,j))
       if (maxval(ext) <= 0.) cycle  ! all its particles at one point
       ia = maxloc(ext,1)
       call grow_cells(st,st%ncell + 2)
       st%axis(c)   = ia
       st%xsplit(c) = real(0.5*(bnd(ia,j) - bnd(ia+3,j)))
       st%left(c)   = st%ncell + 1
       do k = 1,2
          i = st%ncell + k
          st%left(i)  = 0
          st%depth(i) = st%depth(c) + 1
          st%key(i)   = st%key(c)
          if (k == 2) st%key(i) = ibset(st%key(c),nbits_total - 1 - st%depth(c))
          nnew = nnew + 1
          newcand(nnew) = i
       enddo
       st%ncell = st%ncell + 2
    enddo
    if (nnew == 0) exit refine

    ! the particles of the cells that were split move to their children
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(np,icell,st,xyzh) &
    !$omp private(i,c)
    do i = 1,np
       c = icell(i)
       if (c <= 0) cycle
       if (st%left(c) > 0) then
          if (xyzh(st%axis(c),i) >= st%xsplit(c)) then
             icell(i) = st%left(c) + 1
          else
             icell(i) = st%left(c)
          endif
       endif
    enddo
    !$omp end parallel do
    deallocate(sums,bnd,cand)
    call move_alloc(newcand,cand)
    nc = nnew
 enddo refine

 ! the sections are the leaves, in depth first (key) order
 call list_sections(st)
 allocate(st%secw(st%nsec),st%secn(st%nsec),st%owner(st%nsec))
 do k = 1,st%nsec
    st%secw(k) = cellw(st%seccell(k))
    st%secn(k) = celln(st%seccell(k))
 enddo

 ! cut the sections into one run per task, keeping room in the particle arrays
 npmax = frac_maxp*size(xyzh,2)
 call cut_sections(st%nsec,st%secw,st%secn,nprocs,npmax,st%owner,share)

end subroutine build_sections

!--------------------------------------------------------------------------
!+
!  cut n sections, in order, into runs, one per task (owner 0 to nprocs-1),
!  so that the work of the busiest task is as small as possible, with at
!  most npmax particles on a task (unless that cannot be done at all).
!  The least work per run for which fill_runs succeeds is found by
!  bisection, between the largest of a section's work and the mean work
!  per task, and all the work
!+
!--------------------------------------------------------------------------
subroutine cut_sections(n,work,npart,nprocs,npmax,owner,share)
 use io, only:fatal
 integer, intent(in)  :: n,nprocs
 real(8), intent(in)  :: work(n),npart(n),npmax
 integer, intent(out) :: owner(n)
 real,    intent(in), optional :: share(:)
 integer, parameter :: maxiter = 60
 real(8) :: wlo,whi,wmid,ncap
 real    :: s(nprocs)
 integer :: iter

 ! the work of task k is limited to s(k) times the work limit of the cut
 s = 1.
 if (present(share)) then
    if (size(share) == nprocs) s = share
 endif

 ! drop the cap on particles if no cut can meet it
 ncap = npmax
 if (.not.fill_runs(n,work,npart,nprocs,ncap,s,sum(work)/minval(s),owner)) ncap = huge(ncap)

 wlo = max(maxval(work)/maxval(s),sum(work)/sum(s))
 whi = sum(work)/minval(s)
 if (fill_runs(n,work,npart,nprocs,ncap,s,wlo,owner)) return

 do iter = 1,maxiter
    wmid = 0.5*(wlo + whi)
    if (fill_runs(n,work,npart,nprocs,ncap,s,wmid,owner)) then
       whi = wmid
    else
       wlo = wmid
    endif
 enddo

 if (.not.fill_runs(n,work,npart,nprocs,ncap,s,whi,owner)) call fatal('cut_sections','no cut found')

end subroutine cut_sections

!--------------------------------------------------------------------------
!+
!  fill the runs in turn with sections, in order, starting a new run when
!  the next section would take it over s(task) times wmax work or ncap
!  particles, or when only as many sections are left as tasks without one.
!  Returns .false. if that needs more than nprocs runs
!+
!--------------------------------------------------------------------------
logical function fill_runs(n,work,npart,nprocs,ncap,s,wmax,owner)
 integer, intent(in)  :: n,nprocs
 real(8), intent(in)  :: work(n),npart(n),ncap,wmax
 real,    intent(in)  :: s(nprocs)
 integer, intent(out) :: owner(n)
 real(8) :: wrun,nrun
 integer :: k,itask
 logical :: too_much_work,too_many_particles,need_section_per_task

 fill_runs = .false.
 itask = 0
 wrun  = 0.
 nrun  = 0.
 do k = 1,n
    if (k > 1) then
       too_much_work         = (wrun + work(k) > s(itask+1)*wmax)
       too_many_particles    = (nrun + npart(k) > ncap)
       need_section_per_task = (n - k + 1 <= nprocs - 1 - itask)
       if (too_much_work .or. too_many_particles .or. need_section_per_task) then
          itask = itask + 1
          if (itask > nprocs - 1) return
          wrun = 0.
          nrun = 0.
       endif
    endif
    owner(k) = itask
    wrun = wrun + work(k)
    nrun = nrun + npart(k)
 enddo
 fill_runs = .true.

end function fill_runs

!--------------------------------------------------------------------------
!+
!  the sections (leaves of the tree) in depth first order, left first,
!  which is the order of their keys
!+
!--------------------------------------------------------------------------
subroutine list_sections(st)
 type(sectiontree), intent(inout) :: st
 integer :: istack(2*nbits_total+2),nstack,c

 allocate(st%seccell(st%ncell))
 st%isec(1:st%ncell) = 0
 st%nsec = 0
 nstack  = 1
 istack(1) = 1
 do while (nstack > 0)
    c = istack(nstack)
    nstack = nstack - 1
    if (st%left(c) > 0) then
       istack(nstack+1) = st%left(c) + 1
       istack(nstack+2) = st%left(c)
       nstack = nstack + 2
    else
       st%nsec = st%nsec + 1
       st%seccell(st%nsec) = c
       st%isec(c) = st%nsec
    endif
 enddo

end subroutine list_sections

!--------------------------------------------------------------------------
!+
!  the section containing the position x: down the tree from the root,
!  on the side of each split that x is on
!+
!--------------------------------------------------------------------------
pure integer function locate_section(st,x)
 type(sectiontree), intent(in) :: st
 real,              intent(in) :: x(3)
 integer :: c

 c = 1
 do while (st%left(c) > 0)
    if (x(st%axis(c)) >= st%xsplit(c)) then
       c = st%left(c) + 1
    else
       c = st%left(c)
    endif
 enddo
 locate_section = st%isec(c)

end function locate_section

!--------------------------------------------------------------------------
!+
!  the task each particle belongs to (ibelong), from the section it is in.
!  The sections are chosen for the work weights wt on the first call, and
!  again when rebalance is set and either the work of the tasks under the
!  current sections (over their shares) or their busy time since the last
!  cut, if given, is out of balance; otherwise they are kept. A new cut
!  first corrects the shares from the busy times. Dead or accreted
!  particles stay on this task. Optionally the sink particles (xyzmh_ptmass,
!  held by every task) are given the owner of their section, or -1 if they
!  have been accreted, the same on every task
!+
!--------------------------------------------------------------------------
subroutine assign_section_owners(np,xyzh,wt,rebalance,ibelong,nptmass,xyzmh_ptmass,ibelong_sink,busy)
 use io, only:id,nprocs,master,iprint
 integer,           intent(in)    :: np
 real,              intent(in)    :: xyzh(:,:),wt(:)
 logical,           intent(in)    :: rebalance
 integer,           intent(inout) :: ibelong(:)
 integer, optional, intent(in)    :: nptmass
 real,    optional, intent(in)    :: xyzmh_ptmass(:,:)
 integer, optional, intent(out)   :: ibelong_sink(:)
 real,    optional, intent(in)    :: busy(:)  ! busy time of each task since the last check
 integer, parameter :: nchecks_min = 2  ! checks before the busy times are trusted
 integer :: i,nround
 real    :: imb_work,imb_busy
 logical :: new

 if (.not.allocated(share)) then
    allocate(share(nprocs),busy_cut(nprocs))
    share    = 1.
    busy_cut = 0.
 endif
 new = (cur%nsec == 0)
 if (.not.new) then
    call set_owners(np,xyzh,ibelong)
    ! when asked to rebalance, keep the sections unless the work or the
    ! measured busy time is out of balance
    if (rebalance) then
       nchecks  = nchecks + 1
       imb_work = work_imbalance(np,xyzh,wt,ibelong)
       imb_busy = 1.
       if (present(busy)) then
          busy_cut = busy_cut + busy
          if (nchecks >= nchecks_min) imb_busy = busy_imbalance(busy_cut)
       endif
       new = (imb_work > tol_balance .or. imb_busy > tol_balance)
       if (new) then
          if (present(busy) .and. nchecks >= nchecks_min) call update_shares(busy_cut)
          if (id == master) write(iprint,"(a,f6.3,a,f6.3,a,2(1x,f6.3))") ' MPI sections renewed: work',&
             imb_work,' busy',imb_busy,' shares min/max',minval(share),maxval(share)
       endif
    endif
 endif
 if (new) then
    nchecks    = 0
    busy_cut   = 0.
    call build_sections(np,xyzh,wt,cur,nround)
    nown       = count(cur%owner == id)
    kfirst_own = max(findloc(cur%owner,id,1),1)
    call set_owners(np,xyzh,ibelong)
 endif

 if (present(nptmass) .and. present(xyzmh_ptmass) .and. present(ibelong_sink)) then
    do i = 1,nptmass
       if (xyzmh_ptmass(4,i) > 0.) then
          ibelong_sink(i) = cur%owner(locate_section(cur,xyzmh_ptmass(1:3,i)))
       else
          ibelong_sink(i) = -1
       endif
    enddo
 endif

end subroutine assign_section_owners

!--------------------------------------------------------------------------
!+
!  the task each particle belongs to, from the current sections. Dead or
!  accreted particles stay on this task
!+
!--------------------------------------------------------------------------
subroutine set_owners(np,xyzh,ibelong)
 use io, only:id
 integer, intent(in)    :: np
 real,    intent(in)    :: xyzh(:,:)
 integer, intent(inout) :: ibelong(:)
 integer :: i

 !$omp parallel do default(none) schedule(static) &
 !$omp shared(np,xyzh,cur,ibelong,id) &
 !$omp private(i)
 do i = 1,np
    if (xyzh(4,i) > 0.) then
       ibelong(i) = cur%owner(locate_section(cur,xyzh(1:3,i)))
    else
       ibelong(i) = id
    endif
 enddo
 !$omp end parallel do

end subroutine set_owners

!--------------------------------------------------------------------------
!+
!  the largest work of a task over its share of the mean, if the particles
!  moved to the tasks given by ibelong
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
 work_imbalance = real(maxval(work/share)/max(sum(work)/nprocs,tiny(1.d0)))

end function work_imbalance

!--------------------------------------------------------------------------
!+
!  the largest busy time of a task, over the mean
!+
!--------------------------------------------------------------------------
real function busy_imbalance(busy)
 real, intent(in) :: busy(:)

 busy_imbalance = maxval(busy)/max(sum(busy)/size(busy),tiny(1.))

end function busy_imbalance

!--------------------------------------------------------------------------
!+
!  correct the shares of the work from the busy time of each task under
!  the current cut: a task busier than the mean gets a smaller share. The
!  change per cut is limited, as the busy time also depends on the
!  boundaries the cut moves
!+
!--------------------------------------------------------------------------
subroutine update_shares(busy)
 real, intent(in) :: busy(:)
 real, parameter :: fmin = 0.8, fmax = 1.25  ! limits on the change per cut
 real :: bmean

 if (size(busy) /= size(share) .or. minval(busy) <= 0.) return
 bmean = sum(busy)/size(busy)
 share = share*min(max(bmean/busy,fmin),fmax)
 share = share*size(share)/sum(share)

end subroutine update_shares

!--------------------------------------------------------------------------
!+
!  the work weight of each particle: its cost per step, its number of
!  neighbours (cost), times how often it is active
!  (2^ibin with individual timesteps). Particles without a cost yet (before
!  the first density calculation, or new ones) take the mean of the others
!+
!--------------------------------------------------------------------------
subroutine get_work_weights(np,cost,ind_timesteps,ibin,wt)
 integer,         intent(in)  :: np
 real(kind=4),    intent(in)  :: cost(:)
 logical,         intent(in)  :: ind_timesteps
 integer(kind=1), intent(in)  :: ibin(:)
 real,            intent(out) :: wt(:)
 real    :: cost_mean
 integer :: i,nknown

 nknown = count(cost(1:np) > 0.)
 cost_mean = 1.
 if (nknown > 0) cost_mean = real(sum(cost(1:np),mask=(cost(1:np) > 0.)))/nknown

 !$omp parallel do default(none) schedule(static) &
 !$omp shared(np,cost,ind_timesteps,ibin,wt,cost_mean) &
 !$omp private(i)
 do i = 1,np
    wt(i) = cost(i)
    if (wt(i) <= 0.) wt(i) = cost_mean
    if (ind_timesteps) wt(i) = wt(i)*2.**ibin(i)
 enddo
 !$omp end parallel do

end subroutine get_work_weights

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
 integer :: i

 nsecown = nown
 if (nown <= 1) then
    isec(1:n) = 1
    return
 endif
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(n,xyz,isec,cur,kfirst_own,nown) &
 !$omp private(i)
 do i = 1,n
    isec(i) = min(max(locate_section(cur,xyz(1:3,i)) - kfirst_own + 1,1),nown)
 enddo
 !$omp end parallel do

end subroutine local_section_index

!--------------------------------------------------------------------------
!+
!  the key of each section this task owns, in key order (from the last
!  call to assign_section_owners)
!+
!--------------------------------------------------------------------------
subroutine owned_section_keys(keyfirst)
 integer(8), allocatable, intent(out) :: keyfirst(:)

 allocate(keyfirst(max(nown,1)))
 keyfirst = 0
 if (nown > 0) keyfirst(1:nown) = cur%key(cur%seccell(kfirst_own:kfirst_own+nown-1))

end subroutine owned_section_keys

!--------------------------------------------------------------------------
!+
!  the current tree of sections (from the last call to
!  assign_section_owners): for each cell its left child (the right one is
!  next to it; 0 for a section) and, for a section, its index; the owner
!  of each section; and the index of the first section owned by this task
!  (it owns nown in a run from there)
!+
!--------------------------------------------------------------------------
subroutine get_section_tree(ncell,left,isec,nsec,owner,kfirst,nowned)
 integer,              intent(out) :: ncell,nsec,kfirst,nowned
 integer, allocatable, intent(out) :: left(:),isec(:),owner(:)

 ncell  = cur%ncell
 nsec   = cur%nsec
 kfirst = kfirst_own
 nowned = nown
 allocate(left(ncell),isec(ncell),owner(nsec))
 left  = cur%left(1:ncell)
 isec  = cur%isec(1:ncell)
 owner = cur%owner(1:nsec)

end subroutine get_section_tree

!--------------------------------------------------------------------------
!+
!  room for at least n cells in the tree of sections
!+
!--------------------------------------------------------------------------
subroutine grow_cells(st,n)
 type(sectiontree), intent(inout) :: st
 integer,           intent(in)    :: n
 real,       allocatable :: rtmp(:)
 integer(8), allocatable :: ktmp(:)
 integer :: nold,nnew

 nold = 0
 if (allocated(st%axis)) nold = size(st%axis)
 if (n <= nold) return
 nnew = max(n,2*nold)
 call grow_int(st%axis,nold,nnew)
 call grow_int(st%left,nold,nnew)
 call grow_int(st%depth,nold,nnew)
 call grow_int(st%isec,nold,nnew)
 allocate(rtmp(nnew))
 if (nold > 0) rtmp(1:nold) = st%xsplit(1:nold)
 call move_alloc(rtmp,st%xsplit)
 allocate(ktmp(nnew))
 if (nold > 0) ktmp(1:nold) = st%key(1:nold)
 call move_alloc(ktmp,st%key)

end subroutine grow_cells

!--------------------------------------------------------------------------
!+
!  grow an integer array from nold to nnew entries, keeping its values
!+
!--------------------------------------------------------------------------
subroutine grow_int(a,nold,nnew)
 integer, allocatable, intent(inout) :: a(:)
 integer,              intent(in)    :: nold,nnew
 integer, allocatable :: itmp(:)

 allocate(itmp(nnew))
 itmp = 0
 if (nold > 0) itmp(1:nold) = a(1:nold)
 call move_alloc(itmp,a)

end subroutine grow_int

!--------------------------------------------------------------------------
!+
!  room for at least n entries in two real(8) arrays
!+
!--------------------------------------------------------------------------
subroutine grow_real8(a,b,n)
 real(8), allocatable, intent(inout) :: a(:),b(:)
 integer,              intent(in)    :: n
 real(8), allocatable :: tmp(:)
 integer :: nold

 nold = size(a)
 allocate(tmp(max(n,2*nold)))
 tmp(1:nold) = a
 call move_alloc(tmp,a)
 allocate(tmp(max(n,2*nold)))
 tmp(1:nold) = b
 call move_alloc(tmp,b)

end subroutine grow_real8

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
 type(sectiontree) :: st
 real(8), allocatable :: pred(:,:),now(:,:)
 integer :: k,nround

 call build_sections(np,xyzh,wt,st,nround)

 ! work and particles per task: as cut, and as now
 allocate(pred(2,nprocs),now(2,nprocs))
 pred = 0.
 do k = 1,st%nsec
    pred(1,st%owner(k)+1) = pred(1,st%owner(k)+1) + st%secw(k)
    pred(2,st%owner(k)+1) = pred(2,st%owner(k)+1) + st%secn(k)
 enddo
 now = 0.
 now(1,id+1) = sum(wt(1:np),mask=(xyzh(4,1:np) > 0.))
 now(2,id+1) = count(xyzh(4,1:np) > 0.)
 call reduce_in_place_mpi('+',now)

 if (id == master) then
    write(iprint,"(a,i8,a,i3,a,i3)") ' MPI sections: ',st%nsec,' in ',nround, &
         ' rounds, deepest level ',maxval(st%depth(st%seccell(1:st%nsec)))
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
