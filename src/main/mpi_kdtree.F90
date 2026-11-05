!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpikdtree
!
! The tree across MPI tasks: the domain decomposition into sections of the
! geometric tree, the global tree above the sections, and the local tree
! extended with the nodes of other tasks (the combined tree): the nodes copied
! into the global tree (extend_local_tree, icall=0), then the locally
! essential tree imported for the force calculation (icall=1)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, io, kdtree,
!   kernel, mpibalance, mpidomain, mpisections, mpitree, mpiutils,
!   omputils, part, sortutils, timing
!
 use dim,          only:ncellsmax,gravity
 use dtypekdtree,  only:kdnode,node,leaf_is_active,active_leaves,nactive_leaves,ncells
 use kdtree,       only:irootnode,maxdepth,tree_accuracy,maketree,push_pair
 implicit none

 integer, parameter :: maxsec_per_task    = 64    ! room for the top nodes, per task
 integer, parameter :: nmirror_min        = 256   ! fewest nodes a task copies into the global tree
 integer, parameter :: nactive_per_mirror = 4     ! active particles per copied node
 integer, parameter :: nlet_split         = 8192  ! pairs each task's LET descent is split into
 integer, parameter :: nprofile           = 24    ! entries in the size profile of a global leaf
 integer, parameter :: nlet_record = 9 + merge(16,0,gravity) ! reals per imported node

 ! the global tree, the same on every task: the top nodes (1 to ntop) above the sections,
 ! then each task's block of nodes copied from its local tree
 integer,      allocatable, public :: cellatid(:)     ! owner+1 of each leaf, 0 otherwise
 type(kdnode), allocatable, public :: nodeglobal(:)
 integer,      allocatable :: global_owner(:)         ! task owning each node, -1 above the sections
 integer,      allocatable :: global_local_index(:)   ! the owner's local index (negative for a leaf)

 ! the sections, the leaves of the top nodes
 integer, public :: ntop = 0,nsec_top = 0
 integer, allocatable, public :: top_leaf(:),top_owner(:),top_root(:) ! global index, owner, owner's root
 integer, allocatable :: own_lo(:),own_hi(:)           ! range of owners of the sections below each top node
 integer, allocatable :: section_of(:)                 ! section of each top node, 0 if not a leaf
 integer, allocatable :: sec_start(:),sec_by_owner(:)  ! the sections of each task, in task order
 integer, allocatable :: own_roots(:)                  ! the roots of this task's sections
 integer :: nown_roots = 0

 ! the local nodes this task copies into the global tree
 integer, allocatable :: nodemap(:)                       ! index of each in the global tree, else 0
 integer, allocatable :: mirror_list(:)                   ! in the order of this task's block
 integer, allocatable :: nmirror_task(:),mirror_displ(:)  ! size and offset of each task's block
 integer :: nmirror = 0

 ! the combined tree: the local tree, then the global tree (node g at ncells+g, except this
 ! task's nodes, which are its local nodes), then the nodes imported from other tasks (LET)
 integer, public :: ncomb = 0                              ! last node, 0 if none
 integer, public :: irootcomb = 1                          ! root
 integer, allocatable, public :: remote_leaf_owner(:)      ! owner+1 of other tasks' leaves, else 0
 integer, allocatable, public :: remote_local_index(:)     ! the owner's local index (negative for a leaf)
 integer :: ncomb_tree = 0                                 ! last node before the imported nodes
 integer :: ncells_appended = 0                            ! where the global tree was appended
 integer, allocatable :: let_parent(:)                     ! global leaves given imported children
 integer :: nlet_parent = 0
 integer, allocatable :: let_local(:)                      ! global leaves whose particles were imported
 integer :: nlet_local = 0
 real    :: reach_ref = 0.                                 ! typical reach of the active leaves, all tasks
 integer, allocatable :: let_base(:)                       ! position before the nodes from each task

 ! the selection of the nodes to send each task
 type intlist
    integer, allocatable :: x(:)
    integer :: n = 0
 end type intlist
 type(intlist), allocatable :: let_opened(:,:)   ! nodes the walks of each task could open, per thread
 integer :: let_round = 0                        ! marks in let_opened: -(let_round + task + 1)
 integer, allocatable :: let_mark(:)             ! nodes opened for the current receiver (per thread)
 integer :: let_stamp = 0
 integer, allocatable :: own_gleaves(:)          ! this task's leaves of the global tree
 integer :: nown_gleaves = 0
 integer, allocatable :: branch_mark(:)          ! local nodes above active leaves
 integer :: branch_stamp = 0

 ! for each global leaf of another task, the largest size and hmax of a node its walks can
 ! pair with, by the steps taken below the leaf: profile_len(g) values from profile_start(g)
 ! in profile_pool (0 if no active particles, -1 if not received)
 real,    allocatable :: profile_pool(:)
 integer, allocatable :: profile_start(:),profile_len(:)
 logical, allocatable :: global_active(:)          ! a global node with active leaves below
 integer, allocatable :: prof_set(:),inactive_set(:) ! changed by the last call, to undo
 integer :: nprof_set = 0,ninactive_set = 0

 ! statistics of the import since the last summary
 integer(kind=8) :: nlet_sent = 0,nlet_recv = 0
 integer(kind=8) :: nimp_leaves_recv = 0,nimp_part_recv = 0

!$omp threadprivate(let_mark,let_stamp)

 public :: maketreeglobal,extend_local_tree,sync_hmax_mpi,print_let_stats
 public :: allocate_mpikdtree,deallocate_mpikdtree

 private

contains

!-----------------------------------------------------------------------
!+
!  extend the local tree with the nodes of other tasks. icall=0, after a
!  build: copy the top of the local tree into the global tree and append
!  the global tree to the local tree (the combined tree). icall=1, before
!  force: import the nodes and leaves of other tasks that the dual tree
!  walks of this task could open(the locally essential tree)
!+
!-----------------------------------------------------------------------
subroutine extend_local_tree(icall,nactive)
 use io,     only:fatal
 use timing, only:get_timings,increment_timer,itimer_tree_graft
 integer, intent(in) :: icall
 integer, intent(in), optional :: nactive  ! active particles on this task (icall=0)
 real(kind=4) :: t1,t2,tcpu1,tcpu2

 select case(icall)
 case(0)
    if (.not.present(nactive)) call fatal('extend_local_tree','nactive needed with icall=0')
    call get_timings(t1,tcpu1)
    call set_top_owners()
    call graft_local_trees(nactive)
    call combine_trees()
    call get_timings(t2,tcpu2)
    call increment_timer(itimer_tree_graft,t2-t1,tcpu2-tcpu1)
 case(1)
    call let_exchange()
 case default
    call fatal('extend_local_tree','unknown icall',ival=icall)
 end select

end subroutine extend_local_tree

!-----------------------------------------------------------------------
!+
!  allocate the global tree and the arrays of the combined tree. On
!  input nnodes and ncache are the sizes of the tree arrays and walk
!  caches for the local tree; with more than one MPI task they are
!  increased to make room for the combined tree
!+
!-----------------------------------------------------------------------
subroutine allocate_mpikdtree(nnodes,ncache)
 use allocutils, only:allocate_array
 use io,         only:nprocs
 integer, intent(inout) :: nnodes,ncache
 integer(kind=8) :: nglobal

 ! the top nodes, then each task's block of nodes
 nglobal = ntop_max() + int(nprocs,kind=8)*mirror_budget()

 ! room for the combined tree: the global tree appended after the local tree, then as
 ! many nodes again imported from other tasks. The walk caches only the nodes on the
 ! branches of local leaves: the local nodes and the top levels of the global tree
 if (nprocs > 1) then
    nnodes = nnodes + 2*int(nglobal)
    ncache = ncache + ntop_max()
 endif

 call allocate_array('cellatid',          cellatid,          nglobal)
 call allocate_array('nodeglobal',        nodeglobal,        nglobal)
 call allocate_array('nodemap',           nodemap,           ncellsmax+1)
 call allocate_array('remote_leaf_owner', remote_leaf_owner, nnodes)
 call allocate_array('remote_local_index',remote_local_index,nnodes)
 call allocate_array('global_local_index',global_local_index,nglobal)
 call allocate_array('global_owner',      global_owner,      nglobal)
 if (nprocs > 1) then
    call allocate_array('profile_start',profile_start,int(nglobal))
    call allocate_array('profile_len',  profile_len,  int(nglobal))
    allocate(global_active(nglobal))
    call allocate_array('branch_mark',  branch_mark,  ncellsmax+1)
    branch_mark = 0
 endif

 nodemap = 0
 nmirror = 0
 remote_leaf_owner  = 0
 remote_local_index = 0

end subroutine allocate_mpikdtree

!-----------------------------------------------------------------------
!+
!  deallocate the global tree and the arrays of the combined tree
!+
!-----------------------------------------------------------------------
subroutine deallocate_mpikdtree()

 if (allocated(cellatid))           deallocate(cellatid)
 if (allocated(nodeglobal))         deallocate(nodeglobal)
 if (allocated(nodemap))            deallocate(nodemap)
 if (allocated(remote_leaf_owner))  deallocate(remote_leaf_owner)
 if (allocated(remote_local_index)) deallocate(remote_local_index)
 if (allocated(global_local_index)) deallocate(global_local_index)
 if (allocated(global_owner))       deallocate(global_owner)
 if (allocated(profile_pool))       deallocate(profile_pool)
 if (allocated(profile_start))      deallocate(profile_start,profile_len,global_active)
 if (allocated(branch_mark))        deallocate(branch_mark)
 if (allocated(mirror_list))        deallocate(mirror_list)
 if (allocated(nmirror_task))       deallocate(nmirror_task,mirror_displ)
 if (allocated(top_leaf))           deallocate(top_leaf,top_owner,top_root,own_lo,own_hi,section_of,&
                                               sec_start,sec_by_owner)
 if (allocated(own_roots))          deallocate(own_roots)
 if (allocated(own_gleaves))        deallocate(own_gleaves)
 if (allocated(let_parent))         deallocate(let_parent)
 if (allocated(let_local))          deallocate(let_local)
 if (allocated(let_base))           deallocate(let_base)
 if (allocated(let_opened))         deallocate(let_opened)
 if (allocated(prof_set))           deallocate(prof_set)
 if (allocated(inactive_set))       deallocate(inactive_set)
!$omp parallel
 if (allocated(let_mark)) deallocate(let_mark)
!$omp end parallel

 ! no tree until the next build
 ntop          = 0
 nsec_top      = 0
 nown_roots    = 0
 nmirror       = 0
 ncomb         = 0
 ncomb_tree    = 0
 nlet_parent   = 0
 nlet_local    = 0
 nown_gleaves  = 0
 nprof_set     = 0
 ninactive_set = 0

end subroutine deallocate_mpikdtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the global tree (MPI): the particles are moved to the
!  task that owns their section of the geometric tree, each task builds its
!  local tree from its sections, and the cells of the geometric tree above
!  the sections are built from the section roots of all tasks, identically
!  on every task. On substeps with few active particles, the local trees are
!  instead revised for the new positions, keeping the owners (revised=.true.)
!+
!-------------------------------------------------------------------------------
subroutine maketreeglobal(xyzh,np,nactive,apr_tree,revised,nptmass,xyzmh_ptmass)
 use mpibalance,   only:balancedomains
 use mpisections,  only:assign_section_owners
 use mpiutils,     only:reduceall_mpi
 use part,         only:ibelong,ibin,isdead_or_accreted
 use boundary,     only:cross_boundary
 use mpidomain,    only:isperiodic
 use timing,       only:increment_timer,get_timings,itimer_balance,itimer_tree_global,itimer_tree_revise
 use kdtree,       only:revtree
 use dim,          only:ind_timesteps,periodic,maxpsph

 integer, intent(inout) :: np
 integer, intent(in)    :: nactive  ! number of active particles on this task
 real,    intent(inout) :: xyzh(:,:)
 logical, intent(in)    :: apr_tree
 logical, intent(out)   :: revised
 integer, intent(in),    optional :: nptmass
 real,    intent(inout), optional :: xyzmh_ptmass(:,:)
 real,    allocatable :: wt(:)
 real    :: ncount(2)  ! active and live particles on all tasks
 logical :: rebalance
 integer, allocatable :: isinkowner(:)
 integer      :: i,ncross
 logical      :: sinktree
 real(kind=4) :: t1,t2,tcpu1,tcpu2,tg1,tg2,tcpug1,tcpug2,trev,tcpurev

 call get_timings(tg1,tcpug1)
 sinktree = .false.
 if (present(nptmass).and.present(xyzmh_ptmass)) sinktree=.true.

 ! particles that crossed a periodic boundary are wrapped first, so that their key is
 ! that of where they are (construct_root_node would wrap them after the move)
 if (periodic) then
    ncross = 0
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(np,xyzh,isperiodic) &
    !$omp private(i) &
    !$omp reduction(+:ncross)
    do i = 1,np
       if (.not.isdead_or_accreted(xyzh(4,i))) call cross_boundary(isperiodic,xyzh(:,i),ncross)
    enddo
    !$omp end parallel do
 endif

 ! the owner of each particle, from new sections only on full steps (all particles
 ! active), weighting particles by how often they are active
 ncount = reduceall_mpi('+',real([nactive,count(xyzh(4,1:np) > 0.)]))
 rebalance = (ncount(1) >= ncount(2))

 ! with few active particles, the local trees are revised for the new positions, keeping
 ! the owners, and the top of the global tree is rebuilt from their section roots
 revised = .false.
 trev = 0.
 tcpurev = 0.
 if (.not.apr_tree) then
    call get_timings(t1,tcpu1)
    call revtree(node,xyzh,np,leaf_is_active,ncells,revised,nactivetot=ncount(1),nlivetot=ncount(2))
    call get_timings(t2,tcpu2)
    trev = t2-t1
    tcpurev = tcpu2-tcpu1
    call increment_timer(itimer_tree_revise,trev,tcpurev)
 endif
 if (revised) then
    call build_section_top()
    call set_exact_top_sizes()
    call get_timings(tg2,tcpug2)
    ! the global time excludes the revision, which has its own timer
    call increment_timer(itimer_tree_global,(tg2-tg1)-trev,(tcpug2-tcpug1)-tcpurev)
    return
 endif
 allocate(wt(np))
 if (ind_timesteps) then
    wt = 2.**ibin(1:np)
 else
    wt = 1.
 endif

 ! the owner of each particle, and of each sink
 if (sinktree) then
    allocate(isinkowner(max(nptmass,1)))
    call assign_section_owners(np,xyzh,wt,rebalance,ibelong,nptmass,xyzmh_ptmass,isinkowner)
    ibelong(maxpsph+1:maxpsph+nptmass) = isinkowner(1:nptmass)
 else
    call assign_section_owners(np,xyzh,wt,rebalance,ibelong)
 endif

 deallocate(wt)

 ! move particles to where they belong
 call get_timings(t1,tcpu1)
 call balancedomains(np)
 call get_timings(t2,tcpu2)
 call increment_timer(itimer_balance,t2-t1,tcpu2-tcpu1)

 ! local tree
 if (sinktree) then
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,nptmass=nptmass,xyzmh_ptmass=xyzmh_ptmass,&
                  by_section=.true.)
 else
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,by_section=.true.)
 endif

 call get_timings(t1,tcpu1)
 call build_section_top()
 call set_exact_top_sizes()
 call get_timings(tg2,tcpug2)

 ! the global time excludes the local tree and any rejected revision, which have their own timers
 call increment_timer(itimer_tree_global,(tg2-tg1)-(t1-t2)-trev,(tcpug2-tcpug1)-(tcpu1-tcpu2)-tcpurev)

end subroutine maketreeglobal

!--------------------------------------------------------------------------------
!+
!  the top of the global tree from the sections: the cells of the tree of
!  sections, numbered top down, with the sections that hold particles as
!  leaves (cells with one such child are skipped). Each leaf is the root of
!  its section on its owner, gathered from all tasks, and the cells above
!  are built from their children in the same order on every task, so that
!  all tasks hold identical top nodes. The leaves, their owners and the
!  owners' local index of their roots are left in top_leaf, top_owner and
!  top_root
!+
!--------------------------------------------------------------------------------
subroutine build_section_top()
 use io,          only:nprocs
 use kdtree,      only:secroot,nsecroot,combine_child_nodes
 use mpisections, only:get_section_tree
 use mpitree,     only:allgatherv_nodes
 use mpiutils,    only:allgather_mpi,allgatherv_mpi,get_displs_mpi
 integer,      allocatable :: left(:),isec(:),owner(:),idrec(:),idall(:),rep(:),gidx(:)
 integer,      allocatable :: rootof(:),irec(:),istack(:)
 type(kdnode), allocatable :: block(:),blockall(:)
 integer :: ncell,nsec,kfirst,nowned,nrec,ntot,j,k,c,cl,cr,g,nstack,counts(nprocs),displs(nprocs)

 ! the roots of this task's sections with particles, with the index of each section
 call get_section_tree(ncell,left,isec,nsec,owner,kfirst,nowned)
 nrec = count(secroot(1:min(nowned,nsecroot)) > 0)
 allocate(block(max(nrec,1)),idrec(2*max(nrec,1)))

 ! copy each root, with the index of its section and its local index
 nrec = 0
 do j = 1,min(nowned,nsecroot)
    if (secroot(j) <= 0) cycle
    nrec = nrec + 1
    block(nrec) = node(secroot(j))
    idrec(2*nrec-1:2*nrec) = [kfirst+j-1,secroot(j)]
 enddo

 ! gather the roots of the sections of all tasks, with their indices
 call allgather_mpi(nrec,counts)
 displs = get_displs_mpi(counts)
 ntot = sum(counts)

 allocate(blockall(ntot),idall(2*ntot))
 call allgatherv_nodes(block(1:nrec),blockall,counts,displs)
 call allgatherv_mpi(idrec(1:2*nrec),idall,2*counts,2*displs)

 ! for each section with particles, its record and the local index of its root
 allocate(irec(nsec),rootof(nsec))
 irec = 0
 do j = 1,ntot
    irec(idall(2*j-1))   = j
    rootof(idall(2*j-1)) = idall(2*j)
 enddo

 ! the cell standing for each cell: itself, its one child holding particles, or none
 allocate(rep(ncell))
 do c = ncell,1,-1
    if (left(c) == 0) then
       rep(c) = 0
       if (irec(isec(c)) > 0) rep(c) = c
    else
       cl = rep(left(c))
       cr = rep(left(c)+1)
       rep(c) = c
       if (cl == 0) rep(c) = cr
       if (cr == 0) rep(c) = cl
    endif
 enddo

 ! number the cells that stand for themselves, top down (children after their parents)
 allocate(gidx(ncell),istack(2*ncell+2))
 gidx = 0
 ntop = 0
 nstack = 1
 istack(1) = rep(1)
 do while (nstack > 0)
    c = istack(nstack)
    nstack = nstack - 1
    if (c == 0) cycle
    ntop = ntop + 1
    gidx(c) = ntop
    if (left(c) > 0) then
       istack(nstack+1) = rep(left(c)+1)
       istack(nstack+2) = rep(left(c))
       nstack = nstack + 2
    endif
 enddo

 ! the leaves from the roots of the sections, the cells above from their children
 call allocate_top(ntop,ntot)
 nsec_top = 0
 do c = ncell,1,-1
    g = gidx(c)
    if (g == 0) cycle
    if (left(c) == 0) then
       k = isec(c)
       nodeglobal(g) = blockall(irec(k))
       nodeglobal(g)%leftchild  = 0
       nodeglobal(g)%rightchild = 0
       nsec_top = nsec_top + 1
       top_leaf(nsec_top)  = g
       top_owner(nsec_top) = owner(k)
       top_root(nsec_top)  = rootof(k)
    else
       cl = gidx(rep(left(c)))
       cr = gidx(rep(left(c)+1))
       call combine_child_nodes(nodeglobal(cl),nodeglobal(cr),nodeglobal(g))
       nodeglobal(g)%leftchild  = cl
       nodeglobal(g)%rightchild = cr
       nodeglobal(cl)%parent = g
       nodeglobal(cr)%parent = g
    endif
    nodeglobal(g)%level = 0
 enddo

 nodeglobal(1)%parent = 0

end subroutine build_section_top

!--------------------------------------------------------------------------------
!+
!  the exact size of the cells of the global tree above the sections: the
!  largest distance from each cell's centre of mass to the particles below
!  it, over all tasks. Combining the children only bounds it, which opens
!  these cells more often than needed. Each task measures the particles of
!  its own sections; every task then holds the same sizes
!+
!--------------------------------------------------------------------------------
subroutine set_exact_top_sizes()
 use io,       only:id
 use mpiutils, only:reduceall_mpi
 use kdtree,   only:inoderange
 use part,     only:treecache
 integer, parameter :: maxanc = 64
 real,    allocatable :: r2max(:)
 real    :: r2anc(maxanc),dx(3)
 integer :: ianc(maxanc),nanc,k,g,n,ip,j

 ! the largest squared distance from each cell to the particles of this task's sections below it
 allocate(r2max(ntop))
 r2max = 0.
 do k = 1,nsec_top
    if (top_owner(k) /= id) cycle

    ! the cells above this section
    nanc = 0
    g = nodeglobal(top_leaf(k))%parent
    do while (g > 0 .and. nanc < maxanc)
       nanc = nanc + 1
       ianc(nanc) = g
       g = nodeglobal(g)%parent
    enddo
    if (nanc == 0) cycle

    ! the farthest particle of the section from each of their centres
    n = top_root(k)
    r2anc = 0.
    !$omp parallel do default(none) schedule(static) &
    !$omp shared(n,inoderange,treecache,nodeglobal,ianc,nanc) &
    !$omp private(ip,j,dx) &
    !$omp reduction(max:r2anc)
    do ip = inoderange(1,n),inoderange(2,n)
       do j = 1,nanc
          dx = treecache(1:3,ip) - nodeglobal(ianc(j))%xcen
          r2anc(j) = max(r2anc(j),dot_product(dx,dx))
       enddo
    enddo
    !$omp end parallel do

    ! the largest over this task's sections
    do j = 1,nanc
       r2max(ianc(j)) = max(r2max(ianc(j)),r2anc(j))
    enddo
 enddo

 ! the largest over all tasks gives the size of each cell above the sections
 r2max = reduceall_mpi('max',r2max)
 do g = 1,ntop
    if (nodeglobal(g)%leftchild > 0) nodeglobal(g)%size = sqrt(r2max(g)) + epsilon(r2max)
 enddo

end subroutine set_exact_top_sizes

!-----------------------------------------------------------------------
!+
!  copy the top levels of the local tree of every task into the global
!  tree, below the node of its domain
!+
!-----------------------------------------------------------------------
subroutine graft_local_trees(nactive)
 use io,       only:nprocs,id,fatal
 use mpiutils, only:allgather_mpi,get_displs_mpi
 use mpitree,  only:allgatherv_nodes
 use mpiutils, only:allgatherv_mpi
 integer, intent(in) :: nactive  ! number of active particles on this task
 type(kdnode), allocatable :: block(:)
 integer,      allocatable :: ilocal(:)
 integer :: ntot,ifirst,itask,ibase,p,g,k,ileaf

 ! forget the nodes copied from the previous tree
 if (nown_roots > 0) nodemap(own_roots(1:nown_roots)) = 0
 if (nmirror > 0) nodemap(mirror_list(1:nmirror)) = 0
 call list_own_roots()

 ! the roots of this task's sections, as the (negative) top leaf that their copied children
 ! are attached to (pack_mirrored_nodes)
 do k = 1,nsec_top
    if (top_owner(k) == id) nodemap(top_root(k)) = -top_leaf(k)
 enddo

 call select_mirrored_nodes(nactive)
 allocate(block(nmirror))
 call pack_mirrored_nodes(block)

 if (.not.allocated(nmirror_task)) allocate(nmirror_task(nprocs),mirror_displ(nprocs))
 call allgather_mpi(nmirror,nmirror_task)
 mirror_displ = get_displs_mpi(nmirror_task)
 ntot = sum(nmirror_task)
 ifirst = ntop + 1     ! first node after the top nodes
 if (ifirst + ntot - 1 > size(nodeglobal)) call fatal('graft_local_trees','global tree exceeds its node budget')
 call allgatherv_nodes(block,nodeglobal(ifirst:),nmirror_task,mirror_displ)

 ! the local index of each copied node on its owner, negative for a leaf of its local tree
 allocate(ilocal(nmirror))
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nmirror,ilocal,mirror_list) &
 !$omp private(p)
 do p = 1,nmirror
    ilocal(p) = local_index_and_leaf(mirror_list(p))
 enddo
 !$omp end parallel do

 call allgatherv_mpi(ilocal,global_local_index(ifirst:),nmirror_task,mirror_displ)

 ! the top nodes: the leaves (sections) get their children from the blocks below
 cellatid = 0
 global_owner(1:ntop) = -1
 do k = 1,nsec_top
    ileaf = top_leaf(k)
    global_owner(ileaf) = top_owner(k)
    global_local_index(ileaf) = top_root(k)
    nodeglobal(ileaf)%leftchild  = 0
    nodeglobal(ileaf)%rightchild = 0
 enddo

 ! the blocks are linked within themselves: shift them to their place in the global tree,
 ! with the nodes copied from the root of a section attached to that section's leaf
 !$omp parallel do default(none) schedule(dynamic,1) &
 !$omp shared(nprocs,ifirst,mirror_displ,nmirror_task,nodeglobal,cellatid,global_owner) &
 !$omp private(itask,ibase,p,g,ileaf)
 do itask = 1,nprocs
    ibase = ifirst - 1 + mirror_displ(itask)
    do p = 1,nmirror_task(itask)
       g = ibase + p
       if (nodeglobal(g)%leftchild > 0)  nodeglobal(g)%leftchild  = nodeglobal(g)%leftchild  + ibase
       if (nodeglobal(g)%rightchild > 0) nodeglobal(g)%rightchild = nodeglobal(g)%rightchild + ibase
       if (nodeglobal(g)%parent > 0) then
          nodeglobal(g)%parent = nodeglobal(g)%parent + ibase
       else
          ! the left child first: the children of a section root come in that order
          ileaf = -nodeglobal(g)%parent
          nodeglobal(g)%parent = ileaf
          if (nodeglobal(ileaf)%leftchild == 0) then
             nodeglobal(ileaf)%leftchild = g
          else
             nodeglobal(ileaf)%rightchild = g
          endif
       endif
       if (nodeglobal(g)%leftchild == 0 .and. nodeglobal(g)%rightchild == 0) cellatid(g) = itask
    enddo
    global_owner(ibase+1:ibase+nmirror_task(itask)) = itask - 1
 enddo
 !$omp end parallel do

 ! a section with no copied children is a leaf of the global tree, and of its owner's tree
 do k = 1,nsec_top
    ileaf = top_leaf(k)
    if (nodeglobal(ileaf)%leftchild == 0) then
       cellatid(ileaf) = top_owner(k) + 1
       global_local_index(ileaf) = -top_root(k)
    endif
 enddo

 ! where the copied local nodes are in the global tree; the roots of this task's sections
 ! are their leaves
 ibase = ifirst - 1 + mirror_displ(id+1)
 do k = 1,nsec_top
    if (top_owner(k) == id) nodemap(top_root(k)) = top_leaf(k)
 enddo

 ! and the nodes this task copied into its block
 do p = 1,nmirror
    nodemap(mirror_list(p)) = ibase + p
 enddo

end subroutine graft_local_trees

!-----------------------------------------------------------------------
!+
!  the most top nodes of the global tree: two per section
!+
!-----------------------------------------------------------------------
integer function ntop_max()
 use io, only:nprocs

 ntop_max = 2*maxsec_per_task*nprocs

end function ntop_max

!-----------------------------------------------------------------------
!+
!  room for ntop top nodes and nsec sections
!+
!-----------------------------------------------------------------------
subroutine allocate_top(ntopnew,nsec)
 use io, only:nprocs,fatal
 integer, intent(in) :: ntopnew,nsec

 if (ntopnew > ntop_max()) call fatal('allocate_top','too many sections for the global tree')
 ntop = ntopnew
 nsec_top = nsec
 if (allocated(top_leaf)) deallocate(top_leaf,top_owner,top_root,own_lo,own_hi,section_of,sec_start,sec_by_owner)
 allocate(top_leaf(nsec),top_owner(nsec),top_root(nsec))
 allocate(own_lo(ntop),own_hi(ntop),section_of(ntop),sec_start(nprocs+1),sec_by_owner(nsec))

end subroutine allocate_top

!-----------------------------------------------------------------------
!+
!  from the sections and their owners: the section of each top node, the
!  range of tasks owning the sections below each top node (the children
!  of a top node come after it), and the sections of each task
!+
!-----------------------------------------------------------------------
subroutine set_top_owners()
 use io, only:nprocs
 integer :: g,k,il,ir,r
 integer :: next(nprocs)

 ! each section's leaf, and its owner as the range of tasks below it
 section_of = 0
 own_lo = huge(1)
 own_hi = -1
 do k = 1,nsec_top
    g = top_leaf(k)
    section_of(g) = k
    own_lo(g) = top_owner(k)
    own_hi(g) = top_owner(k)
 enddo

 ! the range of each cell above the sections from its children, bottom up
 do g = ntop,1,-1
    if (section_of(g) > 0) cycle
    il = nodeglobal(g)%leftchild
    ir = nodeglobal(g)%rightchild
    if (il > 0) then
       own_lo(g) = min(own_lo(g),own_lo(il))
       own_hi(g) = max(own_hi(g),own_hi(il))
    endif
    if (ir > 0) then
       own_lo(g) = min(own_lo(g),own_lo(ir))
       own_hi(g) = max(own_hi(g),own_hi(ir))
    endif
 enddo

 ! the sections of each task, in the order of their indices
 sec_start = 0
 do k = 1,nsec_top
    sec_start(top_owner(k)+2) = sec_start(top_owner(k)+2) + 1
 enddo

 ! the first section of each task, from the counts
 sec_start(1) = 1
 do r = 2,nprocs+1
    sec_start(r) = sec_start(r) + sec_start(r-1)
 enddo

 ! the sections of each task in that order
 next = sec_start(1:nprocs)
 do k = 1,nsec_top
    r = top_owner(k) + 1
    sec_by_owner(next(r)) = k
    next(r) = next(r) + 1
 enddo

end subroutine set_top_owners

!-----------------------------------------------------------------------
!+
!  the roots of the sections owned by this task
!+
!-----------------------------------------------------------------------
subroutine list_own_roots()
 use io, only:id
 integer :: k

 if (allocated(own_roots)) deallocate(own_roots)

 ! the local roots of the sections whose owner is this task
 allocate(own_roots(max(nsec_top,1)))
 nown_roots = 0
 do k = 1,nsec_top
    if (top_owner(k) == id) then
       nown_roots = nown_roots + 1
       own_roots(nown_roots) = top_root(k)
    endif
 enddo

end subroutine list_own_roots

!-----------------------------------------------------------------------
!+
!  build the combined tree: append the global tree after the local tree,
!  except for this task's sections and block, which are the roots of its
!  sections and the local nodes they were copied from. The root of each section of
!  this task takes the data of its leaf, which the other tasks use, so that
!  every task sees exactly the same node there
!+
!-----------------------------------------------------------------------
subroutine combine_trees()
 use io, only:id,fatal
 type(kdnode) :: rootlinks
 integer :: g,n,ng,k

 ng = ntop + sum(nmirror_task)
 if (int(ncells) + ng > size(node)) call fatal('combine_trees','no room to append the global tree')

 ! forget the nodes of other tasks appended to the previous tree, which may
 ! have started at a different place if the local tree has changed size
 if (ncomb > ncells_appended) remote_leaf_owner(ncells_appended+1:ncomb) = 0
 ncells_appended = int(ncells)

 ! the root of each of this task's sections takes the data of its leaf, links kept
 do k = 1,nsec_top
    if (top_owner(k) /= id) cycle
    n = top_root(k)
    rootlinks = node(n)
    node(n) = nodeglobal(top_leaf(k))
    node(n)%leftchild  = rootlinks%leftchild
    node(n)%rightchild = rootlinks%rightchild
    node(n)%level      = rootlinks%level
    node(n)%parent     = combined_index(nodeglobal(top_leaf(k))%parent)
 enddo

 ! every other node of the global tree is appended after the local tree, with its links
 ! moved to the combined tree
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(ng,ncells,node,nodeglobal,cellatid,global_local_index) &
 !$omp shared(remote_leaf_owner,remote_local_index,leaf_is_active) &
 !$omp private(g,n)
 do g = 1,ng
    if (is_this_task(g)) cycle
    n = int(ncells) + g
    node(n) = nodeglobal(g)
    node(n)%leftchild  = combined_index(nodeglobal(g)%leftchild)
    node(n)%rightchild = combined_index(nodeglobal(g)%rightchild)
    node(n)%parent     = combined_index(nodeglobal(g)%parent)
    remote_leaf_owner(n) = cellatid(g)
    remote_local_index(n) = global_local_index(g)
    leaf_is_active(n)    = min(cellatid(g),1)
 enddo
 !$omp end parallel do

 irootcomb = combined_index(1)
 ncomb = int(ncells) + ng
 ncomb_tree = ncomb
 nlet_parent = 0
 nlet_local  = 0

 ! no profiles received yet: every node may have active leaves below it
 profile_len(1:ng)   = -1
 global_active(1:ng) = .true.
 nprof_set     = 0
 ninactive_set = 0

end subroutine combine_trees

!-----------------------------------------------------------------------
!+
!  forget the nodes imported by the last call to let_exchange: the global
!  leaves they were attached to become leaves again, as combine_trees left
!  them. The rest of the combined tree is unchanged since the tree was built
!+
!-----------------------------------------------------------------------
subroutine forget_let_nodes()
 integer :: k,n,g

 ! the global leaves given imported children: restore the links and leaf status they had
 do k = 1,nlet_parent
    n = let_parent(k)
    g = n - int(ncells)
    node(n)%leftchild  = combined_index(nodeglobal(g)%leftchild)
    node(n)%rightchild = combined_index(nodeglobal(g)%rightchild)
    remote_leaf_owner(n) = cellatid(g)
    leaf_is_active(n)    = min(cellatid(g),1)
 enddo

 ! the global leaves whose particles were imported: leaves of their owners again
 do k = 1,nlet_local
    n = let_local(k)
    remote_leaf_owner(n) = cellatid(n - int(ncells))
 enddo
 nlet_local = 0

 ! and drop the imported nodes
 nlet_parent = 0
 if (ncomb > ncomb_tree) remote_leaf_owner(ncomb_tree+1:ncomb) = 0
 ncomb = ncomb_tree

end subroutine forget_let_nodes

!-----------------------------------------------------------------------
!+
!  whether a node of the global tree is a section of this task or one
!  of the local nodes it copied into the global tree
!+
!-----------------------------------------------------------------------
logical function is_this_task(g)
 use io, only:id
 integer, intent(in) :: g
 integer :: ibase

 ibase = ntop + mirror_displ(id+1)
 if (g <= ntop) then
    is_this_task = (section_of(g) > 0)
    if (is_this_task) is_this_task = (top_owner(section_of(g)) == id)
 else
    is_this_task = (g > ibase .and. g <= ibase + nmirror)
 endif

end function is_this_task

!-----------------------------------------------------------------------
!+
!  index in the combined tree of a node of the global tree (0 for none)
!+
!-----------------------------------------------------------------------
integer function combined_index(g)
 use io, only:id
 integer, intent(in) :: g
 integer :: ibase,k

 ibase = ntop + mirror_displ(id+1)
 combined_index = int(ncells) + g
 if (g <= 0) then
    combined_index = 0
 elseif (g <= ntop) then
    k = section_of(g)
    if (k > 0) then
       if (top_owner(k) == id) combined_index = top_root(k)
    endif
 elseif (g > ibase .and. g <= ibase + nmirror) then
    combined_index = mirror_list(g - ibase)
 endif

end function combined_index

!-----------------------------------------------------------------------
!+
!  local index of a node, negative if it is a leaf of the local tree
!+
!-----------------------------------------------------------------------
integer function local_index_and_leaf(n)
 integer, intent(in) :: n

 local_index_and_leaf = n
 if (node(n)%leftchild <= 0 .and. node(n)%rightchild <= 0) local_index_and_leaf = -n

end function local_index_and_leaf

!-----------------------------------------------------------------------
!+
!  exchange the locally essential tree: each task sends every other task
!  the local nodes below its global leaves that the other task's walks
!  could open, and attaches the nodes it receives below the global leaves
!  of their owners in the combined tree. Called after sync_hmax_mpi, so
!  that the nodes carry the smoothing lengths of the density calculation
!+
!-----------------------------------------------------------------------
subroutine let_exchange()
 use io,       only:nprocs,id
 use mpiutils, only:alltoall_mpi,alltoallv_mpi
 use dim,      only:maximport_leaves,ncellsmax
 use omputils, only:omp_thread_num
 type :: realbuf
    real, allocatable :: x(:)
 end type realbuf
 type(realbuf)        :: records(nprocs)
 real,    allocatable :: sendbuf(:),recvbuf(:)
 integer, allocatable :: list(:),iparent(:),iside(:),implist(:,:,:)
 integer, allocatable :: leaf(:),loc(:),idest(:),ioff(:),itemr(:,:,:),nitemr(:),items(:,:)
 integer :: nsend(nprocs),nrecv(nprocs),nimp(nprocs),nnode(nprocs),nleafreal(nprocs)
 integer :: leafseg(2,nprocs),irecv,nlist,nbuf,k,iproc,nleaves,nitems,ith
 logical :: reset

 call forget_let_nodes()
 call exchange_size_profiles()

 ! select the nodes and leaves to send to each task. The descent for each task is split
 ! into pairs near the top of the tree, shared by the threads; the nodes found are then
 ! listed for each task, one task per thread at a time
 call set_reach_ref()
 call start_let_marks(reset)
 allocate(implist(2,maximport_leaves,nprocs),itemr(3,nlet_split+3,nprocs),nitemr(nprocs))
 nnode  = 0
 nimp   = 0
 nitemr = 0
 !$omp parallel default(none) &
 !$omp shared(nprocs,id,ncomb,ncellsmax,records,implist,nimp,nnode,reset) &
 !$omp shared(itemr,nitemr,items,nitems,let_opened,let_round) &
 !$omp private(iproc,list,iparent,iside,nlist,k,ith)

 ! each thread's marks of the nodes opened, cleared when the marks start a new round
 if (allocated(let_mark)) then
    if (size(let_mark) < ncellsmax+1) deallocate(let_mark)
 endif
 if (.not.allocated(let_mark)) then
    allocate(let_mark(ncellsmax+1))
    let_mark = 0
 endif
 if (reset) let_mark = 0
 ith = omp_thread_num() + 1

 ! the first steps of the descent for each task, breadth first, until it is split into pairs
 !$omp do schedule(dynamic,1)
 do iproc = 0,nprocs-1
    if (iproc == id) cycle
    call split_let_descent(iproc,let_opened(ith,iproc+1),itemr(:,:,iproc+1),nitemr(iproc+1))
 enddo
 !$omp enddo

 ! the pairs of all tasks in one list
 !$omp single
 nitems = sum(nitemr)
 allocate(items(4,max(nitems,1)))
 nitems = 0
 do iproc = 0,nprocs-1
    do k = 1,nitemr(iproc+1)
       nitems = nitems + 1
       items(:,nitems) = [iproc,itemr(:,k,iproc+1)]
    enddo
 enddo
 !$omp end single

 ! the descent below each pair, shared by the threads
 !$omp do schedule(dynamic,1)
 do k = 1,nitems
    iproc = items(1,k)
    call descend_let_pair(iproc,items(2:4,k),let_opened(ith,iproc+1))
 enddo
 !$omp enddo

 ! the nodes and leaves to send to each task, and the records of the nodes
 allocate(list(ncomb),iparent(ncomb),iside(ncomb))
 !$omp do schedule(dynamic,1)
 do iproc = 0,nprocs-1
    if (iproc == id) cycle
    call select_let_nodes(let_opened(:,iproc+1),list,iparent,iside,nlist)
    call select_import_leaves(list,nlist,implist(:,:,iproc+1),nimp(iproc+1))
    allocate(records(iproc+1)%x(nlist*nlet_record))

    ! the record of each node
    do k = 1,nlist
       call put_let_record(records(iproc+1)%x((k-1)*nlet_record+1:k*nlet_record),list(k),iparent(k),iside(k))
    enddo
    nnode(iproc+1) = nlist*nlet_record
 enddo
 !$omp enddo
 deallocate(list,iparent,iside)
 !$omp end parallel

 ! one message to each task: the number of reals of node records, the node
 ! records, then the particles of the leaves of large reach selected for it
 call leaf_layout(implist,nimp,leaf,loc,idest,ioff,nleafreal,nleaves)
 nsend = 0
 do iproc = 0,nprocs-1
    if (nnode(iproc+1) + nleafreal(iproc+1) > 0) nsend(iproc+1) = 1 + nnode(iproc+1) + nleafreal(iproc+1)
 enddo

 ! the node records of each task in the buffer, leaving room for its leaves
 allocate(sendbuf(max(sum(nsend),1)))
 nbuf = 0
 do iproc = 0,nprocs-1
    if (nsend(iproc+1) == 0) cycle
    sendbuf(nbuf+1) = nnode(iproc+1)
    if (nnode(iproc+1) > 0) sendbuf(nbuf+2:nbuf+1+nnode(iproc+1)) = records(iproc+1)%x
    leafseg(1,iproc+1) = nbuf + 1 + nnode(iproc+1)  ! where its leaves start
    nbuf = nbuf + nsend(iproc+1)
 enddo

 ! then the particles of the leaves
 call pack_leaf_particles(nleaves,leaf,loc,idest,ioff,leafseg(1,:),sendbuf)
 nlet_sent = nlet_sent + sum(nnode)/nlet_record

 ! exchange the messages
 call alltoall_mpi(nsend,nrecv)

 allocate(recvbuf(max(sum(nrecv),1)))
 call alltoallv_mpi(sendbuf,nsend,recvbuf,nrecv)

 ! attach the nodes from each task below its global leaves, and find where its leaves are
 if (.not.allocated(let_base)) allocate(let_base(nprocs))
 irecv = 0
 do iproc = 0,nprocs-1
    let_base(iproc+1) = ncomb
    leafseg(:,iproc+1) = irecv
    if (nrecv(iproc+1) > 0) then
       nnode(iproc+1) = nint(recvbuf(irecv+1))
       if (nnode(iproc+1) > 0) call attach_let_nodes(iproc,recvbuf(irecv+2:irecv+1+nnode(iproc+1)))
       leafseg(:,iproc+1) = [irecv + 1 + nnode(iproc+1),irecv + nrecv(iproc+1)]
       nlet_recv = nlet_recv + nnode(iproc+1)/nlet_record
    endif
    irecv = irecv + nrecv(iproc+1)
 enddo

 ! the particles of the leaves of large reach, once the nodes are attached
 call unpack_leaf_particles(recvbuf,leafseg)

end subroutine let_exchange

!-----------------------------------------------------------------------
!+
!  the leaves of this task whose particles are sent to task iproc with
!  the nodes in list. The candidates are the true leaves iproc will have
!  in its combined tree (the leaves in list and this task's leaves in the
!  global tree) that the walks of iproc's active cells could pair with
!  (marked by let_pair_step). Of these, those with a reach (size +
!  radkern*hmax) more than reach_import_factor times the typical reach of
!  the active leaves of all tasks (set_reach_ref) are sent, largest first:
!  large leaves in sparse regions, which many of iproc's cells meet. The
!  reference is neither the median of the candidates (a distant receiver
!  may meet only a few leaves, all of them large) nor of this task's own
!  leaves (a task in a sparse region may hold only large leaves). Only
!  leaves whose particles are all active are sent: the force calculation
!  can wake inactive neighbours, which would then have to be returned.
!  Each leaf is listed in implist with its locator on iproc: its position
!  in list, or minus its index in the global tree
!+
!-----------------------------------------------------------------------
subroutine select_import_leaves(list,nlist,implist,nleaves)
 use dim,      only:maximport_leaves,maximport_part,gr,maxpsph,reach_import_factor
 use kernel,   only:radkern
 use kdtree,   only:inoderange,inodeparts
 use sortutils,only:indexx
 integer, intent(in)  :: list(:),nlist
 integer, intent(out) :: implist(:,:)
 integer, intent(out) :: nleaves
 integer, allocatable :: cand(:),loc(:),indx(:)
 real,    allocatable :: reach(:)
 real    :: reachmin
 integer :: k,n,nc,np,npsent

 nleaves = 0
 if (gr) return   ! the metric of imported particles is not set up

 ! the candidates among the nodes to send
 allocate(cand(nlist+nown_gleaves),loc(nlist+nown_gleaves))
 nc = 0
 do k = 1,nlist
    if (.not.is_import_candidate(list(k))) cycle
    nc = nc + 1
    cand(nc) = list(k)
    loc(nc)  = k
 enddo

 ! and among this task's leaves in the global tree
 do k = 1,nown_gleaves
    n = own_gleaves(k)
    if (.not.is_import_candidate(n)) cycle
    nc = nc + 1
    cand(nc) = n
    loc(nc)  = -nodemap(n)
 enddo

 if (nc == 0) return

 ! their reach, sorted
 allocate(reach(nc),indx(nc))
 reach = node(cand(1:nc))%size + radkern*node(cand(1:nc))%hmax
 call indexx(nc,reach,indx)
 reachmin = reach_import_factor*reach_ref

 ! largest reach first, within the budgets of leaves and particles
 npsent = 0
 do k = nc,1,-1
    if (nleaves >= maximport_leaves .or. reach(indx(k)) <= reachmin) exit
    n  = cand(indx(k))
    np = inoderange(2,n) - inoderange(1,n) + 1
    if (npsent + np > maximport_part) cycle

    ! all particles active (inodeparts > 0), and none of them a sink in the tree
    if (any(inodeparts(inoderange(1,n):inoderange(2,n)) <= 0) .or. &
        any(inodeparts(inoderange(1,n):inoderange(2,n)) > maxpsph)) cycle
    nleaves = nleaves + 1
    implist(:,nleaves) = [n,loc(indx(k))]
    npsent = npsent + np
 enddo

end subroutine select_import_leaves

!-----------------------------------------------------------------------
!+
!  the typical reach (size + radkern*hmax) of the active leaves, the
!  reference for the reach of the leaves sent with their particles: the
!  median over the tasks of the median of each task's active leaves
!+
!-----------------------------------------------------------------------
subroutine set_reach_ref()
 use dim,       only:reach_import_factor
 use io,        only:id,nprocs
 use kernel,    only:radkern
 use mpiutils,  only:reduceall_mpi
 use sortutils, only:indexx
 real,    allocatable :: reach(:)
 integer, allocatable :: indx(:)
 real    :: medtask(nprocs)
 integer :: n,nt

 ! the median reach of this task's active leaves
 n = nactive_leaves
 medtask = 0.
 if (n > 0) then
    allocate(reach(n),indx(n))
    reach = node(active_leaves(1:n))%size + radkern*node(active_leaves(1:n))%hmax
    call indexx(n,reach,indx)
    medtask(id+1) = reach(indx(n/2+1))
    deallocate(reach,indx)
 endif

 ! the median over the tasks with active leaves, the same on every task
 medtask = reduceall_mpi('+',medtask)
 nt = count(medtask > 0.)
 ! no leaf is sent if there are no active leaves
 reach_ref = huge(reach_ref)/reach_import_factor
 if (nt == 0) return
 allocate(reach(nt),indx(nt))
 reach = pack(medtask,medtask > 0.)
 call indexx(nt,reach,indx)
 reach_ref = reach(indx(nt/2+1))

end subroutine set_reach_ref

!-----------------------------------------------------------------------
!+
!  whether local node n is a leaf with particles that the walks of the
!  current receiver could pair with
!+
!-----------------------------------------------------------------------
logical function is_import_candidate(n)
 integer, intent(in) :: n

 is_import_candidate = (local_index_and_leaf(n) < 0 .and. leaf_is_active(n) /= 0 .and. let_mark(n) == let_stamp)

end function is_import_candidate

!-----------------------------------------------------------------------
!+
!  the leaves of large reach selected for each task, in task order: each
!  is sent as a header (locator, number of particles), then its particles.
!  ioff is the position of each leaf in its task's (idest) part of the message
!+
!-----------------------------------------------------------------------
subroutine leaf_layout(implist,nimp,leaf,loc,idest,ioff,nleafreal,nleaves)
 use io,     only:nprocs
 use kdtree, only:inoderange
 integer,              intent(in)  :: implist(:,:,:),nimp(:)
 integer, allocatable, intent(out) :: leaf(:),loc(:),idest(:),ioff(:)
 integer,              intent(out) :: nleafreal(:),nleaves
 integer :: iproc,k,m,np,nrec

 nrec = nreals_per_particle()
 nleaves = sum(nimp)
 allocate(leaf(nleaves),loc(nleaves),idest(nleaves),ioff(nleaves))

 ! the leaves of each task in turn, each after the previous one in that task's part
 m = 0
 nleafreal = 0
 do iproc = 0,nprocs-1
    do k = 1,nimp(iproc+1)
       m = m + 1
       leaf(m)  = implist(1,k,iproc+1)
       loc(m)   = implist(2,k,iproc+1)
       idest(m) = iproc
       ioff(m)  = nleafreal(iproc+1)
       np = inoderange(2,leaf(m)) - inoderange(1,leaf(m)) + 1
       nleafreal(iproc+1) = nleafreal(iproc+1) + 2 + np*nrec
    enddo
 enddo

end subroutine leaf_layout

!-----------------------------------------------------------------------
!+
!  pack the particles of the leaves into the message buffer, where each
!  task's leaves start at istart (one leaf per thread at a time)
!+
!-----------------------------------------------------------------------
subroutine pack_leaf_particles(nleaves,leaf,loc,idest,ioff,istart,sendbuf)
 use part,   only:fill_sendbuf,ipartbufsize
 use kdtree, only:inoderange,inodeparts
 integer, intent(in)    :: nleaves,leaf(:),loc(:),idest(:),ioff(:),istart(:)
 real,    intent(inout) :: sendbuf(:)
 integer :: k,m,n,j,np
 real    :: xtemp(ipartbufsize)

 ! each leaf: its header, then its particles
 !$omp parallel do default(none) schedule(dynamic,1) &
 !$omp shared(nleaves,leaf,loc,idest,ioff,istart) &
 !$omp shared(sendbuf,inoderange,inodeparts) &
 !$omp private(m,n,j,k,xtemp,np)
 do m = 1,nleaves
    n = leaf(m)
    k = istart(idest(m)+1) + ioff(m)
    sendbuf(k+1) = loc(m)
    sendbuf(k+2) = inoderange(2,n) - inoderange(1,n) + 1
    k = k + 2

    ! the particles of the leaf
    do j = inoderange(1,n),inoderange(2,n)
       call fill_sendbuf(inodeparts(j),xtemp,np)   ! np: reals for this particle
       sendbuf(k+1:k+np) = xtemp(1:np)
       k = k + np
    enddo
 enddo
 !$omp end parallel do

end subroutine pack_leaf_particles

!-----------------------------------------------------------------------
!+
!  place the particles of the leaves received after npart in the particle
!  arrays: each such leaf of the combined tree becomes a local leaf, whose
!  particles the walk takes as neighbours. Leaves that do not fit stay
!  leaves of another task. The leaves from task iproc are in
!  recvbuf(leafseg(1,iproc+1)+1:leafseg(2,iproc+1)); each is a node
!  imported from it (positive locator) or one of its leaves in the global tree
!+
!-----------------------------------------------------------------------
subroutine unpack_leaf_particles(recvbuf,leafseg)
 use io,     only:nprocs
 use dim,    only:maxpsph
 use part,   only:npart,unfill_buffer,xyzh,get_pmass,maxphase,maxp,treecache
 use kdtree, only:inoderange,inodeparts
 real,    intent(in) :: recvbuf(:)
 integer, intent(in) :: leafseg(:,:)
 integer, allocatable :: newrange(:,:),leaf(:),ioff(:),ifirst(:)
 integer :: iproc,m,n,j,k,np,nrec,irecv,islot,islot0,nleaves
 logical :: use_gas

 use_gas = .not.(maxphase==maxp)  ! all particles are gas unless types are stored
 nrec = nreals_per_particle()

 ! the leaves of the combined tree need a range of particles in the tree arrays
 if (size(inoderange,2) < size(node)) then
    allocate(newrange(2,size(node)))
    newrange(:,1:size(inoderange,2)) = inoderange
    call move_alloc(newrange,inoderange)
 endif

 ! the leaves that fit, with the position of their particles and their first slot
 n = sum(leafseg(2,:) - leafseg(1,:))/(2 + nrec) + 1
 allocate(leaf(n),ioff(n),ifirst(n))
 nleaves = 0
 ! the particles go after the particle arrays and after the local tree list (which
 ! ends with this task's sinks in the tree)
 islot0 = max(npart,inoderange(2,irootnode))
 islot  = islot0
 do iproc = 0,nprocs-1
    irecv = leafseg(1,iproc+1)
    do while (irecv < leafseg(2,iproc+1))

       ! the node of the combined tree from the leaf's locator
       if (nint(recvbuf(irecv+1)) > 0) then
          n = let_base(iproc+1) + nint(recvbuf(irecv+1))
       else
          n = int(ncells) - nint(recvbuf(irecv+1))
       endif

       ! the leaf becomes local if its particles fit after npart
       np = nint(recvbuf(irecv+2))
       if (islot + np <= maxpsph .and. remote_leaf_owner(n) == iproc + 1) then
          nleaves = nleaves + 1
          leaf(nleaves)   = n
          ioff(nleaves)   = irecv + 2
          ifirst(nleaves) = islot + 1
          inoderange(:,n) = [islot+1,islot+np]
          remote_leaf_owner(n) = 0
          if (n <= ncomb_tree) call add_to_list(let_local,nlet_local,n)
          islot = islot + np
       endif
       irecv = irecv + 2 + np*nrec
    enddo
 enddo

 ! the particles of each leaf, with their tree cache entries
 !$omp parallel do default(none) schedule(dynamic,1) &
 !$omp shared(nleaves,leaf,ioff,ifirst,nrec,recvbuf,inoderange) &
 !$omp shared(inodeparts,treecache,xyzh,use_gas) &
 !$omp private(m,j,k)
 do m = 1,nleaves
    do j = 0,inoderange(2,leaf(m)) - inoderange(1,leaf(m))
       k = ifirst(m) + j
       call unfill_buffer(k,recvbuf(ioff(m)+j*nrec+1:ioff(m)+(j+1)*nrec))
       inodeparts(k)    = k
       treecache(1:4,k) = xyzh(1:4,k)
       treecache(5,k)   = get_pmass(k,use_gas)
    enddo
 enddo
 !$omp end parallel do

 nimp_leaves_recv = nimp_leaves_recv + nleaves
 nimp_part_recv   = nimp_part_recv + (islot - islot0)

end subroutine unpack_leaf_particles

!-----------------------------------------------------------------------
!+
!  the number of reals sent per particle, which depends only on the build
!+
!-----------------------------------------------------------------------
integer function nreals_per_particle()
 use part, only:fill_sendbuf,ipartbufsize
 real :: xtemp(ipartbufsize)

 call fill_sendbuf(1,xtemp,nreals_per_particle)

end function nreals_per_particle

!-----------------------------------------------------------------------
!+
!  start the marks of the descents of this call: room for the opened
!  nodes found by each thread for each task, and a new round of marks
!  (reset, with the marks of every thread, if the round would overflow)
!+
!-----------------------------------------------------------------------
subroutine start_let_marks(reset)
 use io, only:nprocs
!$ use omp_lib, only:omp_get_max_threads
 logical, intent(out) :: reset
 integer :: i,j,nthreads

 ! a list of opened nodes for each thread and task, emptied
 nthreads = 1
!$ nthreads = omp_get_max_threads()
 if (allocated(let_opened)) then
    if (size(let_opened,1) /= nthreads .or. size(let_opened,2) /= nprocs) deallocate(let_opened)
 endif
 if (.not.allocated(let_opened)) allocate(let_opened(nthreads,nprocs))

 ! empty every list
 do j = 1,nprocs
    do i = 1,nthreads
       let_opened(i,j)%n = 0
    enddo
 enddo

 ! a new round of marks, starting again when the values would overflow
 reset = (let_round > huge(let_round) - 2*nprocs)
 if (reset) let_round = 0
 let_round = let_round + nprocs

end subroutine start_let_marks

!-----------------------------------------------------------------------
!+
!  start the descent for task iproc: the walks on iproc could open the
!  nodes of this task that are met in pairs (node on iproc's side, node
!  on this task's side) of the combined tree descending together from the
!  root, as in the unsplit dual tree traversal. Once iproc's side reaches
!  one of its global leaves, whose nodes below are unknown here, a region
!  bound stands in for them and only this task's side descends. The pairs
!  are expanded breadth first until there are nlet_split of them, which
!  are returned as (node on iproc's side, node on this side, steps below
!  a global leaf of iproc) for the threads to share
!+
!-----------------------------------------------------------------------
subroutine split_let_descent(iproc,opened,item,nitem)
 integer,       intent(in)    :: iproc
 type(intlist), intent(inout) :: opened
 integer,       intent(out)   :: item(:,:),nitem
 integer, allocatable :: iqueue(:,:),istack(:,:)
 integer :: ihead,ntail,nstack,k

 ! expand the pairs from the root breadth first, one step each, until there are nlet_split
 allocate(iqueue(3,4*nlet_split),istack(3,4))
 ntail = 0
 call push_pair(iqueue,ntail,irootcomb,irootcomb,0)
 ihead = 1
 do while (ihead <= ntail .and. ntail - ihead + 1 < nlet_split)
    nstack = 1
    istack(:,1) = iqueue(:,ihead)
    ihead = ihead + 1
    call let_pair_step(istack,nstack,iproc,opened)

    ! the pairs below it go to the back of the queue
    do k = 1,nstack
       call push_pair(iqueue,ntail,istack(1,k),istack(2,k),istack(3,k))
    enddo
 enddo

 ! the pairs left in the queue are the items shared by the threads
 nitem = ntail - ihead + 1
 item(:,1:nitem) = iqueue(:,ihead:ntail)

end subroutine split_let_descent

!-----------------------------------------------------------------------
!+
!  the descent below one pair of the split, for task iproc
!+
!-----------------------------------------------------------------------
subroutine descend_let_pair(iproc,pair,opened)
 integer,       intent(in)    :: iproc,pair(3)
 type(intlist), intent(inout) :: opened
 integer, allocatable :: istack(:,:)
 integer :: nstack

 ! every step of the descent below the pair
 allocate(istack(3,1024))
 nstack = 1
 istack(:,1) = pair
 do while (nstack > 0)
    call let_pair_step(istack,nstack,iproc,opened)
 enddo

end subroutine descend_let_pair

!-----------------------------------------------------------------------
!+
!  the local nodes to send to a task, given the nodes its walks could
!  open(found by each thread): the children of the opened nodes below
!  this task's global leaves, in pre-order, with the parent of each
!  (negative: a global node, by its index in the global tree; positive:
!  a position in the list) and which child it is (1 left, 2 right)
!+
!-----------------------------------------------------------------------
subroutine select_let_nodes(opened,list,iparent,iside,nlist)
 type(intlist), intent(in)  :: opened(:)
 integer,       intent(out) :: list(:),iparent(:),iside(:),nlist
 integer, allocatable :: istack(:,:)
 integer :: nstack,ileaf,n,ichild(2),k,i,ipos

 ! mark the nodes opened for this task, as found by every thread
 let_stamp = let_stamp + 1
 do i = 1,size(opened)
    do k = 1,opened(i)%n
       let_mark(opened(i)%x(k)) = let_stamp
    enddo
 enddo

 ! list the children of the marked nodes in pre-order: a node marked below a global
 ! leaf is itself in the list, so its parent comes first
 allocate(istack(3,1024))
 nlist = 0

 ! start from this task's global leaves that are marked
 nstack = 0
 do k = 1,nown_gleaves
    ileaf = own_gleaves(k)
    if (let_mark(ileaf) == let_stamp) call push_pair(istack,nstack,ileaf,-nodemap(ileaf),0)
 enddo

 ! then list the children of each, descending into those that are marked
 do while (nstack > 0)
    n    = istack(1,nstack)
    ipos = istack(2,nstack)
    nstack = nstack - 1

    ! list each child, and descend into it if it is marked
    ichild = [node(n)%leftchild,node(n)%rightchild]
    do k = 1,2
       if (ichild(k) <= 0) cycle
       nlist = nlist + 1
       list(nlist)    = ichild(k)
       iparent(nlist) = ipos
       iside(nlist)   = k
       if (let_mark(ichild(k)) == let_stamp) call push_pair(istack,nstack,ichild(k),nlist,0)
    enddo
 enddo

end subroutine select_let_nodes

!-----------------------------------------------------------------------
!+
!  one step of the descent for task iproc: test the pair on top of
!  the stack, a node on iproc's side against a node on this task's side
!  (both in the combined tree), and push the pairs below it
!+
!-----------------------------------------------------------------------
subroutine let_pair_step(istack,nstack,iproc,opened)
 use io, only:id
 integer, allocatable, intent(inout) :: istack(:,:)
 integer,              intent(inout) :: nstack
 integer,              intent(in)    :: iproc
 type(intlist),        intent(inout) :: opened
 integer :: ia,ib,istep,ka(2),kb(2),i,j
 real    :: sizemax,hmaxmax
 logical :: may_open

 ia    = istack(1,nstack)
 ib    = istack(2,nstack)
 istep = istack(3,nstack)
 nstack = nstack - 1
 if (.not.(on_side(ia,iproc) .and. on_side(ib,id))) return
 ka = [node(ia)%leftchild,node(ia)%rightchild]
 kb = [node(ib)%leftchild,node(ib)%rightchild]

 if (remote_leaf_owner(ia) == iproc + 1) then
    ! a global leaf of iproc, istep steps of the walk above the pair. Below it,
    ! iproc's nodes are unknown here: bound them by its size profile. The leaf
    ! itself, or a leaf of iproc's tree, is known exactly
    if (.not.global_active(ia - int(ncells))) return  ! none of its leaves walk
    if (istep == 0 .or. remote_local_index(ia) < 0) then
       may_open = .not.well_separated(node(ia),node(ib))
    else
       call profile_bounds(ia - int(ncells),istep,sizemax,hmaxmax)
       may_open = profile_may_open(node(ia),sizemax,hmaxmax,node(ib))
    endif
    if (.not.may_open) return

    ! a leaf here can be paired with iproc's leaves (its particles may be imported)
    call mark_opened(ib,-(let_round+iproc+1),opened)

    ! only this task's side descends
    if (all(kb <= 0)) return
    do j = 1,2
       if (kb(j) > 0) call push_pair(istack,nstack,ia,kb(j),istep+1)
    enddo
 elseif (ia == ib .or. .not.well_separated(node(ia),node(ib))) then
    if (ia > ncells) then
       if (global_owner(ia - int(ncells)) == iproc) then
          if (.not.global_active(ia - int(ncells))) return  ! no active leaves below it on iproc
       endif
    endif
    if (all(kb <= 0)) then
       ! leaf lowering: only iproc's side descends
       do i = 1,2
          if (ka(i) > 0) call push_pair(istack,nstack,ka(i),ib,0)
       enddo
    else
       call mark_opened(ib,-(let_round+iproc+1),opened)

       ! both sides descend
       do i = 1,2
          do j = 1,2
             if (ka(i) > 0 .and. kb(j) > 0) call push_pair(istack,nstack,ka(i),kb(j),0)
          enddo
       enddo
    endif
 endif

end subroutine let_pair_step

!-----------------------------------------------------------------------
!+
!  add a local node to the nodes opened for a task, if its children are
!  not already in the global tree. The thread's marks (with a value for
!  the task) skip most repeats
!+
!-----------------------------------------------------------------------
subroutine mark_opened(n,mark,opened)
 use io, only:id
 integer,       intent(in)    :: n,mark
 type(intlist), intent(inout) :: opened

 if (n > ncells) return
 if (nodemap(n) /= 0) then
    if (cellatid(nodemap(n)) /= id + 1) return  ! an internal node of the global tree
 endif
 if (let_mark(n) == mark) return
 let_mark(n) = mark
 call add_to_list(opened%x,opened%n,n)

end subroutine mark_opened

!-----------------------------------------------------------------------
!+
!  whether a node of the combined tree belongs to the given task, or is
!  a node of the global tree above its domain
!+
!-----------------------------------------------------------------------
logical function on_side(n,itask)
 use io, only:id
 integer, intent(in) :: n,itask
 integer :: g

 if (n <= ncells) then
    on_side = (itask == id)
    return
 endif
 g = n - int(ncells)
 if (g <= ntop) then
    ! a top node: does it hold sections of itask?
    on_side = (itask >= own_lo(g) .and. itask <= own_hi(g))
 else
    on_side = (global_owner(g) == itask)
 endif

end function on_side

!-----------------------------------------------------------------------
!+
!  whether two nodes are well separated: the acceptance criterion of
!  the dual tree walk (node_interaction)
!+
!-----------------------------------------------------------------------
logical function well_separated(node1,node2)
 use kernel, only:radkern
 type(kdnode), intent(in) :: node1,node2
 real :: dx,dy,dz,xoff,yoff,zoff,r2,rcut

 call get_sep(node1%xcen,node2%xcen,dx,dy,dz,xoff,yoff,zoff,r2)
 rcut = radkern*max(node1%hmax,node2%hmax)
 well_separated = (r2 > (node1%size + node2%size + rcut)**2) .and. &
                  (tree_accuracy**2*r2 > (node1%size + node2%size)**2)

end function well_separated

!-----------------------------------------------------------------------
!+
!  whether a node below the global leaf noded, of size at most sizemax and
!  hmax at most hmaxmax (from its profile), could fail the acceptance
!  criterion against node n. Its centre is within the leaf's sphere (the
!  centre of mass of particles in a ball is in the ball)
!+
!-----------------------------------------------------------------------
logical function profile_may_open(noded,sizemax,hmaxmax,n_node)
 use kernel, only:radkern
 type(kdnode), intent(in) :: noded,n_node
 real,         intent(in) :: sizemax,hmaxmax
 real :: dx,dy,dz,xoff,yoff,zoff,r2,d,sizes,rcut

 call get_sep(noded%xcen,n_node%xcen,dx,dy,dz,xoff,yoff,zoff,r2)
 d     = sqrt(r2) - noded%size
 sizes = sizemax + n_node%size
 rcut  = radkern*max(hmaxmax,n_node%hmax)
 profile_may_open = (d < sizes + rcut) .or. (tree_accuracy*d < sizes)

end function profile_may_open

!-----------------------------------------------------------------------
!+
!  bounds on the size and hmax of a node of another task's global leaf g
!  that a walk pairs istep steps below it, from its profile (or, if none
!  was received, any node in its sphere: at most twice its radius, and
!  no larger hmax than the leaf)
!+
!-----------------------------------------------------------------------
subroutine profile_bounds(g,istep,sizemax,hmaxmax)
 integer, intent(in)  :: g,istep
 real,    intent(out) :: sizemax,hmaxmax
 integer :: k

 if (profile_len(g) > 0) then
    k = profile_start(g) + 2*min(istep,profile_len(g)-1)
    sizemax = profile_pool(k)
    hmaxmax = profile_pool(k+1)
 else
    sizemax = 2.*nodeglobal(g)%size
    hmaxmax = nodeglobal(g)%hmax
 endif

end subroutine profile_bounds

!-----------------------------------------------------------------------
!+
!  send the size profiles of this task's global leaves to the tasks whose
!  walks they can bound (those whose domains could interact with the leaf),
!  receive theirs, and mark which global nodes of other tasks have active
!  leaves below them. A task missing a profile uses a looser bound, so this
!  affects only how many nodes are imported, never which are needed
!+
!-----------------------------------------------------------------------
subroutine exchange_size_profiles()
 use io,       only:nprocs,id
 use mpiutils, only:alltoall_mpi,alltoallv_mpi
 real,    allocatable :: sendbuf(:),recvbuf(:),prof(:,:,:)
 integer, allocatable :: nlen(:)
 logical, allocatable :: interacts(:,:)
 integer :: nsend(nprocs),nrecv(nprocs),iproc,k,nbuf,irecv,g,n,np

 call mark_active_branches()

 ! this task's global leaves: the roots of its sections that are leaves of the global
 ! tree, and the leaves of its block
 if (allocated(own_gleaves)) deallocate(own_gleaves)
 allocate(own_gleaves(nmirror+nown_roots))
 nown_gleaves = 0
 do k = 1,nown_roots
    if (cellatid(nodemap(own_roots(k))) == id + 1) then
       nown_gleaves = nown_gleaves + 1
       own_gleaves(nown_gleaves) = own_roots(k)
    endif
 enddo

 ! and the copied nodes that are leaves of the global tree
 do k = 1,nmirror
    if (cellatid(nodemap(mirror_list(k))) == id + 1) then
       nown_gleaves = nown_gleaves + 1
       own_gleaves(nown_gleaves) = mirror_list(k)
    endif
 enddo

 ! the tasks each leaf can interact with, and the leaf's profile if there are any
 allocate(prof(2,nprofile,nown_gleaves),nlen(nown_gleaves),interacts(nprocs,nown_gleaves))
 !$omp parallel do default(none) schedule(dynamic,16) &
 !$omp shared(nown_gleaves,own_gleaves,nprocs,id,interacts,prof,nlen) &
 !$omp private(k,iproc)
 do k = 1,nown_gleaves
    do iproc = 0,nprocs-1
       interacts(iproc+1,k) = (iproc /= id)
       if (interacts(iproc+1,k)) interacts(iproc+1,k) = leaf_may_interact(own_gleaves(k),iproc)
    enddo
    nlen(k) = 0
    if (any(interacts(:,k))) call leaf_size_profile(own_gleaves(k),prof(:,:,k),nlen(k))
 enddo
 !$omp end parallel do

 ! records of (global index, length, values) for each leaf, to each task it can interact with
 allocate(sendbuf(max(1,nown_gleaves*(2*nprofile+2)*min(nprocs-1,4))))
 nbuf  = 0
 nsend = 0
 do iproc = 0,nprocs-1
    do k = 1,nown_gleaves
       if (.not.interacts(iproc+1,k)) cycle
       n  = own_gleaves(k)
       np = 2*nlen(k)
       if (nbuf + np + 2 > size(sendbuf)) call grow_real(sendbuf,nbuf + np + 2)
       sendbuf(nbuf+1) = nodemap(n)
       sendbuf(nbuf+2) = nlen(k)
       sendbuf(nbuf+3:nbuf+2+np) = reshape(prof(:,1:nlen(k),k),[np])
       nbuf = nbuf + np + 2
       nsend(iproc+1) = nsend(iproc+1) + np + 2
    enddo
 enddo

 ! exchange the records
 call alltoall_mpi(nsend,nrecv)

 allocate(recvbuf(max(sum(nrecv),1)))
 call alltoallv_mpi(sendbuf,nsend,recvbuf,nrecv)

 ! store the profiles received, after undoing those of the last call
 do k = 1,nprof_set
    profile_len(prof_set(k)) = -1
 enddo

 ! and the nodes marked inactive
 do k = 1,ninactive_set
    global_active(inactive_set(k)) = .true.
 enddo
 nprof_set     = 0
 ninactive_set = 0

 ! room for the profiles received
 if (allocated(profile_pool)) then
    if (size(profile_pool) < sum(nrecv)) deallocate(profile_pool)
 endif
 if (.not.allocated(profile_pool)) allocate(profile_pool(max(sum(nrecv),1)))

 ! each record received: the global leaf, then its profile
 irecv = 0
 do while (irecv < sum(nrecv))
    g  = nint(recvbuf(irecv+1))
    np = 2*nint(recvbuf(irecv+2))
    profile_len(g)   = np/2
    profile_start(g) = irecv + 3
    profile_pool(irecv+3:irecv+2+np) = recvbuf(irecv+3:irecv+2+np)
    irecv = irecv + np + 2
    call add_to_list(prof_set,nprof_set,g)
    if (np == 0) call mark_global_inactive(g)
 enddo

end subroutine exchange_size_profiles

!-----------------------------------------------------------------------
!+
!  whether a local node could interact with any node in the domain of
!  another task: the domain region bound, conservative
!+
!-----------------------------------------------------------------------
logical function leaf_may_interact(n,iproc)
 use kernel, only:radkern
 integer, intent(in) :: n,iproc
 integer :: istack(2*maxdepth+2),nstack,g
 real :: dx,dy,dz,xoff,yoff,zoff,r2,d,sizes,rcut

 ! descend the top nodes holding sections of iproc
 leaf_may_interact = .true.
 nstack = 1
 istack(1) = 1
 do while (nstack > 0)
    g = istack(nstack)
    nstack = nstack - 1
    if (iproc < own_lo(g) .or. iproc > own_hi(g)) cycle

    ! skip a top node whose nodes cannot interact with the leaf: their centres
    ! are inside its sphere and they are at most twice its size
    call get_sep(nodeglobal(g)%xcen,node(n)%xcen,dx,dy,dz,xoff,yoff,zoff,r2)
    d     = sqrt(r2) - nodeglobal(g)%size
    sizes = 2.*nodeglobal(g)%size + node(n)%size
    rcut  = radkern*max(nodeglobal(g)%hmax,node(n)%hmax)
    if (d >= sizes + rcut .and. tree_accuracy*d >= sizes) cycle

    ! a section of iproc in reach
    if (section_of(g) > 0) return

    ! otherwise, the top nodes below it
    istack(nstack+1) = nodeglobal(g)%rightchild
    istack(nstack+2) = nodeglobal(g)%leftchild
    nstack = nstack + 2
 enddo

 leaf_may_interact = .false.

end function leaf_may_interact

!-----------------------------------------------------------------------
!+
!  mark the local nodes on the branches of active leaves
!+
!-----------------------------------------------------------------------
subroutine mark_active_branches()
 integer :: k,n

 ! from each active leaf up to the root of its section, stopping at a node already marked
 branch_stamp = branch_stamp + 1
 do k = 1,nactive_leaves
    n = active_leaves(k)
    do while (branch_mark(n) /= branch_stamp)
       branch_mark(n) = branch_stamp

       ! stop at the root of a section: above it are nodes of the global tree
       if (node(n)%parent <= 0 .or. node(n)%parent > ncells) exit
       n = node(n)%parent
    enddo
 enddo

end subroutine mark_active_branches

!-----------------------------------------------------------------------
!+
!  a leaf of another task with no active particles (from its profile):
!  mark it, and the global nodes above it with no other active leaves
!  below them, as inactive. A leaf with no profile received (too far away
!  to matter) stays active
!+
!-----------------------------------------------------------------------
subroutine mark_global_inactive(g0)
 integer, intent(in) :: g0
 integer :: g,il,ir

 ! up from the leaf, while the other child of the parent is inactive too
 g = g0
 do
    global_active(g) = .false.
    call add_to_list(inactive_set,ninactive_set,g)
    g = nodeglobal(g)%parent
    if (g <= 0) exit
    il = nodeglobal(g)%leftchild
    ir = nodeglobal(g)%rightchild
    if (il > 0) then
       if (global_active(il)) exit
    endif
    if (ir > 0) then
       if (global_active(ir)) exit
    endif
 enddo

end subroutine mark_global_inactive

!-----------------------------------------------------------------------
!+
!  append a value to a list that grows as needed
!+
!-----------------------------------------------------------------------
subroutine add_to_list(list,n,ival)
 integer, allocatable, intent(inout) :: list(:)
 integer,              intent(inout) :: n
 integer,              intent(in)    :: ival

 if (.not.allocated(list)) allocate(list(1024))
 if (n >= size(list)) call grow_int(list,2*size(list),n)
 n = n + 1
 list(n) = ival

end subroutine add_to_list

!-----------------------------------------------------------------------
!+
!  profile of the active branches of the local subtree below n0: entry k
!  bounds the size and hmax of a node the walk can reach k-1 steps below
!  n0, i.e. a node at that depth on an active branch, or an active leaf
!  above it (which the walk keeps while the other side descends). The
!  last entry covers all deeper steps. nlen is 0 with no active leaves
!+
!-----------------------------------------------------------------------
subroutine leaf_size_profile(n0,prof,nlen)
 integer, intent(in)  :: n0
 real,    intent(out) :: prof(2,nprofile)   ! (size, hmax) for each step
 integer, intent(out) :: nlen
 integer, allocatable :: istack(:,:)
 integer :: nstack,n,idepth,k,j,ichild(2)
 real    :: sdepth(2,nprofile),sleaf(2,nprofile),nodeval(2)

 nlen = 0
 prof = 0.
 if (branch_mark(n0) /= branch_stamp) return

 ! the largest size and hmax of the nodes at each depth on the active branches below n0,
 ! and of the active leaves
 sdepth = 0.
 sleaf  = 0.
 allocate(istack(3,256))
 nstack = 0
 call push_pair(istack,nstack,n0,0,0)
 do while (nstack > 0)
    n      = istack(1,nstack)
    idepth = istack(2,nstack)
    nstack = nstack - 1
    k = min(idepth+1,nprofile)
    nlen = max(nlen,k)
    nodeval = [node(n)%size,node(n)%hmax]
    sdepth(:,k) = max(sdepth(:,k),nodeval)
    ichild = [node(n)%leftchild,node(n)%rightchild]
    if (all(ichild <= 0)) then
       sleaf(:,k) = max(sleaf(:,k),nodeval)
    else
       do j = 1,2
          if (ichild(j) > 0) then
             if (branch_mark(ichild(j)) == branch_stamp) call push_pair(istack,nstack,ichild(j),idepth+1,0)
          endif
       enddo
    endif
 enddo

 ! an active leaf found at a shallower depth stays in the walk at all later steps
 do k = 2,nprofile
    sleaf(:,k) = max(sleaf(:,k),sleaf(:,k-1))
 enddo

 ! the profile: the nodes at each depth, and the leaves found above it
 prof(:,1) = sdepth(:,1)
 do k = 2,nlen
    prof(:,k) = max(sdepth(:,k),sleaf(:,k-1))
 enddo

 ! beyond the deepest active node only those leaves remain
 if (nlen < nprofile) then
    nlen = nlen + 1
    prof(:,nlen) = sleaf(:,nlen-1)
 endif

end subroutine leaf_size_profile

!-----------------------------------------------------------------------
!+
!  write one imported-node record: node data, parent, side, local index
!+
!-----------------------------------------------------------------------
subroutine put_let_record(rec,n,iparent,iside)
 real,    intent(out) :: rec(nlet_record)
 integer, intent(in)  :: n,iparent,iside

 rec(1:3) = node(n)%xcen
 rec(4)   = node(n)%size
 rec(5)   = node(n)%hmax
 rec(6)   = node(n)%mass
 rec(7)   = iparent
 rec(8)   = iside
 rec(9)   = local_index_and_leaf(n)
#ifdef GRAVITY
 rec(10:15) = node(n)%quads
 rec(16:25) = node(n)%octs
#endif

end subroutine put_let_record

!-----------------------------------------------------------------------
!+
!  append the nodes received from task iproc to the combined tree,
!  below its global leaves
!+
!-----------------------------------------------------------------------
subroutine attach_let_nodes(iproc,buf)
 integer, intent(in) :: iproc
 real,    intent(in) :: buf(:)
 integer :: nrec,k,n,nbase,ipar,ip

 ! room for the nodes after the combined tree
 nrec  = size(buf)/nlet_record
 nbase = ncomb
 if (nbase + nrec > size(node)) call grow_combined_tree(nbase + nrec)

 ! the global leaves of iproc given children are no longer leaves of the combined tree
 do k = 1,nrec
    ipar = nint(buf((k-1)*nlet_record+7))
    if (ipar >= 0) cycle
    ipar = int(ncells) - ipar
    if (leaf_is_active(ipar) == 0) cycle  ! already listed
    call add_to_list(let_parent,nlet_parent,ipar)
    leaf_is_active(ipar)    = 0
    remote_leaf_owner(ipar) = 0
 enddo

 ! the data of each node, and its parent: a node earlier in the message, or a global leaf of iproc
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nrec,nbase,buf,node,ncells,remote_local_index) &
 !$omp private(k,ip,n,ipar)
 do k = 1,nrec
    ip = (k-1)*nlet_record
    n  = nbase + k
    node(n)%xcen = buf(ip+1:ip+3)
    node(n)%size = buf(ip+4)
    node(n)%hmax = buf(ip+5)
    node(n)%mass = buf(ip+6)
#ifdef GRAVITY
    node(n)%quads = buf(ip+10:ip+15)
    node(n)%octs  = buf(ip+16:ip+25)
#endif
    node(n)%leftchild  = 0
    node(n)%rightchild = 0
    node(n)%level      = 0
    ipar = nint(buf(ip+7))
    if (ipar < 0) then
       ipar = int(ncells) - ipar
    else
       ipar = nbase + ipar
    endif
    node(n)%parent = ipar
    remote_local_index(n) = nint(buf(ip+9))
 enddo
 !$omp end parallel do

 ! link each node to its parent (the two children of a parent set different links)
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nrec,nbase,buf,node) &
 !$omp private(k,n,ipar)
 do k = 1,nrec
    n    = nbase + k
    ipar = node(n)%parent
    if (nint(buf((k-1)*nlet_record+8)) == 1) then
       node(ipar)%leftchild = n
    else
       node(ipar)%rightchild = n
    endif
 enddo
 !$omp end parallel do

 ! the nodes with no children are leaves of iproc
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nrec,nbase,node,leaf_is_active,remote_leaf_owner,iproc) &
 !$omp private(k,n)
 do k = 1,nrec
    n = nbase + k
    if (node(n)%leftchild > 0 .or. node(n)%rightchild > 0) then
       leaf_is_active(n)    = 0
       remote_leaf_owner(n) = 0
    else
       leaf_is_active(n)    = 1
       remote_leaf_owner(n) = iproc + 1
    endif
 enddo
 !$omp end parallel do

 ncomb = nbase + nrec

end subroutine attach_let_nodes

!-----------------------------------------------------------------------
!+
!  make room for at least nmin nodes in the combined tree, keeping the
!  nodes already there (called outside parallel regions)
!+
!-----------------------------------------------------------------------
subroutine grow_combined_tree(nmin)
 integer, intent(in) :: nmin
 type(kdnode), allocatable :: newnode(:)
 integer :: nnew

 nnew = max(nmin,int(1.5*size(node)))
 allocate(newnode(nnew))
 newnode(1:ncomb) = node(1:ncomb)
 call move_alloc(newnode,node)
 call grow_int(leaf_is_active,nnew,ncomb)
 call grow_int(remote_leaf_owner,nnew,ncomb)
 call grow_int(remote_local_index,nnew,ncomb)

end subroutine grow_combined_tree

!-----------------------------------------------------------------------
!+
!  grow an integer array to size nnew, keeping its first nkeep entries
!  and setting the rest to zero
!+
!-----------------------------------------------------------------------
subroutine grow_int(x,nnew,nkeep)
 integer, allocatable, intent(inout) :: x(:)
 integer,              intent(in)    :: nnew,nkeep
 integer, allocatable :: tmp(:)

 allocate(tmp(nnew))
 tmp = 0
 tmp(1:nkeep) = x(1:nkeep)
 call move_alloc(tmp,x)

end subroutine grow_int

!-----------------------------------------------------------------------
!+
!  grow a real array, keeping its contents
!+
!-----------------------------------------------------------------------
subroutine grow_real(x,nmin)
 real, allocatable, intent(inout) :: x(:)
 integer,           intent(in)    :: nmin
 real, allocatable :: tmp(:)

 allocate(tmp(max(nmin,2*size(x))))
 tmp(1:size(x)) = x
 call move_alloc(tmp,x)

end subroutine grow_real

!-----------------------------------------------------------------------
!+
!  print the statistics of the locally essential tree exchange
!+
!-----------------------------------------------------------------------
subroutine print_let_stats(iunit)
 use io,       only:id,master,nprocs
 use mpiutils, only:reduceall_mpi
 integer, intent(in) :: iunit
 integer(kind=8) :: nsent,nrecv,nrecvmax,nimpl,nimpp

 nsent  = reduceall_mpi('+',nlet_sent)
 nrecv  = reduceall_mpi('+',nlet_recv)
 nrecvmax = reduceall_mpi('max',nlet_recv)
 if (id == master .and. nsent > 0) then
    write(iunit,"(a,i0,a,f10.3,a,f10.3)") ' MPI LET nodes imported (all tasks): ',nrecv,&
          '  MB per task (mean,max): ',1.e-6*8*nlet_record*real(nrecv)/nprocs,&
          ' ',1.e-6*8*nlet_record*real(nrecvmax)
 endif

 ! the leaves imported for force
 nimpl = reduceall_mpi('+',nimp_leaves_recv)
 nimpp = reduceall_mpi('+',nimp_part_recv)
 if (id == master .and. nsent > 0) write(iunit,"(a,i0,a,i0,a)") ' MPI leaves imported for force (all tasks): ',&
       nimpl,' with ',nimpp,' particles'
 nlet_sent = 0
 nlet_recv = 0
 nimp_leaves_recv = 0
 nimp_part_recv = 0

end subroutine print_let_stats

!-----------------------------------------------------------------------
!+
!  the most nodes of its local tree a task copies into the global tree
!+
!-----------------------------------------------------------------------
integer function mirror_budget()
 use io,  only:nprocs
 use dim, only:maxmirror_task,maxmirror_total

 mirror_budget = max(2,min(maxmirror_task,maxmirror_total/nprocs))

end function mirror_budget

!-----------------------------------------------------------------------
!+
!  choose the nodes of the local tree to copy into the global tree.
!  Starting from the children of the root, the node reaching furthest
!  (size plus kernel radius) is split next, if it can interact with the
!  domain of another task, until the node budget is used. Other tasks
!  then see this task finely only where they need to: near their domains,
!  and around particles with large smoothing lengths. The nodes are
!  listed with the children of the root first
!+
!-----------------------------------------------------------------------
subroutine select_mirrored_nodes(nactive)
 integer, intent(in) :: nactive
 integer, allocatable :: heapnode(:)
 real,    allocatable :: heapkey(:)
 integer :: nbudget,nheap,n,k

 ! copy a number of nodes proportional to the active particles: gathering and syncing
 ! them then costs a small fraction of the work on those particles, however many tasks
 ! there are, and deep in the timestep hierarchy (few active, the tree rebuilt often,
 ! few exports) only a few nodes are copied
 nbudget = min(mirror_budget(),max(nmirror_min,nactive/nactive_per_mirror,2*nown_roots))
 if (allocated(mirror_list)) then
    if (size(mirror_list) < nbudget) deallocate(mirror_list)
 endif
 if (.not.allocated(mirror_list)) allocate(mirror_list(nbudget))
 allocate(heapnode(nbudget),heapkey(nbudget))

 ! start from the children of the roots of this task's sections
 nmirror = 0
 nheap = 0
 do k = 1,nown_roots
    call add_children(own_roots(k),heapnode,heapkey,nheap)
 enddo

 ! then copy the children of the nodes of largest reach that reach another task's domain
 do while (nheap > 0 .and. nmirror + 2 <= nbudget)
    call heap_pop(heapnode,heapkey,nheap,n)
    if (node(n)%leftchild <= 0 .and. node(n)%rightchild <= 0) cycle  ! a leaf of the local tree
    if (reaches_other_domain(n)) call add_children(n,heapnode,heapkey,nheap)
 enddo

end subroutine select_mirrored_nodes

!-----------------------------------------------------------------------
!+
!  add the children of a node to the list of nodes to copy, and to the
!  heap of nodes that may be split further
!+
!-----------------------------------------------------------------------
subroutine add_children(n,heapnode,heapkey,nheap)
 use kernel, only:radkern
 integer, intent(in)    :: n
 integer, intent(inout) :: heapnode(:),nheap
 real,    intent(inout) :: heapkey(:)
 integer :: ichild(2),k,nc

 ! copy each child, and put it on the heap by its reach
 ichild = [node(n)%leftchild,node(n)%rightchild]
 do k = 1,2
    nc = ichild(k)
    if (nc <= 0) cycle
    nmirror = nmirror + 1
    mirror_list(nmirror) = nc
    call heap_push(heapnode,heapkey,nheap,nc,node(nc)%size + radkern*node(nc)%hmax)
 enddo

end subroutine add_children

!-----------------------------------------------------------------------
!+
!  whether a cell of another task can reach a local node, either as a
!  neighbour (within the kernel radius of either) or by opening it in
!  the gravity walk. The domains of the other tasks are bounded by
!  their nodes in the top levels of the global tree
!+
!-----------------------------------------------------------------------
logical function reaches_other_domain(n)
 use io,     only:id
 use kernel, only:radkern
 integer, intent(in) :: n
 integer :: istack(2*maxdepth+2),nstack,g
 real    :: dx,dy,dz,xoff,yoff,zoff,r,sizen,sized,reachn,reachd

 ! descend the top nodes holding sections of other tasks
 reaches_other_domain = .true.
 sizen  = node(n)%size
 reachn = sizen + radkern*node(n)%hmax
 nstack = 1
 istack(1) = 1
 do while (nstack > 0)
    g = istack(nstack)
    nstack = nstack - 1
    if (own_lo(g) == id .and. own_hi(g) == id) cycle  ! only sections of this task below

    ! skip a top node whose cells cannot reach the node: they are at least
    ! r - sized away from it and at most twice sized
    call get_sep(node(n)%xcen,nodeglobal(g)%xcen,dx,dy,dz,xoff,yoff,zoff)
    r      = sqrt(dx*dx + dy*dy + dz*dz)
    sized  = nodeglobal(g)%size
    reachd = sized + radkern*nodeglobal(g)%hmax
    if (r >= reachn + reachd .and. .not.(gravity .and. tree_accuracy*(r - sized) < sizen + 2.*sized)) cycle

    ! a section of another task in reach
    if (section_of(g) > 0) return

    ! otherwise, the top nodes below it
    istack(nstack+1) = nodeglobal(g)%rightchild
    istack(nstack+2) = nodeglobal(g)%leftchild
    nstack = nstack + 2
 enddo

 reaches_other_domain = .false.

end function reaches_other_domain

!-----------------------------------------------------------------------
!+
!  push a node onto a max-heap ordered by key
!+
!-----------------------------------------------------------------------
pure subroutine heap_push(heapnode,heapkey,nheap,n,key)
 integer, intent(inout) :: heapnode(:),nheap
 real,    intent(inout) :: heapkey(:)
 integer, intent(in)    :: n
 real,    intent(in)    :: key
 integer :: i,iparent

 ! move the new entry up while its parent has a smaller key
 nheap = nheap + 1
 i = nheap
 do while (i > 1)
    iparent = i/2
    if (heapkey(iparent) >= key) exit
    heapnode(i) = heapnode(iparent)
    heapkey(i)  = heapkey(iparent)
    i = iparent
 enddo

 heapnode(i) = n
 heapkey(i)  = key

end subroutine heap_push

!-----------------------------------------------------------------------
!+
!  pop the node with the largest key from a max-heap
!+
!-----------------------------------------------------------------------
pure subroutine heap_pop(heapnode,heapkey,nheap,n)
 integer, intent(inout) :: heapnode(:),nheap
 real,    intent(inout) :: heapkey(:)
 integer, intent(out)   :: n
 integer :: i,ichild,nlast
 real    :: keylast

 ! take the top, then move the last entry down from the top while a child has a larger key
 n = heapnode(1)
 nlast   = heapnode(nheap)
 keylast = heapkey(nheap)
 nheap = nheap - 1
 i = 1
 do
    ichild = 2*i
    if (ichild > nheap) exit
    if (ichild < nheap) then
       if (heapkey(ichild+1) > heapkey(ichild)) ichild = ichild + 1
    endif
    if (keylast >= heapkey(ichild)) exit
    heapnode(i) = heapnode(ichild)
    heapkey(i)  = heapkey(ichild)
    i = ichild
 enddo

 if (nheap > 0) then
    heapnode(i) = nlast
    heapkey(i)  = keylast
 endif

end subroutine heap_pop

!-----------------------------------------------------------------------
!+
!  copy the chosen local nodes into a block, with the links between
!  them given as positions in the block (0 for nodes outside it)
!+
!-----------------------------------------------------------------------
subroutine pack_mirrored_nodes(block)
 type(kdnode), intent(out) :: block(:)
 integer :: p,n

 ! position in the block, temporarily stored in nodemap
 !$omp parallel default(none) &
 !$omp shared(nmirror,nodemap,mirror_list,block,node) &
 !$omp private(p,n)
 !$omp do schedule(static)
 do p = 1,nmirror
    nodemap(mirror_list(p)) = p
 enddo
 !$omp enddo

 ! each copied node, with its links as positions in the block
 !$omp do schedule(static)
 do p = 1,nmirror
    n = mirror_list(p)
    block(p) = node(n)
    block(p)%leftchild  = position_in_block(node(n)%leftchild)
    block(p)%rightchild = position_in_block(node(n)%rightchild)
    block(p)%parent     = position_in_block(node(n)%parent)
 enddo
 !$omp enddo
 !$omp end parallel

end subroutine pack_mirrored_nodes

!-----------------------------------------------------------------------
!+
!  position of a local node in the block of copied nodes, 0 if not in it
!+
!-----------------------------------------------------------------------
integer function position_in_block(n)
 integer, intent(in) :: n

 position_in_block = 0
 if (n > 0) position_in_block = nodemap(n)

end function position_in_block

!-----------------------------------------------------------------------
!+
!  pass the hmax of the copied local nodes to the global tree on all
!  MPI tasks, then raise the hmax of the domain nodes and the levels
!  above them to cover their children
!+
!-----------------------------------------------------------------------
subroutine sync_hmax_mpi
 use mpiutils, only:allgatherv_mpi
 use io,       only:id
 use io,       only:nprocs
 use kdtree,   only:joinnode,njoinnode
 real,    allocatable :: hmirror(:),hall(:),hroot(:),hsec(:)
 integer, allocatable :: nsec_task(:)
 integer :: p,k,n,g,il,ir
 real    :: hchild

 ! the hmax of the nodes each task copied into the global tree
 allocate(hmirror(nmirror),hall(sum(nmirror_task)))
 do p = 1,nmirror
    hmirror(p) = node(mirror_list(p))%hmax
 enddo

 ! into the global tree, after the top nodes
 call allgatherv_mpi(hmirror,hall,nmirror_task,mirror_displ)
 do k = 1,size(hall)
    nodeglobal(ntop+k)%hmax = hall(k)
 enddo

 ! the roots of the sections (some are leaves of the global tree, with no copied children)
 allocate(hroot(nown_roots),hsec(nsec_top),nsec_task(nprocs))
 do k = 1,nown_roots
    hroot(k) = node(own_roots(k))%hmax
 enddo

 ! gathered in the order of their owners, into the sections' leaves
 nsec_task = sec_start(2:nprocs+1) - sec_start(1:nprocs)
 call allgatherv_mpi(hroot,hsec,nsec_task,sec_start(1:nprocs)-1)
 do k = 1,nsec_top
    g = top_leaf(sec_by_owner(k))
    nodeglobal(g)%hmax = max(nodeglobal(g)%hmax,hsec(k))
 enddo

 ! raise the hmax of each section and the top nodes above it to cover its children
 do k = 1,nsec_top
    g  = top_leaf(k)
    il = nodeglobal(g)%leftchild
    ir = nodeglobal(g)%rightchild
    hchild = nodeglobal(g)%hmax
    if (il > 0) hchild = max(hchild,nodeglobal(il)%hmax)
    if (ir > 0) hchild = max(hchild,nodeglobal(ir)%hmax)

    ! the section, and the top nodes above it
    n = g
    nodeglobal(n)%hmax = max(nodeglobal(n)%hmax,hchild)
    do while (nodeglobal(n)%parent /= 0)
       n = nodeglobal(n)%parent
       nodeglobal(n)%hmax = max(nodeglobal(n)%hmax,hchild)
    enddo
 enddo

 ! the join nodes of the local tree above the roots of the sections, which the local
 ! walks pass through (listed children first)
 do k = 1,njoinnode
    n  = joinnode(k)
    il = node(n)%leftchild
    ir = node(n)%rightchild
    node(n)%hmax = max(node(n)%hmax,node(il)%hmax,node(ir)%hmax)
 enddo

 ! the same values in the combined tree
 do k = 1,ncomb_tree - int(ncells)  ! the global tree only, not the imported nodes
    if (.not.is_this_task(k)) node(int(ncells) + k)%hmax = nodeglobal(k)%hmax
 enddo

 ! and the roots of this task's sections, which stand for their leaves
 do k = 1,nsec_top
    if (top_owner(k) == id) node(top_root(k))%hmax = nodeglobal(top_leaf(k))%hmax
 enddo

end subroutine sync_hmax_mpi

!-----------------------------------------------------------
!+
!  get the separation in 3D between two nodes of the tree
!  (mpikdtree and neighkdtree each have a copy, so that it can be
!  inlined in their tree walks without link-time optimisation)
!+
!-----------------------------------------------------------
pure subroutine get_sep(x1,x2,dx,dy,dz,xoffset,yoffset,zoffset,r2)
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound,hdlx,hdly,hdlz
#endif
 real, intent(in)  :: x1(3),x2(3)
 real, intent(out) :: dx,dy,dz,xoffset,yoffset,zoffset
 real, intent(out), optional :: r2

 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 dx = x2(1) - x1(1)
 dy = x2(2) - x1(2)
 dz = x2(3) - x1(3)

#ifdef PERIODIC
 if (abs(dx) > hdlx) then ! mod distances across boundary if periodic BCs
    xoffset = -dxbound*SIGN(1.0,dx)
    dx = dx + xoffset
 endif
 if (abs(dy) > hdly) then
    yoffset = -dybound*SIGN(1.0,dy)
    dy = dy + yoffset
 endif
 if (abs(dz) > hdlz) then
    zoffset = -dzbound*SIGN(1.0,dz)
    dz = dz + zoffset
 endif
#endif

 if (present(r2)) r2 = dx*dx+dy*dy+dz*dz

end subroutine get_sep

end module mpikdtree
