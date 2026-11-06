!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module kdtree
!
! This module implements the k-d tree build
!    and associated tree walking routines
!
! :References:
!    Gafton & Rosswog (2011), MNRAS 418, 770-781
!    Benz, Bowers, Cameron & Press (1990), ApJ 348, 647-667
!    Dehnen (2000), ApJL 536, L39; Dehnen (2002), JCoPh 179, 27
!    Marcello (2017), AJ 154, 92 (angular momentum conservation)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, io, mpidomain,
!   mpisections, mpiutils, part, timing
!
 use dim,         only:maxp,ncellsmax,minpart,use_apr,use_sinktree,maxptmass,maxpsph,gravity
 use io,          only:nprocs
 use dtypekdtree, only:kdnode,lenfgrav
 use part,        only:ll,iphase,treecache,maxphase, &
                       apr_level,aprmassoftype

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)

!
!--tree parameters
!
 integer,          parameter, public :: irootnode    = 1
 character(len=1), parameter, public :: labelax(3)   = (/'x','y','z'/)
 integer,          parameter, public :: maxdepth     = 128
!
!--runtime options for this module
!
 real,    public  :: tree_accuracy    = 0.5
 logical, public  :: use_geosplit     = gravity ! only debug flag / should be on with gravity
 logical, private :: done_init_kdtree = .false.
 logical, private :: already_warned   = .false.
 logical, private :: tree_is_revisable = .false. ! last build was of the particle arrays
 ! revise rather than rebuild the tree if fewer than this fraction of particles are active,
 ! until the particles have had about one step each since the tree was built
 real,    parameter, public :: frac_active_revtree = 0.1
 real,    private :: nactive_since_build = 0.         ! active particles (all tasks) since the last build
 integer, private :: numthreads
 ! scratch space for the parallel partition in build_top_parallel
 real,    allocatable, private :: tcbuf(:,:)
 integer, allocatable, private :: ipbuf(:)

 public :: allocate_kdtree, deallocate_kdtree
 public :: maketree, revtree,kdnode
 public :: empty_tree,combine_child_nodes,push_pair
 integer, public :: maxlevel_indexed, maxlevel

 ! with the domains from sections: the root of the subtree of each section this task owns
 ! (secroot, 0 for a section with no particles here) and the join nodes above them, children first
 integer, allocatable, public :: secroot(:),joinnode(:)
 integer, public :: nsecroot = 0,njoinnode = 0

 ! neighbour cache indices (xyzcache); imported with only: from dens/force
 integer, parameter, public :: ix=1, iy=2, iz=3, ih1=4, im=5

 type kdbuildstack
    integer :: node
    integer :: parent
    integer :: level
    integer :: npnode
    real    :: xmin(3)
    real    :: xmax(3)
 end type kdbuildstack

 private

contains

subroutine allocate_kdtree
 use allocutils, only:allocate_array

 call allocate_array('inoderange', inoderange, 2, ncellsmax+1)
 call allocate_array('inodeparts', inodeparts, maxp)

end subroutine allocate_kdtree

subroutine deallocate_kdtree
 if (allocated(inoderange)) deallocate(inoderange)
 if (allocated(inodeparts)) deallocate(inodeparts)
 if (allocated(tcbuf)) deallocate(tcbuf,ipbuf)

end subroutine deallocate_kdtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the tree from scratch
!
!  Notes/To do:
!  -openMP parallelisation of maketree_stack (done - April 2013)
!  -test centre of mass vs. geometric centre based cell sizes (done - 2013)
!  -need analysis module that times (and checks) build_tree and tree walk
!   for a given dump
!  -test bottom-up vs. top-down neighbour search
!  -should we try to store tree structure with particle arrays?
!  -need to compute centre of mass and moments for each cell on the fly (c.f. revtree?)
!  -need to implement long-range gravitational interaction (done - May 2013)
!  -implement revtree routine to update tree w/out rebuilding (done - Sep 2015)
!+
!-------------------------------------------------------------------------------
subroutine maketree(node, xyzh, np, leaf_is_active, ncells, apr_tree, nptmass,xyzmh_ptmass,by_section)
 use io,   only:fatal,warning,iprint,iverbose
 use timing, only:get_timings,increment_timer,itimer_tree_top,itimer_tree_sub
!$ use omp_lib
 type(kdnode),    intent(out)   :: node(:) !ncellsmax+1)
 integer,         intent(in)    :: np
 real,            intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,         intent(out)   :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 logical,         intent(in)    :: apr_tree
 integer,         intent(in),    optional :: nptmass
 real,            intent(inout), optional :: xyzmh_ptmass(:,:)
 logical,         intent(in),    optional :: by_section  ! one subtree per section owned (MPI)

 integer :: i,j,npnode,il,ir,istack,nl,nr,mymum,ntarget,njoin
 integer :: nnode,minlevel,level,nqueue
 real :: xmini(3),xmaxi(3),xminl(3),xmaxl(3),xminr(3),xmaxr(3)
 integer, parameter :: istacksize = 512, maxqueue = 4096
 type(kdbuildstack), save :: stack(istacksize)
 !$omp threadprivate(stack)
 type(kdbuildstack), allocatable :: queue(:)
 integer, allocatable :: join(:,:)
!$ integer :: threadid
 integer :: npcounter
 logical :: wassplit,finished,sinktree,seeded
 character(len=10) :: string
 real(kind=4) :: t1,t2,tcpu1,tcpu2

 call get_timings(t1,tcpu1)
 ! must initialise: otherwise a garbage .true. takes the sinktree branch
 ! with absent optionals (illegal access of nptmass/xyzmh_ptmass)
 sinktree = .false.
 if (present(nptmass) .and. present(xyzmh_ptmass)) then
    sinktree = .true.
 endif

 leaf_is_active = 0

 ir = 0
 il = 0
 nl = 0
 nr = 0
 wassplit = .false.
 finished = .false.

 ! construct root node, i.e. find bounds of all particles
 if (sinktree) then
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,apr_tree,xyzmh_ptmass,nptmass)
 else
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,apr_tree)
 endif

 if (inoderange(1,irootnode)==0 .or. inoderange(2,irootnode)==0 ) then
    call fatal('maketree','no particles or all particles dead/accreted')
 endif
 ! the tree used to merge particles is not of the particle arrays, so cannot be revised
 tree_is_revisable = .not.apr_tree
 nactive_since_build = 0.

! Put root node on top of stack
 ncells = 1
 maxlevel = 0
 minlevel = maxdepth - 1
 istack = 1

 ! maximum level where 2^k indexing can be used (thus avoiding critical sections)
 ! deeper than this we access cells via a stack as usual
 if (use_geosplit) then
    maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 2
 else
    maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 1
 endif

 ! default number of cells is the size of the `indexed' part of the tree
 ! this can be *increased* by building tree beyond indexed levels
 ! and is decreased afterwards according to the maximum depth actually reached
 ncells = 2**(maxlevel_indexed+1) - 1

 ! need to number of particles in node during build
 ! this is counted above to remove dead/accreted particles
 allocate(queue(maxqueue))
 call push_onto_stack(queue(istack),irootnode,0,0,npcounter,xmini,xmaxi)

 if (.not.done_init_kdtree) then
    ! 1 thread for serial, overwritten when using OpenMP
    numthreads = 1

    ! get number of OpenMP threads
    !$omp parallel default(none) shared(numthreads)
!$  numthreads = omp_get_num_threads()
    !$omp end parallel
    done_init_kdtree = .true.
 endif

 ! with MPI, the subtrees start from the sections this task owns, below
 ! join nodes that are filled in once the subtrees are built
 seeded = .false.
 if (present(by_section)) seeded = by_section .and. .not.apr_tree
 if (seeded) call seed_queue_with_sections(npcounter,queue,istack,join,njoin,seeded,ncells)

 nqueue = numthreads
 ! a subtree for each thread, none with much more than a thread's share of the particles
 ntarget = max(int(1.25*npcounter/nqueue),max(minpart,64))
 ! build the first levels with all threads working on every level (not for the APR
 ! merge tree, whose partition has to leave an even number of particles in each child)
 if (.not.seeded .and. .not.apr_tree .and. nqueue > 1 .and. npcounter > max(minpart,64)) then
    call build_top_parallel(node,queue,istack,ntarget,leaf_is_active)
 endif

 ! then split the largest node in the queue until there is a subtree for each thread
 ! and none has much more than a thread's share of the particles, e.g. a dense core
 ! that build_top_parallel cannot split. The subtrees are then handed out largest first
 over_queue: do while (istack < nqueue .or. maxval(queue(1:max(istack,1))%npnode) > ntarget)
    if (istack + 2 > size(queue)) exit over_queue
    ! if the tree finished while building the queue, then we should just return
    ! only happens for small particle numbers
    if (istack <= 0) then
       finished = .true.
       exit over_queue
    endif

    ! pop the node with most particles off the queue, so the subtrees end up of similar size
    j = maxloc(queue(1:istack)%npnode,1)
    call pop_off_stack(queue(j), istack, nnode, mymum, level, npnode, xmini, xmaxi)

    ! shuffle the rest of the queue forward
    do i=j,istack
       queue(i) = queue(i+1)
    enddo

    ! construct node
    if (sinktree) then
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, apr_tree, xyzmh_ptmass)
    else
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, apr_tree)
    endif

    if (wassplit) then ! add children to back of queue
       istack = istack + 1
       call push_onto_stack(queue(istack),il,nnode,level+1,nl,xminl,xmaxl)
       istack = istack + 1
       call push_onto_stack(queue(istack),ir,nnode,level+1,nr,xminr,xmaxr)
    endif

 enddo over_queue
 ! build_top_parallel completes whole levels, so the queue may hold more than nqueue nodes
 if (.not.finished) nqueue = istack
 if (.not.finished) call sort_queue_by_size(queue,nqueue)

 ! the levels built with all threads, then the subtrees built by one thread each
 call get_timings(t2,tcpu2)
 if (.not.apr_tree) call increment_timer(itimer_tree_top,t2-t1,tcpu2-tcpu1)

 done: if (.not.finished) then

    ! build using a stack which builds depth first
    ! each thread grabs a node from the queue and builds its own subtree

    !$omp parallel default(none) &
    !$omp shared(queue) &
    !$omp shared(ll, leaf_is_active) &
    !$omp shared(xyzmh_ptmass) &
    !$omp shared(np) &
    !$omp shared(node, ncells) &
    !$omp shared(nqueue,apr_tree,sinktree) &
    !$omp private(istack) &
    !$omp private(nnode, mymum, level, npnode, xmini, xmaxi) &
    !$omp private(ir, il, nl, nr) &
    !$omp private(xminr, xmaxr, xminl, xmaxl) &
    !$omp private(threadid) &
    !$omp private(wassplit) &
    !$omp reduction(min:minlevel) &
    !$omp reduction(max:maxlevel)
    !$omp do schedule(dynamic,1)
    do i = 1, nqueue

       stack(1) = queue(i)
       istack = 1

       over_stack: do while(istack > 0)

          ! pop node off top of stack
          call pop_off_stack(stack(istack), istack, nnode, mymum, level, npnode, xmini, xmaxi)

          ! construct node
          if (sinktree) then
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, apr_tree, xyzmh_ptmass)
          else
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, apr_tree)
          endif

          if (wassplit) then ! add children to top of stack
             if (istack+2 > istacksize) call fatal('maketree',&
                                       'stack size exceeded in tree build, increase istacksize and recompile')

             istack = istack + 1
             call push_onto_stack(stack(istack),il,nnode,level+1,nl,xminl,xmaxl)
             istack = istack + 1
             call push_onto_stack(stack(istack),ir,nnode,level+1,nr,xminr,xmaxr)
          endif

       enddo over_stack
    enddo
    !$omp enddo
    !$omp end parallel

 endif done
 if (seeded) call combine_join_nodes(node,join,njoin,leaf_is_active)
 if (present(by_section)) then
    if (allocated(joinnode)) deallocate(joinnode)
    njoinnode = 0
    if (seeded) njoinnode = njoin
    allocate(joinnode(max(njoinnode,1)))
    if (seeded) joinnode(1:njoin) = join(1,njoin:1:-1)  ! children before their parents
 endif
 call get_timings(t1,tcpu1)
 if (.not.apr_tree) call increment_timer(itimer_tree_sub,t1-t2,tcpu1-tcpu2)

 ! store the nodes contiguously, in depth first order: the indexing used while building
 ! leaves most slots of the node array empty, and the walks are faster on a compact tree
 if (apr_tree) then
    ! decrease number of cells if tree is entirely within 2^k indexing limit
    if (maxlevel < maxlevel_indexed) ncells = 2**(maxlevel+1) - 1
 else
    if (present(by_section)) then
       call renumber_depth_first(node,leaf_is_active,ncells,secroot(1:nsecroot),joinnode(1:njoinnode))
    else
       call renumber_depth_first(node,leaf_is_active,ncells)
    endif
 endif

 if (maxlevel > maxlevel_indexed .and. .not.already_warned) then
    write(string,"(i10)") 2**(maxlevel-maxlevel_indexed)
    if (iverbose > 0) call warning('maketree','maxlevel > max_indexed: will run faster if recompiled with '// &
               'NCELLSMAX='//trim(adjustl(string))//'*maxp,')
 endif

 if (iverbose >= 3) then
    write(iprint,"(a,i10,3(a,i2))") ' maketree: nodes = ',ncells,', max level = ',maxlevel,&
       ', min leaf level = ',minlevel,' max level indexed = ',maxlevel_indexed
 endif

end subroutine maketree

!--------------------------------------------------------------------------------
!+
!  renumber the nodes of the tree in depth first order from the root (which
!  stays node 1), so that they fill node(1:ncells) with no gaps, remapping
!  the links between them and the arrays indexed by node
!+
!--------------------------------------------------------------------------------
subroutine renumber_depth_first(node,leaf_is_active,ncells,list1,list2)
 type(kdnode),    intent(inout) :: node(:)
 integer,         intent(inout) :: leaf_is_active(:)
 integer(kind=8), intent(inout) :: ncells
 integer,         intent(inout), optional :: list1(:),list2(:)  ! node indices to renumber (0 kept)
 integer,      allocatable :: order(:),newindex(:)
 integer,      allocatable :: item(:),nitem(:),ioffset(:)
 logical,      allocatable :: is_subtree(:)
 integer :: istack(2*maxdepth+2),nstack,n,k,nnode,nold,nitems,nsmall

 ! the new index of each node, and the size of a subtree listed by one thread
 nold = int(ncells)
 allocate(order(nold),newindex(0:nold),item(nold),is_subtree(nold))
 newindex(0) = 0
 nsmall = max((inoderange(2,irootnode) - inoderange(1,irootnode) + 1)/(8*max(numthreads,1)),1)

 ! the top of the tree, depth first with the left child first, down to the
 ! subtrees with few enough particles to be listed by one thread each
 nitems = 0
 nstack = 1
 istack(1) = irootnode
 do while (nstack > 0)
    n = istack(nstack)
    nstack = nstack - 1
    nitems = nitems + 1
    item(nitems) = n
    is_subtree(nitems) = (node(n)%leftchild <= 0 .or. inoderange(2,n) - inoderange(1,n) + 1 <= nsmall)
    if (.not.is_subtree(nitems)) then
       istack(nstack+1) = node(n)%rightchild
       istack(nstack+2) = node(n)%leftchild
       nstack = nstack + 2
    endif
 enddo

 ! the number of nodes of each item, and where its nodes start in the new order
 allocate(nitem(nitems),ioffset(nitems))
 !$omp parallel do default(none) schedule(dynamic,1) &
 !$omp shared(nitems,item,is_subtree,nitem,node) &
 !$omp private(k)
 do k = 1,nitems
    nitem(k) = 1
    if (is_subtree(k)) nitem(k) = count_subtree(node,item(k))
 enddo
 !$omp end parallel do

 ! where each item's nodes start in the new numbering
 ioffset(1) = 0
 do k = 2,nitems
    ioffset(k) = ioffset(k-1) + nitem(k-1)
 enddo
 nnode = ioffset(nitems) + nitem(nitems)

 ! the order of the nodes, each subtree listed by one thread
 !$omp parallel do default(none) schedule(dynamic,1) &
 !$omp shared(nitems,item,is_subtree,ioffset,node,order,newindex) &
 !$omp private(k)
 do k = 1,nitems
    if (is_subtree(k)) then
       call order_subtree(node,item(k),ioffset(k),order,newindex)
    else
       order(ioffset(k)+1) = item(k)
       newindex(item(k))   = ioffset(k) + 1
    endif
 enddo
 !$omp end parallel do

 ! move the nodes to their new places, with their links renumbered
 call copy_renumbered_nodes(nold,nnode,order,newindex,node,leaf_is_active,inoderange)

 if (present(list1)) list1 = newindex(max(list1,0))
 if (present(list2)) list2 = newindex(max(list2,0))
 ncells = nnode

end subroutine renumber_depth_first

!--------------------------------------------------------------------------------
!+
!  move the nodes to their new indices (order: old index of each new
!  index; newindex: new index of each old index), renumbering the links
!  between them, and clear the nodes left over at the end
!+
!--------------------------------------------------------------------------------
subroutine copy_renumbered_nodes(nold,nnode,order,newindex,node,leaf_is_active,noderange)
 integer,      intent(in)    :: nold,nnode
 integer,      intent(in)    :: order(nnode),newindex(0:nold)
 type(kdnode), intent(inout) :: node(nold)
 integer,      intent(inout) :: leaf_is_active(nold),noderange(2,nold)
 type(kdnode), allocatable :: newnode(:)
 integer,      allocatable :: newactive(:),newrange(:,:)
 integer :: k,n

 ! copy the nodes to their new places, with their links renumbered
 allocate(newnode(nnode),newactive(nnode),newrange(2,nnode))
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nnode,order,node,newnode,newindex,newactive,leaf_is_active,newrange,noderange) &
 !$omp private(k,n)
 do k = 1,nnode
    n = order(k)
    newnode(k) = node(n)
    newnode(k)%leftchild  = newindex(max(node(n)%leftchild,0))
    newnode(k)%rightchild = newindex(max(node(n)%rightchild,0))
    newnode(k)%parent     = newindex(max(node(n)%parent,0))
    newactive(k)  = leaf_is_active(n)
    newrange(:,k) = noderange(:,n)
 enddo
 !$omp end parallel do

 ! copy the renumbered nodes back
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(nold,nnode,node,newnode,leaf_is_active,newactive,noderange,newrange) &
 !$omp private(k)
 do k = 1,nold
    if (k <= nnode) then
       node(k) = newnode(k)
       leaf_is_active(k) = newactive(k)
       noderange(:,k)   = newrange(:,k)
    else
       leaf_is_active(k) = 0
    endif
 enddo
 !$omp end parallel do

end subroutine copy_renumbered_nodes

!--------------------------------------------------------------------------------
!+
!  the number of nodes in the subtree below n0 (n0 included)
!+
!--------------------------------------------------------------------------------
integer function count_subtree(node,n0)
 type(kdnode), intent(in) :: node(:)
 integer,      intent(in) :: n0
 integer :: istack(2*maxdepth+2),nstack,n

 ! visit every node below n0
 count_subtree = 0
 nstack = 1
 istack(1) = n0
 do while (nstack > 0)
    n = istack(nstack)
    nstack = nstack - 1
    count_subtree = count_subtree + 1
    if (node(n)%leftchild > 0) then
       istack(nstack+1) = node(n)%rightchild
       istack(nstack+2) = node(n)%leftchild
       nstack = nstack + 2
    endif
 enddo

end function count_subtree

!--------------------------------------------------------------------------------
!+
!  list the subtree below n0 depth first, left child first, from position
!  ioffset+1 of the new order (order: old index of each new position;
!  newindex: new position of each old index)
!+
!--------------------------------------------------------------------------------
subroutine order_subtree(node,n0,ioffset,order,newindex)
 type(kdnode), intent(in)    :: node(:)
 integer,      intent(in)    :: n0,ioffset
 integer,      intent(inout) :: order(:),newindex(0:)
 integer :: istack(2*maxdepth+2),nstack,n,k

 ! list the nodes below n0 depth first, left child first
 k = ioffset
 nstack = 1
 istack(1) = n0
 do while (nstack > 0)
    n = istack(nstack)
    nstack = nstack - 1
    k = k + 1
    order(k)    = n
    newindex(n) = k
    if (node(n)%leftchild > 0) then
       istack(nstack+1) = node(n)%rightchild
       istack(nstack+2) = node(n)%leftchild
       nstack = nstack + 2
    endif
 enddo

end subroutine order_subtree

!--------------------------------------------------------------------------------
!+
!  sort the queue of subtrees by number of particles, largest first
!  (insertion sort: the queue is short)
!+
!--------------------------------------------------------------------------------
subroutine sort_queue_by_size(queue,nqueue)
 type(kdbuildstack), intent(inout) :: queue(:)
 integer,            intent(in)    :: nqueue
 type(kdbuildstack) :: q
 integer :: i,j

 ! insertion sort, largest first
 do i = 2,nqueue
    q = queue(i)
    j = i - 1
    do while (j >= 1)
       if (queue(j)%npnode >= q%npnode) exit
       queue(j+1) = queue(j)
       j = j - 1
    enddo
    queue(j+1) = q
 enddo

end subroutine sort_queue_by_size

!--------------------------------------------------------------------------------
!+
!  seed the build queue with the sections this task owns (MPI): the
!  particle list is sorted by section, each section with particles
!  is the root of a subtree, and the join nodes above them, from the local
!  root down, are the cells of the geometric tree of the sections: each
!  splits the run of sections where the keys of its first and last first
!  differ. The join nodes are listed in join (node, children, parent and
!  level), parents first. If there is
!  only one section, nothing is done and seeded is returned false
!+
!--------------------------------------------------------------------------------
subroutine seed_queue_with_sections(nproot,queue,istack,join,njoin,seeded,ncells)
 use mpisections,   only:local_section_index,owned_section_keys
 use io,            only:fatal
 integer,              intent(in)    :: nproot
 type(kdbuildstack),   intent(inout) :: queue(:)
 integer,              intent(inout) :: istack
 integer, allocatable, intent(out)   :: join(:,:)
 integer,              intent(out)   :: njoin
 logical,              intent(inout) :: seeded
 integer(kind=8),      intent(inout) :: ncells
 integer,    allocatable :: isec(:),istart(:),ifirst(:),ilast(:),todo(:,:),kown(:)
 integer(kind=8), allocatable :: keyfirst(:),keysec(:)
 integer :: nsecown,nsec,k,ntodo,n,a,b,m,level,mymum,ibit,il,ir
 real    :: xmin(3),xmax(3)

 njoin = 0
 allocate(isec(max(nproot,1)))
 call local_section_index(nproot,treecache(:,1:nproot),isec,nsecown)
 if (allocated(secroot)) deallocate(secroot)
 allocate(secroot(max(nsecown,1)))
 secroot = 0
 nsecroot = max(nsecown,1)
 if (nsecown <= 1) then
    ! one section: the whole local tree
    if (nproot > 0) secroot(1) = irootnode
    seeded = .false.
    return
 endif
 call sort_by_section(nproot,isec,nsecown,istart)

 ! the sections with particles, as ranges of the particle list, and their first keys
 call owned_section_keys(keyfirst)
 allocate(ifirst(nsecown),ilast(nsecown),keysec(nsecown),kown(nsecown))
 nsec = 0
 do k = 1,nsecown
    if (istart(k+1) > istart(k)) then
       nsec = nsec + 1
       ifirst(nsec) = istart(k)
       ilast(nsec)  = istart(k+1) - 1
       keysec(nsec) = keyfirst(k)
       kown(nsec)   = k
    endif
 enddo
 if (nsec > size(queue)) call fatal('maketree','more sections than room in the queue')
 if (nsec <= 1) then
    ! one section with particles: the whole local tree
    if (nsec == 1) secroot(kown(1)) = irootnode
    seeded = .false.
    return
 endif

 ! split the run of sections top-down: (node, first and last section, level, parent)
 allocate(join(5,nsec),todo(5,nsec))
 istack = 0
 ntodo = 1
 todo(:,1) = [irootnode,1,nsec,0,0]
 do while (ntodo > 0)
    n = todo(1,ntodo); a = todo(2,ntodo); b = todo(3,ntodo)
    level = todo(4,ntodo); mymum = todo(5,ntodo)
    ntodo = ntodo - 1
    inoderange(:,n) = [ifirst(a),ilast(b)]
    if (a == b) then
       secroot(kown(a)) = n
       call particle_bounds(ifirst(a),ilast(a),xmin,xmax)
       istack = istack + 1
       call push_onto_stack(queue(istack),n,mymum,level,ilast(a)-ifirst(a)+1,xmin,xmax)
    else
       ! children indexed as in construct_node: 2n and 2n+1 down to maxlevel_indexed
       if (level < maxlevel_indexed) then
          il = 2*n
          ir = il + 1
       else
          ncells = ncells + 2
          ir = int(ncells)
          il = ir - 1
          if (ir > ncellsmax) call fatal('maketree',&
             'number of nodes exceeds array dimensions, increase ncellsmax and recompile',ival=int(ncellsmax))
       endif
       njoin = njoin + 1
       join(:,njoin) = [n,il,ir,mymum,level]

       ! the last section on the left: the last with a 0 at the first bit where
       ! the keys of the first and last sections differ
       ibit = 63 - leadz(ieor(keysec(a),keysec(b)))
       m = a
       do while (m < b - 1 .and. .not.btest(keysec(m+1),ibit))
          m = m + 1
       enddo
       todo(:,ntodo+1) = [il,a,m,level+1,n]
       todo(:,ntodo+2) = [ir,m+1,b,level+1,n]
       ntodo = ntodo + 2
    endif
 enddo

end subroutine seed_queue_with_sections

!--------------------------------------------------------------------------------
!+
!  sort the particle list of the tree (inodeparts, treecache) by the section
!  index isec (1 to nsec), keeping the order within each section. Section k
!  then holds entries istart(k) to istart(k+1)-1
!+
!--------------------------------------------------------------------------------
subroutine sort_by_section(n,isec,nsec,istart)
!$ use omp_lib, only:omp_get_max_threads,omp_get_thread_num
 integer,              intent(in)  :: n,nsec
 integer,              intent(in)  :: isec(:)
 integer, allocatable, intent(out) :: istart(:)
 real,    allocatable :: tcswap(:,:)
 integer, allocatable :: ipswap(:),ncount(:,:)
 integer :: i,k,j,nchunk,ic

 if (allocated(tcbuf)) then
    if (size(tcbuf,2) /= size(treecache,2) .or. size(tcbuf,1) /= size(treecache,1)) deallocate(tcbuf,ipbuf)
 endif
 if (.not.allocated(tcbuf)) allocate(tcbuf(size(treecache,1),size(treecache,2)),ipbuf(size(inodeparts)))

 ! counts per section in each chunk of the list, one chunk per thread
 nchunk = 1
!$ nchunk = omp_get_max_threads()
 allocate(istart(nsec+1),ncount(nsec,0:nchunk-1))
 ncount = 0
 !$omp parallel default(none) &
 !$omp shared(n,isec,nsec,nchunk,ncount,istart,treecache,inodeparts,tcbuf,ipbuf) &
 !$omp private(i,k,j,ic)

 ! each thread counts the entries of each section in its chunk
 ic = 0
!$ ic = omp_get_thread_num()
 do i = int((int(ic,8)*n)/nchunk)+1,int((int(ic+1,8)*n)/nchunk)
    ncount(isec(i),ic) = ncount(isec(i),ic) + 1
 enddo
 !$omp barrier

 ! where each chunk's entries of each section start (one thread)
 !$omp single
 j = 1
 do k = 1,nsec
    istart(k) = j
    do ic = 0,nchunk-1
       i = ncount(k,ic)
       ncount(k,ic) = j
       j = j + i
    enddo
 enddo
 istart(nsec+1) = j
 !$omp end single

 ! each chunk scatters its entries, keeping their order within each section
 ic = 0
!$ ic = omp_get_thread_num()
 do i = int((int(ic,8)*n)/nchunk)+1,int((int(ic+1,8)*n)/nchunk)
    k = isec(i)
    j = ncount(k,ic)
    ncount(k,ic) = j + 1
    tcbuf(:,j) = treecache(:,i)
    ipbuf(j)   = inodeparts(i)
 enddo
 !$omp end parallel

 ! the buffer now holds the particle list in the new order: swap rather than copy back
 call move_alloc(treecache,tcswap)
 call move_alloc(tcbuf,treecache)
 call move_alloc(tcswap,tcbuf)
 call move_alloc(inodeparts,ipswap)
 call move_alloc(ipbuf,inodeparts)
 call move_alloc(ipswap,ipbuf)

end subroutine sort_by_section

!--------------------------------------------------------------------------------
!+
!  the bounding box of entries i1 to i2 of the particle list of the tree
!+
!--------------------------------------------------------------------------------
subroutine particle_bounds(i1,i2,xmin,xmax)
 integer, intent(in)  :: i1,i2
 real,    intent(out) :: xmin(3),xmax(3)
 integer :: i

 ! the bounds of the particles i1 to i2
 xmin =  huge(xmin)
 xmax = -huge(xmax)
 do i = i1,i2
    xmin = min(xmin,treecache(1:3,i))
    xmax = max(xmax,treecache(1:3,i))
 enddo

end subroutine particle_bounds

!--------------------------------------------------------------------------------
!+
!  fill in the join nodes above the sections, listed parents first with
!  their children, once the subtrees have been built: from their
!  particles, as for every other node
!+
!--------------------------------------------------------------------------------
subroutine combine_join_nodes(node,join,njoin,leaf_is_active)
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: join(:,:),njoin
 integer,      intent(inout) :: leaf_is_active(:)
 real    :: x0(3),totmass_node
 integer :: k,n,npnode,mymum,level

 ! the join nodes from their particles, children first
 do k = njoin,1,-1
    n = join(1,k)
    npnode = inoderange(2,n) - inoderange(1,n) + 1
    mymum = join(4,k)
    level = join(5,k)
    call compute_nodes_cofm(npnode,n,x0,totmass_node,.true.)
    call set_nodes_properties(npnode,n,x0,totmass_node,mymum,node(n),level,.true.)
    node(n)%mass       = totmass_node
    node(n)%leftchild  = join(2,k)
    node(n)%rightchild = join(3,k)
    leaf_is_active(n)  = 0
 enddo

end subroutine combine_join_nodes

!--------------------------------------------------------------------------------
!+
!  Build the top levels of the tree a whole level at a time, with all threads
!  working on every level.
!
!  The queue loop in maketree constructs these nodes one at a time, and although
!  the centre of mass and size loops are parallel, the partition that splits each
!  node's particles runs on one thread.  That is log2(nthreads) serial passes over
!  every particle, which grows with the thread count.  Here each node's particles
!  are cut into chunks, the node's share of the threads, and the
!  node construction runs over all chunks of all nodes at once: sums per chunk,
!  combined per node, then a stable two-pass partition into a scratch buffer, which
!  then becomes the particle list (the arrays are swapped rather than copied back).
!  The nodes are those construct_node would give, bar round-off in the sums, but
!  the order of particles within each child differs.
!
!  On entry queue(1:istack) holds one level of the tree (the root).  Returns once
!  istack >= nqueue, or earlier, leaving a consistent queue for the serial loop in
!  maketree to finish, if a level has small nodes or would not fit in the queue.
!+
!--------------------------------------------------------------------------------
subroutine build_top_parallel(node,queue,istack,ntarget,leaf_is_active)
 use part,      only:massoftype,igas,npartoftype
 use dim,       only:maxtypes
 use io,        only:fatal
 type(kdnode),       intent(inout) :: node(:)
 type(kdbuildstack), intent(inout) :: queue(:)
 integer,            intent(inout) :: istack
 integer,            intent(in)    :: ntarget  ! nodes with more particles are split
 integer,            intent(inout) :: leaf_is_active(:)
 integer, allocatable :: cnode(:),chunkl(:),chunkr(:),cnl(:),coffl(:),coffr(:)
 integer, allocatable :: nodecl(:),nodecr(:),jnl(:),nodeax(:),chunk(:),inew(:)
 logical, allocatable :: jdegen(:),jsplit(:)
 real,    allocatable :: psum(:,:),phm(:),pr2(:),pbox(:,:),nodecom(:,:),nodecog(:,:),remain(:)
#ifdef GRAVITY
 real,    allocatable :: pmom(:,:)
 real    :: quads(6),octs(10)
#endif
 type(kdbuildstack), allocatable :: newq(:)
 real,    allocatable :: tcswap(:,:)
 integer, allocatable :: ipswap(:)
 integer :: kmax,cmax,k,j,c,m,ncj,nchunk,i,i1,n,nl,pl,pr,iax,nnode,il,ir,isplit,nnew
 integer(kind=8) :: ntot
 real    :: pmassi,dfac,fac,sx,sy,sz,sm,hm,r2,dx,dy,dz,xpiv,totmass,x0(3),bl(6),br(6),share

 if (allocated(tcbuf)) then
    if (size(tcbuf,2) /= size(treecache,2)) deallocate(tcbuf,ipbuf)
 endif
 if (.not.allocated(tcbuf)) allocate(tcbuf(size(treecache,1),size(treecache,2)),ipbuf(size(inodeparts)))

 ! reference mass, as in construct_node, to keep the centre of mass sums well scaled
 pmassi = massoftype(igas)
 if (pmassi <= 0.) pmassi = massoftype(maxloc(npartoftype(2:maxtypes),1)+1)
 dfac = 1.
 if (pmassi > 0.) dfac = 1./pmassi

 !
 ! work arrays for every level, at their largest: a level has at most size(queue)/2 nodes
 ! (its children must fit in the queue), and at most threads + nodes chunks (the shares
 ! rounded down add up to at most the thread count, plus one for each node whose share
 ! rounds down to zero)
 !
 kmax = size(queue)/2
 cmax = numthreads + kmax
 allocate(nodecl(kmax),nodecr(kmax),jnl(kmax),jdegen(kmax),jsplit(kmax),inew(kmax))
 allocate(nodeax(kmax),nodecom(3,kmax),nodecog(3,kmax))
 allocate(chunk(kmax),remain(kmax),newq(2*kmax))
 allocate(cnode(cmax),chunkl(cmax),chunkr(cmax),cnl(cmax),coffl(cmax),coffr(cmax))
 allocate(psum(4,cmax),phm(cmax),pr2(cmax),pbox(12,cmax))
#ifdef GRAVITY
 allocate(pmom(16,cmax))
#endif

 levels: do
    k = istack
    if (k > kmax) exit levels
    !
    ! split together the nodes of the queue with more than the target, if they can be
    ! split and their children are in the levels numbered 2n/2n+1; the others are carried
    ! over as they are, so that a few small nodes (e.g. around a dense core) do not stop
    ! the rest being split by all threads. Stop when no node is split, or the next queue
    ! would not fit
    !
    nnew = 0
    do j = 1,k
       jsplit(j) = (queue(j)%npnode > ntarget .and. queue(j)%npnode > max(minpart,64) &
                    .and. queue(j)%level < maxlevel_indexed)
       inew(j) = nnew + 1
       nnew = nnew + merge(2,1,jsplit(j))
    enddo
    if (nnew == k .or. nnew > size(queue)) exit levels
    !
    ! each node gets its share of the threads as chunks: the share rounded down (at
    ! least one), then the chunks left over go to the nodes that lost most in the
    ! rounding (largest remainder), so a level with no more nodes than threads has
    ! exactly one chunk per thread.  A chunk never crosses a node boundary.
    !
    ntot = sum(int(queue(1:k)%npnode,kind=8),mask=jsplit(1:k))
    do j = 1,k
       chunk(j)  = 0
       remain(j) = -1.
       if (.not.jsplit(j)) cycle   ! no chunks for a node carried over
       !-- first share estimate (how many threads need a node)
       share     = real(numthreads)*real(queue(j)%npnode)/real(ntot)
       chunk(j)  = max(1,int(share)) !-- cap to 1 as a node should at least have a chunk
       remain(j) = share - int(share)
       if (int(share) == 0) remain(j) = -1. ! already given its minimum of one skip from the remainder share
    enddo
    do while (sum(chunk(1:k)) < numthreads)
       j = maxloc(remain(1:k),1)
       if (remain(j) < 0.) exit
       chunk(j) = chunk(j) + 1
       remain(j) = -1.
    enddo
    nchunk = 0
    do j = 1,k
       ncj = min(chunk(j),queue(j)%npnode)
       nodecl(j) = nchunk + 1 !-- node chunk left
       nchunk    = nchunk + ncj
       nodecr(j) = nchunk     !-- node chunk right
       nodeax(j) = maxloc(queue(j)%xmax - queue(j)%xmin,1)   ! split along the longest axis
       nodecog(:,j) = 0.5*(queue(j)%xmin + queue(j)%xmax)
    enddo
    do j = 1,k
       i1  = inoderange(1,queue(j)%node)
       n   = queue(j)%npnode
       ncj = nodecr(j) - nodecl(j) + 1
       do m = 1,ncj
          c = nodecl(j) + (m-1)
          cnode(c) = j
          chunkl(c)   = i1 + int((int(m-1,8)*n)/ncj)
          chunkr(c)   = i1 + int((int(m,8)*n)/ncj) - 1
       enddo
    enddo

    !$omp parallel default(none) &
    !$omp shared(nchunk,k,chunkl,chunkr,cnode,treecache,dfac,psum,nodecl,nodecr) &
    !$omp shared(nodecom,nodecog,nodeax,jnl,jdegen,node,queue,leaf_is_active,inoderange) &
    !$omp shared(pr2,phm,cnl,coffl,coffr,tcbuf,ipbuf,inodeparts,pbox,newq,use_geosplit) &
    !$omp shared(jsplit,inew) &
#ifdef GRAVITY
    !$omp shared(pmom) &
    !$omp private(quads,octs) &
#endif
    !$omp private(c,j,i,i1,n,nl,pl,pr,iax,nnode,il,ir,isplit) &
    !$omp private(pmassi,fac,sx,sy,sz,sm,hm,r2,dx,dy,dz,xpiv,totmass,x0,bl,br)
    !
    ! mass, centre of mass and hmax of each chunk ...
    !
    !$omp do schedule(static,1)
    do c = 1,nchunk
       sx = 0.; sy = 0.; sz = 0.; sm = 0.;
       do i = chunkl(c),chunkr(c)
          pmassi = treecache(5,i)
          fac = pmassi*dfac
          sm  = sm + pmassi
          sx  = sx + fac*treecache(1,i)
          sy  = sy + fac*treecache(2,i)
          sz  = sz + fac*treecache(3,i)
       enddo
       psum(:,c) = (/sx,sy,sz,sm/)
    enddo
    !$omp enddo
    ! ... combined per node, which also fixes the pivot (the centre of mass)
    !$omp do schedule(static)
    do j = 1,k
       if (.not.jsplit(j)) cycle
       totmass = sum(psum(4,nodecl(j):nodecr(j)))
       if (totmass <= 0.) call fatal('mtree','totmass_node==0',val=totmass)
       nodecom(:,j) = sum(psum(1:3,nodecl(j):nodecr(j)),dim=2)/(totmass*dfac)
       node(queue(j)%node)%mass = totmass
    enddo
    !$omp enddo
    !
    ! size (and multipole moments) about the centre of mass, and the number of
    ! particles left of the pivot, in one pass ...
    !
    !$omp do schedule(static,1)
    do c = 1,nchunk
       iax  = nodeax(cnode(c))
       x0   = nodecom(:,cnode(c))
       if (use_geosplit) then
          xpiv = nodecog(iax,cnode(c))
       else
          xpiv = x0(iax)
       endif
       r2 = 0.
       hm = 0.
       nl = 0
#ifdef GRAVITY
       quads = 0.
       octs  = 0.
#endif
       do i = chunkl(c),chunkr(c)
          dx = treecache(1,i) - x0(1)
          dy = treecache(2,i) - x0(2)
          dz = treecache(3,i) - x0(3)
          hm = max(hm,treecache(4,i))
          r2 = max(r2,dx*dx + dy*dy + dz*dz)
          if (treecache(iax,i) <= xpiv) nl = nl + 1
#ifdef GRAVITY
          pmassi = treecache(5,i)
          call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
       enddo
       phm(c) = hm
       pr2(c) = r2
       cnl(c) = nl
#ifdef GRAVITY
       pmom(1:6,c)  = quads
       pmom(7:16,c) = octs
#endif
    enddo
    !$omp enddo
    ! ... combined per node into the node itself, and each chunk's write positions
    !$omp do schedule(static)
    do j = 1,k
       if (.not.jsplit(j)) cycle
       nnode = queue(j)%node
       i1 = inoderange(1,nnode)
       n  = queue(j)%npnode
       r2 = maxval(pr2(nodecl(j):nodecr(j)))
       node(nnode)%xcen   = nodecom(:,j)
       node(nnode)%size   = sqrt(r2) + epsilon(r2)
       node(nnode)%hmax   = maxval(phm(nodecl(j):nodecr(j)))
       node(nnode)%parent = queue(j)%parent
       node(nnode)%level  = queue(j)%level
#ifdef GRAVITY
       node(nnode)%quads  = sum(pmom(1:6,nodecl(j):nodecr(j)),dim=2)
       node(nnode)%octs   = sum(pmom(7:16,nodecl(j):nodecr(j)),dim=2)
#endif
       il = 2*nnode   ! indexing as per Gafton & Rosswog (2011), as in construct_node
       ir = il + 1
       node(nnode)%leftchild  = il
       node(nnode)%rightchild = ir
       leaf_is_active(nnode)  = 0
       nl = sum(cnl(nodecl(j):nodecr(j)))
       ! all particles on one side: split in half without moving them, as construct_node does
       jdegen(j) = (nl == 0 .or. nl == n)
       if (jdegen(j)) nl = n/2
       jnl(j) = nl
       pl = i1
       pr = i1 + nl
       do c = nodecl(j),nodecr(j)
          coffl(c) = pl
          coffr(c) = pr
          pl = pl + cnl(c)
          pr = pr + (chunkr(c) - chunkl(c) + 1 - cnl(c))
       enddo
    enddo
    !$omp enddo
    !
    ! partition: each chunk scatters its particles into the scratch buffer ...
    !
    !$omp do schedule(static,1)
    do c = 1,nchunk
       if (jdegen(cnode(c))) then
          ! nothing moves, but the buffer becomes the particle list: copy as it is
          tcbuf(:,chunkl(c):chunkr(c)) = treecache(:,chunkl(c):chunkr(c))
          ipbuf(chunkl(c):chunkr(c))   = inodeparts(chunkl(c):chunkr(c))
          cycle
       endif
       iax  = nodeax(cnode(c))
       if (use_geosplit) then
          xpiv = nodecog(iax,cnode(c))
       else
          xpiv = nodecom(iax,cnode(c))
       endif
       pl = coffl(c)
       pr = coffr(c)
       do i = chunkl(c),chunkr(c)
          if (treecache(iax,i) <= xpiv) then
             tcbuf(:,pl) = treecache(:,i)
             ipbuf(pl)   = inodeparts(i)
             pl = pl + 1
          else
             tcbuf(:,pr) = treecache(:,i)
             ipbuf(pr)   = inodeparts(i)
             pr = pr + 1
          endif
       enddo
    enddo
    !$omp enddo
    ! the nodes carried over keep their particles where they are: copy them as they are
    !$omp do schedule(dynamic,1)
    do j = 1,k
       if (jsplit(j)) cycle
       i1 = inoderange(1,queue(j)%node)
       n  = queue(j)%npnode
       tcbuf(:,i1:i1+n-1) = treecache(:,i1:i1+n-1)
       ipbuf(i1:i1+n-1)   = inodeparts(i1:i1+n-1)
    enddo
    !$omp enddo
    ! ... the bounding boxes of the children, from the new order in the buffer
    !$omp do schedule(static,1)
    do c = 1,nchunk
       j  = cnode(c)
       isplit = inoderange(1,queue(j)%node) + jnl(j)   ! first particle of the right child
       bl(1:3) =  huge(1.); bl(4:6) = -huge(1.)
       br = bl
       do i = chunkl(c),chunkr(c)
          if (i < isplit) then
             bl(1:3) = min(bl(1:3),tcbuf(1:3,i))
             bl(4:6) = max(bl(4:6),tcbuf(1:3,i))
          else
             br(1:3) = min(br(1:3),tcbuf(1:3,i))
             br(4:6) = max(br(4:6),tcbuf(1:3,i))
          endif
       enddo
       pbox(1:6,c)  = bl
       pbox(7:12,c) = br
    enddo
    !$omp enddo
    ! the children, or the node carried over, as the next queue
    !$omp do schedule(static)
    do j = 1,k
       if (.not.jsplit(j)) then
          newq(inew(j)) = queue(j)
          cycle
       endif
       nnode = queue(j)%node
       i1 = inoderange(1,nnode)
       n  = queue(j)%npnode
       nl = jnl(j)
       il = 2*nnode
       ir = il + 1
       inoderange(1,il) = i1
       inoderange(2,il) = i1 + nl - 1
       inoderange(1,ir) = i1 + nl
       inoderange(2,ir) = i1 + n - 1
       ! children's boxes from the chunks' partial boxes (into bl/br first, so that no
       ! array temporaries are passed to push_onto_stack)
       bl(1:3) = minval(pbox(1:3,nodecl(j):nodecr(j)),dim=2)
       bl(4:6) = maxval(pbox(4:6,nodecl(j):nodecr(j)),dim=2)
       br(1:3) = minval(pbox(7:9,nodecl(j):nodecr(j)),dim=2)
       br(4:6) = maxval(pbox(10:12,nodecl(j):nodecr(j)),dim=2)
       call push_onto_stack(newq(inew(j)),il,nnode,queue(j)%level+1,nl,bl(1:3),bl(4:6))
       call push_onto_stack(newq(inew(j)+1),ir,nnode,queue(j)%level+1,n-nl,br(1:3),br(4:6))
    enddo
    !$omp enddo
    !$omp end parallel

    ! the buffer now holds the particle list in the new order: swap rather than copy back
    call move_alloc(treecache,tcswap)
    call move_alloc(tcbuf,treecache)
    call move_alloc(tcswap,tcbuf)
    call move_alloc(inodeparts,ipswap)
    call move_alloc(ipbuf,inodeparts)
    call move_alloc(ipswap,ipbuf)

    queue(1:nnew) = newq(1:nnew)
    istack = nnew
 enddo levels

end subroutine build_top_parallel

!----------------------------
!+
! routine to empty the tree
!+
!----------------------------
subroutine empty_tree(node)
 type(kdnode), intent(out) :: node(:)
 integer :: i

!$omp parallel do default(none) shared(node) private(i)
 do i=1,size(node)
    node(i)%xcen = 0.
    node(i)%size = 0.
    node(i)%hmax = 0.
    node(i)%leftchild = 0
    node(i)%rightchild = 0
    node(i)%parent = 0
    node(i)%level  = 0
    node(i)%mass   = 0.
#ifdef GRAVITY
    node(i)%quads = 0.
    node(i)%octs  = 0.
#endif
 enddo
!$omp end parallel do

end subroutine empty_tree

!---------------------------------
!+
! routine to construct root node
!+
!---------------------------------
subroutine construct_root_node(np,nproot,irootnode,xmini,xmaxi,leaf_is_active,xyzh,apr_tree,xyzmh_ptmass,nptmass)
 use, intrinsic :: ieee_arithmetic, only:ieee_is_nan
 use boundary, only:cross_boundary
 use mpidomain,only:isperiodic
 use part, only:iphase,iactive
 use part, only:isdead_or_accreted,ibelong
 use io,   only:fatal,id
 use dim,  only:ind_timesteps,mpi,periodic
 use part, only:isink,massoftype,igas,iamtype,maxphase,maxp,aprmassoftype,apr_level,ihsoft
!$ use omp_lib, only:omp_get_max_threads
 integer, intent(in)    :: np,irootnode
 integer, intent(out)   :: nproot
 real,    intent(out)   :: xmini(3), xmaxi(3)
 integer, intent(inout) :: leaf_is_active(:)
 real,    intent(inout) :: xyzh(:,:)
 logical, intent(in)    :: apr_tree
 real,    intent(inout), optional :: xyzmh_ptmass(:,:)
 integer, intent(in),    optional :: nptmass
 integer :: i,ncross,ic,nchunk,nl
 integer, allocatable :: nlive(:)
 real    :: xminpart,yminpart,zminpart,xmaxpart,ymaxpart,zmaxpart
 real    :: xi, yi, zi

 ! seed with the identity of min/max: an empty task (or a dead first particle)
 ! must not change the box reduced across MPI tasks
 xminpart =  huge(xminpart)
 yminpart =  huge(yminpart)
 zminpart =  huge(zminpart)
 xmaxpart = -huge(xmaxpart)
 ymaxpart = -huge(ymaxpart)
 zmaxpart = -huge(zmaxpart)

 ncross = 0
 nproot = 0
 ! the live particles are also counted per chunk of the arrays in this pass, for the copy below
 nchunk = 1
!$ nchunk = omp_get_max_threads()
 allocate(nlive(0:nchunk))
 nlive = 0
 !$omp parallel default(none) &
 !$omp shared(np,xyzh,nptmass,xyzmh_ptmass) &
 !$omp shared(inodeparts,iphase,treecache,nproot) &
 !$omp shared(id,use_sinktree) &
 !$omp shared(isperiodic,nchunk,nlive) &
 !$omp private(i,ic,nl,xi,yi,zi) &
 !$omp reduction(min:xminpart,yminpart,zminpart) &
 !$omp reduction(max:xmaxpart,ymaxpart,zmaxpart) &
 !$omp reduction(+:ncross)
 !$omp do schedule(static)
 do ic=1,nchunk
    nl = 0
    do i=int((int(ic-1,8)*np)/nchunk)+1,int((int(ic,8)*np)/nchunk)
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          nl = nl + 1
          if (periodic) call cross_boundary(isperiodic,xyzh(:,i),ncross)
          xi = xyzh(1,i)
          yi = xyzh(2,i)
          zi = xyzh(3,i)
          if (ieee_is_nan(xi) .or. ieee_is_nan(yi) .or. ieee_is_nan(zi)) then
             call fatal('maketree','NaN in particle position, likely caused by NaN in force',i,var='x',val=xi)
          endif
          xminpart = min(xminpart,xi)
          yminpart = min(yminpart,yi)
          zminpart = min(zminpart,zi)
          xmaxpart = max(xmaxpart,xi)
          ymaxpart = max(ymaxpart,yi)
          zmaxpart = max(zmaxpart,zi)
       endif
    enddo
    nlive(ic) = nl   ! once per chunk: neighbouring counters share a cache line
 enddo
 !$omp enddo
 !$omp barrier
 if (use_sinktree) then
    if (nptmass>0) then
       !$omp do schedule(guided,1)
       do i=1,nptmass
          if (xyzmh_ptmass(4,i)>0.) then
             if (periodic) call cross_boundary(isperiodic,xyzmh_ptmass(1:3,i),ncross)
             xi = xyzmh_ptmass(1,i)
             yi = xyzmh_ptmass(2,i)
             zi = xyzmh_ptmass(3,i)
             if (ieee_is_nan(xi) .or. ieee_is_nan(yi) .or. ieee_is_nan(zi)) then
                call fatal('maketree','NaN in ptmass position, likely caused by NaN in force',i,var='x',val=xi)
             endif
             xminpart = min(xminpart,xi)
             yminpart = min(yminpart,yi)
             zminpart = min(zminpart,zi)
             xmaxpart = max(xmaxpart,xi)
             ymaxpart = max(ymaxpart,yi)
             zmaxpart = max(zmaxpart,zi)
          endif
       enddo
       !$omp enddo
    endif
 endif
 !$omp end parallel

 !
 ! copy the live particles into the tree's list in parallel: from the counts above, each
 ! chunk of the particle arrays knows where its particles go, and the order is the same
 ! as a serial copy
 !
 do ic=1,nchunk
    nlive(ic) = nlive(ic) + nlive(ic-1)
 enddo
 !$omp parallel do schedule(static) default(none) &
 !$omp shared(np,nchunk,xyzh,nlive,inodeparts,treecache,iphase,massoftype,aprmassoftype,apr_level) &
 !$omp shared(maxp,maxphase,apr_tree) &
 !$omp private(ic,i,nproot)
 do ic=1,nchunk
    nproot = nlive(ic-1)
    do i=int((int(ic-1,8)*np)/nchunk)+1,int((int(ic,8)*np)/nchunk)
       isnotdead: if (.not.isdead_or_accreted(xyzh(4,i))) then
          nproot = nproot + 1

          if (ind_timesteps) then
             if (iactive(iphase(i))) then
                inodeparts(nproot) = i  ! +ve if active
             else
                inodeparts(nproot) = -i ! -ve if inactive
             endif
             ! particles in the tree used for merging are not indexed by their
             ! position in the global arrays, so iphase does not apply to them
             if (apr_tree) inodeparts(nproot) = abs(inodeparts(nproot))
          else
             inodeparts(nproot) = i
          endif
          treecache(1:4,nproot) = xyzh(1:4,i)
          if (apr_tree) then
             ! particles being merged all have the same type and level
             treecache(5,nproot) = 1.
          elseif (maxphase==maxp) then
             if (use_apr) then
                treecache(5,nproot) = aprmassoftype(iamtype(iphase(i)),apr_level(i))
             else
                treecache(5,nproot) = massoftype(iamtype(iphase(i)))
             endif
          elseif (use_apr) then
             treecache(5,nproot) = aprmassoftype(igas,apr_level(i))
          else
             treecache(5,nproot) = massoftype(igas)
          endif
       endif isnotdead
    enddo
 enddo
 !$omp end parallel do
 nproot = nlive(nchunk)
 deallocate(nlive)

 if (use_sinktree) then
    if (nptmass > 0) then
       do i=1,nptmass
          if (mpi) then
             if (ibelong(maxpsph+i) /= id) cycle
          endif
          if (xyzmh_ptmass(4,i)<0.) cycle
          nproot = nproot + 1
          inodeparts(nproot) = (maxpsph) + i
          treecache(1:3,nproot) = xyzmh_ptmass(1:3,i)
          treecache(4,nproot)   = xyzmh_ptmass(ihsoft,i)
          treecache(5,nproot)   = xyzmh_ptmass(4,i)
       enddo
    endif
 endif

 if (nproot /= 0) then
    inoderange(1,irootnode) = 1
    inoderange(2,irootnode) = nproot
 else
    inoderange(:,irootnode) = 0
 endif

 xmini(1) = xminpart
 xmini(2) = yminpart
 xmini(3) = zminpart
 xmaxi(1) = xmaxpart
 xmaxi(2) = ymaxpart
 xmaxi(3) = zmaxpart

end subroutine construct_root_node

! also used for queue push
pure subroutine push_onto_stack(stackentry,node,parent,level,npnode,xmin,xmax)
 type(kdbuildstack), intent(out) :: stackentry
 integer,            intent(in)  :: node,parent,level
 integer,            intent(in)  :: npnode
 real,               intent(in)  :: xmin(3),xmax(3)

 stackentry%node   = node
 stackentry%parent = parent
 stackentry%level  = level
 stackentry%npnode = npnode
 stackentry%xmin   = xmin
 stackentry%xmax   = xmax

end subroutine push_onto_stack

! also used for queue pop
pure subroutine pop_off_stack(stackentry, istack, nnode, mymum, level, npnode, xmini, xmaxi)
 type(kdbuildstack), intent(in)    :: stackentry
 integer,            intent(inout) :: istack
 integer,            intent(out)   :: nnode, mymum, level, npnode
 real,               intent(out)   :: xmini(3), xmaxi(3)

 nnode  = stackentry%node
 mymum  = stackentry%parent
 level  = stackentry%level
 npnode = stackentry%npnode
 xmini  = stackentry%xmin
 xmaxi  = stackentry%xmax
 istack = istack - 1

end subroutine pop_off_stack

subroutine compute_nodes_cofm(npnode,nnode,xyzcofm,totmass_node,doparallel)
 use dim,       only:maxtypes
 use part,      only:massoftype,igas,npartoftype
 integer, intent(in)  :: npnode,nnode
 real,    intent(out) :: xyzcofm(3),totmass_node
 logical, intent(in)  :: doparallel
 real    :: pmassi,fac,dfac
 real    :: xi,yi,zi,xcofm,ycofm,zcofm
 integer :: i1,i
!
! to avoid round off error from repeated multiplication by pmassi (which is small)
! we compute the centre of mass with a factor relative to gas particles
! but only if gas particles are present
!
 pmassi = massoftype(igas)
 fac    = 1.
 totmass_node = 0.
 if (pmassi > 0.) then
    dfac = 1./pmassi
 else
    pmassi = massoftype(maxloc(npartoftype(2:maxtypes),1)+1)
    if (pmassi > 0.) then
       dfac = 1./pmassi
    else
       dfac = 1.
    endif
 endif

 ! note that dfac can be a constant value across all particles even if APR is used
 i1=inoderange(1,nnode)
 xcofm = 0.
 ycofm = 0.
 zcofm = 0.

 ! during initial queue build which is serial, we can parallelise this loop
 if (npnode > 1000 .and. doparallel) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(npnode,dfac) &
    !$omp shared(treecache,i1) &
    !$omp private(i,xi,yi,zi) &
    !$omp firstprivate(pmassi,fac) &
    !$omp reduction(+:xcofm,ycofm,zcofm,totmass_node)
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       pmassi = treecache(5,i)
       fac    = pmassi*dfac ! to avoid round-off error
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       pmassi = treecache(5,i)
       fac    = pmassi*dfac ! to avoid round-off error
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
 endif

 xyzcofm = (/xcofm,ycofm,zcofm/)

 ! if there are no particles in this node, then the cofm will
 ! remain at zero
 if (totmass_node > 0.) then
    xyzcofm(:)   = xyzcofm(:)/(totmass_node*dfac)
 endif

end subroutine compute_nodes_cofm

subroutine set_nodes_properties(npnode,nnode,x0,totmass_node,mymum,nodeentry,level,doparallel)
 type(kdnode), intent(out)   :: nodeentry
 integer,      intent(in)    :: npnode,mymum,level,nnode
 real,         intent(inout) :: totmass_node
 real,         intent(in)    :: x0(3)
 logical,      intent(in)    :: doparallel
 real    :: pmassi
 real    :: dx,dy,dz,dr2,xi,yi,zi,hi
 real    :: hmax,r2max,totmass
 integer :: i1,i
#ifdef GRAVITY
 real    :: quads(6)
 real    :: octs(10)
#endif

 pmassi  = 0.
 totmass = 0.
 r2max = 0.
 hmax  = 0.
#ifdef GRAVITY
 quads(:) = 0.
 octs(:) = 0.
#endif

 i1=inoderange(1,nnode)

 !--compute size of node ! Parallel obsolete when build top tree is on but APR don't have it
 ! parallelise this loop if node is large enough
 ! use !$omp parallel do when doparallel=.true. (not in parallel region)
 ! when doparallel=.false., we're already in a parallel region but can't use nested reductions
 ! so we'll use thread-local accumulators and combine at the end
 if (npnode > 1000 .and. doparallel) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(npnode,treecache,x0,i1,use_geosplit) &
    !$omp private(i,xi,yi,zi,hi,dx,dy,dz,dr2) &
    !$omp firstprivate(pmassi) &
#ifdef GRAVITY
    !$omp reduction(+:totmass,quads,octs) &
#endif
    !$omp reduction(max:r2max,hmax)
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       hi = treecache(4,i)
       dx    = xi - x0(1)
       dy    = yi - x0(2)
       dz    = zi - x0(3)
       dr2   = dx*dx + dy*dy + dz*dz
       r2max = max(r2max,dr2)
       hmax  = max(hmax,hi)
#ifdef GRAVITY
       pmassi = treecache(5,i)
       totmass  = totmass  + pmassi
       call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       hi = treecache(4,i)
       dx    = xi - x0(1)
       dy    = yi - x0(2)
       dz    = zi - x0(3)
       dr2   = dx*dx + dy*dy + dz*dz
       r2max = max(r2max,dr2)
       hmax = max(hmax,hi)
#ifdef GRAVITY
       pmassi = treecache(5,i)
       totmass  = totmass  + pmassi
       call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
    enddo
 endif

 ! assign properties to node
 nodeentry%xcen       = x0(:)
 nodeentry%size       = sqrt(r2max) + epsilon(r2max)
 nodeentry%hmax       = hmax
 nodeentry%parent     = mymum
 nodeentry%level      = level
 nodeentry%mass       = totmass_node
#ifdef GRAVITY
 nodeentry%quads      = quads
 nodeentry%octs       = octs
#endif

end subroutine set_nodes_properties

#ifdef GRAVITY
!----------------------------------------------------------------
!+
!  Accumulate quadrupole and octupole moments of a particle
!  about the node centre of mass (extensive Cartesian form).
!+
!----------------------------------------------------------------
pure subroutine add_node_moments(pmassi,dx,dy,dz,quads,octs)
 real, intent(in)    :: pmassi,dx,dy,dz
 real, intent(inout) :: quads(6),octs(10)
 real :: dx2,dy2,dz2

 dx2 = dx*dx
 dy2 = dy*dy
 dz2 = dz*dz
 quads(1) = quads(1) + pmassi*dx2          ! Q_xx
 quads(2) = quads(2) + pmassi*dx*dy        ! Q_xy
 quads(3) = quads(3) + pmassi*dx*dz        ! Q_xz
 quads(4) = quads(4) + pmassi*dy2          ! Q_yy
 quads(5) = quads(5) + pmassi*dy*dz        ! Q_yz
 quads(6) = quads(6) + pmassi*dz2          ! Q_zz
 octs(1)  = octs(1)  + pmassi*dx2*dx       ! xxx
 octs(2)  = octs(2)  + pmassi*dx2*dy       ! xxy
 octs(3)  = octs(3)  + pmassi*dx2*dz       ! xxz
 octs(4)  = octs(4)  + pmassi*dx*dy2       ! xyy
 octs(5)  = octs(5)  + pmassi*dx*dy*dz     ! xyz
 octs(6)  = octs(6)  + pmassi*dx*dz2       ! xzz
 octs(7)  = octs(7)  + pmassi*dy2*dy       ! yyy
 octs(8)  = octs(8)  + pmassi*dy2*dz       ! yyz
 octs(9)  = octs(9)  + pmassi*dy*dz2       ! yzz
 octs(10) = octs(10) + pmassi*dz2*dz       ! zzz

end subroutine add_node_moments

!----------------------------------------------------------------
!+
!  Add the moments of a child node, about its own centre of mass,
!  to those of its parent, where dx is the offset of the child's
!  centre of mass from the parent's
!+
!----------------------------------------------------------------
pure subroutine add_child_moments(mass,quadsc,octsc,dx,quads,octs)
 real, intent(in)    :: mass,quadsc(6),octsc(10),dx(3)
 real, intent(inout) :: quads(6),octs(10)

 ! the child's own moments, plus those of its mass at its centre of mass
 quads = quads + quadsc
 octs  = octs  + octsc
 call add_node_moments(mass,dx(1),dx(2),dx(3),quads,octs)
 ! cross terms of the child's quadrupole with the offset
 octs(1)  = octs(1)  + 3.*quadsc(1)*dx(1)                                ! xxx
 octs(2)  = octs(2)  + quadsc(1)*dx(2) + 2.*quadsc(2)*dx(1)              ! xxy
 octs(3)  = octs(3)  + quadsc(1)*dx(3) + 2.*quadsc(3)*dx(1)              ! xxz
 octs(4)  = octs(4)  + quadsc(4)*dx(1) + 2.*quadsc(2)*dx(2)              ! xyy
 octs(5)  = octs(5)  + quadsc(2)*dx(3) + quadsc(3)*dx(2) + quadsc(5)*dx(1) ! xyz
 octs(6)  = octs(6)  + quadsc(6)*dx(1) + 2.*quadsc(3)*dx(3)              ! xzz
 octs(7)  = octs(7)  + 3.*quadsc(4)*dx(2)                                ! yyy
 octs(8)  = octs(8)  + quadsc(4)*dx(3) + 2.*quadsc(5)*dx(2)              ! yyz
 octs(9)  = octs(9)  + quadsc(6)*dx(2) + 2.*quadsc(5)*dx(3)              ! yzz
 octs(10) = octs(10) + 3.*quadsc(6)*dx(3)                                ! zzz

end subroutine add_child_moments
#endif

!--------------------------------------------------------------------
!+
!  create all the properties for a given node such as centre of mass,
!  size, max smoothing length, etc
!  will also split the node if necessary, setting wassplit=true
!  returns the left and right child information if split
!+
!--------------------------------------------------------------------
subroutine construct_node(nodeentry, nnode, mymum, level, xmini, xmaxi, npnode, doparallel,&
                          il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                          minlevel, maxlevel, wassplit, apr_tree, &
                          xyzmh_ptmass)
 use dim,       only:maxtypes,ind_timesteps
 use io,        only:fatal,error
 type(kdnode),    intent(out)   :: nodeentry
 integer,         intent(in)    :: nnode, mymum, level
 real,            intent(inout) :: xmini(3), xmaxi(3)
 integer,         intent(in)    :: npnode
 logical,         intent(in)    :: doparallel
 integer,         intent(out)   :: il, ir, nl, nr
 real,            intent(out)   :: xminl(3), xmaxl(3), xminr(3), xmaxr(3)
 integer(kind=8), intent(inout) :: ncells
 integer,         intent(out)   :: leaf_is_active(:)
 integer,         intent(inout) :: maxlevel, minlevel
 logical,         intent(out)   :: wassplit
 logical,         intent(in)    :: apr_tree
 real,            intent(in), optional :: xyzmh_ptmass(:,:)

 integer(kind=8) :: myslot
 real    :: xyzcofm(3)
 real    :: totmass_node
 logical :: nodeisactive
 integer :: i,npcounter,ipart
 real    :: x0(3)
 integer :: iaxis
 real    :: xpivot

 nodeisactive = .false.
 if (inoderange(1,nnode) > 0) then
    checkactive: do i = inoderange(1,nnode),inoderange(2,nnode)
       if (inodeparts(i) > 0) then
          nodeisactive = .true.
          exit checkactive
       endif
    enddo checkactive
    npcounter = inoderange(2,nnode) - inoderange(1,nnode) + 1
 else
    npcounter = 0
 endif

 if (npcounter /= npnode) then
    print*,'constructing node ',nnode,': found ',npcounter,' particles, expected:',npnode,' particles for this node'
    call fatal('maketree', 'expected number of particles in node differed from actual number')
 endif

 ! following lines to avoid compiler warnings on intent(out) variables
 ir = 0
 il = 0
 nl = 0
 nr = 0

 wassplit    = (npnode > minpart)

 xminl = 0.
 xmaxl = 0.
 xminr = 0.
 xmaxr = 0.

 if (npnode  <  1) return ! node has no particles, just quit

 xyzcofm(:) = 0.

 call compute_nodes_cofm(npnode,nnode,xyzcofm,totmass_node,doparallel)
 if (totmass_node<=0. .and. use_apr) call fatal('mtree + apr', &
    'totmass_node==0, something almost certainly wrong with aprmassoftype')
 if (totmass_node<=0.) call fatal('mtree','totmass_node==0',val=totmass_node)

 call set_nodes_properties(npnode,nnode,xyzcofm,totmass_node,mymum,nodeentry,level,doparallel)

 ! for gravity KDtree, we need the geo centre to split the node
 if (use_geosplit) then
    x0 = (xmaxi+xmini)*0.5
 else  !--for default KDtree, we need the split centre to be the centre of mass
    x0 = xyzcofm
 endif

 if (apr_tree)   wassplit = (npnode > 2)

 if (.not. wassplit) then
    nodeentry%leftchild  = 0
    nodeentry%rightchild = 0
    maxlevel = max(level,maxlevel)
    minlevel = min(level,minlevel)

    if (maxlevel > maxdepth) call fatal('maketree','maximum tree depth reached !!')
    ! individual timesteps where we mark leaf node as active/inactive
    if (ind_timesteps) then
       !
       !--mark leaf node as active (contains some active particles)
       !  or inactive by setting the firstincell to +ve (active) or -ve (inactive)
       !
       if (nodeisactive) then
          leaf_is_active(nnode) = 1
       else
          leaf_is_active(nnode) = -1
       endif
    else
       leaf_is_active(nnode) = 1
    endif
 else ! split this node and add children to stack
    iaxis  = maxloc(xmaxi - xmini,1) ! split along longest axis
    xpivot = x0(iaxis)

    maxlevel = max(level,maxlevel)
    if (maxlevel > maxdepth) call fatal('maketree','maximum tree depth reached !!')
    ! create two children nodes and point to them from current node
    if (level < maxlevel_indexed) then
       il = 2*nnode   ! indexing as per Gafton & Rosswog (2011)
       ir = il + 1
    else
       ! no need to lock, we could just atomic the update
       !$omp atomic capture
       ncells = ncells + 2
       myslot = ncells
       !$omp end atomic
       ir = int(myslot)
       il = int(myslot-1)
       if (ir > ncellsmax) call fatal('maketree',&
          'number of nodes exceeds array dimensions, increase ncellsmax and recompile',ival=int(ncellsmax))
    endif
    nodeentry%leftchild  = il
    nodeentry%rightchild = ir

    leaf_is_active(nnode) = 0

    if (npnode > 0) then
       if (apr_tree) then
          ! apr special sort - only used for merging particles
          call special_sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                    inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,treecache,inodeparts,&
                                    npnode)
       else
          ! regular sort
          call sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                  inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,treecache,inodeparts)
       endif

       if (nr + nl  /=  npnode) then
          call error('maketree','number of left + right != parent while splitting (likely cause: NaNs in position arrays)')
       endif

       ! see if all the particles ended up in one node, if so, arbitrarily build 2 cells. This should never happen
       if ((nl==npnode) .or. (nr==npnode)) then
          ! no need to move particles because if they all ended up in one node,
          ! then they are still in the original order
          nl = npnode / 2
          inoderange(1,il) = inoderange(1,nnode)
          inoderange(2,il) = inoderange(1,nnode) + nl - 1
          inoderange(1,ir) = inoderange(1,nnode) + nl
          inoderange(2,ir) = inoderange(2,nnode)
          nr = npnode - nl
       endif

       ! compute min/max with explicit loops for better cache behavior
       xminl(:) =  huge(xminl)
       xmaxl(:) = -huge(xmaxl)
       do ipart=inoderange(1,il),inoderange(2,il)
          xminl(1) = min(xminl(1),treecache(1,ipart))
          xminl(2) = min(xminl(2),treecache(2,ipart))
          xminl(3) = min(xminl(3),treecache(3,ipart))
          xmaxl(1) = max(xmaxl(1),treecache(1,ipart))
          xmaxl(2) = max(xmaxl(2),treecache(2,ipart))
          xmaxl(3) = max(xmaxl(3),treecache(3,ipart))
       enddo

       xminr(:) =  huge(xminr)
       xmaxr(:) = -huge(xmaxr)
       do ipart=inoderange(1,ir),inoderange(2,ir)
          xminr(1) = min(xminr(1),treecache(1,ipart))
          xminr(2) = min(xminr(2),treecache(2,ipart))
          xminr(3) = min(xminr(3),treecache(3,ipart))
          xmaxr(1) = max(xmaxr(1),treecache(1,ipart))
          xmaxr(2) = max(xmaxr(2),treecache(2,ipart))
          xmaxr(3) = max(xmaxr(3),treecache(3,ipart))
       enddo
    else
       nl = 0
       nr = 0
       xminl =  huge(xminl)
       xmaxl = -huge(xmaxl)
       xminr =  huge(xminr)
       xmaxr = -huge(xmaxr)
    endif

 endif

end subroutine construct_node

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis
!+
!----------------------------------------------------------------
subroutine sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,nl,nr,xpivot,&
                                   treecache,inodeparts)
 integer, intent(in)    :: iaxis,imin,imax
 integer, intent(out)   :: min_l,max_l,min_r,max_r,nl,nr
 real,    intent(inout) :: xpivot,treecache(:,:)
 integer, intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot
 integer :: inodeparts_swap,i,j
 real :: xyzh_swap(5)
 real :: xi_coord, xj_coord

 !print*,'nnode ',imin,imax,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 xi_coord = treecache(iaxis,i)
 xj_coord = treecache(iaxis,j)
 i_lt_pivot = xi_coord <= xpivot
 j_lt_pivot = xj_coord <= xpivot
 !  k = 0

 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       xi_coord = treecache(iaxis,i)
       i_lt_pivot = xi_coord <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          xj_coord = treecache(iaxis,j)
          j_lt_pivot = xj_coord <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:5)  = treecache(1:5,i)

          inodeparts(i)   = inodeparts(j)
          treecache(1:5,i) = treecache(1:5,j)

          inodeparts(j)   = inodeparts_swap
          treecache(1:5,j) = xyzh_swap(1:5)

          i = i + 1
          j = j - 1
          xi_coord = treecache(iaxis,i)
          xj_coord = treecache(iaxis,j)
          i_lt_pivot = xi_coord <= xpivot
          j_lt_pivot = xj_coord <= xpivot
          ! k = k + 1
       endif
    endif
 enddo
 if (.not.i_lt_pivot) i = i - 1
 if (j_lt_pivot)      j = j + 1

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine sort_particles_in_cell

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis, but additionally
!  force the cells to have a certain minimum number of particles per cell
!+
!----------------------------------------------------------------
subroutine special_sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,&
                                nl,nr,xpivot,treecache,inodeparts,npnode)
 use io, only:error
 integer, intent(in)    :: iaxis,imin,imax,npnode
 integer, intent(out)   :: min_l,max_l,min_r,max_r,nl,nr
 real,    intent(inout) :: xpivot,treecache(:,:)
 integer, intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot,slide_l,slide_r
 integer :: inodeparts_swap,i,j,nchild_in
 integer :: k,ii,rem_nr,rem_nl
 real :: xyzh_swap(5),dpivot(npnode)

 dpivot = 0.0
 nchild_in = 2

 if (modulo(npnode,nchild_in) > 0) then
    call error('apr sort','number of particles sent in to kdtree is not divisible by 2')
 endif

! print*,'nnode ',imin,imax,npnode,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 i_lt_pivot = treecache(iaxis,i) <= xpivot
 j_lt_pivot = treecache(iaxis,j) <= xpivot
 dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
 !k = 0
 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
       i_lt_pivot = treecache(iaxis,i) <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
          j_lt_pivot = treecache(iaxis,j) <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:5)  = treecache(1:5,i)

          inodeparts(i)   = inodeparts(j)
          treecache(1:5,i) = treecache(1:5,j)

          inodeparts(j)   = inodeparts_swap
          treecache(1:5,j) = xyzh_swap(1:5)

          i = i + 1
          j = j - 1

          dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
          dpivot(j-imin+1) = xpivot - treecache(iaxis,j)

          i_lt_pivot = treecache(iaxis,i) <= xpivot
          j_lt_pivot = treecache(iaxis,j) <= xpivot
       endif
    endif
 enddo

 ! all particles may lie on one side of the pivot, in which case i or j
 ! steps outside the node and there is no distance to record
 if (.not.i_lt_pivot) then
    i = i - 1
    if (i >= imin) dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 endif
 if (j_lt_pivot) then
    j = j + 1
    if (j <= imax) dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
 endif

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

 ! does the pivot need to be adjusted?
 rem_nl = modulo(nl,nchild_in)
 rem_nr = modulo(nr,nchild_in)
 if (rem_nl == 0 .and. rem_nr == 0) return

 ! Decide which direction the pivot needs to go
 if (rem_nl < rem_nr) then
    slide_l = .true.
    slide_r = .false.
 else
    slide_l = .false.
    slide_r = .true.
 endif
 ! Override this if there's less than nchild*2 in the cell
 if (nl < nchild_in) then
    slide_r = .true.
    slide_l = .false.
 elseif (nr < nchild_in) then
    slide_r = .false.
    slide_l = .true.
 endif

 ! Move across particles by distance from xpivot till we get
 ! the right number of particles in each cell
 if (slide_r) then
    do ii = 1,rem_nr
       ! next particle to shift across
       k = minloc(dpivot,dim=1,mask=dpivot > 0.) + imin - 1
       if (k-imin+1==0) k = maxloc(dpivot,dim=1,mask=dpivot < 0.) + imin - 1

       ! swap this with the first particle on the j side
       inodeparts_swap = inodeparts(k)
       xyzh_swap(1:5)  = treecache(1:5,k)

       inodeparts(k)   = inodeparts(j)
       treecache(1:5,k) = treecache(1:5,j)

       inodeparts(j)   = inodeparts_swap
       treecache(1:5,j) = xyzh_swap(1:5)

       ! and now shift to the right
       i = i + 1
       j = j + 1

       ! ditch it, go again
       dpivot(k-imin+1) = huge(k-imin+1)
    enddo
 else
    do ii = 1,rem_nl
       ! next particle to shift across
       k = maxloc(dpivot,dim=1,mask=dpivot < 0.) + imin - 1
       if (k-imin+1==0) k = minloc(dpivot,dim=1,mask=dpivot > 0.) + imin - 1

       ! swap this with the last particle on the i side
       inodeparts_swap = inodeparts(k)
       xyzh_swap(1:5)  = treecache(1:5,k)

       inodeparts(k)   = inodeparts(i)
       treecache(1:5,k) = treecache(1:5,i)

       inodeparts(i)   = inodeparts_swap
       treecache(1:5,i) = xyzh_swap(1:5)

       ! and now shift to the left
       i = i - 1
       j = j - 1

       ! ditch it, go again
       dpivot(k-imin+1) = huge(k-imin+1)

    enddo
 endif

 ! tidy up outputs
 max_l = i
 min_r = j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine special_sort_particles_in_cell

!--------------------------------------------------------------------------------
!+
!  Revise the tree for new particle positions without rebuilding it.
!
!  The particles stay in the leaves of the last build. The leaves are
!  recomputed from their particles, and each parent from its two children,
!  working up the tree a level at a time. A parent's size is then set from
!  its particles, as when built, since the bound from its children would
!  open the node more often than needed.
!
!  Given the numbers of active and live particles (over all MPI tasks), the
!  tree is only revised on substeps with few active particles, until the
!  particles have had about one step each since it was built.
!
!  Returns ok=.false., and the tree must be rebuilt, if it is not to be
!  revised, if the last build was not of these particles, or if particles
!  have been added, accreted or killed since. With more than one MPI task,
!  every task revises its tree or none does.
!+
!--------------------------------------------------------------------------------
subroutine revtree(node,xyzh,np,leaf_is_active,ncells,ok,nactivetot,nlivetot)
 use dim,      only:ind_timesteps,periodic,mpi
 use part,     only:iphase,iactive,isdead_or_accreted,maxphase,igas,massoftype,iamtype
 use mpiutils, only:reduceall_mpi
 type(kdnode),    intent(inout) :: node(:)
 real,            intent(in)    :: xyzh(:,:)
 integer,         intent(in)    :: np
 integer,         intent(inout) :: leaf_is_active(:)
 integer(kind=8), intent(in)    :: ncells
 logical,         intent(out)   :: ok
 real,            intent(in), optional :: nactivetot,nlivetot  ! active and live particles, all tasks
 integer :: nproot,nlive,nbad,ip,i

 ok = .false.
 if (.not.tree_is_revisable .or. use_apr .or. use_sinktree .or. periodic) return
 if (present(nactivetot) .and. present(nlivetot)) then
    nactive_since_build = nactive_since_build + nactivetot
    if (nactivetot >= frac_active_revtree*nlivetot .or. nactive_since_build >= nlivetot) return
 endif
 !
 ! refresh the particle list: positions, masses and active flags
 !
 nproot = inoderange(2,irootnode)
 nlive = 0
 nbad   = 0
 if (inoderange(1,irootnode) /= 1 .or. nproot < 1) then
    nproot = 0
    nbad   = 1
 endif
 !$omp parallel default(none) &
 !$omp shared(np,nproot,xyzh,inodeparts,treecache,iphase,massoftype,maxphase,maxp) &
 !$omp private(ip,i) &
 !$omp reduction(+:nlive,nbad)
 !$omp do schedule(static)
 do i=1,np
    if (.not.isdead_or_accreted(xyzh(4,i))) nlive = nlive + 1
 enddo
 !$omp enddo nowait
 !$omp do schedule(static)
 do ip=1,nproot
    i = abs(inodeparts(ip))
    if (i < 1 .or. i > np) then
       nbad = nbad + 1
       cycle
    endif
    if (isdead_or_accreted(xyzh(4,i))) nbad = nbad + 1
    if (ind_timesteps) then
       if (iactive(iphase(i))) then
          inodeparts(ip) = i  ! +ve if active
       else
          inodeparts(ip) = -i ! -ve if inactive
       endif
    endif
    treecache(1:4,ip) = xyzh(1:4,i)
    if (maxphase==maxp) then
       treecache(5,ip) = massoftype(iamtype(iphase(i)))
    else
       treecache(5,ip) = massoftype(igas)
    endif
 enddo
 !$omp enddo
 !$omp end parallel
 if (nlive /= nproot) nbad = nbad + 1
 ! with more than one MPI task, every task revises its tree or all rebuild theirs
 if (mpi .and. nprocs > 1) nbad = int(reduceall_mpi('+',nbad))
 if (nbad > 0) return
 ok = .true.
 !
 ! recompute the nodes for the new positions, then their sizes from their particles
 !
 call revise_nodes(node,nproot,leaf_is_active)
 call set_exact_parent_sizes(node,nproot)

end subroutine revtree

!--------------------------------------------------------------------------------
!+
!  recompute every node of the tree from its particles (leaves) or its two
!  children (parents). The particles are in tree order, so each thread takes
!  a contiguous range of them and walks down the nodes overlapping it,
!  computing each leaf whose first particle is in its range and each parent
!  within its range once its children are done. The few parents straddling
!  the ranges of several threads are combined afterwards, deepest first
!+
!--------------------------------------------------------------------------------
subroutine revise_nodes(node,nproot,leaf_is_active)
!$ use omp_lib, only:omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: nproot
 integer,      intent(inout) :: leaf_is_active(:)
 integer, parameter   :: maxspan = 2*maxdepth+2
 integer, allocatable :: ispan(:,:),dspan(:,:),nspan(:),list(:)
 integer :: ia,ib,ithread,nthreads,maxthreads,d,k,t,n,j,ifirst(0:maxdepth+1)

 maxthreads = 1
!$ maxthreads = omp_get_max_threads()
 allocate(ispan(maxspan,maxthreads),dspan(maxspan,maxthreads),nspan(maxthreads))
 nspan = 0

 !$omp parallel default(none) &
 !$omp shared(node,nproot,leaf_is_active,ispan,dspan,nspan) &
 !$omp private(ia,ib,ithread,nthreads)
 ithread  = 0
 nthreads = 1
!$ ithread  = omp_get_thread_num()
!$ nthreads = omp_get_num_threads()
 ia = int(int(ithread,kind=8)*nproot/nthreads) + 1
 ib = int(int(ithread+1,kind=8)*nproot/nthreads)
 if (ia <= ib) call revise_walk(node,irootnode,0,ia,ib,leaf_is_active, &
                                ispan(:,ithread+1),dspan(:,ithread+1),nspan(ithread+1))
 !$omp end parallel

 ! the straddling parents, grouped by depth (each may be listed by several threads)
 ifirst = 0
 do t = 1,maxthreads
    do k = 1,nspan(t)
       ifirst(dspan(k,t)+1) = ifirst(dspan(k,t)+1) + 1
    enddo
 enddo
 do d = 1,maxdepth+1
    ifirst(d) = ifirst(d) + ifirst(d-1)
 enddo
 allocate(list(max(ifirst(maxdepth+1),1)))
 do t = 1,maxthreads
    do k = 1,nspan(t)
       d = dspan(k,t)
       ifirst(d) = ifirst(d) + 1
       list(ifirst(d)) = ispan(k,t)
    enddo
 enddo
 ! ifirst(d) is now the last entry of depth d: combine them deepest first, once each
 do d = maxdepth,0,-1
    k = 1
    if (d > 0) k = ifirst(d-1) + 1
    do j = k,ifirst(d)
       n = list(j)
       if (any(list(k:j-1) == n)) cycle
       call combine_child_nodes(node(node(n)%leftchild),node(node(n)%rightchild),node(n))
    enddo
 enddo

end subroutine revise_nodes

!--------------------------------------------------------------------------------
!+
!  the walk of revise_nodes below node n, at depth d, over the particles ia to
!  ib. Parents straddling the range are listed in ispan, with their depth in dspan
!+
!--------------------------------------------------------------------------------
recursive subroutine revise_walk(node,n,d,ia,ib,leaf_is_active,ispan,dspan,nspan)
 use dim, only:ind_timesteps
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: n,d,ia,ib
 integer,      intent(inout) :: leaf_is_active(:)
 integer,      intent(inout) :: ispan(:),dspan(:),nspan
 integer :: np,ip,mymum
 real    :: x0(3),totmass_node
 logical :: nodeisactive

 if (max(ia,inoderange(1,n)) > min(ib,inoderange(2,n))) return

 if (node(n)%leftchild == 0) then
    ! a leaf, computed by the thread holding its first particle
    if (inoderange(1,n) < ia) return
    np = inoderange(2,n) - inoderange(1,n) + 1
    if (inoderange(1,n) <= 0 .or. np < 1) return
    mymum = node(n)%parent
    call compute_nodes_cofm(np,n,x0,totmass_node,.false.)
    call set_nodes_properties(np,n,x0,totmass_node,mymum,node(n),d,.false.)
    node(n)%leftchild  = 0
    node(n)%rightchild = 0
    node(n)%mass       = totmass_node
    if (ind_timesteps) then
       nodeisactive = .false.
       do ip=inoderange(1,n),inoderange(2,n)
          if (inodeparts(ip) > 0) then
             nodeisactive = .true.
             exit
          endif
       enddo
       if (nodeisactive) then
          leaf_is_active(n) = 1
       else
          leaf_is_active(n) = -1
       endif
    else
       leaf_is_active(n) = 1
    endif
    return
 endif

 call revise_walk(node,node(n)%leftchild,d+1,ia,ib,leaf_is_active,ispan,dspan,nspan)
 call revise_walk(node,node(n)%rightchild,d+1,ia,ib,leaf_is_active,ispan,dspan,nspan)
 if (inoderange(1,n) >= ia .and. inoderange(2,n) <= ib) then
    call combine_child_nodes(node(node(n)%leftchild),node(node(n)%rightchild),node(n))
 else
    nspan = nspan + 1
    ispan(nspan) = n
    dspan(nspan) = d
 endif

end subroutine revise_walk

!--------------------------------------------------------------------------------
!+
!  the size of each parent from its particles: the largest distance from
!  its centre, as set when the tree is built. The particles are in tree
!  order, so each thread takes a contiguous range of them and walks down
!  the nodes overlapping it, comparing each particle with all of its
!  ancestors while it is in cache: every particle is read once, by the
!  thread whose memory it is in. Only the nodes that straddle the ranges
!  of several threads, at most a few per level, combine their maxima
!+
!--------------------------------------------------------------------------------
subroutine set_exact_parent_sizes(node,nproot)
!$ use omp_lib, only:omp_get_thread_num,omp_get_num_threads
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: nproot
 real,    allocatable :: r2span(:)
 integer :: ia,ib,ithread,nthreads,k,nspan
 integer :: ispan(2*maxdepth+2)
 real    :: r2(0:maxdepth),xc(3,0:maxdepth),r2s(2*maxdepth+2)

 ! the furthest particle from each straddling node, over the threads
 allocate(r2span(size(node)))

 !$omp parallel default(none) &
 !$omp shared(node,nproot,r2span) &
 !$omp private(ia,ib,ithread,nthreads,k,nspan,ispan,r2,xc,r2s)
 ithread  = 0
 nthreads = 1
!$ ithread  = omp_get_thread_num()
!$ nthreads = omp_get_num_threads()
 ia = int(int(ithread,kind=8)*nproot/nthreads) + 1
 ib = int(int(ithread+1,kind=8)*nproot/nthreads)
 nspan = 0
 if (ia <= ib) call exact_size_walk(node,irootnode,0,ia,ib,r2,xc,ispan,r2s,nspan)

 ! the straddling nodes: set to zero, then the maximum over the threads holding them
 do k = 1,nspan
    r2span(ispan(k)) = 0.
 enddo
 !$omp barrier
 do k = 1,nspan
    !$omp atomic
    r2span(ispan(k)) = max(r2span(ispan(k)),r2s(k))
 enddo
 !$omp barrier
 do k = 1,nspan
    node(ispan(k))%size = sqrt(r2span(ispan(k))) + epsilon(r2span)
 enddo
 !$omp end parallel

end subroutine set_exact_parent_sizes

!--------------------------------------------------------------------------------
!+
!  the walk of set_exact_parent_sizes below node n, at depth d, over the
!  particles ia to ib: r2 and xc hold the running maxima and the centres of
!  the ancestors of n. A parent within the range gets its size here; one
!  straddling the range is listed in ispan, with its maximum over the range in r2s
!+
!--------------------------------------------------------------------------------
recursive subroutine exact_size_walk(node,n,d,ia,ib,r2,xc,ispan,r2s,nspan)
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: n,d,ia,ib
 real,         intent(inout) :: r2(0:maxdepth),xc(3,0:maxdepth)
 integer,      intent(inout) :: ispan(:),nspan
 real,         intent(inout) :: r2s(:)
 integer :: ip,k,i1,i2
 real    :: dx(3)

 i1 = max(ia,inoderange(1,n))
 i2 = min(ib,inoderange(2,n))
 if (i1 > i2) return

 if (node(n)%leftchild <= 0) then
    ! a leaf: each of its particles against every ancestor
    do ip = i1,i2
       do k = 0,d-1
          dx = treecache(1:3,ip) - xc(:,k)
          r2(k) = max(r2(k),dot_product(dx,dx))
       enddo
    enddo
    return
 endif

 r2(d) = 0.
 xc(:,d) = node(n)%xcen
 call exact_size_walk(node,node(n)%leftchild,d+1,ia,ib,r2,xc,ispan,r2s,nspan)
 call exact_size_walk(node,node(n)%rightchild,d+1,ia,ib,r2,xc,ispan,r2s,nspan)
 if (inoderange(1,n) >= ia .and. inoderange(2,n) <= ib) then
    node(n)%size = sqrt(r2(d)) + epsilon(r2)
 else
    ! straddles the range (an ancestor of either end): combined over the threads later
    nspan = nspan + 1
    ispan(nspan) = n
    r2s(nspan) = r2(d)
 endif

end subroutine exact_size_walk

!--------------------------------------------------------------------------------
!+
!  the properties of a parent node from those of its two children: the
!  centre of mass, a sphere bounding theirs, the largest smoothing length
!  and, with self-gravity, the multipole moments about the centre of mass.
!  The links and level of the parent are left unchanged
!+
!--------------------------------------------------------------------------------
pure subroutine combine_child_nodes(left,right,parent)
 type(kdnode), intent(in)    :: left,right
 type(kdnode), intent(inout) :: parent
 real :: ml,mr,x0(3),dxl(3),dxr(3)

 ml = left%mass
 mr = right%mass
 if (ml + mr > 0.) then
    x0 = (ml*left%xcen + mr*right%xcen)/(ml + mr)
 else
    x0 = 0.5*(left%xcen + right%xcen)
 endif
 dxl = left%xcen - x0
 dxr = right%xcen - x0
 parent%xcen = x0
 parent%size = max(sqrt(dot_product(dxl,dxl)) + left%size, &
                   sqrt(dot_product(dxr,dxr)) + right%size)
 parent%hmax = max(left%hmax,right%hmax)
 parent%mass = ml + mr
#ifdef GRAVITY
 parent%quads = 0.
 parent%octs  = 0.
 call add_child_moments(ml,left%quads,left%octs,dxl,parent%quads,parent%octs)
 call add_child_moments(mr,right%quads,right%octs,dxr,parent%quads,parent%octs)
#endif

end subroutine combine_child_nodes

!-----------------------------------------------------------------------
!+
!  push a pair onto a stack that grows as needed
!+
!-----------------------------------------------------------------------
subroutine push_pair(istack,nstack,i1,i2,i3)
 integer, allocatable, intent(inout) :: istack(:,:)
 integer,              intent(inout) :: nstack
 integer,              intent(in)    :: i1,i2,i3
 integer, allocatable :: tmp(:,:)

 if (nstack >= size(istack,2)) then
    allocate(tmp(3,2*size(istack,2)))
    tmp(:,1:nstack) = istack(:,1:nstack)
    call move_alloc(tmp,istack)
 endif
 nstack = nstack + 1
 istack(:,nstack) = [i1,i2,i3]

end subroutine push_pair

end module kdtree
