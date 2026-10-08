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
! :Dependencies: allocutils, boundary, dim, dtypekdtree, io, kernel,
!   mpibalance, mpidomain, mpitree, mpiutils, part, timing
!
 use dim,         only:maxp,ncellsmax,minpart,use_apr,use_sinktree,maxptmass,maxpsph,gravity
 use io,          only:nprocs
 use dtypekdtree, only:kdnode,lenfgrav
 use part,        only:ll,iphase,treecache,maxphase, &
                       apr_level,aprmassoftype

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)
 type(kdnode),     allocatable :: refinementnode(:)
!
!--tree parameters
!
 integer,          parameter, public :: irootnode    = 1
 character(len=1), parameter, public :: labelax(3)   = (/'x','y','z'/)
 integer,          parameter, public :: maxdepth     = 64
!
!--runtime options for this module
!
 real,    public  :: tree_accuracy    = 0.5
 logical, public  :: use_geosplit     = gravity ! only debug flag / should be on with gravity
 logical, private :: done_init_kdtree = .false.
 logical, private :: already_warned   = .false.
 integer, private :: numthreads
 ! scratch space for the parallel partition in build_top_parallel
 real,    allocatable, private :: tcbuf(:,:)
 integer, allocatable, private :: ipbuf(:)

! Index of the last node in the local tree that has been copied to
! the global tree
 integer, public :: irefine

 public :: allocate_kdtree, deallocate_kdtree
 public :: maketree, revtree,kdnode
 public :: maketreeglobal
 public :: empty_tree
 integer, public :: maxlevel_indexed, maxlevel

 ! neighbour cache indices (xyzcache); imported with only: from dens/force
 integer, parameter, public :: ix=1, iy=2, iz=3, ih1=4, im=5, irho=6, izetaomega=7, isoftomega=8

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
 use dim, only:mpi
 use allocutils, only:allocate_array

 call allocate_array('inoderange', inoderange, 2, ncellsmax+1)
 call allocate_array('inodeparts', inodeparts, maxp)
 if (mpi) call allocate_array('refinementnode', refinementnode, ncellsmax+1)

end subroutine allocate_kdtree

subroutine deallocate_kdtree
 use dim, only:mpi
 if (allocated(inoderange)) deallocate(inoderange)
 if (allocated(inodeparts)) deallocate(inodeparts)
 if (mpi .and. allocated(refinementnode)) deallocate(refinementnode)
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
subroutine maketree(node, xyzh, np, leaf_is_active, ncells, apr_tree, refinelevels,nptmass,xyzmh_ptmass)
 use io,   only:fatal,warning,iprint,iverbose
!$ use omp_lib
 type(kdnode),    intent(out)   :: node(:) !ncellsmax+1)
 integer,         intent(in)    :: np
 real,            intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,         intent(out)   :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 logical,         intent(in)    :: apr_tree
 integer,         intent(out),   optional :: refinelevels
 integer,         intent(in),    optional :: nptmass
 real,            intent(inout), optional :: xyzmh_ptmass(:,:)

 integer :: i,npnode,il,ir,istack,nl,nr,mymum
 integer :: nnode,minlevel,level,nqueue
 real :: xmini(3),xmaxi(3),xminl(3),xmaxl(3),xminr(3),xmaxr(3)
 integer, parameter :: istacksize = 512
 type(kdbuildstack), save :: stack(istacksize)
 !$omp threadprivate(stack)
 type(kdbuildstack) :: queue(istacksize)
!$ integer :: threadid
 integer :: npcounter
 logical :: wassplit,finished,sinktree
 character(len=10) :: string

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
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,xyzmh_ptmass,nptmass)
 else
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh)
 endif

 if (inoderange(1,irootnode)==0 .or. inoderange(2,irootnode)==0 ) then
    call fatal('maketree','no particles or all particles dead/accreted')
 endif

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

 nqueue = numthreads
 ! build the first levels with all threads working on every level (not for the APR
 ! merge tree, whose partition has to leave an even number of particles in each child)
 if (.not.apr_tree .and. nqueue > 1 .and. npcounter > max(minpart,64)) then
    call build_top_parallel(node,queue,istack,nqueue,leaf_is_active)
 endif
 ! build using a queue to build level by level until number of nodes = number of threads
 over_queue: do while (istack  <  nqueue)
    ! if the tree finished while building the queue, then we should just return
    ! only happens for small particle numbers
    if (istack <= 0) then
       finished = .true.
       exit over_queue
    endif
    ! pop off front of queue
    call pop_off_stack(queue(1), istack, nnode, mymum, level, npnode, xmini, xmaxi)

    ! shuffle queue forward
    do i=1,istack
       queue(i) = queue(i+1)
    enddo

    ! construct node
    if (sinktree) then
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, .false., apr_tree, xyzmh_ptmass)
    else
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, .false., apr_tree)
    endif

    if (wassplit) then ! add children to back of queue
       if (istack+2 > istacksize) call fatal('maketree',&
                                       'queue size exceeded in tree build, increase istacksize and recompile')

       istack = istack + 1
       call push_onto_stack(queue(istack),il,nnode,level+1,nl,xminl,xmaxl)
       istack = istack + 1
       call push_onto_stack(queue(istack),ir,nnode,level+1,nr,xminr,xmaxr)
    endif

 enddo over_queue
 ! build_top_parallel completes whole levels, so the queue may hold more than nqueue nodes
 if (.not.finished) nqueue = istack

 ! fix the indices

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
                                 minlevel, maxlevel, wassplit, .false., apr_tree, xyzmh_ptmass)
          else
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, .false., apr_tree)
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

 ! decrease number of cells if tree is entirely within 2^k indexing limit
 if (maxlevel < maxlevel_indexed) then
    ncells = 2**(maxlevel+1) - 1
 endif

 if (maxlevel > maxlevel_indexed .and. .not.already_warned) then
    write(string,"(i10)") 2**(maxlevel-maxlevel_indexed)
    if (iverbose > 0) call warning('maketree','maxlevel > max_indexed: will run faster if recompiled with '// &
               'NCELLSMAX='//trim(adjustl(string))//'*maxp,')
 endif

 if (present(refinelevels)) refinelevels = minlevel

 if (iverbose >= 3) then
    write(iprint,"(a,i10,3(a,i2))") ' maketree: nodes = ',ncells,', max level = ',maxlevel,&
       ', min leaf level = ',minlevel,' max level indexed = ',maxlevel_indexed
 endif

end subroutine maketree

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
subroutine build_top_parallel(node,queue,istack,nqueue,leaf_is_active)
 use part,      only:massoftype,igas,npartoftype
 use dim,       only:maxtypes
 use io,        only:fatal
 type(kdnode),       intent(inout) :: node(:)
 type(kdbuildstack), intent(inout) :: queue(:)
 integer,            intent(inout) :: istack
 integer,            intent(in)    :: nqueue
 integer,            intent(inout) :: leaf_is_active(:)
 integer, allocatable :: cnode(:),chunkl(:),chunkr(:),cnl(:),coffl(:),coffr(:)
 integer, allocatable :: nodecl(:),nodecr(:),jnl(:),nodeax(:),chunk(:)
 logical, allocatable :: jdegen(:)
 real,    allocatable :: psum(:,:),phm(:),pr2(:),pbox(:,:),nodecom(:,:),nodecog(:,:),remain(:)
#ifdef GRAVITY
 real,    allocatable :: pmom(:,:)
 real    :: quads(6),octs(10)
#endif
 type(kdbuildstack), allocatable :: newq(:)
 real,    allocatable :: tcswap(:,:)
 integer, allocatable :: ipswap(:)
 integer :: kmax,cmax,k,j,c,m,ncj,nchunk,i,i1,n,nl,pl,pr,iax,nnode,il,ir,isplit
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
 allocate(nodecl(kmax),nodecr(kmax),jnl(kmax),jdegen(kmax))
 allocate(nodeax(kmax),nodecom(3,kmax),nodecog(3,kmax))
 allocate(chunk(kmax),remain(kmax),newq(2*kmax))
 allocate(cnode(cmax),chunkl(cmax),chunkr(cmax),cnl(cmax),coffl(cmax),coffr(cmax))
 allocate(psum(4,cmax),phm(cmax),pr2(cmax),pbox(12,cmax))
#ifdef GRAVITY
 allocate(pmom(16,cmax))
#endif

 levels: do while (istack < nqueue)
    k = istack
    !
    ! the queue holds one whole level (all nodes have the same level); split them all
    ! together, or hand over to the serial loop if the next level would not fit in the
    ! queue, would pass the levels numbered 2n/2n+1, or any node is small
    !
    if (k > kmax .or. queue(1)%level >= maxlevel_indexed) exit levels
    if (any(queue(1:k)%npnode <= max(minpart,64))) exit levels
    !
    ! each node gets its share of the threads as chunks: the share rounded down (at
    ! least one), then the chunks left over go to the nodes that lost most in the
    ! rounding (largest remainder), so a level with no more nodes than threads has
    ! exactly one chunk per thread.  A chunk never crosses a node boundary.
    !
    ntot = sum(int(queue(1:k)%npnode,kind=8))
    do j = 1,k
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

    !$omp parallel default(shared) &
#ifdef GRAVITY
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
       totmass = sum(psum(4,nodecl(j):nodecr(j)))
       if (totmass <= 0.) call fatal('mtree','totmass_node==0',val=totmass)
       nodecom(:,j) = sum(psum(1:3,nodecl(j):nodecr(j)),dim=2)/(totmass*dfac)
#ifdef GRAVITY
       node(queue(j)%node)%mass = totmass
#endif
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
    ! the children, as the next level of the queue
    !$omp do schedule(static)
    do j = 1,k
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
       call push_onto_stack(newq(2*j-1),il,nnode,queue(j)%level+1,nl,bl(1:3),bl(4:6))
       call push_onto_stack(newq(2*j),ir,nnode,queue(j)%level+1,n-nl,br(1:3),br(4:6))
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

    queue(1:2*k) = newq(1:2*k)
    istack = 2*k
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

!$omp parallel do private(i)
 do i=1,size(node)
    node(i)%xcen = 0.
    node(i)%size = 0.
    node(i)%hmax = 0.
    node(i)%leftchild = 0
    node(i)%rightchild = 0
    node(i)%parent = 0
    node(i)%level  = 0
#ifdef GRAVITY
    node(i)%mass  = 0.
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
subroutine construct_root_node(np,nproot,irootnode,xmini,xmaxi,leaf_is_active,xyzh,xyzmh_ptmass,nptmass)
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
 real,    intent(inout), optional :: xyzmh_ptmass(:,:)
 integer, intent(in),    optional :: nptmass
 integer :: i,ncross,ic,nchunk,nl
 integer, allocatable :: nlive(:)
 real    :: xminpart,yminpart,zminpart,xmaxpart,ymaxpart,zmaxpart
 real    :: xi, yi, zi

 xminpart = xyzh(1,1)
 yminpart = xyzh(2,1)
 zminpart = xyzh(3,1)
 xmaxpart = xminpart
 ymaxpart = yminpart
 zmaxpart = zminpart

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
          if (isnan(xi) .or. isnan(yi) .or. isnan(zi)) then
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
             if (isnan(xi) .or. isnan(yi) .or. isnan(zi)) then
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
 !$omp shared(maxp,maxphase) &
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
             if (use_apr) inodeparts(nproot) = abs(inodeparts(nproot))
          else
             inodeparts(nproot) = i
          endif
          treecache(1:4,nproot) = xyzh(1:4,i)
          if (maxphase==maxp) then
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

subroutine set_nodes_properties(npnode,nnode,x0,totmass_node,mymum,nodeentry,xmini,xmaxi,&
                                level,global_build,doparallel)
 use mpitree,   only:reduce_group
 use dim,       only:mpi
 type(kdnode),    intent(out)   :: nodeentry
 integer,         intent(in)    :: npnode,mymum,level,nnode
 real,            intent(inout) :: xmini(3), xmaxi(3), totmass_node
 real,            intent(in)    :: x0(3)
 logical,         intent(in)    :: doparallel,global_build
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

 ! reduce node limits and quads across MPI tasks belonging to this group
 if (mpi .and. global_build) then
    r2max     = reduce_group(r2max,'max',level)
    hmax      = reduce_group(hmax,'max',level)

    xmini(1)  = reduce_group(xmini(1),'min',level)
    xmini(2)  = reduce_group(xmini(2),'min',level)
    xmini(3)  = reduce_group(xmini(3),'min',level)

    xmaxi(1)  = reduce_group(xmaxi(1),'max',level)
    xmaxi(2)  = reduce_group(xmaxi(2),'max',level)
    xmaxi(3)  = reduce_group(xmaxi(3),'max',level)
#ifdef GRAVITY
    quads(1)  = reduce_group(quads(1),'+',level)
    quads(2)  = reduce_group(quads(2),'+',level)
    quads(3)  = reduce_group(quads(3),'+',level)
    quads(4)  = reduce_group(quads(4),'+',level)
    quads(5)  = reduce_group(quads(5),'+',level)
    quads(6)  = reduce_group(quads(6),'+',level)
    octs(1)   = reduce_group(octs(1),'+',level)
    octs(2)   = reduce_group(octs(2),'+',level)
    octs(3)   = reduce_group(octs(3),'+',level)
    octs(4)   = reduce_group(octs(4),'+',level)
    octs(5)   = reduce_group(octs(5),'+',level)
    octs(6)   = reduce_group(octs(6),'+',level)
    octs(7)   = reduce_group(octs(7),'+',level)
    octs(8)   = reduce_group(octs(8),'+',level)
    octs(9)   = reduce_group(octs(9),'+',level)
    octs(10)  = reduce_group(octs(10),'+',level)
#endif
 endif

 ! assign properties to node
 nodeentry%xcen       = x0(:)
 nodeentry%size       = sqrt(r2max) + epsilon(r2max)
 nodeentry%hmax       = hmax
 nodeentry%parent     = mymum
 nodeentry%level      = level
#ifdef GRAVITY
 nodeentry%mass       = totmass_node
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
                          minlevel, maxlevel, wassplit, global_build,apr_tree, &
                          xyzmh_ptmass)
 use dim,       only:maxtypes,mpi,ind_timesteps
 use io,        only:fatal,error
 use mpitree,   only:get_group_cofm,reduce_group
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
 logical,         intent(in)    :: global_build
 logical,         intent(in)    :: apr_tree
 real,            intent(in), optional :: xyzmh_ptmass(:,:)

 integer(kind=8) :: myslot
 real    :: xyzcofm(3)
 real    :: totmass_node
 real    :: xyzcofmg(3)
 real    :: totmassg
 integer :: npnodetot
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

 if (mpi .and. global_build) then
    npnodetot = reduce_group(npnode,'+',level)
 else
    npnodetot = npnode
 endif

 ! following lines to avoid compiler warnings on intent(out) variables
 ir = 0
 il = 0
 nl = 0
 nr = 0

 wassplit    = (npnodetot > minpart)

 xminl = 0.
 xmaxl = 0.
 xminr = 0.
 xmaxr = 0.

 if ((.not. global_build) .and. (npnode  <  1)) return ! node has no particles, just quit

 xyzcofm(:) = 0.

 call compute_nodes_cofm(npnode,nnode,xyzcofm,totmass_node,doparallel)
 ! if this is global node construction, get the cofm and total mass
 ! of all particles in this node (some on other MPI tasks)
 if (mpi .and. global_build) then
    call get_group_cofm(xyzcofm,totmass_node,level,xyzcofmg,totmassg)
    xyzcofm = xyzcofmg
    totmass_node = totmassg
 endif
 ! checks the reduced mass in the case of global maketree
 if (totmass_node<=0. .and. use_apr) call fatal('mtree + apr', &
    'totmass_node==0, something almost certainly wrong with aprmassoftype')
 if (totmass_node<=0.) call fatal('mtree','totmass_node==0',val=totmass_node)

 call set_nodes_properties(npnode,nnode,xyzcofm,totmass_node,mymum,nodeentry,xmini,xmaxi,&
                           level,global_build,doparallel)

 if (use_geosplit) then !--for gravity KDtree, we need the geo centre to split the node
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

    if (maxlevel > maxdepth) call fatal('maketree','maximum tree depth reached !!')
    ! create two children nodes and point to them from current node
    ! always use G&R indexing for global tree
    if (((level < maxlevel_indexed) .or. global_build)) then !.and. (.not. use_geosplit)) then
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
       if ((.not. global_build) .and. ((nl==npnode) .or. (nr==npnode))) then
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
       xminl(1) = treecache(1,inoderange(1,il))
       xminl(2) = treecache(2,inoderange(1,il))
       xminl(3) = treecache(3,inoderange(1,il))
       xmaxl(1) = xminl(1)
       xmaxl(2) = xminl(2)
       xmaxl(3) = xminl(3)
       do ipart=inoderange(1,il)+1,inoderange(2,il)
          xminl(1) = min(xminl(1),treecache(1,ipart))
          xminl(2) = min(xminl(2),treecache(2,ipart))
          xminl(3) = min(xminl(3),treecache(3,ipart))
          xmaxl(1) = max(xmaxl(1),treecache(1,ipart))
          xmaxl(2) = max(xmaxl(2),treecache(2,ipart))
          xmaxl(3) = max(xmaxl(3),treecache(3,ipart))
       enddo

       xminr(1) = treecache(1,inoderange(1,ir))
       xminr(2) = treecache(2,inoderange(1,ir))
       xminr(3) = treecache(3,inoderange(1,ir))
       xmaxr(1) = xminr(1)
       xmaxr(2) = xminr(2)
       xmaxr(3) = xminr(3)
       do ipart=inoderange(1,ir)+1,inoderange(2,ir)
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
       xminl = 0.0
       xmaxl = 0.0
       xminr = 0.0
       xmaxr = 0.0
    endif

    ! Reduce node limits of children across MPI tasks belonging to this group.
    ! The synchronisation needs to happen here, not at the next level, because
    ! the groups will be independent by then.
    if (mpi .and. global_build) then
       xminl(1) = reduce_group(xminl(1),'min',level)
       xminl(2) = reduce_group(xminl(2),'min',level)
       xminl(3) = reduce_group(xminl(3),'min',level)

       xmaxl(1) = reduce_group(xmaxl(1),'max',level)
       xmaxl(2) = reduce_group(xmaxl(2),'max',level)
       xmaxl(3) = reduce_group(xmaxl(3),'max',level)

       xminr(1) = reduce_group(xminr(1),'min',level)
       xminr(2) = reduce_group(xminr(2),'min',level)
       xminr(3) = reduce_group(xminr(3),'min',level)

       xmaxr(1) = reduce_group(xmaxr(1),'max',level)
       xmaxr(2) = reduce_group(xmaxr(2),'max',level)
       xmaxr(3) = reduce_group(xmaxr(3),'max',level)
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

 if (.not.i_lt_pivot) then
    i = i - 1
    dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 endif
 if (j_lt_pivot) then
    j = j + 1
    dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
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

!-----------------------------------------------
!+
!  Routine to update a constructed tree
!  Note: current version ONLY works if
!  tree is built to < maxlevel_indexed
!  That is, it relies on the 2^n style tree
!  indexing to sweep out each level
!+
!-----------------------------------------------
subroutine revtree(node, xyzh, leaf_is_active, ncells)
 use dim,  only:maxp,use_apr,ind_timesteps
 use part, only:maxphase,iphase,igas,massoftype,iamtype,aprmassoftype,&
                apr_level,iactive,treecache,isdead_or_accreted
 use io,   only:fatal
 type(kdnode),    intent(inout) :: node(:) !ncellsmax+1)
 real,            intent(in)    :: xyzh(:,:)
 integer,         intent(inout) :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(in)    :: ncells
 real :: hmax, r2max
 real :: xi, yi, zi, hi
 real :: dx, dy, dz, dr2
#ifdef GRAVITY
 real :: quads(6)
 real :: octs(10)
#endif
 integer :: inode, ipart, ipartidx, i, nptot
 real :: pmassi, totmass
 real :: x0(3)
 real :: xcofm, ycofm, zcofm, fac, dfac
 logical :: nodeisactive

 pmassi = massoftype(igas)

 ! find maximum index in inodeparts that we need to update in treecache
 nptot = 0
 do i=1,int(ncells)
    if (i > 1 .and. node(i)%parent == 0) cycle
    if (inoderange(1,i) > 0 .and. inoderange(2,i) >= inoderange(1,i)) then
       nptot = max(nptot, inoderange(2,i))
    endif
 enddo

 ! update treecache for particles in the tree only
 ! mark dead/accreted particles by setting treecache(4,i) negative
 !$omp parallel default(none) &
 !$omp shared(nptot,inodeparts,xyzh,iphase,apr_level) &
 !$omp shared(massoftype,aprmassoftype,treecache) &
 !$omp shared(maxphase,maxp) &
 !$omp private(i,ipartidx)
 !$omp do schedule(static)
 do i=1,nptot
    if (inodeparts(i) == 0) cycle
    ipartidx = abs(inodeparts(i))
    treecache(1:4,i) = xyzh(1:4,ipartidx)
    ! compute and store mass
    if (maxphase==maxp) then
       if (use_apr) then
          treecache(5,i) = aprmassoftype(iamtype(iphase(ipartidx)),apr_level(ipartidx))
       else
          treecache(5,i) = massoftype(iamtype(iphase(ipartidx)))
       endif
    elseif (use_apr) then
       treecache(5,i) = aprmassoftype(igas,apr_level(ipartidx))
    else
       treecache(5,i) = massoftype(igas)
    endif
 enddo
 !$omp enddo
 !$omp end parallel

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(ncells) &
!$omp shared(node,inoderange,inodeparts,treecache,leaf_is_active) &
!$omp private(hmax,r2max,xi,yi,zi,hi) &
!$omp private(dx,dy,dz,dr2,inode,ipart,x0) &
!$omp private(xcofm,ycofm,zcofm,fac,dfac,nodeisactive) &
#ifdef GRAVITY
!$omp private(quads,octs) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)
!$omp do schedule(guided)
 over_nodes: do inode=1,int(ncells)
    if (inode > 1 .and. node(inode)%parent == 0) cycle
    ! initialize node properties
    node(inode)%xcen(:) = 0.
    node(inode)%size    = 0.
    node(inode)%hmax    = 0.
#ifdef GRAVITY
    node(inode)%mass    = 0.
    node(inode)%quads(:)= 0.
    node(inode)%octs(:) = 0.
#endif
    ! initialize leaf_is_active (will be set for leaf nodes below)
    leaf_is_active(inode) = 0

    ! check if node has particles
    if (inoderange(1,inode) <= 0 .or. inoderange(2,inode) < inoderange(1,inode)) cycle over_nodes

    ! find centre of mass from particle list using same algorithm as maketree
    ! also check for active particles and compute hmax during this loop
    xcofm = 0.
    ycofm = 0.
    zcofm = 0.
    totmass = 0.0
    hmax = 0.
    dfac = 1.
    if (pmassi > 0.) then
       dfac = 1./pmassi
    endif
    nodeisactive = .false.
    do ipart = inoderange(1,inode), inoderange(2,inode)
       if (inodeparts(ipart) == 0) cycle
       xi = treecache(1,ipart)
       yi = treecache(2,ipart)
       zi = treecache(3,ipart)
       hi = treecache(4,ipart)
       ! check condition after loading (dead/accreted particles have hi <= 0)
       if (hi <= 0.) cycle
       hi = abs(hi)
       pmassi = treecache(5,ipart)
       ! check for active particles (for leaf_is_active flag)
       if (ind_timesteps .and. .not. nodeisactive) then
          if (inodeparts(ipart) > 0) nodeisactive = .true.
       endif
       fac = pmassi*dfac
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
       totmass = totmass + pmassi
       hmax = max(hi, hmax)
    enddo
    if (.not. ind_timesteps) nodeisactive = .true.

    if (totmass <= 0.0) cycle over_nodes

    x0(1) = xcofm/(totmass*dfac)
    x0(2) = ycofm/(totmass*dfac)
    x0(3) = zcofm/(totmass*dfac)

    ! update cell size and quads
    r2max = 0.
#ifdef GRAVITY
    quads = 0.
    octs  = 0.
#endif
    do ipart = inoderange(1,inode), inoderange(2,inode)
       ! load all treecache values sequentially (1,2,3,4,5) for cache efficiency
       xi = treecache(1,ipart)
       yi = treecache(2,ipart)
       zi = treecache(3,ipart)
       hi = treecache(4,ipart)
       ! check condition after loading (dead/accreted particles have hi <= 0)
       if (hi <= 0.) cycle
       pmassi = treecache(5,ipart)
       dx = xi - x0(1)
       dy = yi - x0(2)
       dz = zi - x0(3)
       dr2 = dx*dx + dy*dy + dz*dz
       r2max = max(dr2, r2max)
#ifdef GRAVITY
       call add_node_moments(pmassi,dx,dy,dz,quads,octs)
#endif
    enddo

    node(inode)%xcen(1) = x0(1)
    node(inode)%xcen(2) = x0(2)
    node(inode)%xcen(3) = x0(3)
    node(inode)%size = sqrt(r2max) + epsilon(r2max)
    node(inode)%hmax = hmax
#ifdef GRAVITY
    node(inode)%mass = totmass
    node(inode)%quads = quads
    node(inode)%octs  = octs
#endif

    ! set leaf_is_active flag for leaf nodes (matching maketree behavior)
    if (node(inode)%leftchild == 0 .and. node(inode)%rightchild == 0) then
       if (ind_timesteps) then
          if (nodeisactive) then
             leaf_is_active(inode) = 1
          else
             leaf_is_active(inode) = -1
          endif
       else
          leaf_is_active(inode) = 1
       endif
    endif
 enddo over_nodes
!$omp enddo
!$omp end parallel

end subroutine revtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the global level tree
!+
!-------------------------------------------------------------------------------
subroutine maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,&
                          np,cellatid,leaf_is_active,ncells,apr_tree,nptmass,xyzmh_ptmass)
 use io,           only:fatal,warning,id,nprocs,master
 use mpiutils,     only:reduceall_mpi
 use mpibalance,   only:balancedomains
 use mpitree,      only:tree_sync,tree_bcast
 use part,         only:isdead_or_accreted,iactive,ibelong,isink,massoftype,igas,&
                        iamtype,maxphase,maxp,aprmassoftype,apr_level,ihsoft
 use timing,       only:increment_timer,get_timings,itimer_balance
 use dim,          only:ind_timesteps

 type(kdnode),    intent(out)   :: nodeglobal(:)    ! ncellsmax+1
 type(kdnode),    intent(out)   :: node(:)          ! ncellsmax+1
 integer,         intent(out)   :: nodemap(:)       ! ncellsmax+1
 integer,         intent(out)   :: globallevel
 integer,         intent(out)   :: refinelevels
 integer,         intent(inout) :: np
 real,            intent(inout) :: xyzh(:,:)
 integer,         intent(out)   :: cellatid(:)      ! ncellsmax+1
 integer,         intent(out)   :: leaf_is_active(:)  ! ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 logical,         intent(in)    :: apr_tree
 integer,         intent(in),    optional :: nptmass
 real,            intent(inout), optional :: xyzmh_ptmass(:,:)
 real                              :: xmini(3),xmaxi(3)
 real                              :: xminl(3),xmaxl(3)
 real                              :: xminr(3),xmaxr(3)
 integer                           :: minlevel, maxlevel
 integer                           :: idleft, idright
 integer                           :: groupsize,ifirstingroup,groupsplit
 type(kdnode)                      :: mynode(1)
 integer                           :: nl, nr
 integer                           :: il, ir, iself, parent
 integer                           :: level
 integer                           :: nnodestart, nnodeend,locstart,locend
 integer                           :: npcounter
 integer                           :: i, k, offset, roffset, roffset_prev, coffset
 integer                           :: inode
 integer                           :: npnode
 logical                           :: wassplit,sinktree
 real(kind=4)                      :: t1,t2,tcpu1,tcpu2

 sinktree = .false.
 if (present(nptmass).and.present(xyzmh_ptmass)) sinktree=.true.
 parent = 0
 iself = irootnode
 leaf_is_active = 0

 ! root is level 0
 globallevel = int(ceiling(log(real(nprocs)) / log(2.0)))

 minlevel = maxdepth - 1
 maxlevel = 0

 levels: do level = 0, globallevel
    groupsize = 2**(globallevel - level)
    ifirstingroup = (id / groupsize) * groupsize
    if (level == 0) then
       if (sinktree) then
          call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,&
                                   xyzmh_ptmass,nptmass)
       else
          call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh)
       endif
    else
       npcounter = npnode
    endif
    if (sinktree) then
       call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit,.true.,apr_tree,xyzmh_ptmass)
    else
       call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
                        il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                        minlevel, maxlevel, wassplit,.true.,apr_tree)
    endif

    if (.not.wassplit) then
       call fatal('maketreeglobal','insufficient particles for splitting at the global level: '// &
            'use more particles or less MPI threads')
    endif

    ! set which tree child this proc will belong to next
    groupsplit = ifirstingroup + (groupsize / 2)

    ! record parent for next round
    parent = iself

    ! which half of the tree this task is on
    if (id < groupsplit) then
       ! i for the next node we construct
       iself = il
       ! the left and right task IDs
       idleft = id
       idright = id + 2**(globallevel - level - 1)
       xmini = xminl
       xmaxi = xmaxl
    else
       iself = ir
       idleft = id - 2**(globallevel - level - 1)
       idright = id
       xmini = xminr
       xmaxi = xmaxr
    endif
    if (sinktree) then
       if (nptmass>0) then
          ibelong(maxpsph+1:maxpsph+nptmass) = -1
       endif
    endif
    if (npcounter > 0) then
       do i = inoderange(1,il), inoderange(2,il)
          ibelong(abs(inodeparts(i))) = idleft
       enddo
       do i = inoderange(1,ir), inoderange(2,ir)
          ibelong(abs(inodeparts(i))) = idright
       enddo
    endif

    call get_timings(t1,tcpu1)
    ! move particles to where they belong
    call balancedomains(np)
    call get_timings(t2,tcpu2)
    if (sinktree) ibelong(maxpsph+1:maxpsph+nptmass) = int(reduceall_mpi("max", ibelong(maxpsph+1:maxpsph+nptmass)))
    call increment_timer(itimer_balance,t2-t1,tcpu2-tcpu1)
    ! move particles from old array
    ! this is a waste of time, but maintains compatibility
    npnode = 0
    do i=1,np
       npnode = npnode + 1
       !
       ! tag inactive particles with negative index
       ! in the particle list for the node
       !
       if (ind_timesteps) then
          if (iactive(iphase(i))) then
             inodeparts(npnode) = i
          else
             inodeparts(npnode) = -i
          endif
       else
          inodeparts(npnode) = i
       endif
       treecache(1:4,npnode) = xyzh(1:4,i)
       if (maxphase==maxp) then
          if (use_apr) then
             treecache(5,npnode) = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             treecache(5,npnode) = massoftype(iamtype(iphase(i)))
          endif
       elseif (use_apr) then
          treecache(5,npnode) = aprmassoftype(igas,apr_level(i))
       else
          treecache(5,npnode) = massoftype(igas)
       endif
    enddo
    if (sinktree) then
       if (nptmass > 0) then
          do i=1,nptmass
             if (ibelong(maxpsph + i) /= id) cycle
             if (xyzmh_ptmass(4,i) < 0.) cycle ! dead sink particle
             npnode = npnode + 1
             inodeparts(npnode) = maxpsph + i
             treecache(1:3,npnode) = xyzmh_ptmass(1:3,i)
             treecache(4,npnode)   = xyzmh_ptmass(ihsoft,i)
             treecache(5,npnode)   = xyzmh_ptmass(4,i)
          enddo
       endif
    endif

    ! set all particles to belong to this node
    inoderange(1,iself) = 1
    inoderange(2,iself) = npnode

    ! range of newly written tree
    nnodestart = 2**level
    nnodeend = 2**(level + 1) - 1

    ! synchronize tree with other owners if this proc is the first in group
    call tree_sync(mynode,1,nodeglobal(nnodestart:nnodeend),nprocs/groupsize,ifirstingroup,level)

    ! at level 0, tree_sync already 'broadcasts'
    if (level > 0) then
       ! tree broadcast to non-owners
       call tree_bcast(nodeglobal(nnodestart:nnodeend), nnodeend - nnodestart + 1, level)
    endif

 enddo levels

 ! local tree
 if (sinktree) then
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,refinelevels,nptmass,xyzmh_ptmass)
 else
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,refinelevels)
 endif

 ! tree refinement
 refinelevels = int(reduceall_mpi('min',refinelevels),kind=kind(refinelevels))
 roffset_prev = 1

 irefine = 0
 do i = 1,refinelevels
    offset = 2**(globallevel + i)
    roffset = 2**i

    nnodestart = offset
    nnodeend   = 2*nnodestart-1

    if (nnodeend > ncellsmax) call fatal('kdtree', 'global tree refinement has exceeded ncellsmax')

    locstart   = roffset
    locend     = 2*locstart-1

    ! index shift the node to the global level
    do k = roffset,2*roffset-1
       refinementnode(k) = node(k)
       coffset = refinementnode(k)%parent - roffset_prev

       refinementnode(k)%parent = 2**(globallevel + i - 1) + id * roffset_prev + coffset

       if (i /= refinelevels) then
          refinementnode(k)%leftchild  = 2**(globallevel + i + 1) + 2*id*roffset + 2*(k - roffset)
          refinementnode(k)%rightchild = refinementnode(k)%leftchild + 1
       else
          refinementnode(k)%leftchild = 0
          refinementnode(k)%rightchild = 0
       endif
    enddo

    roffset_prev = roffset
    ! sync, replacing level with globallevel, since all procs will get synced
    ! and deeper comms do not exist
    call tree_sync(refinementnode(locstart:locend),roffset, &
                   nodeglobal(nnodestart:nnodeend),nnodestart-nnodeend, &
                   id,globallevel)

    ! get the mapping from the local tree to the global tree, for future hmax updates
    do inode = locstart,locend
       nodemap(inode) = nnodestart + (id * roffset) + (inode - locstart)
    enddo
 enddo
!  The index up to which the local tree is copied to the global tree
 irefine = 2*roffset-1

 ! cellatid is zero by default
 cellatid = 0
 do i = 1,nprocs
    offset = 2**(globallevel+refinelevels)
    roffset = 2**refinelevels
    do k = 1,roffset
       cellatid(offset + (i - 1) * roffset + (k - 1)) = i
    enddo
 enddo

end subroutine maketreeglobal

end module kdtree
