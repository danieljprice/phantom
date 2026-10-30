!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module neighkdtree
!
! This module contains all routines required for
!  tree based neighbour-finding
!
!  THIS VERSION USES A K-D TREE
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters:
!   - tree_accuracy : *tree opening criterion (0.0-1.0)*
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, infile_utils, io,
!   kdtree, kernel, mpitree, mpiutils, part
!
 use dim,          only:ncellsmax,gravity
 use dtypekdtree,  only:kdnode,lenfgrav
 use kdtree,       only:inoderange,inodeparts,irootnode,maxdepth,tree_accuracy,&
                        ih1,im,irho,izetaomega,isoftomega
 implicit none

 ! the global tree, replicated on every MPI task: cellatid is owner+1 for its leaves, 0 otherwise.
 ! Its top levels split the domains of the tasks (the domain of task r is node nprocs+r),
 ! and below them each task's block of nodes copied from its local tree is stored in task order
 integer,      public,  allocatable :: cellatid(:)
 type(kdnode), public,  allocatable :: nodeglobal(:)
 ! index in the global tree of each node of the local tree copied into it, 0 for the others
 integer,               allocatable :: nodemap(:)
 ! the local nodes copied into the global tree, in the order of this task's block
 integer,               allocatable :: mirror_list(:)
 integer                            :: nmirror = 0
 integer,               parameter   :: nmirror_min = 256        ! fewest nodes copied by a task
 integer,               parameter   :: nactive_per_mirror = 4   ! active particles per copied node
 ! the number of nodes in each task's block and the offset of the block
 integer,               allocatable :: nmirror_task(:),mirror_displ(:)
 ! the combined tree: with more than one MPI task, the global tree is appended to the
 ! local tree (global node g at ncells+g), except for the nodes of this task, which are
 ! the local nodes themselves, so that each local leaf's branch runs up to the global root
 integer,      public               :: ncomb = 0        ! last node of the combined tree, 0 if none
 integer,      public               :: irootcomb = 1    ! root of the combined tree
 integer,      public,  allocatable :: remote_leaf_owner(:) ! owner+1 of the leaves of other tasks, else 0
 ! the owner's local index of each node of the global tree, and of each node of other tasks
 ! in the combined tree: negative for a leaf of the owner's local tree. A leaf of another
 ! task with a positive index is truncated: its children were not imported
 integer,               allocatable :: global_local_index(:)
 integer,      public,  allocatable :: remote_local_index(:)
 ! nodes imported from other tasks (the locally essential tree) are appended after ncomb_tree
 integer                            :: ncomb_tree = 0
 integer                            :: ncells_appended = 0 ! where the global tree was appended
 integer,               allocatable :: let_mark(:)       ! local nodes opened for the current receiver
 integer,               allocatable :: global_owner(:)   ! task owning each global node, -1 above the domains
 ! for each global leaf of another task, the largest size of a node that its walks can
 ! pair with a node of this task, by the number of steps the walk has taken below the leaf
 ! (the last entry is for all deeper steps). Only the branches above active leaves count,
 ! since only active leaves walk. Stored in a pool: profile_len(g) values from
 ! profile_start(g); 0 if the leaf has no active particles, -1 if not received
 integer,               parameter   :: nprofile = 24
 real,                  allocatable :: profile_pool(:)
 integer,               allocatable :: profile_start(:),profile_len(:)
 logical,               allocatable :: global_active(:)  ! a global node with active leaves below
 integer,               allocatable :: branch_mark(:)     ! local nodes above active leaves
 integer                            :: let_stamp = 0
 integer(kind=8)                    :: nlet_sent = 0, nlet_recv = 0, nlet_p2p = 0
 ! reals per imported node: xcen, size, hmax, mass, parent, side, local index (+ quads, octs)
 integer,               parameter   :: nlet_record = 9 + merge(16,0,gravity)
 type(kdnode), public,  allocatable :: node(:)
 integer,      public,  allocatable :: leaf_is_active(:) ! : 0 internal node or empty cell, : 1 active cell, :- inactive cell
 integer,      public,  allocatable :: active_leaves(:)  ! the cells with leaf_is_active > 0, in order (set by build_tree)
 integer,      public               :: nactive_leaves = 0
 ! number of active leaves handed to a thread at a time in the density and force loops:
 ! cells nearby in the list share neighbours, so small chunks keep them on one thread,
 ! but with self-gravity the cost per cell varies widely, so the cells are handed out one at a time
 integer,      parameter, public    :: ncellchunk = merge(1,16,gravity)
 integer,      public , allocatable :: listneigh(:)
 integer,      public , allocatable :: listneigh_global(:)
!-- dual tree cache arrays
! cachestate is read without an atomic in the tree walks, which read it for every node visited:
! some compilers (e.g. aocc) put every OpenMP atomic behind one global lock, which serialises
! the walk. This is safe because the state of a node only increases during a walk (0 -> 1 -> 2 -> 3),
! it is written with an atomic and a flush after the data it refers to, and readers flush before
! using cached data: a reader seeing an old value just computes the interaction instead of reusing it.
! It is volatile so that every read is a fresh load from memory
 integer, volatile,     allocatable :: cachestate(:)
 integer,               allocatable :: neighnodecount_branch(:)
 integer,               allocatable :: neighnode_branch(:,:)
 integer,               allocatable :: neighnodecache(:)
 integer,               allocatable :: neighnodecache_start(:)
 integer,               allocatable :: neighnodecache_count(:)
 integer                            :: itail_neigh = 0 ! next free slot in neighnodecache
 real,                  allocatable :: fnodecache(:,:)
 real,                  allocatable :: fnode_branch(:,:)
!$omp threadprivate(fnode_branch,neighnode_branch,neighnodecount_branch)
!$omp threadprivate(listneigh)
 ! leaves of other tasks met in leaf-leaf pairs by the last dual tree walk of this thread
 integer, public, allocatable :: remote_leaves(:)
 integer, public              :: nremote_leaves = 0
!$omp threadprivate(remote_leaves,nremote_leaves)
 integer(kind=8), public            :: ncells
 real,            public            :: dxcell
 real,            public            :: dcellx = 0.,dcelly = 0.,dcellz = 0.
 logical,         public            :: use_dualtree  = .true.
 logical,         public            :: use_dualcache = gravity ! only debug flag / should be on with gravity
 integer,         parameter         :: maxnodecache_local = 512
 integer,         parameter         :: maxneigh_per_node  = 16
 integer,         parameter         :: maxstacksize = 2048
 ! revise rather than rebuild the tree if fewer than this fraction of particles are active,
 ! until the particles have had about one step each since the tree was built
 real,            parameter         :: frac_active_revtree = 0.1
 integer(kind=8)                    :: nactive_since_build = 0
 integer                            :: globallevel

 public :: allocate_neigh, deallocate_neigh
 public :: build_tree, get_neighbour_list, write_options_tree, read_options_tree
 public :: get_distance_from_centre_of_mass, getneigh_pos
 public :: set_hmaxcell,get_hmaxcell
 public :: get_cell_location
 public :: sync_hmax_mpi
 public :: let_exchange,print_let_stats,get_leaf_neighbours
 public :: get_particle_targets
 public :: get_node_node_interaction,expand_fgrav_in_taylor_series

 private

contains

!-----------------------------------------------------------------------
!+
!  allocate memory for the neighbour list
!+
!-----------------------------------------------------------------------
subroutine allocate_neigh
 use allocutils, only:allocate_array
 use kdtree,     only:allocate_kdtree,maxdepth
 use dim,        only:maxp
 use io,         only:nprocs
 integer(kind=8) :: nglobal

 integer :: nnodes,ncache

 ! the top levels, which split the domains of the tasks, then each task's block of nodes
 nglobal = 2*nprocs + int(nprocs,kind=8)*mirror_budget()
 ! room for the combined tree: the global tree appended after the local tree, then as
 ! many nodes again imported from other tasks. The walk caches only the nodes on the
 ! branches of local leaves: the local nodes and the top levels of the global tree
 nnodes = ncellsmax + 1
 ncache = ncellsmax + 1
 if (nprocs > 1) then
    nnodes = nnodes + 2*int(nglobal)
    ncache = ncache + 2*nprocs
 endif
 call allocate_array('cellatid',       cellatid,       nglobal           )
 call allocate_array('leaf_is_active', leaf_is_active, nnodes            )
 call allocate_array('active_leaves',  active_leaves,  ncellsmax+1       )
 call allocate_array('nodeglobal',     nodeglobal,     nglobal           )
 call allocate_array('node',           node,           nnodes            )
 call allocate_array('nodemap',        nodemap,        ncellsmax+1       )
 call allocate_array('remote_leaf_owner',remote_leaf_owner,nnodes        )
 call allocate_array('remote_local_index',remote_local_index,nnodes      )
 call allocate_array('global_local_index',global_local_index,nglobal     )
 call allocate_array('let_mark',       let_mark,       ncellsmax+1       )
 call allocate_array('global_owner',   global_owner,   nglobal           )
 if (nprocs > 1) then
    call allocate_array('profile_start',profile_start,int(nglobal))
    call allocate_array('profile_len',  profile_len,  int(nglobal))
    allocate(global_active(nglobal))
    call allocate_array('branch_mark',  branch_mark,  ncellsmax+1)
    branch_mark = 0
 endif
 nodemap = 0
 nmirror = 0
 remote_leaf_owner = 0
 remote_local_index = 0
 let_mark = 0
 call allocate_kdtree()
 call allocate_array('listneigh_global',listneigh_global,maxp)
 call allocate_array('cachestate', cachestate, ncache)

 if (use_dualcache) then
    call allocate_array('fnodecache', fnodecache, lenfgrav, ncache)
    call allocate_array('neighnodecache',neighnodecache,ncellsmax*maxneigh_per_node)
    call allocate_array('neighnodecache_start',neighnodecache_start,ncache)
    call allocate_array('neighnodecache_count',neighnodecache_count,ncache)
    itail_neigh = 0
 endif

!$omp parallel
 call allocate_array('listneigh',listneigh,maxp)
 call allocate_array('fnode_branch', fnode_branch, lenfgrav, maxdepth+1)
 if (use_dualcache) then
    call allocate_array('neighnodecount_branch',neighnodecount_branch,maxdepth+1)
    call allocate_array('neighnode_branch',neighnode_branch,maxnodecache_local,maxdepth+1)
 endif
!$omp end parallel

end subroutine allocate_neigh

!-----------------------------------------------------------------------
!+
!  deallocate memory for the neighbour list
!+
!-----------------------------------------------------------------------
subroutine deallocate_neigh
 use kdtree,   only:deallocate_kdtree

 if (allocated(cellatid)) deallocate(cellatid)
 if (allocated(leaf_is_active)) deallocate(leaf_is_active)
 if (allocated(active_leaves)) deallocate(active_leaves)
 if (allocated(nodeglobal)) deallocate(nodeglobal)
 if (allocated(node)) deallocate(node)
 if (allocated(nodemap)) deallocate(nodemap)
 if (allocated(remote_leaf_owner)) deallocate(remote_leaf_owner)
 if (allocated(remote_local_index)) deallocate(remote_local_index)
 if (allocated(global_local_index)) deallocate(global_local_index)
 if (allocated(let_mark)) deallocate(let_mark)
 if (allocated(global_owner)) deallocate(global_owner)
 if (allocated(profile_pool)) deallocate(profile_pool)
 if (allocated(profile_start)) deallocate(profile_start,profile_len,global_active)
 if (allocated(branch_mark)) deallocate(branch_mark)
 if (allocated(mirror_list)) deallocate(mirror_list)
 if (allocated(nmirror_task)) deallocate(nmirror_task,mirror_displ)
 if (allocated(cachestate)) deallocate(cachestate)
 if (allocated(fnodecache)) deallocate(fnodecache)
 if (allocated(neighnodecache)) deallocate(neighnodecache)
 if (allocated(neighnodecache_start)) deallocate(neighnodecache_start)
 if (allocated(neighnodecache_count)) deallocate(neighnodecache_count)
!$omp parallel
 if (allocated(neighnode_branch)) deallocate(neighnode_branch)
 if (allocated(neighnodecount_branch)) deallocate(neighnodecount_branch)
 if (allocated(fnode_branch)) deallocate(fnode_branch)
 if (allocated(listneigh)) deallocate(listneigh)
!$omp end parallel
 if (allocated(listneigh_global)) deallocate(listneigh_global)
 call deallocate_kdtree()

end subroutine deallocate_neigh

!-----------------------------------------------------------------------
!+
!  get the hmax value of a cell
!+
!-----------------------------------------------------------------------
subroutine get_hmaxcell(inode,hmaxcell)
 integer, intent(in)  :: inode
 real,    intent(out) :: hmaxcell

 hmaxcell = node(inode)%hmax

end subroutine get_hmaxcell

!-----------------------------------------------------------------------
!+
!  set the hmax value of a cell and propagate the value up the tree
!+
!-----------------------------------------------------------------------
subroutine set_hmaxcell(inode,hmaxcell)
 integer, intent(in) :: inode
 real,    intent(in) :: hmaxcell
 integer :: n
 real    :: hmaxn

 n = inode
 ! never lower it: with individual timesteps the caller's hmaxcell covers only the
 ! cell's active particles (dens.F90 compute_hmax), and the inactive ones in the same
 ! leaf keep their own h, which the symmetric neighbour search in the force pass needs
 ! the leaf's hmax to cover.  The tree is rebuilt every step, so the bound never has
 ! to shrink.
 node(n)%hmax = max(node(n)%hmax, hmaxcell)

 ! walk tree up, stopping at the first ancestor whose hmax already covers hmaxcell:
 ! a node's hmax is never below its children's, so neither is any of its ancestors'.
 ! Other threads update the same ancestors, so both the test and the update are atomic.
 ! Stop at the local root: the global nodes above it are updated by sync_hmax_mpi,
 ! so that they stay the same on every MPI task
 do while (n /= irootnode .and. node(n)%parent /= 0)
    n = node(n)%parent
!$omp atomic read
    hmaxn = node(n)%hmax
    if (hmaxn >= hmaxcell) exit
!$omp atomic
    node(n)%hmax = max(node(n)%hmax, hmaxcell)
 enddo

end subroutine set_hmaxcell

!-----------------------------------------------------------------------
!+
!  Public helper to test compute_M2L (for test and better inlining)
!+
!-----------------------------------------------------------------------
pure subroutine get_node_node_interaction(dx,dy,dz,dr1,q0,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr1,q0,quads(6)
 real, intent(inout) :: fnode(lenfgrav)

 call compute_M2L(dx,dy,dz,dr1,q0,quads,fnode)

end subroutine get_node_node_interaction

!-----------------------------------------------------------------------
!+
!  get the distance from the centre of mass of a cell
!+
!-----------------------------------------------------------------------
subroutine get_distance_from_centre_of_mass(inode,xi,yi,zi,dx,dy,dz)
 integer, intent(in)  :: inode
 real,    intent(in)  :: xi,yi,zi
 real,    intent(out) :: dx,dy,dz
 real :: xdum,ydum,zdum

 call get_sep(node(inode)%xcen,(/xi,yi,zi/),dx,dy,dz,xdum,ydum,zdum)

end subroutine get_distance_from_centre_of_mass

!-----------------------------------------------------------------------
!+
!  build the tree
!+
!-----------------------------------------------------------------------
subroutine build_tree(npart,nactive,xyzh,vxyzu,for_apr)
 use io,           only:nprocs
 use kdtree,       only:maketree,maketreeglobal,revtree
 use dim,          only:mpi,use_sinktree
 use part,         only:nptmass,xyzmh_ptmass,maxp
 use allocutils,   only:allocate_array
 integer, intent(inout) :: npart
 integer, intent(in)    :: nactive
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: vxyzu(:,:)
 logical, intent(in), optional :: for_apr
 logical :: apr_tree,revised

 apr_tree = .false.
 if (present(for_apr)) apr_tree = for_apr

 !
 ! the listneigh array is threadprivate, but if the thread numbers or ids are changed
 ! then the memory might be lost. So the following lines are a failsafe
 ! to ensure that the listneigh array is always allocated for each thread
 !
 !$omp parallel default(none) shared(maxp)
 if (.not. allocated(listneigh)) call allocate_array('listneigh',listneigh,maxp)
 !$omp end parallel

 if (mpi .and. nprocs > 1) then
    if (use_sinktree) then
       call maketreeglobal(nodeglobal,node,globallevel,xyzh,npart,leaf_is_active,ncells,&
                           apr_tree,nptmass,xyzmh_ptmass)
    else
       call maketreeglobal(nodeglobal,node,globallevel,xyzh,npart,leaf_is_active,ncells,&
                           apr_tree)
    endif
    call graft_local_trees(nactive)
    call combine_trees()
 else
    if (use_sinktree) then
       call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree,nptmass=nptmass,xyzmh_ptmass=xyzmh_ptmass)
    else
       ! with few active particles, revise the last tree for the new positions
       ! rather than rebuilding it: the nodes are looser, but only the active
       ! particles have to find their neighbours in it
       revised = .false.
       nactive_since_build = nactive_since_build + nactive
       if (.not.apr_tree .and. nactive < frac_active_revtree*npart .and. nactive_since_build < npart) then
          call revtree(node,xyzh,npart,leaf_is_active,ncells,revised)
       endif
       if (.not.revised) then
          call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree)
          nactive_since_build = 0
       endif
    endif
 endif
 call list_active_leaves()

 ! the dual tree walk cache refers to the old nodes: reset it for the new tree
 if (use_dualcache) then
    cachestate(1:min(max(int(ncells),ncomb),size(cachestate))) = 0 ! 0=untouched, 1=claimed, 2=fnode cached, 3=fully cached
    itail_neigh = 0
 endif

end subroutine build_tree

!-----------------------------------------------------------------------
!+
!  list the cells with leaf_is_active > 0, so that loops over cells
!  (density, force) need not visit every node of the tree: with
!  individual timesteps only a few leaves may be active
!
!  In parallel, in cell order: the active leaves are counted per chunk of
!  cells, and a running total gives each chunk its place in the list.
!+
!-----------------------------------------------------------------------
subroutine list_active_leaves()
!$ use omp_lib, only:omp_get_max_threads
 integer, allocatable :: nlist(:)
 integer :: icell,ic,nchunk,n,k

 n = int(ncells)
 nchunk = 1
!$ nchunk = omp_get_max_threads()
 allocate(nlist(0:nchunk))
 nlist = 0

 !$omp parallel default(none) &
 !$omp shared(n,leaf_is_active,nchunk,nlist,active_leaves,nactive_leaves) &
 !$omp private(icell,ic,k)
 !$omp do schedule(static)
 do ic=1,nchunk
    k = 0
    do icell=int((int(ic-1,8)*n)/nchunk)+1,int((int(ic,8)*n)/nchunk)
       if (leaf_is_active(icell) > 0) k = k + 1
    enddo
    nlist(ic) = k
 enddo
 !$omp enddo
 !$omp single
 do ic=1,nchunk
    nlist(ic) = nlist(ic) + nlist(ic-1)
 enddo
 nactive_leaves = nlist(nchunk)
 !$omp end single
 !$omp do schedule(static)
 do ic=1,nchunk
    k = nlist(ic-1)
    do icell=int((int(ic-1,8)*n)/nchunk)+1,int((int(ic,8)*n)/nchunk)
       if (leaf_is_active(icell) > 0) then
          k = k + 1
          active_leaves(k) = icell
       endif
    enddo
 enddo
 !$omp enddo
 !$omp end parallel
 deallocate(nlist)

end subroutine list_active_leaves

!-----------------------------------------------------------------------
!+
! Using the k-d tree, compiles the neighbour list for the
! current cell (this list is common to all particles in the cell)
!
! the list is returned in 'listneigh' (length nneigh)
!+
!-----------------------------------------------------------------------
subroutine get_neighbour_list(inode,mylistneigh,nneigh,xyzh,xyzcache,ixyzcachesize, &
                              getj,f,remote_export,cell_xpos,cell_xsizei,cell_rcuti)
 use io,       only:nprocs,warning,fatal
 use dim,      only:mpi
 use kernel,   only:radkern
 use part,     only:gravity,periodic
 use boundary, only:dxbound,dybound,dzbound
 integer, intent(in)  :: inode,ixyzcachesize
 integer, intent(out) :: mylistneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: xyzcache(:,:)
 logical, intent(in),  optional :: getj
 real,    intent(out), optional :: f(lenfgrav)
 logical, intent(out), optional :: remote_export(:)
 real,    intent(in),  optional :: cell_xpos(3),cell_xsizei,cell_rcuti
 real :: xpos(3)
 real :: fgrav(lenfgrav),fgrav_global(lenfgrav)
 real :: xsizei,rcuti
 integer :: nrem,ntr,k
 logical :: get_j,global_search,get_f
!
!--retrieve geometric centre of the node and the search radius (e.g. 2*hmax)
!
 if (present(cell_xpos)) then
    xpos = cell_xpos
    xsizei = cell_xsizei
    rcuti = cell_rcuti
 else
    call get_cell_location(inode,xpos,xsizei,rcuti)
 endif

 ! the global tree only exists with more than one MPI task
 global_search = .false.
 if (present(remote_export)) then
    if (nprocs > 1) global_search = .true.
    remote_export = .false.
 endif

 if (periodic) then
    if (rcuti > 0.5*min(dxbound,dybound,dzbound)) then
       call warning('get_neighbour_list', '2h > 0.5*L in periodic neighb. '//&
                'search: USE HIGHER RES, BIGGER BOX or LOWER MINPART IN TREE')
    endif
 endif
 !
 !--perform top-down tree walk to find all particles within radkern*h
 !  and force due to node-node interactions
 !
 get_j = .false.
 if (present(getj)) get_j = getj

 get_f = (gravity .and. present(f))

 if (mpi .and. global_search .and. get_f .and. use_dualtree) then
    ! with more than one MPI task: dual tree walk of the combined tree, with the nodes
    ! imported from other tasks. Leaves of other tasks paired with this cell are listed
    ! in remote_leaves, and the cell is sent to their owners
    call getneigh_dual(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                       leaf_is_active,get_j,get_f,f,inode,nrem,ntr)
    if (ntr > 0) call fatal('get_neighbour_list','dual tree walk met a node of another task whose children were not imported')
    do k = 1,nremote_leaves
       remote_export(remote_leaf_owner(remote_leaves(k))) = .true.
    enddo
    !$omp atomic
    nlet_p2p = nlet_p2p + nrem
    return
 endif

 if (mpi .and. global_search) then ! no sym fmm for now...
    ! Find MPI tasks that have neighbours of this cell, output to remote_export
    call getneigh(nodeglobal,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                  cellatid,get_j,get_f,fgrav_global,remote_export)
 elseif (get_f) then
    ! Set fgrav to zero, which matters if gravity is enabled but global search is not
    fgrav_global = 0.0
 endif

 ! Find neighbours of this cell on this node
 ! the dual tree walk is not yet used with more than one MPI task
 if (get_f .and. .not.(mpi .and. nprocs > 1) .and. use_dualtree) then
    call getneigh_dual(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                          leaf_is_active,get_j,get_f,fgrav,inode)
 else
    call getneigh(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                     leaf_is_active,get_j,get_f,fgrav)
 endif

 if (get_f) f = fgrav + fgrav_global

end subroutine get_neighbour_list

!-----------------------------------------------------------------------
!+
!  get neighbours around an arbitrary position in space
!+
!-----------------------------------------------------------------------
subroutine getneigh_pos(xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,get_j)
 integer, intent(in)  :: ixyzcachesize
 real,    intent(in)  :: xpos(3)
 real,    intent(in)  :: xsizei,rcuti
 integer, intent(out) :: mylistneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(out) :: xyzcache(:,:)
 integer, intent(in)  :: leaf_is_active(:) !ncellsmax+1)
 logical, intent(in), optional :: get_j
 logical :: getj

 getj = .false.
 if (present(get_j)) getj=get_j
 call getneigh(node,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize, &
               leaf_is_active,getj,.false.)

end subroutine getneigh_pos

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_tree(iunit)
 use kdtree,       only:tree_accuracy
 use infile_utils, only:write_inopt
 use part,         only:gravity
 integer, intent(in) :: iunit

 if (gravity) call write_inopt(tree_accuracy,'tree_accuracy','tree opening criterion (0.0-1.0)',iunit)

end subroutine write_options_tree

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_tree(db,nerr)
 use part,         only:gravity
 use kdtree,       only:tree_accuracy
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 if (gravity) call read_inopt(tree_accuracy,'tree_accuracy',db,errcount=nerr,min=0.,max=1.)

end subroutine read_options_tree

!-----------------------------------------------------------------------
!+
!  find the position and size of a tree node
!+
!-----------------------------------------------------------------------
subroutine get_cell_location(inode,xpos,xsizei,rcuti)
 use kernel, only:radkern
 integer, intent(in)  :: inode
 real,    intent(out) :: xpos(3)
 real,    intent(out) :: xsizei
 real,    intent(out) :: rcuti

 xpos    = node(inode)%xcen(1:3)
 xsizei  = node(inode)%size
 rcuti   = radkern*node(inode)%hmax

end subroutine get_cell_location

!-----------------------------------------------------------------------
!+
!  the MPI tasks each particle of an exported cell must be sent to, a
!  subset of the tasks found for the whole cell: those owning a leaf of
!  the global tree within the kernel radius of the particle (rcut) or,
!  with getj, within the kernel radius of the leaf's hmax. Tasks that no
!  particle reaches are removed from remote_export
!+
!-----------------------------------------------------------------------
subroutine get_particle_targets(np,xpos,rcut,getj,remote_export,ptargets)
 integer, intent(in)    :: np
 real,    intent(in)    :: xpos(:,:),rcut(:)
 logical, intent(in)    :: getj
 logical, intent(inout) :: remote_export(:)
 logical, intent(out)   :: ptargets(:,:)
 logical :: pexport(size(remote_export))
 integer :: i,j,listdummy(1),ndummy
 real    :: cachedummy(1,1)

 ptargets = .false.
 do i=1,np
    call getneigh(nodeglobal,xpos(:,i),0.,rcut(i),listdummy,ndummy,cachedummy,0,cellatid,&
                  getj,.false.,remote_export=pexport)
    ptargets(i,:) = pexport .and. remote_export
 enddo
 do j=1,size(remote_export)
    remote_export(j) = any(ptargets(1:np,j))
 enddo

end subroutine get_particle_targets

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
 integer :: ntot,ifirst,itask,ibase,p,g,idomain

 ! forget the nodes copied from the previous tree
 nodemap(1) = 0
 if (nmirror > 0) nodemap(mirror_list(1:nmirror)) = 0

 call select_mirrored_nodes(nactive)
 allocate(block(nmirror))
 call pack_mirrored_nodes(block)

 if (.not.allocated(nmirror_task)) allocate(nmirror_task(nprocs),mirror_displ(nprocs))
 call allgather_mpi(nmirror,nmirror_task)
 mirror_displ = get_displs_mpi(nmirror_task)
 ntot = sum(nmirror_task)
 ifirst = 2*nprocs     ! first node after the top levels
 if (ifirst + ntot - 1 > size(nodeglobal)) call fatal('graft_local_trees','global tree exceeds its node budget')
 call allgatherv_nodes(block,nodeglobal(ifirst:),nmirror_task,mirror_displ)

 ! the local index of each copied node on its owner, negative for a leaf of its local tree
 allocate(ilocal(nmirror))
 do p = 1,nmirror
    ilocal(p) = local_index_and_leaf(mirror_list(p))
 enddo
 call allgatherv_mpi(ilocal,global_local_index(ifirst:),nmirror_task,mirror_displ)

 ! the blocks are linked within themselves: shift them to their place in the global tree
 cellatid = 0
 global_owner(1:nprocs-1) = -1
 do itask = 1,nprocs
    ibase = ifirst - 1 + mirror_displ(itask)
    idomain = nprocs + itask - 1
    do p = 1,nmirror_task(itask)
       g = ibase + p
       if (nodeglobal(g)%leftchild > 0)  nodeglobal(g)%leftchild  = nodeglobal(g)%leftchild  + ibase
       if (nodeglobal(g)%rightchild > 0) nodeglobal(g)%rightchild = nodeglobal(g)%rightchild + ibase
       if (nodeglobal(g)%parent > 0) then
          nodeglobal(g)%parent = nodeglobal(g)%parent + ibase
       else
          nodeglobal(g)%parent = idomain
       endif
       if (nodeglobal(g)%leftchild == 0 .and. nodeglobal(g)%rightchild == 0) cellatid(g) = itask
    enddo
    global_local_index(idomain) = merge(-irootnode,irootnode,nmirror_task(itask) == 0)
    global_owner(idomain) = itask - 1
    global_owner(ibase+1:ibase+nmirror_task(itask)) = itask - 1
    ! the children of the domain node come first in the block
    nodeglobal(idomain)%leftchild = 0
    nodeglobal(idomain)%rightchild = 0
    if (nmirror_task(itask) == 0) then
       cellatid(idomain) = itask
    else
       nodeglobal(idomain)%leftchild = ibase + 1
       if (nmirror_task(itask) > 1) then
          if (nodeglobal(ibase+2)%parent == idomain) nodeglobal(idomain)%rightchild = ibase + 2
       endif
    endif
 enddo

 ! where the copied local nodes are in the global tree; the local root is the domain node
 ibase = ifirst - 1 + mirror_displ(id+1)
 nodemap(irootnode) = nprocs + id
 do p = 1,nmirror
    nodemap(mirror_list(p)) = ibase + p
 enddo

end subroutine graft_local_trees

!-----------------------------------------------------------------------
!+
!  build the combined tree: append the global tree after the local tree,
!  except for this task's domain node and block, which are the local root
!  and the local nodes they were copied from. The local root takes the
!  data of the domain node, which the other tasks use, so that every task
!  sees exactly the same node there
!+
!-----------------------------------------------------------------------
subroutine combine_trees()
 use io, only:nprocs,id,fatal
 type(kdnode) :: rootlinks
 integer :: g,n,ng

 ng = 2*nprocs - 1 + sum(nmirror_task)
 if (int(ncells) + ng > size(node)) call fatal('combine_trees','no room to append the global tree')
 ! forget the nodes of other tasks appended to the previous tree, which may
 ! have started at a different place if the local tree has changed size
 if (ncomb > ncells_appended) remote_leaf_owner(ncells_appended+1:ncomb) = 0
 ncells_appended = int(ncells)

 rootlinks = node(irootnode)
 node(irootnode) = nodeglobal(nprocs + id)
 node(irootnode)%leftchild  = rootlinks%leftchild
 node(irootnode)%rightchild = rootlinks%rightchild
 node(irootnode)%level      = rootlinks%level
 node(irootnode)%parent     = combined_index(nodeglobal(nprocs + id)%parent)

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
 irootcomb = combined_index(1)
 ncomb = int(ncells) + ng
 ncomb_tree = ncomb

end subroutine combine_trees

!-----------------------------------------------------------------------
!+
!  whether a node of the global tree is this task's domain node or one
!  of the local nodes it copied into the global tree
!+
!-----------------------------------------------------------------------
logical function is_this_task(g)
 use io, only:nprocs,id
 integer, intent(in) :: g
 integer :: ibase

 ibase = 2*nprocs - 1 + mirror_displ(id+1)
 is_this_task = (g == nprocs + id) .or. (g > ibase .and. g <= ibase + nmirror)

end function is_this_task

!-----------------------------------------------------------------------
!+
!  index in the combined tree of a node of the global tree (0 for none)
!+
!-----------------------------------------------------------------------
integer function combined_index(g)
 use io, only:nprocs,id
 integer, intent(in) :: g
 integer :: ibase

 ibase = 2*nprocs - 1 + mirror_displ(id+1)
 if (g <= 0) then
    combined_index = 0
 elseif (g == nprocs + id) then
    combined_index = irootnode
 elseif (g > ibase .and. g <= ibase + nmirror) then
    combined_index = mirror_list(g - ibase)
 else
    combined_index = int(ncells) + g
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
 real,    allocatable :: sendbuf(:),recvbuf(:)
 integer, allocatable :: list(:),iparent(:),iside(:)
 integer :: nsend(nprocs),nrecv(nprocs),irecv,nlist,nbuf,k,iproc

 ! forget the nodes imported by a previous call
 call combine_trees()
 call exchange_size_profiles()

 allocate(list(ncomb),iparent(ncomb),iside(ncomb))
 allocate(sendbuf(nlet_record*1024))
 nbuf = 0
 nsend = 0
 do iproc = 0,nprocs-1
    if (iproc == id) cycle
    call select_let_nodes(iproc,list,iparent,iside,nlist)
    if (nbuf + nlist*nlet_record > size(sendbuf)) call grow_real(sendbuf,nbuf + nlist*nlet_record)
    do k = 1,nlist
       call put_let_record(sendbuf(nbuf+1:nbuf+nlet_record),list(k),iparent(k),iside(k))
       nbuf = nbuf + nlet_record
    enddo
    nsend(iproc+1) = nlist*nlet_record
    nlet_sent = nlet_sent + nlist
 enddo

 call alltoall_mpi(nsend,nrecv)
 allocate(recvbuf(max(sum(nrecv),1)))
 call alltoallv_mpi(sendbuf,nsend,recvbuf,nrecv)

 irecv = 0
 do iproc = 0,nprocs-1
    if (nrecv(iproc+1) > 0) call attach_let_nodes(iproc,recvbuf(irecv+1:irecv+nrecv(iproc+1)))
    irecv = irecv + nrecv(iproc+1)
 enddo
 nlet_recv = nlet_recv + sum(nrecv)/nlet_record

 ! the walk caches refer to the nodes and smoothing lengths before the import
 if (use_dualcache) then
    cachestate(1:min(ncomb,size(cachestate))) = 0
    itail_neigh = 0
 endif

end subroutine let_exchange

!-----------------------------------------------------------------------
!+
!  the local nodes to send to task iproc: the children of every node
!  below this task's global leaves that a walk on iproc could open, in
!  pre-order, with the parent of each (negative: a global node, by its
!  index in the global tree; positive: a position in the list) and
!  which child it is (1 left, 2 right)
!+
!-----------------------------------------------------------------------
subroutine select_let_nodes(iproc,list,iparent,iside,nlist)
 use io, only:id
 integer, intent(in)  :: iproc
 integer, intent(out) :: list(:),iparent(:),iside(:),nlist
 integer, allocatable :: istack(:,:)
 integer :: nstack,iglob,ileaf,n,ichild(2),k,ipos

 let_stamp = let_stamp + 1
 !
 ! mark the nodes that a walk on iproc could open, following the walk itself:
 ! pairs (node on iproc's side, node on this task's side) of the combined tree
 ! descend together from the root, as in the unsplit dual tree traversal. Once
 ! iproc's side reaches one of its global leaves, whose nodes below are unknown
 ! here, a region bound stands in for them and only this task's side descends
 !
 allocate(istack(3,1024))
 nstack = 0
 call push_pair(istack,nstack,irootcomb,irootcomb,0)
 do while (nstack > 0)
    call let_pair_step(istack,nstack,iproc)
 enddo
 !
 ! list the children of the marked nodes in pre-order: a node marked
 ! below a global leaf is itself in the list, so its parent comes first
 !
 nlist = 0
 nstack = 0
 do iglob = 1,ncomb_tree - int(ncells)
    if (.not.is_this_task(iglob) .or. cellatid(iglob) /= id + 1) cycle
    ileaf = combined_index(iglob)
    if (let_mark(ileaf) == let_stamp) call push_pair(istack,nstack,ileaf,-iglob,0)
 enddo
 do while (nstack > 0)
    n    = istack(1,nstack)
    ipos = istack(2,nstack)
    nstack = nstack - 1
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
!  one step of the descent in select_let_nodes: test the pair on top of
!  the stack, a node on iproc's side against a node on this task's side
!  (both in the combined tree), and push the pairs below it
!+
!-----------------------------------------------------------------------
subroutine let_pair_step(istack,nstack,iproc)
 use io, only:id
 integer, allocatable, intent(inout) :: istack(:,:)
 integer,              intent(inout) :: nstack
 integer,              intent(in)    :: iproc
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
    if (.not.may_open .or. all(kb <= 0)) return
    call mark_opened(ib)
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
       call mark_opened(ib)
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
!  mark a local node as opened for the current receiver, if its children
!  are not already in the global tree
!+
!-----------------------------------------------------------------------
subroutine mark_opened(n)
 use io, only:id
 integer, intent(in) :: n

 if (n > ncells) return
 if (nodemap(n) /= 0) then
    if (cellatid(nodemap(n)) /= id + 1) return  ! an internal node of the global tree
 endif
 let_mark(n) = let_stamp

end subroutine mark_opened

!-----------------------------------------------------------------------
!+
!  whether a node of the combined tree belongs to the given task, or is
!  a node of the global tree above its domain
!+
!-----------------------------------------------------------------------
logical function on_side(n,itask)
 use io, only:id,nprocs
 integer, intent(in) :: n,itask
 integer :: g,k

 if (n <= ncells) then
    on_side = (itask == id)
    return
 endif
 g = n - int(ncells)
 if (global_owner(g) >= 0) then
    on_side = (global_owner(g) == itask)
 else
    ! a top node (heap indexed): an ancestor of the domain node of itask?
    k = nprocs + itask
    do while (k > g)
       k = k/2
    enddo
    on_side = (k == g)
 endif

end function on_side

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
 well_separated = (tree_accuracy**2*r2 > (node1%size + node2%size)**2) .and. &
                  (r2 > (node1%size + node2%size + rcut)**2)

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
 profile_may_open = (tree_accuracy*d < sizes) .or. (d < sizes + rcut)

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
 real,    allocatable :: sendbuf(:),recvbuf(:)
 integer, allocatable :: leaves(:)
 real    :: prof(2,nprofile)
 integer :: nsend(nprocs),nrecv(nprocs),iproc,k,nleaves,nlen,nbuf,irecv,g,n,ng

 call mark_active_branches()
 ! this task's global leaves: its domain node, if a leaf, or the leaves of its block
 allocate(leaves(nmirror+1))
 nleaves = 0
 if (nmirror == 0) then
    nleaves = 1
    leaves(1) = irootnode
 endif
 do k = 1,nmirror
    if (cellatid(nodemap(mirror_list(k))) == id + 1) then
       nleaves = nleaves + 1
       leaves(nleaves) = mirror_list(k)
    endif
 enddo

 ! records of (global index, length, values) for each leaf, to each task it can interact with
 allocate(sendbuf(max(1,nleaves*(2*nprofile+2)*min(nprocs-1,4))))
 nbuf  = 0
 nsend = 0
 do iproc = 0,nprocs-1
    if (iproc == id) cycle
    do k = 1,nleaves
       n = leaves(k)
       if (.not.leaf_may_interact(n,nprocs + iproc)) cycle
       call leaf_size_profile(n,prof,nlen)
       if (nbuf + 2*nlen + 2 > size(sendbuf)) call grow_real(sendbuf,nbuf + 2*nlen + 2)
       sendbuf(nbuf+1) = nodemap(n)
       sendbuf(nbuf+2) = nlen
       sendbuf(nbuf+3:nbuf+2+2*nlen) = reshape(prof(:,1:nlen),[2*nlen])
       nbuf = nbuf + 2*nlen + 2
       nsend(iproc+1) = nsend(iproc+1) + 2*nlen + 2
    enddo
 enddo
 call alltoall_mpi(nsend,nrecv)
 allocate(recvbuf(max(sum(nrecv),1)))
 call alltoallv_mpi(sendbuf,nsend,recvbuf,nrecv)

 ! store the profiles received
 ng = ncomb_tree - int(ncells)
 profile_len(1:ng) = -1
 if (allocated(profile_pool)) then
    if (size(profile_pool) < sum(nrecv)) deallocate(profile_pool)
 endif
 if (.not.allocated(profile_pool)) allocate(profile_pool(max(sum(nrecv),1)))
 irecv = 0
 do while (irecv < sum(nrecv))
    g    = nint(recvbuf(irecv+1))
    nlen = nint(recvbuf(irecv+2))
    profile_len(g)   = nlen
    profile_start(g) = irecv + 3
    profile_pool(irecv+3:irecv+2+2*nlen) = recvbuf(irecv+3:irecv+2+2*nlen)
    irecv = irecv + 2*nlen + 2
 enddo
 call mark_global_active(ng)

end subroutine exchange_size_profiles

!-----------------------------------------------------------------------
!+
!  whether a local node could interact with any node in the domain of
!  another task: the domain region bound, conservative
!+
!-----------------------------------------------------------------------
logical function leaf_may_interact(n,idomain)
 use kernel, only:radkern
 integer, intent(in) :: n,idomain
 real :: dx,dy,dz,xoff,yoff,zoff,r2,d,sizes,rcut

 call get_sep(nodeglobal(idomain)%xcen,node(n)%xcen,dx,dy,dz,xoff,yoff,zoff,r2)
 d     = sqrt(r2) - nodeglobal(idomain)%size
 sizes = 2.*nodeglobal(idomain)%size + node(n)%size
 rcut  = radkern*max(nodeglobal(idomain)%hmax,node(n)%hmax)
 leaf_may_interact = (tree_accuracy*d < sizes) .or. (d < sizes + rcut)

end function leaf_may_interact

!-----------------------------------------------------------------------
!+
!  mark the local nodes on the branches of active leaves
!+
!-----------------------------------------------------------------------
subroutine mark_active_branches()
 integer :: k,n

 let_stamp = let_stamp + 1
 do k = 1,nactive_leaves
    n = active_leaves(k)
    do while (branch_mark(n) /= let_stamp)
       branch_mark(n) = let_stamp
       if (n == irootnode .or. node(n)%parent <= 0) exit
       n = node(n)%parent
    enddo
 enddo

end subroutine mark_active_branches

!-----------------------------------------------------------------------
!+
!  mark the global nodes with active leaves below them: a leaf of another
!  task is active unless its profile says it has no active particles (or
!  none was received for a leaf too far away to matter, kept as active)
!+
!-----------------------------------------------------------------------
subroutine mark_global_active(ng)
 integer, intent(in) :: ng
 integer :: g,il,ir

 ! children come after their parents in the global tree, except above the
 ! domains, where the heap order also puts children at larger indices
 do g = ng,1,-1
    il = nodeglobal(g)%leftchild
    ir = nodeglobal(g)%rightchild
    if (il <= 0 .and. ir <= 0) then
       global_active(g) = (profile_len(g) /= 0)
    else
       global_active(g) = .false.
       if (il > 0) global_active(g) = global_active(il)
       if (ir > 0) global_active(g) = global_active(g) .or. global_active(ir)
    endif
 enddo

end subroutine mark_global_active

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
 if (branch_mark(n0) /= let_stamp) return
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
             if (branch_mark(ichild(j)) == let_stamp) call push_pair(istack,nstack,ichild(j),idepth+1,0)
          endif
       enddo
    endif
 enddo
 ! an active leaf found at a shallower depth stays in the walk at all later steps
 do k = 2,nprofile
    sleaf(:,k) = max(sleaf(:,k),sleaf(:,k-1))
 enddo
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

 nrec  = size(buf)/nlet_record
 nbase = ncomb
 if (nbase + nrec > size(node)) call grow_combined_tree(nbase + nrec)
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
       ipar = int(ncells) - ipar   ! a global node of iproc
    else
       ipar = nbase + ipar
    endif
    node(n)%parent = ipar
    if (nint(buf(ip+8)) == 1) then
       node(ipar)%leftchild = n
    else
       node(ipar)%rightchild = n
    endif
    ! the parent is no longer a leaf of the combined tree
    leaf_is_active(ipar)    = 0
    remote_leaf_owner(ipar) = 0
    leaf_is_active(n)       = 1
    remote_leaf_owner(n)    = iproc + 1
    remote_local_index(n)   = nint(buf(ip+9))
 enddo
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
 integer(kind=8) :: nsent,nrecv,np2p,nrecvmax

 nsent  = reduceall_mpi('+',nlet_sent)
 nrecv  = reduceall_mpi('+',nlet_recv)
 nrecvmax = reduceall_mpi('max',nlet_recv)
 np2p   = reduceall_mpi('+',nlet_p2p)
 if (id == master .and. nsent > 0) then
    write(iunit,"(a,i0,a,f10.3,a,f10.3)") ' MPI LET nodes imported (all tasks): ',nrecv,&
          '  MB per task (mean,max): ',1.e-6*8*nlet_record*real(nrecv)/nprocs,&
          ' ',1.e-6*8*nlet_record*real(nrecvmax)
    write(iunit,"(a,i0)") ' MPI LET leaf pairs with other tasks: ',np2p
 endif
 nlet_sent = 0
 nlet_recv = 0
 nlet_p2p  = 0

end subroutine print_let_stats

!-----------------------------------------------------------------------
!+
!  the trial neighbours of a cell sent by another task: the particles of
!  the given local leaves, which its dual tree walk paired with the cell
!+
!-----------------------------------------------------------------------
subroutine get_leaf_neighbours(xpos,leaves,nleaves,mylistneigh,nneigh,xyzcache,ixyzcachesize)
 real,    intent(in)  :: xpos(3)
 integer, intent(in)  :: leaves(:),nleaves,ixyzcachesize
 integer, intent(out) :: mylistneigh(:),nneigh
 real,    intent(out) :: xyzcache(:,:)
 integer :: k,maxcache
 real    :: dx,dy,dz,xoffset,yoffset,zoffset

 maxcache = 0
 if (ixyzcachesize > 0) maxcache = size(xyzcache,1)
 nneigh = 0
 do k = 1,nleaves
    call get_sep(xpos,node(leaves(k))%xcen,dx,dy,dz,xoffset,yoffset,zoffset)
    call cache_neighbours(nneigh,leaves(k),ixyzcachesize,maxcache,mylistneigh,xyzcache,xoffset,yoffset,zoffset)
 enddo

end subroutine get_leaf_neighbours

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
 integer, intent(in)  :: nactive
 integer, allocatable :: heapnode(:)
 real,    allocatable :: heapkey(:)
 integer :: nbudget,nheap,n

 ! copy a number of nodes proportional to the active particles: gathering and syncing
 ! them then costs a small fraction of the work on those particles, however many tasks
 ! there are, and deep in the timestep hierarchy (few active, the tree rebuilt often,
 ! few exports) only a few nodes are copied
 nbudget = min(mirror_budget(),max(nmirror_min,nactive/nactive_per_mirror))
 if (allocated(mirror_list)) then
    if (size(mirror_list) < nbudget) deallocate(mirror_list)
 endif
 if (.not.allocated(mirror_list)) allocate(mirror_list(nbudget))
 allocate(heapnode(nbudget),heapkey(nbudget))

 nmirror = 0
 nheap = 0
 call add_children(irootnode,heapnode,heapkey,nheap)
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
 use io,     only:nprocs,id
 use kernel, only:radkern
 integer, intent(in) :: n
 integer :: idomain
 real    :: dx,dy,dz,xoff,yoff,zoff,r,sizen,sized,reachn,reachd

 reaches_other_domain = .true.
 sizen  = node(n)%size
 reachn = sizen + radkern*node(n)%hmax
 do idomain = nprocs,2*nprocs-1
    if (idomain == nprocs + id) cycle
    call get_sep(node(n)%xcen,nodeglobal(idomain)%xcen,dx,dy,dz,xoff,yoff,zoff)
    r = sqrt(dx*dx + dy*dy + dz*dz)
    sized  = nodeglobal(idomain)%size
    reachd = sized + radkern*nodeglobal(idomain)%hmax
    if (r < reachn + reachd) return
    ! a cell of the other domain is at least r - sized away, and no larger than sized
    if (gravity .and. tree_accuracy*(r - sized) < sizen + sized) return
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
 do p = 1,nmirror
    nodemap(mirror_list(p)) = p
 enddo
 do p = 1,nmirror
    n = mirror_list(p)
    block(p) = node(n)
    block(p)%leftchild  = position_in_block(node(n)%leftchild)
    block(p)%rightchild = position_in_block(node(n)%rightchild)
    block(p)%parent     = position_in_block(node(n)%parent)
 enddo

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
 use io,       only:nprocs,id
 use mpiutils, only:allgatherv_mpi
 real, allocatable :: hmirror(:),hall(:)
 integer :: p,k,n,idomain,il,ir
 real    :: hchild

 allocate(hmirror(nmirror),hall(sum(nmirror_task)))
 do p = 1,nmirror
    hmirror(p) = node(mirror_list(p))%hmax
 enddo
 call allgatherv_mpi(hmirror,hall,nmirror_task,mirror_displ)
 do k = 1,size(hall)
    nodeglobal(2*nprocs-1+k)%hmax = hall(k)
 enddo

 do idomain = nprocs,2*nprocs-1
    il = nodeglobal(idomain)%leftchild
    ir = nodeglobal(idomain)%rightchild
    if (il == 0) cycle
    hchild = nodeglobal(il)%hmax
    if (ir > 0) hchild = max(hchild,nodeglobal(ir)%hmax)
    n = idomain
    nodeglobal(n)%hmax = max(nodeglobal(n)%hmax,hchild)
    do while (nodeglobal(n)%parent /= 0)
       n = nodeglobal(n)%parent
       nodeglobal(n)%hmax = max(nodeglobal(n)%hmax,hchild)
    enddo
 enddo

 ! the same values in the combined tree
 do k = 1,ncomb_tree - int(ncells)  ! the global tree only, not the imported nodes
    if (.not.is_this_task(k)) node(int(ncells) + k)%hmax = nodeglobal(k)%hmax
 enddo
 node(irootnode)%hmax = nodeglobal(nprocs + id)%hmax

end subroutine sync_hmax_mpi

!----------------------------------------------------------------!
!+
! - Neighbour finding core routines + gravity far field
!+
!----------------------------------------------------------------!

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
subroutine getneigh(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,&
                    get_hj,get_f,fnode,remote_export,nq)
 use io,       only:fatal,id
 use part,     only:gravity
 use kernel,   only:radkern
 type(kdnode), intent(in)  :: node(:) !ncellsmax+1)
 integer,      intent(in)  :: ixyzcachesize
 real,         intent(in)  :: xpos(3)
 real,         intent(in)  :: xsizei,rcuti
 integer,      intent(out) :: listneigh(:)
 integer,      intent(out) :: nneigh
 real,         intent(out) :: xyzcache(:,:)
 integer,      intent(in)  :: leaf_is_active(:)
 logical,      intent(in)  :: get_hj
 logical,      intent(in)  :: get_f
 real,         intent(out), optional :: fnode(lenfgrav)
 logical,      intent(out), optional :: remote_export(:)
 integer,      intent(in),  optional :: nq
 integer :: maxcache
 integer :: n,istack,il,ir
 integer :: nstack(maxdepth+1)
 real :: dx,dy,dz,xsizej,rcutj
 real :: rcut,rcut2,r2
 real :: xoffset,yoffset,zoffset,tree_acc2
 logical :: open_tree_node
 logical :: global_walk
#ifdef GRAVITY
 real :: quads(6)
 real :: dr,totmass_node
 logical :: count_gravity
#endif
 tree_acc2 = tree_accuracy*tree_accuracy
 if (get_f .and. .not.present(fnode)) then
    call fatal('getneigh','get_f but fnode not passed...')
 endif
 if (present(fnode)) fnode(:) = 0.
 rcut     = rcuti

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache,1)
 else
    maxcache = 0
 endif

 if (present(remote_export)) then
    remote_export = .false.
    global_walk = .true.
 else
    global_walk = .false.
 endif

 nneigh = 0
 istack = 1
 nstack(istack) = irootnode
 open_tree_node = .false.

 over_stack: do while(istack /= 0)
    n = nstack(istack)
    istack = istack - 1
    call get_sep(xpos,node(n)%xcen,dx,dy,dz,xoffset,yoffset,zoffset,r2)
    xsizej  = node(n)%size
    il      = node(n)%leftchild
    ir      = node(n)%rightchild
#ifdef GRAVITY
    totmass_node = node(n)%mass
    quads        = node(n)%quads
#endif

    if (get_hj) then  ! find neighbours within both hi and hj
       rcutj = radkern*node(n)%hmax
       rcut  = max(rcuti,rcutj)
    endif
    rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
    ! only open nodes on the gravity criterion when the walk computes gravity.
    if (get_f) open_tree_node = tree_acc2*r2 < (xsizei + xsizej)**2   ! tree opening criterion for self-gravity
    if_open_node: if ((r2 < rcut2) .or. open_tree_node) then
       if_leaf: if (leaf_is_active(n) /= 0) then ! once we hit a leaf node, retrieve contents into trial neighbour cache
          if_global_walk: if (global_walk) then
             ! id is stored in cellatid (passed through into leaf_is_active) as id + 1
             if (leaf_is_active(n) /= (id + 1)) then
                remote_export(leaf_is_active(n)) = .true.
             endif
          else
             call cache_neighbours(nneigh,n,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
          endif if_global_walk
       else
          if (istack+2 > maxdepth+1) call fatal('getneigh','stack overflow in getneigh')
          if (il /= 0) then
             istack = istack + 1
             nstack(istack) = il
          endif
          if (ir /= 0) then
             istack = istack + 1
             nstack(istack) = ir
          endif
       endif if_leaf
#ifdef GRAVITY
    elseif (get_f) then ! if_open_node
       ! When searching for neighbours of this node, the tree walk may encounter
       ! nodes on the global tree that it does not need to open, so it should
       ! just add the contribution to fnode. However, when walking a different
       ! part of the tree, it may then become necessary to export this node to
       ! a remote task. When it arrives at the remote task, it will then walk
       ! the remote tree.
       !
       ! The complication arises when tree refinment is enabled, which puts part
       ! of the remote tree onto the global tree. fnode will be double counted
       ! if a contribution is made on the global tree and a separate branch
       ! causes it to be sent to a remote task, where that contribution is
       ! counted again.
       !
       ! The solution is to not count the parts of the local tree that have been
       ! added onto the global tree (those with nodemap /= 0). In the global
       ! walk, n indexes the global tree, so nodemap must not be read.
       count_gravity = global_walk
       if (.not.global_walk) count_gravity = (nodemap(n) == 0)
       if_count_gravity: if (count_gravity) then
          !
          !--long range force on node due to distant node, along node centres
          !  along with derivatives in order to perform series expansion
          !
          dr = 1./sqrt(r2)
          call compute_M2L(dx,dy,dz,dr,totmass_node,quads,fnode)

       endif if_count_gravity
#endif

    endif if_open_node
 enddo over_stack

end subroutine getneigh

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search (SFMM version)
!  (all particles within a given h_i and optionally within h_j)
!  A dual tree walk is used to compute
!  every node-node interactions
!+
!----------------------------------------------------------------
subroutine getneigh_dual(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,leaf_is_active,&
                              get_hj,get_f,fnode,icell,nremote,ntrunc)
 type(kdnode), intent(inout) :: node(:) !ncellsmax+1)
 integer,      intent(in)    :: ixyzcachesize
 real,         intent(in)    :: xpos(3)
 real,         intent(in)    :: xsizei,rcuti
 integer,      intent(out)   :: listneigh(:)
 integer,      intent(out)   :: nneigh
 real,         intent(out)   :: xyzcache(:,:)
 integer,      intent(in)    :: leaf_is_active(:)
 logical,      intent(in)    :: get_hj
 logical,      intent(in)    :: get_f
 real,         intent(out)   :: fnode(lenfgrav)
 integer,      intent(in)    :: icell
 integer,      intent(out), optional :: nremote,ntrunc ! leaf pairs with other tasks' leaves, and with truncated ones
 integer :: istack,i,iparent,idstbranch,idst,isrc,maxcache,ibase,nodestate,nrem,ntr
 integer :: branch(maxdepth+1),nparents,stack(3,maxstacksize),startwith(2)
 real    :: dx,dy,dz,xoffset,yoffset,zoffset
 real    :: tree_acc2
 real    :: fnode_acc(lenfgrav)
 logical :: stackit

 tree_acc2 = tree_accuracy*tree_accuracy

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache,1)
 else
    maxcache = 0
 endif

 call get_list_of_parent_nodes(icell,node,branch,nparents,startwith)

 neighnodecount_branch(1:nparents) = 0
 fnode_branch(:,1:nparents) = 0.
 fnode_acc = 0.
 nneigh = 0
 istack = 0
 nrem = 0
 ntr = 0
 nremote_leaves = 0
 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 if (use_dualcache .and. startwith(2) > 0) then
    do i=1,neighnodecache_count(startwith(1))
       isrc = neighnodecache(neighnodecache_start(startwith(1)) + i)
       call open_nodes(stack,istack,node(isrc),isrc,branch,startwith(2),&
                       listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                       maxcache,xoffset,yoffset,zoffset,nrem,ntr)
    enddo
 else
    ! start from the top of the branch: the root of the local tree, or of the
    ! combined tree with more than one MPI task
    istack = istack + 1
    stack(1,istack) = branch(nparents)
    stack(2,istack) = branch(nparents)
    stack(3,istack) = nparents
 endif
!
!-- parallel select algorithm to check every interactions between the tree and the selected branch
!
 do while(istack > 0)
    !-- pop the stack
    idst       = stack(1,istack) ! dest node id
    isrc       = stack(2,istack) ! src node id
    idstbranch = stack(3,istack) ! dest id in branch array
    istack     = istack - 1

    if (idst == isrc) then !-- self interaction ignored (directly push onto stack)
       stackit = .true.
       xoffset = 0.
       yoffset = 0.
       zoffset = 0.
    else
       call node_interaction(idst,node(idst),node(isrc),tree_acc2,fnode_branch(:,idstbranch),stackit,xoffset,yoffset,zoffset)
    endif

    if (stackit) then
       neighnodecount_branch(idstbranch) = neighnodecount_branch(idstbranch) + 1
       !-- if count overflow, we will not cache it during the downward pass
       if (neighnodecount_branch(idstbranch) <= maxnodecache_local) then
          neighnode_branch(neighnodecount_branch(idstbranch),idstbranch) = isrc
       endif

       call open_nodes(stack,istack,node(isrc),isrc,branch,idstbranch,&
                       listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                       maxcache,xoffset,yoffset,zoffset,nrem,ntr)
    endif
 enddo
 if (present(nremote)) nremote = nrem
 if (present(ntrunc))  ntrunc  = ntr

 !
 !-- Downward pass to accumulate on each leaf / Cache and fetch fgrav optimisation
 !
 do i=nparents,2,-1 ! parents(1) is equal to icell
    iparent = branch(i)
    ! -- Cache node if first thread to reach it or fetch fnode in memory
    if (use_dualcache) then
       nodestate = cachestate(iparent) ! no atomic needed, see declaration of cachestate
       if (nodestate == 0) then ! first fence to avoid capture collision
          !$omp atomic capture
          nodestate = cachestate(iparent)
          cachestate(iparent) = max(cachestate(iparent),1)
          !$omp end atomic
          if (nodestate == 0) then ! if still the winner then cache
             !-- winner: publish fnode first ...
             fnodecache(1:lenfgrav,iparent) = fnode_branch(1:lenfgrav,i)
             ! flush: fnodecache must be visible before the state says it is cached
             !$omp flush
             !$omp atomic write
             cachestate(iparent) = 2
             !$omp end atomic

             !-- then store interaction list in the cache array if it fits
             if (neighnodecount_branch(i)> 0 .and. neighnodecount_branch(i) <= maxnodecache_local) then
                !$omp atomic capture
                ibase = itail_neigh
                itail_neigh = itail_neigh + neighnodecount_branch(i)
                !$omp end atomic
                if (ibase+neighnodecount_branch(i) <= size(neighnodecache)) then
                   neighnodecache(ibase+1:ibase+neighnodecount_branch(i)) = neighnode_branch(1:neighnodecount_branch(i),i)
                   neighnodecache_start(iparent) = ibase
                   neighnodecache_count(iparent) = neighnodecount_branch(i)
                   ! flush: interaction list must be visible before the state says it is cached
                   !$omp flush
                   !$omp atomic write
                   cachestate(iparent) = 3
                   !$omp end atomic
                endif
             endif
          endif
       elseif (nodestate>=2) then
          !-- fetch fnode from the cache array
          ! flush: see the fnodecache published before the state was set
          !$omp flush
          fnode_branch(1:lenfgrav,i) = fnodecache(1:lenfgrav,iparent)
       endif
    endif

    call get_sep(node(iparent)%xcen,node(branch(i-1))%xcen,dx,dy,dz,xoffset,yoffset,zoffset)
    fnode = fnode_acc + fnode_branch(:,i)
    call propagate_fnode_to_node(fnode_acc,fnode,dx,dy,dz)
 enddo

 fnode = fnode_acc + fnode_branch(:,1)

end subroutine getneigh_dual

!-----------------------------------------------------------
!+
!  get the separation in 3D between two nodes of the tree
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

!-----------------------------------------------------------
!+
!  get the size and rcut of two interacting nodes
!+
!-----------------------------------------------------------
pure subroutine get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)
 use kernel,   only:radkern
 type(kdnode), intent(in)  :: node_dst,node_src
 real,         intent(out) :: size_src,size_dst
 real,         intent(out) :: rcut_src,rcut_dst

 rcut_src = node_src%hmax*radkern
 rcut_dst = node_dst%hmax*radkern
 size_src = node_src%size
 size_dst = node_dst%size

end subroutine get_node_size

!-----------------------------------------------------------
!+
!  return list of parents of current node
!+
!-----------------------------------------------------------
subroutine get_list_of_parent_nodes(inode,node,parents,nparents,startwith)
 integer,      intent(in)  :: inode
 type(kdnode), intent(in)  :: node(:)
 integer,      intent(out) :: parents(:)
 integer,      intent(out) :: nparents
 integer,      intent(out) :: startwith(2)
 integer :: j,nodestate

 j = inode
 nparents  = 1
 parents   = 0
 startwith = 0
 parents(nparents) = j ! set first elem to inode to use parents for propagation
 do while (node(j)%parent  /=  0)
    j = node(j)%parent
    nparents = nparents + 1
    parents(nparents) = j
    nodestate = cachestate(j) ! no atomic needed, see declaration of cachestate
    if (nodestate==3 .and. startwith(2)==0) then
       ! deepest fully-cached ancestor: candidate pruned start
       startwith(1) = j
       startwith(2) = nparents
    elseif (nodestate<2 .and. startwith(2)/=0) then
       ! if non-cached ancestor node on the pruned branch reset the pruned start
       startwith = 0
    endif
 enddo
 ! flush: the walk will read the cached interaction list published before state 3 was set
 if (startwith(2) > 0) then
    !$omp flush
 endif

end subroutine get_list_of_parent_nodes

!-----------------------------------------------------------
!+
!  Compute node node gravity interactions
!+
!-----------------------------------------------------------
subroutine open_nodes(stack,istack,srcnode,isrc,branch,idstbranch,&
                           listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                           maxcache,xoffset,yoffset,zoffset,nremote,ntrunc)
 use io, only:fatal
 type(kdnode), intent(in)    :: srcnode
 integer,      intent(in)    :: isrc,idstbranch
 integer,      intent(in)    :: branch(:)
 integer,      intent(in)    :: ixyzcachesize,maxcache
 integer,      intent(in)    :: leaf_is_active(:)
 integer,      intent(inout) :: listneigh(:)
 integer,      intent(inout) :: nneigh
 integer,      intent(inout) :: stack(:,:),istack
 real,         intent(inout) :: xyzcache(:,:)
 real,         intent(in)    :: xoffset,yoffset,zoffset
 integer,      intent(inout) :: nremote,ntrunc
 integer :: ir,il,ibranchnext,idstnext
 logical :: isdstleaf

 il = srcnode%leftchild
 ir = srcnode%rightchild

 !-- find the new dst id to push onto the stack
 if (idstbranch-1>0) then !-- if not leaf
    ibranchnext = idstbranch-1
    isdstleaf   = .false.
 else
    ibranchnext = idstbranch ! leaf lowering if upper leaf
    isdstleaf   = .true.
 endif

 idstnext = branch(ibranchnext) ! new dest node id

 is_src_leaf: if (leaf_is_active(isrc) /= 0) then
    ! a leaf of another task with a positive local index is not a leaf on its owner:
    ! its children are missing, so it can be neither lowered against nor paired here
    if (remote_leaf_owner(isrc) /= 0) then
       if (remote_local_index(isrc) > 0) ntrunc = ntrunc + 1
    endif
    is_P2P: if (isdstleaf) then !-- P2P detected should be cached and tagged as neighbours
       if (remote_leaf_owner(isrc) /= 0) then
          ! a leaf of another task: its particles are not here, so the cell is sent to its owner
          nremote = nremote + 1
          call add_remote_leaf(isrc)
       else
          call cache_neighbours(nneigh,isrc,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
       endif
    else ! then you're a leaf -> leaf lowering
       if (istack+1 > maxstacksize) call fatal('getneigh','stack overflow in getneigh')
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = isrc
       stack(3,istack) = ibranchnext
    endif is_P2P
 else
    if (il /= 0) then
       if (istack+1 > maxstacksize) call fatal('getneigh','stack overflow in getneigh')
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = il
       stack(3,istack) = ibranchnext
    endif
    if (ir /= 0) then
       if (istack+1 > maxstacksize) call fatal('getneigh','stack overflow in getneigh')
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = ir
       stack(3,istack) = ibranchnext
    endif
 endif is_src_leaf

end subroutine open_nodes

!-----------------------------------------------------------
!+
!  add a leaf of another task to this thread's list
!+
!-----------------------------------------------------------
subroutine add_remote_leaf(n)
 integer, intent(in) :: n
 integer, allocatable :: tmp(:)

 if (.not.allocated(remote_leaves)) allocate(remote_leaves(256))
 if (nremote_leaves >= size(remote_leaves)) then
    allocate(tmp(2*size(remote_leaves)))
    tmp(1:nremote_leaves) = remote_leaves(1:nremote_leaves)
    call move_alloc(tmp,remote_leaves)
 endif
 nremote_leaves = nremote_leaves + 1
 remote_leaves(nremote_leaves) = n

end subroutine add_remote_leaf

!----------------------------------------------------------------
!+
!  Cache particles within identified neighbour nodes
!+
!----------------------------------------------------------------
subroutine cache_neighbours(nneigh,isrc,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
 use part, only:rho,gradh,treecache
 use dim,  only:igradomega,igradzeta,maxpsph
#ifdef GRAVITY
 use dim,  only:igradsoft
#endif
 integer, intent(in)    :: isrc,ixyzcachesize,maxcache
 real,    intent(in)    :: xoffset,yoffset,zoffset
 integer, intent(inout) :: nneigh
 integer, intent(out)   :: listneigh(:)
 real,    intent(out)   :: xyzcache(:,:)
 integer :: npnode,ipart,num_to_cache,ip,inode

 npnode = inoderange(2,isrc) - inoderange(1,isrc) + 1

 if (nneigh + npnode <= ixyzcachesize) then
    num_to_cache = npnode
 elseif (nneigh < ixyzcachesize) then
    num_to_cache = ixyzcachesize - nneigh
 else
    num_to_cache = 0
 endif

 if (num_to_cache > 0) then
    do ipart=1,num_to_cache
       inode = inoderange(1,isrc)+ipart-1
       ip = abs(inodeparts(inode))
       listneigh(nneigh+ipart)  = ip
       xyzcache(1,nneigh+ipart) = treecache(1,inode) + xoffset
       xyzcache(2,nneigh+ipart) = treecache(2,inode) + yoffset
       xyzcache(3,nneigh+ipart) = treecache(3,inode) + zoffset
       if (maxcache >= 4) then
          xyzcache(ih1,nneigh+ipart) = 1./treecache(4,inode)
       endif
       if (maxcache >= 5) then
          xyzcache(im,nneigh+ipart) = treecache(5,inode)
       endif
       if (maxcache >= 7) then
          if (ip <= maxpsph) then
             xyzcache(irho,nneigh+ipart)        = rho(ip)
             xyzcache(izetaomega,nneigh+ipart) = real(gradh(igradzeta,ip))*real(gradh(igradomega,ip))
          else
             xyzcache(irho,nneigh+ipart)        = 0.
             xyzcache(izetaomega,nneigh+ipart) = 0.
          endif
       endif
#ifdef GRAVITY
       if (maxcache >= 8) then
          if (ip <= maxpsph) then
             xyzcache(isoftomega,nneigh+ipart) = real(gradh(igradsoft,ip))*real(gradh(igradomega,ip))
          else
             xyzcache(isoftomega,nneigh+ipart) = 0.
          endif
       endif
#endif
    enddo
 endif

 if (num_to_cache < npnode) then
    do ipart=num_to_cache+1,npnode
       listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,isrc)+ipart-1))
    enddo
 endif

 nneigh = nneigh + npnode

end subroutine cache_neighbours

!-----------------------------------------------------------
!+
!  Test the separation between the node pair and compute
!  the interaction if needed
!+
!-----------------------------------------------------------
subroutine node_interaction(idst,node_dst,node_src,tree_acc2,fnode,stackit,xoffset,yoffset,zoffset)
 type(kdnode), intent(in)    :: node_dst,node_src
 integer,      intent(in)    :: idst
 real,         intent(in)    :: tree_acc2
 real,         intent(inout) :: fnode(lenfgrav)
 real,         intent(out)   :: xoffset,yoffset,zoffset
 logical,      intent(out)   :: stackit
 real    :: dx,dy,dz,r2
 real    :: rcut_dst,rcut_src,rcut,rcut2
 real    :: size_dst,size_src
 logical :: wellsep
 integer :: dststate
#ifdef GRAVITY
 real    :: dr1
#endif

 call get_sep(node_dst%xcen,node_src%xcen,dx,dy,dz,xoffset,yoffset,zoffset,r2)
 call get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)

 if (use_dualcache) then
    dststate = cachestate(idst) ! no atomic needed, see declaration of cachestate
 else
    dststate = 0
 endif

 rcut  = max(rcut_dst,rcut_src)
 rcut2 = (size_dst+size_src+rcut)**2
 wellsep = (tree_acc2*r2 > (size_dst+size_src)**2) .and. (r2 > rcut2)

 if (wellsep) then
#ifdef GRAVITY
    if (dststate<2) then
       dr1 = 1./sqrt(r2)
       call compute_M2L(dx,dy,dz,dr1,node_src%mass,node_src%quads,fnode)
       call add_torque_correction(dx,dy,dz,dr1,node_dst%mass,node_src%mass, &
                                  node_dst%octs,node_src%octs,fnode)
    endif
#endif
    stackit = .false.
 else
    stackit = .true.
 endif

end subroutine node_interaction

!-----------------------------------------------------------
!+
!  Compute the Taylor expansion coeffs between the node
!  centres using the quadrupole moments (p=3) (Dehnen 2002)
!+
!-----------------------------------------------------------
pure subroutine compute_M2L(dx,dy,dz,dr1,q0,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr1,q0
 real, intent(in)    :: quads(6)
 real, intent(inout) :: fnode(lenfgrav)
 real :: qxx,qxy,qxz,qyy,qyz,qzz,dx2,dx3,dy2,dy3,dz2,dz3
 real :: dr12,D3(10),D2(6),D1(3),g0,g1,g2,g3,g2dx,g2dy,g2dz

! note: dr == 1/sqrt(r2)
 dr12 = dr1*dr1
 dx2  = dx*dx
 dx3  = dx*dx2
 dy2  = dy*dy
 dy3  = dy*dy2
 dz2  = dz*dz
 dz3  = dz*dz2
 ! be careful with the sign of your Green's function, it can mess up everything.
 ! We switched multiple signs here to match the Phantom sign convention
 g0   =  dr1
 g1   = -1.*dr12*g0
 g2   = -3.*dr12*g1
 g3   = -5.*dr12*g2
 g2dx = g2 * dx
 g2dy = g2 * dy
 g2dz = g2 * dz

 !D1, D2, D3 verified and agree with shamrock to float precision
 D3(1)  = 3. * g2dx + g3 * dx3    ! xxx
 D3(2)  = g2dy + g3 * dx2 * dy    ! xxy
 D3(3)  = g2dz + g3 * dx2 * dz    ! xxz
 D3(4)  = g2dx + g3 * dy2 * dx    ! xyy
 D3(5)  = g3 * dx * dy * dz       ! xyz
 D3(6)  = g2dx + g3 * dz2 * dx    ! xzz
 D3(7)  = 3. * g2dy + g3 * dy3    ! yyy
 D3(8)  = g2dz + g3 * dy2 * dz    ! yyz
 D3(9)  = g2dy + g3 * dz2 * dy    ! yzz
 D3(10) = 3. * g2dz + g3 * dz3    ! zzz

 D2(1)  = g1 + g2 * dx2 ! xx
 D2(2)  = g2dx * dy     ! xy
 D2(3)  = g2dx * dz     ! xz
 D2(4)  = g1 + g2 * dy2 ! yy
 D2(5)  = g2dy * dz     ! yz
 D2(6)  = g1 + g2 * dz2 ! zz

 D1(1)  = g1*dx
 D1(2)  = g1*dy
 D1(3)  = g1*dz

 qxx = quads(1)
 qxy = quads(2)
 qxz = quads(3)
 qyy = quads(4)
 qyz = quads(5)
 qzz = quads(6)

 fnode(1)  = fnode(1)  + D1(1)*q0 + 0.5*(D3(1)*qxx + 2.*(D3(2)*qxy + D3(3)*qxz + D3(5)*qyz) + D3(4)*qyy + D3(6)*qzz)    ! C¹_x
 fnode(2)  = fnode(2)  + D1(2)*q0 + 0.5*(D3(2)*qxx + 2.*(D3(4)*qxy + D3(5)*qxz + D3(8)*qyz) + D3(7)*qyy + D3(9)*qzz)    ! C¹_y
 fnode(3)  = fnode(3)  + D1(3)*q0 + 0.5*(D3(3)*qxx + 2.*(D3(5)*qxy + D3(6)*qxz + D3(9)*qyz) + D3(8)*qyy + D3(10)*qzz)   ! C¹_z
 fnode(4)  = fnode(4)  - (D2(1) * q0)  ! C²_xx
 fnode(5)  = fnode(5)  - (D2(2) * q0)  ! C²_xy
 fnode(6)  = fnode(6)  - (D2(3) * q0)  ! C²_xz
 fnode(7)  = fnode(7)  - (D2(4) * q0)  ! C²_yy
 fnode(8)  = fnode(8)  - (D2(5) * q0)  ! C²_yz
 fnode(9)  = fnode(9)  - (D2(6) * q0)  ! C²_zz
 fnode(10) = fnode(10) + D3(1) * q0    ! C³_xxx
 fnode(11) = fnode(11) + D3(2) * q0    ! C³_xxy
 fnode(12) = fnode(12) + D3(3) * q0    ! C³_xxz
 fnode(13) = fnode(13) + D3(4) * q0    ! C³_xyy
 fnode(14) = fnode(14) + D3(5) * q0    ! C³_xyz
 fnode(15) = fnode(15) + D3(6) * q0    ! C³_xzz
 fnode(16) = fnode(16) + D3(7) * q0    ! C³_yyy
 fnode(17) = fnode(17) + D3(8) * q0    ! C³_yyz
 fnode(18) = fnode(18) + D3(9) * q0    ! C³_yzz
 fnode(19) = fnode(19) + D3(10)* q0    ! C³_zzz
 fnode(20) = fnode(20) + g0*q0 + 0.5*(D2(1)*qxx + D2(4)*qyy + D2(6)*qzz + 2*(D2(2)*qxy + D2(3)*qxz + D2(5)*qyz))! C⁰ (potential)

end subroutine compute_M2L

#ifdef GRAVITY
!----------------------------------------------------------------
!+
!  Marcello (2017) TCO torque correction: add a constant
!  acceleration Fc/M_dst to the destination cell so the net
!  cell-cell torque vanishes, while keeping equal-and-opposite
!  forces. Uses the pruned D'_ijkl contraction (his Eq. 19).
!+
!----------------------------------------------------------------
pure subroutine add_torque_correction(dx,dy,dz,dr1,mass_dst,mass_src,octs_dst,octs_src,fnode)
 real, intent(in)    :: dx,dy,dz,dr1,mass_dst,mass_src
 real, intent(in)    :: octs_dst(10),octs_src(10)
 real, intent(inout) :: fnode(lenfgrav)
 real :: s(10)
 real :: sxkk,sykk,szkk,sxrr,syrr,szrr
 real :: r5i,r7i,fac

 if (mass_dst <= 0. .or. mass_src <= 0.) return

 ! S_jkl = M_dst,jkl * M_src - M_dst * M_src,jkl
 s(:) = octs_dst*mass_src - mass_dst*octs_src

 ! traces S_i,kk
 sxkk = s(1) + s(4) + s(6)
 sykk = s(2) + s(7) + s(9)
 szkk = s(3) + s(8) + s(10)

 ! S_iab R_a R_b  (R is the node separation; even in R so dest-src is fine)
 sxrr = s(1)*dx*dx + s(4)*dy*dy + s(6)*dz*dz + 2.*(s(2)*dx*dy + s(3)*dx*dz + s(5)*dy*dz)
 syrr = s(2)*dx*dx + s(7)*dy*dy + s(9)*dz*dz + 2.*(s(4)*dx*dy + s(5)*dx*dz + s(8)*dy*dz)
 szrr = s(3)*dx*dx + s(8)*dy*dy + s(10)*dz*dz + 2.*(s(5)*dx*dy + s(6)*dx*dz + s(9)*dy*dz)

 r5i = dr1**5
 r7i = r5i*dr1*dr1
 ! Appendix Eq. 38 at P=3 uses 1/(n!(P-n)!) = 1/3!, not the 1/2 of Eq. 15.
 ! Combined with the D' contraction this is 3/2 rather than 9/2.
 fac = 1.5

 ! Fc_i = (3/2) (S_ikk/R^5 - 5 S_iab R_a R_b / R^7); add Fc/M_dst
 fnode(1) = fnode(1) - fac*(sxkk*r5i - 5.*sxrr*r7i)/mass_dst
 fnode(2) = fnode(2) - fac*(sykk*r5i - 5.*syrr*r7i)/mass_dst
 fnode(3) = fnode(3) - fac*(szkk*r5i - 5.*szrr*r7i)/mass_dst

end subroutine add_torque_correction
#endif

!-----------------------------------------------------------
!+
!  Taylor expand the contribution from direct parent nodes
!  to the child node centre
!+
!-----------------------------------------------------------
pure subroutine propagate_fnode_to_node(fnode,fnode_sup,dx,dy,dz)
 real, intent(in)  :: fnode_sup(lenfgrav),dx,dy,dz
 real, intent(out) :: fnode(lenfgrav)

 fnode(1)  = fnode_sup(1) + dx*(fnode_sup(4) + 0.5*(dx*fnode_sup(10) + dy*fnode_sup(11) +dz*fnode_sup(12)))& ! xx +0.5(xxx+xxy+xxz)
                          + dy*(fnode_sup(5) + 0.5*(dx*fnode_sup(11) + dy*fnode_sup(13) +dz*fnode_sup(14)))& ! xy +0.5(xxy+xyy+xyz)
                          + dz*(fnode_sup(6) + 0.5*(dx*fnode_sup(12) + dy*fnode_sup(14) +dz*fnode_sup(15)))  ! xz +0.5(xxz+xyz+xzz)
 fnode(2)  = fnode_sup(2) + dx*(fnode_sup(5) + 0.5*(dx*fnode_sup(11) + dy*fnode_sup(13) +dz*fnode_sup(14)))& ! xy +0.5(xxy+xyy+xyz)
                          + dy*(fnode_sup(7) + 0.5*(dx*fnode_sup(13) + dy*fnode_sup(16) +dz*fnode_sup(17)))& ! yy +0.5(xyy+yyy+yyz)
                          + dz*(fnode_sup(8) + 0.5*(dx*fnode_sup(14) + dy*fnode_sup(17) +dz*fnode_sup(18)))  ! yz +0.5(xyz+yyz+yyz)
 fnode(3)  = fnode_sup(3) + dx*(fnode_sup(6) + 0.5*(dx*fnode_sup(12) + dy*fnode_sup(14) +dz*fnode_sup(15)))& ! xz +0.5(xxz+xyz+xzz)
                          + dy*(fnode_sup(8) + 0.5*(dx*fnode_sup(14) + dy*fnode_sup(17) +dz*fnode_sup(18)))& ! yz +0.5(xyz+yyz+yzz)
                          + dz*(fnode_sup(9) + 0.5*(dx*fnode_sup(15) + dy*fnode_sup(18) +dz*fnode_sup(19)))  ! zz +0.5(xzz+yzz+zzz)
 fnode(4)  = fnode_sup(4) + dx*fnode_sup(10) + dy*fnode_sup(11) + dz*fnode_sup(12)                           ! xxx + xxy + xxz
 fnode(5)  = fnode_sup(5) + dx*fnode_sup(11) + dy*fnode_sup(13) + dz*fnode_sup(14)                           ! xxy + xyy + xyz
 fnode(6)  = fnode_sup(6) + dx*fnode_sup(12) + dy*fnode_sup(14) + dz*fnode_sup(15)                           ! xxz + xyz + xzz
 fnode(7)  = fnode_sup(7) + dx*fnode_sup(13) + dy*fnode_sup(16) + dz*fnode_sup(17)                           ! xyy + yyy + yyz
 fnode(8)  = fnode_sup(8) + dx*fnode_sup(14) + dy*fnode_sup(17) + dz*fnode_sup(18)                           ! xyz + yyz + yzz
 fnode(9)  = fnode_sup(9) + dx*fnode_sup(15) + dy*fnode_sup(18) + dz*fnode_sup(19)                           ! xzz + yzz + zzz
 fnode(10) = fnode_sup(10)
 fnode(11) = fnode_sup(11)
 fnode(12) = fnode_sup(12)
 fnode(13) = fnode_sup(13)
 fnode(14) = fnode_sup(14)
 fnode(15) = fnode_sup(15)
 fnode(16) = fnode_sup(16)
 fnode(17) = fnode_sup(17)
 fnode(18) = fnode_sup(18)
 fnode(19) = fnode_sup(19)
 fnode(20) = fnode_sup(20) - dx*(fnode_sup(1)+0.5*(dx*(fnode_sup(4)+(1./3.)*(dx*fnode_sup(10)+dy*fnode_sup(11)+dz*fnode_sup(12)))+&
                                                   dy*(fnode_sup(5)+(1./3.)*(dx*fnode_sup(11)+dy*fnode_sup(13)+dz*fnode_sup(14)))+&
                                                   dz*(fnode_sup(6)+(1./3.)*(dx*fnode_sup(12)+dy*fnode_sup(14)+dz*fnode_sup(15)))))&
                           - dy*(fnode_sup(2)+0.5*(dx*(fnode_sup(5)+(1./3.)*(dx*fnode_sup(11)+dy*fnode_sup(13)+dz*fnode_sup(14)))+&
                                                   dy*(fnode_sup(7)+(1./3.)*(dx*fnode_sup(13)+dy*fnode_sup(16)+dz*fnode_sup(17)))+&
                                                   dz*(fnode_sup(8)+(1./3.)*(dx*fnode_sup(14)+dy*fnode_sup(17)+dz*fnode_sup(18)))))&
                           - dz*(fnode_sup(3)+0.5*(dx*(fnode_sup(6)+(1./3.)*(dx*fnode_sup(12)+dy*fnode_sup(14)+dz*fnode_sup(15)))+&
                                                   dy*(fnode_sup(8)+(1./3.)*(dx*fnode_sup(14)+dy*fnode_sup(17)+dz*fnode_sup(18)))+&
                                                   dz*(fnode_sup(9)+(1./3.)*(dx*fnode_sup(15)+dy*fnode_sup(18)+dz*fnode_sup(19)))))

end subroutine propagate_fnode_to_node
!----------------------------------------------------------------
!+
!  Internal subroutine to compute the Taylor-series expansion
!  of the gravitational force, given the force acting on the
!  centre of the node and its derivatives
!
! INPUT:
!   fnode: array containing force on node due to distant nodes
!          and first derivatives of f (i.e. Jacobian matrix)
!          and second derivatives of f (i.e. Hessian matrix)
!   dx,dy,dz: offset of the particle from the node centre of mass
!
! OUTPUT:
!   fxi,fyi,fzi : gravitational force at the new position
!+
!----------------------------------------------------------------
pure subroutine expand_fgrav_in_taylor_series(fnode,dx,dy,dz,fxi,fyi,fzi,poti)
 real, intent(in)  :: fnode(lenfgrav)
 real, intent(in)  :: dx,dy,dz
 real, intent(out) :: fxi,fyi,fzi,poti
 real :: dfxx,dfxy,dfxz,dfyy,dfyz,dfzz
 real :: d2fxxx,d2fxxy,d2fxxz,d2fxyy,d2fxyz,d2fxzz,d2fyyy,d2fyyz,d2fyzz,d2fzzz

 fxi = fnode(1)
 fyi = fnode(2)
 fzi = fnode(3)
 dfxx = fnode(4)
 dfxy = fnode(5)
 dfxz = fnode(6)
 dfyy = fnode(7)
 dfyz = fnode(8)
 dfzz = fnode(9)
 d2fxxx = fnode(10)
 d2fxxy = fnode(11)
 d2fxxz = fnode(12)
 d2fxyy = fnode(13)
 d2fxyz = fnode(14)
 d2fxzz = fnode(15)
 d2fyyy = fnode(16)
 d2fyyz = fnode(17)
 d2fyzz = fnode(18)
 d2fzzz = fnode(19)
 poti = fnode(20)

 fxi = fxi   + dx*(dfxx + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
             + dy*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
             + dz*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
 fyi = fyi   + dx*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
             + dy*(dfyy + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
             + dz*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
 fzi = fzi   + dx*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
             + dy*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
             + dz*(dfzz + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))
 ! Minus sign here as we are shifted of 1 in the (-1)^k compared to force
 poti = poti - dx*(fnode(1) + 0.5*(dx*dfxx + dy*dfxy + dz*dfxz)) &
             - dy*(fnode(2) + 0.5*(dx*dfxy + dy*dfyy + dz*dfyz)) &
             - dz*(fnode(3) + 0.5*(dx*dfxz + dy*dfyz + dz*dfzz))

end subroutine expand_fgrav_in_taylor_series

end module neighkdtree
