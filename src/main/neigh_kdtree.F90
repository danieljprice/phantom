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
!   kdtree, kernel, mpiutils, part
!
 use dim,          only:ncellsmax,ncellsmaxglobal,gravity
 use dtypekdtree,  only:kdnode,lenfgrav
 use kdtree,       only:inoderange,inodeparts,irootnode,maxdepth,irefine,tree_accuracy,&
                        ih1,im,irho,izetaomega,isoftomega
 implicit none

 integer,               allocatable :: cellatid(:)
 integer,               allocatable :: nodemap(:)
 type(kdnode),          allocatable :: nodeglobal(:)
 type(kdnode), public,  allocatable :: node(:)
 integer,      public,  allocatable :: leaf_is_active(:) ! : 0 internal node or empty cell, : 1 active cell, :- inactive cell
 integer,      public,  allocatable :: active_leaves(:)  ! the cells with leaf_is_active > 0, in order (set by build_tree)
 integer,      public               :: nactive_leaves = 0
 integer,      public , allocatable :: listneigh(:)
 integer,      public , allocatable :: listneigh_global(:)
!-- dual tree cache arrays
 integer,               allocatable :: cachestate(:)
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
 integer(kind=8), public            :: ncells
 real,            public            :: dxcell
 real,            public            :: dcellx = 0.,dcelly = 0.,dcellz = 0.
 logical,         public            :: use_dualtree  = .true.
 logical,         public            :: use_dualcache = gravity ! only debug flag / should be on with gravity
 integer,         parameter         :: maxnodecache_local = 512
 integer,         parameter         :: maxneigh_per_node  = 16
 integer,         parameter         :: maxstacksize = 2048
 integer                            :: globallevel,refinelevels

 public :: allocate_neigh, deallocate_neigh
 public :: build_tree, get_neighbour_list, write_options_tree, read_options_tree
 public :: get_distance_from_centre_of_mass, getneigh_pos
 public :: set_hmaxcell,get_hmaxcell
 public :: get_cell_location
 public :: sync_hmax_mpi
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

 call allocate_array('cellatid',       cellatid,       ncellsmaxglobal+1 )
 call allocate_array('leaf_is_active', leaf_is_active, ncellsmax+1       )
 call allocate_array('active_leaves',  active_leaves,  ncellsmax+1       )
 call allocate_array('nodeglobal',     nodeglobal,     ncellsmaxglobal+1 )
 call allocate_array('node',           node,           ncellsmax+1       )
 call allocate_array('nodemap',        nodemap,        ncellsmax+1       )
 call allocate_kdtree()
 call allocate_array('listneigh_global',listneigh_global,maxp)
 call allocate_array('cachestate', cachestate, ncellsmax+1)

 if (use_dualcache) then
    call allocate_array('fnodecache', fnodecache, lenfgrav, ncellsmax+1)
    call allocate_array('neighnodecache',neighnodecache,ncellsmax*maxneigh_per_node)
    call allocate_array('neighnodecache_start',neighnodecache_start,ncellsmax+1)
    call allocate_array('neighnodecache_count',neighnodecache_count,ncellsmax+1)
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
 node(n)%hmax = hmaxcell

 ! walk tree up, stopping at the first ancestor whose hmax already covers hmaxcell:
 ! a node's hmax is never below its children's, so neither is any of its ancestors'.
 ! Other threads update the same ancestors, so both the test and the update are atomic
 do while (node(n)%parent /= 0)
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
 use kdtree,       only:maketree,maketreeglobal!,revtree
 use dim,          only:mpi,use_sinktree
 use part,         only:nptmass,xyzmh_ptmass,maxp
 use allocutils,   only:allocate_array
 integer, intent(inout) :: npart
 integer, intent(in)    :: nactive
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(in)    :: vxyzu(:,:)
 logical, intent(in), optional :: for_apr
 logical :: apr_tree

 apr_tree = .false.
 if (present(for_apr)) apr_tree = for_apr

 !
 ! the listneigh array is threadprivate, but if the thread numbers or ids are changed
 ! then the memory might be lost. So the following lines are a failsafe
 ! to ensure that the listneigh array is always allocated for each thread
 !
 !$omp parallel
 if (.not. allocated(listneigh)) call allocate_array('listneigh',listneigh,maxp)
 !$omp end parallel

 if (mpi .and. nprocs > 1) then
    if (use_sinktree) then
       call maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,npart,cellatid,leaf_is_active,ncells,&
                           apr_tree,nptmass,xyzmh_ptmass)
    else
       call maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,npart,cellatid,leaf_is_active,ncells,&
                           apr_tree)
    endif
 else
    if (use_sinktree) then
       call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree,nptmass=nptmass,xyzmh_ptmass=xyzmh_ptmass)
    else
       ! use revtree for small numbers of active particles to avoid tree rebuild overhead
       ! threshold: use revtree if < 0.1% of total particles
       !if (npart > 0 .and. nactive < 0.001*npart) then
       !   call revtree(node,xyzh,leaf_is_active,ncells)
       !else
       call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree)
       !endif
    endif
 endif
 call list_active_leaves()

 ! the dual tree walk cache refers to the old nodes: reset it for the new tree
 if (use_dualcache) then
    cachestate(1:ncells) = 0 ! 0=untouched, 1=claimed, 2=fnode cached, 3=fully cached
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
 use io,       only:nprocs,warning
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

 if (present(remote_export)) then
    if (nprocs > 1) global_search = .true.
    remote_export = .false.
 else
    global_search = .false.
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

 if (mpi .and. global_search) then ! no sym fmm for now...
    ! Find MPI tasks that have neighbours of this cell, output to remote_export
    call getneigh(nodeglobal,xpos,xsizei,rcuti,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
                  cellatid,get_j,get_f,fgrav_global,remote_export)
 elseif (get_f) then
    ! Set fgrav to zero, which matters if gravity is enabled but global search is not
    fgrav_global = 0.0
 endif

 ! Find neighbours of this cell on this node
 if (get_f .and. .not.(mpi) .and. use_dualtree) then
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
!  sync the hmax values across all MPI tasks
!+
!-----------------------------------------------------------------------
subroutine sync_hmax_mpi
 use mpiutils,  only:reduceall_mpi
 use io,        only:nprocs
 integer :: i, n
 real    :: hmax(2**(globallevel+refinelevels+1)-1)

 hmax(:) = 0.0
 ! copy hmax values into contiguous array
 do i = 2,2**(refinelevels+1)-1
    hmax(nodemap(i)) = node(i)%hmax
 enddo

 ! reduce across threads
 hmax = reduceall_mpi('max', hmax)

 ! put values back into node
 do i = 2*nprocs,2**(globallevel+refinelevels+1)-1
    nodeglobal(i)%hmax = hmax(i)
 enddo

 ! walk tree up
 do i = 2*nprocs,4*nprocs
    n = i
    do while (nodeglobal(n)%parent /= 0)
       n = nodeglobal(n)%parent
       nodeglobal(n)%hmax = max(nodeglobal(n)%hmax, hmax(i))
    enddo
 enddo

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
    if (gravity) open_tree_node = tree_acc2*r2 < (xsizei + xsizej)**2   ! tree opening criterion for self-gravity
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
       ! added onto the global tree.

       count_gravity: if ( global_walk .or. (n > irefine) ) then
          !
          !--long range force on node due to distant node, along node centres
          !  along with derivatives in order to perform series expansion
          !
          dr = 1./sqrt(r2)
          call compute_M2L(dx,dy,dz,dr,totmass_node,quads,fnode)

       endif count_gravity
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
                              get_hj,get_f,fnode,icell)
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
 integer :: istack,i,iparent,idstbranch,idst,isrc,maxcache,ibase,nodestate
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
 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 if (use_dualcache .and. startwith(2) > 0) then
    do i=1,neighnodecache_count(startwith(1))
       isrc = neighnodecache(neighnodecache_start(startwith(1)) + i)
       call open_nodes(stack,istack,node(isrc),isrc,branch,startwith(2),&
                       listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                       maxcache,xoffset,yoffset,zoffset)
    enddo
 else
    istack = istack + 1
    stack(1,istack) = irootnode
    stack(2,istack) = irootnode
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
                       maxcache,xoffset,yoffset,zoffset)
    endif
 enddo

 !
 !-- Downward pass to accumulate on each leaf / Cache and fetch fgrav optimisation
 !
 do i=nparents,2,-1 ! parents(1) is equal to icell
    iparent = branch(i)
    ! -- Cache node if first thread to reach it or fetch fnode in memory
    if (use_dualcache) then
       ! acquire: if the state says cached, the fnodecache written before it is visible
       !$omp atomic read acquire
       nodestate = cachestate(iparent)
       !$omp end atomic
       if (nodestate == 0) then ! first fence to avoid capture collision
          !$omp atomic capture
          nodestate = cachestate(iparent)
          cachestate(iparent) = max(cachestate(iparent),1)
          !$omp end atomic
          if (nodestate == 0) then ! if still the winner then cache
             !-- winner: publish fnode first ...
             fnodecache(1:lenfgrav,iparent) = fnode_branch(1:lenfgrav,i)
             ! release: fnodecache must be visible before the state says it is cached
             !$omp atomic write release
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
                   ! release: interaction list must be visible before the state says it is cached
                   !$omp atomic write release
                   cachestate(iparent) = 3
                   !$omp end atomic
                endif
             endif
          endif
       elseif (nodestate>=2) then
          !-- fetch fnode from the cache array
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
    ! acquire: state 3 means the walk will read this node's cached interaction list
    !$omp atomic read acquire
    nodestate = cachestate(j)
    !$omp end atomic
    if (nodestate==3 .and. startwith(2)==0) then
       ! deepest fully-cached ancestor: candidate pruned start
       startwith(1) = j
       startwith(2) = nparents
    elseif (nodestate<2 .and. startwith(2)/=0) then
       ! if non-cached ancestor node on the pruned branch reset the pruned start
       startwith = 0
    endif
 enddo

end subroutine get_list_of_parent_nodes

!-----------------------------------------------------------
!+
!  Compute node node gravity interactions
!+
!-----------------------------------------------------------
subroutine open_nodes(stack,istack,srcnode,isrc,branch,idstbranch,&
                           listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                           maxcache,xoffset,yoffset,zoffset)
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
    is_P2P: if (isdstleaf) then !-- P2P detected should be cached and tagged as neighbours
       call cache_neighbours(nneigh,isrc,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
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
    !$omp atomic read
    dststate = cachestate(idst)
    !$omp end atomic
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
