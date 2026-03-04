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
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - tree_accuracy : *tree opening criterion (0.0-1.0)*
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, infile_utils, io,
!   kdtree, kernel, mpiutils, part
!
 use dim,          only:ncellsmax,ncellsmaxglobal
 use dtypekdtree,  only:kdnode
 implicit none

 integer,               allocatable :: cellatid(:)
 integer,               allocatable :: nodemap(:)
 type(kdnode),          allocatable :: nodeglobal(:)
 type(kdnode), public,  allocatable :: node(:)
 integer,      public,  allocatable :: leaf_is_active(:) ! : 0 internal node or empty cell, : 1 active cell, :- inactive cell
 integer,      public , allocatable :: listneigh(:)
 integer,      public , allocatable :: listneigh_global(:)
!$omp threadprivate(listneigh)
 integer(kind=8), public            :: ncells
 real, public                       :: dxcell
 real, public                       :: dcellx = 0.,dcelly = 0.,dcellz = 0.
 logical, public                    :: use_dualtree = .true.
 integer                            :: globallevel,refinelevels

 public :: allocate_neigh, deallocate_neigh
 public :: build_tree, get_neighbour_list, write_options_tree, read_options_tree
 public :: get_distance_from_centre_of_mass, getneigh_pos
 public :: set_hmaxcell,get_hmaxcell
 public :: get_cell_location
 public :: sync_hmax_mpi

 private

contains

!-----------------------------------------------------------------------
!+
!  allocate memory for the neighbour list
!+
!-----------------------------------------------------------------------
subroutine allocate_neigh
 use allocutils, only:allocate_array
 use kdtree,     only:allocate_kdtree
 use dim,        only:maxp

 call allocate_array('cellatid',       cellatid,       ncellsmaxglobal+1 )
 call allocate_array('leaf_is_active', leaf_is_active, ncellsmax+1       )
 call allocate_array('nodeglobal',     nodeglobal,     ncellsmaxglobal+1 )
 call allocate_array('node',           node,           ncellsmax+1       )
 call allocate_array('nodemap',        nodemap,        ncellsmax+1       )
 call allocate_kdtree()
!$omp parallel
 call allocate_array('listneigh',listneigh,maxp)
!$omp end parallel
 call allocate_array('listneigh_global',listneigh_global,maxp)

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
 if (allocated(nodeglobal)) deallocate(nodeglobal)
 if (allocated(node)) deallocate(node)
 if (allocated(nodemap)) deallocate(nodemap)
!$omp parallel
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

 n = inode
 node(n)%hmax = hmaxcell

 ! walk tree up
 do while (node(n)%parent /= 0)
    n = node(n)%parent
!$omp critical (crit_node_hmax)
    node(n)%hmax = max(node(n)%hmax, hmaxcell)
!$omp end critical (crit_node_hmax)
 enddo

end subroutine set_hmaxcell

!-----------------------------------------------------------------------
!+
!  get the distance from the centre of mass of a cell
!+
!-----------------------------------------------------------------------
subroutine get_distance_from_centre_of_mass(inode,xi,yi,zi,dx,dy,dz,xcen)
 integer,   intent(in)           :: inode
 real,      intent(in)           :: xi,yi,zi
 real,      intent(out)          :: dx,dy,dz
 real,      intent(in), optional :: xcen(3)

 if (present(xcen)) then
    dx = xi - xcen(1)
    dy = yi - xcen(2)
    dz = zi - xcen(3)
 else
    dx = xi - node(inode)%xcen(1)
    dy = yi - node(inode)%xcen(2)
    dz = zi - node(inode)%xcen(3)
 endif

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
 integer,           intent(inout) :: npart
 integer,           intent(in)    :: nactive
 real,              intent(inout) :: xyzh(:,:)
 real,              intent(in)    :: vxyzu(:,:)
 logical, optional, intent(in)    :: for_apr
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

end subroutine build_tree

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
 use kdtree,   only:getneigh,getneigh_dual,lenfgrav
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
 use kdtree, only:getneigh
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
 integer,            intent(in)     :: inode
 real,               intent(out)    :: xpos(3)
 real,               intent(out)    :: xsizei
 real,               intent(out)    :: rcuti

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

end module neighkdtree
