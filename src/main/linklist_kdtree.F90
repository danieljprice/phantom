!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module linklist
!
! This module contains all routines required for
!  link-list based neighbour-finding
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
 use part,         only:ll
 use dtypekdtree,  only:kdnode
 implicit none

 integer,               allocatable :: cellatid(:)
 integer,     public,   allocatable :: ifirstincell(:)
 type(kdnode),          allocatable :: nodeglobal(:)
 type(kdnode), public,  allocatable :: node(:)
 integer,               allocatable :: nodemap(:)
 integer,      public ,            allocatable :: listneigh(:)
 integer,      public ,            allocatable :: listneigh_global(:)
!$omp threadprivate(listneigh)
 integer(kind=8), public :: ncells
 real, public            :: dxcell
 real, public :: dcellx = 0.,dcelly = 0.,dcellz = 0.

 integer              :: globallevel,refinelevels

 public :: allocate_linklist, deallocate_linklist
 public :: set_linklist, get_neighbour_list, write_inopts_link, read_inopts_link
 public :: get_distance_from_centre_of_mass, getneigh_pos
 public :: set_hmaxcell,get_hmaxcell,update_hmax_remote
 public :: get_cell_location
 public :: sync_hmax_mpi

 private

contains

subroutine allocate_linklist
 use allocutils, only:allocate_array
 use kdtree,     only:allocate_kdtree
 use dim,        only:maxp

 call allocate_array('cellatid',     cellatid,     ncellsmaxglobal+1 )
 call allocate_array('ifirstincell', ifirstincell, ncellsmax+1       )
 call allocate_array('nodeglobal',   nodeglobal,   ncellsmaxglobal+1 )
 call allocate_array('node',         node,         ncellsmax+1       )
 call allocate_array('nodemap',      nodemap,      ncellsmax+1       )
 call allocate_kdtree()
!$omp parallel
 call allocate_array('listneigh',listneigh,maxp)
!$omp end parallel
 call allocate_array('listneigh_global',listneigh_global,maxp)

end subroutine allocate_linklist

subroutine deallocate_linklist
 use kdtree,   only:deallocate_kdtree

 if (allocated(cellatid)) deallocate(cellatid)
 if (allocated(ifirstincell)) deallocate(ifirstincell)
 if (allocated(nodeglobal)) deallocate(nodeglobal)
 if (allocated(node)) deallocate(node)
 if (allocated(nodemap)) deallocate(nodemap)
!$omp parallel
 if (allocated(listneigh)) deallocate(listneigh)
!$omp end parallel
 if (allocated(listneigh_global)) deallocate(listneigh_global)
 call deallocate_kdtree()

end subroutine deallocate_linklist

subroutine get_hmaxcell(inode,hmaxcell)
 integer, intent(in)  :: inode
 real,    intent(out) :: hmaxcell

 hmaxcell = node(inode)%hmax

end subroutine get_hmaxcell

subroutine set_hmaxcell(inode,hmaxcell)
!!!$ use omputils, only:ipart_omp_lock,nlockgrp
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

subroutine update_hmax_remote(ncells)
 use mpiutils, only:reduceall_mpi
 integer(kind=8), intent(in) :: ncells
 integer :: n,j
 real :: hmaxcell

 ! could do only active cells by checking ifirstincell >= 0
 do n=1,int(ncells)
    if (ifirstincell(n) > 0) then
       ! reduce on leaf nodes
       node(n)%hmax = reduceall_mpi('max',node(n)%hmax)

       ! propagate up local tree from active cells
       hmaxcell = node(n)%hmax
       j = n
       do while (node(j)%parent  /=  0)
          j = node(j)%parent
          node(j)%hmax = max(node(j)%hmax, hmaxcell)
       enddo
    endif
 enddo

end subroutine update_hmax_remote

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

subroutine set_linklist(npart,nactive,xyzh,vxyzu,for_apr)
 use io,           only:nprocs
 use dtypekdtree,  only:ndimtree
 use kdtree,       only:maketree,maketreeglobal
 use dim,          only:mpi,use_sinktree
 use part,         only:nptmass,xyzmh_ptmass
 integer,           intent(inout) :: npart
 integer,           intent(in)    :: nactive
 real,              intent(inout) :: xyzh(:,:)
 real,              intent(in)    :: vxyzu(:,:)
 logical, optional, intent(in)    :: for_apr
 logical :: apr_tree

 apr_tree = .false.
 if (present(for_apr)) apr_tree = for_apr

 if (mpi .and. nprocs > 1) then
    if (use_sinktree) then
       call maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,npart,ndimtree,cellatid,ifirstincell,ncells,&
                           apr_tree,nptmass,xyzmh_ptmass)
    else
       call maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,npart,ndimtree,cellatid,ifirstincell,ncells,&
                           apr_tree)
    endif
 else
    if (use_sinktree) then
       call maketree(node,xyzh,npart,ndimtree,ifirstincell,ncells,apr_tree,nptmass=nptmass,xyzmh_ptmass=xyzmh_ptmass)
    else
       call maketree(node,xyzh,npart,ndimtree,ifirstincell,ncells,apr_tree)
    endif
 endif

end subroutine set_linklist

!-----------------------------------------------------------------------
!+
! Using the k-d tree, compiles the neighbour list for the
! current cell (this list is common to all particles in the cell)
!
! the list is returned in 'listneigh' (length nneigh)
!+
!-----------------------------------------------------------------------
subroutine get_neighbour_list(inode,mylistneigh,nneigh,xyzh,xyzcache,ixyzcachesize, &
                              getj,f,remote_export, &
                              cell_xpos,cell_xsizei,cell_rcuti)
 use io,     only:nprocs
 use dim,    only:mpi
 use kdtree, only:getneigh,lenfgrav
 use kernel, only:radkern
 use part,   only:gravity
#ifdef PERIODIC
 use io,       only:warning
 use boundary, only:dxbound,dybound,dzbound
#endif
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

#ifdef PERIODIC
 if (rcuti > 0.5*min(dxbound,dybound,dzbound)) then
    call warning('get_neighbour_list', '2h > 0.5*L in periodic neighb. '//&
                'search: USE HIGHER RES, BIGGER BOX or LOWER MINPART IN TREE')
 endif
#endif
 !
 !--perform top-down tree walk to find all particles within radkern*h
 !  and force due to node-node interactions
 !
 ! print*,' searching for cell ',icell,' xpos = ',xpos(:),' hmax = ',rcut/radkern
 get_j = .false.
 if (present(getj)) get_j = getj

 get_f = (gravity .and. present(f))

 if (mpi .and. global_search) then
    ! Find MPI tasks that have neighbours of this cell, output to remote_export
    call getneigh(nodeglobal,xpos,xsizei,rcuti,3,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
            cellatid,get_j,get_f,fgrav_global,remote_export)
 elseif (get_f) then
    ! Set fgrav to zero, which matters if gravity is enabled but global search is not
    fgrav_global = 0.0
 endif

 ! Find neighbours of this cell on this node
 call getneigh(node,xpos,xsizei,rcuti,3,mylistneigh,nneigh,xyzcache,ixyzcachesize,&
              ifirstincell,get_j,get_f,fgrav)

 if (get_f) f = fgrav + fgrav_global

end subroutine get_neighbour_list

subroutine getneigh_pos(xpos,xsizei,rcuti,ndim,mylistneigh,nneigh,xyzh,xyzcache,ixyzcachesize,ifirstincell)
 use kdtree, only:getneigh
 integer, intent(in)  :: ndim,ixyzcachesize
 real,    intent(in)  :: xpos(ndim)
 real,    intent(in)  :: xsizei,rcuti
 integer, intent(out) :: mylistneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: xyzcache(:,:)
 integer, intent(in)  :: ifirstincell(:) !ncellsmax+1)

 call getneigh(node,xpos,xsizei,rcuti,ndim,mylistneigh,nneigh,xyzcache,ixyzcachesize, &
               ifirstincell,.false.,.false.)

end subroutine getneigh_pos

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_inopts_link(iunit)
 use kdtree,       only:tree_accuracy
 use infile_utils, only:write_inopt
 use part,         only:gravity
 integer, intent(in) :: iunit

 if (gravity) then
    call write_inopt(tree_accuracy,'tree_accuracy','tree opening criterion (0.0-1.0)',iunit)
 endif

end subroutine write_inopts_link

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_inopts_link(name,valstring,imatch,igotall,ierr)
 use kdtree, only:tree_accuracy
 use part,   only:gravity
 use io,     only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0

 imatch = .false.
 igotall = .true.
 ierr = 0
 if (gravity) then
    imatch = .true.
    igotall = .false.
    select case(trim(name))
    case('tree_accuracy')
       read(valstring,*,iostat=ierr) tree_accuracy
       ngot = ngot + 1
       if ((tree_accuracy < 0. .or. tree_accuracy > 1.0)) &
          call fatal('read_inopts_kdtree','tree accuracy out of range (0.0-1.0)')
    case default
       imatch = .false.
    end select
    if (ngot >= 1) igotall = .true.
 endif

end subroutine read_inopts_link

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

end module linklist
