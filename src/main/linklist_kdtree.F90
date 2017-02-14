!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: linklist
!
!  DESCRIPTION:
!  This module contains all routines required for
!  link-list based neighbour-finding
!
!  THIS VERSION USES A K-D TREE
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    tree_accuracy -- tree opening criterion (0.0-1.0)
!
!  DEPENDENCIES: boundary, dim, infile_utils, io, kdtree, kernel, mpiutils,
!    part
!+
!--------------------------------------------------------------------------
module linklist
 use dim,     only:maxp,ncellsmax
 use part,    only:ll
 implicit none
 character(len=80), parameter, public :: &  ! module version
    modid="$Id$"

 integer, public :: ifirstincell(ncellsmax+1)
 !real(kind=4), dimension(ncellsmax+1), public :: hmaxcell
 integer(kind=8), public :: ncells
 real, public            :: dxcell
 real, public :: dcellx = 0.,dcelly = 0.,dcellz = 0.

 public :: set_linklist, get_neighbour_list, write_inopts_link, read_inopts_link

 public :: set_hmaxcell,get_hmaxcell,update_hmax_remote

private

contains

subroutine get_hmaxcell(inode,hmaxcell)
 use kdtree, only:node
 integer, intent(in)  :: inode
 real,    intent(out) :: hmaxcell

 hmaxcell = node(inode)%hmax

end subroutine get_hmaxcell

subroutine set_hmaxcell(inode,hmaxcell)
!!!$ use omputils, only:ipart_omp_lock,nlockgrp
 use kdtree, only:node
 integer, intent(in) :: inode
 real,    intent(in) :: hmaxcell
 integer :: n

 n = inode
 node(n)%hmax = hmaxcell

 ! walk tree up
 do while (node(n)%parent  /=  0)
   n = node(n)%parent
!$omp critical (hmax)
   node(n)%hmax = max(node(n)%hmax, hmaxcell)
!$omp end critical (hmax)
 enddo

end subroutine set_hmaxcell

subroutine update_hmax_remote(ncells)
 use kdtree,   only:node
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

subroutine set_linklist(npart,nactive,xyzh,vxyzu)
 use kdtree, only:maketree,revtree,ndimtree,maxlevel,maxlevel_indexed
! use timing, only:print_time
 integer, intent(in)    :: npart,nactive
 real,    intent(inout) :: xyzh(4,maxp)
 real,    intent(in)    :: vxyzu(:,:)
 integer, save :: naccum = 0
! real(kind=4) :: t1,t2

! call cpu_time(t1)
 !
 ! make the tree if more than 10% of the particles are active,
 ! or when the accumulated total number of active particles
 ! during revtree calls exceeds the number of particlees
 !
 if (10*nactive > npart .or. naccum==0 .or. naccum > 2*npart &
     .or. maxlevel > maxlevel_indexed) then
    !print*,' MAKETREE ',nactive,naccum
    call maketree(xyzh,vxyzu,npart,ndimtree,ifirstincell,ncells)
    naccum = npart
!    print *, 'maketree'
 else
    naccum = naccum + nactive
!    print *,' REVTREE ',nactive,naccum
    call maketree(xyzh,vxyzu,npart,ndimtree,ifirstincell,ncells)
!    call revtree(xyzh, ifirstincell, ncells)
 endif
! call cpu_time(t2)
! call print_time(t2-t1)

end subroutine set_linklist

!-----------------------------------------------------------------------
!+
! Using the k-d tree, compiles the neighbour list for the
! current cell (this list is common to all particles in the cell)
!
! the list is returned in 'listneigh' (length nneigh)
!+
!-----------------------------------------------------------------------
subroutine get_neighbour_list(inode,listneigh,nneigh,xyzh,xyzcache,ixyzcachesize, &
                              activeonly,hmaxold,inactiveonly,getj,f)
 use kdtree, only:getneigh,node,lenfgrav
 use kernel, only:radkern
#ifdef PERIODIC
 use io,       only:warning
 use boundary, only:dxbound,dybound,dzbound
#endif
 integer, intent(in)  :: inode,ixyzcachesize
 integer, intent(out) :: listneigh(:)
 integer, intent(out) :: nneigh
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: xyzcache(:,:)
 logical, intent(in)  :: activeonly
 logical, intent(in),  optional :: inactiveonly,getj
 real,    intent(in),  optional :: hmaxold
 real,    intent(out), optional :: f(lenfgrav)
 real :: xpos(3)
 real :: fgrav(lenfgrav)
 real :: xsizei,rcuti
 logical :: get_j
!
!--retrieve geometric centre of the node and the search radius (e.g. 2*hmax)
!
 xpos    = node(inode)%xcen(1:3)
 xsizei  = node(inode)%size
 rcuti   = radkern*node(inode)%hmax ! 2h

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

 if (present(f)) then
    call getneigh(xpos,xsizei,rcuti,3,listneigh,nneigh,xyzh,xyzcache,ixyzcachesize,&
               ifirstincell,ll,get_j,fgrav)
    f = fgrav
 else
    call getneigh(xpos,xsizei,rcuti,3,listneigh,nneigh,xyzh,xyzcache,ixyzcachesize,&
               ifirstincell,ll,get_j)
 endif

 return
end subroutine get_neighbour_list

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

end module linklist
