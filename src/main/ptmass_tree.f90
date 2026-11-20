!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ptmass_tree
!
! Module that contains all the routines necessary to build a tree on the ptmass particles.
! This one is used to search efficiently ptmass parts in the accretion routine.
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, dim, dtypekdtree, io
!
 use dtypekdtree, only:ptmasstree
 use dim,        only:maxptmass,nnodeptmassmax
 implicit none

 public :: ptmasstree,build_ptmass_tree,get_ptmass_neigh,nfastacc
 public :: allocate_ptmasstree,deallocate_ptmasstree,ptmasskdtree

 private
 type(ptmasstree)     :: ptmasskdtree
 integer, parameter   :: nmaxleaf   = 2
 integer, parameter   :: iroot      = 1
 integer, parameter   :: istacksize = maxptmass ! enough for the deepest tree (one sink split off per level)

 integer, parameter   :: nfastacc   = 100

contains

subroutine allocate_ptmasstree
 use allocutils, only:allocate_array

 call allocate_array('iptmassnode', ptmasskdtree%iptmassnode, maxptmass)
 call allocate_array('ptmassnode',  ptmasskdtree%nodes,       nnodeptmassmax)

end subroutine allocate_ptmasstree

subroutine deallocate_ptmasstree

 if (allocated(ptmasskdtree%iptmassnode)) deallocate(ptmasskdtree%iptmassnode)
 if (allocated(ptmasskdtree%nodes)) deallocate(ptmasskdtree%nodes)

end subroutine deallocate_ptmasstree
!-----------------------------------------------------------------
!+
! Build ptmass KD-tree (serial: there are at most maxptmass sinks)
!+
!-----------------------------------------------------------------
subroutine build_ptmass_tree(xyzmh_ptmass,nptmass,tree,hacc_min)
 use io,   only:fatal
 use part, only:ihacc
 real,             intent(in)    :: xyzmh_ptmass(:,:)
 integer,          intent(in)    :: nptmass
 type(ptmasstree), intent(inout) :: tree
 real,             intent(in), optional :: hacc_min ! sinks with a smaller hacc cannot accrete (stars)
 integer :: stack(istacksize)
 integer :: i,iaxis,inode,lchild,rchild
 integer :: istart,iend,imed,npnode
 integer :: nnodes,istack
 real    :: xmin(3),xmax(3)
 real    :: dx(3),dr2,r2max
 real    :: xtmp(3),xcen(3)
 real    :: xpivot,mtot

 !-- initialize tree part index array with the alive sinks that can accrete
 npnode = 0
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       if (present(hacc_min)) then
          if (xyzmh_ptmass(ihacc,i) < hacc_min) cycle
       endif
       npnode = npnode + 1
       tree%iptmassnode(npnode) = i
    endif
 enddo

 !-- fill the root inode info
 nnodes               = 1
 tree%nodes(1)%istart = 1
 tree%nodes(1)%iend   = npnode
 tree%nodes(1)%parent = 0
 tree%nodes(1)%lchild = 0
 tree%nodes(1)%rchild = 0
 tree%nodes(1)%xcen   = 0.
 tree%nodes(1)%size   = 0.

 istack = 0
 if (npnode > 0) then
    istack = 1
    stack(istack) = iroot
 endif

 do while (istack > 0)
    inode  = stack(istack)
    istack = istack - 1

    istart = tree%nodes(inode)%istart
    iend   = tree%nodes(inode)%iend
    npnode = iend - istart + 1

    !-- compute bounding box and center for this inode
    xmin =  huge(0.)
    xmax = -huge(0.)
    xcen = 0.
    mtot = 0.

    do i=istart,iend
       xtmp = xyzmh_ptmass(1:3, tree%iptmassnode(i))
       xmin = min(xmin, xtmp)
       xmax = max(xmax, xtmp)
       xcen = xcen + xtmp*xyzmh_ptmass(4, tree%iptmassnode(i))
       mtot = mtot + xyzmh_ptmass(4, tree%iptmassnode(i))
    enddo

    xcen  = xcen/mtot
    r2max = 0.
    do i=istart,iend
       dx    = xyzmh_ptmass(1:3,tree%iptmassnode(i)) - xcen
       dr2   = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)
       r2max = max(r2max,dr2)
    enddo
    tree%nodes(inode)%xcen = xcen
    tree%nodes(inode)%size = sqrt(r2max) + epsilon(r2max)

    !-- leaves are not connected to previous tree builds
    tree%nodes(inode)%lchild = 0
    tree%nodes(inode)%rchild = 0

    !-- split the inode if it holds more than nmaxleaf sinks
    if (npnode > nmaxleaf) then
       !-- choose split dimension: longest dimension
       iaxis  = maxloc(xmax - xmin,1)
       xpivot = xcen(iaxis)

       !-- sort the indices in this inode by chosen dimension
       call sort_tree_ptmass_id(xyzmh_ptmass,tree%iptmassnode,istart,iend,iaxis,xpivot,imed)

       !-- all sinks on one side of the pivot (e.g. same position): split the index range
       if (imed < istart .or. imed >= iend) imed = (istart + iend)/2

       !-- create two child nodes (left: istart..imed, right: imed+1..iend)
       lchild = nnodes + 1
       rchild = nnodes + 2
       nnodes = nnodes + 2
       if (rchild > nnodeptmassmax) call fatal("ptmass_tree","node array overflow...",ival=rchild)

       !-- assign children to this inode
       tree%nodes(inode)%lchild = lchild
       tree%nodes(inode)%rchild = rchild

       !-- fill child metadata
       tree%nodes(lchild)%istart = istart
       tree%nodes(lchild)%iend   = imed
       tree%nodes(lchild)%parent = inode
       tree%nodes(rchild)%istart = imed + 1
       tree%nodes(rchild)%iend   = iend
       tree%nodes(rchild)%parent = inode

       if (istack + 2 > istacksize) call fatal("ptmass_tree","stack overflow...",ival=istack+2)
       stack(istack+1) = lchild
       stack(istack+2) = rchild
       istack = istack + 2
    endif
 enddo
 tree%nnodes = nnodes

end subroutine build_ptmass_tree

!-----------------------------------------------------------------
!+
! Sort a subarray of idx(start:endp) by coordinate.
! A simple quicksort with the pivot as the parents CoM coord along
! the split axis.
!+
!-----------------------------------------------------------------
subroutine sort_tree_ptmass_id(xyzmh_ptmass,iptmassnode,il,ir,iaxis,xpivot,imed)
 use io, only:fatal
 real,    intent(in)    :: xyzmh_ptmass(:, :),xpivot
 integer, intent(inout) :: iptmassnode(:)
 integer, intent(in)    :: il,ir,iaxis
 integer, intent(out)   :: imed
 integer :: i,j,id_swap
 logical :: i_inf_pivot,j_inf_pivot

 i = il
 j = ir
 i_inf_pivot = xyzmh_ptmass(iaxis,iptmassnode(i)) <= xpivot
 j_inf_pivot = xyzmh_ptmass(iaxis,iptmassnode(j)) <= xpivot

 do while (i < j)
    if (i_inf_pivot) then !-- good position if inf to pivot -> cycle i
       i = i + 1
       i_inf_pivot = xyzmh_ptmass(iaxis,iptmassnode(i)) <= xpivot
    else
       if (.not.j_inf_pivot) then !-- good position if sup to pivot -> cycle j
          j = j - 1
          j_inf_pivot = xyzmh_ptmass(iaxis,iptmassnode(j)) <= xpivot
       else !-- if i sup and j inf -> swap i and j and cycle both
          id_swap        = iptmassnode(i)
          iptmassnode(i) = iptmassnode(j)
          iptmassnode(j) = id_swap
          i = i + 1
          j = j - 1
          i_inf_pivot = xyzmh_ptmass(iaxis,iptmassnode(i)) <= xpivot
          j_inf_pivot = xyzmh_ptmass(iaxis,iptmassnode(j)) <= xpivot
       endif
    endif
 enddo

 if (.not.i_inf_pivot) i = i - 1
 if (j_inf_pivot)      j = j + 1

 imed = i

 if ( j /= i+1) call fatal("ptmass_tree","error in sort idx",ival=(j-i))

end subroutine sort_tree_ptmass_id

!-----------------------------------------------------------------
!+
! Neighbour ptmasses search using a top-down algorithm
!+
!-----------------------------------------------------------------
subroutine get_ptmass_neigh(tree,xpos,rsearch,listneigh,nneigh)
 use io, only:fatal
 type(ptmasstree), intent(in)  :: tree
 real,             intent(in)  :: xpos(3)
 real,             intent(in)  :: rsearch
 integer,          intent(out) :: listneigh(:)
 integer,          intent(out) :: nneigh

 integer              :: stack(istacksize)
 integer              :: istack
 integer              :: inode,i,j,k,istart,iend,il,ir,npnode
 real                 :: dr2,dx,dy,dz,size
 logical              :: isleaf,isneigh

 nneigh = 0

 istack = 1
 stack(istack) = iroot

 do while (istack > 0)

    inode = stack(istack)
    istack = istack - 1

    !-- compute squared distance from xpos to inode's max sphere

    dx   = xpos(1) - tree%nodes(inode)%xcen(1)
    dy   = xpos(2) - tree%nodes(inode)%xcen(2)
    dz   = xpos(3) - tree%nodes(inode)%xcen(3)
    dr2  = dx*dx+dy*dy+dz*dz
    size = tree%nodes(inode)%size

    il     = tree%nodes(inode)%lchild
    ir     = tree%nodes(inode)%rchild
    istart = tree%nodes(inode)%istart
    iend   = tree%nodes(inode)%iend

    isneigh = dr2 < (rsearch + size)**2
    isleaf  = (il==0) .and. (ir==0)

    if (isneigh) then
       npnode = iend - istart + 1
       if (npnode > 0) then
          if (isleaf) then
             if (nneigh + npnode > maxptmass) call fatal("ptmass_neigh","search found more sinks than can exist",ival=nneigh+npnode)
             do i = istart, iend
                nneigh = nneigh + 1
                listneigh(nneigh) = tree%iptmassnode(i)
             enddo
          else
             !-- push children (if exist)
             if (istack + 2 <= istacksize) then
                if (il /= 0) then
                   istack = istack + 1
                   stack(istack) = il
                endif
                if (ir /= 0) then
                   istack = istack + 1
                   stack(istack) = ir
                endif
             else
                call fatal("ptmass_neigh","stack overflow in tree search",ival=istack)
             endif
          endif
       endif
    endif
 enddo

 !-- sort by sink index, as in the direct loop over the sinks (the first unconditional match wins)
 do i = 2, nneigh
    j = listneigh(i)
    k = i - 1
    do while (k > 0)
       if (listneigh(k) <= j) exit
       listneigh(k+1) = listneigh(k)
       k = k - 1
    enddo
    listneigh(k+1) = j
 enddo

end subroutine get_ptmass_neigh

end module ptmass_tree

