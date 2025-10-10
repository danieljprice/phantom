!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module ptmass_tree
!
! Module that contains all the routines necessary to build a tree on the ptmass particles.
! This one is used to search efficiently ptmass part in the accretion routine...
!
 use dtypekdtree, only:ptmasstree
 implicit none

 public :: ptmasstree,build_ptmass_tree,get_ptmass_neigh
 public :: allocate_ptmasstree,deallocate_ptmasstree,ptmasskdtree

 private
 type(ptmasstree)     :: ptmasskdtree
 integer, parameter   :: nmaxleaf   = 2
 integer, parameter   :: iroot      = 1
 integer, parameter   :: istacksize = 512

contains

subroutine allocate_ptmasstree
 use dim,        only:maxptmass,nnodeptmassmax
 use allocutils, only:allocate_array

 call allocate_array('iptmassnode', ptmasskdtree%iptmassnode, maxptmass)
 call allocate_array('ptmassnode',  ptmasskdtree%nodes,       nnodeptmassmax)

end subroutine allocate_ptmasstree


subroutine deallocate_ptmasstree

 if (allocated(ptmasskdtree%iptmassnode)) deallocate(ptmasskdtree%iptmassnode)
 if (allocated(ptmasskdtree%nodes)) deallocate(ptmasskdtree%nodes)

end subroutine deallocate_ptmasstree
 !-----------------------------------------------------------------
 ! Build KD-tree (non-recursive construction using a stack)
 ! x(3,N) : coordinates (1..N)
 ! tree    : returned tree (allocated inside)
 !-----------------------------------------------------------------
subroutine build_ptmass_tree(xyzmh_ptmass,nptmass,tree)
 use io, only:fatal
 real,             intent(in)    :: xyzmh_ptmass(:,:)
 integer,          intent(in)    :: nptmass
 type(ptmasstree), intent(inout) :: tree

 integer, allocatable :: stack(:)
 integer :: i,iaxis,inode,lchild,rchild
 integer :: istart,iend,imed,npnode
 integer :: istack
 real    :: xmin(3), xmax(3)
 real    :: dx(3),dr2,r2max
 integer :: nnodes
 real    :: xtmp(3),xcen(3)
 real    :: xpivot

 allocate(stack(istacksize))

 npnode = 0

 !-- initialize tree part index array
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       npnode = npnode + 1
       tree%iptmassnode(npnode) = i
    endif
 enddo

 !-- fill the root inode level info
 nnodes               = 1
 tree%nodes(1)%istart = 1
 tree%nodes(1)%iend   = npnode
 tree%nodes(1)%parent = 0


 !-- push root in the stack
 istack =  1
 stack(istack) = iroot

 do while (istack > 0)
    inode = stack(istack)
    istack = istack - 1

    istart = tree%nodes(inode)%istart
    iend   = tree%nodes(inode)%iend
    npnode = iend - istart + 1

    !-- compute bounding box and center for this inode
    xmin =  huge(0.)
    xmax = -huge(0.)

    do i=istart,iend
       xtmp = xyzmh_ptmass(1:3, tree%iptmassnode(i))
       xmin(1) = min(xmin(1), xtmp(1))
       xmin(2) = min(xmin(2), xtmp(2))
       xmin(3) = min(xmin(3), xtmp(3))
       xmax(1) = max(xmax(1), xtmp(1))
       xmax(2) = max(xmax(2), xtmp(2))
       xmax(3) = max(xmax(3), xtmp(3))
    enddo


    xcen  = 0.5*(xmin + xmax)
    r2max = 0.
    do i=istart,iend
       dx(1) = xyzmh_ptmass(1,tree%iptmassnode(i)) - xcen(1)
       dx(2) = xyzmh_ptmass(2,tree%iptmassnode(i)) - xcen(2)
       dx(3) = xyzmh_ptmass(3,tree%iptmassnode(i)) - xcen(3)
       dr2   = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)
       r2max = max(r2max,dr2)
    enddo
    tree%nodes(inode)%xcen = xcen
    tree%nodes(inode)%size = sqrt(r2max) + epsilon(r2max)

    !-- decide whether to split
    if (npnode > nmaxleaf) then ! inode remains a leaf
       !-- choose split dimension: longest dimension
       iaxis  = maxloc(xmax - xmin,1)
       xpivot = xcen(iaxis)

       !-- sort the indices in this inode by chosen dimension
       call sort_tree_ptmass_id(xyzmh_ptmass,tree%iptmassnode,istart,iend,iaxis,xpivot,imed)

       !-- create two child nodes (left: istart..m, right: m+1..iend)

       lchild  = nnodes + 1
       rchild  = nnodes + 2

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

       nnodes  = nnodes + 2

       !-- push next nodes in the stack
       if ((istack + 2) < istacksize ) then
          istack = istack + 1
          stack(istack) = lchild
          istack = istack + 1
          stack(istack) = rchild
       else
          call fatal("ptmass_tree","stack overflow in tree build")
       endif
    else !-- make sure the leaves are not connected to previous tree builds
       tree%nodes(inode)%lchild = 0
       tree%nodes(inode)%rchild = 0
    endif

 enddo

 tree%nnodes = nnodes

 if (allocated(stack)) deallocate(stack)


end subroutine build_ptmass_tree

 !-----------------------------------------------------------------
 ! Sort a subarray of idx(start:endp) by coordinate x(dim,*)
 ! A simple quicksort (recursive) for clarity.
 !-----------------------------------------------------------------
subroutine sort_tree_ptmass_id(xyzmh_ptmass,iptmassnode,il,ir,iaxis,xpivot,imed)
 use io, only: fatal
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
 ! Radius neighbour search (non-recursive traversal using stack)
 ! x(3,N) : coordinates
 ! tree   : kd-tree
 ! x0(3)  : query position
 ! rsearch: search radius
 ! neigh  : integer array (allocated by caller) to receive point indices
 ! nneigh : number of neighbours found (returned)
 !-----------------------------------------------------------------
subroutine get_ptmass_neigh(tree,xpos,rsearch,listneigh,nneigh)
 use io, only:fatal
 type(ptmasstree), intent(in)  :: tree
 real,             intent(in)  :: xpos(3)
 real,             intent(in)  :: rsearch
 integer,          intent(out) :: listneigh(:)
 integer,          intent(out) :: nneigh

 integer, allocatable :: stack(:)
 integer              :: istack
 integer              :: inode,i,istart,iend,il,ir
 real                 :: r2,dr2,dx,dy,dz,size
 logical              :: isleaf,isneigh

 r2 = rsearch*rsearch
 nneigh = 0

 allocate(stack(istacksize))
 istack = 1
 stack(istack) = iroot

 do while (istack > 0)

    inode = stack(istack)
    istack = istack - 1

    !-- compute squared distance from xpos to inode's bounding cube

    dx   = xpos(1) - tree%nodes(inode)%xcen(1)
    dy   = xpos(2) - tree%nodes(inode)%xcen(2)
    dz   = xpos(3) - tree%nodes(inode)%xcen(3)
    dr2  = dx*dx+dy*dy+dz*dz
    size = tree%nodes(inode)%size

    il = tree%nodes(inode)%lchild
    ir = tree%nodes(inode)%rchild

    isneigh = dr2 < (r2 + size*size)
    isleaf  = (il==0) .and. (ir==0)

    if (isneigh) then
       if (isleaf) then
          istart = tree%nodes(inode)%istart
          iend   = tree%nodes(inode)%iend
          do i = istart, iend
             nneigh = nneigh + 1
             if (nneigh > 1992) print*,inode
             listneigh(nneigh) = tree%iptmassnode(i)
          enddo
       else
          !-- push children (if exist)
          if ((istack + 2) < istacksize) then
             if (il /= 0) then
                istack = istack + 1
                stack(istack) = il
             endif
             if (ir /= 0) then
                istack = istack + 1
                stack(istack) = ir
             endif
          else
             call fatal("ptmass_neigh","stack overflow in tree build",ival=istack)
          endif
       endif
    endif
 enddo

 if (allocated(stack)) deallocate(stack)

end subroutine get_ptmass_neigh

end module ptmass_tree

