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
 integer, parameter   :: istacksize = 512

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
! Build ptmass KD-tree in parallel using omp task and atomic
!+
!-----------------------------------------------------------------
subroutine build_ptmass_tree(xyzmh_ptmass,nptmass,tree)
 use io, only:fatal
 real,             intent(in)    :: xyzmh_ptmass(:,:)
 integer,          intent(in)    :: nptmass
 type(ptmasstree), intent(inout) :: tree
 integer, allocatable :: stack(:)
 integer :: i,k,iaxis,inode,lchild,rchild
 integer :: istart,iend,imed,npnode
 integer :: nlvl,maxlevel_indexed,ilvl
 real    :: xmin(3), xmax(3)
 real    :: dx(3),dr2,r2max
 integer :: nnodes,istack,nnodes_old,myslot,isr,isl
 real    :: xtmp(3),xcen(3)
 real    :: xpivot,mtot
 logical :: buildingtree

 if (.not.allocated(stack)) allocate(stack(istacksize))

 maxlevel_indexed = int(log(real(nnodeptmassmax+1))/log(2.)) - 1

 npnode = 0
 buildingtree = .true.

 !-- initialize tree part index array (can't be parallel)
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

 istack = 0
 ilvl   = 0
 nlvl   = 1
 inode  = 0
 !$omp parallel default(none)&
 !$omp shared(tree,xyzmh_ptmass,nnodeptmassmax)&
 !$omp shared(maxlevel_indexed,ilvl,nnodes,istack,stack)&
 !$omp private(istart,iend,npnode,xmin,xmax,xcen,r2max,dr2,dx,nnodes_old)&
 !$omp private(iaxis,xpivot,imed,lchild,rchild,xtmp,mtot,isr,isl,myslot)&
 !$omp firstprivate(buildingtree,nlvl,inode)
 !$omp single
 do while(buildingtree)
    nnodes_old = nnodes
    do k=1,nlvl
       if (ilvl >= maxlevel_indexed) then
          inode = stack(istack)
          istack = istack - 1
       else
          inode = inode + 1
       endif

       istart = tree%nodes(inode)%istart
       iend   = tree%nodes(inode)%iend
       npnode = iend - istart + 1

       !$omp task firstprivate(inode,istart,iend,npnode)

       !-- compute bounding box and center for this inode
       xmin =  huge(0.)
       xmax = -huge(0.)
       xcen = 0.
       mtot = 0.

       do i=istart,iend
          xtmp = xyzmh_ptmass(1:3, tree%iptmassnode(i))
          xmin(1) = min(xmin(1), xtmp(1))
          xmin(2) = min(xmin(2), xtmp(2))
          xmin(3) = min(xmin(3), xtmp(3))
          xmax(1) = max(xmax(1), xtmp(1))
          xmax(2) = max(xmax(2), xtmp(2))
          xmax(3) = max(xmax(3), xtmp(3))
          xcen    = xcen + xtmp*xyzmh_ptmass(4, tree%iptmassnode(i))
          mtot    = mtot + xyzmh_ptmass(4, tree%iptmassnode(i))
       enddo

       xcen  = xcen/mtot
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

          if (ilvl < maxlevel_indexed) then
             lchild  = inode*2
             rchild  = lchild + 1
          else
             !$omp atomic capture
             myslot  = nnodes
             nnodes  = nnodes + 2
             !$omp end atomic
             lchild  = myslot + 1
             rchild  = myslot + 2
          endif

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

          if ((ilvl>=maxlevel_indexed)) then
             !$omp atomic capture
             myslot  = istack
             istack  = istack + 2
             !$omp end atomic
             isl     = myslot + 1
             isr     = myslot + 2
             if (isr > istacksize) call fatal("ptmass_tree","stack overflow...",ival=isr)
             stack(isl) = lchild
             stack(isr) = rchild
          endif

       else !-- make sure the leaves are not connected to previous tree builds
          tree%nodes(inode)%lchild  = 0
          tree%nodes(inode)%rchild  = 0
       endif
       !$omp end task
    enddo
    !$omp taskwait
    ilvl = ilvl + 1

    if (ilvl >= maxlevel_indexed) then
       nlvl = istack
    else
       nlvl = 2**ilvl
       nnodes = 2**(ilvl)-1
    endif

    buildingtree = nnodes /= nnodes_old
 enddo
 !$omp end single
 !$omp end parallel
 tree%nnodes = nnodes

 if (allocated(stack)) deallocate(stack)

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

 integer, allocatable :: stack(:)
 integer              :: istack
 integer              :: inode,i,istart,iend,il,ir,npnode
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

    isneigh = dr2 < (r2 + size*size)
    isleaf  = (il==0) .and. (ir==0)

    if (isneigh) then
       npnode = (istart - iend) + 1
       if (npnode > 0) then
          if (isleaf) then
             do i = istart, iend
                nneigh = nneigh + 1
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
    endif
 enddo

 if (allocated(stack)) deallocate(stack)

end subroutine get_ptmass_neigh

end module ptmass_tree

