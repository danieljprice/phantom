!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testkdtree
!
! This module performs unit tests of the kdtree module
!   The tests here are specific to the tree, some general
!   tests of neighbour finding are done in test_neigh
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, kdtree, kernel, mpidomain, neighkdtree, part,
!   testutils, timing, unifdis
!
 implicit none
 public :: test_kdtree

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of tree code
!+
!-----------------------------------------------------------------------
subroutine test_kdtree(ntests,npass)
 use dim,         only:maxp,periodic,ind_timesteps
 use io,          only:id,master,iverbose
 use neighkdtree, only:leaf_is_active,ncells,node
 use part,        only:npart,xyzh,hfact,massoftype,igas,maxphase,iphase,isetphase,iactive
 use kernel,      only:hfact_default
 use kdtree,      only:maketree,revtree,kdnode,empty_tree
 use unifdis,     only:set_unifdis
 use testutils,   only:checkvalbuf,checkvalbuf_end,update_test_scores,checkval
 use timing,      only:print_time,getused
 use mpidomain,   only:i_belong
 integer, intent(inout) :: ntests,npass
 logical :: test_revtree, test_all
 integer :: i,nfailed(12),nchecked(12),nfailed_leaf(1),nchecked_leaf(1),ierrmax_leaf(1)
 real    :: psep,tol,errmax(12)
 real(4) :: t2,t1,tmaketree
 type(kdnode), allocatable :: old_tree(:)
 integer, allocatable :: leaf_is_active_saved(:)

 test_all = .true.
 test_revtree = .true.
 iverbose = 2

 if (id==master) write(*,"(a,/)") '--> TESTING KDTREE'

 if (test_revtree .or. test_all) then
    if (id==master) write(*,"(/,a)") '--> testing revtree routine'
    !
    ! set up a random particle distribution
    !
    psep = 1./100.
    hfact = hfact_default
    npart = 0
    call set_unifdis('random',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,&
                     psep,hfact,npart,xyzh,periodic,mask=i_belong)
    massoftype(igas) = 1000./npart
    if (maxphase==maxp) iphase(:) = isetphase(igas,iactive=.true.)

    !
    ! call maketree to build the tree
    !
    call empty_tree(node)
    call cpu_time(t1)
    call maketree(node,xyzh,npart,leaf_is_active,ncells,apr_tree=.false.)
    call cpu_time(t2)
    call print_time(t2-t1,'maketree completed in')
    !
    ! now save the tree structure and leaf_is_active
    !
    allocate(old_tree(int(ncells)))
    old_tree(1:ncells) = node(1:ncells)
    allocate(leaf_is_active_saved(int(ncells)))
    leaf_is_active_saved(1:int(ncells)) = leaf_is_active(1:int(ncells))

    !
    ! erase all information in the existing tree except the structure
    !
    do i=1,int(ncells)
       node(i)%xcen(:) = 0.
       node(i)%size    = 0.
       node(i)%hmax    = 0.
#ifdef GRAVITY
       node(i)%mass    = 0.
       node(i)%quads(:)= 0.
#endif
       leaf_is_active(i) = 0
    enddo

    !
    ! call revtree to rebuild
    !
    tmaketree = t2-t1
    call cpu_time(t1)
    call revtree(node,xyzh,leaf_is_active,ncells)
    call cpu_time(t2)
    call print_time(t2-t1,'revtree completed in')
    if (id==master) print*,' ratio of revtree/maketree: ',(t2-t1)/tmaketree

    !
    ! check that the revised tree matches the tree built
    !
    nfailed(:)  = 0
    nchecked(:) = 0
    errmax(:)   = 0.
    tol = 2.e-11
    do i=1,int(ncells)
       if (i > 1 .and. node(i)%parent == 0) cycle
       ! if (leaf_is_active(i) /= 0) then
       call checkvalbuf(node(i)%xcen(1),old_tree(i)%xcen(1),tol,'x0',nfailed(1),nchecked(1),errmax(1))
       call checkvalbuf(node(i)%xcen(2),old_tree(i)%xcen(2),tol,'y0',nfailed(2),nchecked(2),errmax(2))
       call checkvalbuf(node(i)%xcen(3),old_tree(i)%xcen(3),tol,'z0',nfailed(3),nchecked(3),errmax(3))
!       call checkvalbuf(node(i)%size,old_tree(i)%size,tol,'size',nfailed(4),nchecked(4),errmax(4))
       call checkvalbuf((node(i)%size + tol >= old_tree(i)%size),.true.,'size',nfailed(4),nchecked(4))
       call checkvalbuf(node(i)%hmax,old_tree(i)%hmax,tol,'hmax',nfailed(5),nchecked(5),errmax(5))
#ifdef GRAVITY
       call checkvalbuf(node(i)%mass,old_tree(i)%mass,tol,'mass',nfailed(6),nchecked(6),errmax(6))
       call checkvalbuf(node(i)%quads(1),old_tree(i)%quads(1),tol,'qxx',nfailed(7),nchecked(7),errmax(7))
       call checkvalbuf(node(i)%quads(2),old_tree(i)%quads(2),2.*tol,'qxy',nfailed(8),nchecked(8),errmax(8))
       call checkvalbuf(node(i)%quads(3),old_tree(i)%quads(3),tol,'qxz',nfailed(9),nchecked(9),errmax(9))
       call checkvalbuf(node(i)%quads(4),old_tree(i)%quads(4),tol,'qyy',nfailed(10),nchecked(10),errmax(10))
       call checkvalbuf(node(i)%quads(5),old_tree(i)%quads(5),tol,'qyz',nfailed(11),nchecked(11),errmax(11))
       call checkvalbuf(node(i)%quads(6),old_tree(i)%quads(6),tol,'qzz',nfailed(12),nchecked(12),errmax(12))
#endif
       ! endif
    enddo
    call checkvalbuf_end('x0',nchecked(1),nfailed(1),errmax(1),tol)
    call checkvalbuf_end('y0',nchecked(2),nfailed(2),errmax(2),tol)
    call checkvalbuf_end('z0',nchecked(3),nfailed(3),errmax(3),tol)
    call checkvalbuf_end('size',nchecked(4),nfailed(4),errmax(4),tol)
    call checkvalbuf_end('hmax',nchecked(5),nfailed(5),errmax(5),tol)
#ifdef GRAVITY
    call checkvalbuf_end('mass',nchecked(6),nfailed(6),errmax(6),tol)
    call checkvalbuf_end('qxx',nchecked(7),nfailed(7),errmax(7),tol)
    call checkvalbuf_end('qxy',nchecked(8),nfailed(8),errmax(8),tol)
    call checkvalbuf_end('qxz',nchecked(9),nfailed(9),errmax(9),tol)
    call checkvalbuf_end('qyy',nchecked(10),nfailed(10),errmax(10),tol)
    call checkvalbuf_end('qyz',nchecked(11),nfailed(11),errmax(11),tol)
    call checkvalbuf_end('qzz',nchecked(12),nfailed(12),errmax(12),tol)
#endif
    call update_test_scores(ntests,nfailed,npass)

    !
    ! check that leaf_is_active matches what maketree set
    !
    nfailed_leaf(:) = 0
    nchecked_leaf(:) = 0
    ierrmax_leaf(:) = 0
    do i=1,int(ncells)
       ! only check leaf nodes (non-zero leaf_is_active)
       if (leaf_is_active_saved(i) /= 0) then
          call checkvalbuf(leaf_is_active(i),leaf_is_active_saved(i),0,'leaf_is_active', &
                          nfailed_leaf(1),nchecked_leaf(1),ierrmax_leaf(1))
       endif
    enddo
    if (nchecked_leaf(1) > 0) then
       call checkvalbuf_end('leaf_is_active',nchecked_leaf(1),nfailed_leaf(1),ierrmax_leaf(1),0)
    endif
    call update_test_scores(ntests,nfailed_leaf,npass)

    deallocate(old_tree)
    deallocate(leaf_is_active_saved)
 endif

 if (id==master) write(*,"(/,a,/)") '<-- KDTREE TEST COMPLETE'

end subroutine test_kdtree
end module testkdtree
