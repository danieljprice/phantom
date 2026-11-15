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
! :Dependencies: boundary, dim, io, kdtree, kernel, mpidomain, mpikdtree,
!   mpisections, mpiutils, neighkdtree, part, testutils, timing, unifdis
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
 use dim,         only:maxp,periodic,ind_timesteps,mpi,use_apr,use_sinktree
 use io,          only:id,master,iverbose,nprocs
 use neighkdtree, only:leaf_is_active,ncells,node
 use part,        only:npart,xyzh,hfact,massoftype,igas,maxphase,iphase,isetphase,iactive,&
                       sync_aprmassoftype,init_part
 use kernel,      only:hfact_default
 use kdtree,      only:maketree,revtree,kdnode,empty_tree,inoderange,inodeparts
 use unifdis,     only:set_unifdis
 use testutils,   only:checkvalbuf,checkvalbuf_end,update_test_scores,checkval
 use timing,      only:print_time,getused
 use mpidomain,   only:i_belong
 integer, intent(inout) :: ntests,npass
 logical :: test_revtree, test_all, revised, can_revise
 integer :: i,ip,j,nfailed(22),nchecked(22),nfailed_leaf(1),nchecked_leaf(1),ierrmax_leaf(1)
 real    :: psep,tol,tol_octs,errmax(22),r2max,hmax,dx(3),rand(3)
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
    call init_part()
    psep = 1./100.
    hfact = hfact_default
    npart = 0
    call set_unifdis('random',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,&
                     psep,hfact,npart,xyzh,periodic,mask=i_belong)
    massoftype(igas) = 1000./npart
    call sync_aprmassoftype()
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
       node(i)%octs(:) = 0.
#endif
       leaf_is_active(i) = 0
    enddo

    !
    ! call revtree to rebuild
    !
    tmaketree = t2-t1
    call cpu_time(t1)
    call revtree(node,xyzh,npart,leaf_is_active,ncells,revised)
    call cpu_time(t2)
    call print_time(t2-t1,'revtree completed in')
    if (id==master) print*,' ratio of revtree/maketree: ',(t2-t1)/tmaketree
    ! as in revtree
    can_revise = .not.(periodic .or. use_apr .or. use_sinktree)
    call checkval(revised,can_revise,nfailed(1),'revtree revised the tree')
    call update_test_scores(ntests,nfailed(1:1),npass)
    if (.not.revised) then
       ! the tree was not revised, so restore it for the checks below
       node(1:ncells) = old_tree(1:ncells)
       leaf_is_active(1:int(ncells)) = leaf_is_active_saved(1:int(ncells))
    endif

    !
    ! check that the revised tree matches the tree built
    !
    nfailed(:)  = 0
    nchecked(:) = 0
    errmax(:)   = 0.
    tol = 2.e-11
    ! larger tolerance for the octupole moments, sums of cubes that nearly cancel, which the
    ! order of summation (OpenMP loop order, the particles on each MPI task) changes at round-off
    tol_octs = 3.e-10
    do i=1,int(ncells)
       if (i > 1 .and. node(i)%parent == 0) cycle
       ! if (leaf_is_active(i) /= 0) then
       call checkvalbuf(node(i)%xcen(1),old_tree(i)%xcen(1),tol,'x0',nfailed(1),nchecked(1),errmax(1))
       call checkvalbuf(node(i)%xcen(2),old_tree(i)%xcen(2),tol,'y0',nfailed(2),nchecked(2),errmax(2))
       call checkvalbuf(node(i)%xcen(3),old_tree(i)%xcen(3),tol,'z0',nfailed(3),nchecked(3),errmax(3))
       call checkvalbuf(node(i)%size,old_tree(i)%size,tol,'size',nfailed(4),nchecked(4),errmax(4))
       call checkvalbuf(node(i)%hmax,old_tree(i)%hmax,tol,'hmax',nfailed(5),nchecked(5),errmax(5))
#ifdef GRAVITY
       call checkvalbuf(node(i)%mass,old_tree(i)%mass,tol,'mass',nfailed(6),nchecked(6),errmax(6))
       call checkvalbuf(node(i)%quads(1),old_tree(i)%quads(1),tol,'qxx',nfailed(7),nchecked(7),errmax(7))
       call checkvalbuf(node(i)%quads(2),old_tree(i)%quads(2),2.*tol,'qxy',nfailed(8),nchecked(8),errmax(8))
       call checkvalbuf(node(i)%quads(3),old_tree(i)%quads(3),tol,'qxz',nfailed(9),nchecked(9),errmax(9))
       call checkvalbuf(node(i)%quads(4),old_tree(i)%quads(4),tol,'qyy',nfailed(10),nchecked(10),errmax(10))
       call checkvalbuf(node(i)%quads(5),old_tree(i)%quads(5),tol,'qyz',nfailed(11),nchecked(11),errmax(11))
       call checkvalbuf(node(i)%quads(6),old_tree(i)%quads(6),tol,'qzz',nfailed(12),nchecked(12),errmax(12))
       call checkvalbuf(node(i)%octs(1),old_tree(i)%octs(1),tol_octs,'oxxx',nfailed(13),nchecked(13),errmax(13))
       call checkvalbuf(node(i)%octs(2),old_tree(i)%octs(2),tol_octs,'oxxy',nfailed(14),nchecked(14),errmax(14))
       call checkvalbuf(node(i)%octs(3),old_tree(i)%octs(3),tol_octs,'oxxz',nfailed(15),nchecked(15),errmax(15))
       call checkvalbuf(node(i)%octs(4),old_tree(i)%octs(4),tol_octs,'oxyy',nfailed(16),nchecked(16),errmax(16))
       call checkvalbuf(node(i)%octs(5),old_tree(i)%octs(5),tol_octs,'oxyz',nfailed(17),nchecked(17),errmax(17))
       call checkvalbuf(node(i)%octs(6),old_tree(i)%octs(6),tol_octs,'oxzz',nfailed(18),nchecked(18),errmax(18))
       call checkvalbuf(node(i)%octs(7),old_tree(i)%octs(7),tol_octs,'oyyy',nfailed(19),nchecked(19),errmax(19))
       call checkvalbuf(node(i)%octs(8),old_tree(i)%octs(8),tol_octs,'oyyz',nfailed(20),nchecked(20),errmax(20))
       call checkvalbuf(node(i)%octs(9),old_tree(i)%octs(9),tol_octs,'oyzz',nfailed(21),nchecked(21),errmax(21))
       call checkvalbuf(node(i)%octs(10),old_tree(i)%octs(10),tol_octs,'ozzz',nfailed(22),nchecked(22),errmax(22))
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
    call checkvalbuf_end('oxxx',nchecked(13),nfailed(13),errmax(13),tol_octs)
    call checkvalbuf_end('oxxy',nchecked(14),nfailed(14),errmax(14),tol_octs)
    call checkvalbuf_end('oxxz',nchecked(15),nfailed(15),errmax(15),tol_octs)
    call checkvalbuf_end('oxyy',nchecked(16),nfailed(16),errmax(16),tol_octs)
    call checkvalbuf_end('oxyz',nchecked(17),nfailed(17),errmax(17),tol_octs)
    call checkvalbuf_end('oxzz',nchecked(18),nfailed(18),errmax(18),tol_octs)
    call checkvalbuf_end('oyyy',nchecked(19),nfailed(19),errmax(19),tol_octs)
    call checkvalbuf_end('oyyz',nchecked(20),nfailed(20),errmax(20),tol_octs)
    call checkvalbuf_end('oyzz',nchecked(21),nfailed(21),errmax(21),tol_octs)
    call checkvalbuf_end('ozzz',nchecked(22),nfailed(22),errmax(22),tol_octs)
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

    if (can_revise) then
       !
       ! move the particles, revise the tree, and check that every node
       ! bounds its particles and their smoothing lengths, with the size
       ! of the furthest particle from its centre, as when built
       !
       if (id==master) write(*,"(/,a)") '--> checking revtree for moved particles'
       do i=1,npart
          call random_number(rand)
          xyzh(1:3,i) = xyzh(1:3,i) + 0.5*psep*(rand - 0.5)
       enddo
       call revtree(node,xyzh,npart,leaf_is_active,ncells,revised)
       call checkval(revised,.true.,nfailed(1),'revtree revised the tree')
       call update_test_scores(ntests,nfailed(1:1),npass)
       nfailed(:)  = 0
       nchecked(:) = 0
       errmax(:)   = 0.
       do i=1,int(ncells)
          if (i > 1 .and. node(i)%parent == 0) cycle
          if (inoderange(1,i) <= 0) cycle
          r2max = 0.
          hmax  = 0.
          do ip=inoderange(1,i),inoderange(2,i)
             j  = abs(inodeparts(ip))
             dx = xyzh(1:3,j) - node(i)%xcen
             r2max = max(r2max,dot_product(dx,dx))
             hmax  = max(hmax,xyzh(4,j))
          enddo
          call checkvalbuf((sqrt(r2max) <= node(i)%size*(1. + epsilon(r2max))),.true.,&
                           'size bounds particles',nfailed(1),nchecked(1))
          call checkvalbuf((hmax <= node(i)%hmax),.true.,'hmax bounds particles',nfailed(2),nchecked(2))
          call checkvalbuf(node(i)%size,sqrt(r2max) + epsilon(r2max),tol,'size exact',nfailed(3),nchecked(3),errmax(3))
       enddo
       call checkvalbuf_end('size bounds particles',nchecked(1),nfailed(1))
       call checkvalbuf_end('hmax bounds particles',nchecked(2),nfailed(2))
       call checkvalbuf_end('size exact',nchecked(3),nfailed(3),errmax(3),tol)
       call update_test_scores(ntests,nfailed(1:3),npass)
    endif

    deallocate(old_tree)
    deallocate(leaf_is_active_saved)
 endif

 call test_section_tree(ntests,npass)
 call test_cut_sections(ntests,npass)
 if (mpi .and. nprocs > 1) call test_global_tree(ntests,npass)

 if (id==master) write(*,"(/,a,/)") '<-- KDTREE TEST COMPLETE'

end subroutine test_kdtree

!-----------------------------------------------------------------------
!+
!   Tests of the tree of sections used for the MPI domain decomposition,
!   on the cell centres of a grid in a flat box with two outliers far
!   above and below it: every particle lies in the section it is located
!   in (counts agree), the keys of the sections increase in depth first
!   order, no section cuts through the plane of the grid, and the first
!   split is at the midpoint of the bounds of the particles
!+
!-----------------------------------------------------------------------
subroutine test_section_tree(ntests,npass)
 use mpisections, only:sectiontree,build_sections,locate_section
 use mpiutils,    only:reduceall_mpi
 use testutils,   only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ngrid(3) = [16,16,4], np = ngrid(1)*ngrid(2)*ngrid(3) + 2
 real,    parameter :: boxlen(3) = [4.,4.,0.25]
 type(sectiontree) :: st
 real       :: xyzh(4,np),wt(np),dx(3),xmid
 integer    :: nfailed(4),nbad(4),i,ix,iy,iz,k,nround
 integer, allocatable :: layers(:)
 real,    allocatable :: nsec(:)

 ! cell centres of the grid, plus two particles far above and below it
 dx = boxlen/ngrid
 i = 0
 do iz = 1,ngrid(3)
    do iy = 1,ngrid(2)
       do ix = 1,ngrid(1)
          i = i + 1
          xyzh(1:3,i) = ([ix,iy,iz] - 0.5)*dx
       enddo
    enddo
 enddo
 xyzh(1:3,np-1) = [2.,2., 10.]
 xyzh(1:3,np)   = [2.,2.,-10.]
 xyzh(4,:) = 1.
 wt = 1.
 call build_sections(np,xyzh,wt,st,nround)

 nbad = 0

 ! the particles located in each section, on all tasks, are those it was built from
 allocate(nsec(st%nsec))
 nsec = 0.
 do i = 1,np
    k = locate_section(st,xyzh(1:3,i))
    nsec(k) = nsec(k) + 1.
 enddo
 nsec = reduceall_mpi('+',nsec)
 if (any(abs(nsec - st%secn) > 0.5)) nbad(1) = nbad(1) + 1

 ! keys in increasing order
 do k = 2,st%nsec
    if (st%key(st%seccell(k)) <= st%key(st%seccell(k-1))) nbad(2) = nbad(2) + 1
 enddo

 ! no section cuts through the plane: each holding grid particles has all its layers
 allocate(layers(st%nsec))
 layers = 0
 do i = 1,np-2
    k = locate_section(st,xyzh(1:3,i))
    iz = int(xyzh(3,i)/dx(3)) + 1
    layers(k) = ibset(layers(k),iz-1)
 enddo

 ! a section holds all layers of the grid or none
 do k = 1,st%nsec
    layers(k) = int(reduceall_mpi('max',real(layers(k))))
    if (layers(k) /= 0 .and. layers(k) /= 2**ngrid(3) - 1) nbad(3) = nbad(3) + 1
 enddo

 ! the root is split at the midpoint of the bounds of all particles, along z
 xmid = 0.5*(xyzh(3,np) + xyzh(3,np-1))
 if (st%axis(1) /= 3 .or. abs(st%xsplit(1) - xmid) > epsilon(xmid)) nbad(4) = nbad(4) + 1

 nfailed = 0
 call checkval(nbad(1),0,0,nfailed(1),'particles located in their sections')
 call checkval(nbad(2),0,0,nfailed(2),'section keys in depth first order')
 call checkval(nbad(3),0,0,nfailed(3),'flat distribution not cut through its plane')
 call checkval(nbad(4),0,0,nfailed(4),'split at the midpoint of the particle bounds')

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_section_tree

!-----------------------------------------------------------------------
!+
!   Tests of the cut of the sections into runs, one per task: the
!   largest work of a run is the least possible, every task gets a
!   section, and the cap on particles per task is respected
!+
!-----------------------------------------------------------------------
subroutine test_cut_sections(ntests,npass)
 use mpisections, only:cut_sections
 use testutils,   only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 real    :: w(10),np(10),wrun(4),nrun(4)
 integer :: owner(10),nfailed(3),k

 nfailed = 0

 ! ten sections of equal work on three tasks: the best cut has 4 in the largest run
 w  = 1.
 np = 1.
 call cut_sections(10,w,np,3,huge(1.),owner)
 wrun = 0.
 do k = 1,10
    wrun(owner(k)+1) = wrun(owner(k)+1) + w(k)
 enddo
 call checkval(maxval(wrun(1:3)),4.,0.,nfailed(1),'least largest work of a run')

 ! four sections on four tasks: one each, whatever their work
 w(1:4) = [10.,1.,1.,1.]
 call cut_sections(4,w(1:4),np(1:4),4,huge(1.),owner(1:4))
 call checkval(count(owner(1:4) == [0,1,2,3]),4,0,nfailed(2),'one section for every task')

 ! equal work, but the first five sections hold most of the particles: with at
 ! most 12 particles on a task they cannot share a run
 w  = 1.
 np = [5.,5.,5.,5.,5.,1.,1.,1.,1.,1.]
 call cut_sections(10,w,np,4,12.,owner)
 nrun = 0.
 do k = 1,10
    nrun(owner(k)+1) = nrun(owner(k)+1) + np(k)
 enddo
 call checkval(count(nrun > 12.),0,0,nfailed(3),'particles per task within the cap')

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_cut_sections
!-----------------------------------------------------------------------
!+
!   Tests of the global tree that every MPI task holds a copy of:
!   its links and leaf owners, its hmax after h grows on one task,
!   and that a cell is exported to every task holding a neighbour,
!   also after the trees are revised for moved particles on a
!   substep with few active particles
!+
!-----------------------------------------------------------------------
subroutine test_global_tree(ntests,npass)
 use dim,         only:periodic,maxphase,maxp,use_apr,use_sinktree
 use io,          only:id,master
 use boundary,    only:set_boundary
 use mpiutils,    only:reduceall_mpi
 use neighkdtree, only:build_tree,active_leaves,nactive_leaves,set_hmaxcell
 use mpikdtree,   only:sync_hmax_mpi
 use kdtree,      only:inoderange,inodeparts
 use part,        only:npart,xyzh,vxyzu,hfact,massoftype,igas,iphase,isetphase,init_part
 use kernel,      only:hfact_default
 use unifdis,     only:set_unifdis
 use mpidomain,   only:i_belong
 use testutils,   only:update_test_scores,checkval
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(18),i,ip,icell,k,nactive
 real    :: psep,hmaxcell,rand(3)
 logical :: can_revise,revised

 if (id==master) write(*,"(/,a)") '--> testing global tree'
 call init_part()
 if (periodic) call set_boundary()
 psep = 1./25.
 hfact = hfact_default
 npart = 0
 call set_unifdis('random',id,master,-0.5,0.5,-0.5,0.5,-0.5,0.5,&
                  psep,hfact,npart,xyzh,periodic,mask=i_belong)
 massoftype(igas) = 1./reduceall_mpi('+',npart)
 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,iactive=.true.)

 ! one particle whose kernel reaches into several domains
 if (id==master) xyzh(4,1) = 0.2

 call build_tree(npart,npart,xyzh,vxyzu)
 nfailed = 0
 call check_global_links(nfailed(1:2))
 call check_exports(.false.,nfailed(3))

 ! grow h on one task, as the density calculation does, and pass it to the global tree
 if (id==master) then
    do k=1,nactive_leaves
       icell = active_leaves(k)
       hmaxcell = 0.
       do ip=inoderange(1,icell),inoderange(2,icell)
          i = abs(inodeparts(ip))
          xyzh(4,i) = 1.5*xyzh(4,i)
          hmaxcell = max(hmaxcell,xyzh(4,i))
       enddo
       call set_hmaxcell(icell,hmaxcell)
    enddo
 endif
 call sync_hmax_mpi()
 call check_global_hmax(nfailed(4:5))
 call check_exports(.true.,nfailed(6))
 call check_combined_tree(nfailed(7:10))

 ! a substep with one particle in fifty active, after the particles have moved: the
 ! local trees are revised rather than rebuilt (as revtree allows), and the global
 ! tree from them must still find every neighbour
 nactive = 0
 do i=1,npart
    call random_number(rand)
    xyzh(1:3,i) = xyzh(1:3,i) + 0.25*psep*(rand - 0.5)
    if (maxphase==maxp) iphase(i) = isetphase(igas,iactive=(mod(i,50)==0))
    if (mod(i,50)==0 .or. maxphase /= maxp) nactive = nactive + 1
 enddo
 call build_tree(npart,nactive,xyzh,vxyzu,tree_revised=revised)
 can_revise = maxphase==maxp .and. .not.(periodic .or. use_apr .or. use_sinktree)
 call checkval(revised,can_revise,nfailed(11),'global tree revised for moved particles')
 call check_global_links(nfailed(12:13))
 call check_exports(.true.,nfailed(14))
 call check_combined_tree(nfailed(15:18))
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_global_tree

!-----------------------------------------------------------------------
!+
!   list the nodes of the global tree reached from the root, stopping
!   at leaves and at nodes whose links are broken (counted in nbad)
!+
!-----------------------------------------------------------------------
subroutine walk_global_tree(list,nlist,nbad)
 use mpikdtree,   only:nodeglobal,cellatid
 integer, intent(out) :: list(:),nlist,nbad
 integer :: n,il,ir,i

 nbad = 0

 ! list the nodes reachable from the root, breadth first
 nlist = 1
 list(1) = 1
 i = 0
 do while (i < nlist)
    i = i + 1
    n = list(i)
    if (cellatid(n) /= 0) cycle
    il = nodeglobal(n)%leftchild
    ir = nodeglobal(n)%rightchild
    if (il <= 0 .or. ir <= 0 .or. nlist + 2 > size(list)) then
       nbad = nbad + 1
    elseif (nodeglobal(il)%parent /= n .or. nodeglobal(ir)%parent /= n) then
       nbad = nbad + 1
    else
       list(nlist+1) = il
       list(nlist+2) = ir
       nlist = nlist + 2
    endif
 enddo

end subroutine walk_global_tree

!-----------------------------------------------------------------------
!+
!   check the links of the global tree, that its leaves have no
!   children and one owner, and that the leaves of each task hold its mass
!+
!-----------------------------------------------------------------------
subroutine check_global_links(nfailed)
 use dim,         only:gravity
 use io,          only:id,nprocs
 use mpikdtree,   only:nodeglobal,cellatid
 use part,        only:npart,massoftype,igas
 use testutils,   only:checkval
 integer, intent(out) :: nfailed(2)
 integer, allocatable :: list(:)
 integer :: n,i,nlist,nbad,nleaves,iowner
 real    :: leafmass(nprocs)

 allocate(list(size(nodeglobal)))
 call walk_global_tree(list,nlist,nbad)

 ! each leaf: one owner, no children, and its mass added to its owner's
 nleaves = 0
 leafmass = 0.
 do i=1,nlist
    n = list(i)
    if (cellatid(n) == 0) cycle
    nleaves = nleaves + 1
    iowner = cellatid(n)
    if (iowner < 1 .or. iowner > nprocs .or. nodeglobal(n)%leftchild /= 0 .or. nodeglobal(n)%rightchild /= 0) then
       nbad = nbad + 1
    else
       leafmass(iowner) = leafmass(iowner) + nodeglobal(n)%mass
    endif
 enddo

 ! every node marked as a leaf must be reachable from the root
 if (count(cellatid /= 0) /= nleaves) nbad = nbad + 1
 call checkval(nbad,0,0,nfailed(1),'global tree links and leaf owners')

 ! node masses are only set with self-gravity
 nfailed(2) = 0
 if (gravity) call checkval(leafmass(id+1),npart*massoftype(igas),1.e-12,nfailed(2),'mass in leaves of this task')

end subroutine check_global_links

!-----------------------------------------------------------------------
!+
!   check that every node of the global tree bounds the hmax of its
!   children, and that the leaves of this task bound its particles' h
!+
!-----------------------------------------------------------------------
subroutine check_global_hmax(nfailed)
 use io,          only:id
 use mpikdtree,   only:nodeglobal,cellatid
 use part,        only:npart,xyzh
 use testutils,   only:checkval
 integer, intent(out) :: nfailed(2)
 integer, allocatable :: list(:)
 integer :: n,i,nlist,nbad
 real    :: hmaxleaves

 allocate(list(size(nodeglobal)))
 call walk_global_tree(list,nlist,nbad)

 ! each internal node's hmax covers its children's
 hmaxleaves = 0.
 do i=1,nlist
    n = list(i)
    if (cellatid(n) == id + 1) hmaxleaves = max(hmaxleaves,nodeglobal(n)%hmax)
    if (cellatid(n) /= 0) cycle
    if (nodeglobal(n)%hmax < max(nodeglobal(nodeglobal(n)%leftchild)%hmax,&
                                 nodeglobal(nodeglobal(n)%rightchild)%hmax)) nbad = nbad + 1
 enddo
 call checkval(nbad,0,0,nfailed(1),'global hmax bounds children')
 call checkval(hmaxleaves >= maxval(xyzh(4,1:npart)),.true.,nfailed(2),'global hmax bounds h on this task')

end subroutine check_global_hmax

!-----------------------------------------------------------------------
!+
!   check the combined tree (local tree with the global tree above it):
!   its links, that walking down from its root reaches every local leaf
!   and only leaves of other tasks outside this task, and that every task
!   has exactly the same copy of each section root as its owner, and that
!   every section has one owner, the same on all tasks, with the roots of
!   a task's sections holding exactly the particles of its local tree
!+
!-----------------------------------------------------------------------
subroutine check_combined_tree(nfailed)
 use io,          only:id,nprocs
 use dim,         only:gravity
 use mpiutils,    only:reduceall_mpi
 use neighkdtree, only:node,leaf_is_active,ncells
 use mpikdtree,   only:ncomb,irootcomb,remote_leaf_owner,nsec_top,top_leaf,top_owner,top_root
 use kdtree,      only:irootnode,inoderange
 use testutils,   only:checkval
 integer, intent(out) :: nfailed(4)
 integer, allocatable :: list(:),ownmin(:),ownmax(:)
 real,    allocatable :: secdata(:,:)
 integer :: n,i,k,nlist,nbad,nlocal,nlocal_root,nbadremote,ndiff
 integer :: nbadown,nown,npown,nproot
 real    :: nodedata(6)

 allocate(list(ncomb))

 ! walk from the root of the combined tree
 call walk_tree_from(irootcomb,list,nlist,nbad)
 nlocal = 0
 nbadremote = 0
 do i=1,nlist
    n = list(i)
    if (leaf_is_active(n) == 0) cycle
    if (n <= ncells) then
       nlocal = nlocal + 1
    elseif (remote_leaf_owner(n) < 1 .or. remote_leaf_owner(n) > nprocs .or. remote_leaf_owner(n) == id + 1) then
       nbadremote = nbadremote + 1
    endif
 enddo

 ! the local leaves reached from the local root
 call walk_tree_from(irootnode,list,nlist,ndiff,local_only=.true.)
 nlocal_root = count(leaf_is_active(list(1:nlist)) /= 0)
 call checkval(nbad+ndiff+nbadremote,0,0,nfailed(1),'combined tree links and leaf owners')
 call checkval(nlocal,nlocal_root,0,nfailed(2),'combined tree reaches every local leaf')

 ! each section's root, as its owner has it, against its leaf on the other tasks;
 ! node masses are only set with self-gravity
 allocate(secdata(6,nsec_top))
 secdata = 0.
 do k=1,nsec_top
    if (top_owner(k) /= id) cycle
    n = top_root(k)
    secdata(:,k) = [node(n)%xcen,node(n)%size,node(n)%hmax,0.]
    if (gravity) secdata(6,k) = node(n)%mass
 enddo
 secdata = reduceall_mpi('+',secdata)

 ! the copies of other tasks' section roots match the owners' data
 ndiff = 0
 do k=1,nsec_top
    if (top_owner(k) == id) cycle
    n = int(ncells) + top_leaf(k)
    nodedata = [node(n)%xcen,node(n)%size,node(n)%hmax,0.]
    if (gravity) nodedata(6) = node(n)%mass
    if (any(nodedata /= secdata(:,k))) ndiff = ndiff + 1
 enddo
 call checkval(ndiff,0,0,nfailed(3),'section roots identical on all tasks')

 ! every section has a valid owner, the same on all tasks, and a top leaf of
 ! its own; every task owns a section, and the roots of its sections hold
 ! exactly the particles of its local tree
 nbadown = 0
 ownmin = int(reduceall_mpi('min',top_owner(1:nsec_top)))
 ownmax = int(reduceall_mpi('max',top_owner(1:nsec_top)))
 nown  = 0
 npown = 0
 do k=1,nsec_top
    if (top_owner(k) < 0 .or. top_owner(k) >= nprocs) nbadown = nbadown + 1
    if (ownmin(k) /= ownmax(k)) nbadown = nbadown + 1
    if (count(top_leaf(1:nsec_top) == top_leaf(k)) /= 1) nbadown = nbadown + 1
    if (top_owner(k) == id) then
       nown = nown + 1
       n = top_root(k)
       npown = npown + inoderange(2,n) - inoderange(1,n) + 1
    endif
 enddo
 if (nown == 0) nbadown = nbadown + 1
 nproot = inoderange(2,irootnode) - inoderange(1,irootnode) + 1
 if (npown /= nproot) nbadown = nbadown + 1
 nbadown = int(reduceall_mpi('+',nbadown))
 call checkval(nbadown,0,0,nfailed(4),'every section owned by one task')

end subroutine check_combined_tree

!-----------------------------------------------------------------------
!+
!   list the nodes below (and including) iroot in node, stopping at
!   leaves, counting nodes whose children do not point back to them.
!   With local_only, children whose parent is in the global tree (the
!   roots of the sections, linked to it in the combined tree) are fine
!+
!-----------------------------------------------------------------------
subroutine walk_tree_from(iroot,list,nlist,nbad,local_only)
 use neighkdtree, only:node,leaf_is_active,ncells
 integer, intent(in)  :: iroot
 integer, intent(out) :: list(:),nlist,nbad
 logical, intent(in), optional :: local_only
 integer :: i,n,ichild(2),k,iparent
 logical :: linked

 nbad = 0

 ! list the nodes reachable from iroot, breadth first
 nlist = 1
 list(1) = iroot
 i = 0
 do while (i < nlist)
    i = i + 1
    n = list(i)
    if (leaf_is_active(n) /= 0) cycle

    ! each child links back to its parent
    ichild = [node(n)%leftchild,node(n)%rightchild]
    do k=1,2
       if (ichild(k) <= 0) cycle
       iparent = node(ichild(k))%parent
       linked = (iparent == n)
       if (present(local_only)) linked = linked .or. (local_only .and. iparent > ncells)
       if (.not.linked .or. nlist >= size(list)) then
          nbad = nbad + 1
       else
          nlist = nlist + 1
          list(nlist) = ichild(k)
       endif
    enddo
 enddo

end subroutine walk_tree_from

!-----------------------------------------------------------------------
!+
!   check that the tasks each local cell is exported to include every
!   task holding a particle within the kernel radius, by brute force
!   over all particles on all tasks. Also reports how many exports
!   are needed, as a measure of how fine the global tree is.
!+
!-----------------------------------------------------------------------
subroutine check_exports(getj,nfailed)
 use dim,         only:periodic
 use io,          only:id,master,nprocs
 use boundary,    only:mod_boundary
 use mpiutils,    only:reduceall_mpi,reduce_in_place_mpi
 use neighkdtree, only:get_neighbour_list,active_leaves,nactive_leaves,listneigh
 use kdtree,      only:inoderange,inodeparts
 use kernel,      only:radkern
 use part,        only:npart,xyzh
 use testutils,   only:checkval
 logical, intent(in)  :: getj
 integer, intent(out) :: nfailed
 real,    allocatable :: xall(:,:)
 integer :: nper(nprocs),ntot,ioffset,i,j,ip,k,icell,nneigh,iowner,nmissed
 integer(kind=8) :: nexport,nneeded
 logical :: export(nprocs),needed(nprocs)
 real    :: xyzcache(4,1),dx(3),hij

 ! gather the positions, h and owner of all particles on all tasks
 nper = 0
 nper(id+1) = npart
 nper = int(reduceall_mpi('+',nper))
 ntot = sum(nper)
 ioffset = sum(nper(1:id))
 allocate(xall(5,ntot))
 xall = 0.
 xall(1:4,ioffset+1:ioffset+npart) = xyzh(1:4,1:npart)
 xall(5,ioffset+1:ioffset+npart) = id
 call reduce_in_place_mpi('+',xall)

 ! for each active leaf, the tasks the walk exports to against the tasks with neighbours
 nmissed = 0
 nexport = 0
 nneeded = 0
 do k=1,nactive_leaves
    icell = active_leaves(k)
    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,0,getj=getj,remote_export=export)

    ! the tasks holding a neighbour of any particle of the cell
    needed = .false.
    do ip=inoderange(1,icell),inoderange(2,icell)
       i = abs(inodeparts(ip))
       do j=1,ntot
          iowner = nint(xall(5,j))
          if (iowner == id .or. needed(iowner+1)) cycle
          dx = xyzh(1:3,i) - xall(1:3,j)
          if (periodic) call mod_boundary(dx)
          hij = xyzh(4,i)
          if (getj) hij = max(hij,xall(4,j))
          if (dot_product(dx,dx) < (radkern*hij)**2) needed(iowner+1) = .true.
       enddo
    enddo
    if (any(needed .and. .not.export)) nmissed = nmissed + 1
    nexport = nexport + count(export)
    nneeded = nneeded + count(needed)
 enddo
 nexport = reduceall_mpi('+',nexport)
 nneeded = reduceall_mpi('+',nneeded)
 if (id==master) write(*,"(a,l1,2(a,i0),a,f6.3)") ' getj=',getj,': cells exported to ',nexport,&
    ' tasks, needed ',nneeded,', fraction needed ',real(nneeded)/max(real(nexport),1.)
 call checkval(nmissed,0,0,nfailed,'cells exported to every task with a neighbour')

end subroutine check_exports

end module testkdtree
