!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module kdtree
!
! This module implements the k-d tree build
!    and associated tree walking routines
!
! :References:
!    Gafton & Rosswog (2011), MNRAS 418, 770-781
!    Benz, Bowers, Cameron & Press (1990), ApJ 348, 647-667
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: allocutils, boundary, dim, dtypekdtree, io, kernel,
!   mpibalance, mpidomain, mpitree, mpiutils, part, timing
!
 use dim,         only:maxp,ncellsmax,minpart,use_apr,use_sinktree,maxptmass,maxpsph
 use io,          only:nprocs
 use dtypekdtree, only:kdnode,lenfgrav
 use part,        only:ll,iphase,treecache,maxphase, &
                       apr_level,aprmassoftype

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)
 type(kdnode),     allocatable :: refinementnode(:)
 real,             allocatable :: fnode_branch(:,:)
!$omp threadprivate(fnode_branch)
!
!--tree parameters
!
 integer,          parameter, public :: irootnode  = 1
 character(len=1), parameter, public :: labelax(3) = (/'x','y','z'/)
 integer,          parameter         :: maxdepth   = 32
!
!--runtime options for this module
!
 real,    public  :: tree_accuracy = 0.5
 logical, private :: done_init_kdtree = .false.
 logical, private :: already_warned = .false.
 integer, private :: numthreads

! Index of the last node in the local tree that has been copied to
! the global tree
 integer :: irefine

 public :: allocate_kdtree, deallocate_kdtree
 public :: maketree, revtree, getneigh,getneigh_dual,kdnode,lenfgrav
 public :: maketreeglobal
 public :: empty_tree
 public :: compute_M2L,expand_fgrav_in_taylor_series
 integer, public :: maxlevel_indexed, maxlevel

 type kdbuildstack
    integer :: node
    integer :: parent
    integer :: level
    integer :: npnode
    real    :: xmin(3)
    real    :: xmax(3)
 end type kdbuildstack

 private

contains

subroutine allocate_kdtree
 use dim, only:mpi
 use allocutils, only:allocate_array

 call allocate_array('inoderange', inoderange, 2, ncellsmax+1)
 call allocate_array('inodeparts', inodeparts, maxp)
 if (mpi) call allocate_array('refinementnode', refinementnode, ncellsmax+1)
!$omp parallel
 call allocate_array('fnode_branch', fnode_branch, lenfgrav, maxdepth)
!$omp end parallel

end subroutine allocate_kdtree

subroutine deallocate_kdtree
 use dim, only:mpi
 if (allocated(inoderange)) deallocate(inoderange)
 if (allocated(inodeparts)) deallocate(inodeparts)
 if (mpi .and. allocated(refinementnode)) deallocate(refinementnode)
!$omp parallel
 if (allocated(fnode_branch)) deallocate(fnode_branch)
!$omp end parallel

end subroutine deallocate_kdtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the tree from scratch
!
!  Notes/To do:
!  -openMP parallelisation of maketree_stack (done - April 2013)
!  -test centre of mass vs. geometric centre based cell sizes (done - 2013)
!  -need analysis module that times (and checks) build_tree and tree walk
!   for a given dump
!  -test bottom-up vs. top-down neighbour search
!  -should we try to store tree structure with particle arrays?
!  -need to compute centre of mass and moments for each cell on the fly (c.f. revtree?)
!  -need to implement long-range gravitational interaction (done - May 2013)
!  -implement revtree routine to update tree w/out rebuilding (done - Sep 2015)
!+
!-------------------------------------------------------------------------------
subroutine maketree(node, xyzh, np, leaf_is_active, ncells, apr_tree, refinelevels,nptmass,xyzmh_ptmass)
 use io,   only:fatal,warning,iprint,iverbose
!$ use omp_lib
 type(kdnode),      intent(out)   :: node(:) !ncellsmax+1)
 integer,           intent(in)    :: np
 real,              intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,           intent(out)   :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8),   intent(out)   :: ncells
 logical,           intent(in)    :: apr_tree
 integer, optional, intent(out)   :: refinelevels
 integer, optional, intent(in)    :: nptmass
 real,    optional, intent(inout) :: xyzmh_ptmass(:,:)

 integer :: i,npnode,il,ir,istack,nl,nr,mymum
 integer :: nnode,minlevel,level,nqueue
 real :: xmini(3),xmaxi(3),xminl(3),xmaxl(3),xminr(3),xmaxr(3)
 integer, parameter :: istacksize = 512
 type(kdbuildstack), save :: stack(istacksize)
 !$omp threadprivate(stack)
 type(kdbuildstack) :: queue(istacksize)
!$ integer :: threadid
 integer :: npcounter
 logical :: wassplit,finished,sinktree
 character(len=10) :: string

 if (present(nptmass) .and. present(xyzmh_ptmass)) then
    sinktree = .true.
 endif

 leaf_is_active = 0

 ir = 0
 il = 0
 nl = 0
 nr = 0
 wassplit = .false.
 finished = .false.

 ! construct root node, i.e. find bounds of all particles
 if (sinktree) then
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,xyzmh_ptmass,nptmass)
 else
    call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh)
 endif

 if (inoderange(1,irootnode)==0 .or. inoderange(2,irootnode)==0 ) then
    call fatal('maketree','no particles or all particles dead/accreted')
 endif

! Put root node on top of stack
 ncells = 1
 maxlevel = 0
 minlevel = maxdepth - 1
 istack = 1

 ! maximum level where 2^k indexing can be used (thus avoiding critical sections)
 ! deeper than this we access cells via a stack as usual
 maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 1

 ! default number of cells is the size of the `indexed' part of the tree
 ! this can be *increased* by building tree beyond indexed levels
 ! and is decreased afterwards according to the maximum depth actually reached
 ncells = 2**(maxlevel_indexed+1) - 1

 ! need to number of particles in node during build
 ! this is counted above to remove dead/accreted particles
 call push_onto_stack(queue(istack),irootnode,0,0,npcounter,xmini,xmaxi)

 if (.not.done_init_kdtree) then
    ! 1 thread for serial, overwritten when using OpenMP
    numthreads = 1

    ! get number of OpenMP threads
    !$omp parallel default(none) shared(numthreads)
!$  numthreads = omp_get_num_threads()
    !$omp end parallel
    done_init_kdtree = .true.
 endif

 nqueue = numthreads
 ! build using a queue to build level by level until number of nodes = number of threads
 over_queue: do while (istack  <  nqueue)
    ! if the tree finished while building the queue, then we should just return
    ! only happens for small particle numbers
    if (istack <= 0) then
       finished = .true.
       exit over_queue
    endif
    ! pop off front of queue
    call pop_off_stack(queue(1), istack, nnode, mymum, level, npnode, xmini, xmaxi)

    ! shuffle queue forward
    do i=1,istack
       queue(i) = queue(i+1)
    enddo

    ! construct node
    if (sinktree) then
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, .false.,apr_tree,xyzmh_ptmass)
    else
       call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .true., &  ! construct in parallel
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit, .false.,apr_tree)
    endif

    if (wassplit) then ! add children to back of queue
       if (istack+2 > istacksize) call fatal('maketree',&
                                       'queue size exceeded in tree build, increase istacksize and recompile')

       istack = istack + 1
       call push_onto_stack(queue(istack),il,nnode,level+1,nl,xminl,xmaxl)
       istack = istack + 1
       call push_onto_stack(queue(istack),ir,nnode,level+1,nr,xminr,xmaxr)
    endif

 enddo over_queue

 ! fix the indices

 done: if (.not.finished) then

    ! build using a stack which builds depth first
    ! each thread grabs a node from the queue and builds its own subtree

    !$omp parallel default(none) &
    !$omp shared(queue) &
    !$omp shared(ll, leaf_is_active) &
    !$omp shared(xyzmh_ptmass) &
    !$omp shared(np) &
    !$omp shared(node, ncells) &
    !$omp shared(nqueue,apr_tree,sinktree) &
    !$omp private(istack) &
    !$omp private(nnode, mymum, level, npnode, xmini, xmaxi) &
    !$omp private(ir, il, nl, nr) &
    !$omp private(xminr, xmaxr, xminl, xmaxl) &
    !$omp private(threadid) &
    !$omp private(wassplit) &
    !$omp reduction(min:minlevel) &
    !$omp reduction(max:maxlevel)
    !$omp do schedule(static)
    do i = 1, nqueue

       stack(1) = queue(i)
       istack = 1

       over_stack: do while(istack > 0)

          ! pop node off top of stack
          call pop_off_stack(stack(istack), istack, nnode, mymum, level, npnode, xmini, xmaxi)

          ! construct node
          if (sinktree) then
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, .false.,apr_tree,xyzmh_ptmass)
          else
             call construct_node(node(nnode), nnode, mymum, level, xmini, xmaxi, npnode, .false., &  ! don't construct in parallel
                                 il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr, ncells, leaf_is_active, &
                                 minlevel, maxlevel, wassplit, .false.,apr_tree)
          endif

          if (wassplit) then ! add children to top of stack
             if (istack+2 > istacksize) call fatal('maketree',&
                                       'stack size exceeded in tree build, increase istacksize and recompile')

             istack = istack + 1
             call push_onto_stack(stack(istack),il,nnode,level+1,nl,xminl,xmaxl)
             istack = istack + 1
             call push_onto_stack(stack(istack),ir,nnode,level+1,nr,xminr,xmaxr)
          endif

       enddo over_stack
    enddo
    !$omp enddo
    !$omp end parallel

 endif done

 ! decrease number of cells if tree is entirely within 2^k indexing limit
 if (maxlevel < maxlevel_indexed) then
    ncells = 2**(maxlevel+1) - 1
 endif
 if (maxlevel > maxlevel_indexed .and. .not.already_warned) then
    write(string,"(i10)") 2**(maxlevel-maxlevel_indexed)
    if (iverbose > 0) call warning('maketree','maxlevel > max_indexed: will run faster if recompiled with '// &
               'NCELLSMAX='//trim(adjustl(string))//'*maxp,')
 endif

 if (present(refinelevels)) refinelevels = minlevel

 if (iverbose >= 3) then
    write(iprint,"(a,i10,3(a,i2))") ' maketree: nodes = ',ncells,', max level = ',maxlevel,&
       ', min leaf level = ',minlevel,' max level indexed = ',maxlevel_indexed
 endif

end subroutine maketree

!----------------------------
!+
! routine to empty the tree
!+
!----------------------------
subroutine empty_tree(node)
 type(kdnode), intent(out) :: node(:)
 integer :: i

!$omp parallel do private(i)
 do i=1,size(node)
    node(i)%xcen = 0.
    node(i)%size = 0.
    node(i)%hmax = 0.
    node(i)%leftchild = 0
    node(i)%rightchild = 0
    node(i)%parent = 0
#ifdef GRAVITY
    node(i)%mass  = 0.
    node(i)%quads = 0.
#endif
 enddo
!$omp end parallel do

end subroutine empty_tree

!---------------------------------
!+
! routine to construct root node
!+
!---------------------------------
subroutine construct_root_node(np,nproot,irootnode,xmini,xmaxi,leaf_is_active,xyzh,xyzmh_ptmass,nptmass)
 use boundary, only:cross_boundary
 use mpidomain,only:isperiodic
 use part, only:iphase,iactive
 use part, only:isdead_or_accreted,ibelong
 use io,   only:fatal,id
 use dim,  only:ind_timesteps,mpi,periodic
 use part, only:isink,massoftype,igas,iamtype,maxphase,maxp,aprmassoftype,apr_level,ihsoft
 integer,          intent(in)    :: np,irootnode
 integer,          intent(out)   :: nproot
 real,             intent(out)   :: xmini(3), xmaxi(3)
 integer,          intent(inout) :: leaf_is_active(:)
 real,             intent(inout) :: xyzh(:,:)
 real,   optional, intent(inout) :: xyzmh_ptmass(:,:)
 integer, optional, intent(in)    :: nptmass
 integer :: i,ncross
 real    :: xminpart,yminpart,zminpart,xmaxpart,ymaxpart,zmaxpart
 real    :: xi, yi, zi

 xminpart = xyzh(1,1)
 yminpart = xyzh(2,1)
 zminpart = xyzh(3,1)
 xmaxpart = xminpart
 ymaxpart = yminpart
 zmaxpart = zminpart

 ncross = 0
 nproot = 0
 !$omp parallel default(none) &
 !$omp shared(np,xyzh,nptmass,xyzmh_ptmass) &
 !$omp shared(inodeparts,iphase,treecache,nproot) &
 !$omp shared(id,use_sinktree) &
 !$omp shared(isperiodic) &
 !$omp private(i,xi,yi,zi) &
 !$omp reduction(min:xminpart,yminpart,zminpart) &
 !$omp reduction(max:xmaxpart,ymaxpart,zmaxpart) &
 !$omp reduction(+:ncross)
 !$omp do schedule(guided,1)
 do i=1,np
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       if (periodic) call cross_boundary(isperiodic,xyzh(:,i),ncross)
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       if (isnan(xi) .or. isnan(yi) .or. isnan(zi)) then
          call fatal('maketree','NaN in particle position, likely caused by NaN in force',i,var='x',val=xi)
       endif
       xminpart = min(xminpart,xi)
       yminpart = min(yminpart,yi)
       zminpart = min(zminpart,zi)
       xmaxpart = max(xmaxpart,xi)
       ymaxpart = max(ymaxpart,yi)
       zmaxpart = max(zmaxpart,zi)
    endif
 enddo
 !$omp enddo
 !$omp barrier
 if (use_sinktree) then
    if (nptmass>0) then
       !$omp do schedule(guided,1)
       do i=1,nptmass
          if (xyzmh_ptmass(4,i)>0.) then
             if (periodic) call cross_boundary(isperiodic,xyzmh_ptmass(1:3,i),ncross)
             xi = xyzmh_ptmass(1,i)
             yi = xyzmh_ptmass(2,i)
             zi = xyzmh_ptmass(3,i)
             if (isnan(xi) .or. isnan(yi) .or. isnan(zi)) then
                call fatal('maketree','NaN in ptmass position, likely caused by NaN in force',i,var='x',val=xi)
             endif
             xminpart = min(xminpart,xi)
             yminpart = min(yminpart,yi)
             zminpart = min(zminpart,zi)
             xmaxpart = max(xmaxpart,xi)
             ymaxpart = max(ymaxpart,yi)
             zmaxpart = max(zmaxpart,zi)
          endif
       enddo
       !$omp enddo
    endif
 endif
 !$omp end parallel

 do i=1,np
    isnotdead: if (.not.isdead_or_accreted(xyzh(4,i))) then
       nproot = nproot + 1

       if (ind_timesteps) then
          if (iactive(iphase(i))) then
             inodeparts(nproot) = i  ! +ve if active
          else
             inodeparts(nproot) = -i ! -ve if inactive
          endif
          if (use_apr) inodeparts(nproot) = abs(inodeparts(nproot))
       else
          inodeparts(nproot) = i
       endif
       treecache(1:4,nproot) = xyzh(1:4,i)
       if (maxphase==maxp) then
          if (use_apr) then
             treecache(5,nproot) = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             treecache(5,nproot) = massoftype(iamtype(iphase(i)))
          endif
       elseif (use_apr) then
          treecache(5,nproot) = aprmassoftype(igas,apr_level(i))
       else
          treecache(5,nproot) = massoftype(igas)
       endif
    endif isnotdead
 enddo

 if (use_sinktree) then
    if (nptmass > 0) then
       do i=1,nptmass
          if (mpi) then
             if (ibelong(maxpsph+i) /= id) cycle
          endif
          if (xyzmh_ptmass(4,i)<0.) cycle
          nproot = nproot + 1
          inodeparts(nproot) = (maxpsph) + i
          treecache(1:3,nproot) = xyzmh_ptmass(1:3,i)
          treecache(4,nproot)   = xyzmh_ptmass(ihsoft,i)
          treecache(5,nproot)   = xyzmh_ptmass(4,i)
       enddo
    endif
 endif

 if (nproot /= 0) then
    inoderange(1,irootnode) = 1
    inoderange(2,irootnode) = nproot
 else
    inoderange(:,irootnode) = 0
 endif

 xmini(1) = xminpart
 xmini(2) = yminpart
 xmini(3) = zminpart
 xmaxi(1) = xmaxpart
 xmaxi(2) = ymaxpart
 xmaxi(3) = zmaxpart

end subroutine construct_root_node

! also used for queue push
pure subroutine push_onto_stack(stackentry,node,parent,level,npnode,xmin,xmax)
 type(kdbuildstack), intent(out) :: stackentry
 integer,            intent(in)  :: node,parent,level
 integer,            intent(in)  :: npnode
 real,               intent(in)  :: xmin(3),xmax(3)

 stackentry%node   = node
 stackentry%parent = parent
 stackentry%level  = level
 stackentry%npnode = npnode
 stackentry%xmin   = xmin
 stackentry%xmax   = xmax

end subroutine push_onto_stack

! also used for queue pop
pure subroutine pop_off_stack(stackentry, istack, nnode, mymum, level, npnode, xmini, xmaxi)
 type(kdbuildstack), intent(in)    :: stackentry
 integer,            intent(inout) :: istack
 integer,            intent(out)   :: nnode, mymum, level, npnode
 real,               intent(out)   :: xmini(3), xmaxi(3)

 nnode  = stackentry%node
 mymum  = stackentry%parent
 level  = stackentry%level
 npnode = stackentry%npnode
 xmini  = stackentry%xmin
 xmaxi  = stackentry%xmax
 istack = istack - 1

end subroutine pop_off_stack

!--------------------------------------------------------------------
!+
!  create all the properties for a given node such as centre of mass,
!  size, max smoothing length, etc
!  will also split the node if necessary, setting wassplit=true
!  returns the left and right child information if split
!+
!--------------------------------------------------------------------
subroutine construct_node(nodeentry, nnode, mymum, level, xmini, xmaxi, npnode, doparallel,&
                          il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                          minlevel, maxlevel, wassplit, global_build,apr_tree, &
                          xyzmh_ptmass)
 use dim,       only:maxtypes,mpi,ind_timesteps
 use part,      only:massoftype,igas,iamtype,npartoftype,isink,ihsoft
 use io,        only:fatal,error
 use mpitree,   only:get_group_cofm,reduce_group
 type(kdnode),      intent(out)   :: nodeentry
 integer,           intent(in)    :: nnode, mymum, level
 real,              intent(inout) :: xmini(3), xmaxi(3)
 integer,           intent(in)    :: npnode
 logical,           intent(in)    :: doparallel
 integer,           intent(out)   :: il, ir, nl, nr
 real,              intent(out)   :: xminl(3), xmaxl(3), xminr(3), xmaxr(3)
 integer(kind=8),   intent(inout) :: ncells
 integer,           intent(out)   :: leaf_is_active(:)
 integer,           intent(inout) :: maxlevel, minlevel
 logical,           intent(out)   :: wassplit
 logical,           intent(in)    :: global_build
 logical,           intent(in)    :: apr_tree
 real,    optional, intent(in)    :: xyzmh_ptmass(:,:)

 integer(kind=8) :: myslot
 real    :: xyzcofm(3)
 real    :: totmass_node
 real    :: xyzcofmg(3)
 real    :: totmassg
 integer :: npnodetot

 logical :: nodeisactive,sinktree
 integer :: i,npcounter,i1,ipart
 real    :: xi,yi,zi,hi,dx,dy,dz,dr2
 real    :: r2max, hmax
 real    :: xcofm,ycofm,zcofm,fac,dfac
 real    :: x0(3)
 integer :: iaxis
 real    :: xpivot
#ifdef GRAVITY
 real    :: quads(6)
#endif
 real    :: pmassi

 sinktree = present(xyzmh_ptmass)
 nodeisactive = .false.
 if (inoderange(1,nnode) > 0) then
    checkactive: do i = inoderange(1,nnode),inoderange(2,nnode)
       if (inodeparts(i) > 0) then
          nodeisactive = .true.
          exit checkactive
       endif
    enddo checkactive
    npcounter = inoderange(2,nnode) - inoderange(1,nnode) + 1
 else
    npcounter = 0
 endif

 if (npcounter /= npnode) then
    print*,'constructing node ',nnode,': found ',npcounter,' particles, expected:',npnode,' particles for this node'
    call fatal('maketree', 'expected number of particles in node differed from actual number')
 endif

 ! following lines to avoid compiler warnings on intent(out) variables
 ir = 0
 il = 0
 nl = 0
 nr = 0
 wassplit = .false.
 if ((.not. global_build) .and. (npnode  <  1)) return ! node has no particles, just quit

 r2max = 0.
 hmax  = 0.
 xyzcofm(:) = 0.
 xcofm = 0.
 ycofm = 0.
 zcofm = 0.
!
! to avoid round off error from repeated multiplication by pmassi (which is small)
! we compute the centre of mass with a factor relative to gas particles
! but only if gas particles are present
!
 pmassi = massoftype(igas)
 fac    = 1.
 totmass_node = 0.
 if (pmassi > 0.) then
    dfac = 1./pmassi
 else
    pmassi = massoftype(maxloc(npartoftype(2:maxtypes),1)+1)
    if (pmassi > 0.) then
       dfac = 1./pmassi
    else
       dfac = 1.
    endif
 endif
 ! note that dfac can be a constant value across all particles even if APR is used

 i1=inoderange(1,nnode)
 ! during initial queue build which is serial, we can parallelise this loop
 if (npnode > 1000 .and. doparallel) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(maxp,maxphase,maxpsph,inodeparts) &
    !$omp shared(npnode,massoftype,dfac,aprmassoftype) &
    !$omp shared(treecache,i1) &
    !$omp shared(xyzmh_ptmass,sinktree) &
    !$omp private(i,xi,yi,zi,hi) &
    !$omp firstprivate(pmassi,fac) &
    !$omp reduction(+:xcofm,ycofm,zcofm,totmass_node) &
    !$omp reduction(max:hmax)
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       hi = treecache(4,i)
       pmassi = treecache(5,i)
       fac    = pmassi*dfac ! to avoid round-off error
       hmax  = max(hmax,hi)
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       hi = treecache(4,i)
       pmassi = treecache(5,i)
       fac    = pmassi*dfac ! to avoid round-off error
       hmax  = max(hmax,hi)
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
 endif

 xyzcofm = (/xcofm,ycofm,zcofm/)

 ! if there are no particles in this node, then the cofm will
 ! remain at zero
 if (totmass_node > 0.) then
    xyzcofm(:)   = xyzcofm(:)/(totmass_node*dfac)
 endif

 ! if this is global node construction, get the cofm and total mass
 ! of all particles in this node (some on other MPI tasks)
 if (mpi .and. global_build) then
    call get_group_cofm(xyzcofm,totmass_node,level,xyzcofmg,totmassg)
    xyzcofm = xyzcofmg
    totmass_node = totmassg
 endif

 ! checks the reduced mass in the case of global maketree
 if (totmass_node<=0. .and. use_apr) call fatal('mtree + apr', &
    'totmass_node==0, something almost certainly wrong with aprmassoftype')
 if (totmass_node<=0.) call fatal('mtree','totmass_node==0',val=totmass_node)

!--for gravity, we need the centre of the node to be the centre of mass
 x0(:) = xyzcofm(:)
 r2max = 0.
#ifdef GRAVITY
 quads(:) = 0.
#endif

 !--compute size of node
 ! parallelise this loop if node is large enough
 ! use !$omp parallel do when doparallel=.true. (not in parallel region)
 ! when doparallel=.false., we're already in a parallel region but can't use nested reductions
 ! so we'll use thread-local accumulators and combine at the end
 if (npnode > 1000 .and. doparallel) then
    !$omp parallel do schedule(static) default(none) &
    !$omp shared(npnode,treecache,x0,i1,maxp) &
    !$omp shared(massoftype,sinktree,maxphase,maxpsph,inodeparts) &
    !$omp shared(xyzmh_ptmass,aprmassoftype) &
    !$omp private(i,xi,yi,zi,dx,dy,dz,dr2) &
    !$omp firstprivate(pmassi) &
#ifdef GRAVITY
    !$omp reduction(+:quads) &
#endif
    !$omp reduction(max:r2max)
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       dx    = xi - x0(1)
       dy    = yi - x0(2)
       dz    = zi - x0(3)
       dr2   = dx*dx + dy*dy + dz*dz
       r2max = max(r2max,dr2)
#ifdef GRAVITY
       pmassi = treecache(5,i)
       quads(1) = quads(1) + pmassi*(dx*dx)  ! Q_xx
       quads(2) = quads(2) + pmassi*(dx*dy)  ! Q_xy = Q_yx
       quads(3) = quads(3) + pmassi*(dx*dz)  ! Q_xz = Q_zx
       quads(4) = quads(4) + pmassi*(dy*dy)  ! Q_yy
       quads(5) = quads(5) + pmassi*(dy*dz)  ! Q_yz = Q_zy
       quads(6) = quads(6) + pmassi*(dz*dz)  ! Q_zz
#endif
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = treecache(1,i)
       yi = treecache(2,i)
       zi = treecache(3,i)
       dx    = xi - x0(1)
       dy    = yi - x0(2)
       dz    = zi - x0(3)
       dr2   = dx*dx + dy*dy + dz*dz
       r2max = max(r2max,dr2)
#ifdef GRAVITY
       pmassi = treecache(5,i)
       quads(1) = quads(1) + pmassi*(dx*dx)  ! Q_xx
       quads(2) = quads(2) + pmassi*(dx*dy)  ! Q_xy = Q_yx
       quads(3) = quads(3) + pmassi*(dx*dz)  ! Q_xz = Q_zx
       quads(4) = quads(4) + pmassi*(dy*dy)  ! Q_yy
       quads(5) = quads(5) + pmassi*(dy*dz)  ! Q_yz = Q_zy
       quads(6) = quads(6) + pmassi*(dz*dz)  ! Q_zz
#endif
    enddo
 endif

 ! reduce node limits and quads across MPI tasks belonging to this group
 if (mpi .and. global_build) then
    npnodetot = reduce_group(npnode,'+',level)
    r2max     = reduce_group(r2max,'max',level)
    hmax      = reduce_group(hmax,'max',level)

    xmini(1)  = reduce_group(xmini(1),'min',level)
    xmini(2)  = reduce_group(xmini(2),'min',level)
    xmini(3)  = reduce_group(xmini(3),'min',level)

    xmaxi(1)  = reduce_group(xmaxi(1),'max',level)
    xmaxi(2)  = reduce_group(xmaxi(2),'max',level)
    xmaxi(3)  = reduce_group(xmaxi(3),'max',level)
#ifdef GRAVITY
    quads(1)  = reduce_group(quads(1),'+',level)
    quads(2)  = reduce_group(quads(2),'+',level)
    quads(3)  = reduce_group(quads(3),'+',level)
    quads(4)  = reduce_group(quads(4),'+',level)
    quads(5)  = reduce_group(quads(5),'+',level)
    quads(6)  = reduce_group(quads(6),'+',level)
#endif
 else
    npnodetot = npnode
 endif

 ! assign properties to node
 nodeentry%xcen    = x0(:)
 nodeentry%size    = sqrt(r2max) + epsilon(r2max)
 nodeentry%hmax    = hmax
 nodeentry%parent  = mymum
#ifdef GRAVITY
 nodeentry%mass    = totmass_node
 nodeentry%quads   = quads
#endif

 wassplit = (npnodetot > minpart)
 if (apr_tree) wassplit = (npnode > 2)

 if (.not. wassplit) then
    nodeentry%leftchild  = 0
    nodeentry%rightchild = 0
    maxlevel = max(level,maxlevel)
    minlevel = min(level,minlevel)
    ! individual timesteps where we mark leaf node as active/inactive
    if (ind_timesteps) then
       !
       !--mark leaf node as active (contains some active particles)
       !  or inactive by setting the firstincell to +ve (active) or -ve (inactive)
       !
       if (nodeisactive) then
          leaf_is_active(nnode) = 1
       else
          leaf_is_active(nnode) = -1
       endif
    else
       leaf_is_active(nnode) = 1
    endif
 else ! split this node and add children to stack
    iaxis  = maxloc(xmaxi - xmini,1) ! split along longest axis
    xpivot = xyzcofm(iaxis)          ! split on centre of mass

    ! create two children nodes and point to them from current node
    ! always use G&R indexing for global tree
    if ((level < maxlevel_indexed) .or. global_build) then
       il = 2*nnode   ! indexing as per Gafton & Rosswog (2011)
       ir = il + 1
    else
       ! no need to lock, we could just atomic the update
       !$omp atomic capture
       ncells = ncells + 2
       myslot = ncells
       !$omp end atomic
       ir = int(myslot)
       il = int(myslot-1)
       if (ir > ncellsmax) call fatal('maketree',&
          'number of nodes exceeds array dimensions, increase ncellsmax and recompile',ival=int(ncellsmax))
    endif
    nodeentry%leftchild  = il
    nodeentry%rightchild = ir

    leaf_is_active(nnode) = 0

    if (npnode > 0) then
       if (apr_tree) then
          ! apr special sort - only used for merging particles
          call special_sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                    inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,treecache,inodeparts,&
                                    npnode)
       else
          ! regular sort
          call sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                  inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,treecache,inodeparts)
       endif

       if (nr + nl  /=  npnode) then
          call error('maketree','number of left + right != parent while splitting (likely cause: NaNs in position arrays)')
       endif

       ! see if all the particles ended up in one node, if so, arbitrarily build 2 cells
       if ( (.not. global_build) .and. ((nl==npnode) .or. (nr==npnode)) ) then
          ! no need to move particles because if they all ended up in one node,
          ! then they are still in the original order
          nl = npnode / 2
          inoderange(1,il) = inoderange(1,nnode)
          inoderange(2,il) = inoderange(1,nnode) + nl - 1
          inoderange(1,ir) = inoderange(1,nnode) + nl
          inoderange(2,ir) = inoderange(2,nnode)
          nr = npnode - nl
       endif

       ! compute min/max with explicit loops for better cache behavior
       xminl(1) = treecache(1,inoderange(1,il))
       xminl(2) = treecache(2,inoderange(1,il))
       xminl(3) = treecache(3,inoderange(1,il))
       xmaxl(1) = xminl(1)
       xmaxl(2) = xminl(2)
       xmaxl(3) = xminl(3)
       do ipart=inoderange(1,il)+1,inoderange(2,il)
          xminl(1) = min(xminl(1),treecache(1,ipart))
          xminl(2) = min(xminl(2),treecache(2,ipart))
          xminl(3) = min(xminl(3),treecache(3,ipart))
          xmaxl(1) = max(xmaxl(1),treecache(1,ipart))
          xmaxl(2) = max(xmaxl(2),treecache(2,ipart))
          xmaxl(3) = max(xmaxl(3),treecache(3,ipart))
       enddo

       xminr(1) = treecache(1,inoderange(1,ir))
       xminr(2) = treecache(2,inoderange(1,ir))
       xminr(3) = treecache(3,inoderange(1,ir))
       xmaxr(1) = xminr(1)
       xmaxr(2) = xminr(2)
       xmaxr(3) = xminr(3)
       do ipart=inoderange(1,ir)+1,inoderange(2,ir)
          xminr(1) = min(xminr(1),treecache(1,ipart))
          xminr(2) = min(xminr(2),treecache(2,ipart))
          xminr(3) = min(xminr(3),treecache(3,ipart))
          xmaxr(1) = max(xmaxr(1),treecache(1,ipart))
          xmaxr(2) = max(xmaxr(2),treecache(2,ipart))
          xmaxr(3) = max(xmaxr(3),treecache(3,ipart))
       enddo
    else
       nl = 0
       nr = 0
       xminl = 0.0
       xmaxl = 0.0
       xminr = 0.0
       xmaxr = 0.0
    endif

    ! Reduce node limits of children across MPI tasks belonging to this group.
    ! The synchronisation needs to happen here, not at the next level, because
    ! the groups will be independent by then.
    if (mpi .and. global_build) then
       xminl(1) = reduce_group(xminl(1),'min',level)
       xminl(2) = reduce_group(xminl(2),'min',level)
       xminl(3) = reduce_group(xminl(3),'min',level)

       xmaxl(1) = reduce_group(xmaxl(1),'max',level)
       xmaxl(2) = reduce_group(xmaxl(2),'max',level)
       xmaxl(3) = reduce_group(xmaxl(3),'max',level)

       xminr(1) = reduce_group(xminr(1),'min',level)
       xminr(2) = reduce_group(xminr(2),'min',level)
       xminr(3) = reduce_group(xminr(3),'min',level)

       xmaxr(1) = reduce_group(xmaxr(1),'max',level)
       xmaxr(2) = reduce_group(xmaxr(2),'max',level)
       xmaxr(3) = reduce_group(xmaxr(3),'max',level)
    endif

 endif

end subroutine construct_node

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis
!+
!----------------------------------------------------------------
subroutine sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,nl,nr,xpivot,&
                                   treecache,inodeparts)
 integer, intent(in)  :: iaxis,imin,imax
 integer, intent(out) :: min_l,max_l,min_r,max_r,nl,nr
 real, intent(inout)  :: xpivot,treecache(:,:)
 integer,         intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot
 integer :: inodeparts_swap,i,j
 real :: xyzh_swap(5)
 real :: xi_coord, xj_coord

 !print*,'nnode ',imin,imax,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 xi_coord = treecache(iaxis,i)
 xj_coord = treecache(iaxis,j)
 i_lt_pivot = xi_coord <= xpivot
 j_lt_pivot = xj_coord <= xpivot
 !  k = 0

 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       xi_coord = treecache(iaxis,i)
       i_lt_pivot = xi_coord <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          xj_coord = treecache(iaxis,j)
          j_lt_pivot = xj_coord <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:5)  = treecache(1:5,i)

          inodeparts(i)   = inodeparts(j)
          treecache(1:5,i) = treecache(1:5,j)

          inodeparts(j)   = inodeparts_swap
          treecache(1:5,j) = xyzh_swap(1:5)

          i = i + 1
          j = j - 1
          xi_coord = treecache(iaxis,i)
          xj_coord = treecache(iaxis,j)
          i_lt_pivot = xi_coord <= xpivot
          j_lt_pivot = xj_coord <= xpivot
          ! k = k + 1
       endif
    endif
 enddo
 if (.not.i_lt_pivot) i = i - 1
 if (j_lt_pivot)      j = j + 1

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine sort_particles_in_cell

!----------------------------------------------------------------
!+
!  Categorise particles into daughter nodes by whether they
!  fall to the left or the right of the pivot axis, but additionally
!  force the cells to have a certain minimum number of particles per cell
!+
!----------------------------------------------------------------
subroutine special_sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,&
                                nl,nr,xpivot,treecache,inodeparts,npnode)
 use io, only:error
 integer, intent(in)  :: iaxis,imin,imax,npnode
 integer, intent(out) :: min_l,max_l,min_r,max_r,nl,nr
 real, intent(inout)  :: xpivot,treecache(:,:)
 integer,         intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot,slide_l,slide_r
 integer :: inodeparts_swap,i,j,nchild_in
 integer :: k,ii,rem_nr,rem_nl
 real :: xyzh_swap(5),dpivot(npnode)

 dpivot = 0.0
 nchild_in = 2

 if (modulo(npnode,nchild_in) > 0) then
    call error('apr sort','number of particles sent in to kdtree is not divisible by 2')
 endif

! print*,'nnode ',imin,imax,npnode,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 i_lt_pivot = treecache(iaxis,i) <= xpivot
 j_lt_pivot = treecache(iaxis,j) <= xpivot
 dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
 !k = 0
 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
       i_lt_pivot = treecache(iaxis,i) <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
          j_lt_pivot = treecache(iaxis,j) <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:5)  = treecache(1:5,i)

          inodeparts(i)   = inodeparts(j)
          treecache(1:5,i) = treecache(1:5,j)

          inodeparts(j)   = inodeparts_swap
          treecache(1:5,j) = xyzh_swap(1:5)

          i = i + 1
          j = j - 1

          dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
          dpivot(j-imin+1) = xpivot - treecache(iaxis,j)

          i_lt_pivot = treecache(iaxis,i) <= xpivot
          j_lt_pivot = treecache(iaxis,j) <= xpivot
       endif
    endif
 enddo

 if (.not.i_lt_pivot) then
    i = i - 1
    dpivot(i-imin+1) = xpivot - treecache(iaxis,i)
 endif
 if (j_lt_pivot) then
    j = j + 1
    dpivot(j-imin+1) = xpivot - treecache(iaxis,j)
 endif

 min_l = imin
 max_l = i
 min_r = j
 max_r = imax

 if ( j /= i+1) print*,' ERROR ',i,j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

 ! does the pivot need to be adjusted?
 rem_nl = modulo(nl,nchild_in)
 rem_nr = modulo(nr,nchild_in)
 if (rem_nl == 0 .and. rem_nr == 0) return

 ! Decide which direction the pivot needs to go
 if (rem_nl < rem_nr) then
    slide_l = .true.
    slide_r = .false.
 else
    slide_l = .false.
    slide_r = .true.
 endif
 ! Override this if there's less than nchild*2 in the cell
 if (nl < nchild_in) then
    slide_r = .true.
    slide_l = .false.
 elseif (nr < nchild_in) then
    slide_r = .false.
    slide_l = .true.
 endif

 ! Move across particles by distance from xpivot till we get
 ! the right number of particles in each cell
 if (slide_r) then
    do ii = 1,rem_nr
       ! next particle to shift across
       k = minloc(dpivot,dim=1,mask=dpivot > 0.) + imin - 1
       if (k-imin+1==0) k = maxloc(dpivot,dim=1,mask=dpivot < 0.) + imin - 1

       ! swap this with the first particle on the j side
       inodeparts_swap = inodeparts(k)
       xyzh_swap(1:5)  = treecache(1:5,k)

       inodeparts(k)   = inodeparts(j)
       treecache(1:5,k) = treecache(1:5,j)

       inodeparts(j)   = inodeparts_swap
       treecache(1:5,j) = xyzh_swap(1:5)

       ! and now shift to the right
       i = i + 1
       j = j + 1

       ! ditch it, go again
       dpivot(k-imin+1) = huge(k-imin+1)
    enddo
 else
    do ii = 1,rem_nl
       ! next particle to shift across
       k = maxloc(dpivot,dim=1,mask=dpivot < 0.) + imin - 1
       if (k-imin+1==0) k = minloc(dpivot,dim=1,mask=dpivot > 0.) + imin - 1

       ! swap this with the last particle on the i side
       inodeparts_swap = inodeparts(k)
       xyzh_swap(1:5)  = treecache(1:5,k)

       inodeparts(k)   = inodeparts(i)
       treecache(1:5,k) = treecache(1:5,i)

       inodeparts(i)   = inodeparts_swap
       treecache(1:5,i) = xyzh_swap(1:5)

       ! and now shift to the left
       i = i - 1
       j = j - 1

       ! ditch it, go again
       dpivot(k-imin+1) = huge(k-imin+1)

    enddo
 endif

 ! tidy up outputs
 max_l = i
 min_r = j
 nl = max_l - min_l + 1
 nr = max_r - min_r + 1

end subroutine special_sort_particles_in_cell

!----------------------------------------------------------------
!+
!  Cache particles within identified neighbour nodes
!+
!----------------------------------------------------------------
subroutine cache_neighbours(nneigh,isrc,ixyzcachesize,maxcache,listneigh,xyzcache,xoffset,yoffset,zoffset)
 integer, intent(in)    :: isrc,ixyzcachesize,maxcache
 real,    intent(in)    :: xoffset,yoffset,zoffset
 integer, intent(inout) :: nneigh
 integer, intent(out)   :: listneigh(:)
 real,    intent(out)   :: xyzcache(:,:)
 integer :: npnode,ipart,num_to_cache

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
       listneigh(nneigh+ipart)  = abs(inodeparts(inoderange(1,isrc)+ipart-1))
       xyzcache(1,nneigh+ipart) = treecache(1,inoderange(1,isrc)+ipart-1) + xoffset
       xyzcache(2,nneigh+ipart) = treecache(2,inoderange(1,isrc)+ipart-1) + yoffset
       xyzcache(3,nneigh+ipart) = treecache(3,inoderange(1,isrc)+ipart-1) + zoffset
       if (maxcache >= 4) then
          xyzcache(4,nneigh+ipart) = 1./treecache(4,inoderange(1,isrc)+ipart-1)
       endif
    enddo
 endif

 if (num_to_cache < npnode) then
    do ipart=num_to_cache+1,npnode
       listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,isrc)+ipart-1))
    enddo
 endif

 nneigh = nneigh + npnode

end subroutine cache_neighbours

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
 type(kdnode), intent(in)           :: node(:) !ncellsmax+1)
 integer, intent(in)                :: ixyzcachesize
 real,    intent(in)                :: xpos(3)
 real,    intent(in)                :: xsizei,rcuti
 integer, intent(out)               :: listneigh(:)
 integer, intent(out)               :: nneigh
 real,    intent(out)               :: xyzcache(:,:)
 integer, intent(in)                :: leaf_is_active(:)
 logical, intent(in)                :: get_hj
 logical, intent(in)                :: get_f
 real,    intent(out),    optional  :: fnode(lenfgrav)
 logical, intent(out),    optional  :: remote_export(:)
 integer, intent(in),     optional  :: nq
 integer :: maxcache
 integer :: n,istack,il,ir
 integer :: nstack(maxdepth)
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
          if (istack+2 > ncellsmax+1) call fatal('getneigh','stack overflow in getneigh')
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
 use io,       only:fatal
 type(kdnode), intent(in)   :: node(:) !ncellsmax+1)
 integer,      intent(in)   :: ixyzcachesize
 real,         intent(in)   :: xpos(3)
 real,         intent(in)   :: xsizei,rcuti
 integer,      intent(out)  :: listneigh(:)
 integer,      intent(out)  :: nneigh
 real,         intent(out)  :: xyzcache(:,:)
 integer,      intent(in)   :: leaf_is_active(:)
 logical,      intent(in)   :: get_hj
 logical,      intent(in)   :: get_f
 real,         intent(out)  :: fnode(lenfgrav)
 integer,      intent(in)   :: icell
 integer :: istack,i,idstbranch,idst,isrc,maxcache
 integer :: branch(maxdepth),nparents,stack(3,maxdepth)
 real    :: dx,dy,dz,xoffset,yoffset,zoffset
 real    :: tree_acc2
 logical :: stackit

 tree_acc2 = tree_accuracy*tree_accuracy

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache,1)
 else
    maxcache = 0
 endif

 call get_list_of_parent_nodes(icell,node,branch,nparents)

 fnode_branch = 0.

 nneigh = 0
 istack = 1
 stack(1,istack) = irootnode
 stack(2,istack) = irootnode
 stack(3,istack) = nparents ! root id in the branch

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
       call node_interaction(node(idst),node(isrc),tree_acc2,fnode_branch(:,idstbranch),stackit,xoffset,yoffset,zoffset)
    endif

    if (stackit) then
       call open_nodes(stack,istack,node(isrc),isrc,branch,idstbranch,&
                       listneigh,xyzcache,ixyzcachesize,nneigh,leaf_is_active,&
                       maxcache,xoffset,yoffset,zoffset)
    endif
 enddo

 !
 !-- Downward pass to accumulate on each leaf
 !
 do i=nparents,2,-1 ! parents(1) is equal to icell
    call get_sep(node(branch(i-1))%xcen,node(branch(i))%xcen,dx,dy,dz,xoffset,yoffset,zoffset)
    call propagate_fnode_to_node(fnode_branch(:,i-1),fnode_branch(:,i),dx,dy,dz)
 enddo

 !-- final result is accumulated in the first column of fnode_branch -> store into fnode to be used in force
 fnode = fnode_branch(:,1)

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
 real,           intent(in)  :: x1(3),x2(3)
 real,           intent(out) :: dx,dy,dz,xoffset,yoffset,zoffset
 real, optional, intent(out) :: r2

 xoffset = 0.
 yoffset = 0.
 zoffset = 0.

 dx = x1(1) - x2(1)
 dy = x1(2) - x2(2)
 dz = x1(3) - x2(3)

#ifdef PERIODIC
 if (abs(dx) > hdlx) then ! mod distances across boundary if periodic BCs
    xoffset = dxbound*SIGN(1.0,dx)
    dx = dx - xoffset
 endif
 if (abs(dy) > hdly) then
    yoffset = dybound*SIGN(1.0,dy)
    dy = dy - yoffset
 endif
 if (abs(dz) > hdlz) then
    zoffset = dzbound*SIGN(1.0,dz)
    dz = dz - zoffset
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
 type(kdnode),   intent(in)  :: node_dst,node_src
 real,           intent(out) :: size_src,size_dst
 real,           intent(out) :: rcut_src,rcut_dst

 rcut_src = node_src%hmax*radkern
 rcut_dst = node_dst%hmax*radkern
 size_src = node_src%size
 size_dst = node_dst%size

end subroutine get_node_size

!-----------------------------------------------------------
!+
!  Taylor expand the contribution from direct parent nodes
!  to the child node centre
!+
!-----------------------------------------------------------
pure subroutine propagate_fnode_to_node(fnode,fnode_sup,dx,dy,dz)
 real, intent(in)    :: fnode_sup(lenfgrav),dx,dy,dz
 real, intent(inout) :: fnode(lenfgrav)

 fnode(1)  = fnode(1)  + fnode_sup(1) + dx*(fnode_sup(4) + 0.5*(dx*fnode_sup(10) + dy*fnode_sup(11) +dz*fnode_sup(12)))& ! xx +0.5(xxx+xxy+xxz)
                       + dy*(fnode_sup(5) + 0.5*(dx*fnode_sup(11) + dy*fnode_sup(13) +dz*fnode_sup(14)))& ! xy +0.5(xxy+xyy+xyz)
                       + dz*(fnode_sup(6) + 0.5*(dx*fnode_sup(12) + dy*fnode_sup(14) +dz*fnode_sup(15)))  ! xz +0.5(xxz+xyz+xzz)
 fnode(2)  = fnode(2)  + fnode_sup(2) + dx*(fnode_sup(5) + 0.5*(dx*fnode_sup(11) + dy*fnode_sup(13) +dz*fnode_sup(14)))& ! xy +0.5(xxy+xyy+xyz)
                       + dy*(fnode_sup(7) + 0.5*(dx*fnode_sup(13) + dy*fnode_sup(16) +dz*fnode_sup(17)))& ! yy +0.5(xyy+yyy+yyz)
                       + dz*(fnode_sup(8) + 0.5*(dx*fnode_sup(14) + dy*fnode_sup(17) +dz*fnode_sup(18)))  ! yz +0.5(xyz+yyz+yyz)
 fnode(3)  = fnode(3)  + fnode_sup(3) + dx*(fnode_sup(6) + 0.5*(dx*fnode_sup(12) + dy*fnode_sup(14) +dz*fnode_sup(15)))& ! xz +0.5(xxz+xyz+xzz)
                       + dy*(fnode_sup(8) + 0.5*(dx*fnode_sup(14) + dy*fnode_sup(17) +dz*fnode_sup(18)))& ! yz +0.5(xyz+yyz+yzz)
                       + dz*(fnode_sup(9) + 0.5*(dx*fnode_sup(15) + dy*fnode_sup(18) +dz*fnode_sup(19)))  ! zz +0.5(xzz+yzz+zzz)
 fnode(4)  = fnode(4)  + fnode_sup(4) + dx*fnode_sup(10) + dy*fnode_sup(11) + dz*fnode_sup(12)                           ! xxx + xxy + xxz
 fnode(5)  = fnode(5)  + fnode_sup(5) + dx*fnode_sup(11) + dy*fnode_sup(13) + dz*fnode_sup(14)                           ! xxy + xyy + xyz
 fnode(6)  = fnode(6)  + fnode_sup(6) + dx*fnode_sup(12) + dy*fnode_sup(14) + dz*fnode_sup(15)                           ! xxz + xyz + xzz
 fnode(7)  = fnode(7)  + fnode_sup(7) + dx*fnode_sup(13) + dy*fnode_sup(16) + dz*fnode_sup(17)                           ! xyy + yyy + yyz
 fnode(8)  = fnode(8)  + fnode_sup(8) + dx*fnode_sup(14) + dy*fnode_sup(17) + dz*fnode_sup(18)                           ! xyz + yyz + yzz
 fnode(9)  = fnode(9)  + fnode_sup(9) + dx*fnode_sup(15) + dy*fnode_sup(18) + dz*fnode_sup(19)                           ! xzz + yzz + zzz
 fnode(10) = fnode(10) + fnode_sup(10)
 fnode(11) = fnode(11) + fnode_sup(11)
 fnode(12) = fnode(12) + fnode_sup(12)
 fnode(13) = fnode(13) + fnode_sup(13)
 fnode(14) = fnode(14) + fnode_sup(14)
 fnode(15) = fnode(15) + fnode_sup(15)
 fnode(16) = fnode(16) + fnode_sup(16)
 fnode(17) = fnode(17) + fnode_sup(17)
 fnode(18) = fnode(18) + fnode_sup(18)
 fnode(19) = fnode(19) + fnode_sup(19)
 fnode(20) = fnode(20) + fnode_sup(20) + dx*(fnode_sup(1)+0.5*(dx*fnode_sup(4)+dy*fnode_sup(5)+dz*fnode_sup(6)))&
                                       + dy*(fnode_sup(2)+0.5*(dx*fnode_sup(5)+dy*fnode_sup(7)+dz*fnode_sup(8)))&
                                       + dz*(fnode_sup(3)+0.5*(dx*fnode_sup(6)+dy*fnode_sup(8)+dz*fnode_sup(9)))

end subroutine propagate_fnode_to_node

!-----------------------------------------------------------
!+
!  return list of parents of current node
!+
!-----------------------------------------------------------
pure subroutine get_list_of_parent_nodes(inode,node,parents,nparents)
 integer,      intent(in)  :: inode
 type(kdnode), intent(in)  :: node(:)
 integer,      intent(out) :: parents(:)
 integer,      intent(out) :: nparents
 integer :: j

 j = inode
 nparents = 1
 parents  = 0
 parents(nparents) = j ! set first elem to inode to use parents for propagation
 do while (node(j)%parent  /=  0)
    j = node(j)%parent
    nparents = nparents + 1
    parents(nparents) = j
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
 type(kdnode), intent(in)     :: srcnode
 integer,      intent(in)     :: isrc,idstbranch
 integer,      intent(in)     :: branch(:)
 integer,      intent(in)     :: ixyzcachesize,maxcache
 integer,      intent(in)     :: leaf_is_active(:)
 integer,      intent(inout)  :: listneigh(:)
 integer,      intent(inout)  :: nneigh
 integer,      intent(inout)  :: stack(:,:),istack
 real,         intent(inout)  :: xyzcache(:,:)
 real,         intent(in)     :: xoffset,yoffset,zoffset
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
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = isrc
       stack(3,istack) = ibranchnext
    endif is_P2P
 else
    if (il /= 0) then
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = il
       stack(3,istack) = ibranchnext
    endif
    if (ir /= 0) then
       istack = istack + 1
       stack(1,istack) = idstnext
       stack(2,istack) = ir
       stack(3,istack) = ibranchnext
    endif
 endif is_src_leaf

end subroutine open_nodes

!-----------------------------------------------------------
!+
!  Test the separation between the node pair and compute
!  the interaction if needed
!+
!-----------------------------------------------------------
pure subroutine node_interaction(node_dst,node_src,tree_acc2,fnode,stackit,xoffset,yoffset,zoffset)
 type(kdnode), intent(in)    :: node_dst,node_src
 real,         intent(in)    :: tree_acc2
 real,         intent(inout) :: fnode(lenfgrav)
 real,         intent(out)   :: xoffset,yoffset,zoffset
 logical,      intent(out)   :: stackit
 real    :: dx,dy,dz,r2,dr1
 real    :: rcut_dst,rcut_src,rcut,rcut2
 real    :: size_dst,size_src,mass_src,quads_src(6)
 logical :: wellsep

 call get_sep(node_dst%xcen,node_src%xcen,dx,dy,dz,xoffset,yoffset,zoffset,r2)
 call get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)

 rcut  = max(rcut_dst,rcut_src)
 rcut2 = (size_dst+size_src+rcut)**2
 wellsep = (tree_acc2*r2 > (size_dst+size_src)**2) .and. (r2 > rcut2)

 if (wellsep) then
    dr1 = 1./sqrt(r2)
#ifdef GRAVITY
    mass_src=node_src%mass
    quads_src=node_src%quads
#else
    mass_src=0.
    quads_src=0.
#endif
    call compute_M2L(dx,dy,dz,dr1,mass_src,quads_src,fnode)
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
 g0   = - dr1
 g1   =  1.*dr12*g0
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

 fnode(1)  = fnode(1)  + (D1(1)*q0  +&
                     0.5*(D3(1)*qxx + 2.*(D3(2)*qxy + D3(3)*qxz + D3(5)*qyz) + D3(4)*qyy + D3(6)*qzz ))    ! C¹_x
 fnode(2)  = fnode(2)  + (D1(2)*q0  +&
                     0.5*(D3(2)*qxx + 2.*(D3(4)*qxy + D3(5)*qxz + D3(8)*qyz) + D3(7)*qyy + D3(9)*qzz ))    ! C¹_y
 fnode(3)  = fnode(3)  + (D1(3)*q0  +&
                     0.5*(D3(3)*qxx + 2.*(D3(5)*qxy + D3(6)*qxz + D3(9)*qyz) + D3(8)*qyy + D3(10)*qzz))   ! C¹_z
 fnode(4)  = fnode(4)  + D2(1) * q0    ! C²_xx
 fnode(5)  = fnode(5)  + D2(2) * q0    ! C²_xy
 fnode(6)  = fnode(6)  + D2(3) * q0    ! C²_xz
 fnode(7)  = fnode(7)  + D2(4) * q0    ! C²_yy
 fnode(8)  = fnode(8)  + D2(5) * q0    ! C²_yz
 fnode(9)  = fnode(9)  + D2(6) * q0    ! C²_zz
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
 fnode(20) = fnode(20) + g0*q0 - 0.5*(D2(1)*qxx + D2(4)*qyy + D2(6)*qzz + 2*(D2(2)*qxy + D2(3)*qxz + D2(5)*qyz))! C⁰ (potential)

end subroutine compute_M2L

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

 fxi = fxi + dx*(dfxx + 0.5*(dx*d2fxxx + dy*d2fxxy + dz*d2fxxz)) &
           + dy*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dz*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz))
 fyi = fyi + dx*(dfxy + 0.5*(dx*d2fxxy + dy*d2fxyy + dz*d2fxyz)) &
           + dy*(dfyy + 0.5*(dx*d2fxyy + dy*d2fyyy + dz*d2fyyz)) &
           + dz*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz))
 fzi = fzi + dx*(dfxz + 0.5*(dx*d2fxxz + dy*d2fxyz + dz*d2fxzz)) &
           + dy*(dfyz + 0.5*(dx*d2fxyz + dy*d2fyyz + dz*d2fyzz)) &
           + dz*(dfzz + 0.5*(dx*d2fxzz + dy*d2fyzz + dz*d2fzzz))
 poti = poti - (dx*(fxi - 0.5*(dx*dfxx + dy*dfxy + dz*dfxy)) + &
                dy*(fyi - 0.5*(dx*dfxy + dy*dfyy + dz*dfyz)) + &
                dz*(fzi - 0.5*(dx*dfxz + dy*dfyz + dz*dfzz)))

end subroutine expand_fgrav_in_taylor_series

!-----------------------------------------------
!+
!  Routine to update a constructed tree
!  Note: current version ONLY works if
!  tree is built to < maxlevel_indexed
!  That is, it relies on the 2^n style tree
!  indexing to sweep out each level
!+
!-----------------------------------------------
subroutine revtree(node, xyzh, leaf_is_active, ncells)
 use dim,  only:maxp,use_apr,ind_timesteps
 use part, only:maxphase,iphase,igas,massoftype,iamtype,aprmassoftype,&
                apr_level,iactive,treecache,isdead_or_accreted
 use io,   only:fatal
 type(kdnode), intent(inout) :: node(:) !ncellsmax+1)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(inout) :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(in) :: ncells
 real :: hmax, r2max
 real :: xi, yi, zi, hi
 real :: dx, dy, dz, dr2
#ifdef GRAVITY
 real :: quads(6)
#endif
 integer :: inode, ipart, ipartidx, i, nptot
 real :: pmassi, totmass
 real :: x0(3)
 real :: xcofm, ycofm, zcofm, fac, dfac
 logical :: nodeisactive
 pmassi = massoftype(igas)

 ! find maximum index in inodeparts that we need to update in treecache
 nptot = 0
 do i=1,int(ncells)
    if (i > 1 .and. node(i)%parent == 0) cycle
    if (inoderange(1,i) > 0 .and. inoderange(2,i) >= inoderange(1,i)) then
       nptot = max(nptot, inoderange(2,i))
    endif
 enddo

 ! update treecache for particles in the tree only
 ! mark dead/accreted particles by setting treecache(4,i) negative
 !$omp parallel default(none) &
 !$omp shared(nptot,inodeparts,xyzh,iphase,apr_level) &
 !$omp shared(massoftype,aprmassoftype,treecache) &
 !$omp shared(maxphase,maxp) &
 !$omp private(i,ipartidx)
 !$omp do schedule(static)
 do i=1,nptot
    if (inodeparts(i) == 0) cycle
    ipartidx = abs(inodeparts(i))
    treecache(1:4,i) = xyzh(1:4,ipartidx)
    ! compute and store mass
    if (maxphase==maxp) then
       if (use_apr) then
          treecache(5,i) = aprmassoftype(iamtype(iphase(ipartidx)),apr_level(ipartidx))
       else
          treecache(5,i) = massoftype(iamtype(iphase(ipartidx)))
       endif
    elseif (use_apr) then
       treecache(5,i) = aprmassoftype(igas,apr_level(ipartidx))
    else
       treecache(5,i) = massoftype(igas)
    endif
 enddo
 !$omp enddo
 !$omp end parallel

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(ncells) &
!$omp shared(node,inoderange,inodeparts,treecache,leaf_is_active) &
!$omp private(hmax,r2max,xi,yi,zi,hi) &
!$omp private(dx,dy,dz,dr2,inode,ipart,x0) &
!$omp private(xcofm,ycofm,zcofm,fac,dfac,nodeisactive) &
#ifdef GRAVITY
!$omp private(quads) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)
!$omp do schedule(guided)
 over_nodes: do inode=1,int(ncells)
    if (inode > 1 .and. node(inode)%parent == 0) cycle
    ! initialize node properties
    node(inode)%xcen(:) = 0.
    node(inode)%size    = 0.
    node(inode)%hmax    = 0.
#ifdef GRAVITY
    node(inode)%mass    = 0.
    node(inode)%quads(:)= 0.
#endif
    ! initialize leaf_is_active (will be set for leaf nodes below)
    leaf_is_active(inode) = 0

    ! check if node has particles
    if (inoderange(1,inode) <= 0 .or. inoderange(2,inode) < inoderange(1,inode)) cycle over_nodes

    ! find centre of mass from particle list using same algorithm as maketree
    ! also check for active particles and compute hmax during this loop
    xcofm = 0.
    ycofm = 0.
    zcofm = 0.
    totmass = 0.0
    hmax = 0.
    dfac = 1.
    if (pmassi > 0.) then
       dfac = 1./pmassi
    endif
    nodeisactive = .false.
    do ipart = inoderange(1,inode), inoderange(2,inode)
       if (inodeparts(ipart) == 0) cycle
       xi = treecache(1,ipart)
       yi = treecache(2,ipart)
       zi = treecache(3,ipart)
       hi = treecache(4,ipart)
       ! check condition after loading (dead/accreted particles have hi <= 0)
       if (hi <= 0.) cycle
       hi = abs(hi)
       pmassi = treecache(5,ipart)
       ! check for active particles (for leaf_is_active flag)
       if (ind_timesteps .and. .not. nodeisactive) then
          if (inodeparts(ipart) > 0) nodeisactive = .true.
       endif
       fac = pmassi*dfac
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
       totmass = totmass + pmassi
       hmax = max(hi, hmax)
    enddo
    if (.not. ind_timesteps) nodeisactive = .true.

    if (totmass <= 0.0) cycle over_nodes

    x0(1) = xcofm/(totmass*dfac)
    x0(2) = ycofm/(totmass*dfac)
    x0(3) = zcofm/(totmass*dfac)

    ! update cell size and quads
    r2max = 0.
#ifdef GRAVITY
    quads = 0.
#endif
    do ipart = inoderange(1,inode), inoderange(2,inode)
       ! load all treecache values sequentially (1,2,3,4,5) for cache efficiency
       xi = treecache(1,ipart)
       yi = treecache(2,ipart)
       zi = treecache(3,ipart)
       hi = treecache(4,ipart)
       ! check condition after loading (dead/accreted particles have hi <= 0)
       if (hi <= 0.) cycle
       pmassi = treecache(5,ipart)
       dx = xi - x0(1)
       dy = yi - x0(2)
       dz = zi - x0(3)
       dr2 = dx*dx + dy*dy + dz*dz
       r2max = max(dr2, r2max)
#ifdef GRAVITY
       quads(1) = quads(1) + pmassi*(dx*dx)  ! Q_xx
       quads(2) = quads(2) + pmassi*(dx*dy)  ! Q_xy = Q_yx
       quads(3) = quads(3) + pmassi*(dx*dz)  ! Q_xz = Q_zx
       quads(4) = quads(4) + pmassi*(dy*dy)  ! Q_yy
       quads(5) = quads(5) + pmassi*(dy*dz)  ! Q_yz = Q_zy
       quads(6) = quads(6) + pmassi*(dz*dz)  ! Q_zz
#endif
    enddo

    node(inode)%xcen(1) = x0(1)
    node(inode)%xcen(2) = x0(2)
    node(inode)%xcen(3) = x0(3)
    node(inode)%size = sqrt(r2max) + epsilon(r2max)
    node(inode)%hmax = hmax
#ifdef GRAVITY
    node(inode)%mass = totmass
    node(inode)%quads = quads
#endif

    ! set leaf_is_active flag for leaf nodes (matching maketree behavior)
    if (node(inode)%leftchild == 0 .and. node(inode)%rightchild == 0) then
       if (ind_timesteps) then
          if (nodeisactive) then
             leaf_is_active(inode) = 1
          else
             leaf_is_active(inode) = -1
          endif
       else
          leaf_is_active(inode) = 1
       endif
    endif
 enddo over_nodes
!$omp enddo
!$omp end parallel

end subroutine revtree

!--------------------------------------------------------------------------------
!+
!  Routine to build the global level tree
!+
!-------------------------------------------------------------------------------
subroutine maketreeglobal(nodeglobal,node,nodemap,globallevel,refinelevels,xyzh,&
                          np,cellatid,leaf_is_active,ncells,apr_tree,nptmass,xyzmh_ptmass)
 use io,           only:fatal,warning,id,nprocs,master
 use mpiutils,     only:reduceall_mpi
 use mpibalance,   only:balancedomains
 use mpitree,      only:tree_sync,tree_bcast
 use part,         only:isdead_or_accreted,iactive,ibelong,isink,massoftype,igas,&
                        iamtype,maxphase,maxp,aprmassoftype,apr_level,ihsoft
 use timing,       only:increment_timer,get_timings,itimer_balance
 use dim,          only:ind_timesteps

 type(kdnode),     intent(out)     :: nodeglobal(:)    ! ncellsmax+1
 type(kdnode),     intent(out)     :: node(:)          ! ncellsmax+1
 integer,          intent(out)     :: nodemap(:)       ! ncellsmax+1
 integer,          intent(out)     :: globallevel
 integer,          intent(out)     :: refinelevels
 integer,          intent(inout)   :: np
 real,             intent(inout)   :: xyzh(:,:)
 integer,          intent(out)     :: cellatid(:)      ! ncellsmax+1
 integer,          intent(out)     :: leaf_is_active(:)  ! ncellsmax+1)
 integer(kind=8),  intent(out)     :: ncells
 logical,          intent(in)      :: apr_tree
 integer, optional, intent(in)      :: nptmass
 real,   optional, intent(inout)   :: xyzmh_ptmass(:,:)
 real                              :: xmini(3),xmaxi(3)
 real                              :: xminl(3),xmaxl(3)
 real                              :: xminr(3),xmaxr(3)
 integer                           :: minlevel, maxlevel
 integer                           :: idleft, idright
 integer                           :: groupsize,ifirstingroup,groupsplit
 type(kdnode)                      :: mynode(1)
 integer                           :: nl, nr
 integer                           :: il, ir, iself, parent
 integer                           :: level
 integer                           :: nnodestart, nnodeend,locstart,locend
 integer                           :: npcounter
 integer                           :: i, k, offset, roffset, roffset_prev, coffset
 integer                           :: inode
 integer                           :: npnode
 logical                           :: wassplit,sinktree
 real(kind=4)                      :: t1,t2,tcpu1,tcpu2

 sinktree = .false.
 if (present(nptmass).and.present(xyzmh_ptmass)) sinktree=.true.
 parent = 0
 iself = irootnode
 leaf_is_active = 0

 ! root is level 0
 globallevel = int(ceiling(log(real(nprocs)) / log(2.0)))

 minlevel = maxdepth - 1
 maxlevel = 0

 levels: do level = 0, globallevel
    groupsize = 2**(globallevel - level)
    ifirstingroup = (id / groupsize) * groupsize
    if (level == 0) then
       if (sinktree) then
          call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh,&
                                   xyzmh_ptmass,nptmass)
       else
          call construct_root_node(np,npcounter,irootnode,xmini,xmaxi,leaf_is_active,xyzh)
       endif
    else
       npcounter = npnode
    endif
    if (sinktree) then
       call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
                           il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                           minlevel, maxlevel, wassplit,.true.,apr_tree,xyzmh_ptmass)
    else
       call construct_node(mynode(1), iself, parent, level, xmini, xmaxi, npcounter, .false., &
                        il, ir, nl, nr, xminl, xmaxl, xminr, xmaxr,ncells, leaf_is_active, &
                        minlevel, maxlevel, wassplit,.true.,apr_tree)
    endif

    if (.not.wassplit) then
       call fatal('maketreeglobal','insufficient particles for splitting at the global level: '// &
            'use more particles or less MPI threads')
    endif

    ! set which tree child this proc will belong to next
    groupsplit = ifirstingroup + (groupsize / 2)

    ! record parent for next round
    parent = iself

    ! which half of the tree this task is on
    if (id < groupsplit) then
       ! i for the next node we construct
       iself = il
       ! the left and right task IDs
       idleft = id
       idright = id + 2**(globallevel - level - 1)
       xmini = xminl
       xmaxi = xmaxl
    else
       iself = ir
       idleft = id - 2**(globallevel - level - 1)
       idright = id
       xmini = xminr
       xmaxi = xmaxr
    endif
    if (sinktree) then
       if (nptmass>0) then
          ibelong(maxpsph+1:maxpsph+nptmass) = -1
       endif
    endif
    if (npcounter > 0) then
       do i = inoderange(1,il), inoderange(2,il)
          ibelong(abs(inodeparts(i))) = idleft
       enddo
       do i = inoderange(1,ir), inoderange(2,ir)
          ibelong(abs(inodeparts(i))) = idright
       enddo
    endif

    call get_timings(t1,tcpu1)
    ! move particles to where they belong
    call balancedomains(np)
    call get_timings(t2,tcpu2)
    if (sinktree) ibelong(maxpsph+1:maxpsph+nptmass) = int(reduceall_mpi("max", ibelong(maxpsph+1:maxpsph+nptmass)))
    call increment_timer(itimer_balance,t2-t1,tcpu2-tcpu1)
    ! move particles from old array
    ! this is a waste of time, but maintains compatibility
    npnode = 0
    do i=1,np
       npnode = npnode + 1
       !
       ! tag inactive particles with negative index
       ! in the particle list for the node
       !
       if (ind_timesteps) then
          if (iactive(iphase(i))) then
             inodeparts(npnode) = i
          else
             inodeparts(npnode) = -i
          endif
       else
          inodeparts(npnode) = i
       endif
       treecache(1:4,npnode) = xyzh(1:4,i)
       if (maxphase==maxp) then
          if (use_apr) then
             treecache(5,npnode) = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             treecache(5,npnode) = massoftype(iamtype(iphase(i)))
          endif
       elseif (use_apr) then
          treecache(5,npnode) = aprmassoftype(igas,apr_level(i))
       else
          treecache(5,npnode) = massoftype(igas)
       endif
    enddo
    if (sinktree) then
       if (nptmass > 0) then
          do i=1,nptmass
             if (ibelong(maxpsph + i) /= id) cycle
             if (xyzmh_ptmass(4,i) < 0.) cycle ! dead sink particle
             npnode = npnode + 1
             inodeparts(npnode) = maxpsph + i
             treecache(1:3,npnode) = xyzmh_ptmass(1:3,i)
             treecache(4,npnode)   = xyzmh_ptmass(ihsoft,i)
             treecache(5,npnode)   = xyzmh_ptmass(4,i)
          enddo
       endif
    endif

    ! set all particles to belong to this node
    inoderange(1,iself) = 1
    inoderange(2,iself) = npnode

    ! range of newly written tree
    nnodestart = 2**level
    nnodeend = 2**(level + 1) - 1

    ! synchronize tree with other owners if this proc is the first in group
    call tree_sync(mynode,1,nodeglobal(nnodestart:nnodeend),nprocs/groupsize,ifirstingroup,level)

    ! at level 0, tree_sync already 'broadcasts'
    if (level > 0) then
       ! tree broadcast to non-owners
       call tree_bcast(nodeglobal(nnodestart:nnodeend), nnodeend - nnodestart + 1, level)
    endif

 enddo levels

 ! local tree
 if (sinktree) then
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,refinelevels,nptmass,xyzmh_ptmass)
 else
    call maketree(node,xyzh,np,leaf_is_active,ncells,apr_tree,refinelevels)
 endif

 ! tree refinement
 refinelevels = int(reduceall_mpi('min',refinelevels),kind=kind(refinelevels))
 roffset_prev = 1

 irefine = 0
 do i = 1,refinelevels
    offset = 2**(globallevel + i)
    roffset = 2**i

    nnodestart = offset
    nnodeend   = 2*nnodestart-1

    if (nnodeend > ncellsmax) call fatal('kdtree', 'global tree refinement has exceeded ncellsmax')

    locstart   = roffset
    locend     = 2*locstart-1

    ! index shift the node to the global level
    do k = roffset,2*roffset-1
       refinementnode(k) = node(k)
       coffset = refinementnode(k)%parent - roffset_prev

       refinementnode(k)%parent = 2**(globallevel + i - 1) + id * roffset_prev + coffset

       if (i /= refinelevels) then
          refinementnode(k)%leftchild  = 2**(globallevel + i + 1) + 2*id*roffset + 2*(k - roffset)
          refinementnode(k)%rightchild = refinementnode(k)%leftchild + 1
       else
          refinementnode(k)%leftchild = 0
          refinementnode(k)%rightchild = 0
       endif
    enddo

    roffset_prev = roffset
    ! sync, replacing level with globallevel, since all procs will get synced
    ! and deeper comms do not exist
    call tree_sync(refinementnode(locstart:locend),roffset, &
                   nodeglobal(nnodestart:nnodeend),nnodestart-nnodeend, &
                   id,globallevel)

    ! get the mapping from the local tree to the global tree, for future hmax updates
    do inode = locstart,locend
       nodemap(inode) = nnodestart + (id * roffset) + (inode - locstart)
    enddo
 enddo
!  The index up to which the local tree is copied to the global tree
 irefine = 2*roffset-1

 ! cellatid is zero by default
 cellatid = 0
 do i = 1,nprocs
    offset = 2**(globallevel+refinelevels)
    roffset = 2**refinelevels
    do k = 1,roffset
       cellatid(offset + (i - 1) * roffset + (k - 1)) = i
    enddo
 enddo

end subroutine maketreeglobal

end module kdtree
