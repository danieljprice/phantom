!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
 use part,        only:ll,iphase,xyzh_soa,iphase_soa,maxphase, &
                       apr_level,apr_level_soa,aprmassoftype

 implicit none

 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: inodeparts(:)
 type(kdnode),     allocatable :: refinementnode(:)
 integer,          allocatable :: nstack(:)
 integer,          allocatable :: queue_dual(:)
!$omp threadprivate(queue_dual,nstack)
!
!--tree parameters
!
 integer,          parameter, public :: irootnode = 1
 character(len=1), parameter, public :: labelax(3) = (/'x','y','z'/)
 integer,          parameter         :: maxlevelcrazy = 31
!
!--runtime options for this module
!
 real,    public  :: tree_accuracy = 0.55
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
 call allocate_array('queue_dual', queue_dual, ncellsmax+1)
 call allocate_array('nstack', nstack, ncellsmax+1)
!$omp end parallel

end subroutine allocate_kdtree

subroutine deallocate_kdtree
 use dim, only:mpi
 if (allocated(inoderange)) deallocate(inoderange)
 if (allocated(inodeparts)) deallocate(inodeparts)
 if (mpi .and. allocated(refinementnode)) deallocate(refinementnode)
!$omp parallel
 if (allocated(queue_dual)) deallocate(queue_dual)
 if (allocated(nstack)) deallocate(nstack)
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
 minlevel = 31
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

    ! get number of OpenMPthreads
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
 use part, only:isink
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
 !$omp shared(inodeparts,iphase,xyzh_soa,iphase_soa,nproot,apr_level_soa) &
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
       xyzh_soa(nproot,:) = xyzh(:,i)
       iphase_soa(nproot) = iphase(i)
       if (use_apr) apr_level_soa(nproot) = apr_level(i)
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
          xyzh_soa(nproot,:) = xyzmh_ptmass(1:4,i)
          iphase_soa(nproot) = isink
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
 use part,      only:massoftype,igas,iamtype,maxphase,maxp,npartoftype,isink,ihsoft
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
 integer :: i,npcounter,i1
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
    !$omp shared(xyzh_soa,apr_level_soa,i1,iphase_soa) &
    !$omp shared(xyzmh_ptmass,sinktree) &
    !$omp private(i,xi,yi,zi,hi) &
    !$omp firstprivate(pmassi,fac) &
    !$omp reduction(+:xcofm,ycofm,zcofm,totmass_node) &
    !$omp reduction(max:hmax)
    do i=i1,i1+npnode-1
       xi = xyzh_soa(i,1)
       yi = xyzh_soa(i,2)
       zi = xyzh_soa(i,3)
       hi = xyzh_soa(i,4)
       if (maxphase==maxp) then
          if (sinktree .and. (iamtype(iphase_soa(i)) == isink)) then
             hi = xyzmh_ptmass(ihsoft,inodeparts(i)-maxpsph)
             pmassi = xyzh_soa(i,4)
          elseif (use_apr) then
             pmassi = aprmassoftype(iamtype(iphase_soa(i)),apr_level_soa(i))
          else
             pmassi = massoftype(iamtype(iphase_soa(i)))
          endif
          fac    = pmassi*dfac ! to avoid round-off error
       elseif (use_apr) then
          pmassi = aprmassoftype(igas,apr_level_soa(i))
          fac    = pmassi*dfac ! to avoid round-off error
       endif
       hmax  = max(hmax,hi)
       totmass_node = totmass_node + pmassi
       xcofm = xcofm + fac*xi
       ycofm = ycofm + fac*yi
       zcofm = zcofm + fac*zi
    enddo
    !$omp end parallel do
 else
    do i=i1,i1+npnode-1
       xi = xyzh_soa(i,1)
       yi = xyzh_soa(i,2)
       zi = xyzh_soa(i,3)
       hi = xyzh_soa(i,4)
       if (maxphase==maxp) then
          if (sinktree .and. (iamtype(iphase_soa(i)) == isink)) then
             hi = xyzmh_ptmass(ihsoft,inodeparts(i)-maxpsph)
             pmassi = xyzh_soa(i,4)
          elseif (use_apr) then
             pmassi = aprmassoftype(iamtype(iphase_soa(i)),apr_level_soa(i))
          else
             pmassi = massoftype(iamtype(iphase_soa(i)))
          endif
          fac    = pmassi*dfac ! to avoid round-off error
       elseif (use_apr) then
          pmassi = aprmassoftype(igas,apr_level_soa(i))
          fac    = pmassi*dfac ! to avoid round-off error
       endif
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
 !!$omp parallel do if (npnode > 1000 .and. doparallel) &
 !!$omp default(none) schedule(static) &
 !!$omp shared(npnode,xyzh_soa,x0,i1,apr_level_soa) &
 !!$omp shared(iphase_soa,massoftype,sinktree) &
 !!$omp private(i,xi,yi,zi,dx,dy,dz,dr2,pmassi) &
#ifdef GRAVITY
 !!$omp reduction(+:quads) &
#endif
 !!$omp reduction(max:r2max)
 do i=i1,i1+npnode-1
    xi = xyzh_soa(i,1)
    yi = xyzh_soa(i,2)
    zi = xyzh_soa(i,3)
    dx    = xi - x0(1)
    dy    = yi - x0(2)
    dz    = zi - x0(3)
    dr2   = dx*dx + dy*dy + dz*dz
    r2max = max(r2max,dr2)
#ifdef GRAVITY
    if (maxphase==maxp) then
       if (use_apr) then
          pmassi = aprmassoftype(iamtype(iphase_soa(i)),apr_level_soa(i))
       elseif (sinktree .and. (iamtype(iphase_soa(i)) == isink)) then
          pmassi = xyzh_soa(i,4)
       else
          pmassi = massoftype(iamtype(iphase_soa(i)))
       endif
    endif
    quads(1) = quads(1) + pmassi*(dx*dx)  ! Q_xx
    quads(2) = quads(2) + pmassi*(dx*dy)        ! Q_xy = Q_yx
    quads(3) = quads(3) + pmassi*(dx*dz)        ! Q_xz = Q_zx
    quads(4) = quads(4) + pmassi*(dy*dy)  ! Q_yy
    quads(5) = quads(5) + pmassi*(dy*dz)        ! Q_yz = Q_zy
    quads(6) = quads(6) + pmassi*(dz*dz)  ! Q_zz
#endif
 enddo
 !!$omp end parallel do

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
                                    inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts,&
                                    npnode,apr_level_soa)
       else
          ! regular sort
          call sort_particles_in_cell(iaxis,inoderange(1,nnode),inoderange(2,nnode),inoderange(1,il),inoderange(2,il),&
                                  inoderange(1,ir),inoderange(2,ir),nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts,apr_level_soa)
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

       xminl(1) = minval(xyzh_soa(inoderange(1,il):inoderange(2,il),1))
       xminl(2) = minval(xyzh_soa(inoderange(1,il):inoderange(2,il),2))
       xminl(3) = minval(xyzh_soa(inoderange(1,il):inoderange(2,il),3))

       xmaxl(1) = maxval(xyzh_soa(inoderange(1,il):inoderange(2,il),1))
       xmaxl(2) = maxval(xyzh_soa(inoderange(1,il):inoderange(2,il),2))
       xmaxl(3) = maxval(xyzh_soa(inoderange(1,il):inoderange(2,il),3))

       xminr(1) = minval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),1))
       xminr(2) = minval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),2))
       xminr(3) = minval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),3))

       xmaxr(1) = maxval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),1))
       xmaxr(2) = maxval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),2))
       xmaxr(3) = maxval(xyzh_soa(inoderange(1,ir):inoderange(2,ir),3))
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
subroutine sort_particles_in_cell(iaxis,imin,imax,min_l,max_l,min_r,max_r,nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts,apr_level_soa)
 integer, intent(in)  :: iaxis,imin,imax
 integer, intent(out) :: min_l,max_l,min_r,max_r,nl,nr
 real, intent(inout)  :: xpivot,xyzh_soa(:,:)
 integer(kind=1), intent(inout) :: iphase_soa(:),apr_level_soa(:)
 integer,         intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot
 integer(kind=1) :: iphase_swap,apr_swap
 integer :: inodeparts_swap,i,j
 real :: xyzh_swap(4)

 !print*,'nnode ',imin,imax,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
 j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
 !  k = 0

 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:4)  = xyzh_soa(i,1:4)
          iphase_swap     = iphase_soa(i)
          if (use_apr) apr_swap = apr_level_soa(i)

          inodeparts(i)   = inodeparts(j)
          xyzh_soa(i,1:4) = xyzh_soa(j,1:4)
          iphase_soa(i)   = iphase_soa(j)
          if (use_apr) apr_level_soa(i)= apr_level_soa(j)

          inodeparts(j)   = inodeparts_swap
          xyzh_soa(j,1:4) = xyzh_swap(1:4)
          iphase_soa(j)   = iphase_swap
          if (use_apr) apr_level_soa(j)= apr_swap

          i = i + 1
          j = j - 1
          i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
          j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
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
                                nl,nr,xpivot,xyzh_soa,iphase_soa,inodeparts,npnode,apr_level_soa)
 use io, only:error
 integer, intent(in)  :: iaxis,imin,imax,npnode
 integer, intent(out) :: min_l,max_l,min_r,max_r,nl,nr
 real, intent(inout)  :: xpivot,xyzh_soa(:,:)
 integer(kind=1), intent(inout) :: iphase_soa(:),apr_level_soa(:)
 integer,         intent(inout) :: inodeparts(:)
 logical :: i_lt_pivot,j_lt_pivot,slide_l,slide_r
 integer(kind=1) :: iphase_swap,apr_swap
 integer :: inodeparts_swap,i,j,nchild_in
 integer :: k,ii,rem_nr,rem_nl
 real :: xyzh_swap(4),dpivot(npnode)

 dpivot = 0.0
 nchild_in = 2

 if (modulo(npnode,nchild_in) > 0) then
    call error('apr sort','number of particles sent in to kdtree is not divisible by 2')
 endif

! print*,'nnode ',imin,imax,npnode,' pivot = ',iaxis,xpivot
 i = imin
 j = imax

 i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
 j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
 dpivot(i-imin+1) = xpivot - xyzh_soa(i,iaxis)
 dpivot(j-imin+1) = xpivot - xyzh_soa(j,iaxis)
 !k = 0
 do while(i < j)
    if (i_lt_pivot) then
       i = i + 1
       dpivot(i-imin+1) = xpivot - xyzh_soa(i,iaxis)
       i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
    else
       if (.not.j_lt_pivot) then
          j = j - 1
          dpivot(j-imin+1) = xpivot - xyzh_soa(j,iaxis)
          j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
       else
          ! swap i and j positions in list
          inodeparts_swap = inodeparts(i)
          xyzh_swap(1:4)  = xyzh_soa(i,1:4)
          iphase_swap     = iphase_soa(i)
          apr_swap        = apr_level_soa(i)

          inodeparts(i)   = inodeparts(j)
          xyzh_soa(i,1:4) = xyzh_soa(j,1:4)
          iphase_soa(i)   = iphase_soa(j)
          apr_level_soa(i)= apr_level_soa(j)

          inodeparts(j)   = inodeparts_swap
          xyzh_soa(j,1:4) = xyzh_swap(1:4)
          iphase_soa(j)   = iphase_swap
          apr_level_soa(j)= apr_swap

          i = i + 1
          j = j - 1

          dpivot(i-imin+1) = xpivot - xyzh_soa(i,iaxis)
          dpivot(j-imin+1) = xpivot - xyzh_soa(j,iaxis)

          i_lt_pivot = xyzh_soa(i,iaxis) <= xpivot
          j_lt_pivot = xyzh_soa(j,iaxis) <= xpivot
       endif
    endif
 enddo

 if (.not.i_lt_pivot) then
    i = i - 1
    dpivot(i-imin+1) = xpivot - xyzh_soa(i,iaxis)
 endif
 if (j_lt_pivot) then
    j = j + 1
    dpivot(j-imin+1) = xpivot - xyzh_soa(j,iaxis)
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
       xyzh_swap(1:4)  = xyzh_soa(k,1:4)
       iphase_swap     = iphase_soa(k)

       inodeparts(k)   = inodeparts(j)
       xyzh_soa(k,1:4) = xyzh_soa(j,1:4)
       iphase_soa(k)   = iphase_soa(j)

       inodeparts(j)   = inodeparts_swap
       xyzh_soa(j,1:4) = xyzh_swap(1:4)
       iphase_soa(j)   = iphase_swap

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
       xyzh_swap(1:4)  = xyzh_soa(k,1:4)
       iphase_swap     = iphase_soa(k)

       inodeparts(k)   = inodeparts(i)
       xyzh_soa(k,1:4) = xyzh_soa(i,1:4)
       iphase_soa(k)   = iphase_soa(i)

       inodeparts(i)   = inodeparts_swap
       xyzh_soa(i,1:4) = xyzh_swap(1:4)
       iphase_soa(i)   = iphase_swap

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
 integer :: ipart,n,istack,il,ir,npnode
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
    maxcache = size(xyzcache(1,:))
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

 if (present(nq)) then
    nstack(1:nq) = queue_dual(1:nq)
    istack = nq
 endif

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
             npnode = inoderange(2,n) - inoderange(1,n) + 1
             if_cache_fits: if (nneigh + npnode <= ixyzcachesize) then
                if (maxcache >= 4) then
                   do ipart=1,npnode
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                      xyzcache(nneigh+ipart,4) = 1./xyzh_soa(inoderange(1,n)+ipart-1,4)
                   enddo
                else
                   do ipart=1,npnode
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                   enddo
                endif
             elseif (nneigh < ixyzcachesize) then
                if (maxcache >= 4) then
                   do ipart=1,ixyzcachesize-nneigh
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                      xyzcache(nneigh+ipart,4) = 1./xyzh_soa(inoderange(1,n)+ipart-1,4)
                   enddo
                else
                   do ipart=1,ixyzcachesize-nneigh
                      listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                      xyzcache(nneigh+ipart,1) = xyzh_soa(inoderange(1,n)+ipart-1,1) + xoffset
                      xyzcache(nneigh+ipart,2) = xyzh_soa(inoderange(1,n)+ipart-1,2) + yoffset
                      xyzcache(nneigh+ipart,3) = xyzh_soa(inoderange(1,n)+ipart-1,3) + zoffset
                   enddo
                endif
                do ipart=ixyzcachesize-nneigh+1,npnode
                   listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                enddo
             else
                do ipart=1,npnode
                   listneigh(nneigh+ipart) = abs(inodeparts(inoderange(1,n)+ipart-1))
                enddo
             endif if_cache_fits
             nneigh = nneigh + npnode
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
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
subroutine getneigh_dual(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,ifirstincell,&
                         get_hj,get_f,fnode,remote_export,icell)
 use io,       only:fatal
 type(kdnode), intent(in)           :: node(:) !ncellsmax+1)
 integer, intent(in)                :: ixyzcachesize
 real,    intent(in)                :: xpos(3)
 real,    intent(in)                :: xsizei,rcuti
 integer, intent(out)               :: listneigh(:)
 integer, intent(out)               :: nneigh
 real,    intent(out)               :: xyzcache(:,:)
 integer, intent(in)                :: ifirstincell(:)
 logical, intent(in)                :: get_hj
 logical, intent(in)                :: get_f
 real,    intent(out),    optional  :: fnode(lenfgrav)
 integer, intent(in),     optional  :: icell
 logical, intent(out),    optional  :: remote_export(:)
 integer :: nq
 integer :: istack,inext,i
 integer :: parents(maxlevel+1),nparents,inode
 real    :: fnode_old(lenfgrav),dx,dy,dz,xdum,ydum,zdum
 real    :: tree_acc2
 logical :: open_tree_node

 tree_acc2 = tree_accuracy*tree_accuracy
 if (get_f .and. (.not.present(fnode) .or. .not.present(icell))) then
    call fatal('getneigh_dual','get_f but fnode not passed...')
 endif
 if (present(fnode)) fnode(:) = 0.

 call get_list_of_parent_nodes(icell,node,parents,nparents)

 nneigh = 0
 nq = 1
 istack = 0
 queue_dual(nq) = irootnode
 open_tree_node = .false.
 do i=nparents,2,-1 ! parents(1) is equal to icell
    inode = parents(i)
    call check_node_interactions(node(inode),inode,node,fnode,tree_acc2,queue_dual,nq,nstack,istack)
    fnode_old = fnode
    inext = parents(i-1)
    call get_sep(node(inext)%xcen,node(inode)%xcen,dx,dy,dz,xdum,ydum,zdum)
    call propagate_fnode_to_node(fnode,fnode_old,dx,dy,dz)
    nq = istack
    queue_dual(1:nq) = nstack(1:nq)
    istack = 0
 enddo

 call getneigh(node,xpos,xsizei,rcuti,listneigh,nneigh,xyzcache,ixyzcachesize,ifirstincell,&
               get_hj,get_f,fnode_old,remote_export=remote_export,nq=nq)

 fnode = fnode + fnode_old

end subroutine getneigh_dual

!-----------------------------------------------------------
!+
!  get the separation in 3D between two nodes of the tree
!+
!-----------------------------------------------------------
subroutine get_sep(x1,x2,dx,dy,dz,xoffset,yoffset,zoffset,r2)
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
subroutine get_node_size(node_dst,node_src,size_dst,size_src,rcut_dst,rcut_src)
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
subroutine propagate_fnode_to_node(fnode,fnode_old,dx,dy,dz)
 real, intent(in)    :: fnode_old(lenfgrav),dx,dy,dz
 real, intent(inout) :: fnode(lenfgrav)

 fnode(1)  = fnode_old(1) + dx*(fnode_old(4) + 0.5*(dx*fnode_old(10) + dy*fnode_old(11) +dz*fnode_old(12)))& ! xx +0.5(xxx+xxy+xxz)
                          + dy*(fnode_old(5) + 0.5*(dx*fnode_old(11) + dy*fnode_old(13) +dz*fnode_old(14)))& ! xy +0.5(xxy+xyy+xyz)
                          + dz*(fnode_old(6) + 0.5*(dx*fnode_old(12) + dy*fnode_old(14) +dz*fnode_old(15)))  ! xz +0.5(xxz+xyz+xzz)
 fnode(2)  = fnode_old(2) + dx*(fnode_old(5) + 0.5*(dx*fnode_old(11) + dy*fnode_old(13) +dz*fnode_old(14)))& ! xy +0.5(xxy+xyy+xyz)
                          + dy*(fnode_old(7) + 0.5*(dx*fnode_old(13) + dy*fnode_old(16) +dz*fnode_old(17)))& ! yy +0.5(xyy+yyy+yyz)
                          + dz*(fnode_old(8) + 0.5*(dx*fnode_old(14) + dy*fnode_old(17) +dz*fnode_old(18)))  ! yz +0.5(xyz+yyz+yyz)
 fnode(3)  = fnode_old(3) + dx*(fnode_old(6) + 0.5*(dx*fnode_old(12) + dy*fnode_old(14) +dz*fnode_old(15)))& ! xz +0.5(xxz+xyz+xzz)
                          + dy*(fnode_old(8) + 0.5*(dx*fnode_old(14) + dy*fnode_old(17) +dz*fnode_old(18)))& ! yz +0.5(xyz+yyz+yzz)
                          + dz*(fnode_old(9) + 0.5*(dx*fnode_old(15) + dy*fnode_old(18) +dz*fnode_old(19)))  ! zz +0.5(xzz+yzz+zzz)
 fnode(4)  = fnode_old(4) + dx*fnode_old(10) + dy*fnode_old(11) + dz*fnode_old(12)                           ! xxx + xxy + xxz
 fnode(5)  = fnode_old(5) + dx*fnode_old(11) + dy*fnode_old(13) + dz*fnode_old(14)                           ! xxy + xyy + xyz
 fnode(6)  = fnode_old(6) + dx*fnode_old(12) + dy*fnode_old(14) + dz*fnode_old(15)                           ! xxz + xyz + xzz
 fnode(7)  = fnode_old(7) + dx*fnode_old(13) + dy*fnode_old(16) + dz*fnode_old(17)                           ! xyy + yyy + yyz
 fnode(8)  = fnode_old(8) + dx*fnode_old(14) + dy*fnode_old(17) + dz*fnode_old(18)                           ! xyz + yyz + yzz
 fnode(9)  = fnode_old(9) + dx*fnode_old(15) + dy*fnode_old(18) + dz*fnode_old(19)                           ! xzz + yzz + zzz
 fnode(10) = fnode_old(10)
 fnode(11) = fnode_old(11)
 fnode(12) = fnode_old(12)
 fnode(13) = fnode_old(13)
 fnode(14) = fnode_old(14)
 fnode(15) = fnode_old(15)
 fnode(16) = fnode_old(16)
 fnode(17) = fnode_old(17)
 fnode(18) = fnode_old(18)
 fnode(19) = fnode_old(19)
 fnode(20) = fnode_old(20) + dx*(fnode_old(1)+0.5*(dx*fnode_old(4)+dy*fnode_old(5)+dz*fnode_old(6)))&
                           + dy*(fnode_old(2)+0.5*(dx*fnode_old(5)+dy*fnode_old(7)+dz*fnode_old(8)))&
                           + dz*(fnode_old(3)+0.5*(dx*fnode_old(6)+dy*fnode_old(8)+dz*fnode_old(9)))

end subroutine propagate_fnode_to_node

!-----------------------------------------------------------
!+
!  return list of parents of current node
!+
!-----------------------------------------------------------
subroutine get_list_of_parent_nodes(inode,node,parents,nparents)
 integer, intent(in) :: inode
 type(kdnode), intent(in) :: node(:)
 integer, intent(out) :: parents(:)
 integer, intent(out) :: nparents
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
!  Check the node node separation, compute the interaction
!  if well-separated or stack children if it's a
!  self interaction or not well separated.
!+
!-----------------------------------------------------------
subroutine check_node_interactions(node_dst,inode,node,fnode,tree_acc2,queue,nq,nstack,istack)
 use io, only:fatal
 type(kdnode), intent(in)    :: node(:),node_dst
 real,         intent(in)    :: tree_acc2
 real,         intent(inout) :: fnode(lenfgrav)
 integer,      intent(in)    :: inode
 integer,      intent(inout) :: queue(:),nq
 integer,      intent(out)   :: nstack(:),istack
 integer :: i
 logical :: stackit

 istack = 0

 stack_internode: do while(nq /= 0)
    i = queue(nq)
    if (i == inode) then
       stackit = .true.
    else
       call node_interaction(node_dst,node(i),tree_acc2,fnode,stackit)
    endif

    if (stackit) then
       if (istack+2 > ncellsmax+1) call fatal('check_node_interactions','stack overflow in node interactions')
       call open_nodes(nstack,istack,node(i),i)
    endif
    nq = nq - 1

 enddo stack_internode

end subroutine check_node_interactions

!-----------------------------------------------------------
!+
!  Compute node node gravity interactions
!+
!-----------------------------------------------------------
subroutine open_nodes(nstack,istack,node,inode)
 type(kdnode), intent(in)     :: node
 integer,      intent(in)     :: inode
 integer,      intent(inout)  :: nstack(:),istack
 integer :: ir,il

 il = node%leftchild
 ir = node%rightchild

 if (il /= 0) then
    istack = istack + 1
    nstack(istack) = il
 endif
 if (ir /= 0) then
    istack = istack + 1
    nstack(istack) = ir
 else ! if ir or il == 0 then you're a leaf, should be stack
    istack = istack + 1
    nstack(istack) = inode
 endif

end subroutine open_nodes

!-----------------------------------------------------------
!+
!  Test the separation between the node pair and compute
!  the interaction if needed
!+
!-----------------------------------------------------------
subroutine node_interaction(node_dst,node_src,tree_acc2,fnode,stackit)
 type(kdnode), intent(in)    :: node_dst,node_src
 real,         intent(in)    :: tree_acc2
 real,         intent(inout) :: fnode(lenfgrav)
 logical,      intent(out)   :: stackit
 real    :: dx,dy,dz,r2,dr1
 real    :: rcut_dst,rcut_src,rcut,rcut2
 real    :: size_dst,size_src,mass_src,quads_src(6)
 real    :: xdum,ydum,zdum
 logical :: wellsep

 call get_sep(node_dst%xcen,node_src%xcen,dx,dy,dz,xdum,ydum,zdum,r2)
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
pure subroutine compute_M2L(dx,dy,dz,dr,totmass,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr,totmass
 real, intent(in)    :: quads(6)
 real, intent(inout) :: fnode(lenfgrav)
 real :: dr2,dr3,dr4,dr5,dr3m,dr4m3,rx,ry,rz,qxx,qxy,qxz,qyy,qyz,qzz
 real :: fqx,fqy,fqz,rijQij,Qii

! note: dr == 1/sqrt(r2)
 dr2  = dr*dr
 dr3  = dr2*dr
 dr4  = dr2*dr2
 dr5  = dr4*dr
 dr3m  = totmass*dr3
 dr4m3 = 3.*totmass*dr4
 rx  = dx*dr
 ry  = dy*dr
 rz  = dz*dr
 qxx = quads(1)
 qxy = quads(2)
 qxz = quads(3)
 qyy = quads(4)
 qyz = quads(5)
 qzz = quads(6)
 rijQij = (rx*rx*qxx + 2.*ry*rx*qxy + 2*ry*rz*qyz + ry*ry*qyy + 2*rx*rz*qxz + rz*rz*qzz)
 Qii    = (qxx + qyy + qzz)
 fqx    = dr4*(1.5*(rx*(3*qxx+qyy+qzz) + 2.*ry*qxy + 2.*rz*qxz) - 7.5*rx*rijQij)
 fqy    = dr4*(1.5*(ry*(3*qyy+qxx+qzz) + 2.*rx*qxy + 2.*rz*qyz) - 7.5*ry*rijQij)
 fqz    = dr4*(1.5*(rz*(3*qzz+qyy+qxx) + 2.*ry*qyz + 2.*rx*qxz) - 7.5*rz*rijQij)

 fnode(1)  = fnode(1)  - dr3m*dx + fqx                          ! C¹_x
 fnode(2)  = fnode(2)  - dr3m*dy + fqy                          ! C¹_y
 fnode(3)  = fnode(3)  - dr3m*dz + fqz                          ! C¹_z
 fnode(4)  = fnode(4)  + dr3m*(3.*rx*rx -   1.)                 ! C²_xx
 fnode(5)  = fnode(5)  + dr3m*(3.*rx*ry       )                 ! C²_xy
 fnode(6)  = fnode(6)  + dr3m*(3.*rx*rz       )                 ! C²_xz
 fnode(7)  = fnode(7)  + dr3m*(3.*ry*ry -   1.)                 ! C²_yy
 fnode(8)  = fnode(8)  + dr3m*(3.*ry*rz       )                 ! C²_yz
 fnode(9)  = fnode(9)  + dr3m*(3.*rz*rz -   1.)                 ! C²_zz
 fnode(10) = fnode(10) - dr4m3*(5.*rx*rx*rx - 3.*rx)            ! C³_xxx
 fnode(11) = fnode(11) - dr4m3*(5.*rx*rx*ry - ry)               ! C³_xxy
 fnode(12) = fnode(12) - dr4m3*(5.*rx*rx*rz - rz)               ! C³_xxz
 fnode(13) = fnode(13) - dr4m3*(5.*rx*ry*ry - rx)               ! C³_xyy
 fnode(14) = fnode(14) - dr4m3*(5.*rx*ry*rz)                    ! C³_xyz
 fnode(15) = fnode(15) - dr4m3*(5.*rx*rz*rz - rx)               ! C³_xzz
 fnode(16) = fnode(16) - dr4m3*(5.*ry*ry*ry - 3.*ry)            ! C³_yyy
 fnode(17) = fnode(17) - dr4m3*(5.*ry*ry*rz - rz)               ! C³_yyz
 fnode(18) = fnode(18) - dr4m3*(5.*ry*rz*rz - ry)               ! C³_yzz
 fnode(19) = fnode(19) - dr4m3*(5.*rz*rz*rz - 3.*rz)            ! C³_zzz
 fnode(20) = fnode(20) - (totmass*dr + 0.5*dr3*(3.*rijQij-Qii)) ! C⁰ (potential)

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
 use dim,  only:maxp
 use part, only:maxphase,iphase,igas,massoftype,iamtype
 use io,   only:fatal
 type(kdnode), intent(inout) :: node(:) !ncellsmax+1)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: leaf_is_active(:) !ncellsmax+1)
 integer(kind=8), intent(in) :: ncells
 real :: hmax, r2max
 real :: xi, yi, zi, hi
 real :: dx, dy, dz, dr2
#ifdef GRAVITY
 real :: quads(6)
#endif
 integer :: icell, i, level, il, ir
 real :: pmassi, totmass
 real :: x0(3)
 type(kdnode) :: nodel,noder

 pmassi = massoftype(igas)
 do i=1,int(ncells)
#ifdef GRAVITY
    ! cannot update centre of node without gravity
    ! as it is not at the centre of mass
    node(i)%xcen(:) = 0.
#endif
    node(i)%size    = 0.
    node(i)%hmax    = 0.
#ifdef GRAVITY
    node(i)%mass    = 0.
    node(i)%quads(:)= 0.
#endif
 enddo

!$omp parallel default(none) &
!$omp shared(maxp,maxphase) &
!$omp shared(xyzh, leaf_is_active, ncells, apr_level) &
!$omp shared(node, ll, iphase, massoftype, maxlevel,aprmassoftype) &
!$omp private(hmax, r2max, xi, yi, zi, hi, il, ir, nodel, noder) &
!$omp private(dx, dy, dz, dr2, icell, i, x0) &
#ifdef GRAVITY
!$omp private(quads) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)
!$omp do schedule(guided, 2)
 over_cells: do icell=1,int(ncells)

    i = abs(leaf_is_active(icell))
    if (i==0) cycle over_cells

    ! find centre of mass
    ! this becomes the new node center
    x0 = 0.
    totmass = 0.0
    calc_cofm: do while (i /= 0)
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       if (maxphase==maxp) then
          if (use_apr) then
             pmassi = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             pmassi = massoftype(iamtype(iphase(i)))
          endif
       endif
       x0(1) = x0(1) + pmassi*xi
       x0(2) = x0(2) + pmassi*yi
       x0(3) = x0(3) + pmassi*zi
       totmass = totmass + pmassi

       i = abs(ll(i))
    enddo calc_cofm

    x0 = x0/totmass
#ifdef GRAVITY
    node(icell)%xcen(1) = x0(1)
    node(icell)%xcen(2) = x0(2)
    node(icell)%xcen(3) = x0(3)
#endif

    i = abs(leaf_is_active(icell))

    ! update cell size, hmax
    r2max = 0.
    hmax = 0.
#ifdef GRAVITY
    quads = 0.
#endif
    over_parts: do while (i /= 0)
       xi = xyzh(1,i)
       yi = xyzh(2,i)
       zi = xyzh(3,i)
       hi = xyzh(4,i)
       dx = xi - node(icell)%xcen(1)
       dy = yi - node(icell)%xcen(2)
       dz = zi - node(icell)%xcen(3)
       dr2 = dx*dx + dy*dy + dz*dz
       r2max = max(dr2, r2max)
       hmax  = max(hi, hmax)
#ifdef GRAVITY
       if (maxphase==maxp) then
          if (use_apr) then
             pmassi = aprmassoftype(iamtype(iphase(i)),apr_level(i))
          else
             pmassi = massoftype(iamtype(iphase(i)))
          endif
       endif
       quads(1) = quads(1) + pmassi*(3.*dx*dx - dr2)
       quads(2) = quads(2) + pmassi*(3.*dx*dy)
       quads(3) = quads(3) + pmassi*(3.*dx*dz)
       quads(4) = quads(4) + pmassi*(3.*dy*dy - dr2)
       quads(5) = quads(5) + pmassi*(3.*dy*dz)
       quads(6) = quads(6) + pmassi*(3.*dz*dz - dr2)
#endif
       ! move to next particle in list
       i = abs(ll(i))
    enddo over_parts

    node(icell)%size = sqrt(r2max) + epsilon(r2max)
    node(icell)%hmax = hmax
#ifdef GRAVITY
    node(icell)%mass = totmass
    node(icell)%quads = quads
#endif
 enddo over_cells
!$omp enddo
!
! propagate information to parent nodes
! here we sweep across each level at a time
! and update each node from its two children
!
 do level=maxlevel-1,0,-1
!$omp do
    do i=2**level,2**(level+1)-1
       ! get child nodes
       il = node(i)%leftchild
       ir = node(i)%rightchild
       if (il > 0 .and. ir > 0) then
          nodel = node(il)
          noder = node(ir)
          call add_child_nodes(nodel,noder,node(i))
       else
          if (il > 0 .or. ir > 0) then
             ! should never happen, should have two children or none
             call fatal('revtree','node with only one child during tree revision',var='ir',ival=ir)
          endif
       endif
    enddo
!$omp enddo
 enddo
!$omp end parallel

end subroutine revtree

!-----------------------------------------------------------------
!+
!  Update parent node from the properties of the two child nodes
!  IN:
!    l, r - two child nodes
!  OUT:
!    nodei - updated parent node
!+
!-----------------------------------------------------------------
subroutine add_child_nodes(l,r,nodei)
 type(kdnode), intent(in)  :: l,r
 type(kdnode), intent(out) :: nodei
 real :: xl(3),sl,hl
 real :: xr(3),sr,hr
#ifdef GRAVITY
 real :: ql(6),qr(6),mr,ml,mnode,dm
#endif
 real :: dx,dy,dz,dr2,dr

 xl = l%xcen
 hl = l%hmax
 sl = l%size
#ifdef GRAVITY
 ml = l%mass
 ql = l%quads
#endif

 xr = r%xcen
 hr = r%hmax
 sr = r%size
#ifdef GRAVITY
 mr = r%mass
 qr = r%quads
 mnode = ml + mr
 dm    = 1./mnode
#endif
 dx = xl(1) - xr(1)
 dy = xl(2) - xr(2)
 dz = xl(3) - xr(3)
 dr2 = dx*dx + dy*dy + dz*dz
 dr  = sqrt(dr2)
#ifdef GRAVITY
 ! centre of mass
 nodei%xcen = (xl*ml + xr*mr)*dm
 ! size, formula as in Benz et al. 1990
 ! and from thinking about it...
 nodei%size = max(ml*dm*dr+sr,mr*dm*dr+sl)
#else
 ! distance between left child and node centre
 dx = xl(1) - nodei%xcen(1)
 dy = xl(2) - nodei%xcen(2)
 dz = xl(3) - nodei%xcen(3)
 dr = sqrt(dx*dx + dy*dy + dz*dz)
 nodei%size = dr+sl
 ! distance between right child and node centre
 dx = xr(1) - nodei%xcen(1)
 dy = xr(2) - nodei%xcen(2)
 dz = xr(3) - nodei%xcen(3)
 dr = sqrt(dx*dx + dy*dy + dz*dz)
 nodei%size = max(nodei%size,dr+sr)
#endif
 nodei%hmax = max(hl,hr)
#ifdef GRAVITY
 nodei%mass = mnode
 ! quadrupole moments, see Benz et al. (1990), this is also
 ! the parallel axis theorem
 nodei%quads(1) = ql(1) + qr(1) + ml*mr*(3.*dx*dx - dr2)*dm
 nodei%quads(2) = ql(2) + qr(2) + ml*mr*(3.*dx*dy)*dm
 nodei%quads(3) = ql(3) + qr(3) + ml*mr*(3.*dx*dz)*dm
 nodei%quads(4) = ql(4) + qr(4) + ml*mr*(3.*dy*dy - dr2)*dm
 nodei%quads(5) = ql(5) + qr(5) + ml*mr*(3.*dy*dz)*dm
 nodei%quads(6) = ql(6) + qr(6) + ml*mr*(3.*dz*dz - dr2)*dm
#endif

end subroutine add_child_nodes

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
 use part,         only:isdead_or_accreted,iactive,ibelong,isink
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

 minlevel = 31
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
       xyzh_soa(npnode,:) = xyzh(:,i)
       iphase_soa(npnode) = iphase(i)
       if (use_apr) apr_level_soa(npnode) = apr_level(i)
    enddo
    if (sinktree) then
       if (nptmass > 0) then
          do i=1,nptmass
             if (ibelong(maxpsph + i) /= id) cycle
             if (xyzmh_ptmass(4,i) < 0.) cycle ! dead sink particle
             npnode = npnode + 1
             inodeparts(npnode) = maxpsph + i
             xyzh_soa(npnode,:) = xyzmh_ptmass(1:4,i)
             iphase_soa(npnode) = isink
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
