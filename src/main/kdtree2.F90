module kdtree
 use dtypekdtree, only:kdnode
 use dim, only:ncellsmax,minpart
 implicit none
 integer, parameter :: ndimtree = 3
 integer, private :: maxlevel_indexed
 integer, private :: maxlevel
 integer, public,  allocatable :: inoderange(:,:)
 integer, public,  allocatable :: iorder(:)
 integer, private, allocatable :: iorder_was(:)
 integer, parameter, public :: irootnode = 1
 real, public :: tree_accuracy = 0.5
 integer, parameter, public :: lenfgrav = 20
 real, allocatable :: treecache(:,:)
 logical, allocatable :: i_am_active(:)

 type kdbuildstack
    integer :: node
    integer :: parent
    integer :: level
    integer :: i1
    integer :: i2
    real    :: xmin(ndimtree)
    real    :: xmax(ndimtree)
 end type

contains

subroutine allocate_kdtree
 use allocutils, only:allocate_array


end subroutine allocate_kdtree

subroutine deallocate_kdtree

end subroutine deallocate_kdtree

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

subroutine maketree(node, xyzh, np, ndim, ifirstincell, ncells)
 use timing, only:getused,get_timings
 use part, only:isdead_or_accreted,igas,maxphase,maxp,massoftype,iphase,iamtype,iactive
 type(kdnode),    intent(out)   :: node(ncellsmax+1)
 integer,         intent(in)    :: np,ndim
 real,            intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,         intent(out)   :: ifirstincell(ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 real :: xmini(3),xmaxi(3)
 integer :: i,inode,n
 integer(1) :: iphasei
 real(4) :: t1,t2,t3,t4,t5,tcpu1,tcpu2,tcpu3,tcpu4,tcpu5

 if (allocated(iorder) .and. size(iorder) < np) then
    deallocate(inoderange,iorder,iorder_was,treecache)
    if (allocated(i_am_active)) deallocate(i_am_active)
 endif
 if (.not.allocated(inoderange)) allocate(inoderange(2,ncellsmax+1))
 if (.not.allocated(iorder)) allocate(iorder(np))
 if (.not.allocated(iorder_was)) allocate(iorder_was(np))
 if (.not.allocated(treecache)) allocate(treecache(5,np))
 if (.not.allocated(i_am_active) .and. maxphase==maxp) allocate(i_am_active(np))

 ! maximum level where 2^k indexing can be used (thus avoiding critical sections)
 ! deeper than this we access cells via a stack as usual
 maxlevel_indexed = int(log(real(ncellsmax+1))/log(2.)) - 1
 maxlevel = 0

 ! default number of cells is the size of the `indexed' part of the tree
 ! this can be *increased* by building tree beyond indexed levels
 ! and is decreased afterwards according to the maximum depth actually reached
 ncells = 2**(maxlevel_indexed+1) - 1

 call get_timings(t1,tcpu1)
! call getused(t1)

 ! wipe tree indices
 call empty_tree(node)
 ifirstincell = 0
 inoderange = 0

 ! get bounds of particle distribution, and enforce periodicity
 call get_particle_bounds(np,xyzh,xmini,xmaxi)
! call getused(t2)
 call get_timings(t2,tcpu2)

 n = 0
 do i=1,np
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       n = n + 1
       iorder(n) = i
    endif
 enddo
 print*,' NDEAD = ',np-n
 call get_timings(t3,tcpu3)


! call getused(t3)
 call build_tree_index1(n,node,xyzh,xmini,xmaxi,ifirstincell,iorder,iorder_was,inoderange,ncells)
! call build_tree_index(n,node,xyzh,xmini,xmaxi,ifirstincell,iorder,iorder_was,inoderange,ncells)
 call get_timings(t4,tcpu4)
 !call getused(t4)

! call build_cache(n,iorder,xyzh,treecache)
 !$omp parallel default(none) &
 !$omp shared(maxphase,maxp,treecache,xyzh,n,massoftype,iphase,iorder) &
 !$omp shared(node,ncells,inoderange,ifirstincell,i_am_active) &
 !$omp private(i,inode,iphasei)
 !$omp do
 do i=1,n
    treecache(1:4,i) = xyzh(:,iorder(i))
    if (maxphase==maxp) then
       iphasei    = iphase(iorder(i))
       treecache(5,i) = massoftype(iamtype(iphasei))
       i_am_active(i) = iactive(iphasei)
    else
       treecache(5,i) = massoftype(igas)
    endif
 enddo
 !$omp end do
 !$omp barrier
 !$omp do schedule(dynamic)
 do inode=1,int(ncells)
    call construct_node(node(inode),inode,inoderange(1,inode),inoderange(2,inode),xyzh,iorder,ifirstincell)
 enddo
 !$omp end do
 !$omp end parallel
! call getused(t5)
 call get_timings(t5,tcpu5)


 print*,' TIMINGS: init: ',t3-t2,t2-t1,' index:',t4-t3,' construct:',t5-t4
 print*,' SCALING: init: ',(tcpu3-tcpu2)/(t3-t2),(tcpu2-tcpu1)/(t2-t1),&
        ' index:',(tcpu4-tcpu3)/(t4-t3),' construct:',(tcpu5-tcpu4)/(t5-t4)

end subroutine maketree

subroutine build_tree_index1(np,node,xyzh,xmin,xmax,ifirstincell,iorder,iorder_was,inoderange,ncells)
 use io, only:fatal,iverbose
!$ use omp_lib
 integer, intent(in) :: np
 type(kdnode), intent(inout) :: node(:)
 real,    intent(in)    :: xyzh(:,:),xmin(ndimtree),xmax(ndimtree)
 integer, intent(inout) :: ifirstincell(:)
 integer, intent(inout) :: iorder(:),iorder_was(:),inoderange(:,:)
 integer(kind=8), intent(inout) :: ncells
 integer :: mymum,level,i1,i2,i1l,i2l,i1r,i2r,il,ir
 integer :: i,istack,istack_global,nnode,ndone,npdone
 integer, save :: numthreads
 logical, save :: done_init_kdtree = .false.
 real :: xmini(ndimtree),xmaxi(ndimtree)
 real :: xminl(ndimtree),xmaxl(ndimtree),xminr(ndimtree),xmaxr(ndimtree)
 logical :: wassplit
 integer, parameter :: istacksize = 512
 type(kdbuildstack), save :: stack(istacksize)
 type(kdbuildstack), save :: stack_global(istacksize)
 !$omp threadprivate(stack)

 istack = 1
 wassplit = .false.

 ! need to number of particles in node during build to allocate space for local link list
 ! this is counted above to remove dead/accreted particles

 if (.not.done_init_kdtree) then
    ! 1 thread for serial, overwritten when using OpenMP
    numthreads = 1

    ! get number of OpenMPthreads
    !$omp parallel default(none) shared(numthreads)
    !$ numthreads = omp_get_num_threads()
    !$omp end parallel
    done_init_kdtree = .true.
 endif

 ! build using a stack which builds depth first
 ! each thread grabs a node from the queue and builds its own subtree
 istack_global = 1
 call push_onto_stack(stack_global(istack_global),irootnode,0,0,1,np,xmin,xmax)

 !$omp parallel default(none) &
 !$omp shared(ifirstincell,inoderange,iorder,iorder_was) &
 !$omp shared(xyzh,stack_global,istack_global) &
 !$omp shared(np) &
 !$omp shared(node, ncells) &
 !$omp shared(numthreads,iverbose) &
 !$omp private(istack) &
 !$omp private(nnode, mymum, level, xmini, xmaxi) &
 !$omp private(il,ir,i1,i2,i1l,i2l,i1r,i2r) &
 !$omp private(xminr, xmaxr, xminl, xmaxl) &
 !$omp private(wassplit,ndone,npdone) &
 !$omp reduction(max:maxlevel)
 !$omp do schedule(static)
 do i = 1, numthreads

    istack = 0
    ndone = 0
    npdone = 0
    mymum = 0

    over_stack: do while(istack > 0 .or. istack_global > 0 .or. ndone < 2**6)

       ! pop off stack
       if (istack > 0) then ! pop off local stack
          call pop_off_stack(stack(istack),istack,nnode,mymum,level,i1,i2,xmini,xmaxi)
       else ! pop off global stack
          mymum = 0; nnode = 0; level = 0; i1=0; i2=0 ! to avoid compiler warnings
          !$omp critical(globalstack)
          if (istack_global > 0) then
             call pop_off_stack(stack_global(istack_global),istack_global,nnode,mymum,level,i1,i2,xmini,xmaxi)
          endif
          !$omp end critical(globalstack)
          if (nnode==0) cycle over_stack ! occurs if istack_global decreased before we had time to pop
       endif
       ndone = ndone + 1
       npdone = npdone + i2-i1+1

       ! construct node
       call build_node(node(nnode),nnode,mymum,i1,i2,&
             level, xmini, xmaxi,il,ir,i1l,i2l,i1r,i2r,xminl, xmaxl, xminr, xmaxr, &
             ncells,ifirstincell(nnode),inoderange(:,nnode),xyzh,iorder,iorder_was,wassplit)

       if (wassplit) then
          if (istack+2 > istacksize) call fatal('maketree',&
                                  'local stack size exceeded in tree build, increase istacksize and recompile')
          if (level < 6 .and. istack_global + 2 < istacksize) then
             !$omp critical(globalstack)
             if (istack_global + 2 < istacksize) then
                istack_global = istack_global + 1
                call push_onto_stack(stack_global(istack_global),il,nnode,level+1,i1l,i2l,xminl,xmaxl)
             else
                istack = istack + 1
                call push_onto_stack(stack(istack),il,nnode,level+1,i1l,i2l,xminl,xmaxl)
             endif
             !$omp end critical(globalstack)
          else
             istack = istack + 1
             call push_onto_stack(stack(istack),il,nnode,level+1,i1l,i2l,xminl,xmaxl)
          endif
          istack = istack + 1
          call push_onto_stack(stack(istack),ir,nnode,level+1,i1r,i2r,xminr,xmaxr)
       endif
    enddo over_stack
    if (iverbose >= 0) print*,' thread ',i,' done ',ndone,' npdone = ',npdone
 enddo
!$omp enddo
!$omp end parallel

end subroutine build_tree_index1

! also used for queue push
pure subroutine push_onto_stack(stackentry,node,parent,level,i1,i2,xmin,xmax)
 type(kdbuildstack), intent(out) :: stackentry
 integer,            intent(in)  :: node,parent,level
 integer,            intent(in)  :: i1,i2
 real,               intent(in)  :: xmin(ndimtree),xmax(ndimtree)

 stackentry%node   = node
 stackentry%parent = parent
 stackentry%level  = level
 stackentry%i1     = i1
 stackentry%i2     = i2
 stackentry%xmin   = xmin
 stackentry%xmax   = xmax

end subroutine push_onto_stack

! also used for queue pop
pure subroutine pop_off_stack(stackentry, istack, nnode, mymum, level, i1, i2, xmini, xmaxi)
 type(kdbuildstack), intent(in)    :: stackentry
 integer,            intent(inout) :: istack
 integer,            intent(out)   :: nnode, mymum, level, i1, i2
 real,               intent(out)   :: xmini(ndimtree), xmaxi(ndimtree)

 nnode  = stackentry%node
 mymum  = stackentry%parent
 level  = stackentry%level
 i1     = stackentry%i1
 i2     = stackentry%i2
 xmini  = stackentry%xmin
 xmaxi  = stackentry%xmax
 istack = istack - 1

end subroutine pop_off_stack

subroutine build_node(nodei,inode,mymum,i1,i2,level,xmini,xmaxi,il,ir,&
                      i1l,i2l,i1r,i2r,xminl,xmaxl,xminr,xmaxr,ncells,&
                      ifirstincelli,inoderangei,xyzh,iorder,iorder_was,was_split)
use io, only:fatal
type(kdnode), intent(out) :: nodei
integer, intent(in) :: inode,mymum,i1,i2,level
real, intent(in)    :: xmini(ndimtree),xmaxi(ndimtree)
integer, intent(out) :: il,ir,i1l,i2l,i1r,i2r
real, intent(out)    :: xminl(ndimtree),xmaxl(ndimtree),xminr(ndimtree),xmaxr(ndimtree)
integer(kind=8), intent(inout) :: ncells
integer, intent(out) :: ifirstincelli,inoderangei(2)
real, intent(in)     :: xyzh(:,:)
integer, intent(inout) :: iorder(:),iorder_was(:)
logical, intent(out) :: was_split
integer :: npnode,nl,nr,iaxis,i,j
real :: xpivot

! set start and end of particle list
if (inode > ncellsmax) call fatal('maketree','number of nodes exceeds array size',var='nodes',ival=inode)
inoderangei(1) = i1  ! this works even if node empty
inoderangei(2) = i2  ! in this case i2 = i1-1 and therefore npnode = 0
npnode = i2-i1+1

!maxlevel = max(level,maxlevel)
if (npnode < minpart) then
   ! we are a leaf node, return
   nodei%parent = mymum
   nodei%leftchild = 0
   nodei%rightchild = 0
   ifirstincelli = 1
   il = 0
   ir = 0
   i1l = 0
   i2l = 0
   i1r = 0
   i2r = 0
   was_split = .false.
   return
else
   ifirstincelli = 0
   was_split = .true. ! has children
endif

!iaxis = maxloc(xmaxi - xmini,1) ! cut longest axis
iaxis = mod(level,3) + 1 ! cycle x->y->z as we descend levels

xpivot = 0.5*(xmini(iaxis) + xmaxi(iaxis))
!print*,' pivot = ',xpivot,' on axis ',iaxis,xmin(iaxis),xmax(iaxis)

xminl = xmini
xmaxl = xmaxi
xminr = xmini
xmaxr = xmaxi
xmaxl(iaxis) = xpivot
xminr(iaxis) = xpivot
!
! sort particles into left or right child nodes
!
!print*,inode,i1,i2
!$omp parallel do if (i2-i1 > 10000) private(j)
do j=i1,i2
   iorder_was(j) = iorder(j)
enddo
!$omp end parallel do
nl = 0
nr = 0
do j=i1,i2
   i = iorder_was(j)
   if (xyzh(iaxis,i) < xpivot) then
      nl = nl + 1
      iorder(i1+nl-1) = i
   else
      nr = nr + 1
      iorder(i2-(nr-1)) = i
   endif
enddo
!print*,' splitting, left  = ',nl,i1,i1+nl-1,' iorder=',iorder(i1:i1+10)
!print*,' splitting, right = ',nr,i1+nl,i2,' iorder=',iorder(i1+nl+1:i1+nl+1+10)
if (level < maxlevel_indexed) then
   il = 2*inode
   ir = il + 1
else
   !$omp critical(ncells)
   ncells = ncells + 1
   il = int(ncells)
   ncells = ncells + 1
   ir = int(ncells)
   !$omp end critical(ncells)
endif
nodei%parent = mymum
nodei%leftchild = il
nodei%rightchild = ir
!
! set particle lists for child nodes
!
i1l = i1       ! left node
i2l = i1+nl-1
i1r = i1+nl    ! right node
i2r = i2

end subroutine build_node

!-----------------------------------
!+
!  get particle bounds
!+
!-----------------------------------
subroutine get_particle_bounds(np,xyzh,xmin,xmax)
#ifdef PERIODIC
 use boundary, only:cross_boundary
 use domain,   only:isperiodic
#endif
 integer, intent(in) :: np
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(out)   :: xmin(3),xmax(3)
 integer :: ncross,i

 xmin(:) = huge(0.)
 xmax(:) = -huge(0.)
 ncross = 0
 !$omp parallel do default(none) &
 !$omp shared(np,xyzh) private(i) &
#ifdef PERIODIC
 !$omp shared(isperiodic) &
#endif
 !$omp reduction(min:xmin) &
 !$omp reduction(max:xmax) &
 !$omp reduction(+:ncross)
 do i=1,np
#ifdef PERIODIC
    call cross_boundary(isperiodic,xyzh(:,i),ncross)
#endif
    xmin(1) = min(xmin(1),xyzh(1,i))
    xmax(1) = max(xmax(1),xyzh(1,i))
    xmin(2) = min(xmin(2),xyzh(2,i))
    xmax(2) = max(xmax(2),xyzh(2,i))
    xmin(3) = min(xmin(3),xyzh(3,i))
    xmax(3) = max(xmax(3),xyzh(3,i))
 enddo
 !$omp end parallel do

end subroutine get_particle_bounds
!-----------------------------------
!+
!  compute centre of mass, size
!  hmax and quadrupole moments
!  for a given node
!+
!-----------------------------------
subroutine construct_node(nodei,inode,i1,i2,xyzh,iorder,ifirstincell)
 use part, only:massoftype,npartoftype,maxtypes,igas,iactive
 type(kdnode), intent(inout) :: nodei
 integer, intent(in) :: inode,i1,i2
 real,    intent(in) :: xyzh(:,:)
 integer, intent(in) :: iorder(:)
 integer, intent(inout) :: ifirstincell(:)
 real :: x0(3),pmassi,fac,dfac
 real :: xi,yi,zi,hi,hmax,dx,dy,dz,dr2,r2max,totmass_node
#ifdef GRAVITY
 real :: quads(6)
#endif
 integer :: j,npnode
 logical :: nodeisactive

 if (i2 < i1 .or. i1==0 .or. i2==0) then
    nodei%xcen = 0.
    nodei%size = 0.
    nodei%hmax = 0.
#ifdef GRAVITY
    nodei%mass  = 0.
    nodei%quads = 0.
#endif
    return
 endif
 npnode = i2 - i1 + 1
!
! tag leaf node as inactive if it contains no active particles
!
 if (npnode < minpart .and. allocated(i_am_active)) then
    nodeisactive = .false.
    checkactive: do j=i1,i2
       if (i_am_active(j)) then
!       if (iactive(iphase(iorder(j)))) then
          nodeisactive = .true.
          exit checkactive
       endif
    enddo checkactive
    if (.not.nodeisactive) ifirstincell(inode) = -1
 endif
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
 !
 ! get hmax and centre of mass
 !
 hmax  = 0.
 x0(:) = 0.
 do j=i1,i2
    xi = treecache(1,j)
    yi = treecache(2,j)
    zi = treecache(3,j)
    hi = treecache(4,j)
    pmassi = treecache(5,j)
    hmax  = max(hmax,hi)
    fac    = pmassi*dfac ! to avoid round-off error
    totmass_node = totmass_node + pmassi
    x0(1) = x0(1) + fac*xi
    x0(2) = x0(2) + fac*yi
    x0(3) = x0(3) + fac*zi
 enddo
 ! if we have no particles, then cofm is zero anyway
 if (totmass_node > 0.) then
    x0(:)   = x0(:)/(totmass_node*dfac)
 endif

 !
 ! get node size and quadrupole moments
 ! relative to centre of mass
 !
 r2max = 0.
#ifdef GRAVITY
 quads(:) = 0.
#endif
!--compute size of node
 do j=i1,i2
    dx = treecache(1,j) - x0(1)
    dy = treecache(2,j) - x0(2)
    dz = treecache(3,j) - x0(3)
    dr2   = dx*dx + dy*dy + dz*dz
    r2max = max(r2max,dr2)
#ifdef GRAVITY
    pmassi = treecache(5,j)
    quads(1) = quads(1) + pmassi*(3.*dx*dx - dr2)  ! Q_xx
    quads(2) = quads(2) + pmassi*(3.*dx*dy)        ! Q_xy = Q_yx
    quads(3) = quads(3) + pmassi*(3.*dx*dz)        ! Q_xz = Q_zx
    quads(4) = quads(4) + pmassi*(3.*dy*dy - dr2)  ! Q_yy
    quads(5) = quads(5) + pmassi*(3.*dy*dz)        ! Q_yz = Q_zy
    quads(6) = quads(6) + pmassi*(3.*dz*dz - dr2)  ! Q_zz
#endif
 enddo

 ! assign properties to node
 nodei%xcen    = x0(:)
 nodei%size    = sqrt(r2max) + epsilon(r2max)
 nodei%hmax    = hmax
#ifdef GRAVITY
 nodei%mass    = totmass_node
 nodei%quads   = quads
#endif

end subroutine construct_node

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
subroutine getneigh(node,xpos,xsizei,rcuti,ndim,listneigh,nneigh,xyzh,xyzcache,ixyzcachesize,&
& ifirstincell,get_hj,fnode,remote_export)
 use dim,      only:maxneigh
#ifdef PERIODIC
 use boundary, only:dxbound,dybound,dzbound
#endif
 use io,       only:fatal,id
 use part,     only:gravity
#ifdef FINVSQRT
 use fastmath, only:finvsqrt
#endif
 use kernel,   only:radkern
 type(kdnode), intent(in)           :: node(ncellsmax+1)
 integer, intent(in)                :: ndim,ixyzcachesize
 real,    intent(in)                :: xpos(ndim)
 real,    intent(in)                :: xsizei,rcuti
 integer, intent(out)               :: listneigh(maxneigh)
 integer, intent(out)               :: nneigh
 real,    intent(in)                :: xyzh(:,:)
 real,    intent(out)               :: xyzcache(:,:)
 integer, intent(in)                :: ifirstincell(:)
 logical, intent(in)                :: get_hj
 real,    intent(out),    optional  :: fnode(lenfgrav)
 logical, intent(out),    optional  :: remote_export(:)
 integer, parameter :: istacksize = 300
 integer :: maxcache
 integer :: nstack(istacksize)
 integer :: n,istack,il,ir,npnode,i1,i2,j,jmax,nn
 real :: dx,dy,dz,xsizej,rcutj
 real :: rcut,rcut2,r2,hmax
 real :: xoffset,yoffset,zoffset,tree_acc2
 logical :: open_tree_node
#ifdef GRAVITY
 real :: quads(6)
 real :: dr,totmass_node
#endif
#ifdef PERIODIC
 real :: hdlx,hdly,hdlz

 hdlx = 0.5*dxbound
 hdly = 0.5*dybound
 hdlz = 0.5*dzbound
#endif
 tree_acc2 = tree_accuracy*tree_accuracy
 if (present(fnode)) fnode(:) = 0.
 rcut     = rcuti

 if (ixyzcachesize > 0) then
    maxcache = size(xyzcache(:,1))
 else
    maxcache = 0
 endif

 if (present(remote_export)) remote_export = .false.

 nneigh = 0
 istack = 1
 nstack(istack) = irootnode
 open_tree_node = .false.

 over_stack: do while(istack /= 0)
    n = nstack(istack)
    istack = istack - 1
    dx = xpos(1) - node(n)%xcen(1)      ! distance between node centres
    dy = xpos(2) - node(n)%xcen(2)
#ifndef TREEVIZ
    dz = xpos(3) - node(n)%xcen(3)
#endif
    xsizej       = node(n)%size
    hmax         = node(n)%hmax
    il           = node(n)%leftchild
    ir           = node(n)%rightchild
#ifdef GRAVITY
    totmass_node = node(n)%mass
    quads        = node(n)%quads
#endif
    !call get_child_nodes(n,il,ir)
    xoffset = 0.
    yoffset = 0.
    zoffset = 0.
#ifdef PERIODIC
    if (abs(dx) > hdlx) then            ! mod distances across boundary if periodic BCs
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
    r2    = dx*dx + dy*dy + dz*dz
    if (get_hj) then  ! find neighbours within both hi and hj
       rcutj = radkern*hmax
       rcut  = max(rcuti,rcutj)
    endif
    rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
#ifdef GRAVITY
    open_tree_node = tree_acc2*r2 < xsizej*xsizej    ! tree opening criterion for self-gravity
#endif
    if ((r2 < rcut2) .or. open_tree_node) then
       if_leaf: if (abs(ifirstincell(n)) > 0) then ! once we hit a leaf node, retrieve contents into trial neighbour cache
          if_global_walk: if (present(remote_export)) then
             ! id is stored in ipart as id + 1
             if (ifirstincell(n) /= (id + 1)) then
                remote_export(ifirstincell(n)) = .true.
             endif
          else
             i1=inoderange(1,n)
             i2=inoderange(2,n)
             npnode = i2 - i1 + 1
             do j=1,npnode
                listneigh(nneigh+j) = iorder(i1+j-1)
             enddo
             if (nneigh < ixyzcachesize) then ! not strictly necessary, loop doesn't execute anyway
                nn = min(nneigh+npnode,ixyzcachesize)
                !print*,' got ',nneigh+npnode,' capping to ',nn
                jmax = max(nn-nneigh,0)
                !print*,'cache build ',nneigh+1,'->',nneigh+jmax,'parts ',i1,'->',i1+jmax-1,&
                  !  ' nneigh = ',nneigh+1,'->',nneigh+npnode
                if (maxcache >= 4) then
                   do j=1,jmax
                      xyzcache(1,nneigh+j) = treecache(1,i1+j-1) + xoffset
                      xyzcache(2,nneigh+j) = treecache(2,i1+j-1) + yoffset
                      xyzcache(3,nneigh+j) = treecache(3,i1+j-1) + zoffset
                      xyzcache(4,nneigh+j) = 1./treecache(4,i1+j-1)
                   enddo
                else
                   do j=1,jmax
                      xyzcache(1,nneigh+j) = treecache(1,i1+j-1) + xoffset
                      xyzcache(2,nneigh+j) = treecache(2,i1+j-1) + yoffset
                      xyzcache(3,nneigh+j) = treecache(3,i1+j-1) + zoffset
                   enddo
                endif
             endif
             nneigh = nneigh + npnode
          endif if_global_walk
       else
          if (istack+2 > istacksize) call fatal('getneigh','stack overflow in getneigh')
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
    elseif (present(fnode) .and. ((.not. present(remote_export)) .or. n < 2*nprocs-1)) then
!
!--long range force on node due to distant node, along node centres
!  along with derivatives in order to perform series expansion
!
#ifdef FINVSQRT
       dr = finvsqrt(r2)
#else
       dr = 1./sqrt(r2)
#endif
       call compute_fnode(dx,dy,dz,dr,totmass_node,quads,fnode)
#endif
    endif
 enddo over_stack

end subroutine getneigh


!-----------------------------------------------------------
!+
!  Compute the gravitational force between the node centres
!  along with the derivatives and second derivatives
!  required for the Taylor series expansions.
!+
!-----------------------------------------------------------
pure subroutine compute_fnode(dx,dy,dz,dr,totmass,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr,totmass
 real, intent(in)    :: quads(6)
 real, intent(inout) :: fnode(lenfgrav)
 real :: dr3,dr4,dr5,dr6,dr3m,rx,ry,rz,qxx,qxy,qxz,qyy,qyz,qzz
 real :: dr4m3,rijQij,riQix,riQiy,riQiz,fqx,fqy,fqz
 real :: dfxdxq,dfxdyq,dfxdzq,dfydyq,dfydzq,dfzdzq
 real :: d2fxxxq,d2fxxyq,d2fxxzq,d2fxyyq,d2fxyzq
 real :: d2fxzzq,d2fyyyq,d2fyyzq,d2fyzzq,d2fzzzq

 ! note: dr == 1/sqrt(r2)
 dr3  = dr*dr*dr
 dr4  = dr*dr3
 dr5  = dr*dr4
 dr6  = dr*dr5
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
 rijQij = (rx*rx*qxx + ry*ry*qyy + rz*rz*qzz + 2.*(rx*ry*qxy + rx*rz*qxz + ry*rz*qyz))
 riQix = (rx*qxx + ry*qxy + rz*qxz)
 riQiy = (rx*qxy + ry*qyy + rz*qyz)
 riQiz = (rx*qxz + ry*qyz + rz*qzz)
 fqx = dr4*(riQix - 2.5*rx*rijQij)
 fqy = dr4*(riQiy - 2.5*ry*rijQij)
 fqz = dr4*(riQiz - 2.5*rz*rijQij)
 dfxdxq = dr5*(qxx - 10.*rx*riQix - 2.5*rijQij   + 17.5*rx*rx*rijQij)
 dfxdyq = dr5*(qxy -  5.*ry*riQix - 5.0*rx*riQiy + 17.5*rx*ry*rijQij)
 dfxdzq = dr5*(qxz -  5.*rx*riQiz - 5.0*rz*riQix + 17.5*rx*rz*rijQij)
 dfydyq = dr5*(qyy - 10.*ry*riQiy - 2.5*rijQij   + 17.5*ry*ry*rijQij)
 dfydzq = dr5*(qyz -  5.*ry*riQiz - 5.0*rz*riQiy + 17.5*ry*rz*rijQij)
 dfzdzq = dr5*(qzz - 10.*rz*riQiz - 2.5*rijQij   + 17.5*rz*rz*rijQij)
 d2fxxxq = dr6*(-15.*qxx*rx + 105.*rx*rx*riQix - 15.*riQix - 157.5*rx*rx*rx*rijQij + 52.5*rx*rijQij)
 d2fxxyq = dr6*(35.*rx*rx*riQiy -  5.*qxx*ry - 5.*riQiy + 17.5*ry*rijQij - 157.5*rx*rx*ry*rijQij &
              + 70.*rx*ry*riQix - 10.*qxy*rx)
 d2fxxzq = dr6*(35.*rx*rx*riQiz -  5.*qxx*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*rx*rx*rz*rijQij &
              + 70.*rx*rz*riQix - 10.*qxz*rx)
 d2fxyyq = dr6*(70.*rx*ry*riQiy - 10.*qxy*ry - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*ry*ry*rijQij &
              + 35.*ry*ry*riQix -  5.*qyy*rx)
 d2fxyzq = dr6*(35.*rx*ry*riQiz -  5.*qyz*rx  &
              + 35.*ry*rz*riQix -  5.*qxz*ry  &
              + 35.*rx*rz*riQiy -  5.*qxy*rz                             - 157.5*rx*ry*rz*rijQij)
 d2fxzzq = dr6*(70.*rx*rz*riQiz - 10.*qxz*rz - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*rz*rz*rijQij &
              + 35.*rz*rz*riQix -  5.*qzz*rx)
 d2fyyyq = dr6*(-15.*qyy*ry + 105.*ry*ry*riQiy - 15.*riQiy - 157.5*ry*ry*ry*rijQij + 52.5*ry*rijQij)
 d2fyyzq = dr6*(35.*ry*ry*riQiz -  5.*qyy*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*ry*ry*rz*rijQij &
              + 70.*ry*rz*riQiy - 10.*qyz*ry)
 d2fyzzq = dr6*(70.*ry*rz*riQiz - 10.*qyz*rz - 5.*riQiy + 17.5*ry*rijQij - 157.5*ry*rz*rz*rijQij &
              + 35.*rz*rz*riQiy -  5.*qzz*ry)
 d2fzzzq = dr6*(-15.*qzz*rz + 105.*rz*rz*riQiz - 15.*riQiz - 157.5*rz*rz*rz*rijQij + 52.5*rz*rijQij)

 fnode( 1) = fnode( 1) - dx*dr3m + fqx ! fx
 fnode( 2) = fnode( 2) - dy*dr3m + fqy ! fy
 fnode( 3) = fnode( 3) - dz*dr3m + fqz ! fz
 fnode( 4) = fnode( 4) + dr3m*(3.*rx*rx - 1.) + dfxdxq ! dfx/dx
 fnode( 5) = fnode( 5) + dr3m*(3.*rx*ry)      + dfxdyq ! dfx/dy = dfy/dx
 fnode( 6) = fnode( 6) + dr3m*(3.*rx*rz)      + dfxdzq ! dfx/dz = dfz/dx
 fnode( 7) = fnode( 7) + dr3m*(3.*ry*ry - 1.) + dfydyq ! dfy/dy
 fnode( 8) = fnode( 8) + dr3m*(3.*ry*rz)      + dfydzq ! dfy/dz = dfz/dy
 fnode( 9) = fnode( 9) + dr3m*(3.*rz*rz - 1.) + dfzdzq ! dfz/dz
 fnode(10) = fnode(10) - dr4m3*(5.*rx*rx*rx - 3.*rx) + d2fxxxq ! d2fxdxdx
 fnode(11) = fnode(11) - dr4m3*(5.*rx*rx*ry - ry)    + d2fxxyq ! d2fxdxdy
 fnode(12) = fnode(12) - dr4m3*(5.*rx*rx*rz - rz)    + d2fxxzq ! d2fxdxdz
 fnode(13) = fnode(13) - dr4m3*(5.*rx*ry*ry - rx)    + d2fxyyq ! d2fxdydy
 fnode(14) = fnode(14) - dr4m3*(5.*rx*ry*rz)         + d2fxyzq ! d2fxdydz
 fnode(15) = fnode(15) - dr4m3*(5.*rx*rz*rz - rx)    + d2fxzzq ! d2fxdzdz
 fnode(16) = fnode(16) - dr4m3*(5.*ry*ry*ry - 3.*ry) + d2fyyyq ! d2fydydy
 fnode(17) = fnode(17) - dr4m3*(5.*ry*ry*rz - rz)    + d2fyyzq ! d2fydydz
 fnode(18) = fnode(18) - dr4m3*(5.*ry*rz*rz - ry)    + d2fyzzq ! d2fydzdz
 fnode(19) = fnode(19) - dr4m3*(5.*rz*rz*rz - 3.*rz) + d2fzzzq ! d2fzdzdz
 fnode(20) = fnode(20) - totmass*dr - 0.5*rijQij*dr3   ! potential

end subroutine compute_fnode

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
 poti = poti - (dx*fxi + dy*fyi + dz*fzi)

 return
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
subroutine revtree(node, xyzh, ifirstincell, ncells)
 use dim,  only:maxp
 use part, only:maxphase,iphase,igas,massoftype,iamtype
 use io,   only:fatal
 type(kdnode), intent(inout) :: node(ncellsmax+1)
 real,    intent(in)  :: xyzh(:,:)
 integer, intent(in)  :: ifirstincell(ncellsmax+1)
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
!$omp shared(xyzh, ifirstincell, ncells) &
!$omp shared(node, iphase, massoftype, maxlevel) &
!$omp private(hmax, r2max, xi, yi, zi, hi, il, ir, nodel, noder) &
!$omp private(dx, dy, dz, dr2, icell, i, x0) &
#ifdef GRAVITY
!$omp private(quads) &
#endif
!$omp firstprivate(pmassi) &
!$omp private(totmass)
!$omp do schedule(guided, 2)
 over_cells: do icell=1,int(ncells)

    i = abs(ifirstincell(icell))
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
          pmassi = massoftype(iamtype(iphase(i)))
       endif
       x0(1) = x0(1) + pmassi*xi
       x0(2) = x0(2) + pmassi*yi
       x0(3) = x0(3) + pmassi*zi
       totmass = totmass + pmassi

       !i = abs(ll(i))
    enddo calc_cofm

    x0 = x0/totmass
#ifdef GRAVITY
    node(icell)%xcen(1) = x0(1)
    node(icell)%xcen(2) = x0(2)
    node(icell)%xcen(3) = x0(3)
#endif

    i = abs(ifirstincell(icell))

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
          pmassi = massoftype(iamtype(iphase(i)))
       endif
       quads(1) = quads(1) + pmassi*(3.*dx*dx - dr2)
       quads(2) = quads(2) + pmassi*(3.*dx*dy)
       quads(3) = quads(3) + pmassi*(3.*dx*dz)
       quads(4) = quads(4) + pmassi*(3.*dy*dy - dr2)
       quads(5) = quads(5) + pmassi*(3.*dy*dz)
       quads(6) = quads(6) + pmassi*(3.*dz*dz - dr2)
#endif
       ! move to next particle in list
       !i = abs(ll(i))
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
 real :: ql(6),qr(6)
#endif
 real :: dx,dy,dz,dr2,dr
 real :: mr,ml,mnode,dm

 xl = l%xcen
 hl = l%hmax
 sl = l%size
#ifdef GRAVITY
 ml = l%mass
 ql = l%quads
#else
 ml = 1.
#endif

 xr = r%xcen
 hr = r%hmax
 sr = r%size
#ifdef GRAVITY
 mr = r%mass
 qr = r%quads
#else
 mr = 1.
#endif
 mnode = ml + mr
 if (mnode > 0.) then
    dm = 1./mnode
 else
    dm = 0.
 endif
 dx = xl(1) - xr(1)
 dy = xl(2) - xr(2)
 dz = xl(3) - xr(3)
 dr2 = dx*dx + dy*dy + dz*dz
 dr  = sqrt(dr2)
 ! centre of mass
 nodei%xcen = (xl*ml + xr*mr)*dm
 ! size, formula as in Benz et al. 1990
 ! and from thinking about it...
 nodei%size = max(ml*dm*dr+sr,mr*dm*dr+sl)
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

end module kdtree
