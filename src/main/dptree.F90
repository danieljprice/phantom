module dptree
 use dim, only:ncellsmax,minpart
 use dtypekdtree, only:kdnode
 implicit none
 integer, parameter :: ndimtree = 3
 real, public :: tree_accuracy = 0.5

! type kdnode
!   sequence
!   real :: xcen(ndimtree)
!   real :: size
!   real :: hmax
!#ifdef GRAVITY
!   real :: mass
!   real :: quads(6)
!#endif
! end type kdnode

 public :: maketree1,getneigh1,climb_tree_hmax
 public :: get_node_list
 private

 integer, public :: maxlevel
 integer, allocatable, public :: iorder(:)
 integer, allocatable, public :: inoderange(:,:)
 real, allocatable :: pkey(:)
 integer, parameter :: irootnode = 1

contains

!-------------------------------------
!+
!  Build the tree, from the ground up
!+
!-------------------------------------
subroutine maketree1(node, xyzh, np, ndim, ifirstincell, ncells)
 use kdtree, only:empty_tree
 type(kdnode),    intent(out)   :: node(ncellsmax+1)
 integer,         intent(in)    :: np,ndim
 real,            intent(inout) :: xyzh(:,:)  ! inout because of boundary crossing
 integer,         intent(out)   :: ifirstincell(ncellsmax+1)
 integer(kind=8), intent(out)   :: ncells
 integer :: leaf_level,i
 real :: xmin(3),xmax(3)

 call empty_tree(node)
 ! get sort key
 if (.not.allocated(pkey)) allocate(pkey(np))

 ! get bounds of particle distribution, and enforce periodicity
 call get_particle_bounds(np,xyzh,xmin,xmax)

 ! get Hilbert keys for all particles
 call get_sort_key(np,xyzh,xmin,xmax,pkey)

 if (allocated(iorder) .and. size(iorder) /= np) deallocate(iorder)
 if (.not.allocated(iorder)) allocate(iorder(np))
 if (.not.allocated(inoderange)) allocate(inoderange(2,ncellsmax+1))

 ! sort particles into order based on the sort key
 call parallel_sort(np,pkey,iorder)

 ! divide this list into equal leaf nodes (nearest power of 2)
 ! and get particle lists for all parent nodes
 call build_tree_index1(np,minpart,leaf_level,inoderange,ncells,pkey)

 ! get properties of all leaf nodes (in parallel)
 call build_leaf_nodes(leaf_level,node,xyzh,iorder)

 ! propagate information to parent nodes
 call climb_tree(node,inoderange,leaf_level,xyzh,iorder)

 maxlevel = leaf_level
 return
     write(*,"(a,i10,3(a,i2))") ' maketree1: nodes = ',ncells,', leaf level = ',maxlevel
     block
        real :: sizesum,sizemax,pkeysize
        integer :: i,j,istart,iend,nleaf,level,maxpnode,minpnode,npnode
        print*,' size root node = ',node(1)%size
        do level=leaf_level,0,-1
           call get_node_list(level,istart,iend)
           sizesum = 0.
           sizemax = 0.
           do i=istart,iend
              sizesum = sizesum + node(i)%size
              sizemax = max(sizemax,node(i)%size)
           enddo
           nleaf = iend - istart + 1
           if (level==leaf_level) print*,' mean parts per leaf node = ',np/real(nleaf)
           print*,level,' SCORE = ',node(1)%size*real(nleaf)/sizesum,' mean size = ',&
                        sizesum/real(nleaf),'max=',sizemax,'sum=',sizesum,' nnodes=',nleaf
        enddo
        call get_node_list(leaf_level,istart,iend)
        minpnode = huge(maxpnode)
        maxpnode = 0
        do i=istart,iend
           npnode = inoderange(2,i) - inoderange(1,i) + 1
           minpnode = min(minpnode,npnode)
           maxpnode = max(maxpnode,npnode)
        enddo
        print*,' MAX particles per node = ',maxpnode,' MIN=',minpnode

 open(1,file='pkey.list')
 write(1,"(a)") '# i  nodeid  id  x  y  z  h  pkey  nodesize   pkeysize'
 call get_node_list(leaf_level,istart,iend)
 do j=istart,iend
    pkeysize = pkey(iorder(inoderange(2,j))) - pkey(iorder(inoderange(1,j)))
    do i=inoderange(1,j),inoderange(2,j)
       write(1,*) i,j,iorder(i),xyzh(:,iorder(i)),pkey(iorder(i)),node(j)%size,pkeysize
    enddo
 enddo
 close(1)
 stop

     end block


end subroutine maketree1

!-----------------------------------
!+
!  Build the tree index, i.e. index
!  of first and last particle on
!  every node
!+
!-----------------------------------
subroutine build_tree_index(np,parts_per_leaf,leaf_level,inoderange,ncells,pkey)
 use io, only:fatal
 integer, intent(in)  :: np,parts_per_leaf
 integer, intent(out) :: leaf_level
 integer, intent(out) :: inoderange(:,:)
 real,    intent(in)  :: pkey(:)
 integer(kind=8), intent(out) :: ncells
 real :: ratio,parts_per_leafi,pkeyrangemean,pkeyrangemax,pkeyrange
 integer :: nleaf,inode,index,istart,iend
 integer :: level,il,ir,i1,i2,j

 ratio           = np/real(parts_per_leaf)
 leaf_level      = int(log(ratio)/log(2.)) + 1
 nleaf           = 2**leaf_level
 parts_per_leafi = np/real(nleaf)

 !
 ! divide sorted list into approximately equal length chunks
 !
 call get_node_list(leaf_level,istart,iend)
 ncells = int(iend)
 if (ncells > size(inoderange(1,:))) call fatal('maketree','array size too small for number of levels')
 pkeyrangemean = 0.
 pkeyrangemax = 0.
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(istart,iend,parts_per_leafi,inoderange,pkey,iorder) &
 !$omp private(inode,index,pkeyrange) reduction(max:pkeyrangemax) reduction(+:pkeyrangemean)
 do inode=istart,iend
    index = inode - istart
    inoderange(1,inode) = int(index*parts_per_leafi) + 1
    inoderange(2,inode) = int(((index+1)*parts_per_leafi))
    pkeyrange = pkey(iorder(inoderange(2,inode))) - pkey(iorder(inoderange(1,inode)))
    pkeyrangemean = pkeyrangemean + pkeyrange
    pkeyrangemax = max(pkeyrangemax,pkeyrange)
 enddo
 !$omp end parallel do

 pkeyrangemean = pkeyrangemean/(iend-istart+1)
 print*,' MEAN PKEY RANGE = ',pkeyrangemean,' MAX = ',pkeyrangemax
 do j=1,0
 do inode=istart,iend-1
    i2 = inoderange(2,inode)
    i1 = inoderange(1,inode+1)
    pkeyrange = pkey(iorder(inoderange(2,inode))) - pkey(iorder(inoderange(1,inode)))
    ! move last particle to next bin if it is closer
    ! in key distance to first particle of next leaf
    if (inoderange(2,inode) > inoderange(1,inode) + 1) then
       if (pkeyrange > 10.*pkeyrangemean) then
          !print*,'inode',inode,' parts ',i2,'->',i1,' got ',pkeyrange,' ratio =',pkeyrange/pkeyrangemean
          inoderange(2,inode) = inoderange(2,inode) - 1
          inoderange(1,inode+1) = inoderange(1,inode+1) - 1
      !    print*,'shifting boundary, parts=',inoderange(2,inode) - inoderange(1,inode)
       endif
    endif
 enddo
 !do inode=iend,istart+1
!    i1 = inoderange(1,inode)
!    i2 = inoderange(2,inode-1)
    ! move last particle to next bin if it is closer
    ! in key distance to first particle of next leaf
    !print*,'inode ',i2,i1,' got '
!    if (inoderange(2,inode) > inoderange(1,inode)) then
!       if ((pkey(i1)-pkey(i2)) < (pkey(i1+1) - pkey(i1))) then
!          inoderange(2,inode-1) = inoderange(2,inode-1) + 1
!          inoderange(1,inode) = inoderange(1,inode) + 1
!       endif
!    endif
 !enddo
 enddo
 !
 ! sweep up the levels, building index ranges of parent nodes
 ! based on index ranges of the lower level chunks
 !
 !$omp parallel default(none) &
 !$omp shared(inoderange,leaf_level) &
 !$omp private(level,inode,istart,iend,il,ir)
 do level=leaf_level-1,0,-1
    call get_node_list(level,istart,iend)
    !print*,' *** level ',level,' nodes ',istart,':',iend,' ***'
    !$omp do schedule(static)
    do inode=istart,iend
       call get_child_nodes(inode,il,ir)
       inoderange(1,inode) = inoderange(1,il)
       inoderange(2,inode) = inoderange(2,ir)
       !if (level < 2) print*,'node ',inode,': ',inoderange(1,inode),'->',inoderange(2,inode)
    enddo
    !$omp end do
 enddo
 !$omp end parallel

 ! sanity checks
 if (inoderange(1,irootnode) /= 1) call fatal('maketree','starting index of root node /= 1')
 if (inoderange(2,irootnode) /= np) call fatal('maketree','finishing index of root node /= np')

end subroutine build_tree_index

!-----------------------------------
!+
!  Build the tree index, i.e. index
!  of first and last particle on
!  every node
!+
!-----------------------------------
subroutine build_tree_index1(np,parts_per_leaf,leaf_level,inoderange,ncells,pkey)
 use io, only:fatal
 integer, intent(in)  :: np,parts_per_leaf
 integer, intent(out) :: leaf_level
 integer, intent(out) :: inoderange(:,:)
 real,    intent(in)  :: pkey(:)
 integer(kind=8), intent(out) :: ncells
 real :: ratio,parts_per_leafi,pkeyrangemean,pkeyrangemax,pkeyrange
 integer :: nleaf,inode,index,istart,iend
 integer :: level,il,ir,i1,i2,j,ipivot

 ratio           = np/real(parts_per_leaf)
 leaf_level      = int(log(ratio)/log(2.)) + 1
 nleaf           = 2**leaf_level
 parts_per_leafi = np/real(nleaf)
 !
 ! start with all particles in root node
 !
 inoderange(1,irootnode) = 1
 inoderange(2,irootnode) = np
 !
 ! recursively partition into
 !
 !$omp parallel default(none) &
 !$omp shared(inoderange,leaf_level,iorder,pkey,parts_per_leaf) &
 !$omp private(level,inode,istart,iend,il,ir,i1,i2,ipivot)
 do level=0,leaf_level-1
    call get_node_list(level,istart,iend)
    !$omp do
    do inode=istart,iend
       call get_child_nodes(inode,il,ir)
       i1 = inoderange(1,inode)
       i2 = inoderange(2,inode)
       if (i2 > i1) then
          ipivot = get_median(level,parts_per_leaf,inoderange(1,inode),inoderange(2,inode),iorder,pkey)
          inoderange(1,il) = inoderange(1,inode)
          inoderange(2,il) = ipivot
          inoderange(1,ir) = ipivot+1
          inoderange(2,ir) = inoderange(2,inode)
       else
          inoderange(1,il) = inoderange(1,inode)
          inoderange(2,il) = inoderange(1,inode)-1
          inoderange(1,ir) = inoderange(1,inode)
          inoderange(2,ir) = inoderange(1,inode)-1
       endif
    enddo
    !$omp end do
 enddo
 !$omp end parallel

end subroutine build_tree_index1

integer function get_median(level,parts_per_leaf,i1,i2,iorder,x)
 integer, intent(in) :: level,parts_per_leaf,i1,i2
 integer, intent(in) :: iorder(:)
 real,    intent(in) :: x(:)
 integer :: n,i,minpart_level,mylevel
 real :: xmed,ratio,parts_per_leafi

 !get_median = i1
 !return
 !get_median = (i1 - 1) + (i2 - i1)/2
 !return
 if (i2 <= i1) then
    get_median = i1
    return
 endif

 ratio = (i2-i1)/real(parts_per_leaf)
 mylevel = int(log(ratio)/log(2.)) + 1
 minpart_level = 2**(mylevel+2)
 !parts_per_leafi = (i2-i1)/real(minpart_level)
 !print*,' ALLOWED ',minpart_level, ' ON LEVEL ',level,' with ',parts_per_leafi,' parts per leaf'
 n = 0
 xmed = 0.
 do i=i1,i2
    n = n + 1
    xmed = xmed + x(iorder(i))
 enddo
 xmed = xmed/real(n)
 !print*,' got median ',xmed
 i = i1
 do while (x(iorder(i)) < xmed .and. i < i2)
    i = i + 1
 enddo
 get_median = i
 !return
 !print*,' WANT PIVOT ',i,' MINPART ON LEVEL ',level,' = ',minpart_level
 if (i > i1 + minpart_level) then
    i = i1 + minpart_level
 elseif (i < i2 - minpart_level) then
    i = i2 - minpart_level
 endif
 !print*,' GOT PIVOT ',i
 get_median = i

 return

end function get_median

!-----------------------------------
!+
!  get index of first and last node
!  on a given tree level
!+
!-----------------------------------
pure subroutine get_node_list(level,istart,iend)
 integer, intent(in) :: level
 integer, intent(out) :: istart,iend

 istart = 2**level
 iend   = 2**(level+1) - 1

end subroutine get_node_list

!-----------------------------------
!+
!  get index of two child nodes
!  for a given parent node
!+
!-----------------------------------
pure subroutine get_child_nodes(inode,il,ir)
 integer, intent(in) :: inode
 integer, intent(out) :: il,ir

 il = 2*inode   ! indexing as per Gafton & Rosswog (2011)
 ir = il + 1

end subroutine get_child_nodes

!-----------------------------------
!+
!  get properties of leaf nodes
!  i.e. centre of mass, size, hmax
!  and quadrupole moments
!+
!-----------------------------------
subroutine build_leaf_nodes(leaf_level,node,xyzh,iorder)
 integer,      intent(in)  :: leaf_level
 type(kdnode), intent(out) :: node(:)
 real,         intent(in)  :: xyzh(:,:)
 integer,      intent(in)  :: iorder(:)
 integer :: inode,istart,iend

 call get_node_list(leaf_level,istart,iend)
 !$omp parallel do default(none) &
 !$omp shared(istart,iend,node,inoderange,iorder,xyzh) &
 !$omp private(inode)
 do inode=istart,iend
    call construct_node(node(inode),inoderange(1,inode),inoderange(2,inode),xyzh,iorder)
 enddo
 !$omp end parallel do

end subroutine build_leaf_nodes

!---------------------------------------------
!+
!  propagate information from the leaf nodes
!  back up to the parent nodes
!+
!---------------------------------------------
subroutine climb_tree(node,inoderange,leaf_level,xyzh,iorder)
 use kdtree, only:add_child_nodes
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: inoderange(:,:),leaf_level
 real,         intent(in)    :: xyzh(:,:)
 integer,      intent(in)    :: iorder(:)
 integer :: level,inode,il,ir,istart,iend
 type(kdnode) :: nodel,noder

 ! climb up tree, computing properties of parent nodes from their children
 !$omp parallel default(none) &
 !$omp shared(leaf_level) &
 !$omp shared(node,inoderange,xyzh,iorder) &
 !$omp private(level,inode,istart,iend,il,ir,nodel,noder)
 do level=leaf_level-1,0,-1
    call get_node_list(level,istart,iend)
    !$omp do
    do inode=2**level,2**(level+1)-1    ! nodes stored using Gafton & Rosswog index order
       ! get child nodes
       call get_child_nodes(inode,il,ir)
       nodel = node(il)
       noder = node(ir)
       call add_child_nodes(nodel,noder,node(inode))
       !print*,' approx size =',node(inode)%size
       call recompute_node_size(node(inode),inoderange(1,inode),inoderange(2,inode),xyzh,iorder)
       !print*,' recomputed size =',node(inode)%size
    enddo
    !$omp enddo
 enddo
 !$omp end parallel

end subroutine climb_tree

!---------------------------------------------
!+
!  as for climb tree, but only updating hmax
!+
!---------------------------------------------
subroutine climb_tree_hmax(node,leaf_level)
 type(kdnode), intent(inout) :: node(:)
 integer,      intent(in)    :: leaf_level
 integer :: level,inode,il,ir,istart,iend

 !$omp parallel default(none) &
 !$omp shared(leaf_level) &
 !$omp shared(node) &
 !$omp private(level,inode,istart,iend,il,ir)
 do level=leaf_level-1,0,-1
    call get_node_list(level,istart,iend)
    !$omp do
    do inode=istart,iend
       call get_child_nodes(inode,il,ir)       ! get child nodes
       node(inode)%hmax = max(node(il)%hmax,node(ir)%hmax)
    enddo
    !$omp enddo
 enddo
 !$omp end parallel

end subroutine climb_tree_hmax

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
!  get sort key
!+
!-----------------------------------
subroutine get_sort_key(np,xyzh,xmin,xmax,rkey)
 use hilbert, only:p_to_h
 integer, intent(in) :: np
 real, intent(in) :: xyzh(:,:),xmin(3),xmax(3)
 real, intent(out) :: rkey(np)
 integer :: i
 integer(8) :: ip(3),ngrid
 integer(8) :: iz,m
 real :: dx1(3)

 ngrid = np ! size of discrete grid (np**3)
 m = nint(log(real(ngrid))/log(2.))
 ngrid = 2**m

 dx1 = 1./(xmax - xmin)*ngrid
 !$omp parallel do default(none) &
 !$omp shared(xyzh,xmin,dx1,rkey,np,m) &
 !$omp private(i,ip,iz) &
 !$omp schedule(static)
 do i=1,np
    ip(1) = int((xyzh(1,i) - xmin(1))*dx1(1))
    ip(2) = int((xyzh(2,i) - xmin(2))*dx1(2))
    ip(3) = int((xyzh(3,i) - xmin(3))*dx1(3))
    iz = p_to_h(ip,m)
    if (iz < 0) stop 'exceeded hilbert key precision'
    rkey(i) = real(iz)
 enddo
 !$omp end parallel do

end subroutine get_sort_key

!-----------------------------------
!+
!  sort particles into desired order
!  based on the sort key
!+
!-----------------------------------
subroutine parallel_sort(np,rkey,iorder)
 use sortutils, only:indexx
 integer, intent(in)  :: np
 real, intent(in) :: rkey(np)
 integer, intent(out) :: iorder(np)

 call indexx(np,rkey,iorder)

end subroutine parallel_sort

!-----------------------------------
!+
!  compute centre of mass, size
!  hmax and quadrupole moments
!  for a given node
!+
!-----------------------------------
subroutine construct_node(nodei,i1,i2,xyzh,iorder)
 use part, only:massoftype,iphase,iamtype,npartoftype,maxtypes,maxphase,maxp,igas
 type(kdnode), intent(inout) :: nodei
 integer, intent(in) :: i1,i2
 real,    intent(in) :: xyzh(:,:)
 integer, intent(in) :: iorder(:)
 real :: x0(3),pmassi,fac,dfac
 real :: xi,yi,zi,hi,hmax,dx,dy,dz,dr2,r2max,totmass_node
#ifdef GRAVITY
 real :: quads(6)
#endif
 integer :: i,j

 !if (i2 < i1) return
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
    i = iorder(j)
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    hmax  = max(hmax,hi)
    if (maxphase==maxp) then
       pmassi = massoftype(iamtype(iphase(i)))
       fac    = pmassi*dfac ! to avoid round-off error
    endif
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
    i = iorder(j)
    dx = xyzh(1,i) - x0(1)
    dy = xyzh(2,i) - x0(2)
    dz = xyzh(3,i) - x0(3)
    dr2   = dx*dx + dy*dy + dz*dz
    r2max = max(r2max,dr2)
#ifdef GRAVITY
    if (maxphase==maxp) then
       pmassi = massoftype(iamtype(iphase_soa(i)))
    endif
    quads(1) = quads(1) + pmassi*(3.*dx*dx - dr2)  ! Q_xx
    quads(2) = quads(2) + pmassi*(3.*dx*dy)        ! Q_xy = Q_yx
    quads(3) = quads(3) + pmassi*(3.*dx*dz)        ! Q_xz = Q_zx
    quads(4) = quads(4) + pmassi*(3.*dy*dy - dr2)  ! Q_yy
    quads(5) = quads(5) + pmassi*(3.*dy*dz)        ! Q_yz = Q_zy
    quads(6) = quads(6) + pmassi*(3.*dz*dz - dr2)  ! Q_zz
#endif
 enddo
 if (sqrt(r2max) > 20.) then
    !print*,'LARGE NODE: SIZE=',sqrt(r2max),x0
   !do j=i1,i2
       !print*,' -> part ',iorder(j),xyzh(:,iorder(j)),pkey(iorder(j))
    !enddo
 endif
 ! assign properties to node
 nodei%xcen    = x0(:)
 nodei%size    = sqrt(r2max) + epsilon(r2max)
 nodei%hmax    = hmax
#ifdef GRAVITY
 nodei%mass    = totmass_node
 nodei%quads   = quads
#endif

end subroutine construct_node

!-----------------------------------
!+
!  compute node size exactly
!  from summing over node particles
!+
!-----------------------------------
subroutine recompute_node_size(nodei,i1,i2,xyzh,list)
 type(kdnode), intent(inout) :: nodei
 integer,      intent(in) :: i1,i2
 real,         intent(in) :: xyzh(:,:)
 integer,      intent(in) :: list(:)
 real :: x0,y0,z0,dx,dy,dz,dr2,r2max
 integer :: i,j

 x0 = nodei%xcen(1)
 y0 = nodei%xcen(2)
 z0 = nodei%xcen(3)
 r2max = 0.
 do j=i1,i2
    i = list(j)
    dx = xyzh(1,i) - x0
    dy = xyzh(2,i) - y0
    dz = xyzh(3,i) - z0
    dr2   = dx*dx + dy*dy + dz*dz
    r2max = max(r2max,dr2)
 enddo
 nodei%size = sqrt(r2max) + epsilon(r2max)

end subroutine recompute_node_size

!----------------------------------------------------------------
!+
!  Routine to walk tree for neighbour search
!  (all particles within a given h_i and optionally within h_j)
!+
!----------------------------------------------------------------
subroutine getneigh1(node,xpos,xsizei,rcuti,ndim,listneigh,nneigh,xyzh,xyzcache,ixyzcachesize,&
& get_hj,fnode,remote_export)
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
 use kdtree,   only:lenfgrav
 type(kdnode), intent(in)           :: node(ncellsmax+1)
 integer, intent(in)                :: ndim,ixyzcachesize
 real,    intent(in)                :: xpos(ndim)
 real,    intent(in)                :: xsizei,rcuti
 integer, intent(out)               :: listneigh(maxneigh)
 integer, intent(out)               :: nneigh
 real,    intent(in)                :: xyzh(:,:)
 real,    intent(out)               :: xyzcache(:,:)
 logical, intent(in)                :: get_hj
 real,    intent(out),    optional  :: fnode(lenfgrav)
 logical, intent(out),    optional  :: remote_export(:)
 integer, parameter :: istacksize = 300
 integer :: maxcache
 integer :: nstack(istacksize)
 integer :: n,istack,il,ir,npnode,i1,i2,i,j,jmax,nn
 real :: dx,dy,dz,xsizej,rcutj
 real :: rcut,rcut2,r2
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
#ifdef GRAVITY
    totmass_node = node(n)%mass
    quads        = node(n)%quads
#endif
    call get_child_nodes(n,il,ir)
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
       rcutj = radkern*node(n)%hmax
       rcut  = max(rcuti,rcutj)
    endif
    rcut2 = (xsizei + xsizej + rcut)**2   ! node size + search radius
#ifdef GRAVITY
    open_tree_node = tree_acc2*r2 < xsizej*xsizej    ! tree opening criterion for self-gravity
#endif
    if ((r2 < rcut2) .or. open_tree_node) then
       if_leaf: if (n >= 2**maxlevel) then ! once we hit a leaf node, retrieve contents into trial neighbour cache
          if_global_walk: if (present(remote_export)) then
         !    ! id is stored in ipart as id + 1
         !    if (ifirstincell(n) /= (id + 1)) then
         !       remote_export(ifirstincell(n)) = .true.
         !    endif
          else
             i1=inoderange(1,n)
             i2=inoderange(2,n)
!             if (nneigh > 100000) print*,'NODE ',n,' adding ',npnode,rcuti,xsizei,xsizej
             npnode = i2 - i1 + 1
             listneigh(nneigh+1:nneigh+npnode) = iorder(i1:i2)
             if (nneigh < ixyzcachesize) then ! not strictly necessary, loop doesn't execute anyway
                nn = min(nneigh+npnode,ixyzcachesize)
                !print*,' got ',nneigh+npnode,' capping to ',nn
                jmax = max(nn-nneigh,0)
                !print*,'cache build ',nneigh+1,'->',nneigh+jmax,'parts ',i1,'->',i1+jmax-1,&
                  !  ' nneigh = ',nneigh+1,'->',nneigh+npnode
                if (maxcache >= 4) then
                   do j=1,jmax
                      i = iorder(i1+j-1)
                      xyzcache(1,nneigh+j) = xyzh(1,i) + xoffset
                      xyzcache(2,nneigh+j) = xyzh(2,i) + yoffset
                      xyzcache(3,nneigh+j) = xyzh(3,i) + zoffset
                      xyzcache(4,nneigh+j) = 1./xyzh(4,i)
                   enddo
                else
                   do j=1,jmax
                      i = iorder(i1+j-1)
                      xyzcache(1,nneigh+j) = xyzh(1,i) + xoffset
                      xyzcache(2,nneigh+j) = xyzh(2,i) + yoffset
                      xyzcache(3,nneigh+j) = xyzh(3,i) + zoffset
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

end subroutine getneigh1

end module dptree
