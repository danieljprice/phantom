module dptree
 use dim, only:ncellsmax,minpart
 use dtypekdtree, only:kdnode
 implicit none
 integer, parameter :: ndimtree = 3
 !real, public :: tree_accuracy = 0.5

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

 public :: maketree1
 private

 integer, public :: maxlevel
 integer, allocatable :: iorder(:)
 integer, allocatable :: inoderange(:,:)
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
 integer :: leaf_level
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
 call build_tree_index(np,minpart,leaf_level,inoderange,ncells)

 ! get properties of all leaf nodes (in parallel)
 call build_leaf_nodes(leaf_level,node,xyzh,iorder)

 ! propagate information to parent nodes
 call climb_tree(node,inoderange,leaf_level,xyzh,iorder)

 maxlevel = leaf_level
 !return
    write(*,"(a,i10,3(a,i2))") ' maketree1: nodes = ',ncells,', leaf level = ',maxlevel
    block
       real :: sizesum
       integer :: i,istart,iend,nleaf,level
       print*,' size root node = ',node(1)%size
       do level=leaf_level,0,-1
          call get_node_list(level,istart,iend)
          sizesum = 0.
          do i=istart,iend
             sizesum = sizesum + node(i)%size
          enddo
          nleaf = iend - istart + 1
          if (level==leaf_level) print*,' mean parts per leaf node = ',np/real(nleaf)
          print*,level,' SCORE = ',node(1)%size*real(nleaf)/sizesum,' mean size = ',&
                       sizesum/real(nleaf),'sum=',sizesum,' nnodes=',nleaf
       enddo
    end block

end subroutine maketree1

!-----------------------------------
!+
!  Build the tree index, i.e. index
!  of first and last particle on
!  every node
!+
!-----------------------------------
subroutine build_tree_index(np,parts_per_leaf,leaf_level,inoderange,ncells)
 use io, only:fatal
 integer, intent(in)  :: np,parts_per_leaf
 integer, intent(out) :: leaf_level
 integer, intent(out) :: inoderange(:,:)
 integer(kind=8), intent(out) :: ncells
 real :: ratio,parts_per_leafi
 integer :: nleaf,inode,index,istart,iend
 integer :: level,il,ir

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
 !$omp parallel do default(none) schedule(static) &
 !$omp shared(istart,iend,parts_per_leafi,inoderange) &
 !$omp private(inode,index)
 do inode=istart,iend
    index = inode - istart
    inoderange(1,inode) = int(index*parts_per_leafi) + 1
    inoderange(2,inode) = int(((index+1)*parts_per_leafi))
 enddo
 !$omp end parallel do
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

end module dptree
