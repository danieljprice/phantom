!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for debugging/visualisation of kdtree build
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, giza, io, kdtree, kernel, linklist
!+
!--------------------------------------------------------------------------
module analysis
 use kdtree, only:kdnode
 implicit none
 character(len=20), parameter, public :: analysistype = 'kdtree'
 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use kdtree,   only:node,minpart
 use linklist, only:set_linklist,ncells
 use io,  only:iverbose
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 iverbose= 3
 print*,'> Building 2D kd-tree... minpart = ',minpart
 call set_linklist(npart,npart,xyzh,vxyzu)

 !print*,'> Printing kd-tree '
 !do i=1,numnodes
 !   print*,node(i)
 !enddo

 call viz_kdtree(xyzh,npart,node,int(ncells))

end subroutine do_analysis

subroutine viz_kdtree(xyzh,np,tree,numnodes)
 use kdtree, only:irootnode,ndimtree
 use giza
 integer,      intent(in) :: np,numnodes
 real,         intent(in) :: xyzh(:,:)
 type(kdnode), intent(in) :: tree(:)
 integer :: maxp
 real :: xmin(ndimtree),xmax(ndimtree)

 call giza_open('/eps')
 do maxp=1,ndimtree
    xmin(maxp) = minval(xyzh(maxp,1:np))
    xmax(maxp) = maxval(xyzh(maxp,1:np))
 enddo
 call giza_set_environment(xmin(1),xmax(1),xmin(2),xmax(2),1,0)
 call giza_label('x','y','')

 !--plot particles
 print*,'> Plotting particles '
 call giza_points(np,xyzh(1,1:np),xyzh(2,1:np),1)

#ifdef TREEVIZ
 !--plot lines bounding nodes
 print*,'> Plotting kd tree '
 call giza_set_line_width(1.5)
 call giza_set_colour_index(2)
 call plot_nodes(irootnode,0,ndimtree,tree)
#else
 !--visualise the neighbour finding
 print*,'> Checking neighbour find'
 call check_neighbours(xyzh,tree)
#endif
 call giza_close()

end subroutine viz_kdtree

#ifdef TREEVIZ
recursive subroutine plot_nodes(inode,level,ndim,tree)
 use kdtree, only:labelax
 use giza
 use linklist, only:ifirstincell
 !use part, only:xyzh,ll
 integer,      intent(in) :: inode,level,ndim
 type(kdnode), intent(in) :: tree(:)
 real :: xpivot
 real :: xpts(ndim,2)
 integer :: iaxis
 logical :: plot_pivot_lines_only, plot_node_centres, plot_leaf_node_sizes

 plot_leaf_node_sizes  = .false.
 plot_node_centres     = .false.
 plot_pivot_lines_only = .true.

 call giza_set_fill(giza_fill_hollow)
 !print*,' level = ',level
 if (inode <= 0) return
 if (ifirstincell(inode) /= 0) then
    !print*,' node ',inode,' is on level ',level,' kids = ',tree(inode)%leftchild,tree(inode)%rightchild,&
    !    ' max h = ',tree(inode)%hmax
    if (plot_leaf_node_sizes) then
       call giza_set_colour_index(3)
       call giza_single_point(tree(inode)%xcen(1),tree(inode)%xcen(2),3)
       call giza_circle(tree(inode)%xcen(1),tree(inode)%xcen(2),tree(inode)%size)
    endif
    !print*,'hit leaf on level ',level
    !return
 endif
 !xmaxtemp(1:ndim) = 2.*tree(inode)%xcen(1:ndim) - tree(inode)%xmin(1:ndim)

 iaxis = maxloc(tree(inode)%xmax(1:ndim) - tree(inode)%xmin(1:ndim),1)
 if (iaxis < 1 .or. iaxis > 2) stop 'iaxis out of range'
 xpivot = tree(inode)%xcen(iaxis)
 print*,' node ',inode,' is on level ',level,' kids = ',tree(inode)%leftchild,tree(inode)%rightchild,&
        ' cofm = ',tree(inode)%xcen,' max h = ',tree(inode)%hmax

 xpts(1:ndim,1) = tree(inode)%xmin(1:ndim)
 xpts(1:ndim,2) = tree(inode)%xmax(1:ndim)

 ! plot node centre as a cross hair
 if (plot_node_centres) then
    call giza_set_character_height(3.0)
    call giza_single_point(tree(inode)%xcen(1),tree(inode)%xcen(2),5)
 endif

 ! plot the line where this node is split, or else full rectangle
 if (plot_pivot_lines_only) then
    if (ifirstincell(inode)==0) then  ! only if node is NOT a leaf node
       xpts(iaxis,1:2) = xpivot
       call giza_line(2,xpts(1,:),xpts(2,:))
    endif
 else
    call giza_rectangle(tree(inode)%xmin(1),tree(inode)%xmax(1),tree(inode)%xmin(2),tree(inode)%xmax(2))
 endif
 !print*,'plotting node ',inode
 !print*,labelax(1),' = ',tree(inode)%xmin(1),' -> ',tree(inode)%xmax(1)
 !print*,labelax(2),' = ',tree(inode)%xmin(2),' -> ',tree(inode)%xmax(2)
 !if (xpivot < tree(inode)%xmin(iaxis) .or. xpivot > xmaxtemp(iaxis)) then
 !   print*,'ERROR: pivot out of range ',inode,xpivot,tree(inode)%xmin(iaxis),xmaxtemp(iaxis)
 !endif
 !read*

!--call the left daughter
 !print*,'starting left ',level,tree(inode)%leftchild
 call plot_nodes(tree(inode)%leftchild,level+1,ndim,tree)

!--call the right daughter
 !print*,'starting right ',level,tree(inode)%rightchild
 call plot_nodes(tree(inode)%rightchild,level+1,ndim,tree)

end subroutine plot_nodes
#endif

subroutine check_neighbours(xyzh,tree)
 use dim, only:maxneigh
 use linklist, only:ncells,get_neighbour_list,ifirstincell
 use giza
 use kernel, only:radkern
 real,         intent(in) :: xyzh(:,:)
 type(kdnode), intent(in) :: tree(:)
 integer :: listneigh(maxneigh)
 real    :: xyzcache(3,1)
 integer :: nneigh,inode,i

 over_cells: do inode=1,ncells
    if ((norm2(tree(inode)%xcen(1:3) - (/-0.25,0.,0./)) < 5.e-2) .and. ifirstincell(inode) > 0) then
       call get_neighbour_list(inode,listneigh,nneigh,xyzh,xyzcache,0)

       ! plot all trial neighbours
       call giza_set_character_height(1.0)
       call giza_set_colour_index(4)
       do i=1,nneigh
          call giza_single_point(xyzh(1,listneigh(i)),xyzh(2,listneigh(i)),17)
       enddo

       ! plot search circle
       call giza_set_character_height(2.0)
       call giza_set_line_width(1.5)
       call giza_set_colour_index(2)
       call giza_single_point(tree(inode)%xcen(1),tree(inode)%xcen(2),3)
       call giza_set_fill(giza_fill_hollow)

       call giza_circle(tree(inode)%xcen(1),tree(inode)%xcen(2),&
                        tree(inode)%size + radkern*tree(inode)%hmax)

       !call giza_rectangle(tree(inode)%xmin(1),tree(inode)%xmax(1),tree(inode)%xmin(2),tree(inode)%xmax(2))
!       call giza_circle(tree(inode)%xcen(1),tree(inode)%xcen(2),tree(inode)%size + 0.1)
       print*,'nneigh = ',nneigh

       exit over_cells
    endif
 enddo over_cells

end subroutine check_neighbours

end module
