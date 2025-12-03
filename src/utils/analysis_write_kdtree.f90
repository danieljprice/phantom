!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which computes 3D kd-tree and writes to file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, neighkdtree
!
 use dim, only:gravity
 implicit none
 character(len=20), parameter, public :: analysistype = 'write_kdtree'

 public :: do_analysis,write_kdtree_file,read_kdtree_file

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use neighkdtree, only:build_tree
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real, allocatable :: dumxyzh(:,:)
 integer :: np

 !****************************************
 ! 1. Build kdtree
 ! --> global (shared) neighbour lists for all particles in tree cell
 !****************************************

 print*, 'Building kdtree: '

 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 np = npart
 call build_tree(np,np,dumxyzh,vxyzu)

 print*, '- Done'

 !********************************************
 ! 2. Write kdtree to binary file
 !*******************************************
 call write_kdtree_file(dumpfile)

end subroutine do_analysis

!--------------------------------------------------------------------
!+
! Writes 3D kd-tree to binary file
!+
!--------------------------------------------------------------------
subroutine write_kdtree_file(dumpfile)
 use neighkdtree, only:ncells,node
 character(len=*), intent(in) :: dumpfile
 character(9) :: filetag
 character(100) :: treefile
 integer :: iu

 treefile = 'kdtree_'//trim(dumpfile)
 print'(a,a)', 'Writing kdtree to binary file ', trim(treefile)

 ! Write tag indicating if this is from a run with or without gravity
 if (gravity) then
    filetag = 'gravity'
    print '(a,a,i7)', 'This file contains masses: ', filetag, ncells
 else
    filetag = 'nogravity'
    print '(a,a,i7)', 'This file does not contains masses: ', filetag, ncells
 endif

 open(newunit=iu,file=treefile,form='unformatted')

! Write header data
 write(iu) filetag, ncells
! Now write tree data
 write(iu) node(1:ncells)%xcen(1)
 write(iu) node(1:ncells)%xcen(2)
 write(iu) node(1:ncells)%xcen(3)
 write(iu) node(1:ncells)%size
! Node mass only stored if gravity active
 if (gravity) write(iu) node(1:ncells)%mass
 write(iu) node(1:ncells)%leftchild
 write(iu) node(1:ncells)%rightchild
 write(iu) node(1:ncells)%parent
 close(iu)

 print '(a)', 'kdtree file write complete'

end subroutine write_kdtree_file

!--------------------------------------------------------------------
!+
! Writes 3D kd-tree to binary file
!+
!--------------------------------------------------------------------
subroutine read_kdtree_file(dumpfile)
 use neighkdtree, only:ncells,node
 character(len=*), intent(in) :: dumpfile
 character(7) :: filetag
 character(100) :: treefile
 integer :: iu

 treefile = 'kdtree_'//trim(dumpfile)
 print'(a,a)', 'Reading kdtree from binary file ', trim(treefile)

 open(newunit=iu,file=treefile,form='unformatted')
! Read header
 read(iu) filetag, ncells

 if (filetag=='gravity') then
    print '(a)','Tree from a gravity run'
 else
    print '(a)','Tree from a non-gravity run'
 endif

 print '(a,i7)', 'Tree has ',ncells, ' active cells'

! Now write tree data
 read(iu) node(1:ncells)%xcen(1)
 read(iu) node(1:ncells)%xcen(2)
 read(iu) node(1:ncells)%xcen(3)
 read(iu) node(1:ncells)%size
! Node mass only stored if gravity active
 if (gravity) read(iu) node(1:ncells)%mass
 read(iu) node(1:ncells)%leftchild
 read(iu) node(1:ncells)%rightchild
 read(iu) node(1:ncells)%parent
 close(iu)

 print '(a)', 'kdtree file read complete'

end subroutine read_kdtree_file

end module analysis
