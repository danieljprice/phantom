!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine which computes 3D kd-tree and writes to file
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: kdtree, linklist, part
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'write_kdtree'

 public :: do_analysis

 private

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 use part, only: iphase
 use linklist, only: set_linklist

 implicit none

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time

 real,allocatable,dimension(:,:) :: dumxyzh


 !****************************************
 ! 1. Build kdtree and linklist
 ! --> global (shared) neighbour lists for all particles in tree cell
 !****************************************

 print*, 'Building kdtree and linklist: '

 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

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

 use linklist, only: ncells
 use kdtree, only: node

 implicit none

 character(len=*), intent(in):: dumpfile
 character(7) :: filetag
 character(100) :: treefile
 integer :: icell

 treefile = 'kdtree_'//TRIM(dumpfile)
 print'(a,a)', 'Writing kdtree to binary file ', TRIM(treefile)

 ! Write tag indicating if this is from a run with or without gravity
#ifdef GRAVITY
 filetag = 'gravity'
 print '(a,a,I7)', 'This file contains masses: ', filetag, ncells
#else
 filetag = 'nogravi'
 print '(a,a,I7)', 'This file does not contains masses: ', filetag, ncells
#endif

 open(10,file=treefile, form='unformatted')

! Write header data
 write(10) filetag, ncells

! Now write tree data
 write(10) node(1:ncells)%xcen(1)
 write(10) node(1:ncells)%xcen(2)
 write(10) node(1:ncells)%xcen(3)
 write(10) node(1:ncells)%size

! Node mass only stored if gravity active
#ifdef GRAVITY
 write(10) node(1:ncells)%mass
#endif

 write(10) node(1:ncells)%leftchild
 write(10) node(1:ncells)%rightchild
 write(10) node(1:ncells)%parent

 close(10)

 print '(a)', 'kdtree file write complete'

end subroutine write_kdtree_file

!--------------------------------------------------------------------
!+
! Writes 3D kd-tree to binary file
!+
!--------------------------------------------------------------------
subroutine read_kdtree_file(dumpfile)

 use linklist, only: ncells
 use kdtree, only: node

 implicit none
 character(len=*), intent(in):: dumpfile
 character(7) :: filetag
 character(100) :: treefile

 treefile = 'kdtree_'//TRIM(dumpfile)
 print'(a,a)', 'Reading kdtree from binary file ', TRIM(treefile)

 open(10,file=treefile, form='unformatted')
! Read header
 read(10) filetag, ncells

 if (filetag=='gravity') then
    print'(a)', 'Tree from a gravity run'
 else
    print '(a)', 'Tree from a non-gravity run'
 endif

 print '(a,I7)', 'Tree has ',ncells, ' active cells'

! Now write tree data
 read(10) node(1:ncells)%xcen(1)
 read(10) node(1:ncells)%xcen(2)
 read(10) node(1:ncells)%xcen(3)
 read(10) node(1:ncells)%size

! Node mass only stored if gravity active
#ifdef GRAVITY
 read(10) node(1:ncells)%mass
#endif

 read(10) node(1:ncells)%leftchild
 read(10) node(1:ncells)%rightchild
 read(10) node(1:ncells)%parent

 close(10)

 print '(a)', 'kdtree file read complete'

end subroutine read_kdtree_file

end module
