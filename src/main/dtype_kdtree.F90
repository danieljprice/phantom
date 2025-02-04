!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dtypekdtree
!
! None
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: io, mpi, mpiutils
!
 implicit none

#ifdef TREEVIZ
 integer, parameter :: ndimtree = 2  ! 2D for visualisation/debugging only
#else
 integer, parameter :: ndimtree = 3
#endif

 integer, parameter :: kdnode_bytes = &
                      8*ndimtree &  ! xcen(ndimtree)
                    + 8 &           ! size
                    + 8 &           ! hmax
                    + 4 &           ! leftchild
                    + 4 &           ! rightchild
                    + 4 &           ! parent
#ifdef GRAVITY
                    + 8 &           ! mass
                    + 8*6 &         ! quads(6)
#endif
#ifdef TREEVIZ
 + 8*ndimtree &  ! xmin(ndimtree)
                    + 8*ndimtree &  ! xmax(ndimtree)
#endif
                    + 0

 private
 public :: ndimtree
 public :: kdnode
 public :: kdnode_bytes
 public :: get_mpitype_of_kdnode
 type kdnode
    sequence
    real :: xcen(ndimtree)
    real :: size
    real :: hmax
    real :: dum   ! avoid ifort warning: align on 4-byte boundary
    integer :: leftchild
    integer :: rightchild
    integer :: parent
    integer :: idum ! avoid ifort warning: align on 4-byte boundary
#ifdef GRAVITY
    real :: mass
    real :: quads(6)
#endif
#ifdef TREEVIZ
    real :: xmin(ndimtree)
    real :: xmax(ndimtree)
#endif
 end type kdnode

contains

!----------------------------------------------------------------
!+
!  get the MPI datatype for kdnode
!  Todo:
!  - MPI commit the datatype just once, rather than rebuilding it every time
!+
!----------------------------------------------------------------
subroutine get_mpitype_of_kdnode(dtype)
#ifdef MPI
 use mpi
 use mpiutils, only:mpierr
 use io,       only:error

 integer, parameter              :: ndata = 20
 integer, intent(out)            :: dtype
 integer                         :: nblock, blens(ndata), mpitypes(ndata)
 integer(kind=MPI_ADDRESS_KIND)  :: disp(ndata)
 type(kdnode)                    :: node
 integer(kind=MPI_ADDRESS_KIND)  :: addr,start,lb,extent

 nblock = 0

 call MPI_GET_ADDRESS(node,start,mpierr)

 nblock = nblock + 1
 blens(nblock) = size(node%xcen)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%xcen,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%size,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%hmax,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(node%leftchild,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(node%rightchild,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(node%parent,addr,mpierr)
 disp(nblock) = addr - start

#ifdef GRAVITY
 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%mass,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(node%quads)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%quads,addr,mpierr)
 disp(nblock) = addr - start
#endif
#ifdef TREEVIZ
 nblock = nblock + 1
 blens(nblock) = size(node%xmin)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%xmin,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(node%xmax)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(node%xmax,addr,mpierr)
 disp(nblock) = addr - start
#endif

 call MPI_TYPE_CREATE_STRUCT(nblock,blens(1:nblock),disp(1:nblock),mpitypes(1:nblock),dtype,mpierr)
 call MPI_TYPE_COMMIT(dtype,mpierr)

 ! check extent okay
 call MPI_TYPE_GET_EXTENT(dtype,lb,extent,mpierr)
 if (extent /= sizeof(node)) then
    call error('dtype_kdtree','MPI_TYPE_GET_EXTENT has calculated the extent incorrectly')
 endif

#else
 integer, intent(out) :: dtype
 dtype = 0
#endif
end subroutine get_mpitype_of_kdnode

end module dtypekdtree
