!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: mpidens
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Conrad Chan
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: dim, io, mpi, mpiutils
!+
!--------------------------------------------------------------------------
module mpidens
 use io,       only:nprocs,fatal
 use dim,      only:minpart,maxrhosum,maxprocs,stacksize,maxxpartvecidens
#ifdef MPI
 use mpiutils, only:mpierr
#endif
 implicit none
 private

 public :: celldens
 public :: stackdens
#ifdef MPI
 public :: get_mpitype_of_celldens
#endif

 type celldens
    sequence
    real             :: h(minpart)                             ! don't put this in xpartvec because it is modified inplace
    real             :: h_old(minpart)                         ! original h
    real             :: xpartvec(maxxpartvecidens,minpart)
    real             :: rhosums(maxrhosum,minpart)
    real             :: xpos(3)
    real             :: xsizei
    real             :: rcuti
    real             :: hmax
    integer          :: icell
    integer          :: npcell                                 ! number of particles in here
    integer          :: arr_index(minpart)
    integer          :: owner                                  ! id of the process that owns this
    integer          :: nits                                   ! number of density iterations done so far
    integer          :: nneightry
    integer          :: nneigh(minpart)                        ! number of actual neighbours (diagnostic)
    integer          :: waiting_index
    logical          :: remote_export(maxprocs)                ! remotes we are waiting for
    integer(kind=1)  :: iphase(minpart)
    integer(kind=1)  :: pad(8 - mod(4 * (6 + 2 * minpart + maxprocs) + minpart, 8))
 endtype

 type stackdens
    sequence
    type(celldens), pointer   :: cells(:)
    integer                   :: maxlength = 0
    integer                   :: n = 0
    integer                   :: mem_start
    integer                   :: mem_end
 endtype

contains

#ifdef MPI
subroutine get_mpitype_of_celldens(dtype)
 use mpi

 integer, parameter              :: ndata = 20

 integer, intent(out)            :: dtype
 integer                         :: dtype_old
 integer                         :: nblock, blens(ndata), mpitypes(ndata)
 integer(kind=4)                 :: disp(ndata)

 type(celldens)                 :: cell
 integer(kind=MPI_ADDRESS_KIND)  :: addr,start,lb,extent

 nblock = 0

 call MPI_GET_ADDRESS(cell,start,mpierr)

 nblock = nblock + 1
 blens(nblock) = size(cell%h)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%h,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%h_old)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%h_old,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%xpartvec)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xpartvec,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%rhosums)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%rhosums,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%xpos)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xpos,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xsizei,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%rcuti,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%hmax,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%icell,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%npcell,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%arr_index)
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%arr_index,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%owner,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nits,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nneightry,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%nneigh)
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nneigh,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%waiting_index,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%remote_export)
 mpitypes(nblock) = MPI_LOGICAL
 call MPI_GET_ADDRESS(cell%remote_export,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%iphase)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%iphase,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 8 - mod(4 * (6 + 2 * minpart + maxprocs) + minpart, 8)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%pad,addr,mpierr)
 disp(nblock) = addr - start

 call MPI_TYPE_STRUCT(nblock,blens(1:nblock),disp(1:nblock),mpitypes(1:nblock),dtype,mpierr)
 call MPI_TYPE_COMMIT(dtype,mpierr)

 ! check extent okay
 call MPI_TYPE_GET_EXTENT(dtype,lb,extent,mpierr)
 if (extent /= sizeof(cell)) then
    dtype_old = dtype
    lb = 0
    extent = sizeof(cell)
    call MPI_TYPE_CREATE_RESIZED(dtype_old,lb,extent,dtype,mpierr)
    call MPI_TYPE_COMMIT(dtype,mpierr)
    call MPI_TYPE_FREE(dtype_old,mpierr)
 endif

end subroutine get_mpitype_of_celldens
#endif

end module mpidens
