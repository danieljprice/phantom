!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: mpiforce
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
module mpiforce
 use io,       only:nprocs,fatal
 use dim,      only:minpart,maxfsum,maxprocs,stacksize,maxxpartveciforce
#ifdef MPI
 use mpiutils, only:mpierr
#endif
 implicit none
 private

 public :: cellforce
 public :: stackforce
#ifdef MPI
 public :: get_mpitype_of_cellforce
#endif

 type cellforce
    sequence
    real             :: xpartvec(maxxpartveciforce,minpart)
    real             :: fsums(maxfsum,minpart)
    real             :: fgrav(20)
    real             :: xpos(3)
    real             :: xsizei
    real             :: rcuti
    real             :: dtdrag(minpart)
    real             :: vsigmax(minpart)
    integer          :: icell
    integer          :: npcell                                 ! number of particles in here
    integer          :: arr_index(minpart)
    integer          :: ndrag
    integer          :: nstokes
    integer          :: nsuper
    integer          :: owner                                  ! id of the process that owns this
    integer          :: waiting_index
    logical          :: remote_export(maxprocs)                ! remotes we are waiting for
    integer(kind=1)  :: iphase(minpart)
    integer(kind=1)  :: ibinneigh(minpart)
    integer(kind=1)  :: pad(8 - mod(4 * (7 + minpart + maxprocs) + 2*minpart, 8)) !padding to maintain alignment of elements
 endtype

 type stackforce
    sequence
    type(cellforce), pointer  :: cells(:)
    integer                   :: maxlength = 0
    integer                   :: n = 0
    integer                   :: mem_start
    integer                   :: mem_end
 endtype

contains

#ifdef MPI
subroutine get_mpitype_of_cellforce(dtype)
 use mpi

 integer, parameter              :: ndata = 20

 integer, intent(out)            :: dtype
 integer                         :: dtype_old
 integer                         :: nblock, blens(ndata), mpitypes(ndata)
 integer(kind=4)                 :: disp(ndata)

 type(cellforce)                 :: cell
 integer(kind=MPI_ADDRESS_KIND)  :: addr,start,lb,extent

 nblock = 0

 call MPI_GET_ADDRESS(cell,start,mpierr)

 nblock = nblock + 1
 blens(nblock) = size(cell%xpartvec)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%xpartvec,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%fsums)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%fsums,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%fgrav)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%fgrav,addr,mpierr)
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
 blens(nblock) = size(cell%dtdrag)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%dtdrag,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%vsigmax)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%vsigmax,addr,mpierr)
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
 call MPI_GET_ADDRESS(cell%ndrag,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nstokes,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%nsuper,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = 1
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%owner,addr,mpierr)
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
 blens(nblock) = 8 - mod(4 * (7 + minpart + maxprocs) + minpart, 8)
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

end subroutine get_mpitype_of_cellforce
#endif

end module mpiforce
