!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpidens
!
! None
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, mpi
!
 use io,       only:nprocs,fatal
 use dim,      only:minpart,maxrhosum,maxxpartvecidens

 implicit none
 private

 public :: celldens
 public :: stackdens
 public :: get_mpitype_of_celldens
 public :: free_mpitype_of_celldens

 integer, parameter :: ndata = 19 ! number of elements in the cell (including padding)
 integer, parameter :: nbytes_celldens = 8 * minpart                    + & !  h(minpart)
                                         8 * minpart                    + & !  h_old(minpart)
                                         8 * maxxpartvecidens * minpart + & !  xpartvec(maxxpartvecidens,minpart)
                                         8 * maxrhosum * minpart        + & !  rhosums(maxrhosum,minpart)
                                         8 * 3                          + & !  xpos(3)
                                         8                              + & !  xsizei
                                         8                              + & !  rcuti
                                         8                              + & !  hmax
                                         4                              + & !  icell
                                         4                              + & !  npcell
                                         4 * minpart                    + & !  arr_index(minpart)
                                         4                              + & !  owner
                                         4                              + & !  nits
                                         4                              + & !  nneightry
                                         4 * minpart                    + & !  nneigh(minpart)
                                         4                              + & !  waiting_index
                                         1 * minpart                    + & !  iphase(minpart)
                                         1 * minpart                        !  apr_level

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
    integer(kind=1)  :: iphase(minpart)
    integer(kind=1)  :: apr(minpart)                           ! apr resolution level (not in xpartvec because integer)

    ! pad the array to 8-byte boundaries
    integer(kind=1)  :: pad(8 - mod(nbytes_celldens, 8))
 end type celldens

 type stackdens
    sequence
    type(celldens), pointer   :: cells(:)
    integer                   :: maxlength = 0
    integer                   :: n = 0
    integer                   :: number
    integer                   :: idum     ! to avoid ifort warning
 end type stackdens

contains

subroutine get_mpitype_of_celldens(dtype)
#ifdef MPI
 use mpi
 use io,       only:error
#endif
 integer, intent(out)            :: dtype
#ifdef MPI
 integer                         :: nblock, blens(ndata), mpitypes(ndata)
 integer(kind=MPI_ADDRESS_KIND)  :: disp(ndata)

 type(celldens)                  :: cell
 integer(kind=MPI_ADDRESS_KIND)  :: addr,start,lb,extent
 integer                         :: mpierr

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
 blens(nblock) = size(cell%iphase)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%iphase,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%apr)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%apr,addr,mpierr)
 disp(nblock) = addr - start

 ! padding must come last
 nblock = nblock + 1
 blens(nblock) = 8 - mod(4 * (6 + 2 * minpart) + 2*minpart, 8)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%pad,addr,mpierr)
 disp(nblock) = addr - start

 call MPI_TYPE_CREATE_STRUCT(nblock,blens(1:nblock),disp(1:nblock),mpitypes(1:nblock),dtype,mpierr)
 call MPI_TYPE_COMMIT(dtype,mpierr)

 ! check extent okay
 call MPI_TYPE_GET_EXTENT(dtype,lb,extent,mpierr)
 if (extent /= sizeof(cell)) then
    call fatal('mpi_dens','MPI_TYPE_GET_EXTENT has calculated the extent incorrectly')
 endif

#else
 dtype = 0
#endif
end subroutine get_mpitype_of_celldens

subroutine free_mpitype_of_celldens(dtype)
#ifdef MPI
 use mpi
#endif
 integer, intent(inout) :: dtype
#ifdef MPI
 integer                :: mpierr

 call MPI_Type_free(dtype,mpierr)
#else
 dtype = 0
#endif

end subroutine free_mpitype_of_celldens

end module mpidens
