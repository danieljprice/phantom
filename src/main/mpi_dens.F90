!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
 public :: pack_celldens,unpack_celldens
 public :: put_reals,get_reals

 !
 ! cells are sent as real buffers holding only what the other task needs:
 ! outbound (to the task with neighbours) the particle data, and on return
 ! (to the owner) the partial sums, both for the npcell particles only
 !
 integer, parameter, public :: maxpack_celldens = max(8 + minpart*(maxxpartvecidens + 4), &
                                                      5 + minpart*(maxrhosum + 2))
 ! number of density sums to return, set by densityiterate: the sums at the end
 ! of rhosum are only accumulated with unequal masses, otherwise they are zero
 integer, public :: nrhosum_pack = maxrhosum

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

 type chunkdens
    sequence
    type(celldens), pointer   :: cells(:) => null()
 end type chunkdens

 !
 ! stack of cells, stored in chunks of chunksize cells that are allocated as required,
 ! so the stack can grow from inside a parallel region without moving existing cells
 !
 type stackdens
    sequence
    type(chunkdens), pointer  :: chunks(:) => null()
    integer                   :: chunksize = 0
    integer                   :: maxlength = 0   ! number of cells that can be stored without growing
    integer                   :: n = 0
    integer                   :: ipad = 0        ! pads the type to a multiple of 8 bytes (ifort error #6380)
 end type stackdens

contains

subroutine get_mpitype_of_celldens(dtype)
#ifdef MPI
 use mpi
 use io,       only:error
#endif
 integer, intent(out) :: dtype
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

!-----------------------------------------------------------------------
!+
!  pack a cell into a real buffer: header (npcell, owner, waiting_index),
!  then either the particle data (outbound) or the density sums (return)
!+
!-----------------------------------------------------------------------
subroutine pack_celldens(cell,outbound,buf,n,mask)
 type(celldens), intent(in)  :: cell
 logical,        intent(in)  :: outbound
 real,           intent(out) :: buf(:)
 integer,        intent(out) :: n
 logical,        intent(in), optional :: mask(:)
 integer :: np,i,k,idx(minpart)

 !--particles sent: those selected by mask (outbound only), each with its
 !  position in the owner's cell, so that the results can be added back
 if (outbound .and. present(mask)) then
    np = 0
    do i=1,cell%npcell
       if (mask(i)) then
          np = np + 1
          idx(np) = i
       endif
    enddo
 elseif (outbound) then
    np = cell%npcell
    idx(1:np) = (/(i,i=1,np)/)
 else
    np = cell%npcell
    idx(1:np) = cell%arr_index(1:np)
 endif
 n  = 0
 call put_reals(buf,n,real((/np,cell%owner,cell%waiting_index/)))
 if (outbound) then
    call put_reals(buf,n,(/cell%xpos,cell%xsizei,cell%rcuti/))
    call put_reals(buf,n,cell%h(idx(1:np)))
    do k=1,np
       call put_reals(buf,n,cell%xpartvec(:,idx(k)))
    enddo
    call put_reals(buf,n,real(cell%iphase(idx(1:np))))
    call put_reals(buf,n,real(cell%apr(idx(1:np))))
 else
    call put_reals(buf,n,(/real(nrhosum_pack)/))
    call put_reals(buf,n,reshape(cell%rhosums(1:nrhosum_pack,1:np),(/nrhosum_pack*np/)))
    call put_reals(buf,n,real(cell%nneigh(1:np)))
    call put_reals(buf,n,(/real(cell%nneightry)/))
 endif
 call put_reals(buf,n,real(idx(1:np)))

end subroutine pack_celldens

!-----------------------------------------------------------------------
!+
!  unpack a cell packed by pack_celldens; n returns the number of reals
!  read. Fields that are not sent are set for a cell computed remotely.
!  arr_index returns the position of each particle in the owner's cell
!  (only on the owner is it the index into the tree, set by start_cell)
!+
!-----------------------------------------------------------------------
subroutine unpack_celldens(buf,outbound,cell,n)
 real,           intent(in)    :: buf(:)
 logical,        intent(in)    :: outbound
 type(celldens), intent(inout) :: cell
 integer,        intent(out)   :: n
 real    :: hdr(3),geom(5),tmp(minpart),one(1)
 integer :: np,nx

 n = 0
 call get_reals(buf,n,hdr)
 np = nint(hdr(1))
 cell%npcell        = np
 cell%owner         = nint(hdr(2))
 cell%waiting_index = nint(hdr(3))
 if (outbound) then
    nx = size(cell%xpartvec,1)
    call get_reals(buf,n,geom)
    cell%xpos   = geom(1:3)
    cell%xsizei = geom(4)
    cell%rcuti  = geom(5)
    call get_reals(buf,n,cell%h(1:np))
    cell%xpartvec(:,1:np) = reshape(buf(n+1:n+nx*np),(/nx,np/))
    n = n + nx*np
    call get_reals(buf,n,tmp(1:np))
    cell%iphase(1:np) = int(tmp(1:np),kind=1)
    call get_reals(buf,n,tmp(1:np))
    cell%apr(1:np) = int(tmp(1:np),kind=1)
    cell%nits      = 0
    cell%nneightry = 0
 else
    call get_reals(buf,n,one)
    nx = nint(one(1))
    cell%rhosums(1:nx,1:np) = reshape(buf(n+1:n+nx*np),(/nx,np/))
    cell%rhosums(nx+1:,1:np) = 0.
    n = n + nx*np
    call get_reals(buf,n,tmp(1:np))
    cell%nneigh(1:np) = nint(tmp(1:np))
    call get_reals(buf,n,one)
    cell%nneightry = nint(one(1))
 endif
 call get_reals(buf,n,tmp(1:np))
 cell%arr_index(1:np) = nint(tmp(1:np))

end subroutine unpack_celldens

!-----------------------------------------------------------------------
!+
!  append values to a buffer, advancing the position n
!+
!-----------------------------------------------------------------------
pure subroutine put_reals(buf,n,vals)
 real,    intent(inout) :: buf(:)
 integer, intent(inout) :: n
 real,    intent(in)    :: vals(:)

 buf(n+1:n+size(vals)) = vals
 n = n + size(vals)

end subroutine put_reals

!-----------------------------------------------------------------------
!+
!  read values from a buffer, advancing the position n
!+
!-----------------------------------------------------------------------
pure subroutine get_reals(buf,n,vals)
 real,    intent(in)    :: buf(:)
 integer, intent(inout) :: n
 real,    intent(out)   :: vals(:)

 vals = buf(n+1:n+size(vals))
 n = n + size(vals)

end subroutine get_reals

end module mpidens
