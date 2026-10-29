!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpiforce
!
! None
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, mpi, mpidens
!
 use io,       only:nprocs,fatal
 use dim,      only:minpart,maxfsum,maxxpartveciforce,maxxpartvecGR,gr
 use mpidens,  only:put_reals,get_reals

 implicit none
 private

 public :: cellforce
 public :: stackforce
 public :: get_mpitype_of_cellforce
 public :: free_mpitype_of_cellforce
 public :: pack_cellforce,unpack_cellforce

 !
 ! cells are sent as real buffers holding only what the other task needs:
 ! outbound (to the task with neighbours) the particle data, and on return
 ! (to the owner) the partial sums, both for the npcell particles only
 !
 integer, parameter, public :: maxpack_cellforce = max(8 + minpart*(maxxpartveciforce + 4), &
                                                       26 + minpart*(maxfsum + 4))
 ! particle values sent outbound: the metric values at the end of xpartvec only with GR
 integer, parameter :: nxpartvec_pack = maxxpartveciforce - merge(0,maxxpartvecGR,gr)

 integer, parameter :: ndata = 21 ! number of elements in the cell (including padding)
 integer, parameter :: nbytes_cellforce = 8 * maxxpartveciforce * minpart + &  !  xpartvec(maxxpartveciforce,minpart)
                                          8 * maxfsum * minpart           + &  !  fsums(maxfsum,minpart)
                                          8 * 20                          + &  !  fgrav(20)
                                          8 * 3                           + &  !  xpos(3)
                                          8                               + &  !  xsizei
                                          8                               + &  !  rcuti
                                          8 * minpart                     + &  !  tsmin(minpart)
                                          8 * minpart                     + &  !  vsigmax(minpart)
                                          4                               + &  !  icell
                                          4                               + &  !  npcell
                                          4 * minpart                     + &  !  arr_index(minpart)
                                          4 * minpart                     + &  !  ipos(minpart)
                                          4                               + &  !  ndrag
                                          4                               + &  !  nstokes
                                          4                               + &  !  nsuper
                                          4                               + &  !  owner
                                          4                               + &  !  waiting_index
                                          1 * minpart                     + &  !  iphase(minpart)
                                          1 * minpart                     + &  !  ibinneigh(minpart)
                                          1 * minpart                          !  apr_level

 type cellforce
    sequence
    real             :: xpartvec(maxxpartveciforce,minpart)
    real             :: fsums(maxfsum,minpart)
    real             :: fgrav(20)
    real             :: xpos(3)
    real             :: xsizei
    real             :: rcuti
    real             :: tsmin(minpart)
    real             :: vsigmax(minpart)
    integer          :: icell
    integer          :: npcell                                 ! number of particles in here
    integer          :: arr_index(minpart)
    integer          :: ipos(minpart)                          ! position of each particle in the owner's cell
    integer          :: ndrag
    integer          :: nstokes
    integer          :: nsuper
    integer          :: owner                                  ! id of the process that owns this
    integer          :: waiting_index
    integer(kind=1)  :: iphase(minpart)
    integer(kind=1)  :: ibinneigh(minpart)
    integer(kind=1)  :: apr(minpart)                           ! apr resolution level (not in xpartvec because integer)

    ! pad the array to 8-byte boundaries
    integer(kind=1)  :: pad(8 - mod(nbytes_cellforce, 8)) !padding to maintain alignment of elements
 end type cellforce

 type chunkforce
    sequence
    type(cellforce), pointer  :: cells(:) => null()
 end type chunkforce

 !
 ! stack of cells, stored in chunks of chunksize cells that are allocated as required,
 ! so the stack can grow from inside a parallel region without moving existing cells
 !
 type stackforce
    sequence
    type(chunkforce), pointer :: chunks(:) => null()
    integer                   :: chunksize = 0
    integer                   :: maxlength = 0   ! number of cells that can be stored without growing
    integer                   :: n = 0
    integer                   :: ipad = 0        ! pads the type to a multiple of 8 bytes (ifort error #6380)
 end type stackforce

contains

subroutine get_mpitype_of_cellforce(dtype)
#ifdef MPI
 use mpi
#endif
 integer, intent(out) :: dtype
#ifdef MPI
 integer                         :: nblock, blens(ndata), mpitypes(ndata)
 integer(kind=MPI_ADDRESS_KIND)  :: disp(ndata)

 type(cellforce)                 :: cell
 integer(kind=MPI_ADDRESS_KIND)  :: addr,start,lb,extent
 integer                         :: mpierr

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
 blens(nblock) = size(cell%tsmin)
 mpitypes(nblock) = MPI_REAL8
 call MPI_GET_ADDRESS(cell%tsmin,addr,mpierr)
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
 blens(nblock) = size(cell%ipos)
 mpitypes(nblock) = MPI_INTEGER4
 call MPI_GET_ADDRESS(cell%ipos,addr,mpierr)
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
 blens(nblock) = size(cell%iphase)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%iphase,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%ibinneigh)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%ibinneigh,addr,mpierr)
 disp(nblock) = addr - start

 nblock = nblock + 1
 blens(nblock) = size(cell%apr)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%apr,addr,mpierr)
 disp(nblock) = addr - start

 ! padding must come last
 nblock = nblock + 1
 blens(nblock) = 8 - mod(nbytes_cellforce, 8)
 mpitypes(nblock) = MPI_INTEGER1
 call MPI_GET_ADDRESS(cell%pad,addr,mpierr)
 disp(nblock) = addr - start

 call MPI_TYPE_CREATE_STRUCT(nblock,blens(1:nblock),disp(1:nblock),mpitypes(1:nblock),dtype,mpierr)
 call MPI_TYPE_COMMIT(dtype,mpierr)

 ! check extent okay
 call MPI_TYPE_GET_EXTENT(dtype,lb,extent,mpierr)
 if (extent /= sizeof(cell)) then
    call fatal('mpi_force','MPI_TYPE_GET_EXTENT has calculated the extent incorrectly')
 endif

#else
 dtype = 0
#endif

end subroutine get_mpitype_of_cellforce

subroutine free_mpitype_of_cellforce(dtype)
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
end subroutine free_mpitype_of_cellforce

!-----------------------------------------------------------------------
!+
!  pack a cell into a real buffer: header (npcell, owner, waiting_index),
!  then either the particle data (outbound) or the force sums (return)
!+
!-----------------------------------------------------------------------
subroutine pack_cellforce(cell,outbound,buf,n,mask)
 type(cellforce), intent(in)  :: cell
 logical,         intent(in)  :: outbound
 real,            intent(out) :: buf(:)
 integer,         intent(out) :: n
 logical,         intent(in), optional :: mask(:)
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
    idx(1:np) = cell%ipos(1:np)
 endif
 n  = 0
 call put_reals(buf,n,real((/np,cell%owner,cell%waiting_index/)))
 if (outbound) then
    call put_reals(buf,n,(/cell%xpos,cell%xsizei,cell%rcuti/))
    do k=1,np
       call put_reals(buf,n,cell%xpartvec(1:nxpartvec_pack,idx(k)))
    enddo
    call put_reals(buf,n,real(cell%iphase(idx(1:np))))
    call put_reals(buf,n,real(cell%apr(idx(1:np))))
    call put_reals(buf,n,real(cell%arr_index(idx(1:np))))
 else
    call put_reals(buf,n,reshape(cell%fsums(:,1:np),(/size(cell%fsums,1)*np/)))
    call put_reals(buf,n,cell%fgrav)
    call put_reals(buf,n,cell%tsmin(1:np))
    call put_reals(buf,n,cell%vsigmax(1:np))
    call put_reals(buf,n,real(cell%ibinneigh(1:np)))
    call put_reals(buf,n,real((/cell%ndrag,cell%nstokes,cell%nsuper/)))
 endif
 call put_reals(buf,n,real(idx(1:np)))

end subroutine pack_cellforce

!-----------------------------------------------------------------------
!+
!  unpack a cell packed by pack_cellforce; n returns the number of reals
!  read. Fields that are not sent are set for a cell computed remotely
!  (the force sums, timestep and fgrav are set by the remote calculation).
!+
!-----------------------------------------------------------------------
subroutine unpack_cellforce(buf,outbound,cell,n)
 real,            intent(in)    :: buf(:)
 logical,         intent(in)    :: outbound
 type(cellforce), intent(inout) :: cell
 integer,         intent(out)   :: n
 real    :: hdr(3),geom(5),tmp(minpart),counts(3)
 integer :: np,nx

 n = 0
 call get_reals(buf,n,hdr)
 np = nint(hdr(1))
 cell%npcell        = np
 cell%owner         = nint(hdr(2))
 cell%waiting_index = nint(hdr(3))
 if (outbound) then
    nx = nxpartvec_pack
    call get_reals(buf,n,geom)
    cell%xpos   = geom(1:3)
    cell%xsizei = geom(4)
    cell%rcuti  = geom(5)
    cell%xpartvec(1:nx,1:np) = reshape(buf(n+1:n+nx*np),(/nx,np/))
    n = n + nx*np
    call get_reals(buf,n,tmp(1:np))
    cell%iphase(1:np) = int(tmp(1:np),kind=1)
    call get_reals(buf,n,tmp(1:np))
    cell%apr(1:np) = int(tmp(1:np),kind=1)
    call get_reals(buf,n,tmp(1:np))
    cell%arr_index(1:np) = nint(tmp(1:np))  ! the same on every task for sinks only
    cell%ndrag   = 0
    cell%nstokes = 0
    cell%nsuper  = 0
 else
    nx = size(cell%fsums,1)
    cell%fsums(:,1:np) = reshape(buf(n+1:n+nx*np),(/nx,np/))
    n = n + nx*np
    call get_reals(buf,n,cell%fgrav)
    call get_reals(buf,n,cell%tsmin(1:np))
    call get_reals(buf,n,cell%vsigmax(1:np))
    call get_reals(buf,n,tmp(1:np))
    cell%ibinneigh(1:np) = int(tmp(1:np),kind=1)
    call get_reals(buf,n,counts)
    cell%ndrag   = nint(counts(1))
    cell%nstokes = nint(counts(2))
    cell%nsuper  = nint(counts(3))
 endif
 call get_reals(buf,n,tmp(1:np))
 cell%ipos(1:np) = nint(tmp(1:np))

end subroutine unpack_cellforce

end module mpiforce
