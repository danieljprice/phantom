!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpimemory
!
! Storage of the cells exchanged between MPI tasks during density and force.
! Cells sent by this task and waiting for results are kept as full cells on
! stacks; cells received from other tasks are kept packed (as sent) in a
! pool of reals, and unpacked one at a time when they are computed.
!
! :References: None
!
! :Owner: Conrad Chan
!
! :Runtime parameters: None
!
! :Dependencies: dim, io, mpidens, mpiforce
!
 use io,          only:fatal,iprint
 use mpidens,     only:celldens,stackdens,maxpack_celldens
 use mpiforce,    only:cellforce,stackforce,maxpack_cellforce

 implicit none

 interface allocate_stack
  module procedure allocate_stack_dens,allocate_stack_force
 end interface allocate_stack

 interface deallocate_stack
  module procedure deallocate_stack_dens,deallocate_stack_force
 end interface deallocate_stack

 interface swap_stacks ! force doesn't require a stack swap
  module procedure swap_stacks_dens
 end interface swap_stacks

 interface push_onto_stack
  module procedure push_onto_stack_dens,push_onto_stack_force
 end interface push_onto_stack

 interface get_cell
  module procedure get_cell_dens,get_cell_force
 end interface get_cell

 interface get_cell_ptr
  module procedure get_cell_ptr_dens,get_cell_ptr_force
 end interface get_cell_ptr

 interface write_cell
  module procedure write_cell_dens,write_cell_force
 end interface write_cell

 interface reserve_stack
  module procedure reserve_stack_dens,reserve_stack_force
 end interface reserve_stack

 interface get_packed_cell
  module procedure get_packed_celldens,get_packed_cellforce
 end interface get_packed_cell

 public :: allocate_mpi_memory
 public :: deallocate_mpi_memory
 public :: swap_stacks
 public :: push_onto_stack
 public :: get_cell
 public :: get_cell_ptr
 public :: write_cell
 public :: reserve_stack
 public :: reset_stacks
 public :: push_packed
 public :: get_packed_cell
 public :: clear_pool
 public :: exchange_memory_mb

 ! chunk of a pool of packed cells
 type poolchunk
    real, allocatable :: x(:)
 end type poolchunk

 ! cells received from other tasks, stored packed one after the other in chunks
 ! of reals (a cell never spans two chunks); pos(i) is the start of cell i,
 ! as (chunk - 1)*chunksize + offset
 type packedpool
    type(poolchunk), allocatable :: chunks(:)
    integer(kind=8), allocatable :: pos(:)
    integer         :: chunksize = 0
    integer         :: n = 0        ! number of cells in the pool
    integer         :: maxn = 0     ! largest number of cells held, for the MPI statistics
    integer         :: ichunk = 1   ! chunk being filled
    integer         :: maxchunks = 1 ! largest number of chunks used, for the MPI statistics
    integer         :: fill = 0     ! reals used in that chunk
 end type packedpool

 ! cells received from other tasks during density and force
 type(packedpool), public :: dens_pool, force_pool

 ! stacks to be referenced from density and force routines
 ! (stack 1 of each is kept for the generic interfaces, the received cells are in the pools)
 type(stackdens),  public :: dens_stack_1
 type(stackdens),  public :: dens_stack_2
 type(stackdens),  public :: dens_stack_3
 type(stackforce), public :: force_stack_1
 type(stackforce), public :: force_stack_2

 ! number of cells in each chunk of a stack (the initial size of the stack)
 integer, public :: stacksize

 private

 ! maximum number of chunks per stack
 integer :: maxchunks = 0

contains

!----------------------------------------------------------------
!+
!  allocate the stacks used to exchange cells between MPI tasks
!+
!----------------------------------------------------------------
subroutine allocate_mpi_memory(npart, stacksize_in)
 integer(kind=8), intent(in), optional :: npart
 integer,         intent(in), optional :: stacksize_in

 if (present(stacksize_in)) stacksize = stacksize_in
 if (present(npart)) call calculate_stacksize(npart)

 call allocate_stack(dens_stack_2)
 call allocate_stack(dens_stack_3)
 call allocate_stack(force_stack_2)
 call allocate_pool(dens_pool,maxpack_celldens)
 call allocate_pool(force_pool,maxpack_cellforce)

end subroutine allocate_mpi_memory

subroutine calculate_stacksize(npart)
 use dim, only:mpi,minpart
 use io,  only:nprocs,id,master
 integer(kind=8), intent(in) :: npart ! total number of particles on all tasks
 integer, parameter  :: safety = 2

 ! size of the stack needed for communication,
 ! should be at least the maximum number of cells that need
 ! to be exported to other tasks.
 !
 ! if it is not large enough, it will be automatically expanded

 ! number of particles per cell, divided by number of tasks
 if (mpi .and. nprocs > 1) then
    ! the stacks hold the cells this task sends, so at most its own cells (leaves hold
    ! up to minpart particles), with a safety factor; they grow if more are needed
    stacksize = int((npart / minpart / nprocs) * safety)

    ! enough chunks to hold one cell per particle, an upper limit
    ! on the number of cells that can be imported or exported
    maxchunks = int(npart / max(stacksize,1)) + 2

    if (id == master) then
       write(iprint, *) 'MPI memory stack size = ', stacksize
       write(iprint, *) '  (total number of cells that can be exported by a single task)'
    endif
 else
    stacksize = 0
 endif

end subroutine calculate_stacksize

subroutine deallocate_mpi_memory

 call deallocate_stack(dens_stack_2)
 call deallocate_stack(dens_stack_3)
 call deallocate_stack(force_stack_2)
 call deallocate_pool(dens_pool)
 call deallocate_pool(force_pool)

end subroutine deallocate_mpi_memory

!----------------------------------------------------------------
!+
!  allocate the first chunk of a stack, with stacksize cells
!+
!----------------------------------------------------------------
subroutine allocate_stack_dens(stack)
 type(stackdens), intent(inout) :: stack
 integer :: allocstat

 if (associated(stack%chunks)) then
    if (max(stacksize,1) /= stack%chunksize) call fatal('stack','dens stack already allocated with a different size')
    return
 endif
 stack%chunksize = max(stacksize,1)
 stack%maxlength = stack%chunksize
 allocate(stack%chunks(max(maxchunks,1000)),stat=allocstat)
 if (allocstat == 0) allocate(stack%chunks(1)%cells(stack%chunksize),stat=allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 stack%n = 0

end subroutine allocate_stack_dens

subroutine allocate_stack_force(stack)
 type(stackforce), intent(inout) :: stack
 integer :: allocstat

 if (associated(stack%chunks)) then
    if (max(stacksize,1) /= stack%chunksize) call fatal('stack','force stack already allocated with a different size')
    return
 endif
 stack%chunksize = max(stacksize,1)
 stack%maxlength = stack%chunksize
 allocate(stack%chunks(max(maxchunks,1000)),stat=allocstat)
 if (allocstat == 0) allocate(stack%chunks(1)%cells(stack%chunksize),stat=allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 stack%n = 0

end subroutine allocate_stack_force

subroutine deallocate_stack_dens(stack)
 type(stackdens), intent(inout) :: stack
 integer :: k

 if (.not.associated(stack%chunks)) return
 do k=1,size(stack%chunks)
    if (associated(stack%chunks(k)%cells)) deallocate(stack%chunks(k)%cells)
 enddo
 deallocate(stack%chunks)
 stack%chunksize = 0
 stack%maxlength = 0
 stack%n = 0

end subroutine deallocate_stack_dens

subroutine deallocate_stack_force(stack)
 type(stackforce), intent(inout) :: stack
 integer :: k

 if (.not.associated(stack%chunks)) return
 do k=1,size(stack%chunks)
    if (associated(stack%chunks(k)%cells)) deallocate(stack%chunks(k)%cells)
 enddo
 deallocate(stack%chunks)
 stack%chunksize = 0
 stack%maxlength = 0
 stack%n = 0

end subroutine deallocate_stack_force

!----------------------------------------------------------------
!+
!  allocate the chunk that holds cell i, if not already allocated
!  (threadsafe: called from inside parallel regions)
!+
!----------------------------------------------------------------
subroutine grow_stack_dens(stack,i)
 use io, only:id
 type(stackdens), intent(inout) :: stack
 integer,         intent(in)    :: i
 integer :: k,allocstat

 k = (i-1)/stack%chunksize + 1
 if (k > size(stack%chunks)) call fatal('dens','MPI stack exceeded maximum size')
 allocstat = 0
 !$omp critical (grow_mpi_stack)
 if (.not.associated(stack%chunks(k)%cells)) then
    allocate(stack%chunks(k)%cells(stack%chunksize),stat=allocstat)
    stack%maxlength = max(stack%maxlength,k*stack%chunksize)
    write(iprint,*) 'MPI dens stack on',id,'increased to',stack%maxlength,'cells'
 endif
 !$omp end critical (grow_mpi_stack)
 if (allocstat /= 0) call fatal('dens','error increasing MPI stack size')

end subroutine grow_stack_dens

subroutine grow_stack_force(stack,i)
 use io, only:id
 type(stackforce), intent(inout) :: stack
 integer,          intent(in)    :: i
 integer :: k,allocstat

 k = (i-1)/stack%chunksize + 1
 if (k > size(stack%chunks)) call fatal('force','MPI stack exceeded maximum size')
 allocstat = 0
 !$omp critical (grow_mpi_stack)
 if (.not.associated(stack%chunks(k)%cells)) then
    allocate(stack%chunks(k)%cells(stack%chunksize),stat=allocstat)
    stack%maxlength = max(stack%maxlength,k*stack%chunksize)
    write(iprint,*) 'MPI force stack on',id,'increased to',stack%maxlength,'cells'
 endif
 !$omp end critical (grow_mpi_stack)
 if (allocstat /= 0) call fatal('force','error increasing MPI stack size')

end subroutine grow_stack_force

subroutine swap_stacks_dens(stack_a, stack_b)
 type(stackdens), intent(inout) :: stack_a
 type(stackdens), intent(inout) :: stack_b
 type(stackdens) :: temp

 temp    = stack_a
 stack_a = stack_b
 stack_b = temp

end subroutine swap_stacks_dens

subroutine push_onto_stack_dens(stack,cell)
 type(stackdens), intent(inout) :: stack
 type(celldens),  intent(in)    :: cell
 type(celldens),  pointer       :: p
 integer :: i

 call reserve_stack(stack,i)

 ! no other thread will write to the same position, so it is threadsafe to write without a critical section
 p => get_cell_ptr(stack,i)
 p = cell

end subroutine push_onto_stack_dens

subroutine push_onto_stack_force(stack,cell)
 type(stackforce), intent(inout) :: stack
 type(cellforce),  intent(in)    :: cell
 type(cellforce),  pointer       :: p
 integer :: i

 call reserve_stack(stack,i)

 ! no other thread will write to the same position, so it is threadsafe to write without a critical section
 p => get_cell_ptr(stack,i)
 p = cell

end subroutine push_onto_stack_force

!----------------------------------------------------------------
!+
!  pointer to cell i in the stack
!+
!----------------------------------------------------------------
function get_cell_ptr_dens(stack,i) result(p)
 type(stackdens), intent(in) :: stack
 integer,         intent(in) :: i
 type(celldens),  pointer    :: p
 integer :: k

 k = (i-1)/stack%chunksize
 p => stack%chunks(k+1)%cells(i - k*stack%chunksize)

end function get_cell_ptr_dens

function get_cell_ptr_force(stack,i) result(p)
 type(stackforce), intent(in) :: stack
 integer,          intent(in) :: i
 type(cellforce),  pointer    :: p
 integer :: k

 k = (i-1)/stack%chunksize
 p => stack%chunks(k+1)%cells(i - k*stack%chunksize)

end function get_cell_ptr_force

type(celldens) function get_cell_dens(stack,i)
 type(stackdens), intent(in) :: stack
 integer,         intent(in) :: i
 type(celldens),  pointer    :: p

 if (stack%n < i) call fatal('dens','attempting to read invalid stack address')
 p => get_cell_ptr(stack,i)
 get_cell_dens = p

end function get_cell_dens

type(cellforce) function get_cell_force(stack,i)
 type(stackforce), intent(in) :: stack
 integer,          intent(in) :: i
 type(cellforce),  pointer    :: p

 if (stack%n < i) call fatal('force','attempting to read invalid stack address')
 p => get_cell_ptr(stack,i)
 get_cell_force = p

end function get_cell_force

subroutine write_cell_dens(stack,cell)
 type(stackdens), intent(inout) :: stack
 type(celldens),  intent(inout) :: cell
 type(celldens),  pointer       :: p

 if (cell%waiting_index > stack%n) call fatal('dens','attempting to write to invalid stack address')
 p => get_cell_ptr(stack,cell%waiting_index)
 p = cell

end subroutine write_cell_dens

subroutine write_cell_force(stack,cell)
 type(stackforce), intent(inout) :: stack
 type(cellforce),  intent(inout) :: cell
 type(cellforce),  pointer       :: p

 if (cell%waiting_index > stack%n) call fatal('force','attempting to write to invalid stack address')
 p => get_cell_ptr(stack,cell%waiting_index)
 p = cell

end subroutine write_cell_force

!----------------------------------------------------------------
!+
!  reserve a position on the stack, growing it if required
!+
!----------------------------------------------------------------
subroutine reserve_stack_dens(stack,i)
 type(stackdens), intent(inout) :: stack
 integer,         intent(out)   :: i

 !$omp atomic capture
 stack%n = stack%n + 1
 i = stack%n
 !$omp end atomic

 ! cells beyond the first chunk always go through the critical section in
 ! grow_stack, which ensures their chunk is allocated and visible to this thread
 if (i > stack%chunksize) call grow_stack_dens(stack,i)

end subroutine reserve_stack_dens

subroutine reserve_stack_force(stack,i)
 type(stackforce), intent(inout) :: stack
 integer,          intent(out)   :: i

 !$omp atomic capture
 stack%n = stack%n + 1
 i = stack%n
 !$omp end atomic

 ! cells beyond the first chunk always go through the critical section in
 ! grow_stack, which ensures their chunk is allocated and visible to this thread
 if (i > stack%chunksize) call grow_stack_force(stack,i)

end subroutine reserve_stack_force

subroutine reset_stacks
 dens_stack_1%n=0
 dens_stack_2%n=0
 dens_stack_3%n=0

 force_stack_1%n=0
 force_stack_2%n=0

 ! release the memory of the pools above their first chunk, so that one
 ! step with many imported cells does not hold memory for the rest of the run
 call clear_pool(dens_pool,release=.true.)
 call clear_pool(force_pool,release=.true.)
end subroutine reset_stacks

!----------------------------------------------------------------
!+
!  largest memory (in MB) used by the stacks and pools so far on this task
!+
!----------------------------------------------------------------
real function exchange_memory_mb()
 type(celldens)  :: cd
 type(cellforce) :: cf

 exchange_memory_mb = 1.e-6*(storage_size(cd)/8*real(dens_stack_2%maxlength + dens_stack_3%maxlength) &
                    + storage_size(cf)/8*real(force_stack_2%maxlength) &
                    + storage_size(1.)/8*(real(dens_pool%maxchunks)*dens_pool%chunksize &
                                        + real(force_pool%maxchunks)*force_pool%chunksize))

end function exchange_memory_mb

!----------------------------------------------------------------
!+
!  allocate the first chunk of a pool of packed cells, large enough
!  for a few times the expected number of cells of maxpack reals
!+
!----------------------------------------------------------------
subroutine allocate_pool(pool,maxpack)
 type(packedpool), intent(inout) :: pool
 integer,          intent(in)    :: maxpack
 integer :: allocstat
 integer, parameter :: maxpoolchunks = 100000

 if (allocated(pool%chunks)) return
 ! packed cells are typically a quarter of their maximum size
 pool%chunksize = int(min(max(stacksize,1000)*int(maxpack,kind=8)/4,2_8**30))
 allocate(pool%chunks(maxpoolchunks),pool%pos(max(stacksize,1000)),stat=allocstat)
 if (allocstat == 0) allocate(pool%chunks(1)%x(pool%chunksize),stat=allocstat)
 if (allocstat /= 0) call fatal('pool','memory allocation error')
 pool%n = 0
 pool%ichunk = 1
 pool%fill = 0

end subroutine allocate_pool

subroutine deallocate_pool(pool)
 type(packedpool), intent(inout) :: pool

 if (allocated(pool%chunks)) deallocate(pool%chunks)
 if (allocated(pool%pos)) deallocate(pool%pos)
 pool%n = 0

end subroutine deallocate_pool

!----------------------------------------------------------------
!+
!  empty a pool, optionally freeing all but its first chunk
!+
!----------------------------------------------------------------
subroutine clear_pool(pool,release)
 type(packedpool), intent(inout) :: pool
 logical,          intent(in)    :: release
 integer :: k

 pool%n = 0
 pool%ichunk = 1
 pool%fill = 0
 if (release .and. allocated(pool%chunks)) then
    do k=2,size(pool%chunks)
       if (allocated(pool%chunks(k)%x)) deallocate(pool%chunks(k)%x)
    enddo
 endif

end subroutine clear_pool

!----------------------------------------------------------------
!+
!  add a packed cell (buf(1:nbuf), as received) to a pool
!  (threadsafe: space is reserved in a critical section, then copied)
!+
!----------------------------------------------------------------
subroutine push_packed(pool,buf,nbuf)
 type(packedpool), intent(inout) :: pool
 real,             intent(in)    :: buf(:)
 integer,          intent(in)    :: nbuf
 integer(kind=8), allocatable :: tmp(:)
 integer :: i,k,ioff,allocstat

 if (nbuf > pool%chunksize) call fatal('pool','packed cell larger than a chunk')
 allocstat = 0
 !$omp critical (packed_pool)
 if (pool%fill + nbuf > pool%chunksize) then  ! start the next chunk
    pool%ichunk = pool%ichunk + 1
    pool%fill = 0
 endif
 k = pool%ichunk
 if (k > size(pool%chunks)) call fatal('pool','MPI pool exceeded maximum size')
 pool%maxchunks = max(pool%maxchunks,k)
 if (.not.allocated(pool%chunks(k)%x)) allocate(pool%chunks(k)%x(pool%chunksize),stat=allocstat)
 ioff = pool%fill
 pool%fill = pool%fill + nbuf
 i = pool%n + 1
 if (i > size(pool%pos)) then
    allocate(tmp(2*size(pool%pos)))
    tmp(1:size(pool%pos)) = pool%pos
    call move_alloc(tmp,pool%pos)
 endif
 pool%pos(i) = int(k-1,kind=8)*pool%chunksize + ioff
 pool%n = i
 pool%maxn = max(pool%maxn,i)
 !$omp end critical (packed_pool)
 if (allocstat /= 0) call fatal('pool','error increasing MPI pool size')

 pool%chunks(k)%x(ioff+1:ioff+nbuf) = buf(1:nbuf)

end subroutine push_packed

!----------------------------------------------------------------
!+
!  unpack cell i of a pool of received cells
!+
!----------------------------------------------------------------
subroutine get_packed_celldens(pool,i,cell)
 use mpidens, only:unpack_celldens
 type(packedpool), intent(in)    :: pool
 integer,          intent(in)    :: i
 type(celldens),   intent(inout) :: cell
 integer :: k,ioff,nbuf

 if (i > pool%n) call fatal('dens','attempting to read invalid pool address')
 k    = int(pool%pos(i)/pool%chunksize) + 1
 ioff = int(pool%pos(i) - int(k-1,kind=8)*pool%chunksize)
 call unpack_celldens(pool%chunks(k)%x(ioff+1:),.true.,cell,nbuf)

end subroutine get_packed_celldens

subroutine get_packed_cellforce(pool,i,cell)
 use mpiforce, only:unpack_cellforce
 type(packedpool), intent(in)    :: pool
 integer,          intent(in)    :: i
 type(cellforce),  intent(inout) :: cell
 integer :: k,ioff,nbuf

 if (i > pool%n) call fatal('force','attempting to read invalid pool address')
 k    = int(pool%pos(i)/pool%chunksize) + 1
 ioff = int(pool%pos(i) - int(k-1,kind=8)*pool%chunksize)
 call unpack_cellforce(pool%chunks(k)%x(ioff+1:),.true.,cell,nbuf)

end subroutine get_packed_cellforce

end module mpimemory
