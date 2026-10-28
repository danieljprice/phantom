!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module mpimemory
!
! None
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
 use mpidens,     only:celldens,stackdens
 use mpiforce,    only:cellforce,stackforce

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

 public :: allocate_mpi_memory
 public :: deallocate_mpi_memory
 public :: swap_stacks
 public :: push_onto_stack
 public :: get_cell
 public :: get_cell_ptr
 public :: write_cell
 public :: reserve_stack
 public :: reset_stacks

 ! stacks to be referenced from density and force routines
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

 call allocate_stack(dens_stack_1)
 call allocate_stack(dens_stack_2)
 call allocate_stack(dens_stack_3)
 call allocate_stack(force_stack_1)
 call allocate_stack(force_stack_2)

end subroutine allocate_mpi_memory

subroutine calculate_stacksize(npart)
 use dim, only:mpi,minpart
 use io,  only:nprocs,id,master
 integer(kind=8), intent(in) :: npart ! total number of particles on all tasks
 integer, parameter  :: safety = 8

 ! size of the stack needed for communication,
 ! should be at least the maximum number of cells that need
 ! to be exported to other tasks.
 !
 ! if it is not large enough, it will be automatically expanded

 ! number of particles per cell, divided by number of tasks
 if (mpi .and. nprocs > 1) then
    ! assume that every cell will be exported, with some safety factor
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

 call deallocate_stack(dens_stack_1)
 call deallocate_stack(dens_stack_2)
 call deallocate_stack(dens_stack_3)
 call deallocate_stack(force_stack_1)
 call deallocate_stack(force_stack_2)

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
end subroutine reset_stacks

end module mpimemory
