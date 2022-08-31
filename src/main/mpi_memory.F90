!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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

 interface swap_stacks ! force doesn't require a stack swap
  module procedure swap_stacks_dens
 end interface swap_stacks

 interface push_onto_stack
  module procedure push_onto_stack_dens,push_onto_stack_force
 end interface push_onto_stack

 interface pop_off_stack
  module procedure pop_off_stack_dens,pop_off_stack_force
 end interface pop_off_stack

 interface get_cell
  module procedure get_cell_force
 end interface get_cell

 interface write_cell
  module procedure write_cell_force
 end interface write_cell

 interface reserve_stack
  module procedure reserve_stack_dens,reserve_stack_force
 end interface reserve_stack

 public :: allocate_mpi_memory
 public :: deallocate_mpi_memory
 public :: allocate_stack
 public :: swap_stacks
 public :: push_onto_stack
 public :: pop_off_stack
 public :: get_cell
 public :: write_cell
 public :: reserve_stack
 public :: reset_stacks

 ! stacks to be referenced from density and force routines
 type(stackdens),  public :: dens_stack_1
 type(stackdens),  public :: dens_stack_2
 type(stackdens),  public :: dens_stack_3
 type(stackforce), public :: force_stack_1
 type(stackforce), public :: force_stack_2

 private

 integer :: stacksize_dens
 integer :: stacksize_force
 !$omp threadprivate(stacksize_force)

 ! primary chunk of memory requested using alloc
 type(celldens),  allocatable, target :: dens_cells(:,:)
 type(cellforce), allocatable, target :: force_cells(:,:)

contains

subroutine allocate_mpi_memory(npart, stacksize_in, reallocate)
 integer, optional,  intent(in) :: npart
 integer, optional,  intent(in) :: stacksize_in
 logical, optional,  intent(in) :: reallocate
 integer :: allocstat
 logical :: re_allocate = .false.

 allocstat = 0

 if (present(stacksize_in)) stacksize_dens = stacksize_in
 if (present(stacksize_in)) stacksize_force = stacksize_in
 if (present(npart)) call calculate_stacksize(npart)
 if (present(reallocate)) re_allocate = reallocate

 if (.not. allocated(dens_cells)) allocate(dens_cells(stacksize_dens,3), stat=allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')

 ! If reallocating an existing stack for expanding MPI memory,
 ! give the stack the same address that it previously had. This
 ! may not be the same as the initial order because density stacks
 ! can be swapped. Force stacks do not get swapped.
 if (re_allocate) then
    call allocate_stack(dens_stack_1, dens_stack_1%number)
    call allocate_stack(dens_stack_2, dens_stack_2%number)
    call allocate_stack(dens_stack_3, dens_stack_3%number)
 else
    call allocate_stack(dens_stack_1, 1)
    call allocate_stack(dens_stack_2, 2)
    call allocate_stack(dens_stack_3, 3)
 endif

 if (.not. allocated(force_cells)) allocate(force_cells(stacksize_force,2), stat=allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 call allocate_stack(force_stack_1, 1)
 call allocate_stack(force_stack_2, 2)

end subroutine allocate_mpi_memory

subroutine increase_mpi_memory_dens
 use io, only:id
 real, parameter :: factor = 1.5
 integer         :: stacksize_new
 integer         :: allocstat

 ! temporary memory for increasing stack sizes
 type(celldens),  allocatable, target :: dens_cells_tmp(:,:)

 stacksize_new = int(real(stacksize_dens) * factor)
 write(iprint, *) 'MPI stack exceeded on', id, 'increasing size to', stacksize_new

 ! Expand density
 call move_alloc(dens_cells, dens_cells_tmp)
 allocate(dens_cells(stacksize_new,3), stat=allocstat)
 if (allocstat /= 0) call fatal('stack', 'error increasing dens stack size')
 dens_cells(1:stacksize_dens,:) = dens_cells_tmp(:,:)
 deallocate(dens_cells_tmp)

 stacksize_dens = stacksize_new
 call allocate_stack(dens_stack_1, dens_stack_1%number)
 call allocate_stack(dens_stack_2, dens_stack_2%number)
 call allocate_stack(dens_stack_3, dens_stack_3%number)

end subroutine increase_mpi_memory_dens

subroutine increase_mpi_memory_force
 use io, only:id
 real, parameter :: factor = 1.5
 integer         :: stacksize_new
 integer         :: allocstat

 ! temporary memory for increasing stack sizes
 type(cellforce), allocatable, target :: force_cells_tmp(:,:)

 stacksize_new = int(real(stacksize_force) * factor)
 write(iprint, *) 'MPI stack exceeded on', id, 'increasing size to', stacksize_new

 ! Expand force
 call move_alloc(force_cells, force_cells_tmp)
 allocate(force_cells(stacksize_new,2), stat=allocstat)
 if (allocstat /= 0) call fatal('stack', 'error increasing force stack size')
 force_cells(1:stacksize_force,:) = force_cells_tmp(:,:)
 deallocate(force_cells_tmp)

 ! Set new stacksize value
 ! Reallocate, with memory already containing cells
!  call allocate_mpi_memory(stacksize_in=stacksize_new, reallocate=.true.)
 ! TEMP: while only force has parallel recv, it's easier to just call the allocate here
 ! rather than using allocate_mpi_memory
 stacksize_force = stacksize_new
 call allocate_stack(force_stack_1, 1)
 call allocate_stack(force_stack_2, 2)

end subroutine increase_mpi_memory_force

subroutine calculate_stacksize(npart)
 use dim, only:mpi,minpart
 use io,  only:nprocs,id,master
 integer, intent(in) :: npart
 integer, parameter  :: safety = 4

 ! size of the stack needed for communication,
 ! should be at least the maximum number of cells that need
 ! to be exported to other tasks.
 !
 ! if it is not large enough, it will be automatically expanded

 ! number of particles per cell, divided by number of tasks
 if (mpi .and. nprocs > 1) then
    ! assume that every cell will be exported, with some safety factor
    stacksize_dens = (npart / minpart / nprocs) * safety
    stacksize_force = stacksize_dens

    if (id == master) then
       write(iprint, *) 'MPI memory stack size = ', stacksize_dens
       write(iprint, *) '  (total number of cells that can be exported by a single task)'
    endif
 else
    stacksize_dens = 0
    stacksize_force = 0
 endif

end subroutine calculate_stacksize

subroutine deallocate_mpi_memory
 if (allocated(dens_cells )) deallocate(dens_cells )
 if (allocated(force_cells)) deallocate(force_cells)
end subroutine deallocate_mpi_memory

subroutine allocate_stack_dens(stack, i)
 type(stackdens), intent(inout) :: stack
 integer,         intent(in)    :: i

 stack%number = i
 stack%cells => dens_cells(1:stacksize_dens,stack%number)
 stack%maxlength = stacksize_dens

end subroutine allocate_stack_dens

subroutine allocate_stack_force(stack, i)
 type(stackforce),   intent(inout) :: stack
 integer,            intent(in)    :: i

 stack%number = i
 stack%cells => force_cells(1:stacksize_force,stack%number)
 stack%maxlength = stacksize_force

end subroutine allocate_stack_force

subroutine swap_stacks_dens(stack_a, stack_b)
 type(stackdens),   intent(inout) :: stack_a
 type(stackdens),   intent(inout) :: stack_b

 integer :: temp_n
 integer :: temp_number

 if (stack_a%maxlength /= stack_b%maxlength) call fatal('stack', 'stack swap of unequal size')

 ! counters
 temp_n = stack_a%n
 stack_a%n = stack_b%n
 stack_b%n = temp_n

 ! addresses
 temp_number = stack_a%number
 stack_a%number = stack_b%number
 stack_b%number = temp_number

 ! change pointers
 stack_a%cells => dens_cells(1:stacksize_dens,stack_a%number)
 stack_b%cells => dens_cells(1:stacksize_dens,stack_b%number)

end subroutine swap_stacks_dens

subroutine push_onto_stack_dens(stack,cell)
 type(stackdens),    intent(inout)  :: stack
 type(celldens),     intent(in)     :: cell

 if (stack%n + 1 > stack%maxlength) call increase_mpi_memory_dens
 ! after increasing stack size, cells can be added to because it is just a pointer
 stack%n = stack%n + 1
 stack%cells(stack%n) = cell
end subroutine push_onto_stack_dens

subroutine push_onto_stack_force(stack,cell)
 type(stackforce),   intent(inout)  :: stack
 type(cellforce),    intent(in)     :: cell

 if (stack%n + 1 > stack%maxlength) call increase_mpi_memory_force
 ! after increasing stack size, cells can be added to because it is just a pointer
 stack%n = stack%n + 1
 stack%cells(stack%n) = cell
end subroutine push_onto_stack_force

subroutine pop_off_stack_dens(stack,cell)
 type(stackdens),    intent(inout)  :: stack
 type(celldens),     intent(out)    :: cell

 if (stack%n <= 0) call fatal('density','attempting to pop empty stack')
 cell = stack%cells(stack%n)
 stack%n = stack%n - 1
end subroutine pop_off_stack_dens

subroutine pop_off_stack_force(stack,cell)
 type(stackforce),   intent(inout)  :: stack
 type(cellforce),    intent(out)    :: cell

 if (stack%n <= 0) call fatal('force','attempting to pop empty stack')
 cell = stack%cells(stack%n)
 stack%n = stack%n - 1
end subroutine pop_off_stack_force

type(cellforce) function get_cell_force(stack,i)
 type(stackforce),   intent(in)  :: stack
 integer,            intent(in)  :: i

 if (stack%n < i) call fatal('force','attempting to read invalid stack address')
 get_cell_force = stack%cells(i)
end function get_cell_force

subroutine write_cell_force(stack,cell)
 type(stackforce),   intent(inout)  :: stack
 type(cellforce),    intent(inout)  :: cell

 !$omp critical (crit_stack)
 if (cell%waiting_index > stack%maxlength) call fatal('force','attempting to write to invalid stack address')
 stack%cells(cell%waiting_index) = cell
 !$omp end critical (crit_stack)

end subroutine write_cell_force

subroutine reserve_stack_dens(stack,i)
 type(stackdens),    intent(inout) :: stack
 integer,            intent(out)   :: i

 if (stack%n + 1 > stack%maxlength) call increase_mpi_memory_dens
 stack%n = stack%n + 1
 i = stack%n

end subroutine reserve_stack_dens

subroutine reserve_stack_force(stack,i)
 type(stackforce),   intent(inout) :: stack
 integer,            intent(out)   :: i

 !$omp critical (crit_stack)
 if (stack%n + 1 > stack%maxlength) call increase_mpi_memory_force
 stack%n = stack%n + 1
 i = stack%n
 !$omp end critical (crit_stack)

end subroutine reserve_stack_force

subroutine reset_stacks
 dens_stack_1%n=0
 dens_stack_2%n=0
 dens_stack_3%n=0

 force_stack_1%n=0
 force_stack_2%n=0
end subroutine reset_stacks

end module mpimemory
