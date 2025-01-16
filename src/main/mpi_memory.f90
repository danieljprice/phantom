!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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

 interface swap_stacks ! force doesn't require a stack swap
  module procedure swap_stacks_dens
 end interface swap_stacks

 interface push_onto_stack
  module procedure push_onto_stack_dens,push_onto_stack_force
 end interface push_onto_stack

 interface get_cell
  module procedure get_cell_dens,get_cell_force
 end interface get_cell

 interface write_cell
  module procedure write_cell_dens,write_cell_force
 end interface write_cell

 interface reserve_stack
  module procedure reserve_stack_dens,reserve_stack_force
 end interface reserve_stack

 public :: allocate_mpi_memory
 public :: deallocate_mpi_memory
 public :: allocate_stack
 public :: swap_stacks
 public :: push_onto_stack
 public :: get_cell
 public :: write_cell
 public :: reserve_stack
 public :: reset_stacks
 public :: increase_mpi_memory

 ! stacks to be referenced from density and force routines
 type(stackdens),  public :: dens_stack_1
 type(stackdens),  public :: dens_stack_2
 type(stackdens),  public :: dens_stack_3
 type(stackforce), public :: force_stack_1
 type(stackforce), public :: force_stack_2

 integer, public :: stacksize

 private

 ! primary chunk of memory requested using alloc
 type(celldens),  allocatable, target :: dens_cells(:,:)
 type(cellforce), allocatable, target :: force_cells(:,:)

contains

subroutine allocate_mpi_memory(npart, stacksize_in)
 integer, optional,  intent(in) :: npart
 integer, optional,  intent(in) :: stacksize_in
 integer :: allocstat

 allocstat = 0

 if (present(stacksize_in)) stacksize = stacksize_in
 if (present(npart)) call calculate_stacksize(npart)

 if (allocated(dens_cells)) then
    if (stacksize /= size(dens_cells,1)) then
       call fatal('stack', 'dens_cells already allocated with a different size')
    endif
 endif

 if (allocated(force_cells)) then
    if (stacksize /= size(force_cells,1)) then
       call fatal('stack', 'force_cells already allocated with a different size')
    endif
 endif

 if (.not. allocated(dens_cells)) allocate(dens_cells(stacksize,3), stat=allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 call allocate_stack(dens_stack_1, 1)
 call allocate_stack(dens_stack_2, 2)
 call allocate_stack(dens_stack_3, 3)

 if (.not. allocated(force_cells)) allocate(force_cells(stacksize,2), stat=allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 call allocate_stack(force_stack_1, 1)
 call allocate_stack(force_stack_2, 2)

end subroutine allocate_mpi_memory

subroutine increase_mpi_memory
 use io, only:id
 real, parameter :: factor = 1.5
 integer         :: stacksize_new
 integer         :: allocstat

 ! temporary memory for increasing stack sizes
 type(celldens),  allocatable, target :: dens_cells_tmp(:,:)
 type(cellforce), allocatable, target :: force_cells_tmp(:,:)

 stacksize_new = int(real(stacksize) * factor)
 write(iprint, *) 'MPI dens stack exceeded on', id, 'increasing size to', stacksize_new

 ! Expand density
 call move_alloc(dens_cells, dens_cells_tmp)
 allocate(dens_cells(stacksize_new,3), stat=allocstat)
 if (allocstat /= 0) call fatal('stack', 'error increasing dens stack size')
 dens_cells(1:stacksize,:) = dens_cells_tmp(:,:)
 deallocate(dens_cells_tmp)

 ! Expand force
 call move_alloc(force_cells, force_cells_tmp)
 allocate(force_cells(stacksize_new,2), stat=allocstat)
 if (allocstat /= 0) call fatal('stack', 'error increasing force stack size')
 force_cells(1:stacksize,:) = force_cells_tmp(:,:)
 deallocate(force_cells_tmp)

 stacksize = stacksize_new
 call allocate_stack(force_stack_1, 1)
 call allocate_stack(force_stack_2, 2)
 call allocate_stack(dens_stack_1, dens_stack_1%number)
 call allocate_stack(dens_stack_2, dens_stack_2%number)
 call allocate_stack(dens_stack_3, dens_stack_3%number)
end subroutine increase_mpi_memory

subroutine calculate_stacksize(npart)
 use dim, only:mpi,minpart
 use io,  only:nprocs,id,master
 integer, intent(in) :: npart
 integer, parameter  :: safety = 8

 ! size of the stack needed for communication,
 ! should be at least the maximum number of cells that need
 ! to be exported to other tasks.
 !
 ! if it is not large enough, it will be automatically expanded

 ! number of particles per cell, divided by number of tasks
 if (mpi .and. nprocs > 1) then
    ! assume that every cell will be exported, with some safety factor
    stacksize = (npart / minpart / nprocs) * safety

    if (id == master) then
       write(iprint, *) 'MPI memory stack size = ', stacksize
       write(iprint, *) '  (total number of cells that can be exported by a single task)'
    endif
 else
    stacksize = 0
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
 stack%cells => dens_cells(1:stacksize,stack%number)
 stack%maxlength = stacksize

end subroutine allocate_stack_dens

subroutine allocate_stack_force(stack, i)
 type(stackforce),   intent(inout) :: stack
 integer,            intent(in)    :: i

 stack%number = i
 stack%cells => force_cells(1:stacksize,stack%number)
 stack%maxlength = stacksize

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
 stack_a%cells => dens_cells(1:stacksize,stack_a%number)
 stack_b%cells => dens_cells(1:stacksize,stack_b%number)

end subroutine swap_stacks_dens

subroutine push_onto_stack_dens(stack,cell)
 type(stackdens),    intent(inout)  :: stack
 type(celldens),     intent(in)     :: cell

 integer :: i

 call reserve_stack(stack,i)

 ! no other thread will write to the same position, so it is threadsafe to write without a critical section
 stack%cells(i) = cell
end subroutine push_onto_stack_dens

subroutine push_onto_stack_force(stack,cell)
 type(stackforce),   intent(inout)  :: stack
 type(cellforce),    intent(in)     :: cell

 integer :: i

 call reserve_stack(stack,i)

 ! no other thread will write to the same position, so it is threadsafe to write without a critical section
 stack%cells(i) = cell
end subroutine push_onto_stack_force

type(celldens) function get_cell_dens(stack,i)
 type(stackdens),    intent(in)  :: stack
 integer,            intent(in)  :: i

 if (stack%n < i) call fatal('dens','attempting to read invalid stack address')
 get_cell_dens = stack%cells(i)
end function get_cell_dens

type(cellforce) function get_cell_force(stack,i)
 type(stackforce),   intent(in)  :: stack
 integer,            intent(in)  :: i

 if (stack%n < i) call fatal('force','attempting to read invalid stack address')
 get_cell_force = stack%cells(i)
end function get_cell_force

subroutine write_cell_dens(stack,cell)
 type(stackdens),   intent(inout)  :: stack
 type(celldens),    intent(inout)  :: cell

 if (cell%waiting_index > stack%maxlength) call fatal('dens','attempting to write to invalid stack address')
 stack%cells(cell%waiting_index) = cell

end subroutine write_cell_dens

subroutine write_cell_force(stack,cell)
 type(stackforce),   intent(inout)  :: stack
 type(cellforce),    intent(inout)  :: cell

 if (cell%waiting_index > stack%maxlength) call fatal('force','attempting to write to invalid stack address')
 stack%cells(cell%waiting_index) = cell

end subroutine write_cell_force

subroutine reserve_stack_dens(stack,i)
 type(stackdens),    intent(inout) :: stack
 integer,            intent(out)   :: i

 !$omp atomic capture
 stack%n = stack%n + 1
 i = stack%n
 !$omp end atomic

 if (i > stack%maxlength) call fatal('dens','MPI stack exceeded')

end subroutine reserve_stack_dens

subroutine reserve_stack_force(stack,i)
 type(stackforce),   intent(inout) :: stack
 integer,            intent(out)   :: i

 !$omp atomic capture
 stack%n = stack%n + 1
 i = stack%n
 !$omp end atomic

 if (i > stack%maxlength) call fatal('force','MPI stack exceeded')

end subroutine reserve_stack_force

subroutine reset_stacks
 dens_stack_1%n=0
 dens_stack_2%n=0
 dens_stack_3%n=0

 force_stack_1%n=0
 force_stack_2%n=0
end subroutine reset_stacks

end module mpimemory
