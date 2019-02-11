!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: stack
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
!  DEPENDENCIES: dim, io, mpidens, mpiforce
!+
!--------------------------------------------------------------------------
module stack
#ifdef MPI
 use dim,         only:stacksize
 use io,          only:fatal
 use mpidens,     only:celldens,stackdens
 use mpiforce,    only:cellforce,stackforce

 implicit none

 interface allocate_stack
  module procedure allocate_stack_dens,allocate_stack_force
 end interface

 interface deallocate_stack
  module procedure deallocate_stack_dens,deallocate_stack_force
 end interface

 interface swap_stacks ! force doesn't require a stack swap
  module procedure swap_stacks_dens
 end interface

 interface push_onto_stack
  module procedure push_onto_stack_dens,push_onto_stack_force
 end interface

 interface pop_off_stack
  module procedure pop_off_stack_dens,pop_off_stack_force
 end interface

 interface reserve_stack
  module procedure reserve_stack_dens,reserve_stack_force
 end interface

 public :: init_mpi_memory
 public :: finish_mpi_memory
 public :: allocate_stack
 public :: deallocate_stack
 public :: swap_stacks
 public :: push_onto_stack
 public :: pop_off_stack
 public :: reserve_stack

 ! primary chunk of memory requested using alloc
 integer, parameter  :: n_dens_cells  = stacksize*3
 integer, parameter  :: n_force_cells = stacksize*2
 type(celldens),  allocatable, target :: dens_cells(:)
 type(cellforce), allocatable, target :: force_cells(:)

 ! memory allocation counters
 integer             :: idens
 integer             :: iforce

 ! stacks to be referenced from density and force routines
 type(stackdens)     :: dens_stack_1
 type(stackdens)     :: dens_stack_2
 type(stackdens)     :: dens_stack_3
 type(stackforce)    :: force_stack_1
 type(stackforce)    :: force_stack_2

contains

subroutine init_mpi_memory
 integer :: idens, iforce ! memory allocation counters
 integer :: allocstat

 allocate(dens_cells(n_dens_cells), stat = allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 idens = 1
 call allocate_stack(dens_stack_1, idens)
 call allocate_stack(dens_stack_2, idens)
 call allocate_stack(dens_stack_3, idens)
 if (idens - 1 > n_dens_cells) call fatal('stack','phantom memory allocation error')

 allocate(force_cells(n_force_cells), stat = allocstat)
 if (allocstat /= 0) call fatal('stack','fortran memory allocation error')
 iforce = 1
 call allocate_stack(force_stack_1,iforce)
 call allocate_stack(force_stack_2,iforce)
 if (iforce - 1 > n_force_cells) call fatal('stack','phantom memory allocation error')
end subroutine init_mpi_memory

subroutine finish_mpi_memory
 !
 !--Only called at the end, so deallocation_stack is not strictly necessary.
 !  May be useful in the future. TODO: Allow for unordered deallocation.
 !

 ! call deallocate_stack(dens_stack_3, idens)
 ! call deallocate_stack(dens_stack_2, idens)
 ! call deallocate_stack(dens_stack_1, idens)
 deallocate(dens_cells)

 ! call deallocate_stack(force_stack_2, idens)
 ! call deallocate_stack(force_stack_1, idens)
 deallocate(force_cells)
end subroutine finish_mpi_memory

subroutine allocate_stack_dens(stack, i)
 type(stackdens), intent(inout) :: stack
 integer,         intent(inout) :: i

 stack%mem_start = i
 stack%mem_end   = i + stacksize - 1

 stack%cells => dens_cells(stack%mem_start:stack%mem_end)
 stack%maxlength = stacksize

 i = i + stacksize
end subroutine allocate_stack_dens

subroutine allocate_stack_force(stack, i)
 type(stackforce),   intent(inout) :: stack
 integer,            intent(inout) :: i

 stack%mem_start = i
 stack%mem_end   = i + stacksize - 1

 stack%cells => force_cells(stack%mem_start:stack%mem_end)
 stack%maxlength = stacksize

 i = i + stacksize
end subroutine allocate_stack_force

subroutine deallocate_stack_dens(stack, i)
 type(stackdens),  intent(inout) :: stack
 integer,          intent(inout) :: i

 if (i - 1 /= stack%mem_end) call fatal('stack','memory deallocation not from top of stack')
 i = stack%mem_start - 1
 stack%mem_start = -1
 stack%mem_end   = -1
end subroutine deallocate_stack_dens

subroutine deallocate_stack_force(stack, i)
 type(stackforce),  intent(inout) :: stack
 integer,           intent(inout) :: i

 if (i - 1 /= stack%mem_end) call fatal('stack','memory deallocation not from top of stack')
 i = stack%mem_start - 1
 stack%mem_start = -1
 stack%mem_end   = -1
end subroutine deallocate_stack_force

subroutine swap_stacks_dens(stack_a, stack_b)
 type(stackdens),   intent(inout) :: stack_a
 type(stackdens),   intent(inout) :: stack_b

 integer :: temp_n
 integer :: temp_start
 integer :: temp_end

 if (stack_a%maxlength /= stack_b%maxlength) call fatal('stack', 'stack swap of unequal size')

 ! counters
 temp_n = stack_a%n
 stack_a%n = stack_b%n
 stack_b%n = temp_n

 ! addresses
 temp_start = stack_a%mem_start
 temp_end   = stack_a%mem_end
 stack_a%mem_start = stack_b%mem_start
 stack_a%mem_end   = stack_b%mem_end
 stack_b%mem_start = temp_start
 stack_b%mem_end   = temp_end

 ! change pointers
 stack_a%cells => dens_cells(stack_a%mem_start:stack_a%mem_end)
 stack_b%cells => dens_cells(stack_b%mem_start:stack_b%mem_end)

end subroutine swap_stacks_dens

subroutine push_onto_stack_dens(stack,cell)
 type(stackdens),    intent(inout)  :: stack
 type(celldens),     intent(in)     :: cell

 stack%n = stack%n + 1
 if (stack%n > stack%maxlength) call fatal('density','stack overflow')
 stack%cells(stack%n) = cell
end subroutine push_onto_stack_dens

subroutine push_onto_stack_force(stack,cell)
 type(stackforce),   intent(inout)  :: stack
 type(cellforce),    intent(in)     :: cell

 stack%n = stack%n + 1
 if (stack%n > stack%maxlength) call fatal('force','stack overflow')
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

subroutine reserve_stack_dens(stack,i)
 type(stackdens),    intent(inout) :: stack
 integer,            intent(out)   :: i

 stack%n = stack%n + 1
 i = stack%n
 if (stack%n > stack%maxlength) call fatal('density','stack overflow')
end subroutine reserve_stack_dens

subroutine reserve_stack_force(stack,i)
 type(stackforce),   intent(inout) :: stack
 integer,            intent(out)   :: i

 stack%n = stack%n + 1
 i = stack%n
 if (stack%n > stack%maxlength) call fatal('force','stack overflow')
end subroutine reserve_stack_force
#endif
end module stack
