!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
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
!  DEPENDENCIES: io, mpidens, mpiforce
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

 interface push_onto_stack
  module procedure push_onto_stack_dens,push_onto_stack_force
 end interface

 interface pop_off_stack
  module procedure pop_off_stack_dens,pop_off_stack_force
 end interface

 interface reserve_stack
  module procedure reserve_stack_dens,reserve_stack_force
 end interface


 public :: allocate_stack
 public :: deallocate_stack
 public :: push_onto_stack
 public :: pop_off_stack
 public :: reserve_stack

contains

subroutine allocate_stack_dens(stack)
   type(stackdens),  intent(inout) :: stack

   allocate(stack%cells(stacksize))
   stack%maxlength = stacksize
end subroutine allocate_stack_dens

subroutine allocate_stack_force(stack)
   type(stackforce),  intent(inout) :: stack

   allocate(stack%cells(stacksize))
   stack%maxlength = stacksize
end subroutine allocate_stack_force

subroutine deallocate_stack_dens(stack)
   type(stackdens),  intent(inout) :: stack

   deallocate(stack%cells)
   stack%maxlength = 0
end subroutine deallocate_stack_dens

subroutine deallocate_stack_force(stack)
   type(stackforce),  intent(inout) :: stack

   deallocate(stack%cells)
   stack%maxlength = 0
end subroutine deallocate_stack_force

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
