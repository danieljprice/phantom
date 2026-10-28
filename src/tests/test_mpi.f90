!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testmpi
!
! MPI unit tests
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: io, mpiforce, mpimemory, physcon, testutils, units
!
 use testutils, only:checkval,checkvalbuf,checkvalbuf_end,update_test_scores
 implicit none

 public :: test_mpi

 private

contains

subroutine test_mpi(ntests,npass)
 use io,      only:id,master
 use units,   only:set_units
 use physcon, only:solarm
 integer, intent(inout) :: ntests,npass

 call set_units(mass=1.d6*solarm,G=1.d0,c=1.d0)
 if (id==master) write(*,"(/,a,/)") '--> TESTING MPI'
 call test_grow_mpi_stack(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- MPI TESTS COMPLETE'

end subroutine test_mpi

!-----------------------------------------------------------------------
!+
!  test that the MPI stacks grow automatically when more cells are
!  pushed than fit in the initial allocation, including when pushed
!  concurrently from many OpenMP threads
!+
!-----------------------------------------------------------------------
subroutine test_grow_mpi_stack(ntests,npass)
 use mpimemory, only:allocate_mpi_memory,deallocate_mpi_memory,stacksize,&
                     force_stack_1,push_onto_stack,get_cell,reset_stacks
 use mpiforce,  only:cellforce
 integer, intent(inout) :: ntests,npass
 integer, parameter :: tiny_stacksize = 4, ncells = 1000
 type(cellforce) :: cell
 integer         :: nerr(3), ncheck(3), i, j, stacksize_orig
 logical         :: found(ncells)
 real            :: maxerr(3)

 nerr = 0
 ncheck = 0
 maxerr = 0.

 ! Save original stacksize, assuming they're the same for dens and force
 stacksize_orig = stacksize

 ! Allocate the stacks again at a tiny size, so they must grow
 call deallocate_mpi_memory
 call allocate_mpi_memory(stacksize_in=tiny_stacksize)

 ! Push cells onto the stack in serial
 do i=1,ncells
    cell%xpos = [1.,2.,3.] * i
    call push_onto_stack(force_stack_1, cell)
 enddo
 call checkval(force_stack_1%n,ncells,0,nerr(1),'number of cells on stack after growth')
 call update_test_scores(ntests,nerr,npass)

 ! Check cell data is the same as what was written into cells above
 nerr = 0
 do i=1,ncells
    cell = get_cell(force_stack_1,i)
    call checkvalbuf(cell%xpos(1),1.*i,1.e-15,'xpos(1) after stack growth',nerr(1),ncheck(1),maxerr(1))
    call checkvalbuf(cell%xpos(2),2.*i,1.e-15,'xpos(2) after stack growth',nerr(2),ncheck(2),maxerr(2))
    call checkvalbuf(cell%xpos(3),3.*i,1.e-15,'xpos(3) after stack growth',nerr(3),ncheck(3),maxerr(3))
 enddo
 call checkvalbuf_end('xpos(1) after stack growth',ncheck(1),nerr(1),maxerr(1),1.e-15)
 call checkvalbuf_end('xpos(2) after stack growth',ncheck(2),nerr(2),maxerr(2),1.e-15)
 call checkvalbuf_end('xpos(3) after stack growth',ncheck(3),nerr(3),maxerr(3),1.e-15)
 call update_test_scores(ntests,nerr,npass)

 ! Push cells from many threads at once: each cell must appear exactly once
 call reset_stacks
 !$omp parallel do default(none) shared(force_stack_1) private(i,cell) schedule(dynamic,1)
 do i=1,ncells
    cell%xpos = [1.,2.,3.] * i
    call push_onto_stack(force_stack_1, cell)
 enddo
 !$omp end parallel do

 nerr = 0
 call checkval(force_stack_1%n,ncells,0,nerr(1),'number of cells on stack after threaded growth')
 found = .false.
 do i=1,ncells
    cell = get_cell(force_stack_1,i)
    j = nint(cell%xpos(1))
    if (j >= 1 .and. j <= ncells) then
       if (.not.found(j) .and. abs(cell%xpos(3) - 3.*j) < 1.e-15) found(j) = .true.
    endif
 enddo
 call checkval(count(found),ncells,0,nerr(2),'distinct cells on stack after threaded growth')
 call update_test_scores(ntests,nerr,npass)

 ! Reallocate previous stack
 call deallocate_mpi_memory
 call allocate_mpi_memory(stacksize_in=stacksize_orig)

end subroutine test_grow_mpi_stack

end module testmpi
