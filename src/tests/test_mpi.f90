!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
 integer, intent(inout)   :: ntests,npass

 call set_units(mass=1.d6*solarm,G=1.d0,c=1.d0)
 if (id==master) write(*,"(/,a,/)") '--> TESTING MPI'
 call test_increase_mpi_memory(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- MPI TESTS COMPLETE'

end subroutine test_mpi

subroutine test_increase_mpi_memory(ntests,npass)
 use mpimemory, only:allocate_mpi_memory,increase_mpi_memory,&
                     deallocate_mpi_memory,stacksize,force_stack_1,&
                     push_onto_stack
 use mpiforce,  only:cellforce
 integer, intent(inout) :: ntests,npass
 integer, parameter :: new_stacksize=100
 type(cellforce) :: cell
 integer         :: nerr(3), ncheck(3), i, stacksize_orig
 real            :: maxerr(3)

 nerr = 0
 ncheck = 0
 maxerr = 0.

 ! Save original stacksize, assuming they're the same for dens and force
 stacksize_orig = stacksize

 ! Deallocate existing stack
 call deallocate_mpi_memory

 ! Allocate the stacks again at a smaller size.
 call allocate_mpi_memory(stacksize_in=new_stacksize)

 ! Write some data to each cell
 do i=1,new_stacksize
    cell%xpos = [1.,2.,3.] * i
    call push_onto_stack(force_stack_1, cell)
 enddo

 ! Ensure size of force_stack_1 is what we expect it to be
 call checkval(force_stack_1%n,new_stacksize,0,nerr(1),'stacksize after pushing cells')
 call update_test_scores(ntests,nerr,npass)

 ! Trigger a stacksize increase - if this doesn't segfault, that's a good sign
 call increase_mpi_memory

 ! Ensure stack size hasn't changed
 call checkval(force_stack_1%n,new_stacksize,0,nerr(1),'stacksize after mem increase')
 call update_test_scores(ntests,nerr,npass)

 nerr = 0
 ncheck = 0
 ! Check cell data is the same as what was written into cells above
 do i=1,new_stacksize
    call checkvalbuf(force_stack_1%cells(i)%xpos(1),1.*i,1.e-15,'error in xpos(1) after mem increase',nerr(1),ncheck(1),maxerr(1))
    call checkvalbuf(force_stack_1%cells(i)%xpos(2),2.*i,1.e-15,'error in xpos(2) after mem increase',nerr(2),ncheck(2),maxerr(2))
    call checkvalbuf(force_stack_1%cells(i)%xpos(3),3.*i,1.e-15,'error in xpos(3) after mem increase',nerr(3),ncheck(3),maxerr(3))
 enddo

 call checkvalbuf_end('error in xpos(1) after mem increase asfgd',ncheck(1),nerr(1),maxerr(1),1.e-15)
 call checkvalbuf_end('error in xpos(2) after mem increase asfgd',ncheck(2),nerr(2),maxerr(2),1.e-15)
 call checkvalbuf_end('error in xpos(3) after mem increase asfgd',ncheck(3),nerr(3),maxerr(3),1.e-15)
 call update_test_scores(ntests,nerr,npass)

 ! Reallocate previous stack
 call deallocate_mpi_memory
 call allocate_mpi_memory(stacksize_in=stacksize_orig)

end subroutine test_increase_mpi_memory

end module testmpi
