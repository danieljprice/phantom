!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: io, mpimemory, physcon, testutils, units
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
 use io,        only:id
 use mpimemory, only:allocate_mpi_memory,increase_mpi_memory,&
                     deallocate_mpi_memory,stacksize
 integer, intent(inout) :: ntests,npass
 integer :: nerr(1)

 integer :: stacksize_orig

 nerr = 0

 print*,id,"Original stacksize=",stacksize

 ! Save original stacksize
 stacksize_orig = stacksize

 ! Deallocate existing stack
 call deallocate_mpi_memory

 ! Allocate the stacks again at half the size. Assuming that the original allocation
 ! was fine, increasing it by a factor of 1.5 should be no problem.
 call allocate_mpi_memory(stacksize_in=int(stacksize_orig*0.4))

 print*,id,"Allocated stacksize=",int(stacksize_orig*0.4)

 ! Trigger a stacksize increase - if this doesn't segfault, that's a good sign
 call increase_mpi_memory

 print*,id,"Increased to stacksize=",stacksize

 !  Reallocate previous stack
 call deallocate_mpi_memory
 call allocate_mpi_memory(stacksize_in=stacksize_orig)

 call update_test_scores(ntests,nerr,npass)

end subroutine test_increase_mpi_memory

end module testmpi
