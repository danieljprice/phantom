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
 use mpimemory, only: increase_mpi_memory
 integer, intent(inout) :: ntests,npass
 integer :: nerr(1)

 nerr = 0

 call increase_mpi_memory

 call update_test_scores(ntests,nerr,npass)

end subroutine test_increase_mpi_memory

end module testmpi
