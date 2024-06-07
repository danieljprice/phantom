!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantomtest
!
! Run all unit tests of the Phantom SPH code
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantomtest [no arguments]
!
! :Dependencies: initial, io, memory, mpiutils, test
!
 use memory,          only:allocate_memory
 use mpiutils,        only:init_mpi,finalise_mpi
 use initial,         only:initialise,finalise
 use io,              only:id,nprocs,set_io_unit_numbers
 use test,            only:testsuite
 implicit none
 integer :: nargs,i,ntests,npass,nfail
 character(len=120) :: string
 integer(kind=8), parameter :: maxp_test = 1000000

 ntests = 0
 npass  = 0
 nfail  = 0

 call init_mpi(id,nprocs)
 call set_io_unit_numbers
 !
 ! run the phantom internal test suite
 !
 call initialise()
 call allocate_memory(maxp_test)

 nargs = command_argument_count()
 if (nargs >= 1) then
    !
    ! extract command line arguments to run particular tests
    !
    do i=1,nargs
       call get_command_argument(i,string)
       call testsuite(trim(string),(i==1),(i==nargs),ntests,npass,nfail)
    enddo
 else
    call testsuite('all',.true.,.true.,ntests,npass,nfail)
 endif
 call finalise()
 call finalise_mpi()
 !
 ! stop with an error code if test suite failed
 !
 if (ntests > 0 .and. nfail > 0) stop 666

end program phantomtest
