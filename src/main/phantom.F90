!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom
!
!  DESCRIPTION: The Phantom SPH code, by Daniel Price.
!
!  This code is designed to be an ultra-sleek, ultra-low-memory,
!  code for high resolution SPH simulations
!
!  The requirements mean we need to store as few quantities as possible
!  (aim is to be able to run 10^7 particles in under 1Gb)
!  and to use the fastest possible implementation
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantom infilename
!
!  DEPENDENCIES: dim, evolve, initial, io, mpiderivs, mpiutils, stack, test
!+
!--------------------------------------------------------------------------
program phantom
 use memory,          only:allocate_memory
 use dim,             only:tagline,maxp_hard
 use mpiutils,        only:init_mpi,finalise_mpi
#ifdef MPI
 use mpiderivs,       only:init_tree_comms,finish_tree_comms
 use stack,           only:init_mpi_memory,finish_mpi_memory
#endif
 use initial,         only:initialise,startrun,endrun
 use io,              only:id,master,nprocs,set_io_unit_numbers,die
 use evolve,          only:evol
 use test,            only:testsuite
 implicit none
 integer :: nargs,i,ntests,npass,nfail
 character(len=120) :: infile,logfile,evfile,dumpfile

 id = 0
 ntests = 0
 npass  = 0
 nfail  = 0

 call init_mpi(id,nprocs)
 call set_io_unit_numbers
 !
 ! get name of run from the command line
 !
 nargs = command_argument_count()
 if (nargs < 1) then
    if (id==master) then
       print "(a,/)",trim(tagline)
       print "(a)",' Usage: phantom infilename '
    endif
    call die
 endif
 call get_command_argument(1,infile)
 !
 ! catch error if .setup is on command line instead of .in
 !
 if (index(infile,'.setup') > 0) then
    if (id==master) then
       print "(a,/)",trim(tagline)
       print "(a)",'ERROR: I think you mean ./phantomsetup '//trim(infile)
    endif
    call die
 endif

#ifdef MPI
 call init_tree_comms()
 call init_mpi_memory()
#endif
 if (trim(infile)=='test') then
    !
    ! run the phantom internal test suite
    !
    call initialise()
    call allocate_memory(maxp_hard)
    if (nargs >= 2) then
       do i=2,nargs
          call get_command_argument(i,infile)
          call testsuite(trim(infile),(i==2),(i==nargs),ntests,npass,nfail)
       enddo
    else
       call testsuite('all',.true.,.true.,ntests,npass,nfail)
    endif
 else
    !
    ! perform a simulation
    !
    if (index(infile,'.in')==0) then
       infile = trim(infile)//'.in'
    endif
    call startrun(infile,logfile,evfile,dumpfile)
    call evol(infile,logfile,evfile,dumpfile)
    if (id==master) call endrun()
 endif

#ifdef MPI
 call finish_tree_comms()
 call finish_mpi_memory()
#endif
 call finalise_mpi()

 !
 ! stop with an error code if test suite failed
 !
 if (ntests > 0 .and. nfail > 0) stop 666

end program phantom
