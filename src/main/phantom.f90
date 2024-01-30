!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantom
!
! The Phantom SPH code, by Daniel Price.
!
!  This code is designed to be an ultra-sleek, ultra-low-memory,
!  code for high resolution SPH simulations
!
!  The requirements mean we need to store as few quantities as possible
!  (aim is to be able to run 10^7 particles in under 1Gb)
!  and to use the fastest possible implementation
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantom infilename
!
! :Dependencies: dim, evolve, initial, io, mpiutils
!
 use dim,             only:tagline
 use mpiutils,        only:init_mpi,finalise_mpi
 use initial,         only:initialise,finalise,startrun,endrun
 use io,              only:id,master,nprocs,set_io_unit_numbers,die
 use evolve,          only:evol
 implicit none
 integer            :: nargs
 character(len=120) :: infile,logfile,evfile,dumpfile

 id = 0

 call init_mpi(id,nprocs)
 call set_io_unit_numbers
 !
 ! get name of run from the command line
 !
 nargs = command_argument_count()
 if (nargs < 1) then
    if (id==master) then
       print "(a,/)",trim(tagline)
       print "(a)",' Usage: phantom infilename'
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
 !
 ! perform a simulation
 !
 if (index(infile,'.in')==0) infile = trim(infile)//'.in'
 call startrun(infile,logfile,evfile,dumpfile)
 call evol(infile,logfile,evfile,dumpfile)
 if (id==master) call endrun()

 call finalise_mpi()

end program phantom
