!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantomsetup
!
! Wrapper to routines for setting up Phantom simulations
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantomsetup fileprefix --maxp=10000000 --nprocsfake=1
!
! :Dependencies: boundary, checksetup, dim, eos, fileutils, gravwaveutils,
!   io, krome_interface, memory, mpidomain, mpiutils, options, part,
!   physcon, readwrite_dumps, readwrite_infile, setBfield, setup,
!   setup_params, systemutils, timestep, units
!
 use memory,          only:allocate_memory,deallocate_memory
 use dim,             only:tagline,maxvxyzu,mpi,&
                           ndivcurlv,ndivcurlB,maxp_hard
 use part,            only:xyzh,massoftype,hfact,vxyzu,npart,npartoftype, &
                           Bxyz,Bextx,Bexty,Bextz,rhoh,&
                           isetphase,igas,iamtype,labeltype,mhd,init_part
 use setBfield,       only:set_Bfield
 use eos,             only:polyk,gamma
 use io,              only:set_io_unit_numbers,id,master,nprocs,iwritein,fatal,warning
 use readwrite_dumps, only:init_readwrite_dumps,write_fulldump
 use readwrite_infile,only:write_infile,read_infile
 use options,         only:set_default_options
 use setup,           only:setpart
 use setup_params,    only:ihavesetupB,npart_total
 use checksetup,      only:check_setup
 use timestep,        only:time
 use physcon,         only:pi
 use units,           only:set_units,print_units,c_is_unity
 use mpiutils,        only:init_mpi,finalise_mpi,reduceall_mpi
 use mpidomain,       only:init_domains
 use boundary,        only:set_boundary
 use fileutils,       only:strip_extension
 use gravwaveutils,   only:calc_gravitwaves
 use systemutils,     only:get_command_option
#ifdef KROME
 use krome_interface, only:write_KromeSetupFile
#endif
 implicit none
 integer                     :: nargs,nprocsfake,nerr,nwarn,myid,myid1
 integer(kind=8)             :: ntotal,n_alloc
 integer, parameter          :: lenprefix = 120
 character(len=lenprefix)    :: fileprefix
 character(len=lenprefix+10) :: dumpfile,infile,evfile,logfile
 logical                     :: iexist

 nprocs = 1    ! for MPI, this is not initialised until init_mpi, but an initialised value is required for init_part
 call set_io_unit_numbers
 call set_units
 call set_boundary
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantomsetup fileprefix --maxp=10000000 --nprocsfake=1'
    print "(/,a)",' e.g. "phantomsetup mysim"'
    stop
 endif
 call get_command_argument(1,fileprefix)

 ! strip .in and .setup if present
 call strip_extension(fileprefix,'.in')
 call strip_extension(fileprefix,'.setup')

 if (index(fileprefix,'0000') /= 0) then
    print*,'Error: File prefix should not contain _0000'
    print*,'       (these are assigned automatically)'
    print "(/,a)",' e.g. "phantomsetup mysim"'
    stop
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)

!
!--In general, setup routines do not know the number of particles until they
!  are written. Need to allocate up to the hard limit. Legacy setup routines may
!  also rely on maxp being set to the number of desired particles. Allocate only
!  part, not kdtree or linklist
!
 n_alloc = get_command_option('maxp',default=maxp_hard)
 call allocate_memory(n_alloc, part_only=.true.)

 call set_default_options
!
!--if input file exists, read it
!
 if (iexist) call read_infile(infile,logfile,evfile,dumpfile)

!
!--reset logfile name
!
 logfile = trim(fileprefix)//'01.log'
!
!--setup particles
!
 time = 0.
 call init_part

 if (mpi) then
    call init_mpi(id,nprocs)
    call init_domains(nprocs)
    nprocsfake = 1
 else ! non-mpi
    nprocsfake = int(get_command_option('nprocsfake',default=1))
    nprocs= nprocsfake
    if (nprocs > 1) print*,' nprocs = ',nprocs
    call init_domains(nprocs)
    id = 0
 endif
 do myid=0,nprocsfake-1

    myid1 = myid
    if (mpi) myid1 = id
    call setpart(myid1,npart,npartoftype(:),xyzh,massoftype(:),vxyzu,polyk,gamma,hfact,time,fileprefix)
!
!--setup magnetic field if code compiled with MHD
!
    if (mhd .and. .not.ihavesetupB) then
       call set_Bfield(npart,npartoftype(:),xyzh,massoftype(:),vxyzu,polyk,Bxyz,Bextx,Bexty,Bextz)
    endif
!
!--perform sanity checks on the output of setpart routine
!
    call check_setup(nerr,nwarn)

    if (nwarn > 0) call warning('initial','warnings during particle setup',var='warnings',ival=nwarn)
    if (nerr > 0)  call fatal('initial','errors in particle setup',var='errors',ival=nerr)

    if (nprocsfake > 1) then
       ntotal = npart_total
    else
       ntotal = reduceall_mpi('+',npart)
    endif
    ! calculate gravitational wave strain automatically
    ! if code is run in relativistic units (c=1)
    if (c_is_unity()) calc_gravitwaves = .true.

    if (id==master .and. nerr==0 .and. nwarn==0) call print_units()
!
!--dumpfile name should end in .tmp unless density has been calculated
!  (never true using phantomsetup)
!
    dumpfile = trim(fileprefix)//'_00000.tmp'
    !if (index(dumpfile,'.tmp')==0) dumpfile = trim(dumpfile)//'.tmp'
!
!--write initial conditions to the dump file
!
    call init_readwrite_dumps()
    call write_fulldump(time,dumpfile,ntotal)
!
!--write an input file if it doesn't already exist
!
    if (id==master) then
       call write_infile(infile,logfile,evfile,dumpfile,iwritein,6)
       print "(a,/,/,a)",' To start the calculation, use: ',' ./phantom '//trim(infile)
    endif
 enddo

#ifdef KROME
 inquire(file='krome.setup',exist=iexist)
 if (.not. iexist) call write_KromeSetupFile
#endif

 call finalise_mpi
 call deallocate_memory(part_only=.true.)

end program phantomsetup
