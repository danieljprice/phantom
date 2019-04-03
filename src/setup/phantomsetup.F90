!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantomsetup
!
!  DESCRIPTION: Wrapper to routines for setting up Phantom simulations
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantomsetup fileprefix [nprocsfake]
!
!  DEPENDENCIES: boundary, checksetup, dim, domain, eos, fileutils, io,
!    memory, mpiutils, options, part, physcon, readwrite_dumps,
!    readwrite_infile, setBfield, setup, setup_params, units
!+
!--------------------------------------------------------------------------
program phantomsetup
 use memory,          only:allocate_memory,deallocate_memory
 use dim,             only:tagline,maxp,maxvxyzu,&
                           ndivcurlv,ndivcurlB,maxp_hard
 use part,            only:xyzh,massoftype,hfact,vxyzu,npart,npartoftype, &
                           Bxyz,Bextx,Bexty,Bextz,rhoh,iphase,maxphase,&
                           isetphase,igas,iamtype,labeltype,mhd,init_part
 use setBfield,       only:set_Bfield
 use eos,             only:polyk,gamma,en_from_utherm
 use io,              only:set_io_unit_numbers,id,master,nprocs,iwritein,fatal,warning
 use readwrite_dumps, only:write_fulldump
 use readwrite_infile,only:write_infile,read_infile
 use options,         only:set_default_options
 use setup,           only:setpart
 use setup_params,    only:ihavesetupB,npart_total
 use checksetup,      only:check_setup
 use physcon,         only:pi
 use units,           only:set_units,print_units
 use mpiutils,        only:init_mpi,finalise_mpi,use_mpi,reduceall_mpi
 use domain,          only:init_domains
 use boundary,        only:set_boundary
 use fileutils,       only:strip_extension
#ifdef LIGHTCURVE
 use part,            only:luminosity,maxlum,lightcurve
#endif
 implicit none
 integer                     :: nargs,i,nprocsfake,nerr,nwarn,myid,myid1
 integer(kind=8)             :: ntotal
 integer, parameter          :: lenprefix = 120
 character(len=lenprefix)    :: fileprefix
 character(len=lenprefix+10) :: dumpfile,infile,evfile,logfile,string
 real                        :: time,pmassi
 logical                     :: iexist

 call set_io_unit_numbers
 call set_units
 call set_boundary
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantomsetup fileprefix [nprocsfake]'
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
 elseif (fileprefix=='test') then
    print*,'Error: cannot use ''test'' as the job name, please rename your .setup file'
    stop
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)

 call set_default_options
!
!--if input file exists, read it
!
 if (iexist) call read_infile(infile,logfile,evfile,dumpfile)

!
!--In general, setup routines do not know the number of particles until they
!  are written. Need to allocate up to the hard limit. Legacy setup routines may
!  also rely on maxp being set to the number of desired particles. Allocate only
!  part, not kdtree or linklist
!
 call allocate_memory(maxp_hard, part_only=.true.)

!
!--reset logfile name
!
 logfile = trim(fileprefix)//'01.log'
!
!--setup particles
!
 time = 0.
 call init_part

 if (use_mpi) then
    call init_mpi(id,nprocs)
    call init_domains(nprocs)
    nprocsfake = 1
 else ! non-mpi
    if (nargs >= 3) then
       call get_command_argument(3,string)
       read(string,*) nprocsfake
    else
       nprocsfake = 1
    endif
    nprocs= nprocsfake
    print*,' nprocs = ',nprocs
    call init_domains(nprocs)
    id = 0
 endif

 do myid=0,nprocsfake-1

    myid1 = myid
    if (use_mpi) myid1 = id
    call setpart(myid1,npart,npartoftype(:),xyzh,massoftype(:),vxyzu,polyk,gamma,hfact,time,fileprefix)
!
!--setup magnetic field if code compiled with MHD
!
    if (mhd .and. .not.ihavesetupB) then
       call set_Bfield(npart,npartoftype(:),xyzh,massoftype(:),vxyzu,polyk, &
                       Bxyz,Bextx,Bexty,Bextz)
    endif
!
!--perform sanity checks on the output of setpart routine
!
    call check_setup(nerr,nwarn)
    if (nwarn > 0) call warning('initial','warnings during particle setup',var='warnings',ival=nwarn)
    if (nerr > 0)  call fatal('initial','errors in particle setup',var='errors',ival=nerr)
!
!--setup defines thermal energy: if we are using the entropy then
!  we need to convert this into the entropy variable before writing the dump file
!  (the dump file write converts back to utherm)
!
    if (maxvxyzu==4) then
       pmassi = massoftype(igas)
       do i=1,npart
          if (maxphase==maxp) pmassi = massoftype(iamtype(iphase(i)))
          vxyzu(maxvxyzu,i) = en_from_utherm(vxyzu(maxvxyzu,i),rhoh(xyzh(4,i),pmassi))
       enddo
    endif

    if (nprocsfake > 1) then
       ntotal = npart_total
    else
       ntotal = reduceall_mpi('+',npart)
    endif
    if (id==master) call print_units()
!
!--dumpfile name should end in .tmp unless density has been calculated
!  (never true using phantomsetup)
!
    dumpfile = trim(fileprefix)//'_00000.tmp'
    !if (index(dumpfile,'.tmp')==0) dumpfile = trim(dumpfile)//'.tmp'
!
!--write initial conditions to the dump file
!
    call write_fulldump(time,dumpfile,ntotal)
!
!--write an input file if it doesn't already exist
!
    if (id==master) then
       call write_infile(infile,logfile,evfile,dumpfile,iwritein,6)
       print "(a,/,/,a)",' To start the calculation, use: ',' ./phantom '//trim(infile)
    endif
 enddo

 call finalise_mpi
 call deallocate_memory(part_only=.true.)

end program phantomsetup
