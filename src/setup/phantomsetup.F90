!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
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
!  DEPENDENCIES: boundary, centreofmass, checksetup, dim, domain, eos,
!    fileutils, io, mpiutils, options, part, physcon, readwrite_dumps,
!    readwrite_infile, setBfield, setup, setup_params, sortutils, units
!+
!--------------------------------------------------------------------------
program phantomsetup
 use memory,          only:allocate_memory, deallocate_memory
 use dim,             only:tagline,maxp,maxvxyzu,maxalpha,maxgrav,&
                           ndivcurlv,ndivcurlB
 use part,            only:xyzh,massoftype,hfact,vxyzu,npart,npartoftype, &
                           Bevol,Bxyz,Bextx,Bexty,Bextz,rhoh,iphase,maxphase,isetphase,igas,iamtype, &
                           labeltype,xyzmh_ptmass,vxyz_ptmass,maxp_h2,iHI,abundance,&
                           mhd,alphaind,divcurlv,divcurlB,poten,dustfrac
 use setBfield,       only:set_Bfield
 use eos,             only:polyk,gamma,en_from_utherm
 use io,              only:set_io_unit_numbers,id,master,nprocs,iwritein,fatal,warning
 use readwrite_dumps, only:write_fulldump
 use readwrite_infile,only:write_infile,read_infile
 use options,         only:set_default_options,use_dustfrac
 use setup,           only:setpart
 use setup_params,    only:ihavesetupB,npart_total
 use checksetup,      only:check_setup
 use physcon,         only:pi
 use units,           only:set_units,print_units
 use mpiutils,        only:init_mpi,finalise_mpi,use_mpi,reduceall_mpi
 use domain,          only:init_domains
 use boundary,        only:set_boundary
 use fileutils,       only:strip_extension
#ifdef SORT_RADIUS_INIT
 use sortutils,       only:indexxfunc,r2func_origin,set_r2func_origin
 use centreofmass,    only:get_centreofmass
#endif
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
#ifdef SORT_RADIUS_INIT
 integer :: iorder(maxp)
 real                     :: x0(3),v0(3)
#endif
 logical                  :: iexist

 call allocate_memory

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
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)

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
 npart = 0
 npartoftype(:) = 0
 massoftype(:)  = 0.
!--initialise point mass arrays to zero
 xyzmh_ptmass = 0.
 vxyz_ptmass  = 0.

 ! initialise arrays not passed to setup routine to zero
 if (mhd) Bevol = 0.
 if (maxphase > 0) iphase = 0 ! phases not set
 if (maxalpha==maxp)  alphaind = 0.
 if (ndivcurlv > 0) divcurlv = 0.
 if (ndivcurlB > 0) divcurlB = 0.
 if (maxgrav > 0) poten = 0.
 if (use_dustfrac) dustfrac = 0.
#ifdef LIGHTCURVE
 if (lightcurve) luminosity = 0.
#endif
!
!--initialise chemistry arrays if this has been compiled
!  (these may be altered by the specific setup routine)
!
 if (maxp_h2==maxp) then
    abundance(:,:)   = 0.
    abundance(iHI,:) = 1.  ! assume all atomic hydrogen initially
 endif

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
#ifdef SORT_RADIUS_INIT
    if (id==master) write(*,"(a)",ADVANCE='NO') ' Sorting particles by radius...'
    call get_centreofmass(x0,v0,npart,xyzh,vxyzu)
    if (id==master) print*,' setting origin for sort to ',x0
    call set_r2func_origin(x0(1),x0(2),x0(3))
    call indexxfunc(npart,r2func_origin,xyzh,iorder)
    if (id==master) write(*,"(a)") ' done'
    call write_fulldump(time,dumpfile,ntotal,iorder)
#else
    call write_fulldump(time,dumpfile,ntotal)
#endif
!
!--write an input file if it doesn't already exist
!
    if (id==master) then
       call write_infile(infile,logfile,evfile,dumpfile,iwritein,6)
       print "(a,/,/,a)",' To start the calculation, use: ',' ./phantom '//trim(infile)
    endif
 enddo

 call finalise_mpi()

 call deallocate_memory

end program phantomsetup
