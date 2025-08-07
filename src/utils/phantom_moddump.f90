!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantommoddump
!
! This program is a simple utility for modifying a dump file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: phantom_moddump dumpfilein dumpfileout [time] [outformat] --maxp=50000000
!
! :Dependencies: checkconserved, checksetup, dim, eos, io, memory, moddump,
!   options, part, prompting, readwrite_dumps, readwrite_infile, setBfield,
!   setup_params, systemutils
!
 use dim,             only:tagline,maxp_alloc
 use eos,             only:polyk
 use part,            only:xyzh,hfact,massoftype,vxyzu,npart,npartoftype, &
                           Bxyz,Bextx,Bexty,Bextz,mhd
 use io,              only:set_io_unit_numbers,iprint,idisk1,warning,fatal,iwritein,id,master
 use readwrite_dumps, only:read_dump,write_fulldump,is_not_mhd
 use setBfield,       only:set_Bfield
 use moddump,         only:modify_dump,flags=>moddump_flags
 use readwrite_infile,only:write_infile,read_infile
 use options,         only:set_default_options
 use setup_params,    only:ihavesetupB
 use prompting,       only:prompt
 use checksetup,      only:check_setup
 use checkconserved,  only:get_conserv
 use memory,          only:allocate_memory
 use systemutils,     only:get_command_option
 implicit none
 integer :: nargs, i, nposargs
 character(len=120) :: dumpfilein,dumpfileout,arg,string
 real :: time,timeout
 integer :: ierr,nerr,nwarn,iloc
 logical :: idumpsphNG,iexist,ians
 integer, parameter          :: lenprefix = 120
 character(len=lenprefix)    :: fileprefix
 character(len=lenprefix+10) :: dumpfile,infile,evfile,logfile,progname

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 2) then
    print "(a,/)",trim(tagline)
    call get_command_argument(0,progname)
    print "(a)",' Usage: '//trim(progname)//' dumpfilein dumpfileout [time] [outformat] --maxp=50000000 '//trim(flags)
    stop
 endif

 ! Process positional arguments (those not starting with '-')
 idumpsphNG = .false.
 timeout = -1.
 nposargs = 0
 do i = 1, nargs
    call get_command_argument(i,arg)
    if (arg(1:1) /= '-') then
       nposargs = nposargs + 1
       select case(nposargs)
       case(1)
          dumpfilein = arg
       case(2)
          dumpfileout = arg
       case(3)
          string = arg
          read(string,*,iostat=ierr) timeout
          if (ierr /= 0) then
             if (index(string,'sphNG') /= 0) then
                idumpsphNG = .true.
             else
                stop 'error reading output time from command line'
             endif
          else
             print*,' setting time = ',timeout
          endif
       case(4)
          if (index(arg,'sphNG') /= 0) idumpsphNG = .true.
       end select
    endif
 enddo

 if (nposargs < 2 .or. len_trim(dumpfilein) == 0 .or. len_trim(dumpfileout) == 0) then
    print "(a,/)",trim(tagline)
    call get_command_argument(0,progname)
    print "(a)",' Usage: '//trim(progname)//' dumpfilein dumpfileout [time] [outformat] --maxp=50000000 '//trim(flags)
    stop
 endif

 print "(/,a,/)",' Phantom moddump: pimp my dumpfiles'
!
!--look for an existing input file with name corresponding to the INPUT dump file
!  read this if it exists
!
 call set_default_options
 iloc = index(dumpfilein,'_0')
 if (iloc > 1) then
    fileprefix = trim(dumpfilein(1:iloc-1))
 else
    fileprefix = trim(dumpfilein)
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)
 if (iexist) then
    print "(/,2a,/)",' Reading default values from ', trim(infile)
    call read_infile(infile,logfile,evfile,dumpfile)
 else
    print "(/,64('*'),/,a,/,64('*'))",&
      ' *** WARNING: '//trim(infile)//' NOT FOUND, USING DEFAULT OPTIONS ***'
 endif
!
!--look for an existing input file with name corresponding to the OUTPUT dump file
!  read this if it exists; this will overwrite the values from the input .in file
!
 iloc = index(dumpfileout,'_0')
 if (iloc > 1) then
    fileprefix = trim(dumpfileout(1:iloc-1))
 else
    fileprefix = trim(dumpfileout)
    dumpfileout = trim(dumpfileout)//'_00000'
 endif
 infile = trim(fileprefix)//'.in'
 inquire(file=trim(infile),exist=iexist)
 if (iexist) then
    print "(2a,/)",' Reading revised default values from ', trim(infile)
    call read_infile(infile,logfile,evfile,dumpfile)
 endif
!
!--reset logfile name
!
 logfile = trim(fileprefix)//'01.log'
 evfile  = trim(fileprefix)//'01.ev'
 if (mhd) then
    ihavesetupB = .true.
 else
    ihavesetupB = .false.
 endif
!
!--allocate memory BEFORE reading the first file
!  will be reallocated automatically if npart > maxp_alloc
!  but allows user to manually preset array sizes if necessary
!
 maxp_alloc = get_command_option('maxp',default=int(maxp_alloc))
 call allocate_memory(maxp_alloc)
!
!--read particle setup from dumpfile
!
 call read_dump(trim(dumpfilein),time,hfact,idisk1,iprint,0,1,ierr)
 if (timeout < 0.) timeout = time
 if (mhd .and. ierr==is_not_mhd) then
    ihavesetupB = .false.
 elseif (ierr /= 0) then
    stop 'error reading dumpfile'
 endif
 call check_setup(nerr,nwarn,restart=.true.)
 if (nwarn > 0) call warning('moddump','warnings from original setup',var='warnings',ival=nwarn)
 if (nerr > 0) call warning('moddump','ERRORS in original setup',var='errors',ival=nerr)
!
!--modify the dump file
!
 call modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 get_conserv = 1.
!
!--perform sanity checks on the output of modify_dump routine
!
 call check_setup(nerr,nwarn,restart=.true.)
 if (nwarn > 0) call warning('moddump','warnings from modified setup',var='warnings',ival=nwarn)
 if (nerr > 0)  call fatal('moddump','errors in modified setup',var='errors',ival=nerr)

 if (mhd) then
    ians = .false.
    if (.not.ihavesetupB) call prompt(' add/reset magnetic fields?',ians)
    if (ians) then
       call set_Bfield(npart,npartoftype(:),xyzh,massoftype(:),vxyzu,polyk, &
                       Bxyz,Bextx,Bexty,Bextz)

    endif
 endif

 call write_fulldump(timeout,dumpfileout,sphNG=idumpsphNG)
!
!--write a fresh input file, whether it exists or not
!
 if (id==master .and. .not.idumpsphNG) then
    call write_infile(infile,logfile,evfile,dumpfileout,iwritein,6)
    print "(a,/,/,a)",' To start the calculation, use: ',' ./phantom '//trim(infile)
 endif

 print "(/,a,/)",' Phantom moddump: another happy customer'

end program phantommoddump
