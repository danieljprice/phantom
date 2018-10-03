!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantomanalysis
!
!  DESCRIPTION: This program is a wrapper for post-processing analysis tools
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: phantomanalysis dumpfile(s)
!
!  DEPENDENCIES: analysis, dim, eos, fileutils, infile_utils, io, part,
!    readwrite_dumps
!+
!--------------------------------------------------------------------------
program phantomanalysis
 use dim,             only:tagline
 use part,            only:xyzh,hfact,massoftype,vxyzu,npart !,npartoftype
 use io,              only:set_io_unit_numbers,iprint,idisk1,ievfile,ianalysis
 use readwrite_dumps, only:read_dump,read_smalldump,is_small_dump
 use infile_utils,    only:open_db_from_file,inopts,read_inopt,close_db
 use fileutils,       only:numfromfile,basename
 use analysis,        only:do_analysis,analysistype
 use eos,             only:ieos
 implicit none
 integer            :: nargs,iloc,ierr,iarg,i
 real               :: time
 logical            :: iexist
 character(len=120) :: dumpfile,fileprefix,infile
 type(inopts), dimension(:), allocatable :: db

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a)",trim(tagline)
    print "(a,/)",' Analysis module built for '//trim(analysistype)//' analysis'
    call get_command_argument(0,dumpfile)
    print "(a)",' Usage: '//trim(basename(dumpfile))//' dumpfile(s)'
    stop
 endif

 print "(/,a,/)",' Phantom analysis ('//trim(analysistype)//'): You data, we analyse'

 over_args: do iarg=1,nargs

    call get_command_argument(iarg,dumpfile)
!
!--If the first dumpfile, then read the .in file (if it exists) to obtain the equation of state
!
    if (iarg==1) then

       iloc = index(dumpfile,'_0')

       if (iloc > 1) then
          fileprefix = trim(dumpfile(1:iloc-1))
       else
          fileprefix = trim(dumpfile)
       endif
       infile = trim(fileprefix)//'.in'
       inquire(file=trim(infile),exist=iexist)
       if (iexist) then
          call open_db_from_file(db,infile,ianalysis,ierr)
          call read_inopt(ieos,'ieos',db,ierr)
          call close_db(db)
          close(ianalysis)
       endif
    endif
!
!--read particle setup from dumpfile
!
    if (index(analysistype,'header') /= 0) then
       !--only read the dumpfile header
       call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr,headeronly=.true.)
       if (ierr==0) print "(a,/)",' (finished reading file -- this analysis reads the header only)'
    else
       call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr,dustydisc=.true.)
    endif

    if (ierr==is_small_dump) then
       close(idisk1)
       !--if it is a small dump, look for a file called dumpfile.binary
       !  (this is a workaround to at least obtain what information there is in small dump files)
       open(idisk1,file=trim(dumpfile)//'.binary',form='unformatted',status='old',iostat=ierr)
       if (ierr==0) then
          print "(a)",' reading from '//trim(dumpfile)//'.binary for small dump'
          read(idisk1) time,npart
          print*,'npart = ',npart
          do i=1,npart
             read(idisk1) xyzh(1,i),xyzh(2,i),xyzh(3,i),massoftype(1),xyzh(4,i)
          enddo
          close(idisk1)
          vxyzu(:,:) = 0.
       else
          call read_smalldump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
          vxyzu(:,:) = 0.
       endif
    elseif (ierr /= 0) then
       stop 'error reading dumpfile'
    endif

    if (hfact < tiny(hfact)) then
       print "(a,f6.2,a)",' WARNING! hfact = ',hfact,' from dump file, resetting to 1.2'
       hfact = 1.2
    endif

    call do_analysis(trim(dumpfile),numfromfile(dumpfile),xyzh,vxyzu, &
                     massoftype(1),npart,time,ievfile)
 enddo over_args

 print "(/,a,/)",' Phantom analysis: may your paper be a happy one'

end program phantomanalysis
