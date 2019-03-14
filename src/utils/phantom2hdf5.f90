!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: phantom2hdf5
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  USAGE: phantom2hdf5 dumpfile(s)
!
!  DEPENDENCIES: dim, io, part, readwrite_dumps, readwrite_dumps_hdf5
!+
!--------------------------------------------------------------------------
program phantom2hdf5
 use dim,                  only:tagline
 use part,                 only:hfact
 use io,                   only:set_io_unit_numbers,iprint,idisk1
 use readwrite_dumps,      only:read_dump,read_smalldump,write_fulldump
 use readwrite_dumps_hdf5, only:read_dump_hdf5=>read_dump, write_fulldump_hdf5=>write_fulldump
 use readwrite_dumps_hdf5, only:write_smalldump_hdf5=>write_smalldump
 implicit none
 integer :: nargs,iarg
 character(len=120) :: dumpfile
 real :: time
 integer :: ierr
 logical :: fulldump

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantom2hdf5 dumpfile(s)'
    stop
 endif

 print "(/,a,/)",' Phantom2hdf5: The best conversion in the west...'

 over_args: do iarg=1,nargs
    call get_command_argument(iarg,dumpfile)
   !
   !--read particle setup from dumpfile
   !
    fulldump = .true.
    call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)

    ! Try opening small dump if there is an error opening full dump
    if (ierr /= 0) then
       fulldump = .false.
       call read_smalldump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    endif

    ! If there is still an error, skip to the next file
    if (ierr /= 0) then
       print*,'error reading dumpfile: ',trim(dumpfile)
       print*,'skipping to next one...'
       cycle over_args
    endif

    if (fulldump) then
       call write_fulldump_hdf5(time,trim(dumpfile))
    else
       call write_smalldump_hdf5(time,trim(dumpfile))
    endif

 enddo over_args

 print "(/,a,/)",' Phantom2hdf5: Enjoy a better file format.'

end program phantom2hdf5
