!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program phantom2hdf5
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Usage: phantom2hdf5 dumpfile(s)
!
! :Dependencies: dim, eos, externalforces, io, part,
!   readwrite_dumps_fortran, readwrite_dumps_hdf5
!
 use dim,                     only:tagline
 use part,                    only:hfact,dt_in
 use io,                      only:set_io_unit_numbers,iprint,idisk1
 use readwrite_dumps_fortran, only:read_dump_fortran,read_smalldump_fortran,write_fulldump_fortran
 use readwrite_dumps_hdf5,    only:read_dump_hdf5,write_dump_hdf5
 use eos,                     only:extract_eos_from_hdr
 use externalforces,          only:extract_iextern_from_hdr
 implicit none
 integer :: nargs,iarg
 character(len=120) :: dumpfile
 real :: time
 integer :: ierr
 logical :: fulldump

 extract_eos_from_hdr     = .true.
 extract_iextern_from_hdr = .true.

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
    call read_dump_fortran(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)

    ! Try opening small dump if there is an error opening full dump
    if (ierr /= 0) then
       fulldump = .false.
       call read_smalldump_fortran(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
    endif

    ! If there is still an error, skip to the next file
    if (ierr /= 0) then
       print*,'error reading dumpfile: ',trim(dumpfile)
       print*,'skipping to next one...'
       cycle over_args
    endif

    call write_dump_hdf5(time,trim(dumpfile),fulldump=fulldump,dtind=dt_in)

 enddo over_args

 print "(/,a,/)",' Phantom2hdf5: Enjoy a better file format.'

end program phantom2hdf5
