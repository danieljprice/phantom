program phantom2hdf5
 use dim,                  only:tagline
 use part,                 only:hfact
 use io,                   only:set_io_unit_numbers,iprint,idisk1
 use readwrite_dumps,      only:read_dump,write_fulldump
 use readwrite_dumps_hdf5, only:read_dump_hdf5=>read_dump, write_fulldump_hdf5=>write_fulldump
 implicit none
 integer :: nargs
 character(len=120) :: dumpfile
 real :: time
 integer :: ierr

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs /= 1) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: phantom2hdf5 dumpfile'
    stop
 endif
 call get_command_argument(1,dumpfile)

 print "(/,a,/)",' Phantom2hdf5: The best conversion in the west...'
!
!--read particle setup from dumpfile
!
 call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading dumpfile'

 call write_fulldump_hdf5(time,trim(dumpfile))

 print "(/,a,/)",' Phantom2hdf5: Enjoy a better file format.'

end program phantom2hdf5
