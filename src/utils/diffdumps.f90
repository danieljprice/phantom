!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
program diffdumps
!
! This program is a simple utility for finding
!   any differences between two phantom dumps (e.g. for checking
!   that two calculations have produced the same answer)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Usage: diffdumps firstdumpfilename seconddumpfilename [tolerance]
!
! :Dependencies: dim, io, part, readwrite_dumps, sort_particles, testutils
!
 use dim,     only:maxp,maxvxyzu,tagline
 use part,    only:xyzh,vxyzu,npart,hfact,iorig
 use io,      only:set_io_unit_numbers,iprint,idisk1,real4
 use readwrite_dumps, only:read_dump
 use testutils,       only:checkval
 use sort_particles,  only:sort_part_id
 implicit none
 integer                      :: nargs
 character(len=120)           :: dumpfile,dumpfile2,tolstring
 integer(kind=8), allocatable :: iorig2(:)
 real,            allocatable :: xyzh2(:,:)
 real,            allocatable :: vxyzu2(:,:)
 integer :: ierr,ndiff,ndiffx,ndiffy,ndiffz,ndiffh,ndiffvx,ndiffvy,ndiffvz,ndiffu,ndiffid
 real    :: time,time2,hfact2,tolerance,err(8)
 logical :: idiff

 call set_io_unit_numbers
 iprint = 6
!
!--get name of run from the command line
!
 nargs = command_argument_count()
 if (nargs < 2) then
    print "(a,/)",trim(tagline)
    print "(a)",' Usage: diffdumps firstdumpfilename seconddumpfilename [tolerance]'
    stop
 endif
 call get_command_argument(1,dumpfile)
 call get_command_argument(2,dumpfile2)
 if (nargs >= 3) then
    call get_command_argument(3,tolstring)
    read(tolstring,*)tolerance
 else
    tolerance = tiny(0.)
 endif

 print "(/,a,/)",' diffdumps: we welcome you'

!
!--read particle setup from first dumpfile
!
 call read_dump(trim(dumpfile),time,hfact,idisk1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading first dumpfile'

 allocate (iorig2(maxp),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store iorig'

 allocate (xyzh2(4,maxp),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store positions'

 allocate (vxyzu2(maxvxyzu,maxp),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store velocities'

!  Sort first dump by ID
 call sort_part_id

 iorig2 = iorig
 xyzh2 = xyzh
 vxyzu2 = vxyzu
!
!--read particle setup from second dumpfile
!
 call read_dump(trim(dumpfile2),time2,hfact2,idisk1+1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading second dumpfile'

!  Sort second dump by ID
 call sort_part_id

! call diffarr(hfact2,hfact,0,'hfact',ndiffh)
! call diffarr(time2,time,0,'time',ndiffh)

 idiff = .false.
 ndiffu = 0
 err = 0.
 call checkval(npart,iorig(:),iorig2(:),int(0,kind=8),ndiffid,'id')

 call checkval(npart,xyzh(1,:),xyzh2(1,:),tolerance,ndiffx,'x',rmserr=err(1))
 call checkval(npart,xyzh(2,:),xyzh2(2,:),tolerance,ndiffy,'y',rmserr=err(2))
 call checkval(npart,xyzh(3,:),xyzh2(3,:),tolerance,ndiffz,'z',rmserr=err(3))
 call checkval(npart,xyzh(4,:),xyzh2(4,:),tolerance,ndiffh,'h',rmserr=err(4))

 call checkval(npart,vxyzu(1,:),vxyzu2(1,:),tolerance,ndiffvx,'vx',rmserr=err(5))
 call checkval(npart,vxyzu(2,:),vxyzu2(2,:),tolerance,ndiffvy,'vy',rmserr=err(6))
 call checkval(npart,vxyzu(3,:),vxyzu2(3,:),tolerance,ndiffvz,'vz',rmserr=err(7))
 if (maxvxyzu >= 4) call checkval(npart,vxyzu(4,:),vxyzu2(4,:),tolerance,ndiffu,'u',rmserr=err(8))

 print "(/,a,/)",' diffdumps: we wish you a pleasant journey '
 print *,'      particle IDs differ ',ndiffid,' times'
 print *,'         positions differ ',max(ndiffx,ndiffy,ndiffz),' times'
 print *,' smoothing lengths differ ',ndiffh,' times'
 print *,'        velocities differ ',max(ndiffvx,ndiffvy,ndiffvz),' times'
 if (maxvxyzu >=4) print *,'  thermal energies differ ',ndiffu,' times'

 ndiff = ndiffid + ndiffx + ndiffy + ndiffz + ndiffh + ndiffvx + ndiffvy + ndiffvz + ndiffu
 if (allocated(iorig2)) deallocate(iorig2)
 if (allocated(xyzh2)) deallocate(xyzh2)
 if (allocated(vxyzu2)) deallocate(vxyzu2)

 print "(/a,es10.4)",'MAX RMS ERROR: ',maxval(err)

 if (ndiff > 0) then
    print "(/,a)",' FILES DIFFER'
    if (tolerance <= tiny(0.)) then
       print "(/,a,/,/,a/)",&
        ' To specify a tolerance, use:', &
        '   diffdumps '//trim(dumpfile)//' '//trim(dumpfile2)//' 1.e-15'
    endif
 else
    print "(/,a)",' FILES ARE IDENTICAL '
 endif

 if (ndiff > 0) then
    call exit(1)
 endif

end program diffdumps
