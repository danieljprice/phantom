!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  PROGRAM: diffdumps
!
!  DESCRIPTION: This program is a simple utility for finding
!   any differences between two phantom dumps (e.g. for checking
!   that two calculations have produced the same answer)
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  USAGE: diffdumps firstdumpfilename seconddumpfilename [tolerance]
!
!  DEPENDENCIES: dim, io, part, readwrite_dumps, testutils
!+
!--------------------------------------------------------------------------
program diffdumps
 use dim,     only:maxp,maxvxyzu,tagline
 use part,    only:xyzh,vxyzu,npart,hfact
 use io,      only:set_io_unit_numbers,iprint,idisk1,real4
 use readwrite_dumps, only:read_dump
 use testutils,       only:checkval
 implicit none
 integer                      :: nargs
 character(len=120)           :: dumpfile,dumpfile2,tolstring
 real, allocatable :: xyzh2(:,:)
 real, allocatable :: vxyzu2(:,:)
 integer :: ierr,ndiff,ndiffx,ndiffy,ndiffz,ndiffh,ndiffvx,ndiffvy,ndiffvz,ndiffu
 real    :: time,time2,hfact2,tolerance
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

 allocate (xyzh2(4,maxp),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store positions'

 allocate (vxyzu2(maxvxyzu,maxp),stat=ierr)
 if (ierr /= 0) stop 'error allocating memory to store velocities'

 xyzh2 = xyzh
 vxyzu2 = vxyzu
!
!--read particle setup from second dumpfile
!
 call read_dump(trim(dumpfile2),time2,hfact2,idisk1+1,iprint,0,1,ierr)
 if (ierr /= 0) stop 'error reading second dumpfile'

! call diffarr(hfact2,hfact,0,'hfact',ndiffh)
! call diffarr(time2,time,0,'time',ndiffh)

 idiff = .false.
 ndiffu = 0
 call checkval(npart,xyzh(1,:),xyzh2(1,:),tolerance,ndiffx,'x')
 call checkval(npart,xyzh(2,:),xyzh2(2,:),tolerance,ndiffy,'y')
 call checkval(npart,xyzh(3,:),xyzh2(3,:),tolerance,ndiffz,'z')
 call checkval(npart,xyzh(4,:),xyzh2(4,:),tolerance,ndiffh,'h')

 call checkval(npart,vxyzu(1,:),vxyzu2(1,:),tolerance,ndiffvx,'vx')
 call checkval(npart,vxyzu(2,:),vxyzu2(2,:),tolerance,ndiffvy,'vy')
 call checkval(npart,vxyzu(3,:),vxyzu2(3,:),tolerance,ndiffvz,'vz')
 if (maxvxyzu >= 4) call checkval(npart,vxyzu(4,:),vxyzu2(4,:),tolerance,ndiffu,'u')

 print "(/,a,/)",' diffdumps: we wish you a pleasant journey '
 print *,'         positions differ ',max(ndiffx,ndiffy,ndiffz),' times'
 print *,' smoothing lengths differ ',ndiffh,' times'
 print *,'        velocities differ ',max(ndiffvx,ndiffvy,ndiffvz),' times'
 if (maxvxyzu >=4) print *,'  thermal energies differ ',ndiffu,' times'

 ndiff = ndiffx + ndiffy + ndiffz + ndiffh + ndiffvx + ndiffvy + ndiffvz + ndiffu
 if (allocated(xyzh2)) deallocate(xyzh2)
 if (allocated(vxyzu2)) deallocate(vxyzu2)

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

end program diffdumps

