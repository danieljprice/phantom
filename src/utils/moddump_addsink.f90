!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Moddump to add a sink particle with user-configurable properties.
!
! A parameter file 'addsink.moddump' is read if present. If missing (or if
! parsing fails), a template is written and the program stops, prompting the
! user to edit the file and rerun phantommoddump.
!
! :References: None
!
! :Owner: Taj Jankovič & Aleksej Jurca
!
! :Runtime parameters (addsink.moddump):
!  - use_direct : *use direct specification of position/velocity?*
!  - msink      : *sink mass [code units] (default: 1 Msun in code units)*
!  - rsink      : *sink softening length used as stellar radius [code units]*
!  - racc       : *sink accretion radius [code units] (default: 0 to avoid dead particles)*
!  - z0         : *initial z position for inclination-mode [code units]*
!  - y0         : *initial y (impact parameter) for inclination-mode [code units]*
!  - incl_deg   : *inclination angle in degrees for inclination-mode*
!  - v0         : *speed for inclination-mode [code units]*
!  - x0,y0_dir,z0_dir : *direct initial position [code units]*
!  - vx0,vy0,vz0_dir  : *direct initial velocity [code units]*
!
! :Dependencies: infile_utils, io, part, physcon, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''
 character(len=*), parameter :: paramfile = 'addsink.moddump'

 !--module-level parameters so they can be read/written easily
 logical :: use_direct = .false.
 real    :: msink, rsink, racc
 real    :: z0, y0, incl_deg, v0                ! inclination-mode parameters
 real    :: x0, y0_dir, z0_dir                  ! direct-mode position
 real    :: vx0, vy0, vz0_dir                   ! direct-mode velocity

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,        only: nptmass, maxptmass, xyzmh_ptmass, vxyz_ptmass, ihsoft, ihacc, iJ2, iReff
 use io,          only: fatal, id, master
 use physcon,     only: deg_to_rad
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 logical :: iexist
 integer :: ierr
 real    :: xp(3), vp(3)
 real    :: incl, s

 !--defaults (will be overridden by addsink.moddump if present)
 call set_defaults_addsink()

 !--read parameter file (or write template and stop)
 inquire(file=paramfile,exist=iexist)
 if (iexist) call read_setupfile(paramfile,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(paramfile)
       print*,' Edit '//trim(paramfile)//' and rerun phantommoddump'
    endif
    stop
 endif

 !--bounds check to avoid out-of-bounds access on ptmass arrays
 if (nptmass >= maxptmass) then
    call fatal('moddump_addsink','Too many sink particles: nptmass >= maxptmass')
 endif
 nptmass = nptmass + 1

 !--construct initial position/velocity
 if (use_direct) then
    xp = (/ x0, y0_dir, z0_dir /)
    vp = (/ vx0, vy0, vz0_dir /)
 else
    incl = incl_deg * deg_to_rad
    s = sin(incl)
    if (abs(s) < 1.e-12) then
       call fatal('moddump_addsink','incl_deg too close to 0 or 180 degrees for inclination-mode (sin(incl) ~ 0)')
    endif
    ! Place sink so that (for y0=0) the straight-line trajectory passes through the origin:
    ! x0 = -z0*cot(incl), z0 = +z0, v = v0*(cos(incl), 0, -sin(incl)).
    xp = (/ -z0*cos(incl)/s, y0, z0 /)
    vp = (/ v0*cos(incl), 0.0, -v0*s /)
 endif

 !--assign to sink arrays
 xyzmh_ptmass(1:3,nptmass) = xp
 xyzmh_ptmass(4,nptmass)   = msink
 xyzmh_ptmass(ihacc,nptmass)  = racc
 xyzmh_ptmass(ihsoft,nptmass) = rsink   ! use softening length as Rstar
 xyzmh_ptmass(iJ2,nptmass)    = 0.0
 xyzmh_ptmass(iReff,nptmass)  = 0.0

 vxyz_ptmass(1:3,nptmass) = vp

 if (id==master) then
    print "(a,i0)",        'Added sink particle #', nptmass
    print "(a,es20.10)",   '  msink  = ', msink
    print "(a,es20.10)",   '  rsink  = ', rsink
    print "(a,es20.10)",   '  racc   = ', racc
    print "(a,3es20.10)",  '  x      = ', xp
    print "(a,3es20.10)",  '  v      = ', vp
 endif

 return
end subroutine modify_dump

!-----------------------------------------------------------------------
!+
!  Set default parameters for addsink moddump
!+
!-----------------------------------------------------------------------
subroutine set_defaults_addsink()
 use units,   only: umass
 use physcon, only: solarm
 ! defaults (code units)
 msink    = solarm/umass   ! 1 Msun in code units
 rsink    = 1.0            ! softening length used as "stellar radius"
 racc     = 0.0            ! accretion radius (0 avoids "dead" particles)
 use_direct = .false.

 ! inclination-mode defaults
 z0       = 5.0
 y0       = 0.0
 incl_deg = 90.0
 v0       = 1.0

 ! direct-mode defaults (only used if use_direct = .true.)
 x0       = 0.0
 y0_dir   = 0.0
 z0_dir   = z0
 vx0      = 0.0
 vy0      = 0.0
 vz0_dir  = -v0
end subroutine set_defaults_addsink

!-----------------------------------------------------------------------
!+
!  Write moddump parameter file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only: write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# parameter file for addsink moddump'
 write(iunit,"(a)") '# values are in code units unless otherwise stated'
 write(iunit,"(a)") '# if use_direct=F, (z0,y0,incl_deg,v0) define a straight-line approach'
 write(iunit,"(a)") '# if use_direct=T, (x0,y0_dir,z0_dir,vx0,vy0,vz0_dir) are used directly'

 call write_inopt(use_direct,'use_direct','use direct x0,v0 specification (T) or inclination-mode (F)',iunit)

 write(iunit,"(/,a)") '# sink properties'
 call write_inopt(msink,'msink','sink mass [code units] (default: 1 Msun in code units)',iunit)
 call write_inopt(rsink,'rsink','sink softening length used as stellar radius [code units]',iunit)
 call write_inopt(racc,'racc','sink accretion radius [code units] (default 0)',iunit)

 write(iunit,"(/,a)") '# inclination-mode (used if use_direct=F)'
 call write_inopt(z0,'z0','initial z position [code units]',iunit)
 call write_inopt(y0,'y0','initial y offset / impact parameter [code units]',iunit)
 call write_inopt(incl_deg,'incl_deg','inclination angle [deg]',iunit)
 call write_inopt(v0,'v0','speed [code units]',iunit)

 write(iunit,"(/,a)") '# direct-mode (used if use_direct=T)'
 call write_inopt(x0,'x0','initial x position [code units]',iunit)
 call write_inopt(y0_dir,'y0_dir','initial y position [code units]',iunit)
 call write_inopt(z0_dir,'z0_dir','initial z position [code units]',iunit)
 call write_inopt(vx0,'vx0','initial vx [code units]',iunit)
 call write_inopt(vy0,'vy0','initial vy [code units]',iunit)
 call write_inopt(vz0_dir,'vz0_dir','initial vz [code units]',iunit)

 close(iunit)
end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read moddump parameter file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only: open_db_from_file, inopts, read_inopt, close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 ierr = 0
 print "(a)",' reading moddump options from '//trim(filename)

 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) then
    call close_db(db)
    return
 endif

 call read_inopt(use_direct,'use_direct',db,errcount=nerr)

 call read_inopt(msink,'msink',db,min=0.0,errcount=nerr)
 call read_inopt(rsink,'rsink',db,min=0.0,errcount=nerr)
 call read_inopt(racc,'racc',db,min=0.0,errcount=nerr)

 call read_inopt(z0,'z0',db,errcount=nerr)
 call read_inopt(y0,'y0',db,errcount=nerr)
 call read_inopt(incl_deg,'incl_deg',db,min=0.0,max=180.0,errcount=nerr)
 call read_inopt(v0,'v0',db,errcount=nerr)

 call read_inopt(x0,'x0',db,errcount=nerr)
 call read_inopt(y0_dir,'y0_dir',db,errcount=nerr)
 call read_inopt(z0_dir,'z0_dir',db,errcount=nerr)
 call read_inopt(vx0,'vx0',db,errcount=nerr)
 call read_inopt(vy0,'vy0',db,errcount=nerr)
 call read_inopt(vz0_dir,'vz0_dir',db,errcount=nerr)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of moddump parameter file: re-writing...'
    ierr = nerr
 endif
 call close_db(db)

end subroutine read_setupfile

end module moddump
