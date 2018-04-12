!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION: Setup GR TDE
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS: none
!
!  DEPENDENCIES: infile_utils, io, part, physcon, setbinary, spherical,
!    timestep, units, metric, eos
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 ! real :: m1,m2,rp,semia,hacc1,rstar,norbits
 ! integer :: dumpsperorbit,nr

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas,set_particle_type
 use setbinary, only:set_binary,get_a_from_period
 use spherical, only:set_sphere
 use units,     only:set_units,umass,udist
 use physcon,   only:solarm,au,pi,solarr,km
 use io,        only:master,fatal
 use timestep,  only:tmax,dtmax
 use eos,       only:ieos
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer :: ierr,i,dumpsperorbit,nr
 logical :: iexist
 real    :: mstar,rstar,beta,ecc,rtidal,rp,semia,norbits,mhole
 real    :: vxyzstar(3),xyzstar(3),psep,period,hacc1,hacc2,massr

 !
 !--general parameters
 !
  time  = 0.
  polyk = 1.e-3    ! <== uconst
  gamma = 5./3.
  ieos  = 4

 !
 !--space available for injected gas particles
 !
  npart          = 0
  npartoftype(:) = 0
  xyzh(:,:)      = 0.
  vxyzu(:,:)     = 0.
  nptmass        = 0

!
!--Set central mass to M=1 in code units
!

 mhole = 1.e6*solarm ! (solar masses)
 call set_units(mass=mhole,c=1.,G=1.)

!
!--Default runtime parameters
!
 mstar         = 1.*solarm/umass
 rstar         = 1.*solarr/udist
 beta          = 5.
 ecc           = 0.8

 rtidal        = rstar*(mass1/mstar)**(1./3.)
 rp            = rtidal/beta
 semia         = rp/(1.-ecc)
 norbits       = 1.
 dumpsperorbit = 100
 nr            = 20
 psep          = rstar/nr

 period        = 2.*pi*sqrt(semia**3/mass1)
 tmax          = norbits*period
 dtmax         = period/dumpsperorbit

!
!--Read runtime parameters from setup file
!
 ! if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' White Dwarf tidal disruption'
 ! filename = trim(fileprefix)//'.setup'
 ! inquire(file=filename,exist=iexist)
 ! if (iexist) call read_setupfile(filename,ierr)
 ! if (.not. iexist .or. ierr /= 0) then
 !    if (id==master) then
 !       call write_setupfile(filename)
 !       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
 !    endif
 !    stop
 ! endif

!
!--Set a binary orbit given the desired orbital parameters
!
 hacc1    = rstar/1.e8 ! Something small
 hacc2    = hacc1
 massr    = mstar/mass1
 call set_binary(mass1,massr,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass)
 vxyzstar = vxyz_ptmass(1:3,2)
 xyzstar  = xyzmh_ptmass(1:3,2)

!
!--Delete sinks, replace with potential and a collection of gas particles
!
 nptmass = 0
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,xyz_origin=xyzstar)
 npartoftype(igas) = npart
 massoftype(igas)  = mstar/npart
 do i=1,npart
    call set_particle_type(i,igas)
    vxyzu(1:3,i) = vxyzstar(1:3)
 enddo

 if (npart == 0)   call fatal('setup','no hydro particles setup')
 if (ierr /= 0)    call fatal('setup','ERROR during setup')

end subroutine setpart


! !
! !---Read/write setup file--------------------------------------------------
! !
! subroutine write_setupfile(filename)
!  use infile_utils, only:write_inopt
!  character(len=*), intent(in) :: filename
!  integer, parameter :: iunit = 20
!
!  print "(a)",' writing setup options file '//trim(filename)
!  open(unit=iunit,file=filename,status='replace',form='formatted')
!  write(iunit,"(a)") '# input file for binary setup routines'
!  call write_inopt(m1,           'm1',           'mass of white dwarf (solar mass)',                 iunit)
!  call write_inopt(m2,           'm2',           'mass of asteroid (ceres mass)',                    iunit)
!  call write_inopt(rp,           'rp',           'pericentre distance (solar radii)',                iunit)
!  call write_inopt(semia,        'semia',        'semi-major axis (solar radii)',                    iunit)
!  call write_inopt(hacc1,        'hacc1',        'white dwarf (sink) accretion radius (solar radii)',iunit)
!  call write_inopt(rstar,    'rstar',    'radius of asteroid (km)',                          iunit)
!  call write_inopt(norbits,      'norbits',      'number of orbits',                                 iunit)
!  call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',                        iunit)
!  call write_inopt(nr           ,'nr'           ,'particles per asteroid radius (i.e. resolution)',  iunit)
!  close(iunit)
!
! end subroutine write_setupfile
!
! subroutine read_setupfile(filename,ierr)
!  use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
!  use io,           only:error
!  character(len=*), intent(in)  :: filename
!  integer,          intent(out) :: ierr
!  integer, parameter :: iunit = 21
!  integer :: nerr
!  type(inopts), allocatable :: db(:)
!
!  print "(a)",'reading setup options from '//trim(filename)
!  nerr = 0
!  ierr = 0
!  call open_db_from_file(db,filename,iunit,ierr)
!  call read_inopt(m1,           'm1',           db,min=0.,errcount=nerr)
!  call read_inopt(m2,           'm2',           db,min=0.,errcount=nerr)
!  call read_inopt(rp,           'rp',           db,min=0.,errcount=nerr)
!  call read_inopt(semia,        'semia',        db,min=0.,errcount=nerr)
!  call read_inopt(hacc1,        'hacc1',        db,min=0.,errcount=nerr)
!  call read_inopt(rstar,    'rstar',    db,min=0.,errcount=nerr)
!  call read_inopt(norbits,      'norbits',      db,min=0.,errcount=nerr)
!  call read_inopt(dumpsperorbit,'dumpsperorbit',db,min=0 ,errcount=nerr)
!  call read_inopt(nr,           'nr',           db,min=0 ,errcount=nerr)
!  call close_db(db)
!  if (nerr > 0) then
!     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
!     ierr = nerr
!  endif
!
! end subroutine read_setupfile

end module setup
