!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION: Setup a polytrope in flat space, Minkowski metric (special relativity)
!
!  REFERENCES: None
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!
!  DEPENDENCIES: infile_utils, io, part, physcon, setbinary, spherical,
!    timestep, units, metric, eos
!+
!--------------------------------------------------------------------------
module setup
 implicit none

 public :: setpart

 private

 integer :: nr = 25 !-- Default number of particles in star radius

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,        only:igas,set_particle_type,rhoh
 use spherical,   only:set_sphere
 use units,       only:set_units,umass,udist
 use physcon,     only:solarm,solarr
 use io,          only:master,fatal
 use timestep,    only:tmax,dtmax
 use options,     only:nfulldump
 use eos,         only:ieos
 use rho_profile, only:rho_polytrope
 use prompting,   only:prompt
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename,infile
 integer, parameter :: ntab=5000
 integer :: i,npts,ierr
 real    :: psep
 real    :: rtab(ntab),rhotab(ntab)
 real    :: densi,mstar,rstar
 logical :: iexist

 infile = trim(fileprefix)//'.in'
 iexist = .false.
 inquire(file=trim(infile),exist=iexist)

!-- general parameters
 time  = 0.
 polyk = 1.e-10
 gamma = 5./3.
 ieos  = 2

!-- set tmax and dtmax if no infile found, otherwise we use whatever values it had
 if (.not.iexist) then
    tmax      = 20000.
    dtmax     = 100.
    nfulldump = 1
 endif

 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.

!-- set units
 call set_units(mass=1.e6*solarm,c=1.,G=1.)
 mstar = 1.*solarm/umass
 rstar = 1.*solarr/udist

 !
 !-- Read runtime parameters from setup file
 !
  if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' SR polytrope'
  filename = trim(fileprefix)//'.setup'
  inquire(file=filename,exist=iexist)
  if (iexist) call read_setupfile(filename,ierr)
  if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call prompt('Resolution -- number of radial particles',nr,0)
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
  endif

!-- resolution
 psep  = rstar/nr

!-- polytrope
 call rho_polytrope(gamma,polyk,mstar,rtab,rhotab,npts,set_polyk=.true.,Rstar=rstar)
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,xyz_origin=(/0.,0.,0./),rhotab=rhotab(1:npts),rtab=rtab(1:npts))

!-- mass and number of gas particles
 npartoftype(igas) = npart
 massoftype(igas)  = mstar/npart

!-- set thermal energy from density
 do i=1,npart
    call set_particle_type(i,igas)
    densi        = rhoh(xyzh(4,i),massoftype(igas))
    vxyzu(4,i)   = polyk*densi**(gamma-1.) / (gamma-1.)
 enddo

 if (id==master) print "(/,a,i10,/)",' Number of particles setup = ',npart
 if (id==master) print*,' polyk = ',polyk
 if (id==master) print*,' mstar = ',mstar
 if (npart == 0) call fatal('setup','no particles setup')

end subroutine setpart

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for SR polytrope setup'
 call write_inopt(nr,'nr','resolution (number of radial particles)',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(nr,'nr',db,min=0,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
