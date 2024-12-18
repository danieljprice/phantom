!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup asteroid orbits using data from the IAU Minor Planet Center
!
! :References: https://minorplanetcenter.net/data
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - dumpsperorbit : *number of dumps per orbit*
!   - norbits       : *number of orbits*
!
! :Dependencies: centreofmass, datautils, ephemeris, infile_utils, io, mpc,
!   part, physcon, setbinary, timestep, units
!
 implicit none
 public :: setpart

 real :: norbits
 integer :: dumpsperorbit
 logical :: asteroids
 character(len=20) :: epoch

 private

contains
!----------------------------------------------------------------
!+
!  setup for solar system orbits
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,idust,set_particle_type,&
                        grainsize,graindens,ndustlarge,ndusttypes,ndustsmall,ihacc
 use setbinary,    only:set_binary
 use units,        only:set_units,umass,udist,unit_density
 use physcon,      only:solarm,au,pi,km,solarr
 use io,           only:master,fatal
 use timestep,     only:tmax,dtmax
 use centreofmass, only:reset_centreofmass
 use setbodies,    only:set_minor_planets,add_sun_and_planets,add_body
 use kernel,       only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout) :: massoftype(:)
 real,              intent(inout) :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: ierr
 integer :: values(8),year,month,day
 logical :: iexist
 real    :: period,semia,mtot
 character(len=120) :: filename
!
! default runtime parameters
!
 norbits       = 1000.
 dumpsperorbit = 1
 asteroids = .true.

 call date_and_time(values=values)
 year = values(1); month = values(2); day = values(3)
 write(epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day
!
! read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Superb Solar System Setup'

 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif
!
! set units
!
 call set_units(mass=solarm,dist=au,G=1.d0)
!
! general parameters
!
 time  = 0.
 polyk = 0.
 gamma = 1.
 hfact = hfact_default
!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 semia  = 1.*au/udist  !  Earth
 mtot   = solarm/umass !  mass around which all bodies should orbit

 period = 2.*pi*sqrt(semia**3/mtot)
 tmax   = norbits*period
 dtmax  = period/dumpsperorbit

 if (asteroids) then
    call set_minor_planets(npart,npartoftype,massoftype,xyzh,vxyzu,&
                           mtot,itype=idust,sample_orbits=.false.)
    print*,'npart = ',npart,' npartoftype = ',npartoftype(idust)
 endif
 !
 ! treat minor bodies as km-sized dust particles
 !
 ndustlarge = 1
 ndustsmall = 0
 ndusttypes = 1
 grainsize(ndustlarge) = km/udist         ! assume km-sized bodies
 graindens(ndustlarge) = 2./unit_density  ! 2 g/cm^3
 !
 ! add the planets
 !
 call add_sun_and_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 call add_body('apophis',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 !call add_body('sun',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 !call add_body('jupiter',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (ierr /= 0) call fatal('setup','ERROR during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  write setup parameters to file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 write(iunit,"(a)") '# input file for solar system setup routines'
 call write_inopt(norbits,'norbits','number of orbits',iunit)
 call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',iunit)
 call write_inopt(asteroids,'asteroids','add all distant minor bodies as km-sized dust particles',iunit)
 call write_inopt(epoch,'epoch','epoch to query ephemeris, YYYY-MM-DD, blank = today',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read setup parameters from file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(norbits,      'norbits',      db,min=0.,errcount=nerr)
 call read_inopt(dumpsperorbit,'dumpsperorbit',db,min=0 ,errcount=nerr)
 call read_inopt(asteroids,'asteroids',db,errcount=nerr)
 call read_inopt(epoch,'epoch',db,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
