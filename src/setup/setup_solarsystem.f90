!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
!   - asteroids  : *add distant minor bodies as km-sized dust particles*
!   - dtmax_in   : *time between dumps (e.g. 1 hr)*
!   - epoch      : *epoch to query ephemeris, YYYY-MMM-DD HH:MM:SS.fff, blank = today*
!   - np_apophis : *number of particles used to represent apophis (0=none; 1=sink; n=gas)*
!   - tmax_in    : *end time of simulation (e.g. 3 days)*
!
! :Dependencies: centreofmass, eos_tillotson, infile_utils, io, kernel,
!   options, part, physcon, setbinary, setsolarsystem, setup_params,
!   spherical, timestep, units
!
 implicit none
 public :: setpart

 integer :: np_apophis
 logical :: asteroids
 character(len=20) :: epoch,tmax_in,dtmax_in

 private

contains
!----------------------------------------------------------------
!+
!  setup for solar system orbits
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,         only:nptmass,xyzmh_ptmass,vxyz_ptmass,idust,set_particle_type,&
                        grainsize,graindens,ndustlarge,ndusttypes,ndustsmall,ihacc,igas
 use setbinary,     only:set_binary
 use units,         only:set_units,umass,udist,unit_density,unit_velocity,utime,in_code_units
 use physcon,       only:solarm,au,pi,km,solarr,ceresm,earthm,earthr,days
 use io,            only:master,fatal
 use timestep,      only:tmax,dtmax
 use centreofmass,  only:reset_centreofmass
 use setsolarsystem,only:set_minor_planets,add_sun_and_planets,add_body
 use kernel,        only:hfact_default
 use eos_tillotson, only:rho_0,A
 use spherical,     only:set_sphere
 use options,       only:ieos
 use setup_params,  only:npart_total
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout) :: massoftype(:)
 real,              intent(inout) :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: ierr,i
 !integer :: values(8),year,month,day
 logical :: iexist
 real    :: period,semia,mtot,dx
 real    :: r_apophis,m_apophis,rtidal,spsoundmin
 character(len=120) :: filename
!
! default runtime parameters
!
 tmax_in = '1000 yr'
 dtmax_in = '1 yr'
 asteroids = .true.
 np_apophis = 0
 !call date_and_time(values=values)
 !year = values(1); month = values(2); day = values(3)
 !write(epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day
 epoch='2029-04-10'   ! encounter is on Friday 13th April
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
 npart_total = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 semia  = 1.*au/udist  !  Earth
 mtot   = solarm/umass !  mass around which all bodies should orbit

 period = 2.*pi*sqrt(semia**3/mtot)
 tmax   = in_code_units(tmax_in,ierr,unit_type='time')
 if (ierr /= 0) call fatal('setup_solarsystem',' could not parse tmax')
 dtmax  = in_code_units(dtmax_in,ierr,unit_type='time')
 if (ierr /= 0) call fatal('setup_solarsystem',' could not parse dtmax')

 if (asteroids) then
    call set_minor_planets(npart,npartoftype,massoftype,xyzh,vxyzu,&
                           mtot,itype=idust,sample_orbits=.false.)
    print*,'npart = ',npart,' npartoftype = ',npartoftype(idust)
    !
    ! treat minor bodies as km-sized dust particles
    !
    ndustlarge = 1
    ndustsmall = 0
    ndusttypes = 1
    grainsize(ndustlarge) = km/udist         ! assume km-sized bodies
    graindens(ndustlarge) = 2./unit_density  ! 2 g/cm^3
 endif
 !
 ! add the planets
 !
 call add_sun_and_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)
 !
 ! add the bringer of death
 !
 if (np_apophis > 0) then
    call add_body('apophis',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,epoch)

    r_apophis = xyzmh_ptmass(5,nptmass)
    m_apophis = 4./3.*pi*(rho_0/unit_density)*r_apophis**3
    xyzmh_ptmass(4,nptmass) = m_apophis

    print "(a,2(es10.3,a))",' mass of apophis is ',m_apophis*umass,&
                            ' g or ',m_apophis*umass/ceresm,' ceres masses'
    print "(a,1pg10.3,a)",' density is ',m_apophis/(4./3.*pi*r_apophis**3)*unit_density,' g/cm^3'

    rtidal = r_apophis*(earthm/umass/m_apophis)**(1./3.)
    print "(3(a,1pg10.3),a)",' r_tidal is ',rtidal,' au,',rtidal*udist/km,' km, or ',rtidal*udist/earthr,' earth radii'

    if (np_apophis > 1) then
       !
       ! replace the sink particle with a ball of stuff
       !
       dx = r_apophis/40.
       call set_sphere('closepacked',id,master,0.,r_apophis,dx,hfact,npart,xyzh,npart_total,&
                       xyz_origin=xyzmh_ptmass(1:3,nptmass),exactN=.true.,np_requested=np_apophis)

       do i=1,npart
          vxyzu(1:3,i) = vxyz_ptmass(1:3,nptmass)
       enddo
       massoftype(igas) = m_apophis / npart
       npartoftype(igas) = npart
       nptmass = nptmass - 1
       !
       ! print quantities from the equation of state to give an idea of the timestep
       !
       if (ieos==23) then
          spsoundmin = sqrt(A/rho_0)/unit_velocity
          print "(a,1pg11.4,a)",'     sound speed min = ',spsoundmin*unit_velocity/km,' km/s'
          print "(a,1pg10.3,a)",' sound crossing time = ',(r_apophis/spsoundmin)*utime,' seconds'
       endif
    endif
 endif
 !
 ! set centre of mass as the origin
 !
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
 call write_inopt(tmax_in,'tmax_in','end time of simulation (e.g. 3 days)',iunit)
 call write_inopt(dtmax_in,'dtmax_in','time between dumps (e.g. 1 hr)',iunit)
 call write_inopt(asteroids,'asteroids','add distant minor bodies as km-sized dust particles',iunit)
 call write_inopt(np_apophis,'np_apophis','number of particles used to represent apophis (0=none; 1=sink; n=gas)',iunit)
 call write_inopt(epoch,'epoch','epoch to query ephemeris, YYYY-MMM-DD HH:MM:SS.fff, blank = today',iunit)
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
 call read_inopt(tmax_in, 'tmax_in',db,errcount=nerr)
 call read_inopt(dtmax_in,'dtmax_in',db,errcount=nerr)
 call read_inopt(asteroids,'asteroids',db,errcount=nerr)
 call read_inopt(np_apophis,'np_apophis',db,errcount=nerr)
 call read_inopt(epoch,'epoch',db,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
