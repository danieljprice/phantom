!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
 logical :: use_dem,apophis_only
character(len=256) :: apophis_shape_file

 real :: scale_vel
 real :: scale_pos
 real :: scale_r_apophis
 real :: scale_rho

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
 use units,         only:set_units,umass,udist,unit_density,unit_velocity,utime,in_code_units,in_units
 use physcon,       only:solarm,pi,au,km,solarr,ceresm,earthm,earthr,days
 use io,            only:master,fatal,warning
 use timestep,      only:tmax,dtmax
 use centreofmass,  only:reset_centreofmass
 use setsolarsystem,only:set_minor_planets,add_sun_and_planets,add_body
 use kernel,        only:hfact_default
 use eos_tillotson, only:rho_0,A
 use shape,         only:set_shape
 use options,       only:ieos
 use setup_params,  only:npart_total
 use infile_utils,  only:get_options
 use ptmass,        only:isink_potential
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout) :: massoftype(:)
 real,              intent(inout) :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: ierr,i,nerr
 !integer :: values(8),year,month,day
 real    :: period,semia,mtot,dx
 real    :: r_apophis,m_apophis,rtidal,spsoundmin
!
! default runtime parameters
!
 tmax_in = '1000 yr'
 dtmax_in = '1 yr'
 asteroids = .true.
 np_apophis = 0
 use_dem = .false.
 apophis_only = .false.
 !call date_and_time(values=values)
 !year = values(1); month = values(2); day = values(3)
 !write(epoch,"(i4.4,'-',i2.2,'-',i2.2)") year,month,day
 epoch='2029-04-10'   ! encounter is on Friday 13th April
 scale_vel=1.
 scale_pos=1.
 scale_r_apophis=1.
 scale_rho=1.
 apophis_shape_file='apophis.shape'
!
! read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Superb Solar System Setup'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
!
! set units
!
 call set_units(mass=solarm,dist=km,G=1.d0)
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
 ierr = 0
 call add_sun_and_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
 if (nerr > 0) ierr = ierr + nerr

 if (apophis_only) nptmass = 0
 !
 ! add the bringer of death
 !
 if (np_apophis > 0) then
    call add_body('apophis',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot,nerr,epoch)
    if (nerr > 0) call warning('apophis','missing some information')

    r_apophis = xyzmh_ptmass(5,nptmass) * scale_r_apophis
    xyzmh_ptmass(5,nptmass) = r_apophis
    print "(a,1pg10.3)",' apophis radius scaled by ',scale_r_apophis

    m_apophis = 4./3.*pi*(rho_0*scale_rho/unit_density)*r_apophis**3
    xyzmh_ptmass(4,nptmass) = m_apophis
    print "(a,2(es10.3,a))",' mass of apophis is ',m_apophis*umass,&
                            ' g or ',m_apophis*umass/ceresm,' ceres masses'
    print "(a,1pg10.3,a)",' density is ',m_apophis/(4./3.*pi*r_apophis**3)*unit_density,' g/cm^3'
    print "(a,1pg10.3,a)",' Apophis-Earth relative velocity = ',&
                            in_units(sqrt(sum((vxyz_ptmass(1:3,nptmass)-vxyz_ptmass(1:3,4))**2)),'km/s'),' km/s'

    rtidal = r_apophis*(earthm/umass/m_apophis)**(1./3.)
    print "(3(a,1pg10.3),a)",' r_tidal is ',rtidal,' au,',rtidal*udist/km,' km, or ',rtidal*udist/earthr,' earth radii'

    vxyz_ptmass(1:3,nptmass) = vxyz_ptmass(1:3,nptmass)*scale_vel
    print "(a,1pg10.3)",' velocity of apophis scaled by ',scale_vel

    xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass)*scale_pos
    print "(a,1pg10.3)",' initial position of apophis scaled by ',scale_pos


    if (np_apophis > 1) then
       !
       ! replace the sink particle with a ball of stuff
       !
       call set_shape('closepacked',id,master,np_apophis,xyzmh_ptmass(1:3,nptmass),r_apophis,&
                      hfact,npart,xyzh,npart_total,objfile=apophis_shape_file)
       !call set_sphere('closepacked',id,master,0.,r_apophis,dx,hfact,npart,xyzh,npart_total,&
       !                xyz_origin=xyzmh_ptmass(1:3,nptmass),exactN=.true.,np_requested=np_apophis)

       do i=1,npart
          vxyzu(1:3,i) = vxyz_ptmass(1:3,nptmass)
       enddo
       massoftype(igas) = m_apophis / npart
       npartoftype(igas) = npart
       nptmass = nptmass - 1

       if (use_dem) then
          call replace_gas_with_dem(id,npart,npartoftype(igas),massoftype(igas),&
                                    xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,hfact)
          isink_potential = 2
       endif
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

 if (ierr /= 0) call fatal('setup','ERRORS during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  replace gas with discrete element method particles
!+
!----------------------------------------------------------------
subroutine replace_gas_with_dem(id,npart,ngas,pmass,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass,hfact)
 use part, only:iReff,ihacc
 use units, only:udist
 use physcon, only:km
 use io, only:master
 integer, intent(in)    :: id
 integer, intent(inout) :: npart,ngas,nptmass
 real, intent(inout) :: pmass,xyzh(:,:),vxyzu(:,:)
 real, intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real, intent(in)    :: hfact
 integer :: i,j
 real :: dxij,dyij,dzij,dmin,reff,sep(npart)
 real :: sep_med,sep_min,sep_max

 if (npart < 1) return

 !
 ! DEM sphere radius from cropped lattice geometry (not r_apophis/40).
 ! Use per-particle Reff = 0.5*nearest-neighbour distance so ellipsoid/mesh
 ! surfaces do not share one radius (median) that overlaps close pairs on step 1.
 !
 sep_med = 0.
 sep_min = huge(sep_min)
 sep_max = 0.
 if (npart == 1) then
    sep(1) = xyzh(4,1)/hfact
 else
    do i=1,npart
       dmin = huge(dmin)
       do j=1,npart
          if (j == i) cycle
          dxij = xyzh(1,i) - xyzh(1,j)
          dyij = xyzh(2,i) - xyzh(2,j)
          dzij = xyzh(3,i) - xyzh(3,j)
          dmin = min(dmin,sqrt(dxij*dxij + dyij*dyij + dzij*dzij))
       enddo
       sep(i) = dmin
       sep_min = min(sep_min,dmin)
       sep_max = max(sep_max,dmin)
    enddo
    do i=2,npart
       dmin = sep(i)
       j = i - 1
       do while (j >= 1 .and. sep(j) > dmin)
          sep(j+1) = sep(j)
          j = j - 1
       enddo
       sep(j+1) = dmin
    enddo
    if (mod(npart,2) == 1) then
       sep_med = sep((npart+1)/2)
    else
       sep_med = 0.5*(sep(npart/2) + sep(npart/2+1))
    endif
 endif

 if (id == master) then
    print "(a,1pg12.4,a)",' DEM NN spacing (min/median/max) = ',sep_min*udist/km,&
          ' / ',sep_med*udist/km,' / ',sep_max*udist/km,' km'
 endif

 xyzmh_ptmass(:,nptmass+1:) = 0.
 do i=1,npart
    nptmass = nptmass + 1
    vxyz_ptmass(1:3,nptmass) = vxyzu(1:3,i)
    xyzmh_ptmass(1:3,nptmass) = xyzh(1:3,i)
    xyzmh_ptmass(4,nptmass) = pmass
    reff = 0.5*sep(i)
    xyzmh_ptmass(iReff,nptmass) = reff
    xyzmh_ptmass(ihacc,nptmass) = reff
 enddo
 npart = 0
 ngas = 0
 pmass = 0.

end subroutine replace_gas_with_dem

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

 call write_inopt(use_dem,'use_dem','use the discrete element method for sink-sink interactions',iunit)
 call write_inopt(apophis_only,'apophis_only','only add apophis',iunit)

 call write_inopt(scale_vel,'scale_vel','scaling factor for apophis velocity',iunit)
 call write_inopt(scale_pos,'scale_pos','scaling factor for apophis initial position',iunit)
 call write_inopt(scale_r_apophis,'scale_r_apophis','scaling factor for apophis radius',iunit)
 call write_inopt(scale_rho,'scale_rho','scaling factor for apophis bulk density',iunit)
call write_inopt(apophis_shape_file,'apophis_shape_file','shape config file for lattice cropping',iunit)

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
 call read_inopt(np_apophis,'np_apophis',db,min=0,errcount=nerr)
 call read_inopt(epoch,'epoch',db,errcount=nerr)

 call read_inopt(use_dem,'use_dem',db,errcount=nerr)
 call read_inopt(apophis_only,'apophis_only',db,errcount=nerr)

 call read_inopt(scale_vel,'scale_vel',db,default=1.0,errcount=nerr)
 call read_inopt(scale_pos,'scale_pos',db,default=1.0,errcount=nerr)
 call read_inopt(scale_r_apophis,'scale_r_apophis',db,default=1.0,errcount=nerr)
 call read_inopt(scale_rho,'scale_rho',db,default=1.0,errcount=nerr)
call read_inopt(apophis_shape_file,'apophis_shape_file',db,default='apophis.shape',errcount=nerr)

 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
