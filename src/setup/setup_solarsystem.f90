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
 use mpc,          only:read_mpc,mpc_entry
 use datautils,    only:find_datafile
 use centreofmass, only:reset_centreofmass
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(inout) :: massoftype(:)
 real,              intent(inout) :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=120) :: filename
 integer :: ierr,nbodies,i,j,n,nsample
 logical :: iexist
 real    :: period,semia,mtot,hpart
 integer, parameter :: max_bodies = 2000000
 type(mpc_entry), allocatable :: dat(:)
!
! default runtime parameters
!
 norbits       = 1000.
 dumpsperorbit = 1
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
!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 semia  = 1.  !  Earth
 mtot   = solarm/umass
 hpart  = 10.*au/udist

 period = 2.*pi*sqrt(semia**3/mtot)
 tmax   = norbits*period
 dtmax  = period/dumpsperorbit

 nbodies = 0
 filename = find_datafile('Distant.txt',url='https://www.minorplanetcenter.net/iau/MPCORB/')
 call read_mpc(filename,nbodies,dat=dat)
 print "(a,i0,a)",' read orbital data for ',nbodies,' minor planets'

 nbodies = 0
 n = 0
 nsample = 1  ! can place many particles evenly sampling the orbit if desired
 do i=1,nbodies
    !
    ! for each solar system object get the xyz positions from the orbital parameters
    !
    !print*,i,'aeiOwM=',dat(i)%a,dat(i)%ecc,dat(i)%inc,dat(i)%O,dat(i)%w,dat(i)%M
    do j=1,nsample
       n = n + 1
       if (nsample==1) then
          call set_binary(mtot,epsilon(0.),dat(i)%a,dat(i)%ecc,0.02,1.e-15,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,incl=dat(i)%inc,&
                    arg_peri=dat(i)%w,posang_ascnode=dat(i)%O,&
                    mean_anomaly=dat(i)%M,verbose=.false.)
       else
          call set_binary(mtot,epsilon(0.),dat(i)%a,dat(i)%ecc,0.02,1.e-15,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,incl=dat(i)%inc,&
                    arg_peri=dat(i)%w,posang_ascnode=dat(i)%O,&
                    f=360.*(n-1)/nsample,verbose=.false.)
       endif
       !
       ! now delete the point masses but set a dust particle as the secondary
       !
       nptmass = 0
       xyzh(1:3,n)  = xyzmh_ptmass(1:3,2)
       xyzh(4,n)    = hpart  ! give a random length scale as the smoothing length
       vxyzu(1:3,n) = vxyz_ptmass(1:3,2)
       call set_particle_type(n,idust)
    enddo

 enddo
 !
 ! add the Sun
 !
 nptmass = 1
 xyzmh_ptmass(:,1) = 0.
 xyzmh_ptmass(4,1) = mtot
 xyzmh_ptmass(ihacc,1) = solarr/udist
 !
 ! add the planets
 !
 call add_body('jupiter',nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
 !call set_solarsystem_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
 !
 ! set mass of all the minor bodies equal
 !
 print "(a)"
 npart = nbodies*nsample
 ndustlarge = 1
 ndustsmall = 0
 ndusttypes = 1
 npartoftype(idust) = nbodies*nsample
 massoftype(idust) = 1.e-20
 grainsize(1:ndustlarge) = km/udist         ! assume km-sized bodies
 graindens(1:ndustlarge) = 2./unit_density  ! 2 g/cm^3


 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 hfact = 1.2

 if (ierr /= 0) call fatal('setup','ERROR during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  setup the solar system planets by querying their ephemeris
!  from the JPL server
!+
!----------------------------------------------------------------
subroutine set_solarsystem_planets(nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
 integer, intent(inout) :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,    intent(in)    :: mtot
 integer,          parameter :: nplanets = 9
 character(len=*), parameter :: planet_name(nplanets) = &
     (/'mercury', &
       'venus  ', &
       'earth  ', &
       'mars   ', &
       'jupiter', &
       'saturn ', &
       'uranus ', &
       'neptune', &
       'pluto  '/)  ! for nostalgia's sake
 integer :: i

 do i=1,nplanets
    call add_body(planet_name(i),nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
 enddo

 print*,' nptmass = ',nptmass

end subroutine set_solarsystem_planets

!----------------------------------------------------------------
!+
!  setup a body in the solar system by querying their ephemeris
!  from the JPL server
!+
!----------------------------------------------------------------
subroutine add_body(body_name,nptmass,xyzmh_ptmass,vxyz_ptmass,mtot)
 use ephemeris, only:get_ephemeris,nelem
 use units,     only:umass,udist
 use physcon,   only:gg,km,solarm,earthm,au
 use setbinary, only:set_binary
 character(len=*), intent(in)    :: body_name
 integer,          intent(inout) :: nptmass
 real,             intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,             intent(in)    :: mtot
 real    :: elems(nelem)
 real    :: xyz_tmp(size(xyzmh_ptmass(:,1)),2),vxyz_tmp(3,2),gm_cgs
 real    :: mbody,rbody,a,e,inc,O,w,f
 integer :: ierr,ntmp

 elems = get_ephemeris(body_name,ierr)
 if (ierr > 0) then
    print "(a)",' ERROR: could not read ephemeris data for '//body_name
    !return ! skip if error reading ephemeris file
 endif
 a   = elems(1)*km/udist
 e   = elems(2)
 inc = elems(3)
 O   = elems(4)
 w   = elems(5)
 f   = elems(6)
 gm_cgs = elems(7)*km**3
 mbody  = (gm_cgs/gg)/umass
 rbody = elems(8)*km/udist
 print "(1x,a,1pg10.3)",   '           m/msun = ',mbody*umass/solarm
 print "(1x,a,1pg10.3)",   '         m/mearth = ',mbody*umass/earthm
 print "(1x,a,1pg10.4,a)", '                a = ',a*udist/au,' au'
 print "(1x,a,1pg10.3,a)", '           radius = ',elems(8),' km'
 print "(1x,a,1pg10.3,a)", '          density = ',elems(9),' g/cm^3'
 ntmp = 0
 call set_binary(mtot,0.,a,e,0.01,rbody,&
      xyz_tmp,vxyz_tmp,ntmp,ierr,incl=inc,&
      arg_peri=w,posang_ascnode=O,f=f,verbose=.false.)

 ! add a point mass for each body
 nptmass = nptmass + 1
 xyzmh_ptmass(1:3,nptmass) = xyz_tmp(:,2)
 xyzmh_ptmass(4,nptmass) = mbody
 xyzmh_ptmass(5,nptmass) = rbody
 vxyz_ptmass(1:3,nptmass) = vxyz_tmp(:,2)
 print "(1x,a,3(1x,1pg10.3))",' x = ',xyz_tmp(1:3,2)
 print "(1x,a,3(1x,1pg10.3))",' v = ',vxyz_tmp(1:3,2)

end subroutine add_body

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
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
