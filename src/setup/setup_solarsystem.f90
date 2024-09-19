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
                        grainsize,graindens,ndustlarge,ndusttypes
 use setbinary,    only:set_binary
 use units,        only:set_units,umass,udist,unit_density
 use physcon,      only:solarm,au,pi,km
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
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' solar system'
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

 filename = find_datafile('Distant.txt',url='https://www.minorplanetcenter.net/iau/MPCORB/')
 call read_mpc(filename,nbodies,dat=dat)
 print "(a,i0,a)",' read orbital data for ',nbodies,' minor planets'

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
 ! restore the Sun
 !
 nptmass = 1
 !
 ! set mass of all the minor bodies equal
 !
 npart = nbodies*nsample
 print*,' n = ',n,' npart = ',npart
 ndustlarge = 1
 ndusttypes = 1
 npartoftype(idust) = nbodies*nsample
 massoftype(idust) = 1.e-20
 grainsize(1:ndustlarge) = km/udist         ! assume km-sized bodies
 graindens(1:ndustlarge) = 2./unit_density  ! 2 g/cm^3

 !
 ! add the planets
 !
 call set_solarsystem_planets(nptmass,xyzmh_ptmass,vxyz_ptmass)

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
subroutine set_solarsystem_planets(nptmass,xyzmh_ptmass,vxyz_ptmass)
 use ephemeris, only:get_ephemeris,nelem
 use units,     only:umass,udist
 use physcon,   only:gg,km,solarm,earthm,au
 use setbinary, only:set_binary
 integer, intent(inout) :: nptmass
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
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
 real    :: elems(nelem),xyz_tmp(size(xyzmh_ptmass(:,1)),2),vxyz_tmp(3,2),gm_cgs
 real    :: msun,mplanet,a,e,inc,O,w,f
 integer :: i,ierr,ntmp

 msun = solarm/umass
 do i=1,nplanets
    elems = get_ephemeris(planet_name(i),ierr)
    if (ierr /= 0) then
       print "(a)",' ERROR: could not read ephemeris data for '//planet_name(i)
       cycle  ! skip if error reading ephemeris file
    endif
    gm_cgs  = elems(1)*km**3
    mplanet = (gm_cgs/gg)/umass
    a   = elems(2)*km/udist
    e   = elems(3)
    inc = elems(4)
    O   = elems(5)
    w   = elems(6)
    f   = elems(7)
    print*,' mplanet/mearth = ',mplanet*umass/earthm,' a = ',a*udist/au,' au'
    ntmp = 0
    call set_binary(msun,mplanet,a,e,0.01,0.01,&
                    xyz_tmp,vxyz_tmp,ntmp,ierr,incl=inc,&
                    arg_peri=w,posang_ascnode=O,f=f,verbose=.false.)
    nptmass = nptmass + 1
    xyzmh_ptmass(:,nptmass) = xyz_tmp(:,2)
    vxyz_ptmass(:,nptmass)  = vxyz_tmp(:,2)
 enddo

 print*,' nptmass = ',nptmass

end subroutine set_solarsystem_planets

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
