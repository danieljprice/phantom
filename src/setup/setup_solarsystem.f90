!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup asteroid orbits using data from the IAU Minor Planet Center
!
! :References: https://minorplanetcenter.net/data
!
! :Owner: Not Committed Yet
!
! :Runtime parameters:
!   - dumpsperorbit : *number of dumps per orbit*
!   - norbits       : *number of orbits*
!
! :Dependencies: infile_utils, io, mpc, part, physcon, setbinary, timestep,
!   units
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
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,idust,ndustlarge,set_particle_type,&
                     grainsize,graindens,ndusttypes
 use setbinary, only:set_binary
 use units,     only:set_units,umass,udist,unit_density
 use physcon,   only:solarm,au,pi,km
 use io,        only:master,fatal
 use timestep,  only:tmax,dtmax
 use mpc,       only:read_mpc,mpc_entry
 use datautils, only:find_datafile
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
 integer :: ierr,nbodies,i
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

 semia  = 5.2  !  Jupiter
 mtot   = solarm/umass
 hpart  = 100.*au/udist

 period = 2.*pi*sqrt(semia**3/mtot)
 tmax   = norbits*period
 dtmax  = period/dumpsperorbit
!
! read the orbital data from the data file
!
 filename = find_datafile('Distant.txt',url='https://www.minorplanetcenter.net/iau/MPCORB/')
 call read_mpc(filename,nbodies,dat=dat)
 print "(a,i0,a)",' read orbital data for ',nbodies,' minor planets'

 do i=1,nbodies
    !
    ! for each solar system object get the xyz positions from the orbital parameters
    !
    !print*,i,'aeiOwM=',dat(i)%a,dat(i)%ecc,dat(i)%inc,dat(i)%O,dat(i)%w,dat(i)%M
    call set_binary(mtot,epsilon(0.),dat(i)%a,dat(i)%ecc,1.,1.e-15,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,incl=dat(i)%inc,&
                    arg_peri=dat(i)%w,posang_ascnode=dat(i)%O,&
                    mean_anomaly=dat(i)%M,verbose=.false.)
    !
    ! now delete the point masses but set a dust particle as the secondary
    !
    nptmass = 0
    xyzh(1:3,i)  = xyzmh_ptmass(1:3,2)
    xyzh(4,i)    = hpart  ! give a random length scale as the smoothing length
    vxyzu(1:3,i) = vxyz_ptmass(1:3,2)
    call set_particle_type(i,idust)
 enddo
 !
 ! restore the Sun
 !
 nptmass = 1
 !
 ! set mass of all the minor bodies equal
 !
 npart = nbodies
 ndustlarge = 1
 ndusttypes = 1
 npartoftype(idust) = nbodies
 massoftype(idust) = 1.e-20
 grainsize(1:ndustlarge) = km/udist         ! assume km-sized bodies
 graindens(1:ndustlarge) = 2./unit_density  ! 2 g/cm^3

 hfact = 1.2

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
