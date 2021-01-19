!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup of an asteroid being tidally disrupted by a white dwarf
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters:
!   - dumpsperorbit : *number of dumps per orbit*
!   - hacc1         : *white dwarf (sink) accretion radius (solar radii)*
!   - m1            : *mass of white dwarf (solar mass)*
!   - m2            : *mass of asteroid (ceres mass)*
!   - norbits       : *number of orbits*
!   - nr            : *particles per asteroid radius (i.e. resolution)*
!   - rasteroid     : *radius of asteroid (km)*
!   - rp            : *pericentre distance (solar radii)*
!   - semia         : *semi-major axis (solar radii)*
!
! :Dependencies: infile_utils, io, part, physcon, setbinary, spherical,
!   timestep, units
!
 implicit none
 public :: setpart

 real :: m1,m2,rp,semia,hacc1,rasteroid,norbits
 integer :: dumpsperorbit,nr

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,idust,set_particle_type
 use setbinary, only:set_binary,get_a_from_period
 use spherical, only:set_sphere
 use units,     only:set_units,umass,udist
 use physcon,   only:solarm,au,pi,solarr,ceresm,km
 use io,        only:master,fatal
 use timestep,  only:tmax,dtmax
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
 integer :: ierr,i
 logical :: iexist
 real    :: vbody(3),xyzbody(3),massbody,psep,period,hacc2,ecc

 call set_units(mass=solarm,dist=solarr,G=1.d0)

!
!--Default runtime parameters
!
 m1            = 0.6  ! (solar masses)
 m2            = 0.1  ! (ceres masses)
 rp            = 1.2  ! (solar radii)
 semia         = 24.  ! (solar radii)
 hacc1         = 0.5  ! (solar radii)
 rasteroid     = 100. ! (km)
 norbits       = 1.
 dumpsperorbit = 100
 nr            = 50

!
!--Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' White Dwarf tidal disruption'
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
 !--Convert to code units
 !
 m1    = m1*solarm/umass
 m2    = m2*ceresm/umass
 rp    = rp*solarr/udist
 semia = semia*solarr/udist
 hacc1 = hacc1*solarr/udist
 rasteroid = rasteroid*km/udist

!
!--general parameters
!
 time = 0.
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

 ecc    = 1.-rp/semia
 period = sqrt(4.*pi**2*semia**3/(m1+m2))
 hacc2  = hacc1/1.e10
 tmax   = norbits*period
 dtmax  = period/dumpsperorbit

!
!--Set a binary orbit given the desired orbital parameters
!
 call set_binary(m1,m2,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
 vbody    = vxyz_ptmass(1:3,2)
 xyzbody  = xyzmh_ptmass(1:3,2)
 massbody = xyzmh_ptmass(4,2)

!
!--Delete second sink and replace with collection of dust particles
!
 nptmass = 1
 psep  = rasteroid/nr
 call set_sphere('cubic',id,master,0.,rasteroid,psep,hfact,npart,xyzh,xyz_origin=xyzbody)
 if (id==master) print "(1(/,a,i10,a,/))",' Replaced second sink with ',npart,' dust particles'
 npartoftype(idust) = npart
 massoftype(idust)  = massbody/npart
 do i=1,npart
    call set_particle_type(i,idust)
    vxyzu(1:3,i) = vbody(1:3)
 enddo

 if (nptmass == 0) call fatal('setup','no sink particles setup')
 if (npart == 0)   call fatal('setup','no hydro particles setup')
 if (ierr /= 0)    call fatal('setup','ERROR during setup')

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
 write(iunit,"(a)") '# input file for binary setup routines'
 call write_inopt(m1,           'm1',           'mass of white dwarf (solar mass)',                 iunit)
 call write_inopt(m2,           'm2',           'mass of asteroid (ceres mass)',                    iunit)
 call write_inopt(rp,           'rp',           'pericentre distance (solar radii)',                iunit)
 call write_inopt(semia,        'semia',        'semi-major axis (solar radii)',                    iunit)
 call write_inopt(hacc1,        'hacc1',        'white dwarf (sink) accretion radius (solar radii)',iunit)
 call write_inopt(rasteroid,    'rasteroid',    'radius of asteroid (km)',                          iunit)
 call write_inopt(norbits,      'norbits',      'number of orbits',                                 iunit)
 call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',                        iunit)
 call write_inopt(nr           ,'nr'           ,'particles per asteroid radius (i.e. resolution)',  iunit)
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
 call read_inopt(m1,           'm1',           db,min=0.,errcount=nerr)
 call read_inopt(m2,           'm2',           db,min=0.,errcount=nerr)
 call read_inopt(rp,           'rp',           db,min=0.,errcount=nerr)
 call read_inopt(semia,        'semia',        db,min=0.,errcount=nerr)
 call read_inopt(hacc1,        'hacc1',        db,min=0.,errcount=nerr)
 call read_inopt(rasteroid,    'rasteroid',    db,min=0.,errcount=nerr)
 call read_inopt(norbits,      'norbits',      db,min=0.,errcount=nerr)
 call read_inopt(dumpsperorbit,'dumpsperorbit',db,min=0 ,errcount=nerr)
 call read_inopt(nr,           'nr',           db,min=0 ,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
