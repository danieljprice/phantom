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
!  RUNTIME PARAMETERS:
!     mhole         ---  mass of black hole
!     mstar         ---  mass of star
!     rstar         ---  radius of star
!     beta          ---  penetration factor
!     ecc           ---  eccentricity
!     norbits       ---  number of orbits
!     dumpsperorbit ---  number of dumps per orbit
!     nr            ---  particles per star radius (i.e. resolution)
!
!  DEPENDENCIES: infile_utils, io, part, physcon, setbinary, spherical,
!    timestep, units, metric, eos
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 real    :: mhole,mstar,rstar,beta,ecc,norbits
 integer :: dumpsperorbit,nr

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
 use metric,    only:mass1,a
 use eos,       only:ieos
 use externalforces,only:accradius1,accradius1_hard
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
 real    :: rtidal,rp,semia,psep,period,hacc1,hacc2,massr
 real    :: vxyzstar(3),xyzstar(3)

!
!-- general parameters
!
 time  = 0.
 polyk = 1.e-10    ! <== uconst
 gamma = 5./3.
 ieos  = 4

!
!-- space available for injected gas particles
!
 npart          = 0
 npartoftype(:) = 0
 xyzh(:,:)      = 0.
 vxyzu(:,:)     = 0.
 nptmass        = 0

!
!-- Default runtime parameters
!
!
 mhole         = 1.e6  ! (solar masses)
 mstar         = 1.    ! (solar masses)
 rstar         = 1.    ! (solar radii)
 beta          = 5.
 ecc           = 0.8
 norbits       = 5.
 dumpsperorbit = 100
 nr            = 50

!
!-- Read runtime parameters from setup file
!
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Tidal disruption in GR'
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
!-- Convert to code untis
!
 mhole = mhole*solarm
 call set_units(mass=mhole,c=1.,G=1.) !--Set central mass to M=1 in code units
 mstar         = mstar*solarm/umass
 rstar         = rstar*solarr/udist

 rtidal          = rstar*(mass1/mstar)**(1./3.)
 rp              = rtidal/beta
 semia           = rp/(1.-ecc)
 psep            = rstar/nr
 accradius1_hard = 5.*mass1
 accradius1      = accradius1_hard
 a               = 0.
 period          = 2.*pi*sqrt(semia**3/mass1)
 tmax            = norbits*period
 dtmax           = period/dumpsperorbit

!
!-- Set a binary orbit given the desired orbital parameters
!
 hacc1    = rstar/1.e8    ! Something small so that set_binary doesnt warn about Roche lobe
 hacc2    = hacc1
 massr    = mstar/mass1
 call set_binary(mass1,massr,semia,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass)
 vxyzstar = vxyz_ptmass(1:3,2)
 xyzstar  = xyzmh_ptmass(1:3,2)

!
!-- Delete sinks, replace central sink with Kerr metric and a collection of gas particles
!
 nptmass = 0
 call set_sphere('cubic',id,master,0.,rstar,psep,hfact,npart,xyzh,xyz_origin=xyzstar)
 npartoftype(igas) = npart
 massoftype(igas)  = mstar/npart
 do i=1,npart
    call set_particle_type(i,igas)
    vxyzu(1:3,i) = vxyzstar(1:3)
 enddo

 if (npart == 0)   call fatal('setup','no particles setup')
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
 call write_inopt(mhole,        'mhole',        'mass of black hole (solar mass)',            iunit)
 call write_inopt(mstar,        'mstar',        'mass of star       (solar mass)',            iunit)
 call write_inopt(rstar,        'rstar',        'radius of star     (solar radii)',           iunit)
 call write_inopt(beta,         'beta',         'penetration factor',                         iunit)
 call write_inopt(ecc,          'ecc',          'eccentricity',                               iunit)
 call write_inopt(norbits,      'norbits',      'number of orbits',                           iunit)
 call write_inopt(dumpsperorbit,'dumpsperorbit','number of dumps per orbit',                  iunit)
 call write_inopt(nr           ,'nr'           ,'particles per star radius (i.e. resolution)',iunit)
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
 call read_inopt(mhole,        'mhole',        db,min=0.,errcount=nerr)
 call read_inopt(mstar,        'mstar',        db,min=0.,errcount=nerr)
 call read_inopt(rstar,        'rstar',        db,min=0.,errcount=nerr)
 call read_inopt(beta,         'beta',         db,min=0.,errcount=nerr)
 call read_inopt(ecc,          'ecc',          db,min=0.,errcount=nerr)
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
