!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!    Setup of two sink particles in a binary
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    a            -- semi-major axis (+ve) or period (-ve)
!    eccentricity -- eccentricity
!    hacc1        -- primary accretion radius
!    hacc2        -- secondary accretion radius
!    m1           -- mass of primary
!    massr        -- mass ratio
!
!  DEPENDENCIES: externalforces, infile_utils, io, options, part, physcon,
!    setbinary, units
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 real :: massr,m1,a,hacc1,hacc2,ecc

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft
 use setbinary, only:set_binary,get_a_from_period
 use units,     only:set_units
 use physcon,   only:solarm,au,pi
 use options,   only:iexternalforce
 use externalforces, only:iext_corotate,omega_corotate
 use io,        only:master
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
 integer :: ierr
 logical :: iexist
!
!--units
!
! call set_units(mass=solarm,dist=au,G=1.d0)
 call set_units(mass=solarm,G=1.d0,c=1.d0)
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
 massoftype = 1d-9

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 m1    = 1.
 massr = 1.
 a     = 1. !-2.*pi ! 70 AU
 ecc   = 0.7
 hacc1  = 2. !17 ! 1.7 AU
 hacc2  = 2 !17

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Welcome to the Ultimate Binary Setup'
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
 if (a < 0.) a = get_a_from_period(m1,massr*m1,abs(a))

 if (iexternalforce==iext_corotate) then
    call set_binary(m1,massr,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass,omega_corotate)
 else
    call set_binary(m1,massr,a,ecc,hacc1,hacc2,xyzmh_ptmass,vxyz_ptmass,nptmass)
 endif

end subroutine setpart

subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'
 call write_inopt(m1,'m1','mass of primary',iunit)
 call write_inopt(massr,'massr','mass ratio',iunit)
 call write_inopt(a,'a','semi-major axis (+ve) or period (-ve)',iunit)
 call write_inopt(ecc,'ecc','eccentricity',iunit)
 call write_inopt(hacc1,'hacc1','primary accretion radius',iunit)
 call write_inopt(hacc2,'hacc2','secondary accretion radius',iunit)
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
 call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
 call read_inopt(massr,'massr',db,min=0.,errcount=nerr)
 call read_inopt(a,'a',db,errcount=nerr)
 call read_inopt(ecc,'ecc',db,min=0.,errcount=nerr)
 call read_inopt(hacc1,'hacc1',db,min=0.,errcount=nerr)
 call read_inopt(hacc2,'hacc2',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
