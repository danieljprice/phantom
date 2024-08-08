!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of two stars or sink particles in a binary
!
! :References: None
!
! :Owner: Ana Lourdes Juarez
!
! :Runtime parameters:
!   - a    : *semi-major axis*
!   - hacc : *accretion radius of the companion star*
!   - macc : *mass of the companion star*
!   - mdon : *mass of the donor star*
!
! :Dependencies: centreofmass, eos, extern_corotate, externalforces,
!   infile_utils, io, options, part, setbinary, setunits, timestep
!

 implicit none
 public :: setpart
 real    :: a,mdon,macc,hacc

 private

contains

!----------------------------------------------------------------
!+
!  setup for binary star simulations (with or without gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc
 use setbinary,      only:set_binary,get_period_from_a
 use centreofmass,    only:reset_centreofmass
 use options,        only:iexternalforce
 use externalforces, only:iext_corotate,omega_corotate
 use extern_corotate, only:icompanion_grav,companion_xpos,companion_mass,hsoft
 use io,             only:master,fatal
 use eos,            only:ieos, gmw
 use setunits,       only:mass_unit,dist_unit
 use timestep,       only:tmax,dtmax
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
 real    :: period,ecc,hdon,mass_ratio
!
!--general parameters
!
 dist_unit = 'solarr'
 mass_unit = 'solarm'
 time = 0.
 polyk = 0.
 gamma = 5./3.
!
!--space available for injected gas particles
!  in case only sink particles are used
!
 npart = 0
 npartoftype(:) = 0
 massoftype = 0.

 iexternalforce = iext_corotate
 icompanion_grav = 1
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0
 a    = 266.34
 mdon = 6.97
 macc = 1.41
 hacc = 1.
 ieos = 2
 gmw  = 0.6
 ecc  = 0.
 hdon = 1.

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Ultimate Binary Setup'

 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ieos,polyk,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun phantomsetup'
    endif
    stop
 endif
 !
 !
 !--if a is negative or is given time units, interpret this as a period
 !

 period = get_period_from_a(mdon,macc,a)
 tmax = 10.*period
 dtmax = tmax/200.
 !
 !--now setup orbit using fake sink particles
 !
 call set_binary(mdon,macc,a,ecc,hdon,hacc,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,&
                  verbose=(id==master))

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)


 if (ierr /= 0) call fatal ('setup_binary','error in call to set_binary')

 companion_mass = mdon
 companion_xpos = xyzmh_ptmass(1,1)
 mass_ratio = mdon / macc
 hsoft = 0.1 * 0.49 * mass_ratio**(2./3.) / (0.6*mass_ratio**(2./3.) + &
               log( 1. + mass_ratio**(1./3.) ) ) * a
 !
 !--delete donor sink
 !
 nptmass=1
 xyzmh_ptmass(:,1) = xyzmh_ptmass(:,2)
 vxyz_ptmass(1:3,1) = 0.

 !--restore options
 !


end subroutine setpart

!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setunits,     only:write_options_units
 use eos,          only:write_options_eos
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'

 call write_options_units(iunit)
 call write_options_eos(iunit)

 write(iunit,"(/,a)") '# orbit settings'
 call write_inopt(a,'a','semi-major axis',iunit)
 call write_inopt(mdon,'mdon','mass of the donor star',iunit)
 call write_inopt(macc,'macc','mass of the companion star',iunit)
 call write_inopt(hacc,'hacc','accretion radius of the companion star',iunit)

 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ieos,polyk,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use io,           only:error,fatal
 use setunits,     only:read_options_and_set_units
 character(len=*), intent(in) :: filename
 integer,          intent(inout) :: ieos
 real,             intent(inout) :: polyk
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0

 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr)

 call read_inopt(ieos,'ieos',db,errcount=nerr) ! equation of state
 call read_inopt(a,'a',db,errcount=nerr)
 call read_inopt(mdon,'mdon',db,errcount=nerr)
 call read_inopt(macc,'macc',db,errcount=nerr)
 call read_inopt(hacc,'hacc',db,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
