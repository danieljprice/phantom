!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup of N sink particles in an hierarchical configuration
!
! :References: None
!
! :Owner: Simone Ceppi
!
! :Runtime parameters:
!   - a            : *semi-major axis (+ve) or period (-ve)*
!   - eccentricity : *eccentricity*
!   - hacc1        : *primary accretion radius*
!   - hacc2        : *secondary accretion radius*
!   - m1           : *mass of primary*
!   - m2           : *mass of secondary*
!
! :Dependencies: externalforces, infile_utils, io, options, part, physcon,
!   setbinary, units, sethierarchical
!
 implicit none
 public :: setpart

 !hierarchical system parameters
 character(len=100) :: hierarchy
 integer :: sink_num, hl_num
 character(len=10) :: sink_labels(10), hl_labels(10)
 real :: mass(10),accr(10),a(10),e(10),inc(10),O(10),w(10),f(10)

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
 use sethierarchical, only:set_hierarchical_default_options,set_hierarchical
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

 call set_hierarchical_default_options(hierarchy, sink_num, sink_labels, hl_labels, hl_num, &
                                              mass, accr, a, e, inc, O, w, f)

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Welcome to CHESS (Complete Hierarchical Endless System Setup)'
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


 call set_hierarchical(hierarchy,sink_num, sink_labels, hl_labels, hl_num, mass, accr, a, e, inc, O, w, f, &
                              nptmass, xyzmh_ptmass, vxyz_ptmass, ierr)

end subroutine setpart

subroutine write_setupfile(filename)
  use sethierarchical, only: write_hierarchical_setupfile
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')

 call write_hierarchical_setupfile(hierarchy, iunit, sink_num, sink_labels, hl_labels, hl_num, &
                                          mass, accr, a, e, inc, O, w, f)

 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
  use infile_utils, only:open_db_from_file,inopts,close_db!,read_inopt
    use sethierarchical, only: read_hierarchical_setupfile
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

 call read_hierarchical_setupfile(hierarchy, db, nerr, sink_num, sink_labels, hl_labels, hl_num, &
                                         mass, accr, a, e, inc, O, w, f)
    
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
