!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of N sink particles in an hierarchical configuration
!
! :References: None
!
! :Owner: Simone Ceppi
!
! :Runtime parameters: None
!
! :Dependencies: infile_utils, io, kernel, part, physcon, sethierarchical,
!   units
!

 implicit none
 public :: setpart

 private

contains

!----------------------------------------------------------------
!+
!  setup for sink particle binary simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,gr
 use sethierarchical, only:set_hierarchical_default_options,set_hierarchical,print_chess_logo
 use units,           only:set_units
 use physcon,         only:solarm,au,pi
 use infile_utils,    only:get_options
 use io,              only:master
 use kernel,          only:hfact_default
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer :: ierr
!
!--units
!
 if (gr) then
    call set_units(mass=solarm,G=1.d0,c=1.d0)
 else
    call set_units(mass=solarm,dist=au,G=1.d0)
 endif
!
!--general parameters
!
 time = 0.
 polyk = 0.
 gamma = 1.
 hfact = hfact_default
!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 massoftype = 1d-9

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0

 call set_hierarchical_default_options()
 if (id==master) call print_chess_logo()

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                 read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 call set_hierarchical(fileprefix, nptmass, xyzmh_ptmass, vxyz_ptmass, ierr)

end subroutine setpart

subroutine write_setupfile(filename)
 use sethierarchical, only:write_hierarchical_setupfile
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_hierarchical_setupfile(iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db!,read_inopt
 use sethierarchical, only:read_hierarchical_setupfile
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
 call read_hierarchical_setupfile(db, nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
