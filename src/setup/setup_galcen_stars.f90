!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for simulations of the Galactic Centre
!    Adapted by Daniel Price in collaboration with Jorge Cuadra
!
! :References: Paumard et al. (2006)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - datafile : *filename for star data (m,x,y,z,vx,vy,vz)*
!   - h_sink   : *sink particle radii in arcsec at 8kpc*
!   - m_gas    : *gas mass resolution in solar masses*
!
! :Dependencies: datafiles, dim, eos, infile_utils, io, part, physcon,
!   setup_params, spherical, timestep, units
!
 implicit none
 public :: setpart

 !
 ! setup options and default values for these
 !
 character(len=120) :: datafile = 'stars.m.pos.-vel.txt'
 real :: m_gas = 1.e-6 ! gas mass resolution in Msun
 real :: h_sink = 0.05 ! sink particle radii in arcsec at 8kpc

 private

contains

!----------------------------------------------------------------
!+
!  setup for galactic centre simulation (no gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,ihsoft,igas
 use units,     only:set_units,umass !,udist
 use physcon,   only:solarm,kpc,pi,au
 use io,        only:fatal,master
 use eos,       only:gmw
 use timestep,  only:dtmax
 use spherical, only:set_sphere
 use datafiles, only:find_phantom_datafile
 use setup_params, only:npart_total
 use infile_utils, only:get_options
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=len(datafile)) :: filename
 integer :: ierr,i
 real    :: scale,psep
!
! units (mass = mass of black hole, length = 1 arcsec at 8kpc)
!
 scale = (1./3600.)*(pi/180.)*8.*kpc
 call set_units(mass=3.5e6*solarm,dist=scale,G=1.d0)
!
! general parameters
!
 time = 0.
 hfact = 1.2
 polyk = 0.
 gamma = 5./3.
 gmw = 0.6  ! completely ionized, solar abu; eventually needs to be WR abu
 dtmax = 0.01
 !
 ! read setup parameters from the .setup file
 !
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
!
! space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 massoftype = m_gas*(solarm/umass)  ! mass resolution

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
!
! Set up the black hole at the Galactic centre
!
 nptmass = 1
 xyzmh_ptmass(4,1) = 1.  ! M=1 in code units by definition
 xyzmh_ptmass(ihacc,1)  = 0.1 ! accretion radius
 xyzmh_ptmass(ihsoft,1) = 0.1 ! no softening
!
! Read positions, masses and velocities of stars from file
!
 filename=find_phantom_datafile(datafile,'galcen')
 call read_ptmass_data(filename,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
 do i=2,nptmass
    xyzmh_ptmass(1:3,i)  = xyzmh_ptmass(1:3,i)
    xyzmh_ptmass(4,i)    = xyzmh_ptmass(4,i)
    xyzmh_ptmass(ihacc,i)  = h_sink
    xyzmh_ptmass(ihsoft,i) = h_sink
 enddo
!
! setup initial sphere of particles to prevent initialisation problems
!
 psep = 1.0
 npart_total = 0
 call set_sphere('cubic',id,master,0.,20.,psep,hfact,npart,xyzh,nptot=npart_total)
 vxyzu(4,:) = 5.317e-4
 npartoftype(igas) = npart

 if (nptmass == 0) call fatal('setup','no particles setup')
 if (ierr /= 0) call fatal('setup','ERROR during setup')

end subroutine setpart

!----------------------------------------------------------------
!+
!  read sink particle masses, positions and velocities from file
!+
!----------------------------------------------------------------
subroutine read_ptmass_data(filename,xyzmh_ptmass,vxyz_ptmass,n,ierr)
 use io, only:error
 character(len=*), intent(in) :: filename
 real,    intent(out)   :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: n
 integer, intent(out)   :: ierr
 integer :: iunit,n_input

 n_input = n
 open(newunit=iunit,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) then
    print "(/,2(a,/))",' ERROR opening "'//trim(filename)//'" for read of point mass data', &
                       ' -> this file should contain m,x,y,z,vx,vy,vz for each point mass, one per line'
 endif
 do while(ierr==0)
    n = n + 1
    if (n > size(xyzmh_ptmass(1,:))) then
       ierr = 66
    else
       read(iunit,*,iostat=ierr) xyzmh_ptmass(4,n),xyzmh_ptmass(1:3,n),vxyz_ptmass(1:3,n)
    endif
    if (ierr /= 0) n = n - 1
 enddo
 print "(a,i4,a)",' READ',n - n_input,' point masses from '//trim(filename)
 if (ierr==66) then
    call error('read_ptmass_data','array size exceeded in read_ptmass_data, recompile with MAXPTMASS=n',var='n',ival=n+1)
 endif

 ! end of file error is OK
 if (ierr < 0) ierr = 0

end subroutine read_ptmass_data

!------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use dim,          only:tagline
 character(len=*), intent(in) :: filename
 integer                      :: lu,ierr1,ierr2

 write(*,"(a)") ' Writing '//trim(filename)//' with setup options'
 open(newunit=lu,file=filename,status='replace',form='formatted')
 write(lu,"(a)") '# '//trim(tagline)
 write(lu,"(a)") '# input file for Phantom galactic centre setup'

 write(lu,"(/,a)") '# datafile'
 call write_inopt(datafile,'datafile','filename for star data (m,x,y,z,vx,vy,vz)',lu,ierr1)

 write(lu,"(/,a)") '# resolution'
 call write_inopt(m_gas, 'm_gas','gas mass resolution in solar masses',lu,ierr2)
 call write_inopt(h_sink, 'h_sink','sink particle radii in arcsec at 8kpc',lu,ierr2)
 close(lu)

end subroutine write_setupfile

!------------------------------------------
!+
!  Read setup parameters from input file
!+
!------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 use dim,          only:maxvxyzu
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*,'(1x,2a)') 'Setup_galcen: Reading setup options from ',trim(filename)

 nerr = 0
 call read_inopt(datafile,'datafile',db,errcount=nerr)
 call read_inopt(m_gas,'m_gas',db,errcount=nerr)
 call read_inopt(h_sink,'h_sink',db,errcount=nerr)

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_galcen: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif
 call close_db(db)

end subroutine read_setupfile

end module setup
