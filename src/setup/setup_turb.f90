!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Sets up a calculation of supersonic turbulence in a periodic box.
!  Works for hydro, mhd, and dusty turbulence.
!
! :References:
!    Price & Federrath (2010), MNRAS
!    Tricco, Price & Federrath (2016), MNRAS
!    Tricco, Price & Laibe (2017), MNRAS Letters
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Bz_0     : *initial magnetic field strength*
!   - ilattice : *lattice type (1=cubic, 2=closepacked)*
!   - npartx   : *number of particles in x direction*
!   - polykset : *sound speed in code units (sets polyk)*
!   - rhozero  : *density (gives particle mass)*
!
! :Dependencies: boundary, dim, dust, infile_utils, io, kernel, mpidomain,
!   options, part, physcon, prompting, set_dust, set_dust_options,
!   setup_params, timestep, unifdis, units
!
 use dim,          only:mhd,use_dust,maxdustsmall
 use dust,         only:grainsizecgs,graindenscgs,ilimitdustflux
 use part,         only:ndustsmall
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: polykset
 real,    private :: Bz_0
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:use_dust,maxdustsmall,maxvxyzu,periodic,curlv
 use options,      only:use_dustfrac,nfulldump
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis,latticetype
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd,dustfrac,grainsize,graindens,ndusttypes,ndustsmall,igas
 use physcon,      only:pi,solarm,pc,km
 use units,        only:set_units,udist,umass
 use set_dust,     only:set_dustfrac
 use timestep,     only:dtmax,tmax
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options,infile_exists
 use kernel,       only:hfact_default
 use set_dust_options, only:dustbinfrac,set_dust_grain_distribution,ilimitdustfluxinp,&
                            ndustsmallinp,set_dust_default_options,dust_to_gas
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 integer :: i,ierr
 real :: totmass,deltax
 real :: Bz_0

 !
 !--Default runtime parameters
 !
 npartx = 64
 ilattice = 1
 polykset = 1.
 rhozero = 1.
 call set_dust_default_options()
 Bz_0 = 1.4142e-5

 !
 !--Read runtime parameters from setup file
 !
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Turbulence setup'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 print "(/,a,/)", ' Setup for turbulence in a periodic box'
 !
 !--boundaries
 !
 call set_boundary(0.,1.,0.,1.,0.,1.)
 !
 !--units
 !
 ! Molecular cloud conditions:
 !  L = 3 pc, rho = 1e-20 g/cm^3, c_s = 0.2 km/s
 call set_units(dist=3.0*pc, mass=1e-20*(3.0*pc)**3, time=3.0*pc/(0.2*km))
 !
 !--general parameters
 !
 time = 0.
 hfact = hfact_default
 gamma = 1.
 !
 !--setup particles
 !
 deltax = dxbound/npartx

 if (maxvxyzu < 4) then
    polyk = polykset**2
    print*,' polyk = ',polyk
 else
    polyk = 0.
    polykset = 0.
 endif

 if (use_dust) then
    print *, ''
    print *, 'Setting up dusty turbulence with dust-as-mixture with ',ndustsmallinp,' small grain species'
    print *, ''

    !--currently only works with dustfrac
    use_dustfrac = .true.

    ! set dust grid
    ndustsmall = ndustsmallinp
    call set_dust_grain_distribution(ndusttypes,dustbinfrac,grainsize,graindens,udist,umass)
    ilimitdustflux = ilimitdustfluxinp
 endif

 if (mhd) print "(/,a,/)",' MHD turbulence w/uniform field in z-direction'

 ! setup preferred values of .in file
 if (.not. infile_exists(fileprefix)) then
    tmax         = 1.00   ! run for 20 turbulent crossing times
    dtmax        = 0.0025
    nfulldump    = 5      ! output 4 full dumps per crossing time
    curlv        = .true.
 endif
 npart = 0
 npart_total = 0

 call set_unifdis(latticetype(ilattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print "(a,i0,1x,i0)", ' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 print "(a,1pg12.5)", ' particle mass = ',massoftype(igas)

 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(3,i) = Bz_0
    endif
    !--one fluid dust: set dust fraction on gas particles
    if (use_dust .and. use_dustfrac) then
       if (ndusttypes > 1) then
          dustfrac(1:ndusttypes,i) = dust_to_gas*dustbinfrac(1:ndusttypes)
       else
          call set_dustfrac(dust_to_gas,dustfrac(:,i))
       endif
    endif
 enddo

 if (mhd) ihavesetupB = .true.

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Write setup parameters to .setup file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename)
 use set_dust_options, only:write_dust_setup_options
 use infile_utils,     only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for turbulence setup routine'
 call write_inopt(npartx,'npartx','number of particles in x direction',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 call write_inopt(rhozero,'rhozero','density (gives particle mass)',iunit)
 call write_inopt(polykset,'polykset','sound speed in code units (sets polyk)',iunit)
 if (use_dust) call write_dust_setup_options(iunit,method=1)

 if (mhd) call write_inopt(Bz_0,'Bz_0','initial magnetic field strength',iunit)
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,     only:open_db_from_file,inopts,read_inopt,close_db
 use set_dust_options, only:read_dust_setup_options
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartx,'npartx',db,min=16,errcount=nerr)
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(polykset,'polykset',db,min=0.,errcount=nerr)
 if (use_dust) call read_dust_setup_options(db,nerr,method=1)
 if (mhd) call read_inopt(Bz_0,'Bz_0',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

!-----------------------------------------------------------------------
!+
!  Interactive setup
!+
!-----------------------------------------------------------------------
subroutine setup_interactive()
 use prompting,        only:prompt
 use set_dust_options, only:set_dust_interactive

 call prompt('Enter number of particles in x ',npartx,16)
 call prompt('Select lattice type (1=cubic, 2=closepacked)',ilattice,1,2)
 call prompt('Enter density (gives particle mass)',rhozero,0.)
 call prompt('Enter sound speed in code units (sets polyk)',polykset,0.)
 if (use_dust) call set_dust_interactive(method=1)
 if (mhd) call prompt('Enter initial magnetic field strength ',Bz_0)

end subroutine setup_interactive

end module setup
