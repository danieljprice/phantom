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
!   - Bz_0           : *initial magnetic field strength*
!   - dust_to_gas    : *total dust to gas ratio*
!   - graindenscgs   : *grain density in g/cm^3*
!   - grainsizecgs   : *grain size in cm*
!   - ilattice       : *lattice type (1=cubic, 2=closepacked)*
!   - ilimitdustflux : *limit the dust flux*
!   - ndustsmall     : *number of grain sizes*
!   - npartx         : *number of particles in x direction*
!   - polykset       : *sound speed in code units (sets polyk)*
!   - rhozero        : *density (gives particle mass)*
!   - sindex         : *power-law index, e.g. MRN*
!   - smaxcgs        : *maximum grain size in cm*
!   - smincgs        : *minimum grain size in cm*
!
! :Dependencies: boundary, dim, dust, infile_utils, io, kernel, mpidomain,
!   options, part, physcon, prompting, set_dust, setup_params, table_utils,
!   timestep, unifdis, units
!
 use dim,          only:mhd,use_dust,maxdustsmall
 use dust,         only:grainsizecgs,graindenscgs,ilimitdustflux
 use part,         only:ndustsmall
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: polykset
 real,    private :: dust_to_gas
 real,    private :: smincgs,smaxcgs,sindex
 real,    private :: Bz_0
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:use_dust,maxdustsmall,maxvxyzu,periodic
 use options,      only:use_dustfrac,nfulldump,beta
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis,latticetype
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use part,         only:Bxyz,mhd,dustfrac,grainsize,graindens,ndusttypes,ndustsmall,igas
 use physcon,      only:pi,solarm,pc,km
 use units,        only:set_units,udist,umass
 use prompting,    only:prompt
 use set_dust,     only:set_dustfrac,set_dustbinfrac
 use timestep,     only:dtmax,tmax
 use table_utils,  only:logspace
 use mpidomain,    only:i_belong
 use infile_utils, only:get_options,infile_exists
 use kernel,       only:hfact_default
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
 real :: dustbinfrac(maxdustsmall)

 !
 !--Default runtime parameters
 !
 npartx = 64
 ilattice = 1
 polykset = 1.
 rhozero = 1.
 dust_to_gas = 1.e-2
 ndustsmall = 1
 ilimitdustflux = .false.
 smincgs = 1.e-5
 smaxcgs = 1.
 sindex = 3.5
 grainsizecgs = 1.e-4
 graindenscgs = 3.0
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
    print *, 'Setting up dusty turbulence:'
    print *, ''

    !--currently only works with dustfrac
    use_dustfrac = .true.

    ndusttypes = ndustsmall
    grainsize = 0.
    graindens = 0.

    if (ndusttypes > 1) then
       ! grainsize/mass distribution
       call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
       grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
       ! grain density
       graindens(1:ndusttypes) = graindenscgs/umass*udist**3
    else
       grainsize(1) = grainsizecgs/udist
       graindens(1) = graindenscgs/umass*udist**3
    endif

 endif

 if (mhd) print "(/,a,/)",' MHD turbulence w/uniform field in z-direction'

 ! setup preferred values of .in file
 if (.not. infile_exists(fileprefix)) then
    tmax         = 1.00   ! run for 20 turbulent crossing times
    dtmax        = 0.0025
    nfulldump    = 5      ! output 4 full dumps per crossing time
    beta         = 4      ! legacy from Price & Federrath (2010), haven't checked recently if still required
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
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for turbulence setup routine'
 call write_inopt(npartx,'npartx','number of particles in x direction',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 call write_inopt(rhozero,'rhozero','density (gives particle mass)',iunit)
 call write_inopt(polykset,'polykset','sound speed in code units (sets polyk)',iunit)
 if (use_dust) then
    call write_inopt(dust_to_gas,'dust_to_gas','total dust to gas ratio',iunit)
    call write_inopt(ndustsmall,'ndustsmall','number of grain sizes',iunit)
    call write_inopt(ilimitdustflux,'ilimitdustflux','limit the dust flux',iunit)
    call write_inopt(smincgs,'smincgs','minimum grain size in cm',iunit)
    call write_inopt(smaxcgs,'smaxcgs','maximum grain size in cm',iunit)
    call write_inopt(sindex,'sindex','power-law index, e.g. MRN',iunit)
    call write_inopt(grainsizecgs,'grainsizecgs','grain size in cm',iunit)
    call write_inopt(graindenscgs,'graindenscgs','grain density in g/cm^3',iunit)
 endif
 if (mhd) call write_inopt(Bz_0,'Bz_0','initial magnetic field strength',iunit)
 close(iunit)

end subroutine write_setupfile

!-----------------------------------------------------------------------
!+
!  Read setup parameters from .setup file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
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
 if (use_dust) then
    call read_inopt(dust_to_gas,'dust_to_gas',db,min=0.,errcount=nerr)
    call read_inopt(ndustsmall,'ndustsmall',db,min=1,errcount=nerr)
    call read_inopt(ilimitdustflux,'ilimitdustflux',db,errcount=nerr)
    call read_inopt(smincgs,'smincgs',db,min=0.,errcount=nerr)
    call read_inopt(smaxcgs,'smaxcgs',db,min=0.,errcount=nerr)
    call read_inopt(sindex,'sindex',db,errcount=nerr)
    call read_inopt(grainsizecgs,'grainsizecgs',db,min=0.,errcount=nerr)
    call read_inopt(graindenscgs,'graindenscgs',db,min=0.,errcount=nerr)
 endif
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
 use prompting, only:prompt

 call prompt('Enter number of particles in x ',npartx,16)
 call prompt('Select lattice type (1=cubic, 2=closepacked)',ilattice,1,2)
 call prompt('Enter density (gives particle mass)',rhozero,0.)
 call prompt('Enter sound speed in code units (sets polyk)',polykset,0.)
 if (use_dust) then
    call prompt('Enter total dust to gas ratio',dust_to_gas,0.)
    call prompt('How many grain sizes do you want?',ndustsmall,1,maxdustsmall)
    call prompt('Do you want to limit the dust flux?',ilimitdustflux)
    if (ndustsmall > 1) then
       call prompt('Enter minimum grain size in cm',smincgs,0.)
       call prompt('Enter maximum grain size in cm',smaxcgs,0.)
       call prompt('Enter power-law index, e.g. MRN',sindex)
    else
       call prompt('Enter grain size in cm',grainsizecgs,0.)
    endif
    call prompt('Enter grain density in g/cm^3',graindenscgs,0.)
 endif
 if (mhd) call prompt('Enter initial magnetic field strength ',Bz_0)

end subroutine setup_interactive

end module setup
