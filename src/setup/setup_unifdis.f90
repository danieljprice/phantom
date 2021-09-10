!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for uniform distribution
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Bzero       : *magnetic field strength in code units*
!   - cs0         : *initial sound speed in code units*
!   - dist_unit   : *distance unit (e.g. au)*
!   - dust_to_gas : *dust-to-gas ratio*
!   - ilattice    : *lattice type (1=cubic, 2=closepacked)*
!   - mass_unit   : *mass unit (e.g. solarm)*
!   - nx          : *number of particles in x direction*
!   - rhozero     : *initial density in code units*
!   - xmax        : *xmax boundary*
!   - xmin        : *xmin boundary*
!   - ymax        : *ymax boundary*
!   - ymin        : *ymin boundary*
!   - zmax        : *zmax boundary*
!   - zmin        : *zmin boundary*
!
! :Dependencies: boundary, cooling, dim, domain, eos, h2cooling,
!   infile_utils, io, mpiutils, options, part, physcon, prompting,
!   set_dust, setup_params, timestep, unifdis, units
!
 use dim,          only:use_dust,mhd
 use options,      only:use_dustfrac
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer           :: npartx,ilattice
 real              :: cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,Bzero
 character(len=20) :: dist_unit,mass_unit
 real(kind=8)      :: udist,umass

 !--change default defaults to reproduce the test from Section 5.6.7 of Price+(2018)
 logical :: BalsaraKim = .false.

 !--dust
 real    :: dust_to_gas

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu,h2chemistry,gr
 use setup_params, only:npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use part,         only:Bxyz,periodic,abundance,iHI,dustfrac,ndustsmall,ndusttypes,grainsize,graindens
 use physcon,      only:pi,mass_proton_cgs,kboltz,years,pc,solarm,micron
 use set_dust,     only:set_dustfrac
 use units,        only:set_units,unit_density
 use domain,       only:i_belong
 use eos,          only:gmw
 use options,      only:icooling,alpha,alphau
 use timestep,     only:dtmax,tmax,C_cour,C_force,C_cool,tolv
 use cooling,      only:Tfloor
 use h2cooling,    only:abundc,abundo,abundsi,abunde,dust_to_gas_ratio,iphoto
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma
 real,              intent(inout) :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 character(len=40) :: filename
 real    :: totmass,deltax
 integer :: i,ierr
 logical :: iexist
 !
 !--general parameters
 !
 time = 0.
 if (maxvxyzu < 4) then
    gamma = 1.
 else
    gamma = 5./3.
 endif
 !
 ! default units
 !
 mass_unit = 'solarm'
 dist_unit = 'pc'
 !
 ! set boundaries to default values
 !
 xmini = xmin; xmaxi = xmax
 ymini = ymin; ymaxi = ymax
 zmini = zmin; zmaxi = zmax
 !
 ! set default values for input parameters
 !
 npartx   = 64
 ilattice = 1
 rhozero  = 1.
 if (gr) then
    cs0 = 1.e-4
 else
    cs0 = 1.
 endif
 if (use_dust) then
    use_dustfrac = .true.
    dust_to_gas  = 0.01
    ndustsmall   = 1
    ndusttypes   = 1
    grainsize(1) = 1.*micron/udist
    graindens(1) = 3./unit_density
 endif
 if (BalsaraKim) then
    ! there is a typo in Price+ (2018) in stating the physical density;
    ! this mass_unit yields the correct value of 2.3e-24g/cm^3
    mass_unit = '33982786.25*solarm'
    dist_unit = 'kpc'
    xmini = -0.1; xmaxi = 0.1
    ymini = -0.1; ymaxi = 0.1
    zmini = -0.1; zmaxi = 0.1
    Bzero = 0.056117
    cs0   = sqrt(0.3*gamma)
    ilattice = 2
    filename=trim(fileprefix)//'.in'
    inquire(file=filename,exist=iexist)
    if (.not.iexist) then
       tmax    = 0.035688
       dtmax   = 1.250d-4
       C_cour  = 0.200
       C_force = 0.150
       C_cool  = 0.15
       tolv    = 1.0d3
       alpha   = 1
       alphau  = 0.1
       gmw     = 1.22
       Tfloor  = 3.
       if (h2chemistry) then
          ! flags controlling h2chemistry
          icooling = 1
          abundc  = 2.0d-4
          abundo  = 4.5d-4
          abundsi = 3.0d-5
          abunde  = 2.0d-4
          iphoto  = 0
          dust_to_gas_ratio = 0.010
       else
          icooling = 0
       endif
    endif
 endif
 !
 ! get disc setup parameters from file or interactive setup
 !
 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    !--read from setup file
    call read_setupfile(filename,ierr)
    if (id==master) call write_setupfile(filename)
    if (ierr /= 0) then
       stop
    endif
 elseif (id==master) then
    call setup_interactive(id,polyk)
    call write_setupfile(filename)
    stop 'rerun phantomsetup after editing .setup file'
 else
    stop
 endif
 !
 ! set units and boundaries
 !
 if (gr) then
    call set_units(mass=umass,c=1.d0,G=1.d0)
 else
    call set_units(dist=udist,mass=umass,G=1.d0)
 endif
 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 !
 ! setup particles
 !
 deltax = dxbound/npartx
 npart = 0
 npart_total = 0

 select case(ilattice)
 case(2)
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,&
                     npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 case default
    if (ilattice /= 1) print*,' error: chosen lattice not available, using cubic'
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,&
                     periodic,nptot=npart_total,mask=i_belong)
 end select

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(1)
 if (id==master) print*,' initial sound speed = ',cs0,' pressure = ',cs0**2/gamma

 if (maxvxyzu < 4 .or. gamma <= 1.) then
    polyk = cs0**2
 else
    polyk = 0.
 endif
 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (maxvxyzu >= 4 .and. gamma > 1.) vxyzu(4,i) = cs0**2/(gamma*(gamma-1.))
 enddo

 if (use_dustfrac) then
    do i=1,npart
       call set_dustfrac(dust_to_gas,dustfrac(:,i))
    enddo
 endif

 if (h2chemistry) then
    do i=1,npart
       abundance(:,i)   = 0.
       abundance(iHI,i) = 1.  ! assume all atomic hydrogen initially
    enddo
 endif

 if (mhd .and. balsarakim) then
    Bxyz = 0.
    do i = 1,npart
       Bxyz(1,i) = Bzero
    enddo
    ihavesetupB = .true.
 endif
end subroutine setpart

!------------------------------------------------------------------------
!
! interactive setup
!
!------------------------------------------------------------------------
subroutine setup_interactive(id,polyk)
 use io,        only:master
 use mpiutils,  only:bcast_mpi
 use dim,       only:maxp,maxvxyzu
 use prompting, only:prompt
 use units,     only:select_unit
 integer, intent(in)  :: id
 real,    intent(out) :: polyk
 integer              :: ierr

 if (id==master) then
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter mass unit (e.g. solarm,jupiterm,earthm)',mass_unit)
       call select_unit(mass_unit,umass,ierr)
       if (ierr /= 0) print "(a)",' ERROR: mass unit not recognised'
    enddo
    ierr = 1
    do while (ierr /= 0)
       call prompt('Enter distance unit (e.g. au,pc,kpc,0.1pc)',dist_unit)
       call select_unit(dist_unit,udist,ierr)
       if (ierr /= 0) print "(a)",' ERROR: length unit not recognised'
    enddo

    call prompt('enter xmin boundary',xmini)
    call prompt('enter xmax boundary',xmaxi,xmini)
    call prompt('enter ymin boundary',ymini)
    call prompt('enter ymax boundary',ymaxi,ymini)
    call prompt('enter zmin boundary',zmini)
    call prompt('enter zmax boundary',zmaxi,zmini)
 endif
 !
 ! number of particles
 !
 if (id==master) then
    print*,' uniform setup... (max = ',nint((maxp)**(1/3.)),')'
    call prompt('enter number of particles in x direction ',npartx,1)
 endif
 call bcast_mpi(npartx)
 !
 ! mean density
 !
 if (id==master) call prompt(' enter density (gives particle mass)',rhozero,0.)
 call bcast_mpi(rhozero)
 !
 ! sound speed in code units
 !
 if (id==master) then
    call prompt(' enter sound speed in code units (sets polyk)',cs0,0.)
 endif
 call bcast_mpi(cs0)
 !
 ! dust to gas ratio
 !
 if (use_dustfrac) then
    call prompt('Enter dust to gas ratio',dust_to_gas,0.)
    call bcast_mpi(dust_to_gas)
 endif
 !
 ! magnetic field strength
 if (mhd .and. balsarakim) then
    call prompt('Enter magnetic field strength in code units ',Bzero,0.)
    call bcast_mpi(Bzero)
 endif
 !
 ! type of lattice
 !
 if (id==master) then
    call prompt(' select lattice type (1=cubic, 2=closepacked)',ilattice,1)
 endif
 call bcast_mpi(ilattice)
end subroutine setup_interactive

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for uniform setup routine'

 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 !
 ! boundaries
 !
 write(iunit,"(/,a)") '# boundaries'
 call write_inopt(xmini,'xmin','xmin boundary',iunit)
 call write_inopt(xmaxi,'xmax','xmax boundary',iunit)
 call write_inopt(ymini,'ymin','ymin boundary',iunit)
 call write_inopt(ymaxi,'ymax','ymax boundary',iunit)
 call write_inopt(zmini,'zmin','zmin boundary',iunit)
 call write_inopt(zmaxi,'zmax','zmax boundary',iunit)
 !
 ! other parameters
 !
 write(iunit,"(/,a)") '# setup'
 call write_inopt(npartx,'nx','number of particles in x direction',iunit)
 call write_inopt(rhozero,'rhozero','initial density in code units',iunit)
 call write_inopt(cs0,'cs0','initial sound speed in code units',iunit)
 if (use_dustfrac) then
    call write_inopt(dust_to_gas,'dust_to_gas','dust-to-gas ratio',iunit)
 endif
 if (mhd .and. balsarakim) then
    call write_inopt(Bzero,'Bzero','magnetic field strength in code units',iunit)
 endif
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use units,        only:select_unit
 use io,           only:error
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 ! units
 !
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
 !
 ! boundaries
 !
 call read_inopt(xmini,'xmin',db,errcount=nerr)
 call read_inopt(xmaxi,'xmax',db,min=xmini,errcount=nerr)
 call read_inopt(ymini,'ymin',db,errcount=nerr)
 call read_inopt(ymaxi,'ymax',db,min=ymini,errcount=nerr)
 call read_inopt(zmini,'zmin',db,errcount=nerr)
 call read_inopt(zmaxi,'zmax',db,min=zmini,errcount=nerr)
 !
 ! other parameters
 !
 call read_inopt(npartx,'nx',db,min=8,errcount=nerr)
 call read_inopt(rhozero,'rhozero',db,min=0.,errcount=nerr)
 call read_inopt(cs0,'cs0',db,min=0.,errcount=nerr)
 if (use_dustfrac) then
    call read_inopt(dust_to_gas,'dust_to_gas',db,min=0.,errcount=nerr)
 endif
 if (mhd .and. balsarakim) then
    call read_inopt(Bzero,'Bzero',db,min=0.,errcount=nerr)
 endif
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 call close_db(db)
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_unifdis','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_unifdis','length unit not recognised')
    ierr = ierr + 1
 endif

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
