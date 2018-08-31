!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup routine for uniform distribution
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    cs0       -- initial sound speed in code units
!    dist_unit -- distance unit (e.g. au)
!    ilattice  -- lattice type (1=cubic, 2=closepacked)
!    mass_unit -- mass unit (e.g. solarm)
!    nx        -- number of particles in x direction
!    rhozero   -- initial density in code units
!    xmax      -- xmax boundary
!    xmin      -- xmin boundary
!    ymax      -- ymax boundary
!    ymin      -- ymin boundary
!    zmax      -- zmax boundary
!    zmin      -- zmin boundary
!
!  DEPENDENCIES: boundary, dim, infile_utils, io, mpiutils, options, part,
!    physcon, prompting, set_dust, setup_params, unifdis, units
!+
!--------------------------------------------------------------------------
module setup
 use dim,          only:use_dust
 use options,      only:use_dustfrac
 use setup_params, only:rhozero
 implicit none
 public :: setpart

 integer :: npartx,ilattice,dust_method
 real    :: cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi
 character(len=20) :: dist_unit,mass_unit
 real(kind=8) :: udist,umass
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
 use dim,          only:maxvxyzu,h2chemistry
 use setup_params, only:npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use part,         only:abundance,iHI,dustfrac
 use physcon,      only:pi,mass_proton_cgs,kboltz,years,pc,solarm
 use set_dust,     only:set_dustfrac
 use units,        only:set_units
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
 npartx = 64
 ilattice = 1
 cs0 = 1.
 if (use_dust) then
    use_dustfrac = .true.
    dust_to_gas = 0.01
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
 call set_units(dist=udist,mass=umass,G=1.d0)
 call set_boundary(xmini,xmaxi,xmini,xmaxi,xmini,xmaxi)
 !
 ! setup particles
 !
 deltax = dxbound/npartx
 npart = 0
 npart_total = 0

 select case(ilattice)
 case(1)
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 case(2)
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 case default
    print*,' error: chosen lattice not available, using cubic'
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,hfact,npart,xyzh,nptot=npart_total)
 end select

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(1)
 if (id==master) print*,' initial sound speed = ',cs0,' pressure = ',cs0**2/gamma

 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (maxvxyzu >= 4) vxyzu(4,i) = cs0**2/(gamma*(gamma-1.))
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
 integer :: ierr

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
 npartx = 64
 if (id==master) then
    print*,' uniform setup... (max = ',nint((maxp)**(1/3.)),')'
    call prompt('enter number of particles in x direction ',npartx,1)
 endif
 call bcast_mpi(npartx)
 !
 ! mean density
 !
 rhozero = 1.
 if (id==master) call prompt(' enter density (gives particle mass)',rhozero,0.)
 call bcast_mpi(rhozero)
 !
 ! sound speed in code units
 !
 if (id==master) then
    cs0 = 1.
    call prompt(' enter sound speed in code units (sets polyk)',cs0,0.)
 endif
 call bcast_mpi(cs0)
 if (maxvxyzu < 4) then
    polyk = cs0**2
    print*,' polyk = ',polyk
 else
    polyk = 0.
 endif
 !
 ! dust to gas ratio
 !
 if (use_dustfrac) then
    call prompt('Enter dust to gas ratio',dust_to_gas,0.)
    call bcast_mpi(dust_to_gas)
 endif
 !
 ! type of lattice
 !
 if (id==master) then
    ilattice = 1
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
