!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for uniform distribution
!
! :References: None
!
! :Owner: Sergei Biriukov
!
! :Runtime parameters:
!   - cs0      : *initial sound speed in code units*
!   - ilattice : *lattice type (1=cubic, 2=closepacked)*
!   - nx       : *number of particles in x direction*
!   - rhozero  : *initial density in code units*
!   - xmax     : *xmax boundary*
!   - xmin     : *xmin boundary*
!   - ymax     : *ymax boundary*
!   - ymin     : *ymin boundary*
!   - zmax     : *zmax boundary*
!   - zmin     : *zmin boundary*
!
! :Dependencies: boundary, dim, eos, infile_utils, io, kernel, mpidomain,
!   mpiutils, options, part, physcon, set_dust, setunits, setup_params,
!   timestep, unifdis, units
!
 use setup_params, only:rhozero
 use dim,          only:gr
 implicit none
 public :: setpart

 integer :: npartx,ilattice,iradtype
 real    :: cs0,xmini,xmaxi,ymini,ymaxi,zmini,zmaxi

 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params,  only:npart_total
 use io,            only:master,fatal,warning
 use unifdis,       only:set_unifdis
 use boundary,      only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound,set_boundary
 use physcon,       only:pi,mass_proton_cgs,kboltz,years,pc,solarm,c,Rg,steboltz
 use set_dust,      only:set_dustfrac
 use units,         only:set_units,unit_ergg,unit_velocity,unit_opacity,get_c_code,get_radconst_code
 use part,          only:rhoh,igas,rad,radprop,ithick,iradxi,ikappa,periodic
 use eos,           only:gmw
 use kernel,        only:hfact_default
 use timestep,      only:dtmax,tmax,C_rad
 use options,       only:nfulldump
 use mpidomain,     only:i_belong

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

 real :: a,c_code,cv1,kappa_code,pmassi,Tref,xi0
 real :: rhoi,h0,rho0

 call setup_setdefaults(id,polyk,gamma,xmin,xmax,ymin,ymax,zmin,zmax,npartx,cs0)

 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    !--read from setup file
    call read_setupfile(filename,gamma,ierr)
    if (id==master) call write_setupfile(filename,gamma)
    if (ierr /= 0) then
       stop
    endif
 elseif (id==master) then
    call write_setupfile(filename,gamma)
    stop 'rerun phantomsetup after editing .setup file'
 else
    stop
 endif
 !
 ! set boundaries
 !
 call set_boundary(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi)
 !
 ! setup particles
 !
 deltax = dxbound/npartx
 npart = 0
 npart_total = 0
 hfact = hfact_default

 select case(ilattice)
 case(1)
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                     hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 case(2)
    call set_unifdis('closepacked',id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                     deltax,hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 case default
    print*,' error: chosen lattice not available, using cubic'
    call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,deltax,&
                     hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
 end select

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print*,' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 if (id==master) print*,' particle mass = ',massoftype(igas)
 if (id==master) print*,' initial sound speed = ',cs0,' pressure = ',cs0**2/gamma

 vxyzu(1:3,:) = 0.

 c_code = get_c_code()
 cv1 = (gamma-1.)*gmw/Rg*unit_velocity**2
 a   = get_radconst_code()
 pmassi = massoftype(igas)

 Tref = 100

 ! xi0 = a*Tref**4.0/rho0
 !
 select case(iradtype)
 case(1)

 case(2)
    h0   = xyzh(4,1)
    rho0 = rhoh(h0,pmassi)
    kappa_code = 1.0/unit_opacity

    dtmax = C_rad*h0*h0*rho0*kappa_code/c_code*25
    tmax  = 10*dtmax
    nfulldump = 1

    radprop(ithick,:) = 1.
    do i=1,npart
       rhoi = rhoh(xyzh(4,1),pmassi)
       vxyzu(4,i) = (Tref/cv1)/(unit_ergg)
       xi0 = a*Tref**4.0/rhoi
       radprop(ikappa,i) = kappa_code
       rad(iradxi,i) = xi0*(1 + 1e-1*sin(xyzh(1,i)*2*pi/(xmax-xmin)))
    enddo
 case(3)
    h0   = xyzh(4,1)
    rho0 = rhoh(h0,pmassi)
    kappa_code = 1.0/unit_opacity

    dtmax = C_rad*h0*h0*rho0*kappa_code/c_code/5
    tmax  = 200*dtmax
    nfulldump = 1

    radprop(ithick,:) = 1.
    do i=1,npart
       rhoi = rhoh(xyzh(4,1),pmassi)
       vxyzu(4,i) = (Tref/cv1)/(unit_ergg)
       xi0 = a*Tref**4.0/rhoi
       radprop(ikappa,i) = kappa_code
       rad(iradxi,i) = xi0*exp(-500.*xyzh(1,i)**2)
    enddo
 case default
    call fatal('setup_radiativebox', 'radiation setup is not available')
 end select

end subroutine setpart

!------------------------------------------------------------------------
!
! interactive setup
!
!------------------------------------------------------------------------
subroutine setup_setdefaults(id,polyk,gamma,xmin,xmax,ymin,ymax,zmin,zmax,&
                             npartx,cs0)
 use io,        only:master
 use mpiutils,  only:bcast_mpi
 use options,   only:exchange_radiation_energy,limit_radiation_flux
 use setunits,  only:mass_unit,dist_unit
 integer, intent(in)  :: id
 integer, intent(out) :: npartx
 real,    intent(in)  :: xmin,xmax,ymin,ymax,zmin,zmax
 real,    intent(out) :: polyk,gamma,cs0

 mass_unit = 'solarm'
 dist_unit = 'au'

 xmini = xmin
 xmaxi = xmax
 ymini = ymin
 ymaxi = ymax
 zmini = zmin
 zmaxi = zmax

 npartx = 32
 rhozero = 1.
 gamma = 5./3.
 cs0 = 1.
 polyk = 0.
 ilattice = 2
 iradtype = 3
 exchange_radiation_energy = .false.
 limit_radiation_flux = .false.

 if (id==master) then
    call bcast_mpi(npartx)
    call bcast_mpi(rhozero)
    call bcast_mpi(cs0)
    call bcast_mpi(ilattice)
    call bcast_mpi(iradtype)
 endif

end subroutine setup_setdefaults

!------------------------------------------------------------------------
!
! write setup file
!
!------------------------------------------------------------------------
subroutine write_setupfile(filename,gamma)
 use infile_utils, only:write_inopt
 use options,      only:exchange_radiation_energy,limit_radiation_flux
 use setunits,     only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit
 real, intent(in) :: gamma

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for uniform setup routine'

 call write_options_units(iunit)
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
 call write_inopt(gamma,'gamma','',iunit)
 call write_inopt(ilattice,'ilattice','lattice type (1=cubic, 2=closepacked)',iunit)
 call write_inopt(iradtype,'iradtype',&
 'type of radiation setup (1=uniform,2=sin diffusion,3=gaussian faster-than-light diffusion)',iunit)
 call write_inopt(exchange_radiation_energy,'gas-rad_exchange',&
    'do or do not exchange energy  between gas and radiation',iunit)
 call write_inopt(limit_radiation_flux,'flux_limiter',&
    'do or do not limit radiation  flux',iunit)
 close(iunit)

end subroutine write_setupfile

!------------------------------------------------------------------------
!
! read setup file
!
!------------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use setunits,     only:read_options_and_set_units
 use io,           only:error
 use options,      only:exchange_radiation_energy,limit_radiation_flux
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)
 real, intent(out) :: gamma

 print "(a)",' reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 !
 ! units
 !
 call read_options_and_set_units(db,nerr,gr)
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
 call read_inopt(gamma,'gamma',db,min=0.,errcount=nerr)
 call read_inopt(ilattice,'ilattice',db,min=1,max=2,errcount=nerr)
 call read_inopt(iradtype,'iradtype',db,min=1,max=3,errcount=nerr)
 call read_inopt(exchange_radiation_energy,'gas-rad_exchange',db,errcount=nerr)
 call read_inopt(limit_radiation_flux,'flux_limiter',db,errcount=nerr)

 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
