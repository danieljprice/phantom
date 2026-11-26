!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup of two stars or sink particles in a binary
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - corotate   : *set stars in corotation*
!   - deltat     : *output interval as fraction of binary period*
!   - gwinspiral : *set gravitational wave inspiral*
!   - norbits    : *maximum number of binary orbits*
!
! :Dependencies: centreofmass, dim, eos, extern_gwinspiral, externalforces,
!   infile_utils, io, kernel, mpidomain, options, part, physcon,
!   sethier_utils, sethierarchical, setorbit, setstar, setunits,
!   setup_params, timestep, units
!
 use setstar,       only:star_t
 use setorbit,      only:orbit_t
 use sethier_utils, only:max_hier_levels
 use dim,      only:gr
 use eos,      only:ieos
 implicit none
 public :: setpart

 logical :: relax,write_rho_to_file,corotate,gwinspiral
 integer, parameter :: max_stars = max_hier_levels
 type(star_t)  :: star(max_stars)
 type(orbit_t) :: orbit
 real :: norbits,deltat
 integer :: nstar

 private

contains

!----------------------------------------------------------------
!+
!  setup for binary star simulations (with or without gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,ihacc,eos_vars,rad,nsinkproperties,iJ2,&
                          iReff,ispinx,ispinz
 use setorbit,       only:set_defaults_orbits,set_orbits
 use options,        only:iexternalforce,alphau
 use externalforces, only:iext_corotate,iext_geopot,iext_star,iext_gwinspiral,omega_corotate,mass1,accradius1
 use extern_gwinspiral, only:Nstar_gw
 use io,             only:master,fatal
 use setstar,        only:set_defaults_stars,set_stars,shift_stars
 use eos,            only:X_in,Z_in,use_var_comp
 use setup_params,   only:rhozero,npart_total
 use mpidomain,      only:i_belong
 use centreofmass,   only:reset_centreofmass
 use setunits,       only:mass_unit,dist_unit
 use physcon,        only:deg_to_rad
 use kernel,         only:hfact_default
 use infile_utils,   only:get_options,infile_exists
 use timestep,       only:tmax,dtmax
 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma,hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 real,             intent(out)   :: vxyzu(:,:)
 integer :: ierr,nptmass_in,iextern_prev,i
 logical :: add_spin
 real :: xyzmh_ptmass_in(nsinkproperties,max_stars),vxyz_ptmass_in(3,max_stars),angle
 real :: m(max_stars),hacc(max_stars)
 logical, parameter :: set_oblateness = .false.
!
!--general parameters
!
 dist_unit = 'solarr'
 mass_unit = 'solarm'
 time = 0.
 polyk = 0.
 gamma = 5./3.
 ieos  = 2
 hfact = hfact_default
 norbits = 10.
 deltat = 0.1
!
!--space available for injected gas particles
!  in case only sink particles are used
!
 npart = 0
 npartoftype(:) = 0
 massoftype = 0.

 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 nptmass = 0
 nstar = 2
 call set_defaults_stars(star)
 call set_defaults_orbits(nstar,orbit)
 relax = .true.
 corotate = .false.
 gwinspiral = .false.
 use_var_comp = .false.
 write_rho_to_file = .true.

 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",&
   ' Welcome to the Ultimate Binary & Multiple Setup'

 !
 ! read/write from .setup file
 !
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'

 !
 !--setup and relax stars as needed
 !
 iextern_prev = iexternalforce
 iexternalforce = 0

 call set_stars(id,master,nstar,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,&
                X_in,Z_in,relax,use_var_comp,write_rho_to_file,&
                rhozero,npart_total,i_belong,ierr)

 add_spin = .false.
 xyzmh_ptmass_in(:,:) = 0.
 vxyz_ptmass_in(:,:) = 0.
 nptmass_in = 0
 if (nstar >= 2) then
    ! place stars on orbits: a binary for nstar=2, a hierarchical system for nstar>=3,
    ! with star i on sink i
    do i=1,nstar
       m(i)    = star(i)%m_code
       hacc(i) = star(i)%hacc_code
    enddo
    if (nstar == 2 .and. iextern_prev==iext_corotate) then
       call set_orbits(nstar,m,hacc,orbit,fileprefix,xyzmh_ptmass_in,vxyz_ptmass_in,&
                       nptmass_in,(id==master),ierr,omega_corotate)
    else
       call set_orbits(nstar,m,hacc,orbit,fileprefix,xyzmh_ptmass_in,vxyz_ptmass_in,&
                       nptmass_in,(id==master),ierr)
       add_spin = corotate .and. nstar == 2
    endif
    if (ierr /= 0) call fatal ('setup_binary','error in call to set_orbits')
 endif
 !
 !--place stars into orbit, or add real sink particles if iprofile=0
 !
 call shift_stars(nstar,star,xyzmh_ptmass_in,vxyz_ptmass_in,xyzh,vxyzu,&
                  xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,nptmass,corotate=add_spin)
 !
 !--restore options
 !
 iexternalforce = iextern_prev

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (nstar >= 2 .and. (iexternalforce==iext_geopot .or. iexternalforce==iext_star) .and. star(1)%iprofile == 0) then
    ! delete first sink particle and copy its properties to the central potential
    nptmass = nptmass - 1
    mass1 = xyzmh_ptmass(4,nptmass+1)
    accradius1 = xyzmh_ptmass(ihacc,nptmass+1)
    xyzmh_ptmass(:,nptmass) = xyzmh_ptmass(:,nptmass+1)
    vxyz_ptmass(:,nptmass) = vxyz_ptmass(:,nptmass+1)
 elseif (nstar == 2 .and. (gwinspiral .or. iexternalforce==iext_gwinspiral)) then
    Nstar_gw(1) = star(1)%np
    Nstar_gw(2) = npart - star(1)%np
    iexternalforce = iext_gwinspiral
 elseif (nstar >= 1 .and. set_oblateness .and. star(1)%iprofile == 0) then
    ! set J2 for sink particle 1 to be equal to oblateness of Saturn
    xyzmh_ptmass(iJ2,1) = 0.01629
    angle = 30.*deg_to_rad
    xyzmh_ptmass(ispinx,1) = sin(angle)
    xyzmh_ptmass(ispinz,1) = cos(angle)
    xyzmh_ptmass(iReff,1) = xyzmh_ptmass(ihacc,1)
 endif

 alphau = min(0.1,alphau)
 if (.not.infile_exists(fileprefix) .and. nstar == 2) then
    tmax = norbits * orbit%period
    dtmax = deltat * orbit%period
 endif

end subroutine setpart

!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,    only:write_inopt
 use setstar,         only:write_options_stars
 use setorbit,        only:write_options_orbits
 use setunits,        only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary and multiple star setup routines'

 call write_options_units(iunit,gr)
 call write_options_stars(star,relax,write_rho_to_file,ieos,iunit,nstar)
 if (nstar == 2) then
    call write_inopt(corotate,'corotate','set stars in corotation',iunit)
    call write_inopt(gwinspiral,'gwinspiral','set gravitational wave inspiral',iunit)
 endif
 call write_options_orbits(nstar,orbit,iunit)
 if (nstar == 2) then
    write(iunit,"(/,a)") '# timestepping'
    call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
    call write_inopt(deltat,'deltat','output interval as fraction of binary period',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  read options from .setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils,    only:open_db_from_file,inopts,read_inopt,close_db
 use setstar,         only:read_options_stars
 use setorbit,        only:read_options_orbits
 use setunits,        only:read_options_and_set_units
 use units,           only:in_code_units
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr,ierr1
 type(inopts), allocatable :: db(:)
 real :: m(2)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr,gr)
 call read_options_stars(star,ieos,relax,write_rho_to_file,db,nerr,nstar)
 m = 0.
 if (nstar == 2) then
    call read_inopt(corotate,'corotate',db,errcount=nerr)
    call read_inopt(gwinspiral,'gwinspiral',db,errcount=nerr,default=.false.)
    m(1) = in_code_units(star(1)%m,ierr,unit_type='mass')
    m(2) = in_code_units(star(2)%m,ierr1,unit_type='mass')
    if (ierr /= 0 .or. ierr1 /= 0) then
       print*,' ERROR: error converting masses to code units'
       nerr = nerr + 1
    endif
 endif
 call read_options_orbits(nstar,orbit,m,db,nerr)
 if (nstar == 2) then
    call read_inopt(norbits,'norbits',db,errcount=nerr)
    call read_inopt(deltat,'deltat',db,errcount=nerr,default=0.1)
 endif
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
