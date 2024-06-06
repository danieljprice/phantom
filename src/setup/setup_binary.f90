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
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - corotate : *set stars in corotation*
!
! :Dependencies: centreofmass, dim, eos, externalforces, infile_utils, io,
!   kernel, mpidomain, options, part, physcon, setorbit, setstar, setunits,
!   setup_params
!
 use setstar,  only:star_t
 use setorbit, only:orbit_t
 use dim,     only:gr
 implicit none
 public :: setpart

 logical :: relax,corotate
 type(star_t)  :: star(2)
 type(orbit_t) :: orbit

 private

contains

!----------------------------------------------------------------
!+
!  setup for binary star simulations (with or without gas)
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,&
                   polyk,gamma,hfact,time,fileprefix)
 use part,           only:gr,nptmass,xyzmh_ptmass,vxyz_ptmass,&
                          ihacc,ihsoft,eos_vars,rad,nsinkproperties,iJ2,iReff,ispinx,ispinz
 use setorbit,       only:set_defaults_orbit,set_orbit
 use options,        only:iexternalforce
 use externalforces, only:iext_corotate,iext_geopot,iext_star,omega_corotate,mass1,accradius1
 use io,             only:master,fatal
 use setstar,        only:set_defaults_stars,set_stars,shift_stars
 use eos,            only:X_in,Z_in,ieos
 use setup_params,   only:rhozero,npart_total
 use mpidomain,      only:i_belong
 use centreofmass,   only:reset_centreofmass
 use setunits,       only:mass_unit,dist_unit
 use physcon,        only:deg_to_rad
 use kernel,         only:hfact_default
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
 integer :: ierr,nstar,nptmass_in,iextern_prev
 logical :: iexist,write_profile,use_var_comp,add_spin
 real :: xyzmh_ptmass_in(nsinkproperties,2),vxyz_ptmass_in(3,2),angle
 logical, parameter :: set_oblateness = .false.
!
!--general parameters
!
 dist_unit = 'solarr'
 mass_unit = 'solarm'
 time = 0.
 polyk = 0.
 gamma = 1.
 hfact = hfact_default
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
 call set_defaults_orbit(orbit)
 relax = .true.
 corotate = .false.
 ieos = 2

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
 !--setup and relax stars as needed
 !
 use_var_comp = .false.
 write_profile = .false.
 iextern_prev = iexternalforce
 iexternalforce = 0
 gamma = 5./3.
 call set_stars(id,master,nstar,star,xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                X_in,Z_in,relax,use_var_comp,write_profile,&
                rhozero,npart_total,i_belong,ierr)

 nptmass_in = 0
 if (iexternalforce==iext_corotate) then
    call set_orbit(orbit,star(1)%mstar,star(2)%mstar,star(1)%hacc,star(2)%hacc,&
                   xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,(id==master),ierr,omega_corotate)
    add_spin = .false.
 else
    call set_orbit(orbit,star(1)%mstar,star(2)%mstar,star(1)%hacc,star(2)%hacc,&
                   xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,(id==master),ierr)
    add_spin = corotate
 endif
 if (ierr /= 0) call fatal ('setup_binary','error in call to set_orbit')
 !
 !--place stars into orbit, or add real sink particles if iprofile=0
 !
 call shift_stars(nstar,star,xyzmh_ptmass_in,vxyz_ptmass_in,xyzh,vxyzu,&
                  xyzmh_ptmass,vxyz_ptmass,npart,nptmass,corotate=add_spin)
 !
 !--restore options
 !
 iexternalforce = iextern_prev

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (iexternalforce==iext_geopot .or. iexternalforce==iext_star) then
    ! delete first sink particle and copy its properties to the central potential
    nptmass = nptmass - 1
    mass1 = xyzmh_ptmass(4,nptmass+1)
    accradius1 = xyzmh_ptmass(ihacc,nptmass+1)
    xyzmh_ptmass(:,nptmass) = xyzmh_ptmass(:,nptmass+1)
    vxyz_ptmass(:,nptmass) = vxyz_ptmass(:,nptmass+1)
 elseif (set_oblateness) then
    ! set J2 for sink particle 1 to be equal to oblateness of Saturn
    xyzmh_ptmass(iJ2,1) = 0.01629
    angle = 30.*deg_to_rad
    xyzmh_ptmass(ispinx,1) = sin(angle)
    xyzmh_ptmass(ispinz,1) = cos(angle)
    xyzmh_ptmass(iReff,1) = xyzmh_ptmass(ihacc,1)
 endif

end subroutine setpart

!----------------------------------------------------------------
!+
!  write options to .setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use setstar,      only:write_options_stars
 use setorbit,     only:write_options_orbit
 use setunits,     only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'

 call write_options_units(iunit,gr)
 call write_options_stars(star,relax,iunit)
 call write_inopt(corotate,'corotate','set stars in corotation',iunit)
 call write_options_orbit(orbit,iunit)
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
 use setstar,      only:read_options_stars
 use setorbit,     only:read_options_orbit
 use setunits,     only:read_options_and_set_units
 character(len=*), intent(in)    :: filename
 integer,          intent(inout) :: ieos
 real,             intent(inout) :: polyk
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr,need_iso
 type(inopts), allocatable :: db(:)

 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 call read_options_and_set_units(db,nerr,gr)
 call read_options_stars(star,need_iso,ieos,polyk,relax,db,nerr)
 if (need_iso==1) call fatal('setup_binary','incompatible setup for eos')
 call read_inopt(corotate,'corotate',db,errcount=nerr)
 call read_options_orbit(orbit,db,nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
