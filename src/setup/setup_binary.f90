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
!   - O            : *position angle of ascending node (deg)*
!   - a            : *semi-major axis (e.g. 1 au), period (e.g. 10*days) or rp if e=1*
!   - corotate     : *set stars in corotation*
!   - eccentricity : *eccentricity*
!   - f            : *initial true anomaly (180=apoastron)*
!   - inc          : *inclination (deg)*
!   - relax        : *relax stars into equilibrium*
!   - w            : *argument of periapsis (deg)*
!
! :Dependencies: centreofmass, dim, eos, externalforces, infile_utils, io,
!   mpidomain, options, part, physcon, relaxstar, setbinary, setstar,
!   setunits, setup_params, units
!
 use setstar, only:star_t
 use dim,     only:gr
 implicit none
 public :: setpart

 real    :: a,ecc,inc,O,w,f
 logical :: relax,corotate
 type(star_t) :: star(2)
 character(len=20) :: semi_major_axis

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
 use setbinary,      only:set_binary,get_a_from_period
 use units,          only:is_time_unit,in_code_units,utime
 use physcon,        only:solarm,au,pi,solarr,days
 use options,        only:iexternalforce
 use externalforces, only:iext_corotate,iext_geopot,iext_star,omega_corotate,mass1,accradius1
 use io,             only:master,fatal
 use setstar,        only:set_star,set_defaults_star,shift_star
 use eos,            only:X_in,Z_in,ieos
 use setup_params,   only:rhozero,npart_total
 use mpidomain,      only:i_belong
 use centreofmass,   only:reset_centreofmass
 use setunits,       only:mass_unit,dist_unit
 use physcon,        only:deg_to_rad
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
 integer :: ierr,i,nstar,nptmass_in,iextern_prev
 logical :: iexist,write_profile,use_var_comp,add_spin
 real :: xyzmh_ptmass_in(nsinkproperties,2),vxyz_ptmass_in(3,2),angle
!
!--general parameters
!
 dist_unit = 'solarr'
 mass_unit = 'solarm'
 time = 0.
 polyk = 0.
 gamma = 1.
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
 do i=1,nstar
    call set_defaults_star(star(i))
 enddo
 relax = .true.
 corotate = .false.
 semi_major_axis = '10.'
 a    = 10.
 ecc  = 0.
 inc = 0.
 O = 0.
 w = 270.
 f = 180.
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
 do i=1,nstar
    if (star(i)%iprofile > 0) then
       print "(/,a,i0,a)",' --- STAR ',i,' ---'
       call set_star(id,master,star(i),xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                     massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,&
                     X_in,Z_in,relax,use_var_comp,write_profile,&
                     rhozero,npart_total,i_belong,ierr,itype=i)
    endif
 enddo

 !
 !--if a is negative or is given time units, interpret this as a period
 !
 a = in_code_units(semi_major_axis,ierr)
 if (is_time_unit(semi_major_axis) .and. ierr == 0) then
    a = -abs(a)
    print "(a,g0,a,g0,a)",' Using PERIOD = ',abs(a),' = ',abs(a)*utime/days,' days'
 endif
 if (a < 0.) a = get_a_from_period(star(1)%mstar,star(2)%mstar,abs(a))
 !
 !--now setup orbit using fake sink particles
 !
 nptmass_in = 0
 if (iexternalforce==iext_corotate) then
    call set_binary(star(1)%mstar,star(2)%mstar,a,ecc,star(1)%hacc,star(2)%hacc,&
                    xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,ierr,omega_corotate,&
                    posang_ascnode=O,arg_peri=w,incl=inc,f=f,verbose=(id==master))
    add_spin = .false.  ! already in corotating frame
 else
    call set_binary(star(1)%mstar,star(2)%mstar,a,ecc,star(1)%hacc,star(2)%hacc,&
                    xyzmh_ptmass_in,vxyz_ptmass_in,nptmass_in,ierr,&
                    posang_ascnode=O,arg_peri=w,incl=inc,f=f,verbose=(id==master))
    add_spin = corotate
 endif
 if (ierr /= 0) call fatal ('setup_binary','error in call to set_binary')
 !
 !--place stars into orbit, or add real sink particles if iprofile=0
 !
 do i=1,nstar
    if (star(i)%iprofile > 0) then
       call shift_star(npart,xyzh,vxyzu,x0=xyzmh_ptmass_in(1:3,i),&
                       v0=vxyz_ptmass_in(1:3,i),itype=i,corotate=add_spin)
    else
       nptmass = nptmass + 1
       xyzmh_ptmass(:,nptmass) = xyzmh_ptmass_in(:,i)
       vxyz_ptmass(:,nptmass) = vxyz_ptmass_in(:,i)
    endif
 enddo
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
 else
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
 use setstar,      only:write_options_star
 use relaxstar,    only:write_options_relax
 use setunits,     only:write_options_units
 character(len=*), intent(in) :: filename
 integer :: iunit

 print "(a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for binary setup routines'

 call write_options_units(iunit,gr)
 call write_options_star(star(1),iunit,label='1')
 call write_options_star(star(2),iunit,label='2')

 write(iunit,"(/,a)") '# orbit settings'
 call write_inopt(semi_major_axis,'a','semi-major axis (e.g. 1 au), period (e.g. 10*days) or rp if e=1',iunit)
 call write_inopt(ecc,'ecc','eccentricity',iunit)
 call write_inopt(inc,'inc','inclination (deg)',iunit)
 call write_inopt(O,'O','position angle of ascending node (deg)',iunit)
 call write_inopt(w,'w','argument of periapsis (deg)',iunit)
 call write_inopt(f,'f','initial true anomaly (180=apoastron)',iunit)
 call write_inopt(corotate,'corotate','set stars in corotation',iunit)

 if (any(star(:)%iprofile > 0)) then
    write(iunit,"(/,a)") '# relaxation options'
    call write_inopt(relax,'relax','relax stars into equilibrium',iunit)
    call write_options_relax(iunit)
 endif
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
 use setstar,      only:read_options_star
 use relaxstar,    only:read_options_relax
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
 call read_options_star(star(1),need_iso,ieos,polyk,db,nerr,label='1')
 call read_options_star(star(2),need_iso,ieos,polyk,db,nerr,label='2')
 if (need_iso==1) call fatal('setup_binary','incompatible setup for eos')

 call read_inopt(semi_major_axis,'a',db,errcount=nerr)
 call read_inopt(ecc,'ecc',db,min=0.,errcount=nerr)
 call read_inopt(inc,'inc',db,errcount=nerr)
 call read_inopt(O,'O',db,errcount=nerr)
 call read_inopt(w,'w',db,errcount=nerr)
 call read_inopt(f,'f',db,errcount=nerr)

 call read_inopt(corotate,'corotate',db,errcount=nerr)

 if (any(star(:)%iprofile > 0)) then
    call read_inopt(relax,'relax',db,errcount=nerr)
    call read_options_relax(db,nerr)
 endif

 call close_db(db)
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
