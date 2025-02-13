!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setstar
!
! General routine for setting up a 3D star from a 1D profile
! This is the main functionality from setup_star but in a single routine
! that can also be called from other setups. Also contains
! routines to setup multiple stars
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - EOSopt            : *EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)*
!   - X                 : *hydrogen mass fraction*
!   - gamma             : *Adiabatic index*
!   - ieos              : *1=isothermal,2=adiabatic,10=MESA,12=idealplusrad*
!   - irecomb           : *Species to include in recombination (0: H2+H+He, 1:H+He, 2:He*
!   - metallicity       : *metallicity*
!   - mu                : *mean molecular weight*
!   - nstars            : *number of stars to add (0-'//trim(int_to_string(size(star)))//')*
!   - relax             : *relax stars into equilibrium*
!   - use_var_comp      : *Use variable composition (X, Z, mu)*
!   - write_rho_to_file : *write density profile(s) to file*
!
! :Dependencies: apr, centreofmass, dim, eos, eos_piecewise,
!   extern_densprofile, infile_utils, io, mpiutils, part, physcon,
!   prompting, radiation_utils, relaxstar, setstar_utils, unifdis, units,
!   vectorutils
!
 use setstar_utils, only:ikepler,imesa,ibpwpoly,ipoly,iuniform,ifromfile,ievrard,&
                         need_polyk,need_mu
 implicit none

 !
 ! define a data type with all the options needed
 ! to setup star (these are per-star, not per-simulation options)
 !
 type star_t
    integer :: iprofile
    integer :: isoftcore
    logical :: isinkcore
    integer :: isofteningopt
    integer :: np
    character(len=20) :: m,r,rcore,mcore,lcore,hsoft,hacc
    real :: ui_coef
    real :: polyk
    real :: initialtemp
    character(len=120) :: input_profile,dens_profile,compfile
    character(len=120) :: outputfilename ! outputfilename is the path to the cored profile
    character(len=2) :: label ! used to rename relax_star snapshots to relax1, relax2 etc.
 end type star_t

 public :: star_t
 public :: set_star,set_stars
 public :: set_defaults_star,set_defaults_stars
 public :: shift_star,shift_stars
 public :: write_options_star,write_options_stars
 public :: read_options_star,read_options_stars
 public :: set_stars_interactive
 public :: ikepler,imesa,ibpwpoly,ipoly,iuniform,ifromfile,ievrard
 public :: need_polyk

 integer, parameter :: istar_offset = 3 ! offset for particle type to distinguish particles
 ! placed in stars from other particles in the simulation

 integer, private :: EOSopt = 1

 private

contains
!--------------------------------------------------------------------------
!+
!  default options for a particular star
!  (see also set_defaults_given_profile which selects defaults
!   based on the value of iprofile)
!+
!--------------------------------------------------------------------------
subroutine set_defaults_star(star)
 type(star_t), intent(out) :: star

 star%iprofile    = 2
 star%r           = '1.0*rsun'
 star%m           = '1.0*msun'
 star%ui_coef     = 0.05
 star%polyk       = 0.
 star%initialtemp = 1.0e7
 star%isoftcore   = 0
 star%isinkcore   = .false.
 star%hsoft          = '0.0'
 star%hacc           = '1.0'
 star%rcore          = '0.0'
 star%mcore          = '0.0'
 star%lcore          = '0.*lsun'
 star%isofteningopt  = 1 ! By default, specify rcore
 star%np             = 1000
 star%input_profile  = 'P12_Phantom_Profile.data'
 star%outputfilename = 'mysoftenedstar.dat'
 star%dens_profile   = 'density.profile'
 star%compfile       = 'kepler.comp'
 star%label          = ''

end subroutine set_defaults_star

!--------------------------------------------------------------------------
!+
!  same as above but does it for multiple stars
!+
!--------------------------------------------------------------------------
subroutine set_defaults_stars(stars)
 use eos, only:use_var_comp,gmw,X_in,Z_in
 type(star_t), intent(out) :: stars(:)
 integer :: i

 EOSopt      = 1
 gmw         = 0.5988
 X_in        = 0.74
 Z_in        = 0.02
 use_var_comp = .false.
 do i=1,size(stars)
    call set_defaults_star(stars(i))
 enddo

end subroutine set_defaults_stars

!--------------------------------------------------------------------------
!+
!  utility routine to convert properties to code units while checking
!  for errors
!+
!--------------------------------------------------------------------------
subroutine check_and_convert(var,desc,unit_type,val,nerr)
 use units, only:in_code_units
 use io,    only:error
 character(len=*), intent(in)    :: var,desc,unit_type
 real,             intent(out)   :: val
 integer,          intent(inout) :: nerr
 integer :: ierr

 val = in_code_units(var,ierr,unit_type=unit_type)
 if (ierr /= 0) then
    call error('set_star','error parsing units for '//trim(desc),var=var)
    nerr = nerr + 1
 endif

end subroutine check_and_convert

!--------------------------------------------------------------------------
!+
!  utility routine to convert the stellar properties to code units
!+
!--------------------------------------------------------------------------
subroutine get_star_properties_in_code_units(star,rstar,mstar,rcore,mcore,hsoft,lcore,hacc,nerr)
 type(star_t), intent(in)    :: star
 real,         intent(out)   :: rstar,mstar,rcore,mcore,hsoft,lcore,hacc
 integer,      intent(inout) :: nerr

 call check_and_convert(star%r,'stellar radius','length',rstar,nerr)
 call check_and_convert(star%m,'stellar mass','mass',mstar,nerr)
 call check_and_convert(star%rcore,'rcore','length',rcore,nerr)
 call check_and_convert(star%mcore,'mcore','mass',mcore,nerr)
 call check_and_convert(star%hsoft,'core softening','length',hsoft,nerr)
 call check_and_convert(star%lcore,'core luminosity','luminosity',lcore,nerr)
 call check_and_convert(star%hacc,'accretion radius','length',hacc,nerr)

end subroutine get_star_properties_in_code_units

!--------------------------------------------------------------------------
!+
!  Master routine to setup a star from a specified file or density profile
!+
!--------------------------------------------------------------------------
subroutine set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
                    npart,npartoftype,massoftype,hfact,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
                    relax,use_var_comp,write_rho_to_file,&
                    rhozero,npart_total,mask,ierr,x0,v0,itype,&
                    write_files,density_error,energy_error)
 use centreofmass,       only:reset_centreofmass
 use dim,                only:do_radiation,gr,gravity,maxvxyzu
 use io,                 only:fatal,error,warning
 use eos,                only:eos_outputs_mu,polyk_eos=>polyk
 use setstar_utils,      only:set_stellar_core,read_star_profile,set_star_density, &
                              set_star_composition,set_star_thermalenergy,&
                              write_kepler_comp
 use radiation_utils,    only:set_radiation_and_gas_temperature_equal
 use relaxstar,          only:relax_star
 use part,               only:ihsoft,igas,imu,set_particle_type,ilum
 use extern_densprofile, only:write_rhotab
 use unifdis,            only:mask_prototype
 use physcon,            only:pi
 use units,              only:umass,udist,utime,unit_density
 use mpiutils,           only:reduceall_mpi
 type(star_t), intent(inout)  :: star
 integer,      intent(in)     :: id,master
 integer,      intent(inout)  :: npart,npartoftype(:),nptmass
 real,         intent(inout)  :: xyzh(:,:),vxyzu(:,:),eos_vars(:,:),rad(:,:)
 real,         intent(inout)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,         intent(inout)  :: massoftype(:)
 real,         intent(in)     :: hfact
 logical,      intent(in)     :: relax,use_var_comp,write_rho_to_file
 integer,      intent(in)     :: ieos
 real,         intent(inout)  :: gamma
 real,         intent(in)     :: X_in,Z_in
 real,         intent(out)    :: rhozero
 integer(kind=8), intent(out) :: npart_total
 integer,      intent(out)    :: ierr
 real,         intent(in),  optional :: x0(3),v0(3)
 integer,      intent(in),  optional :: itype
 logical,      intent(in),  optional :: write_files
 real,         intent(out), optional :: density_error,energy_error
 procedure(mask_prototype)      :: mask
 integer                        :: npts,ierr_relax,nerr
 integer                        :: ncols_compo,npart_old,i
 real, allocatable              :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),mu(:)
 real, allocatable              :: composition(:,:)
 real                           :: rmin,rhocentre,rmserr,en_err
 real                           :: rstar,mstar,rcore,mcore,hsoft,lcore,hacc
 character(len=20), allocatable :: comp_label(:)
 character(len=30)              :: lattice  ! The lattice type if stretchmap is used
 logical                        :: use_exactN,composition_exists,write_dumps

 use_exactN = .true.
 composition_exists = .false.
 ierr_relax = 0
 rhozero = 0.
 npart_old = npart
 write_dumps = .true.
 ierr = 0
 if (present(write_files)) write_dumps = write_files
 !
 ! do nothing if iprofile is invalid or zero (sink particle)
 !
 if (star%iprofile <= 0) then
    ierr = 1
    return
 endif
 !
 ! perform unit conversion on input quantities (if necessary)
 !
 nerr = 0
 call get_star_properties_in_code_units(star,rstar,mstar,rcore,mcore,hsoft,lcore,hacc,nerr)
 if (nerr /= 0) then
    ierr = 2
    return
 endif
 !
 ! get the desired tables of density, pressure, temperature and composition
 ! as a function of radius / mass fraction
 !
 call read_star_profile(star%iprofile,ieos,star%input_profile,gamma,star%polyk,&
                        star%ui_coef,r,den,pres,temp,en,mtab,X_in,Z_in,Xfrac,Yfrac,mu,&
                        npts,rmin,rstar,mstar,rhocentre,&
                        star%isoftcore,star%isofteningopt,rcore,mcore,&
                        hsoft,star%outputfilename,composition,&
                        comp_label,ncols_compo)

 !
 ! set up particles to represent the desired stellar profile
 !
 lattice = 'closepacked'
 if (relax) lattice='random'
 if (star%np < 1 .and. npart_old==0) then
    call fatal('set_star','cannot set up a star with zero particles')
    ierr = 2
    return
 endif
 if (mstar < 0.) then
    call fatal('set_star','cannot set up a star with negative mass!')
    ierr = 2
    return
 endif
 call set_star_density(lattice,id,master,rmin,rstar,mstar,hfact,&
                       npts,den,r,npart,npartoftype,massoftype,xyzh,use_exactN,&
                       star%np,rhozero,npart_total,mask)
 !
 ! die if stupid things done with GR
 !
 if (gr) then
    if (rstar < 6.*mstar) call fatal('set_star','R < 6GM/c^2 for star in GR violates weak field assumption')
 endif

 !
 ! add sink particle stellar core
 !
 if (star%isinkcore) call set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,&
                                           mcore,hsoft,ilum,lcore,ierr)
 if (ierr==1) call fatal('set_stellar_core','mcore <= 0')
 if (ierr==2) call fatal('set_stellar_core','hsoft <= 0')
 if (ierr==3) call fatal('set_stellar_core','lcore < 0')
 !
 ! Write the desired profile to file (do this before relaxation)
 !
 if (write_rho_to_file) call write_rhotab(star%dens_profile,&
                                          r,den,npts,star%polyk,gamma,rhocentre,ierr)
 !
 ! mask any existing particles as accreted so they are
 ! excluded from the centre of mass and relax_star calculations
 !
 if (npart_old > 0) xyzh(4,1:npart_old) = -abs(xyzh(4,1:npart_old))
 !
 ! relax the density profile to achieve nice hydrostatic equilibrium
 !
 if (relax) then
    if (reduceall_mpi('+',npart)==npart) then
       polyk_eos = star%polyk
       call relax_star(npts,den,pres,r,npart,xyzh,use_var_comp,Xfrac,Yfrac,&
                       mu,ierr_relax,npin=npart_old,label=star%label,&
                       write_dumps=write_dumps,density_error=rmserr,energy_error=en_err)
       if (present(density_error)) density_error = rmserr
       if (present(energy_error)) energy_error = en_err
    else
       call error('setup_star','cannot run relaxation with MPI setup, please run setup on ONE MPI thread')
    endif
 endif
 !
 ! Reset centre of mass (again)
 !
 call reset_centreofmass(npart,xyzh,vxyzu)

 !
 ! restore previous particles
 !
 if (npart_old > 0) xyzh(4,1:npart_old) = abs(xyzh(4,1:npart_old))

 !
 ! set composition (X,Z,mu, if using variable composition)
 ! of each particle by interpolating from table
 !
 if (use_var_comp .or. eos_outputs_mu(ieos)) then
    call set_star_composition(use_var_comp,eos_outputs_mu(ieos),npart,&
                              xyzh,Xfrac,Yfrac,mu,mtab,mstar,eos_vars,npin=npart_old)
 endif
 !
 ! Write .comp file containing composition of each particle after interpolation
 !
 if (star%iprofile==iKepler) then
    call write_kepler_comp(star%compfile,composition,comp_label,ncols_compo,r,&
                           xyzh,npart,npts,composition_exists,npin=npart_old)
 endif
 !
 ! set the internal energy and temperature
 !
 if (maxvxyzu==4) call set_star_thermalenergy(ieos,den,pres,r,npts,npart,&
                       xyzh,vxyzu,rad,eos_vars,relax,use_var_comp,star%initialtemp,&
                       npin=npart_old)

 if (do_radiation) then
    if (eos_outputs_mu(ieos)) then
       call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,&
                                                    rad,eos_vars(imu,:),npin=npart_old)
    else
       call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad,&
                                                    npin=npart_old)
    endif
 endif

 !
 ! shift star to requested position and velocity
 !
 if (present(x0)) then
    do i=npart_old+1,npart
       xyzh(1:3,i) = xyzh(1:3,i) + x0(:)
    enddo
 endif
 if (present(v0)) then
    do i=npart_old+1,npart
       vxyzu(1:3,i) = vxyzu(1:3,i) + v0(:)
    enddo
 endif
 !
 ! give the particles requested particle type:
 ! this is initially used to tag particles in different stars
 ! so the stars can later be shifted into position
 !
 if (present(itype)) then
    do i=npart_old+1,npart
       call set_particle_type(i,itype+istar_offset)
    enddo
    npartoftype(itype+istar_offset) = npartoftype(itype+istar_offset) + npart - npart_old
    npartoftype(igas) = npartoftype(igas) - (npart - npart_old)
 endif
 !
 ! Print summary to screen
 !
 if (id==master) then
    write(*,"(70('='))")
    if (ieos /= 9) then
       write(*,'(1x,a,f12.5)')       'gamma               = ', gamma
    endif
    if (maxvxyzu <= 4 .and. need_polyk(star%iprofile)) then
       write(*,'(1x,a,f12.5)')       'polyk               = ', star%polyk
       write(*,'(1x,a,f12.6,a)')     'specific int. energ = ', star%polyk*rstar/mstar,' GM/R'
    endif
    call write_mass('particle mass       = ',massoftype(igas),umass)
    call write_dist('Radius              = ',rstar,udist)
    call write_mass('Mass                = ',mstar,umass)
    if (star%iprofile==ipoly) then
       write(*,'(1x,a,g0,a)') 'rho_central         = ', rhocentre*unit_density,' g/cm^3'
    endif
    write(*,'(1x,a,i12)')         'N                   = ', npart_total-npart_old
    write(*,'(1x,a,2(es12.5,a))') 'rho_mean            = ', rhozero*unit_density,  ' g/cm^3 = '&
                                                          , rhozero,               ' code units'
    write(*,'(1x,a,es12.5,a)')    'free fall time      = ', sqrt(3.*pi/(32.*rhozero))*utime,' s'
    write(*,"(70('='))")
 endif

 if ( (star%iprofile==iuniform .and. .not.gravity) .or. star%iprofile==ifromfile) then
    call warning('setup_star','This setup may not be stable')
 endif

 if (ierr_relax /= 0) call warning('setup_star','ERRORS DURING RELAXATION, SEE ABOVE!!')

end subroutine set_star

!--------------------------------------------------------------------------
!+
!  As above but loops over all stars
!+
!--------------------------------------------------------------------------
subroutine set_stars(id,master,nstars,star,xyzh,vxyzu,eos_vars,rad,&
                     npart,npartoftype,massoftype,hfact,&
                     xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,X_in,Z_in,&
                     relax,use_var_comp,write_rho_to_file,&
                     rhozero,npart_total,mask,ierr)
 use unifdis,       only:mask_prototype
 use eos,           only:init_eos,finish_eos
 use eos_piecewise, only:init_eos_piecewise_preset
 use io,            only:error
 type(star_t), intent(inout)  :: star(:)
 integer,      intent(in)     :: id,master,nstars
 integer,      intent(inout)  :: npart,npartoftype(:),nptmass
 real,         intent(inout)  :: xyzh(:,:),vxyzu(:,:),eos_vars(:,:),rad(:,:)
 real,         intent(inout)  :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 real,         intent(inout)  :: massoftype(:)
 real,         intent(in)     :: hfact
 logical,      intent(in)     :: relax,use_var_comp,write_rho_to_file
 integer,      intent(in)     :: ieos
 real,         intent(inout)  :: gamma
 real,         intent(in)     :: X_in,Z_in
 real,         intent(out)    :: rhozero
 integer(kind=8), intent(out) :: npart_total
 integer,      intent(out)    :: ierr
 procedure(mask_prototype)    :: mask
 integer  :: i

 ! initialise piecewise polytropic equation of state if piecewise polytrope used
 if (ieos==9 .or. any(star(:)%iprofile==ibpwpoly)) call init_eos_piecewise_preset(EOSopt)
 if (any(star(:)%iprofile==ibpwpoly)) call init_eos(9,ierr)

 call init_eos(ieos,ierr)
 if (ierr /= 0) then
    call error('setup','could not initialise equation of state')
    return
 endif

 do i=1,min(nstars,size(star))
    if (star(i)%iprofile > 0) then
       print "(/,a,i0,a)",' --- STAR ',i,' ---'
       call set_star(id,master,star(i),xyzh,vxyzu,eos_vars,rad,npart,npartoftype,&
                     massoftype,hfact,xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,gamma,&
                     X_in,Z_in,relax,use_var_comp,write_rho_to_file,&
                     rhozero,npart_total,mask,ierr,itype=i)
    endif
 enddo

 call finish_eos(ieos,ierr)

end subroutine set_stars

!-----------------------------------------------------------------------
!+
!  shift star to the desired position and velocity
!+
!-----------------------------------------------------------------------
subroutine shift_star(npart,npartoftype,xyz,vxyz,x0,v0,itype,corotate)
 use part,        only:get_particle_type,set_particle_type,igas
 use vectorutils, only:cross_product3D
 integer, intent(in)    :: npart
 integer, intent(inout) :: npartoftype(:)
 real, intent(inout) :: xyz(:,:),vxyz(:,:)
 real, intent(in)    :: x0(3),v0(3)
 integer, intent(in), optional :: itype
 logical, intent(in), optional :: corotate
 logical :: add_spin
 integer :: i,mytype
 real    :: omega(3),L(3),lhat(3),rcyl(3),vspin(3)

 !
 !--put star in corotating frame
 !
 add_spin = .false.
 if (present(corotate)) add_spin = corotate
 if (add_spin) then
    call cross_product3D(x0,v0,L)
    lhat   = L/sqrt(dot_product(L,L))
    rcyl   = x0 - dot_product(x0,lhat)*lhat
    omega  = L/dot_product(rcyl,rcyl)
    print*,'Adding spin to star: omega = ',omega
 endif
 if (present(itype)) print "(a,i0,a,2(es10.3,','),es10.3,a)",' MOVING STAR ',itype,' to     (x,y,z) = (',x0(1:3),')'

 over_parts: do i=1,npart
    if (present(itype)) then
       ! get type of current particle
       call get_particle_type(i,mytype)
       ! skip particles that do not match the specified type
       if (mytype /= itype+istar_offset) cycle over_parts
       ! reset type back to gas
       call set_particle_type(i,igas)
       npartoftype(mytype) = npartoftype(mytype) - 1
       npartoftype(igas) = npartoftype(igas) + 1
    endif
    xyz(1:3,i) = xyz(1:3,i) + x0(:)
    vxyz(1:3,i) = vxyz(1:3,i) + v0(:)
    if (add_spin) then
       call cross_product3D(xyz(1:3,i),omega,vspin)
       vxyz(1:3,i) = vxyz(1:3,i) + (vspin(1:3)-v0)
    endif
 enddo over_parts

end subroutine shift_star

!-----------------------------------------------------------------------
!+
!  Shifts all stars to desired positions and velocities
!+
!-----------------------------------------------------------------------
subroutine shift_stars(nstar,star,x0,v0,&
                       xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                       npart,npartoftype,nptmass,corotate)
 use part, only:ihacc,ihsoft
 integer,      intent(in)    :: nstar,npart
 type(star_t), intent(in)    :: star(nstar)
 real,         intent(in)    :: x0(3,nstar),v0(3,nstar)
 real,         intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,         intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer,      intent(inout) :: nptmass,npartoftype(:)
 logical,      intent(in), optional :: corotate
 integer :: i,ierr
 logical :: do_corotate
 real :: rstar,mstar,rcore,mcore,hsoft,lcore,hacc

 do_corotate = .false.
 if (present(corotate)) do_corotate = corotate

 do i=1,nstar
    if (star(i)%iprofile > 0) then
       call shift_star(npart,npartoftype,xyzh,vxyzu,x0=x0(1:3,i),&
                       v0=v0(1:3,i),itype=i,corotate=do_corotate)
    else
       call get_star_properties_in_code_units(star(i),rstar,mstar,rcore,mcore,hsoft,lcore,hacc,ierr)

       nptmass = nptmass + 1
       xyzmh_ptmass(1:3,nptmass)    = x0(1:3,i)
       xyzmh_ptmass(4,nptmass)      = mstar
       xyzmh_ptmass(ihsoft,nptmass) = hsoft
       xyzmh_ptmass(ihacc, nptmass) = hacc
       vxyz_ptmass(1:3,nptmass)     = v0(1:3,i)
    endif
 enddo

end subroutine shift_stars

!-----------------------------------------------------------------------
!+
!  print a distance in both code units and physical units
!+
!-----------------------------------------------------------------------
subroutine write_dist(item_in,dist_in,udist)
 use physcon, only:solarr,km
 real,             intent(in) :: dist_in
 real(kind=8),     intent(in) :: udist
 character(len=*), intent(in) :: item_in

 if ( abs(1.0-km/udist) < 1.0d-4) then
    write(*,'(1x,2(a,es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' km'
 else
    write(*,'(1x,2(a,es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in*udist/solarr,' R_sun'
 endif

end subroutine write_dist
!-----------------------------------------------------------------------
!+
!  print a mass in both code units and physical units
!+
!-----------------------------------------------------------------------
subroutine write_mass(item_in,mass_in,umass)
 use physcon, only:solarm
 real,             intent(in) :: mass_in
 real(kind=8),     intent(in) :: umass
 character(len=*), intent(in) :: item_in

 write(*,'(1x,2(a,es12.5),a)') item_in, mass_in*umass,' g      = ',mass_in*umass/solarm,' M_sun'

end subroutine write_mass

!-----------------------------------------------------------------------
!+
!  Set default options associated with a particular choice of iprofile
!  This routine should not do ANY prompting
!+
!-----------------------------------------------------------------------
subroutine set_defaults_given_profile(iprofile,filename,mstar,polyk)
 integer, intent(in)  :: iprofile
 character(len=120), intent(out) :: filename
 real,    intent(inout) :: polyk
 character(len=*), intent(inout) :: mstar

 select case(iprofile)
 case(ifromfile)
    ! Read the density profile from file (e.g. for neutron star)
    !  Original Author: Mark Bennett
    filename = 'ns-rdensity.tab'
 case(imesa)
    ! sets up a star from a 1D MESA code output
    !  Original Author: Roberto Iaconi
    filename = 'P12_Phantom_Profile.data'
 case(ikepler)
    ! sets up a star from a 1D KEPLER code output
    !  Original Author: Nicole Rodrigues and Megha Sharma
    filename = 'kepler_MS.data'
 case(ibpwpoly)
    !  piecewise polytrope
    !  Original Author: Madeline Marshall & Bernard Field
    !  Supervisors: James Wurster & Paul Lasky
    Mstar = '1.35*msun'
    polyk = 144.
 end select

end subroutine set_defaults_given_profile

!-----------------------------------------------------------------------
!+
!  interactive prompting for setting up a star
!+
!-----------------------------------------------------------------------
subroutine set_star_interactive(star)
 use prompting,     only:prompt
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 use units,         only:in_code_units
 type(star_t), intent(inout) :: star
 integer :: i

 ! Select sphere & set default values
 do i = 0, nprofile_opts
    write(*,"(i2,')',1x,a)") i, profile_opt(i)
 enddo

 call prompt('Enter which density profile to use',star%iprofile,0,nprofile_opts)
 !
 ! set default file output parameters
 !
 write(*,"('Setting up ',a)") trim(profile_opt(star%iprofile))
 call set_defaults_given_profile(star%iprofile,star%input_profile,star%m,star%polyk)

 ! resolution
 if (star%iprofile > 0) then
    star%np = 100000 ! default number of particles
    call prompt('Enter the approximate number of particles in the sphere ',star%np,0)
 endif

 ! star properties
 if (need_inputprofile(star%iprofile)) then
    call prompt('Enter file name containing input profile',star%input_profile)
 else
    call prompt('Enter the mass of the star (e.g. 1*msun)',star%m,noblank=.true.)
    if (need_rstar(star%iprofile)) then
       call prompt('Enter the radius of the star (e.g. 1*rsun)',star%r,noblank=.true.)
    endif
 endif

 select case (star%iprofile)
 case(imesa)
    print*,'Soften the core density profile and add a sink particle core?'
    print "(3(/,a))",'0: Do not soften profile', &
                     '1: Use cubic softened density profile', &
                     '2: Use constant entropy softened profile'
    call prompt('Select option above : ',star%isoftcore,0,2)

    select case(star%isoftcore)
    case(0)
       call prompt('Add a sink particle stellar core?',star%isinkcore)
       if (star%isinkcore) then
          call prompt('Enter mass of the created sink particle core (e.g. 0.1*Msun)',star%mcore)
          call prompt('Enter softening length of the sink particle core (e.g. 0.1*Rsun)',star%hsoft)
          call prompt('Enter sink particle luminosity (e.g. 1*Lsun)',star%lcore)
       endif
    case(1)
       star%isinkcore = .true. ! Create sink particle core automatically
       print*,'Options for core softening:'
       print "(5(/,a))",'1. Specify radius of density softening', &
                        '2. Specify mass of sink particle core (not recommended)', &
                        '3. Specify both radius of density softening length and mass', &
                        '   of sink particle core (if you do not know what you are', &
                        '   doing, you will obtain a poorly softened profile)'
       call prompt('Select option above : ',star%isofteningopt,1,3)

       select case(star%isofteningopt)
       case(1)
          call prompt('Enter core radius (e.g. 0.1*rsun)',star%rcore)
       case(2)
          call prompt('Enter mass of the created sink particle core (e.g. 0.1*msun)',star%mcore)
       case(3)
          call prompt('Enter mass of the created sink particle core (e.g. 0.1*msun)',star%mcore)
          call prompt('Enter core radius (e.g. 1*rsun)',star%rcore)
       end select
       call prompt('Enter sink particle luminosity [e.g. 1*Lsun]',star%lcore)

    case(2)
       star%isinkcore = .true. ! Create sink particle core automatically
       print*,'Specify core radius and initial guess for mass of sink particle core'
       call prompt('Enter core radius (e.g. 0.1*rsun): ',star%rcore)
       call prompt('Enter guess for core mass (e.g. 0.1*Msun): ',star%mcore)
       call prompt('Enter sink particle luminosity (e.g. 1.*lsun',star%lcore)
       call prompt('Enter output file name of cored stellar profile:',star%outputfilename)
    end select
 case(ievrard)
    call prompt('Enter the specific internal energy (units of GM/R) ',star%ui_coef,0.)
 case(:0)
    call prompt('Enter the accretion radius (e.g. 1.0)',star%hacc)
 end select

end subroutine set_star_interactive

!-----------------------------------------------------------------------
!+
!  as above but wrapper routine for multiple stars
!+
!-----------------------------------------------------------------------
subroutine set_stars_interactive(star,ieos,relax,nstar)
 use prompting,    only:prompt
 use infile_utils, only:int_to_string
 type(star_t), intent(inout) :: star(:)
 integer,      intent(inout) :: ieos
 logical,      intent(out)   :: relax
 integer,      intent(out), optional :: nstar
 integer :: i,nstars

 ! optionally ask for number of stars, otherwise fix nstars to the input array size
 if (present(nstar) .and. size(star) > 1) then
    nstar = 1
    call prompt('how many stars to set up (0-'//int_to_string(size(star))//')',nstar,0,size(star))
    nstars = nstar
 else
    nstars = size(star)
 endif

 do i=1,nstars
    print "(/,'------------- STAR ',i0,'-------------')",i
    call set_star_interactive(star(i))
 enddo

 ! prompt for equation of state and relaxation options if any stars made of gas
 if (nstars > 0) then
    if (any(star(1:nstars)%iprofile > 0)) then
       if (any(star(1:nstars)%iprofile==ibpwpoly)) ieos = 9 ! set default eos for piecewise polytropes
       call set_star_eos_interactive(ieos,star)

       relax = .false.
       call prompt('Relax stars automatically during setup?',relax)
    endif
 endif

end subroutine set_stars_interactive

!-----------------------------------------------------------------------
!+
!  write setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine write_options_star(star,iunit,label)
 use infile_utils,  only:write_inopt,get_optstring
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar,need_polyk
 type(star_t),     intent(in) :: star
 integer,          intent(in) :: iunit
 character(len=*), intent(in), optional :: label
 character(len=120) :: string
 character(len=10) :: c

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))

 write(iunit,"(/,a)") '# options for star '//trim(c)
 call get_optstring(nprofile_opts+1,profile_opt,string,4,from_zero=.true.)
 call write_inopt(star%iprofile,'iprofile'//trim(c),trim(string(1:48)),iunit)

 if (star%isoftcore <= 0) then
    if (need_inputprofile(star%iprofile)) then
       call write_inopt(star%input_profile,'input_profile'//trim(c),&
            'Path to input profile',iunit)
    else
       call write_inopt(star%m,'Mstar'//trim(c),'mass of star '//trim(c)//' (code units or e.g. 1*msun)',iunit)
       if (need_rstar(star%iprofile)) &
          call write_inopt(star%r,'Rstar'//trim(c),'radius of star'//trim(c)//' (code units or e.g. 1*rsun)',iunit)
    endif
 endif

 select case(star%iprofile)
 case(imesa)
    call write_inopt(star%isoftcore,'isoftcore'//trim(c),&
                     '0=no core softening, 1=cubic, 2=const. entropy',iunit)

    if (star%isoftcore > 0) then
       call write_inopt(star%input_profile,'input_profile'//trim(c),&
                        'Path to input profile',iunit)
       call write_inopt(star%outputfilename,'outputfilename'//trim(c),&
                        'Output path for softened MESA profile',iunit)

       if (star%isoftcore == 1) then
          call write_inopt(star%isofteningopt,'isofteningopt'//trim(c),&
                           '1=supply rcore, 2=supply mcore, 3=supply both',iunit)
          if ((star%isofteningopt == 1) .or. (star%isofteningopt == 3)) then
             call write_inopt(star%rcore,'rcore'//trim(c),'Radius of core softening',iunit)
          endif
          if ((star%isofteningopt == 2) .or. (star%isofteningopt == 3)) then
             call write_inopt(star%mcore,'mcore'//trim(c),'Mass of point mass stellar core',iunit)
          endif
       elseif (star%isoftcore == 2) then
          call write_inopt(star%rcore,'rcore'//trim(c),'Radius of core softening',iunit)
          call write_inopt(star%mcore,'mcore'//trim(c),&
               'Initial guess for mass of sink particle stellar core',iunit)
       endif
       call write_inopt(star%lcore,'lcore'//trim(c),&
                        'Luminosity of point mass stellar core',iunit)
    else
       call write_inopt(star%isinkcore,'isinkcore'//trim(c),&
               'Add a sink particle stellar core',iunit)
       if (star%isinkcore) then
          call write_inopt(star%mcore,'mcore'//trim(c),&
               'Mass of sink particle stellar core',iunit)
          call write_inopt(star%hsoft,'hsoft'//trim(c),&
               'Softening length of sink particle stellar core',iunit)
       endif
       call write_inopt(star%lcore,'lcore'//trim(c),'Luminosity of sink core particle',iunit)
    endif
 case (ievrard)
    call write_inopt(star%ui_coef,'ui_coef'//trim(c),&
         'specific internal energy (units of GM/R)',iunit)
 case(0)
    call write_inopt(star%hacc,'hacc'//trim(c),'accretion radius for sink'//trim(c),iunit)
 end select

 if (need_polyk(star%iprofile)) call write_inopt(star%polyk,'polyk'//trim(c),'polytropic constant (cs^2 if isothermal)',iunit)

 if (star%iprofile > 0 .and. (len_trim(c)==0 .or. c(1:1)=='1')) then
    call write_inopt(star%np,'np'//trim(c),'number of particles',iunit)
 endif

end subroutine write_options_star

!-----------------------------------------------------------------------
!+
!  read setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine read_options_star(star,db,nerr,label)
 use infile_utils,  only:inopts,read_inopt
 use setstar_utils, only:need_inputprofile,need_rstar,nprofile_opts
 type(star_t),              intent(inout) :: star
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: nerr
 character(len=*),          intent(in), optional :: label
 character(len=10) :: c
 integer :: ierr
 real    :: rstar,mstar,rcore,mcore,hsoft,lcore,hacc

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))
 star%label = trim(c)
 star%dens_profile = 'relax'//trim(c)//'.profile'
 star%compfile = 'relax'//trim(c)//'.comp'

 call read_inopt(star%iprofile,'iprofile'//trim(c),db,errcount=nerr,min=0,max=nprofile_opts)
 call set_defaults_given_profile(star%iprofile,star%input_profile,star%m,star%polyk)

 if (need_inputprofile(star%iprofile)) then
    call read_inopt(star%input_profile,'input_profile'//trim(c),db,errcount=nerr)
 endif

 ! resolution
 if (star%iprofile > 0 .and. (len_trim(c)==0 .or. c(1:1)=='1')) then
    call read_inopt(star%np,'np'//trim(c),db,errcount=nerr,min=10)
    !call read_inopt(use_exactN,'use_exactN',db,errcount=nerr)
 endif

 select case(star%iprofile)
 case(imesa)
    call read_inopt(star%isoftcore,'isoftcore'//trim(c),db,errcount=nerr,min=0)

    if (star%isoftcore <= 0) then ! sink particle core without softening
       call read_inopt(star%isinkcore,'isinkcore'//trim(c),db,errcount=nerr)
       if (star%isinkcore) then
          call read_inopt(star%mcore,'mcore'//trim(c),db,errcount=nerr,err=ierr)
          call read_inopt(star%hsoft,'hsoft'//trim(c),db,errcount=nerr,err=ierr)
       endif
    else
       star%isinkcore = .true.
       call read_inopt(star%outputfilename,'outputfilename'//trim(c),db,errcount=nerr)
       if (star%isoftcore==2) then
          star%isofteningopt=3
       elseif (star%isoftcore==1) then
          call read_inopt(star%isofteningopt,'isofteningopt'//trim(c),db,errcount=nerr,min=0)
       endif

       if ((star%isofteningopt==1) .or. (star%isofteningopt==3)) then
          call read_inopt(star%rcore,'rcore'//trim(c),db,errcount=nerr,err=ierr)
       endif
       if ((star%isofteningopt==2) .or. (star%isofteningopt==3) &
           .or. (star%isoftcore==2)) then
          call read_inopt(star%mcore,'mcore'//trim(c),db,errcount=nerr,err=ierr)
       endif
    endif

    if (star%isinkcore) then
       call read_inopt(star%lcore,'lcore'//trim(c),db,errcount=nerr,err=ierr)
    endif
 case(ievrard)
    call read_inopt(star%ui_coef,'ui_coef'//trim(c),db,errcount=nerr)
 case(:0)
    call read_inopt(star%hacc,'hacc'//trim(c),db,errcount=nerr)
 end select

 if (need_polyk(star%iprofile)) call read_inopt(star%polyk,'polyk'//trim(c),db,errcount=nerr)

 ! star properties
 if (star%isoftcore <= 0) then
    if (need_inputprofile(star%iprofile)) then
       call read_inopt(star%input_profile,'input_profile'//trim(c),db,errcount=nerr)
    else
       call read_inopt(star%m,'Mstar'//trim(c),db,errcount=nerr,err=ierr)
       if (need_rstar(star%iprofile)) then
          call read_inopt(star%r,'Rstar'//trim(c),db,errcount=nerr,err=ierr)
       endif
    endif
 endif

 ! perform a unit conversion, just to check that there are no errors parsing the .setup file
 if (nerr==0) call get_star_properties_in_code_units(star,rstar,mstar,rcore,mcore,hsoft,lcore,hacc,nerr)

end subroutine read_options_star

!-----------------------------------------------------------------------
!+
!  write_options routine that writes options for multiple stars
!+
!-----------------------------------------------------------------------
subroutine write_options_stars(star,relax,write_rho_to_file,ieos,iunit,nstar)
 use relaxstar,    only:write_options_relax
 use infile_utils, only:write_inopt,int_to_string
 use apr,          only:use_apr,write_options_apr
 type(star_t), intent(in) :: star(:)
 integer,      intent(in) :: ieos,iunit
 logical,      intent(in) :: relax,write_rho_to_file
 integer,      intent(in), optional :: nstar
 integer :: i,nstars
 character(len=3) :: label(size(star))

 ! optionally ask for number of stars, otherwise fix nstars to the input array size
 if (present(nstar)) then
    call write_inopt(nstar,'nstars','number of stars to add (0-'//trim(int_to_string(size(star)))//')',iunit)
    nstars = nstar
 else
    nstars = size(star)
 endif

 ! write options for each star
 do i=1,nstars
    label(i) = trim(int_to_string(i))
    call write_options_star(star(i),iunit,label=label(i))
 enddo

 ! write equation of state and relaxation options if any stars made of gas
 if (nstars > 0) then
    if (any(star(1:nstars)%iprofile > 0)) then
       call write_options_stars_eos(nstars,star(1:nstars),label(1:nstars),ieos,iunit)

       write(iunit,"(/,a)") '# relaxation options'
       call write_inopt(relax,'relax','relax stars into equilibrium',iunit)
       call write_options_relax(iunit)
       call write_inopt(write_rho_to_file,'write_rho_to_file','write density profile(s) to file',iunit)
    endif
 endif

 if (use_apr) call write_options_apr(iunit)

end subroutine write_options_stars

!-----------------------------------------------------------------------
!+
!  read_options routine that reads options for multiple stars
!+
!-----------------------------------------------------------------------
subroutine read_options_stars(star,ieos,relax,write_rho_to_file,db,nerr,nstar)
 use relaxstar,    only:read_options_relax
 use infile_utils, only:inopts,read_inopt,int_to_string
 use apr,          only:use_apr,apr_max_in,ref_dir,apr_type,apr_rad,apr_drad
 type(star_t),              intent(inout) :: star(:) ! inout because can set default options manually in calling routine
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: ieos
 logical,                   intent(out)   :: relax,write_rho_to_file
 integer,                   intent(inout) :: nerr
 integer,                   intent(out), optional :: nstar
 integer :: i,nstars
 character(len=3) :: label(size(star))

 ! optionally ask for number of stars
 if (present(nstar)) then
    call read_inopt(nstar,'nstars',db,errcount=nerr,min=0,max=size(star))
    nstars = nstar
 else
    nstars = size(star)
 endif

 ! read options for each star
 do i=1,nstars
    label(i) = trim(int_to_string(i))
    call read_options_star(star(i),db,nerr,label=label(i))
 enddo

 ! equation of state and relaxation options if any stars made of gas
 if (nstars > 0) then
    if (any(star(1:nstars)%iprofile > 0)) then
       if (any(star(1:nstars)%iprofile==ibpwpoly)) ieos = 9 ! set default eos for piecewise polytropes
       ! equation of state options
       call read_options_stars_eos(nstars,star(1:nstars),label(1:nstars),ieos,db,nerr)
       ! relaxation options
       call read_inopt(relax,'relax',db,errcount=nerr)
       call read_options_relax(db,nerr)
       ! option to write density profile to file
       call read_inopt(write_rho_to_file,'write_rho_to_file',db,errcount=nerr)
    endif
 endif

 if (use_apr) then
    call read_inopt(apr_max_in,'apr_max',db,errcount=nerr)
    call read_inopt(ref_dir,'ref_dir',db,errcount=nerr)
    call read_inopt(apr_type,'apr_type',db,errcount=nerr)
    call read_inopt(apr_rad,'apr_rad',db,errcount=nerr)
    call read_inopt(apr_drad,'apr_drad',db,errcount=nerr)
 endif

end subroutine read_options_stars

!-----------------------------------------------------------------------
!+
!  write equation of state options needed to setup stars
!+
!-----------------------------------------------------------------------
subroutine write_options_stars_eos(nstars,star,label,ieos,iunit)
 use eos,          only:use_var_comp,X_in,Z_in,irecomb,gmw,gamma
 use infile_utils, only:write_inopt
 integer, intent(in) :: nstars,ieos,iunit
 type(star_t),     intent(in) :: star(nstars)
 character(len=*), intent(in) :: label(nstars)
 integer :: i

 write(iunit,"(/,a)") '# equation of state used to set the thermal energy profile'
 call write_inopt(ieos,'ieos','1=isothermal,2=adiabatic,10=MESA,12=idealplusrad',iunit)

 if (any(star(:)%iprofile==imesa)) then
    call write_inopt(use_var_comp,'use_var_comp','Use variable composition (X, Z, mu)',iunit)
 endif

 select case(ieos)
 case(9)
    write(iunit,"(/,a)") '# Piecewise Polytrope default options'
    call write_inopt(EOSopt,'EOSopt','EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)',iunit)
 case(2,12)
    call write_inopt(gamma,'gamma','Adiabatic index',iunit)
    if (any(need_mu(star(:)%isoftcore)) .and. (.not. use_var_comp)) call write_inopt(gmw,'mu','mean molecular weight',iunit)
 case(10,20)
    if (ieos==20) call write_inopt(irecomb,'irecomb','Species to include in recombination (0: H2+H+He, 1:H+He, 2:He',iunit)
    if ( (.not. use_var_comp) .and. any(need_mu(star(:)%isoftcore))) then
       call write_inopt(X_in,'X','hydrogen mass fraction',iunit)
       call write_inopt(Z_in,'Z','metallicity',iunit)
    endif
 case(15)
    ! options for setting initial thermal energy (e.g. if degenerate matter eos)
    do i=1,nstars
       if (star(i)%iprofile > 0) then
          call write_inopt(star(i)%initialtemp,'initialtemp'//trim(label(i)),&
                          'initial temperature of star (e.g. if degenerate matter eos)',iunit)
       endif
    enddo
 end select

end subroutine write_options_stars_eos

!-----------------------------------------------------------------------
!+
!  read equation of state options needed to setup stars
!+
!-----------------------------------------------------------------------
subroutine read_options_stars_eos(nstars,star,label,ieos,db,nerr)
 use eos,          only:use_var_comp,X_in,Z_in,irecomb,gamma,gmw
 use infile_utils, only:inopts,read_inopt
 integer,                   intent(in)    :: nstars
 type(star_t),              intent(inout) :: star(nstars)
 character(len=*),          intent(in)    :: label(nstars)
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: ieos
 integer,                   intent(inout) :: nerr
 integer :: i

 ! equation of state
 call read_inopt(ieos,'ieos',db,errcount=nerr)
 if (any(star(:)%iprofile==imesa)) call read_inopt(use_var_comp,'use_var_comp',db,errcount=nerr)

 select case(ieos)
 case(9)
    call read_inopt(EOSopt,'EOSopt',db,min=0,max=4,errcount=nerr)
 case(2,12)
    call read_inopt(gamma,'gamma',db,min=1.,max=7.,errcount=nerr)
    if ((.not. use_var_comp) .and. any(need_mu(star(:)%isoftcore))) then
       call read_inopt(gmw,'mu',db,min=0.,errcount=nerr)
    endif
 case(10,20)
    if (ieos==20) call read_inopt(irecomb,'irecomb',db,errcount=nerr)
    ! if softening stellar core, composition is automatically determined at R/2
    if ((.not. use_var_comp) .and. any(need_mu(star(:)%isoftcore))) then
       call read_inopt(X_in,'X',db,min=0.,max=1.,errcount=nerr)
       call read_inopt(Z_in,'Z',db,min=0.,max=1.,errcount=nerr)
    endif
 case(15)
    do i=1,nstars
       if (star(i)%iprofile > 0) then
          call read_inopt(star(i)%initialtemp,'initialtemp'//trim(label(i)),&
                          db,min=0.,max=1e12,errcount=nerr)
       endif
    enddo
 end select

end subroutine read_options_stars_eos

!------------------------------------------------------------------------
!+
!  interactive prompt for equation of state options needed to setup stars
!+
!------------------------------------------------------------------------
subroutine set_star_eos_interactive(ieos,star)
 use prompting, only:prompt
 use eos,       only:use_var_comp,X_in,Z_in,irecomb,gamma,gmw
 integer,      intent(inout) :: ieos
 type(star_t), intent(in)    :: star(:)

 ! equation of state
 call prompt('Enter the desired EoS (1=isothermal,2=adiabatic,10=MESA,12=idealplusrad)',ieos)
 if (any(star(:)%iprofile==imesa)) call prompt('Use variable composition?',use_var_comp)

 select case(ieos)
 case(9)
    write(*,'(a)') 'EOS options: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)'
    call prompt('Enter equation of state type',EOSopt,1,4)
 case(2,12)
    call prompt('Enter gamma (adiabatic index)',gamma,1.,7.)
    if ( (.not. use_var_comp) .and. any(need_mu(star(:)%isoftcore))) then
       call prompt('Enter mean molecular weight',gmw,0.)
    endif
 case(10,20)
    if (ieos==20) call prompt('Enter irecomb (0: H2+H+He, 1:H+He, 2:He)',irecomb,0)
    if ( (.not. use_var_comp) .and. any(need_mu(star(:)%isoftcore))) then
       call prompt('Enter hydrogen mass fraction (X)',X_in,0.,1.)
       call prompt('Enter metals mass fraction (Z)',Z_in,0.,1.)
    endif
 end select

end subroutine set_star_eos_interactive

end module setstar
