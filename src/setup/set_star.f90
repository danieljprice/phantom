!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setstar
!
! General routine for setting up a 3D star from a 1D profile
! This is the main functionality from setup_star but in a single routine
! that can also be called from other setups. In principle this
! could also be used to setup multiple stars
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: centreofmass, dim, eos, extern_densprofile, infile_utils,
!   io, mpiutils, part, physcon, prompting, radiation_utils, relaxstar,
!   setstar_utils, unifdis, units, vectorutils
!
 use setstar_utils, only:ikepler,imesa,ibpwpoly,ipoly,iuniform,ifromfile,ievrard,&
                         need_polyk
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
    real :: Rstar
    real :: Mstar
    real :: ui_coef
    real :: initialtemp
    real :: rcore
    real :: mcore
    real :: lcore
    real :: hsoft
    real :: hacc   ! accretion radius if star is a sink particle
    character(len=120) :: input_profile,dens_profile
    character(len=120) :: outputfilename ! outputfilename is the path to the cored profile
    character(len=2) :: label ! used to rename relax_star snapshots to relax1, relax2 etc.
 end type star_t

 public :: star_t
 public :: set_star,set_defaults_star,shift_star
 public :: write_options_star,read_options_star,set_star_interactive
 public :: ikepler,imesa,ibpwpoly,ipoly,iuniform,ifromfile,ievrard
 public :: need_polyk

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
 use units,   only:udist,umass
 use physcon, only:solarm,solarr
 type(star_t), intent(out) :: star

 star%iprofile    = 2
 star%rstar       = 1.0*real(solarr/udist)
 star%mstar       = 1.0*real(solarm/umass)
 star%ui_coef     = 0.05
 star%initialtemp = 1.0e7
 star%isoftcore   = 0
 star%isinkcore   = .false.
 star%hsoft          = 0.
 star%hacc           = 1.
 star%rcore          = 0.
 star%mcore          = 0.
 star%lcore          = 0.
 star%isofteningopt  = 1 ! By default, specify rcore
 star%np             = 1000
 star%input_profile  = 'P12_Phantom_Profile.data'
 star%outputfilename = 'mysoftenedstar.dat'
 star%dens_profile   = 'density-profile.tab'
 star%label          = ''

end subroutine set_defaults_star

!--------------------------------------------------------------------------
!+
!  Master routine to setup a star from a specified file or density profile
!+
!--------------------------------------------------------------------------
subroutine set_star(id,master,star,xyzh,vxyzu,eos_vars,rad,&
                    npart,npartoftype,massoftype,hfact,&
                    xyzmh_ptmass,vxyz_ptmass,nptmass,ieos,polyk,gamma,X_in,Z_in,&
                    relax,use_var_comp,write_rho_to_file,&
                    rhozero,npart_total,mask,ierr,x0,v0,itype)
 use centreofmass,       only:reset_centreofmass
 use dim,                only:do_radiation,gr,gravity,maxvxyzu
 use io,                 only:fatal,error,warning
 use eos,                only:eos_outputs_mu
 use setstar_utils,      only:set_stellar_core,read_star_profile,set_star_density, &
                              set_star_composition,set_star_thermalenergy,&
                              write_kepler_comp
 use radiation_utils,    only:set_radiation_and_gas_temperature_equal
 use relaxstar,          only:relax_star
 use part,               only:ihsoft,igas,imu,set_particle_type,ilum
 use extern_densprofile, only:write_rhotab
 use unifdis,            only:mask_prototype
 use physcon,            only:pi
 use units,              only:utime,udist,umass,unit_density
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
 real,         intent(inout)  :: polyk,gamma
 real,         intent(in)     :: X_in,Z_in
 real,         intent(out)    :: rhozero
 integer(kind=8), intent(out) :: npart_total
 integer,      intent(out)    :: ierr
 real,         intent(in), optional :: x0(3),v0(3)
 integer,      intent(in), optional :: itype
 procedure(mask_prototype)    :: mask
 integer                        :: npts,ierr_relax
 integer                        :: ncols_compo,npart_old,i
 real, allocatable              :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),mu(:)
 real, allocatable              :: composition(:,:)
 real                           :: rmin,rhocentre
 character(len=20), allocatable :: comp_label(:)
 character(len=30)              :: lattice  ! The lattice type if stretchmap is used
 logical                        :: use_exactN,composition_exists

 use_exactN = .true.
 composition_exists = .false.
 ierr_relax = 0
 rhozero = 0.
 npart_old = npart
 !
 ! do nothing if iprofile is invalid or zero (sink particle)
 !
 if (star%iprofile <= 0) then
    ierr = 1
    return
 endif
 !
 ! get the desired tables of density, pressure, temperature and composition
 ! as a function of radius / mass fraction
 !
 call read_star_profile(star%iprofile,ieos,star%input_profile,gamma,polyk,&
                        star%ui_coef,r,den,pres,temp,en,mtab,X_in,Z_in,Xfrac,Yfrac,mu,&
                        npts,rmin,star%rstar,star%mstar,rhocentre,&
                        star%isoftcore,star%isofteningopt,star%rcore,star%mcore,&
                        star%hsoft,star%outputfilename,composition,&
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
 if (star%mstar < 0.) then
    call fatal('set_star','cannot set up a star with negative mass!')
    ierr = 2
    return
 endif
 call set_star_density(lattice,id,master,rmin,star%rstar,star%mstar,hfact,&
                       npts,den,r,npart,npartoftype,massoftype,xyzh,use_exactN,&
                       star%np,rhozero,npart_total,mask)
 !
 ! die if stupid things done with GR
 !
 if (gr) then
    if (star%rstar < 6.*star%mstar) call fatal('set_star','R < 6GM/c^2 for star in GR violates weak field assumption')
 endif

 !
 ! add sink particle stellar core
 !
 if (star%isinkcore) call set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,&
                                           star%mcore,star%hsoft,ilum,star%lcore,ierr)
 if (ierr==1) call fatal('set_stellar_core','mcore <= 0')
 if (ierr==2) call fatal('set_stellar_core','hsoft <= 0')
 if (ierr==3) call fatal('set_stellar_core','lcore < 0')
 !
 ! Write the desired profile to file (do this before relaxation)
 !
 if (write_rho_to_file) call write_rhotab(star%dens_profile,&
                                          r,den,npts,polyk,gamma,rhocentre,ierr)
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
       call relax_star(npts,den,pres,r,npart,xyzh,use_var_comp,Xfrac,Yfrac,&
                       mu,ierr_relax,npin=npart_old,label=star%label)
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
                              xyzh,Xfrac,Yfrac,mu,mtab,star%mstar,eos_vars,npin=npart_old)
 endif
 !
 ! Write composition file called kepler.comp containing composition of each particle after interpolation
 !
 if (star%iprofile==iKepler) call write_kepler_comp(composition,comp_label,ncols_compo,r,&
                                  xyzh,npart,npts,composition_exists,npin=npart_old)
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
       call set_particle_type(i,itype)
    enddo
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
       write(*,'(1x,a,f12.5)')       'polyk               = ', polyk
       write(*,'(1x,a,f12.6,a)')     'specific int. energ = ', polyk*star%rstar/star%mstar,' GM/R'
    endif
    call write_mass('particle mass       = ',massoftype(igas),umass)
    call write_dist('Radius              = ',star%rstar,udist)
    call write_mass('Mass                = ',star%mstar,umass)
    if (star%iprofile==ipoly) then
       write(*,'(1x,a,g0,a)') 'rho_central         = ', rhocentre*unit_density,' g/cm^3'
    endif
    write(*,'(1x,a,i12)')         'N                   = ', npart_total-npart_old
    write(*,'(1x,a,2(es12.5,a))') 'rho_mean            = ', rhozero*unit_density,  ' g/cm^3 = '&
                                                          , rhozero,               ' code units'
    write(*,'(1x,a,es12.5,a)')    'free fall time      = ', sqrt(3.*pi/(32.*rhozero))*utime,' s'

    if (composition_exists) then
       write(*,'(a)') 'Composition written to kepler.comp file.'
    endif
    write(*,"(70('='))")
 endif

 if ( (star%iprofile==iuniform .and. .not.gravity) .or. star%iprofile==ifromfile) then
    call warning('setup_star','This setup may not be stable')
 endif

 if (ierr_relax /= 0) call warning('setup_star','ERRORS DURING RELAXATION, SEE ABOVE!!')

end subroutine set_star

!-----------------------------------------------------------------------
!+
!  shift star to the desired position and velocity
!+
!-----------------------------------------------------------------------
subroutine shift_star(npart,xyz,vxyz,x0,v0,itype,corotate)
 use part,        only:get_particle_type,set_particle_type,igas
 use vectorutils, only:cross_product3D
 integer, intent(in) :: npart
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

 over_parts: do i=1,npart
    if (present(itype)) then
       ! get type of current particle
       call get_particle_type(i,mytype)
       ! skip particles that do not match the specified type
       if (mytype /= itype) cycle over_parts
       ! reset type back to gas
       call set_particle_type(i,igas)
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
subroutine set_defaults_given_profile(iprofile,filename,need_iso,ieos,mstar,polyk)
 integer, intent(in)  :: iprofile
 character(len=120), intent(out) :: filename
 integer, intent(out) :: need_iso
 integer, intent(inout) :: ieos
 real,    intent(inout) :: mstar,polyk

 need_iso = 0
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
    ieos         = 9
    !dist_unit    = 'km'
    Mstar        = 1.35
    polyk        = 144.
 case(ievrard)
    need_iso    = -1
 end select

end subroutine set_defaults_given_profile

!-----------------------------------------------------------------------
!+
!  interactive prompting for setting up a star
!+
!-----------------------------------------------------------------------
subroutine set_star_interactive(id,master,star,need_iso,use_var_comp,ieos,polyk)
 use prompting,     only:prompt
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 use units,         only:in_solarm,in_solarr,in_solarl,udist,umass,unit_luminosity
 use physcon,       only:solarr,solarm,solarl
 type(star_t), intent(out)   :: star
 integer,      intent(in)    :: id,master
 logical,      intent(out)   :: use_var_comp
 integer,      intent(out)   :: need_iso
 integer,      intent(inout) :: ieos
 real,         intent(inout) :: polyk
 integer :: i
 real :: mstar_msun,rstar_rsun,rcore_rsun,mcore_msun,lcore_lsun,hsoft_rsun

 ! set defaults
 call set_defaults_star(star)
 mstar_msun = real(in_solarm(star%mstar))
 rstar_rsun = real(in_solarr(star%rstar))
 mcore_msun = real(in_solarm(star%mcore))
 lcore_lsun = real(in_solarl(star%lcore))
 rcore_rsun = real(in_solarr(star%rcore))
 hsoft_rsun = real(in_solarr(star%hsoft))

 ! Select sphere & set default values
 do i = 1, nprofile_opts
    write(*,"(i2,')',1x,a)") i, profile_opt(i)
 enddo

 call prompt('Enter which density profile to use',star%iprofile,1,nprofile_opts)
 !
 ! set default file output parameters
 !
 if (id==master) write(*,"('Setting up ',a)") trim(profile_opt(star%iprofile))
 call set_defaults_given_profile(star%iprofile,star%input_profile,&
                                 need_iso,ieos,star%mstar,polyk)

 ! resolution
 if (star%iprofile > 0) then
    star%np = 100000 ! default number of particles
    call prompt('Enter the approximate number of particles in the sphere ',star%np,0)
 endif

 ! star properties
 if (need_inputprofile(star%iprofile)) then
    call prompt('Enter file name containing input profile',star%input_profile)
 else
    call prompt('Enter the mass of the star (Msun)',mstar_msun,0.)
    star%mstar = mstar_msun*real(solarm/umass)
    if (need_rstar(star%iprofile)) then
       call prompt('Enter the radius of the star (Rsun)',rstar_rsun,0.)
       star%rstar = rstar_rsun*real(solarr/udist)
    endif
 endif

 select case (star%iprofile)
 case(imesa)
    call prompt('Use variable composition?',use_var_comp)

    print*,'Soften the core density profile and add a sink particle core?'
    print "(3(/,a))",'0: Do not soften profile', &
                     '1: Use cubic softened density profile', &
                     '2: Use constant entropy softened profile'
    call prompt('Select option above : ',star%isoftcore,0,2)

    select case(star%isoftcore)
    case(0)
       call prompt('Add a sink particle stellar core?',star%isinkcore)
       if (star%isinkcore) then
          call prompt('Enter mass of the created sink particle core [Msun]',mcore_msun,0.)
          call prompt('Enter softening length of the sink particle core [Rsun]',hsoft_rsun,0.)
          call prompt('Enter sink particle luminosity [Lsun]',lcore_lsun,0.)
          star%mcore = mcore_msun*real(solarm/umass)
          star%hsoft = hsoft_rsun*real(solarr/udist)
          star%lcore = lcore_lsun*real(solarl/unit_luminosity)
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
          call prompt('Enter core radius [Rsun]',rcore_rsun,0.)
          star%rcore = rcore_rsun*real(solarr/udist)
       case(2)
          call prompt('Enter mass of the created sink particle core [Msun]',mcore_msun,0.)
          star%mcore = mcore_msun*real(solarm/umass)
       case(3)
          call prompt('Enter mass of the created sink particle core [Msun]',mcore_msun,0.)
          call prompt('Enter core radius [Rsun]',rcore_rsun,0.)
          star%mcore = mcore_msun*real(solarm/umass)
          star%rcore = rcore_rsun*real(solarr/udist)
       end select
       call prompt('Enter sink particle luminosity [Lsun]',lcore_lsun,0.)
       star%lcore = lcore_lsun*real(solarl/unit_luminosity)

    case(2)
       star%isinkcore = .true. ! Create sink particle core automatically
       print*,'Specify core radius and initial guess for mass of sink particle core'
       call prompt('Enter core radius in Rsun : ',rcore_rsun,0.)
       call prompt('Enter guess for core mass in Msun : ',mcore_msun,0.)
       call prompt('Enter sink particle luminosity [Lsun]',lcore_lsun,0.)
       call prompt('Enter output file name of cored stellar profile:',star%outputfilename)
       star%mcore = mcore_msun*real(solarm/umass)
       star%rcore = rcore_rsun*real(solarr/udist)
       star%lcore = lcore_lsun*real(solarl/unit_luminosity)
    end select
 case(ievrard)
    call prompt('Enter the specific internal energy (units of GM/R) ',star%ui_coef,0.)
 case(:0)
    call prompt('Enter the accretion radius in code units',star%hacc,0.)
 end select

end subroutine set_star_interactive

!-----------------------------------------------------------------------
!+
!  write setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine write_options_star(star,iunit,label)
 use infile_utils,  only:write_inopt,get_optstring
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 use units,         only:in_solarm,in_solarr,in_solarl
 type(star_t),     intent(in) :: star
 integer,          intent(in) :: iunit
 character(len=*), intent(in), optional :: label
 character(len=120) :: string
 character(len=10) :: c

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))

 write(iunit,"(/,a)") '# options for star '//trim(c)
 call get_optstring(nprofile_opts,profile_opt,string,4)
 call write_inopt(star%iprofile,'iprofile'//trim(c),'0=Sink,'//trim(string(1:40)),iunit)

 if (star%isoftcore <= 0) then
    if (need_inputprofile(star%iprofile)) then
       call write_inopt(star%input_profile,'input_profile'//trim(c),&
            'Path to input profile',iunit)
    else
       call write_inopt(in_solarm(star%Mstar),'Mstar'//trim(c),'mass of star '//trim(c)//' [Msun]',iunit)
       if (need_rstar(star%iprofile)) &
          call write_inopt(in_solarr(star%Rstar),'Rstar'//trim(c),'radius of star'//trim(c)//' [Rsun]',iunit)
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
             call write_inopt(in_solarr(star%rcore),'rcore'//trim(c),'Radius of core softening [Rsun]',iunit)
          endif
          if ((star%isofteningopt == 2) .or. (star%isofteningopt == 3)) then
             call write_inopt(in_solarm(star%mcore),'mcore'//trim(c),&
                              'Mass of point mass stellar core [Msun]',iunit)
          endif
       elseif (star%isoftcore == 2) then
          call write_inopt(in_solarr(star%rcore),'rcore'//trim(c),&
               'Radius of core softening [Rsun]',iunit)
          call write_inopt(in_solarm(star%mcore),'mcore'//trim(c),&
               'Initial guess for mass of sink particle stellar core [Msun]',iunit)
       endif
       call write_inopt(in_solarl(star%lcore),'lcore'//trim(c),&
                              'Luminosity of point mass stellar core [Lsun]',iunit)
    else
       call write_inopt(star%isinkcore,'isinkcore'//trim(c),&
               'Add a sink particle stellar core',iunit)
       if (star%isinkcore) then
          call write_inopt(in_solarm(star%mcore),'mcore'//trim(c),&
               'Mass of sink particle stellar core',iunit)
          call write_inopt(in_solarr(star%hsoft),'hsoft'//trim(c),&
               'Softening length of sink particle stellar core [Rsun]',iunit)
       endif
       call write_inopt(in_solarl(star%lcore),'lcore'//trim(c),&
               'Luminosity of sink core particle [Lsun]',iunit)
    endif
 case (ievrard)
    call write_inopt(star%ui_coef,'ui_coef'//trim(c),&
         'specific internal energy (units of GM/R)',iunit)
 case(0)
    call write_inopt(star%hacc,'hacc'//trim(c),'accretion radius for sink'//trim(c),iunit)
 end select

 if (star%iprofile > 0 .and. (len_trim(c)==0 .or. c(1:1)=='1')) then
    call write_inopt(star%np,'np'//trim(c),'number of particles',iunit)
    !call write_inopt(use_exactN,'use_exactN','find closest particle number to np',iunit)
 endif

end subroutine write_options_star

!-----------------------------------------------------------------------
!+
!  read setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine read_options_star(star,need_iso,ieos,polyk,db,nerr,label)
 use infile_utils,  only:inopts,read_inopt
 use setstar_utils, only:need_inputprofile,need_rstar,nprofile_opts
 use units,         only:umass,udist,unit_luminosity
 use physcon,       only:solarm,solarr,solarl
 type(star_t),              intent(out)   :: star
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(out)   :: need_iso
 integer,                   intent(inout) :: ieos
 real,                      intent(inout) :: polyk
 integer,                   intent(inout) :: nerr
 character(len=*),          intent(in), optional :: label
 character(len=10) :: c
 real :: mcore_msun,rcore_rsun,lcore_lsun,mstar_msun,rstar_rsun,hsoft_rsun

 ! set defaults
 call set_defaults_star(star)

 ! append optional label e.g. '1', '2'
 c = ''
 if (present(label)) c = trim(adjustl(label))
 star%label = trim(c)

 call read_inopt(star%iprofile,'iprofile'//trim(c),db,errcount=nerr,min=0,max=nprofile_opts)
 call set_defaults_given_profile(star%iprofile,star%input_profile,&
                                 need_iso,ieos,star%mstar,polyk)

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
    ! core softening options
    call read_inopt(star%isoftcore,'isoftcore'//trim(c),db,errcount=nerr,min=0)
    if (star%isoftcore==2) star%isofteningopt=3

    if (star%isoftcore <= 0) then ! sink particle core without softening
       call read_inopt(star%isinkcore,'isinkcore'//trim(c),db,errcount=nerr)
       if (star%isinkcore) then
          call read_inopt(mcore_msun,'mcore'//trim(c),db,errcount=nerr,min=0.)
          star%mcore = mcore_msun*real(solarm/umass)
          call read_inopt(hsoft_rsun,'hsoft'//trim(c),db,errcount=nerr,min=0.)
          star%hsoft = hsoft_rsun*real(solarr/udist)
       endif
    else
       star%isinkcore = .true.
       call read_inopt(star%input_profile,'input_profile'//trim(c),db,errcount=nerr)
       call read_inopt(star%outputfilename,'outputfilename'//trim(c),db,errcount=nerr)
       if (star%isoftcore==1) call read_inopt(star%isofteningopt,'isofteningopt'//trim(c),&
                                              db,errcount=nerr,min=0)
       if ((star%isofteningopt==1) .or. (star%isofteningopt==3)) then
          call read_inopt(rcore_rsun,'rcore'//trim(c),db,errcount=nerr,min=0.)
          star%rcore = rcore_rsun*real(solarr/udist)
       endif
       if ((star%isofteningopt==2) .or. (star%isofteningopt==3) &
           .or. (star%isoftcore==2)) then
          call read_inopt(mcore_msun,'mcore'//trim(c),db,errcount=nerr,min=0.)
          star%mcore = mcore_msun*real(solarm/umass)
       endif
       call read_inopt(lcore_lsun,'lcore'//trim(c),db,errcount=nerr,min=0.)
       star%lcore = lcore_lsun*real(solarl/unit_luminosity)
    endif
 case(ievrard)
    call read_inopt(star%ui_coef,'ui_coef'//trim(c),db,errcount=nerr,min=0.)
 case(:0)
    call read_inopt(star%hacc,'hacc'//trim(c),db,errcount=nerr,min=0.)
 end select

 ! star properties
 if (star%isoftcore <= 0) then
    if (need_inputprofile(star%iprofile)) then
       call read_inopt(star%input_profile,'input_profile'//trim(c),db,errcount=nerr)
    else
       call read_inopt(mstar_msun,'Mstar'//trim(c),db,errcount=nerr,min=0.)
       star%mstar = mstar_msun*real(solarm/umass)
       if (need_rstar(star%iprofile)) then
          call read_inopt(rstar_rsun,'Rstar'//trim(c),db,errcount=nerr,min=0.)
          star%rstar = rstar_rsun*real(solarr/udist)
       endif
    endif
 endif

end subroutine read_options_star

end module setstar
