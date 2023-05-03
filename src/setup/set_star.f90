!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: eos, eos_piecewise, extern_densprofile, io, part, physcon,
!   radiation_utils, rho_profile, setsoftenedcore, setup_params, sortutils,
!   spherical, table_utils, unifdis, units
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
  real :: hsoft
  character(len=120) :: input_profile,dens_profile
  character(len=120) :: outputfilename ! outputfilename is the path to the cored profile
 end type star_t

 public :: star_t
 public :: set_defaults_star,set_star
 public :: write_star_options,read_star_options,set_star_interactive
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
 type(star_t), intent(out) :: star

 star%iprofile    = 1
 star%rstar       = 1.0
 star%mstar       = 1.0
 star%ui_coef     = 0.05
 star%initialtemp = 1.0e7
 star%isoftcore   = 0
 star%isinkcore   = .false.
 star%hsoft          = 0.
 star%rcore          = 0.
 star%mcore          = 0.
 star%isofteningopt  = 1 ! By default, specify rcore
 star%input_profile  = 'P12_Phantom_Profile.data'
 star%outputfilename = 'mysoftenedstar.dat'
 star%dens_profile   = 'density-profile.tab'

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
                    rhozero,npart_total,mask,ierr)
 use centreofmass,       only:reset_centreofmass
 use dim,                only:do_radiation,gr,gravity,maxvxyzu
 use io,                 only:fatal,error,warning
 use eos,                only:eos_outputs_mu
 use setstar_utils,      only:set_stellar_core,read_star_profile,set_star_density, &
                              set_star_composition,set_star_thermalenergy,&
                              write_kepler_comp
 use radiation_utils,    only:set_radiation_and_gas_temperature_equal
 use relaxstar,          only:relax_star
 use part,               only:ihsoft,igas,imu
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
 real,         intent(out)    :: massoftype(:)
 real,         intent(in)     :: hfact
 logical,      intent(in)     :: relax,use_var_comp,write_rho_to_file
 integer,      intent(in)     :: ieos
 real,         intent(inout)  :: polyk,gamma
 real,         intent(in)     :: X_in,Z_in
 real,         intent(out)    :: rhozero
 integer(kind=8), intent(out) :: npart_total
 integer,      intent(out)    :: ierr
 procedure(mask_prototype)    :: mask
 integer                        :: npts,ierr_relax
 integer                        :: ncols_compo
 real, allocatable              :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),mu(:)
 real, allocatable              :: composition(:,:)
 real                           :: rmin,rhocentre
 character(len=20), allocatable :: comp_label(:)
 character(len=30)              :: lattice  ! The lattice type if stretchmap is used
 logical                        :: use_exactN,composition_exists

 use_exactN = .true.
 composition_exists = .false.
 ierr_relax = 0
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
 call set_star_density(lattice,id,master,rmin,star%rstar,star%mstar,hfact,&
                       npts,den,r,npart,npartoftype,massoftype,xyzh,use_exactN,&
                       star%np,rhozero,npart_total,mask)
 !
 ! add sink particle stellar core
 !
 if (star%isinkcore) call set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,star%mcore,star%hsoft,ierr)
 if (ierr==1) call fatal('set_stellar_core','mcore <= 0')
 if (ierr==2) call fatal('set_stellar_core','hsoft <= 0')
 !
 ! Write the desired profile to file (do this before relaxation)
 !
 if (write_rho_to_file) call write_rhotab(star%dens_profile,r,den,npts,polyk,gamma,rhocentre,ierr)
 !
 ! relax the density profile to achieve nice hydrostatic equilibrium
 !
 if (relax) then
    if (reduceall_mpi('+',npart)==npart) then
       call relax_star(npts,den,pres,r,npart,xyzh,use_var_comp,Xfrac,Yfrac,mu,ierr_relax)
    else
       call error('setup_star','cannot run relaxation with MPI setup, please run setup on ONE MPI thread')
    endif
 endif
 !
 ! set composition (X,Z,mu, if using variable composition) of each particle by interpolating from table
 !
 if (use_var_comp .or. eos_outputs_mu(ieos)) then
    call set_star_composition(use_var_comp,eos_outputs_mu(ieos),npart,&
                              xyzh,Xfrac,Yfrac,mu,mtab,star%mstar,eos_vars)
 endif
 !
 ! Write composition file called kepler.comp contaning composition of each particle after interpolation
 !
 if (star%iprofile==iKepler) call write_kepler_comp(composition,comp_label,ncols_compo,r,&
                                                    xyzh,npart,npts,composition_exists)
 !
 ! set the internal energy and temperature
 !
 if (maxvxyzu==4) call set_star_thermalenergy(ieos,den,pres,r,npts,npart,xyzh,vxyzu,rad,&
                                              eos_vars,relax,use_var_comp,star%initialtemp)

 if (do_radiation) then
    if (eos_outputs_mu(ieos)) then
       call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad,eos_vars(imu,:))
    else
       call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)
    endif
 endif
 !
 ! Reset centre of mass (again)
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (gr) then
    xyzh(1,:)=xyzh(1,:)+10.
    xyzh(2,:)=xyzh(2,:)+10.
 endif
 !
 ! Print summary to screen
 !
 if (id==master) then
    write(*,"(70('='))")
    if (ieos /= 9) then
       write(*,'(1x,a,f12.5)')       'gamma               = ', gamma
    endif
    if (maxvxyzu <= 4) then
       write(*,'(1x,a,f12.5)')       'polyk               = ', polyk
       write(*,'(1x,a,f12.6,a)')     'specific int. energ = ', polyk*star%rstar/star%mstar,' GM/R'
    endif
    call write_mass('particle mass       = ',massoftype(igas),umass)
    call write_dist('Radius              = ',star%rstar,udist)
    call write_mass('Mass                = ',star%mstar,umass)
    if (star%iprofile==ipoly) then
       write(*,'(1x,a,es12.5,a)') 'rho_central         = ', rhocentre*unit_density,' g/cm^3'
    endif
    write(*,'(1x,a,i12)')         'N                   = ', npart_total
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
!  print a distance in both code units and physical units
!+
!-----------------------------------------------------------------------
subroutine write_dist(item_in,dist_in,udist)
 use physcon, only:solarr,km
 real,             intent(in) :: dist_in
 real(kind=8),     intent(in) :: udist
 character(len=*), intent(in) :: item_in

 if ( abs(1.0-solarr/udist) < 1.0d-4) then
    write(*,'(1x,2(a,es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' R_sun'
 elseif ( abs(1.0-km/udist) < 1.0d-4) then
    write(*,'(1x,2(a,es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' km'
 else
    write(*,'(1x,a,es12.5,a)')    item_in, dist_in*udist,' cm'
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

 if ( abs(1.0-solarm/umass) < 1.0d-4) then
    write(*,'(1x,2(a,es12.5),a)') item_in, mass_in*umass,' g      = ',mass_in,' M_sun'
 else
    write(*,'(1x,a,es12.5,a)')    item_in, mass_in*umass,' g'
 endif

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
!  write setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine write_star_options(star,iunit)
 use infile_utils,  only:write_inopt,get_optstring
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 type(star_t), intent(in) :: star
 integer,      intent(in) :: iunit
 character(len=120) :: string

 write(iunit,"(/,a)") '# options for star'
 call get_optstring(nprofile_opts,profile_opt,string,4)
 call write_inopt(star%iprofile,'iprofile',trim(string),iunit)

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(star%np,'np','number of particles',iunit)
 !call write_inopt(use_exactN,'use_exactN','find closest particle number to np',iunit)

 if (star%isoftcore <= 0) then
    write(iunit,"(/,a)") '# star properties'
    if (need_inputprofile(star%iprofile)) then
       call write_inopt(star%input_profile,'input_profile','Path to input profile',iunit)
    else
       call write_inopt(star%Mstar,'Mstar','mass of star',iunit)
       if (need_rstar(star%iprofile)) call write_inopt(star%Rstar,'Rstar','radius of star',iunit)
    endif
 endif

 if (star%iprofile==imesa) then
    write(iunit,"(/,a)") '# core softening and sink stellar core options'
    call write_inopt(star%isoftcore,'isoftcore','0=no core softening, 1=cubic core, 2=constant entropy core',iunit)
    if (star%isoftcore > 0) then
       call write_inopt(star%input_profile,'input_profile','Path to input profile',iunit)
       call write_inopt(star%outputfilename,'outputfilename','Output path for softened MESA profile',iunit)
       if (star%isoftcore == 1) then
          call write_inopt(star%isofteningopt,'isofteningopt','1=supply rcore, 2=supply mcore, 3=supply both',iunit)
          if ((star%isofteningopt == 1) .or. (star%isofteningopt == 3)) then
             call write_inopt(star%rcore,'rcore','Radius of core softening',iunit)
          endif
          if ((star%isofteningopt == 2) .or. (star%isofteningopt == 3)) then
             call write_inopt(star%mcore,'mcore','Mass of sink particle stellar core',iunit)
          endif
       elseif (star%isoftcore == 2) then
          call write_inopt(star%rcore,'rcore','Radius of core softening',iunit)
          call write_inopt(star%mcore,'mcore','Initial guess for mass of sink particle stellar core',iunit)
       endif
    else
       call write_inopt(star%isinkcore,'isinkcore','Add a sink particle stellar core',iunit)
       if (star%isinkcore) then
          call write_inopt(star%mcore,'mcore','Mass of sink particle stellar core',iunit)
          call write_inopt(star%hsoft,'hsoft','Softening length of sink particle stellar core',iunit)
       endif
    endif
 endif

 if (star%iprofile==ievrard) then
    call write_inopt(star%ui_coef,'ui_coef','specific internal energy (units of GM/R)',iunit)
 endif

end subroutine write_star_options

!-----------------------------------------------------------------------
!+
!  interactive prompting for setting up a star
!+
!-----------------------------------------------------------------------
subroutine set_star_interactive(id,master,star,need_iso,use_var_comp,ieos,polyk)
 use prompting,     only:prompt
 use setstar_utils, only:nprofile_opts,profile_opt,need_inputprofile,need_rstar
 type(star_t), intent(out)   :: star
 integer,      intent(in)    :: id,master
 logical,      intent(out)   :: use_var_comp
 integer,      intent(out)   :: need_iso
 integer,      intent(inout) :: ieos
 real,         intent(inout) :: polyk
 integer :: i

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
 star%np = 100000 ! default number of particles
 call prompt('Enter the approximate number of particles in the sphere ',star%np,0)

 ! star properties
 if (need_inputprofile(star%iprofile)) then
    call prompt('Enter file name containing input profile',star%input_profile)
 else
    call prompt('Enter the mass of the star (code units)',star%mstar,0.)
    if (need_rstar(star%iprofile)) call prompt('Enter the radius of the star (code units)',star%rstar,0.)
 endif

 if (star%iprofile==imesa) then
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
          call prompt('Enter mass of the created sink particle core',star%mcore,0.)
          call prompt('Enter softening length of the sink particle core',star%hsoft,0.)
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
          call prompt('Enter core radius',star%rcore,0.)
       case(2)
          call prompt('Enter mass of the created sink particle core',star%mcore,0.)
       case(3)
          call prompt('Enter mass of the created sink particle core',star%mcore,0.)
          call prompt('Enter core radius',star%rcore,0.)
       end select

       call prompt('Enter output file name of cored stellar profile:',star%outputfilename)

    case(2)
       star%isinkcore = .true. ! Create sink particle core automatically
       print*,'Specify core radius and initial guess for mass of sink particle core'
       call prompt('Enter core radius in Rsun : ',star%rcore,0.)
       call prompt('Enter guess for core mass in Msun : ',star%mcore,0.)
       call prompt('Enter output file name of cored stellar profile:',star%outputfilename)
    end select

 endif

 if (star%iprofile==ievrard) then
    call prompt('Enter the specific internal energy (units of GM/R) ',star%ui_coef,0.)
 endif

end subroutine set_star_interactive

!-----------------------------------------------------------------------
!+
!  read setupfile options needed for a star
!+
!-----------------------------------------------------------------------
subroutine read_star_options(star,need_iso,ieos,polyk,db,nerr)
 use infile_utils,  only:inopts,read_inopt
 use setstar_utils, only:need_inputprofile,need_rstar
 type(star_t),              intent(out)   :: star
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(out)   :: need_iso
 integer,                   intent(inout) :: ieos
 real,                      intent(inout) :: polyk
 integer,                   intent(inout) :: nerr

 call read_inopt(star%iprofile,'iprofile',db,errcount=nerr)
 call set_defaults_given_profile(star%iprofile,star%input_profile,&
                                 need_iso,ieos,star%mstar,polyk)

 if (need_inputprofile(star%iprofile)) then
    call read_inopt(star%input_profile,'input_profile',db,errcount=nerr)
 endif

 ! resolution
 call read_inopt(star%np,'np',db,errcount=nerr)
 !call read_inopt(use_exactN,'use_exactN',db,errcount=nerr)

 ! core softening options
 if (star%iprofile==imesa) then
    call read_inopt(star%isoftcore,'isoftcore',db,errcount=nerr)
    if (star%isoftcore==2) star%isofteningopt=3

    if (star%isoftcore <= 0) then ! sink particle core without softening
       call read_inopt(star%isinkcore,'isinkcore',db,errcount=nerr)
       if (star%isinkcore) then
          call read_inopt(star%mcore,'mcore',db,errcount=nerr)
          call read_inopt(star%hsoft,'hsoft',db,errcount=nerr)
       endif
    else
       star%isinkcore = .true.
       call read_inopt(star%input_profile,'input_profile',db,errcount=nerr)
       call read_inopt(star%outputfilename,'outputfilename',db,errcount=nerr)
       if (star%isoftcore==1) call read_inopt(star%isofteningopt,'isofteningopt',db,errcount=nerr)
       if ((star%isofteningopt==1) .or. (star%isofteningopt==3)) call read_inopt(star%rcore,'rcore',db,errcount=nerr)
       if ((star%isofteningopt==2) .or. (star%isofteningopt==3) &
          .or. (star%isoftcore==2)) call read_inopt(star%mcore,'mcore',db,errcount=nerr)
    endif
 endif

 ! star properties
 if (star%isoftcore <= 0) then
    if (need_inputprofile(star%iprofile)) then
       call read_inopt(star%input_profile,'input_profile',db,errcount=nerr)
    else
       call read_inopt(star%mstar,'Mstar',db,errcount=nerr)
       if (need_rstar(star%iprofile)) call read_inopt(star%rstar,'Rstar',db,errcount=nerr)
    endif
 endif

end subroutine read_star_options

end module setstar
