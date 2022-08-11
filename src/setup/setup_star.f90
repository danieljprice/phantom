!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up sphere(s).  There are multiple options,
! as listed in set_sphere.
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - EOSopt            : *EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)*
!   - Mstar             : *mass of star*
!   - Rstar             : *radius of star*
!   - X                 : *hydrogen mass fraction*
!   - dist_unit         : *distance unit (e.g. au)*
!   - gamma             : *Adiabatic index*
!   - hsoft             : *Softening length of sink particle stellar core*
!   - ieos              : *1=isothermal,2=adiabatic,10=MESA,12=idealplusrad*
!   - initialtemp       : *initial temperature of the star*
!   - input_profile     : *Path to input profile*
!   - irecomb           : *Species to include in recombination (0: H2+H+He, 1:H+He, 2:He*
!   - isinkcore         : *Add a sink particle stellar core*
!   - isoftcore         : *0=no core softening, 1=cubic core, 2=constant entropy core*
!   - isofteningopt     : *1=supply rcore, 2=supply mcore, 3=supply both*
!   - mass_unit         : *mass unit (e.g. solarm)*
!   - mcore             : *Mass of sink particle stellar core*
!   - metallicity       : *metallicity*
!   - mu                : *mean molecular weight*
!   - np                : *approx number of particles (in box of size 2R)*
!   - outputfilename    : *Output path for softened MESA profile*
!   - polyk             : *polytropic constant (cs^2 if isothermal)*
!   - rcore             : *Radius of core softening*
!   - relax_star        : *relax star automatically during setup*
!   - ui_coef           : *specific internal energy (units of GM/R)*
!   - use_exactN        : *find closest particle number to np*
!   - use_var_comp      : *Use variable composition (X, Z, mu)*
!   - write_rho_to_file : *write density profile to file*
!
! :Dependencies: centreofmass, dim, eos, eos_gasradrec, eos_piecewise,
!   extern_densprofile, externalforces, infile_utils, io, kernel,
!   mpidomain, mpiutils, options, part, physcon, prompting, relaxstar,
!   setsoftenedcore, setstar, setup_params, timestep, units
!
 use io,             only:fatal,error,master
 use part,           only:gravity,ihsoft
 use physcon,        only:solarm,solarr,km,pi,c,kb_on_mh,radconst
 use options,        only:nfulldump,iexternalforce,calc_erot,use_var_comp
 use timestep,       only:tmax,dtmax
 use eos,            only:ieos
 use externalforces, only:iext_densprofile
 use extern_densprofile, only:nrhotab
 use setsoftenedcore,only:rcore,mcore
 use setstar,        only:iuniform,ipoly,ifromfile,imesa,ikepler,ibpwpoly,ievrard,&
                          nprofile_opts,profile_opt
 implicit none
 !
 ! Input parameters
 !
 integer            :: iprofile,np,EOSopt,isoftcore,isofteningopt
 integer            :: need_iso
 real(kind=8)       :: udist,umass
 real               :: Rstar,Mstar,rhocentre,maxvxyzu,ui_coef,hsoft
 real               :: initialtemp
 logical            :: iexist,input_polyk,isinkcore
 logical            :: use_exactN,relax_star_in_setup,write_rho_to_file
 logical            :: need_inputprofile,need_rstar
 character(len=120) :: input_profile,dens_profile
 character(len=120) :: outputfilename ! outputfilename is the path to the cored profile
 character(len=20)  :: dist_unit,mass_unit
 character(len=30)  :: lattice  ! The lattice type if stretchmap is used

 public             :: setpart
 private

contains

!-----------------------------------------------------------------------
!+
!  Setup routine for stars / spherical collapse calculations
!+
!-----------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use centreofmass,    only:reset_centreofmass
 use units,           only:set_units,select_unit,utime,unit_density
 use kernel,          only:hfact_default
 use extern_densprofile, only:write_rhotab,rhotabfile,read_rhotab_wrapper
 use eos,             only:init_eos,finish_eos,gmw,X_in,Z_in,eos_outputs_mu
 use eos_piecewise,   only:init_eos_piecewise_preset
 use setstar,         only:set_stellar_core,read_star_profile,set_star_density, &
                           set_star_composition,set_star_thermalenergy
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,eos_vars,rad,igas
 use relaxstar,       only:relax_star
 use mpiutils,        only:reduceall_mpi
 use mpidomain,       only:i_belong
 use setup_params,    only:rhozero,npart_total
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer                          :: npts,ierr,ierr_relax
 real                             :: rmin
 real, allocatable                :: r(:),den(:),pres(:),temp(:),en(:),mtab(:),Xfrac(:),Yfrac(:),mu(:)
 logical                          :: setexists
 character(len=120)               :: setupfile,inname
 !
 ! Initialise parameters, including those that will not be included in *.setup
 !
 time         = 0.
 polyk        = 1.0
 gamma        = 5./3.
 hfact        = hfact_default
 maxvxyzu     = size(vxyzu(:,1))
 use_exactN   = .true.
 relax_star_in_setup = .false.
 write_rho_to_file = .false.
 input_polyk  = .false.
 ierr_relax = 0
 !
 ! set default options
 !
 dist_unit   = 'solarr'
 mass_unit   = 'solarm'
 Rstar       = 1.0
 Mstar       = 1.0
 ui_coef     = 1.0
 iprofile    = 1
 EOSopt      = 1
 initialtemp = 1.0e7
 gmw         = 0.5988
 X_in        = 0.74
 Z_in        = 0.02
 isoftcore   = 0
 isinkcore   = .false.
 hsoft          = 0.
 rcore          = 0.
 mcore          = 0.
 isofteningopt  = 1 ! By default, specify rcore
 input_profile  = 'P12_Phantom_Profile.data'
 outputfilename = 'mysoftenedstar.dat'
 dens_profile   = 'density-profile.tab'
 !
 ! defaults needed for error checking
 !
 need_iso = 0       ! -1 = no; 0 = doesn't matter; 1 = yes
 !
 ! determine if an .in file exists
 !
 inname=trim(fileprefix)//'.in'
 inquire(file=inname,exist=iexist)
 if (.not. iexist) then
    tmax      = 100.
    dtmax     = 1.0
    ieos      = 2
 endif
 !
 ! determine if an .setup file exists
 !
 setupfile = trim(fileprefix)//'.setup'
 inquire(file=setupfile,exist=setexists)
 if (setexists) then
    call read_setupfile(setupfile,gamma,polyk,ierr)
    if (ierr /= 0) then
       if (id==master) call write_setupfile(setupfile,gamma,polyk)
       stop 'please rerun phantomsetup with revised .setup file'
    endif
    !--Prompt to get inputs and write to file
 elseif (id==master) then
    print "(a,/)",trim(setupfile)//' not found: using interactive setup'
    call setup_interactive(polyk,gamma,iexist,id,master,ierr)
    call write_setupfile(setupfile,gamma,polyk)
    stop 'please check and edit .setup file and rerun phantomsetup'
 endif

 !
 ! Verify correct pre-processor commands
 !
 if (.not.gravity) then
    iexternalforce = iext_densprofile
    write_rho_to_file = .true.
 endif
 !
 ! set lattice, use closepacked unless relaxation is done automatically
 !
 lattice = 'closepacked'
 if (relax_star_in_setup) lattice='random'

 if (maxvxyzu > 3  .and. need_iso == 1) call fatal('setup','require ISOTHERMAL=yes')
 if (maxvxyzu < 4  .and. need_iso ==-1) call fatal('setup','require ISOTHERMAL=no')
 !
 ! set units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !call set_units(mass=umass, c=1.d0, G=1.d0) ! uncomment if want geometric units
 !
 ! set up particles
 !
 npartoftype(:) = 0
 npart          = 0
 vxyzu          = 0.0
 !
 ! initialise the equation of state
 !
 if (ieos==9) call init_eos_piecewise_preset(EOSopt)
 call init_eos(ieos,ierr)

 !
 ! get the desired tables of density, pressure, temperature and composition
 ! as a function of radius / mass fraction
 !
 call read_star_profile(iprofile,ieos,input_profile,gamma,polyk,ui_coef,r,den,pres,temp,en,mtab,&
                        Xfrac,Yfrac,mu,npts,rmin,Rstar,Mstar,rhocentre,&
                        isoftcore,isofteningopt,rcore,hsoft,outputfilename)
 !
 ! set up particles to represent the desired stellar profile
 !
 call set_star_density(lattice,id,master,rmin,Rstar,Mstar,hfact,&
                       npts,den,r,npart,npartoftype,massoftype,xyzh,use_exactN,np,i_belong)
 !
 ! add sink particle stellar core
 !
 if (isinkcore) call set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,ihsoft,mcore,hsoft)
 !
 ! Write the desired profile to file (do this before relaxation)
 !
 if (write_rho_to_file) call write_rhotab(dens_profile,r,den,npts,polyk,gamma,rhocentre,ierr)
 !
 ! relax the density profile to achieve nice hydrostatic equilibrium
 !
 if (relax_star_in_setup) then
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
    call set_star_composition(use_var_comp,eos_outputs_mu(ieos),npart,xyzh,Xfrac,Yfrac,mu,mtab,Mstar,eos_vars)
 endif
 !
 ! set the internal energy and temperature
 !
 if (maxvxyzu==4) call set_star_thermalenergy(ieos,den,pres,r,npart,xyzh,vxyzu,rad,&
                                              eos_vars,relax_star_in_setup,use_var_comp,initialtemp)

 call finish_eos(ieos,ierr)
 !
 ! Reset centre of mass (again)
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 !
 ! Print summary to screen
 !
 write(*,"(70('='))")
 if (ieos /= 9) then
    write(*,'(1x,a,f12.5)')       'gamma               = ', gamma
 endif
 if (maxvxyzu <= 4) then
    write(*,'(1x,a,f12.5)')       'polyk               = ', polyk
    write(*,'(1x,a,f12.6,a)')     'specific int. energ = ', polyk*Rstar/Mstar,' GM/R'
 endif
 call write_mass('particle mass       = ',massoftype(igas),umass)
 call write_dist('Radius              = ',Rstar,udist)
 call write_mass('Mass                = ',Mstar,umass)
 if (iprofile==ipoly) then
    write(*,'(1x,a,es12.5,a)') 'rho_central         = ', rhocentre*unit_density,' g/cm^3'
 endif
 write(*,'(1x,a,i12)')         'N                   = ', npart_total
 write(*,'(1x,a,2(es12.5,a))')    'rho_mean            = ', rhozero*unit_density,  ' g/cm^3 = '&
                                                       , rhozero,               ' code units'
 write(*,'(1x,a,es12.5,a)')       'free fall time      = ', sqrt(3.*pi/(32.*rhozero))*utime,' s'
 if ( (iprofile==iuniform .and. .not.gravity) .or. iprofile==ifromfile) then
    write(*,'(a)') 'WARNING! This setup may not be stable'
 endif
 write(*,"(70('='))")
 if (ierr_relax /= 0) write(*,"(/,a,/)") ' WARNING: ERRORS DURING RELAXATION, SEE ABOVE!!'

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Ask questions of the user to determine which setup to use
!+
!-----------------------------------------------------------------------
subroutine setup_interactive(polyk,gamma,iexist,id,master,ierr)
 use prompting,     only:prompt
 use units,         only:select_unit
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 real, intent(out)    :: polyk,gamma
 logical, intent(in)  :: iexist
 integer, intent(in)  :: id,master
 integer, intent(out) :: ierr
 integer :: i
 logical :: need_rstar

 ierr = 0
 ! Select sphere & set default values
 do i = 1, nprofile_opts
    write(*,"(i2,')',1x,a)") i, profile_opt(i)
 enddo

 call prompt('Enter which density profile to use',iprofile,1,nprofile_opts)
 !
 ! set default file output parameters
 !
 if (id==master) write(*,"('Setting up ',a)") trim(profile_opt(iprofile))
 call set_defaults_given_profile(iprofile,iexist,need_inputprofile,need_rstar,polyk)

 ! units
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

 ! resolution
 np    = 100000 ! default number of particles
 call prompt('Enter the approximate number of particles in the sphere ',np,0)

 ! equation of state
 call prompt('Enter the desired EoS for setup',ieos)
 select case(ieos)
 case(15) ! Helmholtz
    call prompt('Enter temperature',initialtemp,1.0e3,1.0e11)
 case(9)
    write(*,'(a)') 'EOS options: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)'
    call prompt('Enter equation of state type',EOSopt,1,4)
 case(2)
    call prompt('Enter gamma (adiabatic index)',gamma,1.,7.)
    if (input_polyk) call prompt('Enter polytropic constant',polyk,0.)
 case(1)
    call prompt('Enter polytropic constant (cs^2 if isothermal)',polyk,0.)
 case(20)
    call prompt('Enter irecomb (0: H2+H+He, 1:H+He, 2:He)',irecomb,0)
 end select
 if (iprofile==ievrard) then
    call prompt('Enter the specific internal energy (units of GM/R) ',ui_coef,0.)
 endif

 ! star properties
 if (need_inputprofile) then
    call prompt('Enter file name containing input profile ', input_profile)
 else
    call prompt('Enter the mass of the star (code units) ',  Mstar,0.)
    if (need_rstar) call prompt('Enter the radius of the star (code units) ',Rstar,0.)
 endif

 if (iprofile==imesa) then
    call prompt('Use variable composition?',use_var_comp)

    print*,'Soften the core density profile and add a sink particle core?'
    print "(3(/,a))",'0: Do not soften profile', &
                     '1: Use cubic softened density profile', &
                     '2: Use constant entropy softened profile'
    call prompt('Select option above : ',isoftcore,0,2)

    select case (isoftcore)
    case(0)
       call prompt('Add a sink particle stellar core?',isinkcore)
       if (isinkcore) then
          call prompt('Enter mass of the created sink particle core',mcore,0.)
          call prompt('Enter softening length of the sink particle core',hsoft,0.)
       endif
    case(1)
       isinkcore = .true. ! Create sink particle core automatically
       print*,'Options for core softening:'
       print "(5(/,a))",'1. Specify radius of density softening', &
                        '2. Specify mass of sink particle core (not recommended)', &
                        '3. Specify both radius of density softening length and mass', &
                        '   of sink particle core (if you do not know what you are', &
                        '   doing, you will obtain a poorly softened profile)'
       call prompt('Select option above : ',isofteningopt,1,3)

       select case (isofteningopt)
       case(1)
          call prompt('Enter core radius',rcore,0.)
       case(2)
          call prompt('Enter mass of the created sink particle core',mcore,0.)
       case(3)
          call prompt('Enter mass of the created sink particle core',mcore,0.)
          call prompt('Enter core radius',rcore,0.)
       end select

       call prompt('Enter output file name of cored stellar profile:',outputfilename)

    case(2)
       isinkcore = .true. ! Create sink particle core automatically
       print*,'Specify core radius and initial guess for mass of sink particle core'
       call prompt('Enter core radius in Rsun : ',rcore,0.)
       call prompt('Enter guess for core mass in Msun : ',mcore,0.)
       call prompt('Enter output file name of cored stellar profile:',outputfilename)
    end select

    if ((.not. use_var_comp) .and. (isoftcore<=0)) then
       if ( (ieos==12) .or. (ieos==2) ) call prompt('Enter mean molecular weight',gmw,0.)
       if ( (ieos==10) .or. (ieos==20) ) then
          call prompt('Enter hydrogen mass fraction (X)',X_in,0.,1.)
          call prompt('Enter metals mass fraction (Z)',Z_in,0.,1.)
       endif
    endif
 endif
 call prompt('Relax star automatically during setup?',relax_star_in_setup)

end subroutine setup_interactive

!-----------------------------------------------------------------------
!+
!  Set default options associated with particular star setups
!  This routine should not do ANY prompting
!+
!-----------------------------------------------------------------------
subroutine set_defaults_given_profile(iprofile,iexist,need_inputprofile,need_rstar,polyk)
 integer, intent(in)  :: iprofile
 logical, intent(in)  :: iexist
 logical, intent(out) :: need_inputprofile,need_rstar
 real,    intent(out) :: polyk

 need_rstar = .true.
 need_inputprofile = .false.
 select case(iprofile)
 case(iuniform)
    input_polyk  = .true.
 case(ipoly)
    input_polyk  = .false.
 case(ifromfile)
    ! Read the density profile from file (e.g. for neutron star)
    !  Original Author: Mark Bennett
    input_profile = 'ns-rdensity.tab'
    need_inputprofile = .true.
 case(imesa)
    ! sets up a star from a 1D MESA code output
    !  Original Author: Roberto Iaconi
    input_profile = 'P12_Phantom_Profile.data'
    need_inputprofile = .true.
 case(ikepler)
    ! sets up a star from a 1D KEPLER code output
    !  Original Author: Nicole Rodrigues
    input_profile = 'kepler_MS.data'
    need_inputprofile = .true.
 case(ibpwpoly)
    !  piecewise polytrope
    !  Original Author: Madeline Marshall & Bernard Field
    !  Supervisors: James Wurster & Paul Lasky
    ieos         = 9
    dist_unit    = 'km'
    Mstar        = 1.35
    polyk        = 144.
    calc_erot    = .true.
    need_rstar   = .false.
    input_polyk  = .true.
 case(ievrard)
    ! Evrard Collapse
    if (.not. iexist) then
       tmax      = 3.0
       dtmax     = 0.1
       nfulldump = 1
    endif
    ui_coef     = 0.05
    need_iso    = -1
    input_polyk  = .false.
 end select

end subroutine set_defaults_given_profile

!-----------------------------------------------------------------------
!+
!  Subroutines to write summary to screen
!+
!-----------------------------------------------------------------------
subroutine write_dist(item_in,dist_in,udist)
 real,             intent(in) :: dist_in
 real(kind=8),     intent(in) :: udist
 character(len=*), intent(in) :: item_in

 if ( abs(1.0-solarr/udist) < 1.0d-4) then
    write(*,'(1x,2(a,Es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' R_sun'
 elseif ( abs(1.0-km/udist) < 1.0d-4) then
    write(*,'(1x,2(a,Es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' km'
 else
    write(*,'(1x,a,Es12.5,a)')    item_in, dist_in*udist,' cm'
 endif

end subroutine write_dist

subroutine write_mass(item_in,mass_in,umass)
 real,             intent(in) :: mass_in
 real(kind=8),     intent(in) :: umass
 character(len=*), intent(in) :: item_in

 if ( abs(1.0-solarm/umass) < 1.0d-4) then
    write(*,'(1x,2(a,Es12.5),a)') item_in, mass_in*umass,' g      = ',mass_in,' M_sun'
 else
    write(*,'(1x,a,Es12.5,a)')    item_in, mass_in*umass,' g'
 endif

end subroutine write_mass
!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename,gamma,polyk)
 use infile_utils,  only:write_inopt,get_optstring
 use dim,           only:tagline
 use relaxstar,     only:write_options_relax
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 real,             intent(in) :: gamma,polyk
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20
 character(len=120)           :: string

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(a)") '# input file for Phantom star setup'

 call get_optstring(nprofile_opts,profile_opt,string,4)
 call write_inopt(iprofile,'iprofile',trim(string),iunit)

 write(iunit,"(/,a)") '# units'
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','approx number of particles (in box of size 2R)',iunit)
 call write_inopt(use_exactN,'use_exactN','find closest particle number to np',iunit)

 write(iunit,"(/,a)") '# equation of state'
 call write_inopt(use_var_comp,'use_var_comp','Use variable composition (X, Z, mu)',iunit)
 call write_inopt(ieos,'ieos','1=isothermal,2=adiabatic,10=MESA,12=idealplusrad',iunit)
 select case(ieos)
 case(15) ! Helmholtz
    call write_inopt(initialtemp,'initialtemp','initial temperature of the star',iunit)
 case(9)
    write(iunit,"(/,a)") '# Piecewise Polytrope default options'
    call write_inopt(EOSopt,'EOSopt','EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)',iunit)
 case(2)
    call write_inopt(gamma,'gamma','Adiabatic index',iunit)
    if (input_polyk) call write_inopt(polyk,'polyk','polytropic constant (cs^2 if isothermal)',iunit)
    if ((isoftcore<=0) .and. (.not. use_var_comp)) call write_inopt(gmw,'mu','mean molecular weight',iunit)
 case(1)
    if (input_polyk) call write_inopt(polyk,'polyk','polytropic constant (cs^2 if isothermal)',iunit)
 case(10,20)
    if (ieos==20) call write_inopt(irecomb,'irecomb','Species to include in recombination (0: H2+H+He, 1:H+He, 2:He',iunit)
    if ( (.not. use_var_comp) .and. (isoftcore<=0) ) then
       call write_inopt(X_in,'X','hydrogen mass fraction',iunit)
       call write_inopt(Z_in,'Z','metallicity',iunit)
    endif
 case(12)
    call write_inopt(gamma,'gamma','Adiabatic index',iunit)
    if (.not. use_var_comp) call write_inopt(gmw,'mu','mean molecular weight',iunit)
 end select
 if (iprofile==ievrard) then
    call write_inopt(ui_coef,'ui_coef','specific internal energy (units of GM/R)',iunit)
 endif

 if (isoftcore <= 0) then
    write(iunit,"(/,a)") '# star properties'
    if (need_inputprofile) then
       call write_inopt(input_profile,'input_profile','Path to input profile',iunit)
    else
       call write_inopt(Rstar,'Rstar','radius of star',iunit)
       call write_inopt(Mstar,'Mstar','mass of star',iunit)
    endif
 endif

 if (iprofile==imesa) then
    write(iunit,"(/,a)") '# core softening and sink stellar core options'
    call write_inopt(isoftcore,'isoftcore','0=no core softening, 1=cubic core, 2=constant entropy core',iunit)
    if (isoftcore > 0) then
       call write_inopt(input_profile,'input_profile','Path to input profile',iunit)
       call write_inopt(outputfilename,'outputfilename','Output path for softened MESA profile',iunit)
       if (isoftcore == 1) then
          call write_inopt(isofteningopt,'isofteningopt','1=supply rcore, 2=supply mcore, 3=supply both',iunit)
          if ((isofteningopt == 1) .or. (isofteningopt == 3)) then
             call write_inopt(rcore,'rcore','Radius of core softening',iunit)
          endif
          if ((isofteningopt == 2) .or. (isofteningopt == 3)) then
             call write_inopt(mcore,'mcore','Mass of sink particle stellar core',iunit)
          endif
       elseif (isoftcore == 2) then
          call write_inopt(rcore,'rcore','Radius of core softening',iunit)
          call write_inopt(mcore,'mcore','Initial guess for mass of sink particle stellar core',iunit)
       endif
    else
       call write_inopt(isinkcore,'isinkcore','Add a sink particle stellar core',iunit)
       if (isinkcore) then
          call write_inopt(mcore,'mcore','Mass of sink particle stellar core',iunit)
          call write_inopt(hsoft,'hsoft','Softening length of sink particle stellar core',iunit)
       endif
    endif
 endif

 write(iunit,"(/,a)") '# relaxation options'
 call write_inopt(relax_star_in_setup,'relax_star','relax star automatically during setup',iunit)
 if (relax_star_in_setup) call write_options_relax(iunit)

 call write_inopt(write_rho_to_file,'write_rho_to_file','write density profile to file',iunit)

 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,polyk,ierr)
 use infile_utils,  only:open_db_from_file,inopts,close_db,read_inopt
 use io,            only:error
 use units,         only:select_unit
 use relaxstar,     only:read_options_relax
 use eos,           only:X_in,Z_in,gmw
 use eos_gasradrec, only:irecomb
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma,polyk
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)

 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'setup_star: Reading setup options from ',trim(filename)

 nerr = 0
 call read_inopt(iprofile,'iprofile',db,errcount=nerr)
 call set_defaults_given_profile(iprofile,iexist,need_inputprofile,need_rstar,polyk)
 if (need_inputprofile) call read_inopt(input_profile,'input_profile',db,errcount=nerr)

 ! units
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)

 ! resolution
 call read_inopt(np,'np',db,errcount=nerr)
 call read_inopt(use_exactN,'use_exactN',db,errcount=nerr)

 ! core softening
 if (iprofile==imesa) then
    call read_inopt(use_var_comp,'use_var_comp',db,errcount=nerr)
    call read_inopt(isoftcore,'isoftcore',db,errcount=nerr)
    if (isoftcore==2) isofteningopt=3
 endif

 ! equation of state
 call read_inopt(ieos,'ieos',db,errcount=nerr)
 select case(ieos)
 case(15) ! Helmholtz
    call read_inopt(initialtemp,'initialtemp',db,errcount=nerr)
 case(9)
    call read_inopt(EOSopt,'EOSopt',db,errcount=nerr)
 case(2)
    call read_inopt(gamma,'gamma',db,errcount=nerr)
    if ( (.not. use_var_comp) .and. (isoftcore <= 0)) call read_inopt(gmw,'mu',db,errcount=nerr)
    if (input_polyk) call read_inopt(polyk,'polyk',db,errcount=nerr)
 case(1)
    if (input_polyk) call read_inopt(polyk,'polyk',db,errcount=nerr)
 case(10,20)
    if (ieos==20) call read_inopt(irecomb,'irecomb',db,errcount=nerr)
    ! if softening stellar core, composition is automatically determined at R/2
    if ( (.not. use_var_comp) .and. (isoftcore <= 0)) then
       call read_inopt(X_in,'X',db,errcount=nerr)
       call read_inopt(Z_in,'Z',db,errcount=nerr)
    endif
 case(12)
    ! if softening stellar core, mu is automatically determined at R/2
    call read_inopt(gamma,'gamma',db,errcount=nerr)
    if ( (.not. use_var_comp) .and. (isoftcore <= 0)) call read_inopt(gmw,'mu',db,errcount=nerr)
 end select
 if (iprofile==ievrard) call read_inopt(ui_coef,'ui_coef',db,errcount=nerr)

 ! core softening options
 if (iprofile==imesa) then
    if (isoftcore <= 0) then ! sink particle core without softening
       call read_inopt(isinkcore,'isinkcore',db,errcount=nerr)
       if (isinkcore) then
          call read_inopt(mcore,'mcore',db,errcount=nerr)
          call read_inopt(hsoft,'hsoft',db,errcount=nerr)
       endif
    else
       isinkcore = .true.
       call read_inopt(input_profile,'input_profile',db,errcount=nerr)
       call read_inopt(outputfilename,'outputfilename',db,errcount=nerr)
       if (isoftcore==1) call read_inopt(isofteningopt,'isofteningopt',db,errcount=nerr)
       if ((isofteningopt==1) .or. (isofteningopt==3)) call read_inopt(rcore,'rcore',db,errcount=nerr)
       if ((isofteningopt==2) .or. (isofteningopt==3) .or. (isoftcore==2)) call read_inopt(mcore,'mcore',db,errcount=nerr)
    endif
 endif

 ! star properties
 if (isoftcore <= 0) then
    if (need_inputprofile) then
       call read_inopt(input_profile,'input_profile',db,errcount=nerr)
    else
       call read_inopt(Mstar,'Mstar',db,errcount=nerr)
       if (need_rstar) call read_inopt(Rstar,'Rstar',db,errcount=nerr)
    endif
 endif

 ! relax star options
 call read_inopt(relax_star_in_setup,'relax_star',db,errcount=nerr)
 if (relax_star_in_setup) call read_options_relax(db,nerr)
 if (nerr /= 0) ierr = ierr + 1

 ! option to write density profile to file
 call read_inopt(write_rho_to_file,'write_rho_to_file',db)

 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_star','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_star','length unit not recognised')
    ierr = ierr + 1
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'setup_star: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile

end module setup
