!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!    This module sets up sphere(s).  There are multiple options, incuding
!    1) uniform unit sphere
!    2) single polytrope
!    3) binary polytrope (decommissioned)
!    4) neutron star from file
!    5) red giant (Macquarie)
!    6) neutron star using a piecewise polytrope EOS
!    7) Evrard sphere
!    8) KEPLER star from file
!    9) Helmholtz Equation of state
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    EOSopt             -- EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)
!    Mstar_1            -- mass of star
!    Rstar_1            -- radius of star
!    densityfile        -- File containing data for stellar profile
!    dist_unit          -- distance unit (e.g. au)
!    gamma              -- adiabatic index
!    hdens              -- Radius of core softening
!    hsoft              -- Softening length of sink particle stellar core
!    ieos               -- 1=isothermal,2=adiabatic,10=MESA,12=idealplusrad
!    initialtemp        -- initial temperature of the star
!    isinkcore          -- Add a sink particle stellar core
!    isoftcore          -- Soften the core of an input MESA profile
!    isofteningopt      -- 1=supply hsoft, 2=supply mcore, 3=supply both
!    mass_unit          -- mass unit (e.g. solarm)
!    mcore              -- Mass of sink particle stellar core
!    np                 -- approx number of particles (in box of size 2R)
!    outputfilename     -- Output path for softened MESA profile
!    polyk              -- sound speed .or. constant in EOS
!    relax_star         -- relax star automatically during setup
!    ui_coef            -- specific internal energy (units of GM/R)
!    unsoftened_profile -- Path to MESA profile for softening
!
!  DEPENDENCIES: centreofmass, dim, eos, eos_idealplusrad,
!    extern_densprofile, externalforces, infile_utils, io, kernel, options,
!    part, physcon, prompting, relaxstar, rho_profile, setsoftenedcore,
!    setstellarcore, setup_params, spherical, table_utils, timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use io,             only:fatal,error,master
 use part,           only:gravity
 use physcon,        only:solarm,solarr,km,pi,c,kb_on_mh,radconst
 use options,        only:nfulldump,iexternalforce,calc_erot
 use timestep,       only:tmax,dtmax
 use eos,            only:ieos, p1pwpcgs,gamma1pwp,gamma2pwp,gamma3pwp
 use externalforces, only:iext_densprofile
 use extern_densprofile, only:nrhotab

 implicit none
 !
 ! Input parameters
 !
 integer            :: iprofile,np,EOSopt,isofteningopt
 integer            :: nstar
 integer            :: need_iso, need_temp
 real(kind=8)       :: udist,umass
 real               :: Rstar,Mstar,rhocentre,maxvxyzu,ui_coef
 real               :: initialtemp
 real               :: mcore,hdens,hsoft
 logical            :: iexist,input_polyk,isinkcore,isoftcore
 logical            :: use_exactN,relax_star_in_setup,write_rho_to_file
 logical            :: need_densityfile,need_rstar
 character(len=120) :: unsoftened_profile,densityfile,dens_profile
 character(len=120) :: outputfilename ! outputfilename is the path to the cored profile
 character(len=20)  :: dist_unit,mass_unit
 character(len=30)  :: lattice = 'closepacked'  ! The lattice type if stretchmap is used
 !
 ! Index of setup options
 !
 integer, parameter :: nprofile_opts =  7 ! maximum number of initial configurations
 integer, parameter :: iuniform   = 1
 integer, parameter :: ipoly      = 2
 integer, parameter :: ifromfile  = 3
 integer, parameter :: ikepler    = 4
 integer, parameter :: imesa      = 5
 integer, parameter :: ibpwpoly   = 6
 integer, parameter :: ievrard    = 7

 character(len=*), parameter :: profile_opt(nprofile_opts) = &
    (/'Uniform density             ', &
      'Polytrope                   ', &
      'Density vs r from ascii file', &
      'KEPLER star from file       ', &
      'MESA star from file         ', &
      'Piecewise polytrope         ', &
      'Evrard collapse             '/)

 public             :: setpart
 private

contains

!-----------------------------------------------------------------------
!+
!  Setup routine for stars / spherical collapse calculations
!+
!-----------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params,    only:rhozero,npart_total
 use part,            only:igas,isetphase,iphase,ihsoft
 use spherical,       only:set_sphere
 use centreofmass,    only:reset_centreofmass
 use table_utils,     only:yinterp,interpolator
 use units,           only:set_units,select_unit,utime,unit_density,unit_pressure,unit_ergg
 use kernel,          only:hfact_default
 use rho_profile,     only:rho_uniform,rho_polytrope,rho_piecewise_polytrope, &
                           rho_evrard,read_mesa_file,read_mesa,read_kepler_file, &
                           write_softened_profile
 use extern_densprofile, only:write_rhotab,rhotabfile,read_rhotab_wrapper
 use eos,             only:init_eos,init_eos_9,finish_eos,equationofstate,gmw,X_in,Z_in
 use eos_idealplusrad,only:get_idealplusrad_enfromtemp,get_idealgasplusrad_tempfrompres
 use part,            only:temperature,store_temperature
 use setstellarcore,  only:set_stellar_core
 use setsoftenedcore, only:set_softened_core,find_hsoft_given_mcore,find_mcore_given_hsoft,&
                           check_hsoft_and_mcore
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass
 use relaxstar,       only:relax_star
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer, parameter               :: ng_max = nrhotab
 integer, parameter               :: ng     = 1001
 integer                          :: i,nx,npts,ierr
 real                             :: vol_sphere,psep,rmin,densi,ri,presi
 real                             :: r(ng_max),den(ng_max),pres(ng_max),temp(ng_max),&
                                     enitab(ng_max)
 real, allocatable                :: rho0(:),r0(:),pres0(:),m0(:),ene0(:),temp0(:),&
                                     Xfrac(:),Yfrac(:)
 real                             :: xi,yi,zi,spsoundi,p_on_rhogas,pgas,eni,tempi,Y_in
 logical                          :: calc_polyk,setexists
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
 isoftcore   = .false.
 isinkcore   = .false.
 mcore         = 0.
 hsoft         = 0.
 hdens         = 0.
 isofteningopt = 1 ! By default, specify hdens
 outputfilename = 'mysoftenedstar.dat'
 densityfile  = 'ns-rdensity.tab'
 dens_profile = 'density-profile.tab'
 !
 ! defaults needed for error checking
 !
 need_iso    = 0       ! -1 = no; 0 = doesn't matter; 1 = yes
 need_temp   = 0       ! -1 = no; 0 = doesn't matter; 1 = yes
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

 if (maxvxyzu > 3  .and. need_iso == 1) call fatal('setup','require ISOTHERMAL=yes')
 if (maxvxyzu < 4  .and. need_iso ==-1) call fatal('setup','require ISOTHERMAL=no')
 if (maxvxyzu < 4  .and. need_temp==-1) call fatal('setup','require ISOTHERMAL=no')
 if (need_temp==1 .and. .not. store_temperature) call fatal('setup','require TEMPERATURE=yes')
 !
 ! set units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! setup particles
 !
 npartoftype(:) = 0
 nstar          = 0
 npart          = 0
 npart_total    = 0
 vxyzu          = 0.0
 !
 ! setup tabulated density profile
 !
 calc_polyk = .true.
 if (ieos==9) call init_eos_9(EOSopt)

 print "(/,a,/)",' Using '//trim(profile_opt(iprofile))//' profile'
 select case(iprofile)
 case(ipoly)
    call rho_polytrope(gamma,polyk,Mstar,r,den,npts,rhocentre,calc_polyk,Rstar)
    pres = polyk*den**gamma
 case(ifromfile)
    call read_rhotab_wrapper(trim(densityfile),ng_max,r,den,npts,&
                             polyk,gamma,rhocentre,Mstar,iexist,ierr)
    if (.not.iexist) call fatal('setup','density file does not exist')
    if (ierr > 0)    call fatal('setup','error in reading density file')
    pres = polyk*den**gamma
 case(ibpwpoly)
    call rho_piecewise_polytrope(r,den,rhocentre,Mstar,npts,ierr)
    if (ierr == 1) call fatal('setup','ng_max is too small')
    if (ierr == 2) call fatal('setup','failed to converge to a self-consistent density profile')
    rmin  = r(1)
    Rstar = r(npts)
    pres = polyk*den**gamma
 case(imesa)
    if (isoftcore) then
       !
       ! Read the MESA profile to be softened
       !
       call read_mesa(unsoftened_profile,rho0,r0,pres0,m0,ene0,temp0,Xfrac,Yfrac)
       rmin  = r0(1)
       Rstar = r0(size(r0))
       !
       ! Get mean molecular weight (temporary: Find gmw at R/2)
       !
       call interpolator(r0, 0.5*Rstar, i)
       pgas = pres0(i) - radconst*temp0(i)**4/3. ! Assuming pressure due to ideal gas + radiation
       gmw = (rho0(i) * kb_on_mh * temp0(i)) / pgas
       !
       ! Get representative, fixed X and Z mass fractions (Taken to be at R/2)
       !
       if ((ieos == 10) .or. (ieos == 12)) then
          X_in = Xfrac(i)
          Y_in = Yfrac(i)
          Z_in = 1 - X_in - Y_in
       endif
       !
       ! Get values of hsoft and mcore
       !
       if (isofteningopt == 1) then ! Default
          call find_mcore_given_hsoft(hdens,r0,rho0,m0,mcore,ierr)
          if (ierr==1) call fatal('setup','Cannot find mcore that produces nice profile (mcore/m(h) > 0.98 reached)')
       elseif (isofteningopt == 2) then
          call find_hsoft_given_mcore(mcore,r0,rho0,m0,hdens,ierr)
          if (ierr==1) call fatal('setup','Cannot find softening length that produces nice profile (h/r(mcore) < 1.02 reached)')
       else ! Both hdens and mcore are specified, check if values are sensible
          call check_hsoft_and_mcore(hdens,mcore,r0,rho0,m0,ierr)
          if (ierr==1) call fatal('setup','mcore cannot exceed m(r=h)')
          if (ierr==2) call fatal('setup','softenedrho/rho > tolerance')
          if (ierr==3) call fatal('setup','drho/dr > 0 found in softened profile')
       endif
       hsoft = 0.5*hdens ! This is set by default so that the pressure, energy, and temperature
       ! are same as the original profile for r > hsoft

       call set_softened_core(ieos,gamma,X_in,Y_in,gmw,mcore,hdens,hsoft,rho0,r0,pres0,m0,ene0,temp0,ierr)
       if (ierr==1) call fatal('setup','EoS not one of: adiabatic, ideal gas plus radiation, MESA in set_softened_core')
       !if (ierr==2) call fatal('setup','Xfrac and Yfrac not provided to set_softened_core for ieos=10 (MESA EoS)')
       call set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,mcore,hsoft,ihsoft)
       call write_softened_profile(outputfilename,m0,pres0,temp0,r0,rho0,ene0)
       densityfile = outputfilename ! Have the read_mesa_file subroutine read the softened profile instead
    endif

    call read_mesa_file(trim(densityfile),ng_max,npts,r,den,pres,temp,enitab,Mstar,ierr)
    rmin  = r(1)
    Rstar = r(npts)
    if (ierr==1) call fatal('setup',trim(densityfile)//' does not exist')
    if (ierr==2) call fatal('setup','insufficient data points read from file')
    if (ierr==3) call fatal('setup','too many data points; increase ng')
 case(ikepler)
    call read_kepler_file(trim(densityfile),ng_max,npts,r,den,pres,temp,enitab,Mstar,ierr)
    if (ierr==1) call fatal('setup',trim(densityfile)//' does not exist')
    if (ierr==2) call fatal('setup','insufficient data points read from file')
    if (ierr==3) call fatal('setup','too many data points; increase ng')
    rmin  = r(1)
    Rstar = r(npts)
 case(ievrard)
    call rho_evrard(ng,Mstar,Rstar,r,den)
    npts = ng
    polyk = ui_coef*Mstar/Rstar
    pres = polyk*den**gamma
    print*,' Assuming polyk = ',polyk
 case default  ! set up uniform sphere by default
    call rho_uniform(ng,Mstar,Rstar,r,den) ! use this array for continuity of call to set_sphere
    npts = ng
    pres = polyk*den**gamma
    print*,' Assuming polyk = ',polyk
 end select
 !
 ! place particles in sphere
 !
 vol_sphere  = 4./3.*pi*Rstar**3
 nx          = int(np**(1./3.))
 psep        = vol_sphere**(1./3.)/real(nx)
 call set_sphere(lattice,id,master,rmin,Rstar,psep,hfact,npart,xyzh, &
                 rhotab=den(1:npts),rtab=r(1:npts),nptot=npart_total, &
                 exactN=use_exactN,np_requested=np)
 !
 ! add sink particle stellar core
 !
 if (isinkcore) call set_stellar_core(nptmass,xyzmh_ptmass,vxyz_ptmass,mcore,hsoft,ihsoft)

 nstar = int(npart_total,kind=(kind(nstar)))
 massoftype(igas) = Mstar/nstar
 !
 ! set total particle number (on this MPI thread)
 !
 npart             = nstar
 npartoftype(igas) = npart
 iphase(1:npart)   = isetphase(igas,iactive=.true.)
 !
 ! Write the desired profile to file (do this before relaxation)
 !
 if (write_rho_to_file) call write_rhotab(dens_profile,r,den,npts,polyk,gamma,rhocentre,ierr)
 !
 ! relax the density profile to achieve nice hydrostatic equilibrium
 !
 call init_eos(ieos,ierr)
 if (relax_star_in_setup) then
    if (nstar==npart) then
       call relax_star(npts,den,pres,r,npart,xyzh)
    else
       call error('setup_star','cannot run relaxation with MPI setup, please run setup on ONE MPI thread')
    endif
 endif
 !
 ! set the thermal energy / temperature profile of the star
 !
 do i=1,nstar
    if (maxvxyzu==4) then
       !
       !  Interpolate density and pressure from table
       !
       ri    = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
       densi = yinterp(den(1:npts),r(1:npts),ri)
       presi = yinterp(pres(1:npts),r(1:npts),ri)
       !
       ! Set internal energy of particles given pressure and density depending on EoS
       !
       select case(ieos)
       case(16) ! Shen EoS
          vxyzu(4,i) = initialtemp
       case(15) ! Helmholtz EoS
          xi    = xyzh(1,i)
          yi    = xyzh(2,i)
          zi    = xyzh(3,i)
          tempi = initialtemp
          call equationofstate(ieos,p_on_rhogas,spsoundi,densi,xi,yi,zi,eni,tempi)
          vxyzu(4,i) = eni
          if (store_temperature) temperature(i) = initialtemp
       case(12) ! Ideal gas plus radiation EoS
          call get_idealgasplusrad_tempfrompres(presi*unit_pressure,densi*unit_density,gmw,tempi)
          call get_idealplusrad_enfromtemp(densi*unit_density,tempi,gmw,eni)
          vxyzu(4,i) = eni / unit_ergg
          if (store_temperature) temperature(i) = tempi
       case(10) ! MESA EoS
          vxyzu(4,i) = yinterp(enitab(1:npts),r(1:npts),ri)
          if (store_temperature) temperature(i) = yinterp(temp(1:npts),r(1:npts),ri)
       case default
          if (gamma < 1.00001) then
             vxyzu(4,i) = polyk
          else
             vxyzu(4,i) = presi / ((gamma - 1.) * densi)
          endif
       end select
    endif
 enddo
 call finish_eos(ieos,ierr)
 !
 ! Reset centre of mass (again)
 !
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 !
 ! Print summary to screen
 !
 rhozero = Mstar/vol_sphere
 write(*,"(70('='))")
 if (ieos /= 9) then
    write(*,'(1x,a,f12.5)')       'gamma               = ', gamma
 endif
 if (maxvxyzu <= 4) then
    write(*,'(1x,a,f12.5)')       'polyk               = ', polyk
    write(*,'(1x,a,f10.6,a)')     'specific int. energ = ', polyk*Rstar/Mstar,' GM/R'
 endif
 call write_mass('particle mass       = ',massoftype(igas),umass)
 call write_dist('Radius              = ',Rstar,udist)
 call write_mass('Mass                = ',Mstar,umass)
 if (iprofile==ipoly) then
    write(*,'(1x,a,es12.5,a)') 'rho_central         = ', rhocentre*unit_density,' g/cm^3'
 endif
 write(*,'(1x,a,i12)')         'N                   = ', nstar
 write(*,'(1x,a,2(es12.5,a))')    'rho_mean            = ', rhozero*unit_density,  ' g/cm^3 = '&
                                                       , rhozero,               ' code units'
 write(*,'(1x,a,es12.5,a)')       'free fall time      = ', sqrt(3.*pi/(32.*rhozero))*utime,' s'
 call write_dist('particle separation = ',psep,udist)
 if ( (iprofile==iuniform .and. .not.gravity) .or. iprofile==ifromfile) then
    write(*,'(a)') 'WARNING! This setup may not be stable'
 endif
 write(*,"(70('='))")

end subroutine setpart

!-----------------------------------------------------------------------
!+
!  Ask questions of the user to determine which setup to use
!+
!-----------------------------------------------------------------------
subroutine setup_interactive(polyk,gamma,iexist,id,master,ierr)
 use prompting, only:prompt
 use units,     only:select_unit
 real, intent(out)    :: polyk,gamma
 logical, intent(in)  :: iexist
 integer, intent(in)  :: id,master
 integer, intent(out) :: ierr
 integer :: i
 logical :: need_densityfile,need_rstar

 ierr = 0
 ! Select sphere & set default values
 write(*,*)
 do i = 1, nprofile_opts
    write(*,"(a5,i2,1x,a48)") 'Case ', i, profile_opt(i)
 enddo

 call prompt('Enter which density profile to use',iprofile,1,nprofile_opts)
 !
 ! set default file output parameters
 !
 if (id==master) write(*,"('Setting up ',a)") trim(profile_opt(iprofile))
 call set_defaults_given_profile(iprofile,iexist,need_densityfile,need_rstar,polyk)

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
 end select
 if (iprofile==ievrard) then
    call prompt('Enter the specific internal energy (units of GM/R) ',ui_coef,0.)
 endif

 ! star properties
 if (need_densityfile) then
    call prompt('Enter file name containing density profile ', densityfile)
 else
    call prompt('Enter the mass of the star (code units) ',  Mstar,0.)
    if (need_rstar) call prompt('Enter the radius of the star (code units) ',Rstar,0.)
 endif

 if (iprofile==imesa) then
    call prompt('Soften the core density profile and add a sink particle core?',isoftcore)
    if (.not. isoftcore) then
       call prompt('Add a sink particle stellar core?',isinkcore)
    else
       isinkcore = .true. ! Create sink particle core automatically
       unsoftened_profile = densityfile
       print*,'Options for core softening:'
       print "(5(/,a))",'1. Specify radius of density softening', &
                        '2. Specify mass of sink particle core (not recommended)', &
                        '3. Specify both radius of density softening length and mass', &
                        '   of sink particle core (if you do not know what you are', &
                        '   doing, you will obtain a poorly softened profile)'
       call prompt('Select option above : ',isofteningopt,1,3)
       select case(isofteningopt)
       case(1)
          call prompt('Enter radius of density softening',hdens,0.)
       case(2)
          call prompt('Enter mass of the created sink particle core',mcore,0.)
       case(3)
          call prompt('Enter mass of the created sink particle core',mcore,0.)
          call prompt('Enter radius of density softening',hdens,0.)
       end select
       call prompt('Enter output file name of cored stellar profile:',outputfilename)
    endif
    if (isinkcore .and. (.not. isoftcore)) then
       call prompt('Enter mass of the created sink particle core',mcore,0.)
       call prompt('Enter softening length of the sink particle core',hsoft,0.)
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
subroutine set_defaults_given_profile(iprofile,iexist,need_densityfile,need_rstar,polyk)
 integer, intent(in)  :: iprofile
 logical, intent(in)  :: iexist
 logical, intent(out) :: need_densityfile,need_rstar
 real,    intent(out) :: polyk

 need_rstar = .true.
 need_densityfile = .false.
 select case(iprofile)
 case(iuniform)
    input_polyk  = .true.
 case(ipoly)
    input_polyk  = .false.
 case(ifromfile)
    ! Read the density profile from file (e.g. for neutron star)
    !  Original Author: Mark Bennett
    densityfile = 'ns-rdensity.tab'
    need_densityfile = .true.
 case(imesa)
    ! sets up a star from a 1D MESA code output
    !  Original Author: Roberto Iaconi
    densityfile = 'P12_Phantom_Profile.data'
    need_densityfile = .true.
 case(ikepler)
    ! sets up a star from a 1D KEPLER code output
    !  Original Author: Nicole Rodrigues
    densityfile = 'kepler_MS.data'
    need_densityfile = .true.
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
 use infile_utils, only:write_inopt,get_optstring
 use dim,          only:tagline
 use relaxstar,    only:write_options_relax
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
 case(1)
    if (input_polyk) call write_inopt(polyk,'polyk','polytropic constant (cs^2 if isothermal)',iunit)
 end select
 if (iprofile==ievrard) then
    call write_inopt(ui_coef,'ui_coef','specific internal energy (units of GM/R)',iunit)
 endif

 write(iunit,"(/,a)") '# star properties'
 if (need_densityfile) then
    call write_inopt(densityfile,'densityfile','File containing data for stellar profile',iunit)
 else
    call write_inopt(Rstar,'Rstar','radius of star',iunit)
    call write_inopt(Mstar,'Mstar','mass of star',iunit)
 endif

 if (iprofile==imesa) then
    write(iunit,"(/,a)") '# core softening and sink stellar core options'
    call write_inopt(isoftcore,'isoftcore','Soften the core of an input MESA profile',iunit)
    if (isoftcore) then
       call write_inopt(isofteningopt,'isofteningopt','1=supply hsoft, 2=supply mcore, 3=supply both',iunit)
       call write_inopt(unsoftened_profile,'unsoftened_profile','Path to MESA profile for softening',iunit)
       call write_inopt(outputfilename,'outputfilename','Output path for softened MESA profile',iunit)
       call write_inopt(hdens,'hdens','Radius of core softening',iunit)
       call write_inopt(mcore,'mcore','Mass of sink particle stellar core',iunit)
       call write_inopt(hsoft,'hsoft','Softening length of sink particle stellar core',iunit)
    endif
    call write_inopt(isinkcore,'isinkcore','Add a sink particle stellar core',iunit)
    if (isinkcore .and. (.not. isoftcore)) then
       call write_inopt(mcore,'mcore','Mass of sink particle stellar core',iunit)
       call write_inopt(hsoft,'hsoft','Softening length of sink particle stellar core',iunit)
    endif
 endif

 write(iunit,"(/,a)") '# relaxation options'
 call write_inopt(relax_star_in_setup,'relax_star','relax star automatically during setup',iunit)
 if (relax_star_in_setup) call write_options_relax(iunit)

 call write_inopt(write_rho_to_file,'write_rho_to_file','write density profile to '// &
                  trim(dens_profile),iunit)

 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,polyk,ierr)
 use infile_utils, only:open_db_from_file,inopts,close_db,read_inopt
 use io,           only:error
 use units,        only:select_unit
 use relaxstar,    only:read_options_relax
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
 call set_defaults_given_profile(iprofile,iexist,need_densityfile,need_rstar,polyk)

 ! units
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)

 ! resolution
 call read_inopt(np,'np',db,errcount=nerr)
 call read_inopt(use_exactN,'use_exactN',db,errcount=nerr)
 nstar = np

 ! star properties
 if (need_densityfile) then
    call read_inopt(densityfile,'densityfile',db,errcount=nerr)
 else
    call read_inopt(Mstar,'Mstar',db,errcount=nerr)
    if (need_rstar) call read_inopt(Rstar,'Rstar',db,errcount=nerr)
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
    if (input_polyk) call read_inopt(polyk,'polyk',db,errcount=nerr)
 case(1)
    if (input_polyk) call read_inopt(polyk,'polyk',db,errcount=nerr)
 end select
 if (iprofile==ievrard) then
    call read_inopt(ui_coef,'ui_coef',db,errcount=nerr)
 endif

 ! core softening
 if (iprofile==imesa) then
    call read_inopt(isoftcore,'isoftcore',db,errcount=nerr)
    if (isoftcore) then
       call read_inopt(isofteningopt,'isofteningopt',db,errcount=nerr)
       call read_inopt(unsoftened_profile,'unsoftened_profile',db,errcount=nerr)
       call read_inopt(outputfilename,'outputfilename',db,errcount=nerr)
       if (ieos==2) call read_inopt(gamma,'gamma',db,errcount=nerr)
       if (isofteningopt==1) call read_inopt(hdens,'hdens',db,errcount=nerr)
       if (isofteningopt==2) call read_inopt(mcore,'mcore',db,errcount=nerr)
       if (isofteningopt==3) then
          call read_inopt(mcore,'mcore',db,errcount=nerr)
          call read_inopt(hdens,'hdens',db,errcount=nerr)
       endif
    endif

    ! sink particle core
    call read_inopt(isinkcore,'isinkcore',db,errcount=nerr)
    if (isinkcore .and. (.not. isoftcore)) then
       call read_inopt(mcore,'mcore',db,errcount=nerr)
       call read_inopt(hsoft,'hsoft',db,errcount=nerr)
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
