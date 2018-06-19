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
!    This module sets up sphere(s).  There are multiple options, incuding
!    1) uniform unit sphere
!    2) single polytrope
!    3) binary polytrope (decommissioned)
!    4) neutron star from file
!    5) red giant (Macquarie)
!    6) neutron star using a piecewise polytrope EOS
!    7) Evrard sphere
!    8) KEPLER star from file
!    9) Kelvin-Helmholtz Equation of state
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    EOSopt      -- EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)
!    Mstar_1     -- mass of star
!    Nstar_1     -- particle number in the star
!    Rstar_1     -- radius of star
!    densityfile -- File containing data for stellar profile
!    dist_unit   -- distance unit (e.g. au)
!    gamma       -- Adiabatic index
!    initialtemp -- initial temperature of the star
!    mass_unit   -- mass unit (e.g. solarm)
!    np          -- approx number of particles (in box of size 2R)
!    polyk       -- sound speed .or. constant in EOS
!    ui_coef     -- specific internal energy (units of GM/R)
!
!  DEPENDENCIES: centreofmass, dim, eos, extern_neutronstar,
!    externalforces, infile_utils, io, kernel, options, part, physcon,
!    prompting, rho_profile, setup_params, spherical, table_utils,
!    timestep, units
!+
!--------------------------------------------------------------------------
module setup
 use io,             only: fatal
 use part,           only: gravity
 use physcon,        only: solarm,solarr,km,pi,c
 use options,        only: nfulldump,iexternalforce
 use timestep,       only: tmax,dtmax
 use eos,            only: ieos, p1pwpcgs,gamma1pwp,gamma2pwp,gamma3pwp
 use externalforces, only: iext_neutronstar
 use dim,            only: calc_erot

 implicit none
 !
 ! Input parameters
 !
 integer, parameter :: numEOS       =  4 ! maximum number of piecewise polytrope defaults
 integer, parameter :: numparam     =  4 ! number of parameters governing the piecewise polytrope
 integer            :: isphere,np,EOSopt
 integer            :: nstar
 real(kind=8)       :: udist,umass
 real               :: Rstar,Mstar,rhocentre,maxvxyzu,ui_coef
 real               :: initialtemp
 logical            :: set_vcirc,iexist,input_polyk
 logical            :: use_exactN,use_prompt
 character(len=120) :: densityfile
 character(len=20)  :: dist_unit,mass_unit
 character(len=30)  :: lattice = 'closepacked'  ! The lattice type if stretchmap is used
 !
 ! Index of setup options
 !
 integer, parameter :: nsphere_opts =  9 ! maximum number of initial configurations
 integer, parameter :: iuniform   = 1
 integer, parameter :: ipoly      = 2
 integer, parameter :: ibinpoly   = 3
 integer, parameter :: insfile    = 4
 integer, parameter :: ired       = 5
 integer, parameter :: ibpwpoly   = 6
 integer, parameter :: ievrard    = 7
 integer, parameter :: ikepler    = 8
 integer, parameter :: ihelmholtz = 9

 character(len=48)  :: sphere_opt(nsphere_opts)

 public             :: setpart
 private

contains

!-----------------------------------------------------------------------
!+
!  Setup routine for stars / spherical collapse calculations
!+
!-----------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,             only: master
 use part,           only: igas
 use prompting,      only: prompt
 use spherical,      only: set_sphere
 use centreofmass,   only: reset_centreofmass
 use table_utils,    only: yinterp
 use units,          only: set_units,select_unit,utime,unit_density
 use kernel,         only: hfact_default
 use rho_profile,    only: rho_uniform,rho_polytrope,rho_piecewise_polytrope, &
                           rho_evrard,read_red_giant_file,read_kepler_file
 use extern_neutronstar, only: write_rhotab,rhotabfile,read_rhotab_wrapper
 use eos,            only: init_eos, finish_eos, equationofstate
 use part,           only: rhoh, temperature, store_temperature
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 integer, parameter               :: ng_max = 5000
 integer, parameter               :: ng     = 1024
 integer                          :: i,nx,npts,npmax,ierr
 real                             :: vol_sphere,psep,rmin,densi,ri,polyk_in,presi
 real                             :: r(ng_max),den(ng_max),pres(ng_max),temp(ng_max),enitab(ng_max)
 real                             :: xi, yi, zi, rhoi, spsoundi, p_on_rhogas, eni, tempi
 logical                          :: calc_polyk,write_setup
 character(len=120)               :: setupfile,inname
 !
 !
 ! Initialise parameters, including those that will not be included in *.setup
 !
 time         = 0.
 polyk        = 1.0
 gamma        = 5./3.
 hfact        = hfact_default
 maxvxyzu     = size(vxyzu(:,1))
 write_setup  = .false.
 use_prompt   = .true.   ! allow user input
 use_exactN   = .true.
 call set_option_names()
 !
 ! determine if an .in file exists
 !
 inname=trim(fileprefix)//'.in'
 inquire(file=inname,exist=iexist)
 !
 ! determine if an .setup file exists
 !
 setupfile = trim(fileprefix)//'.setup'
 call read_setupfile(setupfile,gamma,polyk_in,ierr)
 if ( (ierr /= 0 .or. .not.iexist) .and. id==master) then
    ! setup file does not exist or is incomplete
    !
    ! Select sphere & set default values
    call choose_spheres(polyk,iexist,id,master)
    if (isphere==ibinpoly) call fatal('setup','Binary polytope option has been removed.')
    !
    ! Using prompts, determine the parameters the users wishes
    !
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
    npmax = int(2.0/3.0*size(xyzh(1,:))) ! approx max number allowed in sphere given size(xyzh(1,:))
    np    = min(10000,npmax) ! default number of particles
    call prompt('Enter the approximate number of particles in the sphere ',np,0,npmax)
    if (isphere==insfile .or. isphere==ired .or. isphere==ikepler) then
       call prompt('Enter file name containing density profile ', densityfile)
    endif
    if ( use_prompt ) then
       if (isphere==ibpwpoly) then
          write(*,'(a)') 'EOS options: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)'
          call prompt('Enter equation of state type',EOSopt,1,4)
       endif
       call prompt('Enter the radius of the star (code units) ',Rstar,0.)
       call prompt('Enter the mass of the star (code units) ',  Mstar,0.)
       if (isphere==ievrard) then
          call prompt('Enter the specific internal energy (units of GM/R) ',ui_coef,0.)
       endif
       if (isphere/=ibpwpoly) then
          call prompt('Enter the Adiabatic index',gamma,1.)
       endif
       if (input_polyk) then
          call prompt('Enter polyk (sound speed .or. constant in EOS calculation)',polyk,0.)
       endif
       if (isphere==ihelmholtz) then
          call prompt('Enter temperature',initialtemp,1.0e4,1.0e11)
       endif
    endif
    write_setup = .true.
 else
    ! Set default values that are not in the .setup file
    np       = nstar
    if (input_polyk) polyk = polyk_in
 endif
 !
 ! set units
 !
 call set_units(dist=udist,mass=umass,G=1.d0)
 !
 ! set equation of state
 !
 if (isphere==ibpwpoly .and. .not. iexist) call init_ieos9
 !
 ! setup particles
 !
 npartoftype(:) = 0
 nstar          = 0
 npart          = 0
 npart_total    = 0
 vxyzu          = 0.0
 if (isphere==ievrard) polyk = ui_coef*Mstar/Rstar

 !
 ! setup tabulated density profile
 !
  calc_polyk = .true.
 if (isphere==iuniform) then
    call rho_uniform(ng,Mstar,Rstar,r,den) ! use this array for continuity of call to set_sphere
    npts = ng
 elseif (isphere==ipoly .or. isphere==ihelmholtz) then
    call rho_polytrope(gamma,polyk,Mstar,r,den,npts,rhocentre,calc_polyk,Rstar)
 elseif (isphere==insfile) then
    call read_rhotab_wrapper(trim(densityfile),ng_max,r,den,npts,&
                             polyk,gamma,rhocentre,Mstar,iexist,ierr)
    if (.not.iexist) call fatal('setup','density file does not exist')
    if (ierr > 0)    call fatal('setup','error in reading density file')
 elseif (isphere==ibpwpoly) then
    call rho_piecewise_polytrope(r,den,rhocentre,Mstar,npts,ierr)
    if (ierr == 1) call fatal('setup','ng_max is too small')
    if (ierr == 2) call fatal('setup','failed to converge to a self-consistent density profile')
    rmin  = r(1)
    Rstar = r(npts)
 elseif (isphere==ired) then
    call read_red_giant_file(trim(densityfile),ng_max,npts,r,den,pres,temp,enitab,Mstar,ierr)
    if (ierr==1) call fatal('setup',trim(densityfile)//' does not exist')
    if (ierr==2) call fatal('setup','insufficient data points read from file')
    if (ierr==3) call fatal('setup','too many data points; increase ng')
    rmin  = r(1)
    Rstar = r(npts)
 elseif (isphere==ikepler) then
    call read_kepler_file(trim(densityfile),ng_max,npts,r,den,pres,temp,enitab,Mstar,ierr)
    if (ierr==1) call fatal('setup',trim(densityfile)//' does not exist')
    if (ierr==2) call fatal('setup','insufficient data points read from file')
    if (ierr==3) call fatal('setup','too many data points; increase ng')
    rmin  = r(1)
    Rstar = r(npts)
 elseif (isphere==ievrard) then
    call rho_evrard(ng,Mstar,Rstar,r,den)
    npts = ng
 endif
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
 ! reset centre of mass
 !
 nstar = int(npart_total,kind=(kind(nstar)))
 massoftype(igas) = Mstar/nstar
 call reset_centreofmass(nstar,xyzh(:,1:nstar),vxyzu(:,1:nstar))
 !
 ! add energies
 !
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('setup_spheres','error initialising equation of state')
 do i=1,nstar
    if (maxvxyzu==4) then
       if (gamma < 1.00001 .or. isphere==ievrard) then
          vxyzu(4,i) = polyk
       else
          !
          !  Note: Running the polytrope with u stored is not quite
          !  the same as using P = K rho^gamma because we really
          !  should use the actual rho, not rhozero.
          !  Note also that both spheres are currently at the origin
          !
          ri         = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
          densi      = yinterp(den(1:npts),r(1:npts),ri)
          vxyzu(4,i) = polyk*densi**(gamma-1.)/(gamma-1.)
          if ((isphere==ired .or. isphere==ikepler) .and. (ieos==10)) then
             vxyzu(4,i) = yinterp(enitab(1:npts),r(1:npts),ri)
             presi = yinterp(pres(1:npts),r(1:npts),ri)
          else if ((isphere==ired .or. isphere==ikepler) .and. (ieos/=10)) then
             presi = yinterp(pres(1:npts),r(1:npts),ri)
             vxyzu(4,i) = presi / ((gamma - 1.) * densi)
          else if (isphere==ihelmholtz .and. ieos==15) then
             ! set internal energy according to Helmholtz eos
             xi    = xyzh(1,i)
             yi    = xyzh(2,i)
             zi    = xyzh(3,i)
             rhoi  = rhoh(xyzh(4,i),massoftype(igas))
             tempi = initialtemp
             call equationofstate(ieos,p_on_rhogas,spsoundi,rhoi,xi,yi,zi,eni,tempi)
             vxyzu(4,i) = eni
             temperature(i) = initialtemp
          endif
       endif
    endif
 enddo
 call finish_eos(ieos,ierr)
 !
 ! set total particle number
 !
 npart          = nstar
 npartoftype(1) = npart
 !
 ! Reset centre of mass (again)
 !
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! Write the .setup file if necessary; performed here to allow for exact writing of Mstar(2)
 ! (Also rewrite if nstar has changed; may occur if user manually modifies the .setup file)
 !
 if (write_setup) call write_setupfile(setupfile,gamma,polyk)
 !
 ! Write the neutronstar profile to file (if isphere==insfile)
 !
 if (isphere==insfile) then
    call write_rhotab(r,den,npts,polyk,gamma,rhocentre,ierr)
    iexternalforce = iext_neutronstar
 endif
 !
 ! Print summary to screen
 !
 rhozero = Mstar/vol_sphere
 write(*,'(a)') "======================================================================"
 if (isphere/=ibpwpoly) then
    write(*,'(a,F12.5)')       'gamma               = ', gamma
 endif
 if (isphere/=ievrard) then
    write(*,'(a,F12.5)')       'polyk               = ', polyk
 endif
 call write_mass('particle mass       = ',massoftype(igas),umass)
 call write_dist('Radius              = ',Rstar,udist)
 call write_mass('Mass                = ',Mstar,umass)
 if (isphere==ipoly .or. isphere==ihelmholtz) then
    write(*,'(a,Es12.5,a)') 'rho_central         = ', rhocentre*unit_density,' g/cm^3'
 endif
 write(*,'(a,I12)')         'N                   = ', nstar
 write(*,'(a,2(Es12.5,a))')    'rho_mean            = ', rhozero*unit_density,  ' g/cm^3 = '&
                                                       , rhozero,               ' code units'
 write(*,'(a,Es12.5,a)')       'free fall time      = ', sqrt(3.*pi/(32.*rhozero))*utime,' s'
 if (isphere==ievrard) then
    write(*,'(a,F10.6,a)')     'specific internal energy = ', ui_coef,           ' GM/R'
 endif
 call    write_dist('particle separation = ',psep,udist)
 if ( (isphere==iuniform .and. .not.gravity) .or. isphere==insfile) then
    write(*,'(a)') 'WARNING! This setup may not be stable'
 endif
 write(*,'(a)') "======================================================================"
end subroutine setpart
!-----------------------------------------------------------------------
!+
!  The default piecewise polytrope options, as per Read et al (2009)
!  The unlisted values are common to all options; all values are in cgs
!  The array is
!  pw(i,:) = (/ prescrit,gamma1,gamma2,gamma3 /)
!  pw(:,j) = (/ ARP3,SLy,MS1,ENG/)
!+
!-----------------------------------------------------------------------
subroutine init_ieos9
 real :: pw(numEOS,numparam)
 !
 ! Define the default options
 !
 pw(1,:)  = (/ 10**34.392, 3.166, 3.573, 3.281 /)
 pw(2,:)  = (/ 10**34.384, 3.005, 2.988, 2.851 /)
 pw(3,:)  = (/ 10**34.858, 3.224, 3.033, 1.325 /)
 pw(4,:)  = (/ 10**34.437, 3.514, 3.130, 3.168 /)
 !
 ! Choose the default option
 !
 p1pwpcgs  = pw(EOSopt,1)
 gamma1pwp = pw(EOSopt,2)
 gamma2pwp = pw(EOSopt,3)
 gamma3pwp = pw(EOSopt,4)

end subroutine init_ieos9
!-----------------------------------------------------------------------
!+
!  Set the names for the different options
!+
!-----------------------------------------------------------------------
subroutine set_option_names()

 sphere_opt(:)          = 'none'
 sphere_opt(iuniform)   = 'Uniform density sphere'
 sphere_opt(ipoly)      = 'Polytrope'
 sphere_opt(ibinpoly)   = 'Binary polytrope'
 sphere_opt(insfile)    = 'neutron star from file'
 sphere_opt(ired)       = 'Red giant (Macquarie)'
 sphere_opt(ibpwpoly)   = 'neutron star with piecewise polytrope EOS'
 sphere_opt(ievrard)    = 'Evrard collapse'
 sphere_opt(ikepler)    = 'KEPLER star from file'
 sphere_opt(ihelmholtz) = 'Helmholtz free energy eos star'

end subroutine set_option_names
!-----------------------------------------------------------------------
!+
!  Select shock type & set default values
!+
!-----------------------------------------------------------------------
subroutine choose_spheres(polyk,iexist,id,master)
 use prompting,   only: prompt
 use part,        only: store_temperature
 integer, intent(in)  :: id, master
 logical, intent(in)  :: iexist
 real,    intent(out) :: polyk
 integer              :: i,choice,need_grav,need_iso,need_temp

 write(*,*)
 do i = 1, nsphere_opts
    if (trim(sphere_opt(i)) /= 'none') write(*,"(a5,i2,1x,a48)") 'Case ', i, sphere_opt(i)
 enddo

 choice = 1
 call prompt('Enter which setup to use',choice,1,nsphere_opts)
 !
 ! General defaults
 !
 dist_unit   = 'solarr'
 mass_unit   = 'solarm'
 Rstar       = 1.0
 Mstar       = 1.0
 ui_coef     = 1.0
 need_grav   = 1       ! -1 = no; 0 = doesn't matter; 1 = yes
 need_iso    = 0       ! -1 = no; 0 = doesn't matter; 1 = yes
 need_temp   = 0       ! -1 = no; 0 = doesn't matter; 1 = yes
 EOSopt      = 1
 if (.not. iexist) then
    tmax      = 100.
    dtmax     =   1.0
    ieos      =   2
    nfulldump =  10
 endif
 set_vcirc    = .true.
 input_polyk  = .false.
 !
 ! set default file output parameters
 !
 if (id==master) write(*,"('Setting up ',a)") trim(sphere_opt(choice))
 select case (choice)
 case(iuniform)
    ! Uniform density sphere
    polyk       = 0.5
    input_polyk = .true.
    need_grav   = 0 ! to prevent setupfail
 case(ipoly)
    ! Polytrope
    need_iso = 0 ! can be done either with du/dt=P/rho^2 drho/dt or with P=K*rho**gamma
 case(ibinpoly)
    ! Binary Polytrope
    print*, 'This option has been removed.  Please rerun using option ',ipoly, &
            ' and then create the second star using moddump.'
 case(insfile)
    ! Read the density profile from file (for neutron star)
    !  Original Author: Mark Bennett
    !  Note: original densityfile is missing, thus this is the polytrope from
    !        case ipoly; it may *not* be stable
    densityfile = 'ns-rdensity.tab'
    need_grav      = -1
 case(ired)
    ! sets up a star from a 1D code output
    !  Original Author: Roberto Iaconi
    use_exactN  = .false.
    use_prompt  = .false.
    densityfile = 'P12_Phantom_Profile.data'
    call prompt('Enter the desired EoS for setup', ieos)
 case(ikepler)
    ! sets up a star from a 1D code output
    !  Original Author: Nicole Rodrigues
    !  Supervisors: Daniel Price & Alexander Heger
    polyk       = 0.35899
    use_exactN  = .false.
    use_prompt  = .false.
    densityfile = 'kepler_MS.data'
    call prompt('Enter the desired EoS for setup', ieos)
 case(ibpwpoly)
    ! simulates the merger of a pair of neutron stars
    !  Original Author: Madeline Marshall & Bernard Field
    !  Supervisors: James Wurster & Paul Lasky
    if (.not. iexist) then
       tmax           = 5000.0
       dtmax          =    7.5
       ieos           =    9
    endif
    dist_unit    = 'km'
    Mstar        =  1.35
    polyk        = 144.
    calc_erot    = .true.
 case(ievrard)
    ! Evrard Collapse
    if (.not. iexist) then
       tmax      = 3.0
       dtmax     = 0.1
       nfulldump = 1
    endif
    ui_coef     = 0.05
    need_iso    = -1
 case(ihelmholtz)
    ! set up initial polytrope but use Helmholtz free energy eos
    tmax        = 1.0
    dtmax       = 2.0e-4
    nfulldump   = 10
    need_iso    = -1
    ieos        = 15
    initialtemp = 1.0e7
    Rstar       = 0.01
    Mstar       = 0.6
    need_temp   = 1
 end select
 isphere = choice
 !
 ! Verify correct pre-processor commands
 !
 if (      gravity .and. need_grav==-1) call fatal('setup','require GRAVITY=no')
 if (.not. gravity .and. need_grav== 1) call fatal('setup','require GRAVITY=yes')
 if (maxvxyzu > 3  .and. need_iso == 1) call fatal('setup','require ISOTHERMAL=yes')
 if (maxvxyzu < 4  .and. need_iso ==-1) call fatal('setup','require ISOTHERMAL=no')
 if (maxvxyzu < 4  .and. need_temp==-1) call fatal('setup','require ISOTHERMAL=no')
 if (need_temp==1 .and. .not. store_temperature) call fatal('setup','require TEMPERATURE=yes')

end subroutine choose_spheres
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
    write(*,'(2(a,Es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' R_sun'
 else if ( abs(1.0-km/udist) < 1.0d-4) then
    write(*,'(2(a,Es12.5),a)') item_in, dist_in*udist,' cm     = ',dist_in,' km'
 else
    write(*,'(a,Es12.5,a)')    item_in, dist_in*udist,' cm'
 endif

end subroutine write_dist

subroutine write_mass(item_in,mass_in,umass)
 real,             intent(in) :: mass_in
 real(kind=8),     intent(in) :: umass
 character(len=*), intent(in) :: item_in

 if ( abs(1.0-solarm/umass) < 1.0d-4) then
    write(*,'(2(a,Es12.5),a)') item_in, mass_in*umass,' g      = ',mass_in,' M_sun'
 else
    write(*,'(a,Es12.5,a)')    item_in, mass_in*umass,' g'
 endif

end subroutine write_mass
!-----------------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!-----------------------------------------------------------------------
subroutine write_setupfile(filename,gamma,polyk)
 use infile_utils, only: write_inopt,get_optstring
 use dim,          only: tagline
 real,             intent(in) :: gamma,polyk
 character(len=*), intent(in) :: filename
 integer,          parameter  :: iunit = 20
 character(len=120)           :: string

 write(*,"(a)") ' Writing '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# '//trim(tagline)
 write(iunit,"(2a)") '# input file for Phantom spheres setup: ',sphere_opt(isphere)

 call get_optstring(nsphere_opts,sphere_opt,string,4)
 call write_inopt(isphere,'isphere',trim(string),iunit)

 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm)',iunit)

 if (isphere==ibpwpoly) then
    write(iunit,"(/,a)") '# Piecewise Polytrope default options'
    call write_inopt(EOSopt,'EOSopt','EOS: 1=APR3,2=SLy,3=MS1,4=ENG (from Read et al 2009)',iunit)
 endif

 write(iunit,"(/,a)") '# sphere & particle numbers'
 if ( .not. use_exactN ) then
    call write_inopt(np,'np','approx number of particles (in box of size 2R)',iunit)
 endif
 call write_inopt(nstar,'Nstar_1','particle number in the star',iunit)

 write(iunit,"(/,a)") '# sphere properties'
 if (isphere==insfile .or. isphere==ired .or. isphere==ikepler) then
    call write_inopt(densityfile,'densityfile','File containing data for stellar profile',iunit)
 endif
 if (use_prompt) then
    if (isphere/=ibpwpoly) then
       call write_inopt(Rstar,'Rstar_1','radius of star',iunit)
    endif
    call write_inopt(Mstar,'Mstar_1','mass of star',iunit)
    if (isphere/=ibpwpoly) then
       call write_inopt(gamma,'gamma','Adiabatic index',iunit)
    endif
    if (input_polyk) then
       call write_inopt(polyk,'polyk','sound speed .or. constant in EOS',iunit)
    endif
    if (isphere==ievrard) then
       call write_inopt(ui_coef,'ui_coef','specific internal energy (units of GM/R)',iunit)
    endif
    if (isphere==ihelmholtz) then
       call write_inopt(initialtemp,'initialtemp','initial temperature of the star',iunit)
    endif
 endif

 close(iunit)

end subroutine write_setupfile
!-----------------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!-----------------------------------------------------------------------
subroutine read_setupfile(filename,gamma,polyk,ierr)
 use infile_utils, only: open_db_from_file,inopts,close_db,read_inopt
 use io,           only: error
 use units,        only: select_unit
 character(len=*), intent(in)  :: filename
 integer,          parameter   :: lu = 21
 integer,          intent(out) :: ierr
 real,             intent(out) :: gamma,polyk
 integer                       :: nerr
 type(inopts), allocatable     :: db(:)
 !
 call open_db_from_file(db,filename,lu,ierr)
 if (ierr /= 0) return
 write(*, '(1x,2a)') 'Setup_spheres: Reading setup options from ',trim(filename)
 !
 nerr = 0
 call read_inopt(isphere,'isphere',db,errcount=nerr)
 if (isphere==ired .or. isphere==ikepler) then
    use_prompt     = .false.
    use_exactN     = .false.
 endif
 call read_inopt(mass_unit,'mass_unit',db,ierr)
 call read_inopt(dist_unit,'dist_unit',db,ierr)
 if (isphere==insfile .or. isphere==ired .or. isphere==ikepler) then
    call read_inopt(densityfile,'densityfile',db,errcount=nerr)
 endif
 if (use_prompt) then
    if (isphere==ibpwpoly) then
       call read_inopt(EOSopt,'EOSopt',db,errcount=nerr)
    endif
 endif
 if ( .not. use_exactN ) then
    call read_inopt(np,'np',db,errcount=nerr)
 endif
 call read_inopt(nstar,'Nstar_1',db,errcount=nerr)
 if (use_prompt) then
    if (isphere/=ibpwpoly) then
       call read_inopt(Rstar,'Rstar_1',db,errcount=nerr)
    endif
    call read_inopt(Mstar,'Mstar_1',db,errcount=nerr)
    if (isphere/=ibpwpoly) then
       call read_inopt(gamma,'gamma',db,errcount=nerr)
    endif
    if (isphere==iuniform) then
       call read_inopt(polyk,'polyk',db,errcount=nerr)
       input_polyk = .true.
    endif
    if (isphere==ievrard) then
       call read_inopt(ui_coef,'ui_coef',db,errcount=nerr)
    endif
    if (isphere==ihelmholtz) then
       call read_inopt(initialtemp,'initialtemp',db,errcount=nerr)
    endif
 endif
 !
 ! parse units
 !
 call select_unit(mass_unit,umass,nerr)
 if (nerr /= 0) then
    call error('setup_sphereinbox','mass unit not recognised')
    ierr = ierr + 1
 endif
 call select_unit(dist_unit,udist,nerr)
 if (nerr /= 0) then
    call error('setup_sphereinbox','length unit not recognised')
    ierr = ierr + 1
 endif

 if (nerr > 0) then
    print "(1x,a,i2,a)",'Setup_spheres: ',nerr,' error(s) during read of setup file'
    ierr = 1
 endif

 call close_db(db)

end subroutine read_setupfile
!-----------------------------------------------------------------------
end module
