!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! This module sets up accretion discs. The central object(s) can be
!   modelled with sink particles or external potentials. Systems with two
!   sink particles:
!     (i)  in a bound binary can have circumbinary, circumprimary, and
!          circumsecondary discs,
!     (ii) in an unbound binary (i.e. a fly-by) can have circumprimary and
!          circumsecondary discs.
!   In addition to gas, each disc can contain dust, modelled with either the
!   one fluid or two fluid methods. The dust can only grow in the two-fluid
!   method. Embedded planets can be added to single or circumbinary discs.
!
! :References: None
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters:
!   - R_rot          : *Set rotational velocity as Keplerian velocity at R=R_rot*
!   - Ratm_in        : *inner atmosphere radius (planet radii)*
!   - Ratm_out       : *outer atmosphere radius (planet radii)*
!   - Rin_sphere     : *Inner edge of sphere*
!   - Rout_sphere    : *Outer edge of sphere*
!   - T_floor        : *The minimum temperature in the simulation (for any locally isothermal EOS).*
!   - accr1          : *single star accretion radius*
!   - accr1a         : *single star accretion radius*
!   - accr1b         : *single star accretion radius*
!   - accr2          : *secondary accretion radius*
!   - accr2a         : *tight binary primary accretion radius*
!   - accr2b         : *tight binary secondary accretion radius*
!   - add_sphere     : *add sphere around disc?*
!   - add_turbulence : *Add turbulence to the sphere (0=no turbulence, 1=turbulence)*
!   - alphaSS        : *desired alphaSS (0 for minimal needed for shock capturing)*
!   - alpha_z        : *height of transition in tanh vertical temperature profile*
!   - atm_type       : *atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))*
!   - beta_z         : *variation in transition height over radius*
!   - bhspin         : *black hole spin*
!   - bhspinangle    : *black hole spin angle (deg)*
!   - deltat         : *output interval as fraction of orbital period*
!   - discstrat      : *stratify disc? (0=no,1=yes)*
!   - einst_prec     : *include Einstein precession*
!   - eos_file       : *Equation of state file for using lumdisc*
!   - ipotential     : *potential (1=central point mass,*
!   - istrat         : *temperature prescription (0=MAPS, 1=Dartois)*
!   - k              : *Scaling factor of Keplerian rotational velocity*
!   - lumdisc        : *Set qindex from stellar luminosity (ieos=24) (0=no 1=yes)*
!   - m1             : *first hierarchical level primary mass*
!   - m2             : *first hierarchical level secondary mass*
!   - mass_sphere    : *Mass of sphere*
!   - norbits        : *maximum number of orbits at outer disc*
!   - np             : *number of gas particles*
!   - nplanets       : *number of planets*
!   - nsinks         : *number of sinks*
!   - omega_cloud    : *Rotational velocity of the cloud (s^-1)*
!   - q1             : *tight binary 1 mass ratio*
!   - q2             : *tight binary 2 mass ratio*
!   - qatm           : *sound speed power law index of atmosphere*
!   - radkappa       : *constant radiation opacity kappa*
!   - ramp           : *Do you want to ramp up the planet mass slowly?*
!   - rho_core       : *planet core density (cgs units)*
!   - rms_mach       : *RMS Mach number of turbulence*
!   - set_freefall   : *Set the sphere in freefall (0=no freefall, 1=freefall)*
!   - subst          : *star to substitute*
!   - subst1         : *first star to substitute*
!   - subst2         : *second star to substitute*
!   - surface_force  : *model m1 as planet with surface*
!   - temp_atm0      : *atmosphere temperature scaling factor*
!   - temp_mid0      : *midplane temperature scaling factor*
!   - tfact          : *Scale the maximum length scale of the turbulence*
!   - use_mcfost     : *use the mcfost library*
!   - z0             : *z scaling factor*
!
! :Dependencies: centreofmass, datafiles, dim, dust, eos, eos_stamatellos,
!   extern_binary, extern_corotate, extern_lensethirring, externalforces,
!   fileutils, grids_for_setup, growth, infile_utils, io, io_control,
!   kernel, memory, options, orbits, part, partinject, physcon, prompting,
!   radiation_utils, set_dust, set_dust_options, setbinary, setdisc,
!   sethier_utils, sethierarchical, setorbit, setunits, shock_capturing,
!   spherical, systemutils, timestep, units, vectorutils, velfield
!
 use dim,              only:use_dust,maxalpha,use_dustgrowth,maxdusttypes,&
                            maxdustlarge,maxdustsmall,compiled_with_mcfost,gr
 use externalforces,   only:iext_star,iext_binary,iext_lensethirring,&
                            iext_einsteinprec,iext_corot_binary,iext_corotate,&
                            update_externalforce
 use extern_binary,    only:mass2,accradius1,accradius2,ramp,surface_force,eps_soft1
 use fileutils,        only:make_tags_unique
 use growth,           only:ifrag,isnow,rsnow,Tsnow,vfragSI,vfraginSI,vfragoutSI,gsizemincgs,iporosity
 use io,               only:master,warning,error,fatal
 use kernel,           only:hfact_default
 use options,          only:use_dustfrac,iexternalforce,use_hybrid,use_porosity
 use options,          only:use_mcfost,use_mcfost_stellar_parameters
 use part,             only:xyzmh_ptmass,maxvxyzu,vxyz_ptmass,ihacc,ihsoft,&
                            iJ2,ispinx,ispinz,iReff,igas,&
                            idust,iphase,dustprop,dustfrac,ndusttypes,ndustsmall,&
                            ndustlarge,grainsize,graindens,nptmass,iamtype,dustgasprop,&
                            VrelVf,filfac,probastick,rad,radprop,ikappa,iradxi
 use physcon,          only:au,solarm,jupiterm,earthm,pi,twopi,years,hours,deg_to_rad
 use setdisc,          only:scaled_sigma,get_disc_mass,maxbins,get_cs_from_lum
 use set_dust_options, only:set_dust_default_options,dust_method,dust_to_gas,&
                            ndusttypesinp,ndustlargeinp,ndustsmallinp,isetdust,&
                            dustbinfrac,check_dust_method,set_dust_grain_distribution
 use units,            only:umass,udist,utime
 use dim,              only:do_radiation
 use radiation_utils,  only:set_radiation_and_gas_temperature_equal
 use memory,           only:allocate_memory
 use setorbit,         only:orbit_t

 implicit none

 public  :: setpart

 private

 !--resolution
 integer :: np,np_dust(maxdustlarge)

 !--setup filename
 character(len=100) :: filename

 !--hierarchical configuration
 integer :: hl_index
 real :: current_mass, higher_mass
 integer :: higher_disc_index

 !--central objects
 real    :: mcentral
 real    :: m1,m2,m1a,m1b,m2a,m2b,q1,q2,accr1,accr2,accr1a,accr1b,accr2a,accr2b
 type(orbit_t) :: binary,binary1,binary2
 integer :: icentral,ipotential
 integer :: nsinks,subst,subst1,subst2
 real    :: bhspin,bhspinangle
 logical :: einst_prec

 !--stratification
 real    :: temp_atm0,temp_mid0
 real    :: z0_ref

 !--discs
 integer, parameter :: maxdiscs = 10
 real               :: discpos(3),discvel(3)

 character(len=20) :: disclabel
 character(len=*), dimension(maxdiscs), parameter :: disctype(1:4) = &
    (/'binary   ', &
      'primary  ', &
      'secondary', &
      'triple   '/)

 real    :: star_m(maxdiscs)
 real    :: totmass_gas
 real    :: J2star(maxdiscs),spin_period_star(maxdiscs),obliquity_star(maxdiscs)
 real    :: size_star(maxdiscs),kfac_star(maxdiscs)

 integer :: ndiscs
 integer :: onlydisc
 integer :: isetgas(maxdiscs)
 logical :: iuse_disc(maxdiscs)
 integer :: sigmaprofilegas(maxdiscs)
 logical :: sigma_file(maxdiscs)
 logical :: ismoothgas(maxdiscs)
 logical :: itapergas(maxdiscs)
 integer :: itapersetgas(maxdiscs)
 logical :: iwarp(maxdiscs)
 logical :: use_global_iso
 real    :: alphaSS

 real    :: R_in(maxdiscs),R_out(maxdiscs),R_ref(maxdiscs),R_c(maxdiscs)
 real    :: pindex(maxdiscs),disc_m(maxdiscs),sig_ref(maxdiscs),sig_norm(maxdiscs)
 real    :: T_bg,L_star(maxdiscs)
 real    :: qindex(maxdiscs),H_R(maxdiscs),T_floor
 real    :: posangl(maxdiscs),incl(maxdiscs)
 real    :: annulus_m(maxdiscs),R_inann(maxdiscs),R_outann(maxdiscs)
 real    :: R_warp(maxdiscs),H_warp(maxdiscs)
 real    :: Q_min(maxdiscs)

 integer :: sigmaprofiledust(maxdiscs,maxdusttypes)
 logical :: use_sigmadust_file(maxdiscs,maxdusttypes)
 logical :: ismoothdust(maxdiscs,maxdusttypes)
 logical :: itaperdust(maxdiscs,maxdusttypes)
 integer :: itapersetdust(maxdiscs,maxdusttypes)
 real    :: disc_mdust(maxdiscs,maxdusttypes),sig_normdust(maxdiscs,maxdusttypes)
 real    :: R_indust(maxdiscs,maxdusttypes),R_indust_swap(maxdiscs,maxdusttypes)
 real    :: R_outdust(maxdiscs,maxdusttypes),R_outdust_swap(maxdiscs,maxdusttypes)
 real    :: R_c_dust(maxdiscs,maxdusttypes)
 real    :: pindex_dust(maxdiscs,maxdusttypes),qindex_dust(maxdiscs,maxdusttypes)
 real    :: H_R_dust(maxdiscs,maxdusttypes)
 real :: enc_mass(maxbins,maxdiscs)
 real    :: e0(maxdiscs),eindex(maxdiscs),phiperi(maxdiscs)
 integer :: eccprofile(maxdiscs)
 logical :: iecc(maxdiscs)
 !--planets
 integer, parameter :: maxplanets = 9

 character(len=*), dimension(maxplanets), parameter :: num = &
    (/'1','2','3','4','5','6','7','8','9' /)

 logical :: istratify,lumdisc_logi
 integer :: nplanets,discstrat,lumdisc
 real    :: mplanet(maxplanets),rplanet(maxplanets)
 real    :: accrplanet(maxplanets),inclplan(maxplanets)
 real    :: eccplanet(maxplanets)
 real    :: Oplanet(maxplanets),wplanet(maxplanets),fplanet(maxplanets)
 real    :: J2planet(maxplanets),spin_period(maxplanets),obliquity(maxplanets)
 real    :: planet_size(maxplanets),kfac(maxplanets)
 real    :: period_planet_longest

 !--planetary atmosphere
 integer :: npart_planet_atm
 integer :: atm_type
 real    :: rho_core_cgs
 real    :: Ratm_in
 real    :: Ratm_out
 real    :: Natmfrac

 !--sphere of gas around disc
 logical :: add_sphere
 real :: Rin_sphere, Rout_sphere, mass_sphere
 real :: Kep_factor, R_rot, rms_mach, tfact, omega_cloud
 integer :: add_rotation, add_turbulence,set_freefall,dustfrac_method

 !--time
 real    :: tinitial
 real    :: deltat
 integer :: norbits

 !--other
 logical :: ichange_method
 real    :: iradkappa = huge(iradkappa)/10.
contains

!--------------------------------------------------------------------------
!+
! This is the only public subroutine of the module
!+
!--------------------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk
 real,              intent(out)   :: gamma
 real,              intent(out)   :: hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix

 write(*,"(/,65('-'),/,/,5x,a,/,/,65('-'))") 'Welcome to the New Disc Setup'

 !--set default options
 call set_default_options()!-1)

 !--set time
 time = tinitial

 !--get disc setup parameters from file or interactive setup
 call get_setup_parameters(id,fileprefix)

 !--allocate memory
 !nalloc = np
 !if (use_dust) nalloc = nalloc + sum(np_dust)
 !call allocate_memory(nalloc, part_only=.true.)

 !--compute number of discs based on setup options
 call number_of_discs()

 !--setup central object(s), i.e. sink particle(s) or potential
 call setup_central_objects(fileprefix)

 !--setup equation of state
 call equation_of_state(gamma)

 !--set surface density profile based on setup options
 call surface_density_profile()

 !--setup grain size distribution
 call set_dust_grain_distribution(ndusttypes,dustbinfrac,grainsize,graindens,udist,umass)

 !--compute disc mass and surface density normalization
 call calculate_disc_mass()

 !--setup disc(s)
 call setup_discs(id,fileprefix,hfact,gamma,npart,polyk,npartoftype,massoftype,xyzh,vxyzu)

 !--planet atmospheres
 call planet_atmosphere(id,npart,xyzh,vxyzu,npartoftype,gamma,hfact)

 !--initialise dustprop for dust particles only
 call initialise_dustprop(npart)

 !--check dust method for validity
 call check_dust_method(dust_method,ichange_method)

 !--print information about the angular momenta
 call print_angular_momentum(npart,xyzh,vxyzu)

 !--print dust information
 call print_dust()

 !--planets
 call set_planets(npart,massoftype,xyzh)

 !--set sphere of gas around disc
 if (add_sphere) call set_sphere_around_disc(id,npart,xyzh,vxyzu,npartoftype,massoftype,hfact)

 !--reset centre of mass to the origin
 if (any(iecc)) then !Means if eccentricity is present in .setup even if e0=0, it does not reset CM
    print*,'!!!!!!!!! Not resetting CM because one disc is eccentric: CM and ellipse focus do not match !!!!!!!!!'!,e0>0
 else
    call set_centreofmass(npart,xyzh,vxyzu)
 endif

 !--set tmax and dtmax
 call set_tmax_dtmax(fileprefix)

 if (do_radiation) then
    call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)
    radprop(ikappa,1:npart) = iradkappa
 endif

 !--remind user to check for warnings and errors
 write(*,20)
20 format(/, &
   "-----------------------------------------------------------------",/, &
   "",/, &
   "     Check output for warnings and errors",/, &
   "",/, &
   "-----------------------------------------------------------------",/)

 return

end subroutine setpart

!--------------------------------------------------------------------------
!
! Set default options
!
!--------------------------------------------------------------------------
subroutine set_default_options()
 use sethierarchical, only:set_hierarchical_default_options
 use systemutils,     only:get_command_option
 use setorbit,        only:set_defaults_orbit
 use setunits,        only:dist_unit,mass_unit
 use sethier_utils,   only:findloc_local
 integer :: i

 !--time
 tinitial = 0.

 !--units
 dist_unit = 'au'
 mass_unit = 'solarm'

 !--central object(s)
 icentral = 1

 !--external potential
 ipotential = 1

 !--point mass
 iexternalforce = iext_star
 m1    = 1.
 m2    = 1.
 accr1 = 1.
 accr2 = 1.

 !--oblateness of main objects
 J2star = 0.
 spin_period_star = 10.
 obliquity_star = 0.
 size_star = 1.
 kfac_star = 0.205

 !--planetary atmosphere
 surface_force = .false.

 !--sphere around disc
 Rin_sphere = 30.
 Rout_sphere = 500.
 mass_sphere = 0.1
 add_rotation = 0
 Kep_factor = 0.08
 R_rot = 150.
 add_turbulence = 0
 omega_cloud = 1e-11
 dustfrac_method = 0
 set_freefall = 0
 rms_mach = 1.
 tfact = 0.1

 !--spinning black hole (Lense-Thirring)
 einst_prec = .false.
 bhspin      = 1.
 bhspinangle = 0.

 !--sink particle(s)
 nsinks = 1

 !--binary
 call set_defaults_orbit(binary)
 call set_defaults_orbit(binary1)
 call set_defaults_orbit(binary2)

 !--hierarchical
 call set_hierarchical_default_options()!id)

 !--multiple disc options
 iuse_disc = .false.
 iuse_disc(1) = .true.
 ndiscs = 1

 !--eos
 use_global_iso = .false.

 !--dust distribution
 call set_dust_default_options()

 !--gas disc
 R_in         = 1.
 R_out        = 150.
 R_ref        = 10.
 R_c          = 150.
 R_warp       = 0.
 H_warp       = 0.
 isetgas      = 0
 itapergas    = .false.
 itapersetgas = 0
 ismoothgas   = .true.
 iwarp        = .false.
 pindex       = 1.
 qindex       = 0.25
 alphaSS      = 0.005
 posangl      = 0.
 incl         = 0.
 H_R          = 0.05
 disc_m       = 0.05
 sig_norm     = 1.e-02
 sig_ref      = 1.e-02
 Q_min        = 1.0
 annulus_m    = 0.05
 R_inann      = 1.
 R_outann     = 150.
 lumdisc      = 0
 L_star(:)    = 1.
 T_bg         = 5.
 lumdisc_logi = .false.

 !--floor temperature
 T_floor      = 0.0

 !--disc eccentricity
 eccprofile=0
 e0=0.
 eindex=0.
 phiperi=0.

 !--dust disc
 R_indust           = 1.
 R_outdust          = 150.
 pindex_dust        = 1.
 qindex_dust        = 0.25
 H_R_dust           = 0.05
 use_sigmadust_file = .false.
 itaperdust         = .false.
 itapersetdust      = 0
 ismoothdust        = .true.
 R_c_dust           = 150.

 !--dust growth
 ifrag = 1
 isnow = 0
 rsnow = 100.
 Tsnow = 150.
 vfragSI = 15.
 vfraginSI = 5.
 vfragoutSI = 15.
 gsizemincgs = 5.e-3

 !--resolution, default is 1000000 but can be set with --np=N
 !  command line option (used to keep the test suite running fast)
 np = int(get_command_option('np',default=1000000))
 np_dust = np/maxdustlarge/5

 !--planets
 nplanets      = 0
 mplanet       = 1.
 rplanet       = (/ (10.*i, i=1,maxplanets) /)
 accrplanet    = 0.25
 inclplan      = 0.
 J2planet      = 0.
 spin_period   = 10.
 obliquity     = 0.
 planet_size   = 1.
 kfac          = 0.205
 eccplanet     = 0.
 wplanet       = 270.
 Oplanet       = 0.
 fplanet       = 180.

 !--stratification
 istratify     = .false.
 discstrat     = 0
 temp_mid0     = 24.
 temp_atm0     = 63.
 z0_ref        = 9.

 !--simulation time
 deltat  = 0.1
 norbits = 100

 !--planetary atmospheres
 atm_type     = 1
 rho_core_cgs = 5.
 Ratm_in      = 1.
 Ratm_out     = 3.
 Natmfrac     = 0.

end subroutine set_default_options

!--------------------------------------------------------------------------
!
! Get setup parameters from interactive setup or file
!
!--------------------------------------------------------------------------
subroutine get_setup_parameters(id,fileprefix)
 integer,           intent(in) :: id
 character(len=20), intent(in) :: fileprefix
 logical :: iexist,seq_exists
 integer :: j,ierr

 filename=trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) then
    !--read from setup file
    call read_setupfile(filename,ierr)
    if (id==master) call write_setupfile(filename)
    if (ierr /= 0) then
       stop
    endif

    !--setup multiple discs each with different orbital parameters
    !  specified in orbits.dat
    inquire(file='orbits.dat',exist=iexist)
    inquire(file=trim(fileprefix)//'A.setup',exist=seq_exists)
    if (iexist .and. .not.seq_exists) then
       open(unit=23,file='orbits.dat',status='old',iostat=ierr)
       j = 0
       do while(ierr==0)
          binary%input_type = 0
          read(23,*,iostat=ierr) binary%a,binary%e,binary%i,binary%O,binary%w,binary%f
          write(binary%elems%a,"(g0)") binary%a
          binary%elems%a = trim(adjustl(binary%elems%a)) ! convert to string
          if (ierr==0) then
             j = j + 1
             write(filename,"(a)") trim(fileprefix)//achar(j+64)//'.setup'
             if (id==master) call write_setupfile(filename)
          endif
       enddo
       close(unit=23)
       stop
    endif

 elseif (id==master) then

    !--interactive setup
    print "(a,/)",' '//trim(filename)//' not found: using interactive setup'
    call setup_interactive(id)

    !--write setup file from interactive setup
    call write_setupfile(filename)
    print "(/,a)",' >>> please edit '//trim(filename)//' to set parameters for your problem then rerun phantomsetup <<<'
    stop

 else

    stop

 endif

end subroutine get_setup_parameters

!--------------------------------------------------------------------------
!
! Calculate the number of required discs
!
!--------------------------------------------------------------------------
subroutine number_of_discs()
 integer :: i

 ndiscs = max(count(iuse_disc),1)
 !--index of disc (if only one)
 onlydisc = 0
 if (ndiscs==1) then
    do i=1,maxdiscs
       if (iuse_disc(i)) onlydisc = i
    enddo
 endif

end subroutine number_of_discs

!--------------------------------------------------------------------------
!
! Set the equation of state
!
!--------------------------------------------------------------------------
subroutine equation_of_state(gamma)
 use eos,             only:isink,qfacdisc,qfacdisc2,polyk2,beta_z,z0,cs_min,gmw,&
                           ieos,icooling,ipdv_heating,ishock_heating
 use io_control,      only:nfulldump
 use shock_capturing, only:alphau
 use eos_stamatellos, only:init_coolra,read_optab,eos_file
 use physcon, only:rpiontwo,mass_proton_cgs,kboltz
 use units,   only:unit_velocity
 use sethier_utils,   only:findloc_local
 real, intent(out) :: gamma
 real              :: H_R_atm, cs

 logical :: is_isothermal
 integer :: i,ierr

 is_isothermal = (maxvxyzu==3)

 if (compiled_with_mcfost) then
    if (use_mcfost) then
       is_isothermal = .false.
       nfulldump = 1
    else
       is_isothermal = .true.
    endif
 endif

 if (is_isothermal) then

    !--isothermal
    gamma = 1.0
    if (ndiscs /= 1) then
       !--multiple discs
       if (use_global_iso) then
          !--globally isothermal
          ieos = 1
          qindex = 0.
          qfacdisc = qindex(1)
          print "(/,a)",' setting ieos=1 for globally isothermal disc'
          if (iuse_disc(1)) then
             H_R(2) = sqrt(R_ref(2)/R_ref(1)*(m1+m2)/m1) * H_R(1)
             H_R(3) = sqrt(R_ref(3)/R_ref(1)*(m1+m2)/m2) * H_R(1)
             call warning('setup_disc','using circumbinary (H/R)_ref to set global temperature')
          elseif (iuse_disc(2)) then
             H_R(3) = sqrt(R_ref(3)/R_ref(2)*m1/m2) * H_R(2)
             call warning('setup_disc','using circumprimary (H/R)_ref to set global temperature')
          endif
       else
          !--locally isothermal prescription from Farris et al. (2014) for binary system
          if (nsinks>4) then
             ieos = 13
             print "(/,a)",' setting ieos=13 for locally isothermal from generalised Farris et al. (2014) prescription'
             higher_disc_index = findloc_local(iuse_disc, .true.)
             qfacdisc = qindex(higher_disc_index)
             call get_hier_disc_label(higher_disc_index, disclabel)

             call warning('setup_disc','using circum-'//trim(disclabel)//' (H/R)_ref to set global temperature')
          else
             ieos = 14
             print "(/,a)",' setting ieos=14 for locally isothermal from Farris et al. (2014)'
             if (iuse_disc(1)) then
                qfacdisc = qindex(1)
                call warning('setup_disc','using circumbinary (H/R)_ref to set global temperature')
             elseif (iuse_disc(2)) then
                qfacdisc = qindex(2)
                call warning('setup_disc','using circumprimary (H/R)_ref to set global temperature')
             endif
          endif
       endif
    else
       !--single disc
       if (qindex(onlydisc)>= 0.) then
          do i=1,maxdiscs
             !--eos around sink
             if (iuse_disc(i)) isink = i-1
          enddo
          !--locally isothermal
          if (isink /= 0 .and. isink /= 3) then ! isink == 3 special case, to be generalised
             ieos = 6
             print "(/,a)",' setting ieos=6 for locally isothermal disc around sink'
          else
             isink = 0
             if (discstrat > 0) then
                ieos = 7
                print "(/,a)",' setting ieos=7 for locally isothermal disc with stratification'
                call temp_to_HR(temp_mid0,H_R(onlydisc),R_ref(onlydisc),mcentral,cs)
                call temp_to_HR(temp_atm0,H_R_atm,R_ref(onlydisc),mcentral,cs)
                polyk2 = (cs*(1./R_ref(onlydisc))**(-qfacdisc2))**2
                z0 = z0_ref/R_ref(onlydisc)**beta_z
             else
                if (ieos == 6) then
                   ! handle the case where ieos=6 is already set in the .in file; do not override this
                   isink = 1
                   print "(/,a)",' keeping ieos=6 for locally isothermal disc with bright primary'
                else
                   ieos = 3
                   print "(/,a)",' setting ieos=3 for locally isothermal disc around origin'
                endif
             endif
          endif
          qfacdisc = qindex(onlydisc)
       else
          call warning('setup_disc','setting qindex < 0.')
       endif
    endif

 else
    !-- adiabatic
    if (lumdisc > 0) then
       !--for radapprox cooling
       print "(/,a)", ' setting ieos=24 and icooling=9 for radiative cooling approximation'
       lumdisc_logi = .True.
       ieos = 24
       icooling = 9
       gamma = 5./3. ! in case it's needed
       call init_coolra()
       call read_optab(eos_file,ierr)
       if (ndiscs > 1) then
          print *, "We can't set up multiple radapprox discs yet :,("
          stop
       else
          cs = get_cs_from_lum(L_star(1),R_ref(1),T_bg,gamma) / rpiontwo
          H_R(1) = cs * R_ref(1)**0.5 / sqrt(m1) ! single central star, G=1
       endif
    else
       !--adiabatic
       ieos = 2
       gamma = 5./3.
       icooling = 3
    endif

    if (use_mcfost) then
       icooling = 0
       ipdv_heating = 0
       ishock_heating = 0
       alphau = 0
    endif

 endif

 if ( any( ieos==(/3,6,7,13,14/) ) ) then
    print "(/,a)",' Setting floor temperature to ', T_floor, ' K.'
    cs_min =  gmw*T_floor/(mass_proton_cgs/kboltz * unit_velocity**2)
 endif

end subroutine equation_of_state

!--------------------------------------------------------------------------
!
! Set the surface density profile choice
!
!    0 = power law
!    1 = exponentially tapered power law
!    2 = smoothed power law
!    3 = both tapered and smoothed
!    4 = alternative taper
!    5 = alternative taper with smoothing
!
!--------------------------------------------------------------------------
subroutine surface_density_profile()
 integer :: i,j

 !--gas profile
 do i=1,maxdiscs
    sigmaprofilegas(i) = 0
    if (itapergas(i)) sigmaprofilegas(i) = 1
    if (ismoothgas(i)) sigmaprofilegas(i) = 2
    if (itapergas(i) .and. ismoothgas(i)) sigmaprofilegas(i) = 3
    if (itapersetgas(i)==1) sigmaprofilegas(i) = 4
    if (itapersetgas(i)==1 .and. ismoothgas(i)) sigmaprofilegas(i) = 5
    if (sigma_file(i)) sigmaprofilegas(i) = 6
 enddo

 !--dust profile
 if (use_dust) then
    do i=1,maxdiscs
       do j=1,ndusttypes
          sigmaprofiledust(i,j) = 0
          if (itaperdust(i,j)) sigmaprofiledust(i,j) = 1
          if (ismoothdust(i,j)) sigmaprofiledust(i,j) = 2
          if (itaperdust(i,j) .and. ismoothdust(i,j)) sigmaprofiledust(i,j) = 3
          if (itapersetdust(i,j)==1) sigmaprofiledust(i,j) = 4
          if (itapersetdust(i,j)==1 .and. ismoothdust(i,j)) sigmaprofiledust(i,j) = 5
          if (use_sigmadust_file(i,j)) sigmaprofiledust(i,j) = 6
       enddo
       !--swap radii to keep dust profile the same as gas within [R_indust,R_outdust]
       if (isetdust == 2) then
          R_indust_swap(i,:)  = R_in(i)
          R_outdust_swap(i,:) = R_out(i)
       else
          R_indust_swap(i,:)  = R_indust(i,:)
          R_outdust_swap(i,:) = R_outdust(i,:)
       endif
    enddo
 endif

end subroutine surface_density_profile

!--------------------------------------------------------------------------
!
! Set up the central object(s)
!
!--------------------------------------------------------------------------
subroutine setup_central_objects(fileprefix)
 use externalforces,       only:mass1,accradius1
 use extern_lensethirring, only:blackhole_spin,blackhole_spin_angle
 use setbinary,            only:set_binary
 use orbits,               only:convert_flyby_to_elements
 use sethierarchical,      only:set_hierarchical,set_multiple
 use setorbit,             only:set_orbit
 use setunits,             only:dist_unit,mass_unit
 character(len=20), intent(in) :: fileprefix

 integer :: i,ierr

 mcentral = m1
 select case (icentral)
 case (0)
    select case (ipotential)
    case (1)
       print "(/,a)",' Central point mass represented by external force with accretion boundary'
       print "(a,g10.3,a)",'   Object mass:      ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   Accretion Radius: ', accr1, trim(dist_unit)
       mass1      = m1
       accradius1 = accr1
       mcentral   = m1
    case (2)
       print "(/,a)",' Central binary represented by external force with accretion boundary'
       print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
       print "(a,g10.3)",  '   Binary mass ratio:  ', m2/m1
       print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
       print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)
       mass1       = m1
       mass2       = m2
       accradius1  = accr1
       accradius2  = accr2
       if (iexternalforce == iext_corot_binary) then
          mcentral = m1
       else
          mcentral = m1 + m2
       endif
    case (3)
       print "(/,a)",' Central black hole represented by external force with accretion boundary'
       print "(a,g10.3,a)",'   Black hole mass:        ', m1,    trim(mass_unit)
       print "(a,g10.3)",  '   Accretion Radius:       ', accr1
       print "(a,g10.3)",  '   Black hole spin:        ', bhspin
       print "(a,g10.3,a)",'   Black hole spin angle:  ', bhspinangle, 'deg'
       mass1                = m1
       accradius1           = accr1
       blackhole_spin       = bhspin
       blackhole_spin_angle = bhspinangle*deg_to_rad
       mcentral             = m1
    end select
    call update_externalforce(iexternalforce,tinitial,0.)
 case (1)
    select case (nsinks)
    case (1)
       !--single star
       print "(/,a)",' Central object represented by a sink at the system origin'
       print "(a,g10.3,a)",'   Object mass:      ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   Accretion Radius: ', accr1, trim(dist_unit)
       nptmass                      = 1
       xyzmh_ptmass(:,:)            = 0.
       xyzmh_ptmass(1:3,nptmass)    = 0.
       xyzmh_ptmass(4,nptmass)      = m1
       xyzmh_ptmass(ihacc,nptmass)  = accr1
       xyzmh_ptmass(ihsoft,nptmass) = 0.
       vxyz_ptmass                  = 0.
       mcentral                     = m1
       discpos                      = 0.
       discvel                      = 0.
    case (2)
       !--binary
       select case (binary%input_type)
       case (0)
          !--bound
          print "(/,a)",' Central objects represented by two sinks'
          print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
          print "(a,g10.3)",  '   Binary mass ratio:  ', m2/m1
          mcentral = m1 + m2
       case default
          !--unbound (flyby)
          print "(/,a)",' Central object represented by a sink at the system origin with a perturber sink'
          print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
          print "(a,g10.3,a)",'   Perturber mass:     ', m2,    trim(mass_unit)
          mcentral = m1
       end select
       print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
       print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)

       nptmass  = 0
       call set_orbit(binary,m1,m2,accr1,accr2,xyzmh_ptmass,vxyz_ptmass,nptmass,verbose=.true.,ierr=ierr)

       discpos = 0.
       discvel = 0.

    case (5:)
       call set_hierarchical(fileprefix, nptmass, xyzmh_ptmass, vxyz_ptmass, ierr)

    case (3)
       !-- hierarchical triple
       nptmass  = 0
       print "(/,a)",' Central objects represented by a hierarchical triple'
       print "(a,g10.3,a)",'     First hierarchical level primary mass: ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   First hierarchical level secondary mass: ', m2,    trim(mass_unit)
       print "(a,g10.3)",  '                    Wide binary mass ratio: ', m2/m1
       print "(a,g10.3)",  '                   Tight binary mass ratio: ', q2
       print "(a,g10.3)",  '                    Star to be substituted: ', abs(subst)
       print "(a,g10.3,a)",'                        Accretion Radius 1: ', accr1, trim(dist_unit)
       print "(a,g10.3,a)",'                       Accretion Radius 2a: ', accr2a, trim(dist_unit)
       print "(a,g10.3,a)",'                       Accretion Radius 2b: ', accr2b, trim(dist_unit)

       if (subst>0) then
          print "(a,g10.3,a)",'      Tight binary orientation referred to: substituted star orbital plane'
       else
          print "(a,g10.3,a)",'      Tight binary orientation referred to: sky'
       endif

       call set_multiple(m1,m2,semimajoraxis=binary%a,eccentricity=binary%e, &
            posang_ascnode=binary%O,arg_peri=binary%w,incl=binary%i, &
            f=binary%f,accretion_radius1=accr1,accretion_radius2=accr1, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)

       if (iuse_disc(1)) then
          discpos=xyzmh_ptmass(1:3,2)
          discvel=vxyz_ptmass(1:3,2)
       else
          discpos = 0.
          discvel = 0.
       endif

       call set_multiple(m2/(q2+1),m2*q2/(q2+1),semimajoraxis=binary2%a,eccentricity=binary2%e, &
            posang_ascnode=binary2%O,arg_peri=binary2%w,incl=binary2%i, &
            f=binary2%f,accretion_radius1=accr2a,accretion_radius2=accr2b, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass, subst=subst,ierr=ierr)

       mcentral = m2
    case(4)
       !-- hierarchical quadruple
       nptmass  = 0
       print "(/,a)",' Central objects represented by a hierarchical quadruple'
       print "(a,g10.3,a)",'     First hierarchical level primary mass: ', m1,    trim(mass_unit)
       print "(a,g10.3,a)",'   First hierarchical level secondary mass: ', m2,    trim(mass_unit)
       print "(a,g10.3)",  '                    Wide binary mass ratio: ', m2/m1
       print "(a,g10.3)",  '                 Tight binary 1 mass ratio: ', q1
       print "(a,g10.3)",  '                 Tight binary 2 mass ratio: ', q2
       print "(a,g10.3)",  '              First star to be substituted: ', abs(subst1)
       print "(a,g10.3)",  '             Second star to be substituted: ', abs(subst2)
       print "(a,g10.3,a)",'                        Accretion Radius 1: ', accr1, trim(dist_unit)
       print "(a,g10.3,a)",'                       Accretion Radius 1a: ', accr1a, trim(dist_unit)
       print "(a,g10.3,a)",'                       Accretion Radius 1b: ', accr1b, trim(dist_unit)
       print "(a,g10.3,a)",'                       Accretion Radius 2a: ', accr2a, trim(dist_unit)
       print "(a,g10.3,a)",'                       Accretion Radius 2b: ', accr2b, trim(dist_unit)

       if (subst1>0) then
          print "(a,g10.3,a)",'      Tight binary 1 orientation referred to: substituted star orbital plane'
       else
          print "(a,g10.3,a)",'      Tight binary 1 orientation referred to: sky'
       endif
       if (subst2>0) then
          print "(a,g10.3,a)",'      Tight binary 2 orientation referred to: substituted star orbital plane'
       else
          print "(a,g10.3,a)",'      Tight binary 2 orientation referred to: sky'
       endif

       call set_multiple(m1,m2,semimajoraxis=binary%a,eccentricity=binary%e, &
            posang_ascnode=binary%O,arg_peri=binary%w,incl=binary%i, &
            f=binary%f,accretion_radius1=accr1,accretion_radius2=accr1, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)

       discpos=xyzmh_ptmass(1:3,2)
       discvel=vxyz_ptmass(1:3,2)

       call set_multiple(m1/(q1+1),m1*q1/(q1+1),semimajoraxis=binary1%a,eccentricity=binary1%e, &
         posang_ascnode=binary1%O,arg_peri=binary1%w,incl=binary1%i, &
         f=binary1%f,accretion_radius1=accr1a,accretion_radius2=accr1b, &
         xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,subst=subst1,ierr=ierr)

       call set_multiple(m2/(q2+1),m2*q2/(q2+1),semimajoraxis=binary2%a,eccentricity=binary2%e, &
            posang_ascnode=binary2%O,arg_peri=binary2%w,incl=binary2%i, &
            f=binary2%f,accretion_radius1=accr2a,accretion_radius2=accr2b, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,subst=subst2,ierr=ierr)

       mcentral = m2

    end select
 end select

 !--set array of central object masses
 star_m(1:4) = (/mcentral, m1, m2, m1+m2/)
 do i=1,maxdiscs
    if (.not.iuse_disc(i)) star_m(i) = 0.
 enddo

 do i=1,nsinks
    if (abs(J2star(i)) > 0.) then
       call set_sink_oblateness(i,J2star(i),size_star(i),spin_period_star(i),kfac_star(i),obliquity_star(i))
       print "(a)",'# oblateness for object '//num(i)
       call print_oblateness_info(i,spin_period_star(i))
    endif
 enddo

end subroutine setup_central_objects

!--------------------------------------------------------------------------
!
! Calculate the required disc masses
!
!--------------------------------------------------------------------------
subroutine calculate_disc_mass()
 use grids_for_setup, only:init_grid_sigma,deallocate_sigma
 integer :: i,j
 real :: enc_m(maxbins),rad(maxbins)
 real :: Q_mintmp,disc_mtmp,annulus_mtmp
 real :: rgrid_min,rgrid_max,fac

 totmass_gas  = 0.
 disc_mdust = 0.

 do i=1,maxdiscs
    !--initialise the sigma grid file
    if (sigmaprofilegas(i)==6) call init_grid_sigma(R_in(i),R_out(i))

    if (iuse_disc(i)) then
       !
       !--set up a common radial grid for the enclosed mass including gas and dust
       !  even if the gas/dust discs have different radial extents
       !
       rgrid_min = R_in(i)
       rgrid_max = R_out(i)
       if (isetgas(i)==1) then
          rgrid_min = min(rgrid_min,R_inann(i))
          rgrid_max = max(rgrid_max,R_outann(i))
       endif
       if (use_dust) then
          rgrid_min = min(rgrid_min,minval(R_indust_swap(i,1:ndusttypes)))
          rgrid_max = min(rgrid_max,maxval(R_outdust_swap(i,1:ndusttypes)))
       endif
       do j=1,maxbins
          rad(j) = rgrid_min + (j-1) * (rgrid_max-rgrid_min)/real(maxbins-1)
       enddo
       !--gas discs
       select case(isetgas(i))
       case (0)
          !--set sigma normalisation from disc mass
          sig_norm(i) = 1.d0
          call get_disc_mass(disc_mtmp,enc_m,rad,Q_mintmp,sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
          fac = disc_m(i) / disc_mtmp
          sig_norm(i) = sig_norm(i) * fac
          enc_m = enc_m * fac

       case (1)
          !--set disc mass from annulus mass
          sig_norm(i) = 1.d0
          call get_disc_mass(annulus_mtmp,enc_m,rad,Q_mintmp,sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_inann(i),R_outann(i),R_ref(i),R_c(i), &
                             H_R(i))
          sig_norm(i) = sig_norm(i) * annulus_m(i) / annulus_mtmp
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       case (2)
          !--set disc mass from sigma normalisation
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       case (3)
          !--set disc mass from sigma at reference radius
          if (.not.(R_in(i) < R_ref(i)) .and. ismoothgas(i)) call fatal('setup_disc', &
             'if smoothing inner disc and setting disc mass by sigma(R_ref), require R_in < R_ref')
          sig_norm(i) = sig_ref(i) / scaled_sigma(R_ref(i),sigmaprofilegas(i),pindex(i),R_ref(i),R_in(i),R_out(i),R_c(i))
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       case (4)
          !--set disc mass from minimum Toomre Q
          sig_norm(i) = 1.d0
          call get_disc_mass(disc_mtmp,enc_m,rad,Q_mintmp,sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
          fac = Q_mintmp / Q_min(i)
          sig_norm(i) = sig_norm(i) * fac
          !--recompute actual disc mass and Toomre Q
          call get_disc_mass(disc_m(i),enc_m,rad,Q_min(i),sigmaprofilegas(i),sig_norm(i), &
                             star_m(i),pindex(i),qindex(i),R_in(i),R_out(i),R_ref(i),R_c(i), &
                             H_R(i))
       end select

       totmass_gas = totmass_gas + disc_m(i)
       enc_mass(:,i) = enc_m + star_m(i)

       !--dust discs
       if (use_dust) then
          disc_mdust(i,:) = 0.
          do j=1,ndusttypes
             disc_mdust(i,j) = disc_m(i) * dust_to_gas * dustbinfrac(j)
             sig_normdust(i,j) = 1.d0
             call get_disc_mass(disc_mtmp,enc_m,rad,Q_mintmp,sigmaprofiledust(i,j), &
                                sig_normdust(i,j),star_m(i),pindex_dust(i,j),qindex_dust(i,j), &
                                R_indust_swap(i,j),R_outdust_swap(i,j),R_ref(i),R_c_dust(i,j),H_R_dust(i,j))
             fac = disc_mdust(i,j) / disc_mtmp
             sig_normdust(i,j) = sig_normdust(i,j) * fac
             enc_mass(:,i) = enc_mass(:,i) + enc_m(:)*fac
          enddo
       endif
    endif
    !--deallocating datasigma if density profile is read from file,
    !--will be re-initialised with right Rin/out if needed
    if (sigmaprofilegas(i)==6) call deallocate_sigma()
 enddo

end subroutine calculate_disc_mass

!--------------------------------------------------------------------------
!
! Set up the discs
!
!--------------------------------------------------------------------------
subroutine setup_discs(id,fileprefix,hfact,gamma,npart,polyk,&
                       npartoftype,massoftype,xyzh,vxyzu)
 use options,         only:alpha
 use setbinary,       only:Rochelobe_estimate
 use sethierarchical, only:get_hierarchical_level_com,get_hier_level_mass,hs
 use sethier_utils,   only:findloc_local
 use setdisc,         only:set_disc
 use growth,          only:alpha_dg
 integer,           intent(in)    :: id
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: hfact
 real,              intent(in)    :: gamma
 integer,           intent(out)   :: npart
 real,              intent(out)   :: polyk
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(inout) :: xyzh(:,:)
 real,              intent(inout) :: vxyzu(:,:)

 integer            :: i,j,itype
 integer            :: npingasdisc,npindustdisc
 integer            :: iprofilegas,iprofiledust
 real               :: Rochelobe
 real               :: polyk_dust
 real               :: xorigini(3),vorigini(3)
 real               :: alpha_returned(maxdiscs)
 character(len=100) :: prefix
 character(len=100) :: dustprefix(maxdusttypes)

 hfact = hfact_default
 incl    = incl*deg_to_rad
 posangl = posangl*deg_to_rad
 alpha = alphaSS
 alpha_dg = alphaSS
 npart = 0
 npartoftype(:) = 0

 do i=1,maxdiscs
    if (iuse_disc(i)) then
       if (ndiscs > 1) then
          if (nsinks<4) then
             print "(/,a)",'>>> Setting up circum'//trim(disctype(i))//' disc <<<'
             prefix = trim(fileprefix)//'-'//disctype(i)
          else
             call get_hier_disc_label(i, disclabel)
             print "(/,a)",'>>> Setting up circum-'//trim(disclabel)//' disc <<<'
             prefix = trim(fileprefix)//'-'//trim(disclabel)
          endif
       else
          prefix = fileprefix
       endif

       if (nsinks < 4) then
          !--set disc origin
          select case(i)
          case(3)
             !--circumsecondary
             xorigini  = xyzmh_ptmass(1:3,2)
             vorigini  = vxyz_ptmass(1:3,2)
             Rochelobe = Rochelobe_estimate(m1,m2,binary%a)
          case(2)
             !--circumprimary
             xorigini  = xyzmh_ptmass(1:3,1)
             vorigini  = vxyz_ptmass(1:3,1)
             Rochelobe = Rochelobe_estimate(m2,m1,binary%a)
          case default
             !--single disc or circumbinary or circumtriple
             !  centre of mass of binary defined to be zero (see set_binary)
             xorigini  = discpos
             vorigini  = discvel
             Rochelobe = huge(0.)
          end select
       else
          call get_hier_disc_label(i, disclabel)

          m2 = get_hier_level_mass(disclabel)

          if (len(trim(disclabel))>1) then
             m1 = get_hier_level_mass(disclabel(:len(trim(disclabel))-1))-m2

             hl_index = findloc_local(hs%labels%hl, disclabel(:len(trim(disclabel))-1))
             if (hl_index == 0) call fatal('setup_disc','disc level not found in hierarchy')

             Rochelobe = Rochelobe_estimate(m1,m2,hs%levels(hl_index)%a)
          else
             Rochelobe = huge(0.)
          endif

          star_m(i) = m2

          call get_hierarchical_level_com(disclabel, xorigini, vorigini, xyzmh_ptmass, vxyz_ptmass, fileprefix)

       endif

       if ((ndiscs > 1 .and. binary%input_type==0) .and. (R_out(i) > Rochelobe)) then
          call warning('setup_disc', &
             'Outer disc radius for circum-'//trim(disctype(i))//' > Roche lobe of ' &
             //trim(disctype(i)))
       endif

       !--taper gas disc
       iprofilegas = 0
       if (itapergas(i)) iprofilegas = 1
       if (itapersetgas(i) == 1) iprofilegas = 2
       if (sigma_file(i)) iprofilegas = 3

       !--set disc(s)
       if (use_dust .and. use_dustfrac) then

          !--taper dust disc
          iprofiledust = 0
          if (itaperdust(i,1)) iprofiledust = 1
          if (itapersetdust(i,1) == 1) iprofiledust = 2
          if (use_sigmadust_file(i,1)) iprofiledust = 3

          !--gas and dust mixture disc
          npingasdisc = int(disc_m(i)/totmass_gas*np)

          call set_disc(id,master        = master,               &
                        mixture          = .true.,               &
                        npart            = npingasdisc,          &
                        npart_start      = npart + 1,            &
                        rref             = R_ref(i),             &
                        rmin             = R_in(i),              &
                        rmax             = R_out(i),             &
                        rmindust         = R_indust(i,1),        &
                        rmaxdust         = R_outdust(i,1),       &
                        indexprofile     = iprofilegas,          &
                        indexprofiledust = iprofiledust,         &
                        rc               = R_c(i),               &
                        rcdust           = R_c_dust(i,1),        &
                        p_index          = pindex(i),            &
                        p_indexdust      = pindex_dust(i,1),     &
                        q_index          = qindex(i),            &
                        q_indexdust      = qindex_dust(i,1),     &
                        HoverR           = H_R(i),               &
                        HoverRdust       = H_R_dust(i,1),        &
                        disc_mass        = disc_m(i),            &
                        disc_massdust    = sum(disc_mdust(i,:)), &
                        star_mass        = star_m(i),            &
                        gamma            = gamma,                &
                        particle_mass    = massoftype(igas),     &
                        xyz_origin       = xorigini,             &
                        vxyz_origin      = vorigini,             &
                        hfact            = hfact,                &
                        xyzh             = xyzh,                 &
                        vxyzu            = vxyzu,                &
                        polyk            = polyk,                &
                        alpha            = alpha,                &
                        ismooth          = ismoothgas(i),        &
                        position_angle   = posangl(i),           &
                        inclination      = incl(i),              &
                        rwarp            = R_warp(i),            &
                        warp_smoothl     = H_warp(i),            &
                        e0               = e0(i),                &
                        eindex           = eindex(i),            &
                        phiperi          = phiperi(i),           &
                        eccprofile       = eccprofile(i),        &
                        bh_spin          = bhspin,               &
                        enc_mass         = enc_mass(:,i),        &
                        prefix           = prefix,               &
                        lumdisc          = lumdisc_logi,         &
                        L_star           = L_star(1),            &
                        T_bg             = T_bg)

          !--set dustfrac
          call set_dustfrac(i,npart+1,npart+1+npingasdisc,xyzh,xorigini)

          npart = npart + npingasdisc
          npartoftype(igas) = npartoftype(igas) + npingasdisc

          if (use_hybrid) then
             if (ndustlarge > 1) then
                dustprefix = trim(prefix)//'-'
                call make_tags_unique(ndustlarge,dustprefix)
             else
                dustprefix = trim(prefix)
             endif

             !--dust disc(s)
             do j=1,ndustlarge

                npindustdisc = int(disc_mdust(i,j)/sum(disc_mdust(:,j))*np_dust(j))
                itype = idust + j - 1

                !--taper dust disc
                iprofiledust = 0
                if (itaperdust(i,j)) iprofiledust = 1
                if (itapersetdust(i,j) == 1) iprofiledust = 2
                if (use_sigmadust_file(i,j)) iprofiledust = 3

                call set_disc(id,master      = master,             &
                              npart          = npindustdisc,       &
                              npart_start    = npart + 1,          &
                              particle_type  = itype,              &
                              rref           = R_ref(i),           &
                              rmin           = R_indust(i,j),      &
                              rmax           = R_outdust(i,j),     &
                              indexprofile   = iprofiledust,       &
                              rc             = R_c_dust(i,j),      &
                              p_index        = pindex_dust(i,j),   &
                              q_index        = qindex_dust(i,j),   &
                              HoverR         = H_R_dust(i,j),      &
                              disc_mass      = disc_mdust(i,j),    &
                              star_mass      = star_m(i),          &
                              gamma          = gamma,              &
                              particle_mass  = massoftype(itype),  &
                              xyz_origin     = xorigini,           &
                              vxyz_origin    = vorigini,           &
                              hfact          = hfact,              &
                              xyzh           = xyzh,               &
                              vxyzu          = vxyzu,              &
                              polyk          = polyk_dust,         &
                              ismooth        = ismoothdust(i,j),   &
                              position_angle = posangl(i),         &
                              inclination    = incl(i),            &
                              rwarp          = R_warp(i),          &
                              warp_smoothl   = H_warp(i),          &
                              e0             = e0(i),              &
                              eindex         = eindex(i),          &
                              phiperi        = phiperi(i),         &
                              eccprofile     = eccprofile(i),      &
                              bh_spin        = bhspin,             &
                              enc_mass       = enc_mass(:,i),      &
                              prefix         = dustprefix(j),      &
                              lumdisc          = lumdisc_logi,         &
                              L_star           = L_star(1),            &
                              T_bg             = T_bg)

                npart = npart + npindustdisc
                npartoftype(itype) = npartoftype(itype) + npindustdisc

             enddo
          endif
       else

          !--gas disc
          npingasdisc = int(disc_m(i)/totmass_gas*np)
          call set_disc(id,master       = master,             &
                        npart           = npingasdisc,        &
                        npart_start     = npart + 1,          &
                        particle_type   = igas,               &
                        rref            = R_ref(i),           &
                        rmin            = R_in(i),            &
                        rmax            = R_out(i),           &
                        indexprofile    = iprofilegas,        &
                        rc              = R_c(i),             &
                        p_index         = pindex(i),          &
                        q_index         = qindex(i),          &
                        HoverR          = H_R(i),             &
                        disc_mass       = disc_m(i),          &
                        star_mass       = star_m(i),          &
                        gamma           = gamma,              &
                        particle_mass   = massoftype(igas),   &
                        xyz_origin      = xorigini,           &
                        vxyz_origin     = vorigini,           &
                        hfact           = hfact,              &
                        xyzh            = xyzh,               &
                        vxyzu           = vxyzu,              &
                        polyk           = polyk,              &
                        alpha           = alpha,              &
                        ismooth         = ismoothgas(i),      &
                        position_angle  = posangl(i),         &
                        inclination     = incl(i),            &
                        rwarp           = R_warp(i),          &
                        warp_smoothl    = H_warp(i),          &
                        e0              = e0(i),              &
                        eindex          = eindex(i),          &
                        phiperi         = phiperi(i),         &
                        eccprofile      = eccprofile(i),      &
                        bh_spin         = bhspin,             &
                        enc_mass        = enc_mass(:,i),      &
                        prefix          = prefix,             &
                        lumdisc         = lumdisc_logi,       &
                        L_star          = L_star(1),          &
                        T_bg            = T_bg)

          npart = npart + npingasdisc
          npartoftype(igas) = npartoftype(igas) + npingasdisc

          if (use_dust) then

             if (ndustlarge > 1) then
                dustprefix = trim(prefix)//'-'
                call make_tags_unique(ndustlarge,dustprefix)
             else
                dustprefix = trim(prefix)
             endif

             !--dust disc(s)
             do j=1,ndustlarge
                npindustdisc = int(disc_mdust(i,j)/sum(disc_mdust(:,j))*np_dust(j))
                itype = idust + j - 1

                !--taper dust disc
                iprofiledust = 0
                if (itaperdust(i,j)) iprofiledust = 1
                if (itapersetdust(i,j) == 1) iprofiledust = 2
                if (use_sigmadust_file(i,j)) iprofiledust = 3

                call set_disc(id,master      = master,             &
                              npart          = npindustdisc,       &
                              npart_start    = npart + 1,          &
                              particle_type  = itype,              &
                              rref           = R_ref(i),           &
                              rmin           = R_indust(i,j),      &
                              rmax           = R_outdust(i,j),     &
                              indexprofile   = iprofiledust,       &
                              rc             = R_c_dust(i,j),      &
                              p_index        = pindex_dust(i,j),   &
                              q_index        = qindex_dust(i,j),   &
                              HoverR         = H_R_dust(i,j),      &
                              disc_mass      = disc_mdust(i,j),    &
                              star_mass      = star_m(i),          &
                              gamma          = gamma,              &
                              particle_mass  = massoftype(itype),  &
                              xyz_origin     = xorigini,           &
                              vxyz_origin    = vorigini,           &
                              hfact          = hfact,              &
                              xyzh           = xyzh,               &
                              vxyzu          = vxyzu,              &
                              polyk          = polyk_dust,         &
                              ismooth        = ismoothdust(i,j),   &
                              position_angle = posangl(i),         &
                              inclination    = incl(i),            &
                              rwarp          = R_warp(i),          &
                              warp_smoothl   = H_warp(i),          &
                              e0             = e0(i),              &
                              eindex         = eindex(i),          &
                              phiperi        = phiperi(i),         &
                              eccprofile     = eccprofile(i),      &
                              bh_spin        = bhspin,             &
                              enc_mass       = enc_mass(:,i),      &
                              prefix         = dustprefix(j),     &
                              lumdisc        = lumdisc_logi,            &
                              L_star         = L_star(1),             &
                              T_bg           = T_bg)

                npart = npart + npindustdisc
                npartoftype(itype) = npartoftype(itype) + npindustdisc

             enddo
          endif

       endif

       !--reset alpha for each disc
       alpha_returned(i) = alpha

    endif
 enddo

 !--alpha viscosity
 if (ndiscs==1) then
    alpha = alpha_returned(onlydisc)
 else
    call warning('setup_disc', &
       'multiple discs: cannot use alpha_AV for alpha_SS, setting equal to 0.1')
    alpha = 0.1
 endif

end subroutine setup_discs

!--------------------------------------------------------------------------
!
! Set up a planetary atmosphere
!
!--------------------------------------------------------------------------
subroutine planet_atmosphere(id,npart,xyzh,vxyzu,npartoftype,gamma,hfact)
 integer, intent(in)    :: id
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 real,    intent(in)    :: gamma
 real,    intent(in)    :: hfact

 real, parameter :: a0 = 1.

 integer :: npart_disc,itype
 real    :: r_surface
 real    :: udens,rho_core

 if (surface_force) then
    npart_planet_atm = floor(Natmfrac*np)
    npart_disc = npart - npart_planet_atm

    udens = umass/udist**3
    rho_core  = rho_core_cgs/udens
    r_surface = (3./(4.*pi)*m1/rho_core)**(1./3.)
    !--Note the surface of the planet is located at the Plumber softening length
    eps_soft1 = r_surface
    if (eps_soft1 <= 0.) then
       print*,'Something wrong in the surface radius: eps_soft1 =',eps_soft1
    endif
 else
    npart_disc = npart
 endif

 itype = igas

 !--set up an atmosphere around one of the binary masses (i.e. planet)
 if (surface_force .and. npart_planet_atm > 0) then
    call set_planet_atm(id,xyzh,vxyzu,npartoftype,maxvxyzu,itype,a0,R_in(1), &
                        H_R(1),m2,qindex(1),gamma,Ratm_in,Ratm_out,r_surface, &
                        npart,npart_planet_atm,npart_disc,hfact)
 endif

 !--move into the corotating frame with the planet
 if (surface_force .or. iexternalforce == iext_corotate) then
    call make_corotate(xyzh,vxyzu,a0,m2,npart,npart_disc)
 endif

end subroutine planet_atmosphere

!--------------------------------------------------------------------------
!
! Set sphere of particles around one of the binary masses (i.e. planet)
!
!--------------------------------------------------------------------------
subroutine set_planet_atm(id,xyzh,vxyzu,npartoftype,maxvxyzu,itype,a0,R_in, &
                          HoverR,Mstar,q_index,gamma,Ratm_in,Ratm_out,r_surface, &
                          npart,npart_planet_atm,npart_disc,hfact)
 use part,          only:set_particle_type,igas
 use spherical,     only:set_sphere,rho_func
 integer, intent(in)    :: id
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 integer, intent(in)    :: maxvxyzu
 integer, intent(in)    :: itype
 real,    intent(in)    :: a0
 real,    intent(in)    :: R_in
 real,    intent(in)    :: HoverR
 real,    intent(in)    :: Mstar
 real,    intent(in)    :: q_index
 real,    intent(in)    :: gamma
 real,    intent(inout) :: Ratm_in
 real,    intent(inout) :: Ratm_out
 real,    intent(in)    :: r_surface
 integer, intent(inout) :: npart
 integer, intent(inout) :: npart_planet_atm
 integer, intent(inout) :: npart_disc
 real,    intent(in)    :: hfact

 integer(kind=8)    :: nptot
 integer            :: i,nx
 real               :: xyz_orig(3)
 real               :: a_orbit
 real               :: psep,vol_sphere
 real               :: cs0,cs
 procedure(rho_func), pointer :: density_func
 !
 ! place particles in sphere
 !
 Ratm_in   = Ratm_in*r_surface
 Ratm_out  = Ratm_out*r_surface

 if (ramp) then
    xyz_orig(:) = (/a0,0.,0./)
 else
    a_orbit = a0 - mass2/Mstar
    xyz_orig(:) = (/a_orbit,0.,0./)
 endif

 vol_sphere  = 4./3.*pi*Ratm_out**3
 nx          = int(npart_planet_atm**(1./3.))
 psep        = vol_sphere**(1./3.)/real(nx)
 nptot       = npart
 density_func => atm_dens

 call set_sphere('closepacked',id,master,Ratm_in,Ratm_out,psep,hfact,npart,xyzh, &
                 rhofunc=density_func,nptot=nptot, &
                 np_requested=npart_planet_atm,xyz_origin=xyz_orig)

 npart_planet_atm = npart-npart_disc
 npartoftype(igas) = npart
 do i=npart_disc+1,npart
    !--set the particle type for the atmosphere particles
    call set_particle_type(i,itype)
    !--------------------------------------------------------------------------
    !  Set thermal energy
    !  utherm generally should not be stored
    !  for an isothermal equation of state
    if (maxvxyzu >= 4) then
       if (itype==igas) then
          cs0 = sqrt(1.0/R_in)*(HoverR)*sqrt(Mstar)*(1.0/R_in)**(-q_index)
          cs = cs_func(cs0,a0,q_index)
          if (gamma > 1.) then
             vxyzu(4,i) = cs**2/(gamma - 1.)/gamma
          else
             vxyzu(4,i) = 1.5*cs**2
          endif
       else
          vxyzu(4,i) = 0.
       endif
    endif
 enddo

end subroutine set_planet_atm

subroutine set_sphere_around_disc(id,npart,xyzh,vxyzu,npartoftype,massoftype,hfact)
 use partinject,     only:add_or_update_particle
 use velfield,       only:set_velfield_from_cubes
 use datafiles,      only:find_phantom_datafile
 use part,           only:set_particle_type,igas,gravity,dustfrac,ndustsmall
 use set_dust,       only:set_dustfrac,set_dustbinfrac
 use options,        only:use_dustfrac
 use spherical,      only:set_sphere,rho_func
 use dim,            only:maxp
 use physcon,        only:pi
 use eos,            only:get_spsound,gmw,cs_min
 use units,          only:get_kbmh_code,get_G_code,utime
 integer, intent(in)    :: id
 integer, intent(inout) :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 real, intent(inout) :: massoftype(:)
 real,    intent(in)    :: hfact
 integer :: i, ipart
 integer :: itype

 integer :: n_add, np
 integer(kind=8) :: nptot
 real :: delta, pmass, mtot, mdisc, omega, Routmax, Poutmax, ff_in, ff_out
 real :: v_ff_mag, vxi, vyi, vzi, my_vrms, factor, x_pos, y_pos, z_pos
 real :: rhoi, spsound, rms_in, temp, dustfrac_tmp, vol_obj, rpart, rc_in, rc_out, G_code
 integer :: ierr
 real, dimension(:,:), allocatable :: xyzh_add,vxyzu_add
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=120)           :: filex,filey,filez

 itype = igas
 pmass = massoftype(igas)
 mtot = sum(xyzmh_ptmass(4,:))
 mdisc = pmass*npart
 dustfrac_tmp = 0.

 n_add = nint(mass_sphere / pmass)

 if (use_dust) then
    if (use_dustfrac) then
       write(*,*) "Detected one-fluid dust in the simulation, adding smallest dust to cloud."
       if (dustfrac_method == -1) then
          dustfrac_tmp = 0.
       elseif (dustfrac_method == 0) then
          dustfrac_tmp = dust_to_gas
       elseif (dustfrac_method == 1) then
          dustfrac_tmp = sum(dustfrac(1,:npartoftype(igas)))/real(npartoftype(igas))
       endif
       n_add = nint( (mass_sphere * (1. + dustfrac_tmp)) / pmass )
       write(*,*) 'Setting dust-to-gas ratio in the cloud to ',dustfrac_tmp
    endif
 endif

 nptot =  n_add + npartoftype(igas)
 write(*,*) 'Adding ',n_add,' particles to cloud.'

 G_code = get_G_code()
 if (gravity) then
    mtot = mtot + mdisc
 endif

 if (add_rotation == 1) then
    write(*,*) 'Adding Keplerian rotation in the cloud.'
    omega = Kep_factor * sqrt((G_code*mtot)/R_rot**3)
 elseif (add_rotation == 2) then
    write(*,*) 'Adding constant angular velocity in the cloud.'
    omega = omega_cloud*utime
 else
    omega = 0.
 endif

 np = 0

 allocate(xyzh_add(4,n_add),vxyzu_add(4,n_add))

 delta = 1.0
 call set_sphere('random',id,master,Rin_sphere,Rout_sphere,delta,hfact,np,xyzh_add,xyz_origin=(/0., 0., 0./),&
                  np_requested=n_add, nptot=nptot)

 vxyzu_add(1,:) = 0.
 vxyzu_add(2,:) = 0.
 vxyzu_add(3,:) = 0.
 vxyzu_add(4,:) = 5.868e-05 ! T=10K, doesn't seem to be used

 if (add_turbulence==1) then

    filex = find_phantom_datafile(filevx,'velfield')
    filey = find_phantom_datafile(filevy,'velfield')
    filez = find_phantom_datafile(filevz,'velfield')

    call set_velfield_from_cubes(xyzh_add,vxyzu_add,n_add,filex,filey,filez,1.,tfact*Rout_sphere,.false.,ierr)

    if (ierr /= 0) call fatal('setup','error setting up velocity field')

    vol_obj = 4.0/3.0*pi*(Rout_sphere**3 - Rin_sphere**3)

    rhoi = mass_sphere/vol_obj

    if (cs_min > 0.) then
       spsound = cs_min
    else
       write(*,*) 'Warning: Floor temperature not set, assuming T_floor = 10 K'
       temp = 10.
       spsound = sqrt(temp*get_kbmh_code()/gmw)
    endif

    rms_in = spsound*rms_mach

    !--Normalise the energy
    ! rms_curr = sqrt( 1/real(n_add)*sum( (vxyzu_add(1,:)**2 + vxyzu_add(2,:)**2 + vxyzu_add(3,:)**2) ) )

    my_vrms = 0.
    do i=1,n_add
       vxi  = vxyzu_add(1,i)
       vyi  = vxyzu_add(2,i)
       vzi  = vxyzu_add(3,i)
       my_vrms = my_vrms + vxi*vxi + vyi*vyi + vzi*vzi
    enddo

    ! Normalise velocity field
    my_vrms = sqrt(1/real(n_add) * my_vrms)
    factor = rms_in/my_vrms
    do i=1,n_add
       vxyzu_add(1:3,i) = vxyzu_add(1:3,i)*factor
    enddo
 endif

 if (set_freefall == 1) then
    do i=1,n_add
       x_pos = xyzh_add(1,i)
       y_pos = xyzh_add(2,i)
       z_pos = xyzh_add(3,i)

       rpart = sqrt(x_pos*x_pos + y_pos*y_pos + z_pos*z_pos)

       if (rpart > 1.0e-12_8 .and. mtot > 0.0) then
          v_ff_mag = sqrt(2.0 * mtot / rpart)
          vxyzu_add(1,i) = vxyzu_add(1,i) - v_ff_mag * x_pos / rpart
          vxyzu_add(2,i) = vxyzu_add(2,i) - v_ff_mag * y_pos / rpart
          vxyzu_add(3,i) = vxyzu_add(3,i) - v_ff_mag * z_pos / rpart
       endif
    enddo
 endif

 ipart = npart
 do i = 1,n_add
    ipart = ipart + 1
    if (add_rotation == 1) then
       vxyzu_add(1,i) = vxyzu_add(1,i) - omega*xyzh_add(2,i)
       vxyzu_add(2,i) = vxyzu_add(2,i) + omega*xyzh_add(1,i)
    endif
    call add_or_update_particle(igas, xyzh_add(1:3,i), vxyzu_add(1:3,i), xyzh_add(4,i), &
                                 vxyzu_add(4,i), ipart, npart, npartoftype, xyzh, vxyzu)
    if (use_dust) then
       if (use_dustfrac) then
          dustfrac(1, ipart) = dustfrac_tmp
          dustfrac(2:ndustsmall, ipart) = 0.
       endif
    endif
 enddo

 npartoftype(igas) = ipart

 if (npartoftype(igas) > maxp) call fatal('set_sphere_around_disc', &
      'maxp too small, rerun with --maxp=N where N is desired number of particles')

 Routmax = maxval(R_out)
 Poutmax = 2./pi * sqrt(Routmax**3/G_code*mtot)
 write(*,*) 'Maximum disc radius is ',Routmax
 write(*,*) 'Period of outer disc is ',Poutmax*utime/3.15576e7, ' years'

 rc_in   = Rin_sphere**4 * omega**2 / (G_code*mtot)
 rc_out  = Rout_sphere**4 * omega**2 / (G_code*mtot)
 ff_in = sqrt(Rin_sphere**3/(2.*G_code*mtot))
 ff_out = sqrt(Rout_sphere**3/(2.*G_code*mtot))
 write(*,*) 'Centrifugal radius at minimum cloud radius is ', rc_in
 write(*,*) 'Centrifugal radius at maximum cloud radius is ', rc_out
 write(*,*) 'Free-fall time at minimum cloud radius is ', ff_in*utime/3.15576e7, ' years'
 write(*,*) 'which is ', ff_in/Poutmax, ' times the period of the outer disc.'
 write(*,*) 'Free-fall time at maximum cloud radius is ', ff_out*utime/3.15576e7, ' years'
 write(*,*) 'which is ', ff_out/Poutmax, ' times the period of the outer disc.'

end subroutine set_sphere_around_disc

!--------------------------------------------------------------------------
!
! Initialise the dustprop array
!
!--------------------------------------------------------------------------
subroutine initialise_dustprop(npart)
 use physcon,     only:fourpi
 integer, intent(in) :: npart

 integer :: i,iam

 if (use_dustgrowth) then
    do i=1,npart
       iam = iamtype(iphase(i))
       if (iam==idust .or. (use_dustfrac .and. iam==igas)) then
          dustprop(1,i) = fourpi/3.*graindens(1)*grainsize(1)**3
          dustprop(2,i) = graindens(1)
       else
          dustprop(:,i) = 0.
       endif
       filfac(i) = 0.
       probastick(i) = 1.
       dustgasprop(:,i) = 0.
       VrelVf(i)        = 0.
    enddo
 endif

end subroutine initialise_dustprop

!--------------------------------------------------------------------------
!
! Print angular momentum information
!
!--------------------------------------------------------------------------
subroutine print_angular_momentum(npart,xyzh,vxyzu)
 use orbits, only:get_mean_angmom_vector
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 real,    intent(in) :: vxyzu(:,:)
 real :: ldisc(3),lcentral(3)

 if (npart <= 0) return

 ldisc = get_mean_angmom_vector(npart,xyzh,vxyzu)
 print "(a,'(',3(es10.2,1x),')')",' Disc specific angular momentum = ',ldisc
 if (nptmass > 1) then
    lcentral = get_mean_angmom_vector(nptmass,xyzmh_ptmass,vxyz_ptmass)
    print "(a,'(',3(es10.2,1x),')')",' Binary specific angular momentum = ',lcentral
    ! make unit vectors
    lcentral = lcentral/sqrt(dot_product(lcentral,lcentral))
    ldisc    = ldisc/sqrt(dot_product(ldisc,ldisc))
    print "(a,f6.1,a)",' Angle between disc and binary = ',acos(dot_product(lcentral,ldisc))*180./pi,' deg'
 endif

end subroutine print_angular_momentum

!--------------------------------------------------------------------------
!
! Print dust information
!
!--------------------------------------------------------------------------
subroutine print_dust()

 use grids_for_setup, only:init_grid_sigma,deallocate_sigma
 character(len=20) :: duststring(maxdusttypes)
 integer           :: i,j
 real              :: Sigma
 real              :: Sigmadust
 real              :: Stokes(maxdusttypes)
 real              :: R_midpoint

 if (use_dust) then

    print "(/,a)",' --------------------- added dust ---------------------'

    if (use_dustgrowth) then
       print "(a,g10.3,a)", ' initial grain size: ',grainsize(1)*udist,' cm'
    else
       duststring = 'grain size'
       call make_tags_unique(ndusttypes,duststring)
       do i=1,ndusttypes
          print*,adjustr(duststring(i))//' : ',grainsize(i)*udist,' cm'
       enddo
       duststring = 'grain density'
       call make_tags_unique(ndusttypes,duststring)
       do i=1,ndusttypes
          print*,adjustr(duststring(i))//' : ',graindens(i)*umass/udist**3,' g/cm^3'
       enddo
    endif

    do i=1,maxdiscs
       if (iuse_disc(i)) then
          if (sigmaprofilegas(i)==6) call init_grid_sigma(R_in(i),R_out(i))
          R_midpoint = (R_in(i) + R_out(i))/2
          Sigma = sig_norm(i) * &
                  scaled_sigma(R_midpoint,sigmaprofilegas(i),pindex(i),R_ref(i),R_in(i),R_out(i),R_c(i))
          Sigmadust = 0.
          do j=1,ndusttypes
             Sigmadust = Sigmadust + sig_normdust(i,j) * &
                         scaled_sigma(R_midpoint,sigmaprofiledust(i,j),pindex_dust(i,j),&
                                      R_ref(i),R_indust(i,j),R_outdust(i,j),R_c_dust(i,j))
          enddo
          Stokes = 0.5*pi*graindens*grainsize/(Sigma+Sigmadust)
          duststring = 'approx. Stokes'
          call make_tags_unique(ndusttypes,duststring)
          if (ndiscs > 1) then
             print "(/,a,i2,a)",' ---------------------   disc',i,'   ---------------------'
          endif
          do j=1,ndusttypes
             print*,'',adjustr(duststring(j))//' : ',Stokes(j)
          enddo
          if (sigmaprofilegas(i)==6) call deallocate_sigma()
       endif
    enddo
    print "(1x,54('-'),/)"

 else
    print "(/,a,/)",' There is no dust here!'
 endif

end subroutine print_dust

!--------------------------------------------------------------------------
!
! Set up planets
!
!--------------------------------------------------------------------------
subroutine set_planets(npart,massoftype,xyzh)
 use vectorutils, only:rotatevec
 use setbinary, only:set_binary
 use part, only:nsinkproperties
 integer, intent(in) :: npart
 real,    intent(in) :: massoftype(:)
 real,    intent(in) :: xyzh(:,:)

 integer :: i,j,itype,ntmp,ierr
 real    :: dist_bt_sinks,dx(3)
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r
 real    :: Hill(maxplanets)
 real    :: xyz_tmp(nsinkproperties,2),vxyz_tmp(3,2)
 real    :: mtot

 period_planet_longest = 0.
 if (nplanets > 0) then
    print "(a,i2,a)",' --------- added ',nplanets,' planets ------------'
    do i=1,nplanets
       nptmass = nptmass + 1
       phi = 0.
       cosphi = cos(phi*deg_to_rad)
       sinphi = sin(phi*deg_to_rad)
       disc_m_within_r = 0.

       !--disc mass correction
       do j=1,npart
          r2 = xyzh(1,j)**2 + xyzh(2,j)**2 + xyzh(3,j)**2
          if (r2 < rplanet(i)**2) then
             itype = iamtype(iphase(j))
             disc_m_within_r = disc_m_within_r + massoftype(itype)
          endif
       enddo

       !--inner planet mass correction
       if (nplanets>1) then
          do j=1,nplanets
             if (rplanet(j)<rplanet(i)) disc_m_within_r = disc_m_within_r + mplanet(j)*jupiterm/umass
          enddo
       endif

       !--set sink particles
       Hill(i) = (mplanet(i)*jupiterm/umass/(3.*mcentral))**(1./3.) * rplanet(i)
       if (nsinks == 2) then
          dx = xyzmh_ptmass(1:3,1)-xyzmh_ptmass(1:3,2)
          dist_bt_sinks = sqrt(dot_product(dx,dx))
          if (rplanet(i) > dist_bt_sinks) Hill(i) = (mplanet(i)*jupiterm/umass/(3.*m1))**(1./3.) * rplanet(i)
       else
          dist_bt_sinks = 0.
       endif
       mtot = mcentral+disc_m_within_r
       if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) mtot = m1 + disc_m_within_r
       ntmp = 0
       call set_binary(mtot,0.,semimajoraxis=rplanet(i),eccentricity=eccplanet(i),&
                       accretion_radius1=0.0,accretion_radius2=accrplanet(i)*Hill(i),&
                       xyzmh_ptmass=xyz_tmp,vxyz_ptmass=vxyz_tmp,nptmass=ntmp,ierr=ierr,incl=inclplan(i),&
                       arg_peri=wplanet(i),posang_ascnode=Oplanet(i),f=fplanet(i),verbose=.false.)
       xyzmh_ptmass(1:3,nptmass) = xyz_tmp(1:3,2)
       vxyz_ptmass(1:3,nptmass) = vxyz_tmp(1:3,2)

       !xyzmh_ptmass(1:3,nptmass)    = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
       xyzmh_ptmass(4,nptmass)      = mplanet(i)*jupiterm/umass
       xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i)*Hill(i)
       xyzmh_ptmass(ihsoft,nptmass) = 0.
       vphi                         = sqrt((mcentral + disc_m_within_r)/rplanet(i))
       !if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) vphi = sqrt((m1 + disc_m_within_r)/rplanet(i))
       !vxyz_ptmass(1:3,nptmass)     = (/-vphi*sinphi,vphi*cosphi,0./)
       if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) then
          vxyz_ptmass(1:3,nptmass) = vxyz_ptmass(1:3,nptmass) + vxyz_ptmass(1:3,1)
          xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass) + xyzmh_ptmass(1:3,1)
       endif

       !--incline positions and velocities
       !inclplan(i) = inclplan(i)*deg_to_rad
       !u = (/-sin(phi),cos(phi),0./)
       !call rotatevec(xyzmh_ptmass(1:3,nptmass),u,-inclplan(i))
       !call rotatevec(vxyz_ptmass(1:3,nptmass), u,-inclplan(i))

       !--compute obliquity and spin angular momentum
       if (abs(J2planet(i)) > 0.) then
          call set_sink_oblateness(nptmass,J2planet(i),planet_size(i),spin_period(i),kfac(i),obliquity(i))
       endif

       !--print planet information
       omega = vphi/rplanet(i)
       print "(a,i2,a)",             ' >>> planet ',i,' <<<'
       print "(a,g10.3,a)",          ' orbital radius: ',rplanet(i)*udist/au,' au'
       print "(a,g10.3,a,2pf7.3,a)", '          M(<R): ',(disc_m_within_r + mcentral)*umass/solarm, &
                                     ' MSun, disc mass correction is ',disc_m_within_r/mcentral,'%'
       print "(a,g10.3,a)",          '    planet mass: ',mplanet(i),' MJup'
       print "(a,g10.3,a)",          '    planet mass: ',mplanet(i)*jupiterm/earthm,' MEarth'
       print "(a,g10.3,a)",          '    planet mass: ',mplanet(i)*jupiterm/solarm,' MSun'
       print "(a,2(g10.3,a))",       ' orbital period: ',2.*pi*rplanet(i)/vphi*utime/years,' years or ', &
                                                 2*pi*rplanet(i)/vphi,' in code units'
       print "(a,g10.3,a)",          '    Hill radius: ',Hill(i),' au'
       print "(a,g10.3,a,i3,a)",     '   accr. radius: ',xyzmh_ptmass(ihacc,nptmass),' au or ', &
                                                 int(100*xyzmh_ptmass(ihacc,nptmass)/Hill(i)), ' % of Hill radius'
       print "(a,g10.3,a)",          '     resonances:'
       print "(a,g10.3,a)",   '    3:1 : ',(sqrt(mcentral)/(3.*omega))**(2./3.)*udist/au,' au'
       print "(a,g10.3,a)",   '    4:1 : ',(sqrt(mcentral)/(4.*omega))**(2./3.)*udist/au,' au'
       print "(a,g10.3,a)",   '    5:1 : ',(sqrt(mcentral)/(5.*omega))**(2./3.)*udist/au,' au'
       print "(a,g10.3,a)",   '    9:1 : ',(sqrt(mcentral)/(9.*omega))**(2./3.)*udist/au,' au'
       call print_oblateness_info(nptmass,spin_period(i))

       !--check planet accretion radii
       if (accrplanet(i) < 0.05) then
          call warning('setup_disc','accretion radius of planet < 1/20 Hill radius: unnecessarily small')
       elseif (accrplanet(i) > 0.5) then
          call warning('setup_disc','accretion radius of planet > Hill radius: too large')
       elseif (accrplanet(i)*Hill(i) > accr1) then
          call warning('setup_disc','accretion radius of planet > accretion radius of primary star: this is unphysical')
       endif
       if (xyzmh_ptmass(iReff,nptmass) > 0.25*Hill(i)) then
          call warning('setup_disc','planet size exceeds 1/4 of Hill radius: too large')
       endif
       if (xyzmh_ptmass(iReff,nptmass) > max(xyzmh_ptmass(ihacc,nptmass),xyzmh_ptmass(ihsoft,nptmass))) then
          call warning('setup_disc','planet size exceeds accretion radius: too large')
       endif
       print *, ''

       !--determine longest period
       period_planet_longest = max(period_planet_longest, 2.*pi/omega)

    enddo
    print "(1x,45('-'))"
    print *, ''
 endif

end subroutine set_planets

!--------------------------------------------------------------------------
!
!  Set properties needed for geopotential forces from sink particles
!
!--------------------------------------------------------------------------
subroutine set_sink_oblateness(isink,J2,planet_size,spin_period_hrs,kfac,obliquity)
 use physcon, only:jupiterr
 integer, intent(in) :: isink
 real, intent(in) :: J2,planet_size,spin_period_hrs,kfac,obliquity
 real :: spin_am,planet_radius,planet_spin_period

 xyzmh_ptmass(iJ2,isink) = J2
 ! compute spin angular momentum of the body
 planet_radius = planet_size*jupiterr/udist
 planet_spin_period = spin_period_hrs*hours/utime
 spin_am = twopi*kfac*(xyzmh_ptmass(4,isink)*planet_radius**2)/planet_spin_period
 xyzmh_ptmass(ispinx,isink) = spin_am*sin(obliquity*deg_to_rad)
 xyzmh_ptmass(ispinz,isink) = spin_am*cos(obliquity*deg_to_rad)
 xyzmh_ptmass(iReff,isink) = planet_radius

end subroutine set_sink_oblateness

!--------------------------------------------------------------------------
!
!  Reset centre of mass to origin
!
!--------------------------------------------------------------------------
subroutine set_centreofmass(npart,xyzh,vxyzu)
 use centreofmass, only:reset_centreofmass
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)

 integer :: npart_recentre

 if (iexternalforce == iext_corot_binary) then
    npart_recentre = npart - npart_planet_atm
 else
    npart_recentre = npart
 endif
 call reset_centreofmass(npart_recentre,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

end subroutine set_centreofmass

!--------------------------------------------------------------------------
!
!  Set tmax and dtmax for infile
!
!--------------------------------------------------------------------------
subroutine set_tmax_dtmax(fileprefix)
 use orbits,   only:get_T_flyby_hyp,get_T_flyby_par,get_orbital_period,get_time_between_true_anomalies
 use setorbit, only:write_trajectory_to_file
 use timestep, only:tmax,dtmax
 use units,    only:in_code_units
 character(len=*), intent(in) :: fileprefix
 real :: period,period1,period2,mu
 real :: flyby_d
 integer :: ierr

 period2 = 0.
 if (icentral==1 .and. nsinks==2 .and. binary%input_type==0) then
    !--binary orbital period
    period = get_orbital_period(mcentral,binary%a)
 elseif (icentral==1 .and. nsinks==3 .and. binary%input_type==0) then
    !--wide binary orbital period
    period = get_orbital_period(mcentral,binary%a)
    !--tight binary orbital period
    period2 = get_orbital_period(m2,binary2%a)
 elseif (icentral==1 .and. nsinks==4 .and. binary%input_type==0) then
    !--wide binary orbital period
    period = get_orbital_period(mcentral,binary%a)
    !--tight binary 1 orbital period
    period1 = get_orbital_period(m1,binary1%a)
    !--tight binary 2 orbital period
    period2 = get_orbital_period(m2,binary2%a)
 elseif (icentral==1 .and. nsinks==2 .and. binary%input_type>=1) then
    !--time of flyby
    mu = m1+m2
    if (binary%input_type==2) then
       ! for Orbit Reconstructor^TM input, compute time to reach observed separation
       period = get_time_between_true_anomalies(mu,binary%a,binary%e,binary%f,binary%obs%f)
       call write_trajectory_to_file(binary,m1,m2,fileprefix)
    else
       if (binary%e > 1.0) then
          period = get_T_flyby_hyp(mu,binary%e,binary%f,binary%a)
       else
          flyby_d = in_code_units(binary%flyby%d,ierr,unit_type='length')
          period = get_T_flyby_par(mu,binary%a,flyby_d/binary%a)
       endif
    endif
 elseif (nplanets > 0) then
    !--outer planet orbital period
    period = period_planet_longest
 elseif (iwarp(onlydisc)) then
    !--warp period
    period = get_orbital_period(mcentral,R_warp(onlydisc))
 else
    !--outer disc orbital period
    period = get_orbital_period(mcentral,R_out(onlydisc))
 endif

 if (period > 0. .and. nsinks<3) then
    if (deltat > 0.) dtmax = deltat*period
    if (norbits >= 0) tmax = norbits*period
 elseif (period > 0. .and. nsinks==3) then
    if (deltat < 0. .and. period2 > 0.) then
       dtmax = -deltat*period2
    elseif (deltat > 0.) then
       dtmax = deltat*period
    endif
 elseif (nsinks==4) then
    dtmax = deltat*period
    if (norbits < 0 .and. period2 > 0.) then
       tmax = -norbits*period2
    elseif (norbits >= 0) then
       tmax = norbits*period
    endif
 endif

end subroutine set_tmax_dtmax

!--------------------------------------------------------------------------
!
!  Prompt user for desired setup options
!
!--------------------------------------------------------------------------
subroutine setup_interactive(id)
 use prompting,        only:prompt
 use set_dust_options, only:set_dust_interactive
 use sethierarchical,  only:set_hierarchical_default_options,get_hier_level_mass
 use sethierarchical,  only:hs,hierarchy,print_chess_logo,generate_hierarchy_string
 use sethier_utils,    only:findloc_local

 integer, intent(in) :: id
 integer :: i
 real    :: disc_mfac(maxdiscs)

 !--central object(s)
 print "(a)",'==========================='
 print "(a)",'+++  CENTRAL OBJECT(S)  +++'
 print "(a)",'==========================='
 call prompt('Do you want to use sink particles or an external potential?'// &
             new_line('A')//' 0=potential'//new_line('A')//' 1=sinks'// &
             new_line('A'),icentral,0,1)
 select case (icentral)
 case (0)
    !--external potential
    call prompt('Which potential?'//new_line('A')// &
               ' 1=central point mass'//new_line('A')// &
               ' 2=binary potential'//new_line('A')// &
               ' 3=spinning black hole'//new_line('A'),ipotential,1,3)
    select case (ipotential)
    case (1)
       !--point mass
       iexternalforce = iext_star
       m1       = 1.
       accr1    = 1.
    case (2)
       !--fixed binary
       call prompt('Do you want model a surface on one mass (i.e. planet)?',surface_force)
       if (surface_force) then
          iexternalforce = iext_corot_binary
          !--binary
          m1       = 0.001 ! Planet
          m2       = 1.    ! Central star
          accr1    = 0.    ! Surface force, so we don't need an accretion radius
          accr2    = 0.1

          print "(/,a,/)",' NOTE: Using fixed potential, so planet is at R=1 in code units; change dist_unit to "move" the planet'
          call prompt('Enter average core density for the planet (g/cm^3)',rho_core_cgs,0.)
          call prompt('Enter inner atmosphere radius in planet radii',Ratm_in,1.,10.)
          call prompt('Enter outer atmosphere radius in planet radii',Ratm_out,Ratm_in,10.)
          call prompt('Enter atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))',atm_type,1,2)
          call prompt('Enter fraction of particles to be used in planet atmosphere',Natmfrac,0.,1.)
          call prompt('Do you want to ramp up the planet mass slowly?',ramp)
       else
          iexternalforce = iext_binary
          m1       = 1.
          m2       = 1.
          accr1    = 1.
          accr2    = 1.
       endif
    case (3)
       !--spinning black hole (Lense-Thirring)
       iexternalforce = iext_lensethirring
       call prompt('Include Einstein precession?',einst_prec)
       if (einst_prec) iexternalforce = iext_einsteinprec
       m1          = 1.
       accr1       = 30.
       bhspin      = 1.
       bhspinangle = 0.
    end select
 case (1)
    !--sink particle(s)
    call prompt('How many sinks?',nsinks,1)
    select case (nsinks)
    case (1)
       !--single star
       m1       = 1.
       accr1    = 1.
    case (2)
       !--binary
       call prompt('Do you want to specify the binary orbit as bound (elliptic) or'// &
                  ' unbound (parabolic/hyperbolic) [flyby] or as observed dx,dv?'//new_line('A')// &
                  ' 0=bound'//new_line('A')//' 1=unbound'//new_line('A')//' 2=observed dx,dv'//new_line('A'),binary%input_type,0,3)
       select case (binary%input_type)
       case (0)
          !--bound
          m1       = 1.
          m2       = 0.2
          binary%elems%a = '10.'
          accr1    = 1.
          accr2    = 0.5
       case default
          !--unbound (flyby)
          m1       = 1.
          m2       = 1.
          accr1    = 1.
          accr2    = 1.
          !
          ! the following is only if we want to override defaults in set_defaults_orbit
          ! so for input_type=2 or 3 we will just get those defaults
          !
          if (binary%input_type >= 1) then
             binary%flyby%rp = '200.'
             binary%flyby%d = '2000.'
          endif
          binary%e = 2.0
       end select

    case (5:)

       call print_chess_logo()!id)

       binary%input_type = 0

       call generate_hierarchy_string(nsinks)

       call prompt('What is the hierarchy?',hierarchy)
       !call set_hierarchical_interactively()
       call set_hierarchical_default_options()

    case (3)
       !-- hierarchical triple --!
       print "(/,a)",'================================'
       print "(a)",  '+++   HIERARCHICAL TRIPLE    +++'
       print "(a)",  '================================'
       binary%input_type = 0

       !-- Wide binary
       m1       = 1.
       m2       = 0.2
       binary%elems%a = '10.'
       accr1    = 1.

       !-- Tight binary
       subst    = 12
       q2       = 1
       m2a      = m2/(q2+1)
       m2b      = m2*q2/(q2+1)
       binary2%elems%a = '1.'
       accr2a    = 0.1
       accr2b    = 0.1
    case(4)
       !-- hierarchical quadruple --!
       print "(/,a)",'================================'
       print "(a)",  '+++   HIERARCHICAL QUADRUPLE    +++'
       print "(a)",  '================================'
       binary%input_type = 0

       !-- Wide binary
       m1       = 1.
       m2       = 0.2
       binary%elems%a = '10.'
       accr1    = 1.

       !-- Tight binary 1
       subst1    = 11
       q1        = 1
       m1a       = m1/(q1+1)
       m1b       = m1*q1/(q1+1)
       binary1%elems%a = '1.'
       accr1a    = 0.1
       accr1b    = 0.1

       !-- Tight binary 2
       subst2    = 12
       q2       = 1
       m2a      = m2/(q2+1)
       m2b      = m2*q2/(q2+1)
       binary2%elems%a = '1.'
       accr2a    = 0.1
       accr2b    = 0.1

    end select
 end select

 !--multiple disc options
 print "(/,a)",'================='
 print "(a)",  '+++  DISC(S)  +++'
 print "(a)",  '================='
 if ((icentral==1) .and. (nsinks>=2)) then
    !--multiple discs possible
    if (binary%input_type==0 .and. nsinks==2) then
       !--bound binary: circum-binary, -primary, -secondary
       iuse_disc(1) = .true.
       iuse_disc(2) = .false.
       iuse_disc(3) = .false.
       call prompt('Do you want a circumbinary disc?',iuse_disc(1))
       call prompt('Do you want a circumprimary disc?',iuse_disc(2))
       call prompt('Do you want a circumsecondary disc?',iuse_disc(3))
    elseif (binary%input_type>=1) then
       !--unbound binary (flyby): circum-primary, -secondary
       iuse_disc(1) = .false.
       iuse_disc(2) = .true.
       iuse_disc(3) = .false.
       call prompt('Do you want a circumprimary disc?',iuse_disc(2))
       call prompt('Do you want a circumsecondary disc?',iuse_disc(3))
    elseif (nsinks==3) then
       !--bound binary: circum-triple
       iuse_disc(1) = .false.
       iuse_disc(2) = .false.
       iuse_disc(3) = .false.
       iuse_disc(4) = .true.
       call prompt('Do you want a circum-triple disc?',iuse_disc(4))
       if (.not.iuse_disc(4)) then
          call prompt('Do you want a circumbinary disc around the first hierarchical level secondary?',iuse_disc(1))
       else
          print "(/,a)",'Setting circum-triple disc.'
       endif
    elseif (nsinks==4) then
       !--2 bound binaries: circumbinary
       iuse_disc(1) = .true.
       iuse_disc(2) = .false.
       iuse_disc(3) = .false.
       iuse_disc(4) = .false.
       print "(/,a)",'Setting circumbinary disc around the first hierarchical level secondary.'
    elseif (nsinks>=5) then
       !--2 bound binaries: circumbinary
       iuse_disc(:) = .false.

       do i=1,hs%labels%sink_num
          call prompt('Do you want a disc orbiting '//trim(hs%labels%sink(i))//' star?',iuse_disc(i))
       enddo

       do i=1,hs%labels%hl_num
          call prompt('Do you want a disc orbiting '//trim(hs%labels%hl(i))//' hierarchical level?',iuse_disc(i+hs%labels%sink_num))
       enddo

    endif
    if (.not.any(iuse_disc)) iuse_disc(1) = .true.
    !--number of discs
    ndiscs = count(iuse_disc)
    if (ndiscs > 1) then
       use_global_iso = .false.
    endif
 endif

 !--gas disc
 R_in  = accr1
 R_ref = min(10.*R_in,R_out)
 R_c   = R_out
 disc_mfac = 1.
 if (ndiscs > 1) qindex = 0.
 alphaSS = 0.005
 if (surface_force) then
    R_in  = 0.1
    R_out = 3.
    R_ref = 1.
 endif
 if ((icentral==1 .and. nsinks>=2) .and. (binary%input_type==0)) then
    !--don't smooth circumbinary, by default
    ismoothgas(1) = .false.
    !--set appropriate disc radii for bound binary
    R_in(1:4)      = (/2.5*binary%a, accr1, accr2, 2.5*binary%a /)
    R_out(1:4)     = (/5.*R_in(1), 5.*accr1, 5.*accr2, 5.*R_in(1) /)
    R_ref     = R_in
    R_c       = R_out
    disc_mfac(1:4) = (/1., 0.1, 0.01, 1./)
    if (ndiscs > 1) then
       !--set H/R so temperature is globally constant
       call prompt('Do you want a globally isothermal disc (if not Farris et al. 2014)?',use_global_iso)
       !--------------------------------------------------------------------------
       ! N.B. The initializations of multiple discs is not done using the
       ! implementation of the eos a radial profile centred on CM, primary and
       ! secondary is used. The value of H_R used in setpart to set cs0 is the
       ! one of the circumbinary if cb disc is present, otherwise it uses the
       ! circumprimary. The values of H_R used for the other discs are set using
       ! the equations below, however changing them here is not enough. They need
       ! to be changed also in the the equation_of_state function in this module.
       !--------------------------------------------------------------------------
       if (.not. use_global_iso) then
          if (maxvxyzu > 3) then
             call prompt("Do you want to set the disc temperatures from the stellar"// &
             "luminosity? (0=no 1=yes",lumdisc)
          endif
          if (lumdisc > 0) then
             !get luminosity ...
             call prompt("Enter the luminosity of star",L_star(1))
             call prompt("Enter the background temperature e.g. 10 (K)", T_bg)
             qindex(1) = 0.25
             qindex = 0.25
          else
             call prompt('Enter q_index',qindex(1))
             qindex=qindex(1)
          endif
          if (nsinks<5) then
             if (iuse_disc(1)) then
                call prompt('Enter H/R of circumbinary at R_ref',H_R(1))
                H_R(2) = (R_ref(2)/R_ref(1)*(m1+m2)/m1)**(0.5-qindex(1)) * H_R(1)
                H_R(3) = (R_ref(3)/R_ref(1)*(m1+m2)/m2)**(0.5-qindex(1)) * H_R(1)
             else
                if (iuse_disc(2)) then
                   call prompt('Enter H/R of circumprimary at R_ref',H_R(2))
                   H_R(1) = (R_ref(1)/R_ref(2)*m1/(m1+m2))**(0.5-qindex(2)) * H_R(2)
                   H_R(3) = (R_ref(3)/R_ref(2)*m2/m1)**(0.5-qindex(2)) * H_R(2)
                else
                   call prompt('Enter H/R of circumsecondary at R_ref',H_R(3))
                   H_R(1) = sqrt(R_ref(1)/R_ref(3)*m2/(m1+m2))**(0.5-qindex(3)) * H_R(3)
                   H_R(2) = sqrt(R_ref(2)/R_ref(3)*m2/m1)**(0.5-qindex(3)) * H_R(3)
                endif
             endif
             !H_R(2) = nint(H_R(2)*10000.)/10000.
             !H_R(3) = nint(H_R(3)*10000.)/10000.
          else
             higher_disc_index = findloc_local(iuse_disc, .true.)
             call get_hier_disc_label(higher_disc_index, disclabel)
             call prompt('Enter H/R of circum-'//trim(disclabel)//' at R_ref',H_R(higher_disc_index))

             higher_mass = get_hier_level_mass(trim(disclabel))!, mass, sink_num, sink_labels)
             !return
             do i=1,maxdiscs
                if (iuse_disc(i) .and. i /= higher_disc_index) then
                   call get_hier_disc_label(i, disclabel)
                   current_mass = get_hier_level_mass(trim(disclabel))
                   H_R(i) = (R_ref(i)/R_ref(higher_disc_index) * &
                        higher_mass/current_mass)**(0.5-qindex(higher_disc_index)) * &
                        H_R(higher_disc_index)
                endif
             enddo
          endif
          do i=1, maxdiscs
             if (iuse_disc(i)) then
                H_R(i) = nint(H_R(i)*10000.)/10000.
             endif
          enddo
       else
          if (iuse_disc(1)) then
             H_R(2) = sqrt(R_ref(2)/R_ref(1)*(m1+m2)/m1) * H_R(1)
             H_R(3) = sqrt(R_ref(3)/R_ref(1)*(m1+m2)/m2) * H_R(1)
             qindex(2) = qindex(1)
             qindex(3) = qindex(1)
             call warning('setup_disc','using circumbinary (H/R)_ref to set global temperature')
          elseif (iuse_disc(2)) then
             H_R(3) = sqrt(R_ref(3)/R_ref(2)*m1/m2) * H_R(2)
             qindex(3) = qindex(2)
             call warning('setup_disc','using circumprimary (H/R)_ref to set global temperature')
          endif
          H_R(2) = nint(H_R(2)*10000.)/10000.
          H_R(3) = nint(H_R(3)*10000.)/10000.
       endif
    endif
 endif
 do i=1,maxdiscs
    if (iuse_disc(i)) then
       if (ndiscs > 1) then
          if (nsinks<5) then
             print "(/,a)",' >>>  circum'//trim(disctype(i))//' disc  <<<'
          elseif (nsinks>4) then
             call get_hier_disc_label(i, disclabel)
             print "(/,a)",' >>>  circum'//trim(disclabel)//' disc  <<<'
          endif
       endif
       call prompt('How do you want to set the gas disc mass?'//new_line('A')// &
                  ' 0=total disc mass'//new_line('A')// &
                  ' 1=mass within annulus'//new_line('A')// &
                  ' 2=surface density normalisation'//new_line('A')// &
                  ' 3=surface density at reference radius'//new_line('A')// &
                  ' 4=minimum Toomre Q'//new_line('A'),isetgas(i),0,4)
       call prompt('Do you want to exponentially taper the outer gas disc profile?',itapergas(i))
       call prompt('Do you want to read the sigma profile from a grid file?',sigma_file(i))
       call prompt('Do you want to warp the disc?',iwarp(i))
       select case (isetgas(i))
       case (0)
          disc_m(i)    = 0.05   * disc_mfac(i)
       case (1)
          annulus_m(i) = 0.05   * disc_mfac(i)
          R_inann(i)   = R_in(i)
          R_outann(i)  = R_out(i)
       case (2)
          sig_norm(i)  = 1.E-02 * disc_mfac(i)
       case (3)
          sig_ref(i)   = 1.E-02 * disc_mfac(i)
       case (4)
          Q_min(i)     = 1.0
       end select
       if (iwarp(i)) then
          R_warp = 0.5*(R_in + R_out)
          H_warp = 20.
          incl   = 30.
       endif

       call prompt('Do you want the disc to be eccentric?',iecc(i))
       if (iecc(i)) then
          e0(i)=0.1
          eindex(i) = 1.
          phiperi(i) = 0.
          eccprofile(i) = 4
       else
          e0(:)=0.
          eindex(:)=0.
          phiperi(:)=0.
          eccprofile(:)=0.
       endif
    endif
 enddo

 !--dust disc
 if (use_dust) then
    print "(/,a)",'=============='
    print "(a)",  '+++  DUST  +++'
    print "(a)",  '=============='
    !--dust distribution
    call set_dust_interactive()
    !--dust discs
    do i=1,maxdusttypes
       R_indust(:,i)    = R_in
       R_outdust(:,i)   = R_out
       qindex_dust(:,i) = qindex
       H_R_dust(:,i)    = H_R
       ismoothdust(:,i) = ismoothgas
       R_c_dust(:,i)    = R_c
    enddo
    !--dust growth
    if (use_dustgrowth .and. dust_method == 2) then
       print "(/,a)",'================================'
       print "(a)",  '+++  GROWTH & FRAGMENTATION  +++'
       print "(a)",  '================================'
       call prompt('Enter fragmentation model (0=off,1=on,2=Kobayashi)',ifrag,-1,2)
       if (ifrag > 0) then
          call prompt('Enter minimum allowed grain size in cm',gsizemincgs)
          call prompt('Do you want a snow line ? (0=no,1=position based,2=temperature based)',isnow,0,2)
          if (isnow == 0) then
             call prompt('Enter uniform vfrag in m/s',vfragSI,1.)
          else
             if (isnow == 1) call prompt('How far from the star in AU ?',rsnow,0.)
             if (isnow == 2) call prompt('Enter snow line condensation temperature in K',Tsnow,0.)
             call prompt('Enter inward vfragin in m/s',vfraginSI,1.)
             call prompt('Enter outward vfragout in m/s',vfragoutSI,1.)
          endif
       endif
       call prompt('Enter porosity switch (0=off,1=on)',iporosity,0,1)
       if (iporosity == 1) use_porosity = .true.
    endif
 endif

 !--resolution
 if (use_dust .and. .not.use_dustfrac) then
    np_dust = np/ndusttypesinp/5
    !elseif (use_dust .and. use_hybrid) then
    !np_dust = np/ndustlargeinp/5
 else
    np_dust = 0
 endif

 !--planets
 print "(/,a)",'================='
 print "(a)",  '+++  PLANETS  +++'
 print "(a)",  '================='
 call prompt('How many planets?',nplanets,0,maxplanets)

 !--simulation time
 print "(/,a)",'================'
 print "(a)",  '+++  OUTPUT  +++'
 print "(a)",  '================'
 if (nplanets > 0) then
    call prompt('Enter time between dumps as fraction of outer planet period',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 elseif (icentral==1 .and. nsinks==2 .and. binary%input_type==0) then
    call prompt('Enter time between dumps as fraction of binary period',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 elseif (icentral==1 .and. nsinks>=3 .and. binary%input_type==0) then
    call prompt('Enter time between dumps as fraction of binary period'//new_line('A')// &
         '(enter a negative number to refer to the shorter period)',deltat)
    call prompt('Enter number of orbits to simulate'//new_line('A')// &
         '(enter a negative number to refer to the shorter period)',norbits)
 elseif (icentral==1 .and. nsinks==2 .and. binary%input_type > 0) then
    deltat  = 0.01
    norbits = 1
    call prompt('Enter time between dumps as fraction of flyby time',deltat,0.)
 elseif (any(iwarp)) then
    call prompt('Enter time between dumps as fraction of orbital time at warp',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 else
    call prompt('Enter time between dumps as fraction of outer disc orbital time',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 endif

end subroutine setup_interactive

!--------------------------------------------------------------------------
!
! Write setup file
!
!--------------------------------------------------------------------------
subroutine write_setupfile(filename)
 use eos,              only:istrat,alpha_z,beta_z,qfacdisc2
 use infile_utils,     only:write_inopt
 use set_dust_options, only:write_dust_setup_options
 use sethierarchical,  only:write_hierarchical_setupfile,hs
 use setorbit,         only:write_options_orbit
 use setunits,         only:write_options_units
 use eos_stamatellos, only:eos_file
 character(len=*), intent(in) :: filename
 logical            :: done_alpha
 integer            :: i,j,n_possible_discs,iunit
 character(len=20)  :: duststring(maxdusttypes)
 character(len=20)  :: taper_string
 character(len=20)  :: smooth_string
 character(len=40)  :: tmpstr

 done_alpha = .false.
 n_possible_discs = 1
 if ((icentral==1) .and. (nsinks>=2)) n_possible_discs = 3

 print "(/,a)",' writing setup options file '//trim(filename)
 open(newunit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for disc setup routine'
 !--resolution
 write(iunit,"(/,a)") '# resolution'
 call write_inopt(np,'np','number of gas particles',iunit)
 if (use_dust .and. (.not.use_dustfrac .or. use_hybrid)) then
    duststring = 'np_dust'
    call make_tags_unique(ndustlargeinp,duststring)
    do i=1,ndustlargeinp
       call write_inopt(np_dust(i),trim(duststring(i)),'number of large dust particles',iunit)
    enddo
 endif
 !--units
 call write_options_units(iunit,gr=(gr .or. icentral==0 .and. ipotential==3))

 !--central objects(s)/potential
 write(iunit,"(/,a)") '# central object(s)/potential'
 call write_inopt(icentral,'icentral', &
    'use sink particles or external potential (0=potential,1=sinks)',iunit)
 select case (icentral)
 case (0)
    !--external potential
    call write_inopt(ipotential,'ipotential','potential (1=central point mass,'// &
                     '2=binary potential,3=spinning black hole)',iunit)
    select case (ipotential)
    case (1)
       !--point mass
       call write_inopt(m1,'m1','star mass',iunit)
       call write_inopt(accr1,'accr1','star accretion radius',iunit)
    case (2)
       !--fixed binary
       call write_inopt(m1,'m1','primary mass',iunit)
       call write_inopt(m2,'m2','secondary mass',iunit)
       call write_inopt(accr1,'accr1','primary accretion radius',iunit)
       call write_inopt(accr2,'accr2','secondary accretion radius',iunit)

       !--options of planet surface/atmosphere
       write(iunit,"(/,a)") '# options for planet surface/atmosphere'
       call write_inopt(surface_force,'surface_force','model m1 as planet with surface',iunit)
       if (surface_force) then
          call write_inopt(rho_core_cgs,'rho_core','planet core density (cgs units)',iunit)
          call write_inopt(Ratm_in,'Ratm_in','inner atmosphere radius (planet radii)',iunit)
          call write_inopt(Ratm_out,'Ratm_out','outer atmosphere radius (planet radii)',iunit)
          call write_inopt(atm_type,'atm_type','atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))',iunit)
          call write_inopt(Natmfrac,'Natm/Npart','fraction of particles for planet atmosphere',iunit)
          call write_inopt(ramp,'ramp','Do you want to ramp up the planet mass slowly?',iunit)
          if (.not.ramp .and. Natmfrac == 0.) then
             print*,'Warning! Not ramping the mass or initialising an atmosphere will'// &
                     'likely cause explosive collisions between particles'
          elseif (ramp .and. Natmfrac /= 0.) then
             print*,'Warning! The atmosphere will be lost while ramping up the planet mass...'
             print*,'         ...try using one or the other'
          endif
       endif
    case (3)
       !--spinning black hole: Lense-Thirring (+ Einstein precession)
       call write_inopt(einst_prec,'einst_prec','include Einstein precession',iunit)
       call write_inopt(m1,'m1','black hole mass',iunit)
       call write_inopt(accr1,'accr1','black hole accretion radius',iunit)
       call write_inopt(bhspin,'bhspin','black hole spin',iunit)
       call write_inopt(bhspinangle,'bhspinangle','black hole spin angle (deg)',iunit)
    end select
 case (1)
    !--sink particle(s)
    call write_inopt(nsinks,'nsinks','number of sinks',iunit)
    select case (nsinks)
    case (1)
       !--single star
       write(iunit,"(/,a)") '# options for central star'
       call write_inopt(m1,'m1','star mass',iunit)
       call write_inopt(accr1,'accr1','star accretion radius',iunit)
    case (2)
       !--binary
       select case (binary%input_type)
       case (1)
          !--unbound (flyby)
          !--central star
          write(iunit,"(/,a)") '# options for central star'
          call write_inopt(m1,'m1','central star mass',iunit)
          call write_inopt(accr1,'accr1','central star accretion radius',iunit)
          !--perturber
          write(iunit,"(/,a)") '# options for perturber'
          call write_inopt(m2,'m2','perturber mass',iunit)
          call write_inopt(accr2,'accr2','perturber accretion radius',iunit)
       case default
          !--bound
          write(iunit,"(/,a)") '# options for binary'
          call write_inopt(m1,'m1','primary mass',iunit)
          call write_inopt(m2,'m2','secondary mass',iunit)
          call write_inopt(accr1,'accr1','primary accretion radius',iunit)
          call write_inopt(accr2,'accr2','secondary accretion radius',iunit)
       end select
       call write_options_orbit(binary,iunit,prefix='binary')

    case (5:)

       call write_hierarchical_setupfile(iunit)

    case (3)
       !-- hierarchical triple
       write(iunit,"(/,a)") '# options for hierarchical triple'

       !-- masses
       call write_inopt(m1,'m1','first hierarchical level primary mass',iunit)
       call write_inopt(m2,'m2','first hierarchical level secondary mass',iunit)
       call write_inopt(q2,'q2','tight binary mass ratio',iunit)
       call write_inopt(subst,'subst','star to substitute',iunit)

       !-- wide binary parameters
       call write_options_orbit(binary,iunit,prefix='binary',comment_prefix='wide binary',input_type=0)

       !-- tight parameters
       call write_options_orbit(binary2,iunit,prefix='binary2',comment_prefix='tight binary',input_type=0)

       !-- accretion radii
       call write_inopt(accr1,'accr1','single star accretion radius',iunit)
       call write_inopt(accr2a,'accr2a','tight binary primary accretion radius',iunit)
       call write_inopt(accr2b,'accr2b','tight binary secondary accretion radius',iunit)
    case(4)
       !-- hierarchical quadruple
       write(iunit,"(/,a)") '# options for hierarchical quadruple'

       !-- masses
       call write_inopt(m1,'m1','first hierarchical level primary mass',iunit)
       call write_inopt(m2,'m2','first hierarchical level secondary mass',iunit)
       call write_inopt(q1,'q1','tight binary 1 mass ratio',iunit)
       call write_inopt(q2,'q2','tight binary 2 mass ratio',iunit)
       call write_inopt(subst1,'subst1','first star to substitute',iunit)
       call write_inopt(subst2,'subst2','second star to substitute',iunit)

       !-- wide binary parameters
       call write_options_orbit(binary,iunit,prefix='binary',comment_prefix='wide binary',input_type=0)

       !-- tight binary 1 parameters
       call write_options_orbit(binary1,iunit,prefix='binary1',comment_prefix='tight binary 1',input_type=0)

       !-- tight binary 2 parameters
       call write_options_orbit(binary2,iunit,prefix='binary2',comment_prefix='tight binary 2',input_type=0)

       !-- accretion radii
       call write_inopt(accr1a,'accr1a','single star accretion radius',iunit)
       call write_inopt(accr1b,'accr1b','single star accretion radius',iunit)
       call write_inopt(accr2a,'accr2a','tight binary primary accretion radius',iunit)
       call write_inopt(accr2b,'accr2b','tight binary secondary accretion radius',iunit)

    end select

    !--options for oblateness
    write(iunit,"(/,a)") '# oblateness'
    do i=1,nsinks
       call write_oblateness_options(iunit,'_body'//trim(num(i)), &
            J2star(i),size_star(i),spin_period_star(i),kfac_star(i),obliquity_star(i))
    enddo

 end select
 !--multiple disc options
 if (n_possible_discs > 1) then
    if (nsinks == 2) then
       write(iunit,"(/,a)") '# options for multiple discs'
       do i=1,3!maxdiscs-1
          call write_inopt(iuse_disc(i),'use_'//trim(disctype(i))//'disc','setup circum' &
               //trim(disctype(i))//' disc',iunit)
       enddo
    elseif (nsinks == 3) then
       write(iunit,"(/,a)") '# options for multiple discs'
       call write_inopt(iuse_disc(1),'use_'//trim(disctype(1))//'disc','setup circum' &
            //trim(disctype(1))//' disc',iunit)
       call write_inopt(iuse_disc(4),'use_'//trim(disctype(4))//'disc','setup circum' &
            //trim(disctype(4))//' disc',iunit)
    elseif (nsinks == 4) then
       write(iunit,"(/,a)") '# options for multiple discs'
       call write_inopt(iuse_disc(1),'use_'//trim(disctype(1))//'disc','setup circum' &
            //trim(disctype(1))//' disc',iunit)
    elseif (nsinks >= 5) then
       write(iunit,"(/,a)") '# options for multiple discs'

       do i=1,hs%labels%sink_num
          call write_inopt(iuse_disc(i),'use_'//trim(hs%labels%sink(i))//'disc','setup circum-' &
            //trim(hs%labels%sink(i))//' disc',iunit)
       enddo

       do i=1,hs%labels%hl_num
          call write_inopt(iuse_disc(i+hs%labels%sink_num),'use_'//trim(hs%labels%hl(i))//'disc','setup circum-' &
            //trim(hs%labels%hl(i))//' disc',iunit)
       enddo
    endif
    call write_inopt(use_global_iso,'use_global_iso',&
        'globally isothermal or Farris et al. (2014)',iunit)

 endif
 !--individual disc(s)
 do i=1,maxdiscs
    if (iuse_disc(i)) then
       if (n_possible_discs > 1) then
          disclabel = disctype(i)
          if (nsinks > 4) then
             call get_hier_disc_label(i, disclabel)
          endif
       else
          disclabel = ''
       endif
       !--gas disc
       if (n_possible_discs > 1) then
          write(iunit,"(/,a)") '# options for circum-'//trim(disclabel)//' gas disc'
       else
          write(iunit,"(/,a)") '# options for gas accretion disc'
       endif
       call write_inopt(isetgas(i),'isetgas'//trim(disclabel),'how to set gas density profile' // &
          ' (0=total disc mass,1=mass within annulus,2=surface density normalisation,' // &
          '3=surface density at reference radius,4=minimum Toomre Q,5=minimum Toomre Q and Lstar)',iunit)
       call write_inopt(sigma_file(i),'sigma_file'//trim(disclabel), &
           'reading gas profile from file sigma_grid.dat',iunit)
       call write_inopt(itapergas(i),'itapergas'//trim(disclabel), &
          'exponentially taper the outer disc profile',iunit)
       if (itapergas(i)) call write_inopt(itapersetgas(i),'itapersetgas'//trim(disclabel), &
          'how to set taper (0=exp[-(R/R_c)^(2-p)], 1=[1-exp(R-R_out)]',iunit)
       call write_inopt(ismoothgas(i),'ismoothgas'//trim(disclabel),'smooth inner disc',iunit)
       call write_inopt(iwarp(i),'iwarp'//trim(disclabel),'warp disc',iunit)
       call write_inopt(iecc(i),'iecc'//trim(disclabel),'eccentric disc',iunit)
       call write_inopt(R_in(i),'R_in'//trim(disclabel),'inner radius',iunit)
       call write_inopt(R_ref(i),'R_ref'//trim(disclabel),'reference radius',iunit)
       call write_inopt(R_out(i),'R_out'//trim(disclabel),'outer radius',iunit)
       if (itapergas(i) .and. itapersetgas(i)==0) call write_inopt(R_c(i),'R_c'//trim(disclabel), &
          'characteristic radius of the exponential taper',iunit)
       select case (isetgas(i))
       case (0)
          call write_inopt(disc_m(i),'disc_m'//trim(disclabel),'disc mass',iunit)
       case (1)
          call write_inopt(annulus_m(i),'annulus_m'//trim(disclabel),'mass within annulus',iunit)
          call write_inopt(R_inann(i),'R_inann'//trim(disclabel),'inner annulus radius',iunit)
          call write_inopt(R_outann(i),'R_outann'//trim(disclabel),'outer annulus radius',iunit)
       case (2)
          taper_string = ''
          smooth_string = ''
          if (itapergas(i)) then
             if (itapersetgas(i)==0) then
                taper_string = 'exp[-(R/R_c)^(2-p)]'
             elseif (itapersetgas(i)==1) then
                taper_string = '[1-exp(R-R_out)]'
             endif
          endif
          if (ismoothgas(i)) then
             smooth_string = '(1-sqrt(R_in/R))'
          endif
          call write_inopt(sig_norm(i),'sig_norm'//trim(disclabel), &
             'sigma = sig_norm (R/R_ref)^-p '//trim(taper_string)//' '//trim(smooth_string),iunit)
       case (3)
          call write_inopt(sig_ref(i),'sig_ref'//trim(disclabel),'sigma at reference radius',iunit)
       case (4)
          call write_inopt(Q_min(i),'Q_min'//trim(disclabel),'minimum Toomre Q',iunit)

       end select
       call write_inopt(lumdisc,'lumdisc', 'Set qindex from stellar luminosity (ieos=24) (0=no 1=yes)',iunit)
       if (lumdisc > 0) then
          call write_inopt(L_star(i),'L_star'//trim(disclabel),'Stellar luminosity (Lsun)',iunit)
          call write_inopt(T_bg,'T_bg'//trim(disclabel),'background Temperature (K)',iunit)
       endif
       call write_inopt(pindex(i),'pindex'//trim(disclabel),'power law index of surface density sig=sig0*r^-p',iunit)
       if (lumdisc == 0) then
          call write_inopt(qindex(i),'qindex'//trim(disclabel),'power law index of sound speed cs=cs0*r^-q',iunit)
       endif
       call write_inopt(posangl(i),'posangl'//trim(disclabel),'position angle (deg)',iunit)
       call write_inopt(incl(i),'incl'//trim(disclabel),'inclination (deg)',iunit)
       if (discstrat == 0 .and. lumdisc == 0) then
          call write_inopt(H_R(i),'H_R'//trim(disclabel),'H/R at R=R_ref',iunit)
       endif
       if (iwarp(i)) then
          call write_inopt(R_warp(i),'R_warp'//trim(disclabel),'warp radius',iunit)
          call write_inopt(H_warp(i),'H_warp'//trim(disclabel),'warp smoothing length',iunit)
       endif
       if (iecc(i)) then
          call write_inopt(e0(i),'e0'//trim(disclabel),'eccentricity at disc edge',iunit)
          call write_inopt(eindex(i),'eindex'//trim(disclabel),'power of eccentricity profile',iunit)
          call write_inopt(phiperi(i),'phiperi'//trim(disclabel),'longitude of pericentre',iunit)
          call write_inopt(eccprofile(i),'eccprofile'//trim(disclabel),'type of eccentricity profile'// &
                           '(0=const e0,1=power-law,4=from file ecc_grid.dat)',iunit)

       endif
       if (.not.done_alpha) then
          call write_inopt(alphaSS,'alphaSS','desired alphaSS (0 for minimal needed for shock capturing)',iunit)
          done_alpha = .true.
       endif
       !--dust disc
       if (use_dust .and. (isetdust == 1 .or. isetdust == 2)) then
          duststring = 'dust'
          call make_tags_unique(ndusttypesinp,duststring)
          do j=1,ndusttypesinp
             if (n_possible_discs > 1) then
                write(iunit,"(/,a)") '# options for circum'//trim(disclabel)//' '//trim(duststring(j))//' disc'
             else
                write(iunit,"(/,a)") '# options for '//trim(duststring(j))//' accretion disc'
             endif
             tmpstr = trim(duststring(j))//trim(disclabel)
             call write_inopt(itaperdust(i,j),'itaper'//trim(tmpstr), &
                'exponentially taper the outer disc profile',iunit)
             if (itaperdust(i,j)) call write_inopt(itapersetdust(i,j),'itapersetdust'//trim(tmpstr), &
                'how to set taper (0=exp[-(R/R_c)^(2-p)], 1=[1-exp(R-R_out)]',iunit)
             call write_inopt(ismoothdust(i,j),'ismooth'//trim(tmpstr),'smooth inner disc',iunit)
             call write_inopt(R_indust(i,j),'R_in'//trim(tmpstr),'inner radius',iunit)
             call write_inopt(R_outdust(i,j),'R_out'//trim(tmpstr),'outer radius',iunit)
             if (itaperdust(i,j) .and. itapersetdust(i,j)==0) then
                call write_inopt(R_c_dust(i,j),'R_c_'//trim(tmpstr), &
                'characteristic radius of the exponential taper',iunit)
             endif
             call write_inopt(pindex_dust(i,j),'pindex_'//trim(tmpstr),'p index',iunit)
             call write_inopt(qindex_dust(i,j),'qindex_'//trim(tmpstr),'q index',iunit)
             call write_inopt(H_R_dust(i,j),'H_R_'//trim(tmpstr),'H/R at R=R_ref',iunit)
          enddo
       endif
    endif
 enddo
 if (lumdisc > 0) call write_inopt(eos_file,'eos_file','Equation of state file for using lumdisc',iunit)
 !--dust & growth options
 if (use_dust) then
    call write_dust_setup_options(iunit)
 endif
 !-- minimum temperature
 write(iunit,"(/,a)") '# Minimum Temperature in the Simulation'
 call write_inopt(T_floor,'T_floor','The minimum temperature in the simulation (for any locally isothermal EOS).',iunit)
 !--sphere of gas around disc
 write(iunit,"(/,a)") '# set sphere around disc'
 call write_inopt(add_sphere,'add_sphere','add sphere around disc?',iunit)
 if (add_sphere) then
    call write_inopt(mass_sphere,'mass_sphere','Mass of sphere',iunit)
    call write_inopt(Rin_sphere,'Rin_sphere','Inner edge of sphere',iunit)
    call write_inopt(Rout_sphere,'Rout_sphere','Outer edge of sphere',iunit)
    call write_inopt(add_rotation,'add_rotation', &
      'Rotational Velocity of the cloud (0=no rotation, 1=k*(GM/R^3)^0.5, '// &
      '2=Omega (s^-1))',iunit)
    if (add_rotation==1) then
       call write_inopt(Kep_factor,'k','Scaling factor of Keplerian rotational velocity',iunit)
       call write_inopt(R_rot,'R_rot','Set rotational velocity as Keplerian velocity at R=R_rot',iunit)
    elseif (add_rotation==2) then
       call write_inopt(omega_cloud,'omega_cloud','Rotational velocity of the cloud (s^-1)',iunit)
    endif
    call write_inopt(add_turbulence,'add_turbulence','Add turbulence to the sphere (0=no turbulence, 1=turbulence)',iunit)
    if (add_turbulence==1) then
       call write_inopt(rms_mach,'rms_mach','RMS Mach number of turbulence',iunit)
       call write_inopt(tfact,'tfact','Scale the maximum length scale of the turbulence',iunit)
    endif
    call write_inopt(set_freefall,'set_freefall','Set the sphere in freefall (0=no freefall, 1=freefall)',iunit)
    if (use_dust) then
       if (use_dustfrac) then
          call write_inopt(dustfrac_method,'dustfrac_method',&
                        'How to set the dustfrac in the cloud? (-1=no dust, 0=global ratio, 1=bin ratio)',iunit)
       endif
    endif
 endif

 !--planets
 write(iunit,"(/,a)") '# set planets'
 call write_inopt(nplanets,'nplanets','number of planets',iunit)
 if (nplanets > 0) then
    do i=1,nplanets
       write(iunit,"(/,a)") '# planet:'//trim(num(i))
       call write_inopt(mplanet(i),'mplanet'//trim(num(i)),'planet mass (in Jupiter mass)',iunit)
       call write_inopt(rplanet(i),'rplanet'//trim(num(i)),'planet distance from star',iunit)
       call write_inopt(inclplan(i),'inclplanet'//trim(num(i)),'planet orbital inclination (deg)',iunit)
       call write_inopt(accrplanet(i),'accrplanet'//trim(num(i)),'planet accretion radius (in Hill radius)',iunit)
       call write_inopt(eccplanet(i),'eccplanet'//trim(num(i)),'planet eccentricity',iunit)
       call write_inopt(Oplanet(i),'Oplanet'//trim(num(i)),'position angle of ascending node (deg)',iunit)
       call write_inopt(wplanet(i),'wplanet'//trim(num(i)),'argument of periapsis (deg)',iunit)
       call write_inopt(fplanet(i),'fplanet'//trim(num(i)),'true anomaly (deg) 180=apastron',iunit)
       call write_oblateness_options(iunit,'_planet'//trim(num(i)), &
            J2planet(i),planet_size(i),spin_period(i),kfac(i),obliquity(i))
    enddo
 endif
 ! stratification
 write(iunit,"(/,a)") '# thermal stratification'
 call write_inopt(discstrat,'discstrat','stratify disc? (0=no,1=yes)',iunit)
 if (discstrat==1) then
    call write_inopt(istrat,'istrat','temperature prescription (0=MAPS, 1=Dartois)',iunit)
    call write_inopt(z0_ref,'z0', 'z scaling factor',iunit)
    call write_inopt(alpha_z,'alpha_z', 'height of transition in tanh vertical temperature profile',iunit)
    call write_inopt(beta_z,'beta_z', 'variation in transition height over radius',iunit)
    call write_inopt(temp_mid0,'temp_mid0', 'midplane temperature scaling factor',iunit)
    call write_inopt(temp_atm0,'temp_atm0', 'atmosphere temperature scaling factor',iunit)
    call write_inopt(qfacdisc2,'qatm', 'sound speed power law index of atmosphere',iunit)

 endif
 !--timestepping
 write(iunit,"(/,a)") '# timestepping'
 if (nplanets > 0) then
    call write_inopt(norbits,'norbits','maximum number of outer planet orbits',iunit)
 elseif (icentral==1 .and. nsinks>=2 .and. binary%input_type==0) then
    call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
 elseif (icentral==1 .and. nsinks==2 .and. binary%input_type > 0) then
    call write_inopt(norbits,'norbits','maximum number of flyby times',iunit)
 else
    call write_inopt(norbits,'norbits','maximum number of orbits at outer disc',iunit)
 endif
 call write_inopt(deltat,'deltat','output interval as fraction of orbital period',iunit)

 !--mcfost
 if (compiled_with_mcfost) then
    write(iunit,"(/,a)") '# mcfost'
    call write_inopt(use_mcfost,'use_mcfost','use the mcfost library',iunit)
    call write_inopt(use_mcfost_stellar_parameters,'use_mcfost_stars',&
        'Fix the stellar parameters to mcfost values or update using sink mass',iunit)
 endif

 if (do_radiation) call write_inopt(iradkappa,'radkappa','constant radiation opacity kappa',iunit)

 close(iunit)

end subroutine write_setupfile

!--------------------------------------------------------------------------
!
! Read setup file
!
!--------------------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use eos,              only:istrat,alpha_z,beta_z,qfacdisc2
 use dust,             only:ilimitdustflux
 use infile_utils,     only:open_db_from_file,inopts,read_inopt,close_db
 use set_dust_options, only:read_dust_setup_options,ilimitdustfluxinp
 use sethierarchical,  only:read_hierarchical_setupfile,hs
 use setorbit,         only:read_options_orbit
 use setunits,         only:read_options_and_set_units
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr

 type(inopts), allocatable :: db(:)
 integer,      parameter   :: iunit = 21
 integer                   :: nerr,i,j
 character(len=20)         :: duststring(maxdusttypes)
 character(len=40)         :: tmpstr

 print "(a)",' reading setup options from '//trim(filename)

 call open_db_from_file(db,filename,iunit,ierr)

 nerr = 0

 !--units
 call read_options_and_set_units(db,nerr,gr)

 !--central objects(s)/potential
 call read_inopt(icentral,'icentral',db,min=0,max=1,errcount=nerr)
 select case (icentral)
 case (0)
    !--external potential
    call read_inopt(ipotential,'ipotential',db,min=1,max=3,errcount=nerr)
    select case (ipotential)
    case (1)
       !--point mass
       iexternalforce = iext_star
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
    case (2)
       !--fixed binary
       iexternalforce = iext_binary
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
       call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)

       !--options of planet surface/atmosphere
       call read_inopt(surface_force,'surface_force',db,errcount=nerr)
       if (surface_force) then
          iexternalforce = iext_corot_binary
          call read_inopt(rho_core_cgs,'rho_core',db,min=0.,errcount=nerr)
          call read_inopt(Ratm_in,'Ratm_in',db,min=1.,errcount=nerr)
          call read_inopt(Ratm_out,'Ratm_out',db,min=1.,max=10.,errcount=nerr)
          call read_inopt(atm_type,'atm_type',db,min=1,max=2,errcount=nerr)
          call read_inopt(Natmfrac,'Natm/Npart',db,min=0.,max=1.,errcount=nerr)
          call read_inopt(ramp,'ramp',db,errcount=nerr)
       endif
    case (3)
       ! force use of geometric units
       call read_options_and_set_units(db,nerr,gr=.true.)

       !--spinning black hole (Lense-Thirring)
       iexternalforce = iext_lensethirring
       call read_inopt(einst_prec,'einst_prec',db,errcount=nerr)
       if (einst_prec) iexternalforce = iext_einsteinprec
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
       call read_inopt(bhspin,'bhspin',db,min=0.,errcount=nerr)
       call read_inopt(bhspinangle,'bhspinangle',db,min=0.,errcount=nerr)
    end select
 case (1)
    iexternalforce = 0
    !--sink particles
    call read_inopt(nsinks,'nsinks',db,min=1,errcount=nerr)
    select case (nsinks)
    case (1)
       !--single star
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
    case (2)
       !--binary
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
       call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
       call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
       call read_options_orbit(binary,m1,m2,db,nerr,prefix='binary')
    case (5:)

       call read_hierarchical_setupfile(db, nerr)

    case (3)
       !-- hierarchical triple

       !-- masses
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
       call read_inopt(q2,'q2',db,min=0.,max=1.,errcount=nerr)
       call read_inopt(subst,'subst',db,errcount=nerr)

       !-- wide binary parameters
       call read_options_orbit(binary,m1,m2,db,nerr,prefix='binary',input_type=0)

       !-- tight parameters
       call read_options_orbit(binary2,m2/(q2+1.),m2*q2/(q2+1.),db,nerr,prefix='binary2',input_type=0)

       !-- accretion radii
       call read_inopt(accr1,'accr1',db,errcount=nerr)
       call read_inopt(accr2a,'accr2a',db,errcount=nerr)
       call read_inopt(accr2b,'accr2b',db,errcount=nerr)
    case(4)
       !-- hierarchical quadruple

       !-- masses
       call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
       call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
       call read_inopt(q1,'q1',db,min=0.,max=1.,errcount=nerr)
       call read_inopt(q2,'q2',db,min=0.,max=1.,errcount=nerr)
       call read_inopt(subst1,'subst1',db,errcount=nerr)
       call read_inopt(subst2,'subst2',db,errcount=nerr)

       !-- wide binary parameters
       call read_options_orbit(binary,m1,m2,db,nerr,prefix='binary',input_type=0)

       !-- tight binary 1 parameters
       call read_options_orbit(binary1,m1/(q1+1.),m1*q1/(q1+1.),db,nerr,prefix='binary1',input_type=0)

       !-- tight binary 2 parameters
       call read_options_orbit(binary2,m2/(q2+1.),m2*q2/(q2+1.),db,nerr,prefix='binary2',input_type=0)

       !-- accretion radii
       call read_inopt(accr1a,'accr1a',db,errcount=nerr)
       call read_inopt(accr1b,'accr1b',db,errcount=nerr)
       call read_inopt(accr2a,'accr2a',db,errcount=nerr)
       call read_inopt(accr2b,'accr2b',db,errcount=nerr)

    end select
    do i=1,nsinks
       call read_oblateness_options(db,nerr,'_body'//trim(num(i)),&
            J2star(i),size_star(i),spin_period_star(i),kfac_star(i),obliquity_star(i))
    enddo
 end select

 call read_inopt(T_floor,'T_floor',db,errcount=nerr)

 call read_inopt(discstrat,'discstrat',db,errcount=nerr)
 call read_inopt(lumdisc,'lumdisc',db,default=0)

 if (discstrat==1) then
    call read_inopt(istrat,'istrat',db,errcount=nerr)
    call read_inopt(z0_ref,'z0',db,errcount=nerr)
    call read_inopt(alpha_z,'alpha_z',db,errcount=nerr)
    call read_inopt(beta_z,'beta_z',db,errcount=nerr)
    call read_inopt(temp_mid0,'temp_mid0',db,errcount=nerr)
    call read_inopt(temp_atm0,'temp_atm0',db,errcount=nerr)
    call read_inopt(qfacdisc2,'qatm',db,errcount=nerr)
 endif

 !--dust
 if (use_dust) then
    call read_dust_setup_options(db,nerr)
    !--dust method
    select case(dust_method)
    case(1)
       use_dustfrac = .true.
       ilimitdustflux = ilimitdustfluxinp
       ndustsmall = ndusttypesinp
    case(2)
       use_dustfrac = .false.
       ndustlarge = ndusttypesinp
    case(3)
       use_dustfrac   = .true.
       use_hybrid     = .true.
       ndustlarge     = ndustlargeinp
       ndustsmall     = ndustsmallinp
       ilimitdustflux = ilimitdustfluxinp
    end select
    ndusttypes = ndusttypesinp
 endif

 !--resolution
 call read_inopt(np,'np',db,min=0,errcount=nerr)
 if (use_dust .and. (.not.use_dustfrac .or. use_hybrid)) then
    duststring = 'np_dust'
    call make_tags_unique(ndustlargeinp,duststring)
    do i=1,ndustlargeinp
       call read_inopt(np_dust(i),duststring(i),db,min=0,errcount=nerr)
    enddo
 endif

 !--multiple discs
 iuse_disc = .false.
 if ((icentral==1) .and. (nsinks>=2)) then
    if (nsinks==2) then
       if (binary%input_type==0) then
          call read_inopt(iuse_disc(1),'use_binarydisc',db,errcount=nerr)
       endif
       call read_inopt(iuse_disc(2),'use_primarydisc',db,errcount=nerr)
       call read_inopt(iuse_disc(3),'use_secondarydisc',db,errcount=nerr)
    elseif (nsinks == 3) then
       call read_inopt(iuse_disc(4),'use_tripledisc',db,errcount=nerr)
       call read_inopt(iuse_disc(1),'use_binarydisc',db,errcount=nerr)
    elseif (nsinks == 4) then
       call read_inopt(iuse_disc(1),'use_binarydisc',db,errcount=nerr)
    elseif (nsinks >= 5) then
       do i=1,hs%labels%sink_num
          call read_inopt(iuse_disc(i),'use_'//trim(hs%labels%sink(i))//'disc',db,errcount=nerr)
       enddo

       do i=1,hs%labels%hl_num
          call read_inopt(iuse_disc(i+hs%labels%sink_num),'use_'//trim(hs%labels%hl(i))//'disc',db,errcount=nerr)
       enddo
    endif
 else
    iuse_disc(1) = .true.
 endif
 ndiscs = count(iuse_disc)
 if (ndiscs > 1) then
    call read_inopt(use_global_iso,'use_global_iso',db,errcount=nerr)
 endif

 do i=1,maxdiscs
    if (iuse_disc(i)) then
       if (nsinks >= 2) then
          disclabel = disctype(i)
          if (nsinks > 4) then
             call get_hier_disc_label(i, disclabel)
          endif
       else
          disclabel = ''
       endif
       !--gas disc
       call read_inopt(R_in(i),'R_in'//trim(disclabel),db,min=0.,errcount=nerr)
       call read_inopt(R_out(i),'R_out'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       call read_inopt(R_ref(i),'R_ref'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       call read_inopt(sigma_file(i),'sigma_file'//trim(disclabel),db,errcount=nerr)
       call read_inopt(itapergas(i),'itapergas'//trim(disclabel),db,errcount=nerr)
       if (itapergas(i)) call read_inopt(itapersetgas(i),'itapersetgas'//trim(disclabel),db,errcount=nerr)
       call read_inopt(ismoothgas(i),'ismoothgas'//trim(disclabel),db,errcount=nerr)
       call read_inopt(isetgas(i),'isetgas'//trim(disclabel),db,min=0,max=4,errcount=nerr)
       if (itapergas(i)) then
          if (itapersetgas(i)==0) then
             call read_inopt(R_c(i),'R_c'//trim(disclabel),db,min=0.,errcount=nerr)
          endif
       endif
       select case (isetgas(i))
       case (0)
          call read_inopt(disc_m(i),'disc_m'//trim(disclabel),db,min=0.,errcount=nerr)
       case (1)
          call read_inopt(annulus_m(i),'annulus_m'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(R_inann(i),'R_inann'//trim(disclabel),db,min=R_in(i),errcount=nerr)
          call read_inopt(R_outann(i),'R_outann'//trim(disclabel),db,min=R_in(i),errcount=nerr)
       case (2)
          call read_inopt(sig_norm(i),'sig_norm'//trim(disclabel),db,min=0.,errcount=nerr)
       case (3)
          call read_inopt(sig_ref(i),'sig_ref'//trim(disclabel),db,min=0.,errcount=nerr)
       case (4)
          call read_inopt(Q_min(i),'Q_min'//trim(disclabel),db,min=0.,errcount=nerr)
       end select
       call read_inopt(pindex(i),'pindex'//trim(disclabel),db,errcount=nerr)
       if (lumdisc == 0) call read_inopt(qindex(i),'qindex'//trim(disclabel),db,errcount=nerr)
       call read_inopt(posangl(i),'posangl'//trim(disclabel),db,min=-360.,max=360.,errcount=nerr)
       call read_inopt(incl(i),'incl'//trim(disclabel),db,errcount=nerr)
       if (discstrat == 0 .and. lumdisc == 0) then
          call read_inopt(H_R(i),'H_R'//trim(disclabel),db,min=0.,errcount=nerr)
       endif
       call read_inopt(iwarp(i),'iwarp'//trim(disclabel),db,errcount=nerr)
       call read_inopt(iecc(i),'iecc'//trim(disclabel),db,errcount=nerr)
       if (iwarp(i)) then
          call read_inopt(R_warp(i),'R_warp'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(H_warp(i),'H_warp'//trim(disclabel),db,min=0.,errcount=nerr)
       endif
       if (iecc(i)) then
          call read_inopt(e0(i),'e0'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(eindex(i),'eindex'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(phiperi(i),'phiperi'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(eccprofile(i),'eccprofile'//trim(disclabel),db,min=0,max=4,errcount=nerr)
       endif
       !--dust disc
       if (use_dust) then
          call read_inopt(isetdust,'isetdust',db,errcount=nerr)
          duststring = 'dust'
          call make_tags_unique(ndusttypes,duststring)
          do j=1,ndusttypes
             select case (isetdust)
             case (0)
                R_indust(i,j)    = R_in(i)
                R_outdust(i,j)   = R_out(i)
                pindex_dust(i,j) = pindex(i)
                qindex_dust(i,j) = qindex(i)
                H_R_dust(i,j)    = H_R(i)
                use_sigmadust_file(i,j) = sigma_file(i)
                itaperdust(i,j)  = itapergas(i)
                ismoothdust(i,j) = ismoothgas(i)
                R_c_dust(i,j)    = R_c(i)
             case (1,2)
                tmpstr = trim(duststring(j))//trim(disclabel)
                call read_inopt(R_indust(i,j),'R_in'//trim(tmpstr),db,min=R_in(i),err=ierr,errcount=nerr)
                if (ierr /= 0) R_indust(i,j) = R_in(i)

                call read_inopt(R_outdust(i,j),'R_out'//trim(tmpstr),db,min=R_indust(i,j),max=R_out(i),err=ierr,errcount=nerr)
                if (ierr /= 0) R_outdust(i,j) = R_out(i)
                call read_inopt(pindex_dust(i,j),'pindex_'//trim(tmpstr),db,err=ierr,errcount=nerr)
                if (ierr /= 0) pindex_dust(i,j) = pindex(i)
                !call read_inopt(use_sigmadust_file(i,j),'use_sigmadust_file'//trim(tmpstr),db,err=ierr,errcount=nerr)
                call read_inopt(itaperdust(i,j),'itaper'//trim(tmpstr),db,err=ierr,errcount=nerr)
                if (itaperdust(i,j)) call read_inopt(itapersetdust(i,j),'itapersetdust'//trim(tmpstr),db,errcount=nerr)
                call read_inopt(ismoothdust(i,j),'ismooth'//trim(tmpstr),db,err=ierr,errcount=nerr)
                if (itaperdust(i,j)) then
                   if (itapersetdust(i,j)==0) then
                      call read_inopt(R_c_dust(i,j),'R_c_'//trim(tmpstr),db,min=0.,err=ierr,errcount=nerr)
                   endif
                   if (ierr /= 0) R_c_dust(i,j) = R_c(i)
                endif
                call read_inopt(qindex_dust(i,j),'qindex_'//trim(tmpstr),db,min=qindex(i),err=ierr,errcount=nerr)
                if (ierr /= 0) qindex_dust(i,j) = qindex(i)
                call read_inopt(H_R_dust(i,j),'H_R_'//trim(tmpstr),db,min=0.,max=H_R(i),err=ierr,errcount=nerr)
                if (ierr /= 0) H_R_dust(i,j) = H_R(i)
             end select
          enddo
       endif
    endif
 enddo
 if (any(iuse_disc)) call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)

 !--sphere around disc
 call read_inopt(add_sphere,'add_sphere',db,errcount=nerr)
 if (add_sphere) then
    call read_inopt(mass_sphere,'mass_sphere',db,errcount=nerr)
    call read_inopt(Rin_sphere,'Rin_sphere',db,errcount=nerr)
    call read_inopt(Rout_sphere,'Rout_sphere',db,errcount=nerr)
    call read_inopt(add_rotation,'add_rotation',db,errcount=nerr)
    if (add_rotation==1) then
       call read_inopt(Kep_factor,'k',db,errcount=nerr)
       call read_inopt(R_rot,'R_rot',db,errcount=nerr)
    elseif (add_rotation==2) then
       call read_inopt(omega_cloud,'omega_cloud',db,errcount=nerr)
    endif
    call read_inopt(add_turbulence,'add_turbulence',db,errcount=nerr)
    if (add_turbulence==1) then
       call read_inopt(rms_mach,'rms_mach',db,errcount=nerr)
       call read_inopt(tfact,'tfact',db,errcount=nerr)
    endif
    call read_inopt(set_freefall,'set_freefall',db,errcount=nerr)
    if (use_dust) then
       if (use_dustfrac) then
          call read_inopt(dustfrac_method,'dustfrac_method',db,errcount=nerr)
       endif
    endif
 endif

 !--planets
 call read_inopt(nplanets,'nplanets',db,min=0,max=maxplanets,errcount=nerr)
 do i=1,nplanets
    call read_inopt(mplanet(i),'mplanet'//trim(num(i)),db,min=0.,errcount=nerr)
    call read_inopt(rplanet(i),'rplanet'//trim(num(i)),db,min=0.,errcount=nerr)
    call read_inopt(inclplan(i),'inclplanet'//trim(num(i)),db,min=0.,max=180.,errcount=nerr)
    call read_inopt(accrplanet(i),'accrplanet'//trim(num(i)),db,min=0.,errcount=nerr)
    call read_inopt(eccplanet(i),'eccplanet'//trim(num(i)),db,min=0.,errcount=nerr)
    call read_inopt(Oplanet(i),'Oplanet'//trim(num(i)),db,errcount=nerr)
    call read_inopt(wplanet(i),'wplanet'//trim(num(i)),db,errcount=nerr)
    call read_inopt(fplanet(i),'fplanet'//trim(num(i)),db,errcount=nerr)
    call read_oblateness_options(db,nerr,'_planet'//trim(num(i)),&
         J2planet(i),planet_size(i),spin_period(i),kfac(i),obliquity(i))
 enddo
 !--timestepping
 !  following two are optional: not an error if not present
 call read_inopt(norbits,'norbits',db,err=ierr)
 call read_inopt(deltat,'deltat',db,err=ierr)

 !--mcfost
 if (compiled_with_mcfost) then
    call read_inopt(use_mcfost,'use_mcfost',db,err=ierr)
    if (ierr /= 0) use_mcfost = .false. ! no mcfost by default
    call read_inopt(use_mcfost_stellar_parameters,'use_mcfost_stars',db,err=ierr)
    if (ierr /= 0) use_mcfost_stellar_parameters = .false. ! update stellar parameters by default
 endif

 if (do_radiation) call read_inopt(iradkappa,'radkappa',db,err=ierr)

 if (lumdisc > 0) then
    call read_inopt(L_star(1),'L_star',db,min=0.,errcount=nerr)
    call read_inopt(T_bg,'T_bg',db,min=0.,errcount=nerr)
 endif

 call close_db(db)
 ierr = nerr

 if (nerr > 0) print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'

end subroutine read_setupfile

!--------------------------------------------------------------------------
!
! write options needed for oblate sink particles
!
!--------------------------------------------------------------------------
subroutine write_oblateness_options(iunit,label,J2i,sizei,spin_periodi,kfaci,obliquityi)
 use infile_utils, only:write_inopt
 integer,          intent(in) :: iunit
 character(len=*), intent(in) :: label
 real,             intent(in) :: J2i,sizei,spin_periodi,kfaci,obliquityi

 call write_inopt(J2i,'J2'//trim(label),'J2 moment (oblateness)',iunit)
 if (abs(J2i) > 0.) then
    call write_inopt(sizei,'size'//trim(label),'radius (Jupiter radii)',iunit)
    call write_inopt(spin_periodi,'spin_period'//trim(label),'spin period (hrs)',iunit)
    call write_inopt(kfaci,'kfac'//trim(label),'concentration parameter',iunit)
    call write_inopt(obliquityi,'obliquity'//trim(label),'obliquity (degrees)',iunit)
 endif

end subroutine write_oblateness_options

!--------------------------------------------------------------------------
!
! read options needed for oblate sink particles
!
!--------------------------------------------------------------------------
subroutine read_oblateness_options(db,nerr,label,J2i,sizei,spin_periodi,kfaci,obliquityi)
 use infile_utils, only:inopts,read_inopt
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,          intent(inout) :: nerr
 character(len=*), intent(in)    :: label
 real,             intent(inout) :: J2i,sizei,spin_periodi,kfaci,obliquityi

 call read_inopt(J2i,'J2'//trim(label),db,min=-1.0,max=1.0) ! optional, no error if not read
 if (abs(J2i) > 0.) then
    call read_inopt(sizei,'size'//trim(label),db,errcount=nerr)
    call read_inopt(spin_periodi,'spin_period'//trim(label),db,errcount=nerr)
    call read_inopt(kfaci,'kfac'//trim(label),db,min=0.,max=1.,errcount=nerr)
    call read_inopt(obliquityi,'obliquity'//trim(label),db,min=0.,max=180.,errcount=nerr)
 endif

end subroutine read_oblateness_options

!--------------------------------------------------------------------------
!
! print information about oblateness on sink particles
!
!--------------------------------------------------------------------------
subroutine print_oblateness_info(isink,spin_period_hrs)
 use vectorutils, only:unitvec,mag
 use units,       only:unit_angmom
 use physcon,     only:earthr,jupiterr,au
 integer, intent(in) :: isink
 real,    intent(in) :: spin_period_hrs
 real :: u(3)

 if (abs(xyzmh_ptmass(iJ2,isink)) > 0.) then
    print "(a,g10.3)",      '      J2 moment: ',xyzmh_ptmass(iJ2,isink)
    print "(a,g10.3,a)",    '           size: ',xyzmh_ptmass(iReff,isink)*udist/jupiterr,' Jupiter radii'
    print "(a,g10.3,a)",    '           size: ',xyzmh_ptmass(iReff,isink)*udist/earthr,' Earth radii'
    print "(a,g10.3,a)",    '           size: ',xyzmh_ptmass(iReff,isink)*udist/au,' au'
    u = unitvec(xyzmh_ptmass(ispinx:ispinz,isink))
    print "(a,g10.3,a)",    '      obliquity: ',acos(u(3))/deg_to_rad,' degrees to z=0 plane'
    print "(a,g10.3,a)",    '         period: ',spin_period_hrs,' hrs'
    print "(a,3(g10.3,1x))",'       spin vec: ',u
    print "(/,a,g10.3,a)",    '# spin angular momentum =  ',&
           mag(xyzmh_ptmass(ispinx:ispinz,isink))*unit_angmom,' g cm^2 / s'
    print "(/,a,'(',3(es10.2,1x),')')",' specific spin angular momentum = ',&
                   xyzmh_ptmass(ispinx:ispinz,isink)/xyzmh_ptmass(4,isink)
 endif

end subroutine print_oblateness_info

!--------------------------------------------------------------------------
!
! Set dustfrac
!
!--------------------------------------------------------------------------
subroutine set_dustfrac(disc_index,ipart_start,ipart_end,xyzh,xorigini)

 use grids_for_setup, only:init_grid_sigma,deallocate_sigma
 integer, intent(in) :: disc_index
 integer, intent(in) :: ipart_start
 integer, intent(in) :: ipart_end
 real,    intent(in) :: xyzh(:,:)
 real,    intent(in) :: xorigini(3)

 integer :: i,j
 real    :: R,z
 real    :: dust_to_gasi(maxdusttypes)
 real    :: dust_to_gas_disc
 real    :: Hg,Hd
 real    :: sigma_gas,sigma_gas_sum
 real    :: sigma_dust,sigma_dust_sum
 real, parameter :: tol = 1.e-10

 dust_to_gasi   = 0.
 sigma_gas_sum  = 0.
 sigma_dust_sum = 0.
 if (sigmaprofilegas(disc_index)==6) call init_grid_sigma(R_in(disc_index),R_out(disc_index))
 do i=ipart_start,ipart_end

    R = sqrt(dot_product(xyzh(1:2,i)-xorigini(1:2),xyzh(1:2,i)-xorigini(1:2)))
    z = xyzh(3,i) - xorigini(3)

    Hg = get_H(H_R(disc_index)*R_ref(disc_index),qindex(disc_index),R/R_ref(disc_index))
    sigma_gas = sig_norm(disc_index) * scaled_sigma(R,&
                                                    sigmaprofilegas(disc_index),&
                                                    pindex(disc_index),&
                                                    R_ref(disc_index),&
                                                    R_in(disc_index),&
                                                    R_out(disc_index),&
                                                    R_c(disc_index))
    !--Sum the gas masses
    if ((sigma_gas < huge(sigma_gas)) .and. (sigma_gas == sigma_gas)) then
       sigma_gas_sum = sigma_gas_sum + sigma_gas
    endif

    do j=1,ndustsmall
       if (isetdust > 0 .and. (R<R_indust(disc_index,j) .or. R>R_outdust(disc_index,j))) then
          dust_to_gasi(j) = tiny(dust_to_gasi(j))
          sigma_dust = 0.
       else
          Hd = get_H(H_R_dust(disc_index,j)*R_ref(disc_index),qindex_dust(disc_index,j),R/R_ref(disc_index))
          sigma_dust = sig_normdust(disc_index,j) * scaled_sigma(R,&
                                           sigmaprofiledust(disc_index,j),&
                                           pindex_dust(disc_index,j),&
                                           R_ref(disc_index),&
                                           R_indust(disc_index,j),&
                                           R_outdust(disc_index,j),&
                                           R_c_dust(disc_index,j))
          dust_to_gasi(j) = (sigma_dust/sigma_gas) * (Hg/Hd) * exp(-0.5d0*((z/Hd)**2.-(z/Hg)**2.))
       endif
       !--Sum the dust masses
       if ((sigma_dust < huge(sigma_dust)) .and. (sigma_dust == sigma_dust)) then
          sigma_dust_sum = sigma_dust_sum + sigma_dust
       endif
    enddo
    !--Calculate the final dustfrac that will be output to the dump file
    !  Note: dust density and dust fraction have the same dependence on grain size
    !  ===>  dustfrac(:) = sum(dustfrac)*rhodust(:)/sum(rhodust)
    dustfrac(1:ndustsmall,i) = (sum(dust_to_gasi)/(1.+sum(dust_to_gasi)))*dustbinfrac(1:ndustsmall)
 enddo
 !--Check if the total dust-to-gas ratio is equal to the requested ratio in the setup file
 dust_to_gas_disc = sigma_dust_sum/sigma_gas_sum
 if (abs(dust_to_gas_disc-dust_to_gas)/dust_to_gas > tol) then
    write(*,"(a,es15.8)") ' Requested dust-to-gas ratio is ',dust_to_gas
    write(*,"(a,es15.8)") '    Actual dust-to-gas ratio is ',dust_to_gas_disc
    call fatal('setup_disc','dust-to-gas ratio is not correct')
 endif

 if (sigmaprofilegas(disc_index)==6) call deallocate_sigma()

end subroutine set_dustfrac
!--------------------------------------------------------------------------
!
! Scale height as a function of radius
!
!--------------------------------------------------------------------------
real function get_H(h0,qindex,r)
 real, intent(in) :: h0
 real, intent(in) :: qindex
 real, intent(in) :: r

 get_H = h0*(r**(-qindex+1.5))

end function get_H
!--------------------------------------------------------------------------
!
! Spherical density profile as a function of radius
!
!--------------------------------------------------------------------------
real function atm_dens(r)
 use eos, only:gamma
 real, intent(in) :: r

 select case(atm_type)
 case(1)
    atm_dens = r**(-3)
 case(2)
    atm_dens = r**(-1./(gamma - 1.))
 case default
    !atm_dens = exp(-(r-r_planet)/scaleheight)
    stop 'atmosphere not yet implemented...stopping!'
 end select

end function atm_dens

!--------------------------------------------------------------------------
!
! Return the sound speed given the radius
!
!--------------------------------------------------------------------------
pure real function cs_func(cs0,r,q_index)
 real, intent(in) :: cs0
 real, intent(in) :: r
 real, intent(in) :: q_index

 cs_func = cs0*r**(-q_index)

end function cs_func

!--------------------------------------------------------------------------
!
! Convert to a corotating frame around one of the binary masses
!
!--------------------------------------------------------------------------
subroutine make_corotate(xyzh,vxyzu,a0,Mstar,npart,npart_disc)
 use extern_corotate, only:omega_corotate
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real,    intent(in)    :: a0
 real,    intent(in)    :: Mstar
 integer, intent(in)    :: npart
 integer, intent(in)    :: npart_disc

 integer :: i
 real    :: phipart,r
 real    :: v_0(3),vmag,omega0

 !
 !--Change to corotating frame
 !
 ! Calculate velocity at planet
 vmag   = sqrt(Mstar/a0)
 omega0 = sqrt(Mstar/a0**3)

 ! v_phi = v_y at y=0
 ! Obtain the true v_phi at any point (r,phi) via rotation in z axis

 v_0 = (/0.0,vmag,0.0/)

 print *, 'Transforming to corotating frame: angular velocity ', omega0

 do i=1,npart_disc
    r          = sqrt(xyzh(1,i)**2 + xyzh(2,i)**2)
    phipart    = atan2(xyzh(2,i),xyzh(1,i))
    vxyzu(1,i) = vxyzu(1,i) - r*(-omega0)*sin(phipart)
    vxyzu(2,i) = vxyzu(2,i) + r*(-omega0)*cos(phipart)
 enddo
 vxyzu(1:3,npart_disc+1:npart) = 0.

 omega_corotate = omega0

end subroutine make_corotate

subroutine temp_to_HR(temp,H_R,radius,M,cs)
 use units,  only:get_kbmh_code
 use eos,    only:gmw
 real,    intent(in)    :: temp,radius,M
 real,    intent(out)   :: H_R,cs
 real                   :: omega

 cs = sqrt(temp*get_kbmh_code()/gmw)
 omega = sqrt(M/radius**3)
 H_R = cs/(omega*radius)

end subroutine temp_to_HR

subroutine get_hier_disc_label(i, disclabel)
 use sethierarchical, only:hs!sink_num, sink_labels, hl_labels
 character(len=10), intent(out)  :: disclabel
 integer, intent(in) :: i

 if (i <= hs%labels%sink_num) then
    disclabel = trim(hs%labels%sink(i))
 else
    disclabel = trim(hs%labels%hl(i-hs%labels%sink_num))
 endif

end subroutine get_hier_disc_label

end module setup
