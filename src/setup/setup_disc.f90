!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
!   - Ratm_in       : *inner atmosphere radius (planet radii)*
!   - Ratm_out      : *outer atmosphere radius (planet radii)*
!   - accr1         : *single star accretion radius*
!   - accr1a        : *single star accretion radius*
!   - accr1b        : *single star accretion radius*
!   - accr2         : *perturber accretion radius*
!   - accr2a        : *tight binary primary accretion radius*
!   - accr2b        : *tight binary secondary accretion radius*
!   - alphaSS       : *desired alphaSS*
!   - alpha_z       : *height of transition in tanh vertical temperature profile*
!   - atm_type      : *atmosphere type (1:r**(-3); 2:r**(-1./(gamma-1.)))*
!   - beta_z        : *variation in transition height over radius*
!   - bhspin        : *black hole spin*
!   - bhspinangle   : *black hole spin angle (deg)*
!   - binary1_O     : *tight binary 1 Omega, PA of ascending node (deg)*
!   - binary1_a     : *tight binary 1 semi-major axis*
!   - binary1_e     : *tight binary 1 eccentricity*
!   - binary1_f     : *tight binary 1 f, initial true anomaly (deg,180=apastron)*
!   - binary1_i     : *tight binary 1 i, inclination (deg)*
!   - binary1_w     : *tight binary 1 w, argument of periapsis (deg)*
!   - binary2_O     : *tight binary 2 Omega, PA of ascending node (deg)*
!   - binary2_a     : *tight binary 2 semi-major axis*
!   - binary2_e     : *tight binary 2 eccentricity*
!   - binary2_f     : *tight binary 2 f, initial true anomaly (deg,180=apastron)*
!   - binary2_i     : *tight binary 2 i, inclination (deg)*
!   - binary2_w     : *tight binary 2 w, argument of periapsis (deg)*
!   - binary_O      : *wide binary Omega, PA of ascending node (deg)*
!   - binary_a      : *wide binary semi-major axis*
!   - binary_e      : *wide binary eccentricity*
!   - binary_f      : *wide binary f, initial true anomaly (deg,180=apastron)*
!   - binary_i      : *wide binary i, inclination (deg)*
!   - binary_w      : *wide binary w, argument of periapsis (deg)*
!   - deltat        : *output interval as fraction of orbital period*
!   - discstrat     : *stratify disc? (0=no,1=yes)*
!   - dist_unit     : *distance unit (e.g. au,pc,kpc,0.1pc)*
!   - einst_prec    : *include Einstein precession*
!   - flyby_O       : *position angle of ascending node (deg)*
!   - flyby_a       : *distance of minimum approach*
!   - flyby_d       : *initial distance (units of dist. min. approach)*
!   - flyby_i       : *inclination (deg)*
!   - ibinary       : *binary orbit (0=bound,1=unbound [flyby])*
!   - ipotential    : *potential (1=central point mass,*
!   - istrat        : *temperature prescription (0=MAPS, 1=Dartois)*
!   - lumdisc       : *Set qindex from stellar luminosity (ieos=24) (0=no 1=yes)*
!   - m1            : *first hierarchical level primary mass*
!   - m2            : *first hierarchical level secondary mass*
!   - mass_unit     : *mass unit (e.g. solarm,jupiterm,earthm)*
!   - norbits       : *maximum number of orbits at outer disc*
!   - np            : *number of gas particles*
!   - nplanets      : *number of planets*
!   - nsinks        : *number of sinks*
!   - q1            : *tight binary 1 mass ratio*
!   - q2            : *tight binary 2 mass ratio*
!   - qatm          : *sound speed power law index of atmosphere*
!   - radkappa      : *constant radiation opacity kappa*
!   - ramp          : *Do you want to ramp up the planet mass slowly?*
!   - rho_core      : *planet core density (cgs units)*
!   - subst         : *star to substitute*
!   - subst1        : *first star to substitute*
!   - subst2        : *second star to substitute*
!   - surface_force : *model m1 as planet with surface*
!   - temp_atm0     : *atmosphere temperature scaling factor*
!   - temp_mid0     : *midplane temperature scaling factor*
!   - use_mcfost    : *use the mcfost library*
!   - z0            : *z scaling factor*
!
! :Dependencies: centreofmass, dim, dust, eos, eos_stamatellos,
!   extern_binary, extern_corotate, extern_lensethirring, externalforces,
!   fileutils, growth, infile_utils, io, kernel, memory, options, part,
!   physcon, porosity, prompting, radiation_utils, set_dust,
!   set_dust_options, setbinary, setdisc, setflyby, sethierarchical,
!   spherical, timestep, units, vectorutils
!
 use dim,              only:use_dust,maxalpha,use_dustgrowth,maxdusttypes,&
                            maxdustlarge,maxdustsmall,compiled_with_mcfost
 use externalforces,   only:iext_star,iext_binary,iext_lensethirring,&
                            iext_einsteinprec,iext_corot_binary,iext_corotate,&
                            update_externalforce
 use extern_binary,    only:mass2,accradius1,accradius2,ramp,surface_force,eps_soft1
 use fileutils,        only:make_tags_unique
 use growth,           only:ifrag,isnow,rsnow,Tsnow,vfragSI,vfraginSI,vfragoutSI,gsizemincgs
 use porosity,         only:iporosity
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
 use setdisc,          only:scaled_sigma,get_disc_mass,maxbins
 use set_dust_options, only:set_dust_default_options,dust_method,dust_to_gas,&
                            ndusttypesinp,ndustlargeinp,ndustsmallinp,isetdust,&
                            dustbinfrac,check_dust_method
 use units,            only:umass,udist,utime
 use dim,              only:do_radiation
 use radiation_utils,  only:set_radiation_and_gas_temperature_equal
 use memory,           only:allocate_memory

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
 real    :: binary_a,binary_e,binary_i,binary_O,binary_w,binary_f
 real    :: binary1_a,binary1_e,binary1_i,binary1_O,binary1_w,binary1_f
 real    :: binary2_a,binary2_e,binary2_i,binary2_O,binary2_w,binary2_f
 integer :: icentral,ipotential,ibinary
 integer :: nsinks,subst,subst1,subst2
 real    :: flyby_a,flyby_d,flyby_O,flyby_i
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
 logical :: ismoothgas(maxdiscs)
 logical :: itapergas(maxdiscs)
 integer :: itapersetgas(maxdiscs)
 logical :: iwarp(maxdiscs)
 logical :: use_global_iso
 real    :: alphaSS

 real    :: R_in(maxdiscs),R_out(maxdiscs),R_ref(maxdiscs),R_c(maxdiscs)
 real    :: pindex(maxdiscs),disc_m(maxdiscs),sig_ref(maxdiscs),sig_norm(maxdiscs)
 real    :: T_bg,L_star(maxdiscs)
 real    :: qindex(maxdiscs),H_R(maxdiscs)
 real    :: posangl(maxdiscs),incl(maxdiscs)
 real    :: annulus_m(maxdiscs),R_inann(maxdiscs),R_outann(maxdiscs)
 real    :: R_warp(maxdiscs),H_warp(maxdiscs)
 real    :: Q_min(maxdiscs)

 integer :: sigmaprofiledust(maxdiscs,maxdusttypes)
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

 !--planets
 integer, parameter :: maxplanets = 9

 character(len=*), dimension(maxplanets), parameter :: num = &
    (/'1','2','3','4','5','6','7','8','9' /)

 logical :: istratify
 integer :: nplanets,discstrat,lumdisc
 real    :: mplanet(maxplanets),rplanet(maxplanets)
 real    :: accrplanet(maxplanets),inclplan(maxplanets)
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

 !--units
 character(len=20) :: dist_unit,mass_unit

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

 !--setup units
 call setup_units()

 !--compute number of discs based on setup options
 call number_of_discs()

 !--setup central object(s), i.e. sink particle(s) or potential
 call setup_central_objects(fileprefix)

 !--setup equation of state
 call equation_of_state(gamma)

 !--set surface density profile based on setup options
 call surface_density_profile()

 !--setup grain size distribution
 call setup_dust_grain_distribution()

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

 !--reset centre of mass to the origin
 call set_centreofmass(npart,xyzh,vxyzu)

 !--set tmax and dtmax
 call set_tmax_dtmax()

 if (do_radiation) then
    rad(iradxi,1:npart)=0.!call set_radiation_and_gas_temperature_equal(npart,xyzh,vxyzu,massoftype,rad)
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
subroutine set_default_options()!id)
 use sethierarchical, only:set_hierarchical_default_options
!  integer, intent(in) :: id

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

 !--spinning black hole (Lense-Thirring)
 einst_prec = .false.
 bhspin      = 1.
 bhspinangle = 0.

 !--sink particle(s)
 nsinks = 1
 ibinary = 0

 !--binary
 binary_a = 10.
 binary_e = 0.
 binary_i = 0.
 binary_O = 0.
 binary_w = 270.
 binary_f = 180.

 !--hierarchical
 call set_hierarchical_default_options()!id)

 !--flyby
 flyby_a  = 200.
 flyby_d  = 10.
 flyby_O  = 0.
 flyby_i  = 0.

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

 !--dust disc
 R_indust      = 1.
 R_outdust     = 150.
 pindex_dust   = 1.
 qindex_dust   = 0.25
 H_R_dust      = 0.05
 itaperdust    = .false.
 itapersetdust = 0
 ismoothdust   = .true.
 R_c_dust      = 150.

 !--dust growth
 ifrag = 1
 isnow = 0
 rsnow = 100.
 Tsnow = 150.
 vfragSI = 15.
 vfraginSI = 5.
 vfragoutSI = 15.
 gsizemincgs = 5.e-3

 !--resolution
 np = 1000000
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
          read(23,*,iostat=ierr) binary_a,binary_e,binary_i,binary_O,binary_w,binary_f
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
! Set the units
!
!--------------------------------------------------------------------------
subroutine setup_units()
 use units, only:set_units,select_unit

 integer :: ierr

 if (icentral==0 .and. ipotential==3) then
    !--black hole units
    !  note: distance unit not used but (currently) required in the .setup file
    call select_unit(mass_unit,umass,ierr)
    if (ierr /= 0) call error('setup_disc','mass unit not recognised')
    call set_units(mass=umass,c=1.d0)
 else
    !--stellar units
    call select_unit(mass_unit,umass,ierr)
    if (ierr /= 0) call error('setup_disc','mass unit not recognised')
    call select_unit(dist_unit,udist,ierr)
    if (ierr /= 0) call error('setup_disc','length unit not recognised')
    call set_units(dist=udist,mass=umass,G=1.d0)
 endif

end subroutine setup_units

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
 use eos,     only:isink,qfacdisc,qfacdisc2,polyk2,beta_z,z0
 use options, only:ieos,icooling
 use options, only:nfulldump,alphau,ipdv_heating,ishock_heating
 use eos_stamatellos, only:init_coolra
 use physcon, only:rpiontwo
 real, intent(out) :: gamma
 real              :: H_R_atm, cs

 logical :: is_isothermal
 integer :: i

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
             higher_disc_index = findloc(iuse_disc, .true., 1)
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
       if (qindex(onlydisc) > 0.) then
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
       endif
    endif

 else
    !-- adiabatic
    if (lumdisc > 0) then
       !--for radapprox cooling
       print "(/,a)", ' setting ieos=24 and icooling=9 for radiative cooling approximation'
       ieos = 24
       icooling = 9
       gamma = 5./3. ! in case it's needed
       call init_coolra()
       if (ndiscs > 1) then
          print *, "We can't set up multiple radapprox discs yet :,("
          stop
       else
          cs = get_cs_from_lum(L_star(1),R_ref(1)) / rpiontwo
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
 use sethierarchical,      only:set_hierarchical,set_multiple
 use setflyby,             only:set_flyby
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
       print "(a,g10.3,a)",'   Accretion Radius:       ', accr1, trim(dist_unit)
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
       select case (ibinary)
       case (0)
          !--bound
          nptmass  = 0
          print "(/,a)",' Central objects represented by two sinks'
          print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
          print "(a,g10.3)",  '   Binary mass ratio:  ', m2/m1
          print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
          print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)
          call set_binary(m1,m2,semimajoraxis=binary_a,eccentricity=binary_e, &
                          posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
                          f=binary_f,accretion_radius1=accr1,accretion_radius2=accr2, &
                          xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)
          mcentral = m1 + m2
       case (1)
          !--unbound (flyby)
          print "(/,a)",' Central object represented by a sink at the system origin with a perturber sink'
          print "(a,g10.3,a)",'   Primary mass:       ', m1,    trim(mass_unit)
          print "(a,g10.3,a)",'   Perturber mass:     ', m2,    trim(mass_unit)
          print "(a,g10.3,a)",'   Accretion Radius 1: ', accr1, trim(dist_unit)
          print "(a,g10.3,a)",'   Accretion Radius 2: ', accr2, trim(dist_unit)
          call set_flyby(m1,m2,minimum_approach=flyby_a, &
                         initial_dist=flyby_d,posang_ascnode=flyby_O,inclination=flyby_i, &
                         accretion_radius1=accr1,accretion_radius2=accr2, &
                         xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)
          mcentral = m1
       end select
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



       call set_multiple(m1,m2,semimajoraxis=binary_a,eccentricity=binary_e, &
            posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
            f=binary_f,accretion_radius1=accr1,accretion_radius2=accr1, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)

       if (iuse_disc(1)) then
          discpos=xyzmh_ptmass(1:3,2)
          discvel=vxyz_ptmass(1:3,2)
       else
          discpos = 0.
          discvel = 0.
       endif

       call set_multiple(m2/(q2+1),m2*q2/(q2+1),semimajoraxis=binary2_a,eccentricity=binary2_e, &
            posang_ascnode=binary2_O,arg_peri=binary2_w,incl=binary2_i, &
            f=binary2_f,accretion_radius1=accr2a,accretion_radius2=accr2b, &
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

       call set_multiple(m1,m2,semimajoraxis=binary_a,eccentricity=binary_e, &
            posang_ascnode=binary_O,arg_peri=binary_w,incl=binary_i, &
            f=binary_f,accretion_radius1=accr1,accretion_radius2=accr1, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)

       discpos=xyzmh_ptmass(1:3,2)
       discvel=vxyz_ptmass(1:3,2)

       call set_multiple(m1/(q1+1),m1*q1/(q1+1),semimajoraxis=binary1_a,eccentricity=binary1_e, &
         posang_ascnode=binary1_O,arg_peri=binary1_w,incl=binary1_i, &
         f=binary1_f,accretion_radius1=accr1a,accretion_radius2=accr1b, &
         xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass, subst=subst1,ierr=ierr)

       call set_multiple(m2/(q2+1),m2*q2/(q2+1),semimajoraxis=binary2_a,eccentricity=binary2_e, &
            posang_ascnode=binary2_O,arg_peri=binary2_w,incl=binary2_i, &
            f=binary2_f,accretion_radius1=accr2a,accretion_radius2=accr2b, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass, subst=subst2,ierr=ierr)

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
! Set the grain size distribution
!
!--------------------------------------------------------------------------
subroutine setup_dust_grain_distribution()
 use dust,             only:grainsizecgs,graindenscgs
 use set_dust,         only:set_dustbinfrac
 use set_dust_options, only:grainsizeinp,graindensinp,igrainsize,igraindens,&
                            smincgs,smaxcgs,sindex

 if (use_dust) then
    grainsize = 0.
    graindens = 0.
    if (ndusttypes > 1) then
       select case(igrainsize)
       case(0)
          call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
          grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
       case(1)
          grainsize(1:ndusttypes) = grainsizeinp(1:ndusttypes)/udist
       end select
       select case(igraindens)
       case(0)
          graindens(1:ndusttypes) = graindensinp(1)/umass*udist**3
       case(1)
          graindens(1:ndusttypes) = graindensinp(1:ndusttypes)/umass*udist**3
       end select
    else
       grainsize(1) = grainsizeinp(1)/udist
       graindens(1) = graindensinp(1)/umass*udist**3
       grainsizecgs = grainsizeinp(1)
       graindenscgs = graindensinp(1)
    endif
 endif

end subroutine setup_dust_grain_distribution

!--------------------------------------------------------------------------
!
! Calculate the required disc masses
!
!--------------------------------------------------------------------------
subroutine calculate_disc_mass()

 integer :: i,j
 real :: enc_m(maxbins),rad(maxbins)
 real :: Q_mintmp,disc_mtmp,annulus_mtmp
 real :: rgrid_min,rgrid_max,fac

 totmass_gas  = 0.
 disc_mdust = 0.

 do i=1,maxdiscs
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
       print*,'dust'
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
 enddo

end subroutine calculate_disc_mass

!--------------------------------------------------------------------------
!
! Set up the discs
!
!--------------------------------------------------------------------------
subroutine setup_discs(id,fileprefix,hfact,gamma,npart,polyk,&
                       npartoftype,massoftype,xyzh,vxyzu)
 use options,   only:alpha
 use setbinary, only:Rochelobe_estimate
 use sethierarchical, only:get_hierarchical_level_com, get_hier_level_mass
 !use sethierarchical, only:hl_labels, a
 use sethierarchical, only:hs
 use setdisc,   only:set_disc
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
 if (maxalpha==0) alpha = alphaSS
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
             Rochelobe = Rochelobe_estimate(m1,m2,binary_a)
          case(2)
             !--circumprimary
             xorigini  = xyzmh_ptmass(1:3,1)
             vorigini  = vxyz_ptmass(1:3,1)
             Rochelobe = Rochelobe_estimate(m2,m1,binary_a)
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

             hl_index = findloc(hs%labels%hl, disclabel(:len(trim(disclabel))-1), 1)
             Rochelobe = Rochelobe_estimate(m1,m2,hs%levels(hl_index)%a)
          else
             Rochelobe = huge(0.)
          endif


          star_m(i) = m2

          call get_hierarchical_level_com(disclabel, xorigini, vorigini, xyzmh_ptmass, vxyz_ptmass, fileprefix)

       endif

       if ((ndiscs > 1 .and. ibinary==0) .and. (R_out(i) > Rochelobe)) then
          call warning('setup_disc', &
             'Outer disc radius for circum-'//trim(disctype(i))//' > Roche lobe of ' &
             //trim(disctype(i)))
       endif

       !--taper gas disc
       iprofilegas = 0
       if (itapergas(i)) iprofilegas = 1
       if (itapersetgas(i) == 1) iprofilegas = 2

       !--set disc(s)
       if (use_dust .and. use_dustfrac) then

          !--taper dust disc
          iprofiledust = 0
          if (itaperdust(i,1)) iprofiledust = 1
          if (itapersetdust(i,1) == 1) iprofiledust = 2

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
                        bh_spin          = bhspin,               &
                        enc_mass         = enc_mass(:,i),        &
                        prefix           = prefix)

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
                              bh_spin        = bhspin,             &
                              enc_mass       = enc_mass(:,i),      &
                              prefix         = dustprefix(j))

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
                        bh_spin         = bhspin,             &
                        enc_mass        = enc_mass(:,i),      &
                        prefix          = prefix)

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
                              bh_spin        = bhspin,             &
                              enc_mass       = enc_mass(:,i),      &
                              prefix         = dustprefix(j))

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
 use setbinary, only:get_mean_angmom_vector
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)
 real,    intent(in) :: vxyzu(:,:)

!!!real :: ldisc(maxdiscs),lcentral(maxdiscs)
 real :: ldisc(3),lcentral(3)

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
 integer, intent(in) :: npart
 real,    intent(in) :: massoftype(:)
 real,    intent(in) :: xyzh(:,:)

 integer :: i,j,itype
 real    :: dist_bt_sinks
 real    :: phi,vphi,sinphi,cosphi,omega,r2,disc_m_within_r
 real    :: Hill(maxplanets)
 real    :: u(3)

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
          dist_bt_sinks = sqrt(dot_product(xyzmh_ptmass(1:3,1),xyzmh_ptmass(1:3,2)))
          if (rplanet(i) > dist_bt_sinks) Hill(i) = (mplanet(i)*jupiterm/umass/(3.*m1))**(1./3.) * rplanet(i)
       endif
       xyzmh_ptmass(1:3,nptmass)    = (/rplanet(i)*cosphi,rplanet(i)*sinphi,0./)
       xyzmh_ptmass(4,nptmass)      = mplanet(i)*jupiterm/umass
       xyzmh_ptmass(ihacc,nptmass)  = accrplanet(i)*Hill(i)
       xyzmh_ptmass(ihsoft,nptmass) = 0.
       vphi                         = sqrt((mcentral + disc_m_within_r)/rplanet(i))
       if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) vphi = sqrt((m1 + disc_m_within_r)/rplanet(i))
       vxyz_ptmass(1:3,nptmass)     = (/-vphi*sinphi,vphi*cosphi,0./)
       if (nsinks == 2 .and. rplanet(i) < dist_bt_sinks) then
          vxyz_ptmass(1:3,nptmass) = vxyz_ptmass(1:3,nptmass) + vxyz_ptmass(1:3,1)
          xyzmh_ptmass(1:3,nptmass) = xyzmh_ptmass(1:3,nptmass) + xyzmh_ptmass(1:3,1)
       endif

       !--incline positions and velocities
       inclplan(i) = inclplan(i)*deg_to_rad
       u = (/-sin(phi),cos(phi),0./)
       call rotatevec(xyzmh_ptmass(1:3,nptmass),u,-inclplan(i))
       call rotatevec(vxyz_ptmass(1:3,nptmass), u,-inclplan(i))

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
subroutine set_tmax_dtmax()
 use setflyby, only:get_T_flyby
 use timestep, only:tmax,dtmax

 real :: period, period1, period2

 period2 = 0.
 if (icentral==1 .and. nsinks==2 .and. ibinary==0) then
    !--binary orbital period
    period = sqrt(4.*pi**2*binary_a**3/mcentral)
 elseif (icentral==1 .and. nsinks==3 .and. ibinary==0) then
    !--wide binary orbital period
    period = sqrt(4.*pi**2*binary_a**3/mcentral)
    !--tight binary orbital period
    period2 = sqrt(4.*pi**2*binary2_a**3/m2)
 elseif (icentral==1 .and. nsinks==4 .and. ibinary==0) then
    !--wide binary orbital period
    period = sqrt(4.*pi**2*binary_a**3/mcentral)
    !--tight binary 1 orbital period
    period1 = sqrt(4.*pi**2*binary1_a**3/m1)
    !--tight binary 2 orbital period
    period2 = sqrt(4.*pi**2*binary2_a**3/m2)
 elseif (icentral==1 .and. nsinks==2 .and. ibinary==1) then
    !--time of flyby
    period = get_T_flyby(m1,m2,flyby_a,flyby_d)
 elseif (nplanets > 0) then
    !--outer planet orbital period
    period = period_planet_longest
 elseif (iwarp(onlydisc)) then
    !--warp period
    period = sqrt(4.*pi**2*R_warp(onlydisc)**3/mcentral)
 else
    !--outer disc orbital period
    period = sqrt(4.*pi**2*R_out(onlydisc)**3/mcentral)
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
 use set_dust_options, only:set_dust_interactively
 use sethierarchical, only:set_hierarchical_default_options, get_hier_level_mass
 use sethierarchical, only:hs, hierarchy, print_chess_logo, generate_hierarchy_string!sink_num, hl_num, sink_labels, hl_labels

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

          call prompt('Enter orbital radius of the planet (e.g. 5.2au)',dist_unit)
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
       call prompt('Do you want the binary orbit to be bound (elliptic) or'// &
                  ' unbound (parabolic/hyperbolic) [flyby]?'//new_line('A')// &
                  ' 0=bound'//new_line('A')//' 1=unbound'//new_line('A'),ibinary,0,1)
       select case (ibinary)
       case (0)
          !--bound
          m1       = 1.
          m2       = 0.2
          binary_a = 10.
          binary_e = 0.
          binary_i = 0.
          binary_O = 0.
          binary_w = 270.
          binary_f = 180.
          accr1    = 1.
          accr2    = 0.5
       case (1)
          !--unbound (flyby)
          m1       = 1.
          m2       = 1.
          accr1    = 1.
          accr2    = 1.
          flyby_a  = 200.
          flyby_d  = 10.
          flyby_O  = 0.
          flyby_i  = 0.
       end select

    case (5:)

       call print_chess_logo()!id)

       ibinary = 0

       call generate_hierarchy_string(nsinks)

       call prompt('What is the hierarchy?',hierarchy)
       !call set_hierarchical_interactively()
       call set_hierarchical_default_options()

    case (3)
       !-- hierarchical triple --!
       print "(/,a)",'================================'
       print "(a)",  '+++   HIERARCHICAL TRIPLE    +++'
       print "(a)",  '================================'
       ibinary = 0

       !-- Wide binary
       m1       = 1.
       m2       = 0.2
       binary_a = 10.
       binary_e = 0.
       binary_i = 0.
       binary_O = 0.
       binary_w = 270.
       binary_f = 180.
       accr1    = 1.

       !-- Tight binary
       subst    = 12
       q2       = 1
       m2a      = m2/(q2+1)
       m2b      = m2*q2/(q2+1)
       binary2_a = 1.
       binary2_e = 0.
       binary2_i = 0.
       binary2_O = 0.
       binary2_w = 270.
       binary2_f = 180.
       accr2a    = 0.1
       accr2b    = 0.1
    case(4)
       !-- hierarchical quadruple --!
       print "(/,a)",'================================'
       print "(a)",  '+++   HIERARCHICAL QUADRUPLE    +++'
       print "(a)",  '================================'
       ibinary = 0

       !-- Wide binary
       m1       = 1.
       m2       = 0.2
       binary_a = 10.
       binary_e = 0.
       binary_i = 0.
       binary_O = 0.
       binary_w = 270.
       binary_f = 180.
       accr1    = 1.

       !-- Tight binary 1
       subst1    = 11
       q1        = 1
       m1a       = m1/(q1+1)
       m1b       = m1*q1/(q1+1)
       binary1_a = 1.
       binary1_e = 0.
       binary1_i = 0.
       binary1_O = 0.
       binary1_w = 270.
       binary1_f = 180.
       accr1a    = 0.1
       accr1b    = 0.1

       !-- Tight binary 2
       subst2    = 12
       q2       = 1
       m2a      = m2/(q2+1)
       m2b      = m2*q2/(q2+1)
       binary2_a = 1.
       binary2_e = 0.
       binary2_i = 0.
       binary2_O = 0.
       binary2_w = 270.
       binary2_f = 180.
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
    if (ibinary==0 .and. nsinks==2) then
       !--bound binary: circum-binary, -primary, -secondary
       iuse_disc(1) = .true.
       iuse_disc(2) = .false.
       iuse_disc(3) = .false.
       call prompt('Do you want a circumbinary disc?',iuse_disc(1))
       call prompt('Do you want a circumprimary disc?',iuse_disc(2))
       call prompt('Do you want a circumsecondary disc?',iuse_disc(3))
    elseif (ibinary==1) then
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
 if (maxalpha==0) alphaSS = 0.005
 if (surface_force) then
    R_in  = 0.1
    R_out = 3.
    R_ref = 1.
 endif
 if ((icentral==1 .and. nsinks>=2) .and. (ibinary==0)) then
    !--don't smooth circumbinary, by default
    ismoothgas(1) = .false.
    !--set appropriate disc radii for bound binary
    R_in(1:4)      = (/2.5*binary_a, accr1, accr2, 2.5*binary_a /)
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
       ! to be changed also in the the setpart function.
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
             higher_disc_index = findloc(iuse_disc, .true., 1)
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
    endif
 enddo

 !--dust disc
 if (use_dust) then
    print "(/,a)",'=============='
    print "(a)",  '+++  DUST  +++'
    print "(a)",  '=============='
    !--dust distribution
    call set_dust_interactively()
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
 elseif (icentral==1 .and. nsinks==2 .and. ibinary==0) then
    call prompt('Enter time between dumps as fraction of binary period',deltat,0.)
    call prompt('Enter number of orbits to simulate',norbits,0)
 elseif (icentral==1 .and. nsinks>=3 .and. ibinary==0) then
    call prompt('Enter time between dumps as fraction of binary period'//new_line('A')// &
         '(enter a negative number to refer to the shorter period)',deltat)
    call prompt('Enter number of orbits to simulate'//new_line('A')// &
         '(enter a negative number to refer to the shorter period)',norbits)
 elseif (icentral==1 .and. nsinks==2 .and. ibinary==1) then
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
 use sethierarchical, only:write_hierarchical_setupfile
 use sethierarchical, only:hs!sink_num, hl_num, sink_labels, hl_labels
 character(len=*), intent(in) :: filename

 integer, parameter :: iunit = 20
 logical            :: done_alpha
 integer            :: i,j,n_possible_discs
 character(len=20)  :: duststring(maxdusttypes)
 character(len=20)  :: taper_string
 character(len=20)  :: smooth_string
 character(len=40)  :: tmpstr

 done_alpha = .false.
 n_possible_discs = 1
 if ((icentral==1) .and. (nsinks>=2)) n_possible_discs = 3

 print "(/,a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
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
 write(iunit,"(/,a)") '# units'
 call write_inopt(dist_unit,'dist_unit','distance unit (e.g. au,pc,kpc,0.1pc)',iunit)
 call write_inopt(mass_unit,'mass_unit','mass unit (e.g. solarm,jupiterm,earthm)',iunit)
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
       call write_inopt(ibinary,'ibinary','binary orbit (0=bound,1=unbound [flyby])',iunit)
       select case (ibinary)
       case (0)
          !--bound
          write(iunit,"(/,a)") '# options for binary'
          call write_inopt(m1,'m1','primary mass',iunit)
          call write_inopt(m2,'m2','secondary mass',iunit)
          call write_inopt(binary_a,'binary_a','binary semi-major axis',iunit)
          call write_inopt(binary_e,'binary_e','binary eccentricity',iunit)
          call write_inopt(binary_i,'binary_i','i, inclination (deg)',iunit)
          call write_inopt(binary_O,'binary_O','Omega, PA of ascending node (deg)',iunit)
          call write_inopt(binary_w,'binary_w','w, argument of periapsis (deg)',iunit)
          call write_inopt(binary_f,'binary_f','f, initial true anomaly (deg,180=apastron)',iunit)
          call write_inopt(accr1,'accr1','primary accretion radius',iunit)
          call write_inopt(accr2,'accr2','secondary accretion radius',iunit)
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
          call write_inopt(flyby_a,'flyby_a','distance of minimum approach',iunit)
          call write_inopt(flyby_d,'flyby_d','initial distance (units of dist. min. approach)',iunit)
          call write_inopt(flyby_O,'flyby_O','position angle of ascending node (deg)',iunit)
          call write_inopt(flyby_i,'flyby_i','inclination (deg)',iunit)
       end select

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
       call write_inopt(binary_a,'binary_a','wide binary semi-major axis',iunit)
       call write_inopt(binary_e,'binary_e','wide binary eccentricity',iunit)
       call write_inopt(binary_i,'binary_i','wide binary i, inclination (deg)',iunit)
       call write_inopt(binary_O,'binary_O','wide binary Omega, PA of ascending node (deg)',iunit)
       call write_inopt(binary_w,'binary_w','wide binary w, argument of periapsis (deg)',iunit)
       call write_inopt(binary_f,'binary_f','wide binary f, initial true anomaly (deg,180=apastron)',iunit)

       !-- tight parameters
       call write_inopt(binary2_a,'binary2_a','tight binary semi-major axis',iunit)
       call write_inopt(binary2_e,'binary2_e','tight binary eccentricity',iunit)
       call write_inopt(binary2_i,'binary2_i','tight binary i, inclination (deg)',iunit)
       call write_inopt(binary2_O,'binary2_O','tight binary Omega, PA of ascending node (deg)',iunit)
       call write_inopt(binary2_w,'binary2_w','tight binary w, argument of periapsis (deg)',iunit)
       call write_inopt(binary2_f,'binary2_f','tight binary f, initial true anomaly (deg,180=apastron)',iunit)

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
       call write_inopt(binary_a,'binary_a','wide binary semi-major axis',iunit)
       call write_inopt(binary_e,'binary_e','wide binary eccentricity',iunit)
       call write_inopt(binary_i,'binary_i','wide binary i, inclination (deg)',iunit)
       call write_inopt(binary_O,'binary_O','wide binary Omega, PA of ascending node (deg)',iunit)
       call write_inopt(binary_w,'binary_w','wide binary w, argument of periapsis (deg)',iunit)
       call write_inopt(binary_f,'binary_f','wide binary f, initial true anomaly (deg,180=apastron)',iunit)

       !-- tight binary 1 parameters
       call write_inopt(binary1_a,'binary1_a','tight binary 1 semi-major axis',iunit)
       call write_inopt(binary1_e,'binary1_e','tight binary 1 eccentricity',iunit)
       call write_inopt(binary1_i,'binary1_i','tight binary 1 i, inclination (deg)',iunit)
       call write_inopt(binary1_O,'binary1_O','tight binary 1 Omega, PA of ascending node (deg)',iunit)
       call write_inopt(binary1_w,'binary1_w','tight binary 1 w, argument of periapsis (deg)',iunit)
       call write_inopt(binary1_f,'binary1_f','tight binary 1 f, initial true anomaly (deg,180=apastron)',iunit)

       !-- tight binary 2 parameters
       call write_inopt(binary2_a,'binary2_a','tight binary 2 semi-major axis',iunit)
       call write_inopt(binary2_e,'binary2_e','tight binary 2 eccentricity',iunit)
       call write_inopt(binary2_i,'binary2_i','tight binary 2 i, inclination (deg)',iunit)
       call write_inopt(binary2_O,'binary2_O','tight binary 2 Omega, PA of ascending node (deg)',iunit)
       call write_inopt(binary2_w,'binary2_w','tight binary 2 w, argument of periapsis (deg)',iunit)
       call write_inopt(binary2_f,'binary2_f','tight binary 2 f, initial true anomaly (deg,180=apastron)',iunit)

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
       call write_inopt(itapergas(i),'itapergas'//trim(disclabel), &
          'exponentially taper the outer disc profile',iunit)
       if (itapergas(i)) call write_inopt(itapersetgas(i),'itapersetgas'//trim(disclabel), &
          'how to set taper (0=exp[-(R/R_c)^(2-p)], 1=[1-exp(R-R_out)]',iunit)
       call write_inopt(ismoothgas(i),'ismoothgas'//trim(disclabel),'smooth inner disc',iunit)
       call write_inopt(iwarp(i),'iwarp'//trim(disclabel),'warp disc',iunit)
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
       if (.not.done_alpha) then
          if (maxalpha==0) call write_inopt(alphaSS,'alphaSS','desired alphaSS',iunit)
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
 !--dust & growth options
 if (use_dust) then
    call write_dust_setup_options(iunit)
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
 elseif (icentral==1 .and. nsinks>=2 .and. ibinary==0) then
    call write_inopt(norbits,'norbits','maximum number of binary orbits',iunit)
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
 use sethierarchical, only:read_hierarchical_setupfile
 use sethierarchical, only:hs!sink_num, hl_num, sink_labels, hl_labels
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
 call read_inopt(mass_unit,'mass_unit',db,errcount=nerr)
 call read_inopt(dist_unit,'dist_unit',db,errcount=nerr)
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
       call read_inopt(ibinary,'ibinary',db,min=0,max=1,errcount=nerr)
       select case (ibinary)
       case (0)
          !--bound
          call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
          call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
          call read_inopt(binary_a,'binary_a',db,errcount=nerr)
          call read_inopt(binary_e,'binary_e',db,min=0.,errcount=nerr)
          call read_inopt(binary_i,'binary_i',db,errcount=nerr)
          call read_inopt(binary_O,'binary_O',db,errcount=nerr)
          call read_inopt(binary_w,'binary_w',db,errcount=nerr)
          call read_inopt(binary_f,'binary_f',db,errcount=nerr)
          call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
          call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
       case (1)
          !--unbound (flyby)
          !--central star
          call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
          call read_inopt(accr1,'accr1',db,min=0.,errcount=nerr)
          !--perturber
          call read_inopt(m2,'m2',db,min=0.,errcount=nerr)
          call read_inopt(accr2,'accr2',db,min=0.,errcount=nerr)
          call read_inopt(flyby_a,'flyby_a',db,min=0.,errcount=nerr)
          call read_inopt(flyby_d,'flyby_d',db,min=0.,errcount=nerr)
          call read_inopt(flyby_O,'flyby_O',db,min=0.,errcount=nerr)
          call read_inopt(flyby_i,'flyby_i',db,min=0.,errcount=nerr)
       end select
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
       call read_inopt(binary_a,'binary_a',db,errcount=nerr)
       call read_inopt(binary_e,'binary_e',db,errcount=nerr)
       call read_inopt(binary_i,'binary_i',db,errcount=nerr)
       call read_inopt(binary_O,'binary_O',db,errcount=nerr)
       call read_inopt(binary_w,'binary_w',db,errcount=nerr)
       call read_inopt(binary_f,'binary_f',db,errcount=nerr)

       !-- tight parameters
       call read_inopt(binary2_a,'binary2_a',db,errcount=nerr)
       call read_inopt(binary2_e,'binary2_e',db,errcount=nerr)
       call read_inopt(binary2_i,'binary2_i',db,errcount=nerr)
       call read_inopt(binary2_O,'binary2_O',db,errcount=nerr)
       call read_inopt(binary2_w,'binary2_w',db,errcount=nerr)
       call read_inopt(binary2_f,'binary2_f',db,errcount=nerr)

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
       call read_inopt(binary_a,'binary_a',db,errcount=nerr)
       call read_inopt(binary_e,'binary_e',db,errcount=nerr)
       call read_inopt(binary_i,'binary_i',db,errcount=nerr)
       call read_inopt(binary_O,'binary_O',db,errcount=nerr)
       call read_inopt(binary_w,'binary_w',db,errcount=nerr)
       call read_inopt(binary_f,'binary_f',db,errcount=nerr)

       !-- tight binary 1 parameters
       call read_inopt(binary1_a,'binary1_a',db,errcount=nerr)
       call read_inopt(binary1_e,'binary1_e',db,errcount=nerr)
       call read_inopt(binary1_i,'binary1_i',db,errcount=nerr)
       call read_inopt(binary1_O,'binary1_O',db,errcount=nerr)
       call read_inopt(binary1_w,'binary1_w',db,errcount=nerr)
       call read_inopt(binary1_f,'binary1_f',db,errcount=nerr)

       !-- tight binary 2 parameters
       call read_inopt(binary2_a,'binary2_a',db,errcount=nerr)
       call read_inopt(binary2_e,'binary2_e',db,errcount=nerr)
       call read_inopt(binary2_i,'binary2_i',db,errcount=nerr)
       call read_inopt(binary2_O,'binary2_O',db,errcount=nerr)
       call read_inopt(binary2_w,'binary2_w',db,errcount=nerr)
       call read_inopt(binary2_f,'binary2_f',db,errcount=nerr)

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

 call read_inopt(discstrat,'discstrat',db,errcount=nerr)
 call read_inopt(lumdisc,'lumdisc',db,errcount=nerr)

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
       if (ibinary==0) then
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
       call read_inopt(posangl(i),'posangl'//trim(disclabel),db,min=0.,max=360.,errcount=nerr)
       call read_inopt(incl(i),'incl'//trim(disclabel),db,min=0.,max=180.,errcount=nerr)
       if (discstrat == 0 .and. lumdisc == 0) then
          call read_inopt(H_R(i),'H_R'//trim(disclabel),db,min=0.,errcount=nerr)
       endif
       call read_inopt(iwarp(i),'iwarp'//trim(disclabel),db,errcount=nerr)
       if (iwarp(i)) then
          call read_inopt(R_warp(i),'R_warp'//trim(disclabel),db,min=0.,errcount=nerr)
          call read_inopt(H_warp(i),'H_warp'//trim(disclabel),db,min=0.,errcount=nerr)
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
 if (maxalpha==0 .and. any(iuse_disc)) call read_inopt(alphaSS,'alphaSS',db,min=0.,errcount=nerr)
 !--planets
 call read_inopt(nplanets,'nplanets',db,min=0,max=maxplanets,errcount=nerr)
 do i=1,nplanets
    call read_inopt(mplanet(i),'mplanet'//trim(num(i)),db,min=0.,errcount=nerr)
    call read_inopt(rplanet(i),'rplanet'//trim(num(i)),db,min=0.,errcount=nerr)
    call read_inopt(inclplan(i),'inclplanet'//trim(num(i)),db,min=0.,max=180.,errcount=nerr)
    call read_inopt(accrplanet(i),'accrplanet'//trim(num(i)),db,min=0.,errcount=nerr)
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
 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
 endif

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

real function get_cs_from_lum(L_star,r)
 use physcon, only:kb_on_mh,steboltz,solarl,fourpi
 use units,   only:udist,unit_velocity
 real,intent(in) :: L_star,r
 real :: mu

 mu = 2.381 !mean molecular mass
 get_cs_from_lum = sqrt(kb_on_mh/mu) * ( (L_star*solarl/(fourpi*steboltz))**0.125 / &
               (r*udist)**0.25 + sqrt(T_bg) )
 get_cs_from_lum = get_cs_from_lum/unit_velocity
end function get_cs_from_lum

end module setup
