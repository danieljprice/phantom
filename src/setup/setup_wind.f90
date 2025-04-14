!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! initial conditions for binary wind accretion / AGB star wind injection
!
! :References:
!   Siess et al. 2022, A&A, 667, 75
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - Reff2a            : *tight binary primary effective radius (au)*
!   - Reff2b            : *tight binary secondary effective radius (au)*
!   - T_wind            : *wind temperature (K)*
!   - Teff2a            : *tight binary primary effective temperature (K)*
!   - Teff2b            : *tight binary secondary effective temperature (K)*
!   - binary2_a         : *tight binary semi-major axis*
!   - binary2_e         : *tight binary eccentricity*
!   - eccentricity      : *eccentricity of the binary system*
!   - icompanion_star   : *set to 1 for a binary system, 2 for a triple system*
!   - inclination       : *inclination of the tight binary system w.r.t. outer binary (deg)*
!   - lum2a             : *tight binary primary luminosity (Lsun)*
!   - lum2b             : *tight binary secondary luminosity (Lsun)*
!   - mass_of_particles : *particle mass (Msun, overwritten if iwind_resolution <>0)*
!   - primary_Reff      : *primary star effective radius (au)*
!   - primary_Teff      : *primary star effective temperature (K)*
!   - primary_lum       : *primary star luminosity (Lsun)*
!   - primary_mass      : *primary star mass (Msun)*
!   - primary_racc      : *primary star accretion radius (au)*
!   - q2                : *tight binary mass ratio*
!   - racc2a            : *tight binary primary accretion radius*
!   - racc2b            : *tight binary secondary accretion radius*
!   - secondary_Reff    : *secondary star effective radius (au)*
!   - secondary_Teff    : *secondary star effective temperature)*
!   - secondary_lum     : *secondary star luminosity (Lsun)*
!   - secondary_mass    : *secondary star mass (Msun)*
!   - secondary_racc    : *secondary star accretion radius (au)*
!   - semi_major_axis   : *semi-major axis of the binary system (au)*
!   - subst             : *star to substitute*
!   - temp_exponent     : *temperature profile T(r) = T_wind*(r/Reff)^(-temp_exponent)*
!   - wind_gamma        : *adiabatic index (initial if Krome chemistry used)*
!
! :Dependencies: dim, eos, infile_utils, inject, io, part, physcon,
!   prompting, setbinary, sethierarchical, spherical, units
!
 use dim, only:isothermal
 implicit none
 public :: setpart

 private
 real, public :: wind_gamma
 real, public :: T_wind
 real :: temp_exponent
 integer :: icompanion_star,iwind
 real :: semi_major_axis,semi_major_axis_au,eccentricity
 real :: default_particle_mass
 real :: primary_lum_lsun,primary_mass_msun,primary_Reff_au,primary_racc_au
 real :: secondary_lum_lsun,secondary_mass_msun,secondary_Reff_au,secondary_racc_au
 real :: lum2a_lsun,lum2b_lsun,Teff2a,Teff2b,Reff2a_au,Reff2b_au
 real :: binary2_a_au,racc2a_au,racc2b_au,binary2_i,q2
 real :: primary_Reff,primary_Teff,primary_lum,primary_mass,primary_racc
 real :: secondary_Reff,secondary_Teff,secondary_lum,secondary_mass,secondary_racc
 real :: Reff2a,Reff2b
 real :: racc2a,racc2b
 real :: lum2a,lum2b
 real :: binary2_a
 real :: binary2_e
 integer :: subst

contains
!----------------------------------------------------------------
!+
!  default parameter choices
!+
!----------------------------------------------------------------
subroutine set_default_parameters_wind()

 wind_gamma    = 5./3.
 if (isothermal) then
    T_wind                = 100000.
    temp_exponent         = 0.5
    ! primary_racc_au       = 0.465
    ! primary_mass_msun     = 1.5
    ! primary_lum_lsun      = 0.
    ! primary_Reff_au       = 0.465240177008 !100 Rsun
 else
    T_wind = 3000.
    !primary_racc_au       = 1.
    !primary_mass_msun     = 1.5
    !primary_lum_lsun      = 20000.
    !primary_Reff_au       = 0.
 endif
 icompanion_star = 0
 semi_major_axis       = 4.0
 eccentricity          = 0.
 primary_Teff          = 3000.
 secondary_Teff        = 0.
 semi_major_axis_au    = 4.0
 default_particle_mass = 1.e-11
 primary_lum_lsun      = 5315.
 primary_mass_msun     = 1.5
 primary_Reff_au       = 1.
 primary_racc_au       = 1.
 secondary_lum_lsun    = 0.
 secondary_mass_msun   = 1.0
 secondary_Reff_au     = 0.
 secondary_racc_au     = 0.1
 lum2a_lsun            = 0.
 lum2b_lsun            = 0.
 Teff2a                = 0.
 Teff2b                = 0.
 Reff2a_au             = 0.
 Reff2b_au             = 0.
 binary2_a_au          = 0.3
 racc2a_au             = 0.1
 racc2b_au             = 0.1
 binary2_i             = 0.

end subroutine set_default_parameters_wind

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only: xyzmh_ptmass, vxyz_ptmass, nptmass, igas, iTeff, iLum, iReff
 use physcon,   only: au, solarm, mass_proton_cgs, kboltz, solarl
 use units,     only: umass,set_units,unit_velocity,utime,unit_energ,udist
 use inject,    only: set_default_options_inject
 use setbinary, only: set_binary
 use sethierarchical, only: set_multiple
 use io,        only: master
 use eos,       only: gmw,ieos,isink,qfacdisc
 use spherical, only: set_sphere
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 character(len=len(fileprefix)+6) :: filename
 integer :: ierr,k
 logical :: iexist

 call set_units(dist=au,mass=solarm,G=1.)
 call set_default_parameters_wind()
 filename = trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) call set_default_options_inject()

!--general parameters
!
 time = 0.
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call setup_interactive()
       call write_setupfile(filename)
    endif
 endif

!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 if (icompanion_star == 1) then
    call set_binary(primary_mass, &
                    secondary_mass, &
                    semi_major_axis, &
                    eccentricity, &
                    primary_racc, &
                    secondary_racc, &
                    xyzmh_ptmass, vxyz_ptmass, nptmass, ierr)
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
    xyzmh_ptmass(iTeff,2) = secondary_Teff
    xyzmh_ptmass(iReff,2) = secondary_Reff
    xyzmh_ptmass(iLum,2)  = secondary_lum
 elseif (icompanion_star == 2) then
    !-- hierarchical triple
    nptmass  = 0
    print "(/,a)",'----------- Hierarchical triple -----------'
    print "(a,g10.3,a)",'     First hierarchical level primary mass: ', primary_mass_msun
    print "(a,g10.3,a)",'   First hierarchical level secondary mass: ', secondary_mass_msun
    print "(a,g10.3)",  '                    Wide binary mass ratio: ', secondary_mass/primary_mass
    print "(a,g10.3)",  '                   Tight binary mass ratio: ', q2
    print "(a,g10.3)",  '                    Star to be substituted: ', abs(subst)
!        print "(a,g10.3,a)",'                        Accretion Radius 1: ', primary_racc!, trim(dist_unit)
!        print "(a,g10.3,a)",'                       Accretion Radius 2a: ', racc2a!, trim(dist_unit)
!        print "(a,g10.3,a)",'                       Accretion Radius 2b: ', racc2b!, trim(dist_unit)

    if (subst>0) then
       print "(a,g10.3,a)",'      Tight binary orientation referred to: substituted star orbital plane'
    else
       print "(a,g10.3,a)",'      Tight binary orientation referred to: sky'
    endif


    call set_multiple(primary_mass,secondary_mass,semimajoraxis=semi_major_axis,eccentricity=eccentricity, &
            accretion_radius1=primary_racc,accretion_radius2=secondary_racc, &
            xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,ierr=ierr)

    if (subst == 12) then
       call set_multiple(secondary_mass/(q2+1),secondary_mass*q2/(q2+1),semimajoraxis=binary2_a,eccentricity=binary2_e, &
                accretion_radius1=racc2a,accretion_radius2=racc2b, &
                xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,&
                posang_ascnode=0.,arg_peri=0.,incl=binary2_i,subst=subst,ierr=ierr)

       xyzmh_ptmass(iTeff,1) = primary_Teff
       xyzmh_ptmass(iReff,1) = primary_Reff
       xyzmh_ptmass(iLum,1)  = primary_lum
       xyzmh_ptmass(iTeff,2) = Teff2a
       xyzmh_ptmass(iReff,2) = Reff2a
       xyzmh_ptmass(iLum,2)  = lum2a
       xyzmh_ptmass(iTeff,3) = Teff2b
       xyzmh_ptmass(iReff,3) = Reff2b
       xyzmh_ptmass(iLum,3)  = lum2b

    elseif (subst == 11) then
       call set_multiple(primary_mass*q2/(q2+1),primary_mass/(q2+1),semimajoraxis=binary2_a,eccentricity=binary2_e, &
                accretion_radius1=racc2b,accretion_radius2=primary_racc, &
                xyzmh_ptmass=xyzmh_ptmass,vxyz_ptmass=vxyz_ptmass,nptmass=nptmass,&
                posang_ascnode=0.,arg_peri=0.,incl=binary2_i,subst=subst,ierr=ierr)

       xyzmh_ptmass(iTeff,1) = primary_Teff
       xyzmh_ptmass(iReff,1) = primary_Reff
       xyzmh_ptmass(iLum,1)  = primary_lum
       xyzmh_ptmass(iTeff,2) = secondary_Teff
       xyzmh_ptmass(iReff,2) = secondary_Reff
       xyzmh_ptmass(iLum,2)  = secondary_lum
       xyzmh_ptmass(iTeff,3) = Teff2b
       xyzmh_ptmass(iReff,3) = Reff2b
       xyzmh_ptmass(iLum,3)  = lum2b
    endif

    print *,'Sink particles summary'
    print *,'  #    mass       racc      lum         Reff'
    do k=1,nptmass
       print '(i4,2(2x,f9.5),2(2x,es10.3))',k,xyzmh_ptmass(4:5,k),xyzmh_ptmass(iLum,k)/(solarl*utime/unit_energ),&
               xyzmh_ptmass(iReff,k)*udist/au
    enddo
    print *,''

 else
    nptmass = 1
    xyzmh_ptmass(4,1)     = primary_mass
    xyzmh_ptmass(5,1)     = primary_racc
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
 endif

 !
 ! for binary wind simulations the particle mass is IRRELEVANT
 ! since it will be over-written on the first call to init_inject
 !
 massoftype(igas) = default_particle_mass * (solarm / umass)

 if (isothermal) then
    gamma = 1.
    if (iwind == 3) then
       ieos = 6
       qfacdisc = 0.5*temp_exponent
       isink = 1
       T_wind = primary_Teff
    else
       isink = 1
       ieos = 1
    endif
 else
    T_wind = 0.
    gamma = wind_gamma
 endif
 polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)

end subroutine setpart

!----------------------------------------------------------------
!+
!  determine which problem to set up interactively
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use physcon,   only:au,solarm
 use units,     only:umass,udist
 use io,        only:fatal
 integer :: ichoice

 if (isothermal) then
    iwind = 2
 else
    iwind = 1
    call prompt('Type of wind:  1=adia, 2=isoT, 3=T(r)',iwind,1,3)
    if (iwind == 2 .or. iwind == 3) then
       call fatal('setup','If you choose options 2 or 3, the code must be compiled with SETUP=isowind')
    endif
    if (iwind == 3) T_wind = primary_Teff
 endif

 icompanion_star = 0
 call prompt('Add binary?',icompanion_star,0,2)

 !Hierarchical triple system
 if (icompanion_star == 2) then
    !select the tight binary
    ichoice = 1
    print "(a)",'Star to be substituted by a tight binary'
    print "(a)",' 1: primary (2+1)' ,' 2: companion (1+2)'
    call prompt('Select star to be substituted',ichoice,1,2)
    subst = ichoice+10

    !select orbital parameters for outer binary
    semi_major_axis_au = 15.
    eccentricity       = 0.
    ichoice = 1
    print "(a)",'Orbital parameters first hierarchical level binary'
    print "(a)",' 1: semi-axis = 15 au, eccentricity = 0',' 0: custom'
    call prompt('select semi-major axis and ecccentricity',ichoice,0,1)
    if (ichoice == 0) then
       call prompt('enter semi-major axis in au',semi_major_axis_au,0.,100.)
       call prompt('enter eccentricity',eccentricity,0.)
    endif
    semi_major_axis = semi_major_axis_au * au / udist
    ichoice = 1

    !replace companion by tight binary system : 1+2
    if (subst == 12) then
       print "(a)",'Primary star parameters (the single wind launching central star)'
       print "(a)",' 2: Mass = 1.2 Msun, accretion radius = 0.2568 au',&
        ' 1: Mass = 1.5 Msun, accretion radius = 1.2568 au', &
        ' 0: custom'
       call prompt('select mass and radius of primary',ichoice,0,2)
       select case(ichoice)
       case(2)
          primary_mass_msun = 1.2
          primary_racc_au   = 0.2568
       case(1)
          primary_mass_msun = 1.5
          primary_racc_au   = 1.2568
       case default
          primary_mass_msun = 1.5
          primary_racc_au   = 1.
          call prompt('enter primary mass',primary_mass_msun,0.,100.)
          call prompt('enter accretion radius in au ',primary_racc_au,0.)
       end select
       primary_mass = primary_mass_msun * (solarm / umass)
       primary_racc = primary_racc_au * (au / udist)

       ichoice = 1
       print "(a)",'Total mass of tight binary system (1+2)'
       print "(a)",' 1: Total mass tight binary = 1.0 Msun',' 0: custom'
       secondary_mass_msun = 1.
       call prompt('select mass',ichoice,0,1)
       select case(ichoice)
       case(0)
          call prompt('enter total mass tigh binary',secondary_mass_msun,0.,100.)
       end select
       secondary_mass = secondary_mass_msun * (solarm / umass)

       ichoice = 1
       print "(a)",'Mass ratio and accretion radii of stars in tight orbit:'
       print "(a)",' 1: mass ratio m2b/m2a = 1, accretion radius a = 0.01 au, accretion radius b = 0.01 au',' 0: custom'
       call prompt('select mass ratio and accretion radii of tight binary',ichoice,0,1)
       select case(ichoice)
       case(1)
          q2 = 1.
          racc2a_au = 0.1
          racc2b_au = 0.1
       case default
          q2 = 1.
          racc2a_au = 0.1
          racc2b_au = 0.1
          call prompt('enter tight binary mass ratio',q2,0.)
          call prompt('enter accretion radius a in au ',racc2a_au,0.)
          call prompt('enter accretion radius b in au ',racc2b_au,0.)
       end select
       racc2a = racc2a_au * (au / udist)
       racc2b = racc2b_au * (au / udist)
       secondary_racc = racc2a !needs to be /=0 otherwise NaNs in set_multiple

       !replace primary by tight binary system : 2+1
    elseif (subst == 11) then
       print "(a)",'Stellar parameters of the remote single star (2+1)'
       print "(a)",' 1: Mass = 1.0 Msun, accretion radius = 0.1 au',' 0: custom'
       call prompt('select mass and radius of remote single star',ichoice,0,1)
       select case(ichoice)
       case(1)
          secondary_mass_msun = 1.
          secondary_racc_au   = 0.1
       case default
          secondary_mass_msun = 1.
          secondary_racc_au   = 0.1
          call prompt('enter mass of remote single star',secondary_mass_msun,0.,100.)
          call prompt('enter accretion radius in au ',secondary_racc_au,0.)
       end select
       secondary_mass = secondary_mass_msun * (solarm / umass)
       secondary_racc = secondary_racc_au * (au / udist)

       ichoice = 1
       print "(a)",'wind-launching star accretion radius in tigh orbit (called primary)'
       print "(a)",' 2: accretion radius primary = 0.2568 au',&
        ' 1: accretion radius primary = 1.2568 au', &
        ' 0: custom'
       call prompt('select accretion radius of wind launching star',ichoice,0,2)
       select case(ichoice)
       case(2)
          primary_racc_au = 0.2568
       case(1)
          primary_racc_au = 1.2568
       case default
          primary_racc_au = 1.
          call prompt('enter accretion radius in au ',primary_racc_au,0.)
       end select
       primary_racc = primary_racc_au * (au / udist)

       ichoice = 1
       print "(a)",'Total mass of the tight binary system (2+1):'
       print "(a)",' 2: Total mass tight binary = 1.2 Msun',&
        ' 1: Total mass tight binary = 1.5 Msun', &
        ' 0: custom'
       call prompt('select total mass tight binary',ichoice,0,2)
       select case(ichoice)
       case(2)
          primary_mass_msun = 1.2
       case(1)
          primary_mass_msun = 1.5
       case default
          primary_mass_msun = 1.5
          call prompt('enter primary mass',primary_mass_msun,0.,100.)
       end select
       primary_mass = primary_mass_msun * (solarm / umass)

       ichoice = 1
       print "(a)",'Mass ratio and accretion radius of secondary in tight orbit:'
       print "(a)",' 1: mass ratio m1b/m1a = 0.3, accretion radius b = 0.01 au',' 0: custom'
       call prompt('select mass ratio and accretion radius of tight binary',ichoice,0,1)
       select case(ichoice)
       case(1)
          q2 = 0.3
          racc2b_au = 0.1
       case default
          q2 = 0.3
          racc2b_au = 0.1
          call prompt('enter tight binary mass ratio',q2,0.)
          call prompt('enter accretion radius b in au ',racc2b_au,0.)
       end select
       racc2b = racc2b_au * (au / udist)
    endif

    ichoice = 1
    print "(a)",'Orbital parameters of tight system:'
    print "(a)",' 1: semi-axis = 4 au, eccentricity = 0',' 0: custom'
    call prompt('select tight binary semi-major axis and eccentricity',ichoice,0,1)
    select case(ichoice)
    case(1)
       binary2_a_au = 4.
       binary2_e    = 0.
    case default
       binary2_a_au = 4.
       binary2_e    = 0.
       call prompt('enter semi-major axis in au',binary2_a_au,0.,semi_major_axis_au)
       call prompt('enter eccentricity',binary2_e,0.)
    end select
    binary2_a = binary2_a_au * au / udist

    ichoice = 1
    print "(a)",'inclination of orbit tight binary w.r.t. outer binary:'
    print "(a)",' 1: inclination = 0 deg',' 0: custom'
    call prompt('select inclination',ichoice,0,1)
    select case(ichoice)
    case(1)
       binary2_i = 0.
    case default
       binary2_i = 0.
       call prompt('enter inclination',binary2_i,0.,90.)
    end select

    !binary or single star case
 else
    if (icompanion_star == 1) then
       print "(a)",'Primary star parameters'
    else
       print "(a)",'Stellar parameters'
    endif
    ichoice = 2
    print "(a)",' 3: Mass = 1.2 Msun, accretion radius = 1. au (trans-sonic)',&
         ' 2: Mass = 1.2 Msun, accretion radius = 0.2568 au',&
         ' 1: Mass = 1.0 Msun, accretion radius = 1.2568 au', &
        ' 0: custom'
    call prompt('select mass and radius of primary',ichoice,0,3)
    select case(ichoice)
    case(3)
       primary_lum_lsun  = 2.d4
       primary_Teff      = 5.d4
       primary_mass_msun = 1.2
       primary_racc_au   = 1.
       wind_gamma = 1.4
    case(2)
       primary_mass_msun = 1.2
       primary_racc_au   = 0.2568
    case(1)
       primary_mass_msun = 1.
       primary_racc_au   = 1.2568
    case default
       primary_mass_msun = 1.
       primary_racc_au   = 1.
       call prompt('enter primary mass',primary_mass_msun,0.,100.)
       call prompt('enter accretion radius in au ',primary_racc_au,0.)
    end select
    primary_mass = primary_mass_msun * (solarm / umass)
    primary_racc = primary_racc_au * (au / udist)

    if (icompanion_star == 1) then
       ichoice = 1
       print "(a)",'Secondary star parameters'
       print "(a)",' 1: Mass = 1.0 Msun, accretion radius = 0.1 au',' 0: custom'
       call prompt('select mass and radius of secondary',ichoice,0,1)
       select case(ichoice)
       case(1)
          secondary_mass_msun = 1.
          secondary_racc_au   = 0.1
       case default
          secondary_mass_msun = 1.
          secondary_racc_au   = 0.1
          call prompt('enter secondary mass',secondary_mass_msun,0.,100.)
          call prompt('enter accretion radius in au ',secondary_racc_au,0.)
       end select
       secondary_mass = secondary_mass_msun * (solarm / umass)
       secondary_racc = secondary_racc_au * (au / udist)

       ichoice = 1
       print "(a)",'Orbital parameters'
       print "(a)",' 1: semi-axis = 3.7 au, eccentricity = 0',' 0: custom'
       call prompt('select semi-major axis and ecccentricity',ichoice,0,1)
       select case(ichoice)
       case(1)
          semi_major_axis_au = 3.7
          eccentricity       = 0.
       case default
          semi_major_axis_au = 1.
          eccentricity       = 0.
          call prompt('enter semi-major axis in au',semi_major_axis_au,0.,100.)
          call prompt('enter eccentricity',eccentricity,0.)
       end select
       semi_major_axis = semi_major_axis_au * au / udist
    endif
 endif

end subroutine setup_interactive

!----------------------------------------------------------------
!+
!  get luminosity and effective radius in code units
!  from various combinations of L, Teff and Reff in physical inuts
!+
!----------------------------------------------------------------
subroutine get_lum_and_Reff(lum_lsun,reff_au,Teff,lum,Reff)
 use physcon, only:au,steboltz,solarl,pi
 use units,   only:udist,unit_luminosity
 real, intent(inout) :: lum_lsun,reff_au,Teff
 real, intent(out)   :: lum,Reff

 if (Teff <= tiny(0.) .and. lum_lsun > 0. .and. Reff_au > 0.) then
    primary_Teff = (lum_lsun*solarl/(4.*pi*steboltz*(Reff_au*au)**2))**0.25
 elseif (Reff_au <= 0. .and. lum_lsun > 0. .and. Teff > 0.) then
    Reff_au = sqrt(lum_lsun*solarl/(4.*pi*steboltz*Teff**4))/au
 elseif (Reff_au > 0. .and. lum_lsun <= 0. .and. Teff > 0.) then
    lum_lsun = 4.*pi*steboltz*Teff**4*(primary_Reff_au*au)**2/solarl
 endif

 lum  = lum_lsun*(solarl/unit_luminosity)
 Reff = Reff_au*(au/udist)

end subroutine get_lum_and_Reff

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wind setup routine'

 call get_lum_and_Reff(primary_lum_lsun,primary_Reff_au,primary_Teff,primary_lum,primary_Reff)

 if (icompanion_star == 2) then
    call get_lum_and_Reff(secondary_lum_lsun,secondary_Reff_au,secondary_Teff,secondary_lum,secondary_Reff)

    call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system, 2 for a triple system',iunit)
    !-- hierarchical triple
    write(iunit,"(/,a)") '# options for hierarchical triple'
    call write_inopt(subst,'subst','star to substitute',iunit)
    write(iunit,"(/,a)") '# input of primary (wind launching star)'
    if (subst == 12) then
       call write_inopt(primary_mass_msun,'primary_mass','primary star mass (Msun)',iunit)
       call write_inopt(primary_racc_au,'primary_racc','primary star accretion radius (au)',iunit)
       call write_inopt(primary_lum_lsun,'primary_lum','primary star luminosity (Lsun)',iunit)
       call write_inopt(primary_Teff,'primary_Teff','primary star effective temperature (K)',iunit)
       call write_inopt(primary_Reff_au,'primary_Reff','primary star effective radius (au)',iunit)
       call write_inopt(semi_major_axis_au,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
       call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
       write(iunit,"(/,a)") '# input secondary to be replaced by tight binary'
       call write_inopt(secondary_mass_msun,'secondary_mass','total mass of secondary tight binary (Msun)',iunit)
       call write_inopt(q2,'q2','tight binary mass ratio',iunit)
       !-- tight orbital parameters
       call write_inopt(binary2_a,'binary2_a','tight binary semi-major axis',iunit)
       call write_inopt(binary2_e,'binary2_e','tight binary eccentricity',iunit)
       !-- accretion radii, luminosity, radii
       call write_inopt(racc2a_au,'racc2a','tight binary primary accretion radius',iunit)
       call write_inopt(racc2b_au,'racc2b','tight binary secondary accretion radius',iunit)
       call write_inopt(lum2a_lsun,'lum2a','tight binary primary luminosity (Lsun)',iunit)
       call write_inopt(lum2b_lsun,'lum2b','tight binary secondary luminosity (Lsun)',iunit)
       call write_inopt(Teff2a,'Teff2a','tight binary primary effective temperature (K)',iunit)
       call write_inopt(Teff2b,'Teff2b','tight binary secondary effective temperature (K)',iunit)
       call write_inopt(Reff2a_au,'Reff2a','tight binary primary effective radius (au)',iunit)
       call write_inopt(Reff2b_au,'Reff2b','tight binary secondary effective radius (au)',iunit)
    elseif (subst == 11) then
       call write_inopt(primary_racc_au,'primary_racc','primary star accretion radius (au)',iunit)
       call write_inopt(primary_lum_lsun,'primary_lum','primary star luminosity (Lsun)',iunit)
       call write_inopt(primary_Teff,'primary_Teff','primary star effective temperature (K)',iunit)
       call write_inopt(primary_Reff_au,'primary_Reff','primary star effective radius (au)',iunit)
       write(iunit,"(/,a)") '# input tight binary to create close companion'
       call write_inopt(primary_mass_msun,'primary_mass','primary star mass (Msun)',iunit)
       call write_inopt(q2,'q2','tight binary mass ratio',iunit)
       !-- tight orbital parameters
       call write_inopt(binary2_a,'binary2_a','tight binary semi-major axis',iunit)
       call write_inopt(binary2_e,'binary2_e','tight binary eccentricity',iunit)
       !-- accretion radii
       call write_inopt(racc2b_au,'racc2b','tight binary secondary accretion radius',iunit)
       call write_inopt(lum2b_lsun,'lum2b','tight binary secondary luminosity (Lsun)',iunit)
       call write_inopt(Teff2b,'Teff2b','tight binary secondary effective temperature (K)',iunit)
       call write_inopt(Reff2b_au,'Reff2b','tight binary secondary effective radius (au)',iunit)
       write(iunit,"(/,a)") '# input of secondary, outer binary'
       call write_inopt(secondary_mass_msun,'secondary_mass','secondary star mass (Msun)',iunit)
       call write_inopt(secondary_racc_au,'secondary_racc','secondary star accretion radius (au)',iunit)
       call write_inopt(secondary_lum_lsun,'secondary_lum','secondary star luminosity (Lsun)',iunit)
       call write_inopt(secondary_Teff,'secondary_Teff','secondary star effective temperature)',iunit)
       call write_inopt(secondary_Reff_au,'secondary_Reff','secondary star effective radius (au)',iunit)
       call write_inopt(semi_major_axis_au,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
       call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
    endif
    call write_inopt(binary2_i,'inclination','inclination of the tight binary system w.r.t. outer binary (deg)',iunit)
    !binary or single star
 else
    call write_inopt(primary_mass_msun,'primary_mass','primary star mass (Msun)',iunit)
    call write_inopt(primary_racc_au,'primary_racc','primary star accretion radius (au)',iunit)
    call write_inopt(primary_lum_lsun,'primary_lum','primary star luminosity (Lsun)',iunit)
    call write_inopt(primary_Teff,'primary_Teff','primary star effective temperature (K)',iunit)
    call write_inopt(primary_Reff_au,'primary_Reff','primary star effective radius (au)',iunit)
    call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system, 2 for a triple system',iunit)
    if (icompanion_star == 1) then
       call get_lum_and_Reff(secondary_lum_lsun,secondary_Reff_au,secondary_Teff,secondary_lum,secondary_Reff)

       call write_inopt(secondary_mass_msun,'secondary_mass','secondary star mass (Msun)',iunit)
       call write_inopt(secondary_racc_au,'secondary_racc','secondary star accretion radius (au)',iunit)
       call write_inopt(secondary_lum_lsun,'secondary_lum','secondary star luminosity (Lsun)',iunit)
       call write_inopt(secondary_Teff,'secondary_Teff','secondary star effective temperature)',iunit)
       call write_inopt(secondary_Reff_au,'secondary_Reff','secondary star effective radius (au)',iunit)
       call write_inopt(semi_major_axis_au,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
       call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
    endif
 endif

 call write_inopt(default_particle_mass,'mass_of_particles','particle mass (Msun, overwritten if iwind_resolution <>0)',iunit)

 if (isothermal) then
    wind_gamma = 1.
    if (iwind == 3) then
       call write_inopt(primary_Teff,'T_wind','wind temperature at injection radius (K)',iunit)
       call write_inopt(temp_exponent,'temp_exponent','temperature profile T(r) = T_wind*(r/Reff)^(-temp_exponent)',iunit)
    else
       call write_inopt(T_wind,'T_wind','wind temperature (K)',iunit)
    endif
 else
    call write_inopt(wind_gamma,'wind_gamma','adiabatic index (initial if Krome chemistry used)',iunit)
 endif
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use physcon,      only:au,steboltz,solarl,solarm,pi
 use units,        only:udist,umass,utime,unit_energ
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr,ichange

 nerr = 0
 ichange = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(primary_mass_msun,'primary_mass',db,min=0.,max=1000.,errcount=nerr)
 primary_mass = primary_mass_msun * (solarm / umass)
 call read_inopt(primary_lum_lsun,'primary_lum',db,min=0.,max=1e7,errcount=nerr)
 primary_lum = primary_lum_lsun * (solarl * utime / unit_energ)
 call read_inopt(primary_Teff,'primary_Teff',db,min=0.,errcount=nerr)
 call read_inopt(primary_Reff_au,'primary_Reff',db,min=0.,errcount=nerr)
 primary_Reff = primary_Reff_au * au / udist
 call read_inopt(primary_racc_au,'primary_racc',db,min=0.,errcount=nerr)
 primary_racc = primary_racc_au * au / udist
 if (primary_racc < tiny(0.)) then
    print *,'ERROR: primary accretion radius not defined'
    nerr = nerr+1
 endif

 call read_inopt(icompanion_star,'icompanion_star',db,min=0,errcount=nerr)
 if (icompanion_star == 1) then
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(secondary_lum_lsun,'secondary_lum',db,min=0.,max=1e7,errcount=nerr)
    secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
    call read_inopt(secondary_Teff,'secondary_Teff',db,min=0.,errcount=nerr)
    call read_inopt(secondary_Reff_au,'secondary_Reff',db,min=0.,errcount=nerr)
    secondary_Reff = secondary_Reff_au * au / udist
    call read_inopt(secondary_racc_au,'secondary_racc',db,min=0.,errcount=nerr)
    secondary_racc = secondary_racc_au * au / udist
    if (secondary_racc < tiny(0.)) then
       print *,'ERROR: secondary accretion radius not defined'
       nerr = nerr+1
    endif
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
 elseif (icompanion_star == 2) then
    !-- hierarchical triple
    call read_inopt(subst,'subst',db,errcount=nerr)
    !replace primary by tight binary system : 2+1
    if (subst == 11) then
       call read_inopt(secondary_lum_lsun,'secondary_lum',db,min=0.,max=1000.,errcount=nerr)
       secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
       call read_inopt(secondary_Teff,'secondary_Teff',db,min=0.,max=1000.,errcount=nerr)
       call read_inopt(secondary_Reff_au,'secondary_Reff',db,min=0.,max=1000.,errcount=nerr)
       secondary_Reff = secondary_Reff_au * au / udist
       call read_inopt(secondary_racc_au,'secondary_racc',db,min=0.,max=1000.,errcount=nerr)
       secondary_racc = secondary_racc_au * au / udist
    elseif (subst == 12) then
       call read_inopt(lum2a_lsun,'lum2a',db,errcount=nerr)
       lum2a = lum2a_lsun * (solarl * utime / unit_energ)
       !secondary_lum_lsun = lum2a_lsun
       call read_inopt(Teff2a,'Teff2a',db,errcount=nerr)
       call read_inopt(Reff2a_au,'Reff2a',db,errcount=nerr)
       Reff2a = Reff2a_au * au / udist
       !secondary_Reff =  Reff2a
       call read_inopt(racc2a_au,'racc2a',db,errcount=nerr)
       racc2a = racc2a_au * au / udist
    endif
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
    !-- masses
    call read_inopt(q2,'q2',db,min=0.,max=1.,errcount=nerr)
    !-- tight parameters
    call read_inopt(binary2_a_au,'binary2_a',db,errcount=nerr)
    binary2_a = binary2_a_au * au / udist
    call read_inopt(binary2_e,'binary2_e',db,errcount=nerr)
    !-- accretion radii,...
    call read_inopt(racc2b_au,'racc2b',db,errcount=nerr)
    racc2b = racc2b_au * au / udist
    if (racc2b < tiny(0.)) then
       print *,'WARNING: secondary accretion radius not defined'
       !nerr = nerr+1
    endif
    call read_inopt(lum2b_lsun,'lum2b',db,errcount=nerr)
    lum2b = lum2b_lsun * (solarl * utime / unit_energ)
    call read_inopt(Teff2b,'Teff2b',db,errcount=nerr)
    call read_inopt(Reff2b_au,'Reff2b',db,errcount=nerr)
    Reff2b = Reff2b_au * au / udist
    call read_inopt(binary2_i,'inclination',db,errcount=nerr)
 endif

 call read_inopt(default_particle_mass,'mass_of_particles',db,min=0.,errcount=nerr)

 if (isothermal) then
    wind_gamma = 1.
    call read_inopt(T_wind,'T_wind',db,min=0.,errcount=nerr)
    if (iwind == 3) call read_inopt(temp_exponent,'temp_exponent',db,min=0.,max=5.,errcount=nerr)
 else
    call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)
 endif
 call close_db(db)
 ierr = nerr
 call write_setupfile(filename)

end subroutine read_setupfile

end module setup
