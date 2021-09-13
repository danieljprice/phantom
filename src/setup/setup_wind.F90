!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setup
!
! this module does setup
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - T_wind            : *temperature (K)*
!   - eccentricity      : *eccentricity of the binary system*                      !    = binary_e      : *wide binary eccentricity*
!   - icompanion_star   : *set to 1 for a binary system, 2 for a triple system*
!   - mass_of_particles : *mass resolution (Msun)*
!   - primary_Reff      : *primary star effective radius (au)*
!   - primary_Teff      : *primary star effective temperature (K)*
!   - primary_lum       : *primary star luminosity (Lsun)*
!   - primary_mass      : *primary star mass (Msun)*                               !    = m1            : *first hierarchical level primary mass*
!   - primary_racc      : *primary star accretion radius (au)*                     !    = accr1         : *primary star accretion radius*
!   - secondary_Reff    : *secondary star effective radius (au)*
!   - secondary_Teff    : *secondary star effective temperature)*
!   - secondary_lum     : *secondary star luminosity (Lsun)*
!   - secondary_mass    : *secondary star mass (Msun)*                             !    = m2            : *first hierarchical level secondary mass*
!   - secondary_racc    : *secondary star accretion radius (au)*                   !    = accr2         : *perturber accretion radius*
!   - semi_major_axis   : *semi-major axis of the binary system (au)*              !    = binary_a      : *wide binary semi-major axis*
!   - temp_exponent     : *temperature profile T = R^-p (0 = isothermal)*
!   - wind_gamma        : *adiabatic index (initial if Krome chemistry used)*
!   - accr2a            : *tight binary primary accretion radius*
!   - accr2b            : *tight binary secondary accretion radius*
!   - lum2a             : *tight binary luminosity (Lsun)*
!   - lum2b             : *tight binary luminosity (Lsun)*
!   - Teff2a            : *effective temperature (K)*
!   - Teff2b            : *effective temperature (K)*
!   - Reff2a            : *effective radius (au)*
!   - Reff2b            : *effective radius (au)*
!   - binary2_a         : *tight binary semi-major axis*
!   - binary2_e         : *tight binary eccentricity*
!   - q2                : *tight binary mass ratio*                                !    m2a = m2/(q2+1), m2b = m2*q2/(q2+1)
!   - subst             : *star to substitute, or sky if subst=0*


!   incl,posang_ascnode, arg_peri, omega_corotate, f, verbose not needed but may be included

!
! :Dependencies: eos, infile_utils, inject, io, part, physcon, prompting,
!   setbinary, spherical, units
!

 implicit none
 public :: setpart

 private
 real, public :: wind_gamma = 5./3.
#ifdef ISOTHERMAL
 real, public :: T_wind = 50000.
#else
 real, public :: T_wind = 3000.
#endif
 integer, public :: icompanion_star = 0
 real :: semi_major_axis
 real :: eccentricity = 0.

 real :: primary_Teff = 3000.
 real :: primary_Reff
 real :: primary_lum
 real :: primary_mass
 real :: primary_racc
 real :: secondary_Teff = 0.
 real :: secondary_Reff
 real :: secondary_lum
 real :: secondary_mass
 real :: secondary_racc
 real :: semi_major_axis_au = 4.0
 real :: default_particle_mass = 1.e-11
 real :: primary_lum_lsun = 5315.
 real :: primary_mass_msun = 1.5
 real :: primary_Reff_au = 1.
 real :: primary_racc_au = 1.
 real :: secondary_lum_lsun = 0.
 real :: secondary_mass_msun= 1.0
 real :: secondary_Reff_au = 0.
 real :: secondary_racc_au = 0.1
 real :: temp_exponent = 0.5
 real :: q2
 real :: accr2a
 real :: accr2b
 real :: lum2a = 0.
 real :: lum2b = 0.
 real :: Teff2a = 0.
 real :: Teff2b = 0.
 real :: Reff2a = 0.
 real :: Reff2b = 0.
 real :: binary2_a
 real :: binary2_e
 real :: binary2_a_au = 0.3
 real :: racc2a_au = 0.1
 real :: racc2b_au = 0.1
 real :: binary2_i = 0.
 integer :: subst
 
 
contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only: xyzmh_ptmass, vxyz_ptmass, nptmass, igas, iTeff, iLum, iReff
 use physcon,   only: au, solarm, mass_proton_cgs, kboltz
 use units,     only: umass, set_units,unit_velocity
 use inject,    only: init_inject
 use setbinary, only: set_binary,set_multiple
 use io,        only: master
 use eos,       only: gmw,ieos,isink,qfacdisc
 use spherical, only:set_sphere
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
 integer :: ierr
 logical :: iexist

 call set_units(dist=au,mass=solarm,G=1.)
!
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
 else if (icompanion_star > 1) then
    !-- hierarchical triple
       nptmass  = 0
       print "(/,a)",'----------- Hierarchical triple -----------'
       print "(a,g10.3,a)",'     First hierarchical level primary mass: ', primary_mass!,      trim(mass_unit)
       print "(a,g10.3,a)",'   First hierarchical level secondary mass: ', secondary_mass!,    trim(mass_unit)
       print "(a,g10.3)",  '                    Wide binary mass ratio: ', secondary_mass/primary_mass
       print "(a,g10.3)",  '                   Tight binary mass ratio: ', q2
       print "(a,g10.3)",  '                    Star to be substituted: ', abs(subst)
!        print "(a,g10.3,a)",'                        Accretion Radius 1: ', primary_racc!, trim(dist_unit)
!        print "(a,g10.3,a)",'                       Accretion Radius 2a: ', accr2a!, trim(dist_unit)
!        print "(a,g10.3,a)",'                       Accretion Radius 2b: ', accr2b!, trim(dist_unit)

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
                accretion_radius1=accr2a,accretion_radius2=accr2b, &
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
            
       else if (subst == 11) then     !AGB as secondary, because that is the one from whcih the wind is launched
            call set_multiple(primary_mass/(q2+1),primary_mass*q2/(q2+1),semimajoraxis=binary2_a,eccentricity=binary2_e, &
                accretion_radius1=primary_racc,accretion_radius2=accr2b, &
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
       
       
 else
    nptmass = 1
    xyzmh_ptmass(4,1) = primary_mass
    xyzmh_ptmass(5,1) = primary_racc
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
 endif

 !
 ! for binary wind simulations this mass is IRRELEVANT
 ! since it will be over-written on the first call to init_inject
 !
 massoftype(igas) = default_particle_mass * (solarm / umass)

#ifdef ISOTHERMAL
 gamma = 1.
#else
 gamma = wind_gamma
#endif
 if (gamma < 1.0001) then
    if (temp_exponent > 0.) then
       ieos = 6
       qfacdisc = 0.5*temp_exponent
       isink = 1
       polyk = kboltz*primary_Teff/(mass_proton_cgs * gmw * unit_velocity**2)
    else
       ieos = 1
       polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)
    endif
 else
    polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)
 endif


!
! add low density background medium
!
! sphere of 20^3 particles, density will be determined when mass is set
! call set_sphere('cubic',id,master,primary_racc,20.*primary_racc,primary_racc,hfact,npart,xyzh,rhofunc=rhor)
! npartoftype(igas) = npart

! contains
!  real function rhor(r)
!    real, intent(in) :: r
!    rhor = 1./r**2
!  end function rhor

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
 integer :: ichoice,iwind

 ichoice = 1
 iwind = 2
 call prompt('Type of wind: 1=isoT, 2=adia, 3=T(r)',iwind,1,3)
 if (iwind == 1 .or. iwind == 3) wind_gamma = 1.
 if (iwind == 3) temp_exponent = 0.5
#ifndef ISOTHERMAL
 if (iwind == 1 .or. iwind == 3) then
    call fatal('setup','If you choose options 1 or 3, the code must be compiled with SETUP=isowind')
 endif
#endif

 call prompt('Add binary?',icompanion_star,0,2)
 if (icompanion_star > 1) then  !Hierarchical triple case
    ichoice = 1
    print "(a)",'Star to be substituted by tight binary'
    print "(a)",' 1: primary' ,' 0: companion'
    call prompt('Select star to be substituted',ichoice,0,1)
    select case(ichoice)
    case(1)
        subst   = 11
    case default
        subst   = 12
    end select   
    ichoice = 1
    print "(a)",'Orbital parameters first hierarchical level binary'
    print "(a)",' 1: semi-axis: 15 au, eccentricity: 0',' 0: custom'
    call prompt('select semi-major axis and ecccentricity',ichoice,0,1)
    select case(ichoice)
    case(1)
       semi_major_axis = 15. * au / udist
       eccentricity = 0.
    case default
       semi_major_axis_au = 15.
       eccentricity = 0.
       call prompt('enter semi-major axis in au',semi_major_axis_au,0.,100.)
       call prompt('enter eccentricity',eccentricity,0.)
       semi_major_axis = semi_major_axis_au * au / udist
    end select
    semi_major_axis_au = semi_major_axis * udist / au
    if (subst == 12) then
        print "(a)",'Primary star parameters'
        print "(a)",' 2: Mass: 1.2 Msun, accretion radius: 0.2568 au',&
        ' 1: Mass: 1.5 Msun, accretion radius: 1.2568 au', &
        ' 0: custom'
        call prompt('select mass and radius of primary',ichoice,0,2)
        select case(ichoice)
        case(2)
            primary_mass = 1.2 * (solarm / umass)
            primary_racc = 0.2568 * (au / udist)
        case(1)
            primary_mass = 1.5 * (solarm / umass)
            primary_racc = 1.2568 * (au / udist)
        case default
            primary_mass_msun = 1.5
            primary_racc_au = 1.
            call prompt('enter primary mass',primary_mass_msun,0.,100.)
            call prompt('enter accretion radius in au ',primary_racc_au,0.)
            primary_mass = primary_mass_msun * (solarm / umass)
            primary_racc = primary_racc_au * (au / udist)
        end select
        primary_racc_au   = primary_racc*udist/au
        primary_mass_msun = primary_mass * (umass /solarm)
        ichoice = 1
 
        print "(a)",'Secondary, tight binary parameters'
        print "(a)",' 1: Total Mass: 1.0 Msun',' 0: custom'
        call prompt('select mass',ichoice,0,1)
        select case(ichoice)
        case(1)
            secondary_mass = 1. * (solarm / umass)
        case default
            secondary_mass_msun = 1.
            call prompt('enter total mass tigh binary',secondary_mass_msun,0.,100.)
            secondary_mass = secondary_mass_msun * (solarm / umass)
        end select
        secondary_mass_msun = secondary_mass * (umass /solarm)
        ichoice = 1
        print "(a)",' 1: mass ratio m2b/m2a: 1, accretion radius a: 0.01 au, accretion radius b: 0.01 au',' 0: custom'
        call prompt('select mass ratio and accretion radii of tight binary',ichoice,0,1)
        select case(ichoice)
        case(1)
            q2 = 1.
            accr2a = 0.1 * (au / udist)
            accr2b = 0.1 * (au / udist)
        case default
            q2 = 1.
            racc2a_au = 0.1
            racc2b_au = 0.1
            call prompt('enter tight binary mass ratio',q2,0.)
            call prompt('enter accretion radius a in au ',racc2a_au,0.)
            call prompt('enter accretion radius b in au ',racc2b_au,0.)
            accr2a = racc2a_au * (au / udist)
            accr2b = racc2b_au * (au / udist)
        end select
        racc2a_au = accr2a*udist/au
        racc2b_au = accr2b*udist/au
        ichoice = 1
    else if (subst == 11) then
        print "(a)",'Secondary star parameters'
        print "(a)",' 1: Mass: 1.0 Msun, accretion radius: 0.1 au',' 0: custom'
        call prompt('select mass and radius of secondary',ichoice,0,1)
        select case(ichoice)
        case(1)
            secondary_mass = 1. * (solarm / umass)
            secondary_racc = 0.1 * (au / udist)
        case default
            secondary_mass_msun = 1.
            secondary_racc_au = 0.1
            call prompt('enter secondary mass',secondary_mass_msun,0.,100.)
            call prompt('enter accretion radius in au ',secondary_racc_au,0.)
            secondary_mass = secondary_mass_msun * (solarm / umass)
            secondary_racc = secondary_racc_au * (au / udist)
        end select
        secondary_racc_au   = secondary_racc*udist/au
        secondary_mass_msun = secondary_mass * (umass /solarm)
        ichoice = 1
        
        print "(a)",'Primary (wind-launching) star accretion radius'
        print "(a)",' 2: accretion radius primary: 0.2568 au',&
        ' 1: accretion radius primary: 1.2568 au', &
        ' 0: custom'
        call prompt('select radius of primary',ichoice,0,2)
        select case(ichoice)
        case(2)
            primary_racc = 0.2568 * (au / udist)
        case(1)
            primary_racc = 1.2568 * (au / udist)
        case default
            primary_racc_au = 1.
            call prompt('enter accretion radius in au ',primary_racc_au,0.)
            primary_racc = primary_racc_au * (au / udist)
        end select
        primary_racc_au   = primary_racc*udist/au
        ichoice = 1
        print "(a)",'Tight binary parameters:'
        print "(a)",'Total mass:'
        print "(a)",' 2: Total mass tight binary: 1.2 Msun',&
        ' 1: Total mass tight binary: 1.5 Msun', &
        ' 0: custom'
        call prompt('select total mass tight binary',ichoice,0,2)
        select case(ichoice)
        case(2)
            primary_mass = 1.2 * (solarm / umass)
        case(1)
            primary_mass = 1.5 * (solarm / umass)
        case default
            primary_mass_msun = 1.5
            call prompt('enter primary mass',primary_mass_msun,0.,100.)
            primary_mass = primary_mass_msun * (solarm / umass)
        end select
        primary_mass_msun = primary_mass * (umass /solarm)
        ichoice = 1
        print "(a)",'Mass ratio and accretion radius of secondary:'
        print "(a)",' 1: mass ratio m2b/m2a: 0.3, accretion radius b: 0.01 au',' 0: custom'
        call prompt('select mass ratio and accretion radius of tight binary',ichoice,0,1)
        select case(ichoice)
        case(1)
            q2 = 0.3
            accr2b = 0.1 * (au / udist)
        case default
            q2 = 0.3
            racc2b_au = 0.1
            call prompt('enter tight binary mass ratio',q2,0.)
            call prompt('enter accretion radius b in au ',racc2b_au,0.)
            accr2b = racc2b_au * (au / udist)
        end select
        racc2b_au = accr2b*udist/au
        ichoice = 1
    endif
    
    print "(a)",'Orbital parameters tight binary parameters:'
    print "(a)",' 1: semi-axis: 4 au, eccentricity: 0',' 0: custom'
    call prompt('select tight binary semi-major axis and ecccentricity',ichoice,0,1)
    select case(ichoice)
    case(1)
       binary2_a = 4. * au / udist
       binary2_e = 0.
    case default
       binary2_a_au = 4.
       binary2_e = 0.
       call prompt('enter semi-major axis in au',binary2_a_au,0.,semi_major_axis_au)
       call prompt('enter eccentricity',binary2_e,0.)
       binary2_a = binary2_a_au * au / udist
    end select
    binary2_a_au = binary2_a * udist / au
    ichoice = 1
    print "(a)",'inclination of orbit tight binary w.r.t. outer binary:'
    print "(a)",' 1: inclination: 0 deg',' 0: custom'
    call prompt('select inclination',ichoice,0,1)
    select case(ichoice)
    case(1)
        binary2_i = 0.
    case default
        call prompt('enter inclination',binary2_i,0.,90.)
    end select
 else  !binary or single star case
    if (icompanion_star > 0) then  
        print "(a)",'Primary star parameters'
    else
        print "(a)",'Stellar parameters'
    endif
    print "(a)",' 2: Mass: 1.2 Msun, accretion radius: 0.2568 au',&
        ' 1: Mass: 1.0 Msun, accretion radius: 1.2568 au', &
        ' 0: custom'
    call prompt('select mass and radius of primary',ichoice,0,2)
    select case(ichoice)
    case(2)
        primary_mass = 1.2 * (solarm / umass)
        primary_racc = 0.2568 * (au / udist)
    case(1)
        primary_mass = 1. * (solarm / umass)
        primary_racc = 1.2568 * (au / udist)
    case default
        primary_mass_msun = 1.
        primary_racc_au = 1.
        call prompt('enter primary mass',primary_mass_msun,0.,100.)
        call prompt('enter accretion radius in au ',primary_racc_au,0.)
        primary_mass = primary_mass_msun * (solarm / umass)
        primary_racc = primary_racc_au * (au / udist)
    end select
    primary_racc_au   = primary_racc*udist/au
    primary_mass_msun = primary_mass * (umass /solarm)
    ichoice = 1
    if (icompanion_star > 0) then
        print "(a)",'Secondary star parameters'
        print "(a)",' 1: Mass: 1.0 Msun, accretion radius: 0.1 au',' 0: custom'
        call prompt('select mass and radius of secondary',ichoice,0,1)
        select case(ichoice)
        case(1)
            secondary_mass = 1. * (solarm / umass)
            secondary_racc = 0.1 * (au / udist)
        case default
            secondary_mass_msun = 1.
            secondary_racc_au = 0.1
            call prompt('enter secondary mass',secondary_mass_msun,0.,100.)
            call prompt('enter accretion radius in au ',secondary_racc_au,0.)
            secondary_mass = secondary_mass_msun * (solarm / umass)
            secondary_racc = secondary_racc_au * (au / udist)
        end select
        secondary_racc_au   = secondary_racc*udist/au
        secondary_mass_msun = secondary_mass * (umass /solarm)
        ichoice = 1
        print "(a)",'Orbital parameters'
        print "(a)",' 1: semi-axis: 3.7 au, eccentricity: 0',' 0: custom'
        call prompt('select semi-major axis and ecccentricity',ichoice,0,1)
        select case(ichoice)
        case(1)
            semi_major_axis = 3.7 * au / udist
            eccentricity = 0.
        case default
            semi_major_axis_au = 1.
            eccentricity = 0.
            call prompt('enter semi-major axis in au',semi_major_axis_au,0.,100.)
            call prompt('enter eccentricity',eccentricity,0.)
            semi_major_axis = semi_major_axis_au * au / udist
        end select
        semi_major_axis_au = semi_major_axis * udist / au
    endif
 endif
 
end subroutine setup_interactive

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 use physcon,      only:au,steboltz,solarl,solarm,pi
 use units,        only:utime,unit_energ
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wind setup routine'
 if (primary_Teff == 0. .and. primary_lum_lsun > 0. .and. primary_Reff_au > 0.) &
      primary_Teff = (primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Reff_au*au))**0.25
 if (primary_Reff_au == 0. .and. primary_lum_lsun > 0. .and. primary_Teff > 0.) &
      primary_Reff_au = sqrt(primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Teff**4))/au
 if (primary_Reff_au > 0. .and. primary_lum_lsun == 0. .and. primary_Teff > 0.) &
      primary_lum_lsun = 4.*pi*steboltz*primary_Teff**4*(primary_Reff_au*au)**2/solarl
 primary_lum = primary_lum_lsun * (solarl * utime / unit_energ)
 
 if (icompanion_star > 1) then
    if (secondary_Teff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Reff_au > 0.) &
        secondary_Teff = (secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Reff_au*au))**0.25
    if (secondary_Reff_au == 0. .and. secondary_lum_lsun > 0. .and. secondary_Teff > 0.) &
        secondary_Reff_au = sqrt(secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Teff**4))/au
    if (secondary_Reff_au > 0. .and. secondary_lum_lsun == 0. .and. secondary_Teff > 0.) &
        secondary_lum_lsun = 4.*pi*steboltz*secondary_Teff**4*(secondary_Reff_au*au)**2/solarl
    
    secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
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
        call write_inopt(accr2a,'accr2a','tight binary primary accretion radius',iunit)
        call write_inopt(accr2b,'accr2b','tight binary secondary accretion radius',iunit)
        call write_inopt(lum2a,'lum2a','tight binary primary luminosity (Lsun)',iunit)
        call write_inopt(lum2b,'lum2b','tight binary secondary luminosity (Lsun)',iunit) 
        call write_inopt(Teff2a,'Teff2a','tight binary primary effective temperature (K)',iunit)
        call write_inopt(Teff2b,'Teff2b','tight binary secondary effective temperature (K)',iunit)    
        call write_inopt(Reff2a,'Reff2a','tight binary primary effective radius (au)',iunit)
        call write_inopt(Reff2b,'Reff2b','tight binary secondary effective radius (au)',iunit)      
    else if (subst == 11) then
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
        call write_inopt(accr2b,'accr2b','tight binary secondary accretion radius',iunit)        
        call write_inopt(lum2b,'lum2b','tight binary secondary luminosity (Lsun)',iunit) 
        call write_inopt(Teff2b,'Teff2b','tight binary secondary effective temperature (K)',iunit)    
        call write_inopt(Reff2b,'Reff2b','tight binary secondary effective radius (au)',iunit)      
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
 else  !binary or single star
    call write_inopt(primary_mass_msun,'primary_mass','primary star mass (Msun)',iunit)
    call write_inopt(primary_racc_au,'primary_racc','primary star accretion radius (au)',iunit)
    call write_inopt(primary_lum_lsun,'primary_lum','primary star luminosity (Lsun)',iunit)
    call write_inopt(primary_Teff,'primary_Teff','primary star effective temperature (K)',iunit)
    call write_inopt(primary_Reff_au,'primary_Reff','primary star effective radius (au)',iunit)
    call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system, 2 for a triple system',iunit)
    if (icompanion_star == 1) then
        if (secondary_Teff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Reff_au > 0.) &
            secondary_Teff = (secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Reff_au*au))**0.25
        if (secondary_Reff_au == 0. .and. secondary_lum_lsun > 0. .and. secondary_Teff > 0.) &
            secondary_Reff_au = sqrt(secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Teff**4))/au
        if (secondary_Reff_au > 0. .and. secondary_lum_lsun == 0. .and. secondary_Teff > 0.) &
            secondary_lum_lsun = 4.*pi*steboltz*secondary_Teff**4*(secondary_Reff_au*au)**2/solarl
        secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
        call write_inopt(secondary_mass_msun,'secondary_mass','secondary star mass (Msun)',iunit)
        call write_inopt(secondary_racc_au,'secondary_racc','secondary star accretion radius (au)',iunit)
        call write_inopt(secondary_lum_lsun,'secondary_lum','secondary star luminosity (Lsun)',iunit)
        call write_inopt(secondary_Teff,'secondary_Teff','secondary star effective temperature)',iunit)
        call write_inopt(secondary_Reff_au,'secondary_Reff','secondary star effective radius (au)',iunit)
        call write_inopt(semi_major_axis_au,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
        call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
    endif
 endif
 call write_inopt(default_particle_mass,'mass_of_particles','mass resolution (Msun)',iunit)
 
#ifdef ISOTHERMAL
 wind_gamma = 1.
#else
 call write_inopt(wind_gamma,'wind_gamma','adiabatic index (initial if Krome chemistry used)',iunit)
#endif
 if ( wind_gamma < 1.0001) then
    call write_inopt(temp_exponent,'temp_exponent','temperature profile T = R^-p (0 = isothermal)',iunit)
    if (abs(temp_exponent) < tiny(0.)) call write_inopt(T_wind,'T_wind','temperature (K)',iunit)
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
 call read_inopt(icompanion_star,'icompanion_star',db,min=0,errcount=nerr)
 if (icompanion_star == 1) then
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(secondary_lum_lsun,'secondary_lum',db,min=0.,max=1000.,errcount=nerr)
    secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
    call read_inopt(secondary_Teff,'secondary_Teff',db,min=0.,max=1000.,errcount=nerr)
    call read_inopt(secondary_Reff_au,'secondary_Reff',db,min=0.,max=1000.,errcount=nerr)
    secondary_Reff = secondary_Reff_au * au / udist
    call read_inopt(secondary_racc_au,'secondary_racc',db,min=0.,max=1000.,errcount=nerr)
    secondary_racc = secondary_racc_au * au / udist
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
 else if (icompanion_star > 1) then
    !-- hierarchical triple
    call read_inopt(subst,'subst',db,errcount=nerr)
    if (subst == 11) then
        call read_inopt(secondary_lum_lsun,'secondary_lum',db,min=0.,max=1000.,errcount=nerr)
        secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
        call read_inopt(secondary_Teff,'secondary_Teff',db,min=0.,max=1000.,errcount=nerr)
        call read_inopt(secondary_Reff_au,'secondary_Reff',db,min=0.,max=1000.,errcount=nerr)
        secondary_Reff = secondary_Reff_au * au / udist
        call read_inopt(secondary_racc_au,'secondary_racc',db,min=0.,max=1000.,errcount=nerr)
        secondary_racc = secondary_racc_au * au / udist
    else if (subst == 12) then
        call read_inopt(accr2a,'accr2a',db,errcount=nerr)
        call read_inopt(lum2a,'lum2a',db,errcount=nerr)
        call read_inopt(Teff2a,'Teff2a',db,errcount=nerr)
        call read_inopt(Reff2a,'Reff2a',db,errcount=nerr)
    endif
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
        secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
    !-- masses
    call read_inopt(q2,'q2',db,min=0.,max=1.,errcount=nerr)
    !-- tight parameters
    call read_inopt(binary2_a,'binary2_a',db,errcount=nerr)
    call read_inopt(binary2_e,'binary2_e',db,errcount=nerr)
    !-- accretion radii,...
    call read_inopt(accr2b,'accr2b',db,errcount=nerr)
    call read_inopt(lum2b,'lum2b',db,errcount=nerr)
    call read_inopt(Teff2b,'Teff2b',db,errcount=nerr)
    call read_inopt(Reff2b,'Reff2b',db,errcount=nerr)
    call read_inopt(binary2_i,'inclination',db,errcount=nerr)    

 endif
 call read_inopt(default_particle_mass,'mass_of_particles',db,min=0.,errcount=nerr)
#ifdef ISOTHERMAL
 wind_gamma = 1.
#else
 call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)
#endif
 if ( wind_gamma < 1.0001) then
    call read_inopt(temp_exponent,'temp_exponent',db,min=0.,max=5.,errcount=nerr)
    if (abs(temp_exponent) < tiny(0.)) call read_inopt(T_wind,'T_wind',db,min=0.,errcount=nerr)
 endif
 call close_db(db)
 if (primary_Teff == 0. .and. primary_lum_lsun > 0. .and. primary_Reff > 0.) then
    ichange = ichange+1
    primary_Teff = (primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Reff*udist))**0.25
 endif
 if (primary_Reff == 0. .and. primary_lum_lsun > 0. .and. primary_Teff > 0.) then
    ichange = ichange+1
    primary_Reff = sqrt(primary_lum_lsun*solarl/(4.*pi*steboltz*primary_Teff**4))/udist
    primary_Reff_au = primary_Reff * udist / au
 endif
 if (primary_Reff > 0.  .and. primary_lum_lsun == 0. .and. primary_Teff > 0.) then
    ichange = ichange+1
    primary_lum_lsun = 4.*pi*steboltz*primary_Teff**4*(primary_Reff_au*au)**2/solarl
 endif
 if (icompanion_star > 0) then
    if (secondary_Teff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Reff > 0.) then
       ichange = ichange+1
       secondary_Teff = (secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Reff*udist))**0.25
    endif
    if (secondary_Reff == 0. .and. secondary_lum_lsun > 0. .and. secondary_Teff > 0.) then
       ichange = ichange+1
       secondary_Reff = sqrt(secondary_lum_lsun*solarl/(4.*pi*steboltz*secondary_Teff**4))/udist
       secondary_Reff_au = secondary_Reff * udist / au
    endif
 endif
 ierr = nerr
 call write_setupfile(filename)
 
end subroutine read_setupfile
end module setup
