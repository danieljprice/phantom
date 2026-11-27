!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! initial conditions for binary wind accretion / AGB star wind injection,
! from a single star, a binary or a hierarchical system of stars
!
! :References:
!   Siess et al. 2022, A&A, 667, 75
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - T_wind              : *wind temperature (K)*
!   - iresolution         : *particle mass from 0=mass,1=particles_per_shell,2=sonic point*
!   - iwind               : *isothermal wind: 2=constant temperature, 3=T(r) profile*
!   - mass_of_particles   : *particle mass (Msun)*
!   - nstar               : *number of stars (1=single, 2=binary, >=3 hierarchical system)*
!   - nzones_sonic_point  : *particle spacings from injection radius to sonic point*
!   - particles_per_shell : *number of particles per shell of the first star with a wind*
!   - temp_exponent       : *temperature profile T(r) = T_wind*(r/Reff)^(-temp_exponent)*
!   - wind_gamma          : *adiabatic index (initial if Krome chemistry used)*
!
! :Dependencies: dim, eos, infile_utils, inject, io, kernel, part, physcon,
!   prompting, ptmass_radiation, sethier_utils, sethierarchical, setorbit,
!   units, wind
!
 use dim,           only:isothermal
 use setorbit,      only:orbit_t
 use sethier_utils, only:max_hier_levels
 implicit none
 public :: setpart

 private
 real, public  :: wind_gamma
 real, public  :: T_wind
 real :: temp_exponent
 integer :: nstar,iwind
 type(orbit_t) :: orbit
 real :: default_particle_mass
 integer :: iresolution,particles_per_shell,nzones_sonic_point

 ! types of star used to set the defaults
 integer, parameter :: itype_agb = 1, itype_transsonic = 2, itype_ms = 3

 ! which of L, Teff and Reff is computed from the other two
 integer, parameter :: iderived_none = 0, &
                       iderived_lum  = 1, &
                       iderived_Teff = 2, &
                       iderived_Reff = 3

 ! properties of each star, in the units used in the .setup file
 type wind_star_t
    real :: mass_msun,racc_au,lum_lsun,Teff,Reff_au
    character(len=20) :: rinject  ! wind injection radius, with units (e.g. 2.5 au)
    real :: mdot_msun_yr,vwind_km_s,wind_temp,alpha
    real :: veq_km_s,spin(3)
    integer :: iderived
 end type wind_star_t

 integer, parameter :: maxstars = max_hier_levels
 type(wind_star_t) :: star(maxstars)

 ! label used for star i in the .setup file
 character(len=*), parameter :: star_label(maxstars) = &
    (/'primary  ','secondary','star3    ','star4    ','star5    ', &
      'star6    ','star7    ','star8    ','star9    ','star10   '/)

contains
!----------------------------------------------------------------
!+
!  default parameter choices
!+
!----------------------------------------------------------------
subroutine set_default_parameters_wind()
 use setorbit, only:set_defaults_orbit
 integer :: i,ipair

 wind_gamma    = 5./3.
 if (isothermal) then
    T_wind        = 100000.
    temp_exponent = 0.5
    iwind         = 2
 else
    T_wind = 3000.
    iwind  = 1
 endif
 nstar = 1
 default_particle_mass = 1.e-10
 iresolution           = 0
 particles_per_shell   = 1000
 nzones_sonic_point    = 8
 !
 ! an AGB star with a wind, and main-sequence companions without a wind
 !
 call set_star_type(star(1),itype_agb,ipair)
 do i=2,maxstars
    call set_star_type(star(i),itype_ms,ipair)
 enddo
 !
 ! binary orbit, with the periapsis along the x-axis and starting at apastron
 !
 call set_defaults_orbit(orbit)
 orbit%elems%a = '4.0'
 orbit%e = 0.
 orbit%w = 0.
 orbit%f = 180.

end subroutine set_default_parameters_wind

!----------------------------------------------------------------
!+
!  default properties of a star, with no wind and no spin
!+
!----------------------------------------------------------------
subroutine set_default_star(s,mass_msun,racc_au,lum_lsun,Reff_au)
 type(wind_star_t), intent(out) :: s
 real,              intent(in)  :: mass_msun,racc_au,lum_lsun,Reff_au

 s%mass_msun    = mass_msun
 s%racc_au      = racc_au
 s%lum_lsun     = lum_lsun
 s%Teff         = 0.
 s%Reff_au      = Reff_au
 s%rinject      = '0'
 s%mdot_msun_yr = 0.
 s%vwind_km_s   = 0.
 s%wind_temp    = 0.
 s%alpha        = 0.
 s%veq_km_s     = 0.
 s%spin         = 0.
 s%spin(3)      = 1.  ! spin along z-axis
 s%iderived     = iderived_none

end subroutine set_default_star

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,             only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas,ivwind,imloss,iLum,iReff
 use physcon,          only:au,solarm,mass_proton_cgs,kboltz,solarl
 use units,            only:umass,set_units,unit_velocity,utime,unit_energ,udist,get_nearest_unit
 use inject,           only:set_default_options_inject,print_wind_resolution
 use wind,             only:get_wind_particle_mass
 use setorbit,         only:set_orbits
 use io,               only:master,fatal
 use eos,              only:gmw,ieos,isink,qfacdisc,init_eos
 use infile_utils,     only:get_options,infile_exists
 use kernel,           only:hfact_default
 use ptmass_radiation, only:isink_radiation
 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: vxyzu(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma,hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 integer :: ierr,i,wind_type,iref
 real :: m(maxstars),racc(maxstars)

 hfact = hfact_default
 call set_units(dist=au,mass=solarm,G=1.d0)
 call set_default_parameters_wind()
 if (.not. infile_exists(fileprefix)) call set_default_options_inject()

!--general parameters
!
 time = 0.
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Wind setup'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
      read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
!
!--the wind alpha parameters are ignored unless radiation pressure from sinks is switched on
!
 if (any(star(1:nstar)%alpha > 0.) .and. isink_radiation == 0) then
    isink_radiation = 1
    if (id==master) print "(a)",' wind alpha > 0: setting isink_radiation = 1 in the .in file'
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
!
!--place the stars: star i is sink particle i
!
 do i=1,nstar
    m(i)    = star(i)%mass_msun * (solarm / umass)
    racc(i) = star(i)%racc_au * (au / udist)
 enddo
 nptmass = 0
 call set_orbits(nstar,m,racc,orbit,fileprefix,xyzmh_ptmass,vxyz_ptmass,nptmass,id==master,ierr)
 if (ierr /= 0) call fatal('setup','error placing stars on their orbits')
 if (nstar == 2 .and. id==master) print "(2x,a)",'binary period (physical units) : '//&
    trim(get_nearest_unit(orbit%period,'time'))
 do i=1,nstar
    call set_sink_properties(star(i),xyzmh_ptmass(:,i))
 enddo
 if (nstar >= 3 .and. id==master) then
    print *,'Sink particles summary'
    print *,'  #    mass       racc      lum         Reff'
    do i=1,nptmass
       print '(i4,2(2x,f9.5),2(2x,es10.3))',i,xyzmh_ptmass(4:5,i),xyzmh_ptmass(iLum,i)/(solarl*utime/unit_energ),&
             xyzmh_ptmass(iReff,i)*udist/au
    enddo
    print *,''
 endif

 if (isothermal) then
    gamma = 1.
    if (iwind == 3) then
       ieos = 6
       qfacdisc = 0.5*temp_exponent
       isink = 1
       T_wind = star(1)%Teff
    else
       isink = 1
       ieos = 1
    endif
 else
    T_wind = 0.
    gamma = wind_gamma
 endif
 polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)
 ! the 1D wind solution used to set the particle mass needs the equation of state
 call init_eos(ieos,ierr)
 if (ierr /= 0) call fatal('setup','error initialising equation of state')

 !
 ! set the particle mass, either directly or from the desired resolution
 ! of the wind from the first star with a wind (iref). This is the only
 ! place the particle mass is set.
 !
 iref = 0
 do i=nstar,1,-1
    if (xyzmh_ptmass(imloss,i) > 0.) iref = i
 enddo
 if (iref == 0) call fatal('setup','no winds in the simulation: set mdot > 0 for at least one star')
 ! trans-sonic wind (wind_type = 1) if resolving the sonic point, or if no wind
 ! velocity is given so that it comes from the 1D wind solution
 wind_type = 0
 if (iref > 0) then
    if (.not.isothermal .and. (iresolution == 2 .or. xyzmh_ptmass(ivwind,iref) <= 0.)) wind_type = 1
 endif
 if (wind_type == 1) call set_default_options_inject(1)
 select case(iresolution)
 case(1)
    massoftype(igas) = get_wind_particle_mass(xyzmh_ptmass(:,iref),wind_type,particles_per_shell,0)
 case(2)
    massoftype(igas) = get_wind_particle_mass(xyzmh_ptmass(:,iref),wind_type,0,nzones_sonic_point)
 case default
    massoftype(igas) = default_particle_mass * (solarm / umass)
 end select
 call print_wind_resolution(xyzmh_ptmass,nptmass,massoftype(igas),gamma)

end subroutine setpart

!----------------------------------------------------------------
!+
!  set the stellar, wind and spin properties of a sink particle
!  from the properties of a star
!+
!----------------------------------------------------------------
subroutine set_sink_properties(s,xyzmh_ptmassi)
 use part,    only:iTeff,iLum,iReff,ispinx,ispiny,ispinz,ivwind,imloss,iTwind,iwalpha
 use physcon, only:solarm,km,years
 use units,   only:umass,utime,unit_velocity,check_and_convert
 type(wind_star_t), intent(inout) :: s
 real,              intent(inout) :: xyzmh_ptmassi(:)
 real :: lum,Reff,veq,rinject
 integer :: nerr

 call get_lum_and_Reff(s%lum_lsun,s%Reff_au,s%Teff,lum,Reff,s%iderived)
 ! the wind is injected at the stellar photosphere unless an injection radius is given
 nerr = 0
 call check_and_convert(s%rinject,'Rinj','length',rinject,nerr)
 if (rinject <= 0.) rinject = Reff
 xyzmh_ptmassi(iTeff)   = s%Teff
 xyzmh_ptmassi(iReff)   = rinject
 xyzmh_ptmassi(iLum)    = lum
 xyzmh_ptmassi(imloss)  = s%mdot_msun_yr * (solarm / umass) * (utime / years)
 xyzmh_ptmassi(ivwind)  = s%vwind_km_s * (km / unit_velocity)
 xyzmh_ptmassi(iTwind)  = s%wind_temp
 xyzmh_ptmassi(iwalpha) = s%alpha
 veq = s%veq_km_s * (km / unit_velocity)
 xyzmh_ptmassi(ispinx)  = Reff**2*s%spin(1)*veq
 xyzmh_ptmassi(ispiny)  = Reff**2*s%spin(2)*veq
 xyzmh_ptmassi(ispinz)  = Reff**2*s%spin(3)*veq

end subroutine set_sink_properties

!----------------------------------------------------------------
!+
!  set up interactively: the number of stars, then the properties
!  of each star, with defaults that depend on the choice of mass
!  and radius (e.g. trans-sonic), then the orbit and resolution
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use io,        only:fatal
 integer :: i

 nstar = 1
 call prompt('Number of stars (1=single, 2=binary, >=3 hierarchical system)',nstar,1,maxstars)
 if (isothermal) then
    call prompt('Type of wind:  2=isoT, 3=T(r)',iwind,2,3)
 else
    ! the other types need an isothermal build, so this tells the user how to get them
    iwind = 1
    call prompt('Type of wind:  1=adia, 2=isoT, 3=T(r)',iwind,1,3)
    if (iwind == 2 .or. iwind == 3) then
       call fatal('setup','If you choose options 2 or 3, the code must be compiled with SETUP=isowind')
    endif
 endif

 do i=1,nstar
    call get_star(star(i),trim(star_label(i)),i==1)
 enddo
 if (iwind == 3) T_wind = star(1)%Teff

 call get_orbits()
 call get_resolution()

end subroutine setup_interactive

!----------------------------------------------------------------
!+
!  prompt for the properties of a star. The type of star sets the
!  defaults for the values asked for afterwards
!+
!----------------------------------------------------------------
subroutine get_star(s,label,primary)
 use prompting, only:prompt
 type(wind_star_t), intent(inout) :: s
 character(len=*),  intent(in)    :: label
 logical,           intent(in)    :: primary
 integer :: itype,ipair

 print "(/,a)",label//' star parameters'
 itype = itype_agb
 if (.not.primary) itype = itype_ms
 print "(a)",' 1: AGB star (1.2 Msun, racc=0.2568 au, Teff=3000 K, Reff=1 au, Mdot=1e-7 Msun/yr)',&
             ' 2: hot star with a trans-sonic wind (1.2 Msun, racc=1 au, L=2e4 Lsun, Reff=2 au, Mdot=1e-5 Msun/yr)',&
             ' 3: main-sequence star without a wind (1 Msun, racc=0.1 au, L=1 Lsun, Reff=1 Rsun)'
 call prompt('select type of star',itype,1,3)
 call set_star_type(s,itype,ipair)
 call prompt('enter mass (Msun)',s%mass_msun,0.)
 call prompt('enter accretion radius (au)',s%racc_au,0.)
 !
 ! two of L, Teff and Reff: the third is derived from them
 !
 print "(a)",' 1: Teff and Reff',' 2: L and Teff',' 3: L and Reff'
 call prompt('select which two of L, Teff and Reff to set',ipair,1,3)
 call get_lum_Teff_Reff(s,ipair)
 !
 ! wind and rotation
 !
 call prompt('enter wind mass loss rate (Msun/yr, 0=no wind)',s%mdot_msun_yr,0.,1.)
 if (s%mdot_msun_yr > 0.) then
    call prompt('enter wind velocity (km/s, 0=trans-sonic)',s%vwind_km_s,0.,10000.)
    call prompt('enter wind temperature (K, 0=use Teff)',s%wind_temp,0.,1.e7)
    call prompt('enter wind acceleration parameter alpha (0.0=no accel,1.0=free wind)',s%alpha,0.,10.)
    call get_sink_spin(s)
 endif

end subroutine get_star

!----------------------------------------------------------------
!+
!  prompt for the two of L, Teff and Reff chosen by ipair
!  (1 = Teff and Reff, 2 = L and Teff, 3 = L and Reff), with the
!  default values derived from those of the star
!+
!----------------------------------------------------------------
subroutine get_lum_Teff_Reff(s,ipair)
 use prompting, only:prompt
 type(wind_star_t), intent(inout) :: s
 integer,           intent(in)    :: ipair
 real :: lum,Reff

 ! fill in the third value so that every default is consistent
 call get_lum_and_Reff(s%lum_lsun,s%Reff_au,s%Teff,lum,Reff,s%iderived)
 s%iderived = iderived_none
 if (ipair /= 1) then
    call prompt('enter luminosity (Lsun)',s%lum_lsun,0.)
 else
    s%lum_lsun = 0.
 endif
 if (ipair /= 3) then
    call prompt('enter effective temperature (K)',s%Teff,0.)
 else
    s%Teff = 0.
 endif
 if (ipair /= 2) then
    call prompt('enter effective radius (au)',s%Reff_au,0.)
 else
    s%Reff_au = 0.
 endif

end subroutine get_lum_Teff_Reff

!----------------------------------------------------------------
!+
!  prompt for the rotation of a star: equatorial velocity and,
!  if it rotates, the direction of the spin axis
!+
!----------------------------------------------------------------
subroutine get_sink_spin(s)
 use prompting, only:prompt
 type(wind_star_t), intent(inout) :: s
 real :: spin

 call prompt('enter rotation velocity of star at equator (km/s)',s%veq_km_s,0.,1000.)
 if (s%veq_km_s > 0.) then
    call prompt('enter x-component of spin axis',s%spin(1),-1.,1.)
    call prompt('enter y-component of spin axis',s%spin(2),-1.,1.)
    call prompt('enter z-component of spin axis',s%spin(3),-1.,1.)
    ! renormalize vector, using the z-axis if a zero vector is given
    spin = sqrt(s%spin(1)**2+s%spin(2)**2+s%spin(3)**2)
    if (spin > 0.) then
       s%spin = s%spin/spin
    else
       print "(a)",' WARNING: zero spin axis given, using the z-axis'
       s%spin = (/0.,0.,1./)
    endif
 endif

end subroutine get_sink_spin

!----------------------------------------------------------------
!+
!  set the properties of a star of a given type: an AGB star, a hot
!  star with a trans-sonic wind or a main-sequence star without a wind.
!  ipair gives which two of L, Teff and Reff are set
!  (1 = Teff and Reff, 2 = L and Teff, 3 = L and Reff)
!+
!----------------------------------------------------------------
subroutine set_star_type(s,itype,ipair)
 type(wind_star_t), intent(out) :: s
 integer,           intent(in)  :: itype
 integer,           intent(out) :: ipair

 select case(itype)
 case(itype_transsonic)
    call set_default_star(s,mass_msun=1.2,racc_au=1.,lum_lsun=2.e4,Reff_au=2.)
    call set_star_wind(s,mdot_msun_yr=1.e-5,vwind_km_s=0.,wind_temp=5.e4)
    wind_gamma = 1.4
    ipair = 3
 case(itype_ms)
    call set_default_star(s,mass_msun=1.0,racc_au=0.1,lum_lsun=1.,Reff_au=0.00465)
    ipair = 3
 case default
    call set_default_star(s,mass_msun=1.2,racc_au=0.2568,lum_lsun=0.,Reff_au=1.)
    s%Teff = 3000.
    call set_star_wind(s,mdot_msun_yr=1.e-7,vwind_km_s=10.,wind_temp=3000.)
    ipair = 1
 end select

end subroutine set_star_type

!----------------------------------------------------------------
!+
!  set the wind of a star (wind_temp = 0 means use Teff)
!+
!----------------------------------------------------------------
subroutine set_star_wind(s,mdot_msun_yr,vwind_km_s,wind_temp)
 type(wind_star_t), intent(inout) :: s
 real,              intent(in)    :: mdot_msun_yr,vwind_km_s,wind_temp

 s%mdot_msun_yr = mdot_msun_yr
 s%vwind_km_s   = vwind_km_s
 s%wind_temp    = wind_temp

end subroutine set_star_wind

!----------------------------------------------------------------
!+
!  prompt for the orbit of a binary; a hierarchical system gets
!  default orbits, to be edited in the .setup file
!+
!----------------------------------------------------------------
subroutine get_orbits()
 use prompting,       only:prompt
 use sethierarchical, only:generate_hierarchy_string,set_hierarchical_default_options,hs
 integer :: i
 real :: a_au

 if (nstar == 2) then
    print "(/,a)",'Orbital parameters'
    a_au    = 3.7
    orbit%e = 0.
    orbit%f = 180.
    call prompt('enter semi-major axis in au',a_au,0.)
    call prompt('enter eccentricity',orbit%e,0.)
    call prompt('enter true anomaly (180=apastron)',orbit%f,0.,360.)
    write(orbit%elems%a,"(g13.6)") a_au
    orbit%elems%a = adjustl(orbit%elems%a)
 elseif (nstar >= 3) then
    call generate_hierarchy_string(nstar)
    call set_hierarchical_default_options()
    ! each orbit is 8 times wider than the one inside it, the innermost being 3.7 au
    do i=1,hs%labels%hl_num
       hs%levels(i)%a = 3.7*8.**(nstar-1-len_trim(hs%labels%hl(i)))
    enddo
    print "(/,a)",' using default orbits for the hierarchical system: edit the .setup file to change them'
 endif

end subroutine get_orbits

!--------------------------------------------------------
!+
! choose how to set the particle mass: directly, from the
! number of particles per shell, or by resolving the
! sonic point (single trans-sonic star only)
!+
!--------------------------------------------------------
subroutine get_resolution()
 use prompting, only:prompt
 integer :: imax

 imax = 1
 if (sonic_option_allowed()) imax = 2
 print "(/,a)",'Resolution'
 print "(a)",' 0: set particle mass',' 1: set number of particles per shell (first star with a wind)'
 if (imax == 2) print "(a)",' 2: set number of particle spacings between injection radius and sonic point'
 call prompt('select how to set the particle mass',iresolution,0,imax)
 select case(iresolution)
 case(1)
    call prompt('enter number of particles per shell',particles_per_shell,1)
 case(2)
    call prompt('enter number of particle spacings to the sonic point',nzones_sonic_point,1)
 case default
    call prompt('enter particle mass in Msun',default_particle_mass,0.)
 end select

end subroutine get_resolution

!--------------------------------------------------------
!+
! resolving the sonic point is only possible for a
! single, non-isothermal (trans-sonic) wind
!+
!--------------------------------------------------------
logical function sonic_option_allowed()

 sonic_option_allowed = (nstar == 1 .and. .not.isothermal)

end function sonic_option_allowed

!----------------------------------------------------------------
!+
!  get luminosity and effective radius in code units from two of
!  L, Teff and Reff in physical units, computing the third from
!  L = 4 pi sigma Reff^2 Teff^4. The one that was zero is the one
!  computed (stored in iderived so that repeated calls are consistent).
!  If all three are given, L and Reff are kept and Teff is recomputed
!+
!----------------------------------------------------------------
subroutine get_lum_and_Reff(lum_lsun,reff_au,Teff,lum,Reff,ider)
 use physcon,          only:au,steboltz,solarl,pi
 use units,            only:udist,unit_luminosity
 use io,               only:warning
 use ptmass_radiation, only:get_Teff
 real,    intent(inout) :: lum_lsun,reff_au,Teff
 real,    intent(out)   :: lum,Reff
 integer, intent(inout) :: ider
 real :: Teff_LR

 if (ider == iderived_none) then
    if (lum_lsun <= 0.) then
       ider = iderived_lum
    elseif (reff_au <= 0.) then
       ider = iderived_Reff
    else
       ider = iderived_Teff
       if (Teff > 0.) then
          Teff_LR = get_Teff(lum_lsun,reff_au)
          if (abs(Teff_LR-Teff) > 1.e-3*Teff) call warning('setup',&
             'L /= 4*pi*sigma*R^2*T^4: Teff recomputed from L and Reff',var='Teff',val=Teff_LR)
       endif
    endif
 endif

 select case(ider)
 case(iderived_lum)
    lum_lsun = 4.*pi*steboltz*Teff**4*(reff_au*au)**2/solarl
 case(iderived_Reff)
    if (Teff > 0.) reff_au = sqrt(lum_lsun*solarl/(4.*pi*steboltz*Teff**4))/au
 case default
    if (reff_au > 0.) Teff = get_Teff(lum_lsun,reff_au)
 end select

 lum  = lum_lsun*(solarl/unit_luminosity)
 Reff = Reff_au*(au/udist)

end subroutine get_lum_and_Reff

!----------------------------------------------------------------
!+
!  write luminosity, effective temperature and radius of a star
!  to the .setup file, omitting the one computed from the other two
!+
!----------------------------------------------------------------
subroutine write_lum_Teff_Reff(label,lum_lsun,Teff,reff_au,ider,iunit)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: label
 real,             intent(in) :: lum_lsun,Teff,reff_au
 integer,          intent(in) :: ider,iunit

 if (ider /= iderived_lum)  call write_inopt(lum_lsun,label//'_lum',label//' star luminosity (Lsun)',iunit)
 if (ider /= iderived_Teff) call write_inopt(Teff,label//'_Teff',label//' star effective temperature (K)',iunit)
 if (ider /= iderived_Reff) call write_inopt(reff_au,label//'_Reff',label//' star effective radius (au)',iunit)

end subroutine write_lum_Teff_Reff

!----------------------------------------------------------------
!+
!  write the properties of one star to the .setup file
!+
!----------------------------------------------------------------
subroutine write_star(s,label,iunit)
 use infile_utils, only:write_inopt
 type(wind_star_t), intent(inout) :: s
 character(len=*),  intent(in)    :: label
 integer,           intent(in)    :: iunit
 real :: lum,Reff

 ! make L, Teff and Reff consistent, so the derived one is not written
 call get_lum_and_Reff(s%lum_lsun,s%Reff_au,s%Teff,lum,Reff,s%iderived)

 write(iunit,"(/,a)") '# '//label//' star'
 call write_inopt(s%mass_msun,label//'_mass',label//' star mass (Msun)',iunit)
 call write_inopt(s%racc_au,label//'_racc',label//' star accretion radius (au)',iunit)
 call write_lum_Teff_Reff(label,s%lum_lsun,s%Teff,s%Reff_au,s%iderived,iunit)

 call write_inopt(s%mdot_msun_yr,label//'_mdot',label//' wind mass loss rate (in Msun/yr, 0=no wind)',iunit)
 if (s%mdot_msun_yr > 0.) then
    call write_inopt(s%rinject,label//'_Rinj','wind injection radius (e.g. 2.5 au, 0=Reff)',iunit)
    call write_inopt(s%vwind_km_s,label//'_vwind',label//' wind velocity (in km/s)',iunit)
    call write_inopt(s%wind_temp,label//'_wind_temp',label//' wind temperature (K, 0 = use Teff)',iunit)
    call write_inopt(s%alpha,label//'_alpha',label//' alpha parameter',iunit)
    call write_inopt(s%veq_km_s,label//'_veq',label//' equatorial velocity (in km/s)',iunit)
    if (s%veq_km_s /= 0) then
       call write_inopt(s%spin(1),label//'_spinx','x-component of spin direction',iunit)
       call write_inopt(s%spin(2),label//'_spiny','y-component of spin direction',iunit)
       call write_inopt(s%spin(3),label//'_spinz','z-component of spin direction',iunit)
    endif
 endif

end subroutine write_star

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils,    only:write_inopt
 use setorbit,        only:write_options_orbits
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20
 integer :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wind setup routine'
 write(iunit,"(a)") '# set any two of lum, Teff and Reff (third one is derived)'

 write(iunit,"(/,a)") '# stars'
 call write_inopt(nstar,'nstar','number of stars (1=single, 2=binary, >=3 hierarchical system)',iunit)
 do i=1,nstar
    call write_star(star(i),trim(star_label(i)),iunit)
 enddo

 call write_options_orbits(nstar,orbit,iunit,star_label)

 write(iunit,"(/,a)") '# resolution'
 call write_inopt(iresolution,'iresolution','particle mass from 0=mass,1=particles_per_shell,2=sonic point',iunit)
 select case(iresolution)
 case(1)
    call write_inopt(particles_per_shell,'particles_per_shell','number of particles per shell of the first star with a wind',iunit)
 case(2)
    call write_inopt(nzones_sonic_point,'nzones_sonic_point','particle spacings from injection radius to sonic point',iunit)
 case default
    call write_inopt(default_particle_mass,'mass_of_particles','particle mass (Msun)',iunit)
 end select

 write(iunit,"(/,a)") '# equation of state'
 if (isothermal) then
    wind_gamma = 1.
    call write_inopt(iwind,'iwind','isothermal wind: 2=constant temperature, 3=T(r) profile',iunit)
    if (iwind == 3) then
       call write_inopt(star(1)%Teff,'T_wind','wind temperature at injection radius (K)',iunit)
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
!  read the properties of one star from the .setup file
!+
!----------------------------------------------------------------
subroutine read_star(s,label,db,nerr)
 use infile_utils, only:inopts,read_inopt
 use units,        only:check_and_convert
 type(wind_star_t),         intent(inout) :: s
 character(len=*),          intent(in)    :: label
 type(inopts), allocatable, intent(inout) :: db(:)
 integer,                   intent(inout) :: nerr
 real    :: rinject

 s%iderived = iderived_none
 call read_inopt(s%mass_msun,label//'_mass',db,min=0.,max=1000.,errcount=nerr)
 call read_inopt(s%lum_lsun,label//'_lum',db,default=0.,min=0.,max=1e7,errcount=nerr)
 call read_inopt(s%Teff,label//'_Teff',db,default=0.,min=0.,errcount=nerr)
 call read_inopt(s%Reff_au,label//'_Reff',db,default=0.,min=0.,errcount=nerr)
 call read_inopt(s%racc_au,label//'_racc',db,min=0.,errcount=nerr)
 if (s%racc_au < tiny(0.)) then
    print "(a)",' ERROR: '//label//' accretion radius not defined'
    nerr = nerr+1
 endif
 call read_inopt(s%mdot_msun_yr,label//'_mdot',db,min=0.,max=1.,errcount=nerr)
 if (s%mdot_msun_yr > 0.) then
    call read_inopt(s%rinject,label//'_Rinj',db,default='0',errcount=nerr)
    call check_and_convert(s%rinject,label//'_Rinj','length',rinject,nerr)
    call read_inopt(s%vwind_km_s,label//'_vwind',db,min=0.,max=1.e4,errcount=nerr)
    call read_inopt(s%wind_temp,label//'_wind_temp',db,min=0.,max=1.e8,errcount=nerr)
    call read_inopt(s%alpha,label//'_alpha',db,min=0.,max=10.,errcount=nerr)
    call read_inopt(s%veq_km_s,label//'_veq',db,default=0.,min=0.,max=1000.,errcount=nerr)
    if (s%veq_km_s /= 0) then
       call read_inopt(s%spin(1),label//'_spinx',db,min=-1.,max=1.,errcount=nerr)
       call read_inopt(s%spin(2),label//'_spiny',db,min=-1.,max=1.,errcount=nerr)
       call read_inopt(s%spin(3),label//'_spinz',db,min=-1.,max=1.,errcount=nerr)
    endif
 endif

end subroutine read_star

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db,int_to_string
 use physcon,      only:solarm
 use units,        only:umass
 use io,           only:fatal
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr,i,icompanion_star,ierr_mdot
 real :: m(maxstars),mdot

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(nstar,'nstar',db,default=0,min=0,max=maxstars,errcount=nerr)
 if (nstar == 0) then
    ! setup files from before nstar was introduced give icompanion_star
    ! (hierarchical triples set up with icompanion_star = 2 are no longer supported:
    ! stop here so that the file is not overwritten with a binary setup)
    call read_inopt(icompanion_star,'icompanion_star',db,min=0,max=2,errcount=nerr)
    if (icompanion_star == 2) call fatal('setup','icompanion_star = 2 is no longer supported: '//&
       'convert a copy of '//trim(filename)//' to nstar = 3 with the hierarchical options')
    nstar = icompanion_star + 1
    print "(/,a,/)",' WARNING: icompanion_star is deprecated, using nstar = '//trim(int_to_string(nstar))//' instead'
 endif
 do i=1,nstar
    call read_star(star(i),trim(star_label(i)),db,nerr)
    m(i) = star(i)%mass_msun * (solarm / umass)
 enddo
 call read_orbit_options(db,m,nerr)

 ! resolution: setup files without iresolution keep the particle mass they specify
 call read_inopt(iresolution,'iresolution',db,min=0,max=2,default=0,errcount=nerr)
 select case(iresolution)
 case(1)
    call read_inopt(particles_per_shell,'particles_per_shell',db,min=1,errcount=nerr)
 case(2)
    call read_inopt(nzones_sonic_point,'nzones_sonic_point',db,min=1,errcount=nerr)
    if (.not.sonic_option_allowed()) then
       print "(a)",' ERROR: iresolution = 2 requires a single star and a non-isothermal wind'
       nerr = nerr + 1
    endif
 case default
    call read_inopt(default_particle_mass,'mass_of_particles',db,min=0.,errcount=nerr)
 end select

 if (isothermal) then
    wind_gamma = 1.
    call read_inopt(iwind,'iwind',db,min=2,max=3,default=2,errcount=nerr)
    call read_inopt(T_wind,'T_wind',db,min=0.,errcount=nerr)
    if (iwind == 3) then
       call read_inopt(temp_exponent,'temp_exponent',db,min=0.,max=5.,errcount=nerr)
    endif
 else
    call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)
 endif
 ! before the wind was a property of each star, the wind of the primary star
 ! (and the resolution) were set in the .in file
 call read_inopt(mdot,'primary_mdot',db,err=ierr_mdot)
 call close_db(db)
 if (ierr_mdot /= 0) call get_old_wind_options(filename)
 do i=1,nstar
    call check_stellar_props(trim(star_label(i)),star(i)%lum_lsun,star(i)%Teff,star(i)%Reff_au,&
                             i==1 .or. star(i)%mdot_msun_yr > 0.,nerr)
 enddo
 call check_wind_velocities(nerr)
 ierr = nerr

end subroutine read_setupfile

!----------------------------------------------------------------
!+
!  use the wind of the primary star and the particle mass from the
!  old options in the .in file, if any, so that they are written to
!  the .setup file
!+
!----------------------------------------------------------------
subroutine get_old_wind_options(setupfile)
 use infile_utils, only:open_db_from_file,inopts,close_db,infile_exists
 use inject,       only:print_old_wind_options
 character(len=*), intent(in) :: setupfile
 integer, parameter :: iunit = 22
 type(inopts), allocatable :: db(:)
 character(len=len(setupfile)) :: prefix
 real    :: mdot,vwind,wind_temp,rinj_au,alpha,pmass
 integer :: ierr

 prefix = setupfile(1:index(setupfile,'.setup',back=.true.)-1)
 if (.not.infile_exists(prefix)) return
 call open_db_from_file(db,trim(prefix)//'.in',iunit,ierr)
 if (ierr /= 0) return
 call print_old_wind_options(db,mdot,vwind,wind_temp,rinj_au,alpha,pmass)
 call close_db(db)
 if (mdot <= 0.) return

 print "(/,a,/)",' WARNING: using the wind of the primary star and the particle mass from the old options in '//&
                 trim(prefix)//'.in'
 star(1)%mdot_msun_yr = mdot
 star(1)%vwind_km_s   = vwind
 star(1)%wind_temp    = wind_temp
 star(1)%alpha        = alpha
 star(1)%rinject      = '0'
 if (rinj_au > 0.) write(star(1)%rinject,"(f0.4,a)") rinj_au,' au'
 if (pmass > 0.) then
    iresolution = 0
    default_particle_mass = pmass
 endif

end subroutine get_old_wind_options

!----------------------------------------------------------------
!+
!  read the orbit(s) of the stars. Binary setup files from before
!  the generic orbit options give semi_major_axis, eccentricity
!  and true_anomaly instead: these are used to set the orbit (the
!  file is rewritten with the generic options if phantomsetup is
!  run with a command line option, e.g. --nstar=2)
!+
!----------------------------------------------------------------
subroutine read_orbit_options(db,m,nerr)
 use infile_utils, only:inopts,read_inopt
 use setorbit,     only:read_options_orbits,set_defaults_orbit
 type(inopts), allocatable, intent(inout) :: db(:)
 real,                      intent(in)    :: m(:)
 integer,                   intent(inout) :: nerr
 character(len=len(orbit%elems%a)) :: a_string
 integer :: ierr

 ! the semi-major axis is kept as the string given in the file (in au, the code unit);
 ! err is non-zero only if the key is absent, and suppresses the error message
 a_string = ''
 call read_inopt(a_string,'semi_major_axis',db,err=ierr)
 if (nstar == 2 .and. ierr == 0) then
    call set_defaults_orbit(orbit)
    orbit%elems%a = a_string
    call read_inopt(orbit%e,'eccentricity',db,min=0.,errcount=nerr)
    call read_inopt(orbit%f,'true_anomaly',db,min=0.,max=360.,default=180.,errcount=nerr)
    orbit%w = 0.
 else
    call read_options_orbits(nstar,orbit,m,db,nerr)
 endif

end subroutine read_orbit_options

!----------------------------------------------------------------
!+
!  check that two of L, Teff and Reff are given, so that the third
!  can be derived. Stars that are not required may give none of them
!+
!----------------------------------------------------------------
subroutine check_stellar_props(label,lum,Teff,Reff,required,nerr)
 character(len=*), intent(in)    :: label
 real,             intent(in)    :: lum,Teff,Reff
 logical,          intent(in)    :: required
 integer,          intent(inout) :: nerr
 integer :: ngiven

 ngiven = 0
 if (lum  > 0.) ngiven = ngiven + 1
 if (Teff > 0.) ngiven = ngiven + 1
 if (Reff > 0.) ngiven = ngiven + 1
 if (ngiven == 1 .or. (required .and. ngiven == 0)) then
    print "(a)",' ERROR: '//label//' star needs two of '//label//'_lum, '//label//'_Teff and '//label//'_Reff'
    nerr = nerr + 1
 endif

end subroutine check_stellar_props

!----------------------------------------------------------------
!+
!  the wind type (given wind velocity or trans-sonic) is the same
!  for all stars, so either all stars with a wind give a wind
!  velocity or none of them do
!+
!----------------------------------------------------------------
subroutine check_wind_velocities(nerr)
 integer, intent(inout) :: nerr
 integer :: i,nwind,nvwind

 nwind  = 0
 nvwind = 0
 do i=1,nstar
    if (star(i)%mdot_msun_yr > 0.) then
       nwind = nwind + 1
       if (star(i)%vwind_km_s > 0.) nvwind = nvwind + 1
    endif
 enddo
 if (nvwind > 0 .and. nvwind < nwind) then
    print "(a)",' ERROR: either all stars with a wind give *_vwind > 0, or none do (trans-sonic winds)'
    nerr = nerr + 1
 endif

end subroutine check_wind_velocities

end module setup
