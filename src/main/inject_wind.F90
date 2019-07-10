!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION: Handles wind particle injection
!
!  REFERENCES: None
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    bowen_Cprime       -- radiative cooling rate (g.s/cm³)
!    bowen_Tcond        -- dust condensation temperature (K)
!    bowen_delta        -- condensation temperature range (K)
!    bowen_kappa        -- constant gas opacity (cm²/g)
!    bowen_kmax         -- maximum dust opacity (cm²/g)
!    iboundary_spheres  -- number of boundary spheres (integer)
!    iwind_resolution   -- if<>0 set number of particles on the sphere, reset particle mass
!    shift_spheres      -- delay before the ejection of shells
!    sonic_type         -- find transonic solution (1=yes,0=no)
!    star_Lum           -- central star luminosity (Lsun)
!    star_Teff          -- central star effective temperature (K)
!    wind_CO_ratio      -- wind initial C/O ratio
!    wind_alpha         -- fraction of the gravitational acceleration imparted to the gas
!    wind_expT          -- temperature law exponent (if wind_type=2)
!    wind_inject_radius -- radius of injection of the wind (au)
!    wind_mass_rate     -- wind mass loss rate (Msun/yr)
!    wind_osc_period    -- stellar pulsation period (days)
!    wind_osc_vampl     -- velocity amplitude of the pulsations (km/s)
!    wind_shell_spacing -- desired ratio of sphere spacing to particle spacing
!    wind_temperature   -- wind temperature at the injection point (K)
!    wind_type          -- stellar wind (1=isoT,2=T(r),3=adia,4=3+cooling,+10=with dust)
!    wind_velocity      -- velocity at which wind is injected (km/s)
!
!  DEPENDENCIES: bowen_dust, dim, dust_free_wind, dusty_wind, eos,
!    icosahedron, infile_utils, injectutils, io, part, physcon, timestep,
!    units, wind_profile
!+
!--------------------------------------------------------------------------
module inject
 use physcon, only: solarl
 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,radiativeforce
!
!--runtime settings for this module
!
! Read from input file
 integer, public:: iwind_resolution = 0
 integer, public:: iboundary_spheres = 3
 integer, public:: wind_type = 3
 real, public::    wind_expT = 0.5
 real, public::    wind_shell_spacing = 1.
 real, public::    wind_alpha = 0.
 real, public::    shift_spheres = 3.
 real, public::    star_Teff = 2500.
 real, public::    star_Lum = 10000. * solarl
 real, public::    wind_CO_ratio = 2
 real, public::    bowen_kmax = 2.7991
 real, public::    bowen_Tcond = 1500.
 real, public::    bowen_delta = 60.
 real, public::    bowen_kappa = 2.d-4
 real, public::    bowen_Cprime = 3.000d-5
#ifdef BOWEN
 integer, public:: sonic_type = 0
 real, public::    star_Teff = 3000.
 real, public::    star_Lum = 5315. * solarl
 real, public::    bowen_Cprime = 1.000d-5
 real, public::    wind_osc_period_days = 350.
 real, public::    wind_osc_vamplitude_km_s = 3.
 real, public::    wind_velocity_km_s = 0.
 real, public::    wind_mass_rate_Msun_yr = 1.04d-7
 real, public::    wind_temperature = 3000.
 real, public::    wind_injection_radius_au = 1.2568
#endif
#ifdef NUCLEATION
 integer, public:: sonic_type = 1
 real, public::    wind_velocity_km_s = 0.
 real, public::    wind_mass_rate_Msun_yr = 1.d-5
 real, public::    wind_temperature = 2500.
 real, public::    wind_injection_radius_au = 2.37686663
#elseif ISOTHERMAL
 integer, public:: sonic_type = 1
 real, public::    wind_velocity_km_s = 25.
 real, public::    wind_injection_radius_au = 0.465247264
#else
 integer, public:: sonic_type = 0
 real, public::    wind_velocity_km_s = 35.
 real, public::    wind_mass_rate_Msun_yr = 1.00d-8
 real, public::    wind_temperature = 3000.
 real, public::    wind_injection_radius_au = 1.7
#endif

 real, public::    u_to_temperature_ratio
 real, public::    mass_of_particles
! Calculated from the previous parameters
 real, public ::    wind_osc_vamplitude
#ifdef BOWEN
 real, public ::    wind_osc_period,dr3
#endif

 private

 real, private::   dtpulsation = 1.d99
 real, private ::  wind_mass_rate,wind_velocity,mass_of_spheres,time_between_spheres,neighbour_distance,&
      Rstar_cgs,Rstar,Cprime,wind_injection_radius
 integer, private :: particles_per_sphere,nwall_particles,iresolution

 logical, parameter :: wind_verbose = .false.
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3), rho_ini
 character(len=*), parameter :: label = 'inject_wind'

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,           only:fatal
 use timestep,     only:tmax
 use wind_profile, only:init_wind_equations
#ifdef NUCLEATION
 use dusty_wind,     only:setup_dustywind,get_initial_wind_speed
#else
 use dust_free_wind, only:setup_wind,get_initial_wind_speed
#endif
#ifdef BOWEN
 use bowen_dust,   only:init_bowen
#else
 use part,         only:xyzmh_ptmass
#endif
 use physcon,      only:mass_proton_cgs, kboltz, Rg, days, km, au, years, solarm, pi
 use icosahedron,  only:compute_matrices, compute_corners
 use eos,          only:gmw,gamma,polyk
 use units,        only:unit_velocity, umass, utime, udist
 use part,         only:massoftype,igas,iboundary
 use io,           only:iverbose
 use injectutils,  only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 integer, intent(out) :: ierr
 real :: mV_on_MdotR,wind_velocity_cgs,sonic(8)
 real :: dr,dp,mass_of_particles1
 logical :: verbose = .false.

 !
 ! return without error
 !
 ierr = 0
 if (iverbose>0) verbose = .true.
 !
 ! convert input parameters to code units
 !
#ifdef BOWEN
 wind_osc_period        = wind_osc_period_days * (days/utime)
 wind_osc_vamplitude    = wind_osc_vamplitude_km_s * (km / unit_velocity)
 dtpulsation            = wind_osc_period/50.
#else
 wind_osc_vamplitude    = 0.d0
#endif
 wind_velocity          = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate         = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif
 wind_velocity         = max(wind_velocity,wind_osc_vamplitude)
 Cprime = bowen_Cprime / (umass*utime/udist**3)
 wind_injection_radius  = wind_injection_radius_au * au / udist
 Rstar = Rstar_cgs/udist
#ifdef ISOTHERMAL
 wind_temperature = polyk* mass_proton_cgs/kboltz * unit_velocity**2*gmw
#endif

 call init_wind_equations (xyzmh_ptmass(4,wind_emitting_sink), star_Teff, Rstar, &
      wind_expT, u_to_temperature_ratio, wind_type)
#ifdef NUCLEATION
 call setup_dustywind(xyzmh_ptmass(4,wind_emitting_sink), star_Lum, star_Teff, Rstar_cgs, wind_CO_ratio,&
         bowen_Cprime, wind_expT, wind_mass_rate, u_to_temperature_ratio, wind_alpha, wind_type)

#else
 call setup_wind(xyzmh_ptmass(4,wind_emitting_sink), Rstar_cgs,bowen_Cprime, &
         wind_mass_rate, u_to_temperature_ratio, wind_alpha, wind_temperature, wind_type)
#endif

! integrate the wind equation to get the initial velocity required to set the resolution
 if (sonic_type ==  1) then
    ! Temperature profile
    !if (wind_type == 2 .and. wind_temperature /= star_Teff) then
    !   wind_injection_radius_au = Rstar_cgs/au*(star_Teff/wind_temperature)**(1./wind_expT)
    !endif
#ifdef ISOTHERMAL
    wind_injection_radius  = wind_injection_radius_au * au / udist
#else
    wind_injection_radius  = wind_injection_radius_au * au / udist !max(wind_injection_radius_au * au, Rstar_cgs) / udist
#endif
    !Rstar = min(wind_injection_radius_au*au,Rstar_cgs)
    call get_initial_wind_speed(wind_injection_radius*udist,wind_temperature,wind_velocity_cgs,sonic,.true.)
    wind_velocity = wind_velocity_cgs/unit_velocity
 endif


 if (iwind_resolution == 0) then
    !
    ! compute the dimensionless resolution factor m V / (Mdot R)
    ! where m = particle mass and V, Mdot and R are wind parameters
    !
    mass_of_particles = massoftype(igas)
    mV_on_MdotR = mass_of_particles*wind_velocity/(wind_mass_rate*wind_injection_radius)
    !
    ! solve for the integer resolution of the geodesic spheres
    ! gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
    !
    iresolution = get_sphere_resolution(wind_shell_spacing,mV_on_MdotR)
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    print *,'iwind_resolution equivalence = ',iresolution
 else
    if (wind_velocity == 0.)  then
       call fatal(label,'zero input wind velocity')
    endif
    iresolution = iwind_resolution
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    mass_of_particles = wind_shell_spacing*neighbour_distance * wind_injection_radius * &
         wind_mass_rate / (particles_per_sphere * wind_velocity)
    massoftype(igas) = mass_of_particles
    !print *,mass_of_particles,wind_injection_radius,wind_velocity
    print *,'iwind_resolution unchanged = ',iresolution
 endif

 mass_of_spheres = mass_of_particles * particles_per_sphere
 nwall_particles = iboundary_spheres*particles_per_sphere
 time_between_spheres  = mass_of_spheres / wind_mass_rate
 massoftype(iboundary) = mass_of_particles
 if (time_between_spheres > tmax)  then
    print *,'time_between_spheres = ',time_between_spheres,' < tmax = ',tmax
    call fatal(label,'no shell ejection : tmax < time_between_spheres')
 endif

 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 if (iverbose >= 1) then
    mass_of_particles1 = wind_shell_spacing * get_neighb_distance(4) * wind_injection_radius * wind_mass_rate &
                     / (get_parts_per_sphere(4) * wind_velocity)
    print*,'require mass of particle = ',mass_of_particles1,' to get 492 particles per sphere'
    dp = neighbour_distance*wind_injection_radius
    dr = wind_velocity*mass_of_particles*particles_per_sphere/wind_mass_rate
    print*,'particle separation on spheres = ',dp
    print*,'got dr/dp = ',dr/dp,' compared to desired dr on dp = ',wind_shell_spacing
 endif

#ifdef BOWEN
 call setup_bowen(u_to_temperature_ratio,bowen_kappa,bowen_kmax,star_Lum,wind_injection_radius,&
      bowen_Cprime,bowen_Tcond,bowen_delta,star_Teff,wind_osc_vamplitude,wind_osc_period,&
      iboundary_spheres*particles_per_sphere)
 rho_ini = wind_mass_rate / (4.*pi*wind_injection_radius**2*wind_velocity)
 dr3 = 3.*mass_of_spheres/(4.*pi*rho_ini)

 print*,'number of ejected shells per pulsation period (should at least be > 10) ',wind_osc_period/time_between_spheres
 print*,'width of the boundary layer/ R* (should be < 1) = ',1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius
 print*,'radial pulsation amplitude/ R* = ',wind_osc_vamplitude*wind_osc_period/(2.*pi)/wind_injection_radius
 print*,'pulsation period in code units = ',wind_osc_period
 !sanity checks
 ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
 if (wind_osc_period/time_between_spheres < 10. ) print *,'WARNING! only ',wind_osc_period/time_between_spheres,&
      ' shells will be ejected during a pulsation period'
 ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
 if (1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius > 0.2)  print*,'WARNING! the width of the boundary layer = ',&
      1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius,' Rinject'
#endif

 !logging
 print*,'mass_of_particles          = ',mass_of_particles
 print*,'mass_of_spheres            = ',mass_of_spheres
 print*,'distance between spheres   = ',wind_shell_spacing*neighbour_distance
 print*,'particles per sphere       = ',particles_per_sphere
 print*,'time_between_spheres       = ',time_between_spheres
 print*,'wind_injection_radius      = ',wind_injection_radius
 print*,'stellar_radius             = ',Rstar_cgs / udist

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use physcon,      only:pi
 use io,           only:fatal
 use eos,          only:gmw,gamma
#ifdef BOWEN
 use wind_profile, only:pulsating_wind_profile
#elif NUCLEATION
 use dusty_wind,   only:evolve_gail,dusty_wind_profile
 use part,         only:partJstarKmu
#else
 use dust_free_wind, only:stationary_wind_profile
#endif
 use units,        only:udist
 use part,         only:igas,iboundary,nptmass
 use injectutils,  only:inject_geodesic_sphere
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere, inner_sphere, inner_boundary_sphere, first_particle, i, ipart
 real    :: local_time, GM, r, v, u, rho, e, mu, mass_lost, x0(3), v0(3) !, cs2max, dr, dp
#ifdef NUCLEATION
 real :: Jstar, K(0:3), cs
 integer :: j
#elif BOWEN
 real :: surface_radius
#endif
 character(len=*), parameter :: label = 'inject_particles'

 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
    GM = xyzmh_ptmass(4,wind_emitting_sink)
    v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 else
    x0 = 0.
    v0 = 0.
    GM = 0.
 endif

#ifdef NUCLEATION
 !compute the values of K, Jstar and mu for the released particles & update u (if T=cst and mu has changed)
 call evolve_gail(dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,partJstarKmu,npart)
#else
 mu = gmw
#endif

 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1
 inner_sphere = floor(time/time_between_spheres)
 inner_boundary_sphere = inner_sphere + iboundary_spheres

 if (inner_sphere-outer_sphere > iboundary_spheres) call fatal(label,'problem with boundary spheres, timestep likely too large!')
! cs2max = 0.
 if (shift_spheres < 0) then
    ipart = iboundary
 else
    ipart = igas
 endif

 do i=inner_boundary_sphere,outer_sphere,-1
    local_time = time - (i-abs(shift_spheres)) * time_between_spheres

    !compute the radius, velocity, temperature, chemistry of a sphere at the current local time
#ifdef BOWEN
    call pulsating_wind_profile(time,local_time, r, v, u, rho, e, GM, gamma, sphere_number, &
         inner_sphere,inner_boundary_sphere)
#endif
    v = wind_velocity
#ifdef NUCLEATION
    r = Rstar
    call dusty_wind_profile(time,local_time, r, v, u, rho, e, GM, gamma, wind_temperature, Jstar, K, mu, cs)
#else
    r = wind_injection_radius
    call stationary_wind_profile(local_time, r, v, u, rho, e, GM, gamma, mu)
#endif

    if (wind_verbose) then
       print '("sphere ",5(i3),9(1x,es15.8))',i,inner_sphere,iboundary_spheres,outer_sphere,int(shift_spheres),&
            time,local_time,r/xyzmh_ptmass(5,1),v!,v/sqrt(gamma*(gamma-1.)*u),xyzmh_ptmass(5,1)*udist,dtlast
    endif

    if (i > inner_sphere) then
       ! boundary sphere
       first_particle = (iboundary_spheres-i+inner_sphere)*particles_per_sphere+1
       call inject_geodesic_sphere(i, first_particle, iresolution, r, v, u, rho,  geodesic_R, geodesic_V, &
            npart, npartoftype, xyzh, vxyzu, ipart, x0, v0)
#ifdef NUCLEATION
       do j=first_particle,first_particle+particles_per_sphere-1
          partJstarKmu(1,j) = Jstar
          partJstarKmu(2:5,j) = K(0:3)
          partJstarKmu(6,j) = mu
       enddo
#endif
    else
       ! ejected particles
#ifdef NUCLEATION
       do j=npart+1,npart+particles_per_sphere
          partJstarKmu(1,j) = Jstar
          partJstarKmu(2:5,j) = K(0:3)
          partJstarKmu(6,j) = mu
       enddo
#endif
       call inject_geodesic_sphere(i, npart+1, iresolution, r, v, u, rho, geodesic_R, geodesic_V,&
            npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
       ! update the sink particle mass
       if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
          xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_of_spheres
       endif
       print '(" ##### eject sphere ",4(i4),i7,9(1x,es12.5))',i,inner_sphere,iboundary_spheres,outer_sphere,npart,time,local_time,&
            r/xyzmh_ptmass(5,1),v,xyzmh_ptmass(5,1),u
       !stop

    endif
    !cs2max = max(cs2max,gamma*(gamma-1)*u)
 enddo
 ! update sink particle properties
 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_lost
#ifdef BOWEN
    if (wind_osc_vamplitude > 0.) then
       surface_radius = wind_injection_radius + wind_osc_vamplitude*wind_osc_period/(2.*pi)*sin(2.*pi*time/wind_osc_period)
       xyzmh_ptmass(5,wind_emitting_sink) = (surface_radius**3-dr3)**(1./3.)
    endif
#endif
 endif

 !
 ! return timestep constraint to ensure that time between sphere
 ! injections is adequately resolved
 !
 !dr = neighbour_distance*wind_injection_radius
 !dtinject = 0.25*dr/sqrt(cs2max)
 dtinject = min(0.2*time_between_spheres,dtpulsation)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Computes acceleration due to radiation pressure on particles
!+
!-----------------------------------------------------------------------
subroutine radiativeforce(x,y,z,xyzmh_ptmass,fextrad)

 real,    intent(in)  :: x,y,z
 real,    intent(in)  :: xyzmh_ptmass(:,:)
 real,    intent(out) :: fextrad(3)
 real :: dx,dy,dz,dist,fac

 dx = x-xyzmh_ptmass(1,wind_emitting_sink)
 dy = y-xyzmh_ptmass(2,wind_emitting_sink)
 dz = z-xyzmh_ptmass(3,wind_emitting_sink)
 dist = sqrt(dx**2 + dy**2 + dz**2)
 fac = wind_alpha*xyzmh_ptmass(4,wind_emitting_sink)/dist**3
 fextrad(1) = fac*dx
 fextrad(2) = fac*dy
 fextrad(3) = fac*dz

end subroutine radiativeforce

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use dim,          only: maxvxyzu
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_velocity_km_s,'wind_velocity','velocity at which wind is injected (km/s)',iunit)
#if defined (BOWEN)
 call write_inopt(wind_osc_period_days,'wind_osc_period','stellar pulsation period (days)',iunit)
 call write_inopt(wind_osc_vamplitude_km_s,'wind_osc_vampl','velocity amplitude of the pulsations (km/s)',iunit)
#endif
 call write_inopt(wind_injection_radius_au,'wind_inject_radius','radius of injection of the wind (au)',iunit)
#ifndef ISOTHERMAL
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass loss rate (Msun/yr)',iunit)
#endif
 if (maxvxyzu==4) then
    call write_inopt(wind_temperature,'wind_temperature','wind temperature at the injection point (K)',iunit)
 endif
 call write_inopt(iwind_resolution,'iwind_resolution','if<>0 set number of particles on the sphere, reset particle mass',iunit)
 call write_inopt(shift_spheres,'shift_spheres','delay before the ejection of shells',iunit)
 call write_inopt(wind_shell_spacing,'wind_shell_spacing','desired ratio of sphere spacing to particle spacing',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres','number of boundary spheres (integer)',iunit)
 call write_inopt(wind_alpha,'wind_alpha','fraction of the gravitational acceleration imparted to the gas',iunit)
#ifndef ISOTHERMAL
 write(iunit,"(/,a)") '# options controlling dust and cooling'
 call write_inopt(star_Lum/solarl,'star_Lum','central star luminosity (Lsun)',iunit)
 call write_inopt(star_Teff,'star_Teff','central star effective temperature (K)',iunit)
 call write_inopt(wind_type,'wind_type','stellar wind (1=isoT,2=T(r),3=adia,4=3+cooling,+10=with dust)',iunit)
 call write_inopt(wind_expT,'wind_expT','temperature law exponent (if wind_type=2)', iunit)
#ifdef NUCLEATION
 call write_inopt(wind_CO_ratio ,'wind_CO_ratio','wind initial C/O ratio',iunit)
#endif
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
#if defined (BOWEN) || defined(NUCLEATION)
 call write_inopt(bowen_kappa,'bowen_kappa','constant gas opacity (cm²/g)',iunit)
 call write_inopt(bowen_kmax,'bowen_kmax','maximum dust opacity (cm²/g)',iunit)
 call write_inopt(bowen_Tcond,'bowen_Tcond','dust condensation temperature (K)',iunit)
 call write_inopt(bowen_delta,'bowen_delta','condensation temperature range (K)',iunit)
#endif
#endif
 call write_inopt(sonic_type,'sonic_type','find transonic solution (1=yes,0=no)',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal, error, warning
 use physcon, only:pi,steboltz,au
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 integer :: noptions
 logical :: isowind = .true.
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('wind_velocity')
    read(valstring,*,iostat=ierr) wind_velocity_km_s
    ngot = ngot + 1
    if (wind_velocity_km_s < 0.)    call fatal(label,'invalid setting for wind_velocity (<0)')
 case('wind_inject_radius')
    read(valstring,*,iostat=ierr) wind_injection_radius_au
    ngot = ngot + 1
    if (wind_injection_radius_au < 0.) call fatal(label,'invalid setting for wind_inject_radius (<0)')
 case('wind_temperature')
    read(valstring,*,iostat=ierr) wind_temperature
    ngot = ngot + 1
    isowind = .false.
    if (wind_temperature < 0.)    call fatal(label,'invalid setting for wind_temperature (<0)')
 case('iwind_resolution')
    read(valstring,*,iostat=ierr) iwind_resolution
    ngot = ngot + 1
    if (iwind_resolution < 0) call fatal(label,'iwind_resolution must be bigger than zero')
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    ngot = ngot + 1
    if (iboundary_spheres <= 0) call fatal(label,'iboundary_spheres must be > 0')
 case('wind_shell_spacing')
    read(valstring,*,iostat=ierr) wind_shell_spacing
    ngot = ngot + 1
    if (wind_shell_spacing <= 0.) call fatal(label,'wind_shell_spacing must be >=0')
 case('shift_spheres')
    read(valstring,*,iostat=ierr) shift_spheres
    ngot = ngot + 1
 case('sonic_type')
    read(valstring,*,iostat=ierr) sonic_type
    ngot = ngot + 1
    if (sonic_type < 0 .or. sonic_type > 1 ) call fatal(label,'invalid setting for sonic_type ([0,1])')
 case('wind_alpha')
    read(valstring,*,iostat=ierr) wind_alpha
    ngot = ngot + 1
    if (wind_alpha < 0.) call fatal(label,'invalid setting for wind_alpha (must be > 0)')
#ifndef ISOTHERMAL
 case('wind_mass_rate')
    read(valstring,*,iostat=ierr) wind_mass_rate_Msun_yr
    ngot = ngot + 1
    if (wind_mass_rate_Msun_yr < 0.)    call fatal(label,'invalid setting for wind_mass_rate (<0)')
 case('star_Lum')
    read(valstring,*,iostat=ierr) star_Lum
    star_Lum = star_Lum * solarl
    ngot = ngot + 1
    if (star_Lum < 0.) call fatal(label,'invalid setting for star_Lum (<0)')
 case('star_Teff')
    read(valstring,*,iostat=ierr) star_Teff
    ngot = ngot + 1
    if (star_Teff < 0.)    call fatal(label,'invalid setting for star_Teff (<0)')
 case('bowen_kappa')
    read(valstring,*,iostat=ierr) bowen_kappa
    ngot = ngot + 1
    if (bowen_kappa < 0.)    call fatal(label,'invalid setting for bowen_kappa (<0)')
 case('bowen_kmax')
    read(valstring,*,iostat=ierr) bowen_kmax
    ngot = ngot + 1
    if (bowen_kmax < 0.)    call fatal(label,'invalid setting for bowen_kmax (<0)')
 case('bowen_Tcond')
    read(valstring,*,iostat=ierr) bowen_Tcond
    ngot = ngot + 1
    if (bowen_Tcond < 0.) call fatal(label,'invalid setting for bowen_Tcond (<0)')
 case('bowen_delta')
    read(valstring,*,iostat=ierr) bowen_delta
    ngot = ngot + 1
    if (bowen_delta < 0.) call fatal(label,'invalid setting for bowen_delta (<0)')
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
    if (bowen_Cprime < 0.) call fatal(label,'invalid setting for bowen_Cprime (<0)')
 case('wind_type')
    read(valstring,*,iostat=ierr) wind_type
    ngot = ngot + 1
    if (wind_type < 0 .or. wind_type > 4 ) call fatal(label,'invalid setting for wind_type ([0,3])')
 case('wind_expT')
    read(valstring,*,iostat=ierr) wind_expT
    ngot = ngot + 1
#endif
#ifdef BOWEN
 case('wind_osc_period')
    read(valstring,*,iostat=ierr) wind_osc_period_days
    ngot = ngot + 1
    if (wind_osc_period_days < 0.) call fatal(label,'invalid setting for wind_osc_period (<0)')
 case('wind_osc_vampl')
    read(valstring,*,iostat=ierr) wind_osc_vamplitude_km_s
    !wind_velocity_km_s = 0. ! set wind veolicty to zero when pulsating star
    ngot = ngot + 1
    if (wind_osc_vamplitude_km_s <= 0. .and. wind_velocity_km_s <= 0.) call fatal(label,'invalid setting for wind_osc_vamp (<0)')
#endif
#ifdef NUCLEATION
 case('wind_CO_ratio')
    read(valstring,*,iostat=ierr) wind_CO_ratio
    ngot = ngot + 1
    if (wind_CO_ratio < 0.) call fatal(label,'invalid setting for wind_CO_ratio (must be > 0)')
#endif
 case default
    imatch = .false.
 end select
#ifdef BOWEN
 noptions = 18
#else
 noptions = 15
#ifdef NUCLEATION
 noptions = 16
#elif ISOTHERMAL
 noptions = 8
#endif
#endif
!print *,isowind,trim(name),ngot,noptions
 igotall = (ngot >= noptions)
 if (igotall .and. trim(name) /= '') then
#ifdef ISOTHERMAL
    Rstar_cgs = wind_injection_radius_au *au
#else
    Rstar_cgs = sqrt(star_Lum/(4.*pi*steboltz*star_Teff**4))
    !if you launch the wind inside the photosphere, make sure the wind temperature >= star's effective temperature
    if (wind_injection_radius_au <= Rstar_cgs/au) then
       print *,'invalid setting for wind_inject_radius < Rstar (au)',wind_injection_radius_au,Rstar_cgs/au
       !call fatal(label,'invalid setting for wind_inject_radius (< Rstar)')
    endif
    ! if (wind_temperature < star_Teff) then
    !    print *,'Twind < star_Teff',wind_temperature,star_Teff
    !    !call fatal(label,'invalid setting for wind_temperature (Twind < star_Teff)')
    ! endif
#endif
    if (iboundary_spheres > int(shift_spheres)) then
       print *,'shift_spheres too small - imposing shift_spheres = iboundary_spheres = ',iboundary_spheres
       shift_spheres = sign(dble(iboundary_spheres),shift_spheres)
    endif
 endif
 if (trim(name) == '') ngot = 0
end subroutine read_options_inject

end module inject
