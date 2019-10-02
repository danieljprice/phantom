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
!    iboundary_spheres  -- number of boundary spheres (integer)
!    iwind_resolution   -- if<>0 set number of particles on the sphere, reset particle mass
!    outer_boundary     -- kill gas particles outside this radius (au)
!    piston_velocity    -- velocity amplitude of the pulsation (km/s)
!    pulsation_period   -- stellar pulsation period (days)
!    sonic_type         -- find transonic solution (1=yes,0=no)
!    wind_inject_radius -- wind injection radius (au)
!    wind_mass_rate     -- wind mass loss rate (Msun/yr)
!    wind_shell_spacing -- desired ratio of sphere spacing to particle spacing
!    wind_temperature   -- wind temperature at the injection point (K)
!    wind_velocity      -- injection wind velocity (km/s, if sonic_type = 0)
!
!  DEPENDENCIES: dim, eos, icosahedron, infile_utils, injectutils, io,
!    options, part, partinject, physcon, ptmass_radiation, timestep, units,
!    wind, wind_equations
!+
!--------------------------------------------------------------------------
module inject
 use physcon, only: solarl
 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject
!
!--runtime settings for this module
!
! Read from input file
 integer, public:: iboundary_spheres = 5
 integer, public:: iwind_resolution = 0
 real, public::    outer_boundary_au = 10.
 real, public::    wind_shell_spacing = 1.
 real, public::    pulsation_period
 real, public::    pulsation_period_days = 0.
 real, public::    piston_velocity_km_s = 0.
#ifdef NUCLEATION
 integer, public:: sonic_type = 1
 real, public::    wind_velocity_km_s = 0.
 real, public::    wind_mass_rate_Msun_yr = 1.d-5
 real, public::    wind_injection_radius_au = 2.37686663
 real, public::    wind_temperature = 2500.
#elif ISOTHERMAL
 integer, public:: sonic_type = 1
 real, public::    wind_velocity_km_s = 25.
 real, public::    wind_mass_rate_Msun_yr = 1.d-8
 real, public::    wind_injection_radius_au = 0.46524726
 real, public::    wind_temperature
#else
 integer, public:: sonic_type = 0
 real, public::    wind_velocity_km_s = 35.
 real, public::    wind_mass_rate_Msun_yr = 1.d-8
 real, public::    wind_injection_radius_au = 1.7
 real, public::    wind_temperature = 3000.
#endif

 real, public::    Rstar

 private

 real :: dtpulsation = 1.d99
 real :: u_to_temperature_ratio,wind_mass_rate,piston_velocity,wind_velocity,&
      mass_of_spheres,time_between_spheres,neighbour_distance,mass_of_particles,&
      dr3,Rstar_cgs,wind_injection_radius,rho_ini, omega_osc, deltaR_osc, Mstar_cgs
 integer :: particles_per_sphere,nwall_particles,iresolution,nwrite,nreleased

 logical :: pulsating_wind
 logical, parameter :: wind_verbose = .false.
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 character(len=*), parameter :: label = 'inject_wind'

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use options,        only:icooling
 use io,             only:fatal
 use timestep,       only:tmax,dtmax
 use wind_equations, only:init_wind_equations
 use part,           only:xyzmh_ptmass
 use wind,           only:setup_wind
 use physcon,      only:mass_proton_cgs, kboltz, Rg, days, km, au, years, solarm, pi
 use icosahedron,  only:compute_matrices, compute_corners
 use eos,          only:gmw,gamma,polyk
 use units,        only:unit_velocity, umass, utime, udist
 use part,         only:massoftype,igas,iboundary,xyzmh_ptmass,imloss,ilum,iTeff,iReff
 use io,           only:iverbose
 use injectutils,  only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 integer, intent(out) :: ierr
 integer :: ires_min,nzones_per_sonic_point
 real :: mV_on_MdotR,initial_wind_velocity_cgs,sonic(8),rho_inj,dist_to_sonic_point
 real :: dr,dp,mass_of_particles1,Rinject

 !
 ! return without error
 !
#ifdef NUCLEATION
 nwrite = 19
#else
 nwrite = 10
#endif
 if (icooling > 0) nwrite = nwrite+1
 ierr = 0

 pulsating_wind = (pulsation_period_days > 0.) .and. (piston_velocity_km_s > 0.)
 if (pulsating_wind) then
    nreleased = 1
 else
    nreleased = 30
 endif

 !
 ! convert input parameters to code units
 !
 if (pulsating_wind) then
    pulsation_period = pulsation_period_days * (days/utime)
    piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
    dtpulsation      = pulsation_period/50.
    omega_osc        = 2.*pi/pulsation_period
    deltaR_osc       = pulsation_period*piston_velocity/(2.*pi)
    sonic_type       = 0
 else
    piston_velocity  = 0.d0
 endif

 wind_velocity    = wind_velocity_km_s * (km / unit_velocity)
 !wind_velocity    = max(wind_velocity,piston_velocity)
 wind_mass_rate   = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 wind_injection_radius  = wind_injection_radius_au * au / udist
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif

 Rstar_cgs = xyzmh_ptmass(iReff,wind_emitting_sink)*udist
 Mstar_cgs = xyzmh_ptmass(4,wind_emitting_sink)*umass
#ifdef ISOTHERMAL
 Rstar_cgs = wind_injection_radius_au *au
 wind_temperature = polyk* mass_proton_cgs/kboltz * unit_velocity**2*gmw
#endif
 if (wind_injection_radius_au <= Rstar_cgs/au) then
    print *,'invalid setting for wind_inject_radius < Rstar (au)',wind_injection_radius_au,Rstar_cgs/au
    !call fatal(label,'invalid setting for wind_inject_radius (< Rstar)')
 endif
 if (piston_velocity_km_s <= 0. .and. wind_velocity_km_s <= 0. .and. sonic_type == 0) then
    call fatal(label,'invalid setting for the wind velocity parameters (v=0)')
 endif
 if ( .not. pulsating_wind) then
    Rinject = xyzmh_ptmass(iReff,wind_emitting_sink)
    call init_wind_equations (xyzmh_ptmass(4,wind_emitting_sink), Rinject, &
         xyzmh_ptmass(iTeff,wind_emitting_sink), u_to_temperature_ratio)
    call setup_wind(xyzmh_ptmass(4,wind_emitting_sink), Rstar_cgs, wind_mass_rate, &
         u_to_temperature_ratio, wind_temperature)
!if (ieos == 6 .and. wind_temperature /= star_Teff) then
!    wind_injection_radius_au = Rstar_cgs/au*(star_Teff/wind_temperature)**(1./wind_expT)
!    wind_injection_radius  = wind_injection_radius_au * au / udist
!endif
!wind_injection_radius = max(wind_injection_radius_au * au, Rstar_cgs) / udist
!Rstar = min(wind_injection_radius_au*au,Rstar_cgs)

    ! integrate wind equation to get initial velocity required to set the resolution
    initial_wind_velocity_cgs = wind_velocity_km_s*1.e5
    call get_initial_wind_speed(wind_injection_radius*udist,wind_temperature,&
         initial_wind_velocity_cgs,sonic,sonic_type)
    wind_velocity = initial_wind_velocity_cgs/unit_velocity

    !save 1D initial profile for comparison
    call save_windprofile(wind_injection_radius*udist,initial_wind_velocity_cgs,&
         wind_temperature, sonic(4),'windprofile1D.dat')
 endif

 if (iwind_resolution == 0) then
    !
    ! resolution is specified in terms of number of smoothing lengths
    ! per distance to sonic point
    !
    nzones_per_sonic_point = 8
    dist_to_sonic_point = sonic(1)/udist-wind_injection_radius
    dr = abs(dist_to_sonic_point)/nzones_per_sonic_point
    rho_inj = wind_mass_rate/(4.*pi*wind_injection_radius**2*(piston_velocity+wind_velocity))
    mass_of_particles = rho_inj*dr**3
    massoftype(igas) = mass_of_particles
    print*,' suggesting ',mass_of_particles, ' based on desired dr = ',dr,' dist-to-sonic=',dist_to_sonic_point
    !
    ! compute the dimensionless resolution factor m V / (Mdot R)
    ! where m = particle mass and V, Mdot and R are wind parameters
    !
!    mass_of_particles = massoftype(igas)
    mV_on_MdotR = mass_of_particles*(piston_velocity+wind_velocity)/(wind_mass_rate*wind_injection_radius)
    !
    ! solve for the integer resolution of the geodesic spheres
    ! gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
    !
    iresolution = get_sphere_resolution(wind_shell_spacing,mV_on_MdotR)
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    print *,'iwind_resolution equivalence = ',iresolution
 else
    if (wind_velocity+piston_velocity == 0.) call fatal(label,'zero input wind velocity')
    iresolution = iwind_resolution
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    mass_of_particles = wind_shell_spacing*neighbour_distance * wind_injection_radius * &
         wind_mass_rate / (particles_per_sphere * (wind_velocity+piston_velocity))
    massoftype(igas) = mass_of_particles
    print *,'iwind_resolution unchanged = ',iresolution
 endif

 mass_of_spheres = mass_of_particles * particles_per_sphere
 rho_ini = wind_mass_rate / (4.*pi*wind_injection_radius**2*(piston_velocity+wind_velocity))
 dr3 = 3.*mass_of_spheres/(4.*pi*rho_ini)
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
         / (get_parts_per_sphere(4) * (piston_velocity+wind_velocity))
    print*,'required mass of particle  = ',mass_of_particles1,' to get 492 particles per sphere'
    dp = neighbour_distance*wind_injection_radius
    dr = (piston_velocity+wind_velocity)*mass_of_particles*particles_per_sphere/wind_mass_rate
    print*,'particle separation on spheres = ',dp
    print*,'got dr/dp = ',dr/dp,' compared to desired dr on dp = ',wind_shell_spacing
 endif


 !logging
 print*,'mass_of_particles          = ',mass_of_particles
 print*,'particles per sphere       = ',particles_per_sphere
 print*,'distance between spheres   = ',wind_shell_spacing*neighbour_distance
 if (sonic_type == 1) then
    print*,'distance to sonic point    = ',sonic(1)/udist-wind_injection_radius
    print*,'sonic radius               = ',sonic(1)/udist,sonic(1)
    print*,'number of shells to sonic  = ',(sonic(1)/udist-wind_injection_radius)/(wind_shell_spacing*neighbour_distance)
    print*,'time_to_sonic_point        = ',sonic(4)/utime
 endif
 print*,'time_between_spheres       = ',time_between_spheres
 print*,'wind_temperature           = ',wind_temperature
 print*,'wind_injection_radius      = ',wind_injection_radius
 print*,'stellar_radius             = ',Rstar_cgs / udist
 if (pulsating_wind) then
    print*,'number of ejected shells per pulsation period (should at least be > 10) ',pulsation_period/time_between_spheres
    print*,'width of the boundary layer/ R* (should be < 1) = ',1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius
    print*,'radial pulsation amplitude/ R* = ',piston_velocity*pulsation_period/(2.*pi*wind_injection_radius)
    print*,'pulsation period in code units = ',pulsation_period
    !sanity checks
    ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
    if (pulsation_period/time_between_spheres < 10. ) print *,'WARNING! only ',pulsation_period/time_between_spheres,&
         ' shells will be ejected during a pulsation period'
    ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
    if (1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius > 0.2)  print*,'WARNING! the width of the boundary layer = ',&
         1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius,' Rinject'
 else
    !save a few models before the particles reach the sonic point
    if (dtmax > sonic(4)/utime) print *,'WARNING! dtmax > time to sonic point'
    !minimum resolution required so a few shells can be inserted between the injection radius and the sonic point
    ires_min = iboundary_spheres*wind_shell_spacing*0.5257/(sonic(1)/udist/wind_injection_radius-1.)+.5
    if (iwind_resolution < ires_min) print *,'WARNING! resolution too low to pass sonic point : iwind_resolution < ',ires_min
 endif

 xyzmh_ptmass(imloss,wind_emitting_sink) = wind_mass_rate

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use physcon,      only:pi,au
 use io,           only:fatal,iverbose
 use dim,          only:store_dust_temperature
#ifdef NUCLEATION
 use wind,         only:dusty_wind_profile
#else
 use wind,         only:dust_free_wind_profile
#endif
 use part,         only:igas,iTeff,iboundary,nptmass,delete_particles_outside_sphere,dust_temp
 use partinject,   only:add_or_update_particle
 use injectutils,  only:inject_geodesic_sphere
 use units,        only:udist
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere, inner_sphere, inner_boundary_sphere, first_particle, i, ipart, &
            nshell_released, nboundaries
 real    :: local_time, GM, r, v, u, rho, e, mass_lost, x0(3), v0(3), surface_radius
 character(len=*), parameter :: label = 'inject_particles'
 logical, save :: released = .false.
#ifdef NUCLEATION
 real :: JKmuS(7)
#endif

 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
    GM = xyzmh_ptmass(4,wind_emitting_sink)
    v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 else
    x0 = 0.
    v0 = 0.
    GM = 0.
 endif

 !
 ! delete particles that exit the outer boundary
 !
 call delete_particles_outside_sphere(x0,outer_boundary_au*au/udist,npart)

 if (npart > 0) then
    nshell_released = nreleased
    nboundaries = iboundary_spheres
    !release particles and declare inner boundary shells as gas particles so they can exert some pressure
    ipart = igas
    if (.not.released) then
       do i = npart-nshell_released*particles_per_sphere+1,npart
          call add_or_update_particle(igas,xyzh(1:3,i),vxyzu(1:3,i),xyzh(4,i),vxyzu(4,i),i,npart,npartoftype,xyzh,vxyzu)
       enddo
       released = .true.
    endif
 else
    !initialise domain with boundary particles
    ipart = iboundary
    nshell_released = 0
    nboundaries = iboundary_spheres+nreleased
 endif
 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1 + nshell_released
 inner_sphere = floor(time/time_between_spheres)+nshell_released
 inner_boundary_sphere = inner_sphere + nboundaries

 !only one sphere can be ejected at a time
 if (inner_sphere-outer_sphere > nboundaries) call fatal(label,'ejection of more than 1 sphere, timestep likely too large!')

 do i=inner_boundary_sphere,outer_sphere,-1
    local_time = time + (iboundary_spheres+nreleased-i) * time_between_spheres

    !compute the radius, velocity, temperature, chemistry of a sphere at the current local time
    v = wind_velocity
    if (pulsating_wind) then
       call pulsating_wind_profile(time,local_time, r, v, u, rho, e, GM, i, &
            inner_sphere,inner_boundary_sphere,dr3,rho_ini)
    else
#ifdef NUCLEATION
       r = Rstar
       call dusty_wind_profile(time,local_time, r, v, u, rho, e, GM, wind_temperature, JKmuS)
#else
       r = wind_injection_radius
       call dust_free_wind_profile(local_time, r, v, u, rho, e, GM)
#endif
       if (iverbose > 0) print '(" ##### boundary sphere ",i4,3(i4),i7,9(1x,es12.5))',i,&
            inner_sphere,nboundaries,outer_sphere,npart,time,local_time,r,wind_injection_radius,v
    endif

    if (i > inner_sphere) then
       ! boundary sphere
       first_particle = (nboundaries-i+inner_sphere)*particles_per_sphere+1
       !print '(" ##### boundary sphere ",i4,i7,3(i4),i7,9(1x,es12.5))',i,first_particle,inner_sphere,nboundaries,&
       !     outer_sphere,npart,time,local_time,r/xyzmh_ptmass(5,1),v,u,rho
#ifdef NUCLEATION
       call inject_geodesic_sphere(i, first_particle, iresolution, r, v, u, rho,  geodesic_R, geodesic_V, &
            npart, npartoftype, xyzh, vxyzu, ipart, x0, v0, JKmuS)
#else
       call inject_geodesic_sphere(i, first_particle, iresolution, r, v, u, rho,  geodesic_R, geodesic_V, &
            npart, npartoftype, xyzh, vxyzu, ipart, x0, v0)
#endif
       if (store_dust_temperature) dust_temp(first_particle:first_particle+particles_per_sphere-1) = &
            xyzmh_ptmass(iTeff,wind_emitting_sink)
    else
       ! ejected particles
#ifdef NUCLEATION
       call inject_geodesic_sphere(i, npart+1, iresolution, r, v, u, rho, geodesic_R, geodesic_V,&
            npart, npartoftype, xyzh, vxyzu, igas, x0, v0, JKmuS)
#else
       call inject_geodesic_sphere(i, npart+1, iresolution, r, v, u, rho, geodesic_R, geodesic_V,&
            npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
#endif
       !initialize dust temperature to star's effective temperature
       if (store_dust_temperature) dust_temp(npart+1:npart+particles_per_sphere) = xyzmh_ptmass(iTeff,wind_emitting_sink)
       ! update the sink particle mass
       if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
          xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_of_spheres
       endif
       print '(" ##### eject sphere ",4(i4),i7,9(1x,es12.5))',i,inner_sphere,nboundaries,&
            outer_sphere,npart,time,local_time,r/xyzmh_ptmass(5,1),v,u,rho
    endif
    !cs2max = max(cs2max,gamma*(gamma-1)*u)
 enddo
 ! update sink particle properties
 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_lost
    if (pulsating_wind) then
       surface_radius = wind_injection_radius + piston_velocity*pulsation_period/(2.*pi)*sin(2.*pi*time/pulsation_period)
       !v2 xyzmh_ptmass(5,wind_emitting_sink) = (surface_radius**3-dr3)**(1./3.)
       xyzmh_ptmass(5,wind_emitting_sink) = sqrt(xyzh(1,1)**2+xyzh(2,1)**2+xyzh(3,1)**2)
    endif
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
!  Oscillating inner boundary
!+
!-----------------------------------------------------------------------
subroutine pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,sphere_number, &
                                  inner_sphere,inner_boundary_sphere,dr3,rho_ini)
 use physcon,     only:pi
 use eos,         only:gamma
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,GM,dr3,rho_ini
 real,    intent(out) :: r, v, u, rho, e

 integer, parameter :: nrho_index = 10
 integer :: k
 real :: surface_radius,r3
 logical :: verbose = .true.

 v = wind_velocity + piston_velocity* cos(omega_osc*time) !same velocity for all wall particles
 surface_radius = wind_injection_radius + deltaR_osc*sin(omega_osc*time)
 !ejected spheres
 if (sphere_number <= inner_sphere) then
    r = surface_radius
    v = max(piston_velocity,wind_velocity)
 else
    !boundary spheres
    r3 = surface_radius**3-dr3
    do k = 2,sphere_number-inner_sphere
       r3 = r3-dr3*(r3/surface_radius**3)**(nrho_index/3.)
    enddo
    r = r3**(1./3)
 endif
 !r = (surface_radius**3-(sphere_number-inner_sphere)*dr3)**(1./3)
 !rho = rho_ini
 u = wind_temperature * u_to_temperature_ratio
 if (gamma > 1.0001) then
    e = .5*v**2 - GM/r + gamma*u
 else
    e = .5*v**2 - GM/r + u
 endif
 rho = rho_ini*(surface_radius/r)**nrho_index
 if (verbose) then
    if (sphere_number > inner_sphere) then
       print '("boundary, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4)', &
             sphere_number,inner_sphere,surface_radius,r,v,time/pulsation_period
    else
       print '("ejected, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4)', &
             sphere_number,inner_sphere,surface_radius,r,v,time/pulsation_period
    endif
 endif

end subroutine pulsating_wind_profile

!-----------------------------------------------------------------------
!
!  Determine the initial wind speed
!    stype = 0 : do nothing, initial velocity is set
!    stype = 1 : determine the trans-sonic solution
!
!-----------------------------------------------------------------------
subroutine get_initial_wind_speed(r0, T0, v0, sonic, stype)
!all quantities in cgs
 use timestep, only:tmax
 use io,       only:fatal,iverbose
 use units,    only:utime,udist
 use eos,      only:gmw,gamma
 use physcon,  only:Rg,Gg,au,years
 use wind,     only:wind_state,calc_wind_profile
 use ptmass_radiation, only:alpha_rad
 integer, intent(in) :: stype
 real, intent(in)    :: r0, T0
 real, intent(inout) :: v0
 real, intent(out)   :: sonic(8)

 type(wind_state) :: state

 real :: v0min, v0max, v0last, vesc, cs, Rs, alpha_max,vin
 integer, parameter :: ncount_max = 20
 integer :: icount
 character(len=*), parameter :: label = 'get_initial_wind_speed'

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-alpha_rad)/r0)
 cs   = sqrt(gamma*Rg*T0/gmw)
 vin  = cs*(vesc/2./cs)**2*exp(-(vesc/cs)**2/2.+1.5)
 Rs   = Gg*Mstar_cgs*(1.-alpha_rad)/(2.*cs*cs)
 if (iverbose>0) then
    if (vesc> 1.d-50) then
       alpha_max = 1.-(2.*cs/vesc)**2
    else
       alpha_max = 0.
    endif
    print *, "[get_initial_wind_speed] searching for initial velocity."
    print *, ' * stype      = ',stype
    print *, ' * unit(au)   = ',udist/au
    print *, ' * Mstar      = ',Mstar_cgs/1.9891d33
    print *, ' * Twind      = ',T0
#ifndef ISOTHERMAL
    print *, ' * Rstar(au)  = ',Rstar_cgs/1.496d13,Rstar_cgs/69600000000.
    print *, ' * r0(au)     = ',r0/1.496d13,r0/69600000000.
    print *, ' * gamma      = ',gamma
#else
    print *, ' * Rstar (Ro) = ',r0/69600000000.,Rstar_cgs/69600000000.
#endif
    print *, ' * mu         = ',gmw
    print *, ' * cs  (km/s) = ',cs/1e5
    print *, ' * vesc(km/s) = ',vesc/1e5
    if (stype == 1) then
       print *, ' * v0  (km/s) = ',vin/1e5
    else
       print *, ' * v0  (km/s) = ',v0/1e5
    endif
    print *, ' * alpha      = ',alpha_rad
    print *, ' * alpha_max  = ',alpha_max
    print *, ' * tend (s)   = ',tmax*utime,tmax*utime/years
 endif

 !
 ! seach for trans-sonic solution
 !
 if (stype == 1) then

    ! Find lower bound for initial velocity
    v0 = cs*0.991
    v0max = v0
    icount = 0
    do while (icount < ncount_max)
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *,' v0/cs = ',v0/cs
       if (state%spcode == -1) then
          v0min = v0
          exit
       else
          v0max = v0
          v0 = v0 / 2.
       endif
       icount = icount+1
    enddo
    if (iverbose>1) print *, 'Lower bound found for v0/cs :',v0min/cs
    if (icount == ncount_max) call fatal(label,'cannot find v0min, change wind_temperature or wind_injection_radius ?')
    if (v0min/cs > 0.99) call fatal(label,'supersonic wind solution, set sonic_type = 0 and provide wind_velocity')

    ! Find upper bound for initial velocity
    v0 = v0max
    icount = 0
    do while (icount < ncount_max)
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *,' v0/cs = ',v0/cs
       if (state%spcode == 1) then
          v0max = v0
          exit
       else
          v0min = max(v0min, v0)
          v0 = v0 * 1.1
       endif
       icount = icount+1
    enddo
    if (icount == ncount_max) call fatal(label,'cannot find v0max, change wind_temperature or wind_injection_radius ?')
    if (iverbose>1) print *, 'Upper bound found for v0/cs :', v0max/cs

    ! Find sonic point by dichotomy between v0min and v0max
    do
       v0last = v0
       v0 = (v0min+v0max)/2.
       call calc_wind_profile(r0, v0, T0, 0., state)
       if (iverbose>1) print *, 'v0/cs = ',v0/cs
       if (state%spcode == -1) then
          v0min = v0
       elseif (state%spcode == 1) then
          v0max = v0
       else
          exit
       endif
       if (abs(v0-v0last)/v0last < 1.e-5) then
          exit
       endif
    enddo
    !
    !store sonic point properties (location, time to reach, ...)
    !
    sonic(1) = state%r
    sonic(2) = state%v
    sonic(3) = state%c
    sonic(4) = state%time
    sonic(5) = state%Tg
    sonic(6) = state%p
    sonic(7) = state%alpha
    !mdot = 4.*pi*rho*v0*ro*ro

    write (*,'("Sonic point properties  cs (km/s) =",f9.3,", Rs/r0 = ",f7.3,", Rs/R* = ",f7.3,", v0/cs = ",f9.6,", ts =",f8.1)') &
         sonic(2)/1e5,sonic(1)/r0,sonic(1)/Rstar_cgs,v0/sonic(2),sonic(4)/utime

 else
    if (v0 >= cs) then
       print *,' supersonic wind : v0/cs = ',v0/cs
    else
       print *,' sub-sonic wind : v0/cs = ',v0/cs
    endif
 endif

end subroutine get_initial_wind_speed

!-----------------------------------------------------------------------
!
!  Integrate the steady wind equation and save wind profile to a file
!
!-----------------------------------------------------------------------
subroutine save_windprofile(r0, v0, T0, tsonic, filename)
 use units,    only:utime
 use timestep, only:tmax
 use wind,     only:wind_state,wind_step,init_wind
 real, intent(in) :: r0, v0, T0, tsonic
 character(*), intent(in) :: filename
 real, parameter :: Tdust_stop = 1.d0 ! Temperature at outer boundary of wind simulation
 real :: dt_print,time_end
 type(wind_state) :: state
 integer :: n,iter

 write (*,'("Saving 1D model : ")')
 time_end = tmax*utime*2.
 call init_wind(r0, v0, T0, time_end, state)
 open(unit=1337,file=filename)
 call filewrite_header(1337)
 call filewrite_state(1337, state)

 n = 1
 iter = 0
 dt_print = max(min(tsonic/10.,time_end/256.),time_end/5000.)
 do while(state%time < time_end .and. iter < 1000000 .and. state%Tg > Tdust_stop)
    iter = iter+1
    call wind_step(state)
    if (state%time > n*dt_print) then
       n = floor(state%time/dt_print)+1
       call filewrite_state(1337, state)
    endif
 enddo
 if (state%time/time_end < .3) then
    write(*,'(/,"[WARNING] wind integration failed : t/tend = ",f7.5,", dt/tend = ",f7.5," Tgas = ",f6.0," iter = ",i7,/)') &
         state%time/time_end,state%dt/time_end,state%Tg,iter
 endif
 close(1337)
end subroutine save_windprofile

subroutine filewrite_header(iunit)
 use options, only : icooling
 integer, intent(in) :: iunit

#ifdef NUCLEATION
 if (icooling > 0) then
    write(iunit,'("#",11x,a1,19(a20))') 't','r','v','T','c','p','rho','alpha','a',&
         'mu','S','Jstar','K0','K1','K2','K3','tau_lucy','kappa_planck','kappa_ross','Q'
 else
    write(iunit,'("#",11x,a1,18(a20))') 't','r','v','T','c','p','rho','alpha','a',&
         'mu','S','Jstar','K0','K1','K2','K3','tau_lucy','kappa_planck','kappa_ross'
 endif
#else
 if (icooling > 0) then
    write(iunit,'("#",11x,a1,10(a20))') 't','r','v','T','c','p','rho','alpha','a','mu','Q'
 else
    write(iunit,'("#",11x,a1,9(a20))') 't','r','v','T','c','p','rho','alpha','a','mu'
 endif
#endif
end subroutine filewrite_header

subroutine state_to_array(state, array)
  use options, only : icooling
  use wind,     only:wind_state
 type(wind_state), intent(in) :: state
 real, intent(out) :: array(:)
#ifdef NUCLEATION
 real :: f
#endif

 array(1) = state%time
 array(2) = state%r
 array(3) = state%v
 array(4) = state%Tg
 array(5) = state%c
 array(6) = state%p
 array(7) = state%rho
 array(8)  = state%alpha
 array(9) = state%a
#ifdef NUCLEATION
 f = state%r**2 * state%v
 array(10)  = state%JKmuS(6)
 array(11)  = state%JKmuS(7)
 array(12) = state%JKmuS(1)/f
 array(13:16) = state%JKmuS(2:5)/f
 array(17) = state%tau_lucy
 array(18) = state%kappa_planck
 array(19) = state%kappa_ross
#else
 array(10)  = state%mu
#endif
 if (icooling > 0) array(nwrite) = state%Q
end subroutine state_to_array

subroutine filewrite_state(iunit, state)
 use wind,     only:wind_state
 integer, intent(in) :: iunit
 type(wind_state), intent(in) :: state

 real :: array(nwrite)

 call state_to_array(state, array)
 write(iunit, '(20E20.10E3)') array(1:nwrite)
end subroutine filewrite_state

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use dim,          only: maxvxyzu
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_velocity_km_s,'wind_velocity','injection wind velocity (km/s, if sonic_type = 0)',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period','stellar pulsation period (days)',iunit)
 call write_inopt(piston_velocity_km_s,'piston_velocity','velocity amplitude of the pulsation (km/s)',iunit)
 call write_inopt(wind_injection_radius_au,'wind_inject_radius','wind injection radius (au)',iunit)
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass loss rate (Msun/yr)',iunit)
 if (maxvxyzu==4) then
    call write_inopt(wind_temperature,'wind_temperature','wind temperature at the injection point (K)',iunit)
 endif
 call write_inopt(iwind_resolution,'iwind_resolution','if<>0 set number of particles on the sphere, reset particle mass',iunit)
 call write_inopt(wind_shell_spacing,'wind_shell_spacing','desired ratio of sphere spacing to particle spacing',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres','number of boundary spheres (integer)',iunit)
 call write_inopt(sonic_type,'sonic_type','find transonic solution (1=yes,0=no)',iunit)
 call write_inopt(outer_boundary_au,'outer_boundary','kill gas particles outside this radius (au)',iunit)
end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
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
 case('outer_boundary')
    read(valstring,*,iostat=ierr) outer_boundary_au
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
 case('sonic_type')
    read(valstring,*,iostat=ierr) sonic_type
    ngot = ngot + 1
    if (sonic_type < 0 .or. sonic_type > 1 ) call fatal(label,'invalid setting for sonic_type ([0,1])')
 case('wind_mass_rate')
    read(valstring,*,iostat=ierr) wind_mass_rate_Msun_yr
    ngot = ngot + 1
    if (wind_mass_rate_Msun_yr < 0.) call fatal(label,'invalid setting for wind_mass_rate (<0)')
 case('pulsation_period')
    read(valstring,*,iostat=ierr) pulsation_period_days
    ngot = ngot + 1
    if (pulsation_period_days < 0.) call fatal(label,'invalid setting for pulsation_period (<0)')
 case('piston_velocity')
    read(valstring,*,iostat=ierr) piston_velocity_km_s
    !wind_velocity_km_s = 0. ! set wind veolicty to zero when pulsating star
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 noptions = 9
#ifdef NUCLEATION
 noptions = 10
#elif ISOTHERMAL
 noptions = 7
#endif
 !debug
 !print '(a26,i3,i3)',trim(name),ngot,noptions
 igotall = (ngot >= noptions)
 if (trim(name) == '') ngot = 0
end subroutine read_options_inject

end module inject
