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
!    nfill_domain       -- number of spheres used to set the background density profile
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
 private
!
!--runtime settings for this module
!
! Read from input file
#ifdef NUCLEATION
 integer:: sonic_type = 1
 real::    wind_velocity_km_s = 0.
 real::    wind_mass_rate_Msun_yr = 1.d-5
 real::    wind_injection_radius_au = 2.37686663
 real::    wind_temperature = 2500.
#elif ISOTHERMAL
 integer:: sonic_type = 1
 real::    wind_velocity_km_s = 25.
 real::    wind_mass_rate_Msun_yr = 1.d-8
 real::    wind_injection_radius_au = 0.46524726
 real::    wind_temperature
#else
 integer:: sonic_type = 0
 real::    wind_velocity_km_s = 35.
 real::    wind_mass_rate_Msun_yr = 1.d-8
 real::    wind_injection_radius_au = 1.7
 real::    wind_temperature = 3000.
#endif
 integer :: iboundary_spheres = 5
 integer :: iwind_resolution = 0
 integer :: nfill_domain = 30
 real :: outer_boundary_au = 10.
 real :: wind_shell_spacing = 1.
 real :: pulsation_period
 real :: pulsation_period_days = 0.
 real :: piston_velocity_km_s = 0.
 real :: dtpulsation = 1.d99


! global variables
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: u_to_temperature_ratio,wind_mass_rate,piston_velocity,wind_velocity,&
      mass_of_spheres,time_between_spheres,neighbour_distance,mass_of_particles,&
      dr3,Rstar_cgs,Rinject,wind_injection_radius,wind_injection_speed,rho_ini,&
      omega_osc,deltaR_osc,Mstar_cgs
 integer :: particles_per_sphere,nwall_particles,iresolution,nwrite

 logical :: pulsating_wind
 character(len=*), parameter :: label = 'inject_wind'

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use options,        only:icooling,ieos
 use io,             only:fatal,iverbose
 use timestep,       only:tmax,dtmax
 use wind_equations, only:init_wind_equations
 use wind,           only:setup_wind,save_windprofile
 use physcon,        only:mass_proton_cgs, kboltz, Rg, days, km, au, years, solarm, pi, Gg
 use icosahedron,    only:compute_matrices, compute_corners
 use eos,            only:gmw,gamma,polyk
 use units,          only:unit_velocity, umass, utime, udist
 use part,           only:xyzmh_ptmass,massoftype,igas,iboundary,xyzmh_ptmass,imloss,ilum,iTeff,iReff
 use injectutils,    only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 integer, intent(out) :: ierr
 integer :: ires_min,nzones_per_sonic_point
 real :: mV_on_MdotR,initial_wind_velocity_cgs,dist_to_sonic_point
 real :: dr,dp,mass_of_particles1,tcross,tend,vesc,rsonic,tsonic,initial_Rinject


 if (icooling > 0) nwrite = nwrite+1
 ierr = 0

 pulsating_wind = (pulsation_period_days > 0.) .and. (piston_velocity_km_s > 0.)
 if (pulsating_wind .and. ieos == 6) call fatal(label,'cannot use ieos=6 with pulsation')
 !
 ! convert input parameters to code units
 !
 Rstar_cgs        = xyzmh_ptmass(iReff,wind_emitting_sink)*udist
 Mstar_cgs        = xyzmh_ptmass(4,wind_emitting_sink)*umass
 wind_velocity    = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate   = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 if (wind_injection_radius_au == 0.)  wind_injection_radius_au = Rstar_cgs/au
 wind_injection_radius = wind_injection_radius_au * au / udist
 if (pulsating_wind) then
    pulsation_period = pulsation_period_days * (days/utime)
    piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
    dtpulsation      = pulsation_period/50.
    omega_osc        = 2.*pi/pulsation_period
    deltaR_osc       = pulsation_period*piston_velocity/(2.*pi)
    sonic_type       = 1
 else
    deltaR_osc       = 0.d0
    piston_velocity  = 0.d0
 endif
 initial_wind_velocity_cgs = (piston_velocity+wind_velocity)*unit_velocity
 if (initial_wind_velocity_cgs == 0. .and. sonic_type == 0) call fatal(label,'zero input wind velocity')
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif

#ifdef ISOTHERMAL
 !Rstar_cgs = wind_injection_radius_au *au
 wind_temperature = polyk* mass_proton_cgs/kboltz * unit_velocity**2*gmw
#endif
 if (wind_injection_radius_au < Rstar_cgs/au) then
    print *,'WARNING wind_inject_radius < Rstar (au)',wind_injection_radius_au,Rstar_cgs/au
    !call fatal(label,'WARNING wind_inject_radius (< Rstar)')
 endif

 initial_Rinject = wind_injection_radius
 if ( .not. pulsating_wind .or. nfill_domain > 0) then
    if (pulsating_wind) then
       !implement background spheres starting from the smallest radius
       !initial_Rinject = min(initial_Rinject,xyzmh_ptmass(iReff,wind_emitting_sink)-deltaR_osc)
    endif
    if (initial_Rinject < xyzmh_ptmass(5,wind_emitting_sink)) then
       print *,'STOP wind_inject_radius < Racc (au)',wind_injection_radius_au,xyzmh_ptmass(5,wind_emitting_sink)
       call fatal(label,'invalid setting wind_inject_radius < accretion radius')
    endif
    call init_wind_equations (xyzmh_ptmass(4,wind_emitting_sink), &
         xyzmh_ptmass(iTeff,wind_emitting_sink), u_to_temperature_ratio)
!if (ieos == 6 .and. wind_temperature /= star_Teff) then
!    wind_injection_radius_au = Rstar_cgs/au*(star_Teff/wind_temperature)**(1./wind_expT)
!    wind_injection_radius  = wind_injection_radius_au * au / udist
!endif
!wind_injection_radius = max(wind_injection_radius_au * au, Rstar_cgs) / udist
!Rstar = min(wind_injection_radius_au*au,Rstar_cgs)

! integrate wind equation to get initial velocity and radius required to set the resolution
    Rinject = initial_Rinject*udist
    call setup_wind(Mstar_cgs, wind_mass_rate, u_to_temperature_ratio, Rinject,&
         wind_temperature,initial_wind_velocity_cgs,rsonic,tsonic,sonic_type)
    initial_Rinject = Rinject/udist
 endif
 Rinject = initial_Rinject
 wind_injection_speed = initial_wind_velocity_cgs/unit_velocity
 rho_ini = wind_mass_rate/(4.*pi*Rinject**2*wind_injection_speed)

 if (iwind_resolution == 0) then
    !
    ! resolution is specified in terms of number of smoothing lengths
    ! per distance to sonic point
    !
    nzones_per_sonic_point = 8
    dist_to_sonic_point = rsonic/udist-Rinject
    dr = abs(dist_to_sonic_point)/nzones_per_sonic_point
    mass_of_particles = rho_ini*dr**3
    massoftype(igas) = mass_of_particles
    print*,' suggesting ',mass_of_particles, ' based on desired dr = ',dr,' dist-to-sonic=',dist_to_sonic_point
    !
    ! compute the dimensionless resolution factor m V / (Mdot R)
    ! where m = particle mass and V, Mdot and R are wind parameters
    !
!    mass_of_particles = massoftype(igas)
    mV_on_MdotR = mass_of_particles*wind_injection_speed/(wind_mass_rate*Rinject)
    !
    ! solve for the integer resolution of the geodesic spheres
    ! gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
    !
    iresolution = get_sphere_resolution(wind_shell_spacing,mV_on_MdotR)
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    print *,'iwind_resolution equivalence = ',iresolution
 else
    iresolution = iwind_resolution
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    mass_of_particles = wind_shell_spacing*neighbour_distance*Rinject*wind_mass_rate &
         / (particles_per_sphere * wind_injection_speed)
    massoftype(igas) = mass_of_particles
    print *,'iwind_resolution unchanged = ',iresolution
 endif

 mass_of_spheres = mass_of_particles * particles_per_sphere
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

!compute full evolution (to get tcross) and save 1D profile for comparison
 if ( .not. pulsating_wind .or. nfill_domain > 0) then
    tend = max(tmax,(iboundary_spheres+nfill_domain)*time_between_spheres)*utime
    call save_windprofile(Rinject*udist,wind_injection_speed*unit_velocity,&
         wind_temperature, outer_boundary_au*au, tsonic, tend, tcross, 'windprofile1D.dat')
    if ((iboundary_spheres+nfill_domain)*time_between_spheres > tmax) then
       print *,'simulation time < time to reach the last boundary shell'
    endif
    if (tcross < 1.d98) then
       nfill_domain = min(nfill_domain,int(tcross/time_between_spheres)-iboundary_spheres)
       print *,'reduce number of background shells to',nfill_domain
    endif
 endif

 if (iverbose >= 1) then
    mass_of_particles1 = wind_shell_spacing * get_neighb_distance(4) * Rinject * wind_mass_rate &
         / (get_parts_per_sphere(4) * wind_injection_speed)
    print*,'required mass of particle  = ',mass_of_particles1,' to get 492 particles per sphere'
    dp = neighbour_distance*Rinject
    dr = wind_injection_speed*mass_of_particles*particles_per_sphere/wind_mass_rate
    print*,'particle separation on spheres = ',dp
    print*,'got dr/dp = ',dr/dp,' compared to desired dr on dp = ',wind_shell_spacing
 endif


!logging
 vesc = sqrt(2.*Gg*Mstar_cgs/Rstar_cgs)
 print*,'mass_of_particles          = ',mass_of_particles
 print*,'particles per sphere       = ',particles_per_sphere
 print*,'distance between spheres   = ',wind_shell_spacing*neighbour_distance
 print*,'distance between injection = ',time_between_spheres*wind_injection_speed
 print*,'hmax/dist_between_spheres  = ',wind_shell_spacing*neighbour_distance*initial_wind_velocity_cgs**2/&
      (vesc**2-initial_wind_velocity_cgs**2)
 if (sonic_type == 1) then
    print*,'distance to sonic point    = ',rsonic/udist-Rinject
    print*,'sonic radius               = ',rsonic/udist,rsonic
    print*,'number of shells to sonic  = ',(rsonic/udist-Rinject)/(wind_shell_spacing*neighbour_distance)
    print*,'time_to_sonic_point        = ',tsonic/utime
 endif
 print*,'time_between_spheres       = ',time_between_spheres
 print*,'wind_temperature           = ',wind_temperature
 print*,'injection_radius           = ',Rinject
 print*,'stellar_radius             = ',Rstar_cgs / udist
 print*,'rho_ini                    = ',rho_ini
 if (pulsating_wind) then
    print*,'number of ejected shells per pulsation period (should at least be > 10) ',pulsation_period/time_between_spheres
    print*,'width of the boundary layer/ R* (should be < 1) = ',1.-(iboundary_spheres*dr3)**(1./3.)/Rinject
    print*,'radial pulsation amplitude/ R* = ',piston_velocity*pulsation_period/(2.*pi*Rinject),dr3/Rinject
    print*,'pulsation period in code units = ',pulsation_period
    !sanity checks
    ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
    if (pulsation_period/time_between_spheres < 10. ) print *,'WARNING! only ',pulsation_period/time_between_spheres,&
         ' shells will be ejected during a pulsation period'
    ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
    if (1.-(iboundary_spheres*dr3)**(1./3.)/Rinject > 0.2)  print*,'WARNING! the width of the boundary layer = ',&
         1.-(iboundary_spheres*dr3)**(1./3.)/Rinject,' Rinject'
 elseif (rsonic/udist > Rinject) then
    !save a few models before the particles reach the sonic point
    if (dtmax > tsonic/utime) print *,'WARNING! dtmax > time to sonic point'
    !if solution subsonic, minimum resolution required so a few shells can be inserted between the injection radius and the sonic point
    ires_min = int(iboundary_spheres*wind_shell_spacing*0.5257/(rsonic/udist/Rinject-1.)+.5)
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
 use wind,         only:wind_profile
 use part,         only:igas,iTeff,iReff,iboundary,nptmass,delete_particles_outside_sphere,&
      delete_dead_particles_inside_radius,dust_temp,n_nucleation,massoftype

 use partinject,   only:add_or_update_particle
 use injectutils,  only:inject_geodesic_sphere
 use units,        only:udist
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere, inner_sphere, inner_boundary_sphere, first_particle, i, ipart, &
            nreleased, nboundaries
 real    :: local_time, GM, r, v, u, rho, e, mass_lost, x0(3), v0(3), inner_radius, fdone
 character(len=*), parameter :: label = 'inject_particles'
 logical, save :: released = .false.
#ifdef NUCLEATION
 real :: JKmuS(n_nucleation)
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


 if (npart > 0) then
    !
    ! delete particles that exit the outer boundary
    !
    i = npart
    inner_radius = wind_injection_radius + deltaR_osc*sin(omega_osc*time)
    call delete_particles_outside_sphere(x0,outer_boundary_au*au/udist,npart)
    call delete_dead_particles_inside_radius(x0,inner_radius,npart)
    print *,npart,i,inner_radius,dtlast,dtlast*wind_mass_rate/massoftype(igas)
    if (npart.ne.i .and. iverbose > 0) print *,'deleted ',i-npart,'particles, remaining',npart

    nreleased = nfill_domain
    nboundaries = iboundary_spheres
    !release particles and declare inner boundary shells as gas particles so they can exert some pressure
    ipart = igas
    if (.not.released) then
       do i = max(1,npart-nreleased*particles_per_sphere)+1,npart
          call add_or_update_particle(igas,xyzh(1:3,i),vxyzu(1:3,i),xyzh(4,i),vxyzu(4,i),i,npart,npartoftype,xyzh,vxyzu)
       enddo
       released = .true.
    endif
 else
    !initialise domain with boundary particles
    ipart = iboundary
    nreleased = 0
    nboundaries = iboundary_spheres+nfill_domain
 endif
 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1 + nreleased
 inner_sphere = floor(time/time_between_spheres)+nreleased
 inner_boundary_sphere = inner_sphere + nboundaries

 !only one sphere can be ejected at a time
 if (inner_sphere-outer_sphere > nboundaries) call fatal(label,'ejection of more than 1 sphere, timestep likely too large!')

 do i=inner_boundary_sphere,outer_sphere,-1
    local_time = time + (iboundary_spheres+nfill_domain-i) * time_between_spheres

    !compute the radius, velocity, temperature, chemistry of a sphere at the current local time
    v = wind_injection_speed
    r = Rinject
    if (pulsating_wind.and.released) then
       call pulsating_wind_profile(time,local_time, r, v, u, rho, e, GM, i, &
            inner_sphere,inner_boundary_sphere,dr3,rho_ini)
    else
#ifdef NUCLEATION
       call wind_profile(local_time, r, v, u, rho, e, GM, wind_temperature, fdone, JKmuS)
#else
       call wind_profile(local_time, r, v, u, rho, e, GM, wind_temperature, fdone)
#endif
       if (iverbose > 0) print '(" ##### boundary sphere ",i4,3(i4),i7,9(1x,es12.5))',i,&
            inner_sphere,nboundaries,outer_sphere,npart,time,local_time,r,v,fdone
    endif

    if (i > inner_sphere) then
       ! boundary sphere
       first_particle = (nboundaries-i+inner_sphere)*particles_per_sphere+1
       !print '(" ##### boundary sphere ",i4,i7,3(i4),i7,9(1x,es12.5))',i,first_particle,inner_sphere,nboundaries,&
       !     outer_sphere,npart,time,local_time,r/xyzmh_ptmass(iReff,1),v,u,rho
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
            outer_sphere,npart,time,local_time,r/xyzmh_ptmass(iReff,1),v,u,rho
    endif
    !cs2max = max(cs2max,gamma*(gamma-1)*u)
 enddo
 ! update sink particle properties
 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_lost
    if (pulsating_wind) then
       inner_radius = wind_injection_radius + deltaR_osc*sin(omega_osc*time)
       !v2 xyzmh_ptmass(5,wind_emitting_sink) = (inner_radius**3-dr3)**(1./3.) !accretion radius
       xyzmh_ptmass(5,wind_emitting_sink) = inner_radius
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
 use physcon, only:pi
 use eos,     only:gamma
 use units,   only:unit_velocity
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,GM,dr3,rho_ini
 real,    intent(out) :: r, v, u, rho, e

 integer, parameter :: nrho_index = 10
 integer :: k
 real :: inner_radius,r3
 logical :: verbose = .true.

 v = wind_velocity + piston_velocity* cos(omega_osc*time) !same velocity for all wall particles
 inner_radius = wind_injection_radius + deltaR_osc*sin(omega_osc*time)
 !ejected spheres
 if (sphere_number <= inner_sphere) then
    r = inner_radius
    v = max(piston_velocity,wind_velocity)
 else
    !boundary spheres
    r3 = inner_radius**3-dr3
    do k = 2,sphere_number-inner_sphere
       r3 = r3-dr3*(r3/inner_radius**3)**(nrho_index/3.)
    enddo
    r = r3**(1./3)
 endif
 !r = (inner_radius**3-(sphere_number-inner_sphere)*dr3)**(1./3)
 !rho = rho_ini
 u = wind_temperature * u_to_temperature_ratio
 if (gamma > 1.0001) then
    e = .5*v**2 - GM/r + gamma*u
 else
    e = .5*v**2 - GM/r + u
 endif
 rho = rho_ini*(inner_radius/r)**nrho_index
 if (verbose) then
    if (sphere_number > inner_sphere) then
       print '("boundary, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4)', &
             sphere_number,inner_sphere,inner_radius,r,v*unit_velocity,time/pulsation_period
    else
       print '("ejected, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4)', &
             sphere_number,inner_sphere,inner_radius,r,v*unit_velocity,time/pulsation_period
    endif
 endif

end subroutine pulsating_wind_profile

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
 call write_inopt(wind_injection_radius_au,'wind_inject_radius','wind injection radius (au, if 0 take Rstar)',iunit)
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass loss rate (Msun/yr)',iunit)
 if (maxvxyzu==4) then
    call write_inopt(wind_temperature,'wind_temperature','wind temperature at the injection point (K)',iunit)
 endif
 call write_inopt(iwind_resolution,'iwind_resolution','if<>0 set number of particles on the sphere, reset particle mass',iunit)
 call write_inopt(nfill_domain,'nfill_domain','number of spheres used to set the background density profile',iunit)
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
 case('nfill_domain')
    read(valstring,*,iostat=ierr) nfill_domain
    ngot = ngot + 1
    if (nfill_domain < 0) call fatal(label,'nfill_domain must be > 0')
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
 noptions = 10
#ifdef NUCLEATION
 noptions = 11
#elif ISOTHERMAL
 noptions = 8
#endif
 !debug
 !print '(a26,i3,i3)',trim(name),ngot,noptions
 igotall = (ngot >= noptions)
 if (trim(name) == '') ngot = 0
end subroutine read_options_inject

end module inject
