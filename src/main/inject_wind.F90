!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Handles wind injection
!
!  REFERENCES: None
!
!  OWNER: Lionel Siess
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    bowen_Cprime       -- radiative cooling rate (g.s/cm³)
!    star_Lum           -- central star luminosity (Lsun)
!    bowen_Tcond        -- dust condensation temperature (K)
!    star_Teff          -- central star effective temperature (K)
!    bowen_delta        -- condensation temperature range (K)
!    bowen_kappa        -- constant gas opacity (cm²/g)
!    bowen_kmax         -- maximum dust opacity (cm²/g)
!    iboundary_spheres  -- number of boundary spheres (integer)
!    iwind_resolution   -- if<>0 set number of particles on the sphere, reset particle mass
!    shift_spheres      -- delay before the ejection of shells
!    wind_CO_ratio      -- wind initial C/O ratio
!    wind_alpha         -- fraction of the gravitational acceleration imparted to the gas
!    wind_shell_spacing -- desired ratio of sphere spacing to particle spacing
!    wind_expT          -- temperature law exponent (if wind_type=1)
!    wind_inject_radius -- radius of injection of the wind (au)
!    wind_mass_rate     -- wind mass loss rate (Msun/yr)
!    wind_osc_period    -- stellar pulsation period (days)
!    wind_temperature   -- wind temperature at the injection point (K)
!    wind_type          -- stellar wind (1 = no dust, 2 = T(r)+dust, 3 = adia+dust, 4 = 3+cooling)
!    wind_velocity      -- velocity at which wind is injected (km/s)
!
!  DEPENDENCIES: bowen_dust, dim, eos, icosahedron, infile_utils,
!    injectutils, io, part, physcon, units
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
#ifndef KROME
 real, public::    bowen_kmax = 2.7991
 real, public::    bowen_Tcond = 1500.
 real, public::    bowen_delta = 60.
 real, public::    bowen_kappa = 2.d-4
 real, public::    bowen_Cprime = 3.000d-5
 real, public::    wind_CO_ratio = 2
 real, public::    wind_expT = 0.5
 integer, public:: wind_type = 4
#endif
#ifdef BOWEN !pulsating wind
 real, public::    star_Lum = 5315. * solarl
 real, public::    bowen_Cprime = 1.000d-5
 real, public::    wind_osc_period_days = 350.
 real, public::    wind_osc_vamplitude_km_s = 3.
 real, public::    wind_velocity_km_s = 0.
 real, public::    wind_mass_rate_Msun_yr = 1.04d-7
 real, public::    wind_temperature = 3000.
 real, public::    wind_injection_radius_au = 1.2568
#else
 real, public::    wind_velocity_km_s = 35.
 real, public::    wind_mass_rate_Msun_yr = 1.00d-8
 real, public::    wind_temperature = 3000.
 real, public::    wind_injection_radius_au = 1.7
#endif
 integer, public:: iwind_resolution = 0
 integer, public:: iboundary_spheres = 3
 real, public ::   wind_shell_spacing = 1.
 real, public ::   wind_alpha = 0.
 real, public::    shift_spheres = 2.

! Calculated from the previous parameters
 real, public ::    dr3
 real, public ::    wind_osc_vamplitude,wind_osc_period,mass_of_particles

 private

 real, private ::  wind_mass_rate,wind_velocity,mass_of_spheres,time_between_spheres,neighbour_distance,&
  Rstar_cgs,Cprime,wind_injection_radius
 integer, private :: particles_per_sphere,nwall_particles,iresolution

 logical, parameter :: wind_verbose = .false.
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3), u_to_temperature_ratio, rho_ini

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use physcon,     only:Rg, days, km, au, years, solarm,pi
 use icosahedron, only:compute_matrices, compute_corners
 use eos,         only:gmw, gamma
 use units,       only:unit_velocity, umass, utime, udist
 use part,        only:massoftype,igas,iboundary
 use io,          only:iverbose
 use injectutils, only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 integer, intent(out) :: ierr
 real :: mV_on_MdotR
 real :: dr,dp,mass_of_particles1,wind_velocity_max

 !
 ! return without error
 !
 ierr = 0
 !
 ! convert input parameters to code units
 !
#ifdef BOWEN
 wind_osc_period        = wind_osc_period_days * (days/utime)
 wind_osc_vamplitude    = wind_osc_vamplitude_km_s * (km / unit_velocity)
#else
 wind_osc_vamplitude    = 0.d0
#endif
 wind_injection_radius  = wind_injection_radius_au * au / udist
 wind_velocity          = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate         = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif
 wind_velocity_max      = max(wind_velocity,wind_osc_vamplitude)

 if (iresolution == 0) then
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
    iresolution = iwind_resolution
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    mass_of_particles = wind_shell_spacing*neighbour_distance * wind_injection_radius * wind_mass_rate &
       / (particles_per_sphere * wind_velocity)
    massoftype(igas) = mass_of_particles
    print *,'iwind_resolution unchanged = ',iresolution
 endif
 mass_of_spheres = mass_of_particles * particles_per_sphere
 nwall_particles = iboundary_spheres*particles_per_sphere
 time_between_spheres  = mass_of_spheres / wind_mass_rate
 massoftype(iboundary) = mass_of_particles

 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 if (iverbose >= 1) then
    print*,'mass of particle = ',massoftype(igas)
    mass_of_particles1 = wind_shell_spacing * get_neighb_distance(4) * wind_injection_radius * wind_mass_rate &
                     / (get_parts_per_sphere(4) * wind_velocity_max)
    print*,'require mass of particle = ',mass_of_particles1,' to get 492 particles per sphere'
    print*,'particles_per_sphere ',particles_per_sphere
    print*,'neighbour_distance ',neighbour_distance
    dp = neighbour_distance*wind_injection_radius
    dr = wind_velocity_max*mass_of_particles*particles_per_sphere/wind_mass_rate
    print*,'particle separation on spheres = ',dp
    print*,'distance between spheres = ',dr
    print*,'got dr/dp = ',dr/dp,' compared to desired dr on dp = ',wind_shell_spacing
 endif

#ifdef BOWEN
 call init_bowen(u_to_temperature_ratio,bowen_kappa,bowen_kmax,star_Lum,wind_injection_radius,&
      bowen_Cprime,bowen_Tcond,bowen_delta,star_Teff,wind_osc_vamplitude,wind_osc_period,&
      iboundary_spheres*particles_per_sphere)
 rho_ini = wind_mass_rate / (4.*pi*wind_injection_radius**2*wind_velocity_max)
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

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use physcon,     only: pi
 use io,          only:fatal
 use eos,         only:gamma
 use units,       only:udist
 use part,        only:igas,iboundary,nptmass
 use injectutils, only:inject_geodesic_sphere
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere, inner_sphere, inner_boundary_sphere, first_particle, i, ierr, ipart
 real    :: local_time, GM, r, v, u, rho, e, mass_lost, surface_radius, x0(3), v0(3)!, cs2max, dr, dp
 !logical, save :: first_run = .true.
 character(len=*), parameter :: label = 'inject_particles'

 ! Already done in initial.F90
 ! if (first_run) then
 !    call init_inject(ierr)
 !    first_run = .false.
 ! endif

 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
    GM = xyzmh_ptmass(4,wind_emitting_sink)
    v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 else
    x0 = 0.
    v0 = 0.
    GM = 0.
 endif

 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1
 inner_sphere = floor(time/time_between_spheres)
 inner_boundary_sphere = inner_sphere + iboundary_spheres

 if (inner_sphere-outer_sphere > iboundary_spheres) call fatal(label,'problem with boundary spheres, timestep likely too large!')
! cs2max = 0.
  !TESTS
 if (shift_spheres < 0) then
   ipart = iboundary
 else
   ipart = igas
 endif
 do i=inner_boundary_sphere,outer_sphere,-1
    local_time = time - (i-abs(shift_spheres)) * time_between_spheres
    call compute_sphere_properties(time,local_time,gamma,GM,r,v,u,rho,e,i,&
         inner_sphere,inner_boundary_sphere)
    if (wind_verbose) then
    print '("sphere ",5(i3),9(1x,es15.8))',i,inner_sphere,iboundary_spheres,outer_sphere,int(shift_spheres),time,local_time,&
    r/xyzmh_ptmass(5,1),v,v/sqrt(gamma*(gamma-1.)*u),xyzmh_ptmass(5,1)*udist,dtlast
    endif

    if (i > inner_sphere) then
       ! boundary sphere
       first_particle = (iboundary_spheres-i+inner_sphere)*particles_per_sphere+1
       call inject_geodesic_sphere(i, first_particle, iresolution, r, v, u, rho,  geodesic_R, geodesic_V, &
            npart, npartoftype, xyzh, vxyzu, ipart, x0, v0)
    else
      ! ejected particles
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
 ! update the sink particle properties
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
 dtinject = 0.2*time_between_spheres

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Compute the radius, velocity and temperature of a sphere at the current local time
!+
!-----------------------------------------------------------------------
subroutine compute_sphere_properties(time,local_time,gamma,GM,r,v,u,rho,e,sphere_number, &
                                     inner_sphere,inner_boundary_sphere)

!in/out variables in code units (except Jstar,K,mu)
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,gamma,GM
 real,    intent(out) :: r, v, u, rho, e

#ifdef BOWEN
 call bowen_wind_profile(time,local_time,r,v,u,rho,e,GM, gamma,sphere_number, &
                                     inner_sphere,inner_boundary_sphere)
#elif PARKER
 call parker_wind_profile(time,local_time,r,v,u,rho,e,GM, gamma,Jstar,K,mu,cs)
#else
 call stationary_wind_profile(local_time, r, v, u, rho, e, gamma, GM)
#endif

end subroutine compute_sphere_properties


#ifdef BOWEN
!-----------------------------------------------------------------------
!+
!  Oscillating inner boundary : bowen wind
!+
!-----------------------------------------------------------------------
subroutine bowen_wind_profile(time,local_time,r,v,u,rho,e,GM, gamma,sphere_number, &
                                     inner_sphere,inner_boundary_sphere)
 use physcon,     only: pi
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,gamma,GM
 real,    intent(out) :: r, v, u, rho, e

 integer, parameter :: nrho_index = 10
 integer :: k
 real :: surface_radius,r3
 logical :: verbose = .true.

 v = wind_velocity + wind_osc_vamplitude* cos(2.*pi*time/wind_osc_period) !same velocity for all wall particles
 surface_radius = wind_injection_radius + wind_osc_vamplitude*wind_osc_period/(2.*pi)*sin(2.*pi*time/wind_osc_period)
 if (sphere_number <= inner_sphere) then
    r = surface_radius
    v = max(wind_osc_vamplitude,wind_velocity)
 else
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
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)', &
             sphere_number,inner_sphere,surface_radius,r,v,&
             time/wind_osc_period,time_between_spheres/wind_osc_period
    else
       print '("ejected, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)', &
             sphere_number,inner_sphere,surface_radius,r,v,&
             time/wind_osc_period,time_between_spheres/wind_osc_period
    endif
 endif

end subroutine bowen_wind_profile
#endif

!-----------------------------------------------------------------------
!+
!  stationary supersonic wind
!+
!-----------------------------------------------------------------------
subroutine stationary_wind_profile(local_time, r, v, u, rho, e, gamma, GM)
 use physcon,     only: pi
 real, intent(in)  :: local_time, GM, gamma
 real, intent(out) :: r, v, u, rho, e
 real :: dt, rv(2), k1(2), k2(2), k3(2), k4(2), T, new_rv(2)
 integer, parameter :: N = 10000
 integer :: i

 dt = local_time / N
 rv(1) = wind_injection_radius
 rv(2) = wind_velocity
 ! Runge-Kutta iterations
 do i=1,N
    call drv_dt(rv,          k1, GM, gamma)
    new_rv = rv+dt/2.*k1
    call drv_dt(new_rv, k2, GM, gamma)
    new_rv = rv+dt/2.*k2
    call drv_dt(new_rv, k3, GM, gamma)
    new_rv = rv+dt/1.*k3
    call drv_dt(new_rv,    k4, GM, gamma)
    rv = rv + dt/6. * (k1 + 2.*k2 + 2.*k3 + k4)
 enddo
 r = rv(1)
 v = rv(2)
 ! this expression for T is only valid for an adiabatic EOS !
 if (gamma > 1.0001) then
    T = wind_temperature * (wind_injection_radius**2 * wind_velocity / (r**2 * v))**(gamma-1.)
    u = T * u_to_temperature_ratio
    e = .5*v**2 - GM/r + gamma*u
 else
    u = T * u_to_temperature_ratio
    e = .5*v**2 - GM/r + u
 endif
 rho = wind_mass_rate / (4.*pi*r**2*v)

end subroutine stationary_wind_profile

!-----------------------------------------------------------------------
!+
!  Time derivative of r and v, for Runge-Kutta iterations (stationary wind solution)
!+
!-----------------------------------------------------------------------
subroutine drv_dt(rv,drv,GM,gamma)
 use eos, only : polyk
 real, intent(in) :: rv(2),GM,gamma
 real, intent(out) :: drv(2)
 real :: r, v, dr_dt, r2, T, vs2, dv_dr, dv_dt, u

 r = rv(1)
 v = rv(2)
 dr_dt = v
 r2 = r*r
 if (gamma > 1.0001) then
    T = wind_temperature * (wind_injection_radius**2 * wind_velocity / (r2 * v))**(gamma-1.)
    u = T * u_to_temperature_ratio
    vs2 = gamma * (gamma - 1.) * u
 else
    vs2 = polyk
 endif
 dv_dr = (-GM/r2+2.*vs2/r)/(v-vs2/v)
 dv_dt = dv_dr * v

 drv(1) = dr_dt
 drv(2) = dv_dt

end subroutine drv_dt

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
 use units,        only: unit_velocity
 use physcon,      only: mass_proton_cgs, kboltz
 use eos,          only: gmw,polyk
 integer, intent(in) :: iunit

 call write_inopt(wind_velocity_km_s,'wind_velocity','velocity at which wind is injected (km/s)',iunit)
#if defined (BOWEN)
 call write_inopt(wind_osc_period_days,'wind_osc_period','stellar pulsation period (days)',iunit)
 call write_inopt(wind_osc_vamplitude_km_s,'wind_osc_vampl','velocity amplitude of the pulsations (km/s)',iunit)
#endif
 call write_inopt(wind_injection_radius_au,'wind_inject_radius','radius of injection of the wind (au)',iunit)
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass loss rate (Msun/yr)',iunit)
 if (maxvxyzu==4) then
    call write_inopt(wind_temperature,'wind_temperature','wind temperature at the injection point (K)',iunit)
 else
    wind_temperature = polyk* mass_proton_cgs/kboltz * unit_velocity**2*gmw
 endif
 call write_inopt(wind_alpha,'wind_alpha','fraction of the gravitational acceleration imparted to the gas',iunit)
 call write_inopt(iwind_resolution,'iwind_resolution','if<>0 set number of particles on the sphere, reset particle mass',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres','number of boundary spheres (integer)',iunit)
 call write_inopt(wind_shell_spacing,'wind_shell_spacing','desired ratio of sphere spacing to particle spacing',iunit)
 call write_inopt(shift_spheres,'shift_spheres','delay before the ejection of shells',iunit)
#if defined (BOWEN) || defined(PARKER)
 write(iunit,"(/,a)") '# options controlling dust and cooling'

 call write_inopt(star_Lum/solarl,'star_Lum','central star luminosity (Lsun)',iunit)
 call write_inopt(star_Teff,'star_Teff','central star effective temperature (K)',iunit)
 call write_inopt(bowen_kappa,'bowen_kappa','constant gas opacity (cm²/g)',iunit)
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
 call write_inopt(bowen_kmax,'bowen_kmax','maximum dust opacity (cm²/g)',iunit)
 call write_inopt(bowen_Tcond,'bowen_Tcond','dust condensation temperature (K)',iunit)
 call write_inopt(bowen_delta,'bowen_delta','condensation temperature range (K)',iunit)
 call write_inopt(wind_type,'wind_type','stellar wind (1 = no dust, 2 = T(r)+dust, 3 = adia+dust, 4 = 3+cooling)',iunit)
 call write_inopt(wind_expT,'wind_expT','temperature law exponent (if wind_type=1)', iunit)
 call write_inopt(wind_CO_ratio ,'wind_CO_ratio','wind initial C/O ratio',iunit)
#endif

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only: fatal, error, warning
 use physcon, only: pi,steboltz,au
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
    if (wind_velocity < 0.)    call fatal(label,'invalid setting for wind_velocity (<0)')
 case('wind_mass_rate')
    read(valstring,*,iostat=ierr) wind_mass_rate_Msun_yr
    ngot = ngot + 1
    if (wind_mass_rate < 0.)    call fatal(label,'invalid setting for wind_mass_rate (<0)')
 case('wind_temperature')
    read(valstring,*,iostat=ierr) wind_temperature
    ngot = ngot + 1
    isowind = .false.
    if (wind_temperature < 0.)    call fatal(label,'invalid setting for wind_temperature (<0)')
 case('wind_alpha')
     read(valstring,*,iostat=ierr) wind_alpha
     if (wind_alpha < 0.) call fatal(label,'invalid setting for wind_alpha (must be > 0)')
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
 case('wind_inject_radius')
    read(valstring,*,iostat=ierr) wind_injection_radius_au
    ngot = ngot + 1
    if (wind_injection_radius_au < 0.) call fatal(label,'invalid setting for wind_inject_radius (<0)')
#if defined (BOWEN) || defined(PARKER)
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
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
    if (bowen_Cprime < 0.) call fatal(label,'invalid setting for bowen_Cprime (<0)')
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
 case('wind_type')
    read(valstring,*,iostat=ierr) wind_type
    if (wind_type < 0 .or. wind_type > 4 ) call fatal(label,'invalid setting for wind_type ([0,3])')
 case('wind_CO_ratio')
    read(valstring,*,iostat=ierr) wind_CO_ratio
    if (wind_CO_ratio < 0.) call fatal(label,'invalid setting for wind_CO_ratio (must be > 0)')
 case('wind_expT')
    read(valstring,*,iostat=ierr) wind_expT
#endif
 case default
    imatch = .false.
 end select
#ifdef BOWEN
 noptions = 16
#elseif PARKER
 noptions = 14
#else
 noptions = 6
#endif
#if defined (BOWEN) || defined(PARKER)
 Rstar_cgs = sqrt(star_Lum/(4.*pi*steboltz*star_Teff**4))
 !if you launch the wind inside the photosphere, make sure the wind temperature >= star's effective temperature
 if (wind_injection_radius_au <= Rstar_cgs/au .and. wind_temperature < star_Teff) then
    call fatal(label,'invalid setting for wind_temperature (< star_Teff)')
 endif
#endif
 if (isowind) noptions = noptions -1
 igotall = (ngot >= noptions)
 if (iboundary_spheres.gt.int(shift_spheres).and.igotall) then
    print *,'too many boundary shells - imposing iboundary_spheres = shift_spheres',iboundary_spheres,shift_spheres
    iboundary_spheres = int(shift_spheres)
 endif
end subroutine read_options_inject

end module inject
