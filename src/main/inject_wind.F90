!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Handles wind injection
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    bowen_Cprime      -- radiative cooling rate (g.s/cm³)
!    bowen_L           -- central star luminosity (Lsun)
!    bowen_Tcond       -- condensation temperature of dust (K)
!    bowen_Teff        -- effective temperature of the central star (K)
!    bowen_delta       -- condensation temperature range (K)
!    bowen_kappa       -- constant gas opacity (cm²/g)
!    bowen_kmax        -- maximum dust opacity (cm²/g)
!    iboundary_spheres -- number of boundary spheres (integer)
!    wind_dr_on_dp     -- desired ratio of sphere spacing to particle spacing
!    wind_mass_rate    -- wind mass per unit time (Msun/yr)
!    wind_osc_period   -- period of the oscillations (days)
!    wind_temperature  -- wind temperature at the injection point (K)
!
!  DEPENDENCIES: bowen_dust, eos, icosahedron, infile_utils, io, part,
!    partinject, physcon, timestep, units
!+
!--------------------------------------------------------------------------
module inject
 use physcon, only: au, solarm, years, solarl
 implicit none
 real, parameter :: phi = (sqrt(5.)+1.)/2. ! Golden ratio
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: inject_particles, write_options_inject, read_options_inject, wind_init
!
!--runtime settings for this module
!
! Read from input file
#ifdef BOWEN
 real, public ::    bowen_kappa = 2.d-4
 real, public ::    bowen_Teff = 3000.
 real, public ::    bowen_kmax = 2.7991
 real, public ::    bowen_Tcond = 1500.
 real, public ::    bowen_delta = 60.
 real, public ::    bowen_L = 5315. * solarl
 real, public ::    bowen_Cprime = 1.000d-5
 real, public ::    wind_osc_period_days = 350.
 real, public ::    wind_osc_vamplitude_km_s = 3.
 real, public ::    wind_velocity_km_s = 0.
 real, public ::    wind_mass_rate_Msun_yr = 1.04d-7
 real, public ::    wind_temperature = 3000.
! integer, public :: iwind_resolution = 5
! real, public ::    wind_sphdist = 0.2
 real, public ::    shift_spheres = 0.
 integer, public :: iboundary_spheres = 3
 real, public ::    wind_injection_radius = 1.2568 ! code units
#else
 real, public ::    wind_velocity_km_s = 35.
 real, public ::    wind_mass_rate_Msun_yr = 7.65d-7
 real, public ::    wind_temperature = 3000.
 real, public ::    shift_spheres = 2.
 integer, public :: iboundary_spheres = 3
 real, public ::    wind_injection_radius = 1.7 ! code units
 real, parameter :: wind_osc_vamplitude_km_s = 0.
#endif
 real, public ::    wind_dr_on_dp = 1.

! Calculated from the previous parameters
 real, public :: wind_osc_vamplitude,wind_osc_period
 real, private :: wind_mass_rate,wind_velocity,mass_of_spheres, time_between_spheres, neighbour_distance
 integer, private :: particles_per_sphere,iwind_resolution

 private

 logical, parameter :: wind_verbose = .false.
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3), u_to_temperature_ratio, dr3, rho_ini

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine wind_init(setup)
 use physcon,     only:Rg, pi, days, km
 use icosahedron, only:compute_matrices, compute_corners
 use timestep,    only:dtmax
 use eos,         only:gmw, gamma
 use units,       only:unit_velocity, umass, utime
 use part,        only:massoftype,igas
 use io,          only:iverbose
 !use part,        only:xyzmh_ptmass,nptmass
 logical, intent(in) :: setup
 real :: mV_on_MdotR,mass_of_particles
 real :: dr,dp,mass_of_particles1,wind_velocity_max
 real, parameter :: irrational_number_close_to_one = pi/3.
 !
 ! convert input parameters to code units
 !
#ifdef BOWEN
 wind_osc_period  = wind_osc_period_days * (days/utime)
#endif
 wind_osc_vamplitude  = wind_osc_vamplitude_km_s * (km / unit_velocity)
 wind_velocity  = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 wind_velocity_max = max(wind_velocity,wind_osc_vamplitude)
 !
 ! compute the dimensionless resolution factor m V / (Mdot R)
 ! where m = particle mass and V, Mdot and R are wind parameters
 !
 mass_of_particles = massoftype(igas)
 mV_on_MdotR = mass_of_particles*wind_velocity_max/(wind_mass_rate*wind_injection_radius)
 !
 ! solve for the integer resolution of the geodesic spheres
 ! gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
 !
 iwind_resolution     = get_wind_resolution(wind_dr_on_dp,mV_on_MdotR)
 particles_per_sphere = get_parts_per_sphere(iwind_resolution)
 neighbour_distance   = get_neighb_distance(iwind_resolution)

 mass_of_spheres = mass_of_particles * particles_per_sphere
 time_between_spheres = mass_of_spheres / wind_mass_rate

#ifdef BOWEN
 rho_ini = wind_mass_rate / (4.*pi*wind_injection_radius**2*wind_osc_vamplitude)
 dr3 = 3.*mass_of_spheres/(4.*pi*rho_ini)
#endif

 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 if (iverbose >= 1) then
    print*,'mass of particle = ',massoftype(igas)
    mass_of_particles1 = wind_dr_on_dp * get_neighb_distance(4) * wind_injection_radius * wind_mass_rate &
                     / (get_parts_per_sphere(4) * wind_velocity_max)
    print*,'require mass of particle = ',mass_of_particles1,' to get 492 particles per sphere'
    print*,'Mdot*R/(m*V) ',1./mV_on_MdotR
    print*,'wind_resolution ',iwind_resolution
    print*,'particles_per_sphere ',particles_per_sphere
    print*,'neighbour_distance ',neighbour_distance
    dp = neighbour_distance*wind_injection_radius
    dr = wind_velocity_max*mass_of_particles*particles_per_sphere/wind_mass_rate
    print*,'particle separation on spheres = ',dp
    print*,'distance between spheres = ',dr
    print*,'got dr/dp = ',dr/dp,' compared to desired dr on dp = ',wind_dr_on_dp
    print*,'mass_of_particles ',mass_of_particles
    print*,'mass_of_spheres ',mass_of_spheres
    print*,'time_between_spheres ',time_between_spheres
 endif

#ifdef BOWEN
 print*,'number of ejected shells per pulsation period (should at least be > 10) ',wind_osc_period/time_between_spheres
 print*,'width of the boundary layer/ R* (should be < 1) = ',1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius
 print*,'radial amplitude of the pulsations/ R* = ',wind_osc_vamplitude*wind_osc_period/(2.*pi)/wind_injection_radius
 print*,'particles_per_sphere ',particles_per_sphere
 !sanity checks
 ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
 if (wind_osc_period/time_between_spheres < 10. ) print *,'WARNING! only ',wind_osc_period/time_between_spheres,&
      ' shells will be ejected during a pulsation period'
 ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
 if (1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius > 0.2)  print*,'WARNING! the width of the boundary layer = ',&
      1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius,' Rinject'
#endif

 ! adjusting dtmax to avoid uterm < 0 errors
 ! to be removed when the routine provides timestep constraints
 if (setup) then
    dtmax = (.5 * irrational_number_close_to_one * time_between_spheres)
 endif

end subroutine wind_init

!-----------------------------------------------------------------------
!+
!  solve the cubic equation for q, the integer resolution of the spheres
!  given a desired ratio between the particle spacing in the sphere
!  and the radial spacing of spheres
!+
!-----------------------------------------------------------------------
integer function get_wind_resolution(dr_on_dp,mV_on_MdotR)
 real, intent(in) :: dr_on_dp,mV_on_MdotR
 real :: fac,term,f1,f2,res
 !
 ! we solve (20*q^3 - 30*q^2 + 16*q - 3) - 1./fac = 0
 ! the following is the only real solution,
 ! extracted from Mathematica
 !
 fac = mV_on_MdotR*sqrt(2.*(5.+sqrt(5.)))
 term = dr_on_dp/fac
 f1 = 15.**(1./3.)
 f2 = (45.*term + sqrt(15. + 2025.*term**2))**(1./3.)
 res = 0.5*(1. - 1./(f1*f2) + f2/f1**2)
 get_wind_resolution = nint(res) !max(nint(res(1)),4)

end function get_wind_resolution

!-----------------------------------------------------------------------
!+
!  Utility functions for geodesic sphere injection
!+
!-----------------------------------------------------------------------
integer function get_parts_per_sphere(iq)
 integer, intent(in) :: iq

 get_parts_per_sphere = 20 * (2*iq*(iq-1)) + 12

end function get_parts_per_sphere

real function get_neighb_distance(iq)
 integer, intent(in) :: iq

 !unitless: relative to the sphere radius
 get_neighb_distance = 2./((2.*iq-1.)*sqrt(sqrt(5.)*phi))

end function get_neighb_distance

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use io,      only: iprint, fatal
 use physcon, only: au
#ifdef BOWEN
 use bowen_dust, only: bowen_init
#endif
 use units,   only: unit_velocity, unit_density, udist
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: outer_sphere, inner_sphere, inner_handled_sphere, i
 real :: local_time, r, v, u, rho, e, mass_lost
 logical, save :: first_run = .true.
 character(len=*), parameter :: label = 'inject_particles'

 if (first_run) then
    call wind_init(.false.)
#ifdef BOWEN
    call bowen_init(u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,bowen_Cprime,&
         bowen_Tcond,bowen_delta,bowen_Teff,wind_osc_vamplitude,wind_osc_period)
#endif
    first_run = .false.
 endif

 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1
 inner_sphere = floor(time/time_between_spheres)
 inner_handled_sphere = inner_sphere + iboundary_spheres
 if (inner_sphere-outer_sphere > iboundary_spheres) call fatal(label,'problem with handled spheres, timestep likely too large!')
 do i=inner_sphere+iboundary_spheres,outer_sphere,-1
    local_time = time - (i-shift_spheres) * time_between_spheres
    call compute_sphere_properties(time,local_time, r, v, u, rho, e, i, inner_sphere, inner_handled_sphere, xyzmh_ptmass)
    if (wind_verbose) then
       write(iprint,*) '* Sphere:'
       write(iprint,*) '   Number: ', i,(i-shift_spheres)
       write(iprint,*) '   Local Time: ', local_time,time,(i-shift_spheres)*time_between_spheres
       write(iprint,*) '   Radius: ', r, '(', r*(udist/au), ' au)'
       write(iprint,*) '   Expansion velocity: ', v*unit_velocity, '(', v*unit_velocity/1.d5, ' km/s)'
       write(iprint,*) '   Temperature: ', u / u_to_temperature_ratio, ' K'
       write(iprint,*) '   Density: ', rho * unit_density, ' (g/cm³)'
       write(iprint,*) '   Constant e: ', e * unit_velocity**2 , ' (erg/g)'
       write(iprint,*) ''
       !read*
    endif
    if (i > inner_sphere) then
       call inject_geodesic_sphere(i, (iboundary_spheres-i+inner_sphere)*particles_per_sphere+1, r, v, u, rho, &
            npart, npartoftype, xyzh, vxyzu) ! handled sphere
    else
       call inject_geodesic_sphere(i, npart+1, r, v, u, rho, npart, npartoftype, xyzh, vxyzu) ! injected sphere
    endif
 enddo
! print *,'npart',npart,inner_sphere-outer_sphere+1
 mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
 xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_lost

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Time derivative of r and v, for Runge-Kutta iterations (stationary trans-sonic solution)
!+
!-----------------------------------------------------------------------
subroutine drv_dt(rv, drv, GM)
 use eos, only: gamma
 real, intent(in) :: rv(2), GM
 real, intent(out) :: drv(2)

 real :: r, v, dr_dt, r2, T, vs2, dv_dr, dv_dt, u

 r = rv(1)
 v = rv(2)
 dr_dt = v
 r2 = r*r
 T = wind_temperature * (wind_injection_radius**2 * wind_velocity / (r2 * v))**(gamma-1.)
 u = T * u_to_temperature_ratio
 vs2 = gamma * (gamma - 1) * u
 dv_dr = (-GM/r2+2.*vs2/r)/(v-vs2/v)
 dv_dt = dv_dr * v

 drv(1) = dr_dt
 drv(2) = dv_dt

end subroutine

!-----------------------------------------------------------------------
!+
!  Compute the radius, velocity and temperature of a sphere at the current local time
!+
!-----------------------------------------------------------------------
subroutine compute_sphere_properties(time,local_time, r, v, u, rho, e, sphere_number, &
                                     inner_sphere, inner_handled_sphere, xyzmh_ptmass)
#ifdef BOWEN
 !use bowen_dust, only: pulsating_bowen_wind_profile
 use eos,        only: gamma
 use physcon,    only: pi
#endif
 integer, intent(in) :: sphere_number, inner_sphere, inner_handled_sphere
 real, intent(in) :: time,local_time, xyzmh_ptmass(:,:)
 real, intent(out) :: r, v, u, rho, e
#ifdef BOWEN
 real :: surface_radius
 logical :: verbose = .true.

 v = wind_osc_vamplitude* cos(2.*pi*time/wind_osc_period) !same velocity for all wall particles
 surface_radius = wind_injection_radius + wind_osc_vamplitude*wind_osc_period/(2.*pi)*sin(2.*pi*time/wind_osc_period)
 r = (surface_radius**3-(sphere_number-inner_sphere)*dr3)**(1./3)
 u = wind_temperature * u_to_temperature_ratio
 rho = rho_ini
 e = .5*v**2 - xyzmh_ptmass(4,wind_emitting_sink)/r + gamma*u
 ! call pulsating_bowen_wind_profile(time,local_time, r, v, u, rho, e, sphere_number,&
 !       wind_mass_rate,wind_injection_radius,wind_velocity,wind_osc_vamplitude,&
 !       wind_osc_period,shift_spheres,central_star_mass,time_between_spheres,&
 !       wind_temperature,gamma)
 if (verbose) then
  if (sphere_number > inner_sphere) then
   print '("handled, i = ",i5," inner = ",i5," base_r = ",es10.4," r = ",es10.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)'&
          ,sphere_number,inner_sphere,surface_radius,r,v,time/wind_osc_period,time_between_spheres/wind_osc_period
  else
   print '("ejected, i = ",i5," inner = ",i5," base_r = ",es10.4," r = ",es10.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)'&
          ,sphere_number,inner_sphere,surface_radius,r,v,time/wind_osc_period,time_between_spheres/wind_osc_period
  endif
 endif
#else
 call stationary_adiabatic_wind_profile(local_time, r, v, u, rho, e)
#endif

end subroutine

!-----------------------------------------------------------------------
!+
!  stationary adiabatic supersonic wind
!+
!-----------------------------------------------------------------------
subroutine stationary_adiabatic_wind_profile(local_time, r, v, u, rho, e)
 use part,    only: nptmass, xyzmh_ptmass
 use physcon, only: pi
 use eos,     only: gamma
 real, intent(in) :: local_time
 real, intent(out) :: r, v, u, rho, e

 real :: GM
 real :: dt, rv(2), k1(2), k2(2), k3(2), k4(2), T
 integer, parameter :: N = 10000
 integer :: i

 dt = local_time / N
 rv(1) = wind_injection_radius
 rv(2) = wind_velocity
 if (nptmass == 0) then
    GM = 0.
 else
    GM = xyzmh_ptmass(4,wind_emitting_sink)
 endif
 ! Runge-Kutta iterations
 do i=1,N
    call drv_dt(rv,          k1, GM)
    call drv_dt(rv+dt/2.*k1, k2, GM)
    call drv_dt(rv+dt/2.*k2, k3, GM)
    call drv_dt(rv+dt*k3,    k4, GM)
    rv = rv + dt/6. * (k1 + 2.*k2 + 2.*k3 + k4)
 enddo
 r = rv(1)
 v = rv(2)
 ! this expression for T is only valid for an adiabatic EOS !
 T = wind_temperature * (wind_injection_radius**2 * wind_velocity / (r**2 * v))**(gamma-1.)
 u = T * u_to_temperature_ratio
 rho = wind_mass_rate / (4.*pi*r**2*v)
 e = .5*v**2 - GM/r + gamma*u

end subroutine stationary_adiabatic_wind_profile

!-----------------------------------------------------------------------
!+
!  Inject a quasi-spherical distribution of particles.
!+
!-----------------------------------------------------------------------
subroutine inject_geodesic_sphere(sphere_number, first_particle, r, v, u, rho, npart, npartoftype, xyzh, vxyzu)
 use icosahedron, only: pixel2vector
 use partinject,  only: add_or_update_particle
 use part,        only: igas, hrho, xyzmh_ptmass, vxyz_ptmass, nptmass
 integer, intent(in) :: sphere_number, first_particle
 real, intent(in) :: r, v, u, rho
 integer, intent(inout) :: npart, npartoftype(:)
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)

 real :: rotation_angles(3), h_sim
 real :: rotmat(3,3), radial_unit_vector(3), radial_unit_vector_rotated(3)
 real :: particle_position(3), particle_velocity(3)
 integer :: j

 select case (iwind_resolution)
 case(1)
    rotation_angles = (/ 1.28693610288783, 2.97863087745917, 1.03952835451832 /)
 case(2)
    rotation_angles = (/ 1.22718722289660, 2.58239466067315, 1.05360422660344 /)
 case(3)
    rotation_angles = (/ 0.235711384317490, 3.10477287368657, 2.20440220924383 /)
 case(4)
    rotation_angles = (/ 3.05231445647236, 0.397072776282339, 2.27500616856518 /)
 case(5)
    rotation_angles = (/ 0.137429597545199, 1.99860670500403, 1.71609391574493 /)
 case(6)
    rotation_angles = (/ 2.90443293496604, 1.77939686318657, 1.04113050588920 /)
 case(10)
    rotation_angles = (/ 2.40913070927068, 1.91721010369865, 0.899557511636617 /)
 case(15)
    rotation_angles = (/ 1.95605828396746, 0.110825898718538, 1.91174856362170 /)
 case default
    rotation_angles = (/ 1.28693610288783, 2.97863087745917, 1.03952835451832 /)
 end select
 rotation_angles = rotation_angles * sphere_number

 ! Quantities in simulation units
 h_sim = hrho(rho)

 call make_rotation_matrix(rotation_angles, rotmat)
 do j=0,particles_per_sphere-1
    call pixel2vector(j, iwind_resolution, geodesic_R, geodesic_v, radial_unit_vector)
    radial_unit_vector_rotated(1) = radial_unit_vector(1)*rotmat(1,1) &
                                  + radial_unit_vector(2)*rotmat(1,2) &
                                  + radial_unit_vector(3)*rotmat(1,3)
    radial_unit_vector_rotated(2) = radial_unit_vector(1)*rotmat(2,1) &
                                  + radial_unit_vector(2)*rotmat(2,2) &
                                  + radial_unit_vector(3)*rotmat(2,3)
    radial_unit_vector_rotated(3) = radial_unit_vector(1)*rotmat(3,1) &
                                  + radial_unit_vector(2)*rotmat(3,2) &
                                  + radial_unit_vector(3)*rotmat(3,3)
    particle_position = r*radial_unit_vector_rotated
    particle_velocity = v*radial_unit_vector_rotated
    if (nptmass > 0) then
       particle_position = particle_position + xyzmh_ptmass(1:3,wind_emitting_sink)
       particle_velocity = particle_velocity + vxyz_ptmass(1:3,wind_emitting_sink)
    endif
    call add_or_update_particle(igas, particle_position, particle_velocity, h_sim, u, first_particle+j, &
         npart,npartoftype,xyzh,vxyzu)
 enddo
end subroutine

!-----------------------------------------------------------------------
!+
!  Make a 3x3 rotation matrix from three angles.
!+
!-----------------------------------------------------------------------
subroutine make_rotation_matrix(rotation_angles, rot_m)
 real, intent(in)  :: rotation_angles(3)
 real, intent(out) :: rot_m(3,3)

 real :: angle_x, angle_y, angle_z
 real :: c_x, s_x, c_y, s_y, c_z, s_z

 angle_x = rotation_angles(1)
 angle_y = rotation_angles(2)
 angle_z = rotation_angles(3)

 c_x = cos(angle_x)
 s_x = sin(angle_x)
 c_y = cos(angle_y)
 s_y = sin(angle_y)
 c_z = cos(angle_z)
 s_z = sin(angle_z)

 rot_m(1,1) = c_y*c_z
 rot_m(1,2) = -c_y*s_z
 rot_m(1,3) = -s_y
 rot_m(2,1) = -s_x*s_y*c_z + c_x*s_z
 rot_m(2,2) = s_x*s_y*s_z + c_x*c_z
 rot_m(2,3) = -s_x*c_y
 rot_m(3,1) = c_x*s_y*c_z + s_x*s_z
 rot_m(3,2) = -c_x*s_y*s_z + s_x*c_z
 rot_m(3,3) = c_x*c_y

end subroutine make_rotation_matrix

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use physcon,      only: au, solarm, years
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_velocity_km_s,'wind_velocity', &
      'velocity at which wind is injected (km/s)',iunit)
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass per unit time (Msun/yr)',iunit)
 call write_inopt(wind_temperature,'wind_temperature','wind temperature at the injection point (K)',iunit)
 !call write_inopt(shift_spheres,'shift_spheres','shift the spheres of the wind',iunit)
 call write_inopt(wind_injection_radius,'wind_inject_radius', &
      'radius of injection of the wind (code units)',iunit)
 call write_inopt(wind_dr_on_dp,'wind_dr_on_dp','desired ratio of sphere spacing to particle spacing',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres','number of boundary spheres (integer)',iunit)
#ifdef BOWEN
 write(iunit,"(/,a)") '# options controlling bowen dust around central star'

 call write_inopt(bowen_kappa,'bowen_kappa','constant gas opacity (cm²/g)',iunit)
 call write_inopt(bowen_Teff,'bowen_Teff','central star effective temperature (K)',iunit)
 call write_inopt(bowen_kmax,'bowen_kmax','maximum dust opacity (cm²/g)',iunit)
 call write_inopt(bowen_Tcond,'bowen_Tcond','dust condensation temperature (K)',iunit)
 call write_inopt(bowen_delta,'bowen_delta','condensation temperature range (K)',iunit)
 call write_inopt(bowen_L/solarl,'bowen_L','central star luminosity (Lsun)',iunit)
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
 call write_inopt(wind_osc_period_days,'wind_osc_period','stellar pulsation period (days)',iunit)
 call write_inopt(wind_osc_vamplitude_km_s,'wind_osc_vampl', &
      'velocity amplitude of the pulsations (km/s)',iunit)
#endif
end subroutine

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
 real :: Rstar
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('wind_velocity')
    read(valstring,*,iostat=ierr) wind_velocity_km_s
    ngot = ngot + 1
    if (wind_velocity < 0.)    call fatal(label,'invalid setting for wind_velocity (<0)')
    !if (wind_velocity > 1.e10) call error(label,'wind_velocity is huge!!!')
 case('wind_mass_rate')
    read(valstring,*,iostat=ierr) wind_mass_rate_Msun_yr
    ngot = ngot + 1
    if (wind_mass_rate < 0.)    call fatal(label,'invalid setting for wind_mass_rate (<0)')
 case('wind_temperature')
    read(valstring,*,iostat=ierr) wind_temperature
    ngot = ngot + 1
    if (wind_temperature < 0.)    call fatal(label,'invalid setting for wind_temperature (<0)')
 case('iwind_resolution')
    read(valstring,*,iostat=ierr) iwind_resolution
    ngot = ngot + 1
    if (iwind_resolution < 1) call fatal(label,'iwind_resolution must be bigger than zero')
 case('wind_dr_on_dp')
    read(valstring,*,iostat=ierr) wind_dr_on_dp
    ngot = ngot + 1
    if (wind_dr_on_dp <= 0.) call fatal(label,'wind_dr_on_dp must be >=0')
 case('shift_spheres')
    read(valstring,*,iostat=ierr) shift_spheres
    ngot = ngot + 1
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    ngot = ngot + 1
    if (iboundary_spheres <= 0) call fatal(label,'iboundary_spheres must be > 0')
 case('wind_inject_radius')
    read(valstring,*,iostat=ierr) wind_injection_radius
    ngot = ngot + 1
    if (wind_injection_radius < 0.) call fatal(label,'invalid setting for wind_inject_radius (<0)')
#ifdef BOWEN
 case('bowen_kappa')
    read(valstring,*,iostat=ierr) bowen_kappa
    ngot = ngot + 1
    if (bowen_kappa < 0.)    call fatal(label,'invalid setting for bowen_kappa (<0)')
 case('bowen_Teff')
    read(valstring,*,iostat=ierr) bowen_Teff
    ngot = ngot + 1
    if (bowen_Teff < 0.)    call fatal(label,'invalid setting for bowen_Teff (<0)')
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
 case('bowen_L')
    read(valstring,*,iostat=ierr) bowen_L
    bowen_L = bowen_L * solarl
    ngot = ngot + 1
    if (bowen_L < 0.) call fatal(label,'invalid setting for bowen_L (<0)')
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
    if (bowen_Cprime < 0.) call fatal(label,'invalid setting for bowen_Cprime (<0)')
 case('wind_osc_period')
    read(valstring,*,iostat=ierr) wind_osc_period_days
    ngot = ngot + 1
    if (wind_osc_period_days < 0.) call fatal(label,'invalid setting for wind_osc_period (<0)')
 case('wind_osc_vampl')
    read(valstring,*,iostat=ierr) wind_osc_vamplitude_km_s
    wind_velocity_km_s = 0. ! set wind veolicty to zero when pulsating star
    ngot = ngot + 1
    if (wind_osc_vamplitude_km_s <= 0.) call fatal(label,'invalid setting for wind_osc_vamp (<0)')
#endif
 case default
    imatch = .false.
 end select
#ifdef BOWEN
 noptions = 15
 Rstar = sqrt(bowen_L/(4.*pi*steboltz*bowen_Teff**4))/au
 !if you launch the wind inside the photosphere, make sure the wind temperature >= star's effective temperature
 if (wind_injection_radius <= Rstar .and. wind_temperature < bowen_Teff) call fatal(label,'invalid setting for wind_temperature (< bowen_Teff)')
#else
 noptions = 6
#endif
 igotall = (ngot >= noptions)

end subroutine

end module inject
