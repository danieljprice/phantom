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
!    bowen_Tcond       -- dust condensation temperature (K)
!    bowen_Teff        -- central star effective temperature (K)
!    bowen_delta       -- condensation temperature range (K)
!    bowen_kappa       -- constant gas opacity (cm²/g)
!    bowen_kmax        -- maximum dust opacity (cm²/g)
!    iboundary_spheres -- number of boundary spheres (integer)
!    wind_dr_on_dp     -- desired ratio of sphere spacing to particle spacing
!    wind_mass_rate    -- wind mass per unit time (Msun/yr)
!    wind_osc_period   -- stellar pulsation period (days)
!    wind_temperature  -- wind temperature at the injection point (K)
!
!  DEPENDENCIES: bowen_dust, eos, icosahedron, infile_utils, io, part,
!    partinject, physcon, timestep, units
!+
!--------------------------------------------------------------------------
module inject
 use physcon, only: au, solarm, years, solarl, pi
 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject
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
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3), u_to_temperature_ratio
#ifdef BOWEN
 real :: dr3, rho_ini
#endif

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use physcon,     only:Rg, days, km
 use icosahedron, only:compute_matrices, compute_corners
 use eos,         only:gmw, gamma
 use units,       only:unit_velocity, umass, utime
 use part,        only:massoftype,igas,iboundary
 use io,          only:iverbose
 use injectutils, only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 integer, intent(out) :: ierr
 real :: mV_on_MdotR,mass_of_particles
 real :: dr,dp,mass_of_particles1,wind_velocity_max
 real, parameter :: irrational_number_close_to_one = pi/3.
 !
 ! return without error
 !
 ierr = 0
 !
 ! convert input parameters to code units
 !
#ifdef BOWEN
 wind_osc_period        = wind_osc_period_days * (days/utime)
#endif
 wind_osc_vamplitude    = wind_osc_vamplitude_km_s * (km / unit_velocity)
 wind_velocity          = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate         = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 wind_velocity_max      = max(wind_velocity,wind_osc_vamplitude)
 !
 ! compute the dimensionless resolution factor m V / (Mdot R)
 ! where m = particle mass and V, Mdot and R are wind parameters
 !
 mass_of_particles = massoftype(igas)
 massoftype(iboundary) = mass_of_particles
 mV_on_MdotR = mass_of_particles*wind_velocity_max/(wind_mass_rate*wind_injection_radius)
 !
 ! solve for the integer resolution of the geodesic spheres
 ! gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
 !
 iwind_resolution     = get_sphere_resolution(wind_dr_on_dp,mV_on_MdotR)
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
 print*,'radial pulsation amplitude/ R* = ',wind_osc_vamplitude*wind_osc_period/(2.*pi)/wind_injection_radius
 print*,'pulsation period in code units = ',wind_osc_period
 print*,'injection time interval in code units = ',time_between_spheres
 print*,'particles_per_sphere ',particles_per_sphere
 !sanity checks
 ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
 if (wind_osc_period/time_between_spheres < 10. ) print *,'WARNING! only ',wind_osc_period/time_between_spheres,&
      ' shells will be ejected during a pulsation period'
 ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
 if (1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius > 0.2)  print*,'WARNING! the width of the boundary layer = ',&
      1.-(iboundary_spheres*dr3)**(1./3.)/wind_injection_radius,' Rinject'
#endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,          only:iprint,fatal
 use eos,         only:gamma
#ifdef BOWEN
 use bowen_dust,  only:bowen_init
#endif
 use units,       only:unit_velocity,unit_density,udist
 use part,        only:igas,iboundary,nptmass
 use injectutils, only:inject_geodesic_sphere
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real, parameter :: irrational_number_close_to_one = pi/3.
 integer :: outer_sphere, inner_sphere, inner_boundary_sphere, i, ierr
 real    :: local_time, GM, r, v, u, rho, e, mass_lost, x0(3), v0(3)
 logical, save :: first_run = .true.
 character(len=*), parameter :: label = 'inject_particles'

 if (first_run) then
    call init_inject(ierr)
#ifdef BOWEN
    call bowen_init(u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,bowen_Cprime,&
         bowen_Tcond,bowen_delta,bowen_Teff,wind_osc_vamplitude,wind_osc_period,&
         iboundary_spheres*particles_per_sphere)
#endif
    first_run = .false.
 endif

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
 do i=inner_sphere+iboundary_spheres,outer_sphere,-1
    local_time = time - (i-shift_spheres) * time_between_spheres
    call compute_sphere_properties(time,local_time,gamma,GM,r,v,u,rho,e,i,&
         inner_sphere,inner_boundary_sphere)
    if (wind_verbose) then
       write(iprint,*) '   Sphere    : ', i,(i-shift_spheres)
       write(iprint,*) '   Local Time: ', local_time,time,(i-shift_spheres)*time_between_spheres
       write(iprint,*) '   Radius: ', r, '(', r*(udist/au), ' au)'
       write(iprint,*) '   Expansion velocity: ', v*unit_velocity, '(', v*unit_velocity/1.d5, ' km/s)'
       write(iprint,*) '   Density: ', rho * unit_density, ' (g/cm³)'
       write(iprint,*) ''
       !read*
    endif
    if (i > inner_sphere) then
       call inject_geodesic_sphere(i, (iboundary_spheres-i+inner_sphere)*particles_per_sphere+1, &
            iwind_resolution, r, v, u, rho,  geodesic_R, geodesic_V, &
            npart, npartoftype, xyzh, vxyzu, iboundary, x0, v0) ! boundary sphere
    else
       call inject_geodesic_sphere(i, npart+1, iwind_resolution, r, v, u, rho, geodesic_R, geodesic_V,&
            npart, npartoftype, xyzh, vxyzu, igas, x0, v0) ! injected sphere
    endif
 enddo
! print *,'npart',npart,inner_sphere-outer_sphere+1
 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_lost
 endif

 !
 ! return timestep constraint to ensure that time between sphere
 ! injections is adequately resolved
 !
 dtinject = (.5 * irrational_number_close_to_one * time_between_spheres)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Time derivative of r and v, for Runge-Kutta iterations (stationary trans-sonic solution)
!+
!-----------------------------------------------------------------------
subroutine drv_dt(rv,drv,GM,gamma)
 real, intent(in) :: rv(2),GM,gamma
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

end subroutine drv_dt

!-----------------------------------------------------------------------
!+
!  Compute the radius, velocity and temperature of a sphere at the current local time
!+
!-----------------------------------------------------------------------
subroutine compute_sphere_properties(time,local_time,gamma,GM,r,v,u,rho,e,sphere_number, &
                                     inner_sphere,inner_boundary_sphere)
#ifdef BOWEN
 !use bowen_dust, only: pulsating_bowen_wind_profile
#endif
 integer, intent(in)  :: sphere_number, inner_sphere, inner_boundary_sphere
 real,    intent(in)  :: time,local_time,gamma,GM
 real,    intent(out) :: r, v, u, rho, e
#ifdef BOWEN
 integer, parameter :: nrho_index = 10
 integer :: k
 real :: surface_radius,r3
 logical :: verbose = .true.

 v = wind_osc_vamplitude* cos(2.*pi*time/wind_osc_period) !same velocity for all wall particles
 surface_radius = wind_injection_radius + wind_osc_vamplitude*wind_osc_period/(2.*pi)*sin(2.*pi*time/wind_osc_period)
 r3 = surface_radius**3-dr3
 do k = 2,sphere_number-inner_sphere
    r3 = r3-dr3*(r3/surface_radius**3)**(nrho_index/3.)
 enddo
 r = r3**(1./3)
 !r = (surface_radius**3-(sphere_number-inner_sphere)*dr3)**(1./3)
 !rho = rho_ini
 u = wind_temperature * u_to_temperature_ratio
 rho = rho_ini*(surface_radius/r)**nrho_index
 e = .5*v**2 - GM/r + gamma*u
 if (verbose) then
    if (sphere_number > inner_sphere) then
       ! print '("boundary, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
       !       '" r = ",es11.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)', &
       !       sphere_number,inner_sphere,surface_radius,r,v,&
       !       time/wind_osc_period,time_between_spheres/wind_osc_period
    else
       print '("ejected, i = ",i5," inner = ",i5," base_r = ",es11.4,'// &
             '" r = ",es11.4," v = ",es11.4," phase = ",f7.4," feject = ",f4.3)', &
             sphere_number,inner_sphere,surface_radius,r,v,&
             time/wind_osc_period,time_between_spheres/wind_osc_period
    endif
 endif
#else
 call stationary_adiabatic_wind_profile(local_time, r, v, u, rho, e, gamma, GM)
#endif

end subroutine compute_sphere_properties

!-----------------------------------------------------------------------
!+
!  stationary adiabatic supersonic wind
!+
!-----------------------------------------------------------------------
subroutine stationary_adiabatic_wind_profile(local_time, r, v, u, rho, e, gamma, GM)
 real, intent(in)  :: local_time, GM, gamma
 real, intent(out) :: r, v, u, rho, e
 real :: dt, rv(2), k1(2), k2(2), k3(2), k4(2), T
 integer, parameter :: N = 10000
 integer :: i

 dt = local_time / N
 rv(1) = wind_injection_radius
 rv(2) = wind_velocity
 ! Runge-Kutta iterations
 do i=1,N
    call drv_dt(rv,          k1, GM, gamma)
    call drv_dt(rv+dt/2.*k1, k2, GM, gamma)
    call drv_dt(rv+dt/2.*k2, k3, GM, gamma)
    call drv_dt(rv+dt*k3,    k4, GM, gamma)
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
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
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
#ifdef BOWEN
 real :: Rstar
#endif
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
 if (wind_injection_radius <= Rstar .and. wind_temperature < bowen_Teff) then
    call fatal(label,'invalid setting for wind_temperature (< bowen_Teff)')
 endif
#else
 noptions = 6
#endif
 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject
