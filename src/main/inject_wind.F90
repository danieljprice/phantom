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
!    bowen_Cprime     -- radiative cooling rate (g.s/cm³)
!    bowen_L          -- central star luminosity (Lsun)
!    bowen_Tcond      -- condensation temperature of dust (K)
!    bowen_Teff       -- effective temperature of the central star (K)
!    bowen_delta      -- condensation temperature range (K)
!    bowen_kappa      -- constant gas opacity (cm²/g)
!    bowen_kmax       -- maximum dust opacity (cm²/g)
!    ihandled_spheres -- handle inner spheres of the wind (integer)
!    iwind_resolution -- resolution of the wind -- DO NOT CHANGE DURING SIMULATION --
!    shift_spheres    -- shift the spheres of the wind
!    wind_osc_period  -- period of the oscillations (days)
!    wind_sphdist     -- distance between spheres / neighbours -- DO NOT CHANGE DURING SIMULATION --
!    wind_temperature -- wind temperature at the injection point (K)
!
!  DEPENDENCIES: bowen_dust, eos, icosahedron, infile_utils, io, part,
!    partinject, physcon, timestep, units
!+
!--------------------------------------------------------------------------
module inject
 use physcon, only: au, solarm, years, solarl
 implicit none
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
 real, public ::    wind_osc_period = 350. * 86400.
 real, public ::    wind_osc_vamplitude = 2.89d5
 real, public ::    wind_velocity = 8.622d5
 real, public ::    wind_mass_rate = 1.04d-7 * solarm/years
 real, public ::    wind_temperature = 1662.
 integer, public :: iwind_resolution = 5
 real, public ::    wind_sphdist = 0.2
 real, public ::    shift_spheres = 0.
 integer, public :: ihandled_spheres = 2
 real, public ::    wind_injection_radius = 4.786 * au
#else
 real, public ::    wind_velocity = 35.d5
 real, public ::    wind_mass_rate = 7.65d-7 * solarm/years
 real, public ::    wind_temperature = 3000.
 real, public ::    wind_gamma = 5./3.
 integer, public :: iwind_resolution = 4
 real, public ::    wind_sphdist = 10.
 real, public ::    shift_spheres = 3.
 integer, public :: ihandled_spheres = 3
 real, public ::    wind_injection_radius = 1.7 * au
 real, parameter :: wind_osc_vamplitude = 0.
#endif

! Calculated from the previous parameters
 real, public ::    mass_of_particles, mass_of_spheres, time_between_spheres, neighbour_distance
 integer, public :: particles_per_sphere

 private

 logical, parameter :: wind_verbose = .false.
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3), u_to_temperature_ratio

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine wind_init(setup)
 use physcon,     only: Rg, solarm, years, pi
 use icosahedron, only: compute_matrices, compute_corners
 use timestep,    only: dtmax
 use eos,         only: gmw
 use units,       only: utime
 logical, intent(in) :: setup

 real, parameter :: phi = (sqrt(5.)+1.)/2. ! Golden ratio
 real :: irrational_numbre_close_to_one

 u_to_temperature_ratio = Rg/(gmw*(wind_gamma-1.))
 particles_per_sphere = 20 * (2*iwind_resolution*(iwind_resolution-1)) + 12
 neighbour_distance = 2./((2.*iwind_resolution-1.)*sqrt(sqrt(5.)*phi))
 mass_of_particles = wind_sphdist * neighbour_distance * wind_injection_radius * wind_mass_rate &
                    / (particles_per_sphere * max(wind_velocity,wind_osc_vamplitude))
 mass_of_spheres = mass_of_particles * particles_per_sphere
 time_between_spheres = mass_of_spheres / wind_mass_rate
 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 if (wind_verbose) then
    print *,'particles_per_sphere ',particles_per_sphere
    print *,'neighbour_distance ',neighbour_distance
    print *,'mass_of_particles ',mass_of_particles
    print *,'mass_of_spheres ',mass_of_spheres
    print *,'time_between_spheres ',time_between_spheres
 endif

 ! adjusting dtmax to avoid uterm < 0 errors
 if (setup) then
    irrational_numbre_close_to_one = pi / 3.
    dtmax = (.5 * irrational_numbre_close_to_one * time_between_spheres)/utime
 endif
end subroutine

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time_u,dtlast_u,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use io,      only: iprint
 use units,   only: utime, umass
 use physcon, only: au
#ifdef BOWEN
 use bowen_dust, only: bowen_init
#endif
 real,    intent(in)    :: time_u, dtlast_u
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: outer_sphere, inner_sphere, inner_handled_sphere, i
 real :: time, dtlast, local_time, r, v, u, rho, e, mass_lost
 logical, save :: first_run = .true.

 if (first_run) then
    call wind_init(.false.)
#ifdef BOWEN
    call bowen_init(u_to_temperature_ratio,bowen_kappa,bowen_kmax,bowen_L,bowen_Cprime,&
         bowen_Tcond, bowen_delta,bowen_Teff)
#endif
    first_run = .false.
 endif

 time = time_u * utime
 dtlast = dtlast_u * utime
 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1
 inner_sphere = floor(time/time_between_spheres)
 inner_handled_sphere = inner_sphere + ihandled_spheres
 print *,'spheres',inner_sphere,inner_handled_sphere,outer_sphere,npart
 do i=inner_sphere+ihandled_spheres,outer_sphere,-1
    local_time = time - (i-shift_spheres) * time_between_spheres
    call compute_sphere_properties(local_time, r, v, u, rho, e, i, xyzmh_ptmass)
    if (wind_verbose) then
       write(iprint,*) '* Sphere:'
       write(iprint,*) '   Number: ', i
       write(iprint,*) '   Local Time: ', local_time
       write(iprint,*) '   Radius: ', r, '(', r/au, ' au)'
       write(iprint,*) '   Expansion velocity: ', v, '(', v/1.d5, ' km/s)'
       write(iprint,*) '   Temperature: ', u / u_to_temperature_ratio, ' K'
       write(iprint,*) '   Density: ', rho, ' (g/cm³)'
       write(iprint,*) '   Constant e: ', e, ' (J/g)'
       write(iprint,*) ''
    endif
    if (i > inner_sphere) then
       call inject_geodesic_sphere(i, (ihandled_spheres-i+inner_sphere)*particles_per_sphere+1, r, v, u, rho, &
            npart, npartoftype, xyzh, vxyzu) ! handled sphere
       print *,'handled',i,npart
    else
       call inject_geodesic_sphere(i, npart+1, r, v, u, rho, npart, npartoftype, xyzh, vxyzu) ! injected sphere
       print *,'ejected',i,npart
    endif
 enddo
 print *,'npart',npart,inner_sphere-outer_sphere+1
 mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
 xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - mass_lost/umass

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Time derivative of r and v, for Runge-Kutta iterations (stationnary trans-sonic solution)
!+
!-----------------------------------------------------------------------
subroutine drv_dt(rv, drv, GM)
 use physcon, only: Rg
 use eos, only: gmw
 real, intent(in) :: rv(2), GM
 real, intent(out) :: drv(2)

 real :: r, v, dr_dt, r2, T, vs2, dv_dr, dv_dt

 r = rv(1)
 v = rv(2)
 dr_dt = v
 r2 = r*r
 T = wind_temperature * (wind_injection_radius**2 * wind_velocity / (r2 * v))**(wind_gamma-1.)
 vs2 = wind_gamma * Rg * T / gmw
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
subroutine compute_sphere_properties(local_time, r, v, u, rho, e, sphere_number, xyzmh_ptmass)
#ifdef BOWEN
 use bowen_dust, only: pulsating_bowen_wind_profile
#endif
 integer, intent(in) :: sphere_number
 real, intent(in) :: local_time, xyzmh_ptmass(:,:)
 real, intent(out) :: r, v, u, rho, e
#ifdef BOWEN
 real :: central_star_mass, central_star_radius

 central_star_mass   = xyzmh_ptmass(4,wind_emitting_sink)
 central_star_radius = xyzmh_ptmass(5,wind_emitting_sink)
 call pulsating_bowen_wind_profile(local_time, r, v, u, rho, e, sphere_number,&
       wind_mass_rate,wind_injection_radius,wind_velocity,wind_osc_vamplitude,&
       wind_osc_period,shift_spheres,central_star_mass,time_between_spheres,&
       wind_temperature,central_star_radius,wind_gamma)
#else
 call stationnary_adiabatic_wind_profile(local_time, r, v, u, rho, e)
#endif

end subroutine

!-----------------------------------------------------------------------
!+
!  stationnary adiabatic supersonic wind
!+
!-----------------------------------------------------------------------
subroutine stationnary_adiabatic_wind_profile(local_time, r, v, u, rho, e)
 use part,    only: nptmass, xyzmh_ptmass
 use physcon, only: pi, Gg
 use units,   only: umass
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
    GM = Gg * xyzmh_ptmass(4,wind_emitting_sink) * umass
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
 T = wind_temperature * (wind_injection_radius**2 * wind_velocity / (r**2 * v))**(wind_gamma-1.)
 u = T * u_to_temperature_ratio
 rho = wind_mass_rate / (4.*pi*r**2*v)
 e = .5*v**2 - GM/r + wind_gamma*u
end subroutine


!-----------------------------------------------------------------------
!+
!  Inject a quasi-spherical distribution of particles.
!+
!-----------------------------------------------------------------------
subroutine inject_geodesic_sphere(sphere_number, first_particle, r, v, u, rho, npart, npartoftype, xyzh, vxyzu)
 use icosahedron, only: pixel2vector
 use units,       only: udist, utime, umass
 use partinject,  only: add_or_update_particle
 use part,        only: igas, hrho, xyzmh_ptmass, vxyz_ptmass, nptmass
 integer, intent(in) :: sphere_number, first_particle
 real, intent(in) :: r, v, u, rho
 integer, intent(inout) :: npart, npartoftype(:)
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)

 real :: rotation_angles(3), r_sim, v_sim, u_sim, h_sim
 real :: radial_unit_vector(3), rotmat(3,3), radial_unit_vector_rotated(3)
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
 r_sim = r / udist
 v_sim = v / (udist/utime)
 u_sim = u / (udist**2/utime**2)
 h_sim = hrho(rho / (umass/udist**3))

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
    radial_unit_vector = radial_unit_vector_rotated !/ sqrt(dot_product(radial_unit_vector_rotated, radial_unit_vector_rotated))
    particle_position = r_sim*radial_unit_vector
    particle_velocity = v_sim*radial_unit_vector
    if (nptmass > 0) then
       particle_position = particle_position + xyzmh_ptmass(1:3,wind_emitting_sink)
       particle_velocity = particle_velocity + vxyz_ptmass(1:3,wind_emitting_sink)
    endif
    call add_or_update_particle(igas, particle_position, particle_velocity, h_sim, u_sim, first_particle+j, &
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
end subroutine

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use physcon,      only: au, solarm, years
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_velocity,'wind_velocity', &
      'velocity at which wind is injected (cm/s) -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(wind_mass_rate/(solarm/years),'wind_mass_rate', &
      'wind mass per unit time (Msun/yr) -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(wind_temperature,'wind_temperature','wind temperature at the injection point (K)',iunit)
 call write_inopt(iwind_resolution,'iwind_resolution','resolution of the wind -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(wind_sphdist,'wind_sphdist','distance between spheres / neighbours -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(shift_spheres,'shift_spheres','shift the spheres of the wind',iunit)
 call write_inopt(ihandled_spheres,'ihandled_spheres','handle inner spheres of the wind (integer)',iunit)
 call write_inopt(wind_injection_radius/au,'wind_inject_radius', &
      'radius of injection of the wind (au) -- DO NOT CHANGE DURING SIMULATION --',iunit)
#ifdef BOWEN
 write(iunit,"(/,a)") '# options controlling bowen dust around central star'

 call write_inopt(bowen_kappa,'bowen_kappa','constant gas opacity (cm²/g)',iunit)
 call write_inopt(bowen_Teff,'bowen_Teff','effective temperature of the central star (K)',iunit)
 call write_inopt(bowen_kmax,'bowen_kmax','maximum dust opacity (cm²/g)',iunit)
 call write_inopt(bowen_Tcond,'bowen_Tcond','condensation temperature of dust (K)',iunit)
 call write_inopt(bowen_delta,'bowen_delta','condensation temperature range (K)',iunit)
 call write_inopt(bowen_L/solarl,'bowen_L','central star luminosity (Lsun)',iunit)
 call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
 call write_inopt(wind_osc_period/86400.,'wind_osc_period','period of the oscillations (days)',iunit)
 call write_inopt(wind_osc_vamplitude,'wind_osc_vampl', &
      'amplitude of velocity variations during oscillations (cm/s)',iunit)
#endif
end subroutine

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use physcon, only: au, solarm, years
 use io,      only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 integer :: noptions
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('wind_velocity')
    read(valstring,*,iostat=ierr) wind_velocity
    ngot = ngot + 1
    if (wind_velocity < 0.)    call fatal(label,'invalid setting for wind_velocity (<0)')
    if (wind_velocity > 1.e10) call error(label,'wind_velocity is huge!!!')
 case('wind_mass_rate')
    read(valstring,*,iostat=ierr) wind_mass_rate
    wind_mass_rate = wind_mass_rate * (solarm/years)
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
 case('wind_sphdist')
    read(valstring,*,iostat=ierr) wind_sphdist
    ngot = ngot + 1
    if (wind_sphdist <= 0.) call fatal(label,'wind_sphdist must be >=0')
 case('shift_spheres')
    read(valstring,*,iostat=ierr) shift_spheres
    ngot = ngot + 1
 case('ihandled_spheres')
    read(valstring,*,iostat=ierr) ihandled_spheres
    ngot = ngot + 1
    if (ihandled_spheres <= 0) call fatal(label,'ihandled_spheres must be > 0')
 case('wind_inject_radius')
    read(valstring,*,iostat=ierr) wind_injection_radius
    wind_injection_radius = wind_injection_radius * au
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
    read(valstring,*,iostat=ierr) wind_osc_period
    wind_osc_period = wind_osc_period*86400.
    ngot = ngot + 1
    if (wind_osc_period < 0.) call fatal(label,'invalid setting for wind_osc_period (<0)')
 case('wind_osc_vampl')
    read(valstring,*,iostat=ierr) wind_osc_vamplitude
    ngot = ngot + 1
    if (wind_osc_vamplitude < 0.) call fatal(label,'invalid setting for wind_osc_vamp (<0)')
#endif
 case default
    imatch = .false.
 end select
#ifdef BOWEN
 noptions = 17
#else
 noptions = 8
#endif
 igotall = (ngot >= noptions)

end subroutine

end module inject
