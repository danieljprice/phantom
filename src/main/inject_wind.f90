!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles wind particle injection
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - iboundary_spheres  : *number of boundary spheres (integer)*
!   - iwind_resolution   : *if<>0 set number of particles on the sphere, reset particle mass*
!   - nfill_domain       : *number of spheres used to set the background density profile*
!   - outer_boundary     : *delete gas particles outside this radius (au)*
!   - wind_type          : *find transonic solution (1=yes,0=no)*
!   - wind_inject_radius : *wind injection radius (au, if 0 takes Rstar)*
!   - wind_mass_rate     : *wind mass loss rate (Msun/yr)*
!   - wind_shell_spacing : *desired ratio of sphere spacing to particle spacing*
!   - wind_temperature   : *wind temperature at injection radius (K, if 0 takes Teff)*
!   - wind_velocity      : *injection wind velocity (km/s, if wind_type = 0)*
!
! :Dependencies: cooling_molecular, dim, dust_formation, eos, icosahedron,
!   infile_utils, injectutils, io, options, part, partinject, physcon,
!   ptmass_radiation, setbinary, timestep, units, wind, wind_equations
!
 use dim, only:isothermal,nucleation

 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      wind_injection_radius,set_default_options_inject,update_injected_par
 private
!
!--runtime settings for this module
!
! Read from input file
 integer:: wind_type = 0
 integer:: iboundary_spheres = 5
 integer:: iwind_resolution = 5
 integer:: nfill_domain = 0
 real :: wind_velocity_km_s = 10.
 real :: wind_mass_rate_Msun_yr
 real :: wind_mass_rate_alt_Msun_yr = 1e-5
 real :: wind_injection_radius_au
 real :: wind_temperature
 real :: outer_boundary_au = 30.
 real :: wind_shell_spacing = 1.
 real :: pulsation_period
 real :: pulsation_period_days = 0.
 real :: piston_velocity_km_s = 0.
 real :: dtpulsation = huge(0.)
 real :: dt_alt_rate_yr = 0.
 real :: dt_alt_rate
 real :: tbeg_alt_rate_yr = 10.
 real :: tbeg_alt_rate

! global variables
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: u_to_temperature_ratio,wind_mass_rate,piston_velocity,wind_velocity,&
      omega_osc,mass_of_spheres,time_between_spheres,neighbour_distance,&
      dr3,Rstar_cgs,Rinject,wind_injection_radius,wind_injection_speed,rho_ini,&
      deltaR_osc,Mstar_cgs, time_period,orbital_period,mass_of_particles,&
      check_mass,rho_factor,dt_shift,wind_mass_rate_alt,wind_mass_rate_low
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
 use options,           only:icooling,ieos
 use io,                only:fatal,iverbose
 use setbinary,         only:get_eccentricity_vector
 use timestep,          only:tmax,time
 use wind_equations,    only:init_wind_equations
 use wind,              only:setup_wind,save_windprofile
 use physcon,           only:mass_proton_cgs, kboltz, Rg, days, km, au, years, solarm, pi, Gg
 use icosahedron,       only:compute_matrices, compute_corners
 use eos,               only:gmw,gamma,polyk
 use units,             only:unit_velocity, umass, utime, udist
 use part,              only:xyzmh_ptmass,vxyz_ptmass,massoftype,igas,iboundary,imloss,ilum,iTeff,iReff,nptmass,npart
 use injectutils,       only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance

 integer, intent(out) :: ierr
 integer :: nzones_per_sonic_point,new_nfill
 real :: mV_on_MdotR,initial_wind_velocity_cgs,dist_to_sonic_point,semimajoraxis_cgs
 real :: dr,dp,mass_of_particles1,tcross,tend,rsonic,tsonic,initial_Rinject,tboundary
 real :: separation_cgs,ecc(3),eccentricity,Tstar

 if (icooling > 0) nwrite = nwrite+1
 ierr = 0
 dt_shift = 0.
 rho_factor = 1.

 pulsating_wind = (pulsation_period_days > 0.) .and. (piston_velocity_km_s > 0.)
 if (pulsating_wind .and. ieos == 6) call fatal(label,'cannot use ieos=6 with pulsation')

 Tstar              = xyzmh_ptmass(iTeff,wind_emitting_sink)
 Rstar_cgs          = xyzmh_ptmass(iReff,wind_emitting_sink)*udist
 Mstar_cgs          = xyzmh_ptmass(4,wind_emitting_sink)*umass
 !
 ! convert input parameters to code units
 !
 wind_velocity  = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 wind_mass_rate_alt = wind_mass_rate_alt_Msun_yr * (solarm/umass) / (years/utime)
 wind_mass_rate_low = wind_mass_rate
 if (wind_injection_radius_au == 0.)  wind_injection_radius_au = Rstar_cgs/au
 if (wind_temperature == 0.)  wind_temperature = Tstar
 wind_injection_radius = wind_injection_radius_au * au / udist
 if (pulsating_wind) then
    pulsation_period = pulsation_period_days * (days/utime)
    piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
    dtpulsation      = pulsation_period/50.
    omega_osc        = 2.*pi/pulsation_period
    deltaR_osc       = pulsation_period*piston_velocity/(2.*pi)
 else
    omega_osc        = 0.d0
    deltaR_osc       = 0.d0
    piston_velocity  = 0.d0
 endif
 dt_alt_rate = dt_alt_rate_yr * years/utime
 tbeg_alt_rate = tbeg_alt_rate_yr * years/utime
 initial_wind_velocity_cgs = (piston_velocity+wind_velocity)*unit_velocity
 ! if you consider non-trans-sonic wind, provide and input wind velocity /= 0.
 if (initial_wind_velocity_cgs <= 0. .and. wind_type /= 1) call fatal(label,'zero input wind velocity')
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif

 time_period    = 0.
 orbital_period = 1.50
 if (nptmass == 2) then
    separation_cgs = sqrt(sum((xyzmh_ptmass(1:3,2)-xyzmh_ptmass(1:3,1))**2))*udist
    ecc             = get_eccentricity_vector(xyzmh_ptmass,vxyz_ptmass,1,2)
    eccentricity    = sqrt(sum(ecc(1:3)**2))
    !stars initially positionned at apastron (see set_binary)
    semimajoraxis_cgs = separation_cgs/(1.+eccentricity)
    orbital_period  = sqrt(4.*pi**2*semimajoraxis_cgs**3/(Gg*(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))*solarm))    ! cgs
    print *,'eccentricity                 = ',eccentricity
    print *,'orbital_period (days)        = ',orbital_period/86400.
 endif

 if (isothermal) then
    !Rstar_cgs = wind_injection_radius_au *au
    wind_temperature = polyk* mass_proton_cgs/kboltz * unit_velocity**2*gmw
 endif
 if (wind_injection_radius_au < Rstar_cgs/au) then
    print *,'WARNING wind_inject_radius < Rstar (au)',wind_injection_radius_au,Rstar_cgs/au
    !call fatal(label,'WARNING wind_inject_radius (< Rstar)')
 endif

 initial_Rinject = wind_injection_radius
 if (pulsating_wind) then
    !implement background spheres starting from the smallest radius
    !initial_Rinject = min(initial_Rinject,xyzmh_ptmass(iReff,wind_emitting_sink)-deltaR_osc)
 else
    if (initial_Rinject < xyzmh_ptmass(5,wind_emitting_sink)) then
       print *,'stop wind_inject_radius < Racc (au)',wind_injection_radius_au,xyzmh_ptmass(5,wind_emitting_sink)
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
         wind_temperature,initial_wind_velocity_cgs,rsonic,tsonic,wind_type)
    initial_Rinject = Rinject/udist
 endif
 Rinject = initial_Rinject
 wind_injection_speed = initial_wind_velocity_cgs/unit_velocity
 rho_ini = wind_mass_rate/(4.*pi*Rinject**2*wind_injection_speed)
 check_mass = massoftype(igas)

 if (iwind_resolution == 0) then
    !
    ! resolution is specified in terms of number of smoothing lengths
    ! per distance to sonic point (if trans-sonic wind)
    !
    if (wind_type == 1) then
       nzones_per_sonic_point = 8
       dist_to_sonic_point = rsonic/udist-Rinject
       dr = abs(dist_to_sonic_point)/nzones_per_sonic_point
       mass_of_particles = rho_ini*dr**3
       massoftype(igas) = mass_of_particles
       print*,' suggesting ',mass_of_particles, ' based on desired dr = ',dr,' dist-to-sonic=',dist_to_sonic_point
    else
       mass_of_particles = massoftype(igas)
    endif
    !
    ! compute the dimensionless resolution factor m V / (Mdot R)
    ! where m = particle mass and V, Mdot and R are wind parameters
    !
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
 if (npart > 0 .and. abs(log10(check_mass/mass_of_particles)) > 1e-10) then
    print *,'check_mass = ',check_mass
    print *,'particle mass = ',mass_of_particles
    print *,'number of particles = ',npart
    call fatal(label,"you cannot reset the particle's mass")
 endif
 massoftype(iboundary) = mass_of_particles

 mass_of_spheres = mass_of_particles * particles_per_sphere
 nwall_particles = iboundary_spheres*particles_per_sphere
 dr3 = 3.*mass_of_spheres/(4.*pi*rho_ini)
 if (wind_type == 3) then
    time_between_spheres = mass_of_spheres/calc_wind_mass_loss_rate(time)
 else
    time_between_spheres = mass_of_spheres/wind_mass_rate_low
 endif
 if (time_between_spheres > tmax)  then
    call logging(initial_wind_velocity_cgs,rsonic,Tsonic,Tboundary)
    print *,'time_between_spheres = ',time_between_spheres,' < tmax = ',tmax
    call fatal(label,'no shell ejection : tmax < time_between_spheres')
 endif

 ! compute matrices for geodesic
 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 ! compute full evolution : get tcross & save 1D profile
 if ( .not. pulsating_wind .or. nfill_domain > 0) then
    tboundary = (iboundary_spheres+nfill_domain)*time_between_spheres
    tend      = max(tmax,tboundary)*utime
    call save_windprofile(real(Rinject*udist),real(wind_injection_speed*unit_velocity),&
         wind_temperature,real(outer_boundary_au*au), tend, tcross, 'wind_profile1D.dat')
    if (tboundary > tmax) then
       print *,'simulation time < time to reach the last boundary shell'
    endif
    !reduce the value of nfill_domain if it is too large
    if (tcross < 1.d98) then
       if (tcross/time_between_spheres < 1.d4) then
          new_nfill = min(nfill_domain,int(tcross/time_between_spheres)-iboundary_spheres)
          if (new_nfill /= nfill_domain .and. new_nfill > 0) then
             nfill_domain = new_nfill
             print *,'number of background shells set to',nfill_domain
          endif
       endif
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

 xyzmh_ptmass(imloss,wind_emitting_sink) = wind_mass_rate

 ! logging
 call logging(initial_wind_velocity_cgs,rsonic,Tsonic,Tboundary)

end subroutine init_inject


!-----------------------------------------------------------------------

subroutine logging(initial_wind_velocity_cgs,rsonic,Tsonic,Tboundary)

!-----------------------------------------------------------------------

 use physcon,           only:pi,gg,au,km
 use units,             only:udist,unit_velocity,utime
 use timestep,          only:dtmax,tmax
 use ptmass_radiation,  only:alpha_rad
 use part,              only:xyzmh_ptmass, iReff, ispinx, ispiny, ispinz

 real, intent(in) :: initial_wind_velocity_cgs,rsonic,Tsonic,Tboundary
 integer :: ires_min
 real :: vesc,wind_rotation_speed,rotation_speed_crit

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-alpha_rad)/Rstar_cgs)
 write (*,'(/,2(3x,A,es11.4),A,i7)')&
      'mass_of_particles     (cu) : ',mass_of_particles,&
      'time_between_spheres  (cu) : ',time_between_spheres,&
      'particles per sphere       : ',particles_per_sphere
 write (*,'(2(3x,A,es11.4))') &
      'distance between spheres   : ',wind_shell_spacing*neighbour_distance,&
      'distance between injection : ',time_between_spheres*wind_injection_speed
 write (*,'(2(3x,A,es11.4))') &
      'wind_temperature           : ',wind_temperature,&
      'injection_radius (au)      : ',Rinject*au/udist,&
      'stellar_radius (au)        : ',Rstar_cgs/udist, &
      'rho_ini (cgs)              : ',rho_ini
 if (wind_type == 3) then
    write (*,'(2(3x,A,es11.4))') &
         'tmax (cu)                  : ',tmax,&
         'tbeg_alt_rate (cu)         : ',tbeg_alt_rate, &
         'dt_alt_rate (cu)           : ',dt_alt_rate, &
         'ejected particles          : ',dt_alt_rate/time_between_spheres*particles_per_sphere
 endif
 !print*,'hmax/dist_between_spheres  = ',wind_shell_spacing*neighbour_distance*&
 !      initial_wind_velocity_cgs**2/(vesc**2-initial_wind_velocity_cgs**2)
 if (wind_type == 1) then
    write (*,'(2(3x,A,es11.4),2(3x,A,es11.4,3x,i5))') &
         'distance to sonic point    : ',rsonic/udist-Rinject, &
         'sonic radius               : ',rsonic/udist,rsonic, &
         'time_to_sonic_point        : ',tsonic/utime, &
         'time_boundary              : ',tboundary, &
         'number of shells to sonic  : ',(rsonic/udist-Rinject)/(wind_shell_spacing*neighbour_distance)
 endif
 wind_rotation_speed = sqrt(sum(xyzmh_ptmass(ispinx:ispinz,1)**2))/xyzmh_ptmass(iReff,1)**2
 rotation_speed_crit = sqrt(xyzmh_ptmass(4,1)/xyzmh_ptmass(iReff,1))
 if (wind_rotation_speed > 1e-20) then
   ! problem about code units, need to update
    write (*,'(2(3x,A,es11.4))') &
         'rotation speed (km/s)      = ',wind_rotation_speed*unit_velocity/km,&
         'break-up velocity (km/s)   = ',rotation_speed_crit*unit_velocity/km,&
         'rotation_vel/critical_vel  = ',wind_rotation_speed/rotation_speed_crit,&
         'rotation_vel/vinject       = ',wind_rotation_speed/wind_injection_speed
    if (wind_rotation_speed/rotation_speed_crit > 1.) then
      print*,'CAREFUL : rotation velocity exceeding equatorial break-up velocity'
    endif
 endif
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

end subroutine logging


!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use physcon,           only:pi,au,solarm,years
 use io,                only:fatal,iverbose
 use wind,              only:interp_wind_profile !,wind_profile
 use part,              only:igas,iTeff,iReff,iboundary,nptmass,delete_particles_outside_sphere,&
                             delete_dead_particles_inside_radius,dust_temp,n_nucleation
 use partinject,        only:add_or_update_particle
 use injectutils,       only:inject_geodesic_sphere
 use units,             only:udist, utime, umass
 use dust_formation,    only:idust_opacity
 use ptmass_radiation,  only:isink_radiation

 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere, inner_sphere, inner_boundary_sphere, first_particle, i, ipart, &
            nreleased, nboundaries
 real    :: local_time, GM, r, v, u, rho, e, mass_lost, x0(3), v0(3), inner_radius, fdone, dum
 real    :: tolv,dtinject_next
 character(len=*), parameter :: label = 'inject_particles'
 logical, save :: released = .false.
 real :: JKmuS(n_nucleation)

 tolv = 10.
 dum  = 0.

 if (nptmass > 0 .and. wind_emitting_sink <= nptmass) then
    x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
    GM = xyzmh_ptmass(4,wind_emitting_sink)
    v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 else
    x0 = 0.
    v0 = 0.
    GM = 0.
 endif

 time_period  = time_period + dtlast*utime

 if (npart == 0) then
    !initialise domain with boundary particles
    ipart = iboundary
    nreleased = 0
    nboundaries = iboundary_spheres+nfill_domain
 else
    !
    ! delete particles that exit the outer boundary
    !
    inner_radius = wind_injection_radius + deltaR_osc*sin(omega_osc*time)

    if (outer_boundary_au > Rinject) call delete_particles_outside_sphere(x0,real(outer_boundary_au*au/udist),npart)
    call delete_dead_particles_inside_radius(x0,inner_radius,npart)
    if (npart_old /= npart .and. iverbose > 0) print *,'deleted ',npart_old-npart,'particles, remaining',npart
    npart_old = npart

    if (time_period > orbital_period .and. nptmass == 2) time_period = 0.

    nreleased = nfill_domain
    nboundaries = iboundary_spheres
    !release background shells
    ipart = igas
    if (.not.released) then
       if (nfill_domain > 0) print '(/,A,i7,3x,"Mdot : ",es12.4)','release particles :',&
            npart-nreleased*particles_per_sphere+1,npart,&
            wind_mass_rate/(solarm/umass)*(years/utime)
       do i = max(1,npart-nreleased*particles_per_sphere)+1,npart
          if (isothermal) then
             call add_or_update_particle(ipart,xyzh(1:3,i),vxyzu(1:3,i),xyzh(4,i),dum,&
                                         i,npart,npartoftype,xyzh,vxyzu)
          else
             call add_or_update_particle(ipart,xyzh(1:3,i),vxyzu(1:3,i),xyzh(4,i),vxyzu(4,i),&
                                         i,npart,npartoftype,xyzh,vxyzu)
          endif
       enddo
       released = .true.
    endif
 endif

! time varying mass loss rate
 dtinject_next = time_between_spheres
 if (wind_type == 3 .and. dt_alt_rate > 0.) then
    wind_mass_rate = calc_wind_mass_loss_rate(time)
    rho_factor = wind_mass_rate/wind_mass_rate_low !also implemenbt change in velocity
    time_between_spheres = mass_of_spheres / wind_mass_rate
    dtinject_next = mass_of_spheres / calc_wind_mass_loss_rate(time+dtlast)
    if (abs(1.-rho_factor) > 1.d-10) then
       dt_shift = iboundary_spheres*mass_of_spheres*(1./wind_mass_rate_low-1./wind_mass_rate)
    else
       dt_shift = 0.
    endif
    !print *,'!!!!! ',dt_shift,time_between_spheres,wind_mass_rate,wind_mass_rate_low
 endif

 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1 + nreleased
 inner_sphere = floor(time/time_between_spheres)+nreleased
 inner_boundary_sphere = inner_sphere + nboundaries
 !print '("@@@",2(i4),5(1x,es12.4))',inner_boundary_sphere,wind_type,time,dt_alt_rate,&
 !                time_between_spheres,calc_wind_mass_loss_rate(time)

 !only one sphere can be ejected at a time
 if (inner_sphere-outer_sphere > nboundaries) call fatal(label,'ejection of more than 1 sphere, timestep likely too large!')

 ! update properties of boundary particles
 do i=inner_boundary_sphere,outer_sphere,-1
    local_time = time + (iboundary_spheres+nfill_domain-i) * time_between_spheres + dt_shift
    !compute the radius, velocity, temperature, chemistry of a shell at the current local time
    v = wind_injection_speed
    r = Rinject
    if (pulsating_wind.and.released) then
       call pulsating_wind_profile(time,local_time, r, v, u, rho, e, GM, i, inner_sphere,dr3,rho_ini)
    else
       if (idust_opacity == 2) then
          call interp_wind_profile(time, local_time, r, v, u, rho, e, GM, fdone, JKmuS)
       else
          call interp_wind_profile(time, local_time, r, v, u, rho, e, GM, fdone)
       endif
       rho = rho*rho_factor
       if (iverbose > 0) print '(" ## update boundary ",i4,2(i4),i7,8(1x,es12.5))',i,&
            inner_sphere,outer_sphere,npart,time,local_time,r/xyzmh_ptmass(iReff,1),local_time*utime,&
            r*udist,v*udist/utime,wind_mass_rate /(solarm/umass) * (years/utime)
    endif

    if (i > inner_sphere) then
       ! boundary sphere
       first_particle = (nboundaries-i+inner_sphere)*particles_per_sphere+1
       !print '(" @@ update boundary ",i4,2(i4),i7,6(1x,es12.5))',i,inner_sphere,&
       !     outer_sphere,first_particle,time,local_time,r/xyzmh_ptmass(iReff,1),v,u,rho
       if (idust_opacity == 2) then
          call inject_geodesic_sphere(i, first_particle, iresolution, r, v, u, rho,  geodesic_R, geodesic_V, &
               npart, npartoftype, xyzh, vxyzu, ipart, x0, v0, JKmuS)
       else
          call inject_geodesic_sphere(i, first_particle, iresolution, r, v, u, rho,  geodesic_R, geodesic_V, &
               npart, npartoftype, xyzh, vxyzu, ipart, x0, v0)
       endif
       if (isink_radiation > 0) dust_temp(first_particle:first_particle+particles_per_sphere-1) = &
            xyzmh_ptmass(iTeff,wind_emitting_sink)
    else
       ! ejected particles + create new  inner sphere
       if (idust_opacity == 2) then
          call inject_geodesic_sphere(i, npart+1, iresolution, r, v, u, rho, geodesic_R, geodesic_V,&
               npart, npartoftype, xyzh, vxyzu, igas, x0, v0, JKmuS)
       else
          call inject_geodesic_sphere(i, npart+1, iresolution, r, v, u, rho, geodesic_R, geodesic_V,&
               npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
       endif
       !initialize dust temperature to star's effective temperature
       if (isink_radiation > 0) dust_temp(npart+1-particles_per_sphere:npart) = xyzmh_ptmass(iTeff,wind_emitting_sink)
       print '(" ## eject particles [",i7,"-",i7,"] ",3(i4),10(1x,es11.4))',npart+1-particles_per_sphere,npart,i,&
            inner_sphere,outer_sphere,time,local_time,r/xyzmh_ptmass(iReff,1),v,u,rho,wind_mass_rate/(solarm/umass)*(years/utime)
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
 dtinject = min(0.2*time_between_spheres,0.2*dtinject_next,dtpulsation)
 if (time <= 0.) dtinject = 0.01*dtinject

end subroutine inject_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  set the time varying mass loss rate
!+
!-----------------------------------------------------------------------
real function calc_wind_mass_loss_rate(time)

   real, intent(in) :: time

   if (time >= tbeg_alt_rate .and. time < tbeg_alt_rate+dt_alt_rate .and. dt_alt_rate > 0. ) then
      calc_wind_mass_loss_rate = wind_mass_rate_alt
   else
      calc_wind_mass_loss_rate = wind_mass_rate_low
   endif

end function calc_wind_mass_loss_rate


!-----------------------------------------------------------------------
!+
!  Oscillating inner boundary
!+
!-----------------------------------------------------------------------
subroutine pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,sphere_number, &
                                  inner_sphere,dr3,rho_ini)
 use physcon, only:pi
 use eos,     only:gamma
 use units,   only:unit_velocity
 integer, intent(in)  :: sphere_number, inner_sphere
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
! Fit an idealistic density and velocity profile of the stellar wind
!+
!-----------------------------------------------------------------------
subroutine fit_spherical_wind(xyzh,vxyzu,r_sep, r_outer, n_part, n0, m, v_inf)
 use part,   only: rhoh

 ! Data dictionary: Arguments
 real,intent(in)     :: xyzh(:,:), vxyzu(:,:)
 integer, intent(in) :: n_part
 real, intent(out)   :: n0, m, v_inf
 real, intent(in)    :: r_sep, r_outer

 ! Data dictionary: Additional parameters for calculations
 real, dimension(:), allocatable :: allR, allRho, allV
 real, dimension(:), allocatable :: density_fit_compZone, velocity_cutoff
 real, dimension(:), allocatable :: density_log_noCompZone, density_fit_noCompZone
 real, dimension(:), allocatable :: radius_log_noCompZone, radius_fit_noCompZone
 real, parameter                 :: radius_cutoffFactor2 = 9d-1
 integer                         :: n_velocity, n_density_noCompZone, n_density_compZone,i
 real                            :: radius_cutoff1, radius_cutoff2

 allocate(allR(n_part))
 allocate(allRho(n_part))
 allocate(allV(n_part))

 allR   = sqrt(xyzh(1,:)*xyzh(1,:)+xyzh(2,:)*xyzh(2,:)+xyzh(3,:)*xyzh(3,:))
 allV   = sqrt(vxyzu(1,:)*vxyzu(1,:)+vxyzu(2,:)*vxyzu(2,:)+vxyzu(3,:)*vxyzu(3,:))
 do i=1,n_part
    allRho(i) = rhoh(xyzh(4,i),mass_of_particles)
 enddo

 ! Initialisation
 radius_cutoff1 = 1.5 * r_sep
 radius_cutoff2 = radius_cutoffFactor2 * r_outer
 if (radius_cutoff1 > radius_cutoff2) radius_cutoff1 = 0.5 * radius_cutoff2

 ! Fit the exponent m outside compZone
 n_density_noCompZone = count(allR > r_sep)
 allocate(radius_fit_noCompZone(n_density_noCompZone))
 allocate(radius_log_noCompZone(n_density_noCompZone))
 allocate(density_fit_noCompZone(n_density_noCompZone))
 allocate(density_log_noCompZone(n_density_noCompZone))

 radius_fit_noCompZone  = pack(allR, allR > r_sep)
 radius_log_noCompZone  = log(radius_fit_noCompZone)
 density_fit_noCompZone = pack(allRho, allR > r_sep)
 density_log_noCompZone = log(density_fit_noCompZone)

 m = -( sum(radius_log_noCompZone * density_log_noCompZone) - &
          1.D0/n_density_noCompZone * sum(radius_log_noCompZone) * sum(density_log_noCompZone) ) / &
        ( sum(radius_log_noCompZone**2.D0) - 1.D0/n_density_noCompZone * sum(radius_log_noCompZone)**2. )

 ! Fit n0 inside compZone
 n_density_compZone = count(allR <= r_sep)
 allocate(density_fit_compZone(n_density_compZone))

 density_fit_compZone = pack(allRho, allR <= r_sep)
 n0 = sum(density_fit_compZone) / n_density_compZone

 ! Fit velocity profile
 n_velocity = count((radius_cutoff1 <= allR) .and. (allR <= radius_cutoff2))
 allocate(velocity_cutoff(n_velocity))
 velocity_cutoff = pack(allV, (radius_cutoff1 <= allR) .and. (allR <= radius_cutoff2))
 v_inf = sum(velocity_cutoff)/n_velocity

 deallocate(radius_fit_noCompZone)
 deallocate(radius_log_noCompZone)
 deallocate(density_fit_noCompZone)
 deallocate(density_log_noCompZone)
 deallocate(density_fit_compZone)
 deallocate(velocity_cutoff)
 deallocate(allR)
 deallocate(allRho)
 deallocate(allV)

end subroutine fit_spherical_wind


subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
 integer :: icase

 if (present(flag)) then
    icase = flag
 else
    icase = 0
 endif

 if (isothermal) then
    wind_velocity_km_s = 0.
    wind_mass_rate_Msun_yr = 8.2d-8
    wind_injection_radius_au = 0.
 else
    if (icase == 1) then
       !trans-sonic wind
       wind_velocity_km_s = 0.
       wind_mass_rate_Msun_yr = 1.d-5
       wind_injection_radius_au = 2.
       wind_temperature = 50000.
    else
       !super sonic-wind
       wind_velocity_km_s = 20.
       wind_mass_rate_Msun_yr = 1.d-5
       wind_injection_radius_au = 2.
       wind_temperature = 2500.
    endif
 endif

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use dim,          only: maxvxyzu
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_type,'wind_type','wind type 0=std, 1=transonic, 2=pulsation, 3=time dep',iunit)
 call write_inopt(wind_velocity_km_s,'wind_velocity','injection wind velocity (km/s, if wind_type = 0)',iunit)
 if (wind_type == 2) then
    call write_inopt(pulsation_period_days,'pulsation_period','stellar pulsation period (days)',iunit)
    call write_inopt(piston_velocity_km_s,'piston_velocity','velocity amplitude of the pulsation (km/s)',iunit)
 endif
 call write_inopt(wind_injection_radius_au,'wind_inject_radius','wind injection radius (au, if 0 takes Rstar)',iunit)
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass loss rate (Msun/yr)',iunit)
 if (wind_type == 3) then
    call write_inopt(wind_mass_rate_alt_Msun_yr,'wind_mass_rate_alt','high value of the mass loss rate (Msun/yr)',iunit)
    call write_inopt(dt_alt_rate_yr,'dt_alt_rate','duration of mass loss episode (yr, ignored if <=0)',iunit)
    call write_inopt(tbeg_alt_rate_yr,'tbeg_alt_rate','time when high mass loss rate starts (yr)',iunit)
 endif
 if (maxvxyzu==4) then
    call write_inopt(wind_temperature,'wind_temperature','wind temperature at injection radius (K, if 0 takes Teff)',iunit)
 endif
 call write_inopt(iwind_resolution,'iwind_resolution','if<>0 set number of particles on the sphere, reset particle mass',iunit)
 call write_inopt(nfill_domain,'nfill_domain','number of spheres used to set the background density profile',iunit)
 call write_inopt(wind_shell_spacing,'wind_shell_spacing','desired ratio of sphere spacing to particle spacing',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres','number of boundary spheres (integer)',iunit)
 call write_inopt(outer_boundary_au,'outer_boundary','delete gas particles outside this radius (au)',iunit)

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
 logical :: isowind = .true., init_opt = .false.
 character(len=30), parameter :: label = 'read_options_inject'

 if (.not.init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif
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
 case('wind_type')
    read(valstring,*,iostat=ierr) wind_type
    ngot = ngot + 1
    if (wind_type < 0 .or. wind_type > 3 ) call fatal(label,'invalid setting for wind_type ([0,3])')
 case('wind_mass_rate')
    read(valstring,*,iostat=ierr) wind_mass_rate_Msun_yr
    ngot = ngot + 1
    if (wind_mass_rate_Msun_yr < 0.) call fatal(label,'invalid setting for wind_mass_rate (<0)')
 case('wind_mass_rate_alt')
    read(valstring,*,iostat=ierr) wind_mass_rate_alt_Msun_yr
    ngot = ngot + 1
    if (wind_mass_rate_alt_Msun_yr < 0.) call fatal(label,'invalid setting for wind_mass_rate (<0)')
  case('dt_alt_rate')
    read(valstring,*,iostat=ierr) dt_alt_rate_yr
    if (dt_alt_rate_yr <= 0.) dt_alt_rate_yr = huge(0.)
    ngot = ngot + 1
 case('tbeg_alt_rate')
    read(valstring,*,iostat=ierr) tbeg_alt_rate_yr
     if (tbeg_alt_rate_yr < 0.) call fatal(label,'invalid setting for tbeg_alt_rate (<0)')
    ngot = ngot + 1
 case('pulsation_period')
    read(valstring,*,iostat=ierr) pulsation_period_days
    ngot = ngot + 1
    if (pulsation_period_days < 0.) call fatal(label,'invalid setting for pulsation_period (<0)')
 case('piston_velocity')
    read(valstring,*,iostat=ierr) piston_velocity_km_s
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 if (isothermal) then
    noptions = 9
 elseif (wind_type == 2) then
    noptions = 11
 elseif (wind_type == 3) then
    noptions = 12
 else
    noptions = 9
 endif
 igotall = (ngot >= noptions)
 if (trim(name) == '') ngot = 0
 !print '(a26,i3,i3,l1)',trim(name),ngot,noptions,igotall

end subroutine read_options_inject

end module inject
