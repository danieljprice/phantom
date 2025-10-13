!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
!   - B_r                : *radial magnetic field strength (G)*
!   - iboundary_spheres  : *number of boundary spheres (integer)*
!   - iwind_resolution   : *if<>0 set number of particles on the sphere, reset particle mass*
!   - rfill_domain       : *number of spheres used to set the background density profile*
!   - outer_boundary     : *delete gas particles outside this radius (au)*
!   - wind_type          : *find transonic solution (1=yes,0=no)*
!   - wind_inject_radius : *wind injection radius (au, if 0 takes Rstar)*
!   - wind_mass_rate     : *wind mass loss rate (Msun/yr)*
!   - wind_shell_spacing : *desired ratio of sphere spacing to particle spacing*
!   - wind_temperature   : *wind temperature at injection radius (K, if 0 takes Teff)*
!   - wind_velocity      : *injection wind velocity (km/s, if wind_type = 0)*
!
! :Dependencies: cooling_molecular, dim, dust_formation, eos, icosahedron,
!   infile_utils, injectutils, io, options, orbits, part, partinject,
!   physcon, ptmass_radiation, timestep, units, wind, wind_equations
!
 use dim, only:isothermal,nucleation,mhd

 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      wind_injection_radius,set_default_options_inject,update_injected_par

 private
!
!--runtime settings for this module
!
! Read from input file
 integer :: wind_type = 0
 integer :: iboundary_spheres = 5
 integer :: iwind_resolution = 5
 integer :: nfill_domain = 0
 real :: wind_velocity_km_s = 10.
 real :: wind_mass_rate_Msun_yr
 real :: wind_injection_radius_au
 real :: wind_temperature
 real :: outer_boundary_au = 30.
 real :: wind_shell_spacing = 1.
 real :: pulsation_period
 real :: pulsation_period_days = 0.
 real :: rfill_domain_au = 0.
 real :: piston_velocity_km_s = 0.
 real :: dtpulsation = huge(0.)
 real :: jet_edge_velocity = 0.
 real :: jet_opening_angle = 0.
 real :: jet_opening_angle_degree = 0.
 real :: B_r = 0.

! global variables
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: Rstar_cgs,Mstar_cgs,u_to_temperature_ratio,wind_velocity,&
      wind_injection_radius,wind_injection_speed
 real :: omega_puls,deltaR_puls,piston_velocity,dr3 !pulsations
 integer :: particles_per_sphere,nwrite

 logical :: pulsating_wind,onewind
 character(len=*), parameter :: label = 'inject_wind'

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use options,           only:icooling,ieos
 use io,                only:fatal
 use orbits,            only:get_eccentricity_vector,get_orbital_period
 use wind_equations,    only:init_wind_equations
 use wind,              only:setup_wind
 use physcon,           only:mass_proton_cgs, kboltz, Rg, days, km, au, years, solarm, pi, Gg
 use icosahedron,       only:compute_matrices, compute_corners
 use eos,               only:gmw,gamma,polyk
 use units,             only:unit_velocity,umass,udist
 use part,              only:xyzmh_ptmass,vxyz_ptmass,igas,iboundary,imloss,&
                             ilum,iTeff,iReff,nptmass,ivwind,iTwind,nptmass
 use injectutils,       only:get_parts_per_sphere,get_neighb_distance,init_jets
 use ptmass_radiation,  only:isink_radiation

 integer, intent(out) :: ierr
 integer :: isink
 real :: initial_wind_velocity_cgs,semimajoraxis_cgs,orbital_period,time_between_spheres,d_part
 real :: tcross,rinject,rsonic,tsonic,tfill,initial_rinject,tboundary,mu
 real :: separation_cgs,ecc(3),eccentricity

 onewind = abs(xyzmh_ptmass(ivwind,2)) < tiny(0.)
 ! change particle injection method if more than 1 sink is emitting a wind
 if(xyzmh_ptmass(imloss,1) /= 0 .and. xyzmh_ptmass(imloss,2) /= 0) then
    if (wind_type == 3) call init_jets(jet_edge_velocity,jet_opening_angle)
 else
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
 endif

 if (icooling > 0) nwrite = nwrite+1
 ierr = 0

 pulsating_wind = (pulsation_period_days > 0.) .and. (piston_velocity_km_s > 0.)
 if (ieos == 6) call fatal(label,'cannot use ieos=6 with pulsation')

! setup thermo
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif
 if (isothermal) wind_temperature = polyk * mass_proton_cgs/kboltz * unit_velocity**2*gmw

 call init_pulsating_wind(pulsating_wind)

 write(*,'(/,70("-"))')
 do isink = 1, nptmass
    if (abs(xyzmh_ptmass(ivwind,isink)) < tiny(0.)) cycle
    Rstar_cgs = xyzmh_ptmass(iReff,isink)*udist
    Mstar_cgs = xyzmh_ptmass(4,isink)*umass
    wind_temperature = xyzmh_ptmass(iTwind,isink)
    if (wind_temperature < 0.001) wind_temperature = xyzmh_ptmass(iTeff,isink)
    if (isink_radiation == 4)  then
       ! take the initial velocity of the wind to ba a fraction of the terminal velocity ~ 10-20%
       if (xyzmh_ptmass(ivwind,isink) < 1.2e8/unit_velocity) then
          wind_velocity = 0.2 * xyzmh_ptmass(ivwind,isink)
       else
          wind_velocity = 0.1 * xyzmh_ptmass(ivwind,isink)
       endif
    else
       wind_velocity = xyzmh_ptmass(ivwind,isink)
    endif

    if (nptmass < 2)  then
       if (wind_injection_radius_au == 0.) wind_injection_radius_au = Rstar_cgs/au
       wind_injection_radius = wind_injection_radius_au * au / udist
    else
       wind_injection_radius = xyzmh_ptmass(iReff,isink)
       wind_injection_radius_au = wind_injection_radius * udist / au
    endif

    initial_wind_velocity_cgs = (piston_velocity+wind_velocity)*unit_velocity
    ! if you consider non-trans-sonic wind, provide and input wind velocity /= 0.
    if (initial_wind_velocity_cgs <= 0. .and. wind_type /= 1) call fatal(label,'zero input wind velocity')

    if (wind_injection_radius_au < Rstar_cgs/au) then
       print *,'WARNING wind_inject_radius < Rstar (au)',wind_injection_radius_au,Rstar_cgs/au
       !call fatal(label,'WARNING wind_inject_radius (< Rstar)')
    endif

    initial_rinject = wind_injection_radius
    if (pulsating_wind) then
       !implement background spheres starting from the smallest radius
       !initial_rinject = min(initial_rinject,xyzmh_ptmass(iReff,isink)-deltaR_puls)
    else
       if (initial_rinject < xyzmh_ptmass(5,isink)) then
          print *,'stop wind_inject_radius < Racc (au)',wind_injection_radius_au,xyzmh_ptmass(5,isink)
          call fatal(label,'invalid setting wind_inject_radius < accretion radius')
       endif
       call init_wind_equations(xyzmh_ptmass(4,isink),xyzmh_ptmass(iTeff,isink),u_to_temperature_ratio)
       !if (ieos == 6 .and. wind_temperature /= star_Teff) then
       ! wind_injection_radius_au = Rstar_cgs/au*(star_Teff/wind_temperature)**(1./wind_expT)
       ! wind_injection_radius  = wind_injection_radius_au * au / udist
       !endif
       !wind_injection_radius = max(wind_injection_radius_au * au, Rstar_cgs) / udist
       !Rstar = min(wind_injection_radius_au*au,Rstar_cgs)

       ! integrate wind equation to get initial velocity and sonic radius to set resolution
       rinject = initial_rinject*udist
       call setup_wind(isink,Mstar_cgs,xyzmh_ptmass(imloss,isink),u_to_temperature_ratio, rinject,&
            wind_temperature,initial_wind_velocity_cgs,rsonic,tsonic,wind_type)
       initial_rinject = rinject/udist
    endif

    rinject = initial_rinject
    wind_injection_speed = initial_wind_velocity_cgs/unit_velocity

    if (isink == 1) call init_resolution(rinject,rsonic,nptmass)
    call init_sink_resolution(isink,time_between_spheres,d_part,rinject)

    ! compute 1D wind profile to get tcross & save 1D profile
    if (.not. pulsating_wind .or. rfill_domain_au > 0.) then
       call set_1D_wind_profile(isink,d_part,rinject,time_between_spheres,tboundary,tcross,tfill,onewind)
    endif

    ! logging
    if (xyzmh_ptmass(imloss,isink) > 0.) &
         call logging(isink,rinject,time_between_spheres,d_part,rsonic,tsonic,tboundary,tcross,tfill)


 enddo

 if (nptmass == 2) then
    separation_cgs = sqrt(sum((xyzmh_ptmass(1:3,2)-xyzmh_ptmass(1:3,1))**2))*udist
    ecc             = get_eccentricity_vector(xyzmh_ptmass,vxyz_ptmass,1,2)
    eccentricity    = sqrt(sum(ecc(1:3)**2))
    ! stars initially positioned at apastron (see set_binary)
    semimajoraxis_cgs = separation_cgs/(1.+eccentricity)
    mu = Gg*(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))*solarm          ! standard gravitational parameter G*(m1+m2), here in cgs
    orbital_period  = get_orbital_period(mu,semimajoraxis_cgs)    ! period in seconds
    print *,'eccentricity                 = ',eccentricity
    print *,'orbital_period (days)        = ',orbital_period/86400.
    print *,'orbital speed sink 1 (km/s)  = ',sqrt(dot_product(vxyz_ptmass(1:3,1),vxyz_ptmass(1:3,1)))*unit_velocity/km
    print *,'orbital speed sink 2 (km/s)  = ',sqrt(dot_product(vxyz_ptmass(1:3,2),vxyz_ptmass(1:3,2)))*unit_velocity/km
 endif

end subroutine init_inject


!-------------------------------------------------------------------------------
!+
!  set particle mass or resolution depending on iwind_resolution (set by sink 1)
!+
!-------------------------------------------------------------------------------
subroutine init_resolution(rinject,rsonic,nsink)

 use part,         only:xyzmh_ptmass,massoftype,igas,iboundary,imloss,ieject,ivwind,iReff
 use units,        only:udist
 use physcon,      only:pi
 use injectutils,  only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 use io,           only:fatal

 integer, intent(in) :: nsink
 real,    intent(in) :: rinject,rsonic

 integer :: nzones_per_sonic_point,iresolution,isink=1
 real    :: mV_on_MdotR,dr,dist_to_sonic_point,mass_of_particles,neighbour_distance,rho_ini

 if (iwind_resolution == 0 .and. nsink == 1) then
    !
    ! resolution is specified in terms of number of smoothing lengths
    ! per distance to sonic point (if trans-sonic wind)
    !
    if (wind_type == 1) then
       nzones_per_sonic_point = 8
       dist_to_sonic_point = rsonic/udist-rinject
       dr = abs(dist_to_sonic_point)/nzones_per_sonic_point
       rho_ini = xyzmh_ptmass(imloss,isink)/(4.*pi*rinject**2*wind_injection_speed)
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
    mV_on_MdotR = mass_of_particles*xyzmh_ptmass(ivwind,isink)/(xyzmh_ptmass(imloss,isink)*rinject)
    !
    ! solve for the integer resolution of the geodesic spheres
    ! gives number of particles on the sphere via N = 20*(2*q*(q - 1)) + 12
    !
    iresolution = get_sphere_resolution(wind_shell_spacing,mV_on_MdotR)
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    print *,'iwind_resolution equivalence = ',iresolution
 else
    iresolution = abs(iwind_resolution)
    particles_per_sphere = get_parts_per_sphere(iresolution)
    neighbour_distance   = get_neighb_distance(iresolution)
    mass_of_particles    = wind_shell_spacing*neighbour_distance*xyzmh_ptmass(iReff,isink)*&
                           xyzmh_ptmass(imloss,isink)/(particles_per_sphere * xyzmh_ptmass(ivwind,isink))
    massoftype(igas)     = mass_of_particles
 endif
 xyzmh_ptmass(ieject,isink) = particles_per_sphere
 massoftype(iboundary) = mass_of_particles

end subroutine init_resolution


!-----------------------------------------------------------------------
!+
!  set resolution for other wind emitting sink particles
!+
!-----------------------------------------------------------------------
subroutine init_sink_resolution(isink,time_between_spheres,d_part,rinject)

 use io,       only:fatal
 use physcon,  only:pi,solarm,years
 use units,    only:umass,utime
 use timestep, only:tmax
 use part,     only:xyzmh_ptmass,massoftype,igas,iboundary,imloss,&
                            npart,ieject,ivwind,iReff
 use injectutils, only:get_neighb_distance

 integer, intent(in) :: isink
 real, optional, intent(in) :: rinject
 real, intent(out) :: time_between_spheres,d_part
 real :: check_mass,res,mass_of_particles,mass_of_spheres,mdot_save

 mass_of_particles = massoftype(igas)
 mdot_save = xyzmh_ptmass(imloss,isink)

 if (xyzmh_ptmass(imloss,isink) > 0.) then
    res = (sqrt(4.*pi)*wind_shell_spacing*xyzmh_ptmass(iReff,isink)*xyzmh_ptmass(imloss,isink)/&
         (xyzmh_ptmass(ivwind,isink)*mass_of_particles))**(2./3.)
    !print *,res,wind_shell_spacing,xyzmh_ptmass(iReff,isink),xyzmh_ptmass(imloss,isink),xyzmh_ptmass(ivwind,isink),mass_of_particles
    xyzmh_ptmass(ieject,isink)   = nint(res+0.5)
    xyzmh_ptmass(imloss,isink)   = xyzmh_ptmass(imloss,isink) * sqrt(xyzmh_ptmass(ieject,isink)/res)**3
    d_part = get_neighb_distance(int(xyzmh_ptmass(ieject,isink)))
    check_mass = wind_shell_spacing*d_part*xyzmh_ptmass(iReff,isink)*&
         xyzmh_ptmass(imloss,isink)/(xyzmh_ptmass(ieject,isink) * xyzmh_ptmass(ivwind,isink))
    print"(/,' number of particles per sphere for sink',i2,' = ',i7)",isink,nint(res)
    if (abs(log(abs(mdot_save/xyzmh_ptmass(imloss,isink)))) > 1.d-6) &
         print"(' Mdot [Msun/yr] : ',es10.3,' --> ',es10.3)",mdot_save*(umass*years)/(solarm*utime),&
         xyzmh_ptmass(imloss,isink)*(umass*years)/(solarm*utime)
 else
    d_part = 0.
 endif

 if (npart > 0 .and. abs(log10(check_mass/mass_of_particles)) > 1e-10) then
    print *,'check_mass    = ',check_mass
    print *,'particle mass = ',mass_of_particles
    print *,'nbr particles = ',npart
    print *,'shell spacing = ',wind_shell_spacing
    print *,'d_part = ',d_part
    print *,'Reff   = ',xyzmh_ptmass(iReff,isink)
    print *,'mdot   = ',xyzmh_ptmass(imloss,isink)
    print *,'res    = ',xyzmh_ptmass(ieject,isink)
    print *,'vwind  = ',xyzmh_ptmass(ivwind,isink)
    call fatal(label,"you cannot reset the particle's mass")
 endif

 !spheres properties
 mass_of_spheres = massoftype(igas) * xyzmh_ptmass(ieject,isink)
 time_between_spheres = mass_of_spheres / xyzmh_ptmass(imloss,isink)
 if (time_between_spheres > tmax)  then
    if (xyzmh_ptmass(imloss,isink) > 0.) call logging(isink,rinject,time_between_spheres,d_part)
    print *,'time_between_spheres = ',time_between_spheres,' < tmax = ',tmax
    call fatal(label,'no shell ejection : tmax < time_between_spheres')
 endif


end subroutine init_sink_resolution

!-----------------------------------------------------------------------

subroutine logging(isink,rinject,time_between_spheres,neighbour_distance,rsonic,tsonic,&
     tboundary,tcross,tfill)

!-----------------------------------------------------------------------

 use physcon,           only:pi,gg,au,km
 use units,             only:udist,unit_velocity,utime
 use timestep,          only:dtmax
 use ptmass_radiation,  only:alpha_rad
 use part,              only:massoftype,igas,xyzmh_ptmass,iReff,&
                             ispinx,ispiny,ispinz,ieject,imloss

 integer, intent(in) :: isink
 real, intent(in) :: rinject,time_between_spheres,neighbour_distance
 real, optional, intent(in) :: rsonic,tsonic,tboundary,tcross,tfill
 integer :: ires_min
 logical :: lsonic
 real :: vesc,wind_rotation_speed,rotation_speed_crit

 if (.not.present(rsonic)) then
    lsonic = .false.
 else
    lsonic = rsonic/udist > rinject
 endif

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-alpha_rad)/Rstar_cgs)
 write (*,'(/,2(3x,A,es11.4))')&
      'mass_of_particles       : ',massoftype(igas),&
      'time_between_spheres    : ',time_between_spheres
 write (*,'(2(3x,A,es11.4),3x,A,i7)') &
      'dist between spheres    : ',wind_shell_spacing*neighbour_distance,&
      'dist between injection  : ',time_between_spheres*wind_injection_speed,&
      'particles per sphere    : ',nint(xyzmh_ptmass(ieject,isink))
 write (*,'(3(3x,A,es11.4))') &
      'wind_temperature        : ',wind_temperature,&
      'injection_radius   (au) : ',rinject*au/udist
 write (*,'(2(3x,A,es11.4))') &
      'stellar_radius     (au) : ',Rstar_cgs/udist, &
      'rho_ini           (cgs) : ',xyzmh_ptmass(imloss,isink)/(4.*pi*rinject**2*wind_injection_speed)
 if (present(tcross)) then
    write (*,'(3(3x,A,es11.4))') &
      'crossing time      (cu) : ',tcross,&
      'filling time       (cu) : ',tfill/utime,&
      'boundary time      (cu) : ',tboundary
 endif

 !print*,'hmax/dist_between_spheres  = ',wind_shell_spacing*neighbour_distance*&
 !      initial_wind_velocity_cgs**2/(vesc**2-initial_wind_velocity_cgs**2)
 if (wind_type == 1 .and. present(rsonic)) then
    write (*,'(3(3x,A,es11.4),3x,A,i5)') &
         'distance to sonic point    : ',rsonic/udist-rinject, &
         'sonic radius               : ',rsonic/udist, &
         'time_to_sonic_point        : ',tsonic/utime, &
         'number of shells to sonic  : ',(rsonic/udist-rinject)/(wind_shell_spacing*neighbour_distance)
 endif
 wind_rotation_speed = sqrt(sum(xyzmh_ptmass(ispinx:ispinz,isink)**2))/xyzmh_ptmass(iReff,isink)**2
 rotation_speed_crit = sqrt(xyzmh_ptmass(4,isink)/xyzmh_ptmass(iReff,isink))
 if (wind_rotation_speed > 1e-20) then
    write (*,'(4(3x,A,es11.4))') &
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
    print*,'width of the boundary layer/ R* (should be < 1) = ',1.-(iboundary_spheres*dr3)**(1./3.)/rinject
    print*,'radial pulsation amplitude/ R* = ',piston_velocity*pulsation_period/(2.*pi*rinject),dr3/rinject
    print*,'pulsation period in code units = ',pulsation_period
    !sanity checks
    ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
    if (pulsation_period/time_between_spheres < 10. ) print *,'WARNING! only ',pulsation_period/time_between_spheres,&
         ' shells will be ejected during a pulsation period'
    ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
    if (1.-(iboundary_spheres*dr3)**(1./3.)/rinject > 0.2)  print*,'WARNING! the width of the boundary layer = ',&
         1.-(iboundary_spheres*dr3)**(1./3.)/rinject,' rinject'
 elseif (lsonic) then
    !save a few models before the particles reach the sonic point
    if (dtmax > tsonic/utime) print *,'WARNING! dtmax > time to sonic point'
    !if solution subsonic, minimum resolution required so a few shells can be inserted between the injection radius and the sonic point
    ires_min = int(iboundary_spheres*wind_shell_spacing*0.5257/(rsonic/udist/rinject-1.)+.5)
    if (abs(iwind_resolution) < ires_min) print *,'WARNING! resolution too low to pass sonic point : |iwind_resolution| < ',ires_min
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
 use part,              only:massoftype,igas,iReff,iboundary,nptmass,delete_particles_outside_sphere,&
                             delete_dead_particles_inside_radius,n_nucleation,ieject,imloss,ivwind,rhoh,Bevol,Bxyz
 use partinject,        only:add_or_update_particle
 use injectutils,       only:use_icosahedron
 use units,             only:udist,umass,utime
 use dust_formation,    only:idust_opacity
 use ptmass_radiation,  only:isink_radiation

 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere,inner_sphere,inner_boundary_sphere,ifirst,i,ipart,j, &
            nreleased,nboundaries,isink,itype,ires,nfill
 real    :: local_time,GM,r,v,u,rho,e,mass_lost,x0(3),v0(3),inner_radius,fdone,dum
 real    :: mass_of_spheres,time_between_spheres,rinject
 character(len=*), parameter :: label = 'inject_particles'
 logical, save :: released = .false.
 real :: JKmuS(n_nucleation)

 dum  = 0.
 dtinject = huge(dtinject)

 do isink = 1,nptmass

    if (xyzmh_ptmass(imloss,isink) <= 0.) cycle
    if (abs(xyzmh_ptmass(ivwind,isink)) < tiny(0.)) cycle

    x0 = xyzmh_ptmass(1:3,isink)
    GM = xyzmh_ptmass(4,isink)
    v0 = vxyz_ptmass(1:3,isink)

    particles_per_sphere = nint(xyzmh_ptmass(ieject,isink))
    mass_of_spheres =  massoftype(igas) * particles_per_sphere
    dtinject = mass_of_spheres / xyzmh_ptmass(imloss,isink)
    time_between_spheres = dtinject
    if (use_icosahedron) then
       ires = iwind_resolution
    else
       ires = particles_per_sphere
    endif

    wind_injection_radius = xyzmh_ptmass(iReff,isink)
    rinject = wind_injection_radius

    if (isink_radiation == 4)  then
    ! take the initial velocity of the wind as 10% of the terminal velocity
       wind_injection_speed = 0.10 * xyzmh_ptmass(ivwind,isink)
    else
       wind_injection_speed = xyzmh_ptmass(ivwind,isink)
    endif
 !if 2 winds, background particles must originate from sink 2 to keep the memory location of the boundary particles
    if (onewind .or. isink == 2) then
       nfill = nfill_domain
    else
       nfill = 0
    endif

    if (npart > 0 .and. dtlast > 0.) then
    !
    ! delete particles that exit the outer boundary
    !
       inner_radius = wind_injection_radius + deltaR_puls*sin(omega_puls*time)

       if (outer_boundary_au > rinject) call delete_particles_outside_sphere(x0,real(outer_boundary_au*au/udist),npart)
       call delete_dead_particles_inside_radius(x0,inner_radius,npart)
       if (npart_old /= npart .and. iverbose > 0) print *,'deleted ',npart_old-npart,'particles, remaining',npart
       npart_old = npart

       nreleased = nfill_domain
       nboundaries = iboundary_spheres
       !release background shells
       ipart = igas
       if (.not.released .and. (isink == 2 .or. onewind)) then
          if (nfill > 0) print '(/,"release background particles :",i7,3x,"Mdot : ",es12.4,", isink=",i1)',&
               nreleased*particles_per_sphere,xyzmh_ptmass(imloss,isink)/(solarm/umass)*(years/utime),isink
          do i = max(1,npart-nreleased*particles_per_sphere)+1,npart
             if (isothermal) then
                call add_or_update_particle(ipart,xyzh(1:3,i),vxyzu(1:3,i),xyzh(4,i),dum,&
                     i,npart,npartoftype,xyzh,vxyzu)
             else
                call add_or_update_particle(ipart,xyzh(1:3,i),vxyzu(1:3,i),xyzh(4,i),vxyzu(4,i),&
                     i,npart,npartoftype,xyzh,vxyzu)
             endif
             if (mhd) call set_injected_Bfield(xyzmh_ptmass(:,isink),xyzh(:,i),Bevol(:,i),Bxyz(:,i),massoftype(igas))
          enddo
          released = .true.
       endif
    else
    !initialise domain with boundary particles
       ipart = iboundary
       nreleased = 0
       nboundaries = iboundary_spheres+nfill
    endif

    outer_sphere = floor((time-dtlast)/time_between_spheres) + 1 + nreleased
    inner_sphere = floor(time/time_between_spheres)+nreleased
    inner_boundary_sphere = inner_sphere + nboundaries

    if (nptmass < 2) then
! only one sphere can be ejected at a time (not valid for 2 wind-emitting sinks)
       if (inner_sphere-outer_sphere > nboundaries) call fatal(label,'ejection of more than 1 sphere, timestep likely too large!')
    endif

 !
 ! eject particles and update properties of boundary particles
 !
    do i=inner_boundary_sphere,outer_sphere,-1
       local_time = time + (iboundary_spheres+nfill-i) * time_between_spheres
    !compute the radius, velocity, temperature, chemistry of a shell at the current local time
       v = wind_injection_speed
       r = rinject
       if (pulsating_wind.and.released) then
          call pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,i,inner_sphere)
       else
          if (idust_opacity == 2) then
             call interp_wind_profile(time,local_time,r,v,u,rho,e,GM,fdone,isink,JKmuS)
          else
             call interp_wind_profile(time,local_time,r,v,u,rho,e,GM,fdone,isink)
          endif
          if (iverbose > 0) print '(" ## update boundary ",i4,2(i4),i7,i2,8(1x,es12.5))',i,&
               inner_sphere,outer_sphere,npart,isink,time,local_time,r/xyzmh_ptmass(iReff,isink),v*udist/utime,&
               xyzmh_ptmass(imloss,isink) /(solarm/umass) * (years/utime)
       endif

       if (i > inner_sphere) then
       ! boundary sphere
          if (isink == 1) ifirst = (nboundaries-i+inner_sphere)*particles_per_sphere+1
          if (.not. onewind .and. isink ==2) ifirst = (nboundaries-i+inner_sphere)*particles_per_sphere+1 &
               + iboundary_spheres*int(xyzmh_ptmass(ieject,1))
          itype = ipart
          if (iverbose > 0) print '(" @@ update boundary ",i1,i4,2(i4),i7,i7,8(1x,es12.5))',isink,i,inner_sphere,&
               outer_sphere,ipart,ifirst,time,local_time,r/xyzmh_ptmass(iReff,isink),v*udist/utime,u,rho
       else
       ! ejected particles + create new  inner sphere
          ifirst = npart+1
          itype = igas
          if (iverbose > 0) print '(" ## eject particles [",i7,"-",i7,"], sink=",i1,4(i4),7(1x,es11.4))',&
               npart+1-particles_per_sphere,npart,isink,i,inner_sphere,outer_sphere,int(time/time_between_spheres),&
               time,local_time,r/xyzmh_ptmass(iReff,isink),v,u,rho,xyzmh_ptmass(imloss,isink)/(solarm/umass)*(years/utime)
       endif
       if (idust_opacity == 2) then
          call inject_sphere(i,ifirst,ires,r,v,u,rho,npart,npartoftype,xyzh,vxyzu,itype,x0,v0,isink,JKmuS)
       else
          call inject_sphere(i,ifirst,ires,r,v,u,rho,npart,npartoftype,xyzh,vxyzu,itype,x0,v0,isink)
       endif
       if (mhd) then
          do j = ifirst,ifirst+particles_per_sphere-1
             call set_injected_Bfield(xyzmh_ptmass(:,isink),xyzh(:,j),Bevol(:,j),Bxyz(:,j),massoftype(igas))
          enddo
       endif
    enddo

    if (nfill > 0 .and. outer_sphere == 1) print '(3x,"injecting background particles up to r = ",f8.2,&
      &" (au), nparticles = ",i8,", isink=",i1)',r,particles_per_sphere*(nfill_domain+iboundary_spheres),isink

 ! update sink particle properties
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    xyzmh_ptmass(4,isink) = xyzmh_ptmass(4,isink) - mass_lost
    if (pulsating_wind) then
       inner_radius = wind_injection_radius + deltaR_puls*sin(omega_puls*time)
       !v2 xyzh_ptmass(5,isink) = (inner_radius**3-dr3)**(1./3.) !accretion radius
       xyzmh_ptmass(5,isink) = inner_radius
    endif

 !
 ! return timestep constraint to ensure that time between sphere
 ! injections is adequately resolved
 !
 !dr = neighbour_distance*wind_injection_radius
 !dtinject = 0.25*dr/sqrt(cs2max)
    dtinject = min(0.2*time_between_spheres,dtpulsation,dtinject)
 enddo
 if (time <= 0.) dtinject = 0.01*dtinject

end subroutine inject_particles

subroutine set_injected_Bfield(xyzmh_ptmassi,xyzhi,Bevoli,Bxyzi,pmassi)
   use part,  only:rhoh
   use units, only:unit_Bfield
   real, intent(in) :: xyzmh_ptmassi(:),xyzhi(:),pmassi
   real, intent(out) :: Bevoli(:),Bxyzi(:)
   real :: r,r_hat(3),B_r_code,rhoi,dx(3)

   dx(1:3)=xyzhi(1:3)-xyzmh_ptmassi(1:3)
   r = sqrt(dot_product(dx,dx))
   r_hat = xyzhi(1:3)/r
   B_r_code = B_r/unit_Bfield
   rhoi = rhoh(pmassi,xyzhi(4))
   Bevoli(1:3) = B_r_code*r_hat / rhoi
   Bxyzi(1:3) = B_r_code*r_hat

end subroutine set_injected_Bfield


!-----------------------------------------------------------------------
!+
!  inject gas particles and/or reset position of boundary particles
!+
!-----------------------------------------------------------------------
subroutine inject_sphere(i,ifirst,ires,r,v,u,rho,npart,npartoftype,xyzh,vxyzu,itype,x0,v0,isink,JKmuS)

 use ptmass_radiation,  only:isink_radiation
 use part,              only:iTeff,dust_temp,xyzmh_ptmass
 use injectutils,       only:inject_geodesic_sphere

 integer, intent(in) :: i,ifirst,ires,isink,itype
 integer, intent(inout) :: npart
 real,    intent(in), optional :: JKmuS(:)
 real,    intent(in) :: x0(3),v0(3),r,v,u,rho
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)

 if (present(JKmuS)) then
    call inject_geodesic_sphere(i, ifirst, ires, r, v, u, rho,  geodesic_R, geodesic_V, &
         npart, npartoftype, xyzh, vxyzu, itype, x0, v0, isink, JKmuS)
 else
    call inject_geodesic_sphere(i, ifirst, ires, r, v, u, rho,  geodesic_R, geodesic_V, &
         npart, npartoftype, xyzh, vxyzu, itype, x0, v0, isink)
 endif
 if (isink_radiation > 0) dust_temp(ifirst:ifirst+particles_per_sphere-1) = xyzmh_ptmass(iTeff,isink)

end subroutine inject_sphere

!-----------------------------------------------------------------------
!+
!  Updates the injected particles
!+
!-----------------------------------------------------------------------
subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  compute 1D wind profile and define the number of shells to fill the domain
!+
!-----------------------------------------------------------------------
subroutine set_1D_wind_profile (isink,d_part,rinject,time_between_spheres,tboundary,tcross,tfill,onewind)

 use units,     only:udist,utime,unit_velocity
 use physcon,   only:au
 use io,        only:fileprefix
 use timestep,  only:tmax
 use wind,      only:save_windprofile

 integer, intent(in) :: isink
 real, intent(in)    :: rinject,time_between_spheres,d_part
 logical, intent(in) :: onewind
 real, intent(out)   :: tboundary,tcross,tfill
 real :: tend
 character(len=24) :: wfile

 tboundary = (iboundary_spheres+nfill_domain)*time_between_spheres
 tend      = max(tmax,tboundary)*utime
 write(wfile,'(a,i2.2,"_profile.dat")') trim(fileprefix),isink
 call save_windprofile(rinject*udist,wind_injection_speed*unit_velocity,&
      wind_temperature,outer_boundary_au*au,rfill_domain_au*au,tend,tcross,tfill,wfile,isink)
 if (tboundary > tmax) then
    print *,'simulation time < time to reach the last boundary shell'
 endif

!define the number of background shells
 if (tfill < 1.d98 .and. rfill_domain_au > 1e-5 .and. (isink == 2 .or. onewind)) then
    nfill_domain = int(tfill/(utime*time_between_spheres))-iboundary_spheres
    print *,'number of background shells set to',nfill_domain
    if (nfill_domain < 0 ) then
       call logging(isink,rinject,time_between_spheres,d_part,tboundary=tboundary,tcross=tcross,tfill=tfill)
       print *,'[stop set_1D_wind_profile]',tfill,utime*time_between_spheres,iboundary_spheres
       stop
    endif
 endif

end subroutine set_1D_wind_profile

!-----------------------------------------------------------------------
!+
!  initialize oscillating inner boundary
!+
!-----------------------------------------------------------------------
subroutine init_pulsating_wind(pulsating_wind)

 use units,   only:unit_velocity,utime,unit_velocity
 use physcon, only:pi,days,km
 logical, intent(in) :: pulsating_wind

 if (pulsating_wind) then
    pulsation_period = pulsation_period_days * (days/utime)
    dtpulsation      = pulsation_period/50.
    omega_puls       = 2.*pi/pulsation_period
    deltaR_puls      = pulsation_period*piston_velocity/(2.*pi)
    piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
 else
    omega_puls       = 0.d0
    deltaR_puls      = 0.d0
    piston_velocity  = 0.d0
 endif

end subroutine init_pulsating_wind

!-----------------------------------------------------------------------
!+
!  Oscillating inner boundary
!+
!-----------------------------------------------------------------------
subroutine pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,sphere_number,inner_sphere)
 use physcon, only:pi
 use eos,     only:gamma
 use units,   only:unit_velocity
 use part,    only:massoftype,xyzmh_ptmass,ieject,igas,imloss

 integer, intent(in)  :: sphere_number, inner_sphere
 real,    intent(in)  :: time,local_time,GM
 real,  intent(inout) :: r
 real,    intent(out) :: v, u, rho, e

 integer, parameter :: nrho_index = 10
 integer :: k,isink=1
 real :: inner_radius,r3,mass_of_spheres,rho_ini
 logical :: verbose = .true.

 mass_of_spheres = massoftype(igas) * xyzmh_ptmass(ieject,isink)
 rho_ini = xyzmh_ptmass(imloss,isink)/(4.*pi*r**2*wind_injection_speed)
 dr3 = 3.*mass_of_spheres/(4.*pi*rho_ini)
 v = wind_velocity + piston_velocity* cos(omega_puls*time) !same velocity for all wall particles
 inner_radius = wind_injection_radius + deltaR_puls*sin(omega_puls*time)
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
!  Sets default options for the injection module
!+
!-----------------------------------------------------------------------
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
 use dim,          only:maxvxyzu
 use infile_utils, only:write_inopt
 use part,         only:nptmass
 integer, intent(in) :: iunit

 call write_inopt(wind_type,'wind_type','wind type 0=std, 1=transonic, 2=pulsation, 3=jets',iunit)
 if (wind_type == 2) then
    call write_inopt(pulsation_period_days,'pulsation_period','stellar pulsation period (days)',iunit)
    call write_inopt(piston_velocity_km_s,'piston_velocity','velocity amplitude of the pulsation (km/s)',iunit)
 endif
 if (wind_type == 3) then
    call write_inopt(jet_edge_velocity,'jet_edge_velocity','velocity at the edge of the jet (km/s, only for sink1)',iunit)
    call write_inopt(jet_opening_angle_degree,'jet_opening_angle','half opening angle of the jet (degree)',iunit)
 endif
 ! if isink > 1 effective radius should be used as injection radius
 if (nptmass < 2) then
    call write_inopt(wind_injection_radius_au,'wind_inject_radius','wind injection radius (au, if 0 takes Rstar)',iunit)
 endif
 call write_inopt(iwind_resolution,'iwind_resolution','set number of ejected particles, (icosahedron < 0, fibonacci > 0',iunit)
 call write_inopt(rfill_domain_au,'rfill_domain','outer radius of the background density profile',iunit)
 call write_inopt(wind_shell_spacing,'wind_shell_spacing','desired ratio of sphere spacing to particle spacing',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres','number of boundary spheres (integer)',iunit)
 call write_inopt(outer_boundary_au,'outer_boundary','delete gas particles outside this radius (au)',iunit)
 if (mhd) call write_inopt(B_r,'B_r','radial magnetic field strength (G)',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use dim,          only:maxvxyzu
 use physcon,      only:deg_to_rad
 use injectutils,  only:use_icosahedron
 use part,         only:nptmass
 use io,           only:warning
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 logical, save :: init_opt = .false.
 character(len=*), parameter :: label='read_options'

 if (.not.init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif
 call read_inopt(wind_type,'wind_type',db,errcount=nerr,min=0,max=3)
 if (wind_type == 2) then
    call read_inopt(pulsation_period_days,'pulsation_period',db,errcount=nerr,min=0.)
    call read_inopt(piston_velocity_km_s,'piston_velocity',db,errcount=nerr,min=0.)
 endif
 if (wind_type == 3) then
    call read_inopt(jet_edge_velocity,'jet_edge_velocity',db,errcount=nerr,min=0.)
    call read_inopt(jet_opening_angle_degree,'jet_opening_angle',db,errcount=nerr,min=0.,max=90.)
    jet_opening_angle = jet_opening_angle_degree*deg_to_rad
 endif
 ! if isink > 1 effective radius should be used as injection radius
 if (nptmass < 2) then
    call read_inopt(wind_injection_radius_au,'wind_inject_radius',db,errcount=nerr,min=0.)
 endif
 call read_inopt(iwind_resolution,'iwind_resolution',db,errcount=nerr,default=iwind_resolution)
 call read_inopt(rfill_domain_au,'rfill_domain',db,errcount=nerr,min=0.)
 call read_inopt(wind_shell_spacing,'wind_shell_spacing',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(iboundary_spheres,'iboundary_spheres',db,errcount=nerr,min=0)
 call read_inopt(outer_boundary_au,'outer_boundary',db,errcount=nerr,min=0.)
 if (mhd) call read_inopt(B_r,'B_r',db,errcount=nerr)
 if (iwind_resolution < 0) then
    use_icosahedron = .true.
 else
    use_icosahedron = .false.
    if (iwind_resolution < 15) call warning(label,'resolution likely too low for fibonacci (>15)')
 endif

end subroutine read_options_inject

end module inject
