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
!   infile_utils, injectutils, io, options, part, partinject, physcon,
!   ptmass_radiation, setbinary, timestep, units, wind, wind_equations
!
 use dim, only:isothermal,nucleation

 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      wind_injection_radius,set_default_options_inject,update_injected_par,&
      write_headeropts_inject_wind,read_headeropts_inject_wind

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
 real :: wind_injection_radius_au
 real :: wind_temperature
 real :: outer_boundary_au = 30.
 real :: wind_shell_spacing = 1.
 real :: pulsation_period
 real :: pulsation_period_days = 0.
 real :: rfill
 real :: rfill_domain_au = 0.
 real :: piston_velocity_km_s = 0.
 real :: dtpulsation = huge(0.)
 real :: jet_edge_velocity = 0.
 real :: jet_opening_angle = 0.
 real :: jet_opening_angle_degree = 0.

! global variables
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: u_to_temperature_ratio,wind_mass_rate,rho_ini,&
      mass_of_spheres,time_between_spheres,neighbour_distance,wind_velocity,&
      Rstar_cgs,Mstar_cgs,wind_injection_radius,wind_injection_speed,&
      rinject,initial_wind_velocity_cgs
 real :: omega_puls,deltaR_puls,piston_velocity,dr3 !pulsations
 integer :: particles_per_sphere,iresolution,nwrite

 logical :: pulsating_wind
 character(len=*), parameter :: label = 'inject_wind'

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use options,           only:icooling,ieos,iwind_res
 use io,                only:fatal
 use setbinary,         only:get_eccentricity_vector
 use wind_equations,    only:init_wind_equations
 use wind,              only:setup_wind
 use physcon,           only:mass_proton_cgs, kboltz, Rg, days, km, au, years, solarm, pi, Gg
 use icosahedron,       only:compute_matrices, compute_corners
 use eos,               only:gmw,gamma,polyk
 use units,             only:unit_velocity,umass,utime,udist
 use part,              only:xyzmh_ptmass,vxyz_ptmass,igas,imloss,ilum,iTeff,iReff,nptmass
 use injectutils,       only:get_parts_per_sphere,get_neighb_distance,use_icosahedron

 integer, intent(out) :: ierr
 integer :: isink = 1
 real :: semimajoraxis_cgs,orbital_period
 real :: tcross,rsonic,tsonic,tfill,initial_rinject,tboundary
 real :: separation_cgs,ecc(3),eccentricity

! use_icosahedron = .false.
! if (iwind_resolution > 15 .and. use_icosahedron) use_icosahedron = .false.
 use_icosahedron = .true.
 if (use_icosahedron) then
    call compute_matrices(geodesic_R)
    call compute_corners(geodesic_v)
 endif

 if (icooling > 0) nwrite = nwrite+1
 ierr = 0

 pulsating_wind = (pulsation_period_days > 0.) .and. (piston_velocity_km_s > 0.)
 if (pulsating_wind .and. ieos == 6) call fatal(label,'cannot use ieos=6 with pulsation')

! setup thermo
 if (gamma > 1.0001) then
    u_to_temperature_ratio = Rg/(gmw*(gamma-1.)) / unit_velocity**2
 else
    u_to_temperature_ratio = Rg/(gmw*2./3.) / unit_velocity**2
 endif
 if (isothermal) wind_temperature = polyk * mass_proton_cgs/kboltz * unit_velocity**2*gmw

 call init_pulsating_wind(pulsating_wind)

 write(*,'(/,70("-"))')
 Rstar_cgs = xyzmh_ptmass(iReff,isink)*udist
 Mstar_cgs = xyzmh_ptmass(4,isink)*umass
 if (wind_temperature < 0.001)  wind_temperature = xyzmh_ptmass(iTeff,isink)

 wind_velocity  = wind_velocity_km_s * (km / unit_velocity)
 wind_mass_rate = wind_mass_rate_Msun_yr * (solarm/umass) / (years/utime)
 xyzmh_ptmass(imloss,isink) = wind_mass_rate

 if (wind_injection_radius_au == 0.) wind_injection_radius_au = Rstar_cgs/au
 wind_injection_radius = wind_injection_radius_au * au / udist

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
    call init_wind_equations (xyzmh_ptmass(4,isink), &
         xyzmh_ptmass(iTeff,isink), u_to_temperature_ratio)
!if (ieos == 6 .and. wind_temperature /= star_Teff) then
!    wind_injection_radius_au = Rstar_cgs/au*(star_Teff/wind_temperature)**(1./wind_expT)
!    wind_injection_radius  = wind_injection_radius_au * au / udist
!endif
!wind_injection_radius = max(wind_injection_radius_au * au, Rstar_cgs) / udist
!Rstar = min(wind_injection_radius_au*au,Rstar_cgs)

! integrate wind equation to get initial velocity and sonic radius to set resolution
    rinject = initial_rinject*udist
    call setup_wind(Mstar_cgs, wind_mass_rate, u_to_temperature_ratio, rinject,&
         wind_temperature,initial_wind_velocity_cgs,rsonic,tsonic,wind_type)
    initial_rinject = rinject/udist
 endif

 rinject = initial_rinject
 wind_injection_speed = initial_wind_velocity_cgs/unit_velocity
 rho_ini = wind_mass_rate/(4.*pi*rinject**2*wind_injection_speed)

 ! set wind resolution parameters
 if (isink == 1) call init_resolution(rsonic,tsonic)

 ! compute 1D wind profile to get tcross & save 1D profile
 if ( .not. pulsating_wind .or. rfill_domain_au > 0.) then
    call set_1D_wind_profile (tboundary,tcross,tfill)
 else
    tfill = 0.
 endif

! logging
 call logging(rsonic,tsonic,tboundary,tcross,tfill)

 ! if wind resolution has changed, rearrange memory location of boundary particles
 if (iwind_resolution /= iwind_res .and. .not.use_icosahedron) then
    call adjust_boundary_shells(iwind_res)
 endif

 if (nptmass == 2) then
    separation_cgs = sqrt(sum((xyzmh_ptmass(1:3,2)-xyzmh_ptmass(1:3,1))**2))*udist
    ecc             = get_eccentricity_vector(xyzmh_ptmass,vxyz_ptmass,1,2)
    eccentricity    = sqrt(sum(ecc(1:3)**2))
    !stars initially positioned at apastron (see set_binary)
    semimajoraxis_cgs = separation_cgs/(1.+eccentricity)
    orbital_period  = sqrt(4.*pi**2*semimajoraxis_cgs**3/(Gg*(xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2))*solarm))    ! cgs
    print *,'eccentricity                 = ',eccentricity
    print *,'orbital_period (days)        = ',orbital_period/86400.
    print *,'orbital speed sink 1 (km/s)  = ',sqrt(dot_product(vxyz_ptmass(1:3,1),vxyz_ptmass(1:3,1)))*unit_velocity/km
    print *,'orbital speed sink 2 (km/s)  = ',sqrt(dot_product(vxyz_ptmass(1:3,2),vxyz_ptmass(1:3,2)))*unit_velocity/km
 endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  adjust the memory location of the boundary particles
!+
!-----------------------------------------------------------------------
subroutine adjust_boundary_shells(iwind_res)

 use part, only:iphase,xyzh,vxyzu,iorig,fxyzu,fext,npartoftype,npart,iboundary,igas,nucleation
 use dim,  only:do_nucleation
 use io,   only:iverbose

 integer, intent(in) :: iwind_res
 integer :: i,istart,iend,incr,nshift

 !the old boundary particles are stored between [1,iwind_res*iboundary_spheres]
 !the new boundary particles are stored between [1,iwind_resolution*iboundary_spheres]

 if (iwind_resolution > iwind_res) then
    !free space to insert new boundary particles
    nshift = (iwind_resolution-iwind_res)*iboundary_spheres
    istart = npart
    iend   = iwind_res*iboundary_spheres+1
    incr   = -1
 else
    nshift = (iwind_resolution-iwind_res)*iboundary_spheres
    iend   = npart
    istart = iwind_res*iboundary_spheres+1
    incr   = 1
 endif
 if (iverbose > 0) then
    print *,'iboundary_spheres =',iboundary_spheres,', incr=',incr
    print *,'iwind_resolution  =',iwind_res,'-->',iwind_resolution
    print *,'move boundary nold=',iwind_res*iboundary_spheres,'-->',iwind_resolution*iboundary_spheres
    print *,'shifting particles [',istart,',',iend,'] --> [',istart+nshift,',',iend+nshift,'], nshift =',nshift
    print *,npartoftype
 endif
 do i = istart,iend,incr
!print *,i,j,iwind_resolution*iboundary_spheres,npart-iwind_res*iboundary_spheres
    iorig(i+nshift) = iorig(i)
    xyzh(:,i+nshift) = xyzh(:,i)
    vxyzu(:,i+nshift) = vxyzu(:,i)
    fxyzu(:,i+nshift) = fxyzu(:,i)
    fext(:,i+nshift)  = fext(:,i)
    iphase(i+nshift) = iphase(i)
    if (do_nucleation) then
       nucleation(:,i+nshift) = nucleation(:,i)
    endif
 enddo
 if (iwind_resolution*iboundary_spheres > npart .and. incr <0) iphase(npart:iwind_resolution*iboundary_spheres) = 1
 npart = npart+nshift
 npartoftype(igas) = npart
 !npartoftype(iboundary) = iwind_resolution*iboundary_spheres
 !print *,npart,npartoftype

end subroutine adjust_boundary_shells

!-------------------------------------------------------------------------------
!+
!  set particle mass or resolution depending on iwind_resolution (set by sink 1)
!+
!-------------------------------------------------------------------------------
subroutine init_resolution(rsonic,tsonic)

 use part,        only:xyzmh_ptmass,massoftype,igas,iboundary,imloss,npart
 use units,       only:udist
 use timestep,    only:tmax
 use physcon,     only:pi
 use injectutils, only:get_sphere_resolution,get_parts_per_sphere,get_neighb_distance
 use io,          only:fatal

 real, intent(in) :: rsonic,tsonic

 integer :: nzones_per_sonic_point,isink=1
 real :: mV_on_MdotR,check_mass,dr,dist_to_sonic_point,mass_of_particles

 check_mass = massoftype(igas)
 if (iwind_resolution == 0) then
    !
    ! resolution is specified in terms of number of smoothing lengths
    ! per distance to sonic point (if trans-sonic wind)
    !
    if (wind_type == 1) then
       nzones_per_sonic_point = 8
       dist_to_sonic_point = rsonic/udist-rinject
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
    mV_on_MdotR = mass_of_particles*wind_injection_speed/(xyzmh_ptmass(imloss,isink)*rinject)
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
    mass_of_particles    = wind_shell_spacing*neighbour_distance*rinject*&
         xyzmh_ptmass(imloss,isink)/(particles_per_sphere * wind_injection_speed)
    massoftype(igas)     = mass_of_particles
 endif

 !spheres properties
 mass_of_spheres = massoftype(igas) * particles_per_sphere
 time_between_spheres = mass_of_spheres/wind_mass_rate
 dr3 = 3.*mass_of_spheres/(4.*pi*rho_ini)

 if (npart > 0 .and. abs(log10(check_mass/mass_of_particles)) > 1e-10) then
    call logging(rsonic,tsonic)
    print *,'check_mass    = ',check_mass
    print *,'particle mass = ',mass_of_particles
    print *,'number of particles = ',npart,abs(1.-check_mass/mass_of_particles)*100
    call fatal(label,"you cannot reset the particle's mass")
 endif
 massoftype(iboundary) = massoftype(igas)

 if (time_between_spheres > tmax)  then
    call logging(rsonic,tsonic)
    print *,'time_between_spheres = ',time_between_spheres,' < tmax = ',tmax
    call fatal(label,'no shell ejection : tmax < time_between_spheres')
 endif


end subroutine init_resolution


!-----------------------------------------------------------------------

subroutine logging(rsonic,tsonic,tboundary,tcross,tfill)

!-----------------------------------------------------------------------

 use physcon,           only:pi,gg,au,km
 use units,             only:udist,unit_velocity,utime
 use timestep,          only:dtmax
 use ptmass_radiation,  only:alpha_rad
 use part,              only:massoftype,igas,xyzmh_ptmass,iReff,&
                             ispinx,ispiny,ispinz,imloss
 use injectutils,       only:use_icosahedron

 real, intent(in) :: rsonic,tsonic
 real, optional, intent(in) :: tboundary,tcross,tfill
 integer :: ires_min,isink = 1
 integer(kind=8) :: ninject
 real :: vesc,wind_rotation_speed,rotation_speed_crit

 vesc = sqrt(2.*Gg*Mstar_cgs*(1.-alpha_rad)/Rstar_cgs)
 write (*,'(/,2(3x,A,es11.4),3x,A,l1)')&
      'mass_of_particles       : ',massoftype(igas),&
      'time_between_spheres    : ',time_between_spheres,&
      'use_icosahedron         : ',use_icosahedron
 write (*,'(2(3x,A,es11.4),3x,A,i7)') &
      'dist between spheres    : ',wind_shell_spacing*neighbour_distance,&
      'dist between injection  : ',time_between_spheres*wind_injection_speed,&
      'particles per sphere    : ',particles_per_sphere
 write (*,'(3(3x,A,es11.4))') &
      'wind_temperature        : ',wind_temperature,&
      'injection_radius   (au) : ',rinject*au/udist
 write (*,'(2(3x,A,es11.4))') &
      'stellar_radius     (au) : ',Rstar_cgs/udist, &
      'rho_ini           (cgs) : ',rho_ini
 if (present(tcross)) then
    write (*,'(3(3x,A,es11.4))') &
      'crossing time      (cu) : ',tcross,&
      'filling time       (cu) : ',tfill/utime,&
      'boundary time      (cu) : ',tboundary
 endif
 ninject = particles_per_sphere*(nfill_domain+iboundary_spheres)
 if (nfill_domain > 0) print '(3x,"number of background particle = ",i10)',ninject

 !print*,'hmax/dist_between_spheres  = ',wind_shell_spacing*neighbour_distance*&
 !      initial_wind_velocity_cgs**2/(vesc**2-initial_wind_velocity_cgs**2)
 if (wind_type == 1) then
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
 elseif (rsonic/udist > rinject) then
    !save a few models before the particles reach the sonic point
    if (dtmax > tsonic/utime) print *,'WARNING! dtmax > time to sonic point'
    !if solution subsonic, minimum resolution required so a few shells can be inserted between the injection radius and the sonic point
    ires_min = int(iboundary_spheres*wind_shell_spacing*0.5257/(rsonic/udist/rinject-1.)+.5)
    if (iwind_resolution < ires_min) print *,'WARNING! resolution too low to pass sonic point : iwind_resolution < ',ires_min
 endif

end subroutine logging


!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject,is_sink)
 use physcon,           only:pi,au,solarm,years
 use io,                only:fatal,iverbose
 use wind,              only:interp_wind_profile !,wind_profile
 use part,              only:igas,iTeff,iReff,iboundary,nptmass,delete_particles_outside_sphere,&
                             delete_dead_particles_inside_radius,n_nucleation,imloss
 use partinject,        only:add_or_update_particle
 use injectutils,       only:use_icosahedron
 use units,             only:udist,umass,utime
 use dust_formation,    only:idust_opacity
 use options,           only:iwind_res

 integer, intent(in),optional :: is_sink
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere,inner_sphere,inner_boundary_sphere,ifirst,i,ipart, &
            nreleased,nboundaries,isink,itype,ires
 real    :: local_time,GM,r,v,u,rho,e,mass_lost,x0(3),v0(3),inner_radius,fdone,dum
 character(len=*), parameter :: label = 'inject_particles'
 logical, save :: released = .false.
 real :: JKmuS(n_nucleation)

 dum  = 0.
 if (present(is_sink)) then
    isink = is_sink
 else
    isink = 1
 endif

 if (nptmass > 0) then
    x0 = xyzmh_ptmass(1:3,isink)
    GM = xyzmh_ptmass(4,isink)
    v0 = vxyz_ptmass(1:3,isink)
 else
    x0 = 0.
    v0 = 0.
    GM = 0.
 endif

 if (use_icosahedron) then
    ires = iwind_resolution
 else
    ires = particles_per_sphere
 endif

 if (npart > 0 .and. iwind_resolution == iwind_res) then
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
    if (.not.released) then
       if (nfill_domain > 0) print '(/,A,i7,3x,"Mdot : ",es12.4)','background particles :',&
            nreleased*particles_per_sphere,xyzmh_ptmass(imloss,isink)/(solarm/umass)*(years/utime)
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
 else
    !initialise domain with boundary particles
    ipart = iboundary
    nreleased = 0
    nboundaries = iboundary_spheres+nfill_domain
    iwind_res = iwind_resolution
 endif

 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1 + nreleased
 inner_sphere = floor(time/time_between_spheres)+nreleased
 inner_boundary_sphere = inner_sphere + nboundaries
 !print '("@@@",4(i4),5(1x,es12.4))',inner_boundary_sphere,nreleased,nboundaries,wind_type,time,time_between_spheres

 if (nptmass < 2) then
    ! only one sphere can be ejected at a time (not valid for 2 wind-emitting sinks)
    if (inner_sphere-outer_sphere > nboundaries) call fatal(label,'ejection of more than 1 sphere, timestep likely too large!')
 endif

 !
 ! eject particles and update properties of boundary particles
 !
 do i=inner_boundary_sphere,outer_sphere,-1
    local_time = time + (iboundary_spheres+nfill_domain-i) * time_between_spheres
    !compute the radius, velocity, temperature, chemistry of a shell at the current local time
    v = wind_injection_speed
    r = rinject
    if (pulsating_wind.and.released) then
       call pulsating_wind_profile(time,local_time,r,v,u,rho,e,GM,i,inner_sphere,dr3,rho_ini)
    else
       if (idust_opacity == 2) then
          call interp_wind_profile(time,local_time,r,v,u,rho,e,GM,fdone,JKmuS)
       else
          call interp_wind_profile(time,local_time,r,v,u,rho,e,GM,fdone)
       endif
       if (iverbose > 0) print '(" ## update boundary ",i4,2(i4),i7,8(1x,es12.5))',i,&
            inner_sphere,outer_sphere,npart,time,local_time,r/xyzmh_ptmass(iReff,isink),v*udist/utime,&
            xyzmh_ptmass(imloss,isink) /(solarm/umass) * (years/utime)
    endif

    if (i > inner_sphere) then
       ! boundary sphere
       ifirst = (nboundaries-i+inner_sphere)*particles_per_sphere+1
       itype = ipart
       if (iverbose > 0) print '(" @@ update boundary ",i4,3(i4),i7,8(1x,es12.5))',i,inner_sphere,&
            outer_sphere,ipart,ifirst,time,local_time,r/xyzmh_ptmass(iReff,isink),v*udist/utime,u,rho
    else
       ! ejected particles + create new  inner sphere
       ifirst = npart+1
       itype = igas
       if (iverbose > 0) print '(" ## eject particles [",i7,"-",i7,"] ",3(i4),10(1x,es11.4))',&
            npart+1-particles_per_sphere,npart,i,inner_sphere,outer_sphere,time,local_time,&
            r/xyzmh_ptmass(iReff,isink),v,u,rho,xyzmh_ptmass(imloss,isink)/(solarm/umass)*(years/utime)
    endif
    if (idust_opacity == 2) then
       call inject_sphere(i,ifirst,ires,r,v,u,rho,npart,npartoftype,xyzh,vxyzu,itype,x0,v0,isink,JKmuS)
    else
       call inject_sphere(i,ifirst,ires,r,v,u,rho,npart,npartoftype,xyzh,vxyzu,itype,x0,v0,isink)
    endif
    !cs2max = max(cs2max,gamma*(gamma-1)*u)
 enddo

 if (nfill_domain > 0 .and. outer_sphere == 1) print '(3x,"injecting background particles up to r = ",f8.2,&
      &" (au), nparticles = ",i8)',r,particles_per_sphere*(nfill_domain+iboundary_spheres)

 ! update sink particle properties
 if (nptmass > 0) then
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    xyzmh_ptmass(4,isink) = xyzmh_ptmass(4,isink) - mass_lost
    if (pulsating_wind) then
       inner_radius = wind_injection_radius + deltaR_puls*sin(omega_puls*time)
       !v2 xyzh_ptmass(5,isink) = (inner_radius**3-dr3)**(1./3.) !accretion radius
       xyzmh_ptmass(5,isink) = inner_radius
    endif
 endif

 !
 ! return timestep constraint to ensure that time between sphere
 ! injections is adequately resolved
 !
 !dr = neighbour_distance*wind_injection_radius
 !dtinject = 0.25*dr/sqrt(cs2max)
 dtinject = min(0.2*time_between_spheres,dtpulsation)
 if (time <= 0.) dtinject = 0.01*dtinject

end subroutine inject_particles


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

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  compute 1D wind profile and define the number of shells to fill the domain
!+
!-----------------------------------------------------------------------
subroutine set_1D_wind_profile (tboundary,tcross,tfill)

 !use part,      only:xyzmh_ptmass,vxyz_ptmass,massoftype,igas,iboundary,imloss,ilum,iTeff,iReff,nptmass,npart
 use units,     only:udist,utime,unit_velocity
 use physcon,   only:au
 use timestep,  only:tmax
 use wind,      only:save_windprofile

 real, intent(out)   :: tboundary,tcross,tfill
 real :: tend

 tboundary = (iboundary_spheres+nfill_domain)*time_between_spheres
 tend      = max(tmax,tboundary)*utime
 call save_windprofile(rinject*udist,wind_injection_speed*unit_velocity,&
      wind_temperature,outer_boundary_au*au,rfill,tend,tcross,tfill,'wind_profile1D.dat')
 if (tboundary > tmax .and. nfill_domain > 1) then
    print *,'simulation time < time to reach the last boundary shell'
 endif

!define the number of background shells
 if (tfill < 1.d98 .and. rfill_domain_au > 1e-5) then
    nfill_domain = int(tfill/(utime*time_between_spheres))-iboundary_spheres
    print *,'number of background shells set to',nfill_domain
    if (nfill_domain < 0 ) then
       print *,tfill,utime*time_between_spheres,iboundary_spheres
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


subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
 integer :: icase

 if (present(flag)) then
    icase = flag
 else
    icase = 0
 endif

 if (isothermal) then
    wind_type = 1
    wind_velocity_km_s = 0.
    wind_mass_rate_Msun_yr = 8.2d-8
    wind_injection_radius_au = 0.
 else
    if (icase == 1) then
       !trans-sonic wind
       wind_type = 1
       wind_velocity_km_s = 0.
       wind_mass_rate_Msun_yr = 1.d-5
       wind_injection_radius_au = 2.
       wind_temperature = 50000.
    else
       !super sonic-wind
       wind_type = 0
       wind_velocity_km_s = 20.
       wind_mass_rate_Msun_yr = 1.d-5
       wind_injection_radius_au = 2.
       wind_temperature = 2500.
    endif
 endif

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_inject_wind(hdr,ierr)
 use dump_utils,        only:dump_h,add_to_iheader
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out)   :: ierr

 ierr = 0
 call add_to_iheader(iwind_resolution,'iwind_res',hdr,ierr) ! integer

end subroutine write_headeropts_inject_wind

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_inject_wind(hdr,ierr)
 use options,    only:iwind_res
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr

 ierr = 0
 call extract('iwind_res',iwind_res,hdr,ierr) ! integer

end subroutine read_headeropts_inject_wind

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use dim,          only: maxvxyzu
 use infile_utils, only: write_inopt
 use physcon,      only: au
 integer, intent(in) :: iunit

 call write_inopt(wind_type,'wind_type','wind type 0=std, 1=transonic, 2=pulsation',iunit)
 call write_inopt(wind_velocity_km_s,'wind_velocity','injection wind velocity (km/s, if wind_type = 0)',iunit)
 if (wind_type == 2) then
    call write_inopt(pulsation_period_days,'pulsation_period','stellar pulsation period (days)',iunit)
    call write_inopt(piston_velocity_km_s,'piston_velocity','velocity amplitude of the pulsation (km/s)',iunit)
 endif
 call write_inopt(wind_injection_radius_au,'wind_inject_radius','wind injection radius (au, if 0 takes Rstar)',iunit)
 call write_inopt(wind_mass_rate_Msun_yr,'wind_mass_rate','wind mass loss rate (Msun/yr)',iunit)
 if (maxvxyzu==4) then
    call write_inopt(wind_temperature,'wind_temperature','wind temperature at injection radius (K, if 0 takes Teff)',iunit)
 endif
 call write_inopt(iwind_resolution,'iwind_resolution','if<>0 set number of particles on the sphere, reset particle mass',iunit)
 call write_inopt(rfill_domain_au,'rfill_domain','outer radius of the background density profile',iunit)
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
 use physcon, only:au

 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 integer :: noptions
 logical :: init_opt = .false. ! isowind = .true.
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
    if (wind_temperature < 0.)    call fatal(label,'invalid setting for wind_temperature (<0)')
 case('iwind_resolution')
    read(valstring,*,iostat=ierr) iwind_resolution
    ngot = ngot + 1
    if (iwind_resolution < 0) call fatal(label,'iwind_resolution must be bigger than zero')
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    ngot = ngot + 1
    if (iboundary_spheres <= 0) call fatal(label,'iboundary_spheres must be > 0')
 case('rfill_domain')
    read(valstring,*,iostat=ierr) rfill_domain_au
    ngot = ngot + 1
    if (rfill_domain_au < 0) call fatal(label,'rfill_domain must be > 0')
    rfill = rfill_domain_au*au
 case('wind_shell_spacing')
    read(valstring,*,iostat=ierr) wind_shell_spacing
    ngot = ngot + 1
    if (wind_shell_spacing <= 0.) call fatal(label,'wind_shell_spacing must be >=0')
 case('wind_type')
    read(valstring,*,iostat=ierr) wind_type
    ngot = ngot + 1
    if (wind_type < 0 .or. wind_type > 2 ) call fatal(label,'invalid setting for wind_type ([0,2])')
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
