!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
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
!   - jet_edge_velocity  : *velocity at the edge of the jet (km/s, only for sink1)*
!   - jet_opening_angle  : *half opening angle of the jet (degree)*
!   - outer_boundary     : *delete gas particles outside this radius (au)*
!   - piston_velocity    : *velocity amplitude of the pulsation (km/s)*
!   - pulsation_period   : *stellar pulsation period (days)*
!   - rfill_domain       : *outer radius of the background density profile*
!   - wind_shell_spacing : *desired ratio of sphere spacing to particle spacing*
!   - wind_type          : *wind type 0=std, 1=transonic, 2=pulsation, 3=jets*
!
! :Dependencies: dim, dust_formation, eos, infile_utils, injectutils, io,
!   options, orbits, part, partinject, physcon, ptmass_radiation, timestep,
!   units, wind, wind_equations
!
 use dim, only:isothermal,mhd
 use wind, only:rfill_domain_au,wind_shell_spacing

 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par,print_wind_resolution

 private
!
!--runtime settings for this module
!
! Read from input file
 integer :: wind_type = 0
 integer :: iboundary_spheres = 5
 integer :: nfill_domain = 0
 real :: outer_boundary_au = 30.
 real :: pulsation_period
 real :: pulsation_period_days = 0.
 real :: piston_velocity_km_s = 0.
 real :: dtpulsation = huge(0.)
 real :: jet_edge_velocity = 0.
 real :: jet_opening_angle = 0.
 real :: jet_opening_angle_degree = 0.
 real :: B_r = 0.

! global variables
 real :: omega_puls,deltaR_puls,piston_velocity !pulsations
 integer :: nwrite

 logical :: pulsating_wind,onewind
 character(len=*), parameter :: label = 'inject_wind'
 ! fewest particles per shell for the wind to be resolved azimuthally (icosahedron)
 integer, parameter :: min_particles_per_shell = 12

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
 use wind,              only:wind_params,get_params_from_sink,solve_wind
 use physcon,           only:mass_proton_cgs,kboltz,km
 use eos,               only:gmw,polyk
 use units,             only:unit_velocity,udist,get_nearest_unit
 use part,              only:xyzmh_ptmass,vxyz_ptmass,imloss,iReff,ivwind,iTwind,nptmass,&
                              massoftype,igas,iboundary
 use injectutils,       only:init_jets,seed_random

 integer, intent(out) :: ierr
 integer :: isink
 real :: time_between_spheres,d_part
 real :: tcross,rsonic,tsonic,tfill,tboundary,mu,wind_temperature
 real :: dx(3),dv(3),ecc(3)
 integer :: nwinds
 type(wind_params) :: params !reset wind temperature if init_muGamma

 nwinds = count(abs(xyzmh_ptmass(imloss,1:nptmass)) > tiny(0.))
 if (nwinds > 2) call fatal(label,'cannot handle more than 2 emitting sinks (1D wind profile not available)')
 onewind = (nwinds <= 1)
 ! change particle injection method if more than 1 sink is emitting a wind
 if (nwinds > 1) then
    if (wind_type == 3) call init_jets(jet_edge_velocity,jet_opening_angle)
 endif
 if (abs(xyzmh_ptmass(imloss,1)) < tiny(0.)) call fatal(label,'the wind logic imposes that sink 1 is loosing mass')

 if (icooling > 0) nwrite = nwrite+1
 ierr = 0

 seed_random = -1  ! reset seed_random to avoid reproducibility issues

 pulsating_wind = (pulsation_period_days > 0.) .and. (piston_velocity_km_s > 0.)
 if (ieos == 6) call fatal(label,'cannot use ieos=6 with pulsation')

 if (isothermal) wind_temperature = polyk * mass_proton_cgs/kboltz * unit_velocity**2*gmw

 call init_pulsating_wind(pulsating_wind)

 ! the particle mass is set in the setup and is never changed here
 massoftype(iboundary) = massoftype(igas)

 write(*,'(/,70("-"))')
 do isink = 1, nptmass

    ! skip sink if no mass loss
    if (abs(xyzmh_ptmass(imloss,isink)) < tiny(0.)) cycle

    call get_params_from_sink(xyzmh_ptmass(:,isink),wind_type,params)
    ! isothermal: the wind is at the input temperature set by polyk, not Teff
    if (isothermal) params%Twind = wind_temperature

    !implement background spheres starting from the smallest radius if pulsating
    !params%Rinject = min(params%Rinject,(xyzmh_ptmass(iReff,isink)-deltaR_puls)*udist)
    if (.not.pulsating_wind) call solve_wind(xyzmh_ptmass(:,isink),wind_type,params,rsonic,tsonic)

    if (xyzmh_ptmass(ivwind,isink) <= 0.) xyzmh_ptmass(ivwind,isink) = params%vwind/unit_velocity
    if (xyzmh_ptmass(iTwind,isink) <= 0.) xyzmh_ptmass(iTwind,isink) = params%Twind
    xyzmh_ptmass(iReff,isink) = params%Rinject/udist

    call init_sink_resolution(isink,time_between_spheres,d_part)

    ! compute 1D wind profile to get tcross & save 1D profile
    if (.not. pulsating_wind .or. rfill_domain_au > 0.) then
       call set_1D_wind_profile(params,isink,d_part,time_between_spheres,tboundary,tcross,tfill,onewind)
    endif

! logging
    if (xyzmh_ptmass(imloss,isink) > 0.) &
       call logging(params,isink,time_between_spheres,d_part,rsonic,tsonic,tboundary,tcross,tfill)

 enddo

 if (nptmass == 2) then
    ! period from the relative position and velocity, valid at any orbital phase
    dx  = xyzmh_ptmass(1:3,2)-xyzmh_ptmass(1:3,1)
    dv  = vxyz_ptmass(1:3,2)-vxyz_ptmass(1:3,1)
    mu  = xyzmh_ptmass(4,1)+xyzmh_ptmass(4,2)   ! G*(m1+m2) in code units
    ecc = get_eccentricity_vector(xyzmh_ptmass,vxyz_ptmass,1,2)
    print *,'eccentricity                 = ',sqrt(sum(ecc(1:3)**2))
    print *,'orbital period               = ',trim(get_nearest_unit(get_orbital_period(mu,dx,dv),'time'))
    print *,'orbital speed sink 1 (km/s)  = ',sqrt(dot_product(vxyz_ptmass(1:3,1),vxyz_ptmass(1:3,1)))*unit_velocity/km
    print *,'orbital speed sink 2 (km/s)  = ',sqrt(dot_product(vxyz_ptmass(1:3,2),vxyz_ptmass(1:3,2)))*unit_velocity/km
 endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  number of particles per shell (not rounded) needed for the shell
!  spacing to be wind_shell_spacing times the particle spacing on the shell
!+
!-----------------------------------------------------------------------
real function get_particles_per_shell(xyzmh_ptmassi,pmass) result(res)
 use physcon, only:pi
 use part,    only:iReff,imloss,ivwind
 real, intent(in) :: xyzmh_ptmassi(:),pmass

 ! sqrt(4 pi) factor from the neighbour distance on Fibonacci spheres
 res = (sqrt(4.*pi)*wind_shell_spacing*xyzmh_ptmassi(iReff)*xyzmh_ptmassi(imloss)/&
       (xyzmh_ptmassi(ivwind)*pmass))**(2./3.)

end function get_particles_per_shell

!-----------------------------------------------------------------------
!+
!  print the resolution that a given particle mass gives for each
!  wind-emitting sink, so the choice of particle mass can be checked
!+
!-----------------------------------------------------------------------
subroutine print_wind_resolution(xyzmh_ptmass,nptmass,pmass,gamma)
 use physcon, only:solarm
 use units,   only:umass,udist,unit_velocity
 use part,    only:imloss
 use wind,    only:wind_params,get_wind_solution
 real,    intent(in) :: xyzmh_ptmass(:,:),pmass,gamma
 integer, intent(in) :: nptmass
 type(wind_params) :: params
 real :: xyzmh(size(xyzmh_ptmass,1)),rsonic,res,dshell
 integer :: isink,nshell

 print "(/,a,es10.3,a)",' particle mass = ',pmass*umass/solarm,' Msun'
 do isink=1,nptmass
    if (xyzmh_ptmass(imloss,isink) <= 0.) cycle
    call get_wind_solution(xyzmh_ptmass(:,isink),wind_type,xyzmh,params,rsonic)
    res    = get_particles_per_shell(xyzmh,pmass)
    nshell = max(1,nint(res))
    print "(a,i2,a,i8,a,f6.3)",' sink ',isink,': particles per shell = ',nshell,&
          ', shell spacing / particle spacing = ',wind_shell_spacing*(nshell/res)**1.5
    call check_particles_per_shell(nshell)
    call check_wind_unbound(isink,params,gamma)
    if (wind_type == 1) then
       ! distance between shells at injection is v*dt = v*m*N/Mdot
       dshell = real(params%vwind/unit_velocity)*pmass*nshell/xyzmh(imloss)
       print "(a,f8.2)",'         shells between injection radius and sonic point = ',&
             (rsonic-params%Rinject)/udist/dshell
    endif
 enddo

end subroutine print_wind_resolution

!-----------------------------------------------------------------------
!+
!  warn if the particle mass is too large for the wind to be resolved
!  azimuthally, i.e. if there are too few particles per shell
!+
!-----------------------------------------------------------------------
subroutine check_particles_per_shell(nshell)
 use io, only:warning
 integer, intent(in) :: nshell

 if (nshell < min_particles_per_shell) call warning(label,'too few particles per shell to resolve the wind: '//&
    'decrease the particle mass',var='particles per shell',ival=nshell)

end subroutine check_particles_per_shell

!-----------------------------------------------------------------------
!+
!  warn if an adiabatic wind is bound, i.e. if the Bernoulli constant
!  v^2/2 + c_s^2/(gamma-1) - GM(1-alpha)/r is negative at injection.
!  Skipped when cooling, dust or Eddington-factor radiation pressure
!  can change the energy of the wind, since these can unbind it
!+
!-----------------------------------------------------------------------
subroutine check_wind_unbound(isink,params,gamma)
 use io,               only:warning
 use physcon,          only:gg,Rg,km
 use eos,              only:gmw
 use options,          only:icooling
 use dust_formation,   only:idust_opacity
 use ptmass_radiation, only:isink_radiation
 use wind,             only:wind_params
 integer,           intent(in) :: isink
 type(wind_params), intent(in) :: params
 real,              intent(in) :: gamma
 real :: alpha,vesc,bernoulli
 character(len=120) :: msg

 if (gamma <= 1.0001 .or. icooling > 0 .or. idust_opacity > 0 .or. isink_radiation > 1) return

 ! the 1D wind integration ignores the sink alpha unless isink_radiation = 1
 alpha = 0.
 if (isink_radiation == 1) alpha = params%alpha_rad
 vesc      = sqrt(2.*gg*params%Mstar*max(1.-alpha,0.)/params%Rinject)
 bernoulli = 0.5*params%vwind**2 + gamma*Rg*params%Twind/(gmw*(gamma-1.)) - 0.5*vesc**2
 if (bernoulli < 0.) then
    write(msg,"(a,i0,a,f0.1,a,f0.1,a)") 'wind from sink ',isink,' is bound (v_wind = ',params%vwind/km,&
         ' km/s, v_esc = ',vesc/km,' km/s): increase the wind velocity or alpha'
    call warning(label,trim(msg))
 endif

end subroutine check_wind_unbound

!-----------------------------------------------------------------------
!+
!  set the number of particles per shell for a wind-emitting sink.
!  The mass loss rate and particle mass are inputs and are not changed:
!  rounding to an integer number of particles per shell instead changes
!  the spacing between shells slightly
!+
!-----------------------------------------------------------------------
subroutine init_sink_resolution(isink,time_between_spheres,d_part)
 use io,          only:fatal
 use timestep,    only:tmax
 use part,        only:xyzmh_ptmass,massoftype,igas,imloss,ieject
 use injectutils, only:get_neighb_distance
 integer, intent(in)  :: isink
 real,    intent(out) :: time_between_spheres,d_part
 real :: res

 res = get_particles_per_shell(xyzmh_ptmass(:,isink),massoftype(igas))
 xyzmh_ptmass(ieject,isink) = max(1,nint(res))
 d_part = get_neighb_distance(nint(xyzmh_ptmass(ieject,isink)))
 print"(/,' number of particles per sphere for sink',i2,' = ',i7)",isink,nint(xyzmh_ptmass(ieject,isink))
 call check_particles_per_shell(nint(xyzmh_ptmass(ieject,isink)))

 !spheres properties
 time_between_spheres = massoftype(igas)*xyzmh_ptmass(ieject,isink)/xyzmh_ptmass(imloss,isink)
 if (time_between_spheres > tmax)  then
    print *,'time_between_spheres = ',time_between_spheres,' < tmax = ',tmax
    call fatal(label,'no shell ejection : tmax < time_between_spheres')
 endif

end subroutine init_sink_resolution

!-----------------------------------------------------------------------
!+
!  write wind parameters to log
!+
!-----------------------------------------------------------------------
subroutine logging(params,isink,time_between_spheres,neighbour_distance,&
     rsonic,tsonic,tboundary,tcross,tfill)

 use physcon,           only:pi,gg,au,km
 use units,             only:udist,unit_velocity,utime
 use timestep,          only:dtmax
 use wind,              only:wind_params
 use part,              only:massoftype,igas,xyzmh_ptmass,iReff,ivwind,&
                             ispinx,ispinz,ieject,iTwind,iwalpha

 integer,           intent(in) :: isink
 type(wind_params), intent(in) :: params
 real,              intent(in) :: time_between_spheres,neighbour_distance
 real,              intent(in), optional :: rsonic,tsonic,tboundary,tcross,tfill
 logical :: lsonic
 real :: vesc,wind_rotation_speed,rotation_speed_crit,dshell,nshells_sonic

 if (.not.present(rsonic)) then
    lsonic = .false.
 else
    lsonic = rsonic > params%rinject
 endif

 vesc = sqrt(2.*Gg*params%Mstar*(1.-xyzmh_ptmass(iwalpha,isink))/params%Rstar)
 ! distance between injected shells (code units); this differs slightly from
 ! the requested wind_shell_spacing because of the integer number of particles per shell
 dshell = time_between_spheres*real(params%vwind/unit_velocity)
 write (*,'(/,2(3x,A,es11.4))')&
      'mass_of_particles       : ',massoftype(igas),&
      'time_between_spheres    : ',time_between_spheres,&
      'wind velocity (km/s)    : ',params%vwind/km
 write (*,'(2(3x,A,es11.4),3x,A,i7)') &
      'shell spacing requested : ',wind_shell_spacing,&
      'shell spacing actual    : ',dshell/(neighbour_distance*xyzmh_ptmass(iReff,isink)),&
      'particles per sphere    : ',nint(xyzmh_ptmass(ieject,isink))
 write (*,'(2(3x,A,es11.4))') &
      'dist between spheres    : ',dshell
 write (*,'(3(3x,A,es11.4))') &
      'wind_temperature        : ',xyzmh_ptmass(iTwind,isink),&
      'injection_radius   (au) : ',params%rinject/au
 write (*,'(2(3x,A,es11.4))') &
      'stellar_radius     (au) : ',params%Rstar/au, &
      'rho_ini           (cgs) : ',params%Mdot/(4.*pi*params%rinject**2*params%vwind)
 if (present(tcross)) then
    write (*,'(3(3x,A,es11.4))') &
      'crossing time      (cu) : ',tcross,&
      'filling time       (cu) : ',tfill/utime,&
      'boundary time      (cu) : ',tboundary
 endif

 if (wind_type == 1 .and. present(rsonic)) then
    write (*,'(3(3x,A,es11.4),3x,A,i5)') &
         'distance to sonic point : ',(rsonic-params%rinject)/udist, &
         'sonic radius            : ',rsonic/udist, &
         'time_to_sonic_point     : ',tsonic/utime, &
         'number of shells to sonic : ',int((rsonic-params%rinject)/udist/dshell)
 endif
 wind_rotation_speed = sqrt(sum(xyzmh_ptmass(ispinx:ispinz,isink)**2))/xyzmh_ptmass(iReff,isink)**2
 rotation_speed_crit = sqrt(xyzmh_ptmass(4,isink)/xyzmh_ptmass(iReff,isink))
 if (wind_rotation_speed > 1e-20) then
    write (*,'(4(3x,A,es11.4))') &
         'rotation speed (km/s)      = ',wind_rotation_speed*unit_velocity/km,&
         'break-up velocity (km/s)   = ',rotation_speed_crit*unit_velocity/km,&
         'rotation_vel/critical_vel  = ',wind_rotation_speed/rotation_speed_crit,&
         'rotation_vel/vinject       = ',wind_rotation_speed/xyzmh_ptmass(ivwind,isink)
    if (wind_rotation_speed/rotation_speed_crit > 1.) then
       print*,'CAREFUL : rotation velocity exceeding equatorial break-up velocity'
    endif
 endif
 if (pulsating_wind) then
    print*,'number of ejected shells per pulsation period (should at least be > 10) ',pulsation_period/time_between_spheres
    print*,'pulsation period in code units = ',pulsation_period
    !sanity checks
    ! 1 - ensure that a minimum number of shells are ejected during a pulsation period
    if (pulsation_period/time_between_spheres < 10. ) print *,'WARNING! only ',pulsation_period/time_between_spheres,&
         ' shells will be ejected during a pulsation period'
    ! 2 - make sure the size of the boundary layer is not too big (< 0.2 injection_radius)
 elseif (lsonic) then
    !save a few models before the particles reach the sonic point
    if (dtmax > tsonic/utime) print *,'WARNING! dtmax > time to sonic point'
    !if solution subsonic, minimum resolution required so a few shells can be inserted between the injection radius and the sonic point
    nshells_sonic = (rsonic-params%rinject)/udist/dshell
    if (nshells_sonic < iboundary_spheres) print "(a,f6.2,a,i3)",' WARNING! resolution too low to pass sonic point: ',&
       nshells_sonic,' shells between injection radius and sonic point < iboundary_spheres = ',iboundary_spheres
 endif

end subroutine logging

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use physcon,           only:au,solarm,years
 use io,                only:fatal,iverbose,id,master
 use wind,              only:interp_wind_profile
 use part,              only:massoftype,igas,iReff,iboundary,nptmass,delete_particles_outside_sphere,&
                             delete_dead_particles_inside_radius,n_nucleation,ieject,imloss,ivwind,Bevol,Bxyz
 use partinject,        only:add_or_update_particle
 use units,             only:udist,umass,utime
 use dust_formation,    only:idust_opacity
 use ptmass_radiation,  only:isink_radiation

 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(in)    :: rho(:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: outer_sphere,inner_sphere,inner_boundary_sphere,ifirst,i,ipart,j, &
            nreleased,nboundaries,isink,itype,npart_per_sphere,nfill
 real    :: local_time,GM,r,v,u,rhoi,e,mass_lost,x0(3),v0(3),inner_radius,fdone,dum
 real    :: mass_of_spheres,time_between_spheres,rinject,wind_injection_speed
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

    npart_per_sphere = nint(xyzmh_ptmass(ieject,isink))
    mass_of_spheres =  massoftype(igas) * npart_per_sphere
    dtinject = mass_of_spheres / xyzmh_ptmass(imloss,isink)
    time_between_spheres = dtinject

    rinject = xyzmh_ptmass(iReff,isink)

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

    if (dtlast > 0.) then
       !
       ! delete particles that exit the outer boundary
       !
       inner_radius = rinject + deltaR_puls*sin(omega_puls*time)

       if (outer_boundary_au > rinject) call delete_particles_outside_sphere(x0,real(outer_boundary_au*au/udist),npart)
       call delete_dead_particles_inside_radius(x0,inner_radius,npart)
       if (npart_old /= npart .and. iverbose > 0) print *,'deleted ',npart_old-npart,'particles, remaining',npart
       npart_old = npart

       nreleased = nfill_domain
       nboundaries = iboundary_spheres
       !release background shells
       ipart = igas
       if (.not.released .and. (isink == 2 .or. onewind)) then
          if (nfill > 0) print '(/,"release background particles :",i7,3x,", Mdot : ",es12.4,", isink=",i1)',&
               nreleased*npart_per_sphere,xyzmh_ptmass(imloss,isink)/(solarm/umass)*(years/utime),isink
          do i = max(1,npart-nreleased*npart_per_sphere)+1,npart
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
          !call pulsating_wind_profile(time,local_time,r,v,u,rhoi,e,GM,i,inner_sphere)
       else
          if (idust_opacity == 2) then
             call interp_wind_profile(time,local_time,r,v,u,rhoi,e,GM,fdone,isink,JKmuS)
          else
             call interp_wind_profile(time,local_time,r,v,u,rhoi,e,GM,fdone,isink)
          endif
          if (iverbose > 0) print '(" ## update boundary  ",i4,2(i4),i7,i2,5x,8(1x,es12.5))',i,&
               inner_sphere,outer_sphere,npart,isink,time,local_time,r/xyzmh_ptmass(iReff,isink),v*udist/utime,&
               xyzmh_ptmass(imloss,isink) /(solarm/umass) * (years/utime)
       endif

       if (i > inner_sphere) then
          ! boundary sphere
          if (isink == 1) ifirst = (nboundaries-i+inner_sphere)*npart_per_sphere+1
          if (.not. onewind .and. isink ==2) ifirst = (nboundaries-i+inner_sphere)*npart_per_sphere+1 &
               + iboundary_spheres*int(xyzmh_ptmass(ieject,1))
          itype = ipart
          if (iverbose > 0) print '(" @@ update boundary ",i1,i4,2(i4),i7,i7,8(1x,es12.5))',isink,i,inner_sphere,&
               outer_sphere,ipart,ifirst,time,local_time,r/xyzmh_ptmass(iReff,isink),v*udist/utime,u,rhoi
       else
          ! ejected particles + create new  inner sphere
          ifirst = npart+1
          itype = igas
          if (iverbose > 0) print '(" ## eject particles [",i7,"-",i7,"], sink=",i1,4(i4),7(1x,es11.4))',&
               npart+1-npart_per_sphere,npart,isink,i,inner_sphere,outer_sphere,int(time/time_between_spheres),&
               time,local_time,r/xyzmh_ptmass(iReff,isink),v,u,rhoi,xyzmh_ptmass(imloss,isink)/(solarm/umass)*(years/utime)
       endif
       if (idust_opacity == 2) then
          call inject_sphere(i,ifirst,npart_per_sphere,r,v,u,npart,npartoftype,xyzh,vxyzu,rhoi,itype,x0,v0,isink,JKmuS)
       else
          call inject_sphere(i,ifirst,npart_per_sphere,r,v,u,npart,npartoftype,xyzh,vxyzu,rhoi,itype,x0,v0,isink)
       endif
       if (mhd) then
          do j = ifirst,ifirst+npart_per_sphere-1
             call set_injected_Bfield(xyzmh_ptmass(:,isink),xyzh(:,j),Bevol(:,j),Bxyz(:,j),massoftype(igas))
          enddo
       endif
    enddo

    if (nfill > 0 .and. outer_sphere == 1 .and. id==master) &
       print '(3x,"injecting background particles up to r = ",f8.2,&
            " (au), nparticles = ",i8,", isink=",i1)',r,&
            npart_per_sphere*(nfill_domain+iboundary_spheres),isink

    ! update sink particle properties
    mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)
    if (dtlast <= 0. .and. nfill > 0) mass_lost = mass_lost + mass_of_spheres * nfill
    xyzmh_ptmass(4,isink) = xyzmh_ptmass(4,isink) - mass_lost
    if (pulsating_wind) then
       inner_radius = rinject + deltaR_puls*sin(omega_puls*time)
       xyzmh_ptmass(5,isink) = inner_radius
    endif

    !
    ! return timestep constraint to ensure that time between sphere
    ! injections is adequately resolved
    !
    !dr = neighbour_distance*rinject
    !dtinject = 0.25*dr/sqrt(cs2max)
    dtinject = min(0.2*time_between_spheres,dtpulsation,dtinject)
 enddo
 if (time <= 0.) dtinject = 0.01*dtinject

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  set magnetic field
!+
!-----------------------------------------------------------------------
subroutine set_injected_Bfield(xyzmh_ptmassi,xyzhi,Bevoli,Bxyzi,pmassi)
 use part,  only:rhoh
 use units, only:unit_Bfield
 real, intent(in)  :: xyzmh_ptmassi(:),xyzhi(:),pmassi
 real, intent(out) :: Bevoli(:),Bxyzi(:)
 real :: r,r_hat(3),B_r_code,rhoi,dx(3)

 dx(1:3)=xyzhi(1:3)-xyzmh_ptmassi(1:3)
 r = sqrt(dot_product(dx,dx))
 r_hat = xyzhi(1:3)/r
 B_r_code = B_r/unit_Bfield
 ! bootstrap dens for newly injected wind particle (before density sum)
 rhoi = rhoh(xyzhi(4),pmassi)
 Bevoli(1:3) = B_r_code*r_hat / rhoi
 Bxyzi(1:3) = B_r_code*r_hat

end subroutine set_injected_Bfield

!-----------------------------------------------------------------------
!+
!  inject gas particles and/or reset position of boundary particles
!+
!-----------------------------------------------------------------------
subroutine inject_sphere(i,ifirst,ires,r,v,u,npart,npartoftype,xyzh,vxyzu,rho,itype,x0,v0,isink,JKmuS)

 use ptmass_radiation,  only:isink_radiation
 use part,              only:iTeff,dust_temp,xyzmh_ptmass,iReff,ispinx,ispinz,ivwind
 use injectutils,       only:inject_geodesic_sphere

 integer, intent(in)    :: i,ifirst,ires,isink,itype
 integer, intent(inout) :: npart
 real,    intent(in)    :: x0(3),v0(3),r,v,u,rho
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer, intent(inout) :: npartoftype(:)
 real,    intent(in), optional :: JKmuS(:)

 real :: rstar,mstar,omega_vec(3),vwind_terminal

 ! extract sink particle data to pass as parameters
 rstar = xyzmh_ptmass(iReff,isink)
 mstar = xyzmh_ptmass(4,isink)
 omega_vec = xyzmh_ptmass(ispinx:ispinz,isink)
 vwind_terminal = xyzmh_ptmass(ivwind,isink)

 if (present(JKmuS)) then
    call inject_geodesic_sphere(i,ifirst,ires,r,v,u, &
         npart,npartoftype,xyzh,vxyzu,rho,itype,x0,v0,isink,JKmuS, &
         rstar=rstar,mstar=mstar,omega_vec=omega_vec,vwind_terminal=vwind_terminal)
 else
    call inject_geodesic_sphere(i,ifirst,ires,r,v,u, &
         npart,npartoftype,xyzh,vxyzu,rho,itype,x0,v0,isink, &
         rstar=rstar,mstar=mstar,omega_vec=omega_vec,vwind_terminal=vwind_terminal)
 endif
 if (isink_radiation > 0) dust_temp(ifirst:ifirst+ires-1) = xyzmh_ptmass(iTeff,isink)

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
subroutine set_1D_wind_profile(params,isink,d_part,time_between_spheres,tboundary,tcross,tfill,onewind)
 use units,     only:utime
 use physcon,   only:au
 use io,        only:fileprefix,fatal,warning,id,master
 use timestep,  only:tmax
 use wind,      only:save_windprofile,wind_params

 integer,           intent(in)  :: isink
 type(wind_params), intent(in)  :: params
 real,              intent(in)  :: time_between_spheres,d_part
 logical,           intent(in)  :: onewind
 real,              intent(out) :: tboundary,tcross,tfill
 real :: tend,rout,rfill
 character(len=24) :: wfile

 tboundary = (iboundary_spheres+nfill_domain)*time_between_spheres
 tend      = max(tmax,tboundary)*utime
 write(wfile,'(a,i2.2,"_profile.dat")') trim(fileprefix),isink

 rout = outer_boundary_au*au
 rfill = rfill_domain_au*au
 call save_windprofile(params,rout,rfill,tend,tcross,tfill,wfile,isink)
 if (tboundary > tmax) call warning(label,'simulation time < time to reach the last boundary shell')

!define the number of background shells
 if (tfill < 1.d98 .and. rfill_domain_au > 1e-5 .and. (isink == 2 .or. onewind)) then
    nfill_domain = int(tfill/(utime*time_between_spheres))-iboundary_spheres
    if (id==master) print "(a,i0)",' number of background shells set to ',nfill_domain
    if (nfill_domain < 0 ) then
       call logging(params,isink,time_between_spheres,d_part,tboundary=tboundary,tcross=tcross,tfill=tfill)
       call fatal(label,'number of background shells < 0')
    endif
 endif

end subroutine set_1D_wind_profile

!-----------------------------------------------------------------------
!+
!  initialize oscillating inner boundary
!+
!-----------------------------------------------------------------------
subroutine init_pulsating_wind(pulsating_wind)

 use units,   only:unit_velocity,utime
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
!  Sets default options for the injection module
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)

 integer, intent(in), optional :: flag
 integer :: icase

 if (present(flag)) then
    icase = flag
 else
    icase = 0
 endif

! trans-sonic wind
 if (.not. isothermal .and. icase == 1) then
    wind_type = 1
 else
    wind_type = 0
 endif

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
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
 use physcon,      only:deg_to_rad
 use io,           only:warning
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 logical, save :: init_opt = .false.
 integer :: iwind_resolution_old
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
 ! obsolete option: the particle mass is now set in the .setup file (or read from the dump)
 iwind_resolution_old = 0
 call read_inopt(iwind_resolution_old,'iwind_resolution',db,default=0)
 if (iwind_resolution_old > 0) call warning(label,'iwind_resolution is obsolete and ignored: the particle mass '//&
    'is set by phantomsetup; see docs/examples/wind.rst to migrate old input files')
 call read_inopt(rfill_domain_au,'rfill_domain',db,errcount=nerr,min=0.)
 call read_inopt(wind_shell_spacing,'wind_shell_spacing',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(iboundary_spheres,'iboundary_spheres',db,errcount=nerr,min=0)
 call read_inopt(outer_boundary_au,'outer_boundary',db,errcount=nerr,min=0.)
 if (mhd) call read_inopt(B_r,'B_r',db,errcount=nerr)

end subroutine read_options_inject

end module inject
