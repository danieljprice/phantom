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
!   mpiutils, options, orbits, part, partinject, physcon, ptmass_radiation,
!   timestep, units, wind
!
 use dim, only:isothermal,mhd,maxptmass
 use wind, only:rfill_domain_au,wind_shell_spacing

 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,print_old_wind_options,&
      set_default_options_inject,update_injected_par,print_wind_resolution

 private
!
!--runtime settings for this module
!
! Read from input file
 integer :: wind_type = 0
 integer :: iboundary_spheres = 5
 integer :: nfill_domain(maxptmass) = 0  ! number of background shells of each sink
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

 logical :: pulsating_wind
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
 use io,                only:fatal,warning
 use orbits,            only:get_eccentricity_vector,get_orbital_period
 use wind,              only:wind_params,get_params_from_sink,solve_wind
 use physcon,           only:mass_proton_cgs,kboltz,km
 use eos,               only:gmw,polyk
 use units,             only:unit_velocity,udist,get_nearest_unit
 use part,              only:xyzmh_ptmass,vxyz_ptmass,imloss,iReff,ivwind,iTwind,nptmass,&
                              massoftype,igas,iboundary
 use injectutils,       only:init_jets

 integer, intent(out) :: ierr
 integer :: isink
 real :: time_between_spheres,d_part
 real :: tcross,rsonic,tsonic,tfill,tboundary,mu,wind_temperature
 real :: dx(3),dv(3),ecc(3)
 type(wind_params) :: params !reset wind temperature if init_muGamma

 ! polar jets from sink 1
 if (wind_type == 3) call init_jets(jet_edge_velocity,jet_opening_angle)

 if (icooling > 0) nwrite = nwrite+1
 ierr = 0
 nfill_domain = 0

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
       call set_1D_wind_profile(params,isink,d_part,time_between_spheres,tboundary,tcross,tfill)
    endif

! logging
    if (xyzmh_ptmass(imloss,isink) > 0.) &
       call logging(params,isink,time_between_spheres,d_part,rsonic,tsonic,tboundary,tcross,tfill)

 enddo
 if (count(nfill_domain(1:nptmass) > 0) > 1) call warning(label,&
    'background shells of several sinks may overlap: reduce rfill_domain if so')

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
 print "(a)"

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
                             ispinx,ispinz,ieject,iwalpha

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
      'wind_temperature        : ',params%Twind,&
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
!
!  Each wind-emitting sink emits shells of particles, numbered by the
!  time of their emission (shell s is emitted at time s*dt_shell). The
!  number of the next shell to emit is stored on the sink (iwshell), so
!  that no shell is emitted twice, including after a restart.
!
!  The innermost iboundary_spheres shells of each sink are boundary
!  particles, set from the 1D wind profile of the sink. Boundary
!  particles do not need to be stored in any particular order: their
!  shell (and hence their age) is found from their distance to the sink,
!  so this works with MPI and when particles are deleted or reordered.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,       only:iverbose,id,master
 use part,     only:massoftype,igas,iboundary,nptmass,imloss,ieject,iwshell
 use mpiutils, only:reduceall_mpi
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(in)    :: rho(:)
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: isink,ishell_first,ishell_last,nboundary
 real    :: dt_shell

 dtinject = huge(dtinject)
 !
 ! delete particles that leave the domain
 !
 if (dtlast > 0.) call delete_escaped_particles(time,xyzmh_ptmass,nptmass,npart,npart_old)
 nboundary = int(reduceall_mpi('+',npartoftype(iboundary)))

 do isink = 1,nptmass
    if (.not.is_wind_sink(xyzmh_ptmass(:,isink))) cycle

    dt_shell = massoftype(igas)*xyzmh_ptmass(ieject,isink)/xyzmh_ptmass(imloss,isink)
    ishell_last = floor(time/dt_shell)
    ishell_first = get_first_new_shell(isink,nint(xyzmh_ptmass(iwshell,isink)),ishell_last,&
                                       nboundary,time,dtlast)
    !
    ! update boundary shells, then emit the new shells
    !
    call update_boundary_shells(isink,time,dt_shell,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,xyzh,vxyzu)
    call emit_shells(isink,ishell_first,ishell_last,time,dt_shell,xyzmh_ptmass,vxyz_ptmass,&
                     npart,npartoftype,xyzh,vxyzu)
    if (iverbose > 0 .and. id==master .and. ishell_last >= ishell_first) &
       print "(a,i0,a,i0,a,i0)",' sink ',isink,': emitted shells ',ishell_first,' to ',ishell_last
    xyzmh_ptmass(iwshell,isink) = ishell_last + 1
    if (pulsating_wind) xyzmh_ptmass(5,isink) = get_inner_radius(xyzmh_ptmass(:,isink),time)
    !
    ! timestep constraint to ensure that the time between shell
    ! emissions is adequately resolved
    !
    dtinject = min(0.2*dt_shell,dtpulsation,dtinject)
 enddo
 if (time <= 0.) dtinject = 0.01*dtinject

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  whether a sink emits a wind
!+
!-----------------------------------------------------------------------
logical function is_wind_sink(xyzmh_ptmassi)
 use part, only:imloss,ivwind,ieject
 real, intent(in) :: xyzmh_ptmassi(:)

 is_wind_sink = xyzmh_ptmassi(imloss) > 0. .and. abs(xyzmh_ptmassi(ivwind)) > tiny(0.) &
                .and. xyzmh_ptmassi(ieject) >= 1.

end function is_wind_sink

!-----------------------------------------------------------------------
!+
!  number of the first shell to emit, given the number of the next
!  shell stored on the sink (0 if this sink has not emitted any shell)
!  and the number of the last shell due at the current time:
!
!  - on a new calculation (at t=0, with no boundary particles), the
!    boundary shells and the background shells (filling the domain)
!    are all emitted at once
!  - on a restart from a dump written before shells were counted,
!    the existing particles are kept and only new shells emitted
!+
!-----------------------------------------------------------------------
integer function get_first_new_shell(isink,ishell_next,ishell_last,nboundary,time,dtlast) result(ishell_first)
 use io, only:fatal,warning
 integer, intent(in) :: isink,ishell_next,ishell_last,nboundary
 real,    intent(in) :: time,dtlast
 integer :: nnew

 if (ishell_next == 0) then
    if (nboundary > 0 .or. time > 0.) then
       ishell_first = ishell_last + 1
       call warning(label,'number of wind shells not in dump: keeping the existing particles',ival=isink)
    else
       ishell_first = ishell_last - (iboundary_spheres + nfill_domain(isink)) + 1
    endif
 else
    ishell_first = ishell_next
    nnew = ishell_last - ishell_first + 1
    if (nnew < 0 .or. (dtlast <= 0. .and. nnew > 1)) then
       ! starting at a different time, or the shell spacing has changed (e.g. a new
       ! mass loss rate or particle mass): continue from the shell due now
       call warning(label,'wind shell number not consistent with time: restarting the count',ival=isink)
       ishell_first = ishell_last + 1
    elseif (nnew > max(iboundary_spheres,1)) then
       call fatal(label,'emission of more than one shell per step, timestep likely too large!',ival=nnew)
    endif
 endif

end function get_first_new_shell

!-----------------------------------------------------------------------
!+
!  emit shells ishell_first to ishell_last from sink isink at the
!  current time. Shells older than the boundary shells (background
!  shells, on a new calculation) are emitted directly as gas.
!  Each shell is created by one MPI task, but the mass of every shell
!  is taken from the sink on all tasks
!+
!-----------------------------------------------------------------------
subroutine emit_shells(isink,ishell_first,ishell_last,time,dt_shell,xyzmh_ptmass,vxyz_ptmass,&
                       npart,npartoftype,xyzh,vxyzu)
 use io,          only:id,nprocs
 use part,        only:massoftype,igas,iboundary,ieject
 use injectutils, only:get_shell_rotation,get_shell_direction,jets
 integer, intent(in)    :: isink,ishell_first,ishell_last
 real,    intent(in)    :: time,dt_shell
 real,    intent(inout) :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ishell,itype,j,nps
 real    :: age,rotmat(3,3),dir(3)
 logical :: jets_on

 nps = nint(xyzmh_ptmass(ieject,isink))
 jets_on = (jets .and. isink == 1)
 do ishell=ishell_first,ishell_last
    age = time - ishell*dt_shell
    itype = iboundary
    if (age >= iboundary_spheres*dt_shell) itype = igas
    if (modulo(ishell+isink,nprocs) == id) then
       call get_shell_rotation(ishell,jets_on,rotmat)
       do j=0,nps-1
          call get_shell_direction(j,nps,rotmat,jets_on,dir)
          call set_wind_particle(npart+1,itype,isink,dir,age,xyzmh_ptmass(:,isink),vxyz_ptmass(:,isink),&
                                 npart,npartoftype,xyzh,vxyzu)
       enddo
    endif
    xyzmh_ptmass(4,isink) = xyzmh_ptmass(4,isink) - nps*massoftype(igas)
 enddo

end subroutine emit_shells

!-----------------------------------------------------------------------
!+
!  update the boundary particles of sink isink from the 1D wind profile,
!  releasing them as gas once they are older than the boundary shells.
!  The age of a boundary particle is found from its distance to the sink
!  (boundary particles move with constant velocity), rounded to the age
!  of the nearest shell so that all particles of a shell stay together
!+
!-----------------------------------------------------------------------
subroutine update_boundary_shells(isink,time,dt_shell,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype,xyzh,vxyzu)
 use part,       only:iphase,iamtype,iboundary,igas,iunknown,iseed_sink,isdead_or_accreted
 use partinject, only:updated_particle
 use wind,       only:wind_age_at_r
 integer, intent(in)    :: isink
 real,    intent(in)    :: time,dt_shell
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,ishell
 real    :: dx(3),dir(3),r,age

 do i=1,npart
    if (iseed_sink(i) /= isink) cycle
    if (iamtype(iphase(i)) /= iboundary .or. isdead_or_accreted(xyzh(4,i))) cycle
    dx = xyzh(1:3,i) - xyzmh_ptmass(1:3,isink)
    r  = sqrt(dot_product(dx,dx))
    dir = dx/r
    ishell = nint((time - wind_age_at_r(isink,r))/dt_shell)
    age = time - ishell*dt_shell
    if (age < iboundary_spheres*dt_shell) then
       call set_wind_particle(i,iboundary,isink,dir,age,xyzmh_ptmass(:,isink),vxyz_ptmass(:,isink),&
                              npart,npartoftype,xyzh,vxyzu)
    else
       ! release as gas: timestep set in update_injected_particles
       call set_wind_particle(i,igas,isink,dir,age,xyzmh_ptmass(:,isink),vxyz_ptmass(:,isink),&
                              npart,npartoftype,xyzh,vxyzu)
       iphase(i) = int(iunknown,kind=1)
       updated_particle = .true.
    endif
 enddo

end subroutine update_boundary_shells

!-----------------------------------------------------------------------
!+
!  set particle i (or create it if i = npart+1) as a wind particle of
!  sink isink with a given age, in direction dir (unit vector) from the
!  sink, using the 1D wind profile of the sink
!+
!-----------------------------------------------------------------------
subroutine set_wind_particle(i,itype,isink,dir,age,xyzmh_ptmassi,vxyz_ptmassi,npart,npartoftype,xyzh,vxyzu)
 use part,             only:massoftype,igas,hn,n_nucleation,iReff,ispinx,ispinz,ivwind,iTeff,&
                            dust_temp,Bevol,Bxyz,rho
 use partinject,       only:add_or_update_particle
 use injectutils,      only:get_wind_velocity,jets
 use wind,             only:interp_wind_profile
 use dust_formation,   only:idust_opacity
 use ptmass_radiation, only:isink_radiation
 integer, intent(in)    :: i,itype,isink
 real,    intent(in)    :: dir(3),age,xyzmh_ptmassi(:),vxyz_ptmassi(:)
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: r,v,u,rhoi,h,xyz(3),vxyz(3),JKmuS(n_nucleation)

 call interp_wind_profile(isink,age,r,v,u,rhoi,JKmuS)
 call get_wind_velocity(dir,r,v,xyzmh_ptmassi(iReff),xyzmh_ptmassi(4),xyzmh_ptmassi(ispinx:ispinz),&
                        xyzmh_ptmassi(ivwind),(jets .and. isink == 1),xyz,vxyz)
 ! h from the density of the wind, so that boundary particles follow the expansion
 h = hn(rhoi/massoftype(igas))
 xyz  = xyz  + xyzmh_ptmassi(1:3)
 vxyz = vxyz + vxyz_ptmassi(1:3)
 if (idust_opacity == 2) then
    call add_or_update_particle(itype,xyz,vxyz,h,u,i,npart,npartoftype,xyzh,vxyzu,JKmuS,isink)
 else
    call add_or_update_particle(itype,xyz,vxyz,h,u,i,npart,npartoftype,xyzh,vxyzu,isink=isink)
 endif
 ! add_or_update_particle does not reset h or rho of boundary particles
 ! (with APR the stored rho, not h, is the starting guess for the density)
 xyzh(4,i) = h
 rho(i)    = rhoi
 if (isink_radiation > 0) dust_temp(i) = xyzmh_ptmassi(iTeff)
 if (mhd) call set_injected_Bfield(xyzmh_ptmassi,xyzh(:,i),Bevol(:,i),Bxyz(:,i),massoftype(igas))

end subroutine set_wind_particle

!-----------------------------------------------------------------------
!+
!  delete particles outside the outer boundary (centred on the centre
!  of mass of the sinks) and dead particles
!+
!-----------------------------------------------------------------------
subroutine delete_escaped_particles(time,xyzmh_ptmass,nptmass,npart,npart_old)
 use physcon, only:au
 use units,   only:udist
 use io,      only:iverbose
 use part,    only:delete_particles_outside_sphere,delete_dead_particles_inside_radius
 real,    intent(in)    :: time,xyzmh_ptmass(:,:)
 integer, intent(in)    :: nptmass
 integer, intent(inout) :: npart,npart_old
 real    :: xcom(3),mtot
 integer :: isink

 if (outer_boundary_au > 0. .and. nptmass > 0) then
    mtot = sum(xyzmh_ptmass(4,1:nptmass))
    do isink=1,3
       xcom(isink) = dot_product(xyzmh_ptmass(isink,1:nptmass),xyzmh_ptmass(4,1:nptmass))/mtot
    enddo
    call delete_particles_outside_sphere(xcom,real(outer_boundary_au*au/udist),npart)
 endif
 do isink=1,nptmass
    if (.not.is_wind_sink(xyzmh_ptmass(:,isink))) cycle
    call delete_dead_particles_inside_radius(xyzmh_ptmass(1:3,isink),&
         get_inner_radius(xyzmh_ptmass(:,isink),time),npart)
 enddo
 if (npart_old /= npart .and. iverbose > 0) print *,'deleted ',npart_old-npart,'particles, remaining',npart
 npart_old = npart

end subroutine delete_escaped_particles

!-----------------------------------------------------------------------
!+
!  inner radius of the wind of a sink (oscillating if the wind pulsates)
!+
!-----------------------------------------------------------------------
real function get_inner_radius(xyzmh_ptmassi,time) result(rin)
 use part, only:iReff
 real, intent(in) :: xyzmh_ptmassi(:),time

 rin = xyzmh_ptmassi(iReff) + deltaR_puls*sin(omega_puls*time)

end function get_inner_radius

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

 ! radial field from the sink
 dx(1:3) = xyzhi(1:3)-xyzmh_ptmassi(1:3)
 r = sqrt(dot_product(dx,dx))
 r_hat = dx/r
 B_r_code = B_r/unit_Bfield
 ! bootstrap dens for newly injected wind particle (before density sum)
 rhoi = rhoh(xyzhi(4),pmassi)
 Bevoli(1:3) = B_r_code*r_hat / rhoi
 Bxyzi(1:3) = B_r_code*r_hat

end subroutine set_injected_Bfield

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
subroutine set_1D_wind_profile(params,isink,d_part,time_between_spheres,tboundary,tcross,tfill)
 use units,     only:utime
 use physcon,   only:au
 use io,        only:fileprefix,fatal,warning,id,master
 use timestep,  only:tmax
 use wind,      only:save_windprofile,wind_params

 integer,           intent(in)  :: isink
 type(wind_params), intent(in)  :: params
 real,              intent(in)  :: time_between_spheres,d_part
 real,              intent(out) :: tboundary,tcross,tfill
 real :: tend,rout,rfill
 character(len=len(fileprefix)+16) :: wfile ! prefix + NN(NN)_profile.dat

 tboundary = (iboundary_spheres+nfill_domain(isink))*time_between_spheres
 tend      = max(tmax,tboundary)*utime
 write(wfile,'(a,i0.2,"_profile.dat")') trim(fileprefix),isink

 rout = outer_boundary_au*au
 rfill = rfill_domain_au*au
 call save_windprofile(params,rout,rfill,tend,tcross,tfill,wfile,isink)
 if (tboundary > tmax) call warning(label,'simulation time < time to reach the last boundary shell')

!define the number of background shells
 if (tfill < 1.d98 .and. rfill_domain_au > 1e-5) then
    nfill_domain(isink) = int(tfill/(utime*time_between_spheres))-iboundary_spheres
    if (id==master) print "(a,i0)",' number of background shells set to ',nfill_domain(isink)
    if (nfill_domain(isink) < 0 ) then
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
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 logical, save :: init_opt = .false.

 if (.not.init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif
 call read_inopt(wind_type,'wind_type',db,default=0,errcount=nerr,min=0,max=3)
 if (wind_type == 2) then
    call read_inopt(pulsation_period_days,'pulsation_period',db,errcount=nerr,min=0.)
    call read_inopt(piston_velocity_km_s,'piston_velocity',db,errcount=nerr,min=0.)
 endif
 if (wind_type == 3) then
    call read_inopt(jet_edge_velocity,'jet_edge_velocity',db,errcount=nerr,min=0.)
    call read_inopt(jet_opening_angle_degree,'jet_opening_angle',db,errcount=nerr,min=0.,max=90.)
    jet_opening_angle = jet_opening_angle_degree*deg_to_rad
 endif
 ! obsolete options: the wind and the particle mass are now set in the .setup file (or read from the dump)
 call print_old_wind_options(db)
 call read_inopt(rfill_domain_au,'rfill_domain',db,default=0.,errcount=nerr,min=0.)
 call read_inopt(wind_shell_spacing,'wind_shell_spacing',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(iboundary_spheres,'iboundary_spheres',db,errcount=nerr,min=0)
 call read_inopt(outer_boundary_au,'outer_boundary',db,errcount=nerr,min=0.)
 if (mhd) call read_inopt(B_r,'B_r',db,errcount=nerr)

end subroutine read_options_inject

!-----------------------------------------------------------------------
!+
!  .in files from before the wind became a property of each sink give
!  the wind of the primary star and the resolution: these are ignored,
!  so print the options to set in the .setup file instead, or return
!  them if the optional arguments are given (as phantomsetup does).
!  iwind_resolution was the icosahedron resolution q before January
!  2025 (40q(q-1)+12 particles per shell), then the number of particles
!+
!-----------------------------------------------------------------------
subroutine print_old_wind_options(db,mdot_msun_yr,vwind_km_s,wind_temp,rinj_au,alpha,pmass)
 use infile_utils, only:inopts,read_inopt
 use physcon,      only:au,km,years
 type(inopts), intent(inout) :: db(:)
 real, intent(out), optional :: mdot_msun_yr,vwind_km_s,wind_temp,rinj_au,alpha,pmass
 integer :: ierr,isonic,isink_rad,iwindres,nshell
 real    :: mdot,v,temp,r,a,spacing,m
 logical :: old_wind

 mdot = 0.; v = 0.; temp = 0.; r = 0.; a = 0.; spacing = 1.; m = 0.
 isonic = 0; isink_rad = 0; iwindres = 0; nshell = 0
 call read_inopt(iwindres,'iwind_resolution',db,err=ierr)
 call read_inopt(mdot,'wind_mass_rate',db,err=ierr)
 old_wind = (ierr == 0)
 call read_inopt(isonic,'sonic_type',db,err=ierr)
 if (isonic /= 1) call read_inopt(v,'wind_velocity',db,err=ierr)
 call read_inopt(temp,'wind_temperature',db,err=ierr)
 call read_inopt(r,'wind_inject_radius',db,err=ierr)
 call read_inopt(isink_rad,'isink_radiation',db,err=ierr)
 if (isink_rad == 1 .or. isink_rad == 3) call read_inopt(a,'alpha_rad',db,err=ierr)
 call read_inopt(spacing,'wind_shell_spacing',db,err=ierr)
 if (iwindres > 0) then
    nshell = iwindres
    if (iwindres < 12) nshell = 40*iwindres*(iwindres-1) + 12
    ! particle mass of the icosahedral shells: spacing*d*Rinject*Mdot/(N*v), with neighbour
    ! distance d = 2/((2q-1)*sqrt(sqrt(5)*phi)) on a sphere of unit radius
    if (iwindres < 12 .and. mdot*v*r > 0.) &
       m = spacing*2./((2.*iwindres-1.)*sqrt(sqrt(5.)*0.5*(1.+sqrt(5.))))*r*au/(v*km)/years*mdot/nshell
 endif

 if (present(mdot_msun_yr)) then
    mdot_msun_yr = mdot
    if (present(vwind_km_s)) vwind_km_s = v
    if (present(wind_temp))  wind_temp  = temp
    if (present(rinj_au))    rinj_au    = r
    if (present(alpha))      alpha      = a
    if (present(pmass))      pmass      = m
    return
 endif
 if (.not.old_wind .and. iwindres <= 0) return

 print "(/,a)",' WARNING! old wind options in the .in file are ignored: set instead in the .setup file'
 if (old_wind) then
    print "(4x,a,es10.3)",'primary_mdot = ',mdot
    print "(4x,a,f0.2)",'primary_vwind = ',v
    print "(4x,a,f0.1)",'primary_wind_temp = ',temp
    print "(4x,a,f0.3,a)",'primary_Rinj = ',r,' au'
    print "(4x,a,f0.3)",'primary_alpha = ',a
 endif
 if (m > 0.) then
    print "(4x,a,es11.4,a,i0,a)",'iresolution = 0 and mass_of_particles = ',m,&
          ' (the particle mass given by iwind_resolution = ',iwindres,')'
    print "(4x,a,i0,a)",'or iresolution = 1 and particles_per_shell = ',nshell,' (the same number of particles per shell)'
 elseif (nshell > 0) then
    print "(4x,a,i0,a)",'iresolution = 1 and particles_per_shell = ',nshell,' (from iwind_resolution)'
 endif
 print "(a)",''

end subroutine print_old_wind_options

end module inject
