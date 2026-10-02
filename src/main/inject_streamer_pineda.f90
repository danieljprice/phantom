!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection module for "streamer" simulations
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - Mdot         : *mass injection rate, in Msun/yr (peak rate if imdot_func > 0)*
!   - dust_frac    : *Dust fraction in smallest dust bin*
!   - mdot_func    : *functional form of dM/dt(t) (0=const)*
!   - nstreams     : *number of direct Cartesian streamers*
!   - omega        : *angular velocity of cloud stream originates from (s^-1)*
!   - phi0         : *phi0 parameter from the Mendoza+09 streamer, in degrees*
!   - r0           : *r0 parameter from the Mendoza+09 streamer, in au*
!   - r_inj        : *distance from CoM where stream is injected, in au*
!   - stream_width : *radius of injected stream in au*
!   - tend         : *end time of injection in years; negative means infinite*
!   - theta0       : *theta0 parameter from the Mendoza+09 streamer, in degrees*
!   - tstart       : *start time of injection in years*
!   - vr_0         : *radial velocity of cloud stream origin, in km/s*
!
! :Dependencies: dim, infile_utils, options, part, partinject, physcon,
!   random, units, vectorutils
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'streamer'

 public :: inject_particles, write_options_inject, read_options_inject
 public :: init_inject, set_default_options_inject, update_injected_par

 type, private :: streamer_t
    real :: mdot = 1.e-7
    real :: mdot_code = 0.
    real :: width = 10.
    real :: tstart = 0.
    real :: tend = -1.
    real :: tstart_code = 0.
    real :: tend_code = -1.
    real :: position(3) = (/-2414.,-2225.,880./)
    real :: velocity(3) = (/0.35,0.,-0.21/)
    integer :: ninjected = 0
 end type streamer_t

 real, private :: Mdot = 1e-7
 real, private :: Mdotcode = 0.
 integer, private :: nstreams = 1
 type(streamer_t), allocatable, private :: streams(:)
 integer, private :: imdot_func = 0
 integer, private :: stream_model = 0
 integer, private :: sym_stream = 0
 real, private :: stream_width = 10.
 real, private :: r_inj = 100.
 real, private :: phi0 = 0.
 real, private :: theta0 = 90.
 real, private :: r0 = 1500.
 real, private :: omega = 1e-11
 real, private :: vr_0 = 1.5
 real, private :: tstart = 0.
 real, private :: tend = -1.0
 real, private :: tstart_code = 0.
 real, private :: tend_code = -1.0
 real, private :: x_stream = -2414.
 real, private :: y_stream = -2225.
 real, private :: z_stream = 880.
 real, private :: vx_stream = 0.35
 real, private :: vy_stream = 0.0
 real, private :: vz_stream = -0.21
 real, private :: dust_frac = 0.01

contains

subroutine init_inject(ierr)
 use units,   only:umass,utime
 use physcon, only:years,solarm
 integer, intent(out) :: ierr
 integer :: istream,istat
 real :: speed

 ierr = 0
!
!--convert mass injection rate to code units
!
 Mdotcode = Mdot*(solarm/umass)/(years/utime)
 tstart_code = tstart*(years/utime)
 if (tend < 0.) then
    tend_code = tend
 else
    tend_code = tend*(years/utime)
 endif
 if (stream_model == 0) then
    print*,' Mdot is ',Mdot,' Msun/yr, which is ',Mdotcode,' in code units',umass,utime
    print*,' injection time window: ',tstart,' to ',tend,' yr = ',tstart_code,' to ',tend_code,' code units'
 else
    call ensure_stream_array(nstreams,istat)
    if (istat /= 0) then
       print*,' ERROR: could not allocate direct streamers'
       ierr = 1
       return
    endif
    if (imdot_func /= 0 .and. nstreams > 1) then
       print*,' ERROR: multiple direct streamers require mdot_func=0'
       ierr = 1
       return
    endif
    do istream=1,nstreams
       streams(istream)%mdot_code = streams(istream)%mdot*(solarm/umass)/(years/utime)
       streams(istream)%tstart_code = streams(istream)%tstart*(years/utime)
       if (streams(istream)%tend < 0.) then
          streams(istream)%tend_code = streams(istream)%tend
       else
          streams(istream)%tend_code = streams(istream)%tend*(years/utime)
       endif
       speed = sqrt(sum(streams(istream)%velocity**2))
       if (streams(istream)%width <= 0. .or. speed <= tiny(speed) .or. &
           (streams(istream)%tend >= 0. .and. streams(istream)%tend < streams(istream)%tstart)) then
          print*,' ERROR: invalid width, velocity, or time window for streamer ',istream
          ierr = 1
          return
       endif
       print*,' direct streamer ',istream,' Mdot (Msun/yr): ',streams(istream)%mdot
       print*,' direct streamer ',istream,' time window (yr): ', &
               streams(istream)%tstart,streams(istream)%tend
       print*,' direct streamer ',istream,' state (au, km/s): ', &
               streams(istream)%position,streams(istream)%velocity
    enddo
 endif

end subroutine init_inject
!-----------------------------------------------------------------------
!+
!  Main routine handling injection at the L1 point.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,rho,xyzmh_ptmass,vxyz_ptmass,&
           npart,npart_old,npartoftype,dtinject)
 use dim,       only:use_dust,maxdusttypes
 use options,   only:use_dustfrac
 use part,      only:igas,hfact,massoftype,nptmass,gravity,dustfrac,dustevol
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,solarm,years
 use units,     only:udist,umass,utime
 use random,    only:ran2
 use vectorutils, only:make_perp_frame
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(in)    :: rho(:)
 real,    intent(inout) :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real :: mtot,phi0_rad,theta0_rad
 real :: xc,yc,zc,vxc,vyc,vzc,Mdot_now,mass_to_inject
 real :: omega_cu,vr_0_cu,end_time
 integer :: ninject_target,istream

 dtinject = huge(dtinject) ! no timestep constraint from injection

 if (stream_model == 1) then
    if (.not.allocated(streams)) return
    do istream=1,nstreams
       if (imdot_func == 0) then
          ninject_target = cumulative_particle_count(streams(istream),time,massoftype(igas)) - &
                           streams(istream)%ninjected
       else
          if (streams(istream)%tend_code < 0.) then
             end_time = huge(time)
          else
             end_time = streams(istream)%tend_code
          endif
          if (time < streams(istream)%tstart_code .or. time > end_time) cycle
          Mdot_now = Mdotfunc((time-streams(istream)%tstart_code)*utime/years)
          mass_to_inject = Mdot_now*(solarm/umass)/(years/utime)*dtlast
          ninject_target = ceiling(mass_to_inject/massoftype(igas))
       endif
       if (ninject_target <= 0) cycle
       xc = streams(istream)%position(1)
       yc = streams(istream)%position(2)
       zc = streams(istream)%position(3)
       vxc = streams(istream)%velocity(1)*1.0e5*utime/udist
       vyc = streams(istream)%velocity(2)*1.0e5*utime/udist
       vzc = streams(istream)%velocity(3)*1.0e5*utime/udist
       call inject_one_stream(xc,yc,zc,vxc,vyc,vzc,streams(istream)%width,ninject_target)
       streams(istream)%ninjected = streams(istream)%ninjected+ninject_target
    enddo
    return
 endif

 if (tend_code < 0.) then
    end_time = huge(time)
 else
    end_time = tend_code
 endif
 if (time < tstart_code .or. time > end_time) return

 mtot = 0.
 phi0_rad = phi0*pi/180.
 theta0_rad = theta0*pi/180.
 vr_0_cu = vr_0*1.0e5*utime/udist
 omega_cu = omega*utime
 if (gravity) then
    mtot = sum(xyzmh_ptmass(4,1:nptmass)) + npartoftype(igas)*massoftype(igas)
 else
    mtot = sum(xyzmh_ptmass(4,1:nptmass))
 endif
 if (imdot_func > 0) then
    Mdot_now = Mdotfunc((time-tstart_code)*utime/years)
    Mdotcode = Mdot_now*(solarm/umass)/(years/utime)
 endif
 ninject_target = ceiling(Mdotcode*dtlast/massoftype(igas))
 if (ninject_target <= 0) return
 call mendoza_state(mtot,r0,omega_cu,theta0_rad,phi0_rad,vr_0_cu, &
                    r_inj,xc,yc,zc,vxc,vyc,vzc)
 call inject_one_stream(xc,yc,zc,vxc,vyc,vzc,stream_width,ninject_target)

contains
subroutine inject_one_stream(xc,yc,zc,vxc,vyc,vzc,width,ninject_target)
 real,    intent(in) :: xc,yc,zc,vxc,vyc,vzc,width
 integer, intent(in) :: ninject_target
 real :: xyzi(3),vxyz(3),ex(3),ey(3),ez(3)
 real :: h,vt,rand_radius,rand_angle,rrand,theta,dx_loc,dz_loc
 real :: x_si,y_si,z_si
 integer :: ninjected,ipart,iseed

 h = sqrt(hfact*width*width/real(ninject_target))
 vt = sqrt(vxc*vxc+vyc*vyc+vzc*vzc)
 ex = (/vxc,vyc,vzc/)/vt
 call make_perp_frame(ex,ey,ez)
 ninjected = 0
 ipart = npart+1
 iseed = npartoftype(igas)
 do while (ninjected < ninject_target)
    rand_radius = ran2(iseed)
    rand_angle = ran2(iseed)
    iseed = iseed-1
    rrand = width*sqrt(rand_radius)
    theta = 2.*pi*rand_angle
    dx_loc = rrand*cos(theta)
    dz_loc = rrand*sin(theta)
    x_si = xc+dx_loc*ey(1)+dz_loc*ez(1)
    y_si = yc+dx_loc*ey(2)+dz_loc*ez(2)
    z_si = zc+dx_loc*ey(3)+dz_loc*ez(3)
    xyzi = (/x_si,y_si,z_si/)
    vxyz = (/vxc,vyc,vzc/)
    ninjected = ninjected+1
    call add_or_update_particle(igas,xyzi,vxyz,h,rand_radius,ipart, &
                                npart,npartoftype,xyzh,vxyzu)
    call set_injected_dust_properties(ipart,dust_frac)
    ipart = ipart+1
    select case(sym_stream)
    case(1)
       xyzi = (/-x_si,-y_si,z_si/)
       vxyz = (/-vxc,-vyc,vzc/)
    case(2)
       xyzi = (/-x_si,-y_si,-z_si/)
       vxyz = (/-vxc,-vyc,-vzc/)
    case(3)
       xyzi = (/x_si,y_si,-z_si/)
       vxyz = (/vxc,vyc,-vzc/)
    end select
    if (sym_stream > 0) then
       call add_or_update_particle(igas,xyzi,vxyz,h,rand_radius,ipart, &
                                   npart,npartoftype,xyzh,vxyzu)
       call set_injected_dust_properties(ipart,dust_frac)
       ipart = ipart+1
    endif
 enddo
end subroutine inject_one_stream

integer function cumulative_particle_count(stream,time_now,particle_mass)
 type(streamer_t), intent(in) :: stream
 real, intent(in) :: time_now,particle_mass
 real :: active_time,time_stop

 if (time_now <= stream%tstart_code) then
    active_time = 0.
 else
    if (stream%tend_code < 0.) then
       time_stop = time_now
    else
       time_stop = min(time_now,stream%tend_code)
    endif
    active_time = max(0.,time_stop-stream%tstart_code)
 endif
 cumulative_particle_count = floor(stream%mdot_code*active_time/particle_mass)
end function cumulative_particle_count

subroutine set_injected_dust_properties(ipart_in,dust_frac_val)
 integer, intent(in) :: ipart_in
 real,    intent(in) :: dust_frac_val
 real :: dustevol_val

 if (use_dust .and. use_dustfrac) then
    dustfrac(:, ipart_in) = 0.
    dustevol(:, ipart_in) = 0.
    if (maxdusttypes > 0) then
       dustfrac(1, ipart_in) = dust_frac_val
       if (dust_frac_val > 0. .and. dust_frac_val < 1.) then
          dustevol_val = sqrt(dust_frac_val/(1. - dust_frac_val))
          dustevol(1, ipart_in) = dustevol_val
       endif
    endif
 endif

end subroutine set_injected_dust_properties

!-----------------------------------------------------------------------
!+
!  Function to return the total mass injected up to time t
!  by computing the integral \int Mdot dt
!+
!-----------------------------------------------------------------------
real function Mdotfunc(t)
 real, intent(in) :: t

 select case(imdot_func)
 case(1)
    Mdotfunc = Mdotcode*(t/tend)**(-5./3.)*(1.-(t/tend)**(-4./3.))
 case default
    Mdotfunc = Mdotcode
 end select

end function Mdotfunc

end subroutine inject_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file.
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt,int_to_string
 use options, only:use_dustfrac
 use dim, only:use_dust

 integer, intent(in) :: iunit
 integer :: istream,istat
 character(len=12) :: label

 call write_inopt(imdot_func,'mdot_func','functional form of dM/dt(t) (0=const)',iunit)
 call write_inopt(stream_model,'stream_model', &
                  'streamer kinematics (0=PIMS, 1=direct Cartesian states)',iunit)
 call write_inopt(sym_stream,'sym_stream', &
                  'balance angular momentum (0=no, 1=Lx,Ly, 2=Lx,Ly,Lz, 3=Lz)',iunit)

 if (stream_model == 0) then
    call write_inopt(Mdot,'Mdot','mass injection rate, in Msun/yr (peak rate if imdot_func > 0)',iunit)
    call write_inopt(stream_width,'stream_width','radius of injected stream in au',iunit)
    call write_inopt(tstart,'tstart','start time of injection in years',iunit)
    call write_inopt(tend,'tend','end time of injection in years; negative means infinite',iunit)
    call write_inopt(omega,'omega','angular velocity of cloud stream originates from (s^-1)',iunit)
    call write_inopt(r0,'r0','r0 parameter from the Mendoza+09 streamer, in au',iunit)
    call write_inopt(phi0,'phi0','phi0 parameter from the Mendoza+09 streamer, in degrees',iunit)
    call write_inopt(theta0,'theta0','theta0 parameter from the Mendoza+09 streamer, in degrees',iunit)
    call write_inopt(r_inj,'r_inj','distance from CoM where stream is injected, in au',iunit)
    call write_inopt(vr_0,'vr_0','radial velocity of cloud stream origin, in km/s',iunit)
 else
    call ensure_stream_array(nstreams,istat)
    if (istat /= 0) return
    call write_inopt(nstreams,'nstreams','number of direct Cartesian streamers',iunit)
    do istream=1,nstreams
       label = trim(int_to_string(istream))
       call write_inopt(streams(istream)%mdot,'Mdot_stream_'//trim(label), &
                        'mass injection rate in Msun/yr',iunit)
       call write_inopt(streams(istream)%width,'stream_width_'//trim(label), &
                        'stream radius in au',iunit)
       call write_inopt(streams(istream)%tstart,'tstart_stream_'//trim(label), &
                        'injection start time in years',iunit)
       call write_inopt(streams(istream)%tend,'tend_stream_'//trim(label), &
                        'injection end time in years; negative is infinite',iunit)
       call write_inopt(streams(istream)%position(1),'x_stream_'//trim(label), &
                        'x/R.A. injection offset in au',iunit)
       call write_inopt(streams(istream)%position(2),'y_stream_'//trim(label), &
                        'y/Dec. injection offset in au',iunit)
       call write_inopt(streams(istream)%position(3),'z_stream_'//trim(label), &
                        'z/LOS injection offset in au',iunit)
       call write_inopt(streams(istream)%velocity(1),'vx_stream_'//trim(label), &
                        'x/R.A. injection velocity in km/s',iunit)
       call write_inopt(streams(istream)%velocity(2),'vy_stream_'//trim(label), &
                        'y/Dec. injection velocity in km/s',iunit)
       call write_inopt(streams(istream)%velocity(3),'vz_stream_'//trim(label), &
                        'z/LOS injection velocity in km/s',iunit)
       call write_inopt(streams(istream)%ninjected,'ninjected_stream_'//trim(label), &
                        'primary particles injected (restart bookkeeping)',iunit)
    enddo
 endif

 if (use_dust) then
    if (use_dustfrac) then
       call write_inopt(dust_frac,'dust_frac','Dust fraction in smallest dust bin',iunit)
    endif
 endif
end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt,int_to_string
 use options, only:use_dustfrac
 use dim, only:use_dust
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 integer :: istream,istat
 character(len=12) :: label

 call read_inopt(imdot_func,'mdot_func',db,errcount=nerr,min=0)
 call read_inopt(stream_model,'stream_model',db,errcount=nerr,min=0,max=1,default=stream_model)
 call read_inopt(sym_stream,'sym_stream',db,errcount=nerr,min=0,max=3,default=sym_stream)

 ! Read the original scalar options first. They remain the defaults for a
 ! legacy one-stream input file that has no nstreams or indexed keys.
 call read_inopt(omega,'omega',db,errcount=nerr,min=0.,default=omega)
 call read_inopt(r0, 'r0', db,errcount=nerr,min=0.,default=r0)
 call read_inopt(phi0, 'phi0', db,errcount=nerr,default=phi0)
 call read_inopt(theta0, 'theta0', db,errcount=nerr,default=theta0)
 call read_inopt(r_inj,'r_inj',db,errcount=nerr,min=0.,default=r_inj)
 call read_inopt(vr_0,'vr_0',db,errcount=nerr,min=0.,default=vr_0)
 call read_inopt(Mdot,'Mdot',db,errcount=nerr,min=0.,default=Mdot)
 call read_inopt(stream_width,'stream_width',db,errcount=nerr,min=0.,default=stream_width)
 call read_inopt(tstart,'tstart',db,errcount=nerr,default=tstart)
 call read_inopt(tend,'tend',db,errcount=nerr,default=tend)
 call read_inopt(x_stream,'x_stream',db,errcount=nerr,default=x_stream)
 call read_inopt(y_stream,'y_stream',db,errcount=nerr,default=y_stream)
 call read_inopt(z_stream,'z_stream',db,errcount=nerr,default=z_stream)
 call read_inopt(vx_stream,'vx_stream',db,errcount=nerr,default=vx_stream)
 call read_inopt(vy_stream,'vy_stream',db,errcount=nerr,default=vy_stream)
 call read_inopt(vz_stream,'vz_stream',db,errcount=nerr,default=vz_stream)

 if (stream_model == 1) then
    call read_inopt(nstreams,'nstreams',db,errcount=nerr,min=1,default=1)
    call ensure_stream_array(nstreams,istat,reset=.true.)
    if (istat /= 0) then
       nerr = nerr+1
       return
    endif
    do istream=1,nstreams
       label = trim(int_to_string(istream))
       if (nstreams == 1) then
          call read_inopt(streams(istream)%mdot,'Mdot_stream_'//trim(label),db, &
                          errcount=nerr,min=0.,default=streams(istream)%mdot)
          call read_inopt(streams(istream)%position(1),'x_stream_'//trim(label),db, &
                          errcount=nerr,default=streams(istream)%position(1))
          call read_inopt(streams(istream)%position(2),'y_stream_'//trim(label),db, &
                          errcount=nerr,default=streams(istream)%position(2))
          call read_inopt(streams(istream)%position(3),'z_stream_'//trim(label),db, &
                          errcount=nerr,default=streams(istream)%position(3))
          call read_inopt(streams(istream)%velocity(1),'vx_stream_'//trim(label),db, &
                          errcount=nerr,default=streams(istream)%velocity(1))
          call read_inopt(streams(istream)%velocity(2),'vy_stream_'//trim(label),db, &
                          errcount=nerr,default=streams(istream)%velocity(2))
          call read_inopt(streams(istream)%velocity(3),'vz_stream_'//trim(label),db, &
                          errcount=nerr,default=streams(istream)%velocity(3))
       else
          call read_inopt(streams(istream)%mdot,'Mdot_stream_'//trim(label),db, &
                          errcount=nerr,min=0.)
          call read_inopt(streams(istream)%position(1),'x_stream_'//trim(label),db,errcount=nerr)
          call read_inopt(streams(istream)%position(2),'y_stream_'//trim(label),db,errcount=nerr)
          call read_inopt(streams(istream)%position(3),'z_stream_'//trim(label),db,errcount=nerr)
          call read_inopt(streams(istream)%velocity(1),'vx_stream_'//trim(label),db,errcount=nerr)
          call read_inopt(streams(istream)%velocity(2),'vy_stream_'//trim(label),db,errcount=nerr)
          call read_inopt(streams(istream)%velocity(3),'vz_stream_'//trim(label),db,errcount=nerr)
       endif
       call read_inopt(streams(istream)%width,'stream_width_'//trim(label),db, &
                       errcount=nerr,min=0.,default=stream_width)
       call read_inopt(streams(istream)%tstart,'tstart_stream_'//trim(label),db, &
                       errcount=nerr,default=tstart)
       call read_inopt(streams(istream)%tend,'tend_stream_'//trim(label),db, &
                       errcount=nerr,default=tend)
       call read_inopt(streams(istream)%ninjected,'ninjected_stream_'//trim(label),db, &
                       errcount=nerr,min=0,default=0)
    enddo
    if (nstreams == 1) then
       Mdot = streams(1)%mdot
       stream_width = streams(1)%width
       tstart = streams(1)%tstart
       tend = streams(1)%tend
    endif
 endif
 if (use_dust) then
    if (use_dustfrac) then
       call read_inopt(dust_frac,'dust_frac',db,errcount=nerr,default=dust_frac)
    endif
 endif
end subroutine read_options_inject

!-----------------------------------------------------------------------
!+
!  Allocate the direct-stream array. New entries inherit the original
!  scalar values so old one-stream input files remain valid.
!+
!-----------------------------------------------------------------------
subroutine ensure_stream_array(nwanted,istat,reset)
 integer, intent(in) :: nwanted
 integer, intent(out) :: istat
 logical, intent(in), optional :: reset
 logical :: reset_array
 integer :: i

 istat = 0
 reset_array = .false.
 if (present(reset)) reset_array = reset
 if (allocated(streams)) then
    if (size(streams) /= nwanted .or. reset_array) deallocate(streams)
 endif
 if (.not.allocated(streams)) then
    allocate(streams(nwanted),stat=istat)
    if (istat /= 0) return
    do i=1,nwanted
       streams(i)%mdot = Mdot
       streams(i)%width = stream_width
       streams(i)%tstart = tstart
       streams(i)%tend = tend
       streams(i)%position = (/x_stream,y_stream,z_stream/)
       streams(i)%velocity = (/vx_stream,vy_stream,vz_stream/)
    enddo
 endif
end subroutine ensure_stream_array

!-----------------------------------------------------------------------
!+
!  Relevant equations for the Mendoza+09 streamer
!+
!-----------------------------------------------------------------------

subroutine mendoza_state(mstar, r0, omega, theta0m, phi0m,  vr_0,                 &
                         r_inj, x,y,z, vx,vy,vz)
 use units, only:get_G_code
 real, intent(in)  :: mstar, r0, omega, theta0m, phi0m, r_inj, vr_0
 real, intent(out) :: x,y,z, vx,vy,vz
 real :: rc, vk0, mu, nu, eps, ecc, xi0
 real :: theta, phi, vr, vt, vp, r_rc, xi

 ! Based on the PIMS python module by Jess Speedie (https://github.com/jjspeedie/PIMS/blob/main/pims.py)

 call mendonza_invariant_parameters(mstar, r0, omega, theta0m, vr_0,                 &
                         rc, vk0, mu, nu, eps, ecc, xi0)

 ! need a better root finder, so right now r_inj must be r0
 !call theta_at_r(r_inj/rc, theta0m, ecc, xi0, theta)
 theta = theta0m
 phi = phi0m + acos( tan(theta0m) / tan(theta) )   ! Ulrich 1976, eq. (15)

!  Compute the Mondoza+09 streamer velocities
 r_rc = r_inj/rc
 xi   = acos(cos(theta)/cos(theta0m)) + xi0
 vr   = -ecc*sin(theta0m)*sin(xi)/(r_rc*(1.0 - ecc*cos(xi))) * vk0
 vt   =  sin(theta0m)/(sin(theta)*r_rc) *                      &
           sqrt(cos(theta0m)**2 - cos(theta)**2) * vk0
 vp   =  sin(theta0m)**2 /(sin(theta)*r_rc) * vk0

 ! Convert to cartesian coordinates to pass to injection routine
 x  = r_inj*sin(theta)*cos(phi)
 y  = r_inj*sin(theta)*sin(phi)
 z  = r_inj*cos(theta)
 vx = vr*sin(theta)*cos(phi) + vt*cos(theta)*cos(phi) - vp*sin(phi)
 vy = vr*sin(theta)*sin(phi) + vt*cos(theta)*sin(phi) + vp*cos(phi)
 vz = vr*cos(theta)           - vt*sin(theta)

end subroutine mendoza_state

subroutine mendonza_invariant_parameters(mstar, r0, omega, theta0m, vr_0,                 &
                         rc, vk0, mu, nu, eps, ecc, xi0)
 use units, only:get_G_code
 real, intent(in)  :: mstar, r0, omega, theta0m, vr_0
 real, intent(out) :: rc, vk0, mu, nu, eps, ecc, xi0
 real :: G_code

 ! Store the invariants here so we can access them outside of mendoza_state

 G_code = get_G_code()
 rc   = r0**4 * omega**2 / (G_code*mstar)       ! centrifugal radius
 vk0  = sqrt(G_code*mstar/rc)                   ! Keplerian speed at rc, scales other velocities
 mu   = rc/r0
 nu   = (vr_0 * sqrt(rc / (G_code * mstar)))
 eps = nu**2 + mu**2*sin(theta0m)**2 - 2.0*mu
 ecc  = sqrt(1.0 + eps*sin(theta0m)**2)
 xi0 = acos( (1 - mu*sin(theta0m)**2) / ecc ) ! Assuming purely radial motion from the sphere, different from Ulrich 1976

end subroutine mendonza_invariant_parameters

subroutine theta_at_r(r_rc, theta0, ecc, xi0, theta)
! Not used for now
 use physcon, only:pi
 real, intent(in)  :: r_rc, theta0, ecc, xi0
 real, intent(out) :: theta
 real :: a, b, m, f_a, f_m, xi
 integer  :: n
 a = theta0 ;  b = pi/2.0
 do n = 1, 60
    m  = 0.5*(a+b)
    xi = acos(cos(m)/cos(theta0)) + xi0
    f_m = r_rc - sin(theta0)**2 /(1.0 - ecc*cos(xi))
    xi  = acos(cos(a)/cos(theta0)) + xi0
    f_a = r_rc - sin(theta0)**2 /(1.0 - ecc*cos(xi))
    if (f_a*f_m <= 0.0) then
       b = m
    else
       a = m
    endif
    if (abs(b-a) < 1.0e-12) exit
 enddo
 theta = 0.5*(a+b)
 !write(*,*) 'theta = ',theta
end subroutine theta_at_r

subroutine set_default_options_inject(flag)

 integer, intent(in), optional :: flag
end subroutine set_default_options_inject

end module inject
