!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles injection for gas sphere in wind tunnel
!
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters:
!   - BHL_radius       : *radius of the wind cylinder (in code units)*
!   - Rstar            : *radius of sphere where velocities are adjusted (code units)*
!   - handled_layers   : *(integer) number of handled BHL wind layers*
!   - hold_star        : *1: subtract CM velocity of star particles at each timestep*
!   - lattice_type     : *0: cubic distribution, 1: closepacked distribution*
!   - mach             : *mach number of injected particles*
!   - rho_inf          : *ambient density (code units)*
!   - v_inf            : *wind speed (code units)*
!   - wind_injection_x : *x position of the wind injection boundary (in code units)*
!   - wind_length      : *crude wind length (in star radii)*
!
! :Dependencies: dim, eos, infile_utils, io, part, partinject, physcon,
!   units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'windtunnel'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
           set_default_options_inject,update_injected_par
!
!--runtime settings for this module
!
 ! Main parameters: model MS6 from Ruffert & Arnett (1994)
 real,    public :: v_inf = 1.
 real,    public :: rho_inf = 1.
 real,    public :: mach = 13.87

 ! Particle-related parameters
 integer, public :: lattice_type = 1
 integer, public :: handled_layers = 4
 real,    public :: wind_radius = 30.
 real,    public :: wind_injection_x = -10.
 real,    public :: wind_length = 100.

 ! option to fix the motion of a sphere of particles
 integer, public :: hold_star = 0
 real,    public :: Rstar = .1
 integer, public :: nstar  = 0

 private
 real    :: wind_rad,wind_x,psep,distance_between_layers,&
            time_between_layers,h_inf,u_inf
 integer :: max_layers,max_particles,nodd,neven,nstarpart
 logical :: first_run = .true.
 real, allocatable :: layer_even(:,:),layer_odd(:,:)

 logical, parameter :: verbose = .false.

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use eos,        only:gamma
 use part,       only:hfact,massoftype,igas
 use dim,        only:maxp
 use io,         only:fatal
 integer, intent(out) :: ierr
 real :: pmass,element_volume,y,z,cs_inf,pres_inf
 integer :: size_y, size_z, pass, i, j

 ierr = 0

 cs_inf = v_inf/mach
 pres_inf = cs_inf**2*rho_inf/gamma
 !mach = v_inf/cs_inf
 !cs_inf = sqrt(gamma*pres_inf/rho_inf)
 u_inf = pres_inf / (rho_inf*(gamma-1.))
 wind_rad = wind_radius
 wind_x = wind_injection_x
 pmass = massoftype(igas)

 ! Calculate particle separation between layers given rho_inf, depending on lattice type
 element_volume = pmass / rho_inf
 if (lattice_type == 1) then
    psep = (sqrt(2.)*element_volume)**(1./3.)
 elseif (lattice_type == 0) then
    psep = element_volume**(1./3.)
 else
    call fatal("init_inject",'unknown lattice_type (must be 0 or 1)')
 endif

 if (lattice_type == 1) then
    distance_between_layers = psep*sqrt(6.)/3.
    size_y = ceiling(3.*wind_rad/psep)
    size_z = ceiling(3.*wind_rad/(sqrt(3.)*psep/2.))
    do pass=1,2
       if (pass == 2) then
          if (allocated(layer_even)) deallocate(layer_even)
          if (allocated(layer_odd)) deallocate(layer_odd)
          allocate(layer_even(2,neven), layer_odd(2,nodd))
       endif
       neven = 0
       nodd = 0
       do i=1,size_y
          do j=1,size_z
             ! Even layer
             y = -1.5*wind_rad + (i-1)*psep
             z = -1.5*wind_rad + (j-1)*psep*sqrt(3.)/2.
             if (mod(j,2) == 0) y = y + .5*psep
             if (y**2+z**2  <  wind_rad**2) then
                neven = neven + 1
                if (pass == 2) layer_even(:,neven) = (/ y,z /)
             endif
             ! Odd layer
             y = y + psep*.5
             z = z + psep*sqrt(3.)/6.
             if (y**2+z**2  <  wind_rad**2) then
                nodd = nodd + 1
                if (pass == 2) layer_odd(:,nodd) = (/ y,z /)
             endif
          enddo
       enddo
    enddo
 else
    distance_between_layers = psep
    size_y = ceiling(3.*wind_rad/psep)
    size_z = size_y
    do pass=1,2
       if  (pass == 2) allocate(layer_even(2,neven), layer_odd(2,neven))
       neven = 0
       do i=1,size_y
          do j=1,size_z
             y = -1.5*wind_rad+(i-1)*psep
             z = -1.5*wind_rad+(j-1)*psep
             if (y**2+z**2  <  wind_rad**2) then
                neven = neven + 1
                if (pass == 2) layer_even(:,neven) = (/ y,z /)
             endif
          enddo
       enddo
    enddo
    layer_odd(:,:) = layer_even(:,:)
 endif
 h_inf = hfact*(pmass/rho_inf)**(1./3.)
 max_layers = int(wind_length/distance_between_layers)
 max_particles = int(max_layers*(nodd+neven)/2) + nstarpart
 time_between_layers = distance_between_layers/v_inf

 call print_summary(v_inf,cs_inf,rho_inf,pres_inf,mach,pmass,distance_between_layers,&
                    time_between_layers,max_layers,nstarpart,max_particles)

 if (max_particles > maxp) call fatal('windtunnel', 'maxp too small for this simulation, please increase MAXP!')

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use physcon,  only:pi
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 real :: last_time, local_time, x, irrational_number_close_to_one
 integer :: inner_layer, outer_layer, i, i_limited, i_part, np, ierr
 real, allocatable :: xyz(:,:), vxyz(:,:), h(:), u(:)

 if (first_run) then
    call init_inject(ierr)
    first_run = .false.
 endif

 last_time = time-dtlast
 outer_layer = ceiling(last_time/time_between_layers)  ! No. of layers present at t - dt
 inner_layer = ceiling(time/time_between_layers)-1 + handled_layers  ! No. of layers ought to be present at t
 ! Inject layers
 do i=outer_layer,inner_layer  ! loop over layers
    local_time = time - i*time_between_layers  ! time at which layer was injected
    i_limited = mod(i,max_layers)
    i_part = int(i_limited/2)*(nodd+neven)+mod(i_limited,2)*neven
    if (mod(i,2) == 0) then
       allocate(xyz(3,neven), vxyz(3,neven), h(neven), u(neven))
       xyz(2:3,:) = layer_even(:,:)
       np = neven
    else
       allocate(xyz(3,nodd), vxyz(3,nodd), h(nodd), u(nodd))
       xyz(2:3,:) = layer_odd(:,:)
       np = nodd
    endif
    x = wind_x + local_time*v_inf
    xyz(1,:) = x
    vxyz(1,:) = v_inf
    vxyz(2:3,:) = 0.
    h(:) = h_inf
    u(:) = u_inf
    if (verbose) then
       if (i_part  <  npart) then
          if (i  >  max_layers) then
             print *, 'Recycling (i=', i, ', max_layers=', max_layers, ', i_part=', i_part, '):'
          else
             print *, 'Moving:'
          endif
       else
          print *, 'Injecting:'
       endif
       print *, np, ' particles (npart=', npart, '/', max_particles, ')'
    endif
    call inject_or_update_particles(i_part+nstarpart+1, np, xyz, vxyz, h, u, .false.)
    deallocate(xyz, vxyz, h, u)
 enddo

 irrational_number_close_to_one = 3./pi
 dtinject = (irrational_number_close_to_one*time_between_layers)

 if (hold_star > 0) call subtract_star_vcom(nstarpart,xyzh,vxyzu)

end subroutine inject_particles

!
! Inject gas or boundary particles
!
subroutine inject_or_update_particles(ifirst, n, position, velocity, h, u, boundary)
 use part,       only:igas,iboundary,npart,npartoftype,xyzh,vxyzu
 use partinject, only:add_or_update_particle
 implicit none
 integer, intent(in) :: ifirst, n
 double precision, intent(in) :: position(3,n), velocity(3,n), h(n), u(n)
 logical, intent(in) :: boundary

 integer :: i, itype
 real :: position_u(3), velocity_u(3)

 if (boundary) then
    itype = iboundary
 else
    itype = igas
 endif

 do i=1,n
    position_u(:) = position(:,i)
    velocity_u(:) = velocity(:,i)
    call add_or_update_particle(itype,position_u,velocity_u,h(i),u(i),&
     ifirst+i-1,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine inject_or_update_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Subtracts centre-of-mass motion of star particles
!  Assumes star particles have particle IDs 1 to nsphere
!+
!-----------------------------------------------------------------------
subroutine subtract_star_vcom(nsphere,xyzh,vxyzu)
 integer, intent(in) :: nsphere
 real, intent(in)    :: xyzh(:,:)
 real, intent(inout) :: vxyzu(:,:)
 real                :: vstar(3)
 integer             :: i,nbulk

!  vstar = (/ sum(vxyzu(1,1:nsphere)), sum(vxyzu(2,1:nsphere)), sum(vxyzu(3,1:nsphere)) /) / real(nsphere)
 nbulk = 0
 vstar = 0.
 do i=1,nsphere
    if (xyzh(1,i) < 2.*Rstar) then
       vstar = vstar + vxyzu(1:3,i)
       nbulk = nbulk + 1
    endif
 enddo
 vstar = vstar/real(nbulk)

 do i=1,nsphere
    if (xyzh(1,i) < 2.*Rstar) then
       vxyzu(1:3,i) = vxyzu(1:3,i) - vstar
    endif
 enddo

end subroutine subtract_star_vcom


!-----------------------------------------------------------------------
!+
!  Print summary of wind properties (assumes inputs are in code units)
!+
!-----------------------------------------------------------------------
subroutine print_summary(v_inf,cs_inf,rho_inf,pres_inf,mach,pmass,distance_between_layers,&
                         time_between_layers,max_layers,nstar,max_particles)
 use units, only:unit_velocity,unit_pressure,unit_density
 real, intent(in)    :: v_inf,cs_inf,rho_inf,pres_inf,mach,pmass,distance_between_layers,time_between_layers
 integer, intent(in) :: max_layers,nstar,max_particles

 print*, 'wind speed: ',v_inf * unit_velocity / 1e5," km s^-1"
 print*, 'wind cs: ',cs_inf * unit_velocity / 1e5," km s^-1"
 print*, 'wind density: ',rho_inf * unit_density," g cm^-3"
 print*, 'wind pressure: ',pres_inf * unit_pressure," dyn cm^-2"
 print*, 'wind mach number: ', mach

 print*, 'maximum wind layers: ', max_layers
 print*, 'pmass: ',pmass
 print*, 'distance_between_layers: ',distance_between_layers
 print*, 'time_between_layers: ',time_between_layers

 if (hold_star > 0) then
    print*, 'planet crossing time: ',2*Rstar/v_inf
    print*, 'wind impact time: ',(abs(wind_injection_x) - Rstar)/v_inf
    print*, 'nstar: ',nstar
 endif
 print*, 'nstar + max. wind particles: ', max_particles

end subroutine print_summary


!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(v_inf,'v_inf','wind speed (code units)',iunit)
 call write_inopt(mach,'mach','mach number of injected particles',iunit)
 call write_inopt(rho_inf,'rho_inf','ambient density (code units)',iunit)
 call write_inopt(lattice_type,'lattice_type','0: cubic distribution, 1: closepacked distribution',iunit)
 call write_inopt(handled_layers,'handled_layers','(integer) number of handled BHL wind layers',iunit)
 call write_inopt(hold_star,'hold_star','1: subtract CM velocity of star particles at each timestep',iunit)
 if (hold_star > 0) then
    call write_inopt(Rstar,'Rstar','radius of sphere where velocities are adjusted (code units)',iunit)
    call write_inopt(nstar,'nstar','No. of particles that should have their velocity adjusted',iunit)  ! need to write actual no. of particles, not nstar_in
 endif
 call write_inopt(wind_radius,'BHL_radius','radius of the wind cylinder (in code units)',iunit)
 call write_inopt(wind_injection_x,'wind_injection_x','x position of the wind injection boundary (in code units)',iunit)
 call write_inopt(wind_length,'wind_length','crude wind length (in star radii)',iunit)

end subroutine write_options_inject


!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('v_inf')
    read(valstring,*,iostat=ierr) v_inf
    ngot = ngot + 1
    if (v_inf <= 0.)    call fatal(label,'v_inf must be positive')
 case('mach')
    read(valstring,*,iostat=ierr) mach
    ngot = ngot + 1
    if (mach <= 0.) call fatal(label,'mach must be positive')
 case('rho_inf')
    read(valstring,*,iostat=ierr) rho_inf
    ngot = ngot + 1
    if (rho_inf <= 0.) call fatal(label,'rho_inf must be positive')
 case('nstar')
    read(valstring,*,iostat=ierr) nstar
    ngot = ngot + 1
 case('Rstar')
    read(valstring,*,iostat=ierr) Rstar
    ngot = ngot + 1
    if (Rstar <= 0.)    call fatal(label,'invalid setting for Rstar (<=0)')
 case('lattice_type')
    read(valstring,*,iostat=ierr) lattice_type
    ngot = ngot + 1
    if (lattice_type/=0 .and. lattice_type/=1)    call fatal(label,'lattice_type must be 0 or 1')
 case('handled_layers')
    read(valstring,*,iostat=ierr) handled_layers
    ngot = ngot + 1
    if (handled_layers < 0) call fatal(label,'handled_layers must be positive or zero')
 case('BHL_radius')
    read(valstring,*,iostat=ierr) wind_radius
    ngot = ngot + 1
    if (wind_radius <= 0.) call fatal(label,'wind_radius must be >0')
 case('wind_injection_x')
    read(valstring,*,iostat=ierr) wind_injection_x
    ngot = ngot + 1
 case('wind_length')
    read(valstring,*,iostat=ierr) wind_length
    ngot = ngot + 1
    if (wind_length <= 0.) call fatal(label,'wind_length must be positive')
 case('hold_star')
    read(valstring,*,iostat=ierr) hold_star
    ngot = ngot + 1
 end select

 igotall = (ngot >= 9)
end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
