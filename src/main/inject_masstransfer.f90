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
! :Owner: Ana Juarez and Mike Lau
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
 character(len=*), parameter, public :: inject_type = 'masstransfer'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
           set_default_options_inject,update_injected_par
!
!--runtime settings for this module
!
 ! Main parameters: model MS6 from Ruffert & Arnett (1994) for windtunnel
 real,    public :: v_inf = 1.
 real,    public :: rho_inf = 1.
 real,    public :: mach = 13.87

 ! Particle-related parameters
 integer, public :: lattice_type = 0
 integer, public :: handled_layers = 1
 real,    public :: wind_radius = 30.
 real,    public :: wind_injection_x = -10.
 real,    public :: wind_length = 100.

 ! option to fix the motion of a sphere of particles
 integer, public :: hold_star = 0
 real,    public :: Rstar = .1
 integer, public :: nstar  = 0

 !file with mesa values
 character(len=120), public :: filemesa='./test_data.txt'

 private
 real    :: wind_rad,wind_x,psep,time_between_layers,h_inf,u_inf
 logical :: first_run = .true.
 real, allocatable,dimension(:) :: time_mesa, mdot_mesa !Arrays from mesa file

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
 use physcon,    only:pi
 use io,         only:fatal
 use readwrite_mesa,  only:read_masstransferrate
 integer, intent(out) :: ierr
 integer :: nodd,neven,imax_layers,imax_particles
 real :: pmass,cs_inf,pres_inf,rho_inf,distance_between_layers
 real, allocatable :: layer_even(:,:),layer_odd(:,:)

 ierr = 0

 !--Read mesa file
 call read_masstransferrate(filemesa,time_mesa,mdot_mesa,ierr)
 wind_rad = wind_radius
 rho_inf = mdot_mesa(1)/(pi*wind_rad**2*v_inf)
 cs_inf = v_inf/mach
 pres_inf = cs_inf**2*rho_inf/gamma
 u_inf = pres_inf / (rho_inf*(gamma-1.))
 wind_x = wind_injection_x
 pmass = massoftype(igas)

 h_inf = hfact*(pmass/rho_inf)**(1./3.)

 ! Calculate particle separation between layers given rho_inf, depending on lattice type
 call calculate_lattice(lattice_type,rho_inf,pmass,wind_rad,imax_layers,imax_particles, &
      time_between_layers,nodd,neven,layer_even,layer_odd,distance_between_layers)

 call print_summary(v_inf,cs_inf,rho_inf,pres_inf,mach,pmass,distance_between_layers,&
                    time_between_layers,imax_layers,nstar,imax_particles)

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use physcon,          only:pi
 use part,             only:nptmass,delete_dead_particles_inside_radius,delete_particles_outside_sphere,igas
 use extern_corotate,  only: companion_xpos,primarycore_xpos,primarycore_hsoft,hsoft
 use table_utils,      only: find_nearest_index,interp_1d
 use eos,        only:gamma
 use part,       only:hfact,massoftype,igas
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 real :: last_time, local_time, x, irrational_number_close_to_one
 real ::  xyz_acc(3),xyz_don(3),racc,rdon
 integer :: inner_layer, outer_layer, i, i_limited, i_part, np, ierr
 integer :: time_index,t1,t2,nodd,neven
 integer :: imax_layers,imax_particles !LS these variables are not used anymore
 real, allocatable :: xyz(:,:), vxyz(:,:), h(:), u(:), mdot
 real :: pmass,cs_inf,pres_inf,kill_sep,distance_between_layers!,time_between_layers_last
 real, allocatable :: layer_even(:,:),layer_odd(:,:)

   if (first_run) then
      call init_inject(ierr)
      first_run = .false.
   endif

 !Interpolation of the mass transfer rate from the mesa file
 call find_nearest_index(time_mesa,time,time_index)
 t1 = time_index
 t2 = time_index + 1
 mdot = interp_1d(time,time_mesa(t1),time_mesa(t2),mdot_mesa(t1),mdot_mesa(t2))
 wind_rad = wind_radius
 !time_between_layers_last = time_between_layers  !unused
 rho_inf = mdot/(pi*wind_rad**2*v_inf)
 cs_inf = v_inf/mach
 pmass = massoftype(igas)
 pres_inf = cs_inf**2*rho_inf/gamma
 u_inf = pres_inf / (rho_inf*(gamma-1.))
 h_inf = hfact*(pmass/rho_inf)**(1./3.)

 call calculate_lattice(lattice_type,rho_inf,pmass,wind_rad,imax_layers,imax_particles, &
      time_between_layers,nodd,neven,layer_even,layer_odd,distance_between_layers)

 last_time = time-dtlast
 outer_layer = ceiling(last_time/time_between_layers)  ! No. of layers present at t - dt
 inner_layer = ceiling(time/time_between_layers)-1 + handled_layers  ! No. of layers ought to be present at t

 ! enforce ejection of a shell at the start of the simulation otherwise splash complains that the dump is empty
 if (npart == 0) inner_layer = 1

 ! Inject layers
 do i=outer_layer,inner_layer-1  ! loop over layers
    local_time = time - i*time_between_layers  ! time at which layer was injected
    !no need anymore
    !i_limited = mod(i,imax_layers)
    !i_part = int(i_limited/2)*(nodd+neven)+mod(i_limited,2)*neven
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
          if (i  >  imax_layers) then
             print *, 'Recycling (i=', i, ', max_layers=', imax_layers, ', i_part=', i_part, '):'
          else
             print *, 'Moving:'
          endif
       else
          print *, 'Injecting:'
       endif
       print *, np, ' particles (npart=', npart, '/', imax_particles, ')'
    endif
    !print *,'@@ inject ',i,'npart=',npart,'np=',np,neven,nodd,'x=',x,wind_x,local_time,time_between_layers
    !call inject_or_update_particles(i_part+nstar+1, np, xyz, vxyz, h, u, .false.)
    call inject_or_update_particles(np, xyz, vxyz, h, u, .false.)
    deallocate(xyz, vxyz, h, u)
 enddo

!ADD SUBROUTINE DELETE PARTICLES
 if (nptmass == 0) then
    xyz_acc = (/companion_xpos,0.,0./)
    xyz_don = (/primarycore_xpos,0.,0./)
    racc = hsoft
    rdon = primarycore_hsoft
    do i=1,npart
       call delete_particles_inside_sphere(xyz_don,rdon,xyzh(1:3,i),xyzh(4,i))
       call delete_particles_inside_sphere(xyz_acc,racc,xyzh(1:3,i),xyzh(4,i))
    enddo
    kill_sep = (abs(wind_x*2.)+abs(primarycore_xpos))  !LS to be adjusted
    call delete_particles_outside_sphere(xyz_acc,kill_sep,npart)
 endif

 irrational_number_close_to_one = 3./pi
 dtinject = (irrational_number_close_to_one*time_between_layers)

 if (hold_star > 0) call subtract_star_vcom(nstar,xyzh,vxyzu)

end subroutine inject_particles

!
! Inject gas or boundary particles
!
!subroutine inject_or_update_particles(ifirst, n, position, velocity, h, u, boundary)
subroutine inject_or_update_particles(n, position, velocity, h, u, boundary)
 use part,       only:igas,iboundary,npart,npartoftype,xyzh,vxyzu
 use partinject, only:add_or_update_particle
 implicit none
 integer, intent(in) :: n! ifirst
 real,    intent(in) :: position(3,n), velocity(3,n), h(n), u(n)
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
     npart+1,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine inject_or_update_particles

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!----------------------------------------------------------------
!+
!  Calculate the particle positions inside the wind
!  to get the require density (rho). It depends on the lattice type
!  0: cubic distribution, 1: closepacked distribution
!+
!----------------------------------------------------------------

subroutine calculate_lattice(ilattice,rho,pmass,radius,imax_layers,imax_particles, &
     time_between_layers,nodd,neven,positions_layer_even,positions_layer_odd,distance_between_layers)
 use io,         only:fatal
 real, intent(in) :: rho,pmass,radius
 integer, intent(in) :: ilattice
 integer, intent(out) :: imax_layers,imax_particles,nodd,neven
 real, intent(out) :: time_between_layers,distance_between_layers
 real, allocatable, intent(out) :: positions_layer_even(:,:),positions_layer_odd(:,:)
 real, allocatable :: tmp_pos(:,:)
 real :: element_volume,y,z,r2
 real, parameter :: sq3_2=sqrt(3.)/2.,sq3_6=sqrt(3.)/6.
! real :: psep_ini=1.0
 integer :: size_y,size_z,nz,pass,i,j

 element_volume = pmass / rho
 if (ilattice == 1) then
    psep = (sqrt(2.)*element_volume)**(1./3.)
!    psep_ini = (sqrt(2.)*(pmass/1.705e-8))**(1./3.)
 elseif (ilattice == 0) then
    psep = element_volume**(1./3.)
!    psep_ini =  (pmass / 1.705e-8)**(1./3.)
 else
    call fatal("init_inject",'unknown ilattice (must be 0 or 1)')
 endif
 r2 = radius**2
 if (allocated(positions_layer_even)) deallocate(positions_layer_even)
 if (allocated(positions_layer_odd)) deallocate(positions_layer_odd)

 if (ilattice == 1) then
    distance_between_layers = psep*sqrt(6.)/3.
    size_y = ceiling(3.*radius/psep)
    size_z = ceiling(3.*radius/(sqrt(3.)*psep/2.))
    ! do pass=1,2
    !    if (pass == 2) then
    !       if (allocated(positions_layer_even)) deallocate(positions_layer_even)
    !       if (allocated(positions_layer_odd)) deallocate(positions_layer_odd)
    !       allocate(positions_layer_even(2,neven), positions_layer_odd(2,nodd))
    !    endif
    !    neven = 0
    !    nodd = 0
    !    do i=1,size_y
    !       do j=1,size_z
    !          ! Even layer
    !          y = -1.5*radius + (i-1)*psep
    !          z = -1.5*radius + (j-1)*psep*sq3_2
    !          if (mod(j,2) == 0) y = y + .5*psep
    !          if (y**2+z**2  <  r2) then
    !             neven = neven + 1
    !             if (pass == 2) positions_layer_even(:,neven) = (/ y,z /)
    !          endif
    !          ! Odd layer
    !          y = y + psep*.5
    !          z = z + psep*sq3_6
    !          if (y**2+z**2  <  r2) then
    !             nodd = nodd + 1
    !             if (pass == 2) positions_layer_odd(:,nodd) = (/ y,z /)
    !          endif
    !       enddo
    !    enddo
    ! enddo
    nz = size_z*size_z
    allocate(tmp_pos(2,nz))
    neven = 0
    nodd = 0
    do i=1,size_y
       do j=1,size_z
          ! Even layer
          y = -1.5*radius + (i-1)*psep
          z = -1.5*radius + (j-1)*psep*sq3_2
          if (mod(j,2) == 0) y = y + .5*psep
          if (y**2+z**2  <  r2) then
             neven = neven + 1
             tmp_pos(:,neven) = (/ y,z /)
          endif
          ! Odd layer
          y = y + psep*.5
          z = z + psep*sq3_6
          if (y**2+z**2  <  r2) then
             nodd = nodd + 1
             tmp_pos(:,nz-nodd+1) = (/ y,z /)
          endif
       enddo
    enddo
    if (neven+nodd > nz) stop 'memory allocation problem'
    allocate(positions_layer_even(2,neven), positions_layer_odd(2,nodd))
    positions_layer_even(:,1:neven) = tmp_pos(:,1:neven)
    do i = 1,nodd
       positions_layer_odd(:,i) = tmp_pos(:,nz-i+1)
    enddo
    deallocate(tmp_pos)
 else
    distance_between_layers = psep
    size_y = ceiling(3.*radius/psep)
    size_z = size_y
    ! do pass=1,2
    !    if  (pass == 2) allocate(positions_layer_even(2,neven), positions_layer_odd(2,neven))
    !    neven = 0
    !    do i=1,size_y
    !       do j=1,size_z
    !          y = -1.5*radius+(i-1)*psep
    !          z = -1.5*radius+(j-1)*psep
    !          if (y**2+z**2  <  radius**2) then
    !             neven = neven + 1
    !             if (pass == 2) positions_layer_even(:,neven) = (/ y,z /)
    !          endif
    !       enddo
    !    enddo
    ! enddo
    nz = size_y*size_z
    allocate(tmp_pos(2,nz))
    neven = 0
    do i=1,size_y
       do j=1,size_z
          y = -1.5*radius+(i-1)*psep
          z = -1.5*radius+(j-1)*psep
          if (y**2+z**2  <  r2) then
             neven = neven + 1
             tmp_pos(:,neven) = (/ y,z /)
          endif
       enddo
    enddo
    allocate(positions_layer_even(2,neven), positions_layer_odd(2,neven))
    positions_layer_even(:,1:neven) = tmp_pos(:,1:neven)
    deallocate(tmp_pos)
    positions_layer_odd(:,:) = positions_layer_even(:,:)
    nodd = neven
 endif

 time_between_layers = distance_between_layers/v_inf

!these variables are not needed anymore
 imax_particles = int(imax_layers*(nodd+neven)/2) + nstar
 imax_layers = int(wind_length/distance_between_layers)

end subroutine calculate_lattice

!----------------------------------------------------------------
!+
!  Delete particles inside of a defined sphere
!+
!----------------------------------------------------------------
subroutine delete_particles_inside_sphere(center,radius,xyzi,hi)
 real,    intent(in)    :: center(3),radius,xyzi(3)
 real, intent(inout) :: hi
 real    :: r(3), radius_squared

 radius_squared = radius**2

 r = xyzi(1:3) - center
 if (dot_product(r,r) < radius_squared) then
    hi = -abs(hi)
 endif


end subroutine delete_particles_inside_sphere

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
                         time_between_layers,imax_layers,nstar,imax_particles)
 use units, only:unit_velocity,unit_pressure,unit_density
 real, intent(in)    :: v_inf,cs_inf,rho_inf,pres_inf,mach,pmass,distance_between_layers,time_between_layers
 integer, intent(in) :: imax_layers,nstar,imax_particles

 print*, 'wind speed: ',v_inf * unit_velocity / 1e5," km s^-1"
 print*, 'wind cs: ',cs_inf * unit_velocity / 1e5," km s^-1"
 print*, 'wind density: ',rho_inf * unit_density," g cm^-3"
 print*, 'wind pressure: ',pres_inf * unit_pressure," dyn cm^-2"
 print*, 'wind mach number: ', mach

 print*, 'maximum wind layers: ', imax_layers
 print*, 'pmass: ',pmass
 print*, 'distance_between_layers: ',distance_between_layers
 print*, 'time_between_layers: ',time_between_layers

 if (hold_star > 0) then
    print*, 'planet crossing time: ',2*Rstar/v_inf
    print*, 'wind impact time: ',(abs(wind_injection_x) - Rstar)/v_inf
    print*, 'nstar: ',nstar
 endif
 print*, 'nstar + max. wind particles: ', imax_particles

end subroutine print_summary


!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
! use part,     only:nptmass
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
 case('BHL_radius')
    read(valstring,*,iostat=ierr) wind_radius
    ngot = ngot + 1
    if (wind_radius <= 0.) call fatal(label,'wind_radius must be >0')
 case('wind_injection_x')
    read(valstring,*,iostat=ierr) wind_injection_x
    ngot = ngot + 1
 case('hold_star')
    read(valstring,*,iostat=ierr) hold_star
    ngot = ngot + 1

!LS the parameters below can be removed
 case('handled_layers')
    read(valstring,*,iostat=ierr) handled_layers
    ngot = ngot + 1
    if (handled_layers < 0) call fatal(label,'handled_layers must be positive or zero')
 case('wind_length')
    read(valstring,*,iostat=ierr) wind_length
    ngot = ngot + 1
    if (wind_length <= 0.) call fatal(label,'wind_length must be positive')
 end select

 igotall = (ngot >= 9)
end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
