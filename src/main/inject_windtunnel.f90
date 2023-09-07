!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles injection for gas sphere in wind tunnel
!
!
! :Owner: Mike Lau
!
! :Runtime parameters:
!   - lattice_type         : *0: cubic distribution, 1: closepacked distribution*
!   - BHL_handled_layers   : *(integer) number of handled BHL wind layers*
!   - mach_inf             : *BHL wind mach number*
!   - Rstar                : *BHL star radius (in accretion radii)*
!   - BHL_radius           : *radius of the wind cylinder (in star radii)*
!   - BHL_wind_injection_x : *x position of the wind injection boundary (in star radii)*
!   - BHL_wind_length      : *crude wind length (in star radii)*
!
! :Dependencies: dim, eos, infile_utils, io, part, partinject, physcon,
!   units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'windtunnel'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject
!
!--runtime settings for this module
!
 ! Main parameters: model MS6 from Ruffert & Arnett (1994)
 real,    public :: mach_inf = 1.
 real,    public :: rho_inf = 1.
 real,    public :: cs_inf = 1.
 real,    public :: Rstar = .1
 integer, public :: nstar  = 0

 ! Particle-related parameters
 integer, public :: lattice_type = 1
 real,    public :: BHL_handled_layers = 4.
 real,    public :: BHL_wind_cylinder_radius = 30.
 real,    public :: BHL_wind_injection_x = -10.
 real,    public :: BHL_wind_length = 100.

 private
 integer :: handled_layers
 real    :: wind_cylinder_radius,wind_injection_x,psep,distance_between_layers,&
            time_between_layers,h_inf,u_inf,v_inf
 integer :: max_layers,max_particles,nodd,neven
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
 use physcon,    only:gg,pi
 use eos,        only:gamma
 use part,       only:hfact,massoftype,igas
 use dim,        only:maxp
 use io,         only:fatal
 integer, intent(out) :: ierr
 real :: pmass,element_volume,y,z
 integer :: size_y, size_z, pass, i, j

 ierr = 0

 v_inf = mach_inf*cs_inf
 u_inf = cs_inf**2 / (gamma*(gamma-1.))
 handled_layers = int(BHL_handled_layers)
 wind_cylinder_radius = BHL_wind_cylinder_radius * Rstar
 wind_injection_x = BHL_wind_injection_x * Rstar
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
    size_y = ceiling(3.*wind_cylinder_radius/psep)
    size_z = ceiling(3.*wind_cylinder_radius/(sqrt(3.)*psep/2.))
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
             y = -1.5*wind_cylinder_radius + (i-1)*psep
             z = -1.5*wind_cylinder_radius + (j-1)*psep*sqrt(3.)/2.
             if (mod(j,2) == 0) y = y + .5*psep
             if (y**2+z**2  <  wind_cylinder_radius**2) then
                neven = neven + 1
                if (pass == 2) layer_even(:,neven) = (/ y,z /)
             endif
             ! Odd layer
             y = y + psep*.5
             z = z + psep*sqrt(3.)/6.
             if (y**2+z**2  <  wind_cylinder_radius**2) then
                nodd = nodd + 1
                if (pass == 2) layer_odd(:,nodd) = (/ y,z /)
             endif
          enddo
       enddo
    enddo
 else
    distance_between_layers = psep
    size_y = ceiling(3.*wind_cylinder_radius/psep)
    size_z = size_y
    do pass=1,2
       if  (pass == 2) allocate(layer_even(2,neven), layer_odd(2,neven))
       neven = 0
       do i=1,size_y
          do j=1,size_z
             y = -1.5*wind_cylinder_radius+(i-1)*psep
             z = -1.5*wind_cylinder_radius+(j-1)*psep
             if (y**2+z**2  <  wind_cylinder_radius**2) then
                neven = neven + 1
                if (pass == 2) layer_even(:,neven) = (/ y,z /)
             endif
          enddo
       enddo
    enddo
    layer_odd(:,:) = layer_even(:,:)
 endif
 max_layers = int(BHL_wind_length*Rstar/distance_between_layers)
 max_particles = int(max_layers*(nodd+neven)/2) + nstar
 print *, 'BHL maximum layers: ', max_layers
 print *, 'BHL maximum particles: ', max_particles
 print *, 'nstar: ',nstar  
 if (max_particles > maxp) call fatal('BHL', 'maxp too small for this simulation, please increase MAXP!')
 time_between_layers = distance_between_layers/v_inf
 print *, 'distance_between_layers: ',distance_between_layers  
 print *, 'time_between_layers: ',time_between_layers  
 print *, 'pmass: ',pmass
 h_inf = hfact*(pmass/rho_inf)**(1./3.)
 !if (setup) then
!    tmax = (100.*abs(wind_injection_x)/v_inf)/utime
! endif

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use physcon,  only:gg,pi
 use units,    only:utime
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
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
    x = wind_injection_x + local_time*v_inf
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
    call inject_or_update_particles(i_part+nstar+1, np, xyz, vxyz, h, u, .false.)
    deallocate(xyz, vxyz, h, u)
 enddo

 irrational_number_close_to_one = 3./pi
 dtinject = (irrational_number_close_to_one*time_between_layers)/utime

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

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(mach_inf,'mach_inf','BHL wind mach number',iunit)
 call write_inopt(cs_inf,'cs_inf','ambient sound speed',iunit)
 call write_inopt(rho_inf,'rho_inf','ambient density',iunit)
 call write_inopt(Rstar,'Rstar','BHL star radius (in accretion radii)',iunit)
 call write_inopt(nstar,'nstar','No. of particles making up star',iunit)
 call write_inopt(lattice_type,'lattice_type','0: cubic distribution, 1: closepacked distribution',iunit)
 call write_inopt(BHL_handled_layers,'BHL_handled_layers','(integer) number of handled BHL wind layers',iunit)
 call write_inopt(BHL_wind_cylinder_radius,'BHL_radius','radius of the wind cylinder (in star radii)',iunit)
 call write_inopt(BHL_wind_injection_x,'BHL_wind_injection_x','x position of the wind injection boundary (in star radii)',iunit)
 call write_inopt(BHL_wind_length,'BHL_wind_length','crude wind length (in star radii)',iunit)

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
 case('mach_inf')
    read(valstring,*,iostat=ierr) mach_inf
    ngot = ngot + 1
    if (mach_inf <= 0.)    call fatal(label,'invalid setting for mach_inf (<=0)')
 case('cs_inf')
    read(valstring,*,iostat=ierr) cs_inf
    ngot = ngot + 1
    if (cs_inf <= 0.) call fatal(label,'cs_inf must be positive')
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
 case('BHL_handled_layers')
    read(valstring,*,iostat=ierr) BHL_handled_layers
    ngot = ngot + 1
    if (dble(int(BHL_handled_layers))  /=  BHL_handled_layers) call fatal(label,'BHL_handled_layers must be integer')
    if (int(BHL_handled_layers)  <  0) call fatal(label,'BHL_handled_layers must be positive or zero')
 case('BHL_radius')
    read(valstring,*,iostat=ierr) BHL_wind_cylinder_radius
    ngot = ngot + 1
    if (BHL_wind_cylinder_radius <= 0.) call fatal(label,'BHL_wind_cylinder_radius must be >0')
 case('BHL_wind_injection_x')
    read(valstring,*,iostat=ierr) BHL_wind_injection_x
    ngot = ngot + 1
 case('BHL_wind_length')
    read(valstring,*,iostat=ierr) BHL_wind_length
    ngot = ngot + 1
    if (BHL_wind_length <= 0.) call fatal(label,'BHL_wind_length must be positive')
 end select

 igotall = (ngot >= 10)
end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
