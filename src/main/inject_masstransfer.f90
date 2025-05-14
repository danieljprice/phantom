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
!   - wind_radius      : *radius of the wind cylinder (in code units)*
!   - handled_layers   : *(integer) number of handled BHL wind layers*
!   - lattice_type     : *0: cubic distribution, 1: closepacked distribution*
!   - mach             : *mach number of injected particles*
!   - rho_inf          : *ambient density (code units)*
!   - v_inf            : *wind speed (code units)*
!   - wind_injection_x : *x position of the wind injection boundary (in code units)*
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
 real,    public :: v_inf = 1.
 real,    public :: mach = 13.87
 real,    public :: mdot_msun_yr = 1.e-4

 ! Particle-related parameters
 integer, public :: lattice_type = 1
 integer, public :: handled_layers = 4
 real,    public :: wind_radius = 30.
 real,    public :: wind_injection_x = -10.

 ! mdot vs. time data file
 logical, public            :: use_mesa_file=.false.
 character(len=120), public :: filemesa='./test_data.txt'

 private
 logical              :: first_run=.true.,even_layer
 real                 :: mdot,psep,h_inf,u_inf
 logical, parameter   :: verbose=.false.
 integer              :: nlayer_max  ! memory (no. of particles) allocated to y_layer and z_layer arrays
 real, allocatable    :: time_mesa(:),mdot_mesa(:),injection_time(:),y_layer(:,:),z_layer(:,:)
 integer, allocatable :: nlayer(:),ifirst(:)

contains

subroutine init_inject(ierr)
 integer, intent(out) :: ierr

 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject_masstransfer(time,dtlast,ierr)
 use part,            only:massoftype,igas
 use io,              only:fatal
 use readwrite_mesa,  only:read_masstransferrate
 use physcon,         only:solarm,years
 use units,           only:utime,umass
 real, intent(in)     :: time,dtlast
 integer, intent(out) :: ierr
 integer              :: nodd,neven
 real                 :: pmass,cs_inf,pres_inf,rho_inf,distance_between_layers,time_between_layers
 real, allocatable    :: layer_even(:,:),layer_odd(:,:)

 ierr = 0

 if (use_mesa_file) then
    call read_masstransferrate(filemesa,time_mesa,mdot_mesa,ierr)
    call interpolate_mdot(time,time_mesa,mdot_mesa,mdot)
 else
    mdot = mdot_msun_yr * solarm / years * utime / umass
 endif
 call calc_wind_properties(mdot,wind_radius,v_inf,mach,rho_inf,pres_inf,cs_inf,u_inf,h_inf)

 ! Calculate particle separation between layers given rho_inf, depending on lattice type
 pmass = massoftype(igas)
 call calculate_lattice(lattice_type,rho_inf,pmass,wind_radius,&
      time_between_layers,nodd,neven,layer_even,layer_odd,distance_between_layers)

 nlayer_max = 2000
 even_layer = .true.  ! choose first injected layer to be an even layer
 allocate(nlayer(handled_layers),ifirst(handled_layers),injection_time(handled_layers))
 allocate(y_layer(nlayer_max,handled_layers),z_layer(nlayer_max,handled_layers))

 ! initialise arrays
 injection_time = 0.
 nlayer = 0
 y_layer = 0.
 z_layer = 0.
 injection_time(handled_layers) = time - time_between_layers ! ensures first layer is injected after one time step. Do not inject at t=0

 call print_summary(v_inf,cs_inf,rho_inf,pres_inf,mach,pmass,distance_between_layers,&
                    time_between_layers)

 first_run = .false.

end subroutine init_inject_masstransfer

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use io,               only:fatal
 use readwrite_mesa,   only:read_masstransferrate
 use part,             only:nptmass,delete_dead_particles_inside_radius,igas
 use extern_corotate,  only:companion_xpos,primarycore_xpos,primarycore_hsoft,hsoft
 use part,             only:massoftype,igas
 use partinject,       only:add_or_update_particle
 real, intent(in)       :: time,dtlast
 real, intent(inout)    :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real, intent(out)      :: dtinject
 real                   :: irrational_number_close_to_one,xyz_acc(3),xyz_don(3),racc,rdon,xi,xyzi(3)
 real                   :: pmass,cs_inf,rho_inf,pres_inf,kill_rad,time_between_layers,distance_between_layers
 integer                :: i,k,ierr,nodd,neven
 real, allocatable      :: xyz(:,:),layer_even(:,:),layer_odd(:,:)

 if (first_run) call init_inject_masstransfer(time,dtlast,ierr)

 if (use_mesa_file) call interpolate_mdot(time,time_mesa,mdot_mesa,mdot)
 call calc_wind_properties(mdot,wind_radius,v_inf,mach,rho_inf,pres_inf,cs_inf,u_inf,h_inf)

 pmass = massoftype(igas)
 call calculate_lattice(lattice_type,rho_inf,pmass,wind_radius,&
      time_between_layers,nodd,neven,layer_even,layer_odd,distance_between_layers)

 if (max(nodd,neven) > nlayer_max) then
    call fatal('inject_particles','number of particles to be added in layer exceeds array size (nlayer > nlayer_max)')
 endif

 ! inject new layer
 if (time - injection_time(handled_layers) >= time_between_layers) then

    ! shift array entries for old layers
    nlayer(1:handled_layers-1) = nlayer(2:handled_layers)
    injection_time(1:handled_layers-1) = injection_time(2:handled_layers)
    ifirst(1:handled_layers-1) = ifirst(2:handled_layers)
    y_layer(:,1:handled_layers-1) = y_layer(:,2:handled_layers)
    z_layer(:,1:handled_layers-1) = z_layer(:,2:handled_layers)

    if (even_layer) then
       allocate(xyz(3,neven))
       xyz(2:3,:) = layer_even(:,:)
       nlayer(handled_layers) = neven
    else
       allocate(xyz(3,nodd))
       xyz(2:3,:) = layer_odd(:,:)
       nlayer(handled_layers) = nodd
    endif
    
    ifirst(handled_layers) = npart + 1  ! record id of first particle in new layer
    injection_time(handled_layers) = injection_time(handled_layers-1) + time_between_layers  ! record injection time of new layer
    y_layer(1:nlayer(handled_layers),handled_layers) = xyz(2,:)  ! y and z positions of all particles in new layer
    z_layer(1:nlayer(handled_layers),handled_layers) = xyz(3,:)
    deallocate(xyz)

    if (verbose) then
       print*
       print*
       print*,'INJECTING LAYER'
       print*,'t=',time,'time_between_layers=',time_between_layers
       print*,'injection_time=',injection_time
       print*,'iseven=',even_layer,'nlayer=',nlayer
       read*
    endif

    even_layer = .not. even_layer  ! change odd/even-ness for next layer
 endif


 !handle layers
 if (verbose) then
    print*
    print*
    print*,'HANDLING PARTICLES'
    print*,'  injection_time=',injection_time
    print*,'  ifirst=',ifirst
    print*,'  nlayer=',nlayer
    read*
 endif
 do k = 1,handled_layers  ! loop from oldest to the newest layers
    if (nlayer(k) == 0) cycle  ! skip layers that have not yet been injected
    if (verbose) then
       print*
       print*,'Injecting layer ',k
    endif
    xi = wind_injection_x + v_inf * (time-injection_time(k))
    do i = 1,nlayer(k)  ! loop over particles in layer
       xyzi = (/xi, y_layer(i,k), z_layer(i,k)/)
       call add_or_update_particle(igas,xyzi,(/v_inf,0.,0./),h_inf,u_inf,ifirst(k)+i-1,npart,npartoftype,xyzh,vxyzu)
    enddo
 enddo

! delete particles
 if (nptmass == 0) then
    xyz_acc = (/companion_xpos,0.,0./)
    xyz_don = (/primarycore_xpos,0.,0./)
    racc = hsoft
    rdon = primarycore_hsoft
    kill_rad = 2000.+abs(primarycore_xpos)
    do i=1,npart
       call delete_particles_inside_or_outside_sphere(xyz_don,rdon,xyzh(1:3,i),xyzh(4,i))
       call delete_particles_inside_or_outside_sphere(xyz_acc,racc,xyzh(1:3,i),xyzh(4,i))
       call delete_particles_inside_or_outside_sphere(xyz_acc,kill_rad,xyzh(1:3,i),xyzh(4,i),revert=.true.)
    enddo
 endif

 irrational_number_close_to_one = 3./3.1415926536
 dtinject = (irrational_number_close_to_one*time_between_layers)

end subroutine inject_particles


subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par


!----------------------------------------------------------------
!+
!  Calculate wind properties, assuming input is in code units
!+
!----------------------------------------------------------------
subroutine calc_wind_properties(mdot_in,wind_rad_in,vinf_in,mach_in,rho_inf,pres_inf,cs_inf,uinf_out,h_out)
 use physcon, only:pi
 use eos,     only:gamma
 use part,    only:hfact,massoftype,igas
 real, intent(in)  :: mdot_in,wind_rad_in,vinf_in,mach_in
 real, intent(out) :: rho_inf,pres_inf,cs_inf,uinf_out,h_out
 real              :: pmass

 rho_inf = mdot_in/(pi*wind_rad_in**2*vinf_in)
 cs_inf = vinf_in/mach_in
 pmass = massoftype(igas)
 pres_inf = cs_inf**2*rho_inf/gamma
 uinf_out = pres_inf / (rho_inf*(gamma-1.))
 h_out = hfact*(pmass/rho_inf)**(1./3.)

end subroutine calc_wind_properties


!----------------------------------------------------------------
!+
!  Calculate the particle positions inside the wind
!  to get the require density (rho). It depends on the lattice type
!  0: cubic distribution, 1: closepacked distribution
!+
!----------------------------------------------------------------
subroutine calculate_lattice(ilattice,rho,pmass,radius,time_between_layers,nodd,neven,&
                             positions_layer_even,positions_layer_odd,distance_between_layers)
 use io, only:fatal
 real, intent(in)               :: rho,pmass,radius
 integer, intent(in)            :: ilattice
 integer, intent(out)           :: nodd,neven
 real, intent(out)              :: time_between_layers,distance_between_layers
 real, allocatable, intent(out) :: positions_layer_even(:,:),positions_layer_odd(:,:)
 real, allocatable              :: tmp_pos(:,:)
 real                           :: element_volume,y,z,r2
 real, parameter                :: sq3_2=sqrt(3.)/2.,sq3_6=sqrt(3.)/6.
 integer                        :: size_y,size_z,nz,i,j
! real :: psep_ini=1.0

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

end subroutine calculate_lattice


!----------------------------------------------------------------
!+
!  Delete particles inside of a defined sphere
!+
!----------------------------------------------------------------
subroutine delete_particles_inside_or_outside_sphere(center,radius,xyzi,hi,revert)
 real, intent(in)              :: center(3),radius,xyzi(3)
 real, intent(inout)           :: hi
 logical, intent(in), optional :: revert
 real                          :: r(3), radius_squared
 logical                       :: revert_local

 if (present(revert)) then
    revert_local = revert
 else
    revert_local = .false.
 endif

 radius_squared = radius**2

 r = xyzi(1:3) - center
 if (revert_local) then
    if (dot_product(r,r) > radius_squared) then
       hi = -abs(hi)
    endif
 else
   if (dot_product(r,r) < radius_squared) then
       hi = -abs(hi)
    endif
 endif

end subroutine delete_particles_inside_or_outside_sphere


!----------------------------------------------------------------
!+
!  Interpolation of the mass transfer rate from the mesa file
!+
!----------------------------------------------------------------
subroutine interpolate_mdot(time,t_arr,mdot_arr,mdoti)
 use table_utils, only: find_nearest_index,interp_1d
 real, intent(in)  :: time,t_arr(:),mdot_arr(:)
 real, intent(out) :: mdoti
 integer :: t1,t2,time_index

 call find_nearest_index(t_arr,time,time_index)
 t1 = time_index
 t2 = time_index + 1
 mdoti = interp_1d(time,t_arr(t1),t_arr(t2),mdot_arr(t1),mdot_arr(t2))

end subroutine interpolate_mdot


!-----------------------------------------------------------------------
!+
!  Print summary of wind properties (assumes inputs are in code units)
!+
!-----------------------------------------------------------------------
subroutine print_summary(vinf,cs_inf,rho_inf,pres_inf,mach_num,pmass,distance_between_layers,&
                         time_between_layers)
 use units, only:unit_velocity,unit_pressure,unit_density
 real, intent(in)    :: vinf,cs_inf,rho_inf,pres_inf,mach_num,pmass,distance_between_layers,time_between_layers

 print*, 'wind speed: ',v_inf * unit_velocity / 1e5," km s^-1"
 print*, 'wind cs: ',cs_inf * unit_velocity / 1e5," km s^-1"
 print*, 'wind density: ',rho_inf * unit_density," g cm^-3"
 print*, 'wind pressure: ',pres_inf * unit_pressure," dyn cm^-2"
 print*, 'wind mach number: ', mach
 print*, 'pmass: ',pmass
 print*, 'distance_between_layers: ',distance_between_layers
 print*, 'time_between_layers: ',time_between_layers

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
 call write_inopt(use_mesa_file,'use_mesa_file','use mesa data file to specify mdot',iunit)
 if (use_mesa_file) then
    call write_inopt(filemesa,'filemesa','mesa file path',iunit)
 else
    call write_inopt(mdot_msun_yr,'mdot','mass transfer rate in solar mass / yr',iunit)
 endif
 call write_inopt(lattice_type,'lattice_type','0: cubic distribution, 1: closepacked distribution',iunit)
 call write_inopt(handled_layers,'handled_layers','(integer) number of handled BHL wind layers',iunit)
 call write_inopt(wind_radius,'wind_radius','radius of the wind cylinder (in code units)',iunit)
 call write_inopt(wind_injection_x,'wind_injection_x','x position of the wind injection boundary (in code units)',iunit)

end subroutine write_options_inject


!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal,error,warning
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
 case('lattice_type')
    read(valstring,*,iostat=ierr) lattice_type
    ngot = ngot + 1
    if (lattice_type/=0 .and. lattice_type/=1)    call fatal(label,'lattice_type must be 0 or 1')
 case('wind_radius')
    read(valstring,*,iostat=ierr) wind_radius
    ngot = ngot + 1
    if (wind_radius <= 0.) call fatal(label,'wind_radius must be >0')
 case('wind_injection_x')
    read(valstring,*,iostat=ierr) wind_injection_x
    ngot = ngot + 1
 case('handled_layers')
    read(valstring,*,iostat=ierr) handled_layers
    ngot = ngot + 1
    if ((mod(handled_layers,2)> 0) .or. (handled_layers==0)) call fatal(label,'handled_layers must be non-zero and even')
 case('mdot')
    read(valstring,*,iostat=ierr) mdot_msun_yr
    ngot = ngot + 1
    if (mdot_msun_yr <= 0.) call fatal(label,'mdot must be positive')
 case('use_mesa_file')
    read(valstring,*,iostat=ierr) use_mesa_file
    ngot = ngot + 1
 case('filemesa')
    read(valstring,*,iostat=ierr) filemesa
    ngot = ngot + 1

 end select
 igotall = (ngot >= 8)

end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
