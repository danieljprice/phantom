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
! :Owner: Ana Lourdes Juarez
!
! :Runtime parameters:
!   - BHL_radius       : *radius of the wind cylinder (in star radii)*
!   - handled_layers   : *(integer) number of handled BHL wind layers*
!   - wind_injection_x : *x position of the wind injection boundary (in star radii)*
!   - wind_length      : *crude wind length (in star radii)*
!
! :Dependencies: dim, eos, infile_utils, io, part, partinject, physcon,
!   setbinary, units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'windtunnel'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
      set_default_options_inject,update_injected_par
!
!--runtime settings for this module
!

 ! Particle-related parameters
 integer, public :: handled_layers = 4
 real,    public :: wind_radius = 30.
 real,    public :: wind_injection_x = -10.
 real,    public :: wind_length = 100.

 real, private :: Mdot = 1.0e-5
 real, private :: gastemp = 3000.

 private
 real    :: wind_rad,wind_x,psep,distance_between_layers,&
            time_between_layers,h_inf,u_inf
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
 use partinject,only:add_or_update_particle
 use setbinary, only:L1_point
 use physcon,   only:pi,twopi,solarm,years,gg,kboltz,mass_proton_cgs
 use units,     only:udist, umass, utime
 use eos,        only:gamma
 use part,       only:hfact,massoftype,igas
 use dim,        only:maxp
 use io,         only:fatal
 integer, intent(out) :: ierr
 real :: pmass,element_volume,y,z,cs,mach
 integer :: size_y, size_z, pass, i, j

 ierr = 0

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npart_old,npartoftype,dtinject)
 use physcon,  only:pi,twopi,solarm,years,gg,kboltz,mass_proton_cgs
 use units,    only:utime,udist,umass
 use partinject,only:add_or_update_particle
 use setbinary, only:L1_point
 use eos,        only:gamma,gmw
 use part,       only:hfact,massoftype,igas
 use dim,        only:maxp
 use io,         only:fatal
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 real :: last_time, local_time, x, irrational_number_close_to_one
 integer :: inner_layer, outer_layer, i, i_limited, i_part, np, ierr
 real, allocatable :: xyz(:,:), vxyz(:,:), h(:), u(:)
 real :: m1,m2,q,radL1,theta_s,A,mu,theta_rand,r_rand,dNdt_code,Porb,r12,r2L1,r0L1,smag
 real :: eps, spd_inject, phizzs, phinns, lm12, lm1, lm32, sw_chi, sw_gamma, XL1, U1, lsutime,rad_inj
 real :: xyzL1(3),xyzi(3),dr(3),x1(3),x2(3),x0(3),dxyz(3),vxyzL1(3),v1(3),v2(3),xyzinj(3),s(3)
 real :: rptmass(3),rptmass_uni(3),pmass,element_volume,y,z,mach
 real :: rho_noz,vol_noz,cs,u_part,pr_noz
 integer :: size_y, size_z, pass,j

 x1 = xyzmh_ptmass(1:3,1)
 x2 = xyzmh_ptmass(1:3,2)
 v1 = vxyz_ptmass(1:3,1)
 v2 = vxyz_ptmass(1:3,2)
 dr = x2 - x1
 r12 = dist(x2,x1)
 m1 = xyzmh_ptmass(4,1)
 m2 = xyzmh_ptmass(4,2)
 x0 = (m1*x1 + m2*x2)/(m1 + m2)
 q  = m2/m1
 mu = 1./(1 + q)
 radL1      = L1_point(1./q)                     ! find L1 point given binary mass ratio
 XL1        = radL1-(1-mu)
 lsutime = ((m1+m2)/(r12**3))**(-1./2)*utime

!
!--quantities related to the gas injection at/near L1
!
 Porb      = twopi * sqrt( (r12*udist)**3 / (gg*(m1+m2)*umass) )
 eps = Porb/(twopi*r12*udist) * (gastemp*kboltz/(gmw*mass_proton_cgs))**0.5
 A  = mu / abs(XL1 - 1. + mu)**3 + (1. - mu)/abs(XL1 + mu)**3! See Lubow & Shu 1975, eq 13
 theta_s = -acos( -4./(3.*A)+(1-8./(9.*A))**0.5)/2.              ! See Lubow & Shu 1975, eq 24
 xyzL1(1:3) = XL1*dr(:)   ! set as vector position
 r0L1 = dist(xyzL1, (/0., 0., 0./))               ! distance from L1 to center of mass
 r2L1 = dist(xyzL1, x2)                           ! ... and from the mass donor's center
 s =  (/1.,0.,0./)*r2L1*eps/200.0 !(/cos(theta_s),sin(theta_s),0.0/)*r2L1*eps/200.0   ! last factor still a "magic number". Fix. ANA:WHAT??
 smag = sqrt(dot_product(s,s))
 xyzinj(1:3) = xyzL1 + s
 vxyzL1 = v1*dist(xyzL1,x0)/dist(x0, x1) ! orbital motion of L1 point
 U1 = -3*A*sin(2*theta_s)/(4*lsutime/utime)
 spd_inject = U1*dist(xyzinj,xyzL1)  !L&S75 eq 23b
 lm12 = (A - 2. + sqrt(A*(9.*A - 8.)))/2.                        ! See Heerlein+99, eq A8
 lm32 = lm12 - A + 2.                                            ! See Heerlein+99, eq A15
 lm1  = sqrt(A)
 sw_chi = (1e8/udist)**(-2)
 sw_gamma = (1e8/udist)**(-2)
 rad_inj = sw_chi**(-0.5)
 vol_noz = (4./3.)*pi*rad_inj**3
 rho_noz = mdot*dtlast/vol_noz
 u_part = 3.*(kboltz*gastemp/(mu*mass_proton_cgs))/2.
 cs = (u_part*gamma*(gamma-1.))
 mach = spd_inject/cs
 pr_noz = rho_noz*cs**2

 if (first_run) then
    call init_inject(ierr)
    ! Calculate particle separation between layers given rho_inf, depending on lattice type
    element_volume = pmass / rho_noz
    psep = element_volume**(1./3.)

    distance_between_layers = psep
    size_y = ceiling(3.*rad_inj/psep)
    size_z = size_y
    do pass=1,2
       if  (pass == 2) allocate(layer_even(2,neven), layer_odd(2,neven))
       neven = 0
       do i=1,size_y
          do j=1,size_z
             y = -1.5*rad_inj+(i-1)*psep
             z = -1.5*rad_inj+(j-1)*psep
             if (y**2+z**2  <  rad_inj**2) then
                neven = neven + 1
                if (pass == 2) layer_even(:,neven) = (/ y,z /)
             endif
          enddo
       enddo
    enddo
    layer_odd(:,:) = layer_even(:,:)

    h_inf = hfact*(pmass/rho_noz)**(1./3.)
    max_layers = int(rad_inj/distance_between_layers)
    max_particles = int(max_layers*(nodd+neven)/2)
    time_between_layers = distance_between_layers/spd_inject

    call print_summary(spd_inject,cs,rho_noz,pr_noz,mach,pmass,distance_between_layers,&
                       time_between_layers,max_layers,max_particles)

    if (max_particles > maxp) call fatal('windtunnel', 'maxp too small for this simulation, please increase MAXP!')

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
    x = XL1 + local_time*spd_inject
    xyz(1,:) = x
    vxyz(1,:) = spd_inject
    vxyz(2:3,:) = 0.
    h(:) = h_inf
    u(:) = 3.*(kboltz*gastemp/(mu*mass_proton_cgs))/2.
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
    call inject_or_update_particles(i_part+1, np, xyz, vxyz, h, u, .false.)
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

subroutine update_injected_par
 ! -- placeholder function
 ! -- does not do anything and will never be used
end subroutine update_injected_par

!
! Function to get the distance between two points
!
real function dist(x1,x2)
 real, intent(in)  :: x1(3), x2(3)
 real :: dr(3)

 dr = x1 - x2
 dist = sqrt(dot_product(dr, dr))
 return
end function dist

!-----------------------------------------------------------------------
!+
!  Print summary of wind properties (assumes inputs are in code units)
!+
!-----------------------------------------------------------------------
subroutine print_summary(spd_inject,cs,rho_noz,pr_noz,mach,pmass,distance_between_layers,&
                         time_between_layers,max_layers,max_particles)
 use units, only:unit_velocity,unit_pressure,unit_density
 real, intent(in)    :: spd_inject,cs,rho_noz,pr_noz,mach,pmass,distance_between_layers,time_between_layers
 integer, intent(in) :: max_layers,max_particles

 print*, 'wind speed: ',spd_inject * unit_velocity / 1e5," km s^-1"
 print*, 'wind cs: ',cs * unit_velocity / 1e5," km s^-1"
 print*, 'wind density: ',rho_noz * unit_density," g cm^-3"
 print*, 'wind pressure: ',pr_noz * unit_pressure," dyn cm^-2"
 print*, 'wind mach number: ', mach

 print*, 'maximum wind layers: ', max_layers
 print*, 'pmass: ',pmass
 print*, 'max. wind particles: ', max_particles
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

 call write_inopt(handled_layers,'handled_layers','(integer) number of handled BHL wind layers',iunit)
 call write_inopt(wind_radius,'BHL_radius','radius of the wind cylinder (in star radii)',iunit)
 call write_inopt(wind_injection_x,'wind_injection_x','x position of the wind injection boundary (in star radii)',iunit)
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
 end select

 igotall = (ngot >= 10)
end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
