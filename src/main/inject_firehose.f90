!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection module for "firehose" simulations, used for numerical experiments
!  in Liptai et al. paper on tidal disruption events
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Mdot         : *mass injection rate, in Msun/yr (peak rate if imdot_func > 0)*
!   - N            : *number of particles per stream width*
!   - mach         : *Mach number of injected stream*
!   - mdot_func    : *functional form of dM/dt(t) (0=const)*
!   - stream_width : *width of injected stream in Rsun*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon, units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'firehose'

 public :: inject_particles, write_options_inject, read_options_inject
 public :: init_inject, set_default_options_inject, update_injected_par

 real, private :: Mdot = 0.
 real, private :: Mdotcode = 0.
 real, private :: mach = 2.
 integer, private :: imdot_func = 0
 integer, private :: N = 8
 real, private :: stream_width = 1.

contains

subroutine init_inject(ierr)
 use units,   only:umass,utime
 use physcon, only:years,solarm
 integer, intent(out) :: ierr

 ierr = 0
!
!--convert mass injection rate to code units
!
 Mdotcode = Mdot*(solarm/umass)/(years/utime)
 print*,' Mdot is ',Mdot,' Msun/yr, which is ',Mdotcode,' in code units',umass,utime

end subroutine init_inject
!-----------------------------------------------------------------------
!+
!  Main routine handling injection at the L1 point.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass, &
           npart,npart_old,npartoftype,dtinject)
 use part,      only:igas,hfact,massoftype,nptmass
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,solarr,au,solarm,years
 use units,     only:udist,umass,utime
 use eos,       only:gamma
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real :: Rp,Rtidal,Rstar,beta,dt_walls
 real :: xyzi(3),x0(3),vxyz(3),cs,Mdot_now,tb,acrit,mbh_on_mstar,stream_radius
 real :: delta_r,delta_x,h,u,vinject,time_between_walls,local_time,ymin,zmin,rcyl,rcyl2
 integer :: i,iy,iz,i_part,part_type,boundary_walls
 integer :: outer_wall, inner_wall, inner_boundary_wall, particles_per_wall

 ! get the location of the injection point
 ! the following are just to base this on parameters used in our TDE simulation
 ! but Rp in the end is just a location along the x axis
 beta = 1.    ! penetration factor
 Rstar = solarr/udist
 mbh_on_mstar = 1e6
 Rtidal = (solarr/udist)*mbh_on_mstar**(1./3.)
 tb = 100. !2.*pi*sqrt(acrit**3)  ! Mbh = 1 in code units
 Rp = (Rtidal/beta)
 acrit = mbh_on_mstar**(1./3.)*Rtidal/2.
 print*,' Rtidal = ',Rtidal,' = ',Rtidal*udist/au,' au'
 print*,' Rp = ',Rp,' acrit = ',acrit,' time (most bound) = ',tb,' code = ',tb*utime/years,' years'
 print*,' solarr = ',Rstar
 print*,'  1 au = ',au/udist
 print*,' Rschwarzschild = ',2.*udist/au,' au, or ',2*udist/solarr,' Rsun'
 stream_radius = stream_width*Rstar
 print*,' stream radius = ',stream_radius,', or ',stream_radius/Rstar,' Rsun'

 ! specify the injection point
 if (nptmass > 0) then
    x0(1:3) = (/105.,-150.,0./)
 else
    x0(1:3) = (/100.*Rstar,0.,0./)
 endif

 ! give the injection velocity in terms of the velocity of a parabolic orbit at Rp
 ! then use this to specify the sound speed from the desired Mach number
 vinject = 0.5*sqrt(2./Rp)
 cs = vinject/mach
 u = cs**2/(gamma*(gamma-1))

 ! geometric properties of the injected cylinder
 rcyl = stream_radius !0.01*Rp  ! radius of injection cylinder
 rcyl2 = rcyl*rcyl
 ymin = -rcyl       ! to ensure flow is centred around injection point
 zmin = -rcyl
 !N = 8
 delta_r = 2.*rcyl/N  ! particle separation in cylinder
 particles_per_wall = int(0.25*pi*N**2)  ! cross section of cylinder

 ! work out resolution based on Mdot(t)
 if (imdot_func > 0) then
!     Mass_to_inject = dm(time) - dm(time-dtlast)
!     particles_to_inject = Mass_to_inject/
!     walls_to_inject = int(Minject/massoftype(igas)) - total_particles_injected
    Mdot_now = Mdotfunc(time)
    dt_walls = (particles_per_wall*massoftype(igas))/Mdot_now
    delta_x = dt_walls*vinject
    print*,' particle mass = ',massoftype(igas)*umass/solarm
    print*,' particles_per_wall = ',particles_per_wall,' mass per wall (Msun) = ',particles_per_wall*massoftype(igas)*umass/solarm
    print*,' dm/dt (code ) = ',Mdotcode,umass,utime
    print*,' dm/dt (Msun/yr) = ',Mdotcode*(umass/solarm)/(utime/years)
    print*,' dt_walls = ',dt_walls,'dr =',delta_r,' dx = ',delta_x, ' dx/dr = ',delta_x/delta_r
    !read*
 else
    delta_x = delta_r
 endif
 boundary_walls = 4
 h = hfact*max(delta_x,delta_r)

 print*,'injecting at R=',x0(1),' with  v = ',vinject,' Mach # = ',vinject/cs
 print*,' dr =',delta_r,' dx = ',delta_x, ' dx/dr = ',delta_x/delta_r
 !read*
!
!--inject material
!
!  deltat = time - dtlast
!  Minject = deltat*Mdotcode
!  ninject = int(Minject/massoftype(igas)) - 1
!  deltatp = deltat/real(ninject - 1)

 time_between_walls = delta_x/vinject
 outer_wall = ceiling((time-dtlast)/time_between_walls)
 inner_wall = ceiling(time/time_between_walls)-1
 inner_boundary_wall = inner_wall+boundary_walls

 print *, "t = ", time
 print *, "dt last = ", dtlast
 print *, "delta t = ", time_between_walls
 print *, "Injecting wall ", inner_wall, " to ", outer_wall
 print *, "Boundary wall ", inner_boundary_wall, " to ", inner_wall-1
 print *, ' v = ', vinject
 print *, '*** ', time, dtlast, time_between_walls, inner_wall, outer_wall

 do i=inner_boundary_wall,outer_wall,-1
    local_time = time - i*time_between_walls
    if (i  >  inner_wall) then
       ! Boundary layer
       i_part = (inner_boundary_wall-i)*particles_per_wall
       part_type = igas
    else
       ! Live layer
       i_part = npart
       part_type = igas
    endif
    print *, '==== ', i, xyzi(1)
    do iy = 1,N
       do iz = 1,N
          xyzi(2) = local_time * vinject
          xyzi(1) = ymin + (iy-.5)*delta_r
          xyzi(3) = zmin + (iz-.5)*delta_r
          ! crop to cylinder
          if (xyzi(1)**2 + xyzi(3)**2 < rcyl2) then
             ! give injection velocity
             vxyz = (/ 0., vinject, 0. /)
             ! add position offset
             xyzi = xyzi + x0
             i_part = i_part + 1
             ! Another brick in the wall
             call add_or_update_particle(part_type, xyzi, vxyz, h, u, i_part, npart, npartoftype, xyzh, vxyzu)
          endif
       enddo
    enddo
 enddo
 dtinject = huge(dtinject) ! no timestep constraint from injection

contains
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
    Mdotfunc = Mdotcode*(t/tb)**(-5./3.)*(1.-(t/tb)**(-4./3.))
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
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(imdot_func,'mdot_func','functional form of dM/dt(t) (0=const)',iunit)
 call write_inopt(Mdot,'Mdot','mass injection rate, in Msun/yr (peak rate if imdot_func > 0)',iunit)
 call write_inopt(mach,'mach','Mach number of injected stream',iunit)
 call write_inopt(stream_width,'stream_width','width of injected stream in Rsun',iunit)
 call write_inopt(N,'N','number of particles per stream width',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal,error
 use physcon, only:solarm,years
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 select case(trim(name))
 case('Mdot')
    read(valstring,*,iostat=ierr) Mdot
    ngot = ngot + 1
    if (Mdot  <  0.) call fatal(label,'Mdot < 0 in input options')
 case('mach')
    read(valstring,*,iostat=ierr) mach
    ngot = ngot + 1
    if (mach <=  0.) call fatal(label,'mach number <= 0 in input options')
 case('mdot_func')
    read(valstring,*,iostat=ierr) imdot_func
    ngot = ngot + 1
    if (imdot_func <  0) call fatal(label,'imdot_func < 0 in input options')
 case('stream_width')
    read(valstring,*,iostat=ierr) stream_width
    ngot = ngot + 1
    if (stream_width <= 0.) call fatal(label,'stream_width < 0 in input options')
 case('N')
    read(valstring,*,iostat=ierr) N
    ngot = ngot + 1
    if (N <= 1) call fatal(label,'N < 1 in input options')

 case default
    imatch = .false.
 end select

 igotall = (ngot >= 5)

end subroutine read_options_inject

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
