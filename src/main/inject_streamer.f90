!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Injection module for "streamer" simulations
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - Mdot         : *mass injection rate, in Msun/yr (peak rate if imdot_func > 0)*
!   - mdot_func    : *functional form of dM/dt(t) (0=const)*
!   - stream_width : *width of injected stream in Rsun*
!   - r_inj        : *distance from CoM stream is injected*
!   - phi0         : *phi0 parameter from the Mendoza+09 streamer*
!   - theta0       : *theta0 parameter from the Mendoza+09 streamer*
!   - r0           : *r0 parameter from the Mendoza+09 streamer*
!   - omega        : *angular velocity of cloud stream originates from (s^-1)*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon, units
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'streamer'

 public :: inject_particles, write_options_inject, read_options_inject
 public :: init_inject, set_default_options_inject, update_injected_par

 real, private :: Mdot = 1e-7
 real, private :: Mdotcode = 0.
 integer, private :: imdot_func = 0
 integer, private :: sym_stream = 0
 real, private :: stream_width = 10.
 real, private :: r_inj = 100.
 real, private :: phi0 = 0.
 real, private :: theta0 = 0.
 real, private :: r0 = 1500.
 real, private :: omega = 1e-11
 real, private :: vr_0 = 1.5
 real, private :: tstart = 0.
 real, private :: tend = -1.0

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
 use part,      only:igas,hfact,massoftype,nptmass,gravity
 use partinject,only:add_or_update_particle
 use physcon,   only:pi,solarr,au,solarm,years
 use units,     only:udist,umass,utime,get_G_code
 use random,    only:ran2
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 real :: mtot, v,u, phi0_rad, theta0_rad
 real :: xyzi(3),vxyz(3),Mdot_now,stream_radius
 real :: h,ymin,zmin,rcyl,rcyl2
 real :: xc,yc,zc,vxc,vyc,vzc,x_si,y_si,z_si
 real :: ex(3),ey(3),ez(3),vt
 real :: mass_to_inject, omega_cu, vr_0_cu
 real :: rrand, theta, dx_loc, dz_loc
 real :: rc, vk0, mu, nu, eps, ecc, xi0, G_code, end_time
 integer :: ninject_target, ninjected, ipart, iseed, nstreams

if (tend < 0.) end_time = huge(time)
if (time < tstart .or. time > end_time) return

mtot = 0.0

G_code = get_G_code()

phi0_rad = phi0*pi/180.
theta0_rad = theta0*pi/180.
vr_0_cu = (vr_0*1.0e5) * utime / udist 

omega_cu = omega*utime ! unit is s^-1 in input file, convert to code units
 
 if (gravity) then
  write(*,*), "Disc self-gravity is on. Including disc mass in cloud orbit calculation."
  mtot=sum(xyzmh_ptmass(4,:)) + npartoftype(igas)*massoftype(igas)
else
  mtot=sum(xyzmh_ptmass(4,:))
endif

 stream_radius = stream_width
 ! geometric properties of the injected cylinder
 rcyl = stream_radius   ! radius of injection cylinder
 rcyl2 = rcyl*rcyl
 ymin = -rcyl       ! to ensure flow is centred around injection point
 zmin = -rcyl

 ! work out resolution based on Mdot(t)
 if (imdot_func > 0) then
    Mdot_now = Mdotfunc(time-tstart)
    Mdotcode = Mdot_now*(solarm/umass)/(years/utime)
 endif


 mass_to_inject = Mdotcode * dtlast ! (time - dtlast)
ninject_target = ceiling( mass_to_inject / massoftype(igas) )
h = hfact*rcyl2/ninject_target

ninjected = 0



   call mendoza_state( mtot, r0, omega_cu, theta0_rad, phi0_rad, vr_0_cu, &
                     r_inj, xc,yc,zc, vxc,vyc,vzc )

   vt  = sqrt(vxc*vxc + vyc*vyc + vzc*vzc)
   ex  = (/ vxc, vyc, vzc /)/vt             
   call make_perp_frame(ex, ey, ez)   
   

ipart = npart
iseed = npartoftype(igas)
do while (ninjected < ninject_target)
   u = ran2(iseed)                   
   v = ran2(iseed)
   iseed = iseed - 1
   !rrand  = rcyl * sqrt(u)
   rrand  = rcyl * u
   theta  = 2.0*pi*v
   dx_loc = rrand * cos(theta)
   dz_loc = rrand * sin(theta)

   x_si = xc + dx_loc*ey(1) + dz_loc*ez(1)
   y_si = yc + dx_loc*ey(2) + dz_loc*ez(2)
   z_si = zc + dx_loc*ey(3) + dz_loc*ez(3)

   xyzi = (/ x_si, y_si, z_si /) 
   vxyz = (/ vxc, vyc, vzc /)                

   ninjected = ninjected + 1
   call add_or_update_particle( igas, xyzi, vxyz, h, u, ipart, &
                                npart, npartoftype, xyzh, vxyzu )
   ipart = ipart + 1
   if (sym_stream == 1) then
      xyzi = (/ -x_si, -y_si, z_si /) 
      vxyz = (/ -vxc, -vyc, vzc /)    
      call add_or_update_particle( igas, xyzi, vxyz, h, u, ipart, &
                                npart, npartoftype, xyzh, vxyzu )
      ipart = ipart + 1
   elseif (sym_stream == 2) then
      xyzi = (/ -x_si, -y_si, -z_si /) 
      vxyz = (/ -vxc, -vyc, -vzc /)    
      call add_or_update_particle( igas, xyzi, vxyz, h, u, ipart, &
                                npart, npartoftype, xyzh, vxyzu )
      ipart = ipart + 1
   elseif (sym_stream == 3) then
      xyzi = (/ x_si, y_si, -z_si /) 
      vxyz = (/ vxc, vyc, -vzc /)    
      call add_or_update_particle( igas, xyzi, vxyz, h, u, ipart, &
                                npart, npartoftype, xyzh, vxyzu )
      ipart = ipart + 1
   end if
end do

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
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(imdot_func,'mdot_func','functional form of dM/dt(t) (0=const)',iunit)
 call write_inopt(omega,'omega','angular velocity of cloud stream originates from (s^-1)',iunit)
 call write_inopt(r0, 'r0', 'r0 parameter from the Mendoza+09 streamer',iunit)
 call write_inopt(phi0, 'phi0', 'phi0 parameter from the Mendoza+09 streamer',iunit)
 call write_inopt(theta0, 'theta0', 'theta0 parameter from the Mendoza+09 streamer',iunit)
 call write_inopt(r_inj,'r_inj','distance from CoM stream is injected',iunit)
 call write_inopt(vr_0,'vr_0','radial velocity of cloud stream originates from (km/s)',iunit)
 call write_inopt(Mdot,'Mdot','mass injection rate, in Msun/yr (peak rate if imdot_func > 0)',iunit)
 call write_inopt(stream_width,'stream_width','width of injected stream in au',iunit)
 call write_inopt(sym_stream,'sym_stream','balance angular momentum (0=no, 1=Lz, 2=Lx,Ly)',iunit)
 call write_inopt(tstart,'tstart','start time of injection (in years)',iunit)
 call write_inopt(tend,'tend','end time of injection (negative for inf, in years)',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal,error
 use physcon, only:solarm,years,pi
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
 case('mdot_func')
    read(valstring,*,iostat=ierr) imdot_func
    ngot = ngot + 1
    if (imdot_func <  0) call fatal(label,'imdot_func < 0 in input options')
 case('stream_width')
    read(valstring,*,iostat=ierr) stream_width
    ngot = ngot + 1
    if (stream_width <= 0.) call fatal(label,'stream_width < 0 in input options')
 case('r_inj')
    read(valstring,*,iostat=ierr) r_inj
    ngot = ngot + 1
    if (r_inj <= 0.) call fatal(label,'r_inj < 0 in input options')
 case('phi0')
    read(valstring,*,iostat=ierr) phi0
    ngot = ngot + 1
 case('theta0')
    read(valstring,*,iostat=ierr) theta0
    ngot = ngot + 1
 case('r0')
    read(valstring,*,iostat=ierr) r0
    ngot = ngot + 1
    if (r0 <= 0.) call fatal(label,'r0 < 0 in input options')
 case('omega')
    read(valstring,*,iostat=ierr) omega
    ngot = ngot + 1
    if (omega <= 0.) call fatal(label,'omega < 0 in input options')
 case('vr_0')
    read(valstring,*,iostat=ierr) vr_0
    ngot = ngot + 1
    if (vr_0 <= 0.) call fatal(label,'vr_0 < 0 in input options')
 case('sym_stream')
    read(valstring,*,iostat=ierr) sym_stream
    ngot = ngot + 1
    if (sym_stream < 0 .or. sym_stream > 2) call fatal(label,'sym_stream < 0 or > 2 in input options')
 case('tstart')
    read(valstring,*,iostat=ierr) tstart
    ngot = ngot + 1
 case('tend')
    read(valstring,*,iostat=ierr) tend
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 igotall = (ngot >= 8)

end subroutine read_options_inject


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

   call mendonza_invariant_parameters(mstar, r0, omega, theta0m, vr_0,                 &
                         rc, vk0, mu, nu, eps, ecc, xi0)

   ! --- θ at requested r ---------------------------------------------
   !call theta_at_r(r_inj/rc, theta0m, ecc, xi0, theta)
   theta = theta0m
   phi = phi0m + acos( tan(theta0m) / tan(theta) )   ! Ulrich 1976, eq. (15)


   ! --- velocities in spherical coords -------------------------------
   r_rc = r_inj/rc
   xi   = acos(cos(theta)/cos(theta0m)) + xi0
   vr   = -ecc*sin(theta0m)*sin(xi)/(r_rc*(1.0 - ecc*cos(xi))) * vk0
   vt   =  sin(theta0m)/(sin(theta)*r_rc) *                      &
           sqrt(cos(theta0m)**2 - cos(theta)**2) * vk0
   vp   =  sin(theta0m)**2 /(sin(theta)*r_rc) * vk0


   ! --- spherical → Cartesian (disk frame) ---------------------------
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
   xi0 = acos( (1 - mu*sin(theta0m)**2) / ecc ) ! Assuming purely radial motion from the sphere

end subroutine mendonza_invariant_parameters

!-----------------------------------------------------------------------
subroutine theta_at_r(r_rc, theta0, ecc, xi0, theta)
! 1-D bisection: solve r/rc = sin(theta_0)^2 /(1-e cos(xi))
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
      end if
      if (abs(b-a) < 1.0e-12) exit
   end do
   theta = 0.5*(a+b)
end subroutine theta_at_r

subroutine make_perp_frame(a, b, c)
  real, intent(in)  :: a(3)
  real, intent(out) :: b(3), c(3)
  if (abs(a(1)) < 0.9d0) then
     b = (/ 0.d0, -a(3),  a(2) /)
  else
     b = (/ -a(2), a(1),  0.d0 /)
  end if
  b = b / sqrt(sum(b*b))
  c = (/ a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1) /)
end subroutine make_perp_frame

subroutine set_default_options_inject(flag)

 integer, optional, intent(in) :: flag
end subroutine set_default_options_inject

end module inject
