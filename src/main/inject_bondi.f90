!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: inject
!
!  DESCRIPTION:
!  Handles wind injection for the generalised Bondi accretion solution in the Schwarzschild metric
!
!  REFERENCES: Hawley, Smarr and Wilson (1984) ApJ 277, 296
!
!  OWNER: David Liptai
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    isol               -- choice of solution (1 = geodesic flow, 2 = sonic point flow)
!    rcrit              -- critical point
!    wind_gamma         -- adiabatic index gamma
!    iwind_resolution   -- resolution of the wind
!    wind_sphdist       -- distance between spheres / neighbours
!    ihandled_spheres   -- number of handled inner spheres of the wind (integer)
!    wind_inject_radius -- radius of injection of the wind
!
!  DEPENDENCIES: eos, icosahedron, infile_utils, io, part, partinject,
!    physcon, timestep, metric
!+
!--------------------------------------------------------------------------
module inject
 use metric, only: mass1
 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: inject_particles, write_options_inject, read_options_inject, wind_init

!
!--runtime settings for this module
!

!--- Read from input file--------------------------
 integer, public  :: isol               = 2
 integer, private :: iwind_resolution   = 4 !4
 real,    private :: wind_sphdist       = 4. !0.1 !10.
 integer, private :: ihandled_spheres   = 4 !7 !   100 !4
 real, public     :: wind_inject_radius = 5.    ! Injection radius (in units of central mass M)
 real, public     :: wind_gamma         = 5./3.
 !-- Choice of constants for geodesic flow solution
 real, public     :: den0  = 1.    !  12.
 real, public     :: en0   = 1.e-9 !  0.000297118
 !-- Choice of constant for sonic point flow solution
 real, public     :: rcrit = 8.     !--The critical point (in units of central mass M)

! Calculated from the previous parameters
 real,    public :: mass_of_particles

 private

 logical, parameter :: wind_verbose = .false.

 real    :: geodesic_R(0:19,3,3), geodesic_v(0:11,3), u_to_temperature_ratio
 real    :: C1,C2,Tc,n
 real    :: wind_injection_rho, wind_injection_velocity, wind_injection_utherm, wind_mass_rate
 real    :: mass_of_spheres, time_between_spheres, neighbour_distance
 integer :: particles_per_sphere

contains

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables
!+
!-----------------------------------------------------------------------
subroutine wind_init(setup)
 use physcon,     only: Rg, pi
 use icosahedron, only: compute_matrices, compute_corners
 use timestep,    only: dtmax
 use eos,         only: gmw
 logical, intent(in) :: setup
 real, parameter :: phi = (sqrt(5.)+1.)/2. ! Golden ratio
 real :: irrational_numbre_close_to_one

 if (isol==2) call compute_constants(rcrit)

 call get_solution(wind_injection_rho,wind_injection_velocity,wind_injection_utherm,wind_inject_radius)

 u_to_temperature_ratio = Rg/(gmw*(wind_gamma-1.))
 particles_per_sphere   = 20 * (2*iwind_resolution*(iwind_resolution-1)) + 12
 neighbour_distance     = 2./((2.*iwind_resolution-1.)*sqrt(sqrt(5.)*phi))

 wind_mass_rate         = 4.*pi*wind_inject_radius**2*wind_injection_rho*wind_injection_velocity
 mass_of_particles      = wind_sphdist * neighbour_distance * wind_inject_radius * wind_mass_rate &
                           / (particles_per_sphere * wind_injection_velocity)
 mass_of_spheres        = mass_of_particles * particles_per_sphere
 time_between_spheres   = mass_of_spheres / wind_mass_rate

! Dave's alternate attempt
 ! mass_of_particles      = wind_injection_rho*4./3.*pi*neighbour_distance**3 !??? No idea....
 ! mass_of_spheres        = mass_of_particles * particles_per_sphere
 ! time_between_spheres   = wind_inject_radius/wind_injection_velocity  !??? No idea....

 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 ! adjusting dtmax to avoid uterm < 0 errors
 if (setup) then
    irrational_numbre_close_to_one = pi / 3.
    dtmax = (.5 * irrational_numbre_close_to_one * time_between_spheres)
 endif

 print*,'Particles per sphere :',particles_per_sphere
 print*,'Nieghbour distance   :',neighbour_distance
 print*,'Mass of particles    :',mass_of_particles
 print*,'Mass of spheres      :',mass_of_spheres
 print*,'time between spheres ;',time_between_spheres
 print*,'============================================'

end subroutine

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time_u,dtlast_u,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use io,      only: iprint
 real,    intent(in)    :: time_u, dtlast_u
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: outer_sphere, inner_sphere, inner_handled_sphere, i
 real :: time, dtlast, local_time, r, v, u, rho, mass_lost
 logical, save :: first_run = .true.

 if (first_run) then
    call wind_init(.false.)
    first_run = .false.
 endif

 time = time_u
 dtlast = dtlast_u
 outer_sphere = floor((time-dtlast)/time_between_spheres) + 1
 inner_sphere = floor(time/time_between_spheres)
 inner_handled_sphere = inner_sphere + ihandled_spheres
 if (wind_verbose) then
    write(iprint,*) '   time between spheres',time_between_spheres
    write(iprint,*) '   Inner sphere  :', inner_sphere                 ,'Outer sphere:',outer_sphere
    write(iprint,*) '   Loop goes from:', inner_sphere+ihandled_spheres,'to:           ',outer_sphere
endif
 do i=inner_sphere+ihandled_spheres,outer_sphere,-1
    local_time = time - (i-ihandled_spheres) * time_between_spheres
    call compute_sphere_properties(local_time, r, v, u, rho)
    if (wind_verbose) then
       write(iprint,*) '* Sphere:'
       if (i > inner_sphere) then
          write(iprint,*) 'HANDLED'
        else
          write(iprint,*) 'INJECTED'
       endif
       write(iprint,*) '   Number i      : ', i
       write(iprint,*) '   Local Time    : ', local_time
       write(iprint,*) '   Radius        : ', r
       write(iprint,*) '   Velocity      : ', v
       write(iprint,*) '   Temperature   : ', u / u_to_temperature_ratio, ' K'
       write(iprint,*) '   Density (rho) : ', rho
    endif
    if (i > inner_sphere) then
       call inject_geodesic_sphere(i, (ihandled_spheres-i+inner_sphere)*particles_per_sphere+1, r, v, u, rho, &
        npart, npartoftype, xyzh, vxyzu) ! handled sphere
    else
       call inject_geodesic_sphere(i, npart+1, r, v, u, rho, npart, npartoftype, xyzh, vxyzu) ! injected sphere
    endif
 enddo
 mass_lost = mass_of_spheres * (inner_sphere-outer_sphere+1)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Compute the radius, velocity and temperature of a sphere in function of its local time
!+
!-----------------------------------------------------------------------
subroutine compute_sphere_properties(local_time, r, v, u, rho)
 real, intent(in) :: local_time
 real, intent(out) :: r, v, u, rho

 real :: dt
 integer, parameter :: N = 10000
 integer :: i
 real :: v1,v2,v3,v4

 dt = local_time / N
 r = wind_inject_radius
 v = wind_injection_velocity

!--- Note: I don't know if RK4 is necessary here, perhaps simple Euler would do, but speed isn't really affected it seems.
 ! iterations
 do i=1,N
    ! r = r + dt*v
    ! call get_solution(rho,v,u,r)
    call get_solution(rho,v1,u,r)
    call get_solution(rho,v2,u,r+dt/2.*v1)
    call get_solution(rho,v3,u,r+dt/2.*v2)
    call get_solution(rho,v4,u,r+dt*v3)
    r = r + dt/6. * (v1 + 2.*v2 + 2.*v3 + v4)
    call get_solution(rho,v,u,r)
 enddo
end subroutine

!-----------------------------------------------------------------------
!+
!  Inject a quasi-spherical distribution of particles.
!+
!-----------------------------------------------------------------------
subroutine inject_geodesic_sphere(sphere_number, first_particle, r, v, u, rho, npart, npartoftype, xyzh, vxyzu)
 use icosahedron, only: pixel2vector
 use partinject,  only: add_or_update_particle
 use part,        only: igas, hrho
 integer, intent(in) :: sphere_number, first_particle
 real, intent(in) :: r, v, u, rho
 integer, intent(inout) :: npart, npartoftype(:)
 real, intent(inout) :: xyzh(:,:), vxyzu(:,:)

 real :: rotation_angles(3), r_sim, v_sim, u_sim, h_sim
 real :: radial_unit_vector(3), rotmat(3,3), radial_unit_vector_rotated(3)
 real :: particle_position(3), particle_velocity(3)
 integer :: j

 select case (iwind_resolution)
 case(1)
    rotation_angles = (/ 1.28693610288783, 2.97863087745917, 1.03952835451832 /)
 case(2)
    rotation_angles = (/ 1.22718722289660, 2.58239466067315, 1.05360422660344 /)
 case(3)
    rotation_angles = (/ 0.235711384317490, 3.10477287368657, 2.20440220924383 /)
 case(4)
    rotation_angles = (/ 3.05231445647236, 0.397072776282339, 2.27500616856518 /)
 case(5)
    rotation_angles = (/ 0.137429597545199, 1.99860670500403, 1.71609391574493 /)
 case(6)
    rotation_angles = (/ 2.90443293496604, 1.77939686318657, 1.04113050588920 /)
 case(10)
    rotation_angles = (/ 2.40913070927068, 1.91721010369865, 0.899557511636617 /)
 case(15)
    rotation_angles = (/ 1.95605828396746, 0.110825898718538, 1.91174856362170 /)
 case default
    rotation_angles = (/ 1.28693610288783, 2.97863087745917, 1.03952835451832 /)
 end select
 rotation_angles = rotation_angles * sphere_number

 r_sim = r
 v_sim = v
 u_sim = u
 h_sim = hrho(rho)

 call make_rotation_matrix(rotation_angles, rotmat)
 do j=0,particles_per_sphere-1
    call pixel2vector(j, iwind_resolution, geodesic_R, geodesic_v, radial_unit_vector)
    radial_unit_vector_rotated(1) = radial_unit_vector(1)*rotmat(1,1) &
                                  + radial_unit_vector(2)*rotmat(1,2) &
                                  + radial_unit_vector(3)*rotmat(1,3)
    radial_unit_vector_rotated(2) = radial_unit_vector(1)*rotmat(2,1) &
                                  + radial_unit_vector(2)*rotmat(2,2) &
                                  + radial_unit_vector(3)*rotmat(2,3)
    radial_unit_vector_rotated(3) = radial_unit_vector(1)*rotmat(3,1) &
                                  + radial_unit_vector(2)*rotmat(3,2) &
                                  + radial_unit_vector(3)*rotmat(3,3)
    radial_unit_vector = radial_unit_vector_rotated !/ sqrt(dot_product(radial_unit_vector_rotated, radial_unit_vector_rotated))
    particle_position = r_sim*radial_unit_vector
    particle_velocity = v_sim*radial_unit_vector
    call add_or_update_particle(igas, particle_position, particle_velocity, h_sim, u_sim, first_particle+j, &
         npart,npartoftype,xyzh,vxyzu)
 enddo
end subroutine

!-----------------------------------------------------------------------
!+
!  Make a 3x3 rotation matrix from three angles.
!+
!-----------------------------------------------------------------------
subroutine make_rotation_matrix(rotation_angles, rot_m)
 real, intent(in)  :: rotation_angles(3)
 real, intent(out) :: rot_m(3,3)

 real :: angle_x, angle_y, angle_z
 real :: c_x, s_x, c_y, s_y, c_z, s_z

 angle_x = rotation_angles(1)
 angle_y = rotation_angles(2)
 angle_z = rotation_angles(3)

 c_x = cos(angle_x)
 s_x = sin(angle_x)
 c_y = cos(angle_y)
 s_y = sin(angle_y)
 c_z = cos(angle_z)
 s_z = sin(angle_z)

 rot_m(1,1) = c_y*c_z
 rot_m(1,2) = -c_y*s_z
 rot_m(1,3) = -s_y
 rot_m(2,1) = -s_x*s_y*c_z + c_x*s_z
 rot_m(2,2) = s_x*s_y*s_z + c_x*c_z
 rot_m(2,3) = -s_x*c_y
 rot_m(3,1) = c_x*s_y*c_z + s_x*s_z
 rot_m(3,2) = -c_x*s_y*s_z + s_x*c_z
 rot_m(3,3) = c_x*c_y
end subroutine

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 call write_inopt(isol,'isol','choice of solution (1 = geodesic flow, 2 = sonic point flow) --',iunit)
 call write_inopt(rcrit,'rcrit','critical point -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(wind_gamma,'wind_gamma','polytropic indice of the wind',iunit)
 call write_inopt(iwind_resolution,'iwind_resolution','resolution of the wind -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(wind_sphdist,'wind_sphdist','distance between spheres / neighbours -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(ihandled_spheres,'ihandled_spheres','handle inner spheres of the wind (integer)',iunit)
 call write_inopt(wind_inject_radius,'wind_inject_radius', &
      'radius of injection of the wind -- DO NOT CHANGE DURING SIMULATION --',iunit)
end subroutine

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,      only: fatal, error, warning
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 integer :: noptions
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('isol')
    read(valstring,*,iostat=ierr) isol
    ngot = ngot + 1
    if (isol /= 1 .and. isol /= 2) call fatal(label,'invalid setting for isol')
 case('rcrit')
    read(valstring,*,iostat=ierr) rcrit
    ngot = ngot + 1
    if (rcrit < 0.)    call fatal(label,'invalid setting for rcrit (<0)')
 case('wind_gamma')
    read(valstring,*,iostat=ierr) wind_gamma
    ngot = ngot + 1
    if (wind_gamma < 0.)       call fatal(label,'invalid setting for wind_gamma (<0)')
 case('iwind_resolution')
    read(valstring,*,iostat=ierr) iwind_resolution
    ngot = ngot + 1
    if (iwind_resolution < 1)  call fatal(label,'iwind_resolution must be bigger than zero')
 case('wind_sphdist')
    read(valstring,*,iostat=ierr) wind_sphdist
    ngot = ngot + 1
    if (wind_sphdist <= 0.)    call fatal(label,'wind_sphdist must be >=0')
 case('ihandled_spheres')
    read(valstring,*,iostat=ierr) ihandled_spheres
    ngot = ngot + 1
    if (ihandled_spheres <= 0) call fatal(label,'ihandled_spheres must be > 0')
 case('wind_inject_radius')
    read(valstring,*,iostat=ierr) wind_inject_radius
    ngot = ngot + 1
    if (wind_inject_radius < 2.*mass1) call fatal(label,'invalid setting for wind_inject_radius (<2M)')
 case default
    imatch = .false.
 end select

 noptions = 7

 igotall = (ngot >= noptions)

end subroutine

!--- Wrapper for the two different solutions
subroutine get_solution(rho,v,u,r)
   real, intent(out) :: rho,v,u
   real, intent(in)  :: r

   select case(isol)
   case(1)
      call get_solution_geodesic(rho,v,u,r)
   case(2)
      call get_solution_sonicpoint(rho,v,u,r)
   end select

end subroutine get_solution

!------- Geodesic Flow Solution -----------------
subroutine get_solution_geodesic(rho,v,u,r)
   real, intent(out) :: rho,v,u
   real, intent(in)  :: r
   real :: Dr,sqrtg,alpha,Er  !,U0,dens

   Dr    = den0/(r**2*sqrt(2.*mass1/r*(1.- 2.*mass1/r)))
   sqrtg = 1.
   alpha = sqrt(1. - 2.*mass1/r)
   rho   = sqrtg*Dr/alpha
   ! U0    = 1./sqrt(1. - 2.*mass1/r - v**2/(1.-2.*mass1/r))
   ! U0    = 1./(1. - 2.*mass1/r)
   !dens  = rho/(sqrtg*U0)

   Er  = en0/((sqrt(2.*mass1/r)*r**2)**wind_gamma * (1.- 2.*mass1/r)**((wind_gamma + 1.)/4.))
   u   = Er/Dr
   v   = sqrt(2.*mass1/r)*(1. - 2.*mass1/r)
end subroutine get_solution_geodesic

!------- Sonic Point Flow Solution -----------------
subroutine get_solution_sonicpoint(rho,v,u,r)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r
 real, parameter :: adiabat = 1.
 real :: T,uvel,term,u0,sqrtg,dens

 ! Given an r, solve eq 76 for T numerically
 call Tsolve(T,r)

 uvel = C1/(r**2 * T**n)
 dens = adiabat*T**n
 u = T*n

 !get u0 at r
 term = 1./(1.-2.*mass1/r)
 u0  = sqrt(term*(1.+term*uvel**2))
 v   = uvel/u0

 sqrtg = 1. !???? FIX
 rho = sqrtg*u0*dens

end subroutine get_solution_sonicpoint

subroutine compute_constants(rcrit)
 real, intent(in) :: rcrit
 real :: uc2,vc2

 n   = 1./(wind_gamma-1.)
 uc2 = mass1/(2.*rcrit)
 vc2 = uc2/(1.-3.*uc2)
 Tc  = vc2*n/(1.+n-vc2*n*(1.+n))

 C1  = sqrt(uc2) * Tc**n * rcrit**2
 C2  = (1. + (1.+n)*Tc)**2 * (1. - 2.*mass1/rcrit + C1**2/(rcrit**4*Tc**(2.*n)))
 print*,'Constants have been set'
 print*,'C1 = ',C1
 print*,'C2 = ',C2

end subroutine compute_constants

! Newton Raphson
subroutine Tsolve(T,r)
   use io, only: warning
   real, intent(in) :: r
   real, intent(out) :: T
   real :: Tnew, diff
   logical :: converged
   integer :: its
   integer, parameter :: itsmax = 100
   real, parameter :: tol = 1.e-5

   T = Tc !Guess

   converged = .false.
   its = 0
   do while (.not.converged .and. its<itsmax)
      Tnew = T - ffunc(T,r)/df(T,r)
      diff = abs(Tnew - T)/abs(T)
      converged = diff < tol
      T = Tnew
      its = its+1
   enddo

   if (.not. converged) call warning('inject bondi','solution not converged for r =',val=r)

   ! print*,'Found T to be:',T,' with ',its,' iterations'

end subroutine Tsolve

real function ffunc(T,r)
  real, intent(in) :: T,r
  ffunc = (1. + (1. + n)*T)**2*(1. - (2.*mass1)/r + C1**2/(r**4*T**(2.*n))) - C2
end function ffunc

real function df(T,r)
   real, intent(in) :: T,r
   df = (2.*(1. + T + n*T)*((1. + n)*r**3*(-2.*mass1 + r) - C1**2*T**(-1. - 2.*n)*(n + (-1. + n**2)*T)))/r**4
end function df

end module inject
