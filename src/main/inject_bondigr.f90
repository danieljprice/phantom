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
!    gammawind       -- adiabatic index gamma
!    iwindres        -- resolution of the wind
!    drwindsph       -- distance between spheres / neighbours
!    isphereshandled -- number of handled inner spheres of the wind (integer)
!    rin             -- radius of injection of the wind
!
!  DEPENDENCIES: eos, icosahedron, infile_utils, io, part, partinject,
!    physcon, timestep, metric
!+
!--------------------------------------------------------------------------
module inject
 use bondiexact, only:get_bondi_solution,isol
 use metric,     only:mass1
 implicit none
 character(len=*), parameter, public :: inject_type = 'wind'

 public :: inject_particles,write_options_inject,read_options_inject,wind_init

!
!--runtime settings for this module
!

!--- Read from input file--------------------------
 integer, private :: iwindres        = 4
 real,    private :: drwindsph       = 4.
 integer, private :: isphereshandled = 4
 real,    public  :: rin             = 5.    ! Injection radius (in units of central mass M)
 real,    public  :: gammawind       = 5./3.

! Calculated from the previous parameters
 real,    public :: masspart

 private

 logical, parameter :: wind_verbose = .false.
 logical, parameter :: inflow       = .false.

 real    :: geodesic_R(0:19,3,3),geodesic_v(0:11,3),uthermontemp
 real    :: rhoin,vin,uthermin,mdot
 real    :: speed
 real    :: masssphere,dtsphere,neighdist
 integer :: npsphere

contains

! Wrapper
subroutine get_solution(rho,v,u,r)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r
 call get_bondi_solution(rho,v,u,r,mass1,gammawind)
 ! Direction of wind
 if (inflow) v = -v
end subroutine get_solution

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

 rin   = rin*mass1
 call get_solution(rhoin,vin,uthermin,rin)
 speed = abs(vin)

 uthermontemp = Rg/(gmw*(gammawind-1.))
 npsphere     = 20 * (2*iwindres*(iwindres-1)) + 12
 neighdist    = 2./((2.*iwindres-1.)*sqrt(sqrt(5.)*phi))

 mdot       = 4.*pi*rin**2*rhoin*speed
 masspart   = drwindsph * neighdist * rin * mdot / (npsphere * speed)
 masssphere = masspart * npsphere
 dtsphere   = masssphere / mdot

 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)

 ! adjusting dtmax to avoid uterm < 0 errors
 if (setup) then
    irrational_numbre_close_to_one = pi/3.
    dtmax = 0.5 * irrational_numbre_close_to_one * dtsphere
 endif

 print*,'========= GR Bondi Wind Injection =========='
 print*,'Particles per sphere :',npsphere
 print*,'Nieghbour distance   :',neighdist
 print*,'Mass of particles    :',masspart
 print*,'Mass of spheres      :',masssphere
 print*,'time between spheres ;',dtsphere
 print*,'============================================'

end subroutine

!-----------------------------------------------------------------------
!+
!  Main routine handling wind injection.
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use io,  only: iprint
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 logical, save :: first_run = .true.
 logical :: outerbound
 integer :: isphereouter,isphereinner,ihandledinner,i
 real :: tlocal,r,v,u,rho,mass_lost
 real :: r2,outer_boundary

 if (first_run) then
    call wind_init(.false.)
    first_run = .false.
 endif

 isphereouter  = floor((time-dtlast)/dtsphere) + 1
 isphereinner  = floor(time/dtsphere)
 ihandledinner = isphereinner + isphereshandled
 if (wind_verbose) then
    write(iprint,*) '   time between spheres',dtsphere
    write(iprint,*) '   Inner sphere  :'     ,isphereinner ,'Outer sphere:',isphereouter
    write(iprint,*) '   Loop goes from:'     ,ihandledinner,'to:          ',isphereouter
 endif
 do i=ihandledinner,isphereouter,-1
    tlocal = time - (i-isphereshandled) * dtsphere
    call compute_sphere_properties(tlocal,r,v,u,rho)
    if (wind_verbose) then
       write(iprint,*) '* Sphere:'
       if (i > isphereinner) then
          write(iprint,*) 'HANDLED'
       else
          write(iprint,*) 'INJECTED'
       endif
       write(iprint,*) '   Number i      : ',i
       write(iprint,*) '   Local Time    : ',tlocal
       write(iprint,*) '   Radius        : ',r
       write(iprint,*) '   Velocity      : ',v
       write(iprint,*) '   Temperature   : ',u / uthermontemp,' K'
       write(iprint,*) '   Density (rho) : ',rho
    endif
    if (i > isphereinner) then
       call inject_geodesic_sphere(i,(ihandledinner-i)*npsphere+1,r,v,u,rho,npart,npartoftype,xyzh,vxyzu) ! handled sphere
    else
       call inject_geodesic_sphere(i,npart+1,r,v,u,rho,npart,npartoftype,xyzh,vxyzu) ! injected sphere
    endif
 enddo
 mass_lost = masssphere * (isphereinner-isphereouter+1)

 !--- "Accrete" particles after reaching some outer boundary
 outerbound = .false.
 if (outerbound) then
    outer_boundary = 20.*mass1
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,outer_boundary) &
    !$omp private(i,r2)
    do i=1,npart
       r2 = xyzh(1,i)**2 + xyzh(2,i)**2 + xyzh(3,i)**2
       if (r2 >outer_boundary**2) xyzh(4,i) = -abs(xyzh(4,i))
    enddo
    !$omp end parallel do
 endif

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Compute the radius, velocity and temperature of a sphere in function of its local time
!+
!-----------------------------------------------------------------------
subroutine compute_sphere_properties(tlocal,r,v,u,rho)
 real, intent(in)   :: tlocal
 real, intent(out)  :: r,v,u,rho
 integer, parameter :: N = 10000
 integer :: i
 real    :: dt,v1,v2,v3,v4

 dt = tlocal / N
 r = rin
 v = vin

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
subroutine inject_geodesic_sphere(isphere,ifirst,r,v,u,rho,npart,npartoftype,xyzh,vxyzu)
 use icosahedron, only: pixel2vector
 use partinject,  only: add_or_update_particle
 use part,        only: igas, hrho
 integer, intent(in)    :: isphere,ifirst
 real,    intent(in)    :: r,v,u,rho
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: rotation_angles(3),h
 real :: runi(3),rotmat(3,3),runi_rotated(3)
 real :: xyzi(3),veli(3)
 integer :: i

 select case (iwindres)
 case(1)
    rotation_angles = (/ 1.28693610288783,  2.97863087745917,  1.03952835451832 /)
 case(2)
    rotation_angles = (/ 1.22718722289660,  2.58239466067315,  1.05360422660344 /)
 case(3)
    rotation_angles = (/ 0.235711384317490, 3.10477287368657,  2.20440220924383 /)
 case(4)
    rotation_angles = (/ 3.05231445647236,  0.397072776282339, 2.27500616856518 /)
 case(5)
    rotation_angles = (/ 0.137429597545199, 1.99860670500403,  1.71609391574493 /)
 case(6)
    rotation_angles = (/ 2.90443293496604,  1.77939686318657,  1.04113050588920 /)
 case(10)
    rotation_angles = (/ 2.40913070927068,  1.91721010369865,  0.899557511636617 /)
 case(15)
    rotation_angles = (/ 1.95605828396746,  0.110825898718538, 1.91174856362170 /)
 case default
    rotation_angles = (/ 1.28693610288783,  2.97863087745917,  1.03952835451832 /)
 end select
 rotation_angles = rotation_angles * isphere

 h = hrho(rho)

 call make_rotation_matrix(rotation_angles, rotmat)
 do i=0,npsphere-1
    call pixel2vector(i, iwindres, geodesic_R, geodesic_v, runi)
    runi_rotated(1) = runi(1)*rotmat(1,1) + runi(2)*rotmat(1,2) + runi(3)*rotmat(1,3)
    runi_rotated(2) = runi(1)*rotmat(2,1) + runi(2)*rotmat(2,2) + runi(3)*rotmat(2,3)
    runi_rotated(3) = runi(1)*rotmat(3,1) + runi(2)*rotmat(3,2) + runi(3)*rotmat(3,3)
    runi = runi_rotated
    xyzi = r*runi
    veli = v*runi
    call add_or_update_particle(igas,xyzi,veli,h,u,ifirst+i,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine

!-----------------------------------------------------------------------
!+
!  Make a 3x3 rotation matrix from three angles.
!+
!-----------------------------------------------------------------------
subroutine make_rotation_matrix(rotation_angles,rot_m)
 real, intent(in)  :: rotation_angles(3)
 real, intent(out) :: rot_m(3,3)

 real :: angle_x,angle_y,angle_z
 real :: c_x,s_x,c_y,s_y,c_z,s_z

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

 call write_inopt(gammawind,'gammawind','polytropic index of the wind',iunit)
 call write_inopt(iwindres,'iwindres',&
      'resolution of the wind (1-6,10,15)-- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(drwindsph,'drwindsph','distance between spheres / neighbours -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(isphereshandled,'isphereshandled','handle inner spheres of the wind (integer)',iunit)
 call write_inopt(rin,'rin', &
      'radius of injection of the wind (in units of the central mass) -- DO NOT CHANGE DURING SIMULATION --',iunit)
 call write_inopt(isol,'isol','solution type (1 = geodesic flow  |  2 = sonic point flow)',iunit)
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
 case('gammawind')
    read(valstring,*,iostat=ierr) gammawind
    ngot = ngot + 1
    if (gammawind < 0.)       call fatal(label,'invalid setting for gammawind (<0)')
 case('iwindres')
    read(valstring,*,iostat=ierr) iwindres
    ngot = ngot + 1
    if (iwindres < 1)  call fatal(label,'iwindres must be bigger than zero')
 case('drwindsph')
    read(valstring,*,iostat=ierr) drwindsph
    ngot = ngot + 1
    if (drwindsph <= 0.)    call fatal(label,'drwindsph must be >=0')
 case('isphereshandled')
    read(valstring,*,iostat=ierr) isphereshandled
    ngot = ngot + 1
    if (isphereshandled <= 0) call fatal(label,'isphereshandled must be > 0')
 case('rin')
    read(valstring,*,iostat=ierr) rin
    ngot = ngot + 1
    if (rin < 2.*mass1) call fatal(label,'invalid setting for rin (<2M)')
 case('isol')
    read(valstring,*,iostat=ierr) isol
    ngot = ngot + 1
    if (isol /= 1 .and. isol /= 2) call fatal(label,'invalid setting for isol')
 case default
    imatch = .false.
 end select

 noptions = 6

 igotall = (ngot >= noptions)

end subroutine

end module inject
