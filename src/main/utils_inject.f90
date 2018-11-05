!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: injectutils
!
!  DESCRIPTION:
!   Utility routines for geodesic sphere injection
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: icosahedron, part, partinject
!+
!--------------------------------------------------------------------------
module injectutils
 implicit none
 real, parameter :: phi = (sqrt(5.)+1.)/2. ! Golden ratio
 real, parameter :: pi = 4.*atan(1.)

contains

!-----------------------------------------------------------------------
!+
!  solve the cubic equation for q, the integer resolution of the spheres
!  given a desired ratio between the particle spacing in the sphere
!  and the radial spacing of spheres
!+
!-----------------------------------------------------------------------
integer function get_sphere_resolution(dr_on_dp,mV_on_MdotR)
 real, intent(in) :: dr_on_dp,mV_on_MdotR
 real :: fac,term,f1,f2,res
 !
 ! we solve (20*q^3 - 30*q^2 + 16*q - 3) - 1./fac = 0
 ! the following is the only real solution,
 ! extracted from Mathematica
 !
 fac = mV_on_MdotR*sqrt(2.*(5.+sqrt(5.)))
 term = dr_on_dp/fac
 f1 = 15.**(1./3.)
 f2 = (45.*term + sqrt(15. + 2025.*term**2))**(1./3.)
 res = 0.5*(1. - 1./(f1*f2) + f2/f1**2)
 get_sphere_resolution = nint(res) !max(nint(res(1)),4)

end function get_sphere_resolution

!-----------------------------------------------------------------------
!+
!  return number of particles per sphere given q, the integer
!  resolution of the spheres
!+
!-----------------------------------------------------------------------
integer function get_parts_per_sphere(ires)
 integer, intent(in) :: ires

 get_parts_per_sphere = 20 * (2*ires*(ires-1)) + 12

end function get_parts_per_sphere

!-----------------------------------------------------------------------
!+
!  return distance between neighbours on the sphere in units of
!  the sphere radius, given integer resolution of spheres
!+
!-----------------------------------------------------------------------
real function get_neighb_distance(ires)
 integer, intent(in) :: ires

 !unitless: relative to the sphere radius
 get_neighb_distance = 2./((2.*ires-1.)*sqrt(sqrt(5.)*phi))

end function get_neighb_distance

!-----------------------------------------------------------------------
!+
!  Inject a quasi-spherical distribution of particles.
!+
!-----------------------------------------------------------------------
subroutine inject_geodesic_sphere(sphere_number, first_particle, ires, r, v, u, rho, &
           geodesic_R, geodesic_v, npart, npartoftype, xyzh, vxyzu, itype, x0, v0)
 use icosahedron, only:pixel2vector
 use partinject,  only:add_or_update_particle
 use part,        only:hrho
 integer, intent(in) :: sphere_number, first_particle, ires, itype
 real,    intent(in) :: r, v, u, rho, geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real,    intent(in), optional :: x0(3), v0(3)
 integer, intent(inout) :: npart, npartoftype(:)
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)

 real :: rotation_angles(3), h_sim
 real :: rotmat(3,3), radial_unit_vector(3), radial_unit_vector_rotated(3)
 real :: particle_position(3), particle_velocity(3)
 integer :: j, particles_per_sphere

 select case (ires)
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

 ! Quantities in simulation units
 h_sim = hrho(rho)
 particles_per_sphere = get_parts_per_sphere(ires)

 call make_rotation_matrix(rotation_angles, rotmat)
 do j=0,particles_per_sphere-1
    call pixel2vector(j, ires, geodesic_R, geodesic_v, radial_unit_vector)
    radial_unit_vector_rotated(1) = radial_unit_vector(1)*rotmat(1,1) &
                                  + radial_unit_vector(2)*rotmat(1,2) &
                                  + radial_unit_vector(3)*rotmat(1,3)
    radial_unit_vector_rotated(2) = radial_unit_vector(1)*rotmat(2,1) &
                                  + radial_unit_vector(2)*rotmat(2,2) &
                                  + radial_unit_vector(3)*rotmat(2,3)
    radial_unit_vector_rotated(3) = radial_unit_vector(1)*rotmat(3,1) &
                                  + radial_unit_vector(2)*rotmat(3,2) &
                                  + radial_unit_vector(3)*rotmat(3,3)
    particle_position = r*radial_unit_vector_rotated
    particle_velocity = v*radial_unit_vector_rotated
    if (present(x0)) particle_position = particle_position + x0
    if (present(v0)) particle_velocity = particle_velocity + v0
    call add_or_update_particle(itype,particle_position,particle_velocity, &
         h_sim,u,first_particle+j,npart,npartoftype,xyzh,vxyzu)
 enddo

end subroutine inject_geodesic_sphere

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

end subroutine make_rotation_matrix

end module injectutils
