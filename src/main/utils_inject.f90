!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module injectutils
!
! Utility routines for geodesic sphere injection
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: icosahedron, part, partinject, geometry, units, physcon, vector_utils
!
 implicit none
 real, parameter :: phi = (sqrt(5.)+1.)/2. ! Golden ratio
 real, parameter :: pi = 4.*atan(1.)

 logical, public :: use_icosahedron = .false.

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

 if (use_icosahedron) then
    get_parts_per_sphere = 20 * (2*ires*(ires-1)) + 12
 else
    get_parts_per_sphere = ires
 endif

end function get_parts_per_sphere

!-----------------------------------------------------------------------
!+
!  return distance between neighbours on the sphere in units of
!  the sphere radius, given integer resolution of spheres
!+
!-----------------------------------------------------------------------
real function get_neighb_distance(ires)
 integer, intent(in) :: ires

 if (use_icosahedron) then
   ! unitless: relative to the sphere radius
    get_neighb_distance = 2./((2.*ires-1.)*sqrt(sqrt(5.)*phi))
 else
    ! working up to a certain degree, a more precise formulation could exist,
    ! especially when the injection routine only use trigonometric functions

    ! unitless: relative to the sphere radius
    get_neighb_distance = sqrt(4*pi/ires)
 endif

end function get_neighb_distance

!-----------------------------------------------------------------------
!+
!  Inject a quasi-spherical distribution of particles.
!+
!-----------------------------------------------------------------------
subroutine inject_geodesic_sphere(sphere_number, first_particle, ires, r, v, u, rho, &
             geodesic_R, geodesic_v, npart, npartoftype, xyzh, vxyzu, itype, x0, v0, &
             isink,JKmuS)
 use icosahedron, only:pixel2vector,fibonacci_sphere
 use partinject,  only:add_or_update_particle
 use part,        only:hrho,xyzmh_ptmass, iReff, ispinx, ispiny, ispinz, imloss
 use geometry,    only:vector_transform, coord_transform
 use units,       only:unit_velocity
 use physcon,     only:gg,au,solarm,km
 use vectorutils, only:cross_product3D
 integer, intent(in) :: sphere_number, first_particle, ires, itype, isink
 real,    intent(in) :: r,v,u,rho,geodesic_R(0:19,3,3),geodesic_v(0:11,3),x0(3),v0(3)
 real,    intent(in), optional :: JKmuS(:)
 integer, intent(inout) :: npart, npartoftype(:)
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)

 real :: rotation_angles(3), h_sim
 real :: rotmat(3,3), radial_unit_vector(3), radial_unit_vector_rotated(3)
 real :: radial_unit_vector_rotated_new(3), z_axis(3), omega_axis(3)
 real :: particle_position(3), particle_velocity(3), particle_position_proj(3)
 real :: velocity_out(3), position_out(3), vomega_spin(3)
 real :: rot_particle_position(3), rot_particle_velocity(3)
 integer :: j, particles_per_sphere

 real :: omega, rotation_speed_crit, wind_rotation_speed
 integer, parameter :: ndim = 3, igeom = 3

 ! change the method of injection if more than 1 sink emits a wind
 if (xyzmh_ptmass(imloss,2) /= 0 .or. xyzmh_ptmass(imloss,3) /= 0) then
   use_icosahedron = .false.
 endif

 ! Quantities in simulation units
 h_sim = hrho(rho)
 particles_per_sphere = get_parts_per_sphere(ires)

 omega = 0.
 wind_rotation_speed = sqrt(sum(xyzmh_ptmass(ispinx:ispinz,1)**2))/xyzmh_ptmass(iReff,1)**2
 if (wind_rotation_speed > 0.001*km/unit_velocity) then
    rotation_speed_crit = sqrt((gg*xyzmh_ptmass(4,1)*solarm)/(xyzmh_ptmass(iReff,1)*au))/unit_velocity
    omega = wind_rotation_speed/rotation_speed_crit

    ! if the rotation axis is not the z-axis, need to rotate the coordinate system
    omega_axis = xyzmh_ptmass(ispinx:ispinz,1)/(wind_rotation_speed*xyzmh_ptmass(iReff,1)**2)
    ! make sure it is a unit vector
    omega_axis = merge(omega_axis, 0.0, abs(omega_axis) > 1e-5)
 endif

 call optimal_rot_angles(ires,rotation_angles)
 rotation_angles = rotation_angles * sphere_number
 call make_rotation_matrix(rotation_angles, rotmat)

 do j=0,particles_per_sphere-1
    if (use_icosahedron) then
       call pixel2vector(j, ires, geodesic_R, geodesic_v, radial_unit_vector)
    else
       call fibonacci_sphere(j,particles_per_sphere,radial_unit_vector)
    endif

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

    ! get the projected position of particles to compute the centrifugal force
    particle_position_proj = xyzmh_ptmass(iReff,1)*radial_unit_vector_rotated

    if (omega > 0.) then

       if (omega_axis(3) /= 1.) then
          z_axis = [0., 0., 1.0]
          ! rotate the frame to make the omega_axis coincide with z_axis (needed to switch to spherical)
          call rotation_frame(radial_unit_vector_rotated,omega_axis,z_axis,radial_unit_vector_rotated_new)

          rot_particle_position = r*radial_unit_vector_rotated_new
          rot_particle_velocity = v*radial_unit_vector_rotated_new
       else
          rot_particle_position = particle_position
          rot_particle_velocity = particle_velocity
       endif

       ! include velocity component due to spin of the sink particle (Centrifugal force)
       call cross_product3D(xyzmh_ptmass(ispinx:ispinz,1),particle_position_proj,vomega_spin)
       vomega_spin =  vomega_spin/xyzmh_ptmass(iReff,1)**3

       ! geometry subroutines to switch from Cartesian to spherical coordinates
       ! and impose Dwarkadas & Owocki (2002) velocity profile on the radial component

       ! input is cartesian, output spherical polar
       call coord_transform(rot_particle_position, ndim, 1, position_out, ndim, igeom)
       call vector_transform(rot_particle_position, rot_particle_velocity, ndim, 1, velocity_out, ndim, igeom)
       velocity_out(1) = velocity_out(1) * sqrt(1. - omega**2 * sin(position_out(3))**2)

       ! input is spherical polars, output cartesian
       call vector_transform(position_out, velocity_out, ndim, igeom, rot_particle_velocity, ndim, 1)

       if (omega_axis(3) /= 1.) then
          ! rotate the frame to make the z_axis coincide with omega_axis
          call rotation_frame(rot_particle_position,z_axis,omega_axis,particle_position)
          call rotation_frame(rot_particle_velocity,z_axis,omega_axis,particle_velocity)
       else
          particle_position = rot_particle_position
          particle_velocity = rot_particle_velocity
       endif
       particle_position = particle_position + x0
       particle_velocity = particle_velocity + v0 + vomega_spin

    else
       particle_position = particle_position + x0
       particle_velocity = particle_velocity + v0
    endif

    call add_or_update_particle(itype,particle_position,particle_velocity, &
         h_sim,u,first_particle+j,npart,npartoftype,xyzh,vxyzu,JKmuS,isink)
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

!-----------------------------------------------------------------------
!+
!  Rotate the reference frame using Rodrigues' rotation formula
!  to make axis_A coincide with axis_B (unit vectors)
!+
!-----------------------------------------------------------------------
subroutine rotation_frame(vector,axis_A,axis_B,rot_vector)
 use vectorutils, only:cross_product3D
 real, intent(in)  :: vector(3), axis_A(3), axis_B(3)
 real, intent(out) :: rot_vector(3)

 real :: rotation_matrix(3,3), cross_vec(3), cos_theta, sin_theta

 call cross_product3D(axis_A,axis_B,cross_vec)

 cos_theta = dot_product(axis_A,axis_B)
 sin_theta = sqrt(dot_product(cross_vec,cross_vec))

 ! Rodrigues' formula in matrix notation
 rotation_matrix(1,1) = cos_theta + (cross_vec(1)**2)*(1-cos_theta)
 rotation_matrix(1,2) = -cross_vec(3)*sin_theta + cross_vec(1)*cross_vec(2)*(1-cos_theta)
 rotation_matrix(1,3) = cross_vec(2)*sin_theta + cross_vec(1)*cross_vec(3)*(1-cos_theta)
 rotation_matrix(2,1) = cross_vec(3)*sin_theta + cross_vec(1)*cross_vec(2)*(1-cos_theta)
 rotation_matrix(2,2) = cos_theta + (cross_vec(2)**2)*(1-cos_theta)
 rotation_matrix(2,3) = -cross_vec(1)*sin_theta + cross_vec(2)*cross_vec(3)*(1-cos_theta)
 rotation_matrix(3,1) = -cross_vec(2)*sin_theta + cross_vec(1)*cross_vec(3)*(1-cos_theta)
 rotation_matrix(3,2) = cross_vec(1)*sin_theta + cross_vec(2)*cross_vec(3)*(1-cos_theta)
 rotation_matrix(3,3) = cos_theta + (cross_vec(3)**2)*(1-cos_theta)

 rot_vector(1) = dot_product(rotation_matrix(1,1:3),vector(1:3))
 rot_vector(2) = dot_product(rotation_matrix(2,1:3),vector(1:3))
 rot_vector(3) = dot_product(rotation_matrix(3,1:3),vector(1:3))

end subroutine rotation_frame

!-----------------------------------------------------------------------
!+
!  Return the optimal rotation angles given the resolution
! (optimal = which gives the most homogeneous distribution
!  of particles for five consecutives rotations)
!+
!-----------------------------------------------------------------------
subroutine optimal_rot_angles(ires,rotation_angles)

 integer, intent(in) :: ires
 real, intent(out)   :: rotation_angles(3)

 real :: particles_per_sphere

 if (use_icosahedron) then
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
 else
    particles_per_sphere = get_parts_per_sphere(ires)
    if     (particles_per_sphere == 32) then
       rotation_angles = (/ 2.904402425695322, 1.914380562631172, 2.905973222022094 /)
    elseif (particles_per_sphere == 64) then
       rotation_angles = (/ 2.984513038361865, 0.769523415130170, 3.037745563096148 /)
    elseif (particles_per_sphere == 128) then
       rotation_angles = (/ 0.436110737995888, 1.727856784714537, 0.717960533010026 /)
    elseif (particles_per_sphere == 256) then
       rotation_angles = (/ 0.442841078702736, 1.492275685215000, 0.442649331104250 /)
    elseif (particles_per_sphere == 512) then
       rotation_angles = (/ 0.444566807089107, 1.492285272594924, 0.444260010931530 /)
    else
       rotation_angles = (/ 0.00000000000000, 0.00000000000000, 0.00000000000000 /)
    endif
 endif

end subroutine optimal_rot_angles

end module injectutils
