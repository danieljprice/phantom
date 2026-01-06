!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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
! :Dependencies: icosahedron, part, partinject, geometry, units,
!                physcon, vector_utils, random
!
 use physcon, only:pi

 implicit none
 real, parameter :: phi = (sqrt(5.)+1.)/2. ! Golden ratio

 logical :: jets = .false.
 integer :: seed_random = -1
 real    :: edge_velocity, opening_angle

contains


!-----------------------------------------------------------------------
!+
!  return mean distance between neighbours on the sphere in units of
!  the sphere radius, given the total number of particles
!+
!-----------------------------------------------------------------------
real function get_neighb_distance(nps)
 integer, intent(in) :: nps

! unitless: relative to the sphere radius
 get_neighb_distance = sqrt(4.*pi/nps)

end function get_neighb_distance

!-----------------------------------------------------------------------
!+
!  Inject a quasi-spherical distribution of particles.
!+
!-----------------------------------------------------------------------
subroutine inject_geodesic_sphere(sphere_number, first_particle, ires, r, v, u, rho, &
             geodesic_R, geodesic_v, npart, npartoftype, xyzh, vxyzu, itype, x0, v0, &
             isink,JKmuS,rstar,mstar,omega_vec,vwind_terminal)
 use icosahedron, only:pixel2vector,fibonacci_sphere,fibonacci_jets
 use partinject,  only:add_or_update_particle
 use part,        only:hrho
 use units,       only:unit_velocity
 use physcon,     only:km
 integer, intent(in) :: sphere_number,first_particle,ires,itype,isink
 real,    intent(in) :: r,v,u,rho,geodesic_R(0:19,3,3),geodesic_v(0:11,3),x0(3),v0(3)
 real,    intent(in), optional :: rstar,mstar,omega_vec(3),vwind_terminal
 real,    intent(in), optional :: JKmuS(:)
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 integer, intent(inout) :: npart, npartoftype(:)

 real :: omega,rotation_speed_crit,wind_rotation_speed,h_sim
 real :: radial_unit_vector(3),radial_unit_vector_rotated(3),omega_axis(3)
 real :: particle_position(3),particle_velocity(3),rotation_angles(3),rotmat(3,3)
 integer :: j,particles_per_sphere

 ! Quantities in simulation units
 h_sim = hrho(rho)
 particles_per_sphere = ires

 ! check if the wind emitting sink is rotating
 wind_rotation_speed = 0.
 if (present(rstar) .and. present(omega_vec)) then
    if (rstar > 0.) wind_rotation_speed = sqrt(sum(omega_vec**2))/rstar**2
 endif

 omega = 0.
 if (wind_rotation_speed > 0.001*km/unit_velocity .and. present(rstar) .and. present(mstar)) then
    rotation_speed_crit = sqrt(mstar/rstar)
    omega = wind_rotation_speed/rotation_speed_crit
    ! if the rotation axis is not the z-axis, need to rotate the coordinate system
    omega_axis = omega_vec/(wind_rotation_speed*rstar**2)
    ! make sure the rotation axis is a unit vector
    omega_axis = merge(omega_axis, 0., abs(omega_axis) > 1e-5)
 endif

 if (jets .and. isink == 1) then
    call optimal_rot_angles(ires,rotation_angles,.true.)
 else
    call optimal_rot_angles(ires,rotation_angles,.false.)
 endif

 rotation_angles = rotation_angles * sphere_number
 call make_rotation_matrix(rotation_angles, rotmat)

 do j=0,particles_per_sphere-1
    if (jets .and. isink == 1) then
       call fibonacci_jets(j,particles_per_sphere,radial_unit_vector)
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

    if (jets .and. isink == 1 .and. present(vwind_terminal)) then
       call velocity_jets(radial_unit_vector_rotated,r,v,edge_velocity,opening_angle,&
                          particle_position,particle_velocity,vwind_terminal)
    elseif (omega > 0.01 .and. present(rstar) .and. present(omega_vec)) then
       call velocity_rotating_sink(radial_unit_vector_rotated,r,v,omega,omega_axis,&
                                   particle_position,particle_velocity,rstar,omega_vec)
    else
       particle_position = r*radial_unit_vector_rotated
       particle_velocity = v*radial_unit_vector_rotated
    endif
    particle_position = particle_position + x0
    particle_velocity = particle_velocity + v0

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
subroutine rotation_frame(vector,axis_A,axis_B)
 use vectorutils, only:cross_product3D
 real, intent(in)    :: axis_A(3),axis_B(3)
 real, intent(inout) :: vector(3)

 real :: rotation_matrix(3,3),cross_vec(3),cos_theta,sin_theta,temp_vec(3)

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

 temp_vec(1) = dot_product(rotation_matrix(1,1:3),vector(1:3))
 temp_vec(2) = dot_product(rotation_matrix(2,1:3),vector(1:3))
 temp_vec(3) = dot_product(rotation_matrix(3,1:3),vector(1:3))

 vector = temp_vec

end subroutine rotation_frame

!-----------------------------------------------------------------------
!+
!  Impose line-driven wind velocity profile for a rotating star
!  (Dwarkadas & Owocki 2002)
!+
!-----------------------------------------------------------------------
subroutine velocity_rotating_sink(radial_unit_vector,r,v,omega,omega_axis,&
                                  particle_position,particle_velocity,rstar,omega_vec)
 use vectorutils, only:cross_product3D
 use geometry,    only:vector_transform, coord_transform
 real, intent(in)  :: r,v,omega,omega_axis(3),rstar,omega_vec(3)
 real, intent(out) :: particle_position(3),particle_velocity(3)

 real :: radial_unit_vector(3),z_axis(3),particle_position_proj(3)
 real :: velocity_spherical(3),position_spherical(3),vomega_spin(3)
 integer, parameter :: ndim = 3, igeom = 3

 ! get the projected position of particles on the effective stellar surface
 particle_position_proj = rstar*radial_unit_vector

 ! compute centrifugal force due to spin of the sink particle
 call cross_product3D(omega_vec,particle_position_proj,vomega_spin)
 vomega_spin =  vomega_spin/rstar**3

 if (omega_axis(3) /= 1.) then
    !rotate frame to make the omega_axis coincide with z_axis
    z_axis = [0., 0., 1.0]
    call rotation_frame(radial_unit_vector,omega_axis,z_axis)
 endif

 particle_position = r * radial_unit_vector
 particle_velocity = v * radial_unit_vector

 ! switch from cartesian to spherical polar coordinates
 call coord_transform(particle_position, ndim, 1, position_spherical, ndim, igeom)
 call vector_transform(particle_position, particle_velocity, ndim, 1, velocity_spherical, ndim, igeom)

 ! Dwarkadas & Owocki 2002 eq. (8)
 velocity_spherical(1) = velocity_spherical(1) * sqrt(1. - omega**2 * sin(position_spherical(3))**2)

 ! switch from spherical polar to cartesian coordinates
 call vector_transform(position_spherical, velocity_spherical, ndim, igeom, particle_velocity, ndim, 1)

 if (omega_axis(3) /= 1.) then
    ! rotate the frame to make the z_axis coincide with omega_axis
    call rotation_frame(particle_position,z_axis,omega_axis)
    call rotation_frame(particle_velocity,z_axis,omega_axis)
 endif
 particle_velocity = particle_velocity + vomega_spin

end subroutine velocity_rotating_sink

!-----------------------------------------------------------------------
!+
!  Initialize reusable variables for jets properties
!+
!-----------------------------------------------------------------------
subroutine init_jets(jet_edge_vel,jet_opening_angle)
 use physcon, only:km
 use units,   only:unit_velocity
 real, intent(in) :: jet_edge_vel,jet_opening_angle

 jets = .true.
 edge_velocity = jet_edge_vel * (km / unit_velocity)
 opening_angle = jet_opening_angle

end subroutine init_jets

!---------------------------------------------------------------------------
!+
!  Impose polar jets velocity profile (Thomas et al. 2013, MNRAS, 430, 1230)
!+
!---------------------------------------------------------------------------
subroutine velocity_jets(radial_unit_vector,r,v,edge_velocity,opening_angle,&
                         particle_position,particle_velocity,vwind_terminal)
 use geometry, only:vector_transform,coord_transform
 real, intent(out)  :: particle_position(3),particle_velocity(3)
 integer, parameter :: ndim = 3, igeom = 3, p = 2
 real :: r,v,edge_velocity,opening_angle,vwind_terminal
 real :: position_spherical(3),velocity_spherical(3),radial_unit_vector(3)

 particle_position = r * radial_unit_vector
 particle_velocity = v * radial_unit_vector

 ! switch from cartesian to spherical polar coordinates
 call coord_transform(particle_position, ndim, 1, position_spherical, ndim, igeom)
 call vector_transform(particle_position, particle_velocity, ndim, 1, velocity_spherical, ndim, igeom)

 ! Thomas et al. 2013, eq. (1)
 if (position_spherical(3) < pi/2.) then
    velocity_spherical(1) = vwind_terminal + (edge_velocity - vwind_terminal) &
                            * (position_spherical(3)/opening_angle)**p
 else
    velocity_spherical(1) = vwind_terminal + (edge_velocity - vwind_terminal) &
                            * (abs(position_spherical(3)-pi)/opening_angle)**p
 endif

 ! switch from spherical polar to cartesian coordinates
 call vector_transform(position_spherical, velocity_spherical, ndim, igeom, particle_velocity, ndim, 1)

end subroutine velocity_jets

!-----------------------------------------------------------------------
!+
!  Return the optimal rotation angles given the resolution
!  (optimal = which gives the most homogeneous distribution
!  of particles)
!+
!-----------------------------------------------------------------------
subroutine optimal_rot_angles(ires,rotation_angles,rot_jets)
 use random, only:ran2
 logical, optional, intent(in) :: rot_jets
 integer, intent(in) :: ires
 real, intent(out)   :: rotation_angles(3)

! for the fibonacci sphere injection method, the rotation angles are pseudo-random
! (actually follow a sequence of length cycle 2.30584E+18), results are reproducible
 rotation_angles(3) = ran2(seed_random) * 2.* pi
 if (rot_jets) then
! rotate only about the z-axis for the polar jets
    rotation_angles(1) = 0.
    rotation_angles(2) = 0.
 else
    rotation_angles(1) = ran2(seed_random) * 2.* pi
    rotation_angles(2) = ran2(seed_random) * 2.* pi
 endif

end subroutine optimal_rot_angles

end module injectutils
