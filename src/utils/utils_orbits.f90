!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module orbits
!
! orbits module contains functions to calculate various orbital parameters and
! to refine velocity using gradient descent method for a geodesic orbit.
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: physcon, vectorutils
!
 implicit none

 public :: escape, eccentricity_vector, eccentricity_star
 public :: semimajor_axis, period_star, orbital_angles
 public :: isco_kerr, refine_velocity

 private

contains

 !----------------------------------------------------------------
 !+
 !  This function tells if the star escaped the black hole or not
 !+
 !----------------------------------------------------------------
logical function escape(velocity_on_orbit,central_obj_m,position_of_obj)

 ! we require the position of object on orbit/ wrt central object, velocity on
 ! orbit is velocity of object wrt central obj, central_obh_m is mass of central
 ! object

 use physcon,only : gg

 real, intent(in) :: velocity_on_orbit,central_obj_m,position_of_obj
 real             :: escape_vel


 escape_vel = sqrt((2.*gg*central_obj_m)/(position_of_obj))
 if (velocity_on_orbit > escape_vel) then
    print*,'star has escaped',(velocity_on_orbit)/escape_vel,'vel bh/escape vel'
    print*, velocity_on_orbit/(1e5),"vel of star in Km/s"
    escape = .true.
 else
    escape = .false.
    print*, velocity_on_orbit/(1e5),"in Km/s"
 endif
end function escape

 !----------------------------------------------------------------
 !+
 !  This function calculates h vector.
 !+
 !----------------------------------------------------------------
function hvector(pos_vec,vel_vec)
 use vectorutils,     only : cross_product3D

 real,intent(in) :: pos_vec(3),vel_vec(3)
 real,dimension(3) :: hvector

 call cross_product3D(pos_vec,vel_vec,hvector)

end function hvector

 !----------------------------------------------------------------
 !+
 !  This function calculates h vector.
 !+
 !----------------------------------------------------------------
function vcrossh(pos_vec,vel_vec)
 use vectorutils,     only : cross_product3D

 real,intent(in) :: vel_vec(3),pos_vec(3)
 real,dimension(3) :: h_vector,vcrossh

 h_vector = hvector(pos_vec,vel_vec)
 call cross_product3D(vel_vec,h_vector,vcrossh)

end function vcrossh

 !----------------------------------------------------------------
 !+
 !  This function calculates eccentricity vector of the orbit.
 !+
 !----------------------------------------------------------------
function  eccentricity_vector(mass1,mass2,pos_vec,vel_vec)
 use physcon,only : gg,pi

 real, intent(in)  :: pos_vec(3),vel_vec(3)
 real, intent(in)  :: mass1,mass2
 real, dimension(3) :: eccentricity_vector
 real, dimension(3) :: vcrossh_vec
 real               :: pos_vec_mag

 vcrossh_vec = vcrossh(pos_vec,vel_vec)
 pos_vec_mag = sqrt(dot_product(pos_vec,pos_vec))

 ! eccentricity vector = (rdot cross h )/(G(m1+m2)) - rhat
 eccentricity_vector(:) = (vcrossh_vec(:)/(gg*(mass1+mass2))) - (pos_vec/pos_vec_mag)

end function eccentricity_vector

 !----------------------------------------------------------------
 !+
 !  This function calculates eccentricity of the orbit.
 !+
 !----------------------------------------------------------------
real function eccentricity_star(mass1,mass2,pos_vec,vel_vec)

 real, intent(in)  :: pos_vec(3),vel_vec(3)
 real, intent(in)  :: mass1,mass2
 real, dimension(3) :: ecc_vector
 ecc_vector = eccentricity_vector(mass1,mass2,pos_vec,vel_vec)

 eccentricity_star = sqrt(dot_product(ecc_vector,ecc_vector))

end function eccentricity_star

 !----------------------------------------------------------------
 !+
 !  This function calculates semimajor axis of the orbit in cm.
 !+
 !----------------------------------------------------------------
real function semimajor_axis(mass1,mass2,pos_vec,vel_vec)

 use physcon,only : gg,pi

 real, intent(in) :: mass1,mass2
 real, intent(in) :: pos_vec(3),vel_vec(3)
 real :: h2,eccentricity_value,h_vector(3)

 eccentricity_value = eccentricity_star(mass1,mass2,pos_vec,vel_vec)
 h_vector = hvector(pos_vec,vel_vec)
 h2 = dot_product(h_vector,h_vector)

 !formula used is a = h^2/(G(M*+M_BH)*(1-e^2))
 semimajor_axis = h2/((gg*(mass1+mass2))*(1-eccentricity_value**2))


end function semimajor_axis

 !----------------------------------------------------------------
 !+
 !  This function calculates period of the orbit in seconds.
 !+
 !----------------------------------------------------------------
real function period_star(mass1,mass2,pos_vec,vel_vec)

 use physcon,only : gg,pi

 real, intent(in) :: mass1,mass2
 real, intent(in) :: pos_vec(3),vel_vec(3)
 real :: semimajor

 semimajor = semimajor_axis(mass1,mass2,pos_vec,vel_vec)

 ! using Kepler's 3rd law to calculated period
 period_star = sqrt((4*pi**2*abs(semimajor)**3)/(gg*(mass1+mass2)))

end function period_star

 !----------------------------------------------------------------
 !+
 !  This routine calculated the inclination angle,
 !  argument of pariestron, longitude ascending node.
 !+
 !----------------------------------------------------------------
subroutine orbital_angles(mass1,mass2,pos_vec,vel_vec,&
                            inclination_angle,argument_of_periestron,longitude_ascending_node)

 use vectorutils, only : cross_product3D

 real, intent(in)   :: mass1,mass2
 real, intent(in)   :: pos_vec(3),vel_vec(3)
 real               :: i_vector(3),j_vector(3),k_vector(3),n_vector(3)
 real               :: n_vector_mag,h_val,h_vector(3),ecc_vec(3),ecc_val,ecc_hat(3),n_hat(3),h_hat(3)
 real, intent(out)  :: inclination_angle,argument_of_periestron,longitude_ascending_node

 h_vector = hvector(pos_vec,vel_vec)
 h_val = sqrt(dot_product(h_vector,h_vector))
 h_hat = h_vector(:)/h_val

 ecc_vec = eccentricity_vector(mass1,mass2,pos_vec,vel_vec)
 ecc_val = eccentricity_star(mass1,mass2,pos_vec,vel_vec)
 ecc_hat = ecc_vec(:)/ecc_val

 i_vector = (/1,0,0/)
 j_vector = (/0,1,0/)
 k_vector = (/0,0,1/)

 ! n vector = k cross h
 call cross_product3D(k_vector,h_vector,n_vector)
 n_vector_mag = sqrt(dot_product(n_vector,n_vector))
 n_hat = n_vector/n_vector_mag

 ! calculate inclination angle
 inclination_angle = acos(dot_product(k_vector,h_hat))
 argument_of_periestron = acos(dot_product(ecc_hat,n_hat))
 longitude_ascending_node = acos(dot_product(i_vector,n_hat))


end subroutine orbital_angles

 !----------------------------------------------------------------
 !+
 !  This subroutine calculates the kerr metric's Innermost stable
 !  circular orbit. These formulas were obtained from
 !  Bardeen & Teukolsky 1972 Astrophysical Journal, Vol. 178, pp. 347-370
 !+
 !----------------------------------------------------------------

subroutine isco_kerr(a,mass_bh,r_isco)
 real, intent(in) :: a,mass_bh
 real, intent(out):: r_isco
 real  :: z1,z2

 !Modified the formulas for Z1 and Z2 so that it uses the spin parameter
 !implemented in PHANTOM, which is a/M.

 z1 = 1 + (1-a**2.)**(1./3.)*((1+a)**(1./3.) + (1-a)**(1./3.))
 z2 = ((3*a**2.) + z1**2.)**(1./2.)

 !now we check if the value of a is +ve or -ve
 !+ve implies prograde while -ve implies retrograde

 if (a>=0.) then
    print*, "prograde rotation of BH wrt orbit"
    r_isco = mass_bh*(3. + z2 - sqrt((3.-z1)*(3.+z1+2.*z2)))
 else
    r_isco = mass_bh*(3. + z2 + sqrt((3.-z1)*(3.+z1+2.*z2)))
    print*, "retrograde rotation of BH wrt orbit"
 endif
 print*,"----------------------------------------------------------"
 print*, "ISCO of KERR metric with spin ",a," is: ",r_isco
 print*,"----------------------------------------------------------"
end subroutine isco_kerr
!--------------------------------------------------------------------------
!+
!  Obtain velocity in Boyer–Lindquist coordinates using gradient descent method. Iterations are
!  performed until the target specific energy epsilon is reached
!  References: Tejeda et al. (2017), MNRAS 469, 4483–4503
!  Owner: Mario Aguilar Faúndez
!+
!--------------------------------------------------------------------------
subroutine refine_velocity(x, y, z, vx, vy, vz, M_h, a, r, epsilon_target, alpha, delta_v, tol, max_iters)
 real, intent(in) :: x, y, z, M_h, a, r, epsilon_target, alpha, delta_v, tol
 integer, intent(in) :: max_iters
 real, intent(inout) :: vx, vy, vz
 real :: epsilon_0, epsilon_x, epsilon_y
 real :: d_eps_dx, d_eps_dy
 real :: temp_vx, temp_vy
 real :: sign_epsilon
 integer :: iter


 print*, 'Initial velocities: vx = ', vx, ', vy = ', vy, ', vz = ', vz, 'velocity magnitude = ', sqrt(vx**2 + vy**2 + vz**2)
 iter = 0
 do while (iter < max_iters)
    ! Compute epsilon at current velocity
    call compute_epsilon(x, y, z, vx, vy, vz, M_h, a, r, epsilon_0)

    ! Check convergence
    if (abs(epsilon_0 - epsilon_target) < tol) exit

    ! Determine sign of epsilon - epsilon_target
    sign_epsilon = sign(1.0, epsilon_0 - epsilon_target)

    ! We only iterate in the x-y plane.
    ! Compute epsilon_x by changing vx by a small delta
    temp_vx = vx + delta_v
    call compute_epsilon(x, y, z, temp_vx, vy, vz, M_h, a, r, epsilon_x)
    d_eps_dx = (epsilon_x - epsilon_0) / delta_v

    ! Compute epsilon_y by changing vy by a small delta
    temp_vy = vy + delta_v
    call compute_epsilon(x, y, z, vx, temp_vy, vz, M_h, a, r, epsilon_y)
    d_eps_dy = (epsilon_y - epsilon_0) / delta_v

    ! Update velocities using gradient descent with sign adjustment
    vx = vx - alpha * sign_epsilon * d_eps_dx
    vy = vy - alpha * sign_epsilon * d_eps_dy

    iter = iter + 1
 enddo
 print*, 'Total number of iterations for refining velocity: ', iter
 print*, 'Updated velocities: vx = ', vx, ', vy = ', vy, ', vz = ', vz, 'velocity magnitude = ', sqrt(vx**2 + vy**2 + vz**2)
 call compute_epsilon(x, y, z, vx, vy, vz, M_h, a, r, epsilon_0)
 if (iter == max_iters) then
    print *, 'Warning: Gradient descent did not converge after ', max_iters, ' iterations.'
 endif
end subroutine refine_velocity

!--------------------------------------------------------------------------
!+
!  Compute GR total specific energy
!  References: Tejeda et al. (2017), MNRAS 469, 4483–4503
!+
!--------------------------------------------------------------------------
subroutine compute_epsilon(x, y, z, vx, vy, vz, M_h, a, r, epsilon_0)
 real, intent(in) :: x, y, z, vx, vy, vz, M_h, a, r
 real, intent(out) :: epsilon_0
 real :: xy_term, r_dot_term, Gamma_flat, Gamma_curved, Gamma_0
 real :: rho_squared, Delta

 rho_squared = r**2 + (a**2 * z**2) / r**2
 Delta = r**2 - 2 * M_h * r + a**2
 xy_term = x * vy - y * vx
 r_dot_term = r**2 * (x * vx + y * vy) + (r**2 + a**2) * z * vz
 Gamma_flat = 1.0d0 - vx**2 - vy**2 - vz**2
 Gamma_curved = (2.0d0 * M_h * r / rho_squared) * ((1.0d0 - a * xy_term / (r**2 + a**2))**2 + &
                        r_dot_term**2 / (r**2 * Delta * (r**2 + a**2)))
 Gamma_0 = 1.0d0 / sqrt(Gamma_flat - Gamma_curved)

 epsilon_0 = Gamma_0 * (1.0d0 - (2.0d0 * M_h * r / rho_squared) * (1.0d0 - a * xy_term / (r**2 + a**2)))
end subroutine compute_epsilon
end module orbits
