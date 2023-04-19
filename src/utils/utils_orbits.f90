!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module orbits_data
!
! orbits_data
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: physcon, units, vectorutils
!

 implicit none

 public :: escape, orbital_parameters, eccentricity_vector_v, eccentricity_star
 public :: semimajor_axis, period_star, orbital_angles
 public :: isco_kerr

 private

contains

 !----------------------------------------------------------------
 !+
 !  This function tells if the star escaped the black hole or not
 !+
 !----------------------------------------------------------------
logical function escape(velocity_bh,bh_mass,position_bh)

 use units , only : umass
 use physcon,only : gg

 real, intent(in) :: velocity_bh, position_bh, bh_mass
 real             :: escape_vel


 escape_vel = sqrt((2.*gg*bh_mass*umass)/position_bh)

 if (velocity_bh > escape_vel) then
    print*,'star has escaped',velocity_bh/escape_vel,'vel bh/escape vel'
    escape = .true.
 else
    print*,'star is bound',velocity_bh/escape_vel,'vel bh/escape vel'
    escape = .false.
 endif
end function escape
 !----------------------------------------------------------------
 !+
 !  This routine calculates orbital parameters.
 !+
 !----------------------------------------------------------------
subroutine orbital_parameters(angular_momentum_h,bh_mass,mass_star,com_star,position_bh,velocity_wrt_bh)

 use vectorutils,     only : cross_product3D
 use physcon,         only : years,pc,pi

 real, intent(in)   :: angular_momentum_h(3),com_star(3),velocity_wrt_bh(3)
 real, intent(in)   :: bh_mass, mass_star, position_bh
 real               :: eccentricity_value,semimajor_value,period_value
 real               :: h_value
 real               :: vcrossh(3),eccentricity_vector(3)
 real               :: inclination_angle,argument_of_periestron,longitude_ascending_node

 !h_value is the magnitude of angular_momentum_h vector squared.
 h_value = dot_product(angular_momentum_h,angular_momentum_h) !square of the magnitude of h vector
 print*,acos(abs(angular_momentum_h(1))/h_value),'inclination x',abs(angular_momentum_h(1))/h_value
 print*,acos(abs(angular_momentum_h(2))/h_value),'inclination y',abs(angular_momentum_h(2))/h_value
 print*,acos(abs(angular_momentum_h(3))/h_value),'inclination z',abs(angular_momentum_h(3))/h_value
 print*,'h_value',h_value,angular_momentum_h(1),angular_momentum_h(2),angular_momentum_h(3)

 ! semimajor = h_value/((gg*(mass(ngrid-1)+bh_mass)*umass)*(1-eccentricity_star**2))
 ! period_star = sqrt((4*pi**2*abs(semimajor)**3)/(gg*(mass(ngrid-1)+bh_mass)*umass))
 ! print*,semimajor/pc,'semimajor axis in Parsec',period_star/years,'period in years'

 call cross_product3D(velocity_wrt_bh,angular_momentum_h,vcrossh)

 call eccentricity_vector_v(vcrossh, mass_star, bh_mass, com_star, position_bh, eccentricity_vector)

 eccentricity_value = eccentricity_star(eccentricity_vector)
 print*, eccentricity_value, "eccentricity"

 semimajor_value = semimajor_axis(h_value,mass_star,bh_mass,eccentricity_value)
 print*, semimajor_value/pc, "value of semi-major axis `a` in Parsec"

 period_value = period_star(semimajor_value,mass_star,bh_mass)
 print*, period_value/years, "Period of the orbit in years"

 call orbital_angles(angular_momentum_h,eccentricity_vector,eccentricity_value,h_value,&
                              inclination_angle,argument_of_periestron,longitude_ascending_node)

 print*,inclination_angle*(180/pi),'i',argument_of_periestron*(180/pi),'w',longitude_ascending_node*(180/pi),'omega'

end subroutine orbital_parameters

 !----------------------------------------------------------------
 !+
 !  This function calculates eccentricity vector of the orbit.
 !+
 !----------------------------------------------------------------
subroutine eccentricity_vector_v(vcrossh, mass_star, bh_mass, com_star, position_bh,eccentricity_vector)
 use units , only : umass
 use physcon,only : gg,pi

 real, intent(in)  :: vcrossh(3), com_star(3)
 real, intent(in)  :: mass_star, bh_mass, position_bh
 real, intent(out) :: eccentricity_vector(3)

 !eccentricity vector = (rdot cross h )/(G(m1+m2)) - rhat
 eccentricity_vector(:) = (vcrossh(:)/(gg*(mass_star+bh_mass)*umass)) - (com_star(:)/position_bh)

end subroutine eccentricity_vector_v
 !----------------------------------------------------------------
 !+
 !  This function calculates eccentricity of the orbit.
 !+
 !----------------------------------------------------------------
real function eccentricity_star(eccentricity_vector)

 real, intent(in) :: eccentricity_vector(3)

 eccentricity_star = sqrt(dot_product(eccentricity_vector,eccentricity_vector))

end function eccentricity_star

 !----------------------------------------------------------------
 !+
 !  This function calculates semimajor axis of the orbit.
 !+
 !----------------------------------------------------------------
real function semimajor_axis(h2_value,mass_star,bh_mass,eccentricity_value)

 use units , only : umass
 use physcon,only : gg,pi

 real, intent(in) :: h2_value,mass_star,bh_mass,eccentricity_value

 !formula used is a = h^2/(G(M*+M_BH)*(1-e^2))
 semimajor_axis = h2_value/((gg*(mass_star+bh_mass)*umass)*(1-eccentricity_value**2))

end function semimajor_axis

 !----------------------------------------------------------------
 !+
 !  This function calculates period of the orbit.
 !+
 !----------------------------------------------------------------
real function period_star(semimajor_value,mass_star,bh_mass)

 use units , only : umass
 use physcon,only : gg,pi

 real, intent(in) :: semimajor_value,mass_star,bh_mass


 period_star = sqrt((4*pi**2*abs(semimajor_value)**3)/(gg*(mass_star+bh_mass)*umass))

end function period_star

 !----------------------------------------------------------------
 !+
 !  This routine calculated the inclination angle,
 !  argument of pariestron, longitude ascending node.
 !+
 !----------------------------------------------------------------
subroutine orbital_angles(angular_momentum_h,eccentricity_vector,eccentricity_value,h_value,&
                            inclination_angle,argument_of_periestron,longitude_ascending_node)

 use vectorutils, only : cross_product3D

 real, intent(in)   :: angular_momentum_h(3),eccentricity_vector(3)
 real , intent(in)  :: h_value,eccentricity_value
 real               :: i_vector(3),j_vector(3),k_vector(3),n_vector(3)
 real               :: n_vector_mag
 real, intent(out)  :: inclination_angle,argument_of_periestron,longitude_ascending_node


 i_vector = (/1,0,0/)
 j_vector = (/0,1,0/)
 k_vector = (/0,0,1/)

 call cross_product3D(k_vector,angular_momentum_h,n_vector)
 n_vector_mag = sqrt(dot_product(n_vector,n_vector))
 inclination_angle = acos(dot_product(k_vector,angular_momentum_h/sqrt(h_value)))
 argument_of_periestron = acos(dot_product(eccentricity_vector/eccentricity_value,n_vector/n_vector_mag))
 longitude_ascending_node = acos(dot_product(j_vector,n_vector/n_vector_mag))


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

end module orbits_data
