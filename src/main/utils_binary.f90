!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module binaryutils
!
! Utility routines for binary setup, also used in asteroid injection
! Main thing is routines to compute eccentric anomaly from mean anomaly
! by solving Kepler's equation
!
! :References: None
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 real, parameter :: pi = 4.*atan(1.)

 public :: get_E,get_E_from_mean_anomaly
 public :: get_orbit_bits

contains

!---------------------------------------------------------------
!+
! Get eccentric anomaly given time since pericentre
!+
!---------------------------------------------------------------
subroutine get_E(period,ecc,deltat,E)
 real, intent(in)  :: period,ecc,deltat
 real, intent(out) :: E
 real :: mu,M_ref

 mu = 2.*pi/period
 M_ref = mu*deltat ! mean anomaly

 E = get_E_from_mean_anomaly(M_ref,ecc)

end subroutine get_E

!---------------------------------------------------------------
!+
! Get eccentric anomaly from mean anomaly (this function uses
! bisection to guarantee convergence, which is not guaranteed for
! small M or E)
!+
!---------------------------------------------------------------
real function get_E_from_mean_anomaly(M_ref,ecc) result(E)
 real, intent(in) :: M_ref,ecc
 real :: E_left,E_right,E_guess,M_guess
 real, parameter :: tol = 1.e-10

 ! first guess
 E_left = 0.
 E_right = 2.*pi
 E_guess = pi
 M_guess = M_ref - 2.*tol

 do while (abs(M_ref - M_guess) > tol)
    if (ecc < 1.) then     ! eccentric
       M_guess = E_guess - ecc*sin(E_guess)
    elseif (ecc > 1.) then ! hyperbolic
       M_guess = ecc*sinh(E_guess) - E_guess
    else                   ! parabolic
       M_guess = E_guess + 1./3.*E_guess**3
    endif
    if (M_guess > M_ref) then
       E_right = E_guess
    else
       E_left = E_guess
    endif
    E_guess = 0.5*(E_left + E_right)
 enddo

 E = E_guess

end function get_E_from_mean_anomaly

!---------------------------------------------------------------
!+
!  Get eccentric (or parabolic/hyperbolic) anomaly from true anomaly
!  https://space.stackexchange.com/questions/23128/design-of-an-elliptical-transfer-orbit/23130#23130
!+
!---------------------------------------------------------------
real function get_E_from_true_anomaly(theta,ecc) result(E)
 real, intent(in) :: theta  ! true anomaly in radians
 real, intent(in) :: ecc    ! eccentricity

 if (ecc < 1.) then
    E = atan2(sqrt(1. - ecc**2)*sin(theta),(ecc + cos(theta)))
 elseif (ecc > 1.) then ! hyperbolic
    !E = atanh(sqrt(ecc**2 - 1.)*sin(theta)/(ecc + cos(theta)))
    E = 2.*atanh(sqrt((ecc - 1.)/(ecc + 1.))*tan(0.5*theta))
 else ! parabolic
    E = tan(0.5*theta)
 endif

end function get_E_from_true_anomaly

!-----------------------------------------------------------------------
!+
!  Calculate semi-major axis, ecc, ra and rp from radius(3), velocity(3)
!  mass of central object and iexternalforce (for LT corrections)
!+
!-----------------------------------------------------------------------
subroutine get_orbit_bits(vel,rad,m1,iexternalforce,semia,ecc,ra,rp)
 real, intent(in)    :: m1, vel(3), rad(3)
 integer, intent(in) :: iexternalforce
 real, intent(out)   :: semia, ecc, ra, rp
 real                :: speed, r, L_mag
 real                :: spec_energy,L(3),term

 L(1) = rad(2)*vel(3) - rad(3)*vel(2)
 L(2) = rad(3)*vel(1) - rad(1)*vel(3)
 L(3) = rad(1)*vel(2) - rad(2)*vel(1)
 L_mag = sqrt(dot_product(L,L))

 speed = sqrt(dot_product(vel,vel))
 r = sqrt(dot_product(rad,rad))

 spec_energy = 0.5*speed**2 - (1.0*m1/r)
 term = 2.*spec_energy*L_mag**2/(m1**2)

 if (iexternalforce == 11) then
    spec_energy = spec_energy - (3.*m1/(r**2))
    term = 2.*spec_energy*(L_mag**2 - 6.*m1**2)/(m1**2)
 endif

 semia     = -m1/(2.0*spec_energy)
 ecc = sqrt(1.0 + term)

 ra = semia*(1. + ecc)
 rp = semia*(1. - ecc)

end subroutine get_orbit_bits

end module binaryutils
