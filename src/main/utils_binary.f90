!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: binaryutils
!
!  DESCRIPTION:
!   Utility routines for binary setup, also used in asteroid injection
!
!  REFERENCES: None
!
!  OWNER: Bec Nealon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module binaryutils
 implicit none
 real, parameter :: pi = 4.*atan(1.)

contains

!---------------------------------------------------------------
!+
! Get eccentric anomaly (this function uses bisection
! to guarantee convergence, which is not guaranteed for
! small M or E)
!+
!---------------------------------------------------------------
subroutine get_E(period,ecc,deltat,E)
 real, intent(in)  :: period,ecc,deltat
 real, intent(out) :: E
 real :: mu,M_ref,M_guess
 real :: E_left,E_right,E_guess
 real, parameter :: tol = 1.e-10

 mu = 2.*pi/period
 M_ref = mu*deltat ! mean anomaly

 ! first guess
 E_left = 0.
 E_right = 2.*pi
 E_guess = pi
 M_guess = M_ref - 2.*tol

 do while (abs(M_ref - M_guess) > tol)
   M_guess = E_guess - ecc*sin(E_guess)
   if (M_guess > M_ref) then
      E_right = E_guess
   else
      E_left = E_guess
   endif
   E_guess = 0.5*(E_left + E_right)
 enddo

 E = E_guess

end subroutine get_E

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
