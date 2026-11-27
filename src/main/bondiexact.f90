!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module bondiexact
!
! Exact solution for Bondi accretion/wind problem in non-relativistic limit
!
! :References: Barry, Parlange & Li (2000)
!              Cranmer (2004), Am. J. Phys 72, 1397-1403
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: mathfunc, physcon
!
 use physcon,        only:pi
 use mathfunc,       only:lambertw
 implicit none

 public :: get_bondi_solution

 real, public  :: rcrit   = 5.
 real, private :: rhocrit = 1.

 logical, public :: iswind = .false.
 integer, public :: isol   = 0

 private

contains

subroutine get_bondi_solution(rho,v,u,r,mass,gamma)
!
! compute the exact solution for non-relativistic Bondi flow
!
! :Arguments:
!   - rho : density
!   - v   : velocity
!   - u   : internal energy
!   - r   : radius
!   - mass : mass
!   - gamma : adiabatic index
!
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass,gamma
 real :: cs2,vr,mdot

 cs2  = mass/(2.*rcrit)
 mdot = rhocrit*4.*pi*rcrit**2*sqrt(cs2)

 if (r>=rcrit) then
    vr = sqrt(-cs2*lambertw(-func(r,rcrit),0))
 else
    vr = sqrt(-cs2*lambertw(-func(r,rcrit),-1))
 endif

 v   = vr
 rho = mdot/(4.*pi*abs(v)*r**2)
 u   = 0.

end subroutine get_bondi_solution

!--------------------------------------------------------------------------------------------------
!
!---Functions for non-relativistic solution
!
!--------------------------------------------------------------------------------------------------
!--- Function used in the non-rel solution of velocity [D(r)] (See: Cranmer 2004)
real function func(r,rc)
 real, intent(in) :: r,rc
 func = (rc/r)**4 * exp(4.*(1.-rc/r) - 1.)
end function func

end module bondiexact
