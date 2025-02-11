!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module bondiexact
!
! None
!
! :References: None
!
! :Owner: David Liptai
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 use physcon,        only:pi
 implicit none

 public :: get_bondi_solution

 real, public  :: rcrit   = 5.
 real, private :: rhocrit = 1.

 logical, public :: iswind = .false.
 integer, public :: isol   = 0

 private

contains

subroutine get_bondi_solution(rho,v,u,r,mass,gamma)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass,gamma
 real :: cs2,vr,mdot

 cs2  = mass/(2.*rcrit)
 mdot = rhocrit*4.*pi*rcrit**2*sqrt(cs2)

 if (r>=rcrit) then
    vr = sqrt(-cs2*lambertw_0(-func(r,rcrit)))
 else
    vr = sqrt(-cs2*lambertw_neg1(-func(r,rcrit)))
 endif

 v   = vr
 rho = mdot/(4.*pi*abs(v)*r**2)
 u   = 0.

end subroutine get_bondi_solution

!=====----------------------------------------------------------------------------------------------=====
!
!---Functions for non-relativistic solution
!
!--------------------------------------------------------------------------------------------------------
! See Barry, Parlange & Li 2000, for the analytic approximations used for the Lambert W function.
!
!--- Lambert W function for the principal branch (k=0) approaching from the negative
real function lambertw_0(x)
 real, intent(in) :: x
 real, parameter  :: exp1 = exp(1.)
 real :: eta, N1, N2
 eta = 2. + 2.*exp1*x
 N2  = 6. + 3.*Sqrt(2.) - ((-5764. - 4108.*Sqrt(2.) + (2237. + 1457.*Sqrt(2.))*exp1)*eta)/&
      (-796. - 430.*Sqrt(2.) + (215. + 199.*Sqrt(2.))*exp1)
 N1  = (1. - 1./Sqrt(2.))*(Sqrt(2.) + N2)
 lambertw_0 = -1 + Sqrt(eta)/(1 + (Sqrt(eta)*N1)/(Sqrt(eta) + N2))
end function lambertw_0

!--- Lambert W function for the k=-1 branch
real function lambertw_neg1(x)
 real, intent(in) :: x
 real, parameter  :: M1 = 0.3361, M2 = -0.0042, M3 = -0.0201
 real :: sigma
 sigma = -1. - Log(-x)
 lambertw_neg1 = -1. - sigma - (2.*(1. - 1./(1. + (M1*Sqrt(sigma)*(1. + exp(M3*Sqrt(sigma))*M2*sigma))/Sqrt(2.))))/M1
end function lambertw_neg1

!--- Function used in the non-rel solution of velocity [D(r)] (See: Cranmer 2004)
real function func(r,rc)
 real, intent(in) :: r,rc
 func = (rc/r)**4 * exp(4.*(1.-rc/r) - 1.)
end function func

end module bondiexact
