module bondiexact

 implicit none

 public :: get_bondi_solution

 real, private :: den0 = 1.    !  12.
 real, private :: en0  = 1.e-9 !  0.000297118

 ! This is not used. It only exists for compatability reasons
 real, public :: rcrit

 private

contains

subroutine get_bondi_solution(rho,v,u,r,mass1,gamma)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass1,gamma
 real :: sqrtg,alpha,dfunc,efunc

 dfunc = den0/(r**2*sqrt(2.*mass1/r*(1.- 2.*mass1/r)))
 efunc = en0/((sqrt(2.*mass1/r)*r**2)**gamma * (1.- 2.*mass1/r)**((gamma + 1.)/4.))

 sqrtg = 1.
 alpha = sqrt(1. - 2.*mass1/r)
 rho = sqrtg/alpha*dfunc
 v   = sqrt(2.*mass1/r)*(1. - 2.*mass1/r)
 u   = efunc/dfunc

end subroutine get_bondi_solution

end module bondiexact
