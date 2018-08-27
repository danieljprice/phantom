module bondiexact

 implicit none

 public :: get_bondi_solution

! Constants for geodesic  solution
 real, private :: den0 = 1.
 real, private :: en0  = 1.e-9

! Constants for sonic point flow solution
 real, public  :: rcrit   = 8.
 real, private :: adiabat = 1.

 logical, public :: iswind = .false.
 integer, public :: isol   = 1

 private

 real :: c1,c2,Tc,n,mass1,gamma

contains

subroutine get_bondi_solution(rho,v,u,r,mass,gamma_eos)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass,gamma_eos

 select case(isol)
 case(1)
    call bondigr_geodesic(rho,v,u,r,mass,gamma_eos)
 case(2)
    call bondigr_sonicpoint(rho,v,u,r,mass,gamma_eos)
 end select

 if (.not.iswind) v = -v

end subroutine get_bondi_solution

subroutine bondigr_geodesic(rho,v,u,r,mass,gam)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass,gam
 real :: sqrtg,alpha,dfunc,efunc

 dfunc = den0/(r**2*sqrt(2.*mass/r*(1.- 2.*mass/r)))
 efunc = en0/((sqrt(2.*mass/r)*r**2)**gam * (1.- 2.*mass/r)**((gam + 1.)/4.))

 sqrtg = 1.
 alpha = sqrt(1. - 2.*mass/r)
 rho = sqrtg/alpha*dfunc
 v   = sqrt(2.*mass/r)*(1. - 2.*mass/r)
 u   = efunc/dfunc

end subroutine bondigr_geodesic

subroutine bondigr_sonicpoint(rho,v,u,r,mass,gamma_eos)
 real, intent(out) :: rho,v,u
 real, intent(in)  :: r,mass,gamma_eos
 real :: T,uvel,term,u0,dens,sqrtg

 mass1 = mass
 gamma = gamma_eos
 rcrit = rcrit*mass1
 call compute_constants(rcrit)

 ! Given an r, solve eq 76 for T numerically
 call Tsolve(T,r)

 uvel = C1/(r**2 * T**n)
 dens = adiabat*T**n
 u = T*n

 !get u0 at r
 term = 1./(1.-2.*mass1/r)
 u0  = sqrt(term*(1.+term*uvel**2))
 v   = uvel/u0

 sqrtg = 1. !???? FIX
 rho = sqrtg*u0*dens

end subroutine bondigr_sonicpoint

subroutine compute_constants(rcrit)
 real, intent(in) :: rcrit
 real :: uc2,vc2

 n   = 1./(gamma-1.)
 uc2 = mass1/(2.*rcrit)
 vc2 = uc2/(1.-3.*uc2)
 Tc  = vc2*n/(1.+n-vc2*n*(1.+n))

 C1  = sqrt(uc2) * Tc**n * rcrit**2
 C2  = (1. + (1.+n)*Tc)**2 * (1. - 2.*mass1/rcrit + C1**2/(rcrit**4*Tc**(2.*n)))

end subroutine compute_constants

! Newton Raphson
subroutine Tsolve(T,r)
 real, intent(in)  :: r
 real, intent(out) :: T
 real    :: Tnew,diff
 logical :: converged
 integer :: its
 integer, parameter :: itsmax = 100
 real, parameter    :: tol    = 1.e-5

 ! These guess values may need to be adjusted for values of rcrit other than rcrit=8M
 if ((iswind .and. r>=rcrit) .or. (.not.iswind .and. r<rcrit)) then
    T = 0.760326*r**(-1.307)/2.   ! This guess is calibrated for rcrit=8M, and works ok up to r ~ 10^7 M
 elseif ((iswind .and. r<rcrit) .or. (.not.iswind .and. r>=rcrit)) then
    T = 100.
 endif

 converged = .false.
 its = 0
 do while (.not.converged .and. its<itsmax)
    Tnew = T - ffunc(T,r)/df(T,r)
    diff = abs(Tnew - T)/abs(T)
    converged = diff < tol
    T = Tnew
    its = its+1
 enddo

 if(.not.converged) print*,'Bondi exact solution not converged at r = ',r

end subroutine Tsolve

real function ffunc(T,r)
 real, intent(in) :: T,r
 ffunc = (1. + (1. + n)*T)**2*(1. - (2.*mass1)/r + C1**2/(r**4*T**(2.*n))) - C2
end function ffunc

real function df(T,r)
 real, intent(in) :: T,r
 df = (2.*(1. + T + n*T)*((1. + n)*r**3*(-2.*mass1 + r) - C1**2*T**(-1. - 2.*n)*(n + (-1. + n**2)*T)))/r**4
end function df

end module bondiexact
