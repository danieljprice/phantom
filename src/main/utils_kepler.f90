!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_kepler
!
! utils_kepler
!
! :References: None
!
! :Owner: Yrisch
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 use physcon,only: pi
 implicit none

contains
subroutine Espec(v2,r,mu,B)
 real, intent(in) :: v2,r,mu
 real, intent(out) ::  B

 B = 0.5*v2 - mu/r

end subroutine Espec

subroutine extract_a(r,mu,v2,aij)
 real, intent(in) :: r,mu,v2
 real, intent(out):: aij
 aij = (r*mu)/(2.*mu-r*v2)

end subroutine extract_a

subroutine extract_a_dot(r2,r,mu,v2,v,acc,adot)
 real, intent(in)    :: r2,r,mu,v2,v,acc
 real, intent(out)   :: adot
 real :: mu2
 mu2 = mu**2
 adot = 2.*(mu2*v+r2*v*acc)/((2.*mu-r*v2)**2)
end subroutine extract_a_dot

subroutine extract_e(x,y,z,vx,vy,vz,mu,r,eij)
 real, intent(in) :: x,y,z,vx,vy,vz,mu,r
 real, intent(out):: eij
 real :: eijx,eijy,eijz
 real :: hx,hy,hz

 hx = y*vz-z*vy
 hy = z*vx-x*vz
 hz = x*vy-y*vx

 eijx = (vy*hz-vz*hy)/mu - x/r
 eijy = (vz*hx-vx*hz)/mu - y/r
 eijz = (vx*hy-hx*vy)/mu - z/r

 eij = sqrt(eijx**2+eijy**2+eijz**2)

end subroutine extract_e

subroutine extract_ea(x,y,z,vx,vy,vz,mu,aij,eij)
 real, intent(in) :: x,y,z,vx,vy,vz,mu,aij
 real, intent(out):: eij
 real :: hx,hy,hz,h2,neg_e

 hx = y*vz-z*vy
 hy = z*vx-x*vz
 hz = x*vy-y*vx

 h2 = hx**2+hy**2+hz**2

 neg_e = h2/(mu*aij)
 print*,neg_e
 if (neg_e>=1) then
    eij = 0.
 else
    eij = sqrt(1-neg_e)
 endif

end subroutine extract_ea

subroutine extract_kep_elmt(x,y,z,vx,vy,vz,mu,r,a,e,i,argp,longi,M)
 real, intent(in) :: x,y,z,vx,vy,vz,mu,r
 real, intent(out):: a,e,i,argp,longi,M
 real :: hx,hy,hz,ex,ey,ez,v2,h,anoE,nu
 real :: rdote,n,ndote

 v2 = vx**2+vy**2+vz**2

 a = (r*mu)/(2*mu-r*v2)

 hx = y*vz-z*vy
 hy = z*vx-x*vz
 hz = x*vy-y*vx

 h = sqrt(hx*2+hy**2+hz**2)
 i = acos(hz/h)

 ex = (vy*hz-vz*hy)/mu - x/r
 ey = (vz*hx-vx*hz)/mu - y/r
 ez = (vx*hy-hx*vy)/mu - z/r

 e = sqrt(ex**2+ey**2+ez**2)

 rdote = x*ex+y*ey+z*ez

 if (x*vx+y*vy+z*vz>=0) then
    nu = acos(rdote/(e*r))
 else
    nu = 2*pi - acos(rdote/(e*r))
 endif
 anoE = tan(nu*0.5)/sqrt((1+e)/(1-e))
 anoE = 2*atan(anoE)

 M = E-e*sin(E)

 n = sqrt(hy**2+hx**2)
 if (hx>=0) then
    longi = acos(-hy/n)
 else
    longi = 2*pi - acos(-hy/n)
 endif

 ndote = -hy*ex + hx*ey
 if (ez>=0) then
    argp = acos(ndote/(n*e))
 else
    argp = 2*pi - acos(ndote/(n*e))
 endif

end subroutine extract_kep_elmt




end module utils_kepler
