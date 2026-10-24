!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module stellar_hcp
!
! Radially mapped hexagonal close packing with density-dependent masses.
! Construction adapted from StarSmasher, Jamie Lombardi and collaborators
! (GPL-3.0); retain attribution when redistributing this implementation.
! Follows the parent-star construction by Jamie Lombardi in StarSmasher:
! outer support boundary, paired HCP sites, natural cubic interpolation,
! radial mapping and m(r) proportional to rho(r)**(1-q).
!
! :References: Gaburov, Lombardi & Portegies Zwart (2010), MNRAS 402, 105
!   Single-exponent generator: StarSmasher revision 140a5553f31c37c2d0d879728296c9a9afa82146
!
! :Owner: Francesco Radica
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 private
 public :: generate_stellar_hcp,natural_spline,spline_value

contains

subroutine generate_stellar_hcp(r,den,menc,rstar,mstar,np,q,nnopt,surface_factor,hfloor, &
                               xyz,h,pmass,rhoi,redge,a1,ierr,surface_buffer)
 use, intrinsic :: ieee_arithmetic, only:ieee_is_finite
 real, intent(in) :: r(:),den(:),menc(:),rstar,mstar,q,surface_factor,hfloor
 integer, intent(in) :: np,nnopt
 real, allocatable, intent(out) :: xyz(:,:),h(:),pmass(:),rhoi(:)
 real, intent(out) :: redge,a1
 integer, intent(out) :: ierr
 real, intent(in), optional :: surface_buffer
 real, parameter :: pi = 3.1415926535897932384626d0
 real :: hc,redge1,redge2,rhomax,integral,rpreedge,cellvolume,mcentral
 real :: xtry,ytry,ztry,rtry,rpost,ri,amtot,anumden
 real :: edge_factor,rhomid,wmid
 real :: rpre(size(r)),rho2(size(r)),rmap2(size(r))
 integer :: n,i,ix,iy,iz,ixmax,iymax,izmax,nmin,naccepted,ip,pass

 ierr = 1
 if (np <= 0 .or. nnopt <= 0 .or. size(r) < 3) return
 if (size(den) /= size(r) .or. size(menc) /= size(r)) return
 if (.not.all(ieee_is_finite(r)) .or. .not.all(ieee_is_finite(den))) return
 if (.not.all(ieee_is_finite(menc)) .or. any(den < 0.)) return
 if (menc(1) < 0. .or. any(menc(2:) < menc(:size(r)-1))) return
 if (any(r(2:) <= r(:size(r)-1))) return
 if (r(1) < 0.) return
 if (.not.all(ieee_is_finite([rstar,mstar,q,surface_factor,hfloor]))) return
 if (q < 0. .or. q > 1. .or. mstar <= 0. .or. rstar <= 0.) return
 if (surface_factor <= 0. .or. hfloor < 0.) return
 if (present(surface_buffer)) then
    if (.not.ieee_is_finite(surface_buffer)) return
    if (surface_buffer /= -1. .and. surface_buffer <= 0.) return
 endif
 edge_factor = 2.
 if (present(surface_buffer)) then
    if (surface_buffer > 0.) edge_factor = surface_buffer
 endif
 n = np
 rhomax = den(1)
 if (rhomax <= 0.) return

 redge1 = rstar
 do i=1,size(r)
    if (menc(i) > 0.85*mstar) then
       redge1 = r(i)
       exit
    endif
 enddo
 hc = 0.5*rstar*(surface_factor*nnopt/n)**(1./3.)
 redge2 = rstar - edge_factor*hc
 redge = max(redge1,redge2)
 nmin = max(156,int(12.*(0.5*redge/hc)**3))
 n = max(n,nmin)

 call natural_spline(r,den,rho2)
 if (q /= 0.) then
    integral = 0.
    i = 2
    do while (i <= size(r))
       if (r(i) >= redge) exit
       rhomid = 0.5*(den(i)+den(i-1))
       wmid = (rhomid/rhomax)**q
       integral = integral + pi*(r(i)+r(i-1))**2*(r(i)-r(i-1))*wmid
       i = i + 1
    enddo
    rpre(1) = r(1)
    rpreedge = rpre(1)
    do i=2,size(r)
       rhomid = 0.5*(den(i)+den(i-1))
       wmid = (rhomid/rhomax)**q
       if (rpre(i-1) > 0.) then
          rpre(i) = rpre(i-1) + (r(i)-r(i-1))*wmid*(r(i-1)/rpre(i-1))**2
       else
          ! Central finite-density limit for analytic profiles starting at r=0.
          rpre(i) = r(i)*wmid**(1./3.)
       endif
       if (r(i) <= redge) rpreedge = rpre(i)
    enddo
 else
    integral = 4./3.*pi*redge**3
    rpre = r
    rpreedge = redge
 endif
 if (any(rpre(2:) <= rpre(:size(r)-1))) return
 if (integral <= 0.) return
 mcentral = rhomax*integral/n
 ! A point core is a distinct physical model; do not silently omit it.
 if (4.*mcentral >= mstar) then
    ierr = 2
    return
 endif

 call natural_spline(rpre,r,rmap2)
 cellvolume = 2.*integral/n
 a1 = (cellvolume/sqrt(2.))**(1./3.)
 ixmax = int(rpreedge/a1) + 2
 iymax = int(rpreedge/(sqrt(3.)/2.*a1)) + 2
 izmax = int(rpreedge/(0.5*sqrt(8./3.)*a1)) + 2

 ! Count and construct in the same order, including antipodal pairs.
 do pass=1,2
    ip = 0
    do ix=-ixmax,ixmax
       do iy=-iymax,iymax
          do iz=1,izmax
             xtry = (real(ix)-0.5)*a1 + mod(abs(iy),2)*0.5*a1
             ytry = iy*sqrt(3.)/2.*a1 - (mod(abs(iz),2)-0.5)/sqrt(3.)*a1
             ztry = (iz-0.5)*0.5*sqrt(8./3.)*a1
             rtry = sqrt(xtry**2+ytry**2+ztry**2)
             ! Never extrapolate the inverse map into spurious accepted sites.
             if (rtry > rpre(size(r))) cycle
             rpost = spline_value(rpre,r,rmap2,rtry)
             if (rpost > 0. .and. rpost < redge) then
                if (pass == 2) then
                   xyz(:,ip+1) = (/xtry,ytry,ztry/)*rpost/rtry
                   xyz(:,ip+2) = -xyz(:,ip+1)
                endif
                ip = ip + 2
             endif
          enddo
       enddo
    enddo
    if (pass == 1) then
       naccepted = ip
       if (naccepted == 0) return
       allocate(xyz(3,naccepted),h(naccepted),pmass(naccepted),rhoi(naccepted))
    endif
 enddo

 amtot = 0.
 do i=1,naccepted
    ri = sqrt(sum(xyz(:,i)**2))
    rhoi(i) = spline_value(r,den,rho2,ri)
    if (rhoi(i) <= 0.) then
       ierr = 3
       return
    endif
    pmass(i) = integral/naccepted*rhoi(i)**(1.-q)*rhomax**q
    amtot = amtot + pmass(i)
 enddo
 do i=1,naccepted
    pmass(i) = pmass(i)/amtot*mstar
    anumden = rhoi(i)/pmass(i)
    h(i) = (3./32./pi*1.9*nnopt/anumden)**(1./3.) + hfloor
 enddo
 ierr = 0

end subroutine generate_stellar_hcp

subroutine natural_spline(x,y,y2)
 real, intent(in) :: x(:),y(:)
 real, intent(out) :: y2(:)
 real :: u(size(x)),sig,p
 integer :: i,k,n

 n = size(x)
 y2(1) = 0.
 u(1) = 0.
 do i=2,n-1
    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    p = sig*y2(i-1) + 2.
    y2(i) = (sig-1.)/p
    u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i)) - &
            (y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 enddo
 y2(n) = 0.
 do k=n-1,1,-1
    y2(k) = y2(k)*y2(k+1) + u(k)
 enddo

end subroutine natural_spline

pure real function spline_value(x,y,y2,xi) result(yi)
 real, intent(in) :: x(:),y(:),y2(:),xi
 real :: a,b,h
 integer :: lo,hi,k

 lo = 1
 hi = size(x)
 do while (hi-lo > 1)
    k = (hi+lo)/2
    if (x(k) > xi) then
       hi = k
    else
       lo = k
    endif
 enddo
 h = x(hi)-x(lo)
 a = (x(hi)-xi)/h
 b = (xi-x(lo))/h
 yi = a*y(lo)+b*y(hi)+((a**3-a)*y2(lo)+(b**3-b)*y2(hi))*h**2/6.

end function spline_value

end module stellar_hcp
