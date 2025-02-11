!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setcubiccore
!
! This module softens the core of a MESA stellar profile with a cubic
!   density profile, given a core radius core mass, in preparation
!   for adding a sink particle core.
!
! :References: For cubic spline softening of sink gravity, see Price &
!              Monaghan (2007)
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: io, kernel, physcon, table_utils
!
 use physcon,     only:solarm,solarr
 use table_utils, only:interpolator,diff

 implicit none

 ! msoft: Softened mass (mass at softening length minus mass of core
 !        particle)
 ! hsoft: Softening length for the point particle potential, defined in
 !        Price & Monaghan (2007). Set to be 0.5*rcore.

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the cubic core profile
!+
!-----------------------------------------------------------------------
subroutine set_cubic_core(mcore,rcore,rho,r,pres,m)
 use io,          only:fatal
 real, intent(inout):: r(:),rho(:),m(:),pres(:)
 real, allocatable  :: phi(:)
 real, intent(in)   :: mcore,rcore
 real               :: mc,rc,hsoft_cm,msoft
 integer            :: icore,i

 rc       = rcore * solarr      ! Convert to cm
 hsoft_cm = 0.5*rc              ! Convert to cm
 mc       = mcore * solarm      ! Convert to g
 call interpolator(r,rc,icore)  ! Find index in r closest to rc
 msoft    = m(icore) - mc

 call calc_rho_and_m(rho, m, r, mc, rc)
 call calc_phi(r, mc, m-mc, hsoft_cm, phi)

 ! Calculate pressure
 do i = icore,2,-1
    pres(i-1) = pres(i) + rho(i) * (phi(i) - phi(i-1))
 enddo

end subroutine set_cubic_core


!----------------------------------------------------------------
!+
!  Iteratively look for a value of mcore for given rcore that
!  produces a nice softened density profile
!+
!----------------------------------------------------------------
subroutine find_mcore_given_rcore(rcore,r,rho0,m0,mcore,ierr)
 real,    intent(in)  :: r(:),rho0(:),m0(:),rcore
 real,    intent(out) :: mcore
 integer, intent(out) :: ierr
 real,    allocatable :: rho(:),drho(:),m(:)
 real                 :: rc,mc,tolerance
 integer              :: icore,counter=0

 rc = rcore * solarr ! Convert to cm
 call interpolator(r,rc,icore) ! Find index in r closest to rcore
 mc = 0.7*m0(icore) ! Initialise profile to have very large softened mass
 tolerance = 1.3 ! How much we allow the softened density to exceed the original profile by
 ierr = 0

 allocate(rho(size(rho0)))
 allocate(m(size(m0)))
 allocate(drho(size(rho0)-1))
 do
    rho = rho0   ! Reset density
    m   = m0     ! Reset mass
    call calc_rho_and_m(rho, m, r, mc, rc)
    call diff(rho,drho)
    if (all(rho/rho0 < tolerance) .and. all(drho(1:icore) < 0)) exit
    if (mc > 0.999*m0(icore)) then
       ierr = 1
       exit
    endif
    mc = mc + 0.0001*m0(icore) ! Increase mcore/m(rc) by 1 percent
    counter = counter+1
 enddo
 ! Output mcore in solar masses
 mcore = mc / solarm
end subroutine find_mcore_given_rcore

!----------------------------------------------------------------
!+
!  Iteratively look for a value of rcore for given mcore that
!  produces a nice softened density profile
!+
!----------------------------------------------------------------
subroutine find_rcore_given_mcore(mcore,r,rho0,m0,rcore,ierr)
 real,    intent(in)  :: r(:),rho0(:),m0(:),mcore
 real,    intent(out) :: rcore
 integer, intent(out) :: ierr
 real,    allocatable :: rho(:),drho(:),m(:)
 real                 :: rc,mc,mh,tolerance
 integer              :: icore

 mc = mcore * solarm ! Convert to g
 mh = mc / 0.7 ! Initialise rcore such that m(rcore) to be much larger than mcore
 tolerance = 1.3 ! How much we allow the softened density to exceed the original profile by
 ierr = 0
 allocate(rho(size(rho0)))
 allocate(m(size(m0)))
 allocate(drho(size(rho0)-1))
 do
    call interpolator(m0, mh, icore)
    rc   = r(icore)
    rho = rho0   ! Reset density
    m   = m0     ! Reset mass
    call calc_rho_and_m(rho, m, r, mc, rc)
    call diff(rho,drho)
    if (all(rho/rho0 < tolerance) .and. all(drho(1:icore) < 0)) exit
    if (mc > 0.98*m0(icore)) then
       ierr = 1
       exit
    endif
    call interpolator(m0, 1.3*mc, icore)
    mh = 1./(1./mh + 0.01/mc) ! Increase mcore/m(h) by 1 percent
 enddo
 ! Write out hsoft in solar radii
 rcore = rc / solarr
end subroutine find_rcore_given_mcore

!----------------------------------------------------------------
!+
!  Check for sensible values of rcore and mcore
!+
!----------------------------------------------------------------
subroutine check_rcore_and_mcore(rcore,mcore,r,rho0,m0,ierr)
 real,    intent(in)  :: r(:),rho0(:),m0(:),mcore,rcore
 integer, intent(out) :: ierr
 real,    allocatable :: rho(:),drho(:),m(:)
 real                 :: rc,mc,msoft,tolerance
 integer              :: icore
 ierr = 0
 tolerance = 1.5 ! How much we allow the softened density to exceed the original profile by
 rc = rcore * solarr ! Convert to cm
 mc = mcore * solarm ! Convert to g
 call interpolator(r,rc,icore) ! Find index in r closest to h
 msoft = m0(icore) - mc

 if (msoft < 0.) ierr = 1

 allocate(rho(size(rho0)))
 allocate(m(size(m0)))
 allocate(drho(size(rho0)-1))
 rho = rho0
 m = m0
 call calc_rho_and_m(rho, m, r, mc, rc)
 if (any(rho/rho0 > tolerance)) ierr = 2

 call diff(rho, drho)
 if (any(drho(1:icore) > 0)) ierr = 3
end subroutine check_rcore_and_mcore


subroutine calc_rho_and_m(rho,m,r,mc,rc)
 use physcon, only:pi
 real, intent(in)    :: r(:)
 real, intent(inout) :: rho(:),m(:)
 real, intent(in)    :: mc,rc
 real                :: a,b,d,msoft,drhodr_h
 integer             :: icore

 call interpolator(r,rc,icore) ! Find index in r closest to rc
 msoft    = m(icore) - mc
 drhodr_h = (rho(icore+1) - rho(icore)) / (r(icore+1) - r(icore)) ! drho/dr at r = rc

 ! a, b, d: Coefficients of cubic density profile defined by rho(r) = ar**3 + br**2 + d
 a = 2./rc**2 * drhodr_h - 10./rc**3 * rho(icore) + 7.5/pi/rc**6 * msoft
 b = 0.5*drhodr_h/rc - 1.5*a*rc
 d = rho(icore) - 0.5*rc*drhodr_h + 0.5*a*rc**3

 rho(1:icore) = a*r(1:icore)**3 + b*r(1:icore)**2 + d

 ! Mass is then given by m(r) = mcore + 4*pi (1/6 a r^6 + 1/5 b r^5 + 1/3 d r^3)
 m(1:icore) = mc + 4.*pi * (1./6. * a * r(1:icore)**6 + 0.2 * b * r(1:icore)**5 + &
                           1./3. * d * r(1:icore)**3)
end subroutine calc_rho_and_m


subroutine calc_phi(r,mc,mgas,hsoft,phi)
 use kernel, only:kernel_softening
 use physcon, only:gg
 real, intent(in)               :: r(:),mgas(:),mc,hsoft
 real, allocatable              :: q(:),q2(:),phi_core(:),phi_gas(:)
 real, allocatable, intent(out) :: phi(:)
 real                           :: dum
 integer                        :: i

 ! The gravitational potential is needed to integrate the pressure profile using the
 ! equation of hydrostatic equilibrium. First calculate gravitational potential due
 ! to point mass core, using softening kernel selected in Phantom. Then calculate the
 ! gravitational potential due to the softened gas.

 ! Gravitational potential due to primary core
 allocate(phi(size(r)), phi_core(size(r)), phi_gas(size(r)), q(size(r)), q2(size(r)))
 q  = r / hsoft
 q2 = q**2
 do i = 1,size(r)
    call kernel_softening(q2(i),q(i),phi_core(i),dum)
 enddo
 deallocate(q,q2)
 phi_core = gg * phi_core * mc / hsoft

 ! Gravitational potential due to softened gas
 phi_gas(size(r)) = - gg * mgas(size(r)) / r(size(r)) ! Surface boundary condition for phi
 do i = 1,size(r)-1
    phi_gas(size(r)-i) = phi_gas(size(r)-i+1) - gg * mgas(size(r)-i) / r(size(r)-i)**2. &
                                               * (r(size(r)-i+1) - r(size(r)-i))
 enddo

 phi = phi_gas + phi_core

end subroutine calc_phi

end module setcubiccore
