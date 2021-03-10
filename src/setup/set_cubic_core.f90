!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
! :Dependencies: eos, io, kernel, physcon, table_utils
!
 use physcon,          only:pi,gg,solarm,solarr,kb_on_mh
 use table_utils,      only:interpolator,diff

 implicit none
 real    :: rcore,msoft,mcore
 integer :: icore

 ! rcore: Radius below which we replace the original profile with a
 !        softened profile.
 ! mcore: Mass of core particle
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
subroutine set_cubic_core(mcore,rcore,rho,r,pres,m,ene,temp,ierr)
 use eos,         only:calc_temp_and_ene
 use table_utils, only:flip_array
 use io,          only:fatal
 real, intent(inout)        :: r(:),rho(:),m(:),pres(:),ene(:),temp(:)
 real, allocatable          :: phi(:)
 real, intent(in)           :: mcore,rcore
 integer, intent(out)       :: ierr
 real                       :: mc,rc,hsoft_cm,eneguess
 logical                    :: isort_decreasing,iexclude_core_mass
 integer                    :: i

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true.     ! Needs to be true if to be read by Phantom
 !
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true.   ! Needs to be true if to be read by Phantom

 rc       = rcore * solarr      ! Convert to cm
 hsoft_cm = 0.5*rc             ! Convert to cm
 mc       = mcore * solarm      ! Convert to g
 call interpolator(r,rc,icore)   ! Find index in r closest to rc
 msoft = m(icore) - mc

 call calc_rho_and_m(rho, m, r, mc, rc)
 call calc_phi(r, mc, m-mc, hsoft_cm, phi)
 call calc_pres(r, rho, phi, pres)

 call calc_temp_and_ene(rho(1),pres(1),ene(1),temp(1),ierr)
 if (ierr /= 0) call fatal('set_cubic_core','EoS not one of: adiabatic, ideal gas plus radiation, MESA in set_softened_core')
 do i = 2,size(rho)-1
    eneguess = ene(i-1)
    call calc_temp_and_ene(rho(i),pres(i),ene(i),temp(i),ierr,eneguess)
 enddo
 ene(size(rho))  = 0. ! Zero surface internal energy
 temp(size(rho)) = 0. ! Zero surface temperature

 ! Reverse arrays so that data is sorted from stellar surface to stellar centre.
 if (isort_decreasing) then
    call flip_array(m)
    call flip_array(pres)
    call flip_array(temp)
    call flip_array(r)
    call flip_array(rho)
    call flip_array(ene)
 endif

 if (iexclude_core_mass) then
    m = m - mc
 endif

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
 call interpolator(r, rc, icore) ! Find index in r closest to rcore
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

!----------------------------------------------------------------
!+
!  Calculate pressure by integrating the equation of hydrostatic
!  equilibrium given the gravitational potential and the density
!  profile
!+
!----------------------------------------------------------------
subroutine calc_pres(r,rho,phi,pres)
 real, intent(in)  :: rho(:),phi(:),r(:)
 real, intent(out) :: pres(:)
 integer           :: i

 pres(size(r)) = 0 ! Set boundary condition of zero pressure at stellar surface
 do i = 1,size(r)-1
    ! Reverse Euler
    pres(size(r)-i) = pres(size(r)-i+1) + rho(size(r)-i+1) * (phi(size(r)-i+1) - phi(size(r)-i))
 enddo

end subroutine calc_pres

end module setcubiccore
