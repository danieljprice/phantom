!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setsoftenedcore
!
!  DESCRIPTION:
!   This module softens the core of a MESA stellar profile with a cubic
!   density profile, given a softening length and core mass, in preparation
!   for adding a sink particle core.
!
!  REFERENCES: For cubic spline softening of sink gravity, see Price &
!              Monaghan (2007)
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos_idealplusrad, physcon
!+
!--------------------------------------------------------------------------
module setsoftenedcore
 use physcon,          only:pi,gg,solarm,solarr,kb_on_mh
 use eos_idealplusrad, only:get_idealgasplusrad_tempfrompres,&
                            get_idealplusrad_enfromtemp

 implicit none
 real    :: hsoft,msoft,mcore
 integer :: hidx

 ! hsoft: Radius below which we replace the original profile with a
 !        softened profile. Note: This is called 'hdens' in
 !        setup_star.f90
 ! mcore: Mass of core particle
 ! msoft: Softened mass (mass at softening length minus mass of core
 !        particle)
 ! hphi:  Softening length for the point particle potential, defined in
 !        Price & Monaghan (2007). Set to be 0.5*hsoft. Note: This is
 !        called 'hsoft' in setup_star.f90

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calls the subroutines to calculate the cubic
!  core profile
!+
!-----------------------------------------------------------------------
subroutine set_softened_core(ieos,gamma,constX,constY,mu,mcore,hsoft,hphi,&
                             rho,r,pres,m,ene,temp,ierr,Xfrac,Yfrac)
 real, intent(inout)        :: r(:),rho(:),m(:),pres(:),ene(:),temp(:)
 real, allocatable          :: phi(:)
 real, intent(in), optional :: Xfrac(:),Yfrac(:)
 real, intent(in)           :: mu,mcore,hsoft,hphi,gamma,constX,constY
 integer, intent(in)        :: ieos
 integer, intent(out)       :: ierr
 real                       :: mc,h,hphi_cm
 logical                    :: isort_decreasing,iexclude_core_mass

 ! Xfrac, Yfrac:    H and He mass fractions from the MESA profile
 ! constX, constY:  H and He mass fractions at r = R/2
 ! mu:              Mean molecular weight at r = R/2
 ! gamma:           Adiabatic exponent for adiabatic EoS

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true.     ! Needs to be true if to be read by Phantom
 !
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true.   ! Needs to be true if to be read by Phantom

 h       = hsoft * solarr      ! Convert to cm
 hphi_cm = hphi  * solarr      ! Convert to cm
 mc      = mcore * solarm      ! Convert to g
 call interpolator(r,h,hidx)   ! Find index in r closest to h
 msoft = m(hidx) - mc

 call calc_rho_and_m(rho, m, r, mc, h)
 call calc_phi(r, mc, m-mc, hphi_cm, phi)
 call calc_pres(r, rho, phi, pres)

 if (present(Xfrac) .and. present(Yfrac)) then ! This case is temporarily forbidden (until MESA EoS with variable composition is implemented)
    call calc_temp_and_ene(ieos,hidx,mu,constX,constY,gamma,rho,pres,ene,temp,ierr,Xfrac,Yfrac)
 else
    call calc_temp_and_ene(ieos,hidx,mu,constX,constY,gamma,rho,pres,ene,temp,ierr)
 endif

 ! Reverse arrays so that data is sorted from stellar surface to stellar centre.
 if (isort_decreasing) then
    call flip_array(m)
    call flip_array(pres)
    call flip_array(temp)
    call flip_array(r)
    call flip_array(rho)
    call flip_array(ene)
    call flip_array(phi)
 endif

 if (iexclude_core_mass) then
    m = m - mc
 endif

end subroutine set_softened_core


!----------------------------------------------------------------
!+
!  Iteratively look for a value of mcore for given hsoft that
!  produces a nice softened density profile
!+
!----------------------------------------------------------------
subroutine find_mcore_given_hsoft(hsoft,r,rho0,m0,mcore,ierr)
 real,    intent(in)  :: r(:),rho0(:),m0(:),hsoft
 real,    intent(out) :: mcore
 integer, intent(out) :: ierr
 real,    allocatable :: rho(:),drho(:),m(:)
 real                 :: h,mc,tolerance
 integer              :: hidx,counter=0

 h = hsoft * solarr ! Convert to cm
 call interpolator(r, h, hidx) ! Find index in r closest to h
 mc = 0.7*m0(hidx) ! Initialise profile to have very large softened mass
 tolerance = 1.3 ! How much we allow the softened density to exceed the original profile by
 ierr = 0

 allocate(rho( size(rho0) ))
 allocate(m( size(m0) ))
 allocate(drho( size(rho0)-1 ))
 do
    rho = rho0   ! Reset density
    m   = m0     ! Reset mass
    call calc_rho_and_m(rho, m, r, mc, h)
    call diff(rho, drho)
    if (all(rho/rho0 < tolerance) .and. all(drho(1:hidx) < 0)) exit
    if (mc > 0.98*m0(hidx)) then
       ierr = 1
       exit
    endif
    mc = mc + 0.01*m0(hidx) ! Increase mcore/m(h) by 1 percent
    counter = counter+1
 enddo
 ! Output mcore in solar masses
 mcore = mc / solarm
end subroutine find_mcore_given_hsoft

!----------------------------------------------------------------
!+
!  Iteratively look for a value of hsoft for given mcore that
!  produces a nice softened density profile
!+
!----------------------------------------------------------------
subroutine find_hsoft_given_mcore(mcore,r,rho0,m0,hsoft,ierr)
 real,    intent(in)  :: r(:),rho0(:),m0(:),mcore
 real,    intent(out) :: hsoft
 integer, intent(out) :: ierr
 real,    allocatable :: rho(:),drho(:),m(:)
 real                 :: h,mc,mh,tolerance
 integer              :: hidx

 mc = mcore * solarm ! Convert to g
 mh = mc / 0.7 ! Initialise h such that m(h) to be much larger than mcore
 tolerance = 1.3 ! How much we allow the softened density to exceed the original profile by
 ierr = 0
 allocate(rho( size(rho0) ))
 allocate(m( size(m0) ))
 allocate(drho( size(rho0)-1 ))
 do
    call interpolator(m0, mh, hidx)
    h   = r(hidx)
    rho = rho0   ! Reset density
    m   = m0     ! Reset mass
    call calc_rho_and_m(rho, m, r, mc, h)
    call diff(rho, drho)
    if (all(rho/rho0 < tolerance) .and. all(drho(1:hidx) < 0)) exit
    if (mc > 0.98*m0(hidx)) then
       ierr = 1
       exit
    endif
    call interpolator(m0, 1.3*mc, hidx)
    mh = 1./(1./mh + 0.01/mc) ! Increase mcore/m(h) by 1 percent
 enddo
 ! Write out hsoft in solar radii
 hsoft = h / solarr
end subroutine find_hsoft_given_mcore

!----------------------------------------------------------------
!+
!  Check for sensible values of hsoft and mcore, returning errors.
!  Called in setup_star.f90 when both hsoft and mcore are chosen
!  by user
!+
!----------------------------------------------------------------
subroutine check_hsoft_and_mcore(hsoft,mcore,r,rho0,m0,ierr)
 real,    intent(in)  :: r(:),rho0(:),m0(:),mcore
 real,    intent(out) :: hsoft
 integer, intent(out) :: ierr
 real,    allocatable :: rho(:),drho(:),m(:)
 real                 :: h,mc,msoft,tolerance
 integer              :: hidx
 ierr = 0
 tolerance = 1.5 ! How much we allow the softened density to exceed the original profile by
 h = hsoft * solarr ! Convert to cm
 mc = mcore * solarm ! Convert to g
 call interpolator(r, h, hidx) ! Find index in r closest to h
 msoft = m0(hidx) - mc

 if (msoft < 0.) ierr = 1

 allocate(rho( size(rho0) ))
 allocate(m( size(m0) ))
 allocate(drho( size(rho0)-1 ))
 rho = rho0
 m   = m0
 call calc_rho_and_m(rho, m, r, mc, h)
 if (any(rho/rho0 > tolerance)) ierr = 2

 call diff(rho, drho)
 if (any(drho(1:hidx) > 0)) ierr = 3
end subroutine check_hsoft_and_mcore


subroutine calc_rho_and_m(rho,m,r,mc,h)
 real, intent(in)    :: r(:)
 real, intent(inout) :: rho(:),m(:)
 real, intent(in)    :: mc,h
 real                :: a,b,d,msoft,drhodr_h
 integer             :: hidx

 call interpolator(r,h,hidx) ! Find index in r closest to h
 msoft    = m(hidx) - mc
 drhodr_h = (rho(hidx+1) - rho(hidx)) / (r(hidx+1) - r(hidx)) ! drho/dr at r = h

 ! a, b, d: Coefficients of cubic density profile defined by rho(r) = ar**3 + br**2 + d
 a = 2./h**2 * drhodr_h - 10./h**3 * rho(hidx) + 7.5/pi/h**6 * msoft
 b = 0.5*drhodr_h/h - 1.5*a*h
 d = rho(hidx) - 0.5*h*drhodr_h + 0.5*a*h**3

 rho(1:hidx) = a*r(1:hidx)**3 + b*r(1:hidx)**2 + d

 ! Mass is then given by m(r) = mcore + 4*pi (1/6 a r^6 + 1/5 b r^5 + 1/3 d r^3)
 m(1:hidx) = mc + 4.*pi * (1./6. * a * r(1:hidx)**6 + 0.2 * b * r(1:hidx)**5 + &
                           1./3. * d * r(1:hidx)**3)
end subroutine calc_rho_and_m


subroutine calc_phi(r,mc,mgas,hphi,phi)
 real, intent(in)               :: r(:),mgas(:),mc,hphi
 real, allocatable              :: q(:),phi_core(:),phi_gas(:)
 real, allocatable, intent(out) :: phi(:)
 integer                        :: idx2hphi,idxhphi,i
 ! The gravitational potential is needed to integrate the pressure profile using the
 ! equation of hydrostatic equilibrium. First calculate gravitational potential due
 ! to point mass core, according to the cubic spline kernel in Price & Monaghan (2007)
 ! to be consistent with Phantom, then calculate the gravitational potential due to the
 ! softened gas.
 allocate(phi(size(r)), phi_core(size(r)), phi_gas(size(r)))
 ! (i) Gravitational potential due to core particle (cubic spline softening)
 ! For 0 <= r/hphi < 1
 call interpolator(r, hphi, idxhphi) ! Find index corresponding to r = 2*hphi
 allocate(q(1:idxhphi-1))
 q = r(1:idxhphi-1) / hphi
 phi_core(1:idxhphi-1) = gg*mc/hphi * (2./3.*q**2 - 0.3*q**4 + 0.1*q**5 &
                                   - 7./5.)
 deallocate(q)

 ! For 1 <= r/hphi < 2
 call interpolator(r, 2*hphi, idx2hphi) ! Find index corresponding to r = 2*hphi
 allocate(q(idxhphi:idx2hphi-1))
 q = r(idxhphi:idx2hphi-1) / hphi
 phi_core(idxhphi:idx2hphi-1) = gg*mc/hphi * (4./3.*q**2 - q**3 + 0.3*q**4 &
                                        - 1./30.*q**5 - 1.6 + 1./15./q)
 deallocate(q)

 ! For 2 <= r/hphi
 phi_core(idx2hphi:size(r)) = - gg * mc / r(idx2hphi:size(r))

 ! (ii) Gravitational potential due to softened gas
 phi_gas(size(r)) = - gg * mgas(size(r)) / r(size(r)) ! Surface boundary condition for phi
 do i = 1,size(r)-1
    phi_gas(size(r)-i) = phi_gas(size(r)-i+1) - gg * mgas(size(r)-i) / r(size(r)-i)**2. &
                                               * (r(size(r)-i+1) - r(size(r)-i))
 enddo

 ! (iii) Add the potentials
 phi = phi_gas + phi_core
end subroutine calc_phi


subroutine calc_pres(r, rho, phi, pres)
 ! Calculates pressure by integrating the equation of hydrostatic equilibrium
 ! given the gravitational potential and the density profile
 real, intent(in)  :: rho(:),phi(:),r(:)
 real, intent(out) :: pres(:)
 integer           :: i

 pres(size(r)) = 0 ! Set boundary condition of zero pressure at stellar surface
 do i = 1,size(r)-1
    ! Reverse Euler
    pres(size(r)-i) = pres(size(r)-i+1) + rho(size(r)-i+1) * (phi(size(r)-i+1) - phi(size(r)-i))
 enddo
end subroutine calc_pres


subroutine calc_temp_and_ene(ieos,hidx,mu,constX,constY,gamma,rho,pres,ene,temp,ierr,Xfrac,Yfrac)
 real, intent(in)           :: rho(:),pres(:)
 real, intent(in), optional :: Xfrac(:),Yfrac(:)
 real, intent(in)           :: mu,gamma,constX,constY
 real, intent(inout)        :: ene(:),temp(:)
 integer, intent(in)        :: ieos,hidx
 integer, intent(out)       :: ierr
 real                       :: e_rec,tempi,muprofile(size(rho)),Xprofile(size(rho)),Yprofile(size(rho))
 integer                    :: i,endidx

 ! Calculate temperature and internal energy profile depending on EoS
 ierr = 0
 if (ieos == 2) then ! Adiabatic/polytropic EoS
    temp = pres / (rho * kb_on_mh) * mu
    ene = pres / ( (gamma-1.) * rho)

 elseif ( (ieos == 10) .or. (ieos == 12) ) then ! MESA or ideal gas plus rad. EoS
    if (present(Xfrac) .and. present(Yfrac)) then ! This case is temporarily forbidden (until MESA EoS with variable composition is implemented)
       ! Calculate mean molecular weight profile from X, Y
       muprofile = 4. / (6.*Xfrac + Yfrac + 2.)
       Xprofile  = Xfrac
       Yprofile  = Yfrac
       endidx    = hidx ! For r > h, we just use the original MESA profile
    else
       ! Use the fixed mu, X, and Y calculated at r = R/2 for entire profile
       do i = 1,size(muprofile)
          muprofile(i) = mu
          Xprofile(i)  = constX
          Yprofile(i)  = constY
       enddo
       endidx = size(rho) ! For consistency, T and eint have to be recalculated for the entire star
    endif

    ! Calculate temperature using ideal gas plus radiation EoS
    if (ieos == 12) temp(size(rho)) = 0. ! Zero surface temperature
    do i = 1,endidx
       tempi = temp(i)
       call get_idealgasplusrad_tempfrompres(pres(i),rho(i),muprofile(i),tempi)
       temp(i) = tempi
    enddo

    ! Calculate internal energy
    do i = 1,endidx
       ! Contribution from pressure and gas
       call get_idealplusrad_enfromtemp(rho(i),temp(i),muprofile(i),ene(i))
       if (ieos == 10) then
          ! Contribution from ionisation energy for MESA EoS
          call calc_rec_ene(Xprofile(i),Yprofile(i),e_rec)
          ene(i) = ene(i) + e_rec
       endif
    enddo

 else
    ierr = 1
 endif

end subroutine calc_temp_and_ene

!----------------------------------------------------------------
!+
!  Get recombination energy (per unit mass) complete ionisation
!+
!----------------------------------------------------------------
subroutine calc_rec_ene(XX,YY,e_rec)
 real, intent(in)  :: XX, YY
 real, intent(out) :: e_rec
 real              :: e_H2,e_HI,e_HeI,e_HeII
 real, parameter   :: e_ion_H2   = 1.312d13, & ! ionisation energies in erg/mol
                      e_ion_HI   = 4.36d12, &
                      e_ion_HeI  = 2.3723d13, &
                      e_ion_HeII = 5.2505d13

 ! XX     : Hydrogen mass fraction
 ! YY     : Helium mass fraction
 ! e_rec  : Total ionisation energy due to H2, HI, HeI, and HeII

 e_H2   = 0.5 * XX * e_ion_H2
 e_HI   = XX * e_ion_HI
 e_HeI  = 0.25 * YY * e_ion_HeI
 e_HeII = 0.25 * YY * e_ion_HeII
 e_rec  = e_H2 + e_HI + e_HeI + e_HeII

end subroutine calc_rec_ene

!----------------------------------------------------------------
!+
!  Finds index of the array value closest to a given value for an
!  ordered array
!+
!----------------------------------------------------------------
subroutine interpolator(array, value, valueidx)
 real, intent(in)     :: array(:)
 real, intent(in)     :: value
 integer, intent(out) :: valueidx

 valueidx = minloc(abs(array - value), dim = 1)
end subroutine interpolator

!----------------------------------------------------------------
!+
!  Reverses the elements of a 1-d array
!+
!----------------------------------------------------------------
subroutine flip_array(array)
 real, intent(inout) :: array(:)
 real, allocatable   :: flipped_array(:)
 integer             :: i
 allocate(flipped_array(size(array)))
 do i = 1, size(array)
    flipped_array(i) = array(size(array) - i + 1)
 enddo
 array = flipped_array
end subroutine flip_array

!----------------------------------------------------------------
!+
!  Takes a n-dim array and produces a (n-1)-dim array with the
!  ith element being the (i+1)th element minus the ith element of
!  the original array
!+
!----------------------------------------------------------------
subroutine diff(array, darray)
 real, intent(in)               :: array(:)
 real, allocatable, intent(out) :: darray(:)
 integer                        :: i

 allocate(darray(size(array)-1))
 do i = 1, size(array)-1
    darray(i) = array(i+1) - array(i)
 enddo
end subroutine diff

end module setsoftenedcore
