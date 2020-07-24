!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: setfixedentropycore
!
!  DESCRIPTION:
!   This module softens the core of a MESA stellar profile with a constant
!   entropy profile, given a softening length and core mass, in preparation
!   for adding a sink particle core.
!
!  REFERENCES:
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: eos, kernel, physcon, table_utils
!+
!--------------------------------------------------------------------------
module setfixedentropycore
 implicit none

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the constant entropy softened profile
!+
!-----------------------------------------------------------------------
subroutine set_fixedS_softened_core(mcore,hsoft,hphi,rho,r,pres,m,ene,temp,ierr)
 use eos,         only:calc_temp_and_ene
 use physcon,     only:pi,gg,solarm,solarr,kb_on_mh
 use table_utils, only:interpolator,flip_array
 real, intent(inout)        :: r(:),rho(:),m(:),pres(:),ene(:),temp(:),mcore
 real, allocatable          :: r_alloc(:),rho_alloc(:),pres_alloc(:)
 real, intent(in)           :: hsoft,hphi
 integer, intent(out)       :: ierr
 real                       :: mc,msoft,h,hphi_cm,eneguess
 logical                    :: isort_decreasing,iexclude_core_mass
 integer                    :: i,hidx,iSerr

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true.     ! Needs to be true if to be read by Phantom

 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true.   ! Needs to be true if to be read by Phantom

 h       = hsoft * solarr      ! Convert to cm
 hphi_cm = hphi  * solarr      ! Convert to cm
 mc      = mcore * solarm      ! Convert to g
 call interpolator(r,h,hidx)   ! Find index in r closest to h
 msoft = m(hidx) - mc

 ! Make allocatable copies, see instructions of calc_rho_and_pres
 allocate(r_alloc(0:hidx+1))
 r_alloc(0) = 0.
 r_alloc(1:hidx+1) = r(1:hidx+1)
 allocate(rho_alloc(0:hidx))
 rho_alloc(hidx) = rho(hidx)
 allocate(pres_alloc(0:hidx+1))
 pres_alloc(hidx:hidx+1) = pres(hidx:hidx+1)

 call calc_rho_and_pres(r_alloc,mc,m(hidx),rho_alloc,pres_alloc,iSerr)
 if (iSerr == 1) ierr = 2
 mcore = mc / solarm
 write(*,'(1x,a,f12.5,a)') 'Obtained core mass of ',mcore,' Msun'
 rho(1:hidx)  = rho_alloc(1:hidx)
 pres(1:hidx) = pres_alloc(1:hidx)

 call calc_mass_from_rho(r(1:hidx),rho(1:hidx),m(1:hidx))
 m(1:hidx) = m(1:hidx) + mc

 call calc_temp_and_ene(rho(1),pres(1),ene(1),temp(1),ierr)
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

end subroutine set_fixedS_softened_core


!-----------------------------------------------------------------------
!+
!  Returns softened core profile with fixed entropy
!+
!-----------------------------------------------------------------------
subroutine calc_rho_and_pres(r,mcore,mh,rho,pres,ierr)
 real, allocatable, dimension(:), intent(in)    :: r
 real, intent(in)                               :: mh
 real, intent(inout)                            :: mcore
 real, allocatable, dimension(:), intent(inout) :: rho,pres
 integer, intent(out)                           :: ierr
 integer                                        :: Nmax
 real                                           :: Sc,mass,mold,msoft,fac,Sedge

! Instructions

! input variables should be given in the following format

! r(0:Nmax+1): Array of radial grid to be softened, satisfying r(0)=0 and r(Nmax)=hsoft
! mcore:       Core particle mass, need to provide initial guess
! mh:          Mass coordinate at hsoft, m(r=h)
! rho(0:Nmax): Give rho(Nmax)=(rho at hsoft) as input. Outputs density profile.
! p(0:Nmax+1): Give p(Nmax:Nmax+1)=(p at r(Nmax:Nmax+1)) as input. Outputs pressure profile.

! ierr: Is set to 1 when the fixed entropy value exceeds the outer entropy.
!       Should choose a smaller core mass if this happens.

 ierr  = 0
 msoft = mh - mcore
 Nmax  = size(rho)-1 ! Index corresponding to r = h
 Sedge = entropy(rho(Nmax),pres(Nmax))

 ! Start shooting method
 fac  = 0.05
 mass = msoft
 Sc   = Sedge

 do
    mold = mass
    call one_shot(Sc,r,mcore,msoft,rho,pres,mass)
    if (mass < 0.) then
       mcore = mcore * (1. - fac)
       msoft = mh - mcore
    elseif (mass/msoft < 1d-10) then
       exit ! Happy when m(r=0) is sufficiently close to zero
    else
       mcore = mcore * (1. + fac)
       msoft = mh - mcore
    endif

    if (mold * mass < 0.) fac = fac * 0.5
 end do

 if (Sedge < Sc) ierr = 1
 return

end subroutine calc_rho_and_pres


!-----------------------------------------------------------------------
!+
!  Calculate a hydrostatic structure for a given entropy
!+
!-----------------------------------------------------------------------
subroutine one_shot(Sc,r,mcore,msoft,rho,pres,mass)
 use physcon, only: gg,pi
 real, intent(in)                               :: Sc,mcore,msoft
 real, allocatable, dimension(:), intent(in)    :: r
 real, allocatable, dimension(:), intent(inout) :: rho,pres
 real, intent(out)                              :: mass
 integer                                        :: i,Nmax
 real                                           :: hsoft
 real, allocatable, dimension(:)                :: dr,dvol

 Nmax = size(rho)-1
 allocate(dr(1:Nmax+1),dvol(1:Nmax+1))

 do i = 1,Nmax+1
    dr(i) = r(i)-r(i-1)
    dvol(i) = 4.*pi/3. * (r(i)**3 - r(i-1)**3)
 end do

 hsoft = r(Nmax)
 mass  = msoft

 do i = Nmax, 1, -1
    pres(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
                * rho(i) * gg * (mass/r(i)**2 + mcore * gcore(r(i),hsoft)) &
                + dr(i)**2 * pres(i+1) &
                + ( dr(i+1)**2 - dr(i)**2) * pres(i) ) / dr(i+1)**2
    call get_rho_from_p_s(pres(i-1),Sc,rho(i-1))
    mass = mass - 0.5*(rho(i)+rho(i-1)) * dvol(i)
    if (mass < 0.) return ! Choice of entropy fails to give m(r=0) = 0
 end do

 return

end subroutine one_shot


!-----------------------------------------------------------------------
!+
!  Calculates mass-specific entropy (gas + radiation) up to an additive
!  integration constant, from density and pressure.
!+
!-----------------------------------------------------------------------
function entropy(rho,pres)
 use physcon, only:radconst,kb_on_mh
 use eos,     only:gmw,ieos
 real, intent(in) :: rho,pres
 real :: inv_mu,entropy,temp
 real :: corr
 real, parameter :: eoserr=1d-10

 inv_mu = 1/gmw
 corr = 1d99; temp = 1d3

 ! First solve for temperature given density and pressure using Newton-
 ! Raphson method, assumming ideal gas plus radiation EoS
 do while (abs(corr) > eoserr*temp)
    corr = (pres - (radconst*temp**3/3.+ rho*kb_on_mh*inv_mu)* temp) &
          / (-4.*radconst*temp**3/3. - rho*kb_on_mh*inv_mu)
    temp = temp - corr
 end do

 if (ieos == 2) then
    ! Include only gas entropy for adiabatic EoS
    entropy = kb_on_mh * inv_mu * log(temp**1.5/rho)
 else
    ! Include both gas and radiation entropy for MESA and gas plus rad. EoSs
    entropy = kb_on_mh * inv_mu * log(temp**1.5/rho) + 4.*radconst*temp**3 / (3.*rho)
 endif

end function entropy


!-----------------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_rho_from_p_s(pres,S,rho)
 real, intent(in)  :: pres,S
 real, intent(out) :: rho
 real              :: corr,dSdrho,S_plus_dS,rho_plus_drho
 real, parameter   :: eoserr=1d-9,dfac=1d-12

 rho = 1d-8 ! Initial guess
 corr = 1d99

 do while (abs(corr) > eoserr*rho)
    ! First calculate dS/drho
    rho_plus_drho = rho * (1. + dfac)
    S_plus_dS = entropy(rho_plus_drho, pres)
    dSdrho = (S_plus_dS - entropy(rho,pres)) / (rho_plus_drho - rho)
    corr = ( entropy(rho,pres) - S ) / dSdrho
    rho = rho - corr
 end do

return

end subroutine get_rho_from_p_s


!-----------------------------------------------------------------------
!+
!  Acceleration from softened potential of stellar core
!+
!-----------------------------------------------------------------------
function gcore(r,hsoft)
 use kernel, only:kernel_softening
 real, intent(in) :: r,hsoft
 real             :: gcore,dum
 real             :: hphi,q

 hphi = 0.5 * hsoft
 q = r / hphi
 call kernel_softening(q**2,q,dum,gcore)
 gcore = gcore / hphi**2 ! Note: gcore is not multiplied by G or mcore yet.

end function gcore


subroutine calc_mass_from_rho(r,rho,m)
 use physcon, only:pi
 real, intent(in)    :: r(:),rho(:)
 real, intent(inout) :: m(:)
 real                :: densi
 integer             :: i

 m(1) = 0
 do i = 1,size(r)-1
    densi = 0.5*(rho(i) + rho(i+1))
    m(i+1) = m(i) + 4./3. * pi * (r(i+1)**3 - r(i)**3) * densi
 enddo

end subroutine calc_mass_from_rho

end module setfixedentropycore