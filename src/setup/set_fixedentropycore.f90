!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module setfixedentropycore
!
! This module softens the core of a MESA stellar profile with a constant
!   entropy profile, given a core radius and mass, in preparation
!   for adding a sink particle core.
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, kernel, physcon, table_utils
!
 implicit none
 integer :: ientropy

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the constant entropy softened profile
!+
!-----------------------------------------------------------------------
subroutine set_fixedS_softened_core(mcore,rcore,rho,r,pres,m,ene,temp,ierr)
 use eos,         only:calc_temp_and_ene,ieos
 use physcon,     only:pi,gg,solarm,solarr,kb_on_mh
 use table_utils, only:interpolator,flip_array
 use io,          only:fatal
 real, intent(inout)        :: r(:),rho(:),m(:),pres(:),ene(:),temp(:),mcore
 real, allocatable          :: r_alloc(:),rho_alloc(:),pres_alloc(:)
 real, intent(in)           :: rcore
 integer, intent(out)       :: ierr
 real                       :: mc,msoft,rc,eneguess
 logical                    :: isort_decreasing,iexclude_core_mass
 integer                    :: i,icore

 ! Output data to be sorted from stellar surface to interior?
 isort_decreasing = .true.     ! Needs to be true if to be read by Phantom
 ! Exclude core mass in output mass coordinate?
 iexclude_core_mass = .true.   ! Needs to be true if to be read by Phantom

 rc = rcore * solarr     ! Convert to cm
 mc = mcore * solarm     ! Convert to g
 call interpolator(r,rc,icore)   ! Find index in r closest to rc
 msoft = m(icore) - mc

 select case(ieos)
 case(2)
    ientropy = 1
 case(12)
    ientropy = 2
 case(10)
    ientropy = 2
 case default
    call fatal('setfixedentropycore','ieos not one of 2 (adiabatic), 12 (ideal plus rad.), or 10 (MESA)')
 end select

 ! Make allocatable copies, see instructions of calc_rho_and_pres
 allocate(r_alloc(0:icore+1))
 r_alloc(0) = 0.
 r_alloc(1:icore+1) = r(1:icore+1)
 allocate(rho_alloc(0:icore))
 rho_alloc(icore) = rho(icore)
 allocate(pres_alloc(0:icore+1))
 pres_alloc(icore:icore+1) = pres(icore:icore+1)
 call calc_rho_and_pres(r_alloc,mc,m(icore),rho_alloc,pres_alloc)
 mcore = mc / solarm
 write(*,'(1x,a,f12.5,a)') 'Obtained core mass of ',mcore,' Msun'
 write(*,'(1x,a,f12.5,a)') 'Softened mass is ',m(icore)/solarm-mcore,' Msun'
 rho(1:icore)  = rho_alloc(1:icore)
 pres(1:icore) = pres_alloc(1:icore)

 call calc_mass_from_rho(r(1:icore),rho(1:icore),m(1:icore))
 m(1:icore) = m(1:icore) + mc

 call calc_temp_and_ene(rho(1),pres(1),ene(1),temp(1),ierr)
 if (ierr /= 0) call fatal('setfixedentropycore','EoS not one of: adiabatic, ideal gas plus radiation, MESA in set_softened_core')
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
subroutine calc_rho_and_pres(r,mcore,mh,rho,pres)
 use eos, only:entropy
 real, allocatable, dimension(:), intent(in)    :: r
 real, intent(in)                               :: mh
 real, intent(inout)                            :: mcore
 real, allocatable, dimension(:), intent(inout) :: rho,pres
 integer                                        :: Nmax
 real                                           :: Sc,mass,mold,msoft,fac

! INSTRUCTIONS

! Input variables should be given in the following format:

! r(0:Nmax+1): Array of radial grid to be softened, satisfying r(0)=0 and r(Nmax)=rcore
! mcore:       Core particle mass, need to provide initial guess
! mh:          Mass coordinate at rcore, m(r=rcore)
! rho(0:Nmax): Give rho(Nmax)=(rho at rcore) as input. Outputs density profile.
! p(0:Nmax+1): Give p(Nmax:Nmax+1)=(p at r(Nmax:Nmax+1)) as input. Outputs pressure profile.

 msoft = mh - mcore
 Nmax  = size(rho)-1 ! Index corresponding to r = h
 Sc = entropy(rho(Nmax),pres(Nmax),ientropy)

 ! Start shooting method
 fac  = 0.05
 mass = msoft

 do
    mold = mass
    call one_shot(Sc,r,mcore,msoft,rho,pres,mass) ! returned mass is m(r=0)
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
    if (mold == mass) then
       write(*,'(a,f12.5)') 'Warning: Setting fixed entropy for m(r=0)/msoft = ',mass/msoft
       exit
    endif
 enddo

 return

end subroutine calc_rho_and_pres


!-----------------------------------------------------------------------
!+
!  Calculate a hydrostatic structure for a given entropy
!+
!-----------------------------------------------------------------------
subroutine one_shot(Sc,r,mcore,msoft,rho,pres,mass)
 use physcon, only:gg,pi
 real, intent(in)                               :: Sc,mcore,msoft
 real, allocatable, dimension(:), intent(in)    :: r
 real, allocatable, dimension(:), intent(inout) :: rho,pres
 real, intent(out)                              :: mass
 integer                                        :: i,Nmax
 real                                           :: rcore,rhoguess
 real, allocatable, dimension(:)                :: dr,dvol

 Nmax = size(rho)-1
 allocate(dr(1:Nmax+1),dvol(1:Nmax+1))

 do i = 1,Nmax+1
    dr(i) = r(i)-r(i-1)
    dvol(i) = 4.*pi/3. * (r(i)**3 - r(i-1)**3)
 enddo

 rcore = r(Nmax)
 mass  = msoft

 do i = Nmax, 1, -1
    pres(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
                * rho(i) * gg * (mass/r(i)**2 + mcore * gcore(r(i),rcore)) &
                + dr(i)**2 * pres(i+1) &
                + ( dr(i+1)**2 - dr(i)**2) * pres(i) ) / dr(i+1)**2
    if (i == Nmax) then
       rhoguess = 1.e8
    else
       rhoguess = rho(i)
    endif
    call get_rho_from_p_s(pres(i-1),Sc,rho(i-1))
    mass = mass - 0.5*(rho(i)+rho(i-1)) * dvol(i)
    if (mass < 0.) return ! m(r) < 0 encountered, exit and decrease mcore
 enddo

 return

end subroutine one_shot


!-----------------------------------------------------------------------
!+
!  Calculate density given pressure and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_rho_from_p_s(pres,S,rho)
 use eos, only:entropy
 real, intent(in)  :: pres,S
 real, intent(out) :: rho
 real              :: corr,dSdrho,S_plus_dS,rho_plus_drho
 real, parameter   :: eoserr=1d-9,dfac=1d-12

 rho = 1d-8 ! Initial guess
 corr = huge(corr)

 do while (abs(corr) > eoserr*rho)
    ! First calculate dS/drho
    rho_plus_drho = rho * (1. + dfac)
    S_plus_dS = entropy(rho_plus_drho, pres, ientropy)
    dSdrho = (S_plus_dS - entropy(rho,pres,ientropy)) / (rho_plus_drho - rho)
    corr = ( entropy(rho,pres,ientropy) - S ) / dSdrho
    rho = rho - corr
 enddo

 return

end subroutine get_rho_from_p_s


!-----------------------------------------------------------------------
!+
!  Acceleration from softened potential of stellar core
!+
!-----------------------------------------------------------------------
function gcore(r,rcore)
 use kernel, only:kernel_softening
 real, intent(in) :: r,rcore
 real             :: gcore,dum
 real             :: hsoft,q

 hsoft = 0.5 * rcore
 q = r / hsoft
 call kernel_softening(q**2,q,dum,gcore)
 gcore = gcore / hsoft**2 ! Note: gcore is not multiplied by G or mcore yet.

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
