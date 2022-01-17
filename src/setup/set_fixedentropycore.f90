!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
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
subroutine set_fixedS_softened_core(mcore,rcore,rho,r,pres,m,Xcore,Ycore,ierr)
 use eos,         only:ieos
 use physcon,     only:pi,gg,solarm,solarr
 use table_utils, only:interpolator
 use io,          only:fatal
 real, intent(inout)  :: r(:),rho(:),m(:),pres(:),mcore
 real, allocatable    :: r_alloc(:),rho_alloc(:),pres_alloc(:)
 real, intent(in)     :: rcore,Xcore,Ycore
 integer, intent(out) :: ierr
 real                 :: mc,msoft,rc
 integer              :: icore

 ierr = 0
 rc = rcore*solarr  ! convert to cm
 mc = mcore*solarm  ! convert to g
 call interpolator(r,rc,icore)  ! find index in r closest to rc
 msoft = m(icore) - mc
 if (msoft<0.) call fatal('setup','mcore cannot exceed m(r=h)')

 select case(ieos)
 case(2)
    ientropy = 1
 case(10,12,20)
    ientropy = 2
 case default
    call fatal('setfixedentropycore',&
               'ieos not one of 2 (adiabatic), 12 (ideal plus rad.), 10 (MESA), or 20 (gas+rad+recombination)')
 end select

 ! Make allocatable copies, see instructions of calc_rho_and_pres
 allocate(r_alloc(0:icore+1))
 r_alloc(0) = 0.
 r_alloc(1:icore+1) = r(1:icore+1)
 allocate(rho_alloc(0:icore))
 rho_alloc(icore) = rho(icore)
 allocate(pres_alloc(0:icore+1))
 pres_alloc(icore:icore+1) = pres(icore:icore+1)
 call calc_rho_and_pres(r_alloc,mc,m(icore),rho_alloc,pres_alloc,Xcore,Ycore)
 mcore = mc / solarm
 write(*,'(1x,a,f12.5,a)') 'Obtained core mass of ',mcore,' Msun'
 write(*,'(1x,a,f12.5,a)') 'Softened mass is ',m(icore)/solarm-mcore,' Msun'
 rho(1:icore)  = rho_alloc(1:icore)
 pres(1:icore) = pres_alloc(1:icore)
 call calc_mass_from_rho(r(1:icore),rho(1:icore),m(1:icore))
 m(1:icore) = m(1:icore) + mc

end subroutine set_fixedS_softened_core


!-----------------------------------------------------------------------
!+
!  Returns softened core profile with fixed entropy
!+
!-----------------------------------------------------------------------
subroutine calc_rho_and_pres(r,mcore,mh,rho,pres,Xcore,Ycore)
 use eos, only:entropy,get_mean_molecular_weight
 real, allocatable, dimension(:), intent(in)    :: r
 real, intent(in)                               :: mh,Xcore,Ycore
 real, intent(inout)                            :: mcore
 real, allocatable, dimension(:), intent(inout) :: rho,pres
 integer                                        :: Nmax
 real                                           :: Sc,mass,mold,msoft,fac,mu

! INSTRUCTIONS

! Input variables should be given in the following format:

! r(0:Nmax+1): Array of radial grid to be softened, satisfying r(0)=0 and r(Nmax)=rcore
! mcore:       Core particle mass, need to provide initial guess
! mh:          Mass coordinate at rcore, m(r=rcore)
! rho(0:Nmax): Give rho(Nmax)=(rho at rcore) as input. Outputs density profile.
! p(0:Nmax+1): Give p(Nmax:Nmax+1)=(p at r(Nmax:Nmax+1)) as input. Outputs pressure profile.

 msoft = mh - mcore
 Nmax  = size(rho)-1 ! Index corresponding to r = h
 mu = get_mean_molecular_weight(Xcore,1.-Xcore-Ycore)
 Sc = entropy(rho(Nmax),pres(Nmax),mu,ientropy)

 ! Start shooting method
 fac  = 0.05
 mass = msoft

 do
    mold = mass
    call one_shot(Sc,r,mcore,msoft,mu,rho,pres,mass) ! returned mass is m(r=0)
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
subroutine one_shot(Sc,r,mcore,msoft,mu,rho,pres,mass)
 use physcon, only:gg,pi
 use eos, only:get_rho_from_p_s
 real, intent(in)                               :: Sc,mcore,msoft,mu
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
    rhoguess = rho(i)
    call get_rho_from_p_s(pres(i-1),Sc,rho(i-1),mu,rhoguess,ientropy)
    mass = mass - 0.5*(rho(i)+rho(i-1)) * dvol(i)
    if (mass < 0.) return ! m(r) < 0 encountered, exit and decrease mcore
 enddo

 return

end subroutine one_shot

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
