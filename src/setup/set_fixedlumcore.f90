!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setfixedlumcore
!
! This module softens the core of a MESA stellar profile with a specified
! temperature profile that is in equilibrium with a luminosity function (see
! options in the function "luminosity"). This assumes the softened region is
! radiative and so the temperature gradient provides the flux needed to
! transport this luminosity.
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: dim, eos, io, kernel, physcon, table_utils
!
 implicit none
 public :: set_fixedlum_softened_core

 private
 integer, parameter :: ierr_rho=1,ierr_pres=2,ierr_mass=3

contains

!-----------------------------------------------------------------------
!+
!  Main subroutine that calculates the softened profile
!  Lstar in erg/s
!+
!-----------------------------------------------------------------------
subroutine set_fixedlum_softened_core(eos_type,rcore,Lstar,mcore,rho,r,pres,m,Xcore,Ycore,ierr)
 use eos,                 only:calc_temp_and_ene,get_mean_molecular_weight,iopacity_type
 use io,                  only:fatal
 use physcon,             only:solarm,solarr
 use table_utils,         only:interpolator
 use setfixedentropycore, only:calc_mass_from_rho
 integer, intent(in)  :: eos_type
 real, intent(in)     :: rcore,Lstar,Xcore,Ycore
 real, intent(inout)  :: r(:),rho(:),m(:),pres(:),mcore
 real, allocatable    :: r_alloc(:),rho_alloc(:),pres_alloc(:),T_alloc(:)
 integer, intent(out) :: ierr
 real                 :: mc,msoft,rc,eni,mu
 integer              :: i,icore,iverbose

 ierr = 0
 if (Lstar<=tiny(0.)) then
    print *,'Lstar=',Lstar
    call fatal('setfixedlumcore','Lstar must be positive')
 endif
 if (iopacity_type/=1 .and. iopacity_type/=2) then
    print *,'iopacity_type=',iopacity_type
    call fatal('setfixedlumcore','only iopacity_type = 1,2 are supported')
 endif

 rc = rcore*solarr  ! convert to cm
 mc = mcore*solarm  ! convert to g
 call interpolator(r,rc,icore)  ! find index in r closest to rc
 msoft = m(icore) - mc
 if (msoft<0.) call fatal('setfixedlumcore','mcore cannot exceed m(r=h)')

 ! Make allocatable copies, see instructions of calc_rho_and_pres
 allocate(r_alloc(0:icore+1))
 r_alloc(0) = 0.
 r_alloc(1:icore+1) = r(1:icore+1)
 allocate(rho_alloc(0:icore))
 rho_alloc(icore) = rho(icore)
 allocate(pres_alloc(0:icore+1))
 pres_alloc(icore:icore+1) = pres(icore:icore+1)

 ! Allocate and fill in temperature array
 allocate(T_alloc(0:icore+1))
 do i = icore,icore+1
    mu = get_mean_molecular_weight(Xcore,1.-Xcore-Ycore)
    call calc_temp_and_ene(eos_type,rho(i),pres(i),eni,T_alloc(i),ierr,mu_local=mu, &
                           X_local=Xcore,Z_local=1.-Xcore-Ycore)
 enddo

 iverbose = 0
 call shoot_for_mcore(eos_type,r_alloc,mc,m(icore),Lstar,rho_alloc,pres_alloc,T_alloc,Xcore,Ycore,iverbose)
 mcore = mc / solarm
 write(*,'(1x,a,f8.5,a)') 'Obtained core mass of ',mcore,' Msun'
 write(*,'(1x,a,f8.5,a)') 'Softened mass is ',m(icore)/solarm-mcore,' Msun'
 rho(1:icore)  = rho_alloc(1:icore)
 pres(1:icore) = pres_alloc(1:icore)
 call calc_mass_from_rho(r(1:icore),rho(1:icore),m(1:icore))
 m(1:icore) = m(1:icore) + mc

end subroutine set_fixedlum_softened_core


!-----------------------------------------------------------------------
!+
!  Returns softened core profile
!+
!-----------------------------------------------------------------------
subroutine shoot_for_mcore(eos_type,r,mcore,mh,Lstar,rho,pres,temp,Xcore,Ycore,iverbose)
 use eos,     only:get_mean_molecular_weight
 use physcon, only:solarm
 integer, intent(in)                               :: eos_type,iverbose
 real,    allocatable, dimension(:), intent(in)    :: r
 real,    intent(in)                               :: Lstar,mh,Xcore,Ycore
 real,    intent(inout)                            :: mcore
 real,    allocatable, dimension(:), intent(inout) :: rho,pres,temp
 integer                                           :: Nmax,it,ierr
 real                                              :: mass,mold,msoft,fac,mu,mcore_old

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

 ! Start shooting method
 fac  = 0.0005
 mass = msoft
 it = 0
 do
    mold = mass
    mcore_old = mcore
    ierr = 0
    call one_shot(eos_type,r,mcore,msoft,Lstar,mu,rho,pres,temp,mass,iverbose,ierr) ! returned mass is m(r=0)
    it = it + 1

    if (iverbose > 0) write(*,'(1x,i5,4(2x,a,e15.8),2x,a,i1)') it,'m(r=0) = ',mass/solarm,'mcore_old = ',&
                            mcore_old/solarm,'mcore = ',mcore/solarm,'fac = ',fac,'ierr = ',ierr

    if (mass < 0.) then
       mcore = mcore * (1. - fac)
    elseif (mass/msoft < 1d-10 .and. ierr <= ierr_pres) then  ! m(r=0) sufficiently close to zero
       write(*,'(/,1x,a,i5,a,e12.5)') 'Tolerance on central mass reached on iteration no.',it,', fac =',fac
       if (ierr == ierr_rho) write(*,'(a)') 'WARNING: Profile contains density inversion'
       exit
    else
       mcore = mcore * (1. + fac)
    endif
    msoft = mh - mcore

    if (abs(mold-mass) < tiny(0.)) then
       fac = fac * 1.02
    elseif (mold * mass < 0.) then
       fac = fac * 0.99
    endif

    if (abs(mold-mass) < tiny(0.) .and. ierr <= ierr_rho) then
       write(*,'(/,1x,a,e12.5)') 'WARNING: Converged on mcore without reaching tolerance on zero &
                                 &central mass. m(r=0)/msoft = ',mass/msoft
       if (ierr == ierr_rho) write(*,'(1x,a)') 'WARNING: Profile contains density inversion'
       write(*,'(/,1x,a,i4,a,e12.5)') 'Reached iteration ',it,', fac=',fac
       exit
    endif

    
 enddo

end subroutine shoot_for_mcore


!-----------------------------------------------------------------------
!+
!  One shot: Solve structure for given guess for msoft/mcore
!+
!-----------------------------------------------------------------------
subroutine one_shot(eos_type,r,mcore,msoft,Lstar,mu,rho,pres,T,mass,iverbose,ierr)
 use physcon,             only:gg,pi,radconst,c,solarm
 use eos,                 only:calc_rho_from_PT,iopacity_type
 use radiation_utils,     only:get_opacity
 use setfixedentropycore, only:gcore
 use units,               only:unit_density,unit_opacity
 integer, intent(in)                            :: eos_type,iverbose
 real, intent(in)                               :: mcore,msoft,Lstar,mu
 real, allocatable, dimension(:), intent(in)    :: r
 real, allocatable, dimension(:), intent(inout) :: rho,pres,T
 real, intent(out)                              :: mass
 integer, intent(out)                           :: ierr
 integer                                        :: i,Nmax
 real                                           :: kappai,kappa_code,rcore,mu_local,rho_code
 real, allocatable, dimension(:)                :: dr,dvol,lum

 Nmax = size(rho)-1
 allocate(dr(1:Nmax+1),dvol(1:Nmax+1),lum(1:Nmax))

 ! Pre-fill arrays
 do i = 1,Nmax+1
    dr(i) = r(i)-r(i-1)
    dvol(i) = 4.*pi/3. * (r(i)**3 - r(i-1)**3)
 enddo

 rcore = r(Nmax)
 mass  = msoft
 lum(Nmax) = Lstar
 mu_local = mu
 ierr = 0

 do i = Nmax, 1, -1
    pres(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
                * rho(i) * gg * (mass/r(i)**2 + mcore * gcore(r(i),rcore)) &
                + dr(i)**2 * pres(i+1) &
                + ( dr(i+1)**2 - dr(i)**2) * pres(i) ) / dr(i+1)**2
    rho_code = rho(i) / unit_density
    call get_opacity(iopacity_type,rho_code,T(i),kappa_code)
    kappai = kappa_code * unit_opacity
    T(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
             * 3./(16.*pi*radconst*c) * rho(i)*kappai*lum(i) / (r(i)**2*T(i)**3) &
             + dr(i)**2 * T(i+1) &
             + ( dr(i+1)**2 - dr(i)**2) * T(i) ) / dr(i+1)**2
    call calc_rho_from_PT(eos_type,pres(i-1),T(i-1),rho(i-1),ierr,mu_local)
    mass = mass - rho(i)*dvol(i)
    lum(i-1) = luminosity(mass/msoft,Lstar)
    
    if (iverbose > 2) print*,Nmax-i+1,rho(i-1),mass,pres(i-1),T(i-1),kappai
    if (mass < 0.) then ! m(r) < 0 encountered, exit and decrease mcore
       if (iverbose > 1) print*,'WARNING: Negative mass reached at i = ',i, 'm = ',mass/solarm
       ierr = ierr_mass
       return
    endif
    if (rho(i-1)<rho(i)) then
       if (iverbose > 1) then 
          print*,'WARNING: Density inversion at i = ',i, 'm = ',mass/solarm
          write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.7)') i,rho(i),rho(i-1),mass,kappai
       endif
       ierr = ierr_rho
    endif
    if (pres(i-1)<pres(i)) then
       if (iverbose > 1) then
          print*,'WARNING: Pressure inversion at i = ',i, 'm = ',mass/solarm
          write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.7)') i,pres(i-1),rho(i),mass,kappai
       endif
       ierr = ierr_pres
       return
    endif
 enddo
end subroutine one_shot


!-----------------------------------------------------------------------
!+
!  Normalised luminosity function. q can be m/msoft or r/rcore
!  Note: For point mass heating, q = r/(radkern*hsoft), not r/hsoft, so
!        that luminosity reaches target value at r = radkern*hsoft
!+
!-----------------------------------------------------------------------
function luminosity(q,Lstar,hsoft)
!  use kernel, only:radkern,wkern,cnormk
 real, intent(in)           :: q,Lstar
 real, intent(in), optional :: hsoft
 real                       :: luminosity
 integer                    :: ilum

 ilum = 0

 if (q > 1) then
    luminosity = 1.
 else
    select case(ilum)
    case(1)  ! smooth step
       luminosity = 3.*q**2 - 2.*q**3
   !  case(2)  ! kernel softening
   !     r_on_hsoft = q*radkern
   !     luminosity = cnormk*wkern(r_on_hsoft*r_on_hsoft,r_on_hsoft)/hsoft**3
    case default  ! linear (constant heating rate)
       luminosity = q
    end select
   endif

 luminosity = luminosity * Lstar

end function luminosity


end module setfixedlumcore
