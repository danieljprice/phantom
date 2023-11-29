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
 integer, parameter :: ierr_rho=1,ierr_pres=2,ierr_mass=3,ierr_lum=4

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

 iverbose = 1
 call shoot_for_mcore(eos_type,r_alloc,mc,m(icore),Lstar,rho_alloc,pres_alloc,T_alloc,Xcore,Ycore,iverbose)
 mcore = mc / solarm
 write(*,'(1x,a,es24.16e3,a)') 'Obtained core mass of ',mcore,' Msun'
 write(*,'(1x,a,es24.16e3,a)') 'Softened mass is ',m(icore)/solarm-mcore,' Msun'
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
 integer                                           :: Nmax,it_m,it_l,ierr
 real                                              :: mass,mold,msoft,fac_m,fac_l,mu,mcore_old,&
                                                      eps0,epsold,l,lold,tol_eps,tol_m

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
 fac_m  = 0.005
 mass = msoft
 tol_eps = 1.e-10
 tol_m = 1.e-10

 !---------------------------LOOP-OVER-MCORE-------------------------------------
 ! Vary mcore so that m(0) = 0
 it_m = 0
 loop_over_mcore: do
    l = Lstar
    eps0 = Lstar/msoft  ! initial guess for eps0 / erg/g
    mold = mass
    mcore_old = mcore
    !---------------------------LOOP-OVER-HEATING-FACTOR-------------------------------------
    ! Vary heating factor (eps0) so that central luminosity is zero
    it_l = 0
    fac_l = 0.01
    loop_over_eps0: do
       epsold = eps0
       lold = l
       ierr = 0
       call one_shot(eos_type,r,mcore,msoft,Lstar,eps0,mu,rho,pres,temp,mass,l,iverbose,ierr) ! returned mass is m(r=0)
       it_l = it_l + 1

       if (iverbose > 1) write(*,'(2(1x,i5),8(2x,a,e15.8),2x,a,i1)') &
                            it_m,it_l,'eps0=',eps0,'m(0)=',mass/solarm,&
                            'l(0)=',l,'eps0_old=',epsold,'mcore_old = ',&
                            mcore_old/solarm,'mcore=',mcore/solarm,'fac_l=',fac_l,&
                            'fac_m=',fac_m,'ierr=',ierr

       if (l < 0.) then
          eps0 = eps0 * (1. - fac_l)
       elseif (l/Lstar < tol_eps) then ! l(r=0) sufficiently close to zero
         !  write(*,'(/,1x,a,i5,a,e12.5)') 'Tolerance on luminosity reached on iteration no.',it_l,', fac_l =',fac_l
          exit loop_over_eps0
       else
          eps0 = eps0 * (1. + fac_l)
       endif

       if (abs(epsold-eps0) < tiny(0.)) then
          fac_l = fac_l * 1.05
       elseif (lold * l < 0.) then
          fac_l = fac_l * 0.95
       endif
    enddo loop_over_eps0
    !-----------------------------------------------------------------------------------------
    it_m = it_m + 1
    if (iverbose == 1) write(*,'(2(1x,i5),6(2x,a,e15.8),2x,a,i1)') &
                            it_m,it_l,'eps0=',eps0,'m(0)=',mass/solarm,'mcore_old = ',&
                            mcore_old/solarm,'mcore=',mcore/solarm,'fac_l=',fac_l,&
                            'fac_m=',fac_m,'ierr=',ierr

    if (mass < 0.) then
       mcore = mcore * (1. - fac_m)
    elseif (mass/msoft < tol_m .and. ierr <= ierr_pres) then  ! m(r=0) sufficiently close to zero
       write(*,'(/,1x,a,i5,2(1x,a,es24.16e3))') 'Converged on iteration no.',it_m,', fac_m =',fac_m,'eps0 = ',eps0
       if (ierr == ierr_rho) write(*,'(a)') 'WARNING: Profile contains density inversion'
       exit loop_over_mcore
    else
       mcore = mcore * (1. + fac_m)
    endif
    msoft = mh - mcore

    if (mold * mass < 0.) then
      fac_m = fac_m * 0.95
    else
      fac_m = fac_m * 1.05
    endif
   !  if (abs(mold-mass) < tiny(0.)) then
   !     fac_m = fac_m * 1.02
   !  elseif (mold * mass < 0.) then
   !     fac_m = fac_m * 0.99
   !  endif

    if (abs(mold-mass) < tiny(0.) .and. ierr <= ierr_rho) then
       write(*,'(/,1x,a,e12.5)') 'WARNING: Converged on mcore without reaching tolerance on zero &
                                 &central mass. m(r=0)/msoft = ',mass/msoft
       if (ierr == ierr_rho) write(*,'(1x,a)') 'WARNING: Profile contains density inversion'
       write(*,'(/,1x,a,i4,a,e12.5)') 'Reached iteration ',it_m,', fac_m=',fac_m
       exit loop_over_mcore
    endif

 enddo loop_over_mcore
 !-----------------------------------------------------------------------------------------

end subroutine shoot_for_mcore


!-----------------------------------------------------------------------
!+
!  One shot: Solve structure for given guess for msoft/mcore
!+
!-----------------------------------------------------------------------
subroutine one_shot(eos_type,r,mcore,msoft,Lstar,eps0,mu,rho,pres,T,mass,l,iverbose,ierr)
 use physcon,             only:gg,pi,radconst,c,solarm
 use eos,                 only:calc_rho_from_PT,iopacity_type
 use radiation_utils,     only:get_opacity
 use setfixedentropycore, only:gcore
 use units,               only:unit_density,unit_opacity
 integer, intent(in)                            :: eos_type,iverbose
 real, intent(in)                               :: mcore,msoft,Lstar,eps0,mu
 real, allocatable, dimension(:), intent(in)    :: r
 real, allocatable, dimension(:), intent(inout) :: rho,pres,T
 real, intent(out)                              :: mass,l
 integer, intent(out)                           :: ierr
 integer                                        :: i,Nmax
 real                                           :: kappai,kappa_code,rcore,mu_local,rho_code
 real, allocatable, dimension(:)                :: dr,dvol,lum

 Nmax = size(rho)-1
 allocate(dr(1:Nmax+1),dvol(1:Nmax+1),lum(1:Nmax+1))

 ! Pre-fill arrays
 do i = 1,Nmax+1
    dr(i) = r(i)-r(i-1)
    dvol(i) = 4.*pi/3. * (r(i)**3 - r(i-1)**3)
 enddo

 rcore = r(Nmax)
 mass  = msoft
 lum(Nmax:Nmax+1) = Lstar
 mu_local = mu
 ierr = 0

 do i = Nmax, 1, -1
    pres(i-1) = pres(i) + dr(i) * rho(i) * gg * (mass/r(i)**2 + mcore * gcore(r(i),rcore))
   !  pres(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
   !              * rho(i) * gg * (mass/r(i)**2 + mcore * gcore(r(i),rcore)) &
   !              + dr(i)**2 * pres(i+1) &
   !              + ( dr(i+1)**2 - dr(i)**2) * pres(i) ) / dr(i+1)**2
    rho_code = rho(i) / unit_density
    call get_opacity(iopacity_type,rho_code,T(i),kappa_code)
    kappai = kappa_code * unit_opacity
    T(i-1) = T(i) + dr(i) * 3./(16.*pi*radconst*c) * rho(i)*kappai*lum(i) / (r(i)**2*T(i)**3)
   !  T(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
   !           * 3./(16.*pi*radconst*c) * rho(i)*kappai*lum(i) / (r(i)**2*T(i)**3) &
   !           + dr(i)**2 * T(i+1) &
   !           + ( dr(i+1)**2 - dr(i)**2) * T(i) ) / dr(i+1)**2
    call calc_rho_from_PT(eos_type,pres(i-1),T(i-1),rho(i-1),ierr,mu_local)
    mass = mass - rho(i)*dvol(i)
    lum(i-1) = lum(i) - dr(i)*4.*pi*r(i)**2*rho(i)*eps0*eps_heating(r(i),rcore)
   !  lum(i-1) = ( dr(i) * dr(i+1) * sum(dr(i:i+1)) &
   !           * (-4.*pi)*r(i)**2*rho(i)*eps0*eps_heating(r(i),rcore) &
   !           + dr(i)**2 * lum(i+1) &
   !           + ( dr(i+1)**2 - dr(i)**2) * lum(i) ) / dr(i+1)**2
    l = lum(i-1)
    
    if (iverbose > 3) print*,Nmax-i+1,rho(i-1),mass,pres(i-1),T(i-1),kappai,lum(i-1)
    if (mass < 0.) then ! m(r) < 0 encountered, exit and decrease mcore
       if (iverbose > 2) print*,'WARNING: Negative mass reached at i = ',i, 'm = ',mass/solarm
       ierr = ierr_mass
       return
    endif
    if (l < 0.) then ! l(r) < 0 encountered, exit and increase heating pre-factor eps0
       if (iverbose > 2) print*,'WARNING: Negative luminosity reached at i = ',i, 'm = ',mass/solarm
       ierr = ierr_lum
       return
    endif
    if (rho(i-1)<rho(i)) then
       if (iverbose > 2) then 
          print*,'WARNING: Density inversion at i = ',i, 'm = ',mass/solarm
          write(*,'(i5,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.7)') i,rho(i),rho(i-1),mass,kappai
       endif
       ierr = ierr_rho
    endif
    if (pres(i-1)<pres(i)) then
       if (iverbose > 2) then
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
!  Wrapper function to return heating rate per unit mass (epsilon).
!  Warning: normalisation is arbitrary (see heating_kernel)
!+
!-----------------------------------------------------------------------
function eps_heating(r,rcore)
 use kernel,         only:radkern2
 use ptmass_heating, only:heating_kernel,isink_heating
 real, intent(in) :: r,rcore
 real             :: eps_heating,q2

 isink_heating = 1
 q2 = radkern2 * r**2 / rcore**2
 eps_heating = heating_kernel(q2,isink_heating)

end function eps_heating


end module setfixedlumcore
