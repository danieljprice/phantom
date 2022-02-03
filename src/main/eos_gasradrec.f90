!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module eos_gasradrec
!
! EoS from HORMONE that includes internal energy from ideal gas,
!  radiation, and recombination (H2, H, He) (Hirai et al., 2020)
!
! :References: Appendix C, https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.1154H/abstract
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: ionization_mod, physcon
!
 implicit none
 public :: equationofstate_gasradrec,calc_uT_from_rhoP_gasradrec
 private

contains
!-----------------------------------------------------------------------
!+
!  EoS from HORMONE (Hirai et al., 2020). Note eint is internal energy per unit volume
!+
!-----------------------------------------------------------------------
subroutine equationofstate_gasradrec(d,eint,T,imu,X,Y,p,cf)
 use ionization_mod, only:get_erec_imurec
 use physcon,        only:radconst,Rg
 real, intent(in)    :: d,eint
 real, intent(inout) :: T,imu ! imu is 1/mu, an output
 real, intent(in)    :: X,Y
 real, intent(out)   :: p,cf
 real                :: corr,erec,derecdT,dimurecdT,Tdot,logd,dt,gamma_eff
 real, parameter     :: W4err=1.e-2,eoserr=1.e-13
 integer n

 corr=huge(0.); Tdot=0.; logd=log10(d); dt=0.9
 do n = 1,500
    call get_erec_imurec(logd,T,X,Y,erec,imu,derecdT,dimurecdT)
    if (d*erec>=eint) then ! avoid negative thermal energy
       T = 0.9*T; Tdot=0.;cycle
    endif
    corr = (eint-(radconst*T**3+1.5*Rg*d*imu)*T-d*erec) &
           / ( -4.*radconst*T**3-d*(1.5*Rg*(imu+dimurecdT*T)+derecdT) )
    if (abs(corr) > W4err*T) then
       T = T + Tdot*dt
       Tdot = (1.-2.*dt)*Tdot - dt*corr
    else
       T = T-corr
       Tdot = 0.
    endif
    if (abs(corr)<eoserr*T) exit
    if (n>50) dt=0.5
 enddo
 if (n > 500) then
    print*,'Error in equationofstate_gasradrec'
    print*,'d=',d,'eint=',eint,'mu=',1./imu
    stop
 endif
 p = ( Rg*imu*d + radconst*T**3/3. )*T
 gamma_eff = 1.+p/(eint-d*erec)
 cf = sqrt(gamma_eff*p/d)

end subroutine equationofstate_gasradrec


!-----------------------------------------------------------------------
!+
!  Solve for u and T (and mu) from rho and P (ideal gas + radiation + recombination)
!  Inputs and outputs in cgs units
!+
!-----------------------------------------------------------------------
subroutine calc_uT_from_rhoP_gasradrec(rhoi,presi,X,Y,T,eni,mui,ierr)
 use ionization_mod, only: get_imurec,get_erec
 use physcon,        only: radconst,Rg
 real, intent(in)            :: rhoi,presi,X,Y
 real, intent(inout)         :: T
 real, intent(out)           :: eni
 real, optional, intent(out) :: mui
 integer, intent(out)        :: ierr
 integer                     :: n
 real                        :: logrhoi,imu,dimurecdT,dT,Tdot,corr
 real, parameter             :: W4err=1.e-2,eoserr=1.e-13

 if (T <= 0.) T = min((3.*presi/radconst)**0.25, presi/(rhoi*Rg))  ! initial guess for temperature
 ierr = 0
 corr = huge(0.)
 Tdot = 0.
 dT = 0.9
 logrhoi = log10(rhoi)
 do n = 1,500
    call get_imurec(logrhoi,T,X,Y,imu,dimurecdT)
    corr = ( presi - (rhoi*Rg*imu + radconst*T**3/3.)*T ) / &
           ( -rhoi*Rg * ( imu + T*dimurecdT ) - 4.*radconst*T**3/3. )
    if (abs(corr)>W4err*T) then
       T = T + Tdot*dT
       Tdot = (1.-2.*dT)*Tdot - dT*corr
    else
       T = T - corr
       Tdot = 0.
    endif
    if (abs(corr) < eoserr*T) exit
    if (n > 50) dT = 0.5
 enddo
 if (n > 500) then
    print*,'Error in calc_uT_from_rhoP_hormone'
    print*,'rho=',rhoi,'P=',presi,'mu=',1./imu
    ierr = 1
    return
 endif
 eni = ( 1.5*Rg*imu + radconst*T**3/rhoi )*T + get_erec(logrhoi,T,X,Y)
 mui = 1./imu

end subroutine calc_uT_from_rhoP_gasradrec

end module eos_gasradrec
