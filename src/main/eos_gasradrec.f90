!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module eos_gasradrec
!
!  EoS from HORMONE (Hirai+16) that includes internal energy from ideal gas,
!  radiation, and recombination (H2, H, He)
!
! :References: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.93.083006
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
 implicit none
 public :: equationofstate_gasradrec,calc_uT_from_rhoP_gasradrec
 private

contains
!-----------------------------------------------------------------------
!+
!  EoS from HORMONE (Hirai+16). Note eint is internal energy per unit volume
!+
!-----------------------------------------------------------------------
subroutine equationofstate_gasradrec(d,eint,T,imu,X,Y,p,cf)
! PURPOSE: To calculate pressure from density and internal energy without B field
 use ionization_mod, only:get_erec_imurec
 use physcon, only:radconst,kb_on_mh
 real,intent(in):: d,eint
 real,intent(inout):: T,imu ! imu is 1/mu, an output
 real,intent(in):: X,Y
 real,intent(out) :: p,cf
 real:: corr, erec, derecdT, dimurecdT, Tdot, logd, dt, gamma_eff
 real,parameter:: W4err = 1d-2, eoserr=1d-13
 integer n

 corr=1d99;Tdot=0d0;logd=log10(d);dt=0.9d0
 do n = 1, 500
    call get_erec_imurec(logd,T,X,Y,erec,imu,derecdT,dimurecdT)
    if(d*erec>=eint)then ! avoid negative thermal energy
       T = 0.9d0*T; Tdot=0d0;cycle
    end if
    corr = (eint-(radconst*T**3+1.5*kb_on_mh*d*imu)*T-d*erec) &
       / ( -4d0*radconst*T**3-d*(1.5*kb_on_mh*(imu+dimurecdT*T)+derecdT) )
    if(abs(corr)>W4err*T)then
       T = T + Tdot*dt
       Tdot = (1d0-2d0*dt)*Tdot - dt*corr
    else
       T = T-corr
       Tdot = 0d0
    end if
    if(abs(corr)<eoserr*T)exit
    if(n>50)dt=0.5d0
 end do
 if(n>500)then
    print*,'Error in eos_p'
    print*,'d=',d,'eint=',eint,'mu=',1d0/imu
    stop
 end if
 p = ( kb_on_mh*imu*d + radconst*T**3/3d0 )*T
 gamma_eff = 1d0+p/(eint-d*erec)
 cf = sqrt(gamma_eff*p/d)
end subroutine equationofstate_gasradrec
   
   
!-----------------------------------------------------------------------
!+
!  Solve for u and T (and mu) from rho and P (ideal gas + radiation + recombination)
!  Inputs and outputs in cgs units
!+
!-----------------------------------------------------------------------
subroutine calc_uT_from_rhoP_gasradrec(rhoi,presi,X,Y,T,eni,mui,ierr)
 use ionization_mod, only:get_imurec,get_erec_imurec
 use physcon,        only:radconst,kb_on_mh
 real, intent(in)            :: rhoi,presi,X,Y
 real, intent(inout)         :: T
 real, intent(out)           :: eni
 real, optional, intent(out) :: mui
 integer, intent(out)        :: ierr
 integer                     :: n
 real                        :: logrhoi,imu,dimurecdT,ereconrhoi,dT,Tdot,corr
 real, parameter:: W4err=1.e-2, eoserr=1.e-13

 ierr = 0
 corr = 1.e99
 Tdot = 0.
 dT = 0.9
 logrhoi = log10(rhoi)
 do n = 1,500
    call get_imurec(logrhoi,T,X,Y,imu,dimurecdT)
    corr = ( presi - (rhoi*kb_on_mh*imu + radconst*T**3/3.)*T ) / &
           ( -rhoi*kb_on_mh * ( imu + T*dimurecdT ) - 4.*radconst*T**3/3. )
    if (abs(corr)>W4err*T) then
       T = T + Tdot*dT
       Tdot = (1.-2.*dT)*Tdot - dT*corr
    else
       T = T - corr
       Tdot = 0.
    end if
    if (abs(corr) < eoserr*T) exit
    if (n > 50) dT = 0.5
 end do
 if (n > 500) then
    print*,'Error in calc_uT_from_rhoP_hormone'
    print*,'rho=',rhoi,'P=',presi,'mu=',1./imu
    ierr = 1
    return
 end if
 call get_erec_imurec(logrhoi,T,X,Y,ereconrhoi,imu)
 eni = ( 1.5*kb_on_mh*imu + radconst*T**3/rhoi )*T + rhoi*ereconrhoi
 mui = 1./imu

end subroutine calc_uT_from_rhoP_gasradrec

end module eos_gasradrec