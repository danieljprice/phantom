!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_idealplusrad
!
! Ideal gas equation of state plus radiation pressure, assumes
!               inputs are in cgs units
!
! :References: Stellar Structure and Evolution (2nd Edition) (Kippenhahn,
!              Weigert, Weiss)
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 use physcon,  only:Rg,radconst,mass_proton_cgs,kboltz
 implicit none
 real, parameter :: tolerance = 1.e-15

 public :: get_idealplusrad_temp,get_idealplusrad_pres,get_idealplusrad_spsoundi,&
           get_idealgasplusrad_tempfrompres,get_idealplusrad_tempfromrhoS,&
           get_idealplusrad_enfromtemp,get_idealplusrad_rhofrompresT,&
           egas_from_rhoT,erad_from_rhoT

 private

contains

!----------------------------------------------------------------
!+
!  Solve for temperature as a function of (gas+rad) internal energy
!  per unit mass (eni) and density (rhoi)
!+
!----------------------------------------------------------------
subroutine get_idealplusrad_temp(rhoi,eni,mu,tempi,ierr)
 real,    intent(in)    :: rhoi,eni,mu
 real,    intent(inout) :: tempi
 integer, intent(out)   :: ierr
 real                :: gasfac,imu,numerator,denominator,correction
 integer             :: iter
 integer, parameter  :: iter_max = 1000

 gasfac = 1.5 !this is NOT gamma = cp/cv, it refers to the gas being monoatomic
 imu = 1./mu
 if (tempi <= 0. .or. isnan(tempi)) tempi = eni*mu/(gasfac*Rg)  ! Take gas temperature as initial guess

 ierr = 0
 iter = 0
 correction = huge(0.)
 do while (abs(correction) > tolerance*tempi .and. iter < iter_max)
    numerator = eni*rhoi - gasfac*Rg*tempi*rhoi*imu - radconst*tempi**4
    denominator =  - gasfac*Rg*imu*rhoi - 4.*radconst*tempi**3
    correction = numerator/denominator
    tempi= tempi - correction
    iter = iter + 1
 enddo
 if (iter >= iter_max) ierr = 1

end subroutine get_idealplusrad_temp

subroutine get_idealplusrad_pres(rhoi,tempi,mu,presi)
 real, intent(in)  :: rhoi,tempi,mu
 real, intent(out) :: presi

 presi = (Rg*rhoi/mu + radconst*tempi**3/3.)*tempi ! Eq 13.2 (Kippenhahn et al.)

end subroutine get_idealplusrad_pres

subroutine get_idealplusrad_spsoundi(rhoi,presi,eni,spsoundi,gammai)
 real, intent(in)  :: rhoi,presi,eni
 real, intent(out) :: spsoundi,gammai

 gammai = 1. + presi/(eni*rhoi)
 spsoundi = sqrt(gammai*presi/rhoi)

end subroutine get_idealplusrad_spsoundi

!----------------------------------------------------------------
!+
!  Calculates temperature from pressure and density
!+
!----------------------------------------------------------------
subroutine get_idealgasplusrad_tempfrompres(presi,rhoi,mu,tempi)
 real, intent(in)    :: rhoi,presi,mu
 real, intent(inout) :: tempi
 real                :: imu,numerator,denominator,correction,temp_new
 integer             :: iter
 integer, parameter  :: iter_max = 1000

 iter = 0
 correction = huge(0.)
 imu = 1./mu
 tempi = min((3.*presi/radconst)**0.25, presi*mu/(rhoi*Rg))
 do while (abs(correction) > tolerance*tempi .and. iter < iter_max)
    numerator   = presi - rhoi*Rg*tempi*imu - radconst*tempi**4 /3.
    denominator =  - rhoi*Rg*imu - 4./3.*radconst*tempi**3
    correction  = numerator/denominator
    temp_new = tempi - correction
    if (temp_new > 1.2 * tempi) then
       tempi = 1.2 * tempi
    elseif (temp_new < 0.8 * tempi) then
       tempi = 0.8 * tempi
    else
       tempi = temp_new
    endif
    iter = iter + 1
 enddo

end subroutine get_idealgasplusrad_tempfrompres

!----------------------------------------------------------------
!+
!  Entropy residual and derivative for ideal gas + radiation
!+
!----------------------------------------------------------------
pure subroutine entropy_fdf(tempi,inv_mu_mh,log_rho,coeff_rad,cgss,fi,dfi,t2)
 real, intent(in)  :: tempi,inv_mu_mh,log_rho,coeff_rad,cgss
 real, intent(out) :: fi,dfi,t2

 t2 = tempi*tempi
 fi  = inv_mu_mh*(1.5*log(tempi) - log_rho) + coeff_rad*t2*tempi - cgss
 dfi = 1.5*inv_mu_mh/tempi + 3.*coeff_rad*t2

end subroutine entropy_fdf

!----------------------------------------------------------------
!+
!  Solve for temperature and pressure from density and entropy
!  (ideal gas + radiation, cgs units; entropy in erg/g/K)
!+
!----------------------------------------------------------------
subroutine get_idealplusrad_tempfromrhoS(rho,s,mu,temp,pres,niter_out)
 real,    intent(in)    :: rho,s,mu
 real,    intent(inout) :: temp
 real,    intent(out)   :: pres
 integer, intent(out), optional :: niter_out
 real                :: corr,corr_ic,df,f,f_ic,temp_new,temp_gas,temp_rad,temp_ic
 real                :: inv_mu_mh,coeff_rad,log_rho,cgss_tol,best_f,t2,term
 real, parameter     :: eoserr = 1.e-13
 real, parameter     :: one_third = 1./3.
 integer             :: niter
 integer, parameter  :: nitermax = 1000

 niter = 0
 inv_mu_mh = 1./(mu*mass_proton_cgs)
 coeff_rad = 4.*radconst/(3.*rho*kboltz)
 log_rho   = log(rho)
 cgss_tol  = eoserr*abs(s)
 best_f    = huge(best_f)
 temp_ic   = 0.
 corr_ic   = 0.
 f_ic      = best_f

 ! below we try three possible guesses to accelerate convergence
 ! first guess is the input temperature
 if (temp > 0.) then
    call entropy_fdf(temp,inv_mu_mh,log_rho,coeff_rad,s,f,df,t2)
    corr = f/df
    if (abs(corr) <= eoserr*temp .and. abs(f) <= cgss_tol) then
       pres = Rg*rho/mu*temp + radconst*t2*t2*one_third
       if (present(niter_out)) niter_out = niter
       return
    endif
    temp_ic = temp; corr_ic = corr; f_ic = f
    best_f  = abs(f)
 endif

 ! second guess is the radiation temperature (radiation-dominated limit, s > 0)
 ! we take this if the residual is lower than the best so far
 if (s > 0.) then
    temp_rad = (3.*s*rho*kboltz/(4.*radconst))**one_third
    if (temp_rad > 0.) then
       call entropy_fdf(temp_rad,inv_mu_mh,log_rho,coeff_rad,s,f,df,t2)
       if (abs(f) < best_f) then
          temp_ic = temp_rad; corr_ic = f/df; f_ic = f
          best_f  = abs(f)
       endif
    endif
 endif

 ! third guess is the gas temperature (skip expensive exp if IC already good)
 if (best_f > cgss_tol .or. temp_ic <= 0.) then
    term = mu*s*mass_proton_cgs
    if (term < 100.) then
       temp_gas = (rho*exp(mu*s*mass_proton_cgs))**(2./3.)
       call entropy_fdf(temp_gas,inv_mu_mh,log_rho,coeff_rad,s,f,df,t2)
       if (abs(f) < best_f) then
          temp_ic = temp_gas; corr_ic = f/df; f_ic = f
       endif
    endif
 endif

 temp = temp_ic
 t2   = temp*temp
 f    = f_ic
 corr = corr_ic
 do while ((abs(corr) > eoserr*temp .or. abs(f) > cgss_tol) .and. niter < nitermax)
    temp_new = temp - corr
    if (temp_new > 1.2*temp) then
       temp = 1.2*temp
    elseif (temp_new < 0.8*temp) then
       temp = 0.8*temp
    else
       temp = temp_new
    endif
    niter = niter + 1
    call entropy_fdf(temp,inv_mu_mh,log_rho,coeff_rad,s,f,df,t2)
    corr = f/df
 enddo

 pres = Rg*rho*temp/mu + radconst*t2*t2*one_third
 if (present(niter_out)) niter_out = niter

end subroutine get_idealplusrad_tempfromrhoS

!----------------------------------------------------------------
!+
!  Calculates internal energy per unit mass from density
!  and temperature
!+
!----------------------------------------------------------------
subroutine get_idealplusrad_enfromtemp(densi,tempi,mu,eni)
 real, intent(in)  :: densi,tempi,mu
 real, intent(out) :: eni

 eni = egas_from_rhoT(tempi,mu) + erad_from_rhoT(densi,tempi)

end subroutine get_idealplusrad_enfromtemp

!----------------------------------------------------------------
!+
!  Calculates specific gas energy from density and temperature
!+
!----------------------------------------------------------------
real function egas_from_rhoT(tempi,mu) result(egasi)
 real, intent(in) :: tempi,mu

 egasi = 1.5*Rg*tempi/mu

end function egas_from_rhoT

!----------------------------------------------------------------
!+
!  Calculates specific radiation energy from density and temperature
!+
!----------------------------------------------------------------
real function erad_from_rhoT(densi,tempi) result(eradi)
 real, intent(in) :: densi,tempi

 eradi = radconst*tempi**4/densi

end function erad_from_rhoT

!----------------------------------------------------------------
!+
!  Calculates density from pressure and temperature
!+
!----------------------------------------------------------------
subroutine get_idealplusrad_rhofrompresT(presi,tempi,mu,densi)
 real, intent(in)  :: presi,tempi,mu
 real, intent(out) :: densi

 densi = (presi - radconst*tempi**4 /3.) * mu / (Rg*tempi)

end subroutine get_idealplusrad_rhofrompresT

end module eos_idealplusrad
