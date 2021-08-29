!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module eos_idealplusrad
!
! Ideal gas equation of state plus radiation pressure, assumes
!               inputs are in cgs units
!
! :References: Stellar Structure and Evolution (2nd Edition) (Kippenhahn,
!              Weigert, Weiss)
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 use physcon,  only:kb_on_mh,radconst
 implicit none
 real, parameter :: tolerance = 1d-15

 public :: get_idealplusrad_temp,get_idealplusrad_pres,get_idealplusrad_spsoundi,&
           get_idealgasplusrad_tempfrompres,get_idealplusrad_enfromtemp

 private

contains

!----------------------------------------------------------------
!+
!  Solve for temperature as a function of internal energy per
!  unit mass (eni) and density (rhoi)
!+
!----------------------------------------------------------------
subroutine get_idealplusrad_temp(rhoi,eni,mu,tempi)
 real, intent(in)    :: rhoi,eni,mu
 real, intent(inout) :: tempi
 real                :: numerator,denominator,correction,temp_new
 integer             :: iter
 integer, parameter  :: iter_max = 1000

 if (tempi <= 0.) then
    tempi = eni*mu/(1.5*kb_on_mh)  ! Take gas temperature as initial guess
 endif

 iter = 0
 correction = huge(0.)
 do while (abs(correction) > tolerance*tempi .and. iter < iter_max)
    numerator = eni*rhoi - 1.5*kb_on_mh*tempi*rhoi/mu - radconst*tempi**4
    denominator =  - 1.5*kb_on_mh/mu*rhoi - 4.*radconst*tempi**3
    correction = numerator/denominator
    temp_new = tempi - correction
    if (temp_new > 1.2 * tempi) then
       tempi = 1.2 * tempi
    else if (temp_new < 0.8 * tempi) then
       tempi = 0.8 * tempi
    else
       tempi = temp_new
    endif
    iter = iter + 1
 enddo

end subroutine get_idealplusrad_temp


subroutine get_idealplusrad_pres(rhoi,tempi,mu,presi)
 real, intent(in)    :: rhoi,mu
 real, intent(inout) :: tempi
 real, intent(out)   :: presi

 presi = kb_on_mh*rhoi*tempi/mu + 1./3.*radconst*tempi**4 ! Eq 13.2 (Kippenhahn et al.)

end subroutine get_idealplusrad_pres


subroutine get_idealplusrad_spsoundi(rhoi,presi,eni,spsoundi)
 real, intent(in)  :: rhoi,presi,eni
 real, intent(out) :: spsoundi
 real              :: gamma

 gamma = 1. + presi/(eni*rhoi)
 spsoundi = sqrt(gamma*presi/rhoi)

end subroutine get_idealplusrad_spsoundi

!----------------------------------------------------------------
!+
!  Calculates temperature from pressure and density
!+
!----------------------------------------------------------------
subroutine get_idealgasplusrad_tempfrompres(presi,rhoi,mu,tempi)
 real, intent(in)    :: rhoi,presi,mu
 real, intent(inout) :: tempi
 real                :: numerator,denominator,correction,temp_new
 integer             :: iter
 integer, parameter  :: iter_max = 1000

 iter = 0
 correction = huge(0.)
 do while (abs(correction) > tolerance*tempi .and. iter < iter_max)
    numerator   = presi - rhoi*kb_on_mh*tempi/mu - radconst*tempi**4 /3.
    denominator =  - rhoi*kb_on_mh/mu - 4./3.*radconst*tempi**3
    correction  = numerator/denominator
    temp_new = tempi - correction
    if (temp_new > 1.2 * tempi) then
       tempi = 1.2 * tempi
    else if (temp_new < 0.8 * tempi) then
       tempi = 0.8 * tempi
    else
       tempi = temp_new
    endif
    iter = iter + 1
 enddo

end subroutine get_idealgasplusrad_tempfrompres


!----------------------------------------------------------------
!+
!  Calculates internal energy per unit mass from density
!  and temperature
!+
!----------------------------------------------------------------
subroutine get_idealplusrad_enfromtemp(densi,tempi,mu,eni)
 real, intent(in)  :: densi,tempi,mu
 real, intent(out) :: eni

 eni = 1.5*kb_on_mh*tempi/mu + radconst*tempi**4/densi

end subroutine get_idealplusrad_enfromtemp

end module eos_idealplusrad
