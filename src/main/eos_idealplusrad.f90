!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: eos_idealplusrad
!
!  DESCRIPTION: Ideal gas equation of state plus radiation pressure
!
!  REFERENCES:
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: physcon
!+
!--------------------------------------------------------------------------
module eos_idealplusrad
 use physcon,  only:kb_on_mh,radconst
 implicit none
 real, parameter :: tolerance = 1d-15

 public :: get_idealplusrad_temp,get_idealplusrad_ponrhoi,get_idealplusrad_spsoundi,&
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
 real                :: numerator,denominator,correction

 if (tempi <= 0.) then
    tempi = eni*mu/(1.5*kb_on_mh)  ! Take gas temperature as initial guess
 endif

 correction = huge(0.)
 do while (abs(correction) > tolerance*tempi)
    numerator = eni*rhoi - 1.5*kb_on_mh*tempi*rhoi/mu - radconst*tempi**4
    denominator =  - 1.5*kb_on_mh/mu*rhoi - 4.*radconst*tempi**3
    correction = numerator/denominator
    tempi = tempi - correction
 enddo
end subroutine get_idealplusrad_temp


subroutine get_idealplusrad_ponrhoi(rhoi,tempi,mu,ponrhoi)
 real, intent(in)    :: rhoi,mu
 real, intent(inout) :: tempi
 real, intent(out)   :: ponrhoi

 ponrhoi = kb_on_mh*tempi/mu + 1./3.*radconst*tempi**4/rhoi
end subroutine get_idealplusrad_ponrhoi


subroutine get_idealplusrad_spsoundi(ponrhoi,eni,spsoundi)
 real, intent(in)  :: ponrhoi,eni
 real, intent(out) :: spsoundi
 real              :: gamma

 gamma = 1. + ponrhoi/eni
 spsoundi = sqrt(gamma*ponrhoi)
end subroutine get_idealplusrad_spsoundi

!----------------------------------------------------------------
!+
!  Calculates temperature from pressure and density
!+
!----------------------------------------------------------------
subroutine get_idealgasplusrad_tempfrompres(presi,rhoi,mu,tempi)
 real, intent(in)    :: rhoi,presi,mu
 real, intent(inout) :: tempi
 real                :: numerator,denominator,correction

 do
    numerator = presi - rhoi*kb_on_mh*tempi/mu - radconst*tempi**4/3.
    denominator =  - rhoi*kb_on_mh/mu - 4./3.*radconst*tempi**3
    correction = numerator/denominator
    tempi = tempi - correction
    if (abs(correction) < tolerance*tempi) exit
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
