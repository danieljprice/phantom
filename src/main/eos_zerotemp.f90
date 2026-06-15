!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_zerotemp
!
! Implements zero temperature equation of state, e.g. for white dwarfs. Meant to 
!
! :References: Kippenhahn & Weigert, Stellar Structure and Evolution, section 15.2
!
! :Dependencies: infile_utils, io, units, physcon
!
 use units, only:unit_density,unit_velocity
 use physcon,  only:pi,atomic_mass_unit, mass_electron_cgs, planckh, c, 
 implicit none
 real, parameter :: tolerance = 1.e-15

 public :: get_idealplusrad_temp,get_idealplusrad_pres,get_idealplusrad_spsoundi,&
           get_idealgasplusrad_tempfrompres,get_idealplusrad_enfromtemp,&
           get_idealplusrad_rhofrompresT,egas_from_rhoT,erad_from_rhoT

 private

contains

! ----------------------------------------------------------------
!+
!  Calculates the zero temperature pressure for a fully degenerate electron gas, i.e. a
!  white dwarf. See Kippenhahn & Weigert, Stellar Structure and Evolution, section 15.2
!  Note that this is only the electron degeneracy pressure, so does not include the ion
!+
! ----------------------------------------------------------------

subroutine get_zerotemp_Pressure(rhoi,mu,presi)
 real, intent(in)  :: rhoi,mu
 real, intent(out) :: presi

 ! This is the zero temperature pressure for a fully degenerate electron gas, i.e. a white dwarf. See Kippenhahn & Weigert, Stellar Structure and Evolution, section 15.2
 ! Note that this is only the electron degeneracy pressure, so does not include the ion contribution to the pressure. This is a good approximation for white dwarfs where the electrons are highly degenerate but the ions are not.
 ! Note also that this assumes a fully ionised gas, so mu is the mean molecular weight per free electron (e.g. mu=2 for a pure helium gas)

 real :: x

    ! Convert rho to number density of electrons
    n_e = rhoi / (mu_e * atomic_mass_unit)

    x = cbrt(3*n_e * planckh**3 / (8 * pi * mass_electron_cgs**3 * c**3))
    fx= x*(2*x**2 - 3)*np.sqrt(x**2 + 1) + 3*np.log(x + np.sqrt(x**2 + 1))
    presi = (pi*(mass_electron_cgs)**4)*(c**5) * fx/ (3 * planckh**3) 

end subroutine get_zerotemp_Pressure


!----------------------------------------------------------------
!+
!  Calculates sound speed from density (do we need gamma?)
!+
!----------------------------------------------------------------


subroutine get_zerotemp_spsoundi(rhoi,presi,spsoundi,gammai)
 real, intent(in)  :: rhoi,presi
 real, intent(out) :: spsoundi,gammai

 gammai = 1. + presi/(eni*rhoi) ! what should I put this?
 spsoundi = sqrt(gammai*presi/rhoi)

end subroutine get_zerotemp_spsoundi

!----------------------------------------------------------------
!+
!  Calculates density from pressure (under construction)
!+
!----------------------------------------------------------------
subroutine get_zerotemp_rhofrompres(presi,tempi,mu,densi)
 real, intent(in)  :: presi,tempi,mu
 real, intent(out) :: densi

 densi = (presi - radconst*tempi**4 /3.) * mu / (Rg*tempi)

end subroutine get_zerotemp_rhofrompres

end module eos_zerotemp