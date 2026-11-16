!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_HIIR
!
! eos_HIIR
!
! :References: None
!
! :Owner: Yann Bernard
!
! :Runtime parameters: None
!
! :Dependencies: io, physcon, units
!
 implicit none

 public :: get_eos_HIIR_iso,init_eos_HIIR

 real, public, parameter :: Tion = 10000.
 real, public, parameter :: muion = 0.5
 real, public           :: polykion
 real, public           :: csion
 real, public           :: uIon
 real, public           :: Tcold

 private

contains

 !-----------------------------------------------------------------------
 !+
 !  Init eos routine
 !+
 !-----------------------------------------------------------------------

subroutine init_eos_HIIR(gamma,polyk,gmw,temperature_coef,ierr)
 use physcon, only:kb_on_mh
 use units,   only:unit_velocity

 integer, intent(out) :: ierr
 real,    intent(in)  :: gamma,polyk,gmw,temperature_coef

 polykion  = (kb_on_mh*Tion/muion)/(unit_velocity**2)
 csion     = sqrt(polykion)
 Tcold     = polyk*gmw*temperature_coef

 if (gamma>1.) then
    uIon = (kb_on_mh*Tion/(muion*(gamma-1.)))/(unit_velocity**2)
 else
    uIon = 1.5*polykion
 endif

 ierr = 0

end subroutine init_eos_HIIR

 !-----------------------------------------------------------------------
 !+
 !  Main eos routine (isothermal)
 !
 !  two isothermal phases, P = c_s^2 rho, with c_s^2 = polyk in the cold gas
 !  and polykion in the ionised gas. The ionised fraction phi of the front
 !  particle is read from its temperature T = phi*Tion + (1-phi)*Tcold, set by
 !  the HII routines, and P/rho is mixed in the same way. Anything outside
 !  [Tcold,Tion] (e.g. T not set yet) is clipped to the nearest phase
 !+
 !-----------------------------------------------------------------------
subroutine get_eos_HIIR_iso(polyk,temperature_coef,mui,tempi,ponrhoi,spsoundi)
 real, intent(in)    :: polyk,temperature_coef,mui
 real, intent(out)   :: ponrhoi,spsoundi
 real, intent(inout) :: tempi
 real :: phi

 phi      = min(max((tempi - Tcold)/(Tion - Tcold),0.),1.)
 ponrhoi  = phi*polykion + (1. - phi)*polyk
 spsoundi = sqrt(ponrhoi)
 tempi    = phi*Tion + (1. - phi)*Tcold

end subroutine get_eos_HIIR_iso

end module eos_HIIR

