!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
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

 real, public,parameter :: Tion = 10000.
 real, public,parameter :: muion = 0.5
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

 polykion = (kb_on_mh*Tion/muion)/(unit_velocity**2)
 csion    = sqrt(polykion)
 Tcold    = polyk*gmw*temperature_coef


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
 !+
 !-----------------------------------------------------------------------
subroutine get_eos_HIIR_iso(polyk,temperature_coef,mui,tempi,ponrhoi,spsoundi)
 real, intent(in)    :: polyk,temperature_coef
 real, intent(out)   :: ponrhoi,spsoundi,mui
 real, intent(inout) :: tempi

 !
 !--dual medium isothermal eos
 !
 !  :math:`P = c_s^2 \rho`
 !
 !  where :math:`c_s^2 \equiv K` is a constant stored in the dump file header
 !
 if (tempi + epsilon(tempi) > Tion) then
    ponrhoi  = polykion
    spsoundi = sqrt(ponrhoi)
 else
    ponrhoi  = polyk
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi
 endif


end subroutine get_eos_HIIR_iso



end module eos_HIIR

