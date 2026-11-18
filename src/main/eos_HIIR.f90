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

 public :: get_eos_HIIR,init_eos_HIIR

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
 use io,      only:fatal

 integer, intent(out) :: ierr
 real,    intent(in)  :: gamma,polyk,gmw,temperature_coef

 polykion  = (kb_on_mh*Tion/muion)/(unit_velocity**2)
 csion     = sqrt(polykion)
 Tcold     = polyk*gmw*temperature_coef

 if (gamma <= 1.) call fatal('eos_HIIR','ieos=21 needs gamma > 1 for the gas of the front',var='gamma',val=gamma)
 if (gmw <= muion) call fatal('eos_HIIR','gmw lower than muion',var='gmw',val=gmw)
 uIon = (kb_on_mh*Tion/(muion*(gamma-1.)))/(unit_velocity**2)

 ierr = 0

end subroutine init_eos_HIIR

 !-----------------------------------------------------------------------
 !+
 !  Main eos routine
 !
 !  two isothermal phases, P = c_s^2 rho, with c_s^2 = polyk in the cold gas
 !  (T <= Tcold) and polykion in the ionised gas (T >= Tion). The gas in
 !  between, at the ionisation front, is adiabatic, P = (gamma-1) rho u, its
 !  thermal energy relaxing to its phase (see HIIRegion)
 !+
 !-----------------------------------------------------------------------
subroutine get_eos_HIIR(polyk,temperature_coef,mui,gammai,eni,tempi,ponrhoi,spsoundi)
 real, intent(in)  :: polyk,temperature_coef,mui,gammai,eni
 real, intent(out) :: tempi,ponrhoi,spsoundi

 if (mui > muion) then !- neutral or partially ionised
    tempi = (gammai - 1.)*eni*mui*temperature_coef
    if (tempi > Tcold) then !- if hot then in a shell
       ponrhoi  = (gammai - 1.)*eni
       spsoundi = sqrt(gammai*ponrhoi)
    else !- else cold
       ponrhoi  = polyk
       spsoundi = sqrt(polyk)
       tempi    = Tcold
    endif
 else
    ponrhoi  = polykion
    spsoundi = csion
    tempi    = Tion
 endif

end subroutine get_eos_HIIR

end module eos_HIIR

