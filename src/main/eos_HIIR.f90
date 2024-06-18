!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_HIIR
 !
 ! Implements Two temperature eos for HII region expansion
 !
 ! :References: None
 !
 ! :Owner: Yann Bernard
 !
 ! :Runtime parameters: None
 !
 ! :Dependencies: None
 !
 implicit none

 public :: get_eos_HIIR,init_eos_HIIR

 real, parameter :: Tion = 10000.
 real, parameter :: muioninv = 2.
 real, parameter :: muion = 0.5

 real, public    :: polykion

 private

contains

 !-----------------------------------------------------------------------
 !+
 !  Init eos routine
 !+
 !-----------------------------------------------------------------------

subroutine init_eos_HIIR
 use physcon, only:kb_on_mh
 use units,   only:unit_velocity

 polykion = (muioninv*kb_on_mh*Tion)/(unit_velocity**2)


end subroutine init_eos_HIIR


 !-----------------------------------------------------------------------
 !+
 !  Main eos routine
 !+
 !-----------------------------------------------------------------------
subroutine get_eos_HIIR(polyk,temperature_coef,mui,tempi,ponrhoi,spsoundi,isionisedi)
 real, intent(in)    :: polyk,temperature_coef
 real, intent(out)   :: ponrhoi,spsoundi,mui,tempi
 logical, intent(in) :: isionisedi

 !
 !--dual medium isothermal eos
 !
 !  :math:`P = c_s^2 \rho`
 !
 !  where :math:`c_s^2 \equiv K` is a constant stored in the dump file header
 !
 if(isionisedi) then
    ponrhoi  = polykion
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*muion*ponrhoi
 else
    ponrhoi  = polyk
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi
 endif


end subroutine get_eos_HIIR

end module eos_HIIR

