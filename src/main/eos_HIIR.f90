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

 public :: get_eos_HIIR_iso,get_eos_HIIR_adiab,init_eos_HIIR

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
 !  Main eos routine (isothermal)
 !+
 !-----------------------------------------------------------------------
subroutine get_eos_HIIR_iso(polyk,temperature_coef,mui,tempi,ponrhoi,spsoundi,isionisedi)
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
 if (isionisedi) then
    ponrhoi  = polykion
    spsoundi = sqrt(ponrhoi)
    tempi    = Tion
 else
    ponrhoi  = polyk
    spsoundi = sqrt(ponrhoi)
    tempi    = temperature_coef*mui*ponrhoi
 endif


end subroutine get_eos_HIIR_iso


 !-----------------------------------------------------------------------
 !+
 !  Main eos routine (adiabatic)
 !+
 !-----------------------------------------------------------------------
subroutine get_eos_HIIR_adiab(polyk,temperature_coef,mui,tempi,ponrhoi,rhoi,eni,gammai,spsoundi,isionisedi)
 use io, only:fatal
 real,    intent(in)              :: polyk,temperature_coef,rhoi,gammai
 real,    intent(out)             :: ponrhoi,spsoundi,mui,tempi
 logical, intent(in)              :: isionisedi
 real,    intent(in),    optional :: eni


 if (gammai < tiny(gammai)) call fatal('eos','gamma not set for adiabatic eos',var='gamma',val=gammai)


 if (isionisedi) then
    ponrhoi  = polykion
    spsoundi = sqrt(ponrhoi)
    tempi    = Tion
 else
    if (present(eni)) then
       if (eni < 0.) then
          !write(iprint,'(a,Es18.4,a,4Es18.4)')'Warning: eos: u = ',eni,' < 0 at {x,y,z,rho} = ',xi,yi,zi,rhoi
          call fatal('eos','utherm < 0',var='u',val=eni)
       endif
       if (gammai > 1.0001) then
          ponrhoi = (gammai-1.)*eni   ! use this if en is thermal energy
       else
          ponrhoi = 2./3.*eni ! en is thermal energy and gamma = 1
       endif
    else
       ponrhoi = polyk*rhoi**(gammai-1.)
    endif
    spsoundi = sqrt(gammai*ponrhoi)

    tempi = temperature_coef*mui*ponrhoi
 endif


end subroutine get_eos_HIIR_adiab



end module eos_HIIR

