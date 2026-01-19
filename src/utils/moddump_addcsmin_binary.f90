!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters: None
!
! :Dependencies: eos, physcon, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use eos,     only:gmw,cs_min,ieos
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_velocity
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: T_floor

 T_floor = 10.
 if ( any( ieos==(/3,6,7,13,14/) ) ) then
    print "(/,a)",' Setting floor temperature to ', T_floor, ' K.'
    cs_min =  gmw*T_floor/(mass_proton_cgs/kboltz * unit_velocity**2)
    if (ieos == 3) then
       ieos = 14
    endif
 endif

end subroutine modify_dump

end module moddump
