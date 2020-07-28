!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2020 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! default moddump routine: does not make any modifications
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 print*,' no modifications performed '

 return
end subroutine modify_dump

end module moddump

