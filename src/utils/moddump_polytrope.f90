!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! perturb a star with a radial velocity perturbation to excite
! the fundamental oscillation mode
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: systemutils
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = '--amp=1.e-4'

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use systemutils, only:get_command_option_real
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real    :: amp

 amp = get_command_option_real('amp',default=1.e-4)

 do i=1,npart
    vxyzu(1:3,i) = amp*xyzh(1:3,i)
 enddo

end subroutine modify_dump

end module moddump
