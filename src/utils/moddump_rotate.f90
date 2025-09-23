!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! set solid body rotation for all gas particles about the z-axis given angular frequency
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: systemutils
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = '--omega=7.92e-3'

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use systemutils, only:get_command_option_real
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: omega,vphi,R
 integer                :: i

 ! Assume rotation axis is the z-axis
 omega = get_command_option_real('omega',default=7.92e-3)

 print*,"Adding solid body rotation with omega = ",omega

 do i = 1,npart
    R = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
    vphi = omega*R
    vxyzu(1,i) = -omega*xyzh(2,i)
    vxyzu(2,i) = omega*xyzh(1,i)
    vxyzu(3,i) = 0.
 enddo

end subroutine modify_dump

end module moddump

