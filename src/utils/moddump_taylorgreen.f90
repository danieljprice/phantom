!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Add Taylor-Green velocity field
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use physcon, only:pi
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i
 real :: vzero

 print*,' Adding velocity field for Taylor-Green problem '
 vzero = 0.1
 print*,' vzero = ',vzero
 do i=1,npart
    vxyzu(1:3,i) = 0.
    !--velocity field for the Taylor-Green problem
    vxyzu(1,i) =  vzero*sin(2.*pi*xyzh(1,i))*cos(2.*pi*xyzh(2,i))
    vxyzu(2,i) = -vzero*cos(2.*pi*xyzh(1,i))*sin(2.*pi*xyzh(2,i))
 enddo

 return
end subroutine modify_dump

end module moddump

