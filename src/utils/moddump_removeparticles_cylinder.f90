!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Remove particles outside a cylinder
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: part
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part, only:delete_particles_outside_cylinder
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 !integer :: i
 real, dimension(3) :: center
 real :: radius,zmax

 print*,' Phantommoddump: Remove particles outside a cylinder'
 !
 !--set the center and the radius and the height of the cylinder
 !
 center(:)=0.
 radius=1500.
 zmax=1500.0
 !
 !--removing particles
 !
 print*,'Removing particles outside the cylinder centered in ( ', center,' ), with radius ',radius,' and zmax ',zmax,' : '
 call delete_particles_outside_cylinder(center, radius, zmax, npartoftype)

end subroutine modify_dump

end module moddump

