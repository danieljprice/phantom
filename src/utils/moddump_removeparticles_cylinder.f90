!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part
!+
!--------------------------------------------------------------------------
module moddump

 use part, only:delete_particles_outside_cylinder

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 implicit none
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
 radius=60.
 zmax=5.0
 !
 !--removing particles
 !
 print*,'Removing particles outside the cylinder centered in ( ', center,' ), with radius ',radius,' and zmax ',zmax,' : '
 call delete_particles_outside_cylinder(center, radius, zmax)

 return
end subroutine modify_dump

end module moddump

