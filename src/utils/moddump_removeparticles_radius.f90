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
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: part
!+
!--------------------------------------------------------------------------
module moddump

 use part, only:delete_particles_inside_radius

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 implicit none
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 real, dimension(3) :: center
 real :: radius

 print*,' Phantommoddump: Remove particles inside a particular radius'
 !
 !--set the center and the radius
 !
 center(:)=0.
 radius=10.
 !
 !--removing particles
 !
 print*,'Removing particles inside radius ',radius
 call delete_particles_inside_radius(center,radius,npart,npartoftype)

end subroutine modify_dump

end module moddump

