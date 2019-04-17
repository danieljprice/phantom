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
!  DEPENDENCIES: part, prompting
!+
!--------------------------------------------------------------------------
module moddump

 use part,         only:delete_particles_inside_radius,delete_particles_outside_sphere
 use prompting,    only:prompt

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 implicit none
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 real, dimension(3) :: incenter,outcenter
 real :: inradius,outradius
 logical :: icutinside,icutoutside

 icutinside   = .false.
 icutoutside  = .false.
 incenter(:)  = 0.
 outcenter(:) = 0.
 inradius     = 10.
 outradius    = 200.

 !
 !--set the centers and the radius
 !
 call prompt('Deleting particles inside a given radius ?',icutinside)
 call prompt('Deleting particles outside a given radius ?',icutoutside)
 if (icutinside) then
    call prompt('Enter the inward radius in au',inradius,0.)
    call prompt('Enter the x coordinates of the center of that sphere',incenter(1),0.)
    call prompt('Enter the y coordinates of the center of that sphere',incenter(2),0.)
    call prompt('Enter the z coordinates of the center of that sphere',incenter(3),0.)
 endif
 if (icutoutside) then
    call prompt('Enter the outward radius in au',outradius,0.)
    call prompt('Enter the x coordinates of the center of that sphere',outcenter(1),0.)
    call prompt('Enter the y coordinates of the center of that sphere',outcenter(2),0.)
    call prompt('Enter the z coordinates of the center of that sphere',outcenter(3),0.)
 endif

 if (icutinside) then
    print*,'Phantommoddump: Remove particles inside a particular radius'
    print*,'Removing particles inside radius ',inradius
    call delete_particles_inside_radius(incenter,inradius,npart,npartoftype)
 endif

 if (icutoutside) then
    print*,'Phantommoddump: Remove particles outside a particular radius'
    print*,'Removing particles outside radius ',outradius
    call delete_particles_outside_sphere(outcenter,outradius)
 endif

end subroutine modify_dump

end module moddump

