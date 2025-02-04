!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! None
!
! :References: None
!
! :Owner: Daniel Mentiplay
!
! :Runtime parameters: None
!
! :Dependencies: part, prompting
!

 use part,         only:delete_particles_outside_sphere,igas,idust
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
 integer :: iremoveparttype
 real :: inradius,outradius
 logical :: icutinside,icutoutside


 !
 !--set the centers and the radius
 !
 call prompt('Deleting particles inside a given radius ?',icutinside)
 call prompt('Deleting particles outside a given radius ?',icutoutside)
 if (icutinside) then
    call prompt('Enter inward radius in code units',inradius,0.)
    call prompt('Enter x coordinate of the center of that sphere',incenter(1))
    call prompt('Enter y coordinate of the center of that sphere',incenter(2))
    call prompt('Enter z coordinate of the center of that sphere',incenter(3))
 endif
 if (icutoutside) then
    call prompt('Enter outward radius in code units',outradius,0.)
    call prompt('Enter x coordinate of the center of that sphere',outcenter(1))
    call prompt('Enter y coordinate of the center of that sphere',outcenter(2))
    call prompt('Enter z coordinate of the center of that sphere',outcenter(3))
 endif

 call prompt('Deleting which particles (0=all, 1=gas only, 2=dust only)?', iremoveparttype)
 ! add other types of particles here if needed
 select case (iremoveparttype)
 case (1)
    iremoveparttype = igas
 case (2)
    iremoveparttype = idust
 case default
    iremoveparttype = 0
 end select

 if (icutinside) then
    print*,'Phantommoddump: Remove particles inside a particular radius'
    print*,'Removing particles inside radius ',inradius
    if (iremoveparttype > 0) then
       print*,'Removing particles type ',iremoveparttype
       call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.,mytype=iremoveparttype)
    else
       call delete_particles_outside_sphere(incenter,inradius,npart,revert=.true.)
    endif
 endif

 if (icutoutside) then
    print*,'Phantommoddump: Remove particles outside a particular radius'
    print*,'Removing particles outside radius ',outradius
    if (iremoveparttype > 0) then
       print*,'Removing particles type ',iremoveparttype
       call delete_particles_outside_sphere(outcenter,outradius,npart,mytype=iremoveparttype)
    else
       call delete_particles_outside_sphere(outcenter,outradius,npart)
    endif
 endif

end subroutine modify_dump

end module moddump

