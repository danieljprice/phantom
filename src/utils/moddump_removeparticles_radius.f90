!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Remove particles inside or outside a sphere
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: part, prompting, units
!

 use part,         only:delete_particles_outside_sphere,delete_particles_with_large_h,igas,idust
 use prompting,    only:prompt
 use units,        only:unit_density

 implicit none

 character(len=*), parameter, public :: moddump_flags = ''

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 integer, intent(inout) :: npart
 integer, dimension(:), intent(inout) :: npartoftype
 real, dimension(:), intent(inout) :: massoftype
 real, dimension(:,:), intent(inout) :: xyzh,vxyzu
 real, dimension(3) :: incenter,outcenter
 integer :: iremoveparttype,npart_old
 real :: inradius,outradius
 real :: h_on_r_min,rmax,rho_max
 logical :: icutinside,icutoutside,irmlargeh

 npart_old = npart
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
 irmlargeh = .false.
 if (.not. icutinside .and. .not. icutoutside) then
    call prompt('Deleting particles with large h?',irmlargeh)
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

 if (irmlargeh) then
    h_on_r_min = 0.3
    rho_max = 1.e-14/unit_density
    rmax = 1000.
    outcenter = 0.
    print*,'Phantommoddump: Removing particles with large h, satisfying all of'
    print*,'h/r >',h_on_r_min,', rho [code] <',rho_max,', r [code] >',rmax
    call delete_particles_with_large_h(outcenter,npart,h_on_r_min,rho_max,rmax)
 endif

 print*,'Removed a total of ',npart_old-npart,' particles'

end subroutine modify_dump

end module moddump

