!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Convert boundary particles back into non-rigid, gas particles
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: part, boundarypart
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:iamtype,iamboundary,iphase,igas,iboundary,set_particle_type
 use boundarypart, only:set_boundary_particle_velocity
 integer, intent(inout) :: npart,npartoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),massoftype(:)
 integer                :: i,nconverted
 real                   :: rmin,rmin2,r2,boundary_com(3),xyz_CM(3),vxyz_CM(3)

 rmin = 5.
 rmin2 = rmin**2

 boundary_com = 0.
 call set_boundary_particle_velocity(npart,iphase,xyzh,xyz_CM,vxyz_CM)

 nconverted = 0
 do i=1,npart
    if (iamboundary(iamtype(iphase(i)))) then
       r2 = dot_product(xyzh(1:3,i)-xyz_CM,xyzh(1:3,i)-xyz_CM)
       if (r2 > rmin2) then
          call set_particle_type(i,igas)
          nconverted = nconverted + 1
       endif
    endif
 enddo
 npartoftype(iboundary) = npartoftype(iboundary) - nconverted
 npartoftype(igas) = npartoftype(igas) + nconverted

end subroutine modify_dump

end module moddump
