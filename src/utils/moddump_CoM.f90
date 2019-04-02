!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: moddump
!
!  DESCRIPTION:
!  This will reset the particles to the centre of mass; required for
!  some plotting/analysis
!
!  REFERENCES: None
!
!  OWNER: James Wurster
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, centreofmass, dim, part
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only: nptmass,xyzmh_ptmass,vxyz_ptmass
 use centreofmass, only: reset_centreofmass
 use dim,          only: periodic
 use boundary,     only: xmin,ymin,zmin,dxbound,dybound,dzbound
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i
 !
 !--Reset centre of mass
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)
 !
 !--Ensure we're not going outside the box if periodic
 if (periodic) then
    do i = 1,npart
       if (xyzh(1,i) < xmin)         xyzh(1,i) = xyzh(1,i) + dxbound
       if (xyzh(1,i) > xmin+dxbound) xyzh(1,i) = xyzh(1,i) - dxbound
       if (xyzh(2,i) < ymin)         xyzh(2,i) = xyzh(2,i) + dybound
       if (xyzh(2,i) > ymin+dybound) xyzh(2,i) = xyzh(2,i) - dybound
       if (xyzh(3,i) < zmin)         xyzh(3,i) = xyzh(3,i) + dzbound
       if (xyzh(3,i) > zmin+dzbound) xyzh(3,i) = xyzh(3,i) - dzbound
    enddo
 endif
 !
 return
end subroutine modify_dump

end module moddump

