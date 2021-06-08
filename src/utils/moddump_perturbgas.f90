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
!  Give velocity perturbation to gas particles
!
!  REFERENCES: None
!
!  OWNER: Mike Lau
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module moddump
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,              only:xyzmh_ptmass,vxyz_ptmass,nptmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: perturb_factor,sink_perturb_factor
 integer                :: i
 logical                :: perturb_sink

 perturb_sink = .false.
 perturb_factor = 1.5
 sink_perturb_factor = 1.5

 do i=1,npart
    vxyzu(1:3,i) = perturb_factor * vxyzu(1:3,i)
 enddo

 if (perturb_sink) then
    do i=1,nptmass
       vxyz_ptmass(1:3,i) = sink_perturb_factor * vxyz_ptmass(1:3,i)
    enddo
 endif

 return
end subroutine modify_dump

end module moddump

