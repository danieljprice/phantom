!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Give velocity perturbation to gas particles
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: part
!
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

