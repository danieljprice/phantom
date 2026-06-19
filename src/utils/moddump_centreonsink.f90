!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Center positions and velocities of particles around one given sink
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters: None
!
! :Dependencies: part
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,                 only:xyzmh_ptmass,vxyz_ptmass,nptmass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i,sink_ind

 sink_ind=1

 if (nptmass < sink_ind) then
    print*,'Selected sink index larger than number of sinks'
    return
 endif

 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) - xyzmh_ptmass(1:3,sink_ind)
    vxyzu(1:3,i) = vxyzu(1:3,i) - vxyz_ptmass(1:3,sink_ind)
 enddo

 do i=1,nptmass
    !skip sink_ind, because useless if put to 0 first
    if (i==sink_ind) then
       cycle
    else
       xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,sink_ind)
       vxyz_ptmass(1:3,i) = vxyz_ptmass(1:3,i) - vxyz_ptmass(1:3,sink_ind)
    endif
 enddo

 xyzmh_ptmass(1:3,sink_ind) = xyzmh_ptmass(1:3,sink_ind) - xyzmh_ptmass(1:3,sink_ind)
 vxyz_ptmass(1:3,sink_ind) = vxyz_ptmass(1:3,sink_ind) - vxyz_ptmass(1:3,sink_ind)

end subroutine modify_dump

end module moddump

