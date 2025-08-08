!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! moddump routine: store temperature from cluster (HIIfeedback) setup
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use HIIRegion, only:HII_feedback,initialize_H2R,update_ionrates,HII_feedback_ray,iH2R
 use part,      only:xyzmh_ptmass,vxyz_ptmass,nptmass,eos_vars,itemp,&
                     delete_dead_or_accreted_particles,accrete_particles_outside_sphere
 use ptmass,    only:h_acc
 use linklist,  only:set_linklist
 use io,        only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,isinkdeadhead,n,nsinkdead
 integer :: ll(nptmass)

 ll(:) = 0
 call accrete_particles_outside_sphere(10.)
 isinkdeadhead = -1
 nsinkdead = 0
 iH2R=1

 ! list dead sinks
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) < .001 .or. xyzmh_ptmass(5,i)==h_acc) then
       xyzmh_ptmass(:,i) = 0.
       vxyz_ptmass(:,i) = 0.
       ll(i) = isinkdeadhead
       isinkdeadhead = i
       nsinkdead = nsinkdead + 1
    endif
 enddo

 ! remove dead sinks and shuffle ptmass arrays
 n = nptmass
 do while(isinkdeadhead>0)
    if (isinkdeadhead <= n) then
       if (xyzmh_ptmass(4,n) > 0.)then
          xyzmh_ptmass(:,isinkdeadhead) = xyzmh_ptmass(:,n)
          vxyz_ptmass(:,isinkdeadhead) = vxyz_ptmass(:,n)
          isinkdeadhead = ll(isinkdeadhead)
       endif
       n = n - 1
    else
       isinkdeadhead = ll(isinkdeadhead)
    endif
    if (n < 0) call fatal('shuffle','npart < 0')
 enddo

 nptmass = nptmass-nsinkdead

 print*, "number of dead sink particles :",nsinkdead

 call set_linklist(npart,npart,xyzh,vxyzu)
 call initialize_H2R()
 call update_ionrates(nptmass,xyzmh_ptmass,h_acc)
 call HII_feedback(nptmass,npart,xyzh,xyzmh_ptmass,vxyzu,eos_vars,0.)

 !call delete_dead_or_accreted_particles(npart,npartoftype)

end subroutine modify_dump

end module moddump

