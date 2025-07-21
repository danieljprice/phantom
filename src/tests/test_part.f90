!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testpart
!
! Unit tests of the part module, checking basic functionality
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, io, kernel, mpidomain, part, physcon, testutils,
!   unifdis
!
 use testutils, only:update_test_scores,checkval
 implicit none

 public :: test_part

 private

contains
!----------------------------------------------------------
!+
!  unit tests of equation of state module
!+
!----------------------------------------------------------
subroutine test_part(ntests,npass)
 use io,        only:id,master,stdout
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING PART MODULE'

 call test_accrete_and_kill_routines(ntests, npass)

 if (id==master) write(*,"(/,a)") '<-- PART MODULE TEST COMPLETE'

end subroutine test_part


!----------------------------------------------------------
!+
!  test that the initialisation of all eos works correctly
!+
!----------------------------------------------------------
subroutine test_accrete_and_kill_routines(ntests, npass)
 use io,        only:id,master
 use part,      only:xyzh,npartoftype,npart,periodic,init_part,shuffle_part,&
                     igas,isetphase,iphase,maxphase,maxp
 use part,      only:kill_particle,count_dead_particles,isdead_or_accreted, &
                     accrete_particles_outside_sphere,delete_dead_or_accreted_particles
 use boundary,  only:dxbound,dybound,dzbound,xmin,xmax,ymin,ymax,zmin,zmax
 use mpidomain, only:i_belong
 use kernel,    only:hfact_default
 use unifdis,   only:set_unifdis
 use physcon,   only:pi
 integer, intent(inout) :: ntests,npass
 integer :: nfailed(1)
 integer :: i,nkilled,ncount,npartold,nacc,nacc_check
 real :: rhozero,psep,totmass

 if (id==master) write(*,"(/,a)") '--> testing kill_particle and associated routines'

 !
 !--set up particles in a uniform distribution
 !
 call init_part()
 npartoftype(:) = 0
 npart = 0
 rhozero = 1.
 totmass = rhozero*dxbound*dybound*dzbound
 psep = dxbound/100.
 call set_unifdis('cubic',id,master,xmin,xmax,ymin,ymax,zmin,zmax,psep,&
                  hfact_default,npart,xyzh,periodic,mask=i_belong)
 nfailed = 0
 npartoftype(1) = npart
 if (maxphase==maxp) iphase(1:npart) = isetphase(igas,.true.)
 !
 ! accrete particles outside r=1
 !
 nacc = 0
 call accrete_particles_outside_sphere(xmax)
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) nacc = nacc + 1
 enddo
 !
 ! check number of accreted particles is approximately correct
 !
 nacc_check = int((1. - pi/6.)*npart)  ! ratio of volumes: [ (2 r)^3 - 4/3*pi*r^3 ] / (2r)^3
 call checkval(nacc,nacc_check,400,nfailed(1),'accrete_particles_outside_sphere')
 call update_test_scores(ntests,nfailed,npass)
 !
 ! kill every 10th particle, except those already accreted
 !
 nkilled = 0
 do i=1,npart
    if (mod(i,10)==0 .and. .not.isdead_or_accreted(xyzh(4,i))) then
       call kill_particle(i,npartoftype)
       nkilled = nkilled + 1
    endif
 enddo
 !
 ! check that the number killed matches the actual
 !
 ncount  = count_dead_particles()
 call checkval(ncount,nkilled,0,nfailed(1),'count_dead_particles matches n killed')
 call update_test_scores(ntests,nfailed,npass)
 !
 ! call shuffle_part and check revised particle number is correct
 !
 npartold = npart
 call shuffle_part(npart)
 call checkval(npart,npartold-nkilled,0,nfailed(1),'npart after call to shuffle_part')
 call update_test_scores(ntests,nfailed,npass)
 !
 ! check delete_dead_or_accreted_particles works
 !
 npartold = npart
 call delete_dead_or_accreted_particles(npart,npartoftype)
 call checkval(npart,npartold-nacc,0,nfailed(1),'delete_dead_or_accreted has right n')
 call update_test_scores(ntests,nfailed,npass)
 !
 ! check sum(npartoftype) = npart
 !
 call checkval(sum(npartoftype),npart,0,nfailed(1),'npart = sum(npartoftype)')
 call update_test_scores(ntests,nfailed,npass)

end subroutine test_accrete_and_kill_routines

end module testpart
