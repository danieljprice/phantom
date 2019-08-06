!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: testcorotate
!
!  DESCRIPTION:
!   Unit tests of the coriolis & centrifugal force module
!
!  REFERENCES: Tejeda E., Rosswog S., 2013, MNRAS, 433, 1930
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: extern_corotate, io, testutils
!+
!--------------------------------------------------------------------------
module testcorotate
 implicit none
 public :: test_corotate

 private

contains

subroutine test_corotate(ntests,npass)
 use io,              only:id,master
 use testutils,       only:checkval,update_test_scores
 use extern_corotate, only:get_centrifugal_force,get_coriolis_force,omega_corotate
 integer, intent(inout) :: ntests,npass
 logical                :: test_centrifugal,test_coriolis
 real    :: r(3),v(3),f(3),fx,fy,fz,phi
 integer :: nfailed(4)

 if (id==master) write(*,"(/,a,/)") '--> TESTING COROTATION MODULE'

 test_centrifugal = .true.
 test_coriolis    = .true.
 !
 !--Test: check that centrifugal force for particle in z=0 plane is sensible
 !
 testcentrifugal: if (test_centrifugal) then
    if (id==master) write(*,"(/,a)") '--> testing centrifugal force with rotation about z axis'

    omega_corotate = sqrt(2.) ! some irrational number
    r = (/1.32,2.4313,0.548/) ! some arbitrary numbers
    fx = 0.
    fy = 0.
    fz = 0.
    phi = 0.
    call get_centrifugal_force(r,fx,fy,fz,phi)

    call checkval(fx,omega_corotate**2*r(1),epsilon(fx),nfailed(1),'fx=Omega^2 x')
    call checkval(fy,omega_corotate**2*r(2),epsilon(fy),nfailed(2),'fy=Omega^2 y')
    call checkval(fz,0.,epsilon(fz),nfailed(3),'fz=0')
    call checkval(phi,-0.5*omega_corotate**2*(r(1)**2 + r(2)**2),2.*epsilon(phi),nfailed(4),'phi=1/2 Omega^2 R^2')

    call update_test_scores(ntests,nfailed(1:4),npass)

 endif testcentrifugal

 !
 ! Test coriolis force compared to 2D cartesian expression
 !
 testcoriolis: if (test_coriolis) then
    if (id==master) write(*,"(/,a)") '--> testing coriolis force with rotation about z axis'

    omega_corotate = sqrt(2.) ! some irrational number
    r = (/1.32,2.4313,0.548/) ! some arbitrary numbers
    v = (/6.3,25.023,8.325/)  ! ditto
    phi = 0.
    call get_coriolis_force(r,v,f,phi)

    nfailed = 0
    call checkval(f(1),2.*omega_corotate*v(2),epsilon(fx),nfailed(1),'fx=2 Omega vy')
    call checkval(f(2),-2.*omega_corotate*v(1),epsilon(fy),nfailed(2),'fy=-2 Omega vx')
    call checkval(f(3),0.,epsilon(fz),nfailed(3),'fz=0')
    !call checkval(phi,-omega_corotate*(v(2)*r(1) - v(1)*r(2)),epsilon(phi),nfailed(4),'phi=Omega(vy*x-vx*y)')

    call update_test_scores(ntests,nfailed(1:4),npass)

 endif testcoriolis

 if (id==master) write(*,"(/,a)") '<-- COROTATE TEST COMPLETE'

end subroutine test_corotate

end module testcorotate
