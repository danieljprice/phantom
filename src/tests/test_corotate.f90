!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testcorotate
!
! Unit tests of the coriolis & centrifugal force module
!
! :References: Tejeda E., Rosswog S., 2013, MNRAS, 433, 1930
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: extern_corotate, externalforces, io, part, ptmass,
!   setbinary, testutils
!
 implicit none
 public :: test_corotate

 private

contains

subroutine test_corotate(ntests,npass)
 use io,              only:id,master
 use testutils,       only:checkval,update_test_scores
 use extern_corotate, only:get_centrifugal_force,get_coriolis_force,omega_corotate
 integer, intent(inout) :: ntests,npass
 logical                :: test_centrifugal,test_coriolis,do_test_sinkbinary
 real    :: r(3),v(3),f(3),fx,fy,fz,phi
 integer :: nfailed(4)

 if (id==master) write(*,"(/,a,/)") '--> TESTING COROTATION MODULE'

 test_centrifugal = .true.
 test_coriolis    = .true.
 do_test_sinkbinary  = .true.
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

 if (do_test_sinkbinary) call test_sinkbinary(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- COROTATE TEST COMPLETE'

end subroutine test_corotate

!-------------------------------------------------------------------------
!+
!  Test that a binary orbit for sink particles works in a corotating frame
!+
!--------------------------------------------------------------------------
subroutine test_sinkbinary(ntests,npass)
 use io,              only:id,master
 use testutils,       only:checkval,update_test_scores
 use extern_corotate, only:get_centrifugal_force,get_coriolis_force,omega_corotate
 use part,            only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass,dsdt_ptmass
 use setbinary,       only:set_binary
 use externalforces,  only:iext_corotate
 use ptmass,          only:get_accel_sink_sink
 integer, intent(inout) :: ntests,npass
 real :: m1,m2,a,e,hacc,phitot,dtsinksink,ti
 integer :: ierr,nfailed(5),merge_ij(2),merge_n

 if (id==master) write(*,"(/,a)") '--> testing corotating sink particle binary'

 nptmass = 0
 m1 = 14.272  ! arbitrary masses
 m2 = 11.746
 a = 10.61
 e = 0.
 hacc = 1.
 call set_binary(m1,m2,a,e,hacc,hacc,&
                 xyzmh_ptmass,vxyz_ptmass,nptmass,ierr,omega_corotate,verbose=.false.)
 !
 !--check that initial velocities are zero if a corotating frame is used
 !
 call checkval(ierr,0,0,nfailed(1),'ierr=0 from set_binary')
 call checkval(3,vxyz_ptmass(1:3,1),0.,epsilon(0.),nfailed(2),'v1=0 from set_binary')
 call checkval(3,vxyz_ptmass(1:3,2),0.,epsilon(0.),nfailed(3),'v2=0 from set_binary')
 !
 !--check that sink-sink forces are zero
 !
 ti = 0.
 call get_accel_sink_sink(nptmass,xyzmh_ptmass,fxyz_ptmass,phitot,dtsinksink,&
                          iext_corotate,ti,merge_ij,merge_n,dsdt_ptmass)
 call checkval(3,fxyz_ptmass(1:3,1),0.,epsilon(0.),nfailed(4),'sink-sink force1 = 0')
 call checkval(3,fxyz_ptmass(1:3,2),0.,epsilon(0.),nfailed(5),'sink-sink force2 = 0')

 call update_test_scores(ntests,nfailed,npass)

end subroutine test_sinkbinary

end module testcorotate
