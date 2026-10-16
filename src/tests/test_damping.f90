!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testdamping
!
! Unit tests of damping module
!
! :References: exoALMA comparison project
!
! :Owner: Francesco Radica
!
! :Runtime parameters: None
!
! :Dependencies: damping, io, physcon, testutils, units
!
 implicit none
 public :: test_damping

 private

contains
!-----------------------------------------------------------------------
!+
!   Unit tests of the routines in the damping module
!+
!-----------------------------------------------------------------------
subroutine test_damping(ntests,npass)
 use io, only:id,master
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a,/)") '--> TESTING DAMPING TERMS'

 call test_damping_disc(ntests,npass)
 call test_damping_stellar(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- DAMPING TESTS COMPLETE'

end subroutine test_damping

!-----------------------------------------------------------------------
!+
!   Unit tests of the disc damping boundary conditions
!+
!-----------------------------------------------------------------------
subroutine test_damping_disc(ntests,npass)
 use damping,   only:idamp,damp,get_damp_fac_disc,r1in,r2in,r1out,r2out,calc_damp
 use physcon,   only:pi
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer :: nfail(1),ipos
 real :: xyz(3),v0(3),fac,r,damp_fac
 real :: omega,t_orb,t_damp,time
 real, parameter :: tol = 3.e-16
 real :: rpos(4)
 character(len=5) :: label(4)

 ! select disc damping
 idamp = 3
 time = 0.
 damp = 0.01

 ! set positions of damping zones
 r1in = 0.3; r2in = 2.
 r1out = 2.52; r2out = 3.0
 call calc_damp(time,damp_fac)

 ! set position at inner edge of first damping zone
 rpos  = [r1in,r2out,r1out,r2in]
 label = ['r1in ','r2out','r1out','r2in ']

 do ipos=1,size(rpos)
    r = rpos(ipos)
    xyz = [r,0.,0.]
    fac = damp_fac*get_damp_fac_disc(xyz,v0)

    ! check the damping time is the correct multiple
    ! of the orbital time at this radius
    omega = sqrt(1./r**3)
    t_orb = 2.*pi/omega

    if (ipos==3 .or. ipos==4) then
       call checkval(fac,0.,tol,nfail(1),'damping = 0 at r='//trim(label(ipos)))
    else
       t_damp = 1./fac
       call checkval(t_damp,t_orb/damp,tol,nfail(1),'t_damp = f*t_orb at r='//trim(label(ipos)))
    endif
    call update_test_scores(ntests,nfail,npass)

    call checkval(v0(2),r*omega,tol,nfail(1),'v = v_kep at r='//trim(label(ipos)))
    call update_test_scores(ntests,nfail,npass)
 enddo

end subroutine test_damping_disc

!-----------------------------------------------------------------------
!+
!   Isolated stellar damping: physical units, spin, bulk motion and release.
!+
!-----------------------------------------------------------------------
subroutine test_damping_stellar(ntests,npass)
 use damping,   only:idamp,damp_time_s,damp_stop_s,damp_spin_s,calc_damp,apply_damp,&
                    set_stellar_damping_frame,rigid_rotation_velocity
 use units,     only:utime
 use testutils, only:checkval,update_test_scores
 integer, intent(inout) :: ntests,npass
 integer :: nfail(1),old_idamp,i
 real :: old_time,old_stop,old_spin(3),xyz(3),vxyz(3),force(3),fac,omega(3)
 real :: centre(3),bulk(3),residual(3),masses(2),positions(3,2),velocities(3,2),totalforce(3)
 real(kind=8) :: old_utime
 real, parameter :: tol = max(1.e-12,64.*epsilon(1.))

 old_idamp = idamp
 nfail = 0
 old_utime = utime
 old_time = damp_time_s
 old_stop = damp_stop_s
 old_spin = damp_spin_s
 idamp = 4
 utime = 100.d0
 damp_time_s = 200.
 damp_stop_s = 500.
 damp_spin_s = [0.,0.,0.02]
 centre = [10.,-7.,5.]
 bulk = [3.,-2.,4.]
 xyz = centre + [2.,1.,-3.]
 omega = damp_spin_s*utime
 call set_stellar_damping_frame(centre,bulk)
 call calc_damp(0.,fac)
 call checkval(fac,0.5,tol,nfail(1),'stellar damping uses physical seconds')
 call update_test_scores(ntests,nfail,npass)

 vxyz = bulk + rigid_rotation_velocity(xyz-centre,omega)
 call checkval(vxyz(1),bulk(1)-omega(3),tol,nfail(1),'rigid spin sign matches omega cross r')
 call update_test_scores(ntests,nfail,npass)
 call checkval(vxyz(2),bulk(2)+2.*omega(3),tol,nfail(1),'rigid spin azimuthal velocity')
 call update_test_scores(ntests,nfail,npass)
 force = 0.
 call apply_damp(force(1),force(2),force(3),vxyz,xyz,fac)
 call checkval(sum(abs(force)),0.,tol,nfail(1),'rigidly rotating translated star is not damped')
 call update_test_scores(ntests,nfail,npass)

 residual = [1.,-2.,3.]
 force = [4.,5.,6.]
 call apply_damp(force(1),force(2),force(3),vxyz+residual,xyz,fac)
 call checkval(sum(abs(force-([4.,5.,6.]-fac*residual))),0.,tol,nfail(1),'drag removes only residual motion')
 call update_test_scores(ntests,nfail,npass)

 masses = [1.,2.]
 positions(:,1) = centre + [2.,0.,0.]
 positions(:,2) = centre - [1.,0.,0.]
 velocities(:,1) = bulk + [1.,2.,3.]
 velocities(:,2) = bulk - [0.5,1.,1.5]
 totalforce = 0.
 do i=1,2
    force = 0.
    call apply_damp(force(1),force(2),force(3),velocities(:,i),positions(:,i),fac)
    totalforce = totalforce + masses(i)*force
 enddo
 call checkval(sum(abs(totalforce)),0.,tol,nfail(1),'unequal-mass damping preserves total momentum')
 call update_test_scores(ntests,nfail,npass)

 call calc_damp(5.,fac)
 call checkval(fac,0.,tol,nfail(1),'stellar damping switches off at absolute physical time')
 call update_test_scores(ntests,nfail,npass)
 call calc_damp(6.,fac)
 call checkval(fac,0.,tol,nfail(1),'stellar damping remains off after restart')
 call update_test_scores(ntests,nfail,npass)
 force = 0.
 call apply_damp(force(1),force(2),force(3),vxyz+residual,xyz,fac)
 call checkval(sum(abs(force)),0.,tol,nfail(1),'released star has no drag')
 call update_test_scores(ntests,nfail,npass)

 damp_stop_s = 0.
 damp_spin_s = 0.
 call calc_damp(100.,fac)
 force = 0.
 call apply_damp(force(1),force(2),force(3),bulk,xyz,fac)
 call checkval(sum(abs(force)),0.,tol,nfail(1),'static relaxation preserves bulk velocity')
 call update_test_scores(ntests,nfail,npass)
 call checkval(fac,0.5,tol,nfail(1),'zero switch-off time means continuous damping')
 call update_test_scores(ntests,nfail,npass)

 omega = [1.,2.,3.]
 vxyz = rigid_rotation_velocity([4.,5.,6.],omega)
 call checkval(sum(abs(vxyz-[-3.,6.,-3.])),0.,tol,nfail(1),'rigid spin supports arbitrary axes')
 call update_test_scores(ntests,nfail,npass)

 idamp = old_idamp
 utime = old_utime
 damp_time_s = old_time
 damp_stop_s = old_stop
 damp_spin_s = old_spin
 call set_stellar_damping_frame([0.,0.,0.],[0.,0.,0.])

end subroutine test_damping_stellar

end module testdamping
