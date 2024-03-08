!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testdamping
!
! Unit tests of damping module
!
! :References: exoALMA comparison project
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: damping, io, physcon, testutils
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

end module testdamping
