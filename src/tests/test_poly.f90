!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testpoly
!
! Unit tests of the polynomial solver modules
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, quartic, testutils
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_poly

 private

contains
!--------------------------------------------
!+
!  Unit tests of the polynomial solvers
!+
!--------------------------------------------
subroutine test_poly(ntests,npass)
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a)") '--> TESTING POLYNOMIAL SOLVERS'

 call test_quartic(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- POLYNOMIAL TESTS COMPLETE'

end subroutine test_poly

!--------------------------------------------
!+
!  Unit tests of the quartic solver
!+
!--------------------------------------------
subroutine test_quartic(ntests,npass)
 use quartic, only:quarticsolve,quarticf
 integer, intent(inout) :: ntests,npass
 real :: a(0:3),xold,x
 integer :: ierr,nfail(2)
 logical :: moresweep
 real, parameter :: tol = 1.e-6

 if (id==master) write(*,"(/,a)") '--> checking x^4 = 16 gives x=2'

 a = [-16.,0.,0.,0.] ! a0 + a1*x + a2*x^2 + a3*x^3 + x^4 = 0
 xold = -2e6
 !print*,'a = ',a
 call quarticsolve(a,xold,x,moresweep,ierr)
 call checkval(x,2.,tol,nfail(1),'x=2')
 call checkval(ierr,0,0,nfail(2),'ierr=0')
 call checkval(quarticf(a,x),0.,tol,nfail(1),'f(x)=0 for solution found')
 call update_test_scores(ntests,nfail,npass)

 if (id==master) write(*,"(/,a)") '--> checking arbitrary quartic'
 !a = [-6.,-2.,3.,0.]  ! this one fails: solution is real and negative
 !a = [-6.,2.,3.,0.]
 !a = [-6.,2.,3.,2.]
 a = [-6.,2.,3.,1.]
 !a = [0.,1.,-2.,0.]   ! x = 0 is a solution: this also fails
 xold = 1.
 call quarticsolve(a,xold,x,moresweep,ierr)
 call checkval(quarticf(a,x),0.,tol,nfail(1),'f(x)=0 for solution found')
 call checkval(ierr,0,0,nfail(2),'ierr=0')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_quartic

end module testpoly
