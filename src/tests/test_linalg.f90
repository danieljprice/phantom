!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testlinalg
!
! Unit tests of linear algebra
!
! :References:
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: io, linalg, testutils, vectorutils
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_linalg

 private
 real :: A(3,3)

contains
!--------------------------------------------
!+
!  Unit tests of the polynomial solvers
!+
!--------------------------------------------
subroutine test_linalg(ntests,npass)
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a)") '--> TESTING LINEAR ALGEBRA'

 call test_matrix_inversion(ntests,npass)
 call test_bicgstab_quartic(ntests,npass)

 if (id==master) write(*,"(/,a)") '<-- LINEAR ALGEBRA TESTS COMPLETE'

end subroutine test_linalg

!--------------------------------------------
!+
!  Unit tests of the quartic solver
!+
!--------------------------------------------
subroutine test_matrix_inversion(ntests,npass)
 use linalg, only:solve_bicgstab,get_Ax_interface
 use vectorutils, only:matrixinvert3D
 integer, intent(inout) :: ntests,npass
 integer :: ierr,nfail(2),n
 real, parameter :: tol = 1.e-6
 real :: Ainv(3,3),b(3),x(3),xexact(3)
 procedure(get_Ax_interface), pointer :: get_Ax_routine

 if (id==master) write(*,"(/,a)") '--> checking BiCGSTAB solver'

 A = reshape((/10.,2.,3.,&
               4.,-10.,2.,&
               2.,1.,-10./), (/3,3/))
 b = (/5.,0.,1./)
 x = (/0.,0.,0./) ! initial guess
 n = 3

 ! exact solution
 call matrixinvert3D(A,Ainv,ierr)
 xexact = matmul(Ainv,b)

 get_Ax_routine=>get_Ax_local
 call solve_bicgstab(n,b,get_Ax_routine,x,ierr)
 call checkval(ierr,0,0,nfail(1),'ierr=0')
 call checkval(3,x,xexact,tol,nfail(2),'xexact')

 call update_test_scores(ntests,nfail,npass)

contains

function get_Ax_local(n,x) result(Ax)
 integer, intent(in) :: n
 real, intent(in) :: x(n)
 real :: Ax(n)

 Ax = matmul(A,x)

end function get_Ax_local

end subroutine test_matrix_inversion

!--------------------------------------------
!+
!  BiCGSTAB quartic diagonal tests
!+
!--------------------------------------------
subroutine test_bicgstab_quartic(ntests,npass)
 use linalg,  only:solve_bicgstab,get_Ax_interface
 use quartic, only:quarticsolve
 integer, intent(inout) :: ntests,npass
 integer, parameter :: ndim = 3, max_outer = 10
 real,    parameter :: tol_linear = 1.e-6
 real,    parameter :: tol_fixed  = 1.e-6
 real,    parameter :: tol_exact  = 1.e-6
 real,    parameter :: a_diag(ndim) = (/2., 1.5, 0.5/)
 real,    parameter :: b_diag(ndim) = (/1.0, 1.2, 0.8/)
 real,    parameter :: x_ref(ndim)  = (/0.2, -0.3, 0.5/)
 real,    parameter :: c_diag(ndim) = -a_diag*x_ref**4 - b_diag*x_ref
 real,    parameter :: rhs(ndim)    = -c_diag
 real,    parameter :: x0(ndim)     = (/0.05, -0.05, 0.1/)
 real :: diag_work(ndim), x_linear(ndim)
 real :: x_case1(ndim), x_case2(ndim), x_iter(ndim), x_new(ndim), x_exact(ndim)
 real :: coeffs(0:3), soln_mag, sgn, uold_local
 integer :: ierr_case, ierr_quartic(ndim), iter, diff
 logical :: converged, moresweep
 integer :: nfail_case1(2), nfail_case2(3), nfail_case3(3)
 procedure(get_Ax_interface), pointer :: get_Ax_routine

 if (id==master) write(*,"(/,a)") '--> checking BiCGSTAB quartic diagonal cases'

 get_Ax_routine => apply_diag

 ! Case 1: keep matrix built from old x inside get_Ax
 diag_work = a_diag*x0**3 + b_diag
 x_linear = rhs/diag_work
 x_case1 = 0.
 ierr_case = 0
 nfail_case1 = 0
 call solve_bicgstab(ndim,rhs,get_Ax_routine,x_case1,ierr_case)
 call checkval(ierr_case,0,0,nfail_case1(1),'quartic diag (frozen A) ierr')
 call checkval(ndim,x_case1,x_linear,tol_linear,nfail_case1(2), &
      'quartic diag (frozen A) solution')
 call update_test_scores(ntests,nfail_case1,npass)

 ! Case 2: update matrix between outer fixed-point iterations
 nfail_case2 = 0
 converged = .false.
 x_iter = x0
 x_new = x_iter
 ierr_case = 0
 do iter = 1,max_outer
    diag_work = a_diag*x_iter**3 + b_diag
    call solve_bicgstab(ndim,rhs,get_Ax_routine,x_new,ierr_case)
    if (ierr_case /= 0) exit
    if (maxval(abs(x_new - x_iter)) < tol_fixed) then
       converged = .true.
       x_iter = x_new
       exit
    endif
    x_iter = x_new
 enddo
 x_case2 = x_iter
 call checkval(ierr_case,0,0,nfail_case2(1),'quartic diag (updated A) ierr')
 call checkval(ndim,x_case2,x_ref,tol_fixed,nfail_case2(2), &
      'quartic diag (updated A) solution')
 call checkval(converged,.true.,nfail_case2(3), &
      'quartic diag (updated A) converged')
 call update_test_scores(ntests,nfail_case2,npass)

 ! Case 3: analytic quartic solution on diagonal
 nfail_case3 = 0
 do iter = 1,ndim
    coeffs(0) = c_diag(iter)/a_diag(iter)
    coeffs(2) = 0.
    coeffs(3) = 0.
    sgn = sign(1.0,x_case2(iter))
    if (abs(sgn) <= 0.) sgn = 1.
    coeffs(1) = sgn * b_diag(iter)/a_diag(iter)
    moresweep = .false.
    uold_local = max(abs(x_case2(iter)),1.e-10)
    call quarticsolve(coeffs,uold_local,soln_mag,moresweep,ierr_quartic(iter))
    x_exact(iter) = sgn * soln_mag
    call checkval(ierr_quartic(iter),0,0,diff,'quartic diag analytic ierr')
    nfail_case3(1) = nfail_case3(1) + diff
 enddo

 call checkval(ndim,x_exact,x_ref,tol_exact,nfail_case3(2), &
      'quartic diag analytic solution')
 call checkval(ndim,x_case2,x_exact,tol_exact,nfail_case3(3), &
      'quartic diag iter vs analytic')
 call update_test_scores(ntests,nfail_case3,npass)

 if (id==master) write(*,"(a)") '<-- BiCGSTAB quartic diagonal cases complete'

contains

function apply_diag(n_local,x_local) result(Ax)
 integer, intent(in) :: n_local
 real, intent(in) :: x_local(n_local)
 real :: Ax(n_local)

 Ax = diag_work(1:n_local)*x_local

end function apply_diag

end subroutine test_bicgstab_quartic

end module testlinalg
