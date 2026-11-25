!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module pcg
!
! Matrix-free preconditioned conjugate gradient solver for a symmetric
! positive (semi-)definite system A x = b. The operator, preconditioner
! and stopping test are supplied by the caller as procedures, so the
! solver has no knowledge of particles or neighbour lists.
!
! :References:
!   Shewchuk (1994), An Introduction to the Conjugate Gradient Method
!   Without the Agonizing Pain
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
 public :: pcg_solve

 ! error codes returned by pcg_solve
 integer, parameter, public :: ierr_pcg_breakdown = 1

 abstract interface
  !--y = A x, or z = P^{-1} r for the preconditioner
  subroutine pcg_operator(n,x,y)
   integer, intent(in)  :: n
   real,    intent(in)  :: x(n)
   real,    intent(out) :: y(n)
  end subroutine pcg_operator
  !--true when the residual r is small enough to stop
  logical function pcg_stop(n,r)
   integer, intent(in) :: n
   real,    intent(in) :: r(n)
  end function pcg_stop
 end interface

 private

contains

!-----------------------------------------------------------------------
!+
!  Solve A x = b with x starting at zero. On entry r holds b; on exit
!  it holds the final residual b - A x. z, p and w are work arrays of
!  length n. The stopping test is applied to the initial residual and
!  after every update. ierr is non-zero only on breakdown, where p.Ap
!  is not positive.
!+
!-----------------------------------------------------------------------
subroutine pcg_solve(n,apply_op,apply_prec,converged,maxit,x,r,z,p,w,niter,ierr)
 integer, intent(in)    :: n,maxit
 procedure(pcg_operator) :: apply_op,apply_prec
 procedure(pcg_stop)     :: converged
 real,    intent(out)   :: x(n)
 real,    intent(inout) :: r(n)
 real,    intent(out)   :: z(n),p(n),w(n)
 integer, intent(out)   :: niter,ierr
 integer :: i
 real    :: gamma,gammanew,alpha,beta,denom

 niter = 0
 ierr  = 0
 x = 0.
 if (converged(n,r)) return

 call apply_prec(n,r,z)
 p = z
 gamma = dot(n,r,z)
 if (gamma <= 0.) return

 do while (niter < maxit)
    niter = niter + 1
    call apply_op(n,p,w)
    denom = dot(n,p,w)
    if (denom <= 0.) then
       ierr  = ierr_pcg_breakdown
       niter = niter - 1
       exit
    endif
    alpha = gamma/denom
    !$omp parallel do default(none) private(i) shared(n,x,r,p,w,alpha)
    do i=1,n
       x(i) = x(i) + alpha*p(i)
       r(i) = r(i) - alpha*w(i)
    enddo
    !$omp end parallel do
    if (converged(n,r)) exit

    call apply_prec(n,r,z)
    gammanew = dot(n,r,z)
    if (gammanew <= 0.) exit
    beta = gammanew/gamma
    !$omp parallel do default(none) private(i) shared(n,p,z,beta)
    do i=1,n
       p(i) = z(i) + beta*p(i)
    enddo
    !$omp end parallel do
    gamma = gammanew
 enddo

end subroutine pcg_solve

!-----------------------------------------------------------------------
!+
!  Euclidean dot product of a(1:n) and b(1:n)
!+
!-----------------------------------------------------------------------
real function dot(n,a,b)
 integer, intent(in) :: n
 real,    intent(in) :: a(n),b(n)
 integer :: i
 real    :: acc

 acc = 0.
 !$omp parallel do default(none) private(i) shared(n,a,b) reduction(+:acc)
 do i=1,n
    acc = acc + a(i)*b(i)
 enddo
 !$omp end parallel do
 dot = acc

end function dot

end module pcg
