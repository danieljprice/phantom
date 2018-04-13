!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: solvelinearsystem
!
!  DESCRIPTION:
!  module containing routine(s) to solve AX = B
!  using a partial pivoting algorithm and reduced
!  storage.
!
!  REFERENCES: None
!
!  OWNER: Mark Hutchison
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module solvelinearsystem
 implicit none

contains

!***********************************************************************
!* Solve AX = B using a partial pivoting algorithm and reduced storage *
!* ------------------------------------------------------------------- *
!* SAMPLE RUN:                                                         *
!*                                                                     *
!* System to solve:                                                    *
!*   2.0000  -1.0000   1.0000   7.0000 -12.5400    5.0000              *
!*   1.0000   5.0000  -2.0000  -8.0000 100.0000    1.0000              *
!*   3.0000  -2.0000   3.0000  45.0000  27.3333    3.0000              *
!*  11.0000   0.5500  -2.0000  -4.0000   1.0000    4.0000              *
!*  33.0000   2.0000  -3.0000   5.0000   7.3333  -10.0000              *
!*                                                                     *
!* Solution is:                                                        *
!*   2.11149597961869                                                  *
!*  -25.8290267820056                                                  *
!*   8.17423194407132                                                  *
!*  -2.52146730210577                                                  *
!*   1.24210363401706                                                  *
!*                                                                     *
!* ------------------------------------------------------------------- *
!* Ref.: "Wassyng, A. - Solving Ax = b: A method with reduced storage  *
!*        requirements, SIAM J. Numerical Analysis, vol.19 (1982),     *
!*        pp. 197-204".                                                *
!*                                                                     *
!*                                  F90 Release By J-P Moreau, Paris.  *
!*                                         (www.jpmoreau.fr)           *
!***********************************************************************
subroutine dple(rowk, n, a, b, c, ierr)

! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-16  Time: 12:26:32

! ******************************************************************
!        SOLUTION OF LINEAR EQUATIONS WITH REDUCED STORAGE
! ******************************************************************

! Uses the Henderson-Wassyng partial pivot algorithm.
! Wassyng, A. 'Solving Ax = b: A method with reduced storage requirements',
! SIAM J. Numerical Analysis, vol.19 (1982), pp. 197-204.

! The user must provide a routine ROWK to return the requested row of the
! matrix A.

! N.B. Arguments D and IP have been removed.

implicit none
integer, parameter  :: dp = selected_real_kind(14, 60)

integer,  intent(in)    :: n
real(dp), intent(inout) :: a(n,n)
real(dp), intent(in)    :: b(n)
real(dp), intent(out)   :: c(n)
integer,  intent(out)   :: ierr

! external rowk
interface
  subroutine rowk(n, a, k, r)
    implicit none
    integer, parameter  :: dp = selected_real_kind(14, 60)
    integer,   intent(in)    :: n, k
    real(dp),  intent(inout) :: a(n,n)
    real (dp), intent(out)   :: r(:)
  end subroutine rowk
end interface

! Local variables
real (dp)  :: bk, cj, ck, c1, dkj
real (dp), parameter  :: zero = 0.0_dp
real (dp)  :: wk(n*n/4 + n + 3)
integer    :: i, iflag, ij, ijold, ik, iwk(n), j, k, kjold, km1, kp1,   &
              last, lastm1, lcol, lcolp1, m, maxwk, mjold, nm1, np1

! Set the necessary constants
ierr = 0
maxwk = n * n / 4 + n + 3
np1 = n + 1
k = 1
iflag = -1

! Get the first column of the transposed system
call rowk(n, a, 1, c)
bk = b(1)

if (n <= 1) then
  if (c(1) == zero) GO TO 130
  c(1) = bk / c(1)
  return
end if

! Find the pivot for column 1
m = 1
do  i = 2, n
  if (abs(c(m)) < abs(c(i))) m = i
end do

iwk(1) = m
c1 = c(m)
c(m) = c(1)
c(1) = c1
if (c(1) /= zero) then

! Find the first elementary matrix and store it in d
  do  i = 2, n
    wk(i-1) = -c(i) / c(1)
  end do
  wk(n) = bk / c(1)

! k loop - each k for a new column of the transposed system
  do  k = 2, n
    kp1 = k + 1
    km1 = k - 1

! Get column k
    call rowk(n, a, k, c)
    do  j = 1, km1
      m = iwk(j)
      cj = c(j)
      c(j) = c(m)
      c(m) = cj
    end do
    bk = b(k)

    iflag = -iflag
    lcol = np1 - k
    lcolp1 = lcol + 1
    lastm1 = 1
    last = maxwk - n + k
    if (k /= 2) then

      lastm1 = maxwk - n + km1
      if (iflag < 0) last = last - n + k - 2
      if (iflag > 0) lastm1 = lastm1 - n + k - 3
    end if

! j loop - effect of columns 1 to k-1 of l-inverse
    do  j = 1, km1
      cj = c(j)
      ij = (j-1) * lcolp1
      if (j == km1) ij = lastm1 - 1

! i loop - effect of l-inverse on rows k to n+1
      do  i = k, n
        ij = ij + 1
        c(i) = c(i) + wk(ij) * cj
      end do
      bk = bk - wk(ij+1) * cj
    end do

! k=n case
    m = k
    if (k >= n) then
      if (c(k) == zero) GO TO 130
      wk(last) = bk / c(k)
    else

! Find the pivot
      do  i = kp1, n
        if (abs(c(m)) < abs(c(i))) m = i
      end do

      iwk(k) = m
      ck = c(m)
      c(m) = c(k)
      c(k) = ck
      if (c(k) == zero) GO TO 130

! Find the k-th elementary matrix
      ik = last
      do  i = kp1, n
        wk(ik) = -c(i) / c(k)
        ik = ik + 1
      end do
      wk(ik) = bk / c(k)
    end if

! Form the product of the elementary matrices
    do  j = 1, km1
      kjold = j * lcolp1 + k - np1
      mjold = kjold + m - k
      ij = (j-1) * lcol
      ijold = ij + j
      if (j == km1) then

        kjold = lastm1
        mjold = lastm1 + m - k
        ijold = lastm1
      end if

      ik = last - 1
      dkj = wk(mjold)
      wk(mjold) = wk(kjold)
      do  i = kp1, np1
        ij = ij + 1
        ijold = ijold + 1
        ik = ik + 1
        wk(ij) = wk(ijold) + wk(ik) * dkj
      end do
    end do
  end do

  last = maxwk
  if (iflag < 0) last = maxwk - 2
  wk(n) = wk(last)

! Insert the solution in c

  c(1:n) = wk(1:n)

  nm1 = n - 1
  do  i = 1, nm1
    k = n - i
    m = iwk(k)
    ck = c(k)
    c(k) = c(m)
    c(m) = ck
  end do
  return
end if

! The system is singular
130 ierr = k

return
end subroutine dple


subroutine rowk(n, a, k, r)
 integer, parameter  :: dp = selected_real_kind(14, 60)
 integer,   intent(in)    :: n, k
 real(dp),  intent(inout) :: a(n,n)
 real (dp), intent(out)   :: r(:)

 r(:) = a(k,:)

 return
end subroutine rowk


end module solvelinearsystem

