!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module linalg
!
! linalg
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 implicit none

 public :: get_Ax_interface,solve_bicgstab,inverse

 abstract interface
  function get_Ax_interface(n, x) result(Ax)
   integer, intent(in) :: n
   real, intent(in) :: x(n)
   real :: Ax(n)
  end function get_Ax_interface

 end interface

 private

contains
function inverse(matrix,n)
 integer, intent(in) ::n
 real, dimension(n,n), intent(in) :: matrix
 real, dimension(n,2*n) :: a,temp
 integer ::i,j,k
 real :: ratio,divisor
 real, dimension(n,n) :: inverse

 a(:,:) = 0.
 inverse(:,:) = 0.

 !Augmenting Identity Matrix of Order n
 do i=1,n
    do j=1,n
       a(i,j) = matrix(i,j)
    enddo
 enddo

 !we should do a row swap if the elements on diagonal are 0

 do i=1,n
    do j=1,n
       if (i==j) then
          a(i,j+n) = 1.
       endif
    enddo
 enddo

 !we swap the rows if we enounter 0 as diagonal element
 temp = a(:,:)
 do i=1,n
    if (a(i,i)==0.0) then
       do j=1,n
          if (j  /=  i) then
             if (a(j,i)  /=  0.) then
                a(j,:) = temp(i,:)
                a(i,:) = temp(j,:)
                temp = a(:,:)
                exit
             endif
          endif
       enddo
    endif
 enddo

 !Applying Guass Jordan Elimination
 do i=1,n
    do j=1,n
       if (i  /=  j) then
          ratio = a(j,i)/a(i,i)
          do k=1,2*n
             a(j,k) = a(j,k) - ratio*a(i,k)
          enddo
       endif
    enddo
 enddo

 !dividing by the diagonal elements to get identity matrix
 do i=1,n
    divisor = a(i,i)
    do j=1,2*n
       a(i,j) = a(i,j)/divisor
    enddo
 enddo
 do i=1,n
    do j=1,n
       inverse(i,j)=a(i,j+n)
    enddo
 enddo

end function inverse


subroutine solve_bicgstab(n,b,get_Ax,x,ierr)
 integer, intent(in) :: n
 real, intent(in) :: b(n)
 real, intent(out) :: x(n)
 integer, intent(out) :: ierr
 procedure(get_Ax_interface), pointer, intent(in) :: get_Ax
 real, dimension(n) :: Ax,r,rbar,p,pbar
 integer :: i
 real, dimension(n) :: Ap0,r0,rbar0,Apbar0,p0,pbar0
 real :: a0,b0

 x = 0. ! guess
 Ax = get_Ax(n,x)
 r = b-Ax ! residual
 rbar = r
 p = r
 pbar = rbar

 do i = 1,10
    r0 = r
    rbar0 = rbar
    p0 = p
    pbar0 = pbar

    Ap0 = get_Ax(n,p0)
    a0 = dot_product(rbar0,r0) / dot_product(pbar0,Ap0)
    r = r0 - a0*Ap0
    Apbar0 = get_Ax(n,pbar0)
    rbar = rbar0 - a0*Apbar0
    b0 = dot_product(rbar,r) / dot_product(rbar0,r0)
    p = r + b0*p0
    pbar = rbar + b0*pbar0
    x = x + a0*p0
 enddo


end subroutine solve_bicgstab

!--------------------------------------------------------------------------------
!+
!  BiCGSTAB solver (van der Horst, 1992) for inverting A*x = b
!+
!--------------------------------------------------------------------------------
subroutine bicgstab_step(n,x,get_Ax_i,r,rbar,p,pbar)
 integer, intent(in) :: n
 real, intent(inout), dimension(n) :: x,r,rbar,p,pbar
 procedure(get_Ax_interface), pointer :: get_Ax_i
 real :: a0,b0
 real, dimension(n) :: Ap0,r0,rbar0,Apbar0,p0,pbar0

 r0 = r
 rbar0 = rbar
 p0 = p
 pbar0 = pbar

 Ap0 = get_Ax_i(n,p0)
 a0 = dot_product(rbar0,r0) / dot_product(pbar0,Ap0)
 r = r0 - a0*Ap0
 Apbar0 = get_Ax_i(n,pbar0)
 rbar = rbar0 - a0*Apbar0
 b0 = dot_product(rbar,r) / dot_product(rbar0,r0)
 p = r + b0*p0
 pbar = rbar + b0*pbar0
 x = x + a0*p0

end subroutine bicgstab_step

end module linalg
