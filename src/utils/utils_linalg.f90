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

!-----------------------------------------------------------------------
!+
!  Solve a linear system of equations using the BiCGSTAB algorithm
!+
!-----------------------------------------------------------------------
subroutine solve_bicgstab(n,b,get_Ax,x,ierr)
 integer, intent(in) :: n
 real, intent(in) :: b(n)
 real, intent(out) :: x(n)
 integer, intent(out) :: ierr
 procedure(get_Ax_interface), pointer, intent(in) :: get_Ax
 real, dimension(n) :: Ax,r,r0_tilde,p
 integer :: i
 real :: rho

 ierr = 0
 x = 0. ! initial guess x₀
 Ax = get_Ax(n,x)
 r = b - Ax ! initial residual r₀ = b - Ax₀
 r0_tilde = r ! shadow residual r̃₀ (arbitrary choice, commonly r̃₀ = r₀)
 p = r
 rho = dot_product(r0_tilde,r) ! ρ₀ = 1

 do i = 1,10 ! maximum iterations
    ! Check for breakdown
    if (abs(rho) < 1.e-20) then
       ierr = 2 ! breakdown
       return
    endif
    
    ! Call bicgstab_step to perform one iteration
    call bicgstab_step(n,r0_tilde,x,get_Ax,r,p,rho)

    ! Step 12: Convergence check
    if (sqrt(dot_product(r,r)) < 1.e-12) exit
 enddo

 if (i > 100) ierr = 1 ! convergence failure

end subroutine solve_bicgstab

!--------------------------------------------------------------------------------
!+
!  BiCGSTAB solver (van der Horst, 1992) for inverting A*x = b
!  This implements one step of the BiCGSTAB algorithm
!+
!--------------------------------------------------------------------------------
subroutine bicgstab_step(n,r0_tilde,x,get_Ax_i,r,p,rho)
 integer, intent(in) :: n
 real, intent(in), dimension(n) :: r0_tilde
 real, intent(inout), dimension(n) :: x,r,p
 real, intent(inout) :: rho
 procedure(get_Ax_interface), pointer :: get_Ax_i
 real :: beta,omega,alpha,rho_prev
 real, dimension(n) :: y,z,s,t,v
 
 ! Step 4: Solve y from Ky = pᵢ (preconditioning step)
 ! For now, we'll use y = pᵢ (no preconditioning)
 y = p
 
 ! Step 5: Compute vᵢ = Ay
 v = get_Ax_i(n,y)
 
 ! Step 6: Compute α = ρᵢ / (r̃₀, vᵢ)
 alpha = rho / dot_product(r0_tilde,v)
 
 ! Step 7: Compute s = rᵢ₋₁ - αvᵢ
 s = r - alpha * v
 
 ! Step 8: Solve z from Kz = s (preconditioning step)
 ! For now, we'll use z = s (no preconditioning)
 z = s
 
 ! Step 9: Compute t = Az
 t = get_Ax_i(n,z)
 
 ! Step 10: Compute ωᵢ = (K₁⁻¹t, K₁⁻¹s) / (K₁⁻¹t, K₁⁻¹t)
 ! For now, we'll use ωᵢ = (t, s) / (t, t) (no preconditioning)
 omega = dot_product(t,s) / dot_product(t,t)
 
 ! Step 11: Update xᵢ = xᵢ₋₁ + αy + ωᵢz
 x = x + alpha * y + omega * z
 
 ! Step 13: Update rᵢ = s - ωᵢt
 r = s - omega * t

 ! compute ρᵢ for next iteration
 rho_prev = rho
 rho = dot_product(r0_tilde,r)

 ! Step 2: Compute β = (ρᵢ / ρᵢ₋₁)(α / ωᵢ₋₁)
 beta = (rho/rho_prev) * (alpha/omega)

 ! Step 3: Update pᵢ = rᵢ₋₁ + β(pᵢ₋₁ - ωᵢ₋₁vᵢ₋₁)
 p = r + beta * (p - omega * v)

end subroutine bicgstab_step

end module linalg
