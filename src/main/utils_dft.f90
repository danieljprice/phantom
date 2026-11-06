!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_dft
!
! Discrete Fourier transforms on small uniform periodic grids, and a
! spectral solver for the periodic Helmholtz (or Poisson) equation.
! The transforms are direct (slow) DFTs, cost O(n^4) on an n^3 grid, so
! are intended for small grids (e.g. setting up initial data)
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon
!
 implicit none

 public :: solve_helmholtz_periodic,dft_3d,dft_1d,wavenumber

 private

contains

!-----------------------------------------------------------------------
!+
!  solve del^2 phi - c phi = f with periodic boundary conditions on a
!  uniform grid of spacing dx, by discrete Fourier transform (the zero
!  mode of phi is zero if c = 0)
!+
!-----------------------------------------------------------------------
subroutine solve_helmholtz_periodic(f,dx,c,phi)
 real, intent(in)  :: f(:,:,:),dx(3),c
 real, intent(out) :: phi(:,:,:)
 complex, allocatable :: fk(:,:,:)
 real :: k2,denom
 integer :: n(3),i,j,k

 n = shape(f)
 allocate(fk(n(1),n(2),n(3)))
 fk = cmplx(f,0.,kind=kind(f))
 call dft_3d(fk,-1)
 do k=1,n(3)
    do j=1,n(2)
       do i=1,n(1)
          k2 = wavenumber(i,n(1),dx(1))**2 + wavenumber(j,n(2),dx(2))**2 + wavenumber(k,n(3),dx(3))**2
          denom = -k2 - c
          if (abs(denom) > tiny(denom)) then
             fk(i,j,k) = fk(i,j,k)/denom
          else
             fk(i,j,k) = 0.
          endif
       enddo
    enddo
 enddo
 call dft_3d(fk,1)
 phi = real(fk)/product(n)

end subroutine solve_helmholtz_periodic

!-----------------------------------------------------------------------
!+
!  wavenumber of Fourier index i (1..n) on a periodic grid of spacing dx
!+
!-----------------------------------------------------------------------
pure real function wavenumber(i,n,dx)
 use physcon, only:pi
 integer, intent(in) :: i,n
 real,    intent(in) :: dx
 integer :: m

 m = i - 1
 if (m > n/2) m = m - n
 wavenumber = 2.*pi*m/(n*dx)

end function wavenumber

!-----------------------------------------------------------------------
!+
!  unnormalised discrete Fourier transform along each dimension of a 3D
!  array, with exponent sign isign*2 pi i jk/n (slow DFT, for small grids)
!+
!-----------------------------------------------------------------------
subroutine dft_3d(a,isign)
 complex, intent(inout) :: a(:,:,:)
 integer, intent(in)    :: isign
 integer :: i,j,k,n(3)

 n = shape(a)
 do k=1,n(3)
    do j=1,n(2)
       call dft_1d(a(:,j,k),isign)
    enddo
 enddo
 do k=1,n(3)
    do i=1,n(1)
       call dft_1d(a(i,:,k),isign)
    enddo
 enddo
 do j=1,n(2)
    do i=1,n(1)
       call dft_1d(a(i,j,:),isign)
    enddo
 enddo

end subroutine dft_3d

!-----------------------------------------------------------------------
!+
!  unnormalised discrete Fourier transform of one line
!+
!-----------------------------------------------------------------------
subroutine dft_1d(x,isign)
 use physcon, only:pi
 complex, intent(inout) :: x(:)
 integer, intent(in)    :: isign
 complex :: y(size(x))
 integer :: j,m,n

 n = size(x)
 do m=1,n
    y(m) = 0.
    do j=1,n
       y(m) = y(m) + x(j)*exp(cmplx(0.,isign*2.*pi*(j-1)*(m-1)/real(n),kind=kind(x)))
    enddo
 enddo
 x = y

end subroutine dft_1d

end module utils_dft
