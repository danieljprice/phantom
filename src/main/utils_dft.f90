!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module utils_dft
!
! Discrete Fourier transforms on uniform periodic grids, and a spectral
! solver for the periodic Helmholtz (or Poisson) equation. Lines whose
! length is a power of two use a radix-2 fast Fourier transform, costing
! O(n^3 log n) on an n^3 grid; other lengths fall back to a direct DFT
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 ! pi to full precision (physcon:pi is truncated at 11 digits)
 real, parameter :: pi = acos(-1.)

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
 !$omp parallel do default(none) shared(n,dx,c,fk) private(i,j,k,k2,denom)
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
 !$omp end parallel do
 call dft_3d(fk,1)
 phi = real(fk)/product(n)

end subroutine solve_helmholtz_periodic

!-----------------------------------------------------------------------
!+
!  wavenumber of Fourier index i (1..n) on a periodic grid of spacing dx
!+
!-----------------------------------------------------------------------
pure real function wavenumber(i,n,dx)
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
!  array, with exponent sign isign*2 pi i jk/n. Lines are copied to a
!  contiguous buffer and transformed in parallel
!+
!-----------------------------------------------------------------------
subroutine dft_3d(a,isign)
 complex, intent(inout) :: a(:,:,:)
 integer, intent(in)    :: isign
 complex, allocatable :: line(:)
 integer :: i,j,k,n(3)

 n = shape(a)
 !$omp parallel default(none) shared(a,n,isign) private(i,j,k,line)
 allocate(line(maxval(n)))
 !$omp do collapse(2)
 do k=1,n(3)
    do j=1,n(2)
       line(1:n(1)) = a(:,j,k)
       call dft_1d(line(1:n(1)),isign)
       a(:,j,k) = line(1:n(1))
    enddo
 enddo
 !$omp end do
 !$omp do collapse(2)
 do k=1,n(3)
    do i=1,n(1)
       line(1:n(2)) = a(i,:,k)
       call dft_1d(line(1:n(2)),isign)
       a(i,:,k) = line(1:n(2))
    enddo
 enddo
 !$omp end do
 !$omp do collapse(2)
 do j=1,n(2)
    do i=1,n(1)
       line(1:n(3)) = a(i,j,:)
       call dft_1d(line(1:n(3)),isign)
       a(i,j,:) = line(1:n(3))
    enddo
 enddo
 !$omp end do
 deallocate(line)
 !$omp end parallel

end subroutine dft_3d

!-----------------------------------------------------------------------
!+
!  unnormalised discrete Fourier transform of one line, by FFT if the
!  length is a power of two, otherwise directly
!+
!-----------------------------------------------------------------------
subroutine dft_1d(x,isign)
 complex, intent(inout) :: x(:)
 integer, intent(in)    :: isign
 integer :: n

 n = size(x)
 if (n > 1 .and. iand(n,n-1) == 0) then
    call fft_radix2(x,isign)
 else
    call dft_direct(x,isign)
 endif

end subroutine dft_1d

!-----------------------------------------------------------------------
!+
!  direct (O(n^2)) discrete Fourier transform of one line
!+
!-----------------------------------------------------------------------
subroutine dft_direct(x,isign)
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

end subroutine dft_direct

!-----------------------------------------------------------------------
!+
!  in-place iterative radix-2 Cooley-Tukey FFT of a line whose length
!  is a power of two: bit-reversal permutation, then butterflies of
!  doubling length. Twiddle factors are evaluated directly (not by
!  recurrence) to keep round-off at the level of a single exp
!+
!-----------------------------------------------------------------------
subroutine fft_radix2(x,isign)
 complex, intent(inout) :: x(:)
 integer, intent(in)    :: isign
 complex :: w,t
 integer :: n,i,j,m,len,half,k

 n = size(x)
 ! bit-reversal permutation
 j = 1
 do i=1,n
    if (i < j) then
       t = x(j)
       x(j) = x(i)
       x(i) = t
    endif
    m = n/2
    do while (m >= 1 .and. j > m)
       j = j - m
       m = m/2
    enddo
    j = j + m
 enddo
 ! butterflies
 len = 2
 do while (len <= n)
    half = len/2
    do k=0,half-1
       w = exp(cmplx(0.,isign*2.*pi*k/real(len),kind=kind(x)))
       do i=1,n,len
          t = w*x(i+k+half)
          x(i+k+half) = x(i+k) - t
          x(i+k) = x(i+k) + t
       enddo
    enddo
    len = 2*len
 enddo

end subroutine fft_radix2

end module utils_dft
