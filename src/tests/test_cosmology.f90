!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module testcosmology
!
! Unit tests of the Fourier transforms and the Gaussian random field
! utilities used for cosmological initial conditions
!
! :References: None
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: io, set_grf, testutils, utils_dft
!
 use testutils, only:checkval,update_test_scores
 use io,        only:id,master
 implicit none
 public :: test_cosmology

 real, parameter :: pi = acos(-1.)

 private

contains

!--------------------------------------------
!+
!  Unit tests of the FFT and random fields
!+
!--------------------------------------------
subroutine test_cosmology(ntests,npass)
 integer, intent(inout) :: ntests,npass

 if (id==master) write(*,"(/,a)") '--> TESTING FFT AND COSMOLOGICAL INITIAL CONDITIONS'
 call test_fft(ntests,npass)
 call test_helmholtz(ntests,npass)
 call test_white_noise(ntests,npass)
 call test_zeldovich(ntests,npass)
 if (id==master) write(*,"(/,a)") '<-- COSMOLOGY TESTS COMPLETE'

end subroutine test_cosmology

!--------------------------------------------
!+
!  radix-2 FFT agrees with the direct DFT,
!  and the inverse recovers the input
!+
!--------------------------------------------
subroutine test_fft(ntests,npass)
 use utils_dft, only:dft_1d
 integer, intent(inout) :: ntests,npass
 integer, parameter :: n = 64
 complex :: x(n),y(n),z(n)
 integer :: j,m,nfail(2)

 if (id==master) write(*,"(/,a)") '--> checking FFT against direct DFT'
 ! deterministic, irregular input
 do j=1,n
    x(j) = cmplx(sin(1.3*j) + 0.2*j/n,cos(0.7*j*j),kind=kind(x))
 enddo
 y = x
 call dft_1d(y,-1)
 do m=1,n
    z(m) = 0.
    do j=1,n
       z(m) = z(m) + x(j)*exp(cmplx(0.,-2.*pi*(j-1)*(m-1)/real(n),kind=kind(x)))
    enddo
 enddo
 call checkval(maxval(abs(y - z))/maxval(abs(z)),0.,epsilon(0.)*34.,nfail(1),'fft vs direct dft')
 call dft_1d(y,1)
 call checkval(maxval(abs(y/n - x)),0.,epsilon(0.)*3.1,nfail(2),'fft inverse')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_fft

!--------------------------------------------
!+
!  periodic Helmholtz solve is exact for a
!  Fourier mode
!+
!--------------------------------------------
subroutine test_helmholtz(ntests,npass)
 use utils_dft, only:solve_helmholtz_periodic
 integer, intent(inout) :: ntests,npass
 integer, parameter :: n = 16
 real    :: f(n,n,n),phi(n,n,n),ex(n,n,n),dx(3),c,x,y,z
 integer :: i,j,k,nfail(1)

 if (id==master) write(*,"(/,a)") '--> checking periodic Helmholtz solver'
 dx = 2.*pi/n
 c = 0.7
 do k=1,n
    do j=1,n
       do i=1,n
          x = (i-1)*dx(1)
          y = (j-1)*dx(2)
          z = (k-1)*dx(3)
          ex(i,j,k) = sin(2.*x)*cos(3.*y)*sin(z)
          f(i,j,k)  = -(14. + c)*ex(i,j,k)
       enddo
    enddo
 enddo
 call solve_helmholtz_periodic(f,dx,c,phi)
 call checkval(maxval(abs(phi - ex)),0.,epsilon(0.)*12.,nfail(1),'helmholtz solution')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_helmholtz

!--------------------------------------------
!+
!  white noise has unit variance, and the
!  large scales do not depend on the grid
!  when drawn on the same reference grid
!+
!--------------------------------------------
subroutine test_white_noise(ntests,npass)
 use set_grf, only:white_noise_k
 integer, intent(inout) :: ntests,npass
 integer, parameter :: nref = 8, n1 = 16, n2 = 32
 complex :: w1(n1,n1,n1),w2(n2,n2,n2)
 real    :: var,diff
 integer :: iseed,nfail(2)

 if (id==master) write(*,"(/,a)") '--> checking white noise'
 iseed = -4
 call white_noise_k(n1,n1,iseed,w1)
 var = sum(abs(w1)**2)/real(n1)**6
 call checkval(var,1.,8.5e-3,nfail(1),'white noise variance')
 iseed = -4
 call white_noise_k(n1,nref,iseed,w1)
 iseed = -4
 call white_noise_k(n2,nref,iseed,w2)
 ! mode (2,3,1) of each grid, scaled by the grid normalisation
 diff = abs(w1(3,4,2)/sqrt(real(n1)**3) - w2(3,4,2)/sqrt(real(n2)**3))
 call checkval(diff,0.,epsilon(0.)*2.,nfail(2),'same large scales on different grids')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_white_noise

!--------------------------------------------
!+
!  particles moved by the Zel'dovich
!  displacement trace the density contrast:
!  div psi = -delta for a sine mode
!+
!--------------------------------------------
subroutine test_zeldovich(ntests,npass)
 use set_grf, only:zeldovich_displacement
 integer, intent(inout) :: ntests,npass
 integer, parameter :: n = 16
 real    :: delta(n,n,n),psi(3,n,n,n),dx,x,ampl,psiex,err
 integer :: i,nfail(1)

 if (id==master) write(*,"(/,a)") '--> checking Zeldovich displacement'
 dx = 1./n
 ampl = 1.e-3
 do i=1,n
    delta(i,:,:) = ampl*sin(2.*pi*(i-1)*dx)
 enddo
 call zeldovich_displacement(n,dx,delta,psi)
 ! psi_x = (ampl/k) cos(k x) for delta = ampl sin(k x)
 err = 0.
 do i=1,n
    x = (i-1)*dx
    psiex = ampl/(2.*pi)*cos(2.*pi*x)
    err = max(err,maxval(abs(psi(1,i,:,:) - psiex)) + maxval(abs(psi(2:3,i,:,:))))
 enddo
 call checkval(err/ampl,0.,epsilon(0.),nfail(1),'zeldovich displacement')
 call update_test_scores(ntests,nfail,npass)

end subroutine test_zeldovich

end module testcosmology
