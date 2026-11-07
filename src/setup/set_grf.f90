!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module set_grf
!
! Gaussian random field initial conditions for cosmology on a uniform
! periodic grid: the matter power spectrum (from a table, as written by
! CAMB or CLASS, or the Eisenstein & Hu fit), the linear growing mode of
! an Einstein-de Sitter universe in longitudinal gauge, and Zel'dovich
! displacements of particles from a lattice.
!
! Following Macpherson, Price & Lasky (2019) and FLRWSolver
! (github.com/hayleyjm/FLRWSolver_public), with a = 1 initially:
!   phi(k)   = -delta_sync(k)/(C1 k^2)
!   delta(k) = -(C1 k^2 + 2) phi(k)
!   v^i(k)   =  C3 i k^i phi(k)
! with C1 = 2/(3 H^2), C3 = -2/(3 H), where delta_sync is the density
! contrast in synchronous gauge (the usual output of CAMB/CLASS) and
! H is the initial conformal Hubble rate in code units (c = G = 1).
!
! :References:
!   Macpherson, Lasky & Price (2017), Phys. Rev. D 95, 064028
!   Macpherson, Price & Lasky (2019), Phys. Rev. D 99, 063522
!   Eisenstein & Hu (1998), ApJ 496, 605
!
! :Owner:
!
! :Runtime parameters: None
!
! :Dependencies: random, utils_dft
!
 implicit none

 real, parameter :: pi = acos(-1.)

 public :: read_power_spectrum,pk_table,pk_eisenstein_hu
 public :: white_noise_k,delta_sync_k,phi_k_from_delta_sync
 public :: linear_fields_from_phi_k,zeldovich_displacement

 private

contains

!-----------------------------------------------------------------------
!+
!  read a power spectrum table of k [h/Mpc] and P(k) [(Mpc/h)^3] in two
!  columns, skipping lines that do not parse (e.g. comment headers)
!+
!-----------------------------------------------------------------------
subroutine read_power_spectrum(filename,ktab,pktab,ierr)
 character(len=*),  intent(in)  :: filename
 real, allocatable, intent(out) :: ktab(:),pktab(:)
 integer,           intent(out) :: ierr
 real    :: kval,pval
 integer :: iunit,ios,n,pass

 open(newunit=iunit,file=filename,status='old',action='read',iostat=ierr)
 if (ierr /= 0) return
 ! first pass counts the entries, the second stores them
 do pass=1,2
    n = 0
    rewind(iunit)
    do
       read(iunit,'(a)',iostat=ios)
       if (ios /= 0) exit
       backspace(iunit)
       read(iunit,*,iostat=ios) kval,pval
       if (ios /= 0) cycle
       n = n + 1
       if (pass==2) then
          ktab(n)  = kval
          pktab(n) = pval
       endif
    enddo
    if (pass==1) allocate(ktab(n),pktab(n))
 enddo
 close(iunit)
 if (n < 2) ierr = 1

end subroutine read_power_spectrum

!-----------------------------------------------------------------------
!+
!  P(k) from a table by linear interpolation in log k - log P, zero
!  outside the tabulated range
!+
!-----------------------------------------------------------------------
pure real function pk_table(k,ktab,pktab)
 real, intent(in) :: k,ktab(:),pktab(:)
 integer :: ilo,ihi,imid
 real    :: f

 pk_table = 0.
 if (k < ktab(1) .or. k > ktab(size(ktab))) return
 ! bisection for the bracketing interval
 ilo = 1
 ihi = size(ktab)
 do while (ihi - ilo > 1)
    imid = (ilo + ihi)/2
    if (ktab(imid) > k) then
       ihi = imid
    else
       ilo = imid
    endif
 enddo
 if (pktab(ilo) <= 0. .or. pktab(ihi) <= 0.) then
    f = (k - ktab(ilo))/(ktab(ihi) - ktab(ilo))
    pk_table = (1. - f)*pktab(ilo) + f*pktab(ihi)
 else
    f = log(k/ktab(ilo))/log(ktab(ihi)/ktab(ilo))
    pk_table = exp((1. - f)*log(pktab(ilo)) + f*log(pktab(ihi)))
 endif

end function pk_table

!-----------------------------------------------------------------------
!+
!  linear matter power spectrum P(k) [(Mpc/h)^3] at k [h/Mpc] from the
!  Eisenstein & Hu (1998) zero-baryon-oscillation transfer function
!  (their eqs. 26, 28-31), P = A k^ns T(k)^2, normalised to sigma8 at
!  z = 0 and scaled back to redshift z with the Einstein-de Sitter
!  growth factor D = 1/(1+z). Intended as a dependency-free fallback;
!  use a CAMB/CLASS table for production runs
!+
!-----------------------------------------------------------------------
real function pk_eisenstein_hu(k,z,h,omegam,omegab,ns,sigma8)
 real, intent(in) :: k,z,h,omegam,omegab,ns,sigma8
 real :: amp

 amp = (sigma8/sigma_r(8.,h,omegam,omegab,ns))**2
 pk_eisenstein_hu = amp*pk_eh_unnorm(k,h,omegam,omegab,ns)/(1. + z)**2

end function pk_eisenstein_hu

!-----------------------------------------------------------------------
!+
!  unnormalised k^ns T(k)^2 with the Eisenstein & Hu (1998) no-wiggle
!  transfer function, k in h/Mpc
!+
!-----------------------------------------------------------------------
pure real function pk_eh_unnorm(k,h,omegam,omegab,ns)
 real, intent(in) :: k,h,omegam,omegab,ns
 real, parameter :: theta27 = 2.7255/2.7
 real :: om2,fb,s,alphag,gameff,q,l0,c0,kmpc

 om2  = omegam*h**2
 fb   = omegab/omegam
 kmpc = k*h
 ! sound horizon (eq. 26) and the shape parameter (eqs. 30-31)
 s      = 44.5*log(9.83/om2)/sqrt(1. + 10.*(omegab*h**2)**0.75)
 alphag = 1. - 0.328*log(431.*om2)*fb + 0.38*log(22.3*om2)*fb**2
 gameff = omegam*h*(alphag + (1. - alphag)/(1. + (0.43*kmpc*s)**4))
 q      = k*theta27**2/gameff
 ! eqs. 28-29
 l0 = log(2.*exp(1.) + 1.8*q)
 c0 = 14.2 + 731./(1. + 62.5*q)
 pk_eh_unnorm = k**ns*(l0/(l0 + c0*q**2))**2

end function pk_eh_unnorm

!-----------------------------------------------------------------------
!+
!  rms density contrast in spheres of radius r [Mpc/h] of the
!  unnormalised spectrum, by integration in log k
!+
!-----------------------------------------------------------------------
real function sigma_r(r,h,omegam,omegab,ns)
 real, intent(in) :: r,h,omegam,omegab,ns
 integer, parameter :: nk = 4000
 real, parameter :: lnkmin = log(1.e-5), lnkmax = log(1.e3)
 real :: dlnk,k,x,w,total
 integer :: i

 dlnk = (lnkmax - lnkmin)/nk
 total = 0.
 do i=1,nk
    k = exp(lnkmin + (i - 0.5)*dlnk)
    x = k*r
    ! Fourier transform of the top-hat window
    w = 3.*(sin(x) - x*cos(x))/x**3
    total = total + k**3*pk_eh_unnorm(k,h,omegam,omegab,ns)*w**2*dlnk
 enddo
 sigma_r = sqrt(total/(2.*pi**2))

end function sigma_r

!-----------------------------------------------------------------------
!+
!  Fourier transform of real Gaussian white noise of unit variance on an
!  n^3 grid, so that <|w_k|^2> = n^3 and w_{-k} = conj(w_k), with the zero
!  mode removed. The noise is drawn on a grid of nref^3 points (nref <= n);
!  if nref < n its modes below the Nyquist frequency of that grid are
!  embedded in the n^3 grid, so that runs at different resolutions with the
!  same seed and nref share the same large-scale field
!+
!-----------------------------------------------------------------------
subroutine white_noise_k(n,nref,iseed,wk)
 use random,    only:gauss_random
 use utils_dft, only:dft_3d
 integer, intent(in)    :: n,nref
 integer, intent(inout) :: iseed
 complex, intent(out)   :: wk(n,n,n)
 complex, allocatable   :: wref(:,:,:)
 integer :: i,j,k

 allocate(wref(nref,nref,nref))
 ! serial so that the field depends only on the seed
 do k=1,nref
    do j=1,nref
       do i=1,nref
          wref(i,j,k) = cmplx(gauss_random(iseed),0.,kind=kind(wk))
       enddo
    enddo
 enddo
 call dft_3d(wref,-1)
 wref(1,1,1) = 0.
 if (nref == n) then
    wk = wref
    return
 endif
 ! embed in the larger grid, without the Nyquist modes of the reference
 ! grid (whose sign of k is ambiguous)
 wk = 0.
 do k=1,nref
    do j=1,nref
       do i=1,nref
          if (is_resolved(i,nref) .and. is_resolved(j,nref) .and. is_resolved(k,nref)) &
             wk(fine_index(i,nref,n),fine_index(j,nref,n),fine_index(k,nref,n)) = wref(i,j,k)
       enddo
    enddo
 enddo
 ! the same amplitude per mode for the unnormalised DFT on n^3 points
 wk = wk*sqrt((real(n)/real(nref))**3)

end subroutine white_noise_k

!-----------------------------------------------------------------------
!+
!  true unless Fourier index i (1..n) is the Nyquist mode of an even n
!+
!-----------------------------------------------------------------------
pure logical function is_resolved(i,n)
 integer, intent(in) :: i,n

 is_resolved = (n == 1 .or. 2*(i-1) /= n)

end function is_resolved

!-----------------------------------------------------------------------
!+
!  index on an n-point grid of the mode with index i on an nref-point grid
!+
!-----------------------------------------------------------------------
pure integer function fine_index(i,nref,n)
 integer, intent(in) :: i,nref,n
 integer :: m

 m = i - 1
 if (m > nref/2) m = m - nref
 fine_index = modulo(m,n) + 1

end function fine_index

!-----------------------------------------------------------------------
!+
!  synchronous-gauge density contrast in Fourier space (unnormalised
!  DFT convention) with power spectrum P(k), from white noise wk.
!  lbox is the comoving box length in Mpc/h; modes with wavelength
!  shorter than lambda_min [Mpc/h] are removed (none if lambda_min <= 0)
!+
!-----------------------------------------------------------------------
subroutine delta_sync_k(n,lbox,wk,ktab,pktab,lambda_min,deltak)
 use utils_dft, only:wavenumber
 integer, intent(in)  :: n
 real,    intent(in)  :: lbox,ktab(:),pktab(:),lambda_min
 complex, intent(in)  :: wk(n,n,n)
 complex, intent(out) :: deltak(n,n,n)
 real    :: dx,kmag,kcut,dvol,scalefac
 integer :: i,j,k

 dx = lbox/n
 dvol = dx**3
 ! |delta_k|^2 = P(k) V/dV^2 for the unnormalised DFT; wk carries n^3
 scalefac = (lbox**3/dvol**2)/real(n)**3
 kcut = huge(kcut)
 if (lambda_min > 0.) kcut = 2.*pi/lambda_min
 !$omp parallel do default(none) shared(n,dx,wk,ktab,pktab,kcut,scalefac,deltak) &
 !$omp private(i,j,k,kmag)
 do k=1,n
    do j=1,n
       do i=1,n
          kmag = sqrt(wavenumber(i,n,dx)**2 + wavenumber(j,n,dx)**2 + wavenumber(k,n,dx)**2)
          if (kmag > kcut) then
             deltak(i,j,k) = 0.
          else
             deltak(i,j,k) = wk(i,j,k)*sqrt(pk_table(kmag,ktab,pktab)*scalefac)
          endif
       enddo
    enddo
 enddo
 !$omp end parallel do

end subroutine delta_sync_k

!-----------------------------------------------------------------------
!+
!  growing-mode potential phi(k) = -delta_sync(k)/(C1 k^2), with
!  C1 = 2/(3 hub^2), on a grid of spacing dx in code units
!+
!-----------------------------------------------------------------------
subroutine phi_k_from_delta_sync(n,dx,hub,deltak,phik)
 use utils_dft, only:wavenumber
 integer, intent(in)  :: n
 real,    intent(in)  :: dx,hub
 complex, intent(in)  :: deltak(n,n,n)
 complex, intent(out) :: phik(n,n,n)
 real    :: k2,c1
 integer :: i,j,k

 c1 = 2./(3.*hub**2)
 !$omp parallel do default(none) shared(n,dx,c1,deltak,phik) private(i,j,k,k2)
 do k=1,n
    do j=1,n
       do i=1,n
          k2 = wavenumber(i,n,dx)**2 + wavenumber(j,n,dx)**2 + wavenumber(k,n,dx)**2
          if (k2 > 0.) then
             phik(i,j,k) = -deltak(i,j,k)/(c1*k2)
          else
             phik(i,j,k) = 0.
          endif
       enddo
    enddo
 enddo
 !$omp end parallel do

end subroutine phi_k_from_delta_sync

!-----------------------------------------------------------------------
!+
!  real-space linear growing mode from phi(k): the potential phi, the
!  longitudinal-gauge density contrast delta = -(C1 k^2 + 2) phi and the
!  coordinate velocity v^i = C3 i k^i phi, C3 = -2/(3 hub)
!+
!-----------------------------------------------------------------------
subroutine linear_fields_from_phi_k(n,dx,hub,phik,phi,delta,vel)
 use utils_dft, only:wavenumber
 integer, intent(in)  :: n
 real,    intent(in)  :: dx,hub
 complex, intent(in)  :: phik(n,n,n)
 real,    intent(out) :: phi(n,n,n),delta(n,n,n),vel(3,n,n,n)
 complex, allocatable :: fk(:,:,:)
 real    :: kvec(3),c1,c3
 integer :: i,j,k,idim

 c1 = 2./(3.*hub**2)
 c3 = -2./(3.*hub)
 allocate(fk(n,n,n))
 call inverse_transform(phik,phi)
 !$omp parallel do default(none) shared(n,dx,c1,phik,fk) private(i,j,k,kvec)
 do k=1,n
    do j=1,n
       do i=1,n
          kvec = (/wavenumber(i,n,dx),wavenumber(j,n,dx),wavenumber(k,n,dx)/)
          fk(i,j,k) = -(c1*dot_product(kvec,kvec) + 2.)*phik(i,j,k)
       enddo
    enddo
 enddo
 !$omp end parallel do
 call inverse_transform(fk,delta)
 do idim=1,3
    call gradient_k(n,dx,phik,idim,fk)
    call inverse_transform(fk,vel(idim,:,:,:))
    vel(idim,:,:,:) = c3*vel(idim,:,:,:)
 enddo

end subroutine linear_fields_from_phi_k

!-----------------------------------------------------------------------
!+
!  Zel'dovich displacement psi = -grad (del^-2 deltastar) of a periodic
!  density contrast, so that particles moved from a uniform lattice by
!  psi have density contrast deltastar to first order (div psi = -delta)
!+
!-----------------------------------------------------------------------
subroutine zeldovich_displacement(n,dx,deltastar,psi)
 use utils_dft, only:dft_3d,wavenumber
 integer, intent(in)  :: n
 real,    intent(in)  :: dx,deltastar(n,n,n)
 real,    intent(out) :: psi(3,n,n,n)
 complex, allocatable :: dk(:,:,:),fk(:,:,:)
 real    :: k2
 integer :: i,j,k,idim

 allocate(dk(n,n,n),fk(n,n,n))
 dk = cmplx(deltastar,0.,kind=kind(dk))
 call dft_3d(dk,-1)
 ! potential -del^-2 delta, whose gradient is psi
 !$omp parallel do default(none) shared(n,dx,dk) private(i,j,k,k2)
 do k=1,n
    do j=1,n
       do i=1,n
          k2 = wavenumber(i,n,dx)**2 + wavenumber(j,n,dx)**2 + wavenumber(k,n,dx)**2
          if (k2 > 0.) then
             dk(i,j,k) = dk(i,j,k)/k2
          else
             dk(i,j,k) = 0.
          endif
       enddo
    enddo
 enddo
 !$omp end parallel do
 do idim=1,3
    call gradient_k(n,dx,dk,idim,fk)
    call inverse_transform(fk,psi(idim,:,:,:))
 enddo

end subroutine zeldovich_displacement

!-----------------------------------------------------------------------
!+
!  spectral derivative i k_idim f(k) along dimension idim. The Nyquist
!  mode is dropped so that the derivative of a real field stays real
!+
!-----------------------------------------------------------------------
subroutine gradient_k(n,dx,fk,idim,dfk)
 use utils_dft, only:wavenumber
 integer, intent(in)  :: n,idim
 real,    intent(in)  :: dx
 complex, intent(in)  :: fk(n,n,n)
 complex, intent(out) :: dfk(n,n,n)
 integer :: i,j,k,ijk(3)
 real    :: kd

 !$omp parallel do default(none) shared(n,dx,fk,dfk,idim) private(i,j,k,ijk,kd)
 do k=1,n
    do j=1,n
       do i=1,n
          ijk = (/i,j,k/)
          kd = wavenumber(ijk(idim),n,dx)
          if (2*(ijk(idim)-1) == n) kd = 0.
          dfk(i,j,k) = cmplx(0.,kd,kind=kind(dfk))*fk(i,j,k)
       enddo
    enddo
 enddo
 !$omp end parallel do

end subroutine gradient_k

!-----------------------------------------------------------------------
!+
!  normalised inverse transform of fk to a real field f
!+
!-----------------------------------------------------------------------
subroutine inverse_transform(fk,f)
 use utils_dft, only:dft_3d
 complex, intent(in)  :: fk(:,:,:)
 real,    intent(out) :: f(:,:,:)
 complex, allocatable :: work(:,:,:)

 work = fk
 call dft_3d(work,1)
 f = real(work)/real(size(work))

end subroutine inverse_transform

end module set_grf
