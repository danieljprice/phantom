!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: powerspectrums
!
!  DESCRIPTION: None
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------

!
! contains subroutines for taking power spectrums on particle data
!
module powerspectrums
 implicit none
 real, parameter, private :: pi = 3.141592653589
 real, parameter, private :: twopi = 2.*pi
 public :: powerspectrum

 private

contains

subroutine powerspectrum(npts,x,dat,nfreqpts,freq,power,idisordered)
 implicit none
 integer, intent(in) :: npts, nfreqpts
 real, intent(in), dimension(npts) :: x
 real, intent(in), dimension(npts) :: dat
 real, intent(in), dimension(nfreqpts) :: freq
 real, intent(out), dimension(nfreqpts) :: power
 logical, intent(in) :: idisordered
 integer :: ifreq
 real :: datmean, datvar, omega

 if (.not.idisordered) then
    print*,' evaluating fourier transform'
    do ifreq=1,nfreqpts
       omega = twopi*freq(ifreq)
       !--get power at this frequency
       call power_fourier(npts,x,dat,omega,power(ifreq))
    enddo
 else
    print*,'evaluating lomb periodogram...'
!
!--calculate the mean and variance of the data
!
    call mean_variance(dat,npts,datmean,datvar)
    print*,'data mean = ',datmean,' std. dev = ',sqrt(datvar)
    if (datvar <= 0.) then
       print*,'error: variance = 0'
       power = 0.
       return
    endif
    do ifreq=1,nfreqpts
       omega = twopi*freq(ifreq)
       call power_lomb(npts,x,dat,datmean,datvar,omega,power(ifreq))
    enddo

 endif

end subroutine powerspectrum

!-------------------------------------------------------
! subroutine to compute the power spectrum
! of evenly sampled data via a (slow) fourier transform
!--------------------------------------------------------

subroutine power_fourier(npts,x,dat,omega,power)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x, dat
 real, intent(in) :: omega
 real, intent(out) :: power
 integer :: i
 real :: sum1,sum2

 power = 0.
 sum1 = 0.
 sum2 = 0.
 do i=1,npts
    sum1 = sum1 + dat(i)*cos(-omega*x(i))
    sum2 = sum2 + dat(i)*sin(-omega*x(i))
 enddo
 power= sqrt(sum1**2 + sum2**2)/REAL(npts)

 return
end subroutine power_fourier

!----------------------------------------------------------
! Subroutine to compute the power spectrum (periodogram)
! of unevenly sampled data via the Lomb (1976) method
! (algorithm described in Press et al, Numerical Recipes, sec 13.8, p569)
!
! Given the data (dat) on a set of points (x),
! returns an array of nfreq frequencies (freq) between freqmin and freqmax
! together with the power at each frequency (power)
!----------------------------------------------------------
subroutine power_lomb(npts,x,dat,datmean,datvar,omega,power)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x, dat
 real, intent(in) :: datmean,datvar,omega
 real, intent(out) :: power
 integer :: i
 real :: ddat
 real :: tau, tau_numerator, tau_denominator
 real :: term1_numerator, term1_denominator
 real :: term2_numerator, term2_denominator
 real :: omega_dx, cos_term, sin_term
!
!--calculate tau for this frequency
!
 tau_numerator = 0.
 tau_denominator = 0.
 do i=1,npts
    tau_numerator = tau_numerator + sin(2.*omega*x(i))
    tau_denominator = tau_denominator + cos(2.*omega*x(i))
 enddo
 tau = atan(tau_numerator/tau_denominator)/(2.*omega)
!
!--calculate the terms in the power
!
 term1_numerator = 0.
 term1_denominator = 0.
 term2_numerator = 0.
 term2_denominator = 0.
 do i=1,npts
    ddat = dat(i) - datmean
    omega_dx = omega*(x(i) - tau)
    cos_term = cos(omega_dx)
    sin_term = sin(omega_dx)
    term1_numerator = term1_numerator + ddat*cos_term
    term1_denominator = term1_denominator + cos_term**2
    term2_numerator = term2_numerator + ddat*sin_term
    term2_denominator = term2_denominator + sin_term**2
 enddo
!
!--calculate the power at this frequency
!
 power = 1./(datvar)*(term1_numerator**2/term1_denominator + &
                         term2_numerator**2/term2_denominator)

 return
end subroutine power_lomb

!-------------------------------------------------
! Subroutine to calculate the mean and variance
! of a set of data points
! Mean is trivial but variance uses a special
! formula to reduce round-off error
! see Press et al Numerical Recipes, section 14.2
! this is similar to their subroutine avevar
!-------------------------------------------------
subroutine mean_variance(x,npts,xmean,xvariance)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(npts) :: x
 real, intent(out) :: xmean, xvariance
 real :: roundoff, delta
 integer :: i
!
!--calculate average
!
 xmean = 0.
 do i=1,npts
    xmean = xmean + x(i)
 enddo
 xmean = xmean/real(npts)
!
!--calculate variance using the corrected two-pass formula
!
!    var = 1/(n-1)*( sum (x-\bar{x}) - 1/n * (sum(x-\bar{x}) )^2 )
!
!  where the last term corrects for the roundoff error
!  in the first term
!
 xvariance = 0.
 roundoff = 0.

 do i=1,npts
    delta = x(i) - xmean
    roundoff = roundoff + delta
    xvariance = xvariance + delta*delta
 enddo
 xvariance = (xvariance - roundoff**2/npts)/real(npts-1)

 return
end subroutine mean_variance

end module powerspectrums
