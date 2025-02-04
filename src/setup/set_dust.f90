!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module set_dust
!
! Contains utility routines for setting up dust size distributions
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: table_utils
!

 implicit none

 public :: set_dustfrac
 public :: set_dustbinfrac
 private

contains

!--------------------------------------------------------------------------
!+
!  utility function to set the dust fraction given the
!  dust-to-gas ratio. Equation (57) in Price & Laibe (2015)
!+
!--------------------------------------------------------------------------
subroutine set_dustfrac(dust_to_gas,dustfrac)
 real, intent(in)  :: dust_to_gas
 real, intent(out) :: dustfrac(:)

 dustfrac = dust_to_gas/(1.+dust_to_gas)

end subroutine set_dustfrac

!--------------------------------------------------------------------------
!+
!  utility function to set the percentage of dust mass in each dust bin
!  given the grain sizes and the slope of the number density distribution
!+
!--------------------------------------------------------------------------
subroutine set_dustbinfrac(smin,smax,sindex,dustbinfrac,grainsize)
 use table_utils, only:logspace
 real, intent(in)  :: smin
 real, intent(in)  :: smax
 real, intent(in)  :: sindex
 real, intent(out) :: dustbinfrac(:)
 real, intent(out) :: grainsize(:)
 integer :: i,nbins
 real :: rhodust(size(dustbinfrac))
 real :: grid(size(dustbinfrac)+1)
 real :: exact
 real :: power
 real, parameter :: tol = 1.e-10

 rhodust = 0.
 exact = 0.
 power = 0.

 nbins = size(dustbinfrac)
 call logspace(grid,smin,smax)

 !--Dust density is computed from drhodust ∝ dn*mdust where dn ∝ s**(-p)*ds
 !  and mdust ∝ s**(3). This is then integrated across each cell to account
 !  for mass contributions from unrepresented grain sizes
 do i=1,nbins
    !--Find representative s for each cell (geometric mean)
    grainsize(i)=sqrt(grid(i)*grid(i+1))
    if (sindex == 4.) then
       rhodust(i) = log(grid(i+1)/grid(i))
    else
       power = 4. - sindex
       rhodust(i) = 1./power*(grid(i+1)**power - grid(i)**power)
    endif
 enddo

 !--Calculate the percentage of dust in each size bin (NOT equal to dustfrac)
 dustbinfrac(:) = rhodust(:)/sum(rhodust)

 !--Check to make sure the integral determining the contributions is correct
 if (sindex == 4.) then
    exact = log(grid(nbins+1)/grid(1))
 else
    exact = 1./power*(grid(nbins+1)**power - grid(1)**power)
 endif
 if (abs(sum(rhodust)-exact)/exact > tol) then
    print*,' WARNING: Piecewise integration of MRN distribution not matching the exact solution!'
 endif

end subroutine set_dustbinfrac

end module set_dust
