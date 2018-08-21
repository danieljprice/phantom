!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: set_dust
!
!  DESCRIPTION:
!  Contains utility routines for setting up dust size distributions
!
!  REFERENCES:
!
!  OWNER: Daniel Mentiplay
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!
!  DEPENDENCIES:
!+
!--------------------------------------------------------------------------

module set_dust
 implicit none

 interface set_dustfrac
  module procedure set_dustfrac_single,set_dustfrac_power_law
 end interface set_dustfrac
 public :: set_dustfrac
 public :: set_dustfrac_single
 public :: set_dustfrac_power_law
 private

contains

!--------------------------------------------------------------------------
!+
!  utility function to set the dust fraction given the
!  dust-to-gas ratio. Equation (57) in Price & Laibe (2015)
!+
!--------------------------------------------------------------------------
subroutine set_dustfrac_single(dust_to_gas,dustfrac)
 real, intent(in)  :: dust_to_gas
 real, intent(out) :: dustfrac

 dustfrac = dust_to_gas/(1.+dust_to_gas)

end subroutine set_dustfrac_single

!--------------------------------------------------------------------------
!+
!  utility function to set the dust fraction given the dust-to-gas ratio,
!  grain sizes, and slope of dust number density distribution in size
!+
!--------------------------------------------------------------------------
subroutine set_dustfrac_power_law(dust_to_gas_tot,grainsize,sindex,dustfrac)
 use io, only:warning
 real, intent(in)  :: dust_to_gas_tot,grainsize(:),sindex
 real, intent(out) :: dustfrac(:)
 integer :: i,nbins
 real :: dustfrac_tot
 real :: norm
 real :: rhodtot
 real :: rhodust(size(dustfrac))
 real :: exact
 real :: power
 real, parameter :: tol = 1.e-10

 dustfrac_tot = 0.
 norm = 0.
 rhodtot = 0.
 rhodust = 0.
 exact = 0.
 power = 0.

 if (size(dustfrac) /= size(grainsize)) then
    call warning('set_dust','grainsize and dustfrac arrays should have same size')
 endif

 nbins = size(dustfrac)

 !--Dust density is computed from drhodust ∝ dn*mdust where dn ∝ s**(-p)*ds
 !  and mdust ∝ s**(3). This is then integrated across each cell to account
 !  for mass contributions from unrepresented grain sizes
 do i=1,nbins
    if (sindex == 4.) then
       rhodust(i) = log(grainsize(i+1)/grainsize(i))
    else
       power = 4. - sindex
       rhodust(i) = 1./power*(grainsize(i+1)**power - grainsize(i)**power)
    endif
 enddo

 !--Sum the contributions from each cell to get total relative dust content
 rhodtot = sum(rhodust)

 !--Calculate the total dust fraction from the dust-to-gas ratio
 dustfrac_tot = dust_to_gas_tot/(1.+dust_to_gas_tot)

 !--Calculate the normalisation factor (∝ 1/rhotot) and scale the dust fractions
 !  Note: dust density and dust fraction have the same power-law dependence on s.
 norm        = dustfrac_tot/rhodtot
 dustfrac(:) = norm*rhodust(:)

 !--Check to make sure the integral determining the contributions is correct
 if (sindex == 4.) then
    exact = log(grainsize(nbins+1)/grainsize(1))
 else
    exact = 1./power*(grainsize(nbins+1)**power - grainsize(1)**power)
 endif
 if (abs(rhodtot-exact)/exact > tol) then
    call warning('set_dust','Piecewise integration of MRN distribution not matching the exact solution!')
 endif

end subroutine set_dustfrac_power_law

end module set_dust
