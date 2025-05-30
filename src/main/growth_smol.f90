!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module growth_smol
!
! Interface to library for dust growth and fragmentation
!   using Smoluchowsky solver. This module can be compiled
!   with a special rule in the Makefile to link against the
!   Smoluchowsky library
!
! :References:
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: physcon, smol2other, units
!
 implicit none

 public :: grain_growth_smol
 private

contains
!--------------------------------------------------------------------------
!+
!  Interface to Smoluchowsky library routines. This routine
!  receives everything in code units and phantom-style variables
!  and passes to what the library requires
!+
!--------------------------------------------------------------------------
subroutine grain_growth_smol(ndusttypes,dustfraci,rhoi,grainsize,graindens,dt)
 use smol2other, only:grain_growth   ! smoluchowsky library routine
 use units,      only:utime,umass,unit_density
 use physcon,    only:pi
 integer, intent(in) :: ndusttypes
 real,    intent(in)    :: grainsize(ndusttypes),graindens(ndusttypes)
 real,    intent(inout) :: dustfraci(ndusttypes)
 real,    intent(in)    :: rhoi,dt
 real :: rhodust(ndusttypes),grainmass(ndusttypes)
 !
 ! convert from dust fraction to dust density in g/cm^3
 !
 rhodust = rhoi*dustfraci*unit_density
 !
 ! convert grain sizes to grain mass in g
 !
 grainmass = 4./3.*pi*grainsize**3*graindens*umass
 !
 ! call Smoluchowsky solver
 !
 print*,' calling grain growth ntypes = ',ndusttypes
 print*,' grainmass = ',grainmass
 print*,' rhodust = ',rhodust

 call grain_growth(ndusttypes,rhodust,grainmass,dt,utime)

 print*,' AFTER GROWTH: rhodust = ',rhodust
 !
 ! convert from dust density to dust fraction
 !
 dustfraci = rhodust/rhoi

end subroutine grain_growth_smol

end module growth_smol
