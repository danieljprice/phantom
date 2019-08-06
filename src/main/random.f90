!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: random
!
!  DESCRIPTION:
!  this module contains a motley collection of random number
!  generator routines, used in various particle setups
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
module random
 implicit none
 public :: ran2,get_random,rayleigh_deviate

 private

contains

!-----------------------------------------------
!+
!  Simplified interface to get_random routine
!  always takes the same seed for seed 2
!  and stores the second seed internally for
!  subsequent calls
!+
!-----------------------------------------------
real function ran2(s1)
 integer, intent(inout) :: s1
 integer, save :: s2 = 123456789

 ran2 = get_random(s1,s2)

end function ran2

function get_random ( s1, s2 )

!*****************************************************************************
!
!! R8_UNI returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function generates uniformly distributed pseudorandom numbers
!    between 0 and 1, using the 32-bit generator from figure 3 of
!    the article by L'Ecuyer.
!
!    The cycle length is claimed to be 2.30584E+18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original Pascal original version by Pierre L'Ecuyer
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Pierre LEcuyer,
!    Efficient and Portable Combined Random Number Generators,
!    Communications of the ACM,
!    Volume 31, Number 6, June 1988, pages 742-751.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
!    seed for the sequence.  On first call, the user should initialize
!    S1 to a value between 1 and 2147483562;  S2 should be initialized
!    to a value between 1 and 2147483398.
!
!    Output, real ( kind = 8 ) R8_UNI, the next value in the sequence.
!
 integer  :: k
 real     :: get_random
 integer  :: s1
 integer  :: s2
 integer  :: z

 k = s1 / 53668
 s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
 if ( s1 < 0 ) then
    s1 = s1 + 2147483563
 endif

 k = s2 / 52774
 s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
 if ( s2 < 0 ) then
    s2 = s2 + 2147483399
 endif

 z = s1 - s2
 if ( z < 1 ) then
    z = z + 2147483562
 endif

 get_random = real ( z ) / 2147483563.0D+00

 return
end function get_random

!!-------------------------------------------------------------------------
!!
!! Function returns a random number drawn from a Rayleigh distribution
!! P(r) = r*e^(-r^2/(2*s^2))/s^2
!!
!! Useful for drawing amplitudes from a Gaussian distribution,
!! since the modulus is distributed according to a Rayleigh distribution.
!!
!!-------------------------------------------------------------------------
real function rayleigh_deviate(iseed)
 integer :: iseed

 rayleigh_deviate = sqrt(-log(ran2(iseed)))

end function rayleigh_deviate

end module random
