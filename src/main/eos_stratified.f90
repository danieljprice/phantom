!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module eos_stratified
!
! Implements vertically stratified equation of state
!
! :References: Dartois et al. (2003), Law et al. (2021)
!
! :Owner: Caitlyn Hardiman
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 public :: get_eos_stratified

 private

contains

!-----------------------------------------------------------------------
!+
!  Main eos routine
!+
!-----------------------------------------------------------------------
subroutine get_eos_stratified(istrat,xi,yi,zi,polyk,polyk2,qfacdisc,qfacdisc2,alpha_z,beta_z,z0,ponrhoi,spsoundi)
 integer, intent(in)  :: istrat
 real, intent(in)  :: xi,yi,zi,polyk,polyk2
 real, intent(in)  :: qfacdisc,qfacdisc2,alpha_z,beta_z,z0
 real, intent(out) :: ponrhoi,spsoundi
 real              :: r2,cs2mid,cs2atm,cs2,zq

 real, parameter :: pi = 4.*atan(1.0)

 r2 = xi**2 + yi**2
 cs2mid = polyk * r2**(-qfacdisc)
 cs2atm = polyk2 * r2**(-qfacdisc2)
 zq = z0 * r2**(0.5*beta_z)

 ! if istrat == 0 we use the MAPS prescription
 if (istrat==0) then
    ! modified equation 6 from Law et al. (2021)
    cs2 = (cs2mid**4 + 0.5*(1 + tanh((abs(zi) - alpha_z*zq)/zq))*cs2atm**4)**(1./4.)

    ! if istrat == 1 we use the Dartois et al. (2003) prescription
 elseif (istrat==1) then
    if (zi<zq) then
       cs2 = cs2atm + (cs2mid-cs2atm)*(cos((pi/2)*(zi/zq))**2)
    else
       cs2 = cs2atm
    endif
 endif

 ponrhoi = cs2
 spsoundi = sqrt(cs2)

end subroutine get_eos_stratified

end module eos_stratified
