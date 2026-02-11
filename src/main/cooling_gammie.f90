!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_gammie
!
! Simple beta-cooling prescription used for experiments on gravitational
!  instability in discs
!
! :References:
!   Gammie (2001), ApJ 553, 174-183
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - beta_cool : *beta factor in Gammie (2001) cooling*
!
! :Dependencies: externalforces, infile_utils, part
!
 implicit none
 real, private :: beta_cool  = 3.

contains
!-----------------------------------------------------------------------
!+
!   Gammie (2001) cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_Gammie_explicit(xi,yi,zi,ui,dudti)
 use part,           only:xyzmh_ptmass, nptmass
 use externalforces, only:mass1
 real, intent(in)    :: ui,xi,yi,zi
 real, intent(inout) :: dudti

 real :: omegai,r2,tcool1,m1

 if (nptmass > 0) then
    r2 = (xi-xyzmh_ptmass(1,1))**2 + (yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2
    m1 = xyzmh_ptmass(4,1)
 else
    r2 = xi*xi + yi*yi + zi*zi
    m1 = mass1
 endif

 Omegai = sqrt(m1)*r2**(-0.75)
 tcool1 = Omegai/beta_cool
 dudti  = dudti - ui*tcool1

end subroutine cooling_Gammie_explicit

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_gammie(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)

end subroutine write_options_cooling_gammie

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_gammie(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(beta_cool,'beta_cool',db,errcount=nerr,min=1.)

end subroutine read_options_cooling_gammie

end module cooling_gammie
