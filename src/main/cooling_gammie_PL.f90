!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_gammie_PL
!
! Simple beta-cooling prescription used for experiments on gravitational
!  instability in discs
!
! :References:
!   Gammie (2001), ApJ 553, 174-183
!
! :Owner: Cristiano Longarini
!
! :Runtime parameters:
!   - beta_cool : *beta factor in Gammie (2001) cooling @ R_beta*
!   - eta       : *Power law coefficient of the cooling factor*
!   - r_beta    : *Characteristic radius of the cooling power law profile*
!
! :Dependencies: infile_utils
!
 implicit none
 real, private :: beta_cool = 3., eta = 1.5, r_beta = 50.

contains
!-----------------------------------------------------------------------
!+
!   Gammie PL cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_Gammie_PL_explicit(xi,yi,zi,ui,dudti)

 real, intent(in)    :: ui,xi,yi,zi
 real, intent(inout) :: dudti

 real :: omegai,r2,tcool1

 r2     = xi*xi + yi*yi + zi*zi
 Omegai = r2**(-0.75)
 tcool1 = Omegai/(beta_cool * (r2/r_beta**2)** (-eta/2))
 dudti  = dudti - ui*tcool1

end subroutine cooling_Gammie_PL_explicit

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling_gammie_PL(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling @ R_beta',iunit)
 call write_inopt(eta,'eta','Power law coefficient of the cooling factor',iunit)
 call write_inopt(r_beta,'r_beta','Characteristic radius of the cooling power law profile',iunit)

end subroutine write_options_cooling_gammie_PL

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling_gammie_PL(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(beta_cool,'beta_cool',db,errcount=nerr,min=1.)
 call read_inopt(eta,'eta',db,errcount=nerr,default=eta)
 call read_inopt(r_beta,'r_beta',db,errcount=nerr,min=0.,default=r_beta)

end subroutine read_options_cooling_gammie_PL

end module cooling_gammie_PL
