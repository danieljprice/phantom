!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Dependencies: infile_utils, io
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
subroutine read_options_cooling_gammie_PL(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0

 imatch  = .true.
 igotall = .true. ! none of the cooling options are compulsory
 select case(trim(name))
 case('beta_cool')
    read(valstring,*,iostat=ierr) beta_cool
    ngot = ngot + 1
 case('eta')
    read(valstring,*,iostat=ierr) eta
    ngot = ngot + 1
 case('r_beta')
    read(valstring,*,iostat=ierr) r_beta
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (ngot >= 1) igotall = .true.

end subroutine read_options_cooling_gammie_PL

end module cooling_gammie_PL
