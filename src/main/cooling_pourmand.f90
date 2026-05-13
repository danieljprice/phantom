!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_pourmand
!
! 
!
! :Runtime parameters:
!   - rho_norm : *normalization density for cooling timescale* (default: 1e8, in code units)
!
! :Dependencies: infile_utils, part
!
 implicit none
 real, private :: beta_cool  = 3.
 real, private :: rho_norm = 1e7

contains
!----------------------------------------------------------------
!+
!  Ali: this is our cooling subroutine which since we wanted an operator split, we will call after the step call in the main evolve loop.
! It simply cools the particles according to an exponential cooling law, with a cooling timescale tau_cool. 
!You can modify this to implement your own cooling law, and also add a mask to only cool certain particles (e.g. based on temperature or density).
!+
!----------------------------------------------------------------
subroutine cool_pourmand_analytical_exp(dt)
  use part, only:  npart, xyzh, vxyzu, eos_vars, itemp, rhoh, massoftype, iphase, xyzmh_ptmass, nptmass
  use eos,  only: get_u_from_rhoT, eos_type
  integer, intent(in) :: npart
  real, intent(in) :: dt
  integer :: i
  real :: rhoi, Tnew, Told, tau_cool, rho_norm ! I should add this as an input later if it's useful enough
  real :: omegai,r2

  !$omp parallel do default(none) &
  !$omp shared(npart,xyzh,vxyzu,eos_vars,itemp,dt,tau_cool,massoftype,iphase) &
  !$omp private(i,rhoi,Told,Tnew)

  do i = 1, npart
      if (nptmass > 0) then
         r2     = (xi-xyzmh_ptmass(1,1))**2 + (yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2
      else
         r2     = xi*xi + yi*yi + zi*zi
      endif
      rho_units_to_cgs = 5.9 ! for 1M_sun/R_sun^3=5.9g/cm^3, adjust as needed for your code units
      rho_norm = 1e7/rho_units_to_cgs ! assuming this is in code units, adjust as needed
      Omegai = r2**(-0.75)
      Told = eos_vars(itemp,i)

      rhoi = rhoh(xyzh(4,i), massoftype(iphase(i)))   ! adjust if needed
      tau_cool = Omegai/(1+rho_norm/rhoi)

      Tnew = Told * exp(-dt/tau_cool)

      vxyzu(4,i) = get_u_from_rhoT(rhoi, Tnew, eos_type, vxyzu(4,i))


  end do
  !$omp end parallel do
  
end subroutine cool_pourmand_analytical_exp


!-----------------------------------------------------------------------
!+
!   Gammie (2001) cooling
!+
!-----------------------------------------------------------------------
! subroutine cooling_Gammie_explicit(xi,yi,zi,ui,dudti)
!  use part, only:xyzmh_ptmass, nptmass
!  real, intent(in)    :: ui,xi,yi,zi
!  real, intent(inout) :: dudti

!  real :: omegai,r2,tcool1

!  if (nptmass > 0) then
!     r2     = (xi-xyzmh_ptmass(1,1))**2 + (yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2
!  else
!     r2     = xi*xi + yi*yi + zi*zi
!  endif

!  Omegai = r2**(-0.75)
!  tcool1 = Omegai/beta_cool
!  dudti  = dudti - ui*tcool1

! end subroutine cooling_Gammie_explicit

! !-----------------------------------------------------------------------
! !+
! !  writes input options to the input file
! !+
! !-----------------------------------------------------------------------
! subroutine write_options_cooling_gammie(iunit)
!  use infile_utils, only:write_inopt
!  integer, intent(in) :: iunit

!  call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)

! end subroutine write_options_cooling_gammie

! !-----------------------------------------------------------------------
! !+
! !  reads input options from the input file
! !+
! !-----------------------------------------------------------------------
! subroutine read_options_cooling_gammie(db,nerr)
!  use infile_utils, only:inopts,read_inopt
!  type(inopts), intent(inout) :: db(:)
!  integer,      intent(inout) :: nerr

!  call read_inopt(beta_cool,'beta_cool',db,errcount=nerr,min=1.)

! end subroutine read_options_cooling_gammie

end module cooling_pourmand
