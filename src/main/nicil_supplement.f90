!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module nicil_sup
!
! Contains wrapper routines so that NICIL can be used in Phantom
!
! :References: Wurster (2016)
!              Wurster (2021)
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - C_AD           : *constant coefficient for ambipolar diffusion*
!   - C_HE           : *constant coefficient for the Hall effect (incl. sign)*
!   - C_OR           : *constant coefficient for ohmic resistivity*
!   - Cdt_diff       : *coefficient to control the Ohmic & ambipolar timesteps*
!   - Cdt_hall       : *coefficient to control the Hall timestep*
!   - a0_grain       : *grain radius (cm)*
!   - alpha_AD       : *power law exponent for ambipolar diffusion*
!   - an_grain       : *minimum grain radius (cm)*
!   - ax_grain       : *maximum grain radius (cm)*
!   - eta_const_type : *Coefficient type: (1) phys.cnst+B+rho (2) C_NI+B+rho (3) constant*
!   - eta_constant   : *Use constant coefficients for all non-ideal MHD terms*
!   - fdg            : *dust-to-gas mass ratio*
!   - gamma_AD       : *ion-neutral coupling coefficient for ambipolar diffusion*
!   - hall_lt_zero   : *sign of the hall coefficient (<0 if T)*
!   - n_e_cnst       : *constant electron number density*
!   - rho_bulk       : *bulk grain density (g/cm^3)*
!   - rho_i_cnst     : *ionisation density for ambipolar diffusion*
!   - rho_n_cnst     : *neutral density for ambipolar diffusion*
!   - use_ambi       : *Calculate the coefficient for ambipolar diffusion*
!   - use_hall       : *Calculate the coefficient for the Hall effect*
!   - use_ohm        : *Calculate the coefficient for Ohmic resistivity*
!   - zeta           : *cosmic ray ionisation rate (s^-1)*
!
! :Dependencies: infile_utils, nicil
!
 use nicil, only:use_ohm,use_hall,use_ambi,na, &
                  fdg,rho_bulk,a0_grain,an_grain,ax_grain,zeta_cgs,Cdt_diff,Cdt_hall, &
                  eta_constant,eta_const_type,icnstphys,icnstsemi,icnst,C_OR,C_HE,C_AD, &
                  n_e_cnst,hall_lt_zero,rho_i_cnst,rho_n_cnst,alpha_AD,gamma_AD
 implicit none
 !
 !--Subroutines
 public  :: use_consistent_gmw,write_options_nicil,read_options_nicil

 private

contains

!-----------------------------------------------------------------------
!+
!  Ensures a consistent meanmolecular mass is used
!+
!-----------------------------------------------------------------------
subroutine use_consistent_gmw(ierr,gmw_eos,gmw_nicil)
 use nicil,   only:meanmolmass
 integer, intent(out)   :: ierr
 real,    intent(out)   :: gmw_nicil
 real,    intent(inout) :: gmw_eos

 gmw_nicil = meanmolmass
 if (abs(meanmolmass-gmw_eos) > epsilon(gmw_eos)) then
    ierr    = 1
    gmw_eos = meanmolmass
 endif

end subroutine use_consistent_gmw
!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_nicil(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling non-ideal MHD'
 call write_inopt(use_ohm ,    'use_ohm'     ,'Calculate the coefficient for Ohmic resistivity',iunit)
 call write_inopt(use_hall,    'use_hall'    ,'Calculate the coefficient for the Hall effect',iunit)
 call write_inopt(use_ambi,    'use_ambi'    ,'Calculate the coefficient for ambipolar diffusion',iunit)
 call write_inopt(eta_constant,'eta_constant','Use constant coefficients for all non-ideal MHD terms',iunit)
 if ( eta_constant ) then
    call write_inopt(eta_const_type,'eta_const_type','Coefficient type: (1) phys.cnst+B+rho (2) C_NI+B+rho (3) constant',iunit)
    if ( eta_const_type==1 ) then
       if ( use_ohm .or. use_hall  ) then
          call write_inopt(n_e_cnst,'n_e_cnst' ,'constant electron number density',iunit)
          if ( use_hall ) call write_inopt(hall_lt_zero, 'hall_lt_zero' ,'sign of the hall coefficient (<0 if T)',iunit)
       endif
       if ( use_ambi ) then
          call write_inopt(gamma_AD,    'gamma_AD',  'ion-neutral coupling coefficient for ambipolar diffusion',iunit)
          call write_inopt(rho_i_cnst,  'rho_i_cnst','ionisation density for ambipolar diffusion',iunit)
          call write_inopt(rho_n_cnst,  'rho_n_cnst','neutral density for ambipolar diffusion',iunit)
          call write_inopt(alpha_AD,    'alpha_AD',  'power law exponent for ambipolar diffusion',iunit)
       endif
    elseif ( eta_const_type==2 ) then
       if ( use_ohm  ) call write_inopt(C_OR,'C_OR', 'semi-constant coefficient for ohmic resistivity',iunit)
       if ( use_hall ) call write_inopt(C_HE,'C_HE', 'semi-constant coefficient for the Hall effect (incl. sign)',iunit)
       if ( use_ambi ) call write_inopt(C_AD,'C_AD', 'semi-constant coefficient for ambipolar diffusion',iunit)
    elseif ( eta_const_type==3 ) then
       if ( use_ohm  ) call write_inopt(C_OR,'C_OR', 'constant coefficient for ohmic resistivity',iunit)
       if ( use_hall ) call write_inopt(C_HE,'C_HE', 'constant coefficient for the Hall effect (incl. sign)',iunit)
       if ( use_ambi ) call write_inopt(C_AD,'C_AD', 'constant coefficient for ambipolar diffusion',iunit)
    endif
 endif
 call write_inopt(Cdt_diff,          'Cdt_diff',     'coefficient to control the Ohmic & ambipolar timesteps',iunit)
 call write_inopt(Cdt_hall,          'Cdt_hall',     'coefficient to control the Hall timestep',iunit)
 if ( .not. eta_constant ) then
    write(iunit,"(/,a)") '# options controlling ionisation'
    call write_inopt(fdg,          'fdg',         'dust-to-gas mass ratio',iunit)
    call write_inopt(rho_bulk,     'rho_bulk',    'bulk grain density (g/cm^3)',iunit)
    if ( na==1 ) then
       call write_inopt(a0_grain,   'a0_grain',   'grain radius (cm)',iunit)
    else
       call write_inopt(an_grain,   'an_grain',   'minimum grain radius (cm)',iunit)
       call write_inopt(ax_grain,   'ax_grain',   'maximum grain radius (cm)',iunit)
    endif
    call write_inopt(zeta_cgs,     'zeta',         'cosmic ray ionisation rate (s^-1)',iunit)
 endif

end subroutine write_options_nicil
!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_nicil(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(use_ohm,'use_ohm',db,errcount=nerr)
 call read_inopt(use_hall,'use_hall',db,errcount=nerr)
 call read_inopt(use_ambi,'use_ambi',db,errcount=nerr)
 call read_inopt(eta_constant,'eta_constant',db,errcount=nerr)
 if (eta_constant) then
    call read_inopt(eta_const_type,'eta_const_type',db,errcount=nerr)
    if (eta_const_type==1) then
       if (use_ohm .or. use_hall) call read_inopt(n_e_cnst,'n_e_cnst',db,errcount=nerr)
       if (use_hall) call read_inopt(hall_lt_zero,'hall_lt_zero',db,errcount=nerr)
       if (use_ambi) then
          call read_inopt(gamma_AD,'gamma_AD',db,errcount=nerr)
          call read_inopt(rho_i_cnst,'rho_i_cnst',db,errcount=nerr)
          call read_inopt(rho_n_cnst,'rho_n_cnst',db,errcount=nerr)
          call read_inopt(alpha_AD,'alpha_AD',db,errcount=nerr)
       endif
    elseif (eta_const_type==2 .or. eta_const_type==3) then
       if (use_ohm) call read_inopt(C_OR,'C_OR',db,errcount=nerr)
       if (use_hall) call read_inopt(C_HE,'C_HE',db,errcount=nerr)
       if (use_ambi) call read_inopt(C_AD,'C_AD',db,errcount=nerr)
    endif
 endif
 call read_inopt(Cdt_diff,'Cdt_diff',db,errcount=nerr)
 call read_inopt(Cdt_hall,'Cdt_hall',db,errcount=nerr)
 if (.not. eta_constant) then
    call read_inopt(fdg,'fdg',db,errcount=nerr)
    call read_inopt(rho_bulk,'rho_bulk',db,errcount=nerr)
    if (na==1) then
       call read_inopt(a0_grain,'a0_grain',db,errcount=nerr)
    else
       call read_inopt(an_grain,'an_grain',db,errcount=nerr)
       call read_inopt(ax_grain,'ax_grain',db,errcount=nerr)
    endif
    call read_inopt(zeta_cgs,'zeta',db,errcount=nerr)
 endif

end subroutine read_options_nicil

end module nicil_sup
