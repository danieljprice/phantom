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
! :Dependencies: infile_utils, nicil, physcon
!
 use nicil, only: use_ohm,use_hall,use_ambi,na, &
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
subroutine read_options_nicil(name,valstring,imatch,igotall,ierr)
 use physcon, only:fourpi
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer                       :: ngotmax
 integer, save :: ngot = 0

 !--Initialise parameters
 imatch  = .true.
 igotall = .false.
 !--Number of input parameters
 ngotmax = 6

 !--Read input parameters
 select case(trim(name))
 case('use_ohm')
    read(valstring,*,iostat=ierr) use_ohm
    ngot = ngot + 1
 case('use_hall')
    read(valstring,*,iostat=ierr) use_hall
    ngot = ngot + 1
 case('use_ambi')
    read(valstring,*,iostat=ierr) use_ambi
    ngot = ngot + 1
 case('eta_constant')
    read(valstring,*,iostat=ierr) eta_constant
    ngot = ngot + 1
    if (eta_constant) then
       ngotmax = ngotmax + 1
    else
       ngotmax = ngotmax + 4
       if (na==1)  ngotmax = ngotmax + 1
    endif
 case('eta_const_type')
    read(valstring,*,iostat=ierr) eta_const_type
    ngot = ngot + 1
    if (eta_const_type==1) then
       if (use_ohm ) ngotmax = ngotmax + 1
       if (use_hall) ngotmax = ngotmax + 2
       if (use_ambi) ngotmax = ngotmax + 4
    elseif (eta_const_type==2 .or. eta_const_type==3) then
       if (use_ohm ) ngotmax = ngotmax + 1
       if (use_hall) ngotmax = ngotmax + 1
       if (use_ambi) ngotmax = ngotmax + 1
    endif
 case('C_OR')
    read(valstring,*,iostat=ierr) C_OR
    ngot = ngot + 1
 case('C_HE')
    read(valstring,*,iostat=ierr) C_HE
    ngot = ngot + 1
 case('C_AD')
    read(valstring,*,iostat=ierr) C_AD
    ngot = ngot + 1
 case('n_e_cnst')
    read(valstring,*,iostat=ierr) n_e_cnst
    ngot = ngot + 1
 case('hall_lt_zero')
    read(valstring,*,iostat=ierr) hall_lt_zero
    ngot = ngot + 1
 case('gamma_AD')
    read(valstring,*,iostat=ierr) gamma_AD
    ngot = ngot + 1
 case('rho_i_cnst')
    read(valstring,*,iostat=ierr) rho_i_cnst
    ngot = ngot + 1
 case('rho_n_cnst')
    read(valstring,*,iostat=ierr) rho_n_cnst
    ngot = ngot + 1
 case('alpha_AD')
    read(valstring,*,iostat=ierr) alpha_AD
    ngot = ngot + 1
 case('fdg')
    read(valstring,*,iostat=ierr) fdg
    ngot = ngot + 1
 case('rho_bulk')
    read(valstring,*,iostat=ierr) rho_bulk
    ngot = ngot + 1
 case('a0_grain')
    read(valstring,*,iostat=ierr) a0_grain
    ngot = ngot + 1
 case('an_grain')
    read(valstring,*,iostat=ierr) an_grain
    ngot = ngot + 1
 case('ax_grain')
    read(valstring,*,iostat=ierr) ax_grain
    ngot = ngot + 1
 case('zeta')
    read(valstring,*,iostat=ierr) zeta_cgs
    ngot = ngot + 1
 case('Cdt_diff')
    read(valstring,*,iostat=ierr) Cdt_diff
    ngot = ngot + 1
 case('Cdt_hall')
    read(valstring,*,iostat=ierr) Cdt_hall
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if ( ngot >= ngotmax ) igotall = .true.

end subroutine read_options_nicil

!-----------------------------------------------------------------------
end module nicil_sup
