!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2019 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
!+
!  MODULE: nicil_sup
!
!  DESCRIPTION:
!  Contains wrapper routines so that NICIL can be used in Phantom
!
!  REFERENCES: Wurster (2016)
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    C_AD           -- constant coefficient for ambipolar diffusion
!    C_HE           -- constant coefficient for the Hall effect (incl. sign)
!    C_OR           -- constant coefficient for ohmic resistivity
!    C_nimhd        -- coefficient to control the non-ideal MHD timestep
!    a0_grain       -- grain radius (cm)
!    alpha_AD       -- power law exponent for ambipolar diffusion
!    an_grain       -- minimum grain radius (cm)
!    ax_grain       -- maximum grain radius (cm)
!    eta_const_type -- Coefficient type: phys.cnst+B+rho (1), C_NI+B+rho (2), C_NI (3)
!    eta_constant   -- Use constant coefficients for all non-ideal MHD terms
!    fdg            -- dust-to-gas mass ratio
!    g_cnst         -- Use constant (T) / MRN (F) grain distribution
!    gamma_AD       -- ion-neutral coupling coefficient for ambipolar diffusion
!    hall_lt_zero   -- sign of the hall coefficient (<0 if T)
!    ion_rays       -- Include ionisation from cosmic
!    ion_thermal    -- Include thermal ionisation
!    mass_MionR_mp  -- metallic ion mass (m_proton)
!    massfrac_X     -- Hydrogen mass fraction
!    massfrac_Y     -- Helium mass fraction
!    n_e_cnst       -- constant electron number density
!    rho_bulk       -- bulk grain density (g/cm^3)
!    rho_i_cnst     -- ionisation density for ambipolar diffusion
!    rho_n_cnst     -- neutral density for ambipolar diffusion
!    use_ambi       -- Calculate the coefficient for ambipolar diffusion
!    use_hall       -- Calculate the coefficient for the Hall effect
!    use_massfrac   -- Use mass fraction (T) /abundance (F)
!    use_ohm        -- Calculate the coefficient for Ohmic resistivity
!    zeta           -- cosmic ray ionisation rate (s^-1)
!
!  DEPENDENCIES: eos, infile_utils, nicil, physcon
!+
!--------------------------------------------------------------------------
module nicil_sup
 use nicil, only: use_ohm,use_hall,use_ambi, &
                  ion_rays,ion_thermal,use_massfrac,massfrac_X,massfrac_Y, &
                  g_cnst,fdg,rho_bulk,a0_grain,an_grain,ax_grain, &
                  zeta_cgs,mass_MionR_mp,C_nimhd, &
                  eta_constant,eta_const_type,icnstphys,icnstsemi,icnst,C_OR,C_HE,C_AD, &
                  n_e_cnst,hall_lt_zero,rho_i_cnst,rho_n_cnst,alpha_AD,gamma_AD
 implicit none
 !
 !--Subroutines
 public  ::use_consistent_gmw,write_options_nicil,read_options_nicil
 !
 private

contains
!
!-----------------------------------------------------------------------
!+
!  Ensures a consistent meanmolecular mass is used
!+
!-----------------------------------------------------------------------
subroutine use_consistent_gmw(ierr,gmw_old,gmw_new)
 !
 use nicil,   only:meanmolmass
 use eos,     only:gmw
 integer, intent(out) :: ierr
 real,    intent(out) :: gmw_old,gmw_new
 !
 if (abs(meanmolmass-gmw) > epsilon(gmw)) then
    ierr    = 1
    gmw_old = gmw
    gmw_new = meanmolmass
    gmw     = meanmolmass
 endif
 !
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
    call write_inopt(eta_const_type,'eta_const_type','Coefficient type: phys.cnst+B+rho (1), C_NI+B+rho (2), C_NI (3)',iunit)
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
 call write_inopt(C_nimhd,           'C_nimhd',      'coefficient to control the non-ideal MHD timestep',iunit)
 if ( .not. eta_constant ) then
    write(iunit,"(/,a)") '# options controlling ionisation'
    call write_inopt(ion_rays,       'ion_rays',     'Include ionisation from cosmic',iunit)
    call write_inopt(ion_thermal,    'ion_thermal',  'Include thermal ionisation',iunit)
    call write_inopt(use_massfrac,   'use_massfrac', 'Use mass fraction (T) /abundance (F)',iunit)
    if (use_massfrac) then
       call write_inopt(massfrac_X,   'massfrac_X',  'Hydrogen mass fraction',iunit)
       call write_inopt(massfrac_Y,   'massfrac_Y',  'Helium mass fraction',iunit)
    endif
    if (ion_rays) then
       call write_inopt(g_cnst,       'g_cnst',      'Use constant (T) / MRN (F) grain distribution',iunit)
       call write_inopt(fdg,          'fdg',         'dust-to-gas mass ratio',iunit)
       call write_inopt(rho_bulk,     'rho_bulk',    'bulk grain density (g/cm^3)',iunit)
       if ( g_cnst ) then
          call write_inopt(a0_grain,   'a0_grain',   'grain radius (cm)',iunit)
       else
          call write_inopt(an_grain,   'an_grain',   'minimum grain radius (cm)',iunit)
          call write_inopt(ax_grain,   'ax_grain',   'maximum grain radius (cm)',iunit)
       endif
       call write_inopt(zeta_cgs,     'zeta',         'cosmic ray ionisation rate (s^-1)',iunit)
       call write_inopt(mass_MionR_mp,'mass_MionR_mp','metallic ion mass (m_proton)',iunit)
    endif
 endif
 !
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
 ngotmax = 5
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
       ngotmax = ngotmax + 1 ! for eta_const_calc
    else
       ngotmax = ngotmax + 3
    endif
 case('eta_const_type')
    read(valstring,*,iostat=ierr) eta_const_type
    ngot = ngot + 1
    if (eta_const_type==1) then
       if (use_ohm .or. use_hall) ngotmax = ngotmax + 1
       if (use_hall)              ngotmax = ngotmax + 1
       if (use_ambi)              ngotmax = ngotmax + 4
    else
       if (use_ohm )              ngotmax = ngotmax + 1
       if (use_hall)              ngotmax = ngotmax + 1
       if (use_ambi)              ngotmax = ngotmax + 1
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
 case('ion_rays')
    read(valstring,*,iostat=ierr) ion_rays
    ngot = ngot + 1
    if (ion_rays) ngotmax = ngotmax + 5
 case('ion_thermal')
    read(valstring,*,iostat=ierr) ion_thermal
    ngot = ngot + 1
 case('use_massfrac')
    read(valstring,*,iostat=ierr) use_massfrac
    ngot = ngot + 1
    if (use_massfrac) ngotmax = ngotmax + 2
 case('massfrac_X')
    read(valstring,*,iostat=ierr) massfrac_X
    ngot = ngot + 1
 case('massfrac_Y')
    read(valstring,*,iostat=ierr) massfrac_Y
    ngot = ngot + 1
 case('g_cnst')
    read(valstring,*,iostat=ierr) g_cnst
    ngot = ngot + 1
    if (g_cnst) then
       ngotmax = ngotmax + 1
    else
       ngotmax = ngotmax + 2
    endif
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
 case('mass_MionR_mp')
    read(valstring,*,iostat=ierr) mass_MionR_mp
    ngot = ngot + 1
 case('C_nimhd')
    read(valstring,*,iostat=ierr) C_nimhd
    ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if ( ngot >= ngotmax ) igotall = .true.
 !
end subroutine read_options_nicil

!-----------------------------------------------------------------------

end module nicil_sup
