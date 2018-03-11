!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: h2cooling
!
!  DESCRIPTION:
!  Contains routines for cooling
!  Routines are originally by Simon Glover,
!  Translated to Fortran 90 and adapted
!  for use in Phantom by Daniel Price (2011)
!
!  REFERENCES:
!   Sembach et al. (2000) ApJ 528, 310
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS:
!    abundc            -- Carbon abundance
!    abunde            -- electron abundance
!    abundo            -- Oxygen abundance
!    abundsi           -- Silicon abundance
!    dchem             -- distance for chemistry of HI
!    dlq               -- distance for column density in cooling function
!    dphot             -- photodissociation distance used for CO/H2
!    dphotflag         -- photodissociation distance static or radially adaptive (0/1)
!    dust_to_gas_ratio -- dust to gas ratio
!    iflag_atom        -- Which atomic cooling (1:Gal ISM, 2:Z=0 gas)
!    iphoto            -- Photoelectric heating treatment (0=optically thin, 1=w/extinction)
!
!  DEPENDENCIES: fs_data, infile_utils, io, mol_data, part, physcon,
!    splineutils
!+
!--------------------------------------------------------------------------
module h2cooling
 use physcon, only:kboltz,eV
 implicit none
!
! only publicly visible entries are the
! cool_func and init_h2cooling subroutines
!
 public :: cool_func, init_h2cooling, write_options_h2cooling, read_options_h2cooling
 public :: hchem
!
! required constants and parameters
! (this stuff formerly in cool.h)
!
! He:H ratio by number (=> ratio by mass is 4*abhe)
 real, parameter :: abhe = 0.1d0

! Number of entries in cooling table
 integer, parameter :: nmd = 10000

! Number of cooling / heating rates computed in cooling fn.
 integer, parameter, public :: nrates = 12

! Size of abundance array that is passed to cool_func as input
 integer, parameter, public :: nabn = 10

! Number of different quantities stored in cooling look-up table
 integer, parameter :: ncltab = 54

! These varables are initialised in init_h2cooling
 real :: temptab(nmd)
 real :: cltab(ncltab, nmd),dtcltab(ncltab, nmd)
 real :: dtlog, tmax, tmin

! Parameters and tables used for CO rotational cooling
!
! Total number of combinations of column density, temperature for which we have tabulated rates
 integer, parameter, public :: nTco = 1996

! Number of CO column densities for which we have tabuled data
 integer, parameter, public :: ncdco = 46

 real :: co_temptab(nTco), co_colntab(ncdco)
 real :: co_L0(nTco), dTco_L0(nTco)
 real :: co_lte(ncdco,nTco), co_n05(ncdco,nTco), co_alp(ncdco,nTco)
 real :: dTco_lte(ncdco,nTco), dTco_n05(ncdco,nTco), dTco_alp(ncdco,nTco)

! These variables must be initialised during problem setup
! (in Phantom these appear in the input file when cooling is set,
!  here we give them sensible default values)
!
! Total abundances of C, O, Si: Sembach et al. (2000)
 real, public :: abundc  = 1.4d-4
 real, public :: abundo  = 3.2d-4
 real, public :: abundsi = 1.5d-5
 real, public :: abunde  = 2.d-4

! Strength of UV field (in Habing units)
 real :: uv_field_strength = 1.d0

! Dust temperature (in K)
 real :: tdust = 1.d1

! [At lower metallicities, I [SG] generally assume for
!  simplicity that this scales linearly with metallicity,
!   but it need not do so]
!
 real :: dust_to_gas_ratio = 1.d0
!
! Visual extinction (A_V) per unit column density (in cm^-2)
!
 real, public :: AV_conversion_factor = 5.348d-22
!
! Cosmic ray ionization rate of HI (in s^-1)
!
 real, public :: cosmic_ray_ion_rate = 1.d-17
!
! Flag controlling treatment of photoelectric heating
! iphoto = 0 ==> optically thin gas assumed
! iphoto = 1 ==> approximate treatment of effects of extinction
! [See cool_func.F for details]
!
 integer :: iphoto = 1
!
! Flag controlling which atomic cooling function is used.
! iflag_atom = 1 is the appropriate choice for the Galactic ISM
! [iflag_atom = 2 is used for Z=0 gas]

 integer :: iflag_atom = 1

! Distance measurements needed for chemistry

 real(kind=8), public :: dlq = 3.086d19
 real(kind=8), public :: dphot0 = 1.0801d20
 real(kind=8), public :: dchem = 3.086d20
 integer, public       :: dphotflag = 0
 private

contains

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_h2cooling(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(dlq,'dlq','distance for column density in cooling function',iunit)
 call write_inopt(dphot0,'dphot','photodissociation distance used for CO/H2',iunit)
 call write_inopt(dphotflag,'dphotflag','photodissociation distance static or radially adaptive (0/1)',iunit)
 call write_inopt(dchem,'dchem','distance for chemistry of HI',iunit)
 call write_inopt(abundc,'abundc','Carbon abundance',iunit)
 call write_inopt(abundo,'abundo','Oxygen abundance',iunit)
 call write_inopt(abundsi,'abundsi','Silicon abundance',iunit)
 call write_inopt(abunde,'abunde','electron abundance',iunit)
 call write_inopt(uv_field_strength,'uv_field_strength',&
                  'Strength of UV field (in Habing units)',iunit)
 call write_inopt(dust_to_gas_ratio,'dust_to_gas_ratio','dust to gas ratio',iunit)
 call write_inopt(AV_conversion_factor,'AV_conversion_factor',&
                  'Extinction per unit column density (cm^-2)',iunit)
 call write_inopt(cosmic_ray_ion_rate,'cosmic_ray_ion_rate',&
                  'Cosmic ray ionisation rate of H1 (in s^-1)',iunit)
 call write_inopt(iphoto,'iphoto','Photoelectric heating treatment (0=optically thin, 1=w/extinction)',iunit)
 call write_inopt(iflag_atom,'iflag_atom','Which atomic cooling (1:Gal ISM, 2:Z=0 gas)',iunit)

end subroutine write_options_h2cooling

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_h2cooling(name,valstring,imatch,igotall,ierr)
 use part, only:h2chemistry
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 imatch  = .true.
 igotall = .true. ! none of the cooling options are compulsory
 select case(trim(name))
 case('dlq')
    read(valstring,*,iostat=ierr) dlq
 case('dphot')
    read(valstring,*,iostat=ierr) dphot0
 case('dphotflag')
    read(valstring,*,iostat=ierr) dphotflag
 case('dchem')
    read(valstring,*,iostat=ierr) dchem
 case('abundc')
    read(valstring,*,iostat=ierr) abundc
 case('abundo')
    read(valstring,*,iostat=ierr) abundo
 case('abundsi')
    read(valstring,*,iostat=ierr) abundsi
 case('abunde')
    read(valstring,*,iostat=ierr) abunde
 case('uv_field_strength')
    read(valstring,*,iostat=ierr) uv_field_strength
 case('dust_to_gas_ratio')
    read(valstring,*,iostat=ierr) dust_to_gas_ratio
 case('AV_conversion_factor')
    read(valstring,*,iostat=ierr) AV_conversion_factor
 case('cosmic_ray_ion_rate')
    read(valstring,*,iostat=ierr) cosmic_ray_ion_rate
 case('iphoto')
    read(valstring,*,iostat=ierr) iphoto
 case('iflag_atom')
    read(valstring,*,iostat=ierr) iflag_atom
 case default
    imatch = .false.
 end select

end subroutine read_options_h2cooling

!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////            C O O L _ F U N C              \\\\\\\\\\
!
!=======================================================================
!
subroutine cool_func(temp, yn, dl, divv, abundances, ylam, rates)
!
!    Based on cool_h, written for ZEUS-3D by Michael D. Smith and
!    Georgi Pavlovski (Armagh Observatory, 2003) and substantially
!    modified by S. Glover (AMNH, 2003-2005, AIP 2006).
!
!    PURPOSE:  Compute the cooling function for the gas. Note that
!    the convention used here is that Lambda (the net _cooling_ rate)
!    is positive and has units of erg s^-1 cm^-3, so that for gas at
!    rest the equation for the internal energy density, e, can be
!    written as:
!                 de/dt = - Lambda
!    This means that any heating terms which we include in Lambda
!    (e.g. photoelectric heating) must be NEGATIVE; i.e. heating is
!    treated as negative cooling.
!
!    N.B. When used for cosmological simulations, it is necessary to
!         account for the effects of the CMB on the various heating &
!         cooling rates. In the case of the fine structure coolants,
!         this is done explicitly (and exactly) by including the
!         appropriate stimulated emission & absorption terms in the
!         level population equations. For other species (HD, CO, H2O),
!         this effect is not treated explicitly -- instead, we
!         approximate it by first calling cool_func with the gas
!         temperature and then with the CMB temperature, and then
!         subtracting the cooling rate given by the second call from
!         that given by the first. This approximation is of questionable
!         accuracy near T_cmb -- for instance, it can lead to ~50% errors
!         in the HD cooling rate -- but is probably sufficient for our
!         present purposes, given the other uncertainties
!
!  INPUT VARIABLES: temp       -- temperature
!                   yn         -- density
!                   divv       -- velocity divergence
!                   abundances -- chemical abundances
!
!  OUTPUT VARIABLES: Computed values of ylam & rates
!
!  LOCAL VARIABLES:  rates  -- parts of the cooling function
!                       cl[1-*] -- cooling table parameters
!                     dtcl[1-*] -- their change with T = T + dT
!
!  NB: the dt... terms are TEMPERATURE derivatives, not time derivatives
!
!  EXTERNALS:
!
!-----------------------------------------------------------------------
 use fs_data
 use mol_data
 use io, only:fatal
 real,         intent(in)  :: temp
 real,         intent(in)  :: yn
 real(kind=8), intent(in)  :: dl
 real,         intent(in)  :: divv
 real,         intent(in)  :: abundances(nabn)
 real,         intent(out) :: ylam
 real,         intent(out) :: rates(nrates)

 integer :: itemp
 real    :: dtemp

 real    :: ynh2      , ynh     , yne     , ynhp

 real    :: abh2      , abo       , &
            abcI      , abcII     , absiI   , absiII   , abe  , &
            abhp      , abhI      , abco

 real    :: h2var0    , h2var1

 real    :: oxc01, oxc02, &
            oxc10, oxc20, oxn1  , oxn2  , oxn0  , oxr01 , &
            oxr02, oxr12, oxr10 , oxr20 , oxr21 , oxa   , &
            oxb  , oxc  , oxc21 , oxc12

 real    :: cIc01, cIc02, cIc10 , &
            cIc20, cIn1 , cIn2  , cIn0  , cIr01 , &
            cIr02, cIr12, cIr10 , cIr20 , cIr21 , cIa   , &
            cIb  , cIc  , cIc21 , cIc12

 real    :: siIc01, siIc02,  siIc10,                         &
            siIc20, siIn1 , siIn2  , siIn0 , siIr01 ,        &
            siIr02, siIr12, siIr10 , siIr20 , siIr21 , siIa, &
            siIb  , siIc  , siIc21 , siIc12

 real    :: cIIc10,   cIIc01, cIIn1

 real    :: siIIc10, siIIc01, siIIn1

 real    :: eps       , PEvar0   , PEvar1   , PEvar2, &
            G_dust    , AV       , fdust    , NH_tot

 real    :: maximum_CO_column, N_co_eff, dv, temp_co, co_rot_L0, &
            co_rot_lte, co_rot_alpha, co_rot_n05, neff, co_rot_inv, &
            low_T_adjust, sigma_h2, v_e

 real    :: cl1       , cl2     , cl3     , cl6     , cl14  &
           , cl16      , cl17    , cl18    , cl19    , cl20  &
           , cl21      , cl22    , cl23    , cl24    , cl25  &
           , cl26      , cl27    , cl28    , cl29    , cl30  &
           , cl31      , cl32    , cl33    , cl34    , cl35  &
           , cl36      , cl37    , cl38    , cl39    , cl40  &
           , cl41      , cl42    , cl43    , cl44    , cl45  &
           , cl46      , cl47    , cl48    , cl49    , cl50  &
           , cl51      , cl52    , cl53    , cl54
!
 real    :: dtcl1     , dtcl2   , dtcl3   , dtcl6   , dtcl14  &
          , dtcl16    , dtcl17  , dtcl18  , dtcl19  , dtcl20  &
          , dtcl21    , dtcl22  , dtcl23  , dtcl24  , dtcl25  &
          , dtcl26    , dtcl27  , dtcl28  , dtcl29  , dtcl30  &
          , dtcl31    , dtcl32  , dtcl33  , dtcl34  , dtcl35  &
          , dtcl36    , dtcl37  , dtcl38  , dtcl39  , dtcl40  &
          , dtcl41    , dtcl42  , dtcl43  , dtcl44  , dtcl45  &
          , dtcl46    , dtcl47  , dtcl48  , dtcl49  , dtcl50  &
          , dtcl51    , dtcl52  , dtcl53  , dtcl54
 !
! ---------------------------------------------------------------------
!
! Read out tables.
!
 if (temp  <=  tmin) then
    itemp = 1
    dtemp = 0d0
 elseif (temp  >  tmax) then
    itemp = nmd
    dtemp = 0d0
 else
    itemp = int(log10(temp) / dtlog) + 1
    if (itemp  <=  0 .or. itemp  >  nmd) then
       !print*, 'Fatal error in cool_func.F', itemp, temp
       call fatal('cool_func','bad temperature in cooling function',var='temp',val=temp,ival=itemp)
    endif
    dtemp = temp - temptab(itemp)
 endif
!
 dtcl1  = dtcltab(1,itemp)
 dtcl2  = dtcltab(2,itemp)
 dtcl3  = dtcltab(3,itemp)
!
 dtcl6  = dtcltab(6,itemp)
!
 dtcl14 = dtcltab(14,itemp)
!
 dtcl16 = dtcltab(16,itemp)
 dtcl17 = dtcltab(17,itemp)
 dtcl18 = dtcltab(18,itemp)
 dtcl19 = dtcltab(19,itemp)
 dtcl20 = dtcltab(20,itemp)
 dtcl21 = dtcltab(21,itemp)
 dtcl22 = dtcltab(22,itemp)
 dtcl23 = dtcltab(23,itemp)
 dtcl24 = dtcltab(24,itemp)
 dtcl25 = dtcltab(25,itemp)
 dtcl26 = dtcltab(26,itemp)
 dtcl27 = dtcltab(27,itemp)
 dtcl28 = dtcltab(28,itemp)
 dtcl29 = dtcltab(29,itemp)
 dtcl30 = dtcltab(30,itemp)
 dtcl31 = dtcltab(31,itemp)
 dtcl32 = dtcltab(32,itemp)
 dtcl33 = dtcltab(33,itemp)
 dtcl34 = dtcltab(34,itemp)
 dtcl35 = dtcltab(35,itemp)
 dtcl36 = dtcltab(36,itemp)
 dtcl37 = dtcltab(37,itemp)
 dtcl38 = dtcltab(38,itemp)
 dtcl39 = dtcltab(39,itemp)
 dtcl40 = dtcltab(40,itemp)
 dtcl41 = dtcltab(41,itemp)
 dtcl42 = dtcltab(42,itemp)
 dtcl43 = dtcltab(43,itemp)
 dtcl44 = dtcltab(44,itemp)
 dtcl45 = dtcltab(45,itemp)
 dtcl46 = dtcltab(46,itemp)
 dtcl47 = dtcltab(47,itemp)
 dtcl48 = dtcltab(48,itemp)
 dtcl49 = dtcltab(49,itemp)
 dtcl50 = dtcltab(50,itemp)
 dtcl51 = dtcltab(51,itemp)
 dtcl52 = dtcltab(52,itemp)
 dtcl53 = dtcltab(53,itemp)
 dtcl54 = dtcltab(54,itemp)

 cl1  = cltab(1,itemp) + dtemp * dtcl1
 cl2  = cltab(2,itemp) + dtemp * dtcl2
 cl3  = cltab(3,itemp) + dtemp * dtcl3
!
 cl6  = cltab(6,itemp) + dtemp * dtcl6
!
 cl14 = cltab(14,itemp) + dtemp * dtcl14
!
 cl16 = cltab(16,itemp) + dtemp * dtcl16
 cl17 = cltab(17,itemp) + dtemp * dtcl17
 cl18 = cltab(18,itemp) + dtemp * dtcl18
 cl19 = cltab(19,itemp) + dtemp * dtcl19
 cl20 = cltab(20,itemp) + dtemp * dtcl20
 cl21 = cltab(21,itemp) + dtemp * dtcl21
 cl22 = cltab(22,itemp) + dtemp * dtcl22
 cl23 = cltab(23,itemp) + dtemp * dtcl23
 cl24 = cltab(24,itemp) + dtemp * dtcl24
 cl25 = cltab(25,itemp) + dtemp * dtcl25
 cl26 = cltab(26,itemp) + dtemp * dtcl26
 cl27 = cltab(27,itemp) + dtemp * dtcl27
 cl28 = cltab(28,itemp) + dtemp * dtcl28
 cl29 = cltab(29,itemp) + dtemp * dtcl29
 cl30 = cltab(30,itemp) + dtemp * dtcl30
 cl31 = cltab(31,itemp) + dtemp * dtcl31
 cl32 = cltab(32,itemp) + dtemp * dtcl32
 cl33 = cltab(33,itemp) + dtemp * dtcl33
 cl34 = cltab(34,itemp) + dtemp * dtcl34
 cl35 = cltab(35,itemp) + dtemp * dtcl35
 cl36 = cltab(36,itemp) + dtemp * dtcl36
 cl37 = cltab(37,itemp) + dtemp * dtcl37
 cl38 = cltab(38,itemp) + dtemp * dtcl38
 cl39 = cltab(39,itemp) + dtemp * dtcl39
 cl40 = cltab(40,itemp) + dtemp * dtcl40
 cl41 = cltab(41,itemp) + dtemp * dtcl41
 cl42 = cltab(42,itemp) + dtemp * dtcl42
 cl43 = cltab(43,itemp) + dtemp * dtcl43
 cl44 = cltab(44,itemp) + dtemp * dtcl44
 cl45 = cltab(45,itemp) + dtemp * dtcl45
 cl46 = cltab(46,itemp) + dtemp * dtcl46
 cl47 = cltab(47,itemp) + dtemp * dtcl47
 cl48 = cltab(48,itemp) + dtemp * dtcl48
 cl49 = cltab(49,itemp) + dtemp * dtcl49
 cl50 = cltab(50,itemp) + dtemp * dtcl50
 cl51 = cltab(51,itemp) + dtemp * dtcl51
 cl52 = cltab(52,itemp) + dtemp * dtcl52
 cl53 = cltab(53,itemp) + dtemp * dtcl53
 cl54 = cltab(54,itemp) + dtemp * dtcl54
!
! Set abundances
!
 abh2   = abundances(1)
 abhI   = abundances(2)
 abe    = abundances(3)
 abhp   = abundances(4)
 abo    = abundances(5)
 abcI   = abundances(6)
 abcII  = abundances(7)
 absiI  = abundances(8)
 absiII = abundances(9)
 abco   = abundances(10)
!
! Compute useful auxiliary variables

 ynh2 = abh2 * yn
 ynh  = abHI * yn
 yne  = abe  * yn
 ynhp = abhp * yn
!
! Compute effective column density for CO rotational cooling rate
! (See eq. 4 of Neufeld & Kaufman, 1993, ApJ, 418, 263). Tabulated values
! for the effective column density are in terms of cm^-2 / (km s^-1), but
! divv is in units of cm s^-1, so we need to do a unit conversion here.
!
! If divv is very small (or zero), then we set the effective column
! density to the largest tabulated value

 if (abco  >  1d-4 * abundo) then
! Don't bother to compute CO cooling rate when CO abundance very small
    maximum_CO_column = co_colntab(ncdco)
    if (abs(divv) < tiny(divv)) then
       N_co_eff   = maximum_CO_column
    else
       dv = 1d-5 * abs(divv)
       N_co_eff   = log10(abco * yn / dv)
       if (N_co_eff  >  maximum_CO_column) then
          N_co_eff = maximum_CO_column
       endif
    endif
!
! If current temperature is less than the smallest tabulated value, we
! look up rates for smallest tabulated value and then reduce final
! cooling rate by an appropriate exponential factor
!
    if (temp  <=  co_temptab(1)) then
       temp_co = co_temptab(1)
    else
       temp_co = temp
    endif
!
    call co_cool(temp_co, N_co_eff, co_rot_L0, co_rot_lte, co_rot_alpha, co_rot_n05)
 else
    co_rot_L0 = 0.
    co_rot_lte = 0.
    co_rot_n05 = 0.
    co_rot_alpha = 0.
 endif
!
! (R1) -- gas-grain cooling-heating -- dust:gas ratio already incorporated
!         into rate coefficient in init_h2cooling
!
 rates(1) = cl14 * (temp - tdust) * yn**2
!
!
! (R2) -- H2 (vr) cooling
!
 if (ynh2 == 0d0) then
    rates(2) = 0d0
 else
    h2var0   = ynh * cl2 + ynh2 * cl3
    h2var1   = cl1 + h2var0
    rates(2) = ynh2 * cl1 * h2var0 / h2var1
 endif
!
! (R3) -- atomic cooling
!
 rates(3) = cl6 * yn**2 + cl52 * ynh * yne
!
! (R4) -- cosmic ray heating; independent of gas temperature
!       -- following Goldsmith & Langer (1978), we assume that
!          each ionization deposits 20eV as heat
!
 rates(4) =  -3.2d-11 * (ynh2 + ynh) * cosmic_ray_ion_rate
!
! (R5-6) -- photoelectric heating, PAH recombination cooling
!
! If there's no UV field, or if the electron density is very low (in which
! case the photoheating efficiency will also be very low), then we set the
! rates to zero. Otherwise, we compute the heating rate using the
! Bakes & Tielens (1994) formula (as modified by Wolfire et al, 2003),
! and the recombination cooling using the formula from Wolfire et al 2003.
!
! Treatment of dust attenuation follows Bergin et al (2004, ApJ,  612, 921)
! Calculation of A_V uses the conversion factor that was initialized during
! problem setup. This treatment neglects any variation of the photoelectric
! heating efficiency with increasing A_V.
!
 if (uv_field_strength  <  tiny(0.d0) .or. dust_to_gas_ratio  <  tiny(0.d0)) then
    rates(5) = 0d0
    rates(6) = 0d0
 elseif (yne  <  1d-9 * uv_field_strength * cl16) then
    rates(5) = 0d0
    rates(6) = 0d0
 else
    if (iphoto == 0) then
       G_dust = uv_field_strength
    else
       NH_tot = 0.5d0 * dl * yn * (2d0 * abh2 + abhI + abhp)
       AV     = AV_conversion_factor * dust_to_gas_ratio * NH_tot
       fdust  = exp(-2.5d0 * AV)
       G_dust = uv_field_strength * fdust
    endif

    PEvar0 = G_dust * cl16 / yne
    PEvar1 = (1d0 + 4d-3 * PEvar0**0.73)
    PEvar2 = (1d0 + 2d-4 * PEvar0)

    eps = (4.9d-2 / PEvar1) + (cl17 / PEvar2)

! Photoelectric heating:
    rates(5) = -1.3d-24 * eps * G_dust * yn * dust_to_gas_ratio
! Recombination cooling:
    rates(6) = cl54 * PEvar0**cl53 * yne * yn * dust_to_gas_ratio
 endif
!
! (R7) -- OI fine-structure cooling
!
! Total collisional rates:
!
 oxc10  = cl18 * ynh + cl21 * ynh2 + cl24 * yne + cl27 * ynhp
 oxc20  = cl19 * ynh + cl22 * ynh2 + cl25 * yne + cl28 * ynhp
 oxc21  = cl20 * ynh + cl23 * ynh2 + cl26 * yne + cl29 * ynhp

 if (abo  <=  1d-5 * abundo) then
    rates(7) = 0d0
 elseif (oxc10 == 0d0 .and. oxc20 == 0d0 .and. oxc21 == 0d0) then
    rates(7) = 0d0
 else
    oxa =  exp(-oxe10 / (kboltz * temp))
    oxb =  exp(-oxe20 / (kboltz * temp))
    oxc =  exp(-oxe21 / (kboltz * temp))

    oxc01  = 0.6d0 * oxc10 * oxa
    oxc02  = 0.2d0 * oxc20 * oxb
    oxc12  = (1d0 / 3d0) * oxc21 * oxc
!
! Total transition rates:
!
    oxR01  = oxc01
    oxR02  = oxc02
    oxR12  = oxc12
    oxR10  = oxc10 + oxa10
    oxR20  = oxc20 + oxa20
    oxR21  = oxc21 + oxa21
!
    call three_level_pops(oxR01, oxR02, oxR12, oxR10, oxR20, oxR21, oxn0, oxn1, oxn2)
!
! Total emitted energy:
!
    rates(7) = (oxa10 * oxe10 * oxn1 + (oxa20 * oxe20 + oxa21 * oxe21) * oxn2) * abo * yn
 endif
!
! (R8) --  CI fine-structure cooling
!
! Collisional rates:
!
 cIc10  = cl30 * ynh + cl33 * ynh2 + cl36 * yne + cl39 * ynhp
 cIc20  = cl31 * ynh + cl34 * ynh2 + cl37 * yne + cl40 * ynhp
 cIc21  = cl32 * ynh + cl35 * ynh2 + cl38 * yne + cl41 * ynhp

 if (abcI  <=  1d-5 * abundc) then
    rates(8) = 0d0
 elseif (cIc10 == 0d0 .and. cIc20 == 0d0 .and. cIc21 == 0d0) then
    rates(8) = 0d0
 else
    cIa =  exp(-cIe10 / (kboltz * temp))
    cIb =  exp(-cIe20 / (kboltz * temp))
    cIc =  exp(-cIe21 / (kboltz * temp))

    cIc01  = 3d0 * cIc10 * cIa
    cIc02  = 5d0 * cIc20 * cIb
    cIc12  = (5d0 / 3d0) * cIc21 * cIc
!
! Total transition rates:
!
    cIR01  = cIc01
    cIR02  = cIc02
    cIR12  = cIc12
    cIR10  = cIc10 + cIa10
    cIR20  = cIc20 + cIa20
    cIR21  = cIc21 + cIa21
!
    call three_level_pops(cIR01, cIR02, cIR12, cIR10, cIR20, cIR21, cIn0, cIn1, cIn2)
!
! Total emitted energy:
!
    rates(8) = (cIa10 * cIe10 * cIn1 + (cIa20 * cIe20 + cIa21 * cIe21) * cIn2) * abcI * yn
 endif
!
! (R9) --  SiI fine-structure cooling
!
! Proton rates (from HM89) are constant and so there's no point
! tabulating them in init_h2cooling
!
 siIc10 = cl42 * ynh + 7.2d-9 * ynhp
 siIc20 = cl43 * ynh + 7.2d-9 * ynhp
 siIc21 = cl44 * ynh + 2.2d-8 * ynhp

 if (absiI  <=  1d-5 * abundsi) then
    rates(9) = 0d0
 elseif (siIc10 == 0d0 .and. siIc20 == 0d0 .and. siIc21 == 0d0) then
    rates(9) = 0d0
 else
    siIa =  exp(-siIe10 / (kboltz * temp))
    siIb =  exp(-siIe20 / (kboltz * temp))
    siIc =  exp(-siIe21 / (kboltz * temp))
!
    siIc01  = 3d0 * siIc10 * siIa
    siIc02  = 5d0 * siIc20 * siIb
    siIc12  = (5d0 / 3d0) * siIc21 * siIc
!
! Total transition rates:
!
    siIR01  = siIc01
    siIR02  = siIc02
    siIR12  = siIc12
    siIR10  = siIc10 + siIa10
    siIR20  = siIc20 + siIa20
    siIR21  = siIc21 + siIa21
!
    call three_level_pops(siIR01, siIR02, siIR12, siIR10, siIR20, siIR21, siIn0, siIn1, siIn2)
!
! Total emitted energy:
!
    rates(9) = (siIa10 * siIe10 * siIn1 + (siIa20 * siIe20 + siIa21 * siIe21) * siIn2) * absiI * yn
 endif
!
! (R10) -- CII fine-structure cooling
!
 cIIc10 = cl45 * ynh + cl46 * ynh2 + cl47 * yne
 cIIc01 = cl48 * cIIc10

 if (cIIc10 == 0d0 .or. abcII  <=  1d-5 * abundc) then
    rates(10) = 0d0
 else
!
    cIIn1 = cIIc01 / (cIIc01 + cIIc10 + cIIa10)
!
    rates(10) = cIIa10 * cIIe10 * cIIn1 * abcII * yn
 endif
!
! (R11) -- SiII fine-structure cooling
!
 siIIc10 = cl49 * ynh + cl50 * yne
 siIIc01 = cl51 * siIIc10

 if (siIIc10 == 0d0 .or. absiII  <=  1d-5 * abundsi) then
    rates(11) = 0d0
 else
!
    siIIn1 = siIIc01 / (siIIc10 + siIIa10 + siIIc01)
!
    rates(11) = siIIa10 * siIIe10 * siIIn1 * absiII * yn
 endif
!
! (R12) -- CO rotational cooling
!
 sigma_h2 = 3.3d-16 * (temp / 1d3)**(-0.25d0)
 v_e = 1.03d4 * sqrt(temp)
 neff = ynh2 + dsqrt(2d0) * (2.3d-15 / sigma_h2) * ynh &
      + (1.3d-8 / (sigma_h2 * v_e)) * yne
!
 if (abco  <=  1d-4 * abundo .or. neff == 0d0) then
    rates(12) = 0d0
 else
    co_rot_inv = (1d0 / co_rot_L0) + (neff / co_rot_lte) + (1d0 / co_rot_L0) &
              * (1d0 - co_rot_n05 * co_rot_L0  / co_rot_lte) &
              * (neff / co_rot_n05)**co_rot_alpha
!
    rates(12) = abco * neff * yn / co_rot_inv

    if (temp  <=  co_temptab(1)) then
! We don't have accurate data below 5K, but we don't necessarily want to impose a sharp floor there either.
! As a compromise, we reduce the cooling rate exponentially below 5 K.
!
       low_T_adjust = dexp(-1d1 / temp) / dexp(-1d1 / 5d0)
       rates(12) = rates(12) * low_T_adjust
    endif
 endif
!
! Benchmarking suggests that writing this out explicitly is more efficient
! than using a loop (although this is probably only true if the compiler
! optimization is poor).
!
 ylam = rates(1)  + rates(2) + rates(3) + rates(4) + rates(5) + &
        rates(6)  + rates(7) + rates(8) + rates(9) + rates(10) + &
        rates(11) + rates(12)
!
 return
end subroutine cool_func
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////             C O O L _ F U N C             \\\\\\\\\\
!
!=======================================================================

!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
!
!=======================================================================
!
pure subroutine three_level_pops(r01, r02, r12, r10, r20, r21, n0, n1, n2)
 real, intent(in)  :: r01, r02, r12, r10, r20, r21
 real, intent(out) :: n0, n1, n2
 real :: a1 , a2 , a3 , b1 , b2 , b3
!
! If excitation rates are negligibly small, then we assume that all
! of the atoms are in level 0:
!
 if (r01  <  tiny(0d0) .and. r02  <  tiny(0d0)) then
    n0 = 1d0
    n1 = 0d0
    n2 = 0d0
    return
 endif

 a1 = r01 + r02
 a2 = -r10
 a3 = -r20
 b1 = r01
 b2 = -(r10 + r12)
 b3 = r21
!
 n2 = -a1 * (a1 * b2 - b1 * a2) / ((a1 - a2) *  &
      (a1 * b3 - b1 * a3) - (a1 - a3) *  &
      (a1 * b2 - b1 * a2))
!
 n1 = (a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2
!
 n0 = 1d0 - n1 - n2
!
 return
end subroutine three_level_pops
!=======================================================================
!
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
!    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
!
!=======================================================================

!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              C O O L I N M O              \\\\\\\\\\
!
!=======================================================================
!
subroutine init_h2cooling
!
!    Based on an original routine by G. Suttner (University Wuerzburg, 1995)
!    OI, CI cooling added by M. D. Smith (Armagh Observatory, 2000-2001)
!    Extensively rewritten by S. Glover (AMNH, 2003-2005)
!    -- fine structure cooling from CII, SiI, SiII added
!    -- fine structure cooling treatment is now exact
!    -- new H2 cooling function used (taken from Le Bourlot et al 1999)
!    -- high T metal-free cooling from Sutherland & Dopita (1993)
!    -- many other minor changes
!
!  PURPOSE:  Tabulate cooling function for various ions, atoms and molecules
!            before starting run
!
!  REFERENCES:
!       FL77  -- Flower & Launay, 1977, J. Phys. B, 10, 3673
!       B81   -- Black, 1981, MNRAS, 197, 553
!       K86   -- Keenan et al, 1986, MNRAS, 220, 571
!       JBK87 -- Johnson, Burke & Kingston, 1987, J. Phys B, 20, 2553
!       HM89  -- Hollenbach & McKee, 1989, ApJ, 342, 306
!       RLB90 -- Roueff & Le Bourlot, 1990, A&A, 236, 515
!       P90   -- Pequignot, 1990, A&A, 231, 499
!       R90   -- Roueff, 1990, A&A, 234, 567
!       DK91  -- Dufton & Kingston, 1991, MNRAS, 248, 827
!       S91   -- Schroeder et al, 1991, J. Phys. B, 24, 2487
!       C92   -- Cen, 1992, ApJS, 78, 341
!       SD93  -- Sutherland & Dopita, 1993, ApJS, 88, 253
!       BT94  -- Bakes & Tielens, 1994, ApJ, 427, 822
!       W95   -- Wolfire et al, 1995, ApJ,  443, 152
!       WBV96 -- Warin, Benayoun & Viala, 1996, A&A, 308, 535
!       BBT98 -- Bell, Berrington & Thomas, 1998, MNRAS, 293, L83
!       LPF99 -- Le Bourlot, Pineau des Forets & Flower, 1999, MNRAS, 305, 802
!       WB02  -- Wilson & Bell, 2002, MNRAS, 337, 1027
!       W03   -- Wolfire et al, 2003, ApJ, 587, 278
!       AKD07 -- Abrahamsson, Krems & Dalgarno, 2007, ApJ, 654, 1171
!
 use fs_data,     only:cIe10,cIe20,cIe21,oxe10,oxe20,oxe21
 use mol_data,    only:nh2data,h2_h_rate,h2_h2_rate,h2_lte,h2_temp,nco_temp,nco_column,nco_data, &
                       co_temp, co_column, co_data_L0, co_data_LTE, co_data_n05, co_data_alp
 use splineutils, only:spline_eval
!
 integer :: i, j, itemp
!
 integer, parameter :: natom = 81
 real :: coolatom(natom), coolatom_temp(natom)

 integer, parameter :: natom2 = 76
 real :: ca2(natom2), ca2_temp(natom2)
!
 real :: rate0(nmd), rate1(nmd), rate2(nmd), rate3(nmd)
!
! Raw data -- extracted from DATA tables in mol_data.f90
!
 real :: co_lte_raw(nco_temp), co_n05_raw(nco_temp), co_alp_raw(nco_temp)
!
! Temporary variables used in fitting procedure for CO rotational cooling
!
 real :: co_lte_fit(nTco),   co_alp_fit(nTco),   co_n05_fit(nTco)
 real :: co_lte_fit2(ncdco), co_n05_fit2(ncdco), co_alp_fit2(ncdco)

 real :: co_lte_fxT(nco_column), co_n05_fxT(nco_column), co_alp_fxT(nco_column)

 real :: co_lte_smalltab(nco_column,nTco), co_n05_smalltab(nco_column,nTco), &
         co_alp_smalltab(nco_column,nTco)

 real :: temp    , temp2   , f       , gg       , hh      &
       , dtemp   , tinv    , tau     , tsqrt    , opratio &
       , fortho  , brem    , fpara   , atomic   , tloge   &
       , h2e20   , h2e31   , h2n2    , h2n3     , h2q02   &
       , h2q13   , phi_pah , tinth   , tinq     , tintq   &
       , tisqt   , tfix    , tfintq  , tfinq    , tfinth
!
! atomic cooling:  table of SD93 for zero radiation field, non-equil.,
! Fe = -0.5 (table 10)
!
 data coolatom /-68.,-67.,-66.,-65.,-64.,-63.,-62.,-61.,-60.,-59., &
      -58.,-57.,-56.,-55.,-54.,-53.,-52.,-51.,-50.,-49.,-48.,-47., &
      -46.,-45.,-44.,-43.,-42.,-41.,-40.,-39.,-38.,-37.,-36.,-35., &
      -34.,-33.,-32.,-31.,-30.,-26.67,-23.16,-22.83,-22.53,-22.30, &
      -22.20,-22.11,-22.00,-21.88,-21.77,-21.69,-21.60,-21.54,     &
      -21.53,-21.58,-21.73,-21.91,-22.00,-22.07,-22.23,-22.32,     &
      -22.35,-22.34,-22.43,-22.62,-22.75,-22.80,-22.83,-22.86,     &
      -22.85,-22.84,-22.85,-22.87,-22.87,-22.86,-22.85,-22.83,     &
      -22.79,-22.75,-22.70,-22.65,-22.60/
!
! atomic cooling #2: from SD93 -- zero radiation field, zero metals,
! collisional ionization equlibrium assumed. Does not include effects
! of HI collisional excitation. [Based on table 6 of SD93]
!
! NB. Data runs from log(T) = 4.3 - 8.0 in increments of 0.05.
!     First value corresponds to T = T_min, but is utterly arbitrary
!     -- this is here purely to make the spline fitting somewhat
!     easier, but as the actual entries for T < 10**4.3 K are handled
!     specially below, it should not significantly affect the computed
!     rates.
!
 data ca2 /-60.00, -23.38, -23.17, -23.06, -22.95, -22.89, -22.89, &
           -22.93, -22.93, -22.86, -22.69, -22.47, -22.27, -22.18, &
           -22.19, -22.27, -22.36, -22.46, -22.56, -22.66, -22.74, &
           -22.82, -22.89, -22.95, -23.01, -23.05, -23.09, -23.13, &
           -23.16, -23.18, -23.20, -23.22, -23.23, -23.24, -23.25, &
           -23.25, -23.25, -23.25, -23.25, -23.25, -23.24, -23.23, &
           -23.23, -23.23, -23.21, -23.20, -23.19, -23.18, -23.16, &
           -23.15, -23.13, -23.11, -23.10, -23.08, -23.06, -23.05, &
           -23.03, -23.01, -22.99, -22.97, -22.95, -22.93, -22.91, &
           -22.89, -22.86, -22.84, -22.82, -22.80, -22.78, -22.76, &
           -22.74, -22.71, -22.69, -22.67, -22.65, -22.62/
!
! H2 ortho-para ratio (N.B. Only used for OI fine structure cooling)
!
 opratio = 2.4d0
 fortho  = opratio / (1.d0 + opratio)
 fpara   = 1.d0 / (1.d0 + opratio)
!
! establish temperature table
!
 tmin   = 1.d0
 tmax   = 1.d8
 dtlog  = log10(tmax) / (nmd - 1)
!
 do itemp = 1, nmd
    temptab(itemp) = 1.D1**( (itemp - 1) * dtlog )
 enddo
!
! tabulate cooling functions for each temp
!
 over_temp: do itemp = 1, nmd
    temp  = temptab(itemp)
    temp2 = temp * 1d-2
    tau   = temp * 1d-3 + 1d0
    tinv  = 1.d0 / temp
    tinth = tinv**(1d0/3d0)
    tinq  = tinv**0.25d0
    tintq = tinv**0.75d0
    tsqrt = sqrt(temp)
    tisqt = 1.d0 / tsqrt
    tloge = log(temp)
!
! cl[1-3]: H2 cooling (LPF99)
!
!   Values in range 100 -- 10000K are handled at end.
!
!   Rate at T > 10000K is arbitrarily fixed to be the same as at 10000K
!   (this is incorrect, but simple, and in any case Lyman-alpha will
!    dominate at these temperatures)
!
!   Below 100K, we assume that all cooling comes from the J=2-0, 3-1
!   transitions in the vibrational ground state.
!
!   All rates assume an ortho:para ratio of 3:1
!
! (cl1) -- LTE rate
!
    h2e20 = 510.06d0
    h2e31 = 844.94d0
!
    if (temp  <  5d0) then
       cltab(1,itemp) = tiny(cltab)
    else if (temp  <  1d2) then
       h2n2 = 0.25d0 * 5d0 * exp(-h2e20 / temp)
       h2n3 = 0.75d0 * (7d0 / 3d0) * exp(-h2e31 / temp)
       f    = 2.94d-11 * h2e20 * h2n2 * kboltz + &
              4.76d-10 * h2e31 * h2n3 * kboltz
       cltab(1,itemp) = max(f, tiny(f))
    else if (temp  >  1d4) then
       cltab(1,itemp) = 1d1**(-18.253d0)
    endif
!
! (cl2) -- H rate -- 0->2 rate tweaked slightly to ensure that we match up
!                    properly with the tabulated rates at 100K
!
    h2q02 = 1d1**(-9.121d0 - (3.983d0 / tau) - (1.844d-1 / tau**2)) &
            * 5d0 * exp(-h2e20 / temp)
!
    h2q13 = 1d1**(-9.493d0 - (4.435d0 / tau) + (2.707d-1 / tau**2)) &
            * (7d0 / 3d0) * exp(-h2e31 / temp)
!
    if (temp  <  5d0) then
       cltab(2,itemp) = tiny(cltab)
    else if (temp  <  1d2) then
       cltab(2,itemp) = (fpara  * h2q02 * h2e20 + fortho * h2q13 * h2e31) * kboltz
    else if (temp  >  1d4) then
       cltab(2,itemp) = 1d1**(-21.943d0)
    endif
!
! (cl3) -- H2 rate
!
    h2q02 = 1d1**(-9.946d0 - (2.688d0 / tau) + (2.020d-1 / tau**2)) &
            * 5d0 * exp(-h2e20 / temp)
!
    h2q13 = 1d1**(-9.785d0 - (3.024d0 / tau) + (2.930d-1 / tau**2)) &
            * (7d0 / 3d0) * exp(-h2e31 / temp)
!
    if (temp  <  5d0) then
       cltab(3,itemp) = tiny(cltab)
    else if (temp  <  1d2) then
       cltab(3,itemp) = (fpara  * h2q02 * h2e20  +  fortho * h2q13 * h2e31) * kboltz
    else if (temp  >  1d4) then
       cltab(3,itemp) = 1d1**(-22.758d0)
    endif
!
! (cl4, cl5) -- currently unused
!
! (cl6) --  the atomic cooling function - this is computed by fitting
!           a cubic spline to the data specified in coolatom above, and
!           so is calculated after the main loop is done
!
! (cl7 -- cl13) -- unused in this version of the code
!
! (cl14) -- gas-grain cooling (HM89, eqn 2.15)
!
    gg              = 1d0 - 0.8 * exp(-75d0 * tinv)
!
    cltab(14,itemp) = 3.8d-33 * tsqrt * gg * dust_to_gas_ratio
!
! (cl15) --  unused in this version of the code
!
! (cl[16-17]) -- Photoelectric heating
!
! [phi_pah is an adjustable parameter introduced by W03. They quote
!  a value of 0.5, which we also adopt here. Note that the
!  photoelectric heating rate given by BT94 and W95 corresponds to
!  the case phi_pah = 1.0]
!
    phi_pah = 0.5d0

    cltab(16,itemp) = sqrt(temp) / phi_pah
!
    cltab(17,itemp) = 3.7e-2 * (temp / 1d4)**0.7d0
!
! Fine structure cooling -- three-level atoms
!
! We use the convention that the ground state is level 0, the first
! excited state is level 1, and the second excited state is level 2.
!
! We tabulate here the collisional de-excitation rates for all relevent
! collision partners. Everything else is taken care of in the cooling
! function itself at run time.
!
! (cl[18-29]): OI fine-structure lines
!
! Collisional de-excitation rates:
!
! HI: taken from AKD07 below 1000K, extended to higher temperatures
!     with a simple power-law extrapolation
!
! 1 -> 0
    if (temp  <  5d0) then
       tfix = 5d0
       tfintq = 1d0 / tfix**0.75d0
       cltab(18, itemp) = (5d-11 / 3d0) * exp(4.581  &
                        - 156.118       * tfintq     &
                        + 2679.979      * tfintq**2  &
                        - 78996.962     * tfintq**3  &
                        + 1308323.468   * tfintq**4  &
                        - 13011761.861  * tfintq**5  &
                        + 71010784.971  * tfintq**6  &
                        - 162826621.855 * tfintq**7) &
                        * exp(oxe10 / (kboltz * tfix))
    elseif (temp  <  1d3) then
       cltab(18, itemp) = (5d-11 / 3d0) * exp(4.581  &
                        - 156.118       * tintq      &
                        + 2679.979      * tintq**2   &
                        - 78996.962     * tintq**3   &
                        + 1308323.468   * tintq**4   &
                        - 13011761.861  * tintq**5   &
                        + 71010784.971  * tintq**6   &
                        - 162826621.855 * tintq**7)  &
                        * exp(oxe10 / (kboltz * temp))
    else
       cltab(18, itemp) = 6.81d-11 * temp**0.376d0
    endif
! 2 -> 0
    if (temp  <  5d0) then
       tfix = 5d0
       tfintq = 1d0 / tfix**0.75d0
       cltab(19, itemp) = 5d-11 * exp(3.297          &
                        - 168.382       * tfintq     &
                        + 1844.099      * tfintq**2  &
                        - 68362.889     * tfintq**3  &
                        + 1376864.737   * tfintq**4  &
                        - 17964610.169  * tfintq**5  &
                        + 134374927.808 * tfintq**6  &
                        - 430107587.886 * tfintq**7) &
                        * exp(oxe20 / (kboltz * tfix))
    elseif (temp  <  1d3) then
       cltab(19, itemp) = 5d-11 * exp(3.297         &
                        - 168.382       * tintq     &
                        + 1844.099      * tintq**2  &
                        - 68362.889     * tintq**3  &
                        + 1376864.737   * tintq**4  &
                        - 17964610.169  * tintq**5  &
                        + 134374927.808 * tintq**6  &
                        - 430107587.886 * tintq**7) &
                        * exp(oxe20 / (kboltz * temp))
    else
       cltab(19, itemp) = 6.34d-11 * temp**0.36d0
    endif
!
! 2 -> 1:  Low T extrapolation here is necessary because the AKD07
!          fitting function blows up
!
    if (temp  <  5d1) then
       cltab(20, itemp) = 2.62d-12 * temp**0.74d0
    elseif (temp  <  1d3) then
       cltab(20, itemp) = 3d-11 * exp(3.437     &
                        + 17.443    * tisqt     &
                        - 618.761   * tisqt**2  &
                        + 3757.156  * tisqt**3  &
                        - 12736.468 * tisqt**4  &
                        + 22785.266 * tisqt**5  &
                        - 22759.228 * tisqt**6  &
                        + 12668.261 * tisqt**7) &
                        * exp(oxe21 / (kboltz * temp))
    else
       cltab(20,itemp) = 3.61d-10 * temp**0.158d0
    endif
!
! H2 rates supplied by Flower (priv. comm.)
! (NB. If the H2 rate is based on the data from J92, then strictly
! speaking it is only applicable for T < 1500K; however, the rate
! doesn't misbehave too badly at higher T).
!
! H2 - ortho and para states must be accounted for separately
!
! 1 -> 0
    f               = fortho * 2.70d-11 * (temp**0.362)
    hh              = fpara  * 3.46d-11 * (temp**0.316)
    cltab(21,itemp) = f + hh
! 2 -> 0
    f               = fortho * 5.49d-11 * (temp**0.317)
    hh              = fpara  * 7.07d-11 * (temp**0.268)
    cltab(22,itemp) = f + hh
! 2 -> 1
    f               = fortho * 2.74d-14 * (temp**1.060)
    hh              = fpara  * 3.33d-15 * (temp**1.360)
    cltab(23,itemp) = f + hh
!
! Electron rate -- from my fits to BBT98.
!
! 1 -> 0
    cltab(24,itemp) = 5.12d-10  * (temp**(-0.075))
! 2 -> 0
    cltab(25,itemp) = 4.863d-10 * (temp**(-0.026))
! 2 -> 1
    cltab(26,itemp) = 1.082d-14 * (temp**(0.926))
!
! Proton rate -- from P90
!
! 1 -> 0
    if (temp  <  194) then
       cltab(27,itemp) = 6.38d-11 * (temp**0.40)
    else if (temp  <  3686) then
       cltab(27,itemp) = 7.75d-12 * (temp**0.80)
    else
       cltab(27,itemp) = 2.65d-10 * (temp**0.37)
    endif
! 2 -> 0
    if (temp  <  511) then
       cltab(28,itemp) = 6.10d-13 * (temp**1.10)
    else if (temp  <  7510) then
       cltab(28,itemp) = 2.12d-12 * (temp**0.90)
    else
       cltab(28,itemp) = 4.49d-10 * (temp**0.30)
    endif
! 2 -> 1
    if (temp  <  2090) then
       cltab(29,itemp) = 2.029d-11 * (temp**0.56)
    else
       cltab(29,itemp) = 3.434d-10 * (temp**0.19)
    endif
!
! (cl[30-41]): CI fine-structure lines
!
! Collisional de-excitation rates:
!
! HI: taken from AKD07 below 1000K, extended to higher temperatures
!     with a simple power-law extrapolation
!
! 1 -> 0
    if (temp  <  5d0) then
       tfix = 5d0
       tfinq = 1d0 / tfix**0.25d0
       cltab(30, itemp) =  (1d-11 / 3d0) * exp(3.6593 &
                        + 56.6023    * tfinq           &
                        - 802.9765   * tfinq**2        &
                        + 5025.1882  * tfinq**3        &
                        - 17874.4255 * tfinq**4        &
                        + 38343.6655 * tfinq**5        &
                        - 49249.4895 * tfinq**6        &
                        + 34789.3941 * tfinq**7        &
                        - 10390.9809 * tfinq**8)       &
                        * exp(cIe10 / (kboltz * tfix))
    elseif (temp  <  1d3) then
       cltab(30, itemp) = (1d-11 / 3d0) * exp(3.6593 &
                        + 56.6023    * tinq           &
                        - 802.9765   * tinq**2        &
                        + 5025.1882  * tinq**3        &
                        - 17874.4255 * tinq**4        &
                        + 38343.6655 * tinq**5        &
                        - 49249.4895 * tinq**6        &
                        + 34789.3941 * tinq**7        &
                        - 10390.9809 * tinq**8)       &
                        * exp(cIe10 / (kboltz * temp))
    else
       cltab(30, itemp) = 2.57d-11 * temp**0.31d0
    endif
! 2 -> 0

    if (temp  <  5d0) then
       tfix = 5d0
       tfinth = 1d0 / tfix**(1d0/3d0)
       cltab(31, itemp) = 2.0d-12 * exp(10.8377   &
                        - 173.4153    * tfinth     &
                        + 2024.0272   * tfinth**2  &
                        - 13391.6549  * tfinth**3  &
                        + 52198.5522  * tfinth**4  &
                        - 124518.3586 * tfinth**5  &
                        + 178182.5823 * tfinth**6  &
                        - 140970.6106 * tfinth**7  &
                        + 47504.5861  * tfinth**8) &
                        * exp(cIe20 / (kboltz * tfix))
    elseif (temp  <  1d3) then
       cltab(31, itemp) = 2.0d-12 * exp(10.8377  &
                        - 173.4153    * tinth     &
                        + 2024.0272   * tinth**2  &
                        - 13391.6549  * tinth**3  &
                        + 52198.5522  * tinth**4  &
                        - 124518.3586 * tinth**5  &
                        + 178182.5823 * tinth**6  &
                        - 140970.6106 * tinth**7  &
                        + 47504.5861  * tinth**8) &
                        * exp(cIe20 / (kboltz * temp))
    else
       cltab(31, itemp) = 1.69d-11 * temp**0.35d0
    endif
! 2 -> 1
    if (temp  <  5d0) then
       tfix = 5d0
       tfinq = 1d0 / tfix**0.25d0
       cltab(32, itemp) = 6.0d-12 * exp(15.8996 &
                        - 201.3030   * tfinq     &
                        + 1533.6164  * tfinq**2  &
                        - 6491.0083  * tfinq**3  &
                        + 15921.9239 * tfinq**4  &
                        - 22691.1632 * tfinq**5  &
                        + 17334.7529 * tfinq**6  &
                        - 5517.9360  * tfinq**7) &
                        * exp(cIe21 / (kboltz * tfix))
    elseif (temp  <  1d3) then
       cltab(32, itemp) = 6.0d-12 * exp(15.8996 &
                        - 201.3030   * tinq      &
                        + 1533.6164  * tinq**2   &
                        - 6491.0083  * tinq**3   &
                        + 15921.9239 * tinq**4   &
                        - 22691.1632 * tinq**5   &
                        + 17334.7529 * tinq**6   &
                        - 5517.9360  * tinq**7)  &
                        * exp(cIe21 / (kboltz * temp))
    elseif (temp  <  1d3) then
       cltab(32, itemp) = 6.0d-12 * exp(15.8996 &
                        - 201.3030   * tinq      &
                        + 1533.6164  * tinq**2   &
                        - 6491.0083  * tinq**3   &
                        + 15921.9239 * tinq**4   &
                        - 22691.1632 * tinq**5   &
                        + 17334.7529 * tinq**6   &
                        - 5517.9360  * tinq**7)  &
                        * exp(cIe21 / (kboltz * temp))
    else
       cltab(32,itemp) = 4.95d-11 * temp**0.35
    endif
!
! H2 -- ortho and para states must be accounted for separately
!    -- rates from S91 using the fit from WBV96
!
! 1 -> 0
    f  = (8.7d-11 - 6.6d-11 * exp(-temp / 218.3d0) &
          + 6.6d-11 * exp(-2d0 * temp / 218.3d0) ) * fortho
!
    hh = (7.9d-11 - 8.7d-11 * exp(-temp / 126.4d0) &
          + 1.3d-10 * exp(-2d0 * temp / 126.4d0) ) * fpara
!
    cltab(33,itemp) = f + hh
! 2 -> 0
    f  = (1.2d-10 - 6.1d-11 * exp(-temp / 387.3d0)) * fortho
!
    hh = (1.1d-10 - 8.6d-11 * exp(-temp / 223.0d0) &
          + 8.7d-11 * exp(-2d0 * temp / 223.0d0) ) * fpara
!
    cltab(34,itemp) = f + hh
! 2 -> 1
    f  = (2.9d-10 - 1.9d-10 * exp(-temp / 348.9d0)) * fortho
!
    hh = (2.7d-10 - 2.6d-10 * exp(-temp / 250.7d0) &
         + 1.8d-10 * exp(-2d0 * temp / 250.7d0)) * fpara
!
    cltab(35,itemp) = f + hh
!
! Electrons -- from JBK87. Note that the 'electron' rate given in HM89 is
!              actually the proton rate, as footnote f of table 8 makes
!              clear. Note that fits are only valid for T < 10^4 K --
!              at higher temperatures, we assume that the collision
!              strength is constant (and so the rate scales as T^-1/2)
! 1 -> 0
    if (temp  <  1d3) then
       f = -9.25141d0 - 0.773782d0 * tloge + 0.361184d0 * tloge**2 &
           - 0.150892d-1 * tloge**3 - 0.656325d-3 * tloge**4
    elseif (temp  <=  1d4) then
       f = 4.44600d2 - 2.27913d2 * tloge + 42.5952 * tloge**2 &
           - 3.47620d0 * tloge**3 + 0.105085d0 * tloge**4
    else
       f = -0.990634367510893d0
    endif
!
    cltab(36,itemp) = (1d0 / 3d0) * 8.629d-6 * exp(f) / sqrt(temp)
!
! 2 -> 0
!
    if (temp  <  1d3) then
       f = -7.69735d0 - 1.30743d0 * tloge + 0.697638d0 * tloge**2 &
           - 0.111338 * tloge**3 + 0.705277d-2 * tloge**4
    elseif (temp  <=  1d4) then
       f = 3.50609d2 - 1.87474d2 * tloge + 3.61803d1 * tloge**2 &
           - 3.03283d0 * tloge**3 + 0.938138d-1 * tloge**4
    else
       f = -1.40040241654697d0
    endif
!
    cltab(37,itemp) = 0.2d0 * 8.629d-6 * exp(f) / sqrt(temp)
!
! 2 -> 1
    if (temp  <  1d3) then
       f = -7.4387d0 - 0.57443d0 * tloge + 0.358264d0 * &
           tloge**2  - 0.418166d-1 * tloge**3 + 0.235272d-2 * &
           tloge**4
    elseif (temp  <=  1d4) then
       f = 3.86186d2 - 2.02192d2 * tloge + 3.85049d1 * tloge**2 &
           - 3.19268d0 * tloge**3 + 0.978573d-1 * tloge**4
    else
       f = 0.0198312027880547d0
    endif
!
    cltab(38,itemp) = 0.2d0 * 8.629d-6 * exp(f) / sqrt(temp)
!
! Protons -- from RLB90, fits by WBV96 (T < 5000K), SCOG (T > 5000K).
! Assume rate becomes constant above limit of tabulated data.
!
! 1 -> 0
    if (temp  <  5d3) then
       f = (9.6d-11 - 1.8d-14 * temp + 1.9d-18 * temp**2) * temp**0.45
    elseif (temp  <  2d4) then
       f = 8.873d-10 * temp**0.117
    else
       f = 2.8268d-9
    endif
!
    cltab(39,itemp) = f
!
! 2 -> 0
    if (temp  <  5d3) then
       f = (3.1d-12 - 6.0d-16 * temp + 3.9d-20 * temp**2) * temp
    elseif (temp  <  2d4) then
       f = 2.314d-9 * temp**0.0965
    else
       f = 6.0175e-09
    endif
!
    cltab(40,itemp) = f
!
! 2 -> 1
    if (temp  <  5d3) then
       f = (1.0d-10 - 2.2d-14 * temp + 1.7d-18 * temp**2) * temp**0.70
    elseif (temp  <  2d4) then
       f = 9.198d-9 * temp**0.0535
    else
       f = 1.5624e-08
    endif
!
    cltab(41,itemp) = f
!
! (cl[42-43]) -- SiI fine structure lines
!
! Collisional rates -- HI, protons from HM89
!
! HI:
! 1 -> 0
    cltab(42,itemp) = 3.5d-10 * temp2**(-0.03)
! 2 -> 0
    cltab(43,itemp) = 1.7d-11 * temp2**0.17
! 2 -> 1
    cltab(44,itemp) = 5.0d-10 * temp2**0.17
!
! Proton rates are independent of T, so we don't bother to tabulate
! them here -- instead, they're listed in cool_func.F
!
! No data for H2, electrons.
!
! Fine structure cooling -- two-level atoms
!
! We use the convention that the ground state is level 0 and
! the excited state is level 1.
!
! The cooling rate is simply:
!
! Lambda = A(1->0) * E(1->0) * f(1,LTE) * C(1->0) / ( A(1->0) +
!          C(1->0) * (1 + f(1,LTE)) )
!
! where f(1,LTE) is the level population of level 1 in LTE.
! We tabulate f(1,LTE) and the temperature-dependent bits of C(1->0),
! and compute everything else in the cooling function.
!
! (cl[45-48]) -- CII fine structure lines
!
! Collisional rates:
!
! (cl45) CII - HI (HM89 below 2000K; K86 above 2000K)
!
! (Note that the high T coefficient has been tweaked slightly to ensure that
!  the rate is continuous at 2000K -- the adjustment is well within the
!  errors).
!
    if (temp  <=  2d3) then
       cltab(45,itemp) = 8d-10 * temp2**0.07d0
    else
       cltab(45,itemp) = 3.113619d-10 * temp2**0.385d0
    endif
    !
    ! (cl46) CII - H2 (Below 250K, we use the fit from WBV96 to the data from
    !                  FL77. Above 250K, we assume that the rate scales the
    !                  same as the low-temp. HI rate)
    !
    if (temp  <  250d0) then
       f  = (4.7d-10 + 4.6d-13 * temp) * fortho
       hh = (2.5d-10 * temp**0.12) * fpara
    else
       f  = (5.85d-10 * temp**0.07) * fortho
       hh = (4.85d-10 * temp**0.07) * fpara
    endif
!
    cltab(46,itemp) = f + hh
!
! (cl47) CII - electron (WB02).
!
! (Note that the high T coefficient has been tweaked slightly to ensure that
!  the rate is continuous at 2000K -- the adjustment is well within the
!  errors).
!
    if (temp  <=  2d3) then
       cltab(47,itemp) = 3.86d-7 / sqrt(temp2)
    else
       cltab(47,itemp) = 2.426206d-7 / temp2**0.345d0
    endif
!
! Proton rate negligible below 10^4K, ignorable below 10^5K.
!
! Finally, cl48 holds f(1,LTE)
!
    cltab(48,itemp) = 2d0 * exp(-91.25d0 / temp)
!
! (cl[49-51]) -- SiII fine structure lines
!
! Collisional rates:
!
! (cl49) SiII - HI (fit by SCOG to data from R90).
!
    cltab(49,itemp) = 4.95d-10 * temp2**0.24
!
! (cl50) SiII - electron (DK91 -- extrapolated to T < 4000K,
!                         assuming constant collision strength)
!
    cltab(50,itemp) = 1.2d-6 / sqrt(temp2)
!
! cl51 holds f(1,LTE)
!
    cltab(51, itemp) = 2d0 * exp(-412.24d0 / temp)
!
! cl52 -- HI excitation cooling (aka Lyman-alpha cooling) - for
!         iflag_atom = 1, this is included in the atomic cooling
!         term computed below
!
! cl53, cl54 -- these values are used in the calculation of the
! cooling rate due to electron recombination with PAHs
!
    cltab(53, itemp) = 0.74d0 / temp**0.068
    cltab(54, itemp) = 4.65d-30 * phi_pah * temp**0.94d0
!
 enddo over_temp
!
!  Rates requiring spline fits to tabulated data:
!
! H2 rovibrational cooling -- data in mol_data.h
!
 do i = 1, nh2data
    h2_temp(i) = 1d1**(2d0 + 5d-2 * (I - 1))
 enddo
!
 call spline_eval(nh2data, h2_temp, h2_lte, nmd, temptab, rate1)
 call spline_eval(nh2data, h2_temp, h2_h_rate, nmd, temptab, rate2)
 call spline_eval(nh2data, h2_temp, h2_h2_rate, nmd, temptab,rate3)
!
 do itemp = 1, nmd
    if (temptab(itemp)  >=  1d2 .and. temptab(itemp)  <=  1d4) then
       cltab(1,itemp) = 1d1**rate1(itemp)
       cltab(2,itemp) = 1d1**rate2(itemp)
       cltab(3,itemp) = 1d1**rate3(itemp)
    endif
 enddo
!
! Initialize tables for CO, H2O cooling rates
!
! CO data: 1K bins, 5K -> 2000K
!          0.1dex bins in column density
!
 do i = 1, nTco
    co_temptab(i) = 4d0 + 1d0 * i
 enddo
!
 do i = 1, ncdco
    co_colntab(i) = 14.4d0 + 0.1d0 * i
 enddo
!
! CO rotational cooling
!
! For each quantity:
!   Do spline fits over T range for each N value
!
 call spline_eval(nco_temp, co_temp, co_data_L0, nTco, co_temptab, co_L0)
!
 do i = 1, nco_column
    do j = 1, nco_temp
       co_lte_raw(j) = co_data_LTE(nco_temp * (i-1) + j)
       co_n05_raw(j) = co_data_n05(nco_temp * (i-1) + j)
       co_alp_raw(j) = co_data_alp(nco_temp * (i-1) + j)
    enddo
!
    call spline_eval(nco_temp, co_temp, co_lte_raw, nTco, co_temptab, co_lte_fit)
    call spline_eval(nco_temp, co_temp, co_n05_raw, nTco, co_temptab, co_n05_fit)
    call spline_eval(nco_temp, co_temp, co_alp_raw, nTco, co_temptab, co_alp_fit)
!
    do j = 1, nTco
       co_lte_smalltab(i,j) = co_lte_fit(j)
       co_n05_smalltab(i,j) = co_n05_fit(j)
       co_alp_smalltab(i,j) = co_alp_fit(j)
    enddo
 enddo
!
!   do spline fits over N to fill in table
!
 do j = 1, nTco
    do i = 1, nco_column
       co_lte_fxT(i) = co_lte_smalltab(i,j)
       co_n05_fxT(i) = co_n05_smalltab(i,j)
       co_alp_fxT(i) = co_alp_smalltab(i,j)
    enddo
!
    call spline_eval(nco_column, co_column, co_lte_fxT, ncdco, co_colntab, co_lte_fit2)
    call spline_eval(nco_column, co_column, co_n05_fxT, ncdco, co_colntab, co_n05_fit2)
    call spline_eval(nco_column, co_column, co_alp_fxT, ncdco, co_colntab, co_alp_fit2)
!
    do i = 1, ncdco
       co_lte(i,j) = co_lte_fit2(i)
       co_n05(i,j) = co_n05_fit2(i)
       co_alp(i,j) = co_alp_fit2(i)
    enddo
 enddo
!
 do j = 1, nTco-1
    dTco_L0(j)  = co_L0(j+1)  - co_L0(j)
    do i = 1, ncdco
       dTco_lte(i,j) = co_lte(i,j+1) - co_lte(i,j)
       dTco_n05(i,j) = co_n05(i,j+1) - co_n05(i,j)
       dTco_alp(i,j) = co_alp(i,j+1) - co_alp(i,j)
    enddo
 enddo
!
 dTco_L0(nTco) = dTco_L0(nTco-1)
 do i = 1, ncdco
    dTco_lte(i,nTco) = dTco_lte(i,nTco-1)
    dTco_n05(i,nTco) = dTco_n05(i,nTco-1)
    dTco_alp(i,nTco) = dTco_alp(i,nTco-1)
 enddo
!
! (cl6) --  the atomic cooling function
!
! Calculate the temperatures corresponding to the data in coolatom
!
 if (iflag_atom == 1) then
    do i = 1, natom
       coolatom_temp(i) = 10**(8d0 * (i-1) / 79d0)
    enddo
!
! Compute and evaluate a spline fit to the data in coolatom
!
    call spline_eval(natom, coolatom_temp, coolatom, nmd, temptab, rate0)

    do itemp = 1, nmd
       temp   = temptab(itemp)
       atomic = 1d1**(rate0(itemp))
!
! For temp > 10^4 K, include thermal bremsstrahlung component
!
       if (temp  >  1d4)  then
          brem = 1.42d-27 * sqrt(temp)
       else
          brem = 0d0
       endif

       cltab(6,itemp)  = atomic + brem
!
! HI excitation cooling is included in this rate, so we set its entry to zero
       cltab(52,itemp) = 0d0
    enddo
!
 elseif (iflag_atom == 2) then
!
    ca2_temp(1) = tmin
    do i = 2, natom2
       ca2_temp(i) = 10**(4.25d0 + 0.05d0 * i)
    enddo
!
    call spline_eval(natom2, ca2_temp, ca2, nmd, temptab, rate0)
!
    do itemp = 1, nmd
       temp   = temptab(itemp)
       if (temp  <  10**4.3d0) then
          cltab(6,itemp) = 0d0
       else
          cltab(6,itemp) = 1d1**(rate0(itemp))
       endif
!
! Tabulate HI excitation cooling rate separately.
! This expression is from C92, based on B81.
!
       if (temp  <  1d3) then
          cltab(52,itemp) = 0d0
       else
          cltab(52,itemp) = 7.5d-19 * exp(-1.18348d5 / temp) / (1d0 + sqrt(temp / 1d5))
       endif
!
    enddo
 endif
!
! Approximate temperature derivatives
!
 do itemp = 1,nmd-1
    dtemp           = 1d0 / ( temptab(itemp + 1) - temptab(itemp) )
!
    do j = 1, ncltab
       dtcltab(j,itemp) = ( cltab(j, itemp+1) - cltab(j, itemp  ) ) * dtemp
    enddo
 enddo
!
! Manually add on the final tabulated value at itemp = ntabtemp
!
 do j = 1, ncltab
    dtcltab(j,nmd) = dtcltab(j,nmd-1)
 enddo
!
 return
end subroutine init_h2cooling
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////        I N I T _ H 2 C O O L I N G        \\\\\\\\\\
!
!=======================================================================

!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////                C O _ C O O L              \\\\\\\\\\
!
!=======================================================================
!
subroutine co_cool(temp, N_co_eff, co_rot_L0, co_rot_lte, co_rot_alpha, co_rot_n05)
!
!     written by: Simon Glover, AMNH, 2004-2005, AIP, 2006
!
!  PURPOSE: Compute CO cooling rates based on tabulated Neufeld &
!           Kaufman (1993) and Neufeld, Lepp & Melnick (1995)
!           cooling functions [see cool_func.F]
!
!  INPUT VARIABLES: temp, N_co_eff
!
!  OUTPUT VARIABLES: co_rot_L0, co_rot_lte, co_rot_alpha, co_rot_n05,
!
!--------------------------------------------------------------------
!
 use mol_data
 use splineutils, only:spline_eval
 real, intent(in)  :: temp, N_co_eff
 real, intent(out) :: co_rot_L0, co_rot_lte, co_rot_alpha, co_rot_n05
!
 real :: dtemp_co, dN_co, co_rot_lte_1, co_rot_lte_2, co_rot_alp_1, co_rot_alp_2,&
         co_rot_n05_1, co_rot_n05_2
 integer :: itemp_co, iN_co
!
! CO rotational cooling
!
! No CO cooling below 5K, as no good data
!
 if (temp  <  co_temptab(1)) then
    co_rot_L0    = 0.0
    co_rot_lte   = 0.0
    co_rot_n05   = 1.0
    co_rot_alpha = 1.0
    return
 elseif (temp == co_temptab(1)) then
    itemp_co = 1
    dtemp_co = 0d0
 elseif (temp  >=  co_temptab(nTco)) then
    itemp_co = nTco
    dtemp_co = 0d0
 else
    itemp_co = int(temp) - 4    ! Table currently starts at 5K
    dtemp_co = temp - int(temp)
 endif
!
! For column densities that do not lie within the region covered by the
! NK93 or NLM95 data, we use the smallest or largest of the tabulated
! values, as appropriate.
!
 if (N_co_eff  <=  co_colntab(1)) then
    iN_co = 1
    dN_co = 0d0
 elseif (N_co_eff  >=  co_colntab(ncdco)) then
    iN_co = ncdco
    dN_co = 0d0
 else
    iN_co = int((10 * N_co_eff) - 144)
    dN_co = (N_co_eff - co_colntab(iN_co)) / 0.1d0
 endif
!
 co_rot_L0 = co_L0(itemp_co) + dtemp_co * dTco_L0(itemp_co)
!
 co_rot_lte_1 = co_lte(iN_co,itemp_co) + dtemp_co * dTco_lte(iN_co,itemp_co)
 co_rot_alp_1 = co_alp(iN_co,itemp_co) + dtemp_co * dTco_alp(iN_co,itemp_co)
 co_rot_n05_1 = co_n05(iN_co,itemp_co) + dtemp_co * dTco_n05(iN_co,itemp_co)
 if (iN_co == ncdco) then
    co_rot_lte   = co_rot_lte_1
    co_rot_alpha = co_rot_alp_1
    co_rot_n05   = co_rot_n05_1
 else
    co_rot_lte_2 = co_lte(iN_co+1,itemp_co) + dtemp_co * dTco_lte(iN_co+1,itemp_co)
    co_rot_alp_2 = co_alp(iN_co+1,itemp_co) + dtemp_co * dTco_alp(iN_co+1,itemp_co)
    co_rot_n05_2 = co_n05(iN_co+1,itemp_co) + dtemp_co * dTco_n05(iN_co+1,itemp_co)
!
    co_rot_lte   = co_rot_lte_1 + (co_rot_lte_2 - co_rot_lte_1) * dN_co
    co_rot_alpha = co_rot_alp_1 + (co_rot_alp_2 - co_rot_alp_1) * dN_co
    co_rot_n05   = co_rot_n05_1 + (co_rot_n05_2 - co_rot_n05_1) * dN_co
 endif
!
! Do final conversion to correct units:
!
 co_rot_L0  = 10d0**(-co_rot_L0)
 co_rot_lte = 10d0**(-co_rot_lte)
 co_rot_n05 = 10d0**(co_rot_n05)
!
 return
end subroutine co_cool
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               C O _ C O O L               \\\\\\\\\\
!
!=======================================================================


!------------------------------------------------------
!+
! Written by S. Glover, AIP, September 2006.
!
! This routine computes the creation term C and
! destruction term D appearing in the rate equation
! for the number density of atomic hydrogen, n_H,
! i.e.
!
!  dn_H / dt = C - D n_H
!
! The inputs required are the temperature and density
! of the gas, the total HI column density (which
! controls the charge state of the dust and hence
! the recombination rate on grain surfaces) and the
! fractional abundances of electrons and protons
! (which will not, in general, be equal, as elements
! other than hydrogen may also contribute electrons).
!+
!------------------------------------------------------
pure subroutine hchem(temp, yn, NH, abe, abhp, C, D, sqrttemp)
 real,         intent(in)  :: temp, yn, abe, abhp, sqrttemp
 real(kind=8), intent(in)  :: NH
 real,         intent(out) :: C, D

! Temperature in Kelvin such that kb * temp_1ev = 1eV
! real, parameter :: temp_1ev = eV / kboltz
 real, parameter :: dtemp_1ev = kboltz / eV

 real :: lnte, tinv, var1, var2, var3, var4
 real :: AV, G_dust, phi
 real :: yne, ynhp
 real :: lnte2, lnte4, lnte6

 real :: k_ci, k_rec, k_gr, tinvterm

! Compute rate coefficients
 lnte = log(temp * dtemp_1eV)
 tinv = 1d0 / temp
!
! Collisional ionization of HI by electrons
! From A97; based on data from J87
!
 lnte2 = lnte*lnte
 lnte4 = lnte2*lnte2
 lnte6 = lnte4*lnte2

 k_ci = exp(-32.71396786d0                  &
             + 13.5365560d0  * lnte          &
             - 5.73932875d0  * (lnte2)       &
             + 1.56315498d0  * (lnte2*lnte)  &
             - 0.28770560d0  * (lnte4)       &
             + 3.48255977d-2 * (lnte4*lnte)  &
             - 2.63197617d-3 * (lnte6)       &
             + 1.11954395d-4 * (lnte6*lnte)  &
             - 2.03914985d-6 * (lnte4*lnte4))
!
! Case B gas-phase recombination
! From F92.
!
 ! optimisation by DJP: (315614d0 * tinv)**1.500d0
 ! replace with term*sqrt(term): removes fractional power
 tinvterm = (315614d0 * tinv)
 k_rec = 2.753d-14 * tinvterm*sqrt(tinvterm) / &
         ((1d0 + (115188d0 * tinv)**0.407d0)**2.242d0)
 !k_rec = 2.753d-14 * (315614d0 * tinv)**1.500d0 / &
 !        ((1d0 + (115188d0 * tinv)**0.407d0)**2.242d0)

! Recombination on grain surfaces. Rate from WD01.
!
 if (abe  <  1d-20) then
! We do this to avoid numerical problems at v. small abe
    phi = 1d20
 elseif (iphoto == 0) then
    phi = uv_field_strength * sqrttemp / (yn * abe)
 else
    AV     = AV_conversion_factor * dust_to_gas_ratio * NH
    G_dust = uv_field_strength * exp(-2.5d0 * AV)
    phi    = G_dust * sqrttemp / (yn * abe)
 endif

 var1 = 5.087d2 * temp**1.586d-2
 var2 = - 0.4723d0 - 1.102d-5 * log(temp)

 if (phi == 0d0) then
    k_gr = 1.225d-13 * dust_to_gas_ratio
 else
    var3  = 8.074d-6 * phi**1.378d0
    var4  = (1d0 + var1 * phi**var2)
    k_gr  = 1.225d-13 * dust_to_gas_ratio / (1d0 + var3 * var4)
 endif
!
! Compute creation and destruction rates for hydrogen
!
 yne  = abe  * yn
 ynhp = abhp * yn
!
 C = k_rec * yne * ynhp + k_gr * ynhp * yn
!
 D = cosmic_ray_ion_rate + k_ci * yne
!
 return
end subroutine hchem
!
! REFERENCES:
!
!      J87    -- Janev et al, 1987, 'Elementary Processes in
!                Hydrogen-Helium Plasmas', Springer
!      F92    -- Ferland et al, 1992, ApJ, 387, 95
!      A97    -- Abel et al, 1997, New Astron, 2, 181
!      WD01   -- Weingartner & Draine, 2001, ApJ, 563, 842
!

end module h2cooling
