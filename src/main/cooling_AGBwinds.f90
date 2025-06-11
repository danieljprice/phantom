!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module cooling_AGBwinds
!
!  Contains routines for cooling of AGB outflows
!  Routines are originally by Simon Glover,
!  Translated to Fortran 90 and adapted
!  for use in Phantom by Daniel Price (2011)
!
! :References:
!   Sembach et al. (2000) ApJ 528, 310
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - abundc            : *Carbon abundance*
!   - abunde            : *electron abundance*
!   - abundo            : *Oxygen abundance*
!   - abundsi           : *Silicon abundance*
!   - dchem             : *distance for chemistry of HI*
!   - dlq               : *distance for column density in cooling function*
!   - dphot             : *photodissociation distance used for CO/H2*
!   - dphotflag         : *photodissociation distance static or radially adaptive (0/1)*
!   - dust_to_gas_ratio : *dust to gas ratio*
!   - h2chemistry       : *Calculate H2 chemistry*
!   - iflag_atom        : *Which atomic cooling (1:Gal ISM, 2:Z=0 gas)*
!   - iphoto            : *Photoelectric heating treatment (0=optically thin, 1=w/extinction)*
!
! :Dependencies: dim, fs_data, infile_utils, io, mol_data, part, physcon,
!   splineutils, units
!
  use physcon, only:kboltz, mass_proton_cgs
  use dim,     only:nabundances, nabn_AGB
  use dust_formation,   only:eps

 implicit none
!
! Only publicly visible entries are the
! cool_func and init_cooling_ism subroutines
!
 public :: cool_func, init_cooling_AGB
 public :: energ_cooling_AGB

! Required constants and parameters
! (this stuff formerly in cool.h)
!
! He:H ratio by number (=> ratio by mass is 4*abhe)
 real, parameter :: abhe = 0.1d0

! Number of entries in cooling table
 integer, parameter :: nmd = 10000

! Number of cooling / heating rates computed in cooling fn.
 integer, parameter, public :: nrates = 28


! Number of different quantities stored in cooling look-up table
 integer, parameter :: ncltab = 67

! These variables are initialised in init_cooling_AGB
 real :: temptab(nmd)
 real :: cltab(ncltab, nmd),dtcltab(ncltab, nmd)
 real :: dtlog, tmax, tmin

!
! Parameters and tables used for CO rotational cooling
!
! Total number of combinations of column density, temperature for which we have tabulated rates
 integer, parameter, public :: nTco = 1991

! Number of CO column densities for which we have tabuled data
 integer, parameter, public :: ncdco = 46

 real :: co_temptab(nTco), co_colntab(ncdco)
 real :: co_L0(nTco), dTco_L0(nTco)
 real :: co_lte(ncdco,nTco), co_n05(ncdco,nTco), co_alp(ncdco,nTco)
 real :: dTco_lte(ncdco,nTco), dTco_n05(ncdco,nTco), dTco_alp(ncdco,nTco)

!
! Parameters and tables used for CO vibrational cooling
!
! Total number of combinations of column density, temperature for which we have tabulated rates
 integer, parameter, public :: nTco_vib = 3901

! Number of CO column densities for which we have tabuled data
 integer, parameter, public :: ncdco_vib = 61

 real :: co_vib_temptab(nTco_vib), co_vib_colntab(ncdco_vib)
 real :: co_vib_LTE_final(ncdco_vib, nTco_vib)
 real :: dTco_vib_LTE(ncdco_vib, nTco_vib)

!
! Parameters and tables used for H2O rotational cooling
!
! Total number of combinations of column density, temperature for which we have tabulated rates
 integer, parameter, public :: nTh2o = 3991

! Number of H2O column densities for which we have tabuled data
 integer, parameter, public :: ncdh2o = 91

 real :: h2o_temptab(nTh2o), h2o_colntab(ncdh2o)
 real :: h2o_L0_ortho(nTh2o), dTh2o_L0_ortho(nTh2o), &
         h2o_L0_para(nTh2o),  dTh2o_L0_para(nTh2o)
 real :: h2o_LTE_ortho(ncdh2o,nTh2o), h2o_n05_ortho(ncdh2o,nTh2o), h2o_alp_ortho(ncdh2o,nTh2o), &
         h2o_LTE_para(ncdh2o,nTh2o),  h2o_n05_para(ncdh2o,nTh2o),  h2o_alp_para(ncdh2o,nTh2o)
 real :: dTh2o_LTE_ortho(ncdh2o,nTh2o), dTh2o_n05_ortho(ncdh2o,nTh2o), dTh2o_alp_ortho(ncdh2o,nTh2o), &
         dTh2o_LTE_para(ncdh2o,nTh2o),  dTh2o_n05_para(ncdh2o,nTh2o),  dTh2o_alp_para(ncdh2o,nTh2o)

!
! Parameters and tables used for H2O vibrational cooling
!
! Total number of combinations of column density, temperature for which we have tabulated rates
 integer, parameter, public :: nTh2o_vib = 3901

! Number of H2O column densities for which we have tabuled data
 integer, parameter, public :: ncdh2o_vib = 61

 real :: h2o_vib_temptab(nTh2o_vib), h2o_vib_colntab(ncdh2o_vib)
 real :: h2o_vib_LTE_final(ncdh2o_vib, nTh2o_vib)
 real :: dTh2o_vib_LTE(ncdh2o_vib, nTh2o_vib)


! These variables must be initialised during problem setup
! (in Phantom these appear in the input file when cooling is set,
!  here we give them sensible default values)
 real, public :: abund_default(nabundances) = (/0.,1.,0.,0.,0./)
!
! Strength of UV field (in Habing units)
 real, public :: uv_field_strength = 1.d0

! Dust temperature (in K)
 real :: tdust = 1.d1

! [At lower metallicities, I [SG] generally assume for
!  simplicity that this scales linearly with metallicity,
!   but it need not do so]
!
 real, public :: dust_to_gas_ratio = 0.5d-2
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
 integer, public :: iphoto = 1
!
! Flag controlling which atomic cooling function is used.
! iflag_atom = 1 is the appropriate choice for the Galactic ISM
! [iflag_atom = 2 is used for Z=0 gas]

 integer :: iflag_atom = 1

! Distance measurements needed for chemistry

 real(kind=8), public :: dlq = 3.086d19
 real(kind=8), public :: dphot0 = 1.0801d20
 real(kind=8), public :: dchem = 3.086d20
 integer, public      :: dphotflag = 0
 private

contains

!----------------------------------------------------------------
!+
!  Compute cooling term in energy equation (du/dt)
!+
!----------------------------------------------------------------
subroutine energ_cooling_AGB(T_in,rhoi,divv,gmwvar,abund,dudti,ratesq)
 use units,     only:utime,udist,unit_density,unit_ergg
 use physcon,   only:mass_proton_cgs,Rg
 use dust_formation, only:mass_per_H
 real,         intent(in)    :: T_in,rhoi,gmwvar
 real(kind=4), intent(in)    :: divv
 real,         intent(in)    :: abund(nabn_AGB)
 real,         intent(inout) :: dudti
 real, intent(out), optional :: ratesq(nrates)
 real :: dummy_ratesq(nrates)
 real :: ylamq,divv_cgs,np1,T
!
! Compute temperature and number density
!
 T = max(T_in, 20.d0)
 np1     = rhoi*unit_density/mass_per_H  ! n = (5/7)*(rho/mp), gamma=7/5?
!
! Call cooling function with all abundances
!
 divv_cgs = sign(max(abs(divv) / real(utime, kind=4), real(1.e-40, kind=4)) , divv) ! Ensuring that divv is different from 0
 if (present(ratesq)) then
    call cool_func(T,np1,dlq,divv_cgs,abund,ylamq,ratesq)
 else
    call cool_func(T,np1,dlq,divv_cgs,abund,ylamq,dummy_ratesq)
 endif
!
! Compute change in u from 'ylamq' above.
!

 if (dudti /= dudti) dudti = 0.d0

 dudti = (-1.d0*ylamq/(rhoi*unit_density))   ! Returns dudt in cgs

end subroutine energ_cooling_AGB


!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////            C O O L _ F U N C              \\\\\\\\\\
!
!=======================================================================
!
subroutine cool_func(temp, yn, dl, divv, abundances, ylam, rates)
 use fs_data
 use mol_data
 use dim,     only:nElements
 use io, only:fatal
 real,         intent(in)  :: temp
 real,         intent(in)  :: yn
 real(kind=8), intent(in)  :: dl
 real,         intent(in)  :: divv
 real,         intent(in)  :: abundances(nabn_AGB)
 real,         intent(out) :: ylam
 real,         intent(out) :: rates(nrates)

 integer :: itemp
 real    :: dtemp

 real    :: ynh2      , ynh       , yne     , ynhp     , ynhd   , ynhe

 real    :: abh2      , abo       , aboh    , abh2o    , abhd   , &
            abcI      , abcII     , absiI   , absiII   , abe    , &
            abhp      , abhI      , abco    , abheI    , abheII , &
            abheIII
! Indices for cooling species (same as in dust_formation):
 integer, parameter :: icoolH=1, icoolC=2, icoolO=3, icoolSi=4, icoolH2=5, icoolCO=6, &
                       icoolH2O=7, icoolOH=8, icoolHe=12

 integer, parameter :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10

 real    :: abundo, abundc, abundsi

 real    :: h2var0    , h2var1

 real    :: oxc01, oxc02, oxlam1, oxlam2, oxb10, oxb01, oxc10 , oxb02 , &
            oxb20, oxb21, oxb12,  oxc20,  oxn1,  oxn2,  oxn0  , oxr01 , &
            oxr02, oxr12, oxr10,  oxr20,  oxr21, oxa,   oxb   , oxc ,   &
            oxc21 , oxc12

 real    :: cIc01, cIc02, cIlam1, cIlam2, cIb10, cIb01, cIc10 , cIb02 , &
            cIb20, cIb21, cIb12,  cIc20,  cIn1,  cIn2,  cIn0  , cIr01 , &
            cIr02, cIr12, cIr10,  cIr20,  cIr21, cIa,   cIb  , cIc  ,   &
            cIc21 , cIc12

 real    :: siIc01, siIc02, siIlam1, siIlam2, siIb10, siIb01, siIc10 , &
            siIb02, siIb20, siIb21,  siIb12,  siIc20, siIn1  , siIn2 , &
            siIn0,  siIr01, siIr02,  siIr12,  siIr10, siIr20 , siIr21, &
            siIa ,  siIb  , siIc  ,  siIc21 , siIc12

 real    :: cIIc10,   cIIc01, cIIn1

 real    :: siIIc10, siIIc01, siIIn1

 real    :: PEvar0   , PEvar1   , PEvar2, &
            G_dust    , AV       , fdust    , NH_tot

 real    :: cmb_temp
 real, parameter :: redshift = 0.0d0

 real    :: H2_opacity_correction, N_H2_eff, N_H2_jeans, N_H2_lvg, vth

 real    :: tau, tau_fac, cie_opac

 real, parameter :: f_ortho_h2o =  0.75d0
 real, parameter :: log_fo_h2o  = -0.1249387d0
 real, parameter :: f_para_h2o  =  0.25d0
 real, parameter :: log_fp_h2o  = -0.602060d0

 real    :: maximum_CO_column_rot, maximum_CO_column_vib, &
            N_co_eff, n_c13o_eff, n_co18_eff, N_co_eff_vib, &
            n_c13o_eff_vib, n_co18_eff_vib, &
            N_h2o_eff_ortho, n_h2o18_eff_ortho, N_h2o_eff_para,  &
            n_h2o18_eff_para, N_h2o_eff_vib, n_h2o18_eff_vib, &
            co_rot_L0, co_rot_LTE, co_rot_n05, co_rot_alpha, co_vib_LTE_rate, & 
            c13o_rot_L0, c13o_rot_LTE, c13o_rot_n05, c13o_rot_alpha, c13o_vib_LTE_rate, &
            co18_rot_L0, co18_rot_LTE, co18_rot_n05, co18_rot_alpha, co18_vib_LTE_rate, &
            co_vib_L0_rate, co_vib_inv, co_rot_inv, h2o_rot_inv, &
            neff, dv 

 real    :: h2o_rot_L0_ortho,  h2o_rot_LTE_ortho, h2o_rot_n05_ortho, h2o_rot_alpha_ortho
 real    :: h2o_rot_L0_para,   h2o_rot_LTE_para,  h2o_rot_n05_para,  h2o_rot_alpha_para
 real    :: h2o18_rot_L0_ortho,  h2o18_rot_LTE_ortho, h2o18_rot_n05_ortho, h2o18_rot_alpha_ortho
 real    :: h2o18_rot_L0_para,   h2o18_rot_LTE_para,  h2o18_rot_n05_para,  h2o18_rot_alpha_para
 real    :: h2o_vib_LTE_rate,    h2o18_vib_lte_rate,  h2o_vib_L0_rate,     h2o_vib_inv

 real    :: temp_co

 real    :: cl1       , cl2     , cl3     , cl4     , cl5,  &
            cl6       , cl7     , cl8     , cl9     , cl10,  &
            cl11      , cl12    , cl13    , cl14    , cl15,  &
            cl16      , cl17    , cl18    , cl19    , cl20,  &
            cl21      , cl22    , cl23    , cl24    , cl25,  &
            cl26      , cl27    , cl28    , cl29    , cl30,  &
            cl31      , cl32    , cl33    , cl34    , cl35,  &
            cl36      , cl37    , cl38    , cl39    , cl40,  &
            cl41      , cl42    , cl43    , cl44    , cl45,  &
            cl46      , cl47    , cl48    , cl49    , cl50,  &
            cl51      , cl52    , cl53    , cl54    , cl55,  &
            cl56      , cl57    , cl58    , cl59    , cl60,  &
            cl61      , cl62    , cl63    , cl64    , cl65,  &
            cl66      , cl67

!
 real    :: dtcl1     , dtcl2   , dtcl3   , dtcl4   , dtcl5,  &
            dtcl6     , dtcl7   , dtcl8   , dtcl9   , dtcl10,  &
            dtcl11    , dtcl12  , dtcl13  , dtcl14  , dtcl15,  &
            dtcl16    , dtcl17  , dtcl18  , dtcl19  , dtcl20,  &
            dtcl21    , dtcl22  , dtcl23  , dtcl24  , dtcl25,  &
            dtcl26    , dtcl27  , dtcl28  , dtcl29  , dtcl30,  &
            dtcl31    , dtcl32  , dtcl33  , dtcl34  , dtcl35,  &
            dtcl36    , dtcl37  , dtcl38  , dtcl39  , dtcl40,  &
            dtcl41    , dtcl42  , dtcl43  , dtcl44  , dtcl45,  &
            dtcl46    , dtcl47  , dtcl48  , dtcl49  , dtcl50,  &
            dtcl51    , dtcl52  , dtcl53  , dtcl54  , dtcl55,  &
            dtcl56    , dtcl57  , dtcl58  , dtcl59  , dtcl60,  &
            dtcl61    , dtcl62  , dtcl63  , dtcl64  , dtcl65,  &
            dtcl66    , dtcl67
 
!
! ---------------------------------------------------------------------
!
! Set the atomic abundances
!
 abundo  = eps(iOx)
 abundc  = eps(iC)
 abundsi = eps(iSi)
!
! Read out tables.
!
 if (temp <= tmin) then
    itemp = 1
    dtemp = 0d0
 elseif (temp > tmax) then
    itemp = nmd
    dtemp = temp - temptab(itemp)
 else
    itemp = dlog10(temp) / dtlog + 1
    if (itemp /= itemp .or. itemp <= 0.or. itemp > nmd) then
       print*, 'fatal error in cool_func.f', itemp, temp
       stop
    endif
    dtemp = temp - temptab(itemp)
 endif
!
 dtcl1  = dtcltab(1,itemp)
 dtcl2  = dtcltab(2,itemp)
 dtcl3  = dtcltab(3,itemp)
 dtcl4  = dtcltab(4,itemp)
 dtcl5  = dtcltab(5,itemp)
 dtcl6  = dtcltab(6,itemp)
 dtcl7  = dtcltab(7,itemp)
 dtcl8  = dtcltab(8,itemp)
 dtcl9  = dtcltab(9,itemp)
 dtcl10 = dtcltab(10,itemp)
 dtcl11 = dtcltab(11,itemp)
 dtcl12 = dtcltab(12,itemp)
 dtcl13 = dtcltab(13,itemp)
 dtcl14 = dtcltab(14,itemp)
 dtcl15 = dtcltab(15,itemp)
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
 dtcl55 = dtcltab(55,itemp)
 dtcl56 = dtcltab(56,itemp)
 dtcl57 = dtcltab(57,itemp)
 dtcl58 = dtcltab(58,itemp)
 dtcl59 = dtcltab(59,itemp)
 dtcl60 = dtcltab(60,itemp)
 dtcl61 = dtcltab(61,itemp)
 dtcl62 = dtcltab(62,itemp)
 dtcl63 = dtcltab(63,itemp)
 dtcl64 = dtcltab(64,itemp)
 dtcl65 = dtcltab(65,itemp)
 dtcl66 = dtcltab(66,itemp)
 dtcl67 = dtcltab(67,itemp)
 cl1  = cltab(1,itemp) + dtemp * dtcl1
 cl2  = cltab(2,itemp) + dtemp * dtcl2
 cl3  = cltab(3,itemp) + dtemp * dtcl3
 cl4  = cltab(4,itemp) + dtemp * dtcl4
 cl5  = cltab(5,itemp) + dtemp * dtcl5
 cl6  = cltab(6,itemp) + dtemp * dtcl6
 cl7  = cltab(7,itemp) + dtemp * dtcl7
 cl8  = cltab(8,itemp) + dtemp * dtcl8
 cl9  = cltab(9,itemp) + dtemp * dtcl9
 cl10 = cltab(10,itemp) + dtemp * dtcl10
 cl11 = cltab(11,itemp) + dtemp * dtcl11
 cl12 = cltab(12,itemp) + dtemp * dtcl12
 cl13 = cltab(13,itemp) + dtemp * dtcl13
 cl14 = cltab(14,itemp) + dtemp * dtcl14
 cl15 = cltab(15,itemp) + dtemp * dtcl15
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
 cl55 = cltab(55,itemp) + dtemp * dtcl55
 cl56 = cltab(56,itemp) + dtemp * dtcl56
 cl57 = cltab(57,itemp) + dtemp * dtcl57
 cl58 = cltab(58,itemp) + dtemp * dtcl58
 cl59 = cltab(59,itemp) + dtemp * dtcl59
 cl60 = cltab(60,itemp) + dtemp * dtcl60
 cl61 = cltab(61,itemp) + dtemp * dtcl61
 cl62 = cltab(62,itemp) + dtemp * dtcl62
 cl63 = cltab(63,itemp) + dtemp * dtcl63
 cl64 = cltab(64,itemp) + dtemp * dtcl64
 cl65 = cltab(65,itemp) + dtemp * dtcl65
 cl66 = cltab(66,itemp) + dtemp * dtcl66
 cl67 = cltab(67,itemp) + dtemp * dtcl67
!
! Compute CMB temperature (used for Compton cooling and for adjusting
! the fine structure cooling rates to take account of radiative pumping
! by the CMB)
!
 CMB_temp  = 2.726d0 * (1d0 + redshift)
!
! Set abundances
!
 abh2   = abundances(icoolH2)
 abo    = abundances(icoolO)
 aboh   = abundances(icoolOH)
 abh2o  = abundances(icoolH2O)
 abco   = abundances(icoolCO)
 abcI   = abundances(icoolC)
!  abcIi  = abundances(7)
 absiI  = abundances(icoolSi)
!  absiII = abundances(9)
!  abe    = abundances(10)
!  abhp   = abundances(11)
 abhI   = abundances(icoolH)
!  abhd   = abundances(13)
 abhei  = abundances(icoolHe)
!  abheII = abundances(15)
!  abheIII = abundances(16)

 abe  = 1.0d-4  ! To make it consistent with cooling_ism
 abhp = 1.0d-7
 abheII = 0.0d0
 abheIII = 0.0d0
!
! Compute useful auxiliary variables
!
 ynh2 = abh2 * yn
 ynh  = abhI * yn
 yne  = abe  * yn
 ynhp = abhp * yn
!  ynhd = abhd * yn
 ynhe = abhei * yn

 if (temp < 50) then
   ynh2 = 0.5d0 * yn
   ynh  = 0.0d0 * yn
   ynhe = 0.104d0 * yn
 endif
!
! Compute effective column density for CO rotational cooling rate
! (See eq. 4 of Neufeld & Kaufman, 1993, ApJ, 418, 263). Tabulated values
! for the effective column density are in terms of cm^-2 / (km s^-1), but
! divv is in units of cm s^-1, so we need to do a unit conversion here.
!
! If divv is very small (or zero), then we set the effective column
! density to the largest tabulated value

 if (abco  >  1d-5 * abundo) then
! Don't bother to compute CO cooling rate when CO abundance very small
    maximum_co_column_rot = co_colntab(ncdco)
    maximum_co_column_vib = co_vib_colntab(ncdco_vib)
    if (abs(divv) < tiny(divv)) then
       N_co_eff   = maximum_co_column_rot
       n_c13o_eff = maximum_co_column_rot
       n_co18_eff = maximum_co_column_rot
       N_co_eff_vib   = maximum_co_column_vib
       n_c13o_eff_vib = maximum_co_column_vib
       n_co18_eff_vib = maximum_co_column_vib
    else
       dv = 1d-5 * abs(divv)
       N_co_eff   = dlog10(abco * yn / dv)
       ! isotopic abundance ratios are from nlm95
       n_c13o_eff = N_co_eff - 2d0
       n_co18_eff = N_co_eff - 2.69897d0  ! log10(2d-3)
       N_co_eff_vib   = N_co_eff
       n_c13o_eff_vib = n_c13o_eff
       n_co18_eff_vib = n_co18_eff
       if (n_co18_eff > maximum_co_column_rot) then
         N_co_eff   = maximum_co_column_rot
         n_c13o_eff = maximum_co_column_rot
         n_co18_eff = maximum_co_column_rot
       else if (n_c13o_eff > maximum_co_column_rot) then
         N_co_eff   = maximum_co_column_rot
         n_c13o_eff = maximum_co_column_rot
       else if (N_co_eff > maximum_co_column_rot) then
         N_co_eff   = maximum_co_column_rot
       endif

       if (n_co18_eff_vib > maximum_co_column_vib) then
         N_co_eff_vib   = maximum_co_column_vib
         n_c13o_eff_vib = maximum_co_column_vib
         n_co18_eff_vib = maximum_co_column_vib
       else if (n_c13o_eff_vib > maximum_co_column_vib) then
         N_co_eff_vib   = maximum_co_column_vib
         n_c13o_eff_vib = maximum_co_column_vib
       else if (N_co_eff_vib > maximum_co_column_vib) then
         N_co_eff_vib   = maximum_co_column_vib
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
    call co_cool(temp_co, N_co_eff, N_co_eff_vib, co_rot_L0,co_rot_LTE, &
                   co_rot_alpha, co_rot_n05, co_vib_LTE_rate)
!
! if the column density of our main isotope so small that the gas is
! optically thin, then it must also be optically thin to the other
! isotopes
    if (N_co_eff <= co_colntab(1)) then
       c13o_rot_L0  = co_rot_L0
       co18_rot_L0  = co_rot_L0
       c13o_rot_LTE = co_rot_LTE
       co18_rot_LTE = co_rot_LTE
       c13o_rot_alpha = co_rot_alpha
       co18_rot_alpha = co_rot_alpha
       c13o_rot_n05 = co_rot_n05
       co18_rot_n05 = co_rot_n05
       c13o_vib_LTE_rate = co_vib_LTE_rate
       co18_vib_LTE_rate = co_vib_LTE_rate
    else
       call co_cool(temp, n_c13o_eff, n_c13o_eff_vib, c13o_rot_L0,c13o_rot_LTE, &
                      c13o_rot_alpha, c13o_rot_n05, c13o_vib_LTE_rate)
        !
       call co_cool(temp, n_co18_eff, n_co18_eff_vib, co18_rot_L0,co18_rot_LTE, &
                      co18_rot_alpha, co18_rot_n05, co18_vib_LTE_rate)
    endif
 else
    co_rot_L0 = 0.
    co_rot_LTE = 0.
    co_rot_n05 = 0.
    co_rot_alpha = 0.
    co_vib_LTE_rate = 0.
    c13o_rot_L0 = 0.
    c13o_rot_LTE = 0.
    c13o_rot_n05 = 0.
    c13o_rot_alpha = 0.
    c13o_vib_LTE_rate = 0.
    co18_rot_L0 = 0.
    co18_rot_LTE = 0.
    co18_rot_n05 = 0.
    co18_rot_alpha = 0.
    co18_vib_LTE_rate = 0.
 endif
!
! H2O -- standard isotope, plus h2o18
!
 if (abh2o > 1d-5 * abundo) then
    if (abs(divv) < 1d0) then
       N_h2o_eff_para  = h2o_colntab(ncdh2o)
       N_h2o_eff_ortho = h2o_colntab(ncdh2o)
       n_h2o18_eff_para  = h2o_colntab(ncdh2o)
       n_h2o18_eff_ortho = h2o_colntab(ncdh2o)
       N_h2o_eff_vib   = h2o_vib_colntab(ncdh2o_vib)
       n_h2o18_eff_vib = h2o_vib_colntab(ncdh2o_vib)
    else
       dv = 1d-5 * dabs(divv)
       N_h2o_eff_para  = N_h2o_eff_vib + log_fp_h2o
       N_h2o_eff_ortho = N_h2o_eff_vib + log_fo_h2o
       N_h2o_eff_vib = dlog10(abh2o * yn / dv)
       ! isotopic abundance ratio is from nlm95
       n_h2o18_eff_para  = N_h2o_eff_para  - 2.69897d0  ! log10(2d-3)
       n_h2o18_eff_ortho = N_h2o_eff_ortho - 2.69897d0  ! log10(2d-3)
       n_h2o18_eff_vib   = N_h2o_eff_vib   - 2.69897d0  ! log10(2d-3)
       !
       if (N_h2o_eff_para > h2o_colntab(ncdh2o)) then
           N_h2o_eff_para = h2o_colntab(ncdh2o)
       endif
       if (N_h2o_eff_ortho > h2o_colntab(ncdh2o)) then
           N_h2o_eff_ortho = h2o_colntab(ncdh2o)
       endif
       if (N_h2o_eff_vib > h2o_vib_colntab(ncdh2o_vib)) then
           N_h2o_eff_vib = h2o_vib_colntab(ncdh2o_vib)
       endif
       !
       if (n_h2o18_eff_para > h2o_colntab(ncdh2o)) then
           n_h2o18_eff_para = h2o_colntab(ncdh2o)
       endif
       if (n_h2o18_eff_ortho > h2o_colntab(ncdh2o)) then
           n_h2o18_eff_ortho = h2o_colntab(ncdh2o)
       endif
       if (n_h2o18_eff_vib > h2o_vib_colntab(ncdh2o_vib)) then
           n_h2o18_eff_vib = h2o_vib_colntab(ncdh2o_vib)
       endif
    endif
   
    call h2o_rot_cool(temp, N_h2o_eff_para, h2o_rot_L0_para, h2o_rot_LTE_para, h2o_rot_alpha_para,h2o_rot_n05_para, 0)
   
    call h2o_rot_cool(temp, N_h2o_eff_ortho, h2o_rot_L0_ortho,h2o_rot_LTE_ortho, h2o_rot_alpha_ortho, &
                         h2o_rot_n05_ortho, 1)
    call h2o_rot_cool(temp, n_h2o18_eff_para, h2o18_rot_L0_para,h2o18_rot_LTE_para, h2o18_rot_alpha_para, &
                         h2o18_rot_n05_para, 0)
   
    call h2o_rot_cool(temp, n_h2o18_eff_ortho, h2o18_rot_L0_ortho,h2o18_rot_LTE_ortho, h2o18_rot_alpha_ortho, &
                         h2o18_rot_n05_ortho, 1)
   
    call h2o_vib_cool(temp, N_h2o_eff_vib, h2o_vib_LTE_rate)
   
    call h2o_vib_cool(temp, n_h2o18_eff_vib, h2o18_vib_lte_rate)
   
 else
   h2o_rot_L0_para    = 0.
   h2o_rot_LTE_para   = 0. 
   h2o_rot_alpha_para = 0.
   h2o_rot_n05_para   = 0.
   h2o_rot_L0_ortho    = 0.
   h2o_rot_LTE_ortho   = 0. 
   h2o_rot_alpha_ortho = 0.
   h2o_rot_n05_ortho   = 0.
   h2o_vib_LTE_rate   = 0.
   h2o18_rot_L0_para    = 0.
   h2o18_rot_LTE_para   = 0. 
   h2o18_rot_alpha_para = 0.
   h2o18_rot_n05_para   = 0.
   h2o18_rot_L0_ortho    = 0.
   h2o18_rot_LTE_ortho   = 0. 
   h2o18_rot_alpha_ortho = 0.
   h2o18_rot_n05_ortho   = 0.
   h2o18_vib_lte_rate   = 0.
 endif
!
! (R1) -- gas-grain cooling-heating -- dust:gas ratio already incorporated
!         into rate coefficient in init_cooling_ism
!
 rates(1) = cl14 * (temp - tdust) * yn**2
!
!
! (r2) -- H2 (vr) cooling  ars: here's where the h2 cooling is!
!
 if (ynh2 == 0d0) then
   rates(2) = 0d0
 else
!        h2var0   = ynh * cl2 + ynh2 * cl3
    h2var0   = ynh * cl4 + ynh2 * cl5 + ynhe * cl64 + ynhp * cl65 + yne * cl66
    if (abs(h2var0) < 1d-4 * cl1) then
       rates(2) = h2var0 * ynh2
    else
       rates(2) = ynh2 * cl1 / (1d0 + cl1 / h2var0)
    endif
! Ortho-para conversion heating / cooling
    rates(2) = rates(2) + 4.76d-24 * ynhp * ynh2 *(cl67 * 0.25d0 - 0.75d0)
!
! Optically thick h2 cooling - we assume that this is unimportant at
! low gas densities (i.e. below n = 10^8 cm^-3).
    !
    if (yn < 1d8) then
       H2_opacity_correction = 1d0
    else
! Estimate #1 -- h2 column within local jeans length  ars: here's where the h2 cooling is!
       vth = dsqrt(kboltz * temp / mass_proton_cgs)
       N_H2_jeans = ynh2 / (vth / dl)
! Estimate #2 -- lvg estimate
       N_H2_lvg = ynh2 / dabs(divv)
! Take the smallest value -- our local jeans length estimate ensures that
! don't use an artifically large value in the case that divv is very small
       N_H2_eff = min(N_H2_lvg, N_H2_jeans)
       call compute_h2_opacity(temp, N_H2_eff, H2_opacity_correction)
       rates(2) = rates(2) * H2_opacity_correction
    endif
 endif
!
!
! (R3)  -- atomic cooling, table-based; depending on the value of iflag_atom,
!          may include some or all of r25 - r28 (in which case the coefficients
!          for these rates will have been set to zero, to avoid double counting)
! (R25) -- hi electronic excitation cooling [aka lyman-alpha cooling]
! (R26) -- hei electronic excitation cooling
!          (n.b. note the dependence of the heII metastable rate (cl60) on
!           the heII number density and the _square_ of the electron density)
! (R27) -- heII electronic excitation cooling
! (R28) -- thermal bremsstrahlung
!
 rates(3)  = cl6 * yn**2
 rates(25) = cl52 * ynh  * yne
 rates(26) = cl59 * ynhe * yne + cl60 * abheII * yn * yne**2
 rates(27) = cl61 * abheII * yn * yne
 rates(28) = cl62 * (abhp + abheII) * yn * yne + cl63 * abheIII * yn * yne
!
! (r4)  -- cooling through coll. h20 (r) w. h, h2 and he
! (r23) -- cooling through coll. h2018 (r) w. h, h2 and he
!
 neff = ynh2 + dsqrt(2d0) * ynh + 0.5d0 * ynhe
 if (abh2o <= 1d-5 * abundo .or. neff == 0d0) then
    rates(4)  = 0d0
    rates(23) = 0d0
 else
!
! Standard h2o:
!
    h2o_rot_inv = (1d0 / h2o_rot_L0_ortho) + (neff / h2o_rot_LTE_ortho) + (1d0 / h2o_rot_L0_ortho) * &
     (1d0 - h2o_rot_n05_ortho * h2o_rot_L0_ortho / h2o_rot_LTE_ortho) * (neff / h2o_rot_n05_ortho)** &
    h2o_rot_alpha_ortho
    
    rates(4) = f_ortho_h2o * abh2o * neff * yn / h2o_rot_inv
    
    h2o_rot_inv = (1d0 / h2o_rot_L0_para) + (neff / h2o_rot_LTE_para) + (1d0 / h2o_rot_L0_para) * &
     (1d0 - h2o_rot_n05_para * h2o_rot_L0_para / h2o_rot_LTE_para) * (neff / h2o_rot_n05_para)**h2o_rot_alpha_para
    
    rates(4) = rates(4) + f_para_h2o * abh2o * neff * yn / h2o_rot_inv
!
! Isotopes:
!
    h2o_rot_inv = (1d0 / h2o18_rot_L0_ortho) + (neff / h2o18_rot_LTE_ortho) + (1d0 / h2o18_rot_L0_ortho) * &
     (1d0 - h2o18_rot_n05_ortho * h2o18_rot_L0_ortho / h2o18_rot_LTE_ortho) * (neff / h2o18_rot_n05_ortho)** &
    h2o18_rot_alpha_ortho
    
    rates(23) = 2d-3 * f_ortho_h2o * abh2o * neff * yn / h2o_rot_inv
    
    h2o_rot_inv = (1d0 / h2o18_rot_L0_para) + (neff / h2o18_rot_LTE_para) + (1d0 / h2o18_rot_L0_para) * &
     (1d0 - h2o18_rot_n05_para * h2o18_rot_L0_para / h2o18_rot_LTE_para) *(neff / h2o18_rot_n05_para)** &
    h2o18_rot_alpha_para
    
    rates(23) = rates(23) + 2d-3 * f_para_h2o * abh2o * neff * yn / h2o_rot_inv
 endif
!
! (r5) -- h2o vibrational cooling   (coll. with h, h2)
! (r6) -- h2o18 vibrational cooling (coll. with h, h2)
!
 neff = ynh2 + ynh
 if (abh2o <= 1d-5 * abundo .or. neff == 0d0) then
    rates(5) = 0d0
    rates(6) = 0d0
 else
    h2o_vib_L0_rate = cl8 * abh2 + cl9 * abhI
    h2o_vib_inv = (1d0 / h2o_vib_L0_rate) + (yn / h2o_vib_LTE_rate)
    rates(5) = abh2o * yn * neff * (1d0 / h2o_vib_inv)
    !
    h2o_vib_L0_rate = cl8 * abh2 + cl9 * abhI
    h2o_vib_inv = (1d0 / h2o_vib_L0_rate) +(yn / h2o18_vib_lte_rate)
    rates(6) = 2d-3 * abh2o * yn * neff * (1d0 / h2o_vib_inv)
 endif
!
! (r7)  -- cooling through coll. co (r) h, h2 and he
! (r21) -- cooling from c13o (r) with h, h2 and he
! (r22) -- cooling from co18 (r) with h, h2 and he
!
 neff = ynh2 + dsqrt(2d0) * ynh + 0.5d0 * ynhe
 if (abco <= 1d-5 * abundo .or. neff == 0d0) then
    rates(7)  = 0d0
    rates(21) = 0d0
    rates(22) = 0d0
 else
!
! Standard CO:
!
    co_rot_inv = (1d0 / co_rot_L0) + (neff / co_rot_LTE) + (1d0 / co_rot_L0) * (1d0 - co_rot_n05 * &
    co_rot_L0  / co_rot_LTE) * (neff / co_rot_n05)**co_rot_alpha
    rates(7) = abco * neff * yn / co_rot_inv
!
! Isotopes:
!
    co_rot_inv = (1d0 / c13o_rot_L0) + (neff / c13o_rot_LTE) + (1d0 / c13o_rot_L0) * (1d0 - c13o_rot_n05 * &
    c13o_rot_L0  / c13o_rot_LTE) * (neff / c13o_rot_n05)**c13o_rot_alpha
    rates(21) = 1d-2 * abco * neff * yn / co_rot_inv
    
    co_rot_inv = (1d0 / co18_rot_L0) + (neff / co18_rot_LTE) + (1d0 / co18_rot_L0) * (1d0 - co18_rot_n05 * &
    co18_rot_L0  / co18_rot_LTE) * (neff / co18_rot_n05)**co18_rot_alpha
    
    rates(22) = 2d-3 * abco * neff * yn / co_rot_inv
 endif
!
! (r8) -- co vibrational cooling   (coll. with h, h2)
! (r9) -- c13o vibrational cooling (coll. with h, h2)
! (r24) -- co18 vibrational cooling (coll. with h, h2)
!
 neff = ynh2 + ynh
 if (abco <= 1d-5 * abundo .or. neff == 0d0) then
    rates(8)  = 0d0
    rates(9)  = 0d0
    rates(24) = 0d0
 else
    co_vib_L0_rate = cl12 * abh2 + cl13 * abhI
    co_vib_inv = (1d0 / co_vib_L0_rate) + (yn / co_vib_LTE_rate)
    rates(8) = abco * yn * neff * (1d0 / co_vib_inv)
    
    co_vib_inv = (1d0 / co_vib_L0_rate) + (yn / c13o_vib_LTE_rate)
    rates(9) = 1d-2 * abco * yn * neff * (1d0 / co_vib_inv)
    
    co_vib_inv = (1d0 / co_vib_L0_rate) + (yn / co18_vib_LTE_rate)
    rates(24) = 2d-3 * abco * yn * neff * (1d0 / co_vib_inv)
 endif
!
! (r10) -- oh cooling
!
 rates(10) = cl15 * aboh * neff
!
!
! (r13) -- OI fine-structure cooling
!
! Total collisional rates:
!
 oxc10  = cl18 * ynh + cl21 * ynh2 + cl24 * yne + cl27 * ynhp
 oxc20  = cl19 * ynh + cl22 * ynh2 + cl25 * yne + cl28 * ynhp
 oxc21  = cl20 * ynh + cl23 * ynh2 + cl26 * yne + cl29 * ynhp
 if (abo <= 1d-5 * abundo) then
    rates(13) = 0d0
 elseif (oxc10 == 0d0 .and. oxc20 == 0d0 .and. oxc21 == 0d0) then
    rates(13) = 0d0
 else
    oxa =  dexp(-oxe10 / (kboltz * temp))
    oxb =  dexp(-oxe20 / (kboltz * temp))
    oxc =  dexp(-oxe21 / (kboltz * temp))
    oxc01  = 0.6d0 * oxc10 * oxa
    oxc02  = 0.2d0 * oxc20 * oxb
    oxc12  = (1d0 / 3d0) * oxc21 * oxc
!
! Stimulated emission and absorption:
!
    call compute_stim(oxa10, oxe10, cmb_temp, oxb10)
    oxb01 = 0.6d0 * oxb10
    call compute_stim(oxa20, oxe20, cmb_temp, oxb20)
    oxb02 = 0.2d0 * oxb20
    call compute_stim(oxa21, oxe21, cmb_temp, oxb21)
    oxb12 = (1d0 / 3d0) * oxb21
!
! Total transition rates:
!
    oxr01  = oxc01 + oxb01
    oxr02  = oxc02 + oxb02
    oxr12  = oxc12 + oxb12
    oxr10  = oxc10 + oxa10 + oxb10
    oxr20  = oxc20 + oxa20 + oxb20
    oxr21  = oxc21 + oxa21 + oxb21
   
    call three_level_pops(oxr01, oxr02, oxr12, oxr10, oxr20,oxr21, oxn0, oxn1, oxn2)
! Total emitted energy:
!
    oxlam1 = (oxa10 + oxb10) * oxe10 * oxn1 + ((oxa20 +oxb20) * oxe20 + (oxa21 + oxb21) * oxe21) * oxn2
!
! Total absorbed energy:
!
    oxlam2 = (oxb01 * oxe10 + oxb02 * oxe20) * oxn0 +oxb12 * oxe21 * oxn1
!
! Net cooling rate is emission - absorption; note that if the latter
! term is larger, we have heating
!
    rates(13) = (oxlam1 - oxlam2) * abo * yn
 endif
!
!
! (r14) --  CI fine-structure cooling
!
! Collisional rates:
!
 cic10  = cl30 * ynh + cl33 * ynh2 + cl36 * yne + cl39 * ynhp
 cic20  = cl31 * ynh + cl34 * ynh2 + cl37 * yne + cl40 * ynhp
 cic21  = cl32 * ynh + cl35 * ynh2 + cl38 * yne + cl41 * ynhp
 if (abcI <= 1d-5 * abundc) then
    rates(14) = 0d0
 elseif (cic10 == 0d0 .and. cic20 == 0d0 .and. cic21 == 0d0) then
    rates(14) = 0d0
 else
    cia =  dexp(-cie10 / (kboltz * temp))
    cib =  dexp(-cie20 / (kboltz * temp))
    cic =  dexp(-cie21 / (kboltz * temp))
    cic01  = 3d0 * cic10 * cia
    cic02  = 5d0 * cic20 * cib
    cic12  = (5d0 / 3d0) * cic21 * cic
!
! Stimulated emission and absorption:
!
    call compute_stim(cia10, cie10, cmb_temp, cib10)
    cib01 = 3d0 * cib10
    call compute_stim(cia20, cie20, cmb_temp, cib20)
    cib02 = 5d0 * cib20
    call compute_stim(cia21, cie21, cmb_temp, cib21)
    cib12 = (5d0 / 3d0) * cib21
! Total transition rates:
!
    cir01  = cic01 + cib01
    cir02  = cic02 + cib02
    cir12  = cic12 + cib12
    cir10  = cic10 + cia10 + cib10
    cir20  = cic20 + cia20 + cib20
    cir21  = cic21 + cia21 + cib21
    
    call three_level_pops(cir01, cir02, cir12, cir10, cir20,cir21, cin0, cin1, cin2)
!
! Total emitted energy:
!
    cilam1 = (cia10 + cib10) * cie10 * cin1 + ((cia20 + cib20) * cie20 + (cia21 + cib21) * cie21) * cin2
!
! Total absorbed energy:
!
    cilam2 = (cib01 * cie10 + cib02 * cie20) * cin0 + cib12 * cie21 * cin1
!
! Net cooling rate is emission - absorption; note that if the latter
! term is larger, we have heating
!
    rates(14) = (cilam1 - cilam2) * abcI * yn
 endif
!
! (r15) --  SiI fine-structure cooling
!
! Proton rates (from hm89) are constant and so there's no point
! tabulating them in coolinmo.
!
 siIc10 = cl42 * ynh + 7.2d-9 * ynhp
 siIc20 = cl43 * ynh + 7.2d-9 * ynhp
 siIc21 = cl44 * ynh + 2.2d-8 * ynhp
 if (absiI <= 1d-5 * abundsi) then
    rates(15) = 0d0
 elseif (siIc10 == 0d0 .and. siIc20 == 0d0 .and.siIc21 == 0d0) then
    rates(15) = 0d0
 else
    siIa =  dexp(-siIe10 / (kboltz * temp))
    siIb =  dexp(-siIe20 / (kboltz * temp))
    siIc =  dexp(-siIe21 / (kboltz * temp))
    
    siIc01  = 3d0 * siIc10 * siIa
    siIc02  = 5d0 * siIc20 * siIb
    siIc12  = (5d0 / 3d0) * siIc21 * siIc
!
! Stimulated emission and absorption:
!
    call compute_stim(siIa10, siIe10, cmb_temp, siIb10)
    siIb01 = 3d0 * siIb10
    call compute_stim(siIa20, siIe20, cmb_temp, siIb20)
    siIb02 = 5d0 * siIb20
    call compute_stim(siIa21, siIe21, cmb_temp, siIb21)
    siIb12 = (5d0 / 3d0) * siIb21
!
! Total transition rates:
!
    siIr01  = siIc01 + siIb01
    siIr02  = siIc02 + siIb02
    siIr12  = siIc12 + siIb12
    siIr10  = siIc10 + siIa10 + siIb10
    siIr20  = siIc20 + siIa20 + siIb20
    siIr21  = siIc21 + siIa21 + siIb21
    
    call three_level_pops(siIr01, siIr02, siIr12, siIr10,siIr20, siIr21, siIn0, siIn1, siIn2)
!
! Total emitted energy:
!
    siIlam1 = (siIa10 + siIb10) * siIe10 * siIn1 +((siIa20 + siIb20) * siIe20 + (siIa21 + siIb21) *siIe21) * siIn2
!
! Total absorbed energy:
!
    siIlam2 = (siIb01 * siIe10 + siIb02 * siIe20) * siIn0 +siIb12 * siIe21 * siIn1
!
! Net cooling rate is emission - absorption; note that if the latter
! term is larger, we have heating
!
    rates(15) = (siIlam1 - siIlam2) * absiI * yn
 endif
!
! (r20) -- H2 CIE cooling: assumes molecular fraction ~ 1, although this
!          should always be true at the densities at which cie cooling is
!          significant
!
 cie_opac = 1
 if (ynh2 > 1.e10) then
    tau = (ynh2/7e15)**2.8
    tau_fac = (1 - exp(-tau))/tau
    cie_opac = min(1.,tau_fac)
 endif
 rates(20) = cl58 * yn * (2d0 * ynh2) * cie_opac
 !
 ! Benchmarking suggests that writing this out explicitly is more efficient
 ! than using a loop (although this is probably only true if the compiler
 ! optimization is poor).
 !
 rates(1)  = 0.0d0 ! Forcing dust grain heating cooling to be 0
 rates(11) = 0.0d0
 rates(12) = 0.0d0
 rates(16) = 0.0d0
 rates(17) = 0.0d0
 rates(18) = 0.0d0
 rates(19) = 0.0d0
 ylam = rates(1)  + rates(2)  + rates(3)  + rates(4)  + rates(5)  + rates(6)  + rates(7)  + rates(8)  + rates(9)  + &
        rates(10) + rates(11) + rates(12) + rates(13) + rates(14) + rates(15) + rates(16) + rates(17) + rates(18) + &
        rates(19) + rates(20) + rates(21) + rates(22) + rates(23) + rates(24) + rates(25) + rates(26) + rates(27) + &
        rates(28)

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
!    //////////                C O _ C O O L              \\\\\\\\\\
!
!=======================================================================
!
subroutine co_cool(temp, N_co_eff, N_co_eff_vib, co_rot_L0, co_rot_LTE, co_rot_alpha, co_rot_n05, co_vib_LTE_out)
!
!     written by: Simon Glover, amnh, 2004-2005, aip, 2006
!
!  PURPOSE: compute co cooling rates based on tabulated neufeld &
!           kaufman (1993) and neufeld, lepp & melnick (1995)
!           cooling functions [see coolinmo.f]
!
!  INPUT VARIABLES: temp, N_co_eff
!
!  OUTPUT VARIABLES: co_rot_L0, co_rot_LTE, co_rot_alpha, co_rot_n05,
!                    co_vib_LTE_out
!
!--------------------------------------------------------------------
!
 use mol_data
 use splineutils, only:spline_eval
!
 real, intent(in)  :: temp, N_co_eff
 real, intent(out) :: co_rot_L0, co_rot_LTE, co_rot_alpha, co_rot_n05, co_vib_LTE_out
 integer           :: itemp_co, iN_co
!
 real   ::   co_rot_LTE_1, co_rot_LTE_2, co_rot_alp_1, co_rot_alp_2, &
             co_rot_n05_1, co_rot_n05_2, co_vib_LTE_1, co_vib_LTE_2
 real   ::   dtemp_co, dn_co, N_co_eff_vib
!
! co rotational cooling
!
 if (temp <= co_temptab(1)) then
       itemp_co = 1
       dtemp_co = 0d0
 elseif (temp >= co_temptab(ntco)) then
       itemp_co = ntco
       dtemp_co = 0d0
 else
       itemp_co = int(temp) - 9    ! xxx: table currently starts at 10k
       dtemp_co = temp - int(temp)
 endif
!
! For column densities that do not lie within the region covered by the
! NK93 or NLM95 data, we use the smallest or largest of the tabulated
! values, as appropriate.
!
 if (N_co_eff <= co_colntab(1)) then
       iN_co = 1
       dn_co = 0d0
 elseif (N_co_eff >= co_colntab(ncdco)) then
       iN_co = ncdco
       dn_co = 0d0
 else
       iN_co = int((10 * N_co_eff) - 144)
       dn_co = (N_co_eff - co_colntab(iN_co)) / 0.1d0
 endif
 
 co_rot_L0 = co_l0(itemp_co) + dtemp_co * dtco_l0(itemp_co)
 
 co_rot_LTE_1 =   co_lte(iN_co,itemp_co) + dtemp_co *dtco_lte(iN_co,itemp_co)
 co_rot_alp_1 =   co_alp(iN_co,itemp_co) + dtemp_co *dtco_alp(iN_co,itemp_co)
 co_rot_n05_1 =   co_n05(iN_co,itemp_co) + dtemp_co *dtco_n05(iN_co,itemp_co)
 if (iN_co == ncdco) then
       co_rot_LTE   = co_rot_LTE_1
       co_rot_alpha = co_rot_alp_1
       co_rot_n05   = co_rot_n05_1
 else
       co_rot_LTE_2 =   co_lte(iN_co+1,itemp_co) + dtemp_co *dtco_lte(iN_co+1,itemp_co)
       co_rot_alp_2 =   co_alp(iN_co+1,itemp_co) + dtemp_co *dtco_alp(iN_co+1,itemp_co)
       co_rot_n05_2 =   co_n05(iN_co+1,itemp_co) + dtemp_co *dtco_n05(iN_co+1,itemp_co)
       !
       co_rot_LTE = co_rot_LTE_1 + (co_rot_LTE_2 - co_rot_LTE_1) *dn_co
       co_rot_alpha = co_rot_alp_1 + (co_rot_alp_2 - co_rot_alp_1) *dn_co
       co_rot_n05 = co_rot_n05_1 + (co_rot_n05_2 - co_rot_n05_1) *dn_co
 endif
!
! Do final conversion to correct units:
!
 co_rot_L0  = 10d0**(-co_rot_L0)
 co_rot_LTE = 10d0**(-co_rot_LTE)
 co_rot_n05 = 10d0**(co_rot_n05)
!
! co vibrational cooling
!
 if (temp <= co_vib_temptab(1)) then
       itemp_co = 1
       dtemp_co = 0d0
 elseif (temp >= co_vib_temptab(ntco_vib)) then
       itemp_co = ntco_vib
       dtemp_co = 0d0
 else
       itemp_co = int(temp) - 99    ! table starts at 100k
       dtemp_co = temp - int(temp)
 endif
 
 if (N_co_eff_vib <= co_vib_colntab(1)) then
       iN_co = 1
       dn_co = 0d0
 elseif (N_co_eff_vib >= co_vib_colntab(ncdco_vib)) then
       iN_co = ncdco_vib
       dn_co = 0d0
 else
       iN_co = int((10 * N_co_eff) - 129)
       dn_co = (N_co_eff_vib - co_vib_colntab(iN_co)) / 0.1d0
 endif
 co_vib_LTE_1 =   co_vib_LTE_final(iN_co,itemp_co) + dtemp_co *dtco_vib_LTE(iN_co,itemp_co)
 if (iN_co == ncdco_vib) then
       co_vib_LTE_out = co_vib_LTE_1
 else
       co_vib_LTE_2 =   co_vib_LTE_final(iN_co+1,itemp_co) + dtemp_co *dtco_vib_LTE(iN_co+1,itemp_co)
       co_vib_LTE_out = co_vib_LTE_1 + (co_vib_LTE_2 - co_vib_LTE_1) *dn_co
 endif
 
 co_vib_LTE_out = 10d0**(-co_vib_LTE_out) * exp(-3.08d3 / temp)

 return
end subroutine co_cool
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               C O _ C O O L               \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          H 2 O _ R O T _ C O O L          \\\\\\\\\\
!
!=======================================================================
!
subroutine h2o_rot_cool(temp, N_h2o_eff, h2o_rot_l0, h2o_rot_lte, h2o_rot_alpha, h2o_rot_n05, iop_flag)
!
!     Written by: Simon Glover, amnh, 2004-2005, aip 2006
!
!  PURPOSE: compute H2O rotational cooling rates based on tabulated
!           neufeld & kaufman (1993) and neufeld, lepp & melnick (1995)
!           cooling functions [see coolinmo.f]
!
!  INPUT VARIABLES: temp, N_h2o_eff, iop_flag (0 == para, 1 == ortho)
!
!  OUTPUT VARIABLES: h2o_rot_l0, h2o_rot_lte, h2o_rot_alpha, h2o_rot_n05
!
!--------------------------------------------------------------------
!
 use mol_data
 use splineutils, only:spline_eval
 !
 real    :: temp, dtemp_h2o, dN_h2o, N_h2o_eff
 integer :: itemp_h2o, iN_h2o, iop_flag
 !
 real    :: h2o_rot_l0, h2o_rot_lte, h2o_rot_n05,h2o_rot_alpha, &
            h2o_rot_lte_1, h2o_rot_lte_2,h2o_rot_alp_1, &
            h2o_rot_alp_2, h2o_rot_n05_1,h2o_rot_n05_2
 !
 if (temp <= h2o_temptab(1)) then
       itemp_h2o = 1
       dtemp_h2o = 0d0
 elseif (temp >= h2o_temptab(nth2o)) then
       itemp_h2o = nth2o
       dtemp_h2o = 0d0
 else
       itemp_h2o = int(temp) - 9    ! xxx: table currently starts at 10k
       dtemp_h2o = temp - int(temp)
 endif
 !
 ! For column densities that do not lie within the region covered by the
 ! nk93 or nlm95 data, we use the smallest or largest of the tabulated
 ! values, as appropriate.
 !
 if (N_h2o_eff <= h2o_colntab(1)) then
    iN_h2o = 1
    dN_h2o = 0d0
 elseif (N_h2o_eff >= h2o_colntab(ncdh2o)) then
    iN_h2o = ncdh2o
    dN_h2o = 0d0
 else
    iN_h2o = int((10 * N_h2o_eff) - 144)
    dN_h2o = (N_h2o_eff - h2o_colntab(iN_h2o)) / 0.1d0
 endif
 !
 if (iop_flag == 0) then
    h2o_rot_l0 = h2o_l0_para(itemp_h2o) + dtemp_h2o * dth2o_l0_para(itemp_h2o)
    h2o_rot_lte_1 = h2o_lte_para(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_lte_para(iN_h2o,itemp_h2o)
    h2o_rot_alp_1 = h2o_alp_para(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_alp_para(iN_h2o,itemp_h2o)
    h2o_rot_n05_1 = h2o_n05_para(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_n05_para(iN_h2o,itemp_h2o)
    if (iN_h2o == ncdh2o) then
       h2o_rot_lte   = h2o_rot_lte_1
       h2o_rot_alpha = h2o_rot_alp_1
       h2o_rot_n05   = h2o_rot_n05_1
    else
       h2o_rot_lte_2 = h2o_lte_para(iN_h2o+1,itemp_h2o) +dtemp_h2o * dth2o_lte_para(iN_h2o+1,itemp_h2o)
       h2o_rot_alp_2 = h2o_alp_para(iN_h2o+1,itemp_h2o) +dtemp_h2o * dth2o_alp_para(iN_h2o+1,itemp_h2o)
       h2o_rot_n05_2 = h2o_n05_para(iN_h2o+1,itemp_h2o) +dtemp_h2o * dth2o_n05_para(iN_h2o+1,itemp_h2o)
       
       h2o_rot_lte = h2o_rot_lte_1 + (h2o_rot_lte_2 - h2o_rot_lte_1) * dN_h2o
       h2o_rot_alpha = h2o_rot_alp_1 + (h2o_rot_alp_2 - h2o_rot_alp_1) * dN_h2o
       h2o_rot_n05 = h2o_rot_n05_1 + (h2o_rot_n05_2 - h2o_rot_n05_1) * dN_h2o
    endif
 else
    h2o_rot_l0 = h2o_l0_ortho(itemp_h2o) + dtemp_h2o * dth2o_l0_ortho(itemp_h2o)
    h2o_rot_lte_1 = h2o_lte_ortho(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_lte_ortho(iN_h2o,itemp_h2o)
    h2o_rot_alp_1 = h2o_alp_ortho(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_alp_ortho(iN_h2o,itemp_h2o)
    h2o_rot_n05_1 = h2o_n05_ortho(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_n05_ortho(iN_h2o,itemp_h2o)
    if (iN_h2o == ncdh2o) then
       h2o_rot_lte   = h2o_rot_lte_1
       h2o_rot_alpha = h2o_rot_alp_1
       h2o_rot_n05   = h2o_rot_n05_1
    else
       h2o_rot_lte_2 = h2o_lte_ortho(iN_h2o+1,itemp_h2o) + dtemp_h2o * dth2o_lte_ortho(iN_h2o+1,itemp_h2o)
       h2o_rot_alp_2 = h2o_alp_ortho(iN_h2o+1,itemp_h2o) + dtemp_h2o * dth2o_alp_ortho(iN_h2o+1,itemp_h2o)
       h2o_rot_n05_2 = h2o_n05_ortho(iN_h2o+1,itemp_h2o) + dtemp_h2o * dth2o_n05_ortho(iN_h2o+1,itemp_h2o)
       !
       h2o_rot_lte = h2o_rot_lte_1 + (h2o_rot_lte_2 - h2o_rot_lte_1) * dN_h2o
       h2o_rot_alpha = h2o_rot_alp_1 + (h2o_rot_alp_2 - h2o_rot_alp_1) * dN_h2o
       h2o_rot_n05 = h2o_rot_n05_1 + (h2o_rot_n05_2 - h2o_rot_n05_1) * dN_h2o
    endif
 endif
!
! Do final conversion to correct units:
!
 h2o_rot_l0  = 10d0**(-h2o_rot_l0)
 h2o_rot_lte = 10d0**(-h2o_rot_lte)
 h2o_rot_n05 = 10d0**(h2o_rot_n05)

 return
end subroutine h2o_rot_cool
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////          H 2 O _ R O T _ C O O L          \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          H 2 O _ V I B _ C O O L          \\\\\\\\\\
!
!=======================================================================
!
subroutine h2o_vib_cool(temp, N_h2o_eff, h2o_vib_lte_out)
!
!     Written by: Simon Glover, amnh, 2004-2005, aip 2006
!
!  PURPOSE: compute H2O vibrational cooling rates based on tabulated
!           Neufeld & Kaufman (1993) and Neufeld, Lepp & Melnick (1995)
!           cooling functions [see coolinmo.f]
!
!  INPUT VARIABLES: temp, N_h2o_eff
!
!  OUTPUT VARIABLES: h2o_vib_lte_out
!
!--------------------------------------------------------------------
!
 use mol_data
 use splineutils, only:spline_eval

 real, intent(in)   :: temp, N_h2o_eff
 real, intent(out)  :: h2o_vib_lte_out
 real     :: dtemp_h2o, dN_h2o, h2o_vib_lte_1, h2o_vib_lte_2
 integer  :: itemp_h2o, iN_h2o
!
! H2O vibrational cooling
!
 if (temp <= h2o_vib_temptab(1)) then
    itemp_h2o = 1
    dtemp_h2o = 0d0
 elseif (temp >= h2o_vib_temptab(nth2o_vib)) then
    itemp_h2o = nth2o_vib
    dtemp_h2o = 0d0
 else
    itemp_h2o = int(temp) - 99    ! table starts at 100k
    dtemp_h2o = temp - int(temp)
 endif

 if (N_h2o_eff <= h2o_vib_colntab(1)) then
    iN_h2o = 1
    dN_h2o = 0d0
 elseif (N_h2o_eff >= h2o_vib_colntab(ncdh2o_vib)) then
    iN_h2o = ncdh2o_vib
    dN_h2o = 0d0
 else
    iN_h2o = int((10 * N_h2o_eff) - 129)
    dN_h2o = (N_h2o_eff - h2o_vib_colntab(iN_h2o)) / 0.1d0
 endif
 h2o_vib_lte_1 =   h2o_vib_lte_final(iN_h2o,itemp_h2o) + dtemp_h2o * dth2o_vib_lte(iN_h2o,itemp_h2o)
 if (iN_h2o == ncdh2o_vib) then
    h2o_vib_lte_out = h2o_vib_lte_1
 else
    h2o_vib_lte_2 =   h2o_vib_lte_final(iN_h2o+1,itemp_h2o) + dtemp_h2o * dth2o_vib_lte(iN_h2o+1,itemp_h2o)
    h2o_vib_lte_out = h2o_vib_lte_1 + (h2o_vib_lte_2 -h2o_vib_lte_1) * dN_h2o
 endif

 h2o_vib_lte_out = 10d0**(-h2o_vib_lte_out) * exp(-2.325d3 / temp)

 return
end subroutine h2o_vib_cool
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////          H 2 O _ V I B _ C O O L          \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////      T H R E E _ L E V E L _ P O P S      \\\\\\\\\\
!
!=======================================================================
!


subroutine three_level_pops(r01, r02, r12, r10, r20, r21,n0, n1, n2)

 implicit none

 real :: r01, r02, r12, r10, r20, r21
 real :: n0 , n1 , n2
 real :: a1 , a2 , a3 , b1 , b2 , b3
!
! If excitation rates are negligibly small, then we assume that all
! of the atoms are in level 0:
!
 if (r01 == 0d0 .and. r02 == 0d0) then
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
 
 n2 = -a1 * (a1 * b2 - b1 * a2) / ((a1 - a2) *(a1 * b3 - b1 * a3) - (a1 - a3) *(a1 * b2 - b1 * a2))
 
 n1 = (a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2
 
 n0 = 1d0 - n1 - n2
 
 return 
end subroutine three_level_pops
!=======================================================================
!
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E       //////////
!    //////////      T H R E E _ L E V E L _ P O P S     \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////          C O M P U T E _ S T I M          \\\\\\\\\\
!
!=======================================================================
!
subroutine compute_stim(a10, e10, rad_temp, b10)

 implicit none

 real :: a10, e10, rad_temp, b10
 real :: x

 x = e10 / (kboltz * rad_temp)
 if (x < 5d0) then
     b10 = a10 / (dexp(x) - 1d0)
 else
     b10 = 0d0
 endif

 return 
end subroutine compute_stim
!=======================================================================
!
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E       //////////
!    //////////          C O M P U T E _ S T I M         \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////        L O A D _ H 2 _ T A B L E          \\\\\\\\\\
!
!=======================================================================
!

subroutine load_H2_table
 use mol_data

    implicit none

    integer i, j
    open(12, file='H2-cooling-ratios.dat', status='old')
    do i = 1, nh2op
       do j = 1, nh2op
          read(12,*) h2_opac_temp(i), h2_opac_column(j), h2_opac(i,j)
       enddo
    enddo
    close (12, status='keep')

 return
end subroutine load_H2_table
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////         L O A D _ H 2 _ T A B L E         \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////    C O M P U T E _ H 2 _ O P A C I T Y    \\\\\\\\\\
!
!=======================================================================
!


subroutine compute_h2_opacity(temp, N_H2_eff, opac)
 use mol_data

    implicit none
   
    real :: temp, N_H2_eff, opac
    real :: column_min, column_max, logN, logT, diff, dN, dT
    real :: opac_tmp(2)
    integer :: ii, jj

    column_min = 10**(h2_opac_column(1))
    column_max = 10**(h2_opac_column(nh2op))
    if (N_H2_eff <= column_min) then
       opac = 1d0
       return
    elseif (N_H2_eff >= column_max) then
       opac = 0d0
       return
    else
       logN = log10(N_H2_eff)
       jj   = 1 + int(10 * (logN - 17.0))
       diff = h2_opac_column(jj+1) - h2_opac_column(jj)
       dN   = (logN - h2_opac_column(jj)) / diff
       
       logT = log10(temp)
       if (logT <= h2_opac_temp(1)) then
          ii = 1
          dT = 0d0
       elseif (logT >= h2_opac_temp(nh2op)) then
          ii = nh2op
          dT = 0d0
       else
          ii   = 1 + int((logT - 1.5) / 0.03)
          diff = h2_opac_temp(ii+1) - h2_opac_temp(ii)
          dT   = (logT - h2_opac_temp(ii)) / diff
       endif
       
       if (dT > 0d0) then
          opac_tmp(1) = h2_opac(ii,jj) + dT * (h2_opac(ii+1,jj) - h2_opac(ii,jj))
          opac_tmp(2) = h2_opac(ii,jj+1) + dT * (h2_opac(ii+1,jj+1)- h2_opac(ii,jj+1))
       else
          opac_tmp(1) = h2_opac(ii,jj)
          opac_tmp(2) = h2_opac(ii,jj+1)
       endif
       
       opac = opac_tmp(1) + dN * (opac_tmp(2) - opac_tmp(1))
       opac = 1d1**opac
   endif

 return
end subroutine compute_h2_opacity
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////    C O M P U T E _ H 2 _ O P A C I T Y    \\\\\\\\\\
!
!=======================================================================
!
!=======================================================================
!
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////              C O O L I N M O              \\\\\\\\\\
!
!=======================================================================
!
subroutine init_cooling_AGB
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
!       S78   -- Spitzer, 1978, 'Physical processes in the interstellar medium'
!       HM79  -- Hollenbach and McKee, 1979, ApJS, 41, 555
!       B81   -- Black, 1981, MNRAS, 197, 553
!       DRD83 -- Draine, Roberge & Dalgarno, 1983, ApJ, 264, 485
!       K86   -- Keenan et al, 1986, MNRAS, 220, 571
!       SK87  -- Shapiro & Kang, 1987, ApJ, 318, 32
!       JBK87 -- Johnson, Burke & Kingston, 1987, J. Phys B, 20, 2553
!       HM89  -- Hollenbach & McKee, 1989, ApJ, 342, 306
!       G90   -- Gerlich, 1990, J. Chem. Phys., 92, 2377
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
!  use mol_data,    only:nh2data,h2_h_rate,h2_h2_rate,h2_lte,h2_temp,nco_temp,nco_column,nco_data, &
!                       co_temp, co_column, co_data_L0, co_data_LTE, co_data_n05, co_data_alp
 use mol_data
 use splineutils, only:spline_eval

 integer :: i, j, itemp, idx1, idx2

 integer, parameter :: natom = 81
 real :: coolatom(natom), coolatom_temp(natom)

 integer, parameter :: natom2 = 76
 real :: ca2(natom2), ca2_temp(natom2)

 real :: rate0(nmd), rate1(nmd), rate2(nmd), rate3(nmd)
!
! Raw data -- extracted from DATA tables in mol_data.f90
!
 real :: co_lte_raw(nco_temp), co_n05_raw(nco_temp), co_alp_raw(nco_temp)
 real :: co_vib_LTE_raw(nco_vib_temp)

 real ::  h2o_LTE_raw_ortho(nh2o_temp), &
          h2o_n05_raw_ortho(nh2o_temp), &
          h2o_alp_raw_ortho(nh2o_temp)
 real ::  h2o_LTE_raw_para(nh2o_temp), &
          h2o_n05_raw_para(nh2o_temp), &
          h2o_alp_raw_para(nh2o_temp)
 real ::  h2o_vib_LTE_raw(nh2o_vib_temp)
!
! Temporary variables used in fitting procedure for CO rotational cooling
!
! CO rotational cooling:
!
 real :: co_lte_fit(nTco),   co_alp_fit(nTco),   co_n05_fit(nTco)
 real :: co_lte_fit2(ncdco), co_n05_fit2(ncdco), co_alp_fit2(ncdco)
   
 real :: co_lte_fxT(nco_column), co_n05_fxT(nco_column), co_alp_fxT(nco_column)
   
 real :: co_lte_smalltab(nco_column,nTco), co_n05_smalltab(nco_column,nTco), &
         co_alp_smalltab(nco_column,nTco)
!
! H2O rotational cooling:
!
 real :: h2o_lte_fit_ortho(nth2o), h2o_alp_fit_ortho(nth2o), h2o_n05_fit_ortho(nth2o), &
         h2o_lte_fit_para(nth2o),  h2o_alp_fit_para(nth2o),  h2o_n05_fit_para(nth2o)
 real :: h2o_lte_fit2_ortho(ncdh2o), h2o_n05_fit2_ortho(ncdh2o), h2o_alp_fit2_ortho(ncdh2o), &
         h2o_lte_fit2_para(ncdh2o),  h2o_n05_fit2_para(ncdh2o),  h2o_alp_fit2_para(ncdh2o)
 real :: h2o_lte_fxt_ortho(nh2o_column),h2o_n05_fxt_ortho(nh2o_column),h2o_alp_fxt_ortho(nh2o_column), &
         h2o_lte_fxt_para(nh2o_column), h2o_n05_fxt_para(nh2o_column), h2o_alp_fxt_para(nh2o_column)
 real :: h2o_lte_smalltab_ortho(nh2o_column,nth2o), h2o_n05_smalltab_ortho(nh2o_column,nth2o), &
         h2o_alp_smalltab_ortho(nh2o_column,nth2o), h2o_lte_smalltab_para(nh2o_column,nth2o), &
         h2o_n05_smalltab_para(nh2o_column,nth2o),  h2o_alp_smalltab_para(nh2o_column,nth2o)
!
! CO vibrational cooling:
!
 real :: co_vib_LTE_fit(ntco_vib), co_vib_LTE_fit2(ncdco_vib), co_vib_LTE_fxt(nco_vib_column)
 real :: co_vib_LTE_smalltab(nco_vib_column, ntco_vib)
!
! H2O vibrational cooling:
!
 real :: h2o_vib_lte_fit(nth2o_vib), h2o_vib_lte_fit2(ncdh2o_vib), h2o_vib_lte_fxt(nh2o_vib_column)
 real :: h2o_vib_lte_smalltab(nh2o_vib_column, nth2o_vib)
!
 real :: temp    , temp2   , f       , gg       , hh      &
       , dtemp   , tinv    , tau     , tsqrt    , opratio &
       , fortho  , brem    , fpara   , atomic   , tloge   &
       , h2e20   , h2e31   , h2n2    , h2n3     , h2q02   &
       , h2q13   , phi_pah , tinth   , tinq     , tintq   &
       , tisqt   , tfix    , tfintq  , tfinq    , tfinth  &
       , temp3   , tlog    , t4log   , blwatv   , watin   &
       , blwatr  , ctinv   , blcov   , vth      , coc     &
       , etinv   , blwatvh , t34     , ccr      , sig     &
       , frac    , gff     
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

 do itemp = 1, nmd
    temptab(itemp) = 1.D1**( (itemp - 1) * dtlog )
 enddo
!
! tabulate cooling functions for each temp
!
 over_temp: do itemp = 1, nmd
    temp  = temptab(itemp)
    temp2 = temp * 1d-2
    temp3 = temp * 1d-3
    tau   = temp * 1d-3 + 1d0
    tinv  = 1.d0 / temp
    tinth = tinv**(1d0/3d0)
    tinq  = tinv**0.25d0
    tintq = tinv**0.75d0
    tsqrt = sqrt(temp)
    tisqt = 1.d0 / tsqrt
    tlog  = dlog10(temp)
    t4log = tlog - 4d0
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
    elseif (temp  <  1d2) then
       h2n2 = 0.25d0 * (5d0 * dexp(-h2e20 / temp) /(1d0 + 5d0 * dexp(-h2e20 / temp)))
       h2n3 = 0.75d0 * ((7d0 / 3d0) * dexp(-h2e31 / temp) /(1d0 + (7d0 / 3d0) * dexp(-h2e31 / temp)))
       f    = 2.94d-11 * h2e20 * h2n2 * kboltz + &
              4.76d-10 * h2e31 * h2n3 * kboltz
       cltab(1,itemp) = max(f, tiny(f))
    elseif (temp  >  1d4) then
       cltab(1,itemp) = 1d1**(-18.253d0)
    endif
!
! (cl2) -- H rate -- 0->2 rate tweaked slightly to ensure that we match up
!                    properly with the tabulated rates at 100K
!
    h2q02 = 1d1**(-9.121d0 - (3.983d0 / tau) - (1.844d-1 / tau**2)) &
            * 5d0 * dexp(-h2e20 / temp)
!
    h2q13 = 1d1**(-9.493d0 - (4.435d0 / tau) + (2.707d-1 / tau**2)) &
            * (7d0 / 3d0) * dexp(-h2e31 / temp)
!
    if (temp  <  5d0) then
       cltab(2,itemp) = tiny(cltab)
    elseif (temp  <  1d2) then
       cltab(2,itemp) = (fpara  * h2q02 * h2e20 + fortho * h2q13 * h2e31) * kboltz
    elseif (temp  >  1d4) then
       cltab(2,itemp) = 1d1**(-21.943d0)
    endif
!
! (cl3) -- H2 rate
!
    h2q02 = 1d1**(-9.946d0 - (2.688d0 / tau) + (2.020d-1 / tau**2)) &
            * 5d0 * dexp(-h2e20 / temp)
!
    h2q13 = 1d1**(-9.785d0 - (3.024d0 / tau) + (2.930d-1 / tau**2)) &
            * (7d0 / 3d0) * dexp(-h2e31 / temp)
!
    if (temp  <  5d0) then
       cltab(3,itemp) = tiny(cltab)
    elseif (temp  <  1d2) then
       cltab(3,itemp) = (fpara  * h2q02 * h2e20  +  fortho * h2q13 * h2e31) * kboltz
    elseif (temp  >  1d4) then
       cltab(3,itemp) = 1d1**(-22.758d0)
    endif
!
!
! (cl4) -- H2-H, low density rate: based on WF07, assumes 3:1 o:p ratio
!
! Data spans range 100 < t < 6000k, but is well-behaved up to 10000k;
! at higher T, we assume rate remains constant (but note that at these
! temperatures, atomic cooling dominates)
!
    if (temp < 1d1) then
       cltab(4,itemp) = 0d0
    elseif (temp < 1d2) then
        cltab(4,itemp) = 1d1**(-16.818342d0 + 37.383713d0 * dlog10(temp3) + 58.145166d0 * dlog10(temp3)**2 + &
                         48.656103d0 * dlog10(temp3)**3 + 20.159831d0 * dlog10(temp3)**4 + 3.8479610d0 * dlog10(temp3)**5)
    elseif (temp < 1d3) then
        cltab(4,itemp) = 1d1**(-24.311209d0 + 3.5692468d0 * dlog10(temp3) - 11.332860d0 * dlog10(temp3)**2 - &
                         27.850082d0 * dlog10(temp3)**3 - 21.328264d0 * dlog10(temp3)**4 - 4.2519023d0 * dlog10(temp3)**5)
    elseif (temp < 1d4) then
        cltab(4,itemp) = 1d1**(-24.311209d0 + 4.6450521d0 * dlog10(temp3) - 3.7209846d0 * dlog10(temp3)**2 + &
                         5.9369081d0 * dlog10(temp3)**3 - 5.5108047d0 * dlog10(temp3)**4 + 1.5538288d0 * dlog10(temp3)**5)
    else
        cltab(4,itemp) = 1d1**(-24.311209d0 + 4.6450521d0- 3.7209846d0 + 5.9369081d0- 5.5108047d0 + 1.5538288d0)
    endif
!
! (cl5) -- H2-H2, low density rate: based on LPF99, assumes 3:1 o:p ratio
!
! At t > 10000k, assume rate same as at T=10000k
!
    if (temp < 1d1) then
       cltab(5,itemp) = 0d0
    elseif (temp < 1d4) then
       cltab(5,itemp) = 1d1**(-23.962112d0 + 2.09433740d0 * dlog10(temp3) - 0.77151436d0  * dlog10(temp3)**2 + &
                         0.43693353d0 * dlog10(temp3)**3 - 0.14913216d0 * dlog10(temp3)**4 - 0.033638326d0 * dlog10(temp3)**5)
    else
       cltab(5,itemp) = 1d1**(-23.962112d0 + 2.09433740d0- 0.77151436d0 + 0.43693353d0 - 0.14913216d0 - 0.033638326d0)
    endif
!
! (cl64) -- H2-He, low density rate: based on FRZ98, BFD99; assumes 3:1 o:p ratio
!
    if (temp < 1d1) then
      cltab(64,itemp) = 0d0
    elseif (temp < 1d4) then
       cltab(64,itemp) = 1d1**(-23.689237d0 + 2.1892372d0 * dlog10(temp3) - 0.81520438d0 * dlog10(temp3)**2 + &
                         0.29036281d0 * dlog10(temp3)**3 - 0.16596184d0 * dlog10(temp3)**4 + 0.19191375d0 * dlog10(temp3)**5)
    else
       cltab(64,itemp) = 1d1**(-23.689237d0 + 2.1892372d0 - 0.81520438d0 + 0.29036281d0 - 0.16596184d0 + 0.19191375d0)
    endif
!
! (cl65) -- H2-H+, low density rate: based on G90, K02; assumes 3:1 o:p ratio
!
! Fit accurate to within 5%.
!
    if (temp < 1d1) then
       cltab(65,itemp) = 0d0
    elseif (temp < 1d4) then
       cltab(65,itemp) = 1d1**(-21.716699d0 + 1.3865783d0 * dlog10(temp3) - 0.37915285d0 * dlog10(temp3)**2 + &
                         0.11453688d0 * dlog10(temp3)**3 - 0.23214154d0 * dlog10(temp3)**4 + 0.058538864d0 * dlog10(temp3)**5)
    else
       cltab(65,itemp) = 1d1**(-21.716699d0 + 1.3865783d0 - 0.37915285d0 + 0.11453688d0 - 0.23214154d0 + 0.058538864d0)
    endif
!
! (cl66) -- H2-e, low density rate: based on DRD83, assumes 3:1 o:p ratio
!
    if (temp < 1d1) then
       cltab(66,itemp) = 0d0
    elseif (temp < 2d2) then
       cltab(66,itemp) = 1d1**(-34.286155d0 - 48.537163d0 * dlog10(temp3) - 77.121176d0 * dlog10(temp3)**2 - &
                         51.352459d0 * dlog10(temp3)**3 - 15.169160d0 * dlog10(temp3)**4 - 0.98120322d0 * dlog10(temp3)**5)
    elseif (temp < 1d4) then
       cltab(66,itemp) = 1d1**(-22.190316 + 1.5728955 * dlog10(temp3) - 0.21335100 * dlog10(temp3)**2 + &
                         0.96149759 * dlog10(temp3)**3 - 0.91023195 * dlog10(temp3)**4 + 0.13749749 * dlog10(temp3)**5)
    else
       cltab(66,itemp) = 1d1**(-22.190316 + 1.5728955 - 0.21335100 + 0.96149759 - 0.91023195 + 0.13749749)
    endif
!
! (cl67): Equilibrium J=1/J=0 H2 ratio (used to compute heating/cooling
! from J = 0 <-> J=1 transitions in cool_func
!
    cltab(67,itemp) = 9d0 * exp(-170.5d0 / temp)
!
!
! (cl6) --  the atomic cooling function - this is computed by fitting
!           a cubic spline to the data specified in coolatom above, and
!           so is calculated after the main loop is done
!
! (cl[7-9]) --  h2o cooling
!
! MDS  4.2.2001: completely revised
! MDS 10.6.2001: rerevised. H2O and CO both in equilibrium
!
! For H2O vibrational cooling, we can use the low density, op. thin rate
! with little to no error for densities n < 10^10 cm^-3 -- to see this
! for H2, compare L_0, L_lte numbers in table 5 of NK93; for H, substitute
! L_0 as given by HM89
!
! (r) H2O: for H2 and H collisions
!
    watin          = 1.35d0 - 0.3d0 * dlog10(0.001 * temp)
    blwatr         = 1.32d-23 * (0.001 * temp) ** watin
!
    cltab(7,itemp) = blwatr
!
! (v) H2O: for H2 collisions only (NK93)
!
    etinv          = dexp( -2.325d3 * tinv )
    blwatv         = 1.03d-26 * temp * etinv * dexp(-47.5 * tinth)
!
    cltab(8,itemp) = blwatv
!
! (v) H2O: for H collisions only  (HM89)
!
    blwatvh        = 0.74d-26 * temp * etinv * dexp(-34.5 * tinth)
!
    cltab(9,itemp) = blwatvh
!
! (cl[10-13]) -- CO cooling (Mckee et al 1982 (APJ 259, 647))
!
! (r) CO: for H2 and H collisions
!
    t34             = (0.001d0 * temp)**0.75
    ccr             = 3.3d+06 * t34
    sig             = 3.0d-13 * t34 * tinv
    vth             = tsqrt * 1.03d+04
    coc             = 0.5d0 * kboltz * temp * sig * vth
!
! Modification by AR (after discussion with MDS) based on
! statement of errors in Mckee et al. 1982 about the accuracy
! of these values (between equations 5.4 and 5.5 on p. 655)
!
    if(itemp <= 2800) then
!
! This addition for the full run of A3, which had the segmentation
! fault about half way through, still needed an even more drastic
! change
!               7/3/01
!
       if(itemp <= 1626) then
          frac = 2.014d0 * (t4log) + 5.056d0
          frac = 10.0d0**(frac)
       
          cltab(10,itemp) = coc * frac
       
       else
          frac = 0.4d0*(t4log) + 0.7d0
          frac = 10.0d0**(frac)
       
          cltab(10,itemp) = coc * frac
        
       endif
    else
        
        cltab(10,itemp) = coc
        
    endif
!
    cltab(11,itemp) = ccr
!
! Similar to H2O, but smaller n_crit; the optically thin rate is a good
! approx. only for n < 10^8 cm^-3.
!
!
! (v) CO: for H2 collisions only (NK93)
!
    ctinv           = dexp(-3.080d3 * tinv)
    blcov           = 1.83d-26 * temp * ctinv * dexp(-68.0 * tinth)
!
    cltab(12,itemp) = blcov
!
! (v) CO: for H collisions only
!
! Analytical fit by SCOG to data of BYD02. Fit is accurate to within
! 15% for T > 700 k, to within ~50% at lower T. Note, however, that
! vibrational cooling is unlikely to be important compared to rotational
! at low T, so the larger error there is probably unimportant. If it
! does prove to be significant, it would not be difficult to come up
! with a better (albeit more complicated) fit
!
    if (temp < 7d2) then
       cltab(13,itemp) = 1.049d-35 * temp**3.98 * ctinv
    else
       cltab(13,itemp) = 1.062d-29 * temp**1.87 * ctinv
    endif
!
! NB comparison of results of CP02 for CO-He collisions with those
! of BYD02 for CO-H collisions suggests that collisions with he are
! unimportant compared to those with H; we therefore don't need to
! tabulate a rate for CO-He collisions.
!
!
! (cl14) -- gas-grain cooling (HM89, eqn 2.15)
!
    gg              = 1d0 - 0.8 * exp(-75d0 * tinv)
!
    cltab(14,itemp) = 3.8d-33 * tsqrt * gg * dust_to_gas_ratio
!
! (cl15) --  OH cooling - non-lte
! Like water cooling: HM79 universal cooling law
!
    cltab(15,itemp) = 2.84d-28 * temp * tsqrt
!
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
! Assume constant rates above 10000k to avoid high t problems.
!
 if (temp > 1d4) then
    tfix = 1d4
 else
    tfix = temp
 endif
! 1 -> 0
    f               = fortho * 2.70d-11 * (tfix**0.362)
    hh              = fpara  * 3.46d-11 * (tfix**0.316)
    cltab(21,itemp) = f + hh
! 2 -> 0
    f               = fortho * 5.49d-11 * (tfix**0.317)
    hh              = fpara  * 7.07d-11 * (tfix**0.268)
    cltab(22,itemp) = f + hh
! 2 -> 1
    f               = fortho * 2.74d-14 * (tfix**1.060)
    hh              = fpara  * 3.33d-15 * (tfix**1.360)
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
    elseif (temp  <  3686) then
       cltab(27,itemp) = 7.75d-12 * (temp**0.80)
    else
       cltab(27,itemp) = 2.65d-10 * (temp**0.37)
    endif
! 2 -> 0
    if (temp  <  511) then
       cltab(28,itemp) = 6.10d-13 * (temp**1.10)
    elseif (temp  <  7510) then
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
    if (temp < 1d3) then
       cltab(52,itemp) = 0d0
    else
       cltab(52,itemp) = 7.5d-19 * dexp(-1.18348d5 / temp) /(1d0 + dsqrt(temp / 1d5))
    endif
!
! cl53, cl54 -- these values are used in the calculation of the
! cooling rate due to electron recombination with PAHs
!
    cltab(53, itemp) = 0.74d0 / temp**0.068
    cltab(54, itemp) = 4.65d-30 * phi_pah * temp**0.94d0
!
!
! (cl58) H2 collision-induced emission cooling.
!
! from RA04. assumes that H2-H2 collisions dominate
!
    cltab(58,itemp) = 2.289d-49 * temp**4
!
! (cl59) HeI excitation cooling (from n=1)
!
! My own fit to data from BBFT00. Accurate to within 10% for
! 3.75 < log10(t) < 5.75
!
    if (temp < 2d3) then
       cltab(59,itemp) = 0d0
    else
       cltab(59,itemp) = 1.1d-19 * temp**0.082d0 *exp(-2.3d5 / temp)
    endif
!
! (cl62-63) Thermal bremsstrahlung.
!
! Rates from SK87, based on S78.
!
! cl62: singly charged ions (H+, He+)
!
    if (temp < 3.2d5) then
       gff = 0.79464 + 0.1243 * log10(temp)
    else
       gff = 2.13164 - 0.1240 * log10(temp)
    endif
!
  cltab(62,itemp) = 1.426d-27 * dsqrt(temp) * gff
!
! cl63: doubly charged ions (He++)
!
    if (temp < 1.28d6) then
       gff = 0.79464 + 0.1243 * log10(temp / 4d0)
    else
       gff = 2.13164 - 0.1240 * log10(temp / 4d0)
    endif
!
  cltab(63,itemp) = 4d0 * 1.426d-27 * dsqrt(temp) * gff
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
! CO data: 1K bins, 10K -> 2000K
!          0.1dex bins in column density
!
 do i = 1, nTco
    co_temptab(i) = 9d0 + 1d0 * i
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
 do i = 1, nco_temp
    co_temp(i) = log10(co_temp(i))
 enddo
 do i = 1, ntco
    co_temptab(i) = log10(co_temptab(i))
 enddo

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
 do i = 1, ntco
    co_temptab(i) = 1d1**(co_temptab(i))
 enddo
!
!   Do spline fits over N to fill in table
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
! CO vibrational cooling (LTE rate)
!
! 1k bins, 100K -> 4000K [Below 100K, we assume non-exponential term
!                         is constant, but exponential ==> negligible]
! 0.1dex bins in column density
!
 do I = 1, ntco_vib
    co_vib_temptab(i)  = 99d0 + 1d0 * I
 enddo
!
 do I = 1, ncdco_vib
    co_vib_colntab(i) = 12.9d0 + 0.1d0 * I
 enddo
!
 do I = 1, nco_vib_temp
    co_vib_temp(I) = log10(co_vib_temp(I))
 enddo
 do I = 1, ntco_vib
    co_vib_temptab(I) = log10(co_vib_temptab(I))
 enddo
 do i = 1, nco_vib_column
    do j = 1, nco_vib_temp
       idx1 = nco_vib_temp * (i - 1) + j
       co_vib_LTE_raw(j) = co_vib_LTE(idx1)
    enddo
   
    call spline_eval(nco_vib_temp, co_vib_temp,co_vib_LTE_raw, ntco_vib,co_vib_temptab, co_vib_LTE_fit)
   
    do j = 1, ntco_vib
       co_vib_LTE_smalltab(i,j) = co_vib_LTE_fit(j)
    enddo
 enddo
 do i = 1, ntco_vib
   co_vib_temptab(i) = 1d1**(co_vib_temptab(i))
 enddo
!
!  Do spline fits over n to fill in table
!
 do i = 1, ntco_vib
    do j = 1, nco_vib_column
       co_vib_LTE_fxt(j) = co_vib_LTE_smalltab(j,i)
    enddo
   
    call spline_eval(nco_vib_column, co_vib_column,co_vib_LTE_fxt, ncdco_vib,co_vib_colntab, co_vib_LTE_fit2)
   
    do j = 1, ncdco_vib
       co_vib_LTE_final(j,i) = co_vib_LTE_fit2(j)
    enddo
 enddo

 do i = 1, ntco_vib-1
    do j = 1, ncdco_vib
       dtco_vib_LTE(j,i) = co_vib_LTE_final(j,i+1) -co_vib_LTE_final(j,i)
    enddo
 enddo
 do i = 1, ncdco_vib
    dtco_vib_LTE(i,ntco_vib) = dtco_vib_LTE(i,ntco_vib-1)
 enddo
!
! H2O rotational cooling
!
! H2O data: 1k bins, 10k -> 4000k
!           0.1dex bins in column density
!
 do i = 1, nth2o
    h2o_temptab(i)  = 9d0 + 1d0 * i
 enddo
 !
 do i = 1, ncdh2o
    h2o_colntab(i) = 9.9d0 + 0.1d0 * i
 enddo
!
! First, combine low temperature & high temperature data
!
 do i = 1, nh2o_temp_low
    h2o_data_l0_ortho(i) = h2o_l0_low_ortho(i)
    h2o_data_l0_para(i)  = h2o_l0_low_para(i)
    h2o_temp(i)          = h2o_temp_low(i)
 enddo
 
 do i = nh2o_temp_low+1, nh2o_temp
    h2o_data_l0_ortho(i) = h2o_l0_high(i - nh2o_temp_low)
    h2o_data_l0_para(i)  = h2o_l0_high(i - nh2o_temp_low)
    h2o_temp(i)          = h2o_temp_high(i - nh2o_temp_low)
 enddo
 
 do j = 1, nh2o_column
    do i = 1, nh2o_temp
       idx1 = i + (j - 1) * nh2o_temp
       if (i <= nh2o_temp_low) then
          idx2 = i + (j - 1) * nh2o_temp_low
          h2o_data_lte_ortho(idx1) = h2o_lte_low_ortho(idx2)
          h2o_data_lte_para(idx1)  = h2o_lte_low_para(idx2)
          h2o_data_alp_ortho(idx1) = h2o_alp_low_ortho(idx2)
          h2o_data_alp_para(idx1)  = h2o_alp_low_para(idx2)
          h2o_data_n05_ortho(idx1) = h2o_n05_low_ortho(idx2)
          h2o_data_n05_para(idx1)  = h2o_n05_low_para(idx2)
       else
          idx2 = (i - nh2o_temp_low) + (j - 1) * nh2o_temp_high
          h2o_data_lte_ortho(idx1) = h2o_lte_high(idx2)
          h2o_data_lte_para(idx1)  = h2o_lte_high(idx2)
          h2o_data_alp_ortho(idx1) = h2o_alp_high(idx2)
          h2o_data_alp_para(idx1)  = h2o_alp_high(idx2)
          h2o_data_n05_ortho(idx1) = h2o_n05_high(idx2)
          h2o_data_n05_para(idx1)  = h2o_n05_high(idx2)
       endif
    enddo
 enddo
 !
 do i = 1, nh2o_temp
    h2o_temp(i) = log10(h2o_temp(i))
 enddo
 do i = 1, nth2o
    h2o_temptab(i) = log10(h2o_temptab(i))
 enddo
 call spline_eval(nh2o_temp, h2o_temp, h2o_data_l0_ortho,nth2o, h2o_temptab, h2o_l0_ortho)
 call spline_eval(nh2o_temp, h2o_temp, h2o_data_l0_para,nth2o, h2o_temptab, h2o_l0_para)
 !
 do i = 1, nh2o_column
    do j = 1, nh2o_temp
       idx1 = nh2o_temp * (i - 1) + j
       h2o_lte_raw_ortho(j) = h2o_data_lte_ortho(idx1)
       h2o_n05_raw_ortho(j) = h2o_data_n05_ortho(idx1)
       h2o_alp_raw_ortho(j) = h2o_data_alp_ortho(idx1)
       h2o_lte_raw_para(j)  = h2o_data_lte_para(idx1)
       h2o_n05_raw_para(j)  = h2o_data_n05_para(idx1)
       h2o_alp_raw_para(j)  = h2o_data_alp_para(idx1)
    enddo
    
    call spline_eval(nh2o_temp, h2o_temp, h2o_lte_raw_ortho,nth2o, h2o_temptab, h2o_lte_fit_ortho)
    call spline_eval(nh2o_temp, h2o_temp, h2o_lte_raw_para,nth2o, h2o_temptab, h2o_lte_fit_para)
    
    call spline_eval(nh2o_temp, h2o_temp, h2o_n05_raw_ortho,nth2o, h2o_temptab, h2o_n05_fit_ortho)
    call spline_eval(nh2o_temp, h2o_temp, h2o_n05_raw_para,nth2o, h2o_temptab, h2o_n05_fit_para)
    
    call spline_eval(nh2o_temp, h2o_temp, h2o_alp_raw_ortho,nth2o, h2o_temptab, h2o_alp_fit_ortho)
    call spline_eval(nh2o_temp, h2o_temp, h2o_alp_raw_para,nth2o, h2o_temptab, h2o_alp_fit_para)
    
    do j = 1, nth2o
       h2o_lte_smalltab_ortho(i,j) = h2o_lte_fit_ortho(j)
       h2o_n05_smalltab_ortho(i,j) = h2o_n05_fit_ortho(j)
       h2o_alp_smalltab_ortho(i,j) = h2o_alp_fit_ortho(j)
       h2o_lte_smalltab_para(i,j)  = h2o_lte_fit_para(j)
       h2o_n05_smalltab_para(i,j)  = h2o_n05_fit_para(j)
       h2o_alp_smalltab_para(i,j)  = h2o_alp_fit_para(j)
    enddo
 enddo
 
 do i = 1, nth2o
    h2o_temptab(i) = 1d1**(h2o_temptab(i))
 enddo 
 !
 !  Do spline fits over n to fill in table
 !
 do j = 1, nth2o
    do i = 1, nh2o_column
        h2o_lte_fxt_ortho(i) = h2o_lte_smalltab_ortho(i,j)
        h2o_n05_fxt_ortho(i) = h2o_n05_smalltab_ortho(i,j)
        h2o_alp_fxt_ortho(i) = h2o_alp_smalltab_ortho(i,j)
        h2o_lte_fxt_para(i)  = h2o_lte_smalltab_para(i,j)
        h2o_n05_fxt_para(i)  = h2o_n05_smalltab_para(i,j)
        h2o_alp_fxt_para(i)  = h2o_alp_smalltab_para(i,j)
    enddo
    !
    call spline_eval(nh2o_column, h2o_column, h2o_lte_fxt_ortho,ncdh2o, h2o_colntab, h2o_lte_fit2_ortho)
    call spline_eval(nh2o_column, h2o_column, h2o_lte_fxt_para,ncdh2o, h2o_colntab, h2o_lte_fit2_para)
    !
    call spline_eval(nh2o_column, h2o_column, h2o_n05_fxt_ortho,ncdh2o, h2o_colntab, h2o_n05_fit2_ortho)
    call spline_eval(nh2o_column, h2o_column, h2o_n05_fxt_para,ncdh2o, h2o_colntab, h2o_n05_fit2_para)
    !
    call spline_eval(nh2o_column, h2o_column, h2o_alp_fxt_ortho,ncdh2o, h2o_colntab, h2o_alp_fit2_ortho)
    call spline_eval(nh2o_column, h2o_column, h2o_alp_fxt_para,ncdh2o, h2o_colntab, h2o_alp_fit2_para)
    !
    do i = 1, ncdh2o
       h2o_lte_ortho(i,j) = h2o_lte_fit2_ortho(i)
       h2o_n05_ortho(i,j) = h2o_n05_fit2_ortho(i)
       h2o_alp_ortho(i,j) = h2o_alp_fit2_ortho(i)
       h2o_lte_para(i,j)  = h2o_lte_fit2_para(i)
       h2o_n05_para(i,j)  = h2o_n05_fit2_para(i)
       h2o_alp_para(i,j)  = h2o_alp_fit2_para(i)
    enddo
 enddo
 !
 do j = 1, nth2o-1
    dth2o_l0_ortho(j) = h2o_l0_ortho(j+1) - h2o_l0_ortho(j)
    dth2o_l0_para(j)  = h2o_l0_para(j+1)  - h2o_l0_para(j)
    do i = 1, ncdh2o
       dth2o_lte_ortho(i,j) = h2o_lte_ortho(i,j+1) - h2o_lte_ortho(i,j)
       dth2o_n05_ortho(i,j) = h2o_n05_ortho(i,j+1) - h2o_n05_ortho(i,j)
       dth2o_alp_ortho(i,j) = h2o_alp_ortho(i,j+1) - h2o_alp_ortho(i,j)
       dth2o_lte_para(i,j)  = h2o_lte_para(i,j+1) -  h2o_lte_para(i,j)
       dth2o_n05_para(i,j)  = h2o_n05_para(i,j+1) -  h2o_n05_para(i,j)
       dth2o_alp_para(i,j)  = h2o_alp_para(i,j+1) -  h2o_alp_para(i,j)
    enddo
 enddo
 !
 dth2o_l0_ortho(nth2o) = dth2o_l0_ortho(nth2o-1)
 dth2o_l0_para(nth2o)  = dth2o_l0_para(nth2o-1)
 do i = 1, ncdh2o
    dth2o_lte_ortho(i,nth2o) = dth2o_lte_ortho(i,nth2o-1)
    dth2o_n05_ortho(i,nth2o) = dth2o_n05_ortho(i,nth2o-1)
    dth2o_alp_ortho(i,nth2o) = dth2o_alp_ortho(i,nth2o-1)
    dth2o_lte_para(i,nth2o)  = dth2o_lte_para(i,nth2o-1)
    dth2o_n05_para(i,nth2o)  = dth2o_n05_para(i,nth2o-1)
    dth2o_alp_para(i,nth2o)  = dth2o_alp_para(i,nth2o-1)
 enddo
!
! H2O vibrational cooling (lte rate)
!
! 1k bins, 100K -> 4000K [Below 100K, we assume non-exponential term
!                         is constant, but exponential ==> negligible]
! 0.1dex bins in column density
!
 do i = 1, nth2o_vib
    h2o_vib_temptab(i)  = 99d0 + 1d0 * i
 enddo
 !
 do i = 1, ncdh2o_vib
    h2o_vib_colntab(i) = 12.9d0 + 0.1d0 * i
 enddo
 !
 do i = 1, nh2o_vib_temp
    h2o_vib_temp(i) = log10(h2o_vib_temp(i))
 enddo
 do i = 1, nth2o_vib
    h2o_vib_temptab(i) = log10(h2o_vib_temptab(i))
 enddo
 
 do i = 1, nh2o_vib_column
    do j = 1, nh2o_vib_temp
       idx1 = nh2o_vib_temp * (i - 1) + j
       h2o_vib_lte_raw(j) = h2o_vib_lte(idx1)
    enddo
    
    call spline_eval(nh2o_vib_temp, h2o_vib_temp,h2o_vib_lte_raw, nth2o_vib,h2o_vib_temptab, h2o_vib_lte_fit)
    
    do j = 1, nth2o_vib
       h2o_vib_lte_smalltab(i,j) = h2o_vib_lte_fit(j)
    enddo
 enddo
 do i = 1, nth2o_vib
    h2o_vib_temptab(i) = 1d1**(h2o_vib_temptab(i))
 enddo
 do i = 1, nth2o_vib
    do j = 1, nh2o_vib_column
       h2o_vib_lte_fxt(j) = h2o_vib_lte_smalltab(j,i)
    enddo
    
    call spline_eval(nh2o_vib_column, h2o_vib_column,h2o_vib_lte_fxt, ncdh2o_vib,h2o_vib_colntab, h2o_vib_lte_fit2)
    
    do j = 1, ncdh2o_vib
        h2o_vib_lte_final(j,i) = h2o_vib_lte_fit2(j)
    enddo
 enddo
 !
 do i = 1, nth2o_vib-1
    do j = 1, ncdh2o_vib
        dth2o_vib_lte(j,i) = h2o_vib_lte_final(j,i+1) - h2o_vib_lte_final(j,i)
    enddo
 enddo
 do i = 1, ncdh2o_vib
    dth2o_vib_lte(i,nth2o_vib) = dth2o_vib_lte(i,nth2o_vib-1)
 enddo
!
! (cl6) --  The atomic cooling function
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
       cltab(6,itemp)  = atomic
!
! This rate includes cooling from HI, HeI, HeII excitation, so we set 
! the corresponding individual rates to zero to avoid double counting
!
! HI excitation cooling ('Lyman-alpha' cooling):
!
          cltab(52,itemp) = 0d0
!
! HeI excitation cooling
!
          cltab(59,itemp) = 0d0
          cltab(60,itemp) = 0d0
!
! HeII excitation cooling
!
          cltab(61,itemp) = 0d0
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
! HeI, HeII excitation cooling & bremsstrahlung are included
! in the rate given in the DATA statement above, and so we 
! set the individual rates to zero to avoid double counting
!
       cltab(59,itemp) = 0d0
       cltab(60,itemp) = 0d0
       cltab(61,itemp) = 0d0
       cltab(62,itemp) = 0d0
       cltab(63,itemp) = 0d0
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
       dtcltab(j,itemp) = ( cltab(j, itemp+1) - cltab(j, itemp ) ) * dtemp
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
end subroutine init_cooling_AGB
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////        I N I T _ H 2 C O O L I N G        \\\\\\\\\\
!
!=======================================================================

end module cooling_AGBwinds
