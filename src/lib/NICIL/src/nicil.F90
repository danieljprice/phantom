!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!                                                                      !
! This is a stand-alone library that will calculate ionisation values  !
! and the coefficients for the non-ideal MHD terms: Ohmic resistivity, !
! Hall Effect and Ambipolar diffusion.                                 !
!                                                                      !
!                 Copyright (c) 2015-2023 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
!  MODULE: nicil
!
!  DESCRIPTION:
!  Contains routines to calculation the ionisation rate, and related
!  useful quantities.
!  Copyright (c) 2015-2023 James Wurster
!  References: Wurster (2016) PASA, 33:e041.
!              Wurster (2021) MNRAS, 501:5873-5891.
!  See LICENCE file for usage and distribution conditions
!
!
!  REFERENCES: Asplund et al. (2009)
!              Cox (2000)
!              Draine & Lee (1984)
!              Fujii et al. (2011)
!              Keith & Wardle (2014)
!              Gorti & Hollenbach (2008)
!              Kunz & Mouschovias (2009)
!              Lenzuni, Gail & Henning (1995)
!              Liu et al. (2003)
!              Marchand et al. (2016)
!              Nakano, Nishi & Umebayashi (2002)
!              Osterbrock (1961)
!              Pandey & Wardle (2008)
!              Pinto & Galli (2008)
!              Pollack et al. (1994)
!              Umebayashi & Nakano (1990)
!              Umebayashi & Nakano (2009)
!              Wardle (2007)
!              Wardle & Ng (1999)
!              Wang, Bai & Goodman (2019)
!              Wurster, Price & Ayliffe (2014)
!              Wurster, Bate & Price (2016)
!              Xu et al. (2019)
!              Zhao, Caselli & Li (2018)
!
!  AUTHOR: James Wurster
!
!  PRIMARY HARDCODED PARAMETERS:
!    na             -- number of grain species (using MRN distribution if > 1)
!  PRIMARY RUNTIME PARAMETERS:
!    use_ohm        -- Calcualate and use the coefficient for Ohmic resistivity
!    use_hall       -- Calcualate and use the coefficient for the Hall effect
!    use_ambi       -- Calcualate and use the coefficient for ambipolar diffusion
!    eta_constant   -- Use a constant resistivity
!    warn_approx    -- To print warnings of assumption violation to file
!    use_fdg_in     -- Import dust mass density from the parent code for each point
!    fdg            -- Grain Parameter: gas to dust mass ratio
!    a0_grain       -- Grain Parameter if na = 1 : grain radius for constant grain size
!    an_grain       -- Grain Parameter if na > 1: minimum grain radius for power-law distribution
!    ax_grain       -- Grain Parameter if na > 1: maximum grain radius for power-law distribution
!    grain_slope    -- Grain Parameter if na > 1: slope of the power-law distribution
!    rho_bulk       -- Grain Parameter: bulk grain density
!    zeta_of_rho    -- Use constant or variable cosmic ray ionisation rate
!    zeta_cgs       -- Ionisation rate (if zeta_of_rho=false)
!    zeta_cgs_min   -- Minimum allowed ionisation rate if zeta_of_rho=true
!    eta_const_type -- If eta_constant=true: determines the form of the (semi-) constant coefficients
!    C_OR           -- If eta_constant=true & eta_const_type=icnstsemi: eta_OR = C_OR
!    C_HE           -- If eta_constant=true & eta_const_type=icnstsemi (icnst): eta_HE = C_HE*B     (=C_HE)
!    C_AD           -- If eta_constant=true & eta_const_type=icnstsemi (icnst): eta_AD = C_AD*v_A^2 (=C_AD)
!    n_e_cnst       -- If eta_constant=true & eta_const_type=icnstphys: Constant electron number density
!    rho_i_cnst     -- If eta_constant=true & eta_const_type=icnstphys: Density of ionised gas
!    rho_n_cnst     -- If eta_constant=true & eta_const_type=icnstphys: Density of neutral gas
!    alpha_AD       -- If eta_constant=true & eta_const_type=icnstphys: Exponent of Power-law ion density
!    gamma_AD       -- If eta_constant=true & eta_const_type=icnstphys: Collisional coupling coefficient
!    hall_lt_zero   -- If eta_constant=true & eta_const_type=icnstphys: The sign of the Hall coefficient
!
!----------------------------------------------------------------------!
module nicil
 implicit none
 !--INPUT PAREMETERS
 !--Number of grain sizes
 integer, public, parameter :: na                =  1                ! Hardcoded for storage efficiency ( na >=0 )
 !--Turn on/off individual non-ideal MHD coefficients
 logical, public            :: use_ohm           = .true.           ! Calculate the coefficient for Ohmic resistivity
 logical, public            :: use_hall          = .true.           ! Calculate the coefficient for the Hall effect
 logical, public            :: use_ambi          = .true.           ! Calculate the coefficient for ambipolar diffusion
 !--Use a constant (false) or variable (true) cosmic ray ionisation rate
 logical, public            :: zeta_of_rho       = .false.
 !--Import the dust to gas ratio from the parent code for each point
 logical, public            :: use_fdg_in        = .false.
 !--Use constant resistivity coefficients for all three resistivity terms
 logical, public            :: eta_constant      = .false.
 !--To print warnings to file
 logical, public            :: warn_approx       = .false.
 !--To solve a reordered Jacobian for stability
 logical, private           :: reorder_Jacobian  = .false.

 !--Grain properties
 real,    public            :: fdg               = 0.01             ! gas to dust mass ratio
 real,    public            :: a0_grain          = 1.0d-5           ! grain radius for constant grain size [cm]
 real,    public            :: an_grain          = 3.0d-6           ! minimum grain radius for distribution [cm]
 real,    public            :: ax_grain          = 2.5d-5           ! maximum grain radius for grain distribution [cm]
 real,    public            :: grain_slope       = -3.5             ! slope on the grain distribution (-3.5 for MRN)
 real,    public            :: rho_bulk          = 3.0              ! bulk grain density [g/cm^3]
 real,    private,parameter :: fdg_max           = 0.99             ! the maximum dust to gas ratio permitted (< 1)

 !--Cosmic ray ionisation
 real,    public            :: zeta_cgs          =  1.20d-17        ! cosmic ray ionisation rate [s^-1] (if zeta_of_rho=.false.)
 real,    public            :: zeta_cgs_min      =  1.10d-22        ! radio-nuclide ionisation rate of 40K [s^-1] (if zeta_of_rho=.true.)
 real,    public            :: SigmaCR           = 96.0             ! Ionisation parameter: attenuation depth for cosmic rays [g cm^-2]
 real,    public            :: delta_gn          = 1.3              ! Ionisation parameter: multiplicative factor for sigmavgnbyT
 real,    public            :: pnH2              = 0.804            ! Ionisation parameter: polarizability for H2 [angstroms^3]
 real,    public            :: ProbI             = 1.0              ! Sticking probability of ions onto grains
 real,    public            :: ProbE             = 0.6              ! Sticking probability of electrons onto grains

 !--Resistivity coefficients (if fixed as constants)
 integer, public, parameter :: icnstphys         = 1                ! Index for calculating eta using physical parameters
 integer, public, parameter :: icnstsemi         = 2                ! Index for calculating eta using a semi-constant coefficient
 integer, public, parameter :: icnst             = 3                ! Index for calculating eta using a constant coefficient
 integer, public            :: eta_const_type    = icnstsemi        ! Determines the form of the (semi-) constant coefficients
 real,    public            :: C_OR              =  0.1             ! eta_OR = C_OR                                  if eta_const_type = icnst, icnstsemi
 real,    public            :: C_HE              = -0.5             ! eta_HE = C_HE (*B)                             if eta_const_type = icnst (icnstsemi)
 real,    public            :: C_AD              =  0.01            ! eta_AD = C_AD (*v_A^2)                         if eta_const_type = icnst (icnstsemi)
 real,    public            :: n_e_cnst          = 1.0d19           ! Constant electron number density  [cm^-3]      if eta_const_type = icnstphys
 real,    public            :: rho_i_cnst        = 3.8d-11          ! Density of ionised gas            [g/cm^3]     if eta_const_type = icnstphys
 real,    public            :: rho_n_cnst        = 3.8d-08          ! Density of neutral gas            [g/cm^3]     if eta_const_type = icnstphys
 real,    public            :: alpha_AD          = 0.0              ! Exponent of Power-law ion density              if eta_const_type = icnstphys
 real,    public            :: gamma_AD          = 2.6e13           ! Collisional coupling coefficient  [cm^3/(s g)] if eta_const_type = icnstphys
 logical, public            :: hall_lt_zero      = .true.           ! The sign of the Hall coefficient               if eta_const_type = icnstphys

 !--Threshholds
 real,    public            :: Texp_thresh0      = 0.005            ! Will set exp(-chi/kT) = 0 if T is too low
 real,    public            :: eta_thresh_cgs    = 1.0d25           ! threshold on eta above which a warning will be printed [cm^2/s]

 !--Additional parameters
 integer, public            :: NRctrmax          = 200              ! maximum number of Newton–Raphson iterations
 real,    public            :: NRtol             = 1.0d-8           ! default tolerance on Newton–Raphson iterations
 real,    public            :: NRtol4            = 1.0d-6           ! default tolerance on Newton–Raphson iterations if single precision
 real,    public            :: MPtol             = 1.0d-12          ! default tolerance to zero values due to machine precision
 real,    public            :: MPtol4            = 1.0d-6           ! default tolerance to zero values due to machine precision if single precision
 real,    public            :: MPval             = 1.0d-100         ! a value slightly larger than sqrt (min possible value)
 real,    public            :: MPval4            = 1.0d-35          ! a value slightly larger than sqrt (min possible value)if single precision
 real,    public            :: Cdt_diff          = 0.12             ! Coefficient to control the AD & OR timesteps
 real,    public            :: Cdt_hall          = 0.0795774        ! Coefficient to control the HE timestep (==1/4pi)
 !--END OF INPUT PARAMETERS

 !--Misc. parameters not to be modified
 real,            parameter :: warn_ratio        =  0.1             ! fraction within with two values will be assumed equal for warnings
 real,            parameter :: vrms_min_kms      =  20              ! Minimum rms velocity below which rate coefficients become unreliable [km/s] (Pinto & Galli, 2008)
 real,            parameter :: vrms_max_kms      = 500              ! Maximum rms velocity above which rate coefficients become unreliable [km/s] (Pinto & Galli, 2008)

 !--Physical Constants (CGS)
 real(kind=8),    parameter :: pi                =  3.1415926536d0  !  pi
 real(kind=8),    parameter :: twopi             =  6.2831853072d0  ! 2pi
 real(kind=8),    parameter :: fourpi            = 12.5663706144d0  ! 4pi
 real(kind=8),    parameter :: c                 =  2.997924d10     ! Speed of light [cm/s]
 real(kind=8),    parameter :: qe                =  4.8032068d-10   ! charge on electron [esu == statC == (cm^3 g)^0.5 s^-1]
 real(kind=8),    parameter :: mass_electron_cgs =  9.10938291d-28  ! Electron mass [g]
 real(kind=8),    parameter :: eV                =  1.60217657d-12  ! Electron volts [ergs]
 real(kind=8),    parameter :: mass_proton_cgs   =  1.67262158d-24  ! Proton mass [g]
 real(kind=8),    parameter :: kboltz            =  1.38066d-16     ! Boltzmann constant  [erg/K]
 real(kind=8),    parameter :: Ggrav             =  6.67408d-8      ! Gravitational constant [cm^3 g^-1 s^-2]
 real(kind=8),    parameter :: planckh           =  6.6260755d-27   ! Planck's Constant [erg s]

 !--Indicies for warnings
 integer,         parameter :: ierr_neTconv      =  1               ! fatal error code if n_electronT did not converge
 integer,         parameter :: ierr_neTle0       =  2               ! fatal error code if n_electronT < 0
 integer,         parameter :: ierr_fdg_in       =  3               ! fatal error code if use_fdg_in=.true. and fdg_in not passed in
 integer,         parameter :: ierr_nRconv       =  4               ! warning code if nion_R did not converge or yield valid values (setting n_ionR=n_eR=0)
 integer,         parameter :: ierr_highfdg      =  5               ! warning code if grain fraction is too high
 integer,         parameter :: ierr_alion        =  6               ! warning code if fully ionised
 integer,         parameter :: ierr_T            =  7               ! warning code if T==0 (from input)
 integer,         parameter :: ierr_B            =  8               ! warning code if B==0 (from input)
 integer,         parameter :: ierr_sigH         =  9               ! warning code if sigmaH == 0 to prevent errors from roundoff
 integer,         parameter :: ierr_scN          = 10               ! warning code if strong coupling approximation broke (rho_n ~ rho is false) [if warn_approx=true]
 integer,         parameter :: ierr_scI          = 11               ! warning code if strong coupling approximation broke (rho_i <<rho is false) [if warn_approx=true]
 integer,         parameter :: ierr_drift        = 12               ! warning code if drift velocity is large relative to absolute velocity [if warn_approx=true]
 integer,         parameter :: ierr_vrms         = 13               ! warning code if velocity is outside of the range valid for the rate coefficients [if warn_approx=true]
 integer,         parameter :: ierr_or           = 14               ! warning code if eta_ohm > eta_thresh (suggesting a bug in the calculation) [if warn_approx=true]
 integer,         parameter :: ierr_he           = 15               ! warning code if |eta_hall| > eta_thresh (suggesting a bug in the calculation) [if warn_approx=true]
 integer,         parameter :: ierr_ad           = 16               ! warning code if eta_ambi > eta_thresh (suggesting a bug in the calculation) [if warn_approx=true]
 integer,public,  parameter :: n_warn            = 16               ! number of types of errors including fatal & non-fatal warnings
 integer,public,  parameter :: n_fatal           =  3               ! the first n_fatal errors are fatal
 !--Indicies for ions (iterated & stored)
 integer,         parameter :: iGp               =  1               ! Charged grain
 integer,         parameter :: iGn               =  2               ! Charged grain
 integer,         parameter :: iH3p              =  1+2*na          ! Charged ion with no neutral counterpart
 integer,         parameter :: iHCOp             =  2+2*na          ! Charged ion with no neutral counterpart
 integer,         parameter :: iO2p              =  3+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iHp               =  4+2*na          ! Charged ion with temperature-dependent neutral counterpart
 integer,         parameter :: iHep              =  5+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iCp               =  6+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iOp               =  7+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iSip              =  8+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iSp               =  9+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iMgp              = 10+2*na          ! Charged ion with constant neutral counterpart
 integer,         parameter :: iH2p              = 11+2*na          ! Charged ion with constant neutral counterpart (thermal ionisation)
 integer,         parameter :: iKp               = 12+2*na          ! Charged ion with constant neutral counterpart (thermal ionisation)
 integer,         parameter :: iNap              = 13+2*na          ! Charged ion with constant neutral counterpart (thermal ionisation)
 integer,         parameter :: ie                = 14+2*na          ! Electrons from thermal ionisation
 !--Indicies for neutral species (constant)
 integer,         parameter :: iCO               = iHCOp            ! double usage for storage reasons
 integer,         parameter :: iO2               = iO2p
 integer,         parameter :: iH                = iHp
 integer,         parameter :: iHe               = iHep
 integer,         parameter :: iC                = iCp
 integer,         parameter :: iO                = iOp
 integer,         parameter :: iSi               = iSip
 integer,         parameter :: iS                = iSp
 integer,         parameter :: iMg               = iMgp
 integer,         parameter :: iH2               = iH2p
 integer,         parameter :: iK                = iKp              ! Only for thermal ionisation
 integer,         parameter :: iNa               = iNap             ! Only for thermal ionisation
 !--Important start & end indicies for the various ionisation processes
 integer,         parameter :: iirs              = iH3p             ! index of first ion relevant for cosmic ray ionisation
 integer,         parameter :: iire              = iMgp             ! index of last  ion relevant for cosmic ray ionisation
 integer,         parameter :: iitps             = iH2p             ! index of first ion to save for thermal ionisation only
 integer,         parameter :: iitpe             = iNap             ! index of last  ion to save for thermal ionisation only
 integer,         parameter :: iits              = iH               ! index of first neutral relevant for thermal ionisation
 integer,         parameter :: iite              = iNa              ! index of last  neutral relevant for thermal ionisation
 !--Indices for reactions
 integer,         parameter :: ir                =  1               ! Reaction rate
 integer,         parameter :: in_               =  2               ! Reaction * first number density
 integer,         parameter :: i_n               =  3               ! Reaction * second number density
 integer,         parameter :: inn               =  4               ! Reaction * both number densities
 !--Indicies for the specific k-k reactions (where k = neutral, ion, or electron):  r_kk array
 integer,         parameter :: iH3p_CO           =  1
 integer,         parameter :: iH3p_Mg           =  2
 integer,         parameter :: iH3p_e            =  3
 integer,         parameter :: ie_H3p            =  4
 integer,         parameter :: iHCOp_Mg          =  5
 integer,         parameter :: iHCOp_e           =  6
 integer,         parameter :: iHp_O             =  7
 integer,         parameter :: iHp_O2            =  8
 integer,         parameter :: iHp_Si            =  9
 integer,         parameter :: iHp_S             = 10
 integer,         parameter :: iHp_Mg            = 11
 integer,         parameter :: iHp_e             = 12
 integer,         parameter :: iHep_H2           = 13
 integer,         parameter :: iHep_CO           = 14
 integer,         parameter :: iHep_O2           = 15
 integer,         parameter :: iHep_Si           = 16
 integer,         parameter :: iHep_e            = 17
 integer,         parameter :: iCp_Si            = 18
 integer,         parameter :: iCp_S             = 19
 integer,         parameter :: iCp_e             = 20
 integer,         parameter :: iOp_e             = 21
 integer,         parameter :: iO2p_Si           = 22
 integer,         parameter :: iO2p_S            = 23
 integer,         parameter :: iO2p_e            = 24
 integer,         parameter :: iSip_Mg           = 25
 integer,         parameter :: iSip_e            = 26
 integer,         parameter :: iSp_Si            = 27
 integer,         parameter :: iSp_Mg            = 28
 integer,         parameter :: iSp_e             = 29
 integer,         parameter :: iMgp_e            = 30
 !--Indicies for the specific electron-Grain: r_eG array
 integer,         parameter :: ie_G0             =  1
 integer,         parameter :: ie_Gp             =  2
 !--Indicies for the specific k-Grain: r_kG array
 integer,         parameter :: ik_Gn             =  1
 integer,         parameter :: ik_G0             =  2
 !--Indicies for the specific Grain-Grain: r_GG array
 integer,         parameter :: iGp_Gn            =  1
 integer,         parameter :: iGp_G0            =  2
 integer,         parameter :: iGn_G0            =  3

 !--Array lengths
 integer, public, parameter :: n_data_out        = 22+3*na            ! number of array element the optional array data_out
 integer, public, parameter :: n_nden            = iire+1             ! number of species that require saving (elements from CR ionisation + thermal electrons)
 integer,         parameter :: nspecies          = ie                 ! The maximum number of charged species (internal array)
 integer,         parameter :: neqn              = iire               ! The number of equations to solve for cosmic ray ionisation
 !
 !--Local Variables to the NICIL module that may be required elsewhere in the user's code
 real,    public    :: meanmolmass,unit_eta

 !--Local Variables to the NICIL module
 integer, private   :: iprint,iprintw
 real,    private   :: unit_ndensity,unit_ndensity1
 real,    private   :: csqbyfourpi,threehundred1,small,small2,small21,smallXten,eta_thresh
 real,    private   :: warn_ratio2,warn_ratio_m1,warn_ratio_p1,warn_ratio_m12,warn_ratio_p12
 real,    private   :: mass_proton,mass_proton1,mass_neutral_cold1,one_minus_fdg,qe_c
 real,    private   :: eta_ohm_cnst,eta_hall_cnst,eta_ambi_cnst,alpha_AD_p1
 real,    private   :: zetaCRcoef,Saha_coef_HH2,chij_HH2,chiHH2
 real,    private   :: vrms2_min,vrms2_max
 real,    private   :: n_grain_coef(na)
 real,    private   :: mass_ion_cgs(nspecies),mass_ion_mp(nspecies),mass_ion(nspecies)
 real,    private   :: mass_neu_cgs(nspecies),mass_neu_mp(nspecies),mass_neu(nspecies)
 real,    private   :: abundance(nspecies),gjp1_gj(nspecies),Texp_thresh(nspecies)
 real,    private   :: Saha_coef(nspecies),chij(nspecies)
 real,    private   :: Zj(nspecies),aZj(nspecies),sZj(nspecies)
 real,    private   :: beta_coef(nspecies),nu_jn_coef(nspecies),sigmav(nspecies)
 real,    private   :: sqrt_q2_by_ak(na),q2_by_ak(na),ak_by_q2(na), &
                       q2_by_aak(na,na),aak_by_q2(na,na),sqrt_q2_by_aak(na,na), &
                       r_GG_coef(3,na,na),r_kG_coef(na,iire),r_eG_coef(na)
 character(len=2),private   :: symj(nspecies)

 !--Subroutines
 public  :: nicil_initialise,nicil_update_nimhd,nicil_translate_error
 public  :: nicil_get_dudt_nimhd,nicil_get_dt_nimhd
 public  :: nicil_get_vion,nicil_get_ambidrift,nicil_get_halldrift
 public  :: nimhd_get_jcbcb,nimhd_get_dBdt
 private :: nicil_nimhd_get_eta,nicil_nimhd_get_eta_cnst,nimhd_get_DcrossR
 private :: nicil_initialise_species,nicil_ic_error,nicil_print_summary,nicil_version
 private :: nicil_ion_get_sigma,nicil_ionR_get_n,nicil_get_HH2_ratio
 private :: nicil_set_reaction_rates,nicil_reaction_rates_X_ncharged
 private :: nicil_ionR_calc_Jf,get_dn_viaLU
 private :: nicil_ionT_get_ne,nicil_ionT_get_nion

 private

contains

!======================================================================!
! VERSION CONTROL & HISTORY                                            !
!======================================================================!
!----------------------------------------------------------------------!
!+
! Internal Version Control
! important modifications are listed before version number
!+
!----------------------------------------------------------------------!
pure subroutine nicil_version(version)
 character(len=200), intent(out) :: version
           !  8 Dec  2015: Initial Version
 version = "Version 1.0: 8 Dec 2015: Initial Version" ! commit c2203ba
           ! 27 Jan  2016: If iterations fail to converge, will reset initial guess and try again
           !  3 Mar  2016: Thermal ionisation can doubly ionise atoms
           !  8 Mar  2016: Added light ion for Cosmic ray ionisation, thus there is now heavy and light ions
           ! 11 Mar  2016: Rewrote MRN grain calculations to be sums with characteristic sizes, and not a single average
 version = "Version 1.1: 26 April 2016"
           !  7 June 2016: bug fix: collision rates account for arbitrary charge
           !  7 June 2016: bug fix: added statistical weights to Saha equation
           !  8 June 2016: changed default: mod_beta = .false.
           !               modified rate equations to be from Table 1 of Pinto & Galli (2008)
           !  9 June 2016: added warnings that can be optionally printed to file
           !               added optional subroutine to calculate ion velocity
           !               added dissociation of molecular hydrogen for self-consistent fractions of H and H2
           ! 27 June 2016: Rewrote cosmic ray ionisation algorithm.  Solves for n_g(Z=-1),n_g(Z=0),n_g(Z=+1) rather than \bar{Z}
           !               Solving for \bar{Z} remains as a contingency if convergence is not obtained with the above method
 version = "Version 1.2: 27 June 2016"
           !  6 July 2016: bug fix in nimhd_get_dudt
           !  1 Aug  2016: will exit cleanly if B=0 is passed in
           !               bug fix when calculating k_ei(1); temperature dependence now included
 version = "Version 1.2.1: 2 Aug 2016"      ! commit 8220ec8.  Version used in Wurster (2016)
           !  6 Sept 2016: bug fix in v_ion
           ! 12 Sept 2016: Updated references
           ! 22 Sept 2016: bug fix in eta_ambi if eta_constant=true & eta_calc_const=false
           !  8 Nov  2016: added a third option when choosing how to calculate (semi-)constant resistivities
           !  6 Jan  2017: corrected how grain charges are calculated if using the average-Z method
           ! 13 Mar  2017: bug fix in calculation of rho_n
           !               mass_neutral is dynamically calculated based upon Temperature to account for correct fractions of H, H2
 version = "Version 1.2.2: 14 March 2017"
           ! 20 Mar  2017: added option for density-dependent cosmic ionisation rate (off by default)
           !  3 July 2017: added option to pass in the dust-to-gas ratio; will (optionally) pass back ratio of charged dust
           !               corrected calculation of n_n in the Jacobian equations
           !               added option to re-order the Jacobian for stability (off by default due to expense)
           !               general cleanup
           !  4 July 2017: added subroutines to calculate the ion and hall drift velocities; augmented the current v_ion routine
 version = "Version 1.2.3: 5 July 2017"     ! commit 324833c
           ! 10 Apr  2018: Ambipolar diffusion now uses a subtraction rather than a double loop
 version = "Version 1.2.4: 11 April 2018"   ! commit b0d825b
           ! 19 Oct  2018: if rho_n < 0, will scale the neutral components from thermal and rays by their respective n_electron
           !               rho_is_rhogas = .false. by default
           !               n_gas rather than n_total is used for calculation thermal ionisation (if rho_is_rhogas = .false.)
 version = "Version 1.2.5: 19 October 2018" ! commit 28a817d
           !  4 Feb  2019: bug fix when using MRN grain distribution and size(n_R) .ne. 2+2*na
 version = "Version 1.2.6: 28 August 2019"  ! commit 679b501
           ! CLEANING: removing obsolete features in advance of a major rewrite of the chemistry
           ! 19 Sept 2019: removing the option to use the modified \beta
           !               removing the option to use mass fractions of H,He; must use abundances
           !               user *must* set size(n_R) == nimass + 2na in host code
           !               re-wrote the error bookkeeping
           !               Removing option of doubly ionising particles in thermal ionisation
           ! 20 Sept 2019: merged nicil_get_ion_n and nicil_get_eta into nicil_update_nimhd
 version = "Version 1.3: 20 September 2019"  ! commit 5178486
           ! 23 Jan  2020: Rewrote the chemistry for cosmic ray ionisation
           !               Number of grains is now controlled by hardcoding na
           ! 27 Jan  2020: Restructured density calculations: calculates thermal ionisation, then the CR ionisation with the remaining neutrals
           !  9 Feb  2020: Modified the grain reactions to match Kunz & Mouschovias (2009)
           !               Added dust evaporation between 725 < T/K < 1700; above 1700K, there is no dust
           !               Modified cosmic ray chemistry for stability and to exit with nden_ionR = 0 if it cannot converge
           ! 15 Feb  2020: Cleaning up code to remove obsolete variables and subroutines, including nimass, nelements_max and average Z failsafe
           !               Removed ionisation from radionuclide decay
           !               Removed ion_rays and ion_thermal options such that they are always on
           !               Removed option to pass back fraction of charged dust
           !               Replaced C_nimhd with Cdt_diff for ambipolar diffusion and Ohmic resistivity Cdt_hall = 1/4pi for the Hall effect
           ! 20 Feb  2020: Optional argument of input dust-to-gas ratios is now an na-array
           !               Can optionally pass in an na-array of grain radii in to the initialisation routine
           !               removed rho_is_rhogas
           ! 22 Feb  2020: Replaced loops with array operations for enhanced performance (assumed results)
           ! 23 Feb  2020: Replaced array operations with loops for enhanced performance (after time-tests)
 version = "Version 2.0: 23 February 2020" ! commit 16c498a
           ! 27 May  2020: Bug fix for grain coefficients when inputting dust-to-gas ratio but not initial grain sizes
           !               Slope of grain distribution (if use_fdg_in = .false. and not passing in grain sizes) is an input parameter
           ! 28 May  2020: Bug fix for consistent results with the legal case of use_fdg_in = .false. but fdg passed in
           !               Renamed subroutines & altered public/private status
           !               Ion velocity includes the effect of both ambipolar diffusion and the Hall effect
           ! 31 May  2020: Updated error bookkeeping: can now keep a running total of errors & non-fatal are hardcoded to be < 0
           ! 17 Jun  2020: Added optional initial guess of nden_save to be calculated in the initialisation routine
           ! 22 July 2020: Added minimum ionisation rate when using zeta_of_rho = .true.
 version = "Version 2.0.1: 30 July 2020" ! commit a63fdec
           ! 13 Nov  2020: Added chemical species Silicon and Sulfur
           !               Reduced loop-dependence on chemical species tags
 version = "Version 2.1: 13 November 2020" ! commit 201dc39.  Version used in Wurster (2021)
           ! 10 Oct  2022: Added additional warnings if eta is too large
           !               Will zero densities if they are likely noise
           !               Updated initial guesses when solving chemical network
           !  2 Dec  2022: Updated initial guesses when solving chemical network again (possibly unnecessary due to a bug in a parent code)
           !               Added optional out to track which guess number for the chemical network lead to converging results
           !               After convergence of chemical network, will accept results if negative values are clearly due to machine roundoff (values are zeroed)
           !               Updated error message for non-convergence of chemical network
           !               When passing in dust fractions, will zero fractions that are << 1 to prevent machine underflow
           !               Corrected zeroing of densities in Saha solver
 version = "Version 2.1.1: 2 December 2022" ! commit a90c270
           ! 30 Mar  2023: data_out(5) is total number density rather than an incorrect neutral number density


end subroutine nicil_version
!======================================================================!
! INITIALISATION ROUTINES                                              !
! Routines that state the chemical elements, their properties, and     !
! all reactions.  These can be used to understand the ionisation       !
! occurring within NICIL.                                              !
!======================================================================!
!----------------------------------------------------------------------!
!+
! Define the properties of the species
!+
!----------------------------------------------------------------------!
subroutine nicil_initialise_species
 integer :: k

 !--Chemical Symbols
 symj(iCO) = "CO"
 symj(iO2) = "O2"
 symj(iH ) = "H"
 symj(iHe) = "He"
 symj(iC ) = "C"
 symj(iO ) = "O"
 symj(iSi) = "Si"
 symj(iS ) = "S"
 symj(iMg) = "Mg"
 symj(iH2) = "H2"
 symj(iK ) = "K"
 symj(iNa) = "Na"

 !--Neutral abundances (as a fraction of n_H); assumes all hydrogen is molecular
 !  NOTE: This means that, for any given species, there may be more ionised gas than neutral gas.
 !        We might need to instead set these as neutral+ionised abundances, but that would add extra
 !        cost with little benefit; although the Jacobian would remain the same size, there would be
 !        additional elements since neutral gas = (total - ionised).  Preliminary tests (Fall 2022)
 !        show that a partial implementation has a negligible effect at low ionisations and a small
 !        effect at high ionisations; given the smaller values of eta at high ionisations, this
 !        assumption that these are the neutral and not neutral+ion abundances will have minimal effect
 !        on a simulation.
 abundance      = 0.0
 abundance(iH2) = 0.921 * 0.5
 abundance(iHe) = 0.0784
 abundance(iC ) = 7.74d-5
 abundance(iO ) = 3.78d-6
 abundance(iO2) = 4.42d-6
 abundance(iCO) = 3.00d-5
 abundance(iS ) = 2.80d-5
 abundance(iSi) = 1.70d-6
 abundance(iMg) = 1.56d-6
 abundance(iK ) = 2.20d-10
 abundance(iNa) = 3.10d-9

 !--Masses of neutral particle (amu)
 mass_neu_mp      =  0.0
 mass_neu_mp(iH ) =  1.00784
 mass_neu_mp(iHe) =  4.0026
 mass_neu_mp(iC ) = 12.0107
 mass_neu_mp(iO ) = 15.9994
 mass_neu_mp(iO2) = mass_neu_mp(iO)*2.0
 mass_neu_mp(iSi) = 28.0855
 mass_neu_mp(iS ) = 32.065
 mass_neu_mp(iMg) = 24.305
 mass_neu_mp(iCO) = mass_neu_mp(iC) + mass_neu_mp(iO)
 mass_neu_mp(iH2) =  2.01410
 mass_neu_mp(iK ) =  39.0983
 mass_neu_mp(iNa) =  22.9898

 !--Masses of ions (amu)
 mass_ion_mp        = 0.0
 mass_ion_mp(ie   ) = mass_electron_cgs/mass_proton_cgs
 mass_ion_mp(iH3p ) = 3.02293
 mass_ion_mp(iHCOp) = mass_neu_mp(iH) + mass_neu_mp(iC) + mass_neu_mp(iO) - mass_ion_mp(ie)
 do k = iO2p,iNap
    mass_ion_mp(k)  = mass_neu_mp(k ) - mass_ion_mp(ie)
 enddo

 !--Dissociation potential [eV] of H2
 chiHH2    =  4.476

 !--First ionisation potential [eV]
 chij      = 0.0
 chij(iH ) = 13.60
 chij(iHe) = 24.59
 chij(iC ) = 11.26
 chij(iO ) = 13.62
 chij(iS ) = 10.36
 chij(iSi) =  8.15
 chij(iMg) =  7.65
 chij(iH2) = 15.60
 chij(iK ) =  4.34
 chij(iNa) =  5.14

 !--Ratio of statistical weights (first/ground)
 gjp1_gj      = 0.0
 gjp1_gj(iH ) = 1.0/2.0
 gjp1_gj(iHe) = 2.0/1.0
 gjp1_gj(iC ) = 2.0/1.0
 gjp1_gj(iO ) = 2.0/1.0
 gjp1_gj(iS ) = 2.0/1.0
 gjp1_gj(iSi) = 2.0/1.0
 gjp1_gj(iMg) = 2.0/1.0
 gjp1_gj(iH2) = 1.0/2.0
 gjp1_gj(iNa) = 1.0/2.0
 gjp1_gj(iK ) = 1.0/2.0

end subroutine nicil_initialise_species
!----------------------------------------------------------------------!
!+
! A list of all the reactions included in the cosmic ray chemistry
! The reaction chemistry is performed in CGS units from McElroy+ (2013)
! with the coefficients in units of cm^3/s
!+
!----------------------------------------------------------------------!
pure subroutine nicil_set_reaction_rates(T,zeta,nden_neu,r_kk,r_eG,r_kG,r_GG,zeta_n)
 real, intent(in)  :: T,zeta,nden_neu(:)
 real, intent(out) :: r_kk(:,:),r_eG(:,:,:),r_kG(:,:,:,:),r_GG(:,:,:,:),zeta_n(:)
 integer           :: k,j,jj,j1,iGpj,iGnj
 real              :: Tby300,T1,sqrtT,sqrtT1,termA,termB

 !--Initialise temperature variables
 Tby300 = T*threehundred1
 T1     = 1.0/T
 sqrtT  = sqrt(T)
 sqrtT1 = 1./sqrtT

 !--Ion destruction
 r_kk(iH3p_CO, ir) = 1.36d-9 *Tby300**(-0.14)*exp(3.4*T1)    ! H3p  + CO  -> HCOp + H2       (T = 10-  400)
 r_kk(iH3p_Mg, ir) = 1.00d-9                                 ! H3p  + Mg  -> Mgp  + H  + H2  (T = 10-41000)
 r_kk(iH3p_e,  ir) = 4.36d-8 *Tby300**(-0.52)                ! H3p  + e   -> H2   + H        (T = 10- 1000)
 r_kk(ie_H3p,  ir) = 2.34d-8 *Tby300**(-0.52)                ! H3p  + e   -> 3H              (T = 10- 1000)
 r_kk(iHCOp_Mg,ir) = 2.90d-9                                 ! HCOp + Mg  -> Mgp  + HCO      (T = 10-41000)
 r_kk(iHCOp_e, ir) = 2.40d-7 *Tby300**(-0.69)                ! HCOp + e   -> CO   + H        (T = 10-  300)
 r_kk(iHp_O,   ir) = 6.86d-10*Tby300**( 0.26)*exp(-224.3*T1) ! Hp   + O   -> Op   + H        (T = 10-41000)
 r_kk(iHp_O2,  ir) = 2.00d-9                                 ! Hp   + O2  -> O2p  + H        (T = 10-  300)
 r_kk(iHp_Si,  ir) = 9.9d-10                                 ! Hp   + Si  -> Sip  + H        (T = 10-41000)
 r_kk(iHp_S ,  ir) = 1.30d-9                                 ! Hp   + S   -> Sp   + H        (T = 10-41000)
 r_kk(iHp_Mg,  ir) = 1.10d-9                                 ! Hp   + Mg  -> Mgp  + H        (T = 10-41000)
 r_kk(iHp_e,   ir) = 3.50d-12*Tby300**(-0.75)                ! Hp   + e   -> H    + photon   (T = 10-20000)
 r_kk(iHep_H2, ir) = 3.70d-14                *exp(-35.*T1)   ! Hep  + H2  -> Hp   + H  + He  (T = 10-  300)
 r_kk(iHep_CO, ir) = 1.60d-9                                 ! Hep  + CO  -> Cp   + O  + He  (T = 10-41000)
 r_kk(iHep_O2, ir) = 1.10d-9                                 ! Hep  + O2  -> Op   + O  + He  (T = 10-41000)
 r_kk(iHep_Si, ir) = 3.30d-9                                 ! Hep  + Si  -> Sip  + He       (T = 10-41000)
 r_kk(iHep_e,  ir) = 5.36d-12*Tby300**(-0.50)                ! Hep  + e   -> He   + photon   (T = 10- 1000)
 r_kk(iCp_Si,  ir) = 2.10d-9                                 ! Cp   + Si  -> Sip  + C        (T = 10-41000)
 r_kk(iCp_S,   ir) = 5.00d-11                                ! Cp   + S   -> Sp   + C        (T = 10-41000)
 r_kk(iCp_e,   ir) = 2.36d-12*Tby300**(-0.29)*exp(17.6*T1)   ! Cp   + e   -> C    + photon   (T = 10- 1000)
 r_kk(iOp_e,   ir) = 3.24d-12*Tby300**(-0.66)                ! Op   + e   -> O    + photon   (T = 10- 1000)
 r_kk(iO2p_Si, ir) = 1.60d-9                                 ! O2p  + Si  -> Sip  + O2       (T = 10-41000)
 r_kk(iO2p_S , ir) = 5.4d-10                                 ! O2p  + S   -> Sp   + O2       (T = 10-41000)
 r_kk(iO2p_e,  ir) = 1.95d-7 *Tby300**(-0.70)                ! O2p  + e   -> O    + O        (T = 10-  300)
 r_kk(iSip_Mg, ir) = 2.90d-9                                 ! Sip  + Mg  -> Mgp  + Si       (T = 10-41000)
 r_kk(iSip_e,  ir) = 4.26d-12*Tby300**(-0.62)                ! Sip  + e   -> Si   + photon   (T = 10-41000)
 r_kk(iSp_Si,  ir) = 1.60d-9                                 ! Sp   + Si  -> Sip  + S        (T = 10-41000)
 r_kk(iSp_Mg,  ir) = 2.80d-10                                ! Sp   + Mg  -> Mgp  + S        (T = 10-41000)
 r_kk(iSp_e,   ir) = 5.49d-12*Tby300**(-0.59)                ! Sp   + e   -> S    + photon   (T = 10- 1000)
 r_kk(iMgp_e,  ir) = 2.78d-12*Tby300**(-0.68)                ! Mgp  + e   -> Mg   + photon   (T = 10- 1000)

 !--Multiply the reactions by their neutral number densities
 r_kk(iH3p_CO, i_n) = r_kk(iH3p_CO, ir)*nden_neu(iCO)
 r_kk(iH3p_Mg, i_n) = r_kk(iH3p_Mg, ir)*nden_neu(iMg)
 r_kk(iHCOp_Mg,i_n) = r_kk(iHCOp_Mg,ir)*nden_neu(iMg)
 r_kk(iHp_O2,  i_n) = r_kk(iHp_O2,  ir)*nden_neu(iO2)
 r_kk(iHp_Si,  i_n) = r_kk(iHp_Si,  ir)*nden_neu(iSi)
 r_kk(iHp_S ,  i_n) = r_kk(iHp_S ,  ir)*nden_neu(iS)
 r_kk(iHp_O,   i_n) = r_kk(iHp_O,   ir)*nden_neu(iO)
 r_kk(iHp_Mg,  i_n) = r_kk(iHp_Mg,  ir)*nden_neu(iMg)
 r_kk(iHep_Si, i_n) = r_kk(iHep_Si, ir)*nden_neu(iSi)
 r_kk(iHep_CO, i_n) = r_kk(iHep_CO, ir)*nden_neu(iCO)
 r_kk(iHep_O2, i_n) = r_kk(iHep_O2, ir)*nden_neu(iO2)
 r_kk(iHep_H2, i_n) = r_kk(iHep_H2, ir)*nden_neu(iH2)
 r_kk(iCp_Si,  i_n) = r_kk(iCp_Si,  ir)*nden_neu(iSi)
 r_kk(iCp_S,   i_n) = r_kk(iCp_S,   ir)*nden_neu(iS)
 r_kk(iO2p_Si, i_n) = r_kk(iO2p_Si, ir)*nden_neu(iSi)
 r_kk(iO2p_S , i_n) = r_kk(iO2p_S , ir)*nden_neu(iS)
 r_kk(iSip_Mg, i_n) = r_kk(iSip_Mg, ir)*nden_neu(iMg)
 r_kk(iSp_Si,  i_n) = r_kk(iSp_Si,  ir)*nden_neu(iSi)
 r_kk(iSp_Mg,  i_n) = r_kk(iSp_Mg,  ir)*nden_neu(iMg)

 !--Cosmic Ray Equations
 !  these values are relative to the H2 ionisation rate, zeta
 !  the 2 in the iH3p equation is a coefficient for the density, not the CRIR
 zeta_n       = 0.                                         ! initialise to prevent compiler warnings
 zeta_n(iH3p) = 2.0   *nden_neu(iH2)                       ! 2H2 +       cr -> H3p + H + e
 zeta_n(iHp)  = 0.0248*nden_neu(iH2) + 0.498*nden_neu(iH)  ! H2  + 0.0248cr -> Hp  + H + e & H + 0.498cr -> Hp + e
 zeta_n(iHep) = 0.542 *nden_neu(iHe)                       ! He  + 0.542 cr -> Hep     + e
 zeta_n(iCp)  = 1.92  *nden_neu(iC)                        ! C   + 1.92  cr -> Cp      + e
 zeta_n(iOp)  = 2.83  *nden_neu(iO)                        ! O   + 2.83  cr -> Op      + e
 zeta_n       = zeta_n * zeta

 !--Grain Reactions
 r_kG = 0.0
 do j = 1,na
    jj   = 2*(j-1)
    iGpj = iGp+jj
    iGnj = iGn+jj
    ! Grain-grain reactions
    do j1 = 1,na
       r_GG(iGp_Gn,j,j1,ir) = r_GG_coef(iGp_Gn,j,j1)*sqrtT*(1.0+q2_by_aak(j,j1)*T1)*(1.0+sqrt(2.0/(2.0+aak_by_q2(j,j1)*T))) ! Gnj + Gpk -> G0j + G0k
       r_GG(iGp_G0,j,j1,ir) = r_GG_coef(iGp_G0,j,j1)*sqrtT*(1.0+sqrt_q2_by_aak(j,j1)*sqrtT1)                                ! Gpj + G0k -> G0j + Gpk
       r_GG(iGn_G0,j,j1,ir) = r_GG(iGp_G0,j,j1,ir)                                                                          ! Gnj + G0k -> G0j + Gnk
    enddo
    termA = sqrtT*(1.0+sqrt_q2_by_ak(j)*sqrtT1)
    termB = sqrtT*(1.0+q2_by_ak(j)*T1)*(1.0+sqrt(2.0/(2.0+ak_by_q2(j)*T)))
    ! Grain-charged ion reactions
    do k = iirs,iire
       r_kG(ik_G0,j,k,ir) = r_kG_coef(j,k)*termA             ! G0 + i  -> Gp + 0
       r_kG(ik_Gn,j,k,ir) = r_kG_coef(j,k)*termB             ! Gn + i  -> G0 + 0
    enddo
    ! Grain-electron reactions
    r_eG(ie_G0,j,ir) = r_eG_coef(j)*termA                    ! G0 + e  -> Gn
    r_eG(ie_Gp,j,ir) = r_eG_coef(j)*termB                    ! Gp + e  -> G0
 enddo

end subroutine nicil_set_reaction_rates
!======================================================================!
! INITIALISATION & CONTROL ROUTINES                                    !
!======================================================================!
!----------------------------------------------------------------------!
!+
! Initialisation subroutine
! This will initialise all the variable required for NICIL, including
! the frequently used coefficients.
! All coefficients will be converted to code units.
!+
!----------------------------------------------------------------------!
subroutine nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint_in,iprintw_in,nden_nimhd0,a_grain_cgs_in)
 real,              intent(in)  :: utime,umass,udist,unit_Bfield
 integer,           intent(out) :: ierr
 real,    optional, intent(in)  :: a_grain_cgs_in(na)
 integer, optional, intent(in)  :: iprint_in,iprintw_in
 real,    optional, intent(out) :: nden_nimhd0(n_nden)
 integer                        :: j,j1,iGnj,iGpj,iGpj1,k
 integer                        :: ierrlist(n_warn)
 real                           :: unit_density,unit_erg
 real                           :: sigmavgnbyT_coef,sigmav_coef
 real                           :: mass_neutral_cold_cgs
 real                           :: dloga,abundance_sum,rhog_sum
 real                           :: a_grain_sum,mass_grain_red
 real                           :: a_grain_cgs(na),rho_grain(na)
 real                           :: rdummy(6)

 !--Initialise species properties
 call nicil_initialise_species

 !--Verify input parameters are realistic; print error messages for each invalid error
 ierr = 0
 if (fdg  >= 1.0 .or. fdg < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid dust-to-gas fraction','fdg',fdg)
 if ( na==1 ) then
    if (a0_grain  < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid fixed grain radius:','a_grain',a0_grain)
 else
    if (an_grain < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid minimum MRN grain radius:','an_grain',an_grain)
    if (ax_grain < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid maximum MRN grain radius:','ax_grain',ax_grain)
    if (ax_grain < an_grain) call nicil_ic_error(ierr,'nicil_initialise','Invalid radii ordering of MRN grain radii')
    if (abs(grain_slope + 1.0) < epsilon(grain_slope)) then
       call nicil_ic_error(ierr,'nicil_initialise','Invalid grain slope of grain distribution:','grain_slope',grain_slope)
    endif
 endif
 if (rho_bulk < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid bulk grain density:','rho_bulk',rho_bulk)
 if (zeta_cgs < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid ionisation rate:','zeta', zeta_cgs)
 if (NRtol < 1.0e-15 .or. NRtol > 0.1) &
                      call nicil_ic_error(ierr,'nicil_initialise','Poor constraint on Newton–Raphson tolerance','NRtol',NRtol)
 if (NRctrmax <= 10) call nicil_ic_error(ierr,'nicil_initialise','Too few maximum permitted Newton–Raphson iterations' &
 ,'NRctrmax',real(NRctrmax))
 if (eta_constant .and. eta_const_type/=icnstphys .and. eta_const_type/=icnstsemi .and. eta_const_type/=icnst) then
    call nicil_ic_error(ierr,'nicil_initialise','Invalid choice of (semi-) constant eta calculations; correct eta_const_type.')
 endif
 !--Abort setup if errors fatal exist
 if (ierr/=0) return

 !--Unit conversions: cgs -> code
 unit_density      = umass/udist**3
 unit_erg          = umass*(udist/utime)**2
 unit_eta          = udist**2/utime
 unit_ndensity     = 1.0/udist**3
 unit_ndensity1    = 1.0/unit_ndensity

 !--Initialise parameters and constants
 small          = epsilon(small)
 small2         = small*small
 small21        = 1./small2
 smallXten      = 10.*small
 warn_ratio2    = warn_ratio*warn_ratio
 warn_ratio_m1  = 1.0-warn_ratio
 warn_ratio_p1  = 1.0+warn_ratio
 warn_ratio_m12 = warn_ratio_m1*warn_ratio_m1
 warn_ratio_p12 = warn_ratio_p1*warn_ratio_p1
 vrms2_min      = (vrms_min_kms*1.0d5*utime/udist)**2
 vrms2_max      = (vrms_max_kms*1.0d5*utime/udist)**2
 threehundred1  = 1.0/300.
 csqbyfourpi    = c**2/fourpi/unit_eta
 eta_thresh     = eta_thresh_cgs/unit_eta

 !--Calculate the mean molar mass of cold gas; assume all hydrogen is molecular
 meanmolmass   = 0.0
 abundance_sum = 0.0
 do j = 1,iite
    meanmolmass   = meanmolmass   + abundance(j)*mass_neu_mp(j)
    abundance_sum = abundance_sum + abundance(j)
 enddo
 meanmolmass = meanmolmass/abundance_sum
 mass_neutral_cold_cgs  = meanmolmass*mass_proton_cgs

 !--Determine grain distribution; initialise sizes and masses
 a_grain_cgs  = 0.0
 n_grain_coef = 0.0
 rho_grain    = 0.0
 rhog_sum     = 0.0
 do j = 1,na
    ! grain radii
    dloga = (log10(ax_grain) - log10(an_grain))/real(na)
    if (present(a_grain_cgs_in)) then
       a_grain_cgs(j) = a_grain_cgs_in(j)                         ! user's input grain sizes
    else
       if (na==1) then
          a_grain_cgs(j) = a0_grain                               ! constant grain size
       else
          a_grain_cgs(j) = 10**(log10(an_grain) +(j-0.5)*dloga )  ! grain sizes are evenly spaced in log-space
       endif
    endif
    ! grain masses
    mass_ion_mp(iGp+2*(j-1)) = fourpi/3.0*a_grain_cgs(j)**3*rho_bulk / mass_proton_cgs
    mass_ion_mp(iGn+2*(j-1)) = mass_ion_mp(iGp+2*(j-1))
    ! conversion coefficients from total mass density to grain number densities
    if (use_fdg_in) then
       n_grain_coef(j) = 1.0/mass_ion_mp(iGp+2*(j-1))             ! user passes in dust-to-gas mass ratio
    else
       if (na==1 .or. present(a_grain_cgs_in)) then
          n_grain_coef(j) = fdg/(na*mass_ion_mp(iGp+2*(j-1)))     ! if na > 1, this assumes that all dust is present with the same total mass
       else
          ! Assume a power-law grain distribution in number density (e.g. MRN if grain_slope=-3.5)
          n_grain_coef(j) = (( 10**(log10(an_grain) +(j-1)*dloga ) )**(grain_slope+1.0) &
                          - (10**(log10(an_grain) +(j)*dloga ))**(grain_slope+1.0))
       end if
       rho_grain(j) = n_grain_coef(j)*mass_ion_mp(iGp+2*(j-1))
       rhog_sum     = rhog_sum + rho_grain(j)
    end if
 enddo
 if (na > 1 .and. .not. use_fdg_in .and. .not. present(a_grain_cgs_in)) n_grain_coef = n_grain_coef/rhog_sum*fdg
 if (rhog_sum > 0.0) rho_grain = rho_grain/rhog_sum
 if (na > 0) then
    one_minus_fdg = 1.0 - fdg
 else
    one_minus_fdg = 1.0
 endif

 !--Convert masses to additional units
 mass_proton        = mass_proton_cgs           / umass
 mass_proton1       = 1.0/mass_proton
 mass_neutral_cold1 = 1.0/mass_neutral_cold_cgs * umass
 mass_ion_cgs       = mass_ion_mp   * mass_proton_cgs         ! convert from amu to cgs
 mass_neu_cgs       = mass_neu_mp   * mass_proton_cgs         ! convert from amu to cgs
 mass_ion           = mass_ion_cgs  / umass                   ! convert to code units
 mass_neu           = mass_neu_cgs  / umass                   ! convert to code units
 n_grain_coef       = n_grain_coef  * umass / mass_proton_cgs ! convert to code units

 !--Coefficient for the exponential term if density-dependent cosmic ionisation rate
 !  to be multiplied by sqrt (T rho/ m_n); udist converts rho/m_n to cgs, making this term dimensionless
 zetaCRcoef = sqrt(kboltz/(pi*Ggrav*mass_proton*udist**3))/SigmaCR

 !--Rate coefficients for grain reactions
 do j = 1,na
    sqrt_q2_by_ak(j) = sqrt(pi*qe**2/(2.0*a_grain_cgs(j)*kboltz))
    q2_by_ak(j)      = qe**2/(a_grain_cgs(j)*kboltz)
    ak_by_q2(j)      = a_grain_cgs(j)*kboltz/qe**2
    iGpj = iGp+2*(j-1)
    ! Grain-grain reactions
    do j1 = 1,na
       iGpj1 = iGp+2*(j1-1)
       a_grain_sum             = a_grain_cgs(j) + a_grain_cgs(j1)
       mass_grain_red          = mass_ion_cgs(iGpj)*mass_ion_cgs(iGpj1)/( mass_ion_cgs(iGpj) + mass_ion_cgs(iGpj1) )  ! recall mass(iGpj) = mass(iGpn)
       q2_by_aak(j,j1)         = qe**2/(a_grain_sum*kboltz)
       aak_by_q2(j,j1)         = a_grain_sum*kboltz/qe**2
       sqrt_q2_by_aak(j,j1)    = sqrt(pi*qe**2/(2.0*a_grain_sum*kboltz))
       r_GG_coef(iGp_Gn,j,j1) = a_grain_sum**2*sqrt(8.0*pi*kboltz/mass_grain_red)
       r_GG_coef(iGp_G0,j,j1) = r_GG_coef(iGp_Gn,j,j1)*a_grain_cgs(j1)**2/(a_grain_cgs(j)**2+a_grain_cgs(j1)**2)
    enddo
    ! Grain-charged ion reactions
    do k = iirs,iire
       r_kG_coef(j,k) = a_grain_cgs(j)**2*sqrt(8.0*pi*kboltz)/sqrt(mass_ion_cgs(k)  )*ProbI
    enddo
    ! Grain-electron reactions
    r_eG_coef(j)      = a_grain_cgs(j)**2*sqrt(8.0*pi*kboltz)/sqrt(mass_electron_cgs)*ProbE
 enddo

 !--Momentum transfer coefficients between neutrals and ions (CGS: cm^3 s^-1)
 !  use Langevin if there is no temperature-dependent value; temperature-dependent values are listed here for completeness
 sigmav_coef     = 2.81d-9*sqrt(pnH2)
 sigmav          = 0.
 !sigmav(iH3p) = 1.d-9*              (2.693 - (1.238 - (0.664 - 0.089 *logT)*logT)*logT)
 !sigmav(iHCOp)= max(0.0,1.d-9*sqrtT*(1.476 - (1.409 - (0.555 - 0.0775*logT)*logT)*logT))
 sigmav(iO2p)  = sigmav_coef*sqrt( (mass_ion(iO2p )+mass_neu(iH2))/(mass_ion(iO2p )*mass_neu_mp(iH2)))
 !sigmav(iHp)  = 1.d-9*              (1.003 + (0.050 + (0.136 - 0.014 *logT)*logT)*logT)
 sigmav(iHep)  = sigmav_coef*sqrt( (mass_ion(iHep )+mass_neu(iH2))/(mass_ion(iHep )*mass_neu_mp(iH2)))
 sigmav(iCp)   = sigmav_coef*sqrt( (mass_ion(iCp  )+mass_neu(iH2))/(mass_ion(iCp  )*mass_neu_mp(iH2)))
 sigmav(iOp)   = sigmav_coef*sqrt( (mass_ion(iOp  )+mass_neu(iH2))/(mass_ion(iOp  )*mass_neu_mp(iH2)))
 sigmav(iSip)  = sigmav_coef*sqrt( (mass_ion(iSip )+mass_neu(iH2))/(mass_ion(iSip )*mass_neu_mp(iH2)))
 sigmav(iSp)   = sigmav_coef*sqrt( (mass_ion(iSp  )+mass_neu(iH2))/(mass_ion(iSp  )*mass_neu_mp(iH2)))
 sigmav(iMgp)  = sigmav_coef*sqrt( (mass_ion(iMgp )+mass_neu(iH2))/(mass_ion(iMgp )*mass_neu_mp(iH2)))
 sigmav(iH2p)  = sigmav_coef*sqrt( (mass_ion(iH2p )+mass_neu(iH2))/(mass_ion(iH2p )*mass_neu_mp(iH2)))
 sigmav(iKp )  = sigmav_coef*sqrt( (mass_ion(iKp  )+mass_neu(iH2))/(mass_ion(iKp  )*mass_neu_mp(iH2)))
 sigmav(iNap)  = sigmav_coef*sqrt( (mass_ion(iNap )+mass_neu(iH2))/(mass_ion(iNap )*mass_neu_mp(iH2)))
 !sigmav(ie)   = 1.d-9*sqrtT*(0.535 + (0.203 - (0.163 - 0.050 *logT)*logT)*logT)
 sigmavgnbyT_coef  = delta_gn * sqrt(128.0*pi*kboltz/(9.0*mass_proton_cgs))

 !--Collisional Frequencies (CGS = s^-1)
 nu_jn_coef = 1.0/mass_proton_cgs * unit_density
 do j = 1,na
    nu_jn_coef(iGp+2*(j-1)) = sigmavgnbyT_coef*a_grain_cgs(j)**2/mass_proton_cgs * unit_density ! nu_gn
    nu_jn_coef(iGn+2*(j-1)) = nu_jn_coef(iGp+2*(j-1))
 enddo

 !--Fill charge arrays (electric charge,Z; absolute value of Z; sign of Z)
 Zj          =  1.0; aZj       = 1.0; sZj       =  1.0       ! Default since most ions are z = +1
 Zj(ie)      = -1.0; aZj(ie)   = 1.0; sZj(ie)   = -1.0       ! electrons
 do j = 1,na
    iGnj = iGn+2*(j-1)
    Zj(iGnj) = -1.0; aZj(iGnj) = 1.0; sZj(iGnj) = -1.0       ! negatively charged grains
 enddo

 !--Hall parameters (dimensionless after multipiled by B_code)
 beta_coef = aZj*qe/(mass_ion_cgs*c) * unit_Bfield

 !--Conductivity coefficient (code units)
 qe_c      = qe*c / (udist**3 * unit_Bfield)

 !--Coefficient for the Saha equation (Code units)
 Saha_coef_HH2 =             (    pi*mass_proton_cgs  *kboltz/planckh**2 * utime**2*unit_erg/umass )**1.5
 Saha_coef     = 2.0*gjp1_gj*(2.0*pi*mass_electron_cgs*kboltz/planckh**2 * utime**2*unit_erg/umass )**1.5
 chij          = chij*eV/kboltz
 chij_HH2      = chiHH2*eV/kboltz
 Texp_thresh   = chij*Texp_thresh0

 !--Set constant coefficients
 alpha_AD_p1 = 1.0
 if ( eta_const_type==icnst ) then
    eta_ohm_cnst  = C_OR / unit_eta
    eta_hall_cnst = C_HE / unit_eta
    eta_ambi_cnst = C_AD / unit_eta
 else
    if ( eta_const_type==icnstphys ) then
       eta_ohm_cnst  = mass_electron_cgs*c**2/(fourpi*qe**2*n_e_cnst)
       eta_hall_cnst = c/(fourpi*qe*n_e_cnst)
       eta_ambi_cnst = 1.0/(fourpi*gamma_AD*rho_i_cnst/rho_n_cnst**alpha_AD)
       alpha_AD_p1   = alpha_AD + 1.0                     ! alpha_AD /= 0 only makes sense for this case
       if (hall_lt_zero) eta_hall_cnst = -eta_hall_cnst
    else if ( eta_const_type==icnstsemi ) then
       ! the factors of fourpi are included to cancel out the sqrt(fourpi) that should be included in unit_Bfield
       eta_ohm_cnst  = C_OR
       eta_hall_cnst = C_HE / sqrt(fourpi)
       eta_ambi_cnst = C_AD /      fourpi
    endif
    !  Convert units as required
    eta_ohm_cnst  = eta_ohm_cnst                                               / unit_eta
    eta_hall_cnst = eta_hall_cnst *  unit_Bfield                               / unit_eta
    eta_ambi_cnst = eta_ambi_cnst * (unit_Bfield**2/unit_density**alpha_AD_p1) / unit_eta
 endif
 if (.not. use_ohm ) eta_ohm_cnst  = 0.0
 if (.not. use_hall) eta_hall_cnst = 0.0
 if (.not. use_ambi) eta_ambi_cnst = 0.0

 !--Set the location of the printouts
 if (present(iprint_in)) then
    iprint = iprint_in
 else
    iprint = 6
 endif
 if (present(iprintw_in)) then
    iprintw = iprintw_in
 else
    iprintw = 6
 endif

 !--Reset tolerances if required
 if (kind(NRtol)==4) then
    NRtol = NRtol4
    MPtol = MPtol4
    MPval = MPval4
    write(iprint,'(a,Es10.3)')'Resetting NR tolerance to ',NRtol
 endif

 !--Run through with known values to obtain reasonable initial conditions
 !  this will provide a useful initial guess for the regime of high densities and temperatures
 !  (from nicil_ex_disc for rho > 1d-11g/cm^3 & T > 1000K)
 if (present(nden_nimhd0)) then
    ierrlist    = 0
    nden_nimhd0 = 0.
    rdummy      = 1.0d-5
    rdummy(5)   = 1.0d-12/unit_density
    rdummy(6)   = 100.
    call nicil_update_nimhd(1,rdummy(1),rdummy(2),rdummy(3),rdummy(4),rdummy(5),rdummy(6),nden_nimhd0,ierrlist)
 endif

 !--Print Statements to summarise conditions used
 call nicil_print_summary(a_grain_cgs,rho_grain,meanmolmass,present(a_grain_cgs_in))

end subroutine nicil_initialise
!----------------------------------------------------------------------!
!+
! Print Statements to summarise what is being used
!+
!----------------------------------------------------------------------!
subroutine nicil_print_summary(a_grain_cgs,rho_grain,meanmolmass,grain_sizes_in)
 real,    intent(in) :: a_grain_cgs(:),rho_grain(:),meanmolmass
 logical, intent(in) :: grain_sizes_in
 integer             :: j
 real                :: nfrac
 character(len=  2)  :: comma
 character(len=200)  :: ni_terms,version,fmt

 call nicil_version(version)
 write(iprint,'(a)' ) "NICIL:"
 write(iprint,'(2a)') "NICIL: ",trim(version)
 write(iprint,'(a)' ) "NICIL: Copyright (c) 2015-2023 James Wurster"
 write(iprint,'(a)' ) "NICIL: See LICENCE file for usage and distribution conditions"
 write(iprint,'(a)' ) "NICIL: References: Wurster (2016) PASA, 33:e041."
 write(iprint,'(a)' ) "NICIL:             Wurster (2021) MNRAS, 501:5873-5891."
 write(iprint,'(a)' ) "NICIL:"

 ni_terms = ""
 comma    = ""
 if (use_ohm ) then
    write(ni_terms,'(2a)') trim(ni_terms),"Ohmic Resistivity"
    comma = ", "
 endif
 if (use_hall) then
    write(ni_terms,'(3a)') trim(ni_terms),comma,"Hall Effect"
    comma = ", "
 endif
 if (use_ambi) write(ni_terms,'(3a)') trim(ni_terms),comma,"Ambipolar diffusion"

 write(iprint,'(a)') "NICIL: Including ionisation from Cosmic rays & thermal ionisation"
 if(use_ohm .or. use_hall .or. use_ambi) write(iprint,'(2a)')"NICIL: Non-ideal terms used: ",trim(ni_terms)
 if (eta_constant) then
    if (eta_const_type==icnst) then
       write(iprint,'(a)' )                       "NICIL: All resistivity coefficients are constant."
       if (use_ohm ) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ohm  = ", eta_ohm_cnst*unit_eta,  " cm^2 s^{-1}"
       if (use_hall) write(iprint,'(a,Es10.3,a)') "NICIL: eta_Hall = ", eta_hall_cnst*unit_eta, " cm^2 s^{-1}"
       if (use_ambi) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ambi = ", eta_ambi_cnst*unit_eta, " cm^2 s^{-1}"
    else
       write(iprint,'(a)' )                       "NICIL: All resistivity coefficients are semi-constant."
       if (use_ohm ) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ohm  = ", eta_ohm_cnst*unit_eta, "           cm^2 s^{-1}"
       if (use_hall) write(iprint,'(a,Es10.3,a)') "NICIL: eta_Hall = ", eta_hall_cnst*unit_eta,"*|B|       cm^2 s^{-1}"
       if (use_ambi) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ambi = ", eta_ambi_cnst*unit_eta,"*|B|^2/rho cm^2 s^{-1}"
    endif
    write(iprint,'(a,Es10.3,a)') "NICIL: Mean molecular mass:             ",meanmolmass,"m_p"
 else
    if (na==1) then
       write(iprint,'(a)')    "NICIL: Using constant grain size."
    else
       write(iprint,'(a)')    "NICIL: Approximating grain distribution with MRN distribution."
    endif
    write(iprint,'(a)')                          "NICIL: Species: abundance relative to H, assuming only molecular hydrogen"
    j=iH2; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iH;  write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iHe; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iC;  write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iO;  write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iO2; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iMg; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iSi; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iS;  write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iK;  write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iNa; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    j=iCO; write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),":     ",abundance(j)
    write(iprint,'(a,Es10.3,a)')       "NICIL: Mean molecular mass of cold gas: ",meanmolmass,"m_p"
    if (use_fdg_in) then
       write(iprint,'(a         )') "NICIL: Importing local dust-to-gas ratio from parent code"
    else
       write(iprint,'(a,f10.4  )') "NICIL: Dust-to-gas fraction:            ",fdg
       write(iprint,'(a        )') "NICIL: rho_gas = rho_in*(1 - fdg)"
    endif
    if (na==1) then
       write(iprint,'(a         )') "NICIL: There is a single grain population"
       write(iprint,'(a,Es10.3,a)') "NICIL: Grain mass:                      ",mass_ion_mp(1),"m_p"
       write(iprint,'(a,Es10.3,a)') "NICIL: Grain radius:                    ",a_grain_cgs(1),"cm"
    else
       if (use_fdg_in) then
          write(iprint,'(a,I4    ,a)') &
          "NICIL: There are ",na," grain sizes, with radii & masses of"
          do j = 1,na
             write(iprint,'(a,I4,a,2(Es10.3,a))') "NICIL: size ",j,": radius = ",a_grain_cgs(j),"cm & mass = ", &
                                                  mass_ion_mp(j),"m_p"
          enddo
       else
          if (grain_sizes_in) then
             write(iprint,'(a,I4,a)') "NICIL: There are ",na," grain sizes"
          else if (abs(grain_slope + 3.5) < epsilon(grain_slope)) then
             write(iprint,'(a,I4,a)') &
             "NICIL: There are ",na," grain sizes that follow an MRN distribution with a slope of -3.5"
          else
             write(iprint,'(a,I4,a,F6.3)') &
             "NICIL: There are ",na," grain sizes that follow a power-law distribution with a slope of ",grain_slope
          endif
          write(iprint,'(a,I4,a)') &
          "NICIL: The grain radii, masses, relative mass density fractions & relative number density fractions are"
          do j = 1,na
             nfrac = n_grain_coef(j)/sum(n_grain_coef(:))
             if (nfrac > 1.0d-5) then
                fmt = '(a,I4,a,2(Es10.3,a),F8.6,a,1x,F8.6)'
             else
                fmt = '(a,I4,a,2(Es10.3,a),F8.6,a,Es10.3)'
             endif
             write(iprint,fmt) "NICIL: size ",j,": radius = ",a_grain_cgs(j),"cm & mass = ",mass_ion_mp(j), &
                               "m_p & rel. mass density = ", rho_grain(j), " & rel. number density = ", nfrac
          enddo
       endif
    endif
    if (zeta_of_rho) then
       write(iprint,'(a,Es10.3,a)') "NICIL: Unattenuated cosmic ray ionisation rate:     ",zeta_cgs,    " s^{-1}"
       write(iprint,'(a,Es10.3,a)') "NICIL: Minimum allowed ionisation rate:             ",zeta_cgs_min," s^{-1}"
    else
       write(iprint,'(a,Es10.3,a)') "NICIL: Cosmic ray ionisation rate:      ",zeta_cgs," s^{-1}"
    endif
    write(iprint,'(a,f10.4)') "NICIL: Coefficient for the ambipolar & Ohmic timesteps: ",Cdt_diff
    write(iprint,'(a,f10.4)') "NICIL: Coefficient for the Hall timestep:               ",Cdt_hall
 endif

end subroutine nicil_print_summary
!----------------------------------------------------------------------!
!+
! These will print an error message for each error found in NICIL.
! This will NOT end the main programme, but pass the error code to the
! host code such that the code can be properly and cleanly terminated.
! The first subroutine is specifically for initialising NICIL, and
! the second and third subroutines are is for runtime.
!+
!----------------------------------------------------------------------!
!--Initialisation error messages
subroutine nicil_ic_error(num_errors,wherefrom,reason,var,val)
 integer,                    intent(inout) :: num_errors
 character(len=*),           intent(in)    :: wherefrom,reason
 character(len=*), optional, intent(in)    :: var
 real,             optional, intent(in)    :: val

 !--Print cause of error
 if (present(var) .and. present(val)) then
    write(iprint,'(7a,Es16.4)') 'NICIL: ERROR: ',wherefrom,': ',reason,'. ',var,'=',val
 else
    write(iprint,'(4a)') 'NICIL: ERROR: ',wherefrom,': ',reason
 endif
 num_errors = num_errors + 1

end subroutine nicil_ic_error
!----------------------------------------------------------------------!
!--Runtime error messages
subroutine nicil_translate_error(ierrlist,fatal_only,rho,B,T,eta_ohm,eta_hall,eta_ambi,nden_save,fdg_in)
 integer,        intent(in) :: ierrlist(:)
 logical,        intent(in) :: fatal_only
 real, optional, intent(in) :: rho,B,T,eta_ohm,eta_hall,eta_ambi,nden_save(:),fdg_in(:)
 integer                    :: i
 logical                    :: print_n_fdg
 character(len= 16)         :: werr
 character(len= 20)         :: ferr
 character(len= 96)         :: errmsg(n_warn)
 character(len=512)         :: parspace

 !--The fatal error messages
 ferr       = 'NICIL: FATAL ERROR: '
 errmsg(ierr_neTconv) = 'nicil_ionT_get_ne: n_electronT did not converge'
 errmsg(ierr_neTle0 ) = 'nicil_ionT_get_ne: n_electronT < 0'
 errmsg(ierr_fdg_in ) = 'calling error: use_fdg_in=.true., but fdg_in not passed in'
 !--The warning error messages
 werr       = 'NICIL: WARNING: '
 errmsg(ierr_nRconv ) = 'nicil_ionR_get_n: n_{i,e,g} did not converge or yield valid values.  Resetting n_ionR = n_eR = 0'
 errmsg(ierr_highfdg) = 'grain fraction is too high.  Setting eta = 0'
 errmsg(ierr_T      ) = 'input error: T == 0.  Verify your code can legitimately input this.'
 errmsg(ierr_B      ) = 'input error: B == 0.  Verify your code can legitimately input this.'
 errmsg(ierr_sigH   ) = 'sigma_Hall == 0 (would otherwise be dominate by round-off errors)'
 errmsg(ierr_scN    ) = 'strong coupling approximation (rho_n~ rho) invalid since rho_n < 0.1rho'
 errmsg(ierr_scI    ) = 'strong coupling approximation (rho_i<<rho) invalid since rho_i > 0.01rho'
 errmsg(ierr_drift  ) = 'drift velocity ~ 0 approximation invalid since v_d > 0.1v'
 errmsg(ierr_vrms   ) = 'vrms_min < v < vrms_max is not true'
 errmsg(ierr_alion  ) = 'gas is fully ionised'
 errmsg(ierr_or     ) = 'eta_ohm is above the given threshold'
 errmsg(ierr_he     ) = '|eta_hall| is above the given threshold'
 errmsg(ierr_ad     ) = 'eta_ambi is above the given threshold'

 !--Append the message with parameter space, if passed in
 parspace = ""
 if (present(rho)) write (parspace,'(a, Es10.3)') '; rho = ',rho
 if (present(B  )) write (parspace,'(2a,Es10.3)') trim(parspace),'; B = ',B
 if (present(T  )) write (parspace,'(2a,Es10.3)') trim(parspace),'; T = ',T
 if (present(eta_ohm )) write (parspace,'(2a,Es10.3)') trim(parspace),'; eta_ohm = ',eta_ohm
 if (present(eta_hall)) write (parspace,'(2a,Es10.3)') trim(parspace),'; eta_hall = ',eta_hall
 if (present(eta_ambi)) write (parspace,'(2a,Es10.3)') trim(parspace),'; eta_ambi = ',eta_ambi

 !--Print the error messages
 do i = 1,n_fatal
    if (ierrlist(i) > 0) then
       write(iprint, '(3a)') ferr,trim(errmsg(i)),trim(parspace)
       if (iprintw/=iprint) write(iprintw,'(3a)') ferr,trim(errmsg(i)),trim(parspace)
    endif
 enddo
 if (.not. fatal_only) then
    print_n_fdg = .false.
    do i = n_fatal + 1,n_warn
       if (ierrlist(i) < 0) then
          print_n_fdg = .true.
          if (ierrlist(i)==-1 .or. parspace/="") then
             write(iprintw,'(3a)') werr,trim(errmsg(i)),trim(parspace)
          else
             write(iprintw,'(2a,I8,2a)') werr,'Repeated ',-ierrlist(i),' times: ',trim(errmsg(i))
          endif
       endif
    enddo
    !print_n_fdg = .false.  ! Comment this line for additional verbose information for debugging
    if (print_n_fdg) then
       if (present(nden_save) .and. present(fdg_in)) then
          write(iprintw,*) trim(parspace),' n = ',nden_save, ' fdg = ',fdg_in
       elseif (present(nden_save)) then
          write(iprintw,*) trim(parspace),' n = ',nden_save
       elseif (present(fdg_in)) then
          write(iprintw,*) trim(parspace),' fdg = ',fdg_in
       endif
    endif
 endif

end subroutine nicil_translate_error
!======================================================================!
! PRIMARY ROUTINES                                                     !
!======================================================================!
!----------------------------------------------------------------------!
!+
! This is a primary control routine for NICIL to be called by the host
! code.  This single routine will update the ion number densities
! (which is an iterative process) and update the non-ideal coefficients
! (direct calculations).  These processes can be called individually
! or together, where the former might be required depending on where
! properties are updated in the host code.
! Options:
! icall = 0: do everything
! icall = 1: calculate nden_save only (Bfield not required)
! icall = 2: calculate eta only
!+
!----------------------------------------------------------------------!
pure subroutine nicil_update_nimhd(icall,eta_ohm,eta_hall,eta_ambi,Bfield,rho,T,nden_save,ierrlist,data_out,itry,fdg_in)
 integer,          intent(in)    :: icall
 integer,          intent(inout) :: ierrlist(n_warn)
 real,             intent(out)   :: eta_ohm,eta_hall,eta_ambi
 real,             intent(in)    :: Bfield,rho,T
 real,             intent(inout) :: nden_save(:)
 real,   optional, intent(in)    :: fdg_in(:)
 real,   optional, intent(out)   :: data_out(n_data_out)
 integer,optional, intent(out)   :: itry
 integer                         :: j,jj,id2,itry_n0
 real                            :: mass_neutral_mp,nden_electronR,nden_electronT,nden_total,rho_gas,rho_n
 real                            :: zeta,fdg_local
 real                            :: sigmas(8),n_g_tot(na),fdg_inClean(na)
 real                            :: nden_neutral(nspecies)
 real                            :: nden_ion(nspecies),nden_ionT(nspecies)

 itry_n0 = 0
 if (.not.eta_constant) then
    eta_ohm  = 0.0
    eta_hall = 0.0
    eta_ambi = 0.0
    fdg_inClean = 0.0
    if (T > small2 .and. (Bfield > small2 .or. icall==1)) then
       !--Determine actual gas density (i.e. remove the dust component)
       if (.not. use_fdg_in) then
          rho_gas = rho*one_minus_fdg          ! remove dust component using the global dust-to-gas ratio
       else
          if (present(fdg_in)) then
             fdg_inClean = fdg_in
             do j = 1,na
                if (0.0 < fdg_inClean(j) .and. fdg_inClean(j) < MPval) fdg_inClean(j) = 0.  ! reset small, non-zero dust fractions for computational stability
             enddo
             fdg_local = sum(fdg_inClean(:))
             if (fdg_local < fdg_max) then
                rho_gas = rho*(1.0-fdg_local)  ! remove the total dust component from the input ratios
             else
                ! input density is almost all dust & there is not enough gas to perform stable calculations
                if (present(data_out)) data_out = 0.
                ierrlist(ierr_highfdg) = ierrlist(ierr_highfdg) - 1
                nden_save = 0.0
                return
             endif
          else
             !  dust fraction must be passed in if use_fdg_in=true; return and pass fatal warning to parent code
             ierrlist(ierr_fdg_in) = ierrlist(ierr_fdg_in) + 1
             return
          endif
       endif

       !--Calculate the neutral number densities; use the default neutral mass
       nden_total   = rho_gas*mass_neutral_cold1
       nden_neutral = abundance*nden_total

       !--Determine the fractions of molecular and atomic Hydrogen & update the neutral mass
       call nicil_get_HH2_ratio(nden_neutral,T,mass_neutral_mp)

       !--Calculate number densities from thermal ionisation
       !  Always need to do this since icall=1: calculate nden_electronT & icall=2: calculate nden_ionT
       nden_electronT = nden_save(iire+1)
       call nicil_ionT_get_ne(nden_neutral,nden_ionT,nden_electronT,T,ierrlist)
       nden_save(iire+1) = nden_electronT

       if (icall == 0 .or. icall == 1) then
          !--Calculate number densities from cosmic rays (this process is performed in CGS units for stability)
          if (.not.zeta_of_rho) then
             zeta = zeta_cgs
          else
             zeta = zeta_cgs*exp(-zetaCRcoef*sqrt(T*rho_gas/mass_neutral_mp)) + zeta_cgs_min
          endif
          !--Calculate the amount of dust that is burned off
          if (.not.use_fdg_in) then
             n_g_tot = n_grain_coef*rho ! note that fdg is included in n_grain_coef
             if (T > 1700.) then
                n_g_tot = 0.
             else if (T > 725.) then
                n_g_tot = n_g_tot*fraction_unburned(T)
             endif
          else
             n_g_tot = n_grain_coef*rho*fdg_inClean  ! dust burning should be accounted for in the parent code
          endif
          !--inputs are CGS units; outputs are code units
          call nicil_ionR_get_n(nden_save(1:iire),nden_electronR,nden_neutral*unit_ndensity,n_g_tot*unit_ndensity, &
                                T,zeta,itry_n0,ierrlist)
       else
          nden_electronR = nicil_ionR_get_ne(nden_save(1:iire))  ! in this case, need to calculate electron density from ions
          n_g_tot        = 0. ! to prevent compiler warnings
          zeta = 0. ! prevent compiler warnings
       endif

       !--Sum the ion populations from thermal and cosmic ray ionisation
       nden_ion(ie)          = nden_electronR    + nden_electronT
       nden_ion(1    :iire ) = nden_save(1:iire) + nden_ionT(1    :iire )
       nden_ion(iitps:iitpe) =                     nden_ionT(iitps:iitpe)

       if (icall == 0 .or. icall == 2) then
          !--Calculate the conductivities
          call nicil_ion_get_sigma(Bfield,rho_gas,nden_ion,mass_neutral_mp,T,sigmas,rho_n,ierrlist)

          !--Calculate the coefficients
          call nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,sigmas)
       else
          rho_n = 0. ! to prevent compiler warnings
       endif

       !--Copy important data to an output array, if requested
       if (present(data_out)) then
          data_out       = 0.
          ! This array is only complete if icall==0
          if (icall == 0) then
             data_out(   1) = sigmas(1)                     ! Ohmic conductivities
             data_out(   2) = sigmas(5)                     ! Hall conductivities
             data_out(   3) = sigmas(2)                     ! Pedersen conductivities
             data_out(   4) = rho_n                         ! neutral mass density of the gas
             data_out(   5) = nden_total                    ! number density of the gas
             data_out(   6) = nden_neutral(iH2)             ! neutral molecular hydrogen number density
             data_out(   7) = nden_neutral(iH)              ! neutral atomic hydrogen density
             data_out(   8) = nden_electronR+nden_electronT ! electron number density
             data_out(9:21) = nden_ion(iirs:iitpe)          ! ion number densities
             do j = 1,na                                    ! positive, neutral and negative grain number densities
                id2 = 22+3*(j-1)
                jj  = 2*(j-1)
                data_out(id2  ) = nden_ion(iGp+jj)
                data_out(id2+1) = n_g_tot(j)-nden_ion(iGp+jj)-nden_ion(iGn+jj)
                data_out(id2+2) = nden_ion(iGn+jj)
             enddo
             data_out(22+3*na) = zeta                       ! cosmic ray ionisation rate (only useful for density dependent zeta)
          endif
       endif
       !--warning messages of eta is too big (likely indicates a bug; will = 0 if icall=1 thus not trigger)
       if (eta_ohm      > eta_thresh) ierrlist(ierr_or) = ierrlist(ierr_or) - 1
       if (abs(eta_ohm) > eta_thresh) ierrlist(ierr_he) = ierrlist(ierr_he) - 1
       if (eta_ambi     > eta_thresh) ierrlist(ierr_ad) = ierrlist(ierr_ad) - 1
    else
       !--Exit with error message if T = 0 or B = 0
       if ( T      <=small2 ) ierrlist(ierr_T) = ierrlist(ierr_T) - 1
       if ( Bfield <=small2 .and. icall/=1) ierrlist(ierr_B) = ierrlist(ierr_B) - 1
       if (present(data_out)) data_out = 0.
       nden_save = 0.0
    endif
 else
    !--Return constant coefficient version and exit
    call nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
    if (present(data_out)) data_out = 0.
 endif
 if (present(itry)) itry = itry_n0

end subroutine nicil_update_nimhd
!----------------------------------------------------------------------!
!+
!  Calculates the condictivities
!  Terms have been modified such that they require input of sound speed
!  rather than temperature.
!  sigmas(1): Ohmic conductivity: sigma_O
!  sigmas(2): Pedersen conductivity: sigma_P
!  sigmas(3): positive component of Hall conductivity: sigma_H > 0
!  sigmas(4): negative component of Hall conductivity: sigma_H < 0
!  sigmas(5): total Hall conductivity: sigma_H
!  sigmas(6): inverse Ohmic conductivity: 1/sigma_O
!  sigmas(7): inverse square of the perpendicular conductivity:
!             1.0/(sigma_P**2 + sigma_H**2)
!  sigmas(8): rearrangement for the coefficient to Ambipolar:
!             sigma_O*sigma_P - sigma_P*sigma_P - sigma_H*sigma_H
!  if sigmas(8) < 0, then sigmas(8) == 0; commented-out code can be used
!  to calculate this using other methods
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ion_get_sigma(Bmag,rho_gas,nden_ion,mass_neutral_mp,T,sigmas,rho_n,ierrlist)
 real,    intent(out)   :: sigmas(8),rho_n
 real,    intent(in)    :: Bmag,rho_gas,T,mass_neutral_mp
 real,    intent(in)    :: nden_ion(:)
 integer, intent(inout) :: ierrlist(:)
 integer                :: i,k,p
 real                   :: sqrtT,logT,sigma_coef_onB,n_electron
 real                   :: rho_ion
 real                   :: nu_jn(nspecies),ns(nspecies),rho_k(nspecies)
 real                   :: betai(nspecies),beta2p11(nspecies),sigmav_local(nspecies)

 !--Initialise values
 sqrtT          = sqrt(T)
 logT           = log10(T)
 sigma_coef_onB = qe_c/Bmag
 sigmas         = 0.0
 ns             = nden_ion
 n_electron     = nden_ion(ie)
 if (n_electron > 0.0) then               ! normalise to prevent (possible) numerical overflow
    ns             = ns/n_electron
    sigma_coef_onB = sigma_coef_onB*n_electron
 endif
 !--Calculate species and neutral mass densities
 rho_k = nden_ion*mass_ion
 rho_n = rho_gas
 do k = iH3p,ie
    rho_n = rho_n - rho_k(k)
 enddo

 if (warn_approx) then
    !--Ensure that the strong coupling approximation is valid
    rho_ion = rho_gas - rho_n + rho_k(ie)
    if (warn_ratio_m1*rho_gas > rho_n .or. rho_n > warn_ratio_p1*rho_gas) ierrlist(ierr_scN) = ierrlist(ierr_scN) - 1
    if (warn_ratio2*rho_ion > rho_gas) ierrlist(ierr_scI) = ierrlist(ierr_scI) - 1
 endif

 if ( rho_n > 0.0 ) then
    !--Temperature-dependent momentum transfer coefficients (CGS: cm^3 s^-1)
    sigmav_local         = sigmav
    sigmav_local(iH3p)   =     1.d-9*      (2.693 - (1.238 - (0.664 - 0.089 *logT)*logT)*logT)
    sigmav_local(iHCOp)  = max(1.d-9*sqrtT*(1.476 - (1.409 - (0.555 - 0.0775*logT)*logT)*logT),0.0d0)
    sigmav_local(iHp)    =     1.d-9*      (1.003 + (0.050 + (0.136 - 0.014 *logT)*logT)*logT)
    sigmav_local(ie)     =     1.d-9*sqrtT*(0.535 + (0.203 - (0.163 - 0.050 *logT)*logT)*logT)

    !--Collisional Frequencies (CGS: s^-1)
    nu_jn = 0.0
    nu_jn = nu_jn_coef*rho_n/( mass_neutral_mp + mass_ion_mp ) * sigmav_local
    if (na > 0) then
       nu_jn(iGp:iGn+2*(na-1)) = nu_jn_coef(iGp:iGn+2*(na-1))/sqrt(mass_neutral_mp)*sqrtT &
                               * rho_n/( mass_neutral_mp + mass_ion_mp(iGp:iGn+2*(na-1)))
    endif

    !--Hall parameters (dimensionless, but calculated in CGS)
    betai    = 0.0
    do i = 1,nspecies
       if (nu_jn(i) > 0.) betai(i) = beta_coef(i)*Bmag/nu_jn(i)
    enddo
    beta2p11 = 1.0/( 1.0+betai**2 )

    !--Conductivities (code units)
    do i = 1,nspecies
       sigmas(1) =    sigmas(1) + ns(i)*aZj(i)*betai(i)              ! sigma_O
       sigmas(2) =    sigmas(2) + ns(i)*aZj(i)*betai(i)*beta2p11(i)  ! sigma_P
       if (Zj(i) > 0.0) then
          sigmas(3) = sigmas(3) + ns(i)* Zj(i)         *beta2p11(i)  ! sigma_H > 0
       else
          sigmas(4) = sigmas(4) + ns(i)* Zj(i)         *beta2p11(i)  ! sigma_H < 0
       endif
    enddo
    sigmas(5) = sigmas(3)+sigmas(4)                                  ! sigma_H
    if ( abs(sigmas(5)) < smallXten*max(sigmas(3),-sigmas(4))) then
       sigmas(5) = 0.0                                               ! sigma_H (alternative)
       ierrlist(ierr_sigH) = ierrlist(ierr_sigH) - 1
    endif
    sigmas = sigmas*sigma_coef_onB                                   ! scale sigma
    if (sigmas(1) > 0.0) sigmas(6) = 1.0/sigmas(1)                   ! 1/sigma_O
    sigmas(7)   = sigmas(2)*sigmas(2) + sigmas(5)*sigmas(5)          ! perp^2 = P^2 + H^2
    sigmas(8)   = sigmas(1)*sigmas(2) - sigmas(7)                    ! OP - perp^2
    if ( sigmas(8) < smallXten*max(sigmas(1)*sigmas(2),sigmas(7))) then
       sigmas(8) = 0.0
       do i = 1,nspecies
          do p = i+1,nspecies
             sigmas(8) = sigmas(8) + ns(i)*aZj(i)*betai(i)*beta2p11(i) &
                                   * ns(p)*aZj(p)*betai(p)*beta2p11(p) &
                                   * (sZj(i)*betai(i) - sZj(p)*betai(p))**2
          enddo
       enddo
       sigmas(8) = sigmas(8)*sigma_coef_onB**2
    endif
    if ( sigmas(7) > 0.0 ) sigmas(7) = 1.0/sigmas(7)                 ! perp^2 -> 1/perp^2
 else
    !--This is the Ideal MHD regime.  Turn off non-ideal terms.
    ierrlist(ierr_alion) = ierrlist(ierr_alion) - 1
    sigmas = -1.0
    rho_n  =  0.0
 endif

end subroutine nicil_ion_get_sigma
!======================================================================!
! GRAIN-RELATED SUBROUTINES                                            !
!======================================================================!
!----------------------------------------------------------------------!
!+
! This calculates the fraction of dust that is not burned off at medium
! temperatures.  This assumes a dust grain is homogeneously comprised
! of 88.3% Carbon, 11.2% Silicates and 0.5% Aluminium Oxide.  At these
! fractions, this equations roughly reproduces plots 5-7 from
! Lenzuni, Gail & Henning (1995).
!+
!----------------------------------------------------------------------!
pure real function fraction_unburned(T)
 real, intent(in) :: T

 if (T > 1670.) then
     fraction_unburned = 0.005 - 0.0001 * ( T - 1670.0 )
 elseif (T > 1305.) then
     fraction_unburned = 0.005
 elseif (T > 1274.) then
     fraction_unburned = 0.117 - 0.0025 * ( T - 1260.0 )
 elseif (T > 1225.) then
     fraction_unburned = 0.117 - 0.0007 * ( T - 1225.0 )
 elseif (T > 1100.) then
     fraction_unburned = 0.117
 elseif (T > 1009.) then
     fraction_unburned = 0.117 - 0.0005 * ( T - 1100.0 )
 elseif (T > 825.) then
     fraction_unburned = 1.000 - 0.004  * ( T -  800.0 )
 elseif (T > 725.) then
     fraction_unburned = 1.000 - 0.001  * ( T -  725.0 )
 else
     fraction_unburned = 1.000
 endif

end function fraction_unburned
!======================================================================!
! COSMIC RAY IONISATION-RELATED SUBROUTINES                            !
!======================================================================!
!----------------------------------------------------------------------!
!+
! Calculate the electron, ion and grain number densities for
! cosmic ray ionisation, using the Newton–Raphson method in multiple
! dimensions.
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionR_get_n(nden_ionR,n_e,nden_neutral,n_g_tot,T,zeta,itry_n0,ierrlist)
 real,    intent(inout) :: nden_ionR(:)
 real,    intent(in)    :: nden_neutral(:),n_g_tot(:)
 real,    intent(in)    :: T,zeta
 real,    intent(out)   :: n_e
 integer, intent(out)   :: itry_n0
 integer, intent(inout) :: ierrlist(:)
 integer, parameter     :: i_try_max  = 16
 integer, parameter     :: i_try_max2 = i_try_max*i_try_max
 integer                :: i,j,jj,iter,i_try_new_n0
 integer                :: i1_try_new_n0,i2_try_new_n0
 integer                :: feqni(neqn)
 real                   :: r_kk(iMgp_e,inn),r_eG(2,na,inn),r_kG(2,na,iire,inn),r_GG(3,na,na,inn),zeta_n(iire)
 real                   :: nH_nHe,nH_nHe_ten,nH_nHe_small,max_n,fac
 real                   :: n_new(neqn),dn(neqn)
 real                   :: feqn(neqn),Jacob(neqn,neqn)
 logical                :: iterate,negative_n,negative_g,lerr,reset_loop

 !--Convert nden_ionR to CGS; nden_neutral & n_g_tot have been converted upon input
 nden_ionR = nden_ionR*unit_ndensity

 !--Determine the rate coefficients given the current temperature & multiply by the neutral number densities
 call nicil_set_reaction_rates(T,zeta,nden_neutral,r_kk,r_eG,r_kG,r_GG,zeta_n)

 !--Set initial conditions
 iter          = 0
 iterate       = .true.
 negative_n    = .false.
 i_try_new_n0  = 0
 i1_try_new_n0 = 0
 i2_try_new_n0 = 1
 n_new         = 0. ! to avoid compiler warnings
 max_n         = 0. ! to avoid compiler warnings
 nH_nHe        = nden_neutral(iH2) + 2.*nden_neutral(iH) + nden_neutral(iHe)
 nH_nHe_ten    = nH_nHe*10.
 nH_nHe_small  = nH_nHe*small2

 do while (iterate)  ! to turn off cosmic ray ionisation, replace iterate with .false.
    !--Calculate electron & neutral grain densities;
    !  multiply the reaction rates with the charged number densities
    call nicil_reaction_rates_X_ncharged(r_kk,r_eG,r_kG,r_GG,nden_ionR,n_g_tot,n_e)
    !--Update the system of equations and the Jacobian
    call nicil_ionR_calc_Jf(feqni,feqn,Jacob,r_kk,r_eG,r_kG,r_GG,zeta_n)
    !--Calculate dn by solving Jacob*dn = feqn
    call get_dn_viaLU(dn,feqni,feqn,Jacob,lerr)
    if (.not.lerr) then
       !--Calculate the new iterated number densities
       n_new = nden_ionR - dn
       !--Determine if we need to iterate again
       iter       = iter + 1
       iterate    = .false.
       negative_n = .false.
       max_n      = 0.
       do i = 1,neqn
          if (abs(n_new(i)) > 0. .and. abs(n_new(i)) < nH_nHe_small) then
             n_new(i) = 0.0
          elseif ( (abs(nden_ionR(i)-n_new(i)) > abs(nden_ionR(i))*NRtol)) then
             iterate = .true.
          endif
          if (abs(n_new(i)) > nH_nHe_ten) lerr       = .true.  ! ion number densities are growing beyond their neutral number density
          if (n_new(i)      < 0.0       ) negative_n = .true.  ! a flag in case we converge to a negative density
          max_n = max(max_n, n_new(i))
       enddo
       nden_ionR = n_new
       !--Failsafe to prevent an infinite loop
       if (iter > NRctrmax) iterate = .false.
    else
       !--Do not iterate if the Jacobian was not solved
       iterate = .false.
    endif

    !--Not that values are converged, verify they are legitimate
    !  Ensure neutral grain densities are positive
    negative_g = .false.
    do j = 1,na
       jj = 2*(j-1)
       if (n_g_tot(j) - nden_ionR(iGp+jj) - nden_ionR(iGn+jj) < 0.) negative_g = .true.
    enddo

    !--If n_e < 0., then see if this value is likely a result of round-off; if so, then reset them to 0 and continue
    if (n_e < 0. .and. abs(n_e) < MPtol*max_n) n_e = small

    !--If negative_n = true, then see if these values are likely a result of round-off; if so, then reset them to 0 and continue
    if (n_e > 0. .and. .not. negative_g .and. negative_n) then
       negative_n = .false.
       do i = 1,neqn
          if (nden_ionR(i) < 0.0 ) then
             if (-nden_ionR(i) < MPtol*max_n) then
                nden_ionR(i) = 0.
             else
                negative_n = .true.
             endif
          endif
       enddo
    endif

    !--Try again if we converged to invalid options
    if (.not. iterate .and. (negative_n .or. negative_g .or. n_e < 0.0)) lerr = .true.

    !--Iterations unsuccessful: Try again with new default guess
    if (iter > NRctrmax .or. lerr) then
       reset_loop         = .true.
       if (i_try_new_n0==0) then
          nden_ionR       = nH_nHe*1.d-12
          nden_ionR(iHp)  = nden_neutral(iHe)
          nden_ionR(iHep) = nden_neutral(iHe)
       elseif (i_try_new_n0  < i_try_max2+1) then
          i1_try_new_n0 = i1_try_new_n0 + 1
          if (i1_try_new_n0 > i_try_max) then
             i1_try_new_n0 = 1
             i2_try_new_n0 = i2_try_new_n0 + 1
          endif
          fac = 10.0**(-i1_try_new_n0)
          nden_ionR = nden_neutral(1:iire)*fac
          fac = 10.0**(-i2_try_new_n0)
          do j = 1,na
             nden_ionR(iGp+2*(j-1)) = n_g_tot(j)*fac
             nden_ionR(iGn+2*(j-1)) = n_g_tot(j)*fac
          enddo
       elseif (i_try_new_n0 == i_try_max2+1) then
          nden_ionR       = nH_nHe
       elseif (i_try_new_n0 == i_try_max2+2) then
          nden_ionR       = nH_nHe*1.0d-12
       elseif (i_try_new_n0 == i_try_max2+3) then
          nden_ionR       = nH_nHe*1.0d-32
       elseif (i_try_new_n0 == i_try_max2+4) then
          nden_ionR       = nden_neutral(1:iire)
          do j = 1,na
             nden_ionR(iGp+2*(j-1)) = n_g_tot(j)
             nden_ionR(iGn+2*(j-1)) = n_g_tot(j)
          enddo
       elseif (i_try_new_n0 == i_try_max2+5) then
          nden_ionR       = nden_neutral(1:iire)*10.
          do j = 1,na
             nden_ionR(iGp+2*(j-1)) = n_g_tot(j)*10.
             nden_ionR(iGn+2*(j-1)) = n_g_tot(j)*10.
          enddo
       else
          !--All initial conditions failed to converge or converged to a negative number
          reset_loop = .false.
          iterate    = .false.
       endif

       i_try_new_n0 = i_try_new_n0 + 1
       if (reset_loop) then
          iterate      = .true.
          iter         = 0
       else
          ! This may reasonably happen at very high temperatures; non-convergence will be acceptable
          ! since the major contribution to eta will be from thermal ionisation
          ierrlist(ierr_nRconv) = ierrlist(ierr_nRconv) - 1
          nden_ionR = 0.
          n_e       = 0.
       endif
    endif
 enddo
 itry_n0 = i_try_new_n0 - 1

 !--Revert nden_ionR & electron number density to code units
 nden_ionR = nden_ionR*unit_ndensity1
 n_e       = n_e      *unit_ndensity1

end subroutine nicil_ionR_get_n
!----------------------------------------------------------------------!
!+
! Calculate the Jacobian and Rate equations
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionR_calc_Jf(feqni,feqn,Jacob,r_kk,r_eG,r_kG,r_GG,zeta_n)
 integer, intent(out) :: feqni(neqn)
 real,    intent(out) :: feqn(neqn),Jacob(neqn,neqn)
 real,    intent(in)  :: r_kk(:,:),r_eG(:,:,:),r_kG(:,:,:,:),r_GG(:,:,:,:),zeta_n(:)
 integer              :: i,j,j1,jj,jj1,k,p,imax,jmax,feqnitmp(neqn),arrayTrans(neqn,neqn)
 integer              :: iGnj,iGpj,iGnj1,iGpj1
 real                 :: Jacobmax
 real                 :: feqntmp(neqn),Jacobtmp(neqn,neqn),Jacobrow(neqn)

 !--Initialise the values
 feqn   = 0.0
 Jacob  = 0.0

 !--Rate equations & corresponding Jacobian for the ions
 feqn( iH3p) = zeta_n(iH3p) -r_kk(iH3p_CO,inn) -r_kk(iH3p_Mg,inn) -r_kk(iH3p_e,inn) -r_kk(ie_H3p,inn)
 Jacob(iH3p,iH3p )     =    -r_kk(iH3p_CO,i_n) -r_kk(iH3p_Mg,i_n) -r_kk(iH3p_e,i_n) -r_kk(ie_H3p,in_)
 Jacob(iH3p,iirs:iire) = Jacob(iH3p,iirs:iire)                    -r_kk(iH3p_e,in_) -r_kk(ie_H3p,i_n)

 feqn( iHCOp) =            r_kk(iH3p_CO, inn) -r_kk(iHCOp_Mg,inn) -r_kk(iHCOp_e,inn)
 Jacob(iHCOp,iH3p )     =  r_kk(iH3p_CO, i_n)
 Jacob(iHCOp,iHCOp)     =                     -r_kk(iHCOp_Mg,i_n) -r_kk(iHCOp_e,i_n)
 Jacob(iHCOp,iirs:iire) = Jacob(iHCOp,iirs:iire)                  -r_kk(iHCOp_e,in_)

 feqn( iHp) =  zeta_n(iHp) -r_kk(iHp_O, inn) -r_kk(iHp_O2,inn) -r_kk(iHp_Mg,inn) +r_kk(iHep_H2,inn)  &
                           -r_kk(iHp_Si,inn) -r_kk(iHp_S, inn) -r_kk(iHp_e, inn)
 Jacob(iHp,iHp  )     =    -r_kk(iHp_O, i_n) -r_kk(iHp_O2,i_n) -r_kk(iHp_Mg,i_n) &
                           -r_kk(iHp_Si,i_n) -r_kk(iHp_S, i_n) -r_kk(iHp_e, i_n)
 Jacob(iHp,iHep )     =                                                           r_kk(iHep_H2,i_n)
 Jacob(iHp,iirs:iire) = Jacob(iHp,iirs:iire)                   -r_kk(iHp_e, in_)

 feqn( iHep) = zeta_n(iHep) -r_kk(iHep_H2,inn) -r_kk(iHep_O2,inn) -r_kk(iHep_CO,inn) -r_kk(iHep_Si,inn) -r_kk(iHep_e,inn)
 Jacob(iHep,iHep )     =    -r_kk(iHep_H2,i_n) -r_kk(iHep_O2,i_n) -r_kk(iHep_CO,i_n) -r_kk(iHep_Si,i_n) -r_kk(iHep_e,i_n)
 Jacob(iHep,iirs:iire) = Jacob(iHep,iirs:iire)                                       -r_kk(iHep_e, in_)

 feqn( iCp) = zeta_n(iCp) +r_kk(iHep_CO,inn) -r_kk(iCp_Si,inn) -r_kk(iCp_S,inn) -r_kk(iCp_e,inn)
 Jacob(iCp,iHep )     =    r_kk(iHep_CO,i_n)
 Jacob(iCp,iCp  )     =                      -r_kk(iCp_Si,i_n) -r_kk(iCp_S,i_n) -r_kk(iCp_e,i_n)
 Jacob(iCp,iirs:iire) = Jacob(iCp,iirs:iire) -r_kk(iCp_e, in_)

 feqn( iOp) = zeta_n(iOp) +r_kk(iHp_O,inn) +r_kk(iHep_O2,inn) -r_kk(iOp_e,inn)
 Jacob(iOp,iHp  )     =    r_kk(iHp_O,i_n)
 Jacob(iOp,iHep )     =                     r_kk(iHep_O2,i_n)
 Jacob(iOp,iOp  )     =                                       -r_kk(iOp_e,i_n)
 Jacob(iOp,iirs:iire) = Jacob(iOp,iirs:iire)                  -r_kk(iOp_e,in_)

 feqn( iO2p) =            r_kk(iHp_O2,inn)     -r_kk(iO2p_Si,inn) -r_kk(iO2p_S,inn) -r_kk(iO2p_e,inn)
 Jacob(iO2p,iHp  )     =  r_kk(iHp_O2,i_n)
 Jacob(iO2p,iO2p )     =                       -r_kk(iO2p_Si,i_n) -r_kk(iO2p_S,i_n) -r_kk(iO2p_e,i_n)
 Jacob(iO2p,iirs:iire) = Jacob(iO2p,iirs:iire) -r_kk(iO2p_e,in_)

 feqn( iSip) =            r_kk(iHp_Si, inn) +r_kk(iHep_Si,inn) +r_kk(iCp_Si,inn) +r_kk(iO2p_Si,inn) &
                         -r_kk(iSip_Mg,inn) +r_kk(iSp_Si, inn) -r_kk(iSip_e,inn)
 Jacob(iSip,iHp )      =  r_kk(iHp_Si, i_n)
 Jacob(iSip,iHep)      =                     r_kk(iHep_Si,i_n)
 Jacob(iSip,iCp )      =                                        r_kk(iCp_Si,i_n)
 Jacob(iSip,iO2p)      =                                                          r_kk(iO2p_Si,i_n)
 Jacob(iSip,iSip)      = -r_kk(iSip_Mg,i_n)                    -r_kk(iSip_e,i_n)
 Jacob(iSip,iirs:iire) = Jacob(iSip,iirs:iire) -r_kk(iSip_e,in_)

 feqn( iSp) =           r_kk(iHp_S,inn)      +r_kk(iCp_S,inn) +r_kk(iO2p_S,inn) -r_kk(iSp_Mg,inn) -r_kk(iSp_Si,inn) &
                                                              -r_kk(iSp_e, inn)
 Jacob(iSp,iHp )      = r_kk(iHp_S,i_n)
 Jacob(iSp,iCp )      =                       r_kk(iCp_S,i_n)
 Jacob(iSp,iO2p)      =                                       r_kk(iO2p_S, i_n)
 Jacob(iSp,iSp )      =                                      -r_kk(iSp_e,  i_n) -r_kk(iSp_Mg,i_n) -r_kk(iSp_Si,i_n)
 Jacob(iSp,iirs:iire) = Jacob(iSp,iirs:iire)                 -r_kk(iSp_e,  in_)

 feqn( iMgp) =            r_kk(iH3p_Mg,inn) +r_kk(iHCOp_Mg,inn) +r_kk(iHp_Mg,inn) +r_kk(iSip_Mg,inn) +r_kk(iSp_Mg,inn) &
                                                                -r_kk(iMgp_e,inn)
 Jacob(iMgp,iH3p )     =  r_kk(iH3p_Mg,i_n)
 Jacob(iMgp,iHCOp)     =                     r_kk(iHCOp_Mg,i_n)
 Jacob(iMgp,iHp  )     =                                         r_kk(iHp_Mg,i_n)
 Jacob(iMgp,iSip )     =                                                           r_kk(iSip_Mg,i_n)
 Jacob(iMgp,iSp  )     =                                                                              r_kk(iSp_Mg,i_n)
 Jacob(iMgp,iMgp )     =                                        -r_kk(iMgp_e,i_n)
 Jacob(iMgp,iirs:iire) = Jacob(iMgp,iirs:iire)                  -r_kk(iMgp_e,in_)

 !--Supplement the terms with the grain contributions & calculate the grain terms
 do j = 1,na
    jj   = 2*(j-1)
    iGpj = iGp+jj
    iGnj = iGn+jj

    ! The ion reactions that involve electrons (since e = sum(ions) + Gp - Gn; naturally accounted for in feqn)
    Jacob(iH3p, iGpj) = Jacob(iH3p, iGpj) - r_kk(iH3p_e, in_) - r_kk(ie_H3p,i_n)
    Jacob(iH3p, iGnj) = Jacob(iH3p, iGnj) + r_kk(iH3p_e, in_) + r_kk(ie_H3p,i_n)
    Jacob(iHCOp,iGpj) = Jacob(iHCOp,iGpj) - r_kk(iHCOp_e,in_)
    Jacob(iHCOp,iGnj) = Jacob(iHCOp,iGnj) + r_kk(iHCOp_e,in_)
    Jacob(iHp,  iGpj) = Jacob(iHp,  iGpj) - r_kk(iHp_e,  in_)
    Jacob(iHp,  iGnj) = Jacob(iHp,  iGnj) + r_kk(iHp_e,  in_)
    Jacob(iHep, iGpj) = Jacob(iHep, iGpj) - r_kk(iHep_e, in_)
    Jacob(iHep, iGnj) = Jacob(iHep, iGnj) + r_kk(iHep_e, in_)
    Jacob(iCp,  iGpj) = Jacob(iCp,  iGpj) - r_kk(iCp_e,  in_)
    Jacob(iCp,  iGnj) = Jacob(iCp,  iGnj) + r_kk(iCp_e,  in_)
    Jacob(iOp,  iGpj) = Jacob(iOp,  iGpj) - r_kk(iOp_e,  in_)
    Jacob(iOp,  iGnj) = Jacob(iOp,  iGnj) + r_kk(iOp_e,  in_)
    Jacob(iO2p, iGpj) = Jacob(iO2p, iGpj) - r_kk(iO2p_e, in_)
    Jacob(iO2p, iGnj) = Jacob(iO2p, iGnj) + r_kk(iO2p_e, in_)
    Jacob(iSip, iGpj) = Jacob(iSip, iGpj) - r_kk(iSip_e, in_)
    Jacob(iSip, iGnj) = Jacob(iSip, iGnj) + r_kk(iSip_e, in_)
    Jacob(iSp,  iGpj) = Jacob(iSp,  iGpj) - r_kk(iSp_e,  in_)
    Jacob(iSp,  iGnj) = Jacob(iSp,  iGnj) + r_kk(iSp_e,  in_)
    Jacob(iMgp, iGpj) = Jacob(iMgp, iGpj) - r_kk(iMgp_e, in_)
    Jacob(iMgp, iGnj) = Jacob(iMgp, iGnj) + r_kk(iMgp_e, in_)

    ! Grain contribution to charged ions
    do k = iirs,iire
       feqn( k)      = feqn(k)       - r_kG(ik_Gn,j,k,inn) - r_kG(ik_G0,j,k,inn)
       Jacob(k,k   ) = Jacob(k,k   ) - r_kG(ik_Gn,j,k,i_n) - r_kG(ik_G0,j,k,i_n)
       Jacob(k,iGpj) = Jacob(k,iGpj)                       + r_kG(ik_G0,j,k,in_)
       Jacob(k,iGnj) = Jacob(k,iGnj) - r_kG(ik_Gn,j,k,in_) + r_kG(ik_G0,j,k,in_)
    enddo
    ! Ion contribution to grains
    do k = iirs,iire
       feqn(iGpj)       = feqn(iGpj)       + r_kG(ik_G0,j,k,inn)
       feqn(iGnj)       = feqn(iGnj)       - r_kG(ik_Gn,j,k,inn)
       Jacob(iGpj,k   ) = Jacob(iGpj,k   ) + r_kG(ik_G0,j,k,i_n)
       Jacob(iGpj,iGpj) = Jacob(iGpj,iGpj) - r_kG(ik_G0,j,k,in_)
       Jacob(iGpj,iGnj) = Jacob(iGpj,iGnj) - r_kG(ik_G0,j,k,in_)
       Jacob(iGnj,k   ) = Jacob(iGnj,k   ) - r_kG(ik_Gn,j,k,i_n)
       Jacob(iGnj,iGnj) = Jacob(iGnj,iGnj) - r_kG(ik_Gn,j,k,in_)
    enddo
    ! Electron contribution to grains
    feqn(iGpj) = feqn(iGpj) - r_eG(ie_Gp,j,inn)
    feqn(iGnj) = feqn(iGnj) + r_eG(ie_G0,j,inn)
    do k = iirs,iire
       Jacob(iGpj,k) = Jacob(iGpj,k   ) - r_eG(ie_Gp,j,i_n)
       Jacob(iGnj,k) = Jacob(iGnj,k   ) + r_eG(ie_G0,j,i_n)
    enddo
    Jacob(iGpj,iGpj) = Jacob(iGpj,iGpj) - r_eG(ie_Gp,j,in_)
    Jacob(iGnj,iGpj) = Jacob(iGnj,iGpj) - r_eG(ie_G0,j,in_)
    Jacob(iGnj,iGnj) = Jacob(iGnj,iGnj) - r_eG(ie_G0,j,in_)
    do j1 = 1,na
       jj1   = 2*(j1-1)
       iGpj1 = iGp+jj1
       iGnj1 = iGn+jj1
       Jacob(iGpj,iGpj1) = Jacob(iGpj,iGpj1) - r_eG(ie_Gp,j,i_n)
       Jacob(iGpj,iGnj1) = Jacob(iGpj,iGnj1) + r_eG(ie_Gp,j,i_n)
       Jacob(iGnj,iGpj1) = Jacob(iGnj,iGpj1) + r_eG(ie_G0,j,i_n)
       Jacob(iGnj,iGnj1) = Jacob(iGnj,iGnj1) - r_eG(ie_G0,j,i_n)
    enddo
    ! Grain-grain interactions
    if (na==1) then
       feqn(iGpj)       = feqn(iGpj)       - r_GG(iGp_Gn,j,j,inn)
       feqn(iGnj)       = feqn(iGnj)       - r_GG(iGp_Gn,j,j,inn)
       Jacob(iGpj,iGpj) = Jacob(iGpj,iGpj) - r_GG(iGp_Gn,j,j,i_n)
       Jacob(iGpj,iGnj) = Jacob(iGpj,iGnj) - r_GG(iGp_Gn,j,j,in_)
       Jacob(iGnj,iGpj) = Jacob(iGnj,iGpj) - r_GG(iGp_Gn,j,j,i_n)
       Jacob(iGnj,iGnj) = Jacob(iGnj,iGnj) - r_GG(iGp_Gn,j,j,in_)
    else
       do j1 = 1,na
          jj1   = 2*(j1-1)
          iGpj1 = iGp+jj1
          iGnj1 = iGn+jj1
          feqn(iGpj)        = feqn(iGpj )       - r_GG(iGp_Gn,j,j1,inn) - r_GG(iGp_G0,j,j1,inn) + r_GG(iGp_G0,j1,j,inn)
          Jacob(iGpj,iGpj ) = Jacob(iGpj,iGpj ) - r_GG(iGp_Gn,j,j1,i_n) - r_GG(iGp_G0,j,j1,i_n) - r_GG(iGp_G0,j1,j,in_)
          Jacob(iGpj,iGnj ) = Jacob(iGpj,iGnj )                                                 - r_GG(iGp_G0,j1,j,in_)
          Jacob(iGpj,iGpj1) = Jacob(iGpj,iGpj1)                         + r_GG(iGp_G0,j,j1,in_) + r_GG(iGp_G0,j1,j,i_n)
          Jacob(iGpj,iGnj1) = Jacob(iGpj,iGnj1) - r_GG(iGp_Gn,j,j1,in_) + r_GG(iGp_G0,j,j1,in_)
          feqn(iGnj)        = feqn(iGnj )       - r_GG(iGp_Gn,j1,j,inn) - r_GG(iGn_G0,j,j1,inn) + r_GG(iGn_G0,j1,j,inn)
          Jacob(iGnj,iGpj ) = Jacob(iGnj,iGpj )                                                 - r_GG(iGn_G0,j1,j,in_)
          Jacob(iGnj,iGnj ) = Jacob(iGnj,iGnj ) - r_GG(iGp_Gn,j1,j,in_) - r_GG(iGn_G0,j,j1,i_n) - r_GG(iGn_G0,j1,j,in_)
          Jacob(iGnj,iGpj1) = Jacob(iGnj,iGpj1) - r_GG(iGp_Gn,j1,j,i_n) + r_GG(iGn_G0,j,j1,in_)
          Jacob(iGnj,iGnj1) = Jacob(iGnj,iGnj1)                         + r_GG(iGn_G0,j,j1,in_) + r_GG(iGn_G0,j1,j,i_n)
       enddo
    endif
 enddo

!--Rearrange the Jacobian so that the largest numbers are on the diagonal (stability uncertain)
 if (reorder_Jacobian) then
    imax = 1
    jmax = 1
    do i = 1,neqn
       feqni(i) = i
    enddo
    do k = 1,neqn-1
       !--Find the largest entry & put it in the k'th row
       Jacobmax = 0.0
       do i = k,neqn
          do j = k,neqn
             if (abs(Jacob(i,j)) > Jacobmax) then
                imax = i
                jmax = j
                Jacobmax = abs(Jacob(i,j))
             endif
          enddo
       enddo
       if (imax > k) then
          Jacobrow      = Jacob(k,:)
          Jacob(k,:)    = Jacob(imax,:)
          Jacob(imax,:) = Jacobrow
          feqntmp       = feqn
          feqn(k)       = feqntmp(imax)
          feqn(imax)    = feqntmp(k)
       endif
       !--Create and use a transformation matrix to put the Jacobmax in the k'th column
       if (jmax > k) then
          arrayTrans = 0
          do i=1,neqn
             arrayTrans(i,i) = 1
          enddo
          arrayTrans(k,   k   ) = 0
          arrayTrans(jmax,jmax) = 0
          arrayTrans(k,   jmax) = 1
          arrayTrans(jmax,k   ) = 1
          Jacobtmp              = Jacob
          Jacob                 = 0.
          feqnitmp              = feqni
          feqni                 = 0
          do i = 1,neqn
             do j = 1,neqn
                do p = 1,neqn
                   Jacob(i,j) = Jacob(i,j) + arrayTrans(j,p)*Jacobtmp(i,p)
                   if (i==1) feqni(j) = feqni(j) + arrayTrans(j,p)*feqnitmp(p)
                enddo
             enddo
          enddo
       endif
    enddo
 endif

end subroutine nicil_ionR_calc_Jf

!----------------------------------------------------------------------!
!+
! Use LU decomposition of the Jacobian to solve dn in Jacob*dn=feqn
!+
!----------------------------------------------------------------------!
pure subroutine get_dn_viaLU(dn,feqni,feqn,Jacob,lerr)
 real,    intent(in)  :: feqn(neqn),Jacob(neqn,neqn)
 integer, intent(in)  :: feqni(neqn)
 real,    intent(out) :: dn(neqn)
 logical, intent(out) :: lerr
 real                 :: LU(neqn,neqn),Ux(neqn)
 integer              :: i,j,k,p

 lerr = .false.
 !--Decompose into Upper & Lower & overwrite the Jacobian at the same time;
 !  this uses the same array since we define diag(L)==1
 LU = Jacob
 do k = 1,neqn                                ! Track the diagonal
    do j = k,neqn                             ! Track across the rows to solve for U
       do p = 1,k-1
          LU(k,j) = LU(k,j) - LU(k,p)*LU(p,j)
       enddo
    enddo
    if (abs(LU(k,k)) < small2 .or. abs(LU(k,k)) > small21 ) lerr = .true.
    do i = k+1,neqn                           ! Track down the columns to solve for L
       do p = 1,k-1
          LU(i,k) = LU(i,k) - LU(i,p)*LU(p,k)
       enddo
       if (.not.lerr) LU(i,k) = LU(i,k)/LU(k,k)
    enddo
 enddo
 if (lerr) then                               ! exit cleanly since dividing by zero
    dn = 0.0
    LU = 0.0
    return
 endif

 !--Solve L*(Ux) = feqn for Ux == y
 Ux = feqn
 do i = 2,neqn
    do j = 1,i-1
       Ux(i) = Ux(i) - LU(i,j)*Ux(j)
    enddo
 enddo

 !--Solve U*dn = y for dn
 dn = Ux
 do i = neqn,1,-1
    do j = i+1,neqn
       dn(i) = dn(i) - LU(i,j)*dn(j)
    enddo
    dn(i) = dn(i)/LU(i,i)
 enddo

 !--Return the dn array to the correct order
 if (reorder_Jacobian) then
    Ux = dn
    do i = 1,neqn
       dn(feqni(i)) = Ux(i)
    enddo
 endif

end subroutine get_dn_viaLU
!----------------------------------------------------------------------!
!+
! Multiply the reaction rates with the ion, electron and grain
! number densities
!+
!----------------------------------------------------------------------!
pure subroutine nicil_reaction_rates_X_ncharged(r_kk,r_eG,r_kG,r_GG,nden_ionR,n_g_tot,n_e)
 real,    intent(inout) :: r_kk(:,:),r_eG(:,:,:),r_kG(:,:,:,:),r_GG(:,:,:,:)
 real,    intent(in)    :: nden_ionR(:),n_g_tot(:)
 real,    intent(out)   :: n_e
 real                   :: n_G0(na)
 integer                :: j,jj,j1,jj1,k,iGpj,iGnj,iGpj1,iGnj1

 !--Calculate the electron and neutral grain densities
 n_e = 0.
 do k = iirs,iire
    n_e = n_e + nden_ionR(k)
 enddo
 do j = 1,na
    jj      = 2*(j-1)
    n_e     = n_e        + nden_ionR(iGp+jj) - nden_ionR(iGn+jj)
    n_G0(j) = n_g_tot(j) - nden_ionR(iGp+jj) - nden_ionR(iGn+jj)
 enddo

 !--Multiply the k-k reactions by ion & electron number densities
 r_kk(iH3p_CO, in_)=r_kk(iH3p_CO, ir)*nden_ionR(iH3p);  r_kk(iH3p_CO, inn)=r_kk(iH3p_CO, i_n)*nden_ionR(iH3p)
 r_kk(iH3p_Mg, in_)=r_kk(iH3p_Mg, ir)*nden_ionR(iH3p);  r_kk(iH3p_Mg, inn)=r_kk(iH3p_Mg, i_n)*nden_ionR(iH3p)
 r_kk(iHCOp_Mg,in_)=r_kk(iHCOp_Mg,ir)*nden_ionR(iHCOp); r_kk(iHCOp_Mg,inn)=r_kk(iHCOp_Mg,i_n)*nden_ionR(iHCOp)
 r_kk(iHp_O,   in_)=r_kk(iHp_O,   ir)*nden_ionR(iHp);   r_kk(iHp_O,   inn)=r_kk(iHp_O,   i_n)*nden_ionR(iHp)
 r_kk(iHp_O2,  in_)=r_kk(iHp_O2,  ir)*nden_ionR(iHp);   r_kk(iHp_O2,  inn)=r_kk(iHp_O2,  i_n)*nden_ionR(iHp)
 r_kk(iHp_Si,  in_)=r_kk(iHp_Si,  ir)*nden_ionR(iHp);   r_kk(iHp_Si,  inn)=r_kk(iHp_Si,  i_n)*nden_ionR(iHp)
 r_kk(iHp_S ,  in_)=r_kk(iHp_S ,  ir)*nden_ionR(iHp);   r_kk(iHp_S ,  inn)=r_kk(iHp_S ,  i_n)*nden_ionR(iHp)
 r_kk(iHp_Mg,  in_)=r_kk(iHp_Mg,  ir)*nden_ionR(iHp);   r_kk(iHp_Mg,  inn)=r_kk(iHp_Mg,  i_n)*nden_ionR(iHp)
 r_kk(iHep_H2, in_)=r_kk(iHep_H2, ir)*nden_ionR(iHep);  r_kk(iHep_H2, inn)=r_kk(iHep_H2, i_n)*nden_ionR(iHep)
 r_kk(iHep_O2, in_)=r_kk(iHep_O2, ir)*nden_ionR(iHep);  r_kk(iHep_O2, inn)=r_kk(iHep_O2, i_n)*nden_ionR(iHep)
 r_kk(iHep_Si, in_)=r_kk(iHep_Si, ir)*nden_ionR(iHep);  r_kk(iHep_Si, inn)=r_kk(iHep_Si, i_n)*nden_ionR(iHep)
 r_kk(iHep_CO, in_)=r_kk(iHep_CO, ir)*nden_ionR(iHep);  r_kk(iHep_CO, inn)=r_kk(iHep_CO, i_n)*nden_ionR(iHep)
 r_kk(iCp_Si,  in_)=r_kk(iCp_Si,  ir)*nden_ionR(iCp);   r_kk(iCp_Si,  inn)=r_kk(iCp_Si,  i_n)*nden_ionR(iCp)
 r_kk(iCp_S,   in_)=r_kk(iCp_S,   ir)*nden_ionR(iCp);   r_kk(iCp_S,   inn)=r_kk(iCp_S,   i_n)*nden_ionR(iCp)
 r_kk(iO2p_Si, in_)=r_kk(iO2p_Si, ir)*nden_ionR(iO2p);  r_kk(iO2p_Si, inn)=r_kk(iO2p_Si, i_n)*nden_ionR(iO2p)
 r_kk(iO2p_S , in_)=r_kk(iO2p_S , ir)*nden_ionR(iO2p);  r_kk(iO2p_S , inn)=r_kk(iO2p_S , i_n)*nden_ionR(iO2p)
 r_kk(iSip_Mg, in_)=r_kk(iSip_Mg, ir)*nden_ionR(iSip);  r_kk(iSip_Mg, inn)=r_kk(iSip_Mg, i_n)*nden_ionR(iSip)
 r_kk(iSp_Si,  in_)=r_kk(iSp_Si,  ir)*nden_ionR(iSp);   r_kk(iSp_Si,  inn)=r_kk(iSp_Si,  i_n)*nden_ionR(iSp)
 r_kk(iSp_Mg,  in_)= r_kk(iSp_Mg, ir)*nden_ionR(iSp);   r_kk(iSp_Mg,  inn)=r_kk(iSp_Mg,  i_n)*nden_ionR(iSp)

 r_kk(iH3p_e,  in_)=r_kk(iH3p_e,  ir)*nden_ionR(iH3p)
 r_kk(iH3p_e,  i_n)=r_kk(iH3p_e,  ir)*n_e;              r_kk(iH3p_e,  inn)=r_kk(iH3p_e, in_)*n_e
 r_kk(ie_H3p,  i_n)=r_kk(ie_H3p,  ir)*nden_ionR(iH3p)
 r_kk(ie_H3p,  in_)=r_kk(ie_H3p,  ir)*n_e;              r_kk(ie_H3p,  inn)=r_kk(ie_H3p, i_n)*n_e
 r_kk(iHCOp_e, in_)=r_kk(iHCOp_e, ir)*nden_ionR(iHCOp)
 r_kk(iHCOp_e, i_n)=r_kk(iHCOp_e, ir)*n_e;              r_kk(iHCOp_e, inn)=r_kk(iHCOp_e,in_)*n_e
 r_kk(iO2p_e,  in_)=r_kk(iO2p_e,  ir)*nden_ionR(iO2p)
 r_kk(iO2p_e,  i_n)=r_kk(iO2p_e,  ir)*n_e;              r_kk(iO2p_e,  inn)=r_kk(iO2p_e, in_)*n_e
 r_kk(iHp_e,   in_)=r_kk(iHp_e,   ir)*nden_ionR(iHp)
 r_kk(iHp_e,   i_n)=r_kk(iHp_e,   ir)*n_e;              r_kk(iHp_e,   inn)=r_kk(iHp_e,  in_)*n_e
 r_kk(iHep_e,  in_)=r_kk(iHep_e,  ir)*nden_ionR(iHep)
 r_kk(iHep_e,  i_n)=r_kk(iHep_e,  ir)*n_e;              r_kk(iHep_e,  inn)=r_kk(iHep_e, in_)*n_e
 r_kk(iCp_e,   in_)=r_kk(iCp_e,   ir)*nden_ionR(iCp)
 r_kk(iCp_e,   i_n)=r_kk(iCp_e,   ir)*n_e;              r_kk(iCp_e,   inn)=r_kk(iCp_e,  in_)*n_e
 r_kk(iOp_e,   in_)=r_kk(iOp_e,   ir)*nden_ionR(iOp)
 r_kk(iOp_e,   i_n)=r_kk(iOp_e,   ir)*n_e;              r_kk(iOp_e,   inn)=r_kk(iOp_e,  in_)*n_e
 r_kk(iSip_e,  in_)=r_kk(iSip_e,  ir)*nden_ionR(iSip)
 r_kk(iSip_e,  i_n)=r_kk(iSip_e,  ir)*n_e;              r_kk(iSip_e,  inn)=r_kk(iSip_e, in_)*n_e
 r_kk(iSp_e,   in_)=r_kk(iSp_e,   ir)*nden_ionR(iSp)
 r_kk(iSp_e,   i_n)=r_kk(iSp_e,   ir)*n_e;              r_kk(iSp_e,   inn)=r_kk(iSp_e,  in_)*n_e
 r_kk(iMgp_e,  in_)=r_kk(iMgp_e,  ir)*nden_ionR(iMgp)
 r_kk(iMgp_e,  i_n)=r_kk(iMgp_e,  ir)*n_e;              r_kk(iMgp_e,  inn)=r_kk(iMgp_e, in_)*n_e

 !--Multiply the grain reactions by the relevant number densities
 do j = 1,na
    jj   = 2*(j-1)
    iGpj = iGp+jj
    iGnj = iGn+jj
    ! Grain-grain reactions
    do j1 = 1,na
       jj1   = 2*(j1-1)
       iGpj1 = iGp+jj1
       iGnj1 = iGn+jj1
       r_GG(iGp_Gn,j,j1,in_) = r_GG(iGp_Gn,j,j1,ir )*nden_ionR(iGpj )
       r_GG(iGp_Gn,j,j1,i_n) = r_GG(iGp_Gn,j,j1,ir )*nden_ionR(iGnj1)
       r_GG(iGp_Gn,j,j1,inn) = r_GG(iGp_Gn,j,j1,in_)*nden_ionR(iGnj1)

       r_GG(iGp_G0,j,j1,in_) = r_GG(iGp_G0,j,j1,ir )*nden_ionR(iGpj )
       r_GG(iGp_G0,j,j1,i_n) = r_GG(iGp_G0,j,j1,ir )*n_G0(j1)
       r_GG(iGp_G0,j,j1,inn) = r_GG(iGp_G0,j,j1,in_)*n_G0(j1)

       r_GG(iGn_G0,j,j1,in_) = r_GG(iGn_G0,j,j1,ir )*nden_ionR(iGnj )
       r_GG(iGn_G0,j,j1,i_n) = r_GG(iGn_G0,j,j1,ir )*n_G0(j1)
       r_GG(iGn_G0,j,j1,inn) = r_GG(iGn_G0,j,j1,in_)*n_G0(j1)
    enddo
    ! Grain-charged ion reactions
    do k = iirs,iire
       r_kG(ik_G0,j,k,in_) = r_kG(ik_G0,j,k,ir )*nden_ionR(k)
       r_kG(ik_G0,j,k,i_n) = r_kG(ik_G0,j,k,ir )*n_G0(j)
       r_kG(ik_G0,j,k,inn) = r_kG(ik_G0,j,k,in_)*n_G0(j)

       r_kG(ik_Gn,j,k,in_) = r_kG(ik_Gn,j,k,ir )*nden_ionR(k)
       r_kG(ik_Gn,j,k,i_n) = r_kG(ik_Gn,j,k,ir )*nden_ionR(iGnj)
       r_kG(ik_Gn,j,k,inn) = r_kG(ik_Gn,j,k,in_)*nden_ionR(iGnj)
    enddo
    ! Grain-electron reactions
    r_eG(ie_G0,j,in_) = r_eG(ie_G0,j,ir )*n_e
    r_eG(ie_G0,j,i_n) = r_eG(ie_G0,j,ir )*n_G0(j)
    r_eG(ie_G0,j,inn) = r_eG(ie_G0,j,in_)*n_G0(j)

    r_eG(ie_Gp,j,in_) = r_eG(ie_Gp,j,ir )*n_e
    r_eG(ie_Gp,j,i_n) = r_eG(ie_Gp,j,ir )*nden_ionR(iGpj)
    r_eG(ie_Gp,j,inn) = r_eG(ie_Gp,j,in_)*nden_ionR(iGpj)
 enddo

end subroutine nicil_reaction_rates_X_ncharged
!----------------------------------------------------------------------!
!+
! Calculate the electron number density, assuming charge neutrality
!+
!----------------------------------------------------------------------!
pure real function nicil_ionR_get_ne(ndens_ionR)
 real, intent(in) :: ndens_ionR(:)
 integer          :: j,k

 nicil_ionR_get_ne = 0.0
 do k = iirs,iire
    nicil_ionR_get_ne = nicil_ionR_get_ne + ndens_ionR(k)
 enddo
 do j = 1,na
    nicil_ionR_get_ne = nicil_ionR_get_ne + ndens_ionR(iGp+2*(j-1))
    nicil_ionR_get_ne = nicil_ionR_get_ne - ndens_ionR(iGn+2*(j-1))
 enddo

end function nicil_ionR_get_ne
!======================================================================!
! THERMAL IONISATION-RELATED SUBROUTINES                               !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  This will solve the Saha equation to determine the electron number
!  density.  The species of interest and thier properties are given
!  at the beginnig of this module.
!  Each element can only be singly ionised.
!  Once this is complete, the neutral number density is recalculate to
!  account for the ions.
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_ne(nden_neutral,nden_ionT,nden_electronT,T,ierrlist)
 real,    intent(in)    :: T
 real,    intent(out)   :: nden_ionT(:)
 real,    intent(inout) :: nden_neutral(:),nden_electronT
 integer, intent(inout) :: ierrlist(:)
 integer                :: k,iter
 real                   :: fatn,fatndn,ne_rat,ne_old
 real                   :: dne,nden_i
 real                   :: Kjk(nspecies)
 logical                :: iterate,will_thermalise,try_new_ne0

 !--Initialise values
 iter        = 0
 Kjk         = 0.
 ne_rat      = NRtol*2.0
 ne_old      = nden_electronT
 nden_ionT   = 0.
 iterate     = .true.
 try_new_ne0 = .true.
 will_thermalise  = .false.

 !--Calculate the coefficients
 do k = iits,iite
    if (T > Texp_thresh(k) .and. chij(k) > small2) then
       Kjk(k)          = Saha_coef(k)*sqrt(T)**3*exp(-chij(k)/T)
       will_thermalise = .true.  ! comment out this line to turn off thermal ionisation
    endif
 enddo

 if (will_thermalise) then
    do while ( iterate )
       call nicil_ionT_get_nion(nden_neutral,nden_ionT,ne_old,nden_electronT,dne,Kjk)
       fatn   = ne_old - nden_electronT
       if ( abs(fatn) < smallXten*max(ne_old,nden_electronT)) fatn = 0.
       fatndn = 1.0    - dne
       nden_electronT = ne_old - fatn/fatndn
       if (nden_electronT > 0.0) then
          ne_rat = abs( 1.0 - ne_old/nden_electronT )
       else
          ne_rat = 0.0
       endif
       !--Actions if converged
       if (ne_rat < NRtol) iterate = .false.

       !--Actions if errors occurred
       if (iterate .and. (iter > NRctrmax .or. nden_electronT < 0.0)) then
          if (try_new_ne0) then
             !--Try again with default guess rather than value from previous iteration
             try_new_ne0    = .false.
             nden_electronT = nden_neutral(iH) + nden_neutral(iH2) + nden_neutral(iHe)
             iter           = 0
             ne_rat         = NRtol*2.0
          else
             !--New guess failed; trigger fatal warnings
             iterate = .false.
             if (iter >       NRctrmax) ierrlist(ierr_neTconv) = ierrlist(ierr_neTconv) + 1  ! n_electronT did not converge
             if (nden_electronT < 0.0 ) ierrlist(ierr_neTle0 ) = ierrlist(ierr_neTle0 ) + 1  ! n_electronT < 0
          endif
       endif
       iter   = iter + 1
       ne_old = nden_electronT
    enddo

    !--Update the neutral number densities to remove the thermally ionised ions
    do k = iits,iite
       nden_i = max(0.,nden_neutral(k) - nden_ionT(k))
       if (nden_i < nden_neutral(k)*small) nden_i = 0.
       nden_neutral(k) = nden_i
    enddo
 else
    !--Zero values and exit
    nden_ionT      = 0.
    nden_electronT = 0.
 endif

end subroutine nicil_ionT_get_ne
!----------------------------------------------------------------------!
!+
!  This calculates the ion number density from thermal ionisation
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_nion(nden_neutral,nden_ionT,ne_old,ne,dne,Kjk)
 real,    intent(out) :: ne,dne,nden_ionT(:)
 real,    intent(in)  :: ne_old,nden_neutral(:),Kjk(:)
 integer              :: k
 real                 :: term,term1

 nden_ionT = 0.0
 ne        = 0.0
 dne       = 0.0
 do k = iits,iite
    term = ne_old + Kjk(k)
    if (term > 0.0) then
       term1 = 1.0/term
    else
       term1 = 0.0
    endif
    if (Kjk(k) > 0.0) nden_ionT(k) = nden_neutral(k)*Kjk(k)*term1
    ne  = ne  + nden_ionT(k)
    dne = dne - nden_ionT(k)*term1
 enddo

end subroutine nicil_ionT_get_nion
!----------------------------------------------------------------------!
!+
!  This calculates the number density of H and H2
!+
!----------------------------------------------------------------------!
pure subroutine nicil_get_HH2_ratio(nden_neutral,T,mass_neutral_mp)
 real,    intent(inout) :: nden_neutral(:)
 real,    intent(in)    :: T
 real,    intent(out)   :: mass_neutral_mp
 integer                :: j
 real                   :: KHH2,nH2,frac,nden_neutral_sum

 nH2  = nden_neutral(iH2)
 frac = 0.
 KHH2 = 0.

 if ( T > Texp_thresh0*chij_HH2 ) then
    KHH2 = Saha_coef_HH2*T**1.5 * exp(-chij_HH2/T)
    if (8.0d7*nH2 < KHH2) then
       frac = 1.0
    else
       frac = ( sqrt(KHH2 * (KHH2 + 4.0*nH2)) - KHH2 )/(2.0*nH2)
    endif
 endif
 nden_neutral(iH2) = (1.0-frac)*nH2
 nden_neutral(iH ) =  2.0*frac *nH2

 !--Update the neutral mass using the correct ratios of H & H2
 nden_neutral_sum = 0.0
 mass_neutral_mp  = 0.0
 do j = 1,nspecies
    nden_neutral_sum = nden_neutral_sum + nden_neutral(j)
    mass_neutral_mp  = mass_neutral_mp  + nden_neutral(j)*mass_neu_mp(j)
 enddo
 mass_neutral_mp = mass_neutral_mp/nden_neutral_sum

end subroutine nicil_get_HH2_ratio
!======================================================================!
! NON-IDEAL MHD-RELATED SUBROUTINES                                    !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  Calculates the coefficients for the non-ideal MHD terms:
!    Ohmic Resistivity, Hall effect, Ambipolar diffusion
!+
!----------------------------------------------------------------------!
pure subroutine nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,sigmas)
 real, intent(out) :: eta_ohm,eta_hall,eta_ambi
 real, intent(in)  :: sigmas(8)

 eta_ohm  = 0.0
 eta_hall = 0.0
 eta_ambi = 0.0
 if (sigmas(6) > 0.0) then ! a verification that we have not entered the ideal regime
    if (use_ohm)  eta_ohm  = csqbyfourpi * sigmas(6)
    if (use_hall) eta_hall = csqbyfourpi * sigmas(5) * sigmas(7)
    if (use_ambi) eta_ambi = csqbyfourpi * sigmas(8) * sigmas(6) * sigmas(7)
 endif

end subroutine nicil_nimhd_get_eta
!-----------------------------------------------------------------------
pure subroutine nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
 real, intent(out) :: eta_ohm,eta_hall,eta_ambi
 real, intent(in)  :: Bfield,rho

 !--Set the coefficient to the constant coefficient (==0 if the term is off)
 eta_ohm  = eta_ohm_cnst
 eta_hall = eta_hall_cnst
 eta_ambi = eta_ambi_cnst
 if (eta_const_type==icnst) return
 !--Multiply by the variable, if requested
 eta_hall = eta_hall*Bfield
 eta_ambi = eta_ambi*Bfield**2/rho**alpha_AD_p1

end subroutine nicil_nimhd_get_eta_cnst
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD contributions to energy
!  Note: dudthall==0
!  Note: the sign of (36) in Wurster (2016) is incorrect
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_dudt_nimhd(dudtnonideal,eta_ohm,eta_ambi,rho,J,B)
 real, intent(out) :: dudtnonideal
 real, intent(in)  :: J(3),B(3)
 real, intent(in)  :: eta_ohm,eta_ambi,rho
 real              :: B2i,J2i,BJi,BJBJihat

 B2i          = dot_product(B,B)
 J2i          = dot_product(J,J)
 if (B2i > 0.0) then
    BJi       = dot_product(B,J)
    BJBJihat  = BJi*BJi/B2i
 else
    BJBJihat  = 0.0
 endif
 dudtnonideal = ( eta_ohm*J2i + eta_ambi*(J2i - BJBJihat) )/rho

end subroutine nicil_get_dudt_nimhd
!-----------------------------------------------------------------------
!+
!  Calculates the timesteps for non-ideal MHD
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_dt_nimhd(dtohm,dthall,dtambi,h,eta_ohm,eta_hall,eta_ambi)
 real, intent(in)  :: h,eta_ohm,eta_hall,eta_ambi
 real, intent(out) :: dtohm,dthall,dtambi
 real              :: h2

 dtohm  = huge(dtohm)
 dthall = huge(dthall)
 dtambi = huge(dtambi)
 h2     = h*h

 if (use_ohm  .and.     eta_ohm   > tiny(eta_ohm ) ) dtohm  =      Cdt_diff*h2/eta_ohm
 if (use_ambi .and.     eta_ambi  > tiny(eta_ambi) ) dtambi =      Cdt_diff*h2/eta_ambi
 if (use_hall .and. abs(eta_hall) > tiny(eta_hall) ) dthall = abs( Cdt_hall*h2/eta_hall )

end subroutine nicil_get_dt_nimhd
!-----------------------------------------------------------------------
!+
!  Calculates the Hall drift velocity (electron-ion drift)
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_halldrift(eta_hall,Bx,By,Bz,jcurrent,vdrift)
 real,    intent(in)    :: eta_hall,Bx,By,Bz,jcurrent(3)
 real,    intent(out)   :: vdrift(3)
 real                   :: B1,B2

 B2 = Bx*Bx + By*By + Bz*Bz
 if (B2 > 0.) then
    B1 = 1.0/sqrt(B2)
 else
    B1 = 0.
 endif
 vdrift = -eta_hall*jcurrent*B1

end subroutine nicil_get_halldrift
!-----------------------------------------------------------------------
!+
!  Calculates the ion-neutral drift velocity caused by ambipolar diffusion
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_ambidrift(eta_ambi,Bx,By,Bz,jcurrent,vdrift)
 real,    intent(in)    :: eta_ambi,Bx,By,Bz,jcurrent(3)
 real,    intent(out)   :: vdrift(3)
 real                   :: B2,B21

 B2 = Bx*Bx + By*By + Bz*Bz
 if (B2 > 0.) then
    B21 = 1.0/B2
 else
    B21 = 0.
 endif
 vdrift(1) = eta_ambi*( jcurrent(2)*Bz - jcurrent(3)*By )*B21
 vdrift(2) = eta_ambi*( jcurrent(3)*Bx - jcurrent(1)*Bz )*B21
 vdrift(3) = eta_ambi*( jcurrent(1)*By - jcurrent(2)*Bx )*B21

end subroutine nicil_get_ambidrift
!-----------------------------------------------------------------------
!+
!  Calculates the ion velocity, assuming v_electron = v_gas
!  v_ion = v_gas + v_drift_Hall + v_drift_ambi
!  (signs verified in Zhao+2020 & MacLow+1995)
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_vion(eta_hall,eta_ambi,vx,vy,vz,Bx,By,Bz,jcurrent,vion,ierrlist,vdrift_out)
 integer,           intent(inout) :: ierrlist(:)
 real,              intent(in)    :: eta_hall,eta_ambi,vx,vy,vz,Bx,By,Bz,jcurrent(3)
 real,              intent(out)   :: vion(3)
 real,    optional, intent(out)   :: vdrift_out(3)
 real                             :: vion2,v2
 real                             :: vdriftH(3),vdriftA(3)

 !--The drift velocities
 call nicil_get_halldrift(eta_hall,Bx,By,Bz,jcurrent,vdriftH)
 call nicil_get_ambidrift(eta_ambi,Bx,By,Bz,jcurrent,vdriftA)

 !--The ion velocity
 vion(1) = vx + vdriftH(1) + vdriftA(1)
 vion(2) = vy + vdriftH(2) + vdriftA(2)
 vion(3) = vz + vdriftH(3) + vdriftA(3)

 if (warn_approx) then
    !--Ensure that the drift velocity is small
    vion2 = vion(1)*vion(1) + vion(2)*vion(2) + vion(3)*vion(3)
    v2    = vx*vx + vy*vy + vz*vz
    if (warn_ratio_m12*vion2 > v2 .or. v2 > warn_ratio_p12*vion2) ierrlist(ierr_drift) = ierrlist(ierr_drift) - 1
    if (vrms2_min         > vion2 .or. vion2 > vrms2_max        ) ierrlist(ierr_vrms ) = ierrlist(ierr_vrms ) - 1
 endif

 if (present(vdrift_out)) vdrift_out = vdriftH + vdriftA

end subroutine nicil_get_vion
!-----------------------------------------------------------------------
!+
!  Calculates JxB and (JxB)xB
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_jcbcb(jcbcb,jcb,jcurrent,Bx,By,Bz,B1)
 real, intent(out) :: jcb(3), jcbcb(3)
 real, intent(in)  :: jcurrent(3)
 real, intent(in)  :: Bx, By, Bz, B1

 jcb(1)   = ( jcurrent(2)*Bz - jcurrent(3)*By )*B1
 jcb(2)   = ( jcurrent(3)*Bx - jcurrent(1)*Bz )*B1
 jcb(3)   = ( jcurrent(1)*By - jcurrent(2)*Bx )*B1

 jcbcb(1) = ( jcb(2)*Bz      - jcb(3)*By      )*B1
 jcbcb(2) = ( jcb(3)*Bx      - jcb(1)*Bz      )*B1
 jcbcb(3) = ( jcb(1)*By      - jcb(2)*Bx      )*B1

end subroutine nimhd_get_jcbcb
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD contributions to the magnetic field in SPH
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dBdt(dBnonideal,eta_ohm,eta_hall,eta_ambi,jcurrent,jcb,jcbcb,dxr1,dyr1,dzr1)
 real, intent(out) :: dBnonideal(3)
 real, intent(in)  :: jcurrent(3),jcb(3),jcbcb(3)
 real, intent(in)  :: eta_ohm,eta_hall,eta_ambi,dxr1,dyr1,dzr1
 real              :: dBohm(3),dBhall(3),dBambi(3)

 dBohm  = 0.0
 dBhall = 0.0
 dBambi = 0.0

 if (use_ohm ) call nimhd_get_DcrossR(dBohm ,jcurrent,dxr1,dyr1,dzr1,eta_ohm )
 if (use_hall) call nimhd_get_DcrossR(dBhall,jcb     ,dxr1,dyr1,dzr1,eta_hall)
 if (use_ambi) call nimhd_get_DcrossR(dBambi,jcbcb   ,dxr1,dyr1,dzr1,eta_ambi)
 dBnonideal = dBambi - dBhall - dBohm

end subroutine nimhd_get_dBdt
!-----------------------------------------------------------------------
!+
!  performs simple cross product and multiplies by eta
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_DcrossR(DcrossR,D_in,dx,dy,dz,eta)
 real, intent (in)  :: D_in(3)
 real, intent (out) :: DcrossR(3)
 real, intent (in)  :: dx, dy, dz, eta

 DcrossR(1) = (D_in(2)*dz - D_in(3)*dy)*eta
 DcrossR(2) = (D_in(3)*dx - D_in(1)*dz)*eta
 DcrossR(3) = (D_in(1)*dy - D_in(2)*dx)*eta

end subroutine nimhd_get_DcrossR
!----------------------------------------------------------------------!
end module nicil
