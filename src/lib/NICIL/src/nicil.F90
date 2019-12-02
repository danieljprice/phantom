!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!                                                                      !
! This is a stand-alone library that will calculate ionisation values  !
! and the coefficients for the non-ideal MHD terms: Ohmic resistivity, !
! Hall Effect and Ambipolar diffusion.                                 !
!                                                                      !
!                 Copyright (c) 2015-2019 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
!+
!  MODULE: nicil
!
!  DESCRIPTION:
!  Contains routines to calculation the ionisation rate, and related
!  useful quantities.
!  Copyright (c) 2015-2019 James Wurster
!  Reference: Wurster (2016) PASA, 33:e041.
!  See LICENCE file for usage and distribution conditions
!
!
!  REFERENCES: Asplund et al. (2009)
!              Cox (2000)
!              Draine & Lee (1984)
!              Fujii et al. (2011)
!              Keith & Wardle (2014)
!              Liu et al. (2003)
!              Nakano, Nishi & Umebayashi (2002)
!              Osterbrock (1961)
!              Pandey & Wardle (2008)
!              Pinto & Galli (2008)
!              Pollack et al. (1994)
!              Umebayashi & Nakano (2009)
!              Umebayashi & Nakano (1990)
!              Wardle (2007)
!              Wardle & Ng (1999)
!              Wurster, Price & Ayliffe (2014)
!
!  AUTHOR: James Wurster
!
!  PRIMARY RUNTIME PARAMETERS:
!    use_ohm        -- Calcualate and use the coefficient for Ohmic resistivity
!    use_hall       -- Calcualate and use the coefficient for the Hall effect
!    use_ambi       -- Calcualate and use the coefficient for ambipolar diffusion
!    g_cnst         -- Set constant grain size (true) or approximate an MRN grain distribution (false)
!    ion_rays       -- Include ionisation from cosmic (or x-) rays
!    ion_thermal    -- Include thermal ionisation
!    use_massfrac   -- Set mass fraction of H and He (true) or set abundances of 5 elements (false)
!    eta_constant   -- Use a constant resistivity
!    warn_verbose   -- To print warnings to file
!    use_fdg_in     -- Import dust mass density from the parent code for each point
!    rho_is_rhogas  -- set rho_gas == rho; else rho_gas == rho*(1-fdg)
!    fdg            -- Grain Parameter: gas to dust mass ratio
!    a0_grain       -- Grain Parameter if g_cnst=true : grain radius for constant grain size
!    an_grain       -- Grain Parameter if g_cnst=false: minimum grain radius for MRN distribution
!    ax_grain       -- Grain Parameter if g_cnst=false: maximum grain radius for MRN distribution
!    rho_bulk       -- Grain Parameter: bulk grain density
!    na_max         -- Number of grain sizes if using MRN distribution
!    zeta_of_rho    -- Use constant or variable cosmic ray ionisation rate
!    zeta_cgs       -- Ionisation rate (if zeta_of_rho=false)
!    zeta_CR_cgs    -- Unattenuated cosmic ray ionisation rate (if zeta_of_rho=true)
!    zeta_R_cgs     -- Ionisation rate of decaying radionuclides (if zeta_of_rho=true)
!    mass_MionR_mp  -- Mass of a metallic ion (for cosmic ray ionisation)
!    massfrac_X     -- If use_massfrac=true: Mass fraction of Hydrogen
!    massfrac_Y     -- If use_massfrac=true: Mass fraction of Helium
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
 !--Turn on/off individual non-ideal MHD coefficients
 logical, public            :: use_ohm           = .true.           ! Calculate the coefficient for Ohmic resistivity
 logical, public            :: use_hall          = .true.           ! Calculate the coefficient for the Hall effect
 logical, public            :: use_ambi          = .true.           ! Calculate the coefficient for ambipolar diffusion
 !--Set constant grain size (g_cnst=true) or approximate an MRN grain distribution (g_cnst=false)
 logical, public            :: g_cnst            = .true.
 !--Set ionisation type (can have one or both = true)
 logical, public            :: ion_rays          = .true.           ! Include ionisation from cosmic (or x-) rays (=true)
 logical, public            :: ion_thermal       = .true.           ! Include thermal ionisation (=true)
 !--Use mass fractions (true) or abundances (false)
 logical, public            :: use_massfrac      = .false.
 !--Use a constant (false) or variable (true) cosmic ray ionisation rate
 logical, public            :: zeta_of_rho       = .false.
 !--Import the dust mass density from the parent code for each point
 logical, public            :: use_fdg_in        = .false.
 !--Use rho in calculations (true); else use (1-fdg)*rho
 logical, public            :: rho_is_rhogas     = .false.
 !--Use constant resistivity coefficients for all three resistivity terms
 logical, public            :: eta_constant      = .false.
 !--Use the modified Hall parameters
 logical, public            :: mod_beta          = .false.
 !--To print warnings to file
 logical, public            :: warn_verbose      = .false.
 !--To solve a reordered Jacobian for stability
 logical, private           :: reorder_Jacobian  = .false.
 !
 !--Grain properties
 real,    public            :: fdg               = 0.01             ! gas to dust mass ratio
 real,    public            :: a0_grain          = 1.0d-5           ! grain radius for constant grain size [cm]
 real,    public            :: an_grain          = 5.0d-7           ! minimum grain radius for ~MRN distribution [cm]
 real,    public            :: ax_grain          = 2.5d-5           ! maximum grain radius for ~MRN distribution [cm]
 real,    public            :: rho_bulk          = 3.0              ! bulk grain density [g/cm^3]
 integer,         parameter :: na_max            = 5                ! number of bins of grain size for MRN distribution
 real,    private,parameter :: fdg_min           = 1.0d-14          ! the minimum dust to gas ratio permitted (> 0)
 real,    private,parameter :: fdg_max           = 1.0-fdg_min      ! the maximum dust to gas ratio permitted (< 1)
 !
 !--Cosmic ray ionisation
 integer, public, parameter :: nimass            =  2               ! Number of ion masses for cosmic ray ionisation
 real,    public            :: zeta_cgs          =  1.00d-17        ! ionisation rate [s^-1] (if zeta_of_rho=.false.)
 real,    public            :: zeta_CR_cgs       =  9.24d-18        ! unattenuated cosmic ray ionisation rate [s^-1] (if zeta_of_rho=.true.)
 real,    public            :: zeta_R_cgs        =  7.60d-19        ! ionisation rate of decaying radionuclides [s^-1] (if zeta_of_rho=.true.)
 real,    public            :: mass_MionR_mp     = 24.3             ! mass of ion (default is mass of magnesium) [m_proton]
 !
 real,    public            :: delta_gn          = 1.3              ! Ionisation parameter: multiplicative factor for sigmavgnbyT
 real,    public            :: pnH               = 0.667            ! Ionisation parameter: polarizability for H  [angstroms^3]
 real,    public            :: pnH2              = 0.804            ! Ionisation parameter: polarizability for H2 [angstroms^3]
 real,    public            :: pnHe              = 0.207            ! Ionisation parameter: polarizability for He [angstroms^3]
 real,    public            :: SigmaCR           = 96.0             ! Ionisation parameter: attenuation depth for cosmic rays [g cm^-2]
 real,    public            :: alpha_Mg          = 2.8d-12          ! Recombination coefficient for Mg
 real,    public            :: alpha_H           = 3.5d-12          ! Recombination coefficient for atomic Hydrogen
 real,    public            :: alpha_He          = 4.5d-12          ! Recombination coefficient for Helium
 real,    public            :: alpha_expT_Mg     = -0.86            ! Recombination exponent of temperature for Mg
 real,    public            :: alpha_expT_H      = -0.7             ! Recombination exponent of temperature for atomic Hydrogen
 real,    public            :: alpha_expT_He     = -0.67            ! Recombination exponent of temperature for Helium
 !
 !--Mass Fractions (used if use_massfrac=.true.)
 real,    public            :: massfrac_X        =     0.70         ! Mass fraction of Hydrogen
 real,    public            :: massfrac_Y        =     0.28         ! Mass fraction of Helium
 !--Thermal ionisation
 integer, public, parameter :: nelements_max     =     6            ! Maximum number of elements
 integer, public, parameter :: nlevels           =     2            ! Number of calculated ionisation levels
 integer,         parameter :: n_nuin            =  1000            ! Number of ni_inT_coef values to pre-caculate
 real,            parameter :: se                =    1.0           ! Electron ticking coefficient: se \in (1.0d-3, 1.0)
 real,            parameter :: m_min             =    1.0           ! minimum ion mass (units of m_proton) in the table
 real,            parameter :: m_max             =  101.0           ! maximum ion mass (units of m_proton) in the table
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
 !
 !--Threshholds
 real,    public            :: Texp_thresh0      = 0.005            ! Will set exp(-chi/kT) = 0 if T is too low
 !
 !--Additional parameters
 real,            parameter :: epsilon_coef      = 10.0             ! for subtractions, will assume 0 if abs(a-b)<epsilon_coef*epsilon(a)
 integer,         parameter :: NRctrmax          = 200              ! maximum number of Newton–Raphson before quiting
 real,            parameter :: NRtol             = 1.0e-8           ! tolerance on Newton–Raphson iterations
 real,    public            :: C_nimhd           = 0.1591549        ! Coefficient to control the timestep (==1/2pi)
 !
 !--END OF INPUT PARAMETERS
 !
 !--Misc. parameters not to be modified
 integer,         parameter :: zmax              =  1               ! Maximum grain charge; charge range is (-zmax,zmax)
 real,            parameter :: warn_ratio        =  0.1             ! fraction within with two values will be assumed equal for warnings
 real,            parameter :: vrms_min_kms      =  20              ! Minimum rms velocity below which rate coefficients become unreliable [km/s] (Pinto & Galli, 2008)
 real,            parameter :: vrms_max_kms      = 500              ! Maximum rms velocity above which rate coefficients become unreliable [km/s] (Pinto & Galli, 2008)
 integer, public, parameter :: n_data_out        = 17+nelements_max*nlevels-3+1 ! number of array element the optional array data_out

 !
 !--Physical Constants (CGS)
 real,            parameter :: pi                =  3.1415926536d0  !  pi
 real,            parameter :: twopi             =  6.2831853072d0  ! 2pi
 real,            parameter :: fourpi            = 12.5663706144d0  ! 4pi
 real,            parameter :: c                 =  2.997924d10     ! Speed of light [cm/s]
 real,            parameter :: qe                =  4.8032068d-10   ! charge on electron [esu == statC == (cm^3 g)^0.5 s^-1]
 real,            parameter :: mass_electron_cgs =  9.10938291d-28  ! Electron mass [g]
 real,            parameter :: eV                =  1.60217657d-12  ! Electron volts [ergs]
 real,            parameter :: mass_proton_cgs   =  1.67262158d-24  ! Proton mass [g]
 real,            parameter :: kboltz            =  1.38066d-16     ! Boltzmann constant  [erg/K]
 real,            parameter :: Ggrav             =  6.67408d-8      ! Gravitational constant [cm^3 g^-1 s^-2]
 real,            parameter :: planckh           =  6.6260755d-27   ! Planck's Constant [erg s]
 real,            parameter :: Amrn              =  1.5d-25         ! Constant for MRN grain distribution [cm^2.5]
 !
 !--Indicies for the various species; will allow for cleaner bookkeeping arrays
 integer,         parameter :: nspecies_max      =  5+2*na_max      ! The maximum number of charged species (do not modify)
 integer,         parameter :: ine               =  1               ! index for electron number density
 integer,         parameter :: iniHR             =  2               ! index for ion density of light elements (for CR ionisation)
 integer,         parameter :: iniMR             =  3               ! index for ion density of metallic elements (for CR ionisation)
 integer,         parameter :: inisT             =  4               ! index for ion density of singly ionised atoms (for thermal ionisation)
 integer,         parameter :: inidT             =  5               ! index for ion density of doubly ionised atoms (for thermal ionisation)
 integer,         parameter :: ing               =  6               ! index for the number density of the first grain
 !--Indicies for H,H2 and He
 integer,         parameter :: iH2               =  1               ! index for molecular hydrogen
 integer,         parameter :: iH                =  2               ! index for atomic hydrogen
 integer,         parameter :: iHe               =  3               ! index for Helium
 !--Indicies for warnings
 integer,         parameter :: ierr_nRconv       =       1          ! fatal error code if n_R(ion,e .or. grain) did not converge
 integer,         parameter :: ierr_neTconv      =       2          ! fatal error code if n_electronT did not converge
 integer,         parameter :: ierr_neTle0       =       5          ! fatal error code if n_electronT < 0
 integer,         parameter :: ierr_nRle0        =      10          ! fatal error code if n_R < 0
 integer,         parameter :: ierr_Zge0         =      20          ! fatal error code if Z_g > 0 for average Z_grain calculation
 integer,         parameter :: ierr_fdg_in       =      50          ! fatal error code if use_fdg_in=.true. and fdg_in not passed in
 integer,         parameter :: ierr_T            =     100          ! warning code if T==0 (from input)
 integer,         parameter :: ierr_B            =     200          ! warning code if B==0 (from input)
 integer,         parameter :: ierr_Zave         =     500          ! warning code if using average Z_grain method to calculate n_R
 integer,         parameter :: ierr_LagH2        =    1000          ! warning code if Langevin rate used for molecular Hydrogen
 integer,         parameter :: ierr_LagH         =    2000          ! warning code if Langevin rate used for atomic Hydrogen
 integer,         parameter :: ierr_LagHe        =    5000          ! warning code if Langevin rate used for Helium
 integer,         parameter :: ierr_sigH         =   10000          ! warning code if sigmaH == 0 to prevent errors from roundoff
 integer,         parameter :: ierr_scN          =   20000          ! warning code if strong coupling approximation broke (rho_n ~ rho is false)
 integer,         parameter :: ierr_scI          =   50000          ! warning code if strong coupling approximation broke (rho_i <<rho is false)
 integer,         parameter :: ierr_drift        =  100000          ! warning code if drift velocity is large relative to absolute velocity
 integer,         parameter :: ierr_vrms         =  200000          ! warning code if velocity is outside of the range valid for the rate coefficients
 integer,         parameter :: ierr_nedec        =  500000          ! warning code if n_electronR ~ n_electronT
 integer,         parameter :: ierr_alion        = 1000000          ! warning code if fully ionised
 !
 !--Local Variables to the NICIL module that may be required elsewhere in the user's code
 integer, public    :: nelements
 real,    public    :: meanmolmass,unit_eta
 !--Local Variables to the NICIL module
 integer, private   :: iprint,iprintw,nspecies,na,neqn
 real,    private   :: csqbyfourpi,small,epsilonR,coef_epsilonR,coef2_epsilonR,coef4_epsilonR
 real,    private   :: warn_ratio2,warn_ratio_m1,warn_ratio_p1,warn_ratio_m12,warn_ratio_p12
 real,    private   :: mass_proton,mass_proton1,mass_neutral_cold1,mfrac_by_m_notH_mp,one_minus_fdg
 real,    private   :: nu_ei_coef
 real,    private   :: sigma_coef,sigmavenXH_LR,sigmavenXH2_LR,sigmavenY_LR,sigmavenX_coef,sigmavenY_coef
 real,    private   :: eta_ohm_cnst,eta_hall_cnst,eta_ambi_cnst,alpha_AD_p1
 real,    private   :: zeta0,zetaCR,zetaR,zetaCRcoef,umass0,unit_density0
 real,    private   :: Saha_coef_HH2,chij_HH2,chiHH2
 real,    private   :: vrms2_min,vrms2_max
 real,    private   :: a_grain(na_max),a_grain_cgs(na_max),n_grain_coef(na_max)
 real,    private   :: sigmaviRn(3,nimass)
 real,    private   :: k_fac_coef(-zmax:zmax,na_max)
 real,    private   :: k_ig_coef(nimass,na_max),k_eg_coef(na_max),k_ei_coef(nimass+1)
 real,    private   :: n_grain_dist(na_max),n_R0(nspecies_max)
 real,    private   :: log_abunj(nelements_max),mj_mp(nelements_max),abundancej(nelements_max)
 real,    private   :: mass_frac(nelements_max),massj_mp(nspecies_max),mj_mp1(nelements_max)
 real,    private   :: Saha_coef(nlevels,nelements_max),chij(nlevels,nelements_max)
 real,    private   :: gjp1_gj(nlevels,nelements_max),Texp_thresh(nlevels,nelements_max)
 real,    private   :: Zj(nspecies_max),aZj(nspecies_max),sZj(nspecies_max)
 real,    private   :: massj(nspecies_max),beta_coef(nspecies_max),nu_jn_coef(nspecies_max)
 character(len=2),private   :: symj(nelements_max)
 !
 !--Subroutines
 public  :: nicil_initialise,nicil_get_ion_n,nicil_get_eta,nicil_get_vion,nicil_get_vdrift,nicil_get_halldrift
 public  :: nimhd_get_jcbcb,nimhd_get_dBdt,nimhd_get_dudt,nimhd_get_dt
 public  :: nicil_translate_error
 private :: nicil_initialise_species,nicil_print_summary,nicil_ic_error,nicil_version
 private :: nicil_ion_get_sigma,nicil_ionR_get_n,nicil_ionR_get_n_via_Zave,nicil_get_HH2_ratio
 private :: nicil_ionT_calc_k,nicil_ionT_calc_Jf,get_dn_viaLU
 private :: nicil_ionR_predict_ng,nicil_fill_nRcomplete,nicil_unfill_nRcomplete,nicil_ionR_set_guess
 private :: nicil_ionT_get_ne,nicil_ionT_get_n,nicil_ionT_get_nj_Kjk,nicil_ionT_get_nion
 private :: nicil_nimhd_get_eta,nicil_nimhd_get_eta_cnst,nimhd_get_DcrossR
 !
 private
!
contains
!+
!----------------------------------------------------------------------!
!+
! Internal Version Control
! important modifications are listed before version number
!+
!----------------------------------------------------------------------!
pure subroutine nicil_version(version)
 character(len=200), intent(out) :: version
           !  8 Dec  2015: Initial Version
 version = "Version 1.0: 8 Dec 2015: Initial Version"
           ! 27 Jan  2016: If iterations fail to converge, will reset initial guess and try again
           !  3 Mar  2016: Thermal ionisation can doubly ionise atoms
           !  8 Mar  2016: Added light ion for Cosmic ray ionisation, thus there is now heavy and light ions
           ! 11 Mar  2016: Rewrote MRN grain calculations to be sums with characteristic sizes, and not a single average
 version = "Version 1.1: 26 April 2016"
           !  7 June 2016: bug fix: collision rates account for arbitrary charge
           !  7 June 2016: bug fix: added statistical weights to Saha equation
           !  8 June 2016: changed default: mod_beta = .false.
           !  8 June 2016: modified rate equations to be from Table 1 of Pinto & Galli (2008)
           !  9 June 2016: added warnings that can be optionally printed to file
           !  9 June 2016: added optional subroutine to calculate ion velocity
           !  9 June 2016: added dissociation of molecular hydrogen for self-consistent fractions of H and H2
           ! 27 June 2016: Rewrote cosmic ray ionisation algorithm.  Solves for n_g(Z=-1),n_g(Z=0),n_g(Z=+1) rather than \bar{Z}
           !               Solving for \bar{Z} remains as a contingency if convergence is not obtained with the above method
 version = "Version 1.2: 27 June 2016"
           !  6 July 2016: bug fix in nimhd_get_dudt
           !  1 Aug  2016: will exit cleanly if B=0 is passed in
           !  1 Aug  2016: bug fix when calculating k_ei(1); temperature dependence now included
 version = "Version 1.2.1: 2 Aug 2016"
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
 version = "Version 1.2.3: 5 July 2017"
           ! 10 Apr  2018: Ambipolar diffusion now uses a subtraction rather than a double loop
 version = "Version 1.2.4: 11 April 2018"
           ! 19 Oct  2018: if rho_n < 0, will scale the neutral components from thermal and rays by their respective n_electron
           !               rho_is_rhogas = .false. by default
           !               n_gas rather than n_total is used for calculation thermal ionisation (if rho_is_rhogas = .false.)
 version = "Version 1.2.5: 19 October 2018"
           !  4 Feb  2019: bug fix when using MRN grain distribution and size(n_R) .ne. 2+2*na
 version = "Version 1.2.6: 28 August 2019"

end subroutine nicil_version
!----------------------------------------------------------------------!
!+
! Define the properties of the species that will be used for thermal
! ionisation
!+
!----------------------------------------------------------------------!
subroutine nicil_initialise_species
 integer :: i
 !
 !--Initialise/Zero the values
 log_abunj     = 0.0
 chij          = 0.0                ! Default the ionisation potential to huge
 mj_mp         = 0.0
 gjp1_gj       = 0.0                ! ratio of statistical weights
 !
 !--Hydrogen (molecular)
 i             = iH2
 symj(i)       = "H2"               ! Chemical Symbol
 log_abunj(i)  =  0.0               ! Log abundance
 mj_mp(i)      =  2.02              ! Mass [m_proton]
 chij(1,i)     = 15.60              ! First  ionisation potential [eV]
 chiHH2        =  4.476             ! Dissociation potential [eV]
 gjp1_gj(1,i)  = 1.0/2.0            ! ratio of statistical weights (first/ground)
 !
 if (use_massfrac) then
    nelements     = 3               ! Total number of elements
    !--Hydrogen (atomic)
    i             = iH
    symj(i)       = "H"             ! Chemical Symbol
    mj_mp(i)      =  1.00           ! Mass [m_proton]
    chij(1,i)     = 13.60           ! First  ionisation potential [eV]
    gjp1_gj(1,i)  = 1.0/2.0         ! ratio of statistical weights (first/ground)
    mass_frac(i)  = massfrac_X
    !--Helium
    i             = iHe
    symj(i)       = "He"            ! Chemical Symbol
    mj_mp(i)      =  4.00           ! Mass [m_proton]
    chij(1,i)     = 24.59           ! First  ionisation potential [eV]
    chij(2,i)     = 54.42           ! Second ionisation potential [eV]
    gjp1_gj(1,i)  = 2.0/1.0         ! ratio of statistical weights (first/ground)
    gjp1_gj(2,i)  = 1.0/2.0         ! ratio of statistical weights (second/first)
    mass_frac(i)  = massfrac_Y
 else
    nelements     = 6               ! Total number of elements
    !--Hydrogen (atomic)
    i             = iH
    symj(i)       = "H"             ! Chemical Symbol
    log_abunj(i)  = 12.0            ! Log abundance
    mj_mp(i)      =  1.01           ! Mass [m_proton]
    chij(1,i)     = 13.60           ! First  ionisation potential [eV]
    gjp1_gj(1,i)  = 1.0/2.0         ! ratio of statistical weights (first/ground)
    !--Helium
    i             = iHe
    symj(i)       = "He"            ! Chemical Symbol
    log_abunj(i)  = 10.93           ! Log abundance
    mj_mp(i)      =  4.00           ! Mass [m_proton]
    chij(1,i)     = 24.59           ! First  ionisation potential [eV]
    chij(2,i)     = 54.42           ! Second ionisation potential [eV]
    gjp1_gj(1,i)  = 2.0/1.0         ! ratio of statistical weights (first/ground)
    gjp1_gj(2,i)  = 1.0/2.0         ! ratio of statistical weights (second/first)
    !--Sodium
    i             = i + 1
    symj(i)       = "Na"            ! Chemical Symbol
    log_abunj(i)  =  6.24           ! Log abundance
    mj_mp(i)      = 22.98           ! Mass [m_proton]
    chij(1,i)     =  5.14           ! First  ionisation potential [eV]
    chij(2,i)     = 47.29           ! Second ionisation potential [eV]
    gjp1_gj(1,i)  = 1.0/2.0         ! ratio of statistical weights (first/ground)
    gjp1_gj(2,i)  = 6.0/1.0         ! ratio of statistical weights (second/first)
    !--Magnesium
    i             = i + 1
    symj(i)       = "Mg"            ! Chemical Symbol
    log_abunj(i)  =  7.60           ! Log abundance
    mj_mp(i)      = 24.31           ! Mass [m_proton]
    chij(1,i)     =  7.65           ! First  ionisation potential [eV]
    chij(2,i)     = 15.03           ! Second ionisation potential [eV]
    gjp1_gj(1,i)  = 2.0/1.0         ! ratio of statistical weights (first/ground)
    gjp1_gj(2,i)  = 1.0/2.0         ! ratio of statistical weights (second/first)
    !--Potassium
    i             = i + 1
    symj(i)       = "K"             ! Chemical Symbol
    log_abunj(i)  =  5.03           ! Log abundance
    mj_mp(i)      = 39.10           ! Mass [m_proton]
    chij(1,i)     =  4.34           ! First  ionisation potential [eV]
    chij(2,i)     = 31.62           ! Second ionisation potential [eV]
    gjp1_gj(1,i)  = 1.0/2.0         ! ratio of statistical weights (first/ground)
    gjp1_gj(2,i)  = 6.0/1.0         ! ratio of statistical weights (second/first)
 endif
 !
end subroutine nicil_initialise_species
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
subroutine nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint_in,iprintw_in)
 real,              intent(in)  :: utime,umass,udist,unit_Bfield
 integer,           intent(out) :: ierr
 integer, optional, intent(in)  :: iprint_in,iprintw_in
 integer                        :: j,k,p,z,zn,zp
 real                           :: unit_velocity,unit_density,unit_erg
 real                           :: sigmavgnbyT_coef
 real                           :: mass_proton_mp,mass_neutral_cold_cgs,mass_grain_mp
 real                           :: massj_cgs(nspecies_max),n_grain_coef_cgs(na_max)
 real                           :: a01_grain,dloga,n_grain_coef_tot
 real                           :: mu_iXH,mu_iXH2,mu_iY,mu_eXH,mu_eXH2,mu_eY
 real                           :: mass_total,n_rel,nj_rel(nelements_max)
 !
 !--Initialise species properties for thermal ionisation & Calculate abundances & mass fractions
 !  By construction, Sum (abundancej), Sum(mass_frac) == 1
 call nicil_initialise_species
 !
 !--Verify input parameters are realistic; print error messages for each invalid error
 ierr = 0
 if (fdg  >= 1.0 .or. fdg < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid dust-to-gas fraction','fdg',fdg)
 if ( g_cnst ) then
    if (a0_grain  < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid fixed grain radius:','a_grain',a0_grain)
 else
    if (an_grain < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid minimum MRN grain radius:','an_grain',an_grain)
    if (ax_grain < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid maximum MRN grain radius:','ax_grain',ax_grain)
    if (ax_grain < an_grain) call nicil_ic_error(ierr,'nicil_initialise','Invalid radii ordering of MRN grain radii')
 endif
 if (rho_bulk < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid bulk grain density:','rho_bulk',rho_bulk)
 if (zeta_of_rho) then
    if (zeta_CR_cgs < 0.0)call nicil_ic_error(ierr,'nicil_initialise','Invalid unattenuated ionisation rate:', 'zetaCR',zeta_CR_cgs)
    if (zeta_R_cgs  < 0.0)call nicil_ic_error(ierr,'nicil_initialise','Invalid radionuclides ionisation rate:','zetaR', zeta_R_cgs)
 else
    if (zeta_cgs    < 0.0)call nicil_ic_error(ierr,'nicil_initialise','Invalid ionisation rate:',               'zeta', zeta_cgs)
 endif
 if (zmax /= 1)           call nicil_ic_error(ierr,'nicil_initialise','Current Jacobian decomposition requires z_grain = -1,0,+1')
 if (use_fdg_in .and. .not.g_cnst) &
                          call nicil_ic_error(ierr,'nicil_initialise','Cannot use multiple grain sizes and pass in dust density')
 if (mass_MionR_mp < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid ion mass:'     ,'m_Mion'  ,mass_MionR_mp)
 if (use_massfrac) then
    if (massfrac_X < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid Hydrogen mass fraction:','X',massfrac_X)
    if (massfrac_X < 0.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid helium mass fraction:'  ,'Y',massfrac_Y)
    if (massfrac_X+massfrac_Y>1.0) &
      call nicil_ic_error(ierr,'nicil_initialise','Invalid metalicities:'  ,'massfrac_X + massfrac_Y',massfrac_X+massfrac_Y)
 endif
 if (epsilon_coef <= 1.0) call nicil_ic_error(ierr,'nicil_initialise','Invalid subtraction constraint','epsilon_coef',epsilon_coef)
 if (NRtol < 1.0e-15 .or. NRtol > 0.1) &
                     call nicil_ic_error(ierr,'nicil_initialise','Poor constraint on Newton–Raphson tolerance','NRtol',NRtol)
 if (NRctrmax <= 10) call nicil_ic_error(ierr,'nicil_initialise','Too few maximum permitted Newton–Raphson iterations' &
 ,'NRctrmax',real(NRctrmax))
 if (nelements > nelements_max)  call nicil_ic_error(ierr,'nicil_initialise','nelements > nelements_max')
 if (iniHR > iniMR)  call nicil_ic_error(ierr,'nicil_initialise','indicies for light & heavy ion densities out of order')
 if (inisT > inidT)  call nicil_ic_error(ierr,'nicil_initialise','indicies for singly & doubly ionised densities out of order')
 if (trim(symj(iH))/="H" )  call nicil_ic_error(ierr,'nicil_initialise_species','Hygrogen does not have the correct array number')
 if (trim(symj(iHe))/="He") call nicil_ic_error(ierr,'nicil_initialise_species','Helium does not have the correct array number')
 if (eta_constant .and. eta_const_type/=icnstphys .and. eta_const_type/=icnstsemi .and. eta_const_type/=icnst) then
    call nicil_ic_error(ierr,'nicil_initialise','Invalid choice of (semi-) constant eta calculations; correct eta_const_type.')
 endif
 !--Abort setup if errors fatal exist
 if (ierr/=0) return
 !--Reset logicals depending on user's current selection
 !  Set ionisations false if using constant coefficients or not calculating any non-ideal MHD terms
 if (eta_constant .or. (.not.use_ohm .and. .not.use_hall .and. .not.use_ambi)) then
    ion_rays    = .false.
    ion_thermal = .false.
 endif
 !  Set calculate coefficients false if no ionisation and eta_constant = false
 !  (this is not a fatal error, since the user may require ideal mhd while not removing the nicil algorithm)
 if (.not.eta_constant .and. .not.ion_rays .and. .not.ion_thermal) then
    use_ohm     = .false.
    use_hall    = .false.
    use_ambi    = .false.
 endif
 !
 !--Initialise general parameters
 small          = epsilon(small)*epsilon(small)
 epsilonR       = epsilon(epsilonR)
 coef_epsilonR  = epsilon_coef*epsilonR
 coef2_epsilonR = epsilon_coef**2*epsilonR
 coef4_epsilonR = epsilon_coef**4*epsilonR
 warn_ratio2    = warn_ratio*warn_ratio
 warn_ratio_m1  = 1.0-warn_ratio
 warn_ratio_p1  = 1.0+warn_ratio
 warn_ratio_m12 = warn_ratio_m1*warn_ratio_m1
 warn_ratio_p12 = warn_ratio_p1*warn_ratio_p1
 vrms2_min      = (vrms_min_kms*1.0d5*utime/udist)**2
 vrms2_max      = (vrms_max_kms*1.0d5*utime/udist)**2
 mj_mp1         = 1.0/mj_mp
 !
 !--Calculate mass fractions or abundances, as required
 if (use_massfrac) then
    ! Note: Abundance is approximate since we do not explicitly account for Z
    !       We explicitly assume that the hydrogen is molecular; given how the
    !       abundances are used in NICIL, this abundance must be associated
    !       with H
    mass_total       = 1.0/(massfrac_X*mj_mp1(iH2) + massfrac_Y*mj_mp1(iHe))
    abundancej(2)    = massfrac_X*mass_total*mj_mp1(iH2)
    abundancej(3)    = massfrac_Y*mass_total*mj_mp1(iHe)
    massj_mp(iniHR)  = mass_total                                   ! Light element ion mass
 else
    n_rel       = 0.0
    nj_rel      = 0.0
    mass_total  = 0.0
    do j = 1,nelements
       if (log_abunj(j)>0.0) nj_rel(j) = 10**(log_abunj(j) - 12.0)  ! Number density relative to H
       n_rel     = n_rel + nj_rel(j)
    enddo
    do j = 1,nelements
       abundancej(j) = nj_rel(j)/n_rel                              ! Abundance calculated using H input
       mass_total    = mass_total + mj_mp(j)*abundancej(j)
    enddo
    do j = 1,nelements
       mass_frac(j)  = mj_mp(j)*abundancej(j)/mass_total
    enddo
 endif
 !--Calculate the mean molar mass of cold gas; assumes all hydrogen is molecular
 meanmolmass        = 0.0
 mfrac_by_m_notH_mp = 0.0
 do j = 1,nelements
    if (trim(symj(j))=="H") then
       meanmolmass        = meanmolmass        + mass_frac(j)/(2.0*mj_mp(j))
    else
       meanmolmass        = meanmolmass        + mass_frac(j)/mj_mp(j)
       mfrac_by_m_notH_mp = mfrac_by_m_notH_mp + mass_frac(j)/mj_mp(j)
    endif
    if (j==iHe) massj_mp(iniHR) = 1.0/meanmolmass                                 ! Light element ion mass
 enddo
 meanmolmass = 1.0/meanmolmass
 !
 !--Unit conversions from cgs-> code
 unit_velocity     = udist/utime
 unit_density      = umass/udist**3
 unit_erg          = umass*unit_velocity**2
 unit_eta          = udist**2/utime
 umass0            = umass
 unit_density0     = unit_density
 !
 !--Initialise masses (m_proton units)
 mass_proton_mp    = 1.0
 massj_mp(ine)     = mass_electron_cgs/mass_proton_cgs
 massj_mp(iniMR)   = mass_MionR_mp
 !
 !--Initialise variables (CGS variables)
 massj_cgs              = 0.0                                                ! Set all masses to zero (required for thermally ionised masses)
 massj_cgs(ine)         = mass_electron_cgs                                  ! Electron mass
 massj_cgs(iniHR:iniMR) = massj_mp(iniHR:iniMR)*mass_proton_cgs              ! Ion mass for cosmic ray ionisation
 mass_neutral_cold_cgs  = meanmolmass          *mass_proton_cgs              ! Cold neutral mass
 if (.not. use_massfrac) then
    massfrac_X           = mass_frac(iH2)                                    ! mass fraction of Hydrogen
    massfrac_Y           = mass_frac(iHe)                                    ! mass fraction of Helium
 endif
 !
 !--Determine grain distribution and initialise
 a_grain_cgs      = 0.0
 n_grain_coef_cgs = 0.0
 nspecies         = 5
 n_grain_dist     = 1.0
 mass_grain_mp    = 0.0
 a01_grain        = 0.0
 if (ion_rays) then
    if (g_cnst) then
       na = 1
       a_grain_cgs(1)   = a0_grain
       massj_cgs(ing)   = fourpi/3.0*a_grain_cgs(1)**3*rho_bulk
       massj_cgs(ing+1) = massj_cgs(ing)
       mass_grain_mp    = massj_cgs(ing)
       if (use_fdg_in) then
          n_grain_coef_cgs(1) = mass_neutral_cold_cgs/massj_cgs(ing)
       else
          n_grain_coef_cgs(1) = mass_neutral_cold_cgs*fdg/massj_cgs(ing)
       endif
    else
       na          = na_max
       dloga       = (log10(ax_grain) - log10(an_grain))/float(na)
       n_grain_coef_tot = 0.0
       do j = 1,na
          a_grain_cgs(j)         = 10**(log10(an_grain) +(j-0.5)*dloga )
          massj_cgs(ing+2*(j-1)) = fourpi/3.0*a_grain_cgs(j)**3*rho_bulk
          massj_cgs(ing+2*j-1)   = massj_cgs(ing+2*(j-1))
          n_grain_coef_cgs(j)    = 0.4*Amrn*(( 10**(log10(an_grain) +(j-1)*dloga ) )**(-2.5) &
                                 - (10**(log10(an_grain) +(j)*dloga ))**(-2.5))
          n_grain_coef_tot       = n_grain_coef_tot + n_grain_coef_cgs(j)
          mass_grain_mp          = mass_grain_mp + massj_cgs(ing+j-1)
          a01_grain              = a01_grain + a_grain_cgs(j)
       enddo
       n_grain_dist = n_grain_coef_cgs/n_grain_coef_tot
    endif
    massj_mp(ing:ing+2*na-1) = massj_cgs(ing:ing+2*na-1)/mass_proton_cgs
    nspecies      = nspecies + 2*na
    mass_grain_mp = mass_grain_mp/(mass_proton_cgs*na)
    a01_grain     = a01_grain/na
 endif
 one_minus_fdg = 1.0 - fdg
 !
 !--Initialise variables (Code unit variables)
 csqbyfourpi        = c**2/fourpi               / unit_eta
 zeta0              = zeta_cgs                  * utime
 zetaCR             = zeta_CR_cgs               * utime
 zetaR              = zeta_R_cgs                * utime
 mass_proton        = mass_proton_cgs           / umass
 mass_proton1       = 1.0/mass_proton
 mass_neutral_cold1 = 1.0/mass_neutral_cold_cgs * umass
 massj              = massj_cgs                 / umass
 a_grain            = a_grain_cgs               / udist
 n_grain_coef       = n_grain_coef_cgs
 !
 !--Set default values as a fraction of total number density
 n_R0 = epsilon(small)
 n_R0(1:1+nimass) = 1.0/float(nimass)
 !
 !--Coefficient for the exponential term if density-dependent cosmic ionisation rate
 !  to be multiplied by sqrt (T rho/ m_n); udist converts rho/m_n to cgs, making this term dimensionless
 zetaCRcoef = sqrt(kboltz/(pi*Ggrav*mass_proton*udist**3))/SigmaCR
 !
 !--Charge capture Rates coefficients (Code units and scaled)
 do j = 1,na
    do z = -zmax,zmax
       k_fac_coef(z,j) = float(z)*qe**2/(a_grain_cgs(j)*kboltz)
    enddo
    k_eg_coef(j)       = a_grain_cgs(j)**2*sqrt(8.0*pi*kboltz/mass_electron_cgs)     * utime/udist**3
    do k = 1,nimass
       k_ig_coef(k,j)  = a_grain_cgs(j)**2*sqrt(8.0*pi*kboltz/massj_cgs(k-1+iniHR) ) * utime/udist**3
    enddo
 enddo
 k_ei_coef(1) = massfrac_X*alpha_H /300.0**alpha_expT_H                              * utime/udist**3   ! iniHR (Hydrogen component)
 k_ei_coef(2) = massfrac_Y*alpha_He/300.0**alpha_expT_he                             * utime/udist**3   ! iniHR (Helium component)
 k_ei_coef(3) = alpha_Mg/300.0**alpha_expT_Mg                                        * utime/udist**3   ! iniMR
 neqn         = nimass + 1 + 3*na  ! ions, electrons, (Z=-1,0,+1)
 !
 !--Rate coefficientes (CGS; sigma units are cm^3/s)
 !  *_LR are Langevin rates, and the maximum sigmaenX is used at any given time
 !  Note that this sigmaviRn is for cosmic ray ionisation, thus Zion == 1
 mu_eXH2           = (mj_mp(iH2) * massj_mp(ine))       /(mj_mp(iH2) + massj_mp(ine)      )
 mu_eXH            = (mj_mp(iH ) * massj_mp(ine))       /(mj_mp(iH ) + massj_mp(ine)      )
 mu_eY             = (mj_mp(iHe) * massj_mp(ine))       /(mj_mp(iHe) + massj_mp(ine)      )
 do j = 1,nimass
    mu_iXH2        = (mj_mp(iH2) * massj_mp(j-1+iniHR)) /(mj_mp(iH2) + massj_mp(j-1+iniHR))
    mu_iXH         = (mj_mp(iH ) * massj_mp(j-1+iniHR)) /(mj_mp(iH ) + massj_mp(j-1+iniHR))
    mu_iY          = (mj_mp(iHe) * massj_mp(j-1+iniHR)) /(mj_mp(iHe) + massj_mp(j-1+iniHR))
    sigmaviRn(1,j) = 2.81d-9              * sqrt(pnH2 / mu_iXH2)
    sigmaviRn(2,j) = 2.81d-9              * sqrt(pnH  / mu_iXH )
    sigmaviRn(3,j) = 2.81d-9 * massfrac_Y * sqrt(pnHe / mu_iY  )
 enddo
 sigmavenX_coef    = 1.0d-9
 sigmavenY_coef    = massfrac_Y * 0.428 * 1.0d-9
 sigmavenXH2_LR    =            2.81d-9*sqrt(pnH2 / mu_eXH2)
 sigmavenXH_LR     =            2.81d-9*sqrt(pnH  / mu_eXH )
 sigmavenY_LR      = massfrac_Y*2.81d-9*sqrt(pnHe / mu_eY  )
 sigmavgnbyT_coef  = delta_gn * sqrt(128.0*pi*kboltz/(9.0*mass_proton_cgs))
 !
 !--Collisional Frequencies (CGS = s^-1)
 nu_ei_coef              = 51.0                      / udist**3                          ! nu_ei
 nu_jn_coef(ine)         = 1.0      /mass_proton_cgs * unit_density                      ! nu_en
 nu_jn_coef(iniHR:iniMR) = 1.0      /mass_proton_cgs * unit_density                      ! nu_inR
 nu_jn_coef(inisT:inidT) = 2.81d-9  /mass_proton_cgs * unit_density                      ! nu_inT
 p = ing
 do j = 1,na
    do z = 1,2  ! values are required for Z = -1 and +1
       nu_jn_coef(p) = sigmavgnbyT_coef*a_grain_cgs(j)**2/mass_proton_cgs * unit_density ! nu_gn
       p = p + 1
    enddo
 enddo
 !
 !--Fill charge arrays (electric charge,Z; absolute value of Z; sign of Z)
 Zj(ine)   = -1.0; aZj(ine)   = 1.0; sZj(ine)   = -1.0        ! electrons
 Zj(iniHR) =  1.0; aZj(iniHR) = 1.0; sZj(iniHR) =  1.0        ! ions of light elements (for CR ionisation)
 Zj(iniMR) =  1.0; aZj(iniMR) = 1.0; sZj(iniMR) =  1.0        ! ions of metallic elements (for CR ionisation)
 Zj(inisT) =  1.0; aZj(inisT) = 1.0; sZj(inisT) =  1.0        ! ions of singly ionised atoms (for thermal ionisation)
 Zj(inidT) =  2.0; aZj(inidT) = 2.0; sZj(inidT) =  1.0        ! ions of doubly ionised atoms (for thermal ionisation)
 do j = 1,na
    zn = 1+inidT+2*(j-1)
    zp = 2+inidT+2*(j-1)
    Zj(zn)  = -1.0; aZj(zn)    = 1.0; sZj(zn)    = -1.0        ! negatively charged grains
    Zj(zp)  =  1.0; aZj(zp)    = 1.0; sZj(zp)    =  1.0        ! positively charged grains
 enddo
 !
 !--Hall parameters (dimensionless after multipiled by B_code (and mass_ionT_code))
 beta_coef(ine)          = aZj(ine)        *qe/(massj_cgs(ine)        *c) * unit_Bfield
 beta_coef(iniHR:iniMR)  = aZj(iniHR:iniMR)*qe/(massj_cgs(iniHR:iniMR)*c) * unit_Bfield
 beta_coef(inisT:inidT)  = aZj(inisT:inidT)*qe/(                       c) * unit_Bfield / umass
 p = ing
 do j = 1,na
    do z = 1,2  ! values are required for Z = -1 and _+1
       beta_coef(p)      = aZj(p)          *qe/(massj_cgs(p)          *c) * unit_Bfield
       p = p + 1
    enddo
 enddo
 !
 !--Conductivity coefficient (code units)
 sigma_coef   = qe*c                                                      / (udist**3 * unit_Bfield)
 !
 !--Coefficient for the Saha equation (Code units)
 Saha_coef_HH2 =             (    pi*mass_proton*kboltz/planckh**2 * utime**2*unit_erg )**1.5
 Saha_coef     = 2.0*gjp1_gj*(2.0*pi*massj(ine) *kboltz/planckh**2 * utime**2*unit_erg )**1.5
 chij          = chij*eV/kboltz
 chij_HH2      = chiHH2*eV/kboltz
 Texp_thresh   = chij*Texp_thresh0
 !
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
 !
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
 !
 !--Print Statements to summarise conditions used
 call nicil_print_summary(a01_grain,massj_mp(iniHR),mass_grain_mp,meanmolmass)
 !
end subroutine nicil_initialise
!----------------------------------------------------------------------!
!+
! Print Statements to summarise what is being used here
!+
!----------------------------------------------------------------------!
subroutine nicil_print_summary(a01_grain,mass_HionR_mp,mass_grain_mp,meanmolmass)
 real, intent(in)   :: a01_grain,mass_HionR_mp,mass_grain_mp,meanmolmass
 integer            :: j
 character(len=  2) :: comma
 character(len=200) :: ni_terms,version
 !
 call nicil_version(version)
 write(iprint,'(2a)') "NICIL: ",trim(version)
 write(iprint,'(a)' ) "NICIL: Copyright (c) 2015-2019 James Wurster"
 write(iprint,'(a)' ) "NICIL: See LICENCE file for usage and distribution conditions"
 write(iprint,'(a)' ) "NICIL: Reference: Wurster (2016) PASA, 33:e041."
 !
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
 !
 if (eta_constant) then
    write(iprint,'(2a)')                       "NICIL: Non-ideal terms used: ",trim(ni_terms)
    write(iprint,'(a)' )                       "NICIL: All resistivity coefficients are constant."
    if (eta_const_type==icnst) then
       if (use_ohm ) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ohm  = ", eta_ohm_cnst*unit_eta,  " cm^2 s^{-1}"
       if (use_hall) write(iprint,'(a,Es10.3,a)') "NICIL: eta_Hall = ", eta_hall_cnst*unit_eta, " cm^2 s^{-1}"
       if (use_ambi) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ambi = ", eta_ambi_cnst*unit_eta, " cm^2 s^{-1}"
    else
       if (use_ohm ) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ohm  = ", eta_ohm_cnst*unit_eta, "           cm^2 s^{-1}"
       if (use_hall) write(iprint,'(a,Es10.3,a)') "NICIL: eta_Hall = ", eta_hall_cnst*unit_eta,"*|B|       cm^2 s^{-1}"
       if (use_ambi) write(iprint,'(a,Es10.3,a)') "NICIL: eta_ambi = ", eta_ambi_cnst*unit_eta,"*|B|^2/rho cm^2 s^{-1}"
    endif
    write(iprint,'(a,Es10.3,a)') "NICIL: Mean molecular mass:             ",meanmolmass," m_proton"
 else
    if (ion_rays) then
       if (g_cnst) then
          write(iprint,'(a)')    "NICIL: Using constant grain size."
       else
          write(iprint,'(a)')    "NICIL: Approximating grain distribution with MRN distribution."
       endif
    endif
    write(iprint,'(a)')          "NICIL: Species  abundance   mass fraction"
    do j = iH,nelements
       if (mass_frac(j) > 0.01 .or. mass_frac(j) <= small) then
          write(iprint,'(3a,Es10.3,4x,F6.3)')   "NICIL: ",symj(j),"     ",abundancej(j),mass_frac(j)
       else
          write(iprint,'(3a,Es10.3,4x,Es10.3)') "NICIL: ",symj(j),"     ",abundancej(j),mass_frac(j)
       endif
    enddo
    write(iprint,'(a,Es10.3,a)')       "NICIL: Mean molecular mass of cold gas: ",meanmolmass," m_proton"
    if (ion_rays) then
       write(iprint,'(a,Es10.3,a)')    "NICIL: Light element ion mass for CR's: ",mass_HionR_mp," m_proton"
       write(iprint,'(a,Es10.3,a)')    "NICIL: Metallic ion mass for CR's:      ",mass_MionR_mp," m_proton"
       if (use_fdg_in) then
          write(iprint,'(a         )') "NICIL: Importing grain mass density from parent code"
       else
          write(iprint,'(a,Es10.3  )') "NICIL: Dust-to-gas fraction:            ",fdg
       endif
       if (.not. rho_is_rhogas) then
          write(iprint,'(a         )') "NICIL: rho_gas = rho_in*(1 - fdg)"
       endif
       if (g_cnst) then
          write(iprint,'(a,Es10.3,a)') "NICIL: Grain mass:                      ",mass_grain_mp," m_proton"
          write(iprint,'(a,Es10.3,a)') "NICIL: Grain radius:                    ",a0_grain," cm"
       else
          write(iprint,'(a,Es10.3,a)') "NICIL: Average grain mass:              ",mass_grain_mp/mass_proton_cgs," m_proton"
          write(iprint,'(a,Es10.3,a)') "NICIL: Minimum grain radius:            ",an_grain, " cm"
          write(iprint,'(a,Es10.3,a)') "NICIL: Average grain radius:            ",a01_grain," cm"
          write(iprint,'(a,Es10.3,a)') "NICIL: Maximum grain radius:            ",ax_grain, " cm"
       endif
       if (zeta_of_rho) then
          write(iprint,'(a,Es10.3,a)') "NICIL: Unattenuated cosmic ray ionisation rate:     ",zeta_CR_cgs," s^{-1}"
          write(iprint,'(a,Es10.3,a)') "NICIL: Ionisation rate from decaying radionuclides: ",zeta_R_cgs," s^{-1}"
       else
          write(iprint,'(a,Es10.3,a)') "NICIL: Cosmic ray ionisation rate:      ",zeta_cgs," s^{-1}"
       endif
       write(iprint,'(a)')             "NICIL: Including ionisation from Cosmic rays"
    endif
    if (ion_thermal) then
       write(iprint,'(a)')             "NICIL: Including thermal ionisation"
    endif
    if(use_ohm.or.use_hall .or.use_ambi) then
       write(iprint,'(2a)')            "NICIL: Non-ideal terms used: ",trim(ni_terms)
    endif
    if (.not.ion_thermal .and. .not.ion_rays) then
       write(iprint,'(a)')             "NICIL: WARNING! No ionisation sources included!"
       write(iprint,'(a)')             "NICIL: WARNING! No non-ideal MHD coefficients will be calculated!"
       write(iprint,'(a)')             "NICIL: WARNING! This is Ideal MHD!"
    endif
 endif
 !
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
 !
 !--Print cause of  error
 if (present(var) .and. present(val)) then
    write(iprint,'(7a,Es16.4)') 'NICIL: ERROR: ',wherefrom,': ',reason,'. ',var,'=',val
 else
    write(iprint,'(4a)') 'NICIL: ERROR: ',wherefrom,': ',reason
 endif
 num_errors = num_errors + 1
 !
end subroutine nicil_ic_error
!
!--Runtime error messages
subroutine nicil_translate_error(ierr)
 integer, intent(inout)    :: ierr
 integer                   :: i,ilen
 character(len= 1)         :: id
 character(len=16)         :: cerr,werr
 character(len=20)         :: ferr
 character(len=96)         :: errmsg(19)
 logical                   :: is_fatal,is_warning
 !
 write(cerr,*) ierr
 !--The error messages
 ferr       = 'NICIL: FATAL ERROR: '
 werr       = 'NICIL: WARNING: '
 errmsg( 1) = 'nicil_ionR_get_n: n_{i,e,g} did not converge'
 errmsg( 2) = 'nicil_ionT_get_ne: n_electronT did not converge'
 errmsg( 3) = 'nicil_ionT_get_ne: n_electronT < 0'
 errmsg( 4) = 'nicil_ionR_get_n: unrealistic n_R value (i.e. n_R < 0)'
 errmsg( 5) = 'nicil_ionR_get_n: invalid Z_g > 0 in average Z_grain method'
 errmsg( 6) = 'calling error: use_fdg_in=.true., but fdg_in not passed in'
 errmsg( 7) = 'input error: T == 0.  Verify your code can legitimately input this.'
 errmsg( 8) = 'input error: B == 0.  Verify your code can legitimately input this.'
 errmsg( 9) = 'Using average-grain-charge method to approximate densities for cosmic rays.'
 errmsg(10) = 'nicil_ion_get_sigma: Using Langevin cross section for molecular Hydrogen.'
 errmsg(11) = 'nicil_ion_get_sigma: Using Langevin cross section for atomic Hydrogen.'
 errmsg(12) = 'nicil_ion_get_sigma: Using Langevin cross section for Helium.'
 errmsg(13) = 'sigma_Hall == 0 (would otherwise be dominate by round-off errors)'
 errmsg(14) = 'strong coupling approximation (rho_n~ rho) invalid since rho_n < 0.1rho'
 errmsg(15) = 'strong coupling approximation (rho_i<<rho) invalid since rho_i > 0.01rho'
 errmsg(16) = 'drift velocity ~ 0 approximation invalid since v_d > 0.1v'
 errmsg(17) = 'vrms_min < v < vrms_max is not true'
 errmsg(18) = 'thermal/cosmic ray decoupling approximation invalid since (n_eI ~ n_eR)'
 errmsg(19) = 'gas is fully ionised'

 !--Initialise parameters
 ilen       = len(trim(cerr))
 is_warning = .false.
 is_fatal   = .false.
 !--Determine which errors to print
 do i = ilen,1,-1
    id = cerr(i:i)
    if (id/="" .and. id/="0") then
       if (i >= ilen-1) then
          is_fatal   = .true.
       else
          is_warning = .true.
       endif
       if (i == ilen) then
          if ( is_err(id,ierr_nRconv, 0) ) write(iprint, '(2a)') ferr,errmsg( 1)
          if ( is_err(id,ierr_neTconv,0) ) write(iprint, '(2a)') ferr,errmsg( 2)
          if ( is_err(id,ierr_neTle0, 0) ) write(iprint, '(2a)') ferr,errmsg( 3)
       else if (i == ilen-1) then
          if ( is_err(id,ierr_nRle0,  1) ) write(iprint, '(2a)') ferr,errmsg( 4)
          if ( is_err(id,ierr_Zge0,   1) ) write(iprint, '(2a)') ferr,errmsg( 5)
          if ( is_err(id,ierr_fdg_in, 1) ) write(iprint, '(2a)') ferr,errmsg( 6)
       endif
       if (warn_verbose) then
          if (i == ilen-2) then
             if ( is_err(id,ierr_T,    2) ) write(iprintw,'(2a)') werr,errmsg( 7)
             if ( is_err(id,ierr_B,    2) ) write(iprintw,'(2a)') werr,errmsg( 8)
             if ( is_err(id,ierr_Zave, 2) ) write(iprintw,'(2a)') werr,errmsg( 9)
          else if (i == ilen-3) then
             if ( is_err(id,ierr_LagH2,3) ) write(iprintw,'(2a)') werr,errmsg(10)
             if ( is_err(id,ierr_LagH, 3) ) write(iprintw,'(2a)') werr,errmsg(11)
             if ( is_err(id,ierr_LagHe,3) ) write(iprintw,'(2a)') werr,errmsg(12)
          else if (i == ilen-4) then
             if ( is_err(id,ierr_sigH, 4) ) write(iprintw,'(2a)') werr,errmsg(13)
             if ( is_err(id,ierr_scN,  4) ) write(iprintw,'(2a)') werr,errmsg(14)
             if ( is_err(id,ierr_scI,  4) ) write(iprintw,'(2a)') werr,errmsg(15)
          else if (i == ilen-5) then
             if ( is_err(id,ierr_drift,5) ) write(iprintw,'(2a)') werr,errmsg(16)
             if ( is_err(id,ierr_vrms, 5) ) write(iprintw,'(2a)') werr,errmsg(17)
             if ( is_err(id,ierr_nedec,5) ) write(iprintw,'(2a)') werr,errmsg(18)
          else if (i == ilen-6) then
             if ( is_err(id,ierr_alion,6) ) write(iprintw,'(2a)') werr,errmsg(19)
          endif
       endif
       !
    endif
 enddo
 !
 !--Change ierr < 0 if non-fatal warning
 if (is_warning .and. .not.is_fatal) ierr = -ierr
 !
end subroutine nicil_translate_error
!
! Function to state whether or not an error message should be printed
! Note: nzeros is the number of zeros in the error code; passed in so
! they can be removed and only the first digit dealt with
pure logical function is_err(cid,ierr_in,nzeros)
 integer,         intent(in)   :: ierr_in,nzeros
 character(len=1),intent(in)   :: cid
 integer                       :: id,ierrA,ierrB,ierrC
 !
 read(cid,*) id
 ierrA  = ierr_in/10**nzeros
 is_err = .false.
 if (ierrA==1) then
    ierrB = 2
    ierrC = 5
 elseif (ierrA==2) then
    ierrB = 5
    ierrC = 1
 elseif (ierrA==5) then
    ierrB = 1
    ierrC = 2
 else
    ierrB  = -1
    ierrC  = -1
    is_err = .true.
 endif
 if (id==ierrA               .or.  &
     id==ierrA+ierrB         .or.  &
     id==ierrA+ierrC         .or.  &
     id==ierrA+ierrB+ierrC) is_err = .true.
 !
end function is_err
!----------------------------------------------------------------------!
!+
! This is the primary control routine for NICIL.  It will call the
! required ionisation and non-ideal MHD routines to determine the
! coefficients of the non-ideal MHD terms
!+
!----------------------------------------------------------------------!
pure subroutine nicil_get_eta(eta_ohm,eta_hall,eta_ambi,Bfield,rho,T,n_R,n_electronT,ierr,data_out,fdg_in)
 integer,          intent(out) :: ierr
 real,             intent(out) :: eta_ohm,eta_hall,eta_ambi
 real,             intent(in)  :: Bfield,rho,T,n_electronT,n_R(:)
 real,   optional, intent(in)  :: fdg_in
 real,   optional, intent(out) :: data_out(n_data_out)
 integer                       :: j
 real                          :: mass_neutral_mp,n_total,n_cold,n_electronR,rho_gas
 real                          :: zeta,fdg_local,n_cold_total,n_gas,f_Bdust
 real                          :: n_R_complete(nimass+2*na)
 real                          :: n_ionR(nimass),n_grainR(2*na),n_ionT(nlevels),mass_ionT(nlevels)
 real                          :: sigmas(8),n_densities(7),afrac(2),mfrac(2)
 real                          :: njk(nlevels,nelements)
 logical                       :: get_data_out
 !
 ierr = 0                                  ! initialise error code
 !
 if (.not.eta_constant) then
    if (T > small .and. Bfield > small) then
       !--Determine if we will be returning bookkeeping data
       if (present(data_out)) then
          get_data_out = .true.
       else
          get_data_out = .false.
       endif
       !--Determine actual gas density
       if (.not.use_fdg_in) then
          fdg_local = fdg
       else
          if (.not. present(fdg_in)) then
             ierr = ierr + ierr_fdg_in
             return
          endif
          fdg_local = min(fdg_in,   fdg_max)
          fdg_local = max(fdg_local,fdg_min)
       endif
       if (rho_is_rhogas) then
          rho_gas = rho
       else
          if (.not.use_fdg_in) then
             rho_gas = rho*one_minus_fdg
          else
             rho_gas = rho*(1.0-fdg_local)
          endif
       endif
       !
       !--Determine the fractions of molecular and atomic Hydrogen
       n_cold       = rho_gas*mass_neutral_cold1
       n_cold_total = rho*mass_neutral_cold1
       call nicil_get_HH2_ratio(abundancej(iH),n_cold,T,mass_neutral_mp,afrac,mfrac)
       n_gas        = rho_gas*mass_proton1/mass_neutral_mp
       n_total      = rho*mass_proton1/mass_neutral_mp
       !
       !--Calculate electron number densities from cosmic rays
       if (ion_rays) then
          if (.not.zeta_of_rho) then
             zeta = zeta0
          else
             zeta = zetaCR*exp(-zetaCRcoef*sqrt(T*rho_gas/mass_neutral_mp)) + zetaR
          endif
          if (size(n_R)==nimass+2*na .and. coef2_epsilonR*max(n_R(1+nimass),n_R(2+nimass)) < min(n_R(1),n_R(2))) then
             n_grainR    = n_R(1+nimass:nimass+2*na)
             n_electronR = nicil_ionR_get_ne(n_R)
             n_ionR      = n_R(1:nimass)
          else
             ! must re-calculate grain densities since they were not stored
             call nicil_ionR_predict_ng(n_R_complete,n_R(1:nimass+2))
             call nicil_ionR_get_n(n_R_complete,T,rho_gas,n_cold_total,fdg_local,zeta,n_electronR,f_Bdust,ierr)
             n_ionR      = n_R(1:nimass)
             do j = 1,na
                n_grainR(2*j-1) = n_R_complete(nimass+2*j-1)
                n_grainR(2*j  ) = n_R_complete(nimass+2*j  )
             enddo
          endif
       else
          n_electronR = 0.0
          n_ionR      = 0.0
          n_grainR    = 0.0
       endif
       !
       !--Calculate ion number densities from thermal ionisation
       if (ion_thermal) then
          call nicil_ionT_get_n(n_ionT,mass_ionT,njk,n_electronT,T,n_gas,afrac)
       else
          n_ionT      = 0.0
          mass_ionT   = 0.0
          njk         = 0.0
       endif
       if (warn_verbose) then
          ! ensure that these two processes are sufficiently decoupled
          if (warn_ratio_m1*n_electronR < n_electronT .and. n_electronT < warn_ratio_p1*n_electronR) then
             ierr = ierr + ierr_nedec
          endif
       endif
       !
       !--Calculate the conductivities
       call nicil_ion_get_sigma(Bfield,rho_gas,n_electronR,n_electronT,n_ionR,n_grainR,n_ionT,mass_ionT, &
                                mass_neutral_mp,T,sigmas,mfrac,get_data_out,n_densities,ierr)
       !
       !--Calculate the coefficients
       call nicil_nimhd_get_eta(eta_ohm,eta_hall,eta_ambi,sigmas)
       !
       !--Copy optional arrays to output, if requested
       if (present(data_out)) then
          data_out        = 0.0
          data_out(    1) = sigmas(1)                                             ! Ohmic conductivities
          data_out(    2) = sigmas(5)                                             ! Hall conductivities
          data_out(    3) = sigmas(2)                                             ! Pedersen conductivities
          data_out( 4: 5) = n_densities(1:2)                                      ! rho_neutral, rho_ion (total)
          data_out(    6) = n_electronR + n_electronT
          data_out( 7:11) = n_densities(3:7)                                      ! n_neutral, n_ionR (light,metallic), n_ionT(singly,doubly)
          do j = 1,na
             data_out(12)  = data_out(12) + n_grainR(2*j-1)                       ! negatively charged grains
             if (.not.use_fdg_in) then
                data_out(13) = data_out(13) + n_grain_coef(j)*n_cold_total           ! total grains
             else
                data_out(13) = data_out(13) + n_grain_coef(j)*n_cold_total*fdg_local ! total grains
             endif
             data_out(14)    = data_out(14) + n_grainR(2*j  )                     ! positively  charged grains
          enddo
          data_out(13)    = data_out(13) - data_out(12) - data_out(14)            ! convert total grains to neutral grains
          data_out(15:16) = afrac*n_total                                         ! molecular & atomic hydrogen
          data_out(17:17+nelements-1) = njk(1,:)                                  ! singly ionised elements
          data_out(17+nelements:17+nelements*nlevels-3) = njk(2,iHe:nelements)    ! doubly ionised elements
          if (ion_rays .and. zeta_of_rho) data_out(17+nelements*nlevels-2) = zeta ! density dependent cosmic ray ionisation rate
       endif
   else
       !--Exit with error message if T = 0 or B = 0
       if ( T      <=small ) ierr = ierr + ierr_T
       if ( Bfield <=small ) ierr = ierr + ierr_B
       eta_ohm  = 0.0
       eta_hall = 0.0
       eta_ambi = 0.0
       if (present(data_out)) data_out = 0.0
       return
    endif
 else
    !--Return constant coefficient version and exit
    call nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
    if (present(data_out)) data_out = 0.0  ! bookkeeping has no meaning here
    return
 endif
 !
end subroutine nicil_get_eta
!----------------------------------------------------------------------!
!+
!  calculates the condictivities
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
!  The n_densities arrays is for bookkeeping only, and has no real
!  computational value.
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ion_get_sigma(Bmag,rho_gas,n_electronR,n_electronT,n_ionR,n_grainR,n_ionT,mass_ionT, &
                                    mass_neutral_mp,T,sigmas,mfrac,get_data_out,n_densities,ierr)
 integer,intent(inout) :: ierr
 real,   intent(out)   :: sigmas(8),n_densities(7)
 real,   intent(in)    :: Bmag,rho_gas,T,n_electronR,n_electronT,mass_neutral_mp
 real,   intent(in)    :: n_ionR(:),n_grainR(:),n_ionT(:),mass_ionT(:),mfrac(2)
 logical,intent(in)    :: get_data_out
 integer               :: j,k,p
 real                  :: sqrtT,logT,sigmavenXH2,sigmavenXH,sigmavenY,sigma_coef_onB
 real                  :: rho_n,rho_i,nu_ei,sigmaviTn,n_electron
 real                  :: mu_iXH21,mu_iXH1,mu_iY1
 real                  :: sigmaviRntot(nimass),mass_ionT_mp(nlevels),rho_ion(2)
 real                  :: nu_jn(nspecies_max),ns(nspecies_max),rho_j(nspecies_max)
 real                  :: betaj(nspecies_max),beta2p11(nspecies_max)
 !
 !--Initialise values
 sqrtT          = sqrt(T)
 logT           = log10(T)
 sigma_coef_onB = sigma_coef/Bmag
 sigmas         = 0.0
 !
 !--Number densities
 ns               = 0.0
 n_electron       = n_electronR + n_electronT
 ns(ine)          = n_electron
 ns(iniHR:iniMR)  = n_ionR
 ns(inisT:inidT)  = n_ionT
 do j = 1,na
    p = 1+inidT+2*(j-1)
    ns(p:p+1)      = n_grainR(1+2*(j-1):2+2*(j-1))
 enddo
 if (n_electron > 0.0) then
    ! to prevent (possible) numerical overflow
    ns             = ns/n_electron
    sigma_coef_onB = sigma_coef_onB*n_electron
 endif
 !
 !--Densities
 rho_j(ine)         = n_electron*massj(ine)
 rho_j(iniHR:iniMR) = n_ionR*massj(iniHR:iniMR)
 rho_j(inisT:inidT) = n_ionT*mass_ionT
 rho_ion            = 0.
 do j = 1,nimass
    rho_ion(1)      = rho_ion(1) + rho_j(iniHR+j-1)
 enddo
 do k = 1,nlevels
    rho_ion(2)      = rho_ion(2) + rho_j(inisT+k-1)
 enddo
 rho_n = rho_gas - rho_j(ine) - rho_ion(1) - rho_ion(2)
 !
 if (rho_n < small*rho_gas) then
    ! In this regime, ionisation from both thermal ionisation and cosmic rays are important,
    ! and lead to very small or negative neutral gas densities; Scale contribution to prevent this
    rho_n = rho_gas - rho_j(ine) - (rho_ion(1)*n_electronR + rho_ion(2)*n_electronT)/n_electron
 endif
 !
 if (warn_verbose) then
    ! ensure that the strong coupling approximation is valid
    if (warn_ratio_m1*rho_gas > rho_n .or. rho_n > warn_ratio_p1*rho_gas) ierr = ierr + ierr_scN
    rho_i = rho_gas - rho_n - rho_j(ine)
    if (warn_ratio2*rho_i > rho_gas) ierr = ierr + ierr_scN
 endif
 !
 if ( rho_n > 0.0 ) then
    !--Rate Coefficients; use Langevin rate if high temperature
    sigmavenXH2 = sigmavenX_coef*sqrtT*(0.535 + (0.203 - (0.163 - 0.050*logT)*logT)*logT)
    sigmavenXH  = sigmavenX_coef*sqrtT*(2.841 + (0.093 - (0.245 - 0.089*logT)*logT)*logT)
    sigmavenY   = sigmavenY_coef*sqrtT
    if (sigmavenXH2 > sigmavenXH2_LR) then
       sigmavenXH2 = sigmavenXH2_LR
       if (warn_verbose) ierr = ierr + ierr_LagH2
    endif
    if (sigmavenXH > sigmavenXH_LR) then
       sigmavenXH = sigmavenXH_LR
       if (warn_verbose) ierr = ierr + ierr_LagH
    endif
    if (sigmavenY > sigmavenY_LR) then
       sigmavenY = sigmavenY_LR
       if (warn_verbose) ierr = ierr + ierr_LagHe
    endif
    sigmavenXH2 = mfrac(1)*sigmavenXH2
    sigmavenXH  = mfrac(2)*sigmavenXH
    !
    !--Collisional Frequencies
    sigmaviRntot          = mfrac(1)*sigmaviRn(1,1:nimass) + mfrac(2)*sigmaviRn(2,1:nimass) + sigmaviRn(3,1:nimass)
    nu_ei                 = nu_ei_coef           *n_electron/sqrtT**3
    nu_jn                 = 0.0
    nu_jn(ine:iniMR)      = nu_jn_coef(ine:iniMR)*rho_n/( mass_neutral_mp + massj_mp(ine:iniMR) )
    nu_jn(ine)            = nu_jn(ine)           *(sigmavenXH2 + sigmavenXH + sigmavenY)
    nu_jn(iniHR:iniMR)    = nu_jn(iniHR:iniMR)   *sigmaviRntot
    nu_jn(ing:ing+2*na-1) = nu_jn_coef(ing:ing+2*na-1)/sqrt(mass_neutral_mp)*sqrtT &
                          * rho_n/( mass_neutral_mp + massj_mp(ing:ing+2*na-1) )
    if ( ion_thermal ) then
       mass_ionT_mp = mass_ionT*mass_proton1
       do k = 1,nlevels
          if (mass_ionT(k) > 0.0) then
             mu_iXH21     = (mj_mp(iH2) + mass_ionT_mp(k))/(mj_mp(iH2) * mass_ionT_mp(k))
             mu_iXH1      = (mj_mp(iH ) + mass_ionT_mp(k))/(mj_mp(iH ) * mass_ionT_mp(k))
             mu_iY1       = (mj_mp(iHe) + mass_ionT_mp(k))/(mj_mp(iHe) * mass_ionT_mp(k))
             sigmaviTn    = mfrac(1)*sqrt(pnH2 * mu_iXH21) + mfrac(2)*sqrt(pnH * mu_iXH1) + massfrac_Y*sqrt(pnHe * mu_iY1)
             nu_jn(inisT+k-1) = nu_jn_coef(inisT+k-1)*sigmaviTn*rho_n/( mass_neutral_mp + mass_ionT_mp(k) )
          endif
       enddo
    endif
    !
    !--Hall parameters
    betaj        = 0.0
    beta2p11     = 0.0
    if (.not. mod_beta) then
       betaj(ine) = beta_coef(ine) * Bmag/nu_jn(ine)
    else
       betaj(ine) = beta_coef(ine) * Bmag/(nu_jn(ine)  + nu_ei)
    endif
    betaj(iniHR) = betaIj(beta_coef(iniHR),rho_j(iniHR),rho_j(ine),nu_jn(iniHR),nu_ei,Bmag,ion_rays)
    betaj(iniMR) = betaIj(beta_coef(iniMR),rho_j(iniMR),rho_j(ine),nu_jn(iniMR),nu_ei,Bmag,ion_rays)
    betaj(inisT) = betaIj(beta_coef(inisT),rho_j(inisT),rho_j(ine),nu_jn(inisT),nu_ei,Bmag,ion_thermal,mass_ionT(1))
    betaj(inidT) = betaIj(beta_coef(inidT),rho_j(inidT),rho_j(ine),nu_jn(inidT),nu_ei,Bmag,ion_thermal,mass_ionT(2))
    betaj(ing:ing+2*na-1) = beta_coef(ing:ing+2*na-1) *Bmag/ nu_jn(ing:ing+2*na-1)
    beta2p11     = 1.0/( 1.0+betaj**2 )
    !
    !--Conductivities
    do j = 1,nspecies
       sigmas(1) =    sigmas(1) + ns(j)*aZj(j)*betaj(j)              ! sigma_O
       sigmas(2) =    sigmas(2) + ns(j)*aZj(j)*betaj(j)*beta2p11(j)  ! sigma_P
       if (Zj(j) > 0.0) then
          sigmas(3) = sigmas(3) + ns(j)* Zj(j)         *beta2p11(j)  ! sigma_H > 0
       else
          sigmas(4) = sigmas(4) + ns(j)* Zj(j)         *beta2p11(j)  ! sigma_H < 0
       endif
    enddo
    sigmas(5) = sigmas(3)+sigmas(4)                                  ! sigma_H
    if ( abs(sigmas(5)) < coef_epsilonR*max(sigmas(3),-sigmas(4))) then
       sigmas(5) = 0.0                                               ! sigma_H (alternative)
       if (warn_verbose) ierr = ierr + ierr_sigH
    endif
    sigmas = sigmas*sigma_coef_onB                                   ! scale sigma
    if (sigmas(1) > 0.0) sigmas(6) = 1.0/sigmas(1)                   ! 1/sigma_O
    sigmas(7)   = sigmas(2)*sigmas(2) + sigmas(5)*sigmas(5)          ! perp^2 = P^2 + H^2
    sigmas(8)   = sigmas(1)*sigmas(2) - sigmas(7)                    ! OP - perp^2
    if ( sigmas(8) < coef_epsilonR*max(sigmas(1)*sigmas(2),sigmas(7))) then
       sigmas(8) = 0.0
       do j = 1,nspecies
          do k = j+1,nspecies
             sigmas(8) = sigmas(8) + ns(j)*aZj(j)*betaj(j)*beta2p11(j) &
                                   * ns(k)*aZj(k)*betaj(k)*beta2p11(k) &
                                   * (sZj(j)*betaj(j) - sZj(k)*betaj(k))**2
          enddo
       enddo
       sigmas(8) = sigmas(8)*sigma_coef_onB**2
    endif
    if ( sigmas(7) > 0.0 ) sigmas(7) = 1.0/sigmas(7)                 ! perp^2 -> 1/perp^2
 else
    ! This is the Ideal MHD regime.  Turn off non-ideal terms.
    if (warn_verbose) ierr = ierr + ierr_alion
    sigmas = -1.0
    rho_n  =  0.0
 endif
 !
 !--Number densities (for bookkeeping)
 n_densities = 0.0
 if ( get_data_out ) then
    n_densities(1) = rho_n                                ! neutral density
    n_densities(2) = rho_j(iniHR)+rho_j(iniMR) &
                   + rho_j(inisT)+rho_j(inidT)            ! total ion density
    n_densities(3) = rho_n*mass_proton1/mass_neutral_mp   ! neutral number density
    n_densities(4) = n_ionR(1)                            ! light element ion number density from cosmic rays
    n_densities(5) = n_ionR(2)                            ! metallic ion number density from cosmic rays
    n_densities(6) = n_ionT(1)                            ! singly ionised ion number density from thermal radiation
    n_densities(7) = n_ionT(2)                            ! doubly ion number density from thermal radiation
 endif
 !
end subroutine nicil_ion_get_sigma
!
pure real function betaIj(betaj_coef,rhoj,rhoe,nu_jn,nu_ei,Bmag,calc_beta,mass)
 ! Function to calculate the Hall parameter for ions
 ! Recall: nu_ie = nu_ei*rho_e/rho_i
 real,    intent(in)           :: betaj_coef,rhoj,rhoe,nu_jn,nu_ei,Bmag
 real,    intent(in), optional :: mass
 logical, intent(in)           :: calc_beta
 !
 betaIj = 0.0
 if ( calc_beta .and. rhoj > 0.0 .and. nu_jn > 0.0) then
    if (.not.mod_beta) then
       betaIj = betaj_coef*Bmag/ nu_jn
    else
       betaIj = betaj_coef*Bmag/(nu_jn + nu_ei*rhoe/rhoj)
    endif
    if (present(mass)) then
       if (mass > 0.0) betaIj = betaIj/mass
    endif
 endif
 !
end function betaIj
!----------------------------------------------------------------------!
!+
! This is a control routine for NICIL.  It will update number grain
! densities.  These are calculated iteratively, thus will be done prior
! to calculating resistivities.
!+
!----------------------------------------------------------------------!
pure subroutine nicil_get_ion_n(rho,T,n_R,n_electronT,ierr,fdg_in,f_Bdust_out)
 real,              intent(in)    :: rho,T
 real,              intent(inout) :: n_R(:),n_electronT
 integer,           intent(out)   :: ierr
 real,    optional, intent(in)    :: fdg_in
 real,    optional, intent(out)   :: f_Bdust_out
 real                             :: mass_neutral_mp,n_total,n_cold,n_cold_total,n_electronR,rho_gas
 real                             :: zeta,fdg_local,f_Bdust
 real                             :: afrac(2),n_R_complete(nimass+2*na)
 !
 ierr = 0
 !--Exit if using constant resistivities
 if (eta_constant) return

 if (T > small) then
    !--If requesting dust mass density, ensure it has been properly passed
    if (.not. use_fdg_in) then
       fdg_local = fdg
    else
       if (.not. present(fdg_in)) then
          ierr = ierr + ierr_fdg_in
          return
       endif
       fdg_local = min(fdg_in,   fdg_max)
       fdg_local = max(fdg_local,fdg_min)
    endif
    if (rho_is_rhogas) then
       rho_gas = rho
    else
       if (.not.use_fdg_in) then
          rho_gas = rho*one_minus_fdg
       else
          rho_gas = rho*(1.0-fdg_local)
       endif
    endif
 else
    !--Exit with error message if T = 0
    ierr        = ierr + ierr_T
    n_R         = 0.0
    n_electronT = 0.0
    return
 endif
 !
 !--Determine the ratio of H and H2
 n_cold       = rho_gas*mass_neutral_cold1
 n_cold_total = rho*mass_neutral_cold1
 call nicil_get_HH2_ratio(abundancej(iH),n_cold,T,mass_neutral_mp,afrac)
 n_total      = rho_gas*mass_proton1/mass_neutral_mp
 !
 !--Calculate the number densities from cosmic rays
 if (ion_rays) then
    if (.not.zeta_of_rho) then
       zeta = zeta0
    else
       zeta = zetaCR*exp(-zetaCRcoef*sqrt(T*rho_gas/mass_neutral_mp)) + zetaR
    endif
    call nicil_fill_nRcomplete(n_R_complete,n_R,n_total)
    n_R_complete = epsilon(n_R_complete(1))
    call nicil_ionR_get_n(n_R_complete,T,rho_gas,n_cold_total,fdg_local,zeta,n_electronR,f_Bdust,ierr)
    call nicil_unfill_nRcomplete(n_R_complete,n_R)
    if (present(f_Bdust_out)) f_Bdust_out = f_Bdust
 endif
 !
 !--Calculate the grain charge from thermal ionisaion and update electron number density
 if (ion_thermal) then
    if (n_electronT <= small) n_electronT = small ! special consideration for n_electronT ~ 0
    call nicil_ionT_get_ne(n_electronT,n_total,T,afrac,ierr)
 endif
 !
end subroutine nicil_get_ion_n
!----------------------------------------------------------------------!
!+
! Modify the n_R and n_R_complete arrays
!+
!----------------------------------------------------------------------!
!--Convert n_R -> n_R_complete since they may be different sizes
pure subroutine nicil_fill_nRcomplete(n_R_complete,n_R,n0)
 real, intent(out) :: n_R_complete(:)
 real, intent(in)  :: n0,n_R(:)
 !
 n_R_complete(1:nimass) = n_R(1:nimass)
 if (size(n_R) == nimass+2*na) then
    n_R_complete(1+nimass:nimass+2*na) = n_R(1+nimass:nimass+2*na)
 else
    call nicil_ionR_predict_ng(n_R_complete,n_R)
 endif
 call nicil_ionR_set_guess(n_R_complete,n0,.true.)
 !
end subroutine nicil_fill_nRcomplete
!
!--Convert n_R_complete -> n_R since they may be different sizes
pure subroutine nicil_unfill_nRcomplete(n_R_complete,n_R)
 real, intent(in)  :: n_R_complete(:)
 real, intent(out) :: n_R(:)
 integer           :: j
 !
 n_R(1:nimass) = n_R_complete(1:nimass)
 if (size(n_R) == nimass+2*na) then
    n_R(1+nimass:nimass+2*na) = n_R_complete(1+nimass:nimass+2*na)
 else
    n_R(1+nimass) = 0.0
    n_R(2+nimass) = 0.0
    do j = 1,na
       n_R(1+nimass) = n_R(1+nimass) + n_R_complete(1+nimass+2*(j-1))
       n_R(2+nimass) = n_R(2+nimass) + n_R_complete(2+nimass+2*(j-1))
    enddo
 endif
 !
end subroutine nicil_unfill_nRcomplete
!
!--Sets the initial guesses, either for the entire array, or for selected values
pure subroutine nicil_ionR_set_guess(n_R,n0,conditional)
 real,    intent(inout) :: n_R(:)
 real,    intent(in)    :: n0
 logical, intent(in)    :: conditional
 integer                :: i
 !
 if (conditional) then
    do i = 1,size(n_R)
       if (n_R(i) <= small) n_R(i) = n_R0(i)*n0
    enddo
 else
    n_R = n_R0(1:size(n_R))*n0
 endif
 !
end subroutine nicil_ionR_set_guess
!
! If using the MRN grain distribution, then predict the grain number
! densities, assuming only a single grain density per charge is stored.
pure subroutine nicil_ionR_predict_ng(n_R_complete,n_R)
 real, intent(out) :: n_R_complete(:)
 real, intent(in)  :: n_R(:)
 integer           :: j
 !
 n_R_complete(1:nimass) = n_R(1:nimass)
 do j = 1,na
    n_R_complete(nimass+2*j-1) = n_grain_dist(j)*n_R(1+nimass)
    n_R_complete(nimass+2*j  ) = n_grain_dist(j)*n_R(2+nimass)
 enddo
 !
end subroutine nicil_ionR_predict_ng
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
pure subroutine nicil_ionR_get_n(n_R,T,rho_gas,n_total,fdg_local,zeta,n_e_out,f_Bdust,ierr)
 integer, intent(inout)  :: ierr
 real,    intent(inout)  :: n_R(nimass+2*na)
 real,    intent(in)     :: T,rho_gas,n_total,fdg_local,zeta
 real,    intent(out)    :: n_e_out,f_Bdust
 integer                 :: i,j,iter
 integer                 :: feqni(neqn)
 real                    :: n_e,n_g_charged,n_g_neutral
 real                    :: n_g_tot(na),n_i(nimass),n_g(-zmax:zmax,na)
 real                    :: n_old(neqn),n_new(neqn),dn(neqn)
 real                    :: k_ig(-zmax:zmax,nimass,na),k_eg(-zmax:zmax,na),k_ei(nimass)
 real                    :: feqn(neqn),Jacob(neqn,neqn)
 logical                 :: iterate,try_new_n0,lerr
 !
 !--Determine the rate coefficients given the current temperature
 call nicil_ionT_calc_k(T,k_ig,k_eg,k_ei)
 !
 !--Set initial conditions & define the nold array
 if (.not.use_fdg_in) then
    n_g_tot  = n_grain_coef(1:na)*n_total             ! fdg is included in n_grain_coef
 else
    n_g_tot  = n_grain_coef(1:na)*n_total*fdg_local
 endif
 iter       = 0
 iterate    = .true.
 try_new_n0 = .true.
 n_old(1:nimass) = n_R(1:nimass)
 n_old(1+nimass) = nicil_ionR_get_ne(n_R)
 do j = 1,na
    n_old(2+nimass+3*(j-1)) = n_R(nimass+2*j-1)
    n_old(3+nimass+3*(j-1)) = n_g_tot(j) - (n_R(nimass+2*j-1) + n_R(nimass+2*j))
    n_old(4+nimass+3*(j-1)) = n_R(nimass+2*j  )
 enddo
 !
 !--Perform the iterations
 do while (iterate)
    !
    ! update the number densities
    n_i = n_old(1:nimass)
    n_e = n_old(1+nimass)
    do j = 1,na
       n_g(:,j) = n_old(2+nimass+3*(j-1):4+nimass+3*(j-1))
    enddo
    !
    ! update the system of equations and the Jacobian
    call nicil_ionT_calc_Jf(feqni,feqn,Jacob,rho_gas,n_g_tot,n_e,n_i,n_g,k_ig,k_eg,k_ei,zeta)
    ! determine dn by solving Jacob*dn = feqn
    call get_dn_viaLU(dn,feqni,feqn,Jacob,lerr)
    if (.not.lerr) then
       ! calculate the new number densities
       n_new = n_old - dn
       !
       ! determine if we need to iterate again
       iter    = iter + 1
       iterate = .false.
       do i = 1,neqn
          ! first term: n_rat > NRtol, where n_rat = abs( 1.0 - n_new(i)/n_old(i) ), where n_old >= 0
          if ( (abs(n_old(i)-n_new(i)) > n_old(i)*NRtol) .or. n_new(i) < 0.0) then
             n_new(i) = abs(n_new(i))
             iterate = .true.
          endif
       enddo
       n_old = n_new
       ! Values successfully converged
       if (.not. iterate) then
          n_R(1:nimass) = n_new(1:nimass)
          n_g_charged = 0.0
          n_g_neutral = 0.0
          do j = 1,na
             n_R(nimass+2*j-1) = n_new(2+nimass+3*(j-1))
             n_R(nimass+2*j  ) = n_new(4+nimass+3*(j-1))
             n_g_charged = n_g_charged + n_new(2+nimass+3*(j-1)) + n_new(4+nimass+3*(j-1))
             n_g_neutral = n_g_neutral + n_new(3+nimass+3*(j-1))
          enddo
          f_Bdust = n_g_charged/(n_g_charged+n_g_neutral)
       endif
       !
       ! Failsafe to prevent an infinite loop
       if (iter > NRctrmax) iterate = .false.
    else
       iterate = .false.
    endif
    ! Iterations unsuccessful; take appropriate actions
    if (iter > NRctrmax .or. lerr) then
       if (try_new_n0) then
          ! Try again with default guess rather than value from previous iteration
          call nicil_ionR_set_guess(n_old,n_total,.false.)
          try_new_n0 = .false.
          iterate    = .true.
          iter       = 0
       else
          ! New guess failed; Switch Methods
          iterate  = .false.
          call nicil_ionR_get_n_via_Zave(n_R,n_e,rho_gas*mass_neutral_cold1,n_g_tot,T,zeta,k_ig,k_eg,f_Bdust,ierr)
       endif
    endif
 enddo
 n_e_out = n_e
 !
end subroutine nicil_ionR_get_n
!----------------------------------------------------------------------!
!+
! Calculate the capture Rates coefficients now that temperature is known
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_calc_k(T,k_ig,k_eg,k_ei)
 real, intent(in)  :: T
 real, intent(out) :: k_ig(-zmax:zmax,nimass,na),k_eg(-zmax:zmax,na),k_ei(nimass)
 integer           :: j,k
 real              :: sqrtT
 real              :: k_fac(-zmax:zmax,na)
 !
 k_fac(-zmax:zmax,1:na) = k_fac_coef(-zmax:zmax,1:na)/T
 sqrtT = sqrt(T)
 do j = 1,na
    do k = 1,nimass
       k_ig(-1,k,j) = k_ig_coef(k,j)*sqrtT*(1.0-k_fac(-1,j))
       k_ig( 0,k,j) = k_ig_coef(k,j)*sqrtT
       k_ig( 1,k,j) = k_ig_coef(k,j)*sqrtT*exp(-k_fac( 1,j))
    enddo
    k_eg(-1,j) = k_eg_coef(j)*sqrtT*exp(k_fac(-1,j))
    k_eg( 0,j) = k_eg_coef(j)*sqrtT
    k_eg( 1,j) = k_eg_coef(j)*sqrtT*(1.0+k_fac( 1,j))
 enddo
 k_ei(1) = k_ei_coef(1)*T**alpha_expT_H + k_ei_coef(2)*T**alpha_expT_He
 k_ei(2) = k_ei_coef(3)*T**alpha_expT_Mg
 !
end subroutine nicil_ionT_calc_k
!----------------------------------------------------------------------!
!+
! Calculate the Jacobian and Rate equations
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_calc_Jf(feqni,feqn,Jacob,rho_gas,n_g_tot,n_e,n_i,n_g,k_ig,k_eg,k_ei,zeta)
 integer, intent(out)   :: feqni(neqn)
 real,    intent(out)   :: feqn(neqn),Jacob(neqn,neqn)
 real,    intent(in)    :: rho_gas,n_g_tot(na),n_g(-zmax:zmax,na)
 real,    intent(inout) :: n_e,n_i(nimass)
 real,    intent(in)    :: k_ig(-zmax:zmax,nimass,na),k_eg(-zmax:zmax,na),k_ei(nimass)
 real,    intent(in)    :: zeta
 integer                :: i,j,k,p,q,z,dcol,imax,jmax,feqnitmp(neqn),arrayTrans(neqn,neqn)
 real                   :: n_n,zcharge,Jacobmax
 real                   :: feqntmp(neqn),Jacobtmp(neqn,neqn),Jacobrow(neqn)
 !
 !--Neutral number density
 !  since this is an intermediate step, its acceptable to be < 0
 !  using cold gas mass since this process is dominant at cold temperatures
 n_n = (rho_gas - n_e*massj(ine) - n_i(1)*massj(iniHR) - n_i(2)*massj(iniMR))*mass_neutral_cold1
 !
 !--Zero the arrays
 feqn  = 0.0
 Jacob = 0.0
 !--Rate equations & Jacobian
 !  feqn==0 when all the values are correct; note grains are not couples to grains of different sizes
 !  Jacobian: Jacob(p,q)
 !         p: feqn for n_i(1:nimass), n_e, n_g(-zmax:zmax,1),n_g(-zmax:zmax,2),n_g(-zmax:zmax,3)...
 !         q: derivative of the matching primary variable of feqn
 !
 !  Ions
 do k = 1,nimass
    feqn(k) = zeta*n_n - k_ei(k)*n_e*n_i(k)                                               ! ion-electron
    do j = 1,na
       do z = -zmax,zmax
          feqn(k) = feqn(k) - k_ig(z,k,j)*n_g(z,j)*n_i(k)                                 ! ion-grains
       enddo
    enddo
    Jacob(k,1:nimass) = -zeta*massj(iniHR:iniMR)*mass_neutral_cold1                       ! df(n_i)/dn_i
    Jacob(k,k)        = Jacob(k,k) - k_ei(k)*n_e                                          ! df(n_i)/dn_i
    do j = 1,na
       do z = -zmax,zmax
          Jacob(k,k) = Jacob(k,k) - k_ig(z,k,j)*n_g(z,j)                                  ! df(n_i)/dn_i
       enddo
    enddo
    q = nimass + 1
    Jacob(k,q) = -zeta*massj(ine)*mass_neutral_cold1 - k_ei(k)*n_i(k)                     ! df(n_i)/dn_e
    do j = 1,na
       do z = -zmax,zmax
          q = q + 1
          Jacob(k,q) = -k_ig(z,k,j)*n_i(k)                                                ! df(n_i)/dn_g
       enddo
    enddo
 enddo
 !  Electrons using charge neutrality
 p = 1 + nimass
 feqn(p) = -n_e                                                                           ! charge neutrality
 do k = 1,nimass
    feqn(p) = feqn(p) + n_i(k)                                                            ! charge neutrality
 enddo
 do j = 1,na
    zcharge = n_g(1,j) - n_g(-1,j)
    if (abs (zcharge) > 1.0d12*n_e) then
       if (abs(n_g(1,j) - n_g(-1,j)) < 1.0d-8*max(n_g(-1,j),n_g(1,j))) zcharge = 0.0
    endif
    feqn(p) = feqn(p) + zcharge                                                           ! charge neutrality
 enddo
 Jacob(p,1:nimass) =  1.0                                                                 ! df(charge_neutrality)/dn_i
 Jacob(p,1+nimass) = -1.0                                                                 ! df(charge_neutrality)/dn_e
 do j = 1,na
    Jacob(p,2+nimass+3*(j-1)) = -1.0                                                      ! df(charge_neutrality)/dn_g-
    Jacob(p,3+nimass+3*(j-1)) =  0.0                                                      ! df(charge_neutrality)/dn_g0
    Jacob(p,4+nimass+3*(j-1)) =  1.0                                                      ! df(charge_neutrality)/dn_g+
 enddo
 p = p + 1
 !  Grains (Z = -1,0,+1)
 do j = 1,na
    dcol      = nimass+3*(j-1)
    feqn(p  ) = ( -k_eg(-1,j)*n_g(-1,j) + k_eg(0,j)*n_g(0,j) )*n_e                        ! n_g(Z = -1)
    feqn(p+1) = n_g_tot(j)                                                                ! total grain number density
    feqn(p+2) =   -k_eg( 1,j)*n_g( 1,j)                       *n_e                        ! n_g(Z =  1)
    do k = 1,nimass
       feqn(p  ) = feqn(p  )                            - k_ig(-1,k,j)*n_g(-1,j)  *n_i(k) ! n_g(Z = -1)
       feqn(p+2) = feqn(p+2) + ( k_ig( 0,k,j)*n_g( 0,j) - k_ig( 1,k,j)*n_g( 1,j) )*n_i(k) ! n_g(Z =  1)
    enddo
    do z = -zmax,zmax
       feqn(p+1) = feqn(p+1) - n_g(z,j)                                                   ! total grain number density
    enddo
    ! grains (Z = -1)
    Jacob(p,1:nimass) = -k_ig(-1,1:nimass,j)*n_g(-1,j)                                    ! df(n_g-)/dn_i
    Jacob(p,1+nimass) = -k_eg(-1,j)*n_g(-1,j) + k_eg(0,j)*n_g(0,j)                        ! df(n_g-)/dn_e
    Jacob(p,2+dcol)   = -k_eg(-1,j)*n_e                                                   ! df(n_g-)/dn_g-
    do k = 1,nimass
       Jacob(p,2+dcol) = Jacob(p,2+dcol) - k_ig(-1,k,j)*n_i(k)                            ! df(n_g-)/dn_g-
    enddo
    Jacob(p,3+dcol)   = k_eg(0,j)*n_e                                                     ! df(n_g-)/dn_g0
    Jacob(p,4+dcol)   = 0.0                                                               ! df(n_g-)/dn_g+
    p = p + 1
    ! total grain number
    Jacob(p,1:1+nimass)    =  0.0                                                         ! df(total_grains)/dn_{i,e}
    Jacob(p,2+dcol:4+dcol) = -1.0                                                         ! df(total_grains)/dn_g
    p = p + 1
    ! grains (Z =  1)
    Jacob(p,1:nimass) = k_ig(0,1:nimass,j)*n_g(0,j) - k_ig(1,1:nimass,j)*n_g(1,j)         ! df(n_g+)/dn_i
    Jacob(p,1+nimass) =                             - k_eg(1,j)         *n_g(1,j)         ! df(n_g+)/dn_e
    Jacob(p,2+nimass) = 0.0                                                               ! df(n_g+)/dn_g-
    do k = 1,nimass
       Jacob(p,3+dcol) = Jacob(p,3+dcol) + k_ig( 0,k,j)*n_i(k)                            ! df(n_g+)/dn_g0
       Jacob(p,4+dcol) = Jacob(p,4+dcol) - k_ig( 1,k,j)*n_i(k)                            ! df(n_g+)/dn_g+
    enddo
    Jacob(p,4+dcol)   = Jacob(p,4+dcol) - k_eg( 1,j)*n_e                                  ! df(n_g+)/dn_g+
    p = p + 1
 enddo
!
!--Rearrange the Jacobian so that the largest numbers are on the diagonal
 if (reorder_Jacobian) then
    imax = 1
    jmax = 1
    do i = 1,neqn
       feqni(i) = i
    enddo
    do k = 1,neqn-1
       ! find the largest entry & put it in the k'th row
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
       ! Create and use a transformation matrix to put the Jacobmax in the k'th column
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
 !
end subroutine nicil_ionT_calc_Jf
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
 !
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
    if (abs(LU(k,k)) < small) lerr = .true.
    do i = k+1,neqn                           ! Track down the columns to solve for L
       do p = 1,k-1
          LU(i,k) = LU(i,k) - LU(i,p)*LU(p,k)
       enddo
       if (.not.lerr) LU(i,k) = LU(i,k)/LU(k,k)
    enddo
 enddo
 if (lerr) then                              ! exit cleanly since dividing by zero
    dn = 0.0
    LU = 0.0
    return
 endif
 !
 !--Solve L*(Ux) = feqn for Ux == y
 Ux = feqn
 do i = 2,neqn
    do j = 1,i-1
       Ux(i) = Ux(i) - LU(i,j)*Ux(j)
    enddo
 enddo
 !
 !--Solve U*dn = y for dn
 dn = Ux
 do i = neqn,1,-1
    do j = i+1,neqn
       dn(i) = dn(i) - LU(i,j)*dn(j)
    enddo
    dn(i) = dn(i)/LU(i,i)
 enddo
 !
 !--Return the dn array to the correct order
 if (reorder_Jacobian) then
    Ux = dn
    do i = 1,neqn
       dn(feqni(i)) = Ux(i)
    enddo
 endif
 !
end subroutine get_dn_viaLU
!----------------------------------------------------------------------!
!+
! Calculate the electron number density, assuming charge neutrality
! Note: n_R(1:nimass)     are positive ion densities
!       n_R(nimass+2*j-1) are negative grain densities
!       n_R(nimass+2*j  ) are positive grain densities
!+
!----------------------------------------------------------------------!
pure real function nicil_ionR_get_ne(n_R)
 real, intent(in) :: n_R(:)
 integer          :: j,k
 !
 nicil_ionR_get_ne = 0.0
 do k = 1,nimass
    nicil_ionR_get_ne = nicil_ionR_get_ne + n_R(k)
 enddo
 do j = 1,na
    if (abs(n_R(nimass+2*j-1)-n_R(nimass+2*j)) > coef_epsilonR*max(n_R(nimass+2*j-1),n_R(nimass+2*j))) then
       nicil_ionR_get_ne = nicil_ionR_get_ne - n_R(nimass+2*j-1)
       nicil_ionR_get_ne = nicil_ionR_get_ne + n_R(nimass+2*j  )
    endif
 enddo
 !
end function nicil_ionR_get_ne
!----------------------------------------------------------------------!
!+
! Average Z_grain method to approximate number densities if the Jacobian
! method fails.  This method assumes
! n_n ~ n, k_ei = 0, Z_g -> (< Z_g >) < 0
! It solves for < Z_g >, then determines n_i and n_e.  From these
! it then determines the values for n_g(Z=+/-1)
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionR_get_n_via_Zave(n_R,n_e,n0,n_grain,T,zeta,k_ig,k_eg,fBdust,ierr)
 real,    intent(out)   :: n_R(:),n_e,fBdust
 real,    intent(in)    :: n0,T,zeta
 real,    intent(in)    :: n_grain(:),k_ig(-zmax:zmax,nimass,na),k_eg(-zmax:zmax,na)
 integer, intent(inout) :: ierr
 integer                :: j,k,iter
 real                   :: sqrtT, n_grain_tot
 real                   :: k_fac(na),exp_kfacZ(na),k_facZ(na)
 real                   :: n_i(nimass),n_g(-1:1),numer
 real                   :: Z_coef,Z_rat,Z_old,Z_new,fatZ,fatZdZ
 real                   :: kn_ig_coef_Zave(nimass,na),dkn_ig_coef_Zave(nimass,na),kn_ig_Zave(nimass)
 real                   :: kn_ig_Zave1(nimass),dkn_ig_Zave(nimass)
 real                   :: kn_eg_coef_Zave(na),dkn_eg_coef_Zave(na),kn_eg_Zave,kn_eg_Zave1,dkn_eg_Zave
 real                   :: sum_kn,sum_dkn
 logical                :: iterate,negative_n,negative_ng,positive_z
 !
 !--Initialise global values
 !  For k_fac, k_fac_coef(1,j) is used since that value assumes Z==1, thus it can be safely multiplied by Z_ave
 n_g              = 0.0
 n_R              = 0.0
 n_grain_tot      = 0.
 iter             = 0
 sqrtT            = sqrt(T)
 do j = 1,na
    n_grain_tot    = n_grain_tot + n_grain(j)
 enddo
 k_fac            = k_fac_coef(1,1:na)/T
 do k = 1,nimass
    kn_ig_coef_Zave(k,1:na)  = k_ig_coef(k,1:na)*sqrtT*n_grain
    dkn_ig_coef_Zave(k,1:na) = kn_ig_coef_Zave(k,1:na)*k_fac
 enddo
 kn_eg_coef_Zave  = k_eg_coef(  1:na)*sqrtT*n_grain(1:na)
 dkn_eg_coef_Zave = kn_eg_coef_Zave*k_fac
 Z_old            = -small
 Z_coef           = zeta*n0/n_grain_tot
 iterate          = .true.
 negative_n       = .false.
 negative_ng      = .false.
 positive_z       = .false.
 !
 !--Perform the iterations
 do while ( iterate )
    k_facZ      = k_fac*Z_old
    exp_kfacZ   = exp(k_facZ)
    kn_ig_Zave  =  0.0
    kn_eg_Zave  =  0.0
    dkn_ig_Zave =  0.0
    dkn_eg_Zave =  0.0
    do j = 1,na
       kn_ig_Zave  =  kn_ig_Zave +  kn_ig_coef_Zave(:,j)*(1.0-k_facZ(j))
       kn_eg_Zave  =  kn_eg_Zave +  kn_eg_coef_Zave(  j)*exp_kfacZ(j)
       dkn_ig_Zave = dkn_ig_Zave + dkn_ig_coef_Zave(:,j)
       dkn_eg_Zave = dkn_eg_Zave - dkn_eg_coef_Zave(  j)*exp_kfacZ(j)
    enddo
    !
    do k = 1,nimass
       if (kn_ig_Zave(k) > 0.0) then
          kn_ig_Zave1(k) = 1.0/kn_ig_Zave(k)
       else
          kn_ig_Zave1(k) = 0.0
       endif
    enddo
    if (kn_eg_Zave > 0.0) then
       kn_eg_Zave1 = 1.0/kn_eg_Zave
    else
       kn_eg_Zave1 = 0.0
    endif
    sum_kn  =  kn_eg_Zave1
    sum_dkn = dkn_eg_Zave*kn_eg_Zave1**2
    do k = 1,nimass
       sum_kn  = sum_kn  - kn_ig_Zave1(k)
       sum_dkn = sum_dkn - dkn_ig_Zave(k)*kn_ig_Zave1(k)**2
    enddo
    fatZ   = Z_old - Z_coef*sum_kn
    fatZdZ = 1.0   - Z_coef*sum_dkn
    !
    Z_new  = Z_old - fatZ/fatZdZ
    Z_rat  = abs( 1.0 - Z_new/Z_old )
    !
    Z_old  = Z_new
    iter   = iter + 1
    if (iter >= NRctrmax .or. Z_rat < NRtol .or. Z_new > 0.0) iterate = .false.
 enddo
 !
 !--The approximated ion and electron number densities
 n_i = zeta*n0/kn_ig_Zave
 n_e = zeta*n0/kn_eg_Zave
 n_R(1:nimass) =  n_i
 !
 !--Use n_i and n_e to calculate the grain number densities using their correct charge
 do j = 1,na
    n_g   = 0.0
    numer = 0.0
    do k = 1,nimass
       n_g(-1) = n_g(-1) +                    k_ig(-1,k,j) *n_i(k)
       n_g( 1) = n_g( 1) + (2.0*k_ig(0,k,j) + k_ig( 1,k,j))*n_i(k)
       numer   = numer   + k_ig(0,k,j)*n_i(k)
    enddo
    n_g(-1) = n_g(-1) + ( 2.0*k_eg(0,j) + k_eg(-1,j) )*n_e
    n_g( 1) = n_g( 1) +   k_eg( 1,j)                  *n_e
    if (n_g(-1) > 0.0)  n_g(-1) = (1.0+Z_new)*k_eg(0,j)*n_e*n_grain(j)/n_g(-1)
    if (n_g( 1) > 0.0)  n_g( 1) = (1.0-Z_new)*numer        *n_grain(j)/n_g( 1)
    ! recalculate n_g(1) as an average of the current value, and that calculated using n_g(-1) and charge neutrality
    n_g( 1) = 0.5*n_g(1) + 0.5*( n_e + n_g(-1) - n_i(1) - n_i(2) )
    ! recalculate ng(-1) using charge neutrality
    n_g(-1) = n_i(1) + n_i(2) + n_g(1) - n_e
    !
    !--Verify against absurb answers
    do k = 1,nimass
       if (n_R(k) < 0.0) negative_n  = .true.
    enddo
    do k = -zmax,zmax
       if (n_g(k) < 0.0) negative_ng = .true.
    enddo
    !--If grain number densities are negative, but trivial compared to ions, reset to zero
    !  Using coef4 since test show we need more lee-way than allowed by coef
    if (negative_ng) then
       negative_ng = .false.
       do k = -zmax,zmax
          if (k/=0 .and. abs(n_g(k)) > coef4_epsilonR*max(n_i(1),n_i(2))) negative_ng = .true.
       enddo
    endif
    if (negative_ng) then
       negative_n = .true. ! This really is a negative grain number density
    else
       n_g(-1)    = small
       n_g( 1)    = small
    endif
    if (Z_new > 0.0) positive_z = .true.
    !
    !--Update the output array
    n_R(nimass+2*j-1) = n_g(-1)
    n_R(nimass+2*j  ) = n_g( 1)
 enddo
 fBdust = (n_g(-1)+n_g(1))/n_grain_tot
 !
 !--warnings
 if (warn_verbose)     ierr = ierr + ierr_Zave   ! Trigger non-fatal warning since this subroutine is used
 if (negative_n)       ierr = ierr + ierr_nRle0  ! Trigger fatal warning if n < 0 (unphysical)
 if (positive_z)       ierr = ierr + ierr_Zge0   ! Trigger fatal warning if Z > 0 (violates assumption)
 if (iter >= NRctrmax) ierr = ierr + ierr_nRconv ! Trigger fatal warning if too many iterations
 !
end subroutine nicil_ionR_get_n_via_Zave
!======================================================================!
! THERMAL IONISATION-RELATED SUBROUTINES                               !
!======================================================================!
!----------------------------------------------------------------------!
!+
!  This will solve the Saha equation to determine the electron number
!  density.  The species of interest and thier properties are given
!  at the beginnig of this module.
!  This assumes for each element, there are neutral, singly-ionised and
!  doubly-ionised species; hydrogen can only be singly ionised.
!  We will include the electron number density calculated by the cosmic
!  ray ionisaiton as well.
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_ne(n_electronT,n_total,T,afrac,ierr)
 integer               :: iter
 integer,intent(inout) :: ierr
 real,   intent(inout) :: n_electronT,afrac(2)
 real,   intent(in)    :: n_total,T
 real                  :: mass_ionT_mp(nlevels),n_ionT(nlevels)
 real                  :: fatn,fatndn,nerat,neold,nenew
 real                  :: neS,dneSdne,ne0
 real                  :: Kjk(nlevels,nelements),nj(nelements),njk(nlevels,nelements)
 logical               :: iterate,try_new_ne0
 !
 !--Initialise values
 iter        = 0
 nerat       = NRtol*2.0
 neold       = n_electronT
 ne0         = neold
 iterate     = .true.
 try_new_ne0 = .true.
 call nicil_ionT_get_nj_Kjk(nj,Kjk,T,n_total,afrac)
 !
 !--Calcualte n_e
 do while ( iterate )
    call nicil_ionT_get_nion(n_ionT,neS,dneSdne,mass_ionT_mp,neold,njk,nj,Kjk,iter)
    fatn   = neold - neS
    fatndn = 1.0   - dneSdne
    if ( iterate ) then
       nenew = neold - fatn/fatndn
       if (nenew > 0.0) then
          nerat = abs( 1.0 - neold/nenew )
       else
          nerat = 0.0
       endif
       ! Actions if converged
       if (nerat < NRtol) then
          iterate     = .false.
          n_electronT = nenew
       endif
       ! Actions if errors occurred
       if (iterate .and. (iter >= NRctrmax .or. nenew < 0.0)) then
          if (try_new_ne0) then
             ! Try again with default guess rather than value from previous iteration
             try_new_ne0 = .false.
             nenew       = epsilon(nenew)
             iter        = 0
             nerat       = NRtol*2.0
          else
             ! New guess failed; trigger warnings
             ierr = ierr + ierr_neTconv  ! n_electronT did not converge
             ierr = ierr + ierr_neTle0   ! n_electronT < 0
          endif
       endif
    endif
    iter  = iter + 1
    neold = nenew
 enddo
 !
end subroutine nicil_ionT_get_ne
!----------------------------------------------------------------------!
!+
!  This is a stand-alone routine to calculate n_ionT, mass_ionT & njk
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_n(n_ionT,mass_ionT,njk,n_electronT,T,n_total,afrac)
 real,    intent(out)   :: mass_ionT(:),n_ionT(:),njk(:,:)
 real,    intent(in)    :: n_electronT,T,n_total
 real,    intent(inout) :: afrac(2)
 integer                :: iter
 real                   :: neS,dneSdne
 real                   :: Kjk(nlevels,nelements),nj(nelements),mass_ionT_mp(nlevels)
 !
 iter = 0
 call nicil_ionT_get_nj_Kjk(nj,Kjk,T,n_total,afrac)
 call nicil_ionT_get_nion(n_ionT,neS,dneSdne,mass_ionT_mp,n_electronT,njk,nj,Kjk,iter)
 mass_ionT = mass_ionT_mp*mass_proton
 !
end subroutine nicil_ionT_get_n
!----------------------------------------------------------------------!
!+
!  This is a stand-alone routine to calculate nj & Kjk
!  (done here due to restriction on Kjk)
!  Note: fractions are assuming all hydrogen is atomic; if some is
!  molecular, then we need to re-normalise the number densities
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_nj_Kjk(nj,Kjk,T,n_total,afrac)
 real,    intent(out)   :: nj(:),Kjk(:,:)
 real,    intent(in)    :: T,n_total
 real,    intent(inout) :: afrac(2)
 integer                :: i,j,k
 real                   :: afrac_total

 afrac_total = afrac(1)+afrac(2)
 do i = 3,nelements
    afrac_total = afrac_total + abundancej(i)
 enddo

 nj(iH2)         = afrac(1)*n_total
 nj(iH)          = afrac(2)*n_total
 nj(3:nelements) = abundancej(3:nelements)*n_total
 nj              = nj/afrac_total
 afrac           = afrac/afrac_total  ! for correct densities in the output files
 Kjk             = 0.

 do k = 1,nelements
    do j = 1,nlevels
       if (T > Texp_thresh(j,k) .and. chij(j,k) > small) then
          Kjk(j,k) = Saha_coef(j,k)*sqrt(T)**3*exp(-chij(j,k)/T)
       endif
    enddo
 enddo

end subroutine nicil_ionT_get_nj_Kjk
!----------------------------------------------------------------------!
!+
!  This calculates the total ion number density, and average ion mass
!+
!----------------------------------------------------------------------!
pure subroutine nicil_ionT_get_nion(ni,neS,dneSdne,mi_mp,ne,njk,nj,Kjk,iter)
 integer, intent(in)  :: iter
 real,    intent(out) :: neS,dneSdne,ni(:),mi_mp(:),njk(:,:)
 real,    intent(in)  :: ne,nj(:),Kjk(:,:)
 integer              :: j,k
 real                 :: term,term1,ne1,KKonne,m_total(nlevels)
 !
 njk     = 0.0
 ni      = 0.0
 neS     = 0.0
 dneSdne = 0.0
 mi_mp   = 0.0
 m_total = 0.0
 if ( ne > 0.0 .and. (iter==0 .or. (iter>0 .and. abs(ne-epsilon(ne)) > small ) ) ) then
    ne1 = 1.0/ne
    do k = 1,nelements
       KKonne = Kjk(1,k)*Kjk(2,k)*ne1
       term = ne + Kjk(1,k) + KKonne
       if (term > 0.0) then
          term1 = 1.0/term
       else
          term1 = 0.0
       endif
       if (Kjk(1,k) > 0.0) njk(1,k) = nj(k)*Kjk(1,k)*term1
       if (Kjk(2,k) > 0.0) njk(2,k) = njk(1,k)*Kjk(2,k)*ne1
       ni      = ni      + njk(:,k)
       neS     = neS     + njk(1,k) + 2.0*njk(2,k)
       dneSdne = dneSdne - njk(1,k)*term1*(1.0-KKonne*ne1)   &
                         - 2.0*njk(2,k)*term1*ne1*(2.0*ne+Kjk(1,k))
       m_total(1) = m_total(1) + njk(1,k)*(mj_mp(k)-    massj_mp(ine))
       m_total(2) = m_total(2) + njk(2,k)*(mj_mp(k)-2.0*massj_mp(ine))
    enddo
    do j = 1,nlevels
       if (ni(j) > 0.0) mi_mp(j) = m_total(j)/ni(j)  ! this is consistent with how we calculate m_neutral
    enddo
 endif
 !
end subroutine nicil_ionT_get_nion
!----------------------------------------------------------------------!
!+
!  This calculates the number density of H and H2, and the neutral
!  number density based upon these abundances
!+
!----------------------------------------------------------------------!
pure subroutine nicil_get_HH2_ratio(abund,n_cold,T,mass_neutral_mp,afrac,mfrac)
 real,             intent(in)  :: abund,n_cold,T
 real,             intent(out) :: mass_neutral_mp
 real,             intent(out) :: afrac(2)
 real,    optional,intent(out) :: mfrac(2)
 real                          :: KHH2,nH2,frac,mtotal1
 real                          :: mfrac0(2)
 !
 nH2  = 0.5*abund*n_cold
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

 afrac(1)  = 0.5*(1.0-frac)*abund
 afrac(2)  =          frac *abund

 mtotal1   = 1.0/((1.0-frac)*mj_mp(1) + frac*mj_mp(2))
 mfrac0(1) = mj_mp(1)*(1.0-frac)*mtotal1 * mass_frac(iH)
 mfrac0(2) = mj_mp(2)*     frac *mtotal1 * mass_frac(iH)
 if (present(mfrac)) mfrac = mfrac0
 ! calculate the neutral mass using the correct ratios of H & H2
 mass_neutral_mp = 1.0/( mfrac_by_m_notH_mp + mfrac0(1)*mj_mp1(1) + mfrac0(2)*mj_mp1(2) )
 !
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
 !
 eta_ohm  = 0.0
 eta_hall = 0.0
 eta_ambi = 0.0
 if (sigmas(6) > 0.0) then ! a verification that we have not entered the ideal regime
    if (use_ohm)  eta_ohm  = csqbyfourpi * sigmas(6)
    if (use_hall) eta_hall = csqbyfourpi * sigmas(5) * sigmas(7)
    if (use_ambi) eta_ambi = csqbyfourpi * sigmas(8) * sigmas(6) * sigmas(7)
 endif
 !
end subroutine nicil_nimhd_get_eta
!-----------------------------------------------------------------------
pure subroutine nicil_nimhd_get_eta_cnst(eta_ohm,eta_hall,eta_ambi,Bfield,rho)
 real, intent(out) :: eta_ohm,eta_hall,eta_ambi
 real, intent(in)  :: Bfield,rho
 !
 ! set the coefficient to the constant coefficient (==0 if the term is off)
 eta_ohm  = eta_ohm_cnst
 eta_hall = eta_hall_cnst
 eta_ambi = eta_ambi_cnst
 if (eta_const_type==icnst) return
 ! multiply by the variable, if requested
 eta_hall = eta_hall*Bfield
 eta_ambi = eta_ambi*Bfield**2/rho**alpha_AD_p1
 !
end subroutine nicil_nimhd_get_eta_cnst
!-----------------------------------------------------------------------
!+
!  Calculates JxB and (JxB)xB
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_jcbcb(jcbcb,jcb,jcurrent,Bx,By,Bz,B1)
 real, intent(out) :: jcb(3), jcbcb(3)
 real, intent(in)  :: jcurrent(3)
 real, intent(in)  :: Bx, By, Bz, B1
 !
 jcb(1)   = ( jcurrent(2)*Bz - jcurrent(3)*By )*B1
 jcb(2)   = ( jcurrent(3)*Bx - jcurrent(1)*Bz )*B1
 jcb(3)   = ( jcurrent(1)*By - jcurrent(2)*Bx )*B1
 !
 jcbcb(1) = ( jcb(2)*Bz      - jcb(3)*By      )*B1
 jcbcb(2) = ( jcb(3)*Bx      - jcb(1)*Bz      )*B1
 jcbcb(3) = ( jcb(1)*By      - jcb(2)*Bx      )*B1
 !
end subroutine nimhd_get_jcbcb
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD contributions to the magnetic field in SPH
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dBdt(dBnonideal,eta_ohm,eta_hall,eta_ambi &
                              ,jcurrent,jcb,jcbcb,dxr1,dyr1,dzr1)
 real, intent(out) :: dBnonideal(3)
 real, intent(in)  :: jcurrent(3),jcb(3),jcbcb(3)
 real, intent(in)  :: eta_ohm,eta_hall,eta_ambi,dxr1,dyr1,dzr1
 real              :: dBohm(3),dBhall(3),dBambi(3)
 !
 dBohm  = 0.0
 dBhall = 0.0
 dBambi = 0.0
 !
 if (use_ohm ) call nimhd_get_DcrossR(dBohm ,jcurrent,dxr1,dyr1,dzr1,eta_ohm )
 if (use_hall) call nimhd_get_DcrossR(dBhall,jcb     ,dxr1,dyr1,dzr1,eta_hall)
 if (use_ambi) call nimhd_get_DcrossR(dBambi,jcbcb   ,dxr1,dyr1,dzr1,eta_ambi)
 dBnonideal = dBambi - dBhall - dBohm
 !
end subroutine nimhd_get_dBdt
!-----------------------------------------------------------------------
!+
!  performs simple cross product
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_DcrossR(DcrossR,D_in,dx,dy,dz,eta)
 real, intent (in)  :: D_in(3)
 real, intent (out) :: DcrossR(3)
 real, intent (in)  :: dx, dy, dz, eta
 !
 DcrossR(1) = (D_in(2)*dz - D_in(3)*dy)*eta
 DcrossR(2) = (D_in(3)*dx - D_in(1)*dz)*eta
 DcrossR(3) = (D_in(1)*dy - D_in(2)*dx)*eta
 !
end subroutine nimhd_get_DcrossR
!-----------------------------------------------------------------------
!+
!  Calculates the non-ideal MHD contributions to energy
!  Note: dudthall==0
!  Note: the sign of (36) in Wurster (2016) is incorrect
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dudt(dudtnonideal,eta_ohm,eta_ambi,rho,J,B)
 real, intent(out) :: dudtnonideal
 real, intent(in)  :: J(3),B(3)
 real, intent(in)  :: eta_ohm,eta_ambi,rho
 real              :: B2i,J2i,BJi,BJBJihat
 !
 B2i          = dot_product(B,B)
 J2i          = dot_product(J,J)
 if (B2i > 0.0) then
    BJi        = dot_product(B,J)
    BJBJihat   = BJi*BJi/B2i
 else
    BJBJihat   = 0.0
 endif
 dudtnonideal = ( eta_ohm*J2i + eta_ambi*(J2i - BJBJihat) )/rho
 !
end subroutine nimhd_get_dudt
!-----------------------------------------------------------------------
!+
!  Calculates the timesteps for non-ideal MHD
!+
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dt(dtohm,dthall,dtambi,h,eta_ohm,eta_hall,eta_ambi)
 real, intent(in)  :: h,eta_ohm,eta_hall,eta_ambi
 real, intent(out) :: dtohm,dthall,dtambi
 real              :: h2
 !
 dtohm  = huge(dtohm)
 dthall = huge(dthall)
 dtambi = huge(dtambi)
 h2     = h*h
 !
 if (use_ohm  .and.     eta_ohm   > tiny(eta_ohm ) ) dtohm  =      C_nimhd*h2/eta_ohm
 if (use_hall .and. abs(eta_hall) > tiny(eta_hall) ) dthall = abs( C_nimhd*h2/eta_hall )
 if (use_ambi .and.     eta_ambi  > tiny(eta_ambi) ) dtambi =      C_nimhd*h2/eta_ambi
 !
end subroutine nimhd_get_dt
!-----------------------------------------------------------------------
!+
!  Calculates the Hall drift velocity
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_halldrift(eta_hall,Bx,By,Bz,jcurrent,vdrift)
 real,    intent(in)    :: eta_hall,Bx,By,Bz,jcurrent(3)
 real,    intent(out)   :: vdrift(3)
 real                   :: B1
 !
 B1       = 1.0/sqrt(Bx*Bx + By*By + Bz*Bz)
 ! The Hall drift velocity
 vdrift = -eta_hall*jcurrent*B1
 !
end subroutine nicil_get_halldrift
!-----------------------------------------------------------------------
!+
!  Calculates the ion drift velocity caused by ambipolar diffusion
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_vdrift(eta_ambi,Bx,By,Bz,jcurrent,vdrift)
 real,    intent(in)    :: eta_ambi,Bx,By,Bz,jcurrent(3)
 real,    intent(out)   :: vdrift(3)
 real                   :: B21
 !
 B21       = 1.0/(Bx*Bx + By*By + Bz*Bz)
 ! The ion drift velocity
 vdrift(1) = eta_ambi*( jcurrent(2)*Bz - jcurrent(3)*By )*B21
 vdrift(2) = eta_ambi*( jcurrent(3)*Bx - jcurrent(1)*Bz )*B21
 vdrift(3) = eta_ambi*( jcurrent(1)*By - jcurrent(2)*Bx )*B21
 !
end subroutine nicil_get_vdrift
!-----------------------------------------------------------------------
!+
!  Calculates the ion velocity caused by ambipolar diffusion
!  v_ion = v_gas + v_drift
!+
!-----------------------------------------------------------------------
pure subroutine nicil_get_vion(eta_ambi,vx,vy,vz,Bx,By,Bz,jcurrent,vion,ierr,vdrift_out)
 integer,           intent(inout) :: ierr
 real,              intent(in)    :: eta_ambi,vx,vy,vz,Bx,By,Bz,jcurrent(3)
 real,              intent(out)   :: vion(3)
 real,    optional, intent(out)   :: vdrift_out(3)
 real                             :: vion2,v2
 real                             :: vdrift(3)
 !
 ! The drift velocity
 call nicil_get_vdrift(eta_ambi,Bx,By,Bz,jcurrent,vdrift)
 !
 ! The ion velocity
 vion(1) = vx + vdrift(1)
 vion(2) = vy + vdrift(2)
 vion(3) = vz + vdrift(3)
 !
 if (warn_verbose) then
    ! ensure that the drift velocity is small
    vion2 = vion(1)*vion(1) + vion(2)*vion(2) + vion(3)*vion(3)
    v2    = vx*vx + vy*vy + vz*vz
    if (warn_ratio_m12*vion2 > v2 .or. v2 > warn_ratio_p12*vion2) ierr = ierr + ierr_drift
    if (vrms2_min         > vion2 .or. vion2 > vrms2_max        ) ierr = ierr + ierr_vrms
 endif
 !
 if (present(vdrift_out)) vdrift_out = vdrift
 !
end subroutine nicil_get_vion

!----------------------------------------------------------------------!
end module nicil
