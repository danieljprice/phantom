!----------------------------------------------------------------------!
!                               N I C I L                              !
!           Non-Ideal mhd Coefficient and Ionisation Library           !
!                                                                      !
! This is a stand-alone library that will calculate ionisation values  !
! and the coefficients for the non-ideal MHD terms: Ohmic resistivity, !
! Hall Effect and Ambipolar diffusion.                                 !
!                                                                      !
!                 Copyright (c) 2015-2021 James Wurster                !
!        See LICENCE file for usage and distribution conditions        !
!----------------------------------------------------------------------!
! THIS IS A DUMMY-file that will exit since we nicil needs to be
! compiled with real*8
!----------------------------------------------------------------------!
!+
!  MODULE: nicil
!
!  DESCRIPTION:
!  Contains the public values of nicil, but does nothing!
!  Copyright (c) 2015-2021 James Wurster
!  References: Wurster (2016) PASA, 33:e041.
!              Wurster (2021) MNRAS, 501:5873-5891.
!  See LICENCE file for usage and distribution conditions
!
!
!
!  AUTHOR: James Wurster
!
!----------------------------------------------------------------------!
module nicil
 implicit none
 !--INPUT PAREMETERS
 !--Number of grain sizes
 integer, public, parameter :: na                = 1                ! Hardcoded for storage efficiency ( na >=0 )
 !--Turn on/off individual non-ideal MHD coefficients
 logical, public            :: use_ohm           = .true.           ! Calculate the coefficient for Ohmic resistivity
 logical, public            :: use_hall          = .true.           ! Calculate the coefficient for the Hall effect
 logical, public            :: use_ambi          = .true.           ! Calculate the coefficient for ambipolar diffusion
 !--Use a constant (false) or variable (true) cosmic ray ionisation rate
 logical, public            :: zeta_of_rho       = .false.
 !--Import the dust mass density from the parent code for each point
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

 !--Additional parameters
 integer, public            :: NRctrmax          = 200              ! maximum number of Newton–Raphson iterations
 real,    public            :: NRtol             = 1.0d-8           ! default tolerance on Newton–Raphson iterations
 real,    public            :: NRtol4            = 1.0d-6           ! default tolerance on Newton–Raphson iterations if single precsision
 real,    public            :: Cdt_diff          = 0.12             ! Coefficient to control the AD & OR timesteps
 real,    public            :: Cdt_hall          = 0.0795774        ! Coefficient to control the HE timestep (==1/4pi)
 !--END OF INPUT PARAMETERS
 integer, public, parameter :: n_nden            = 1
 integer, public, parameter :: n_warn            = 1
 integer, public, parameter :: n_data_out        = 1
 real,    public            :: meanmolmass,unit_eta

 !--Subroutines
 public  :: nicil_initialise,nicil_update_nimhd,nicil_translate_error
 public  :: nicil_get_dudt_nimhd,nicil_get_dt_nimhd
 public  :: nicil_get_vion,nicil_get_ambidrift,nicil_get_halldrift
 public  :: nimhd_get_jcbcb,nimhd_get_dBdt

 private

contains

!----------------------------------------------------------------------!
subroutine nicil_initialise(utime,udist,umass,unit_Bfield,ierr,iprint_in,iprintw_in,nden_nimhd0,a_grain_cgs_in)
 real,              intent(in)  :: utime,umass,udist,unit_Bfield
 integer,           intent(out) :: ierr
 real,    optional, intent(in)  :: a_grain_cgs_in(na)
 integer, optional, intent(in)  :: iprint_in,iprintw_in
 real,    optional, intent(out) :: nden_nimhd0(n_nden)

 return

end subroutine nicil_initialise
!----------------------------------------------------------------------!
subroutine nicil_translate_error(ierrlist,fatal_only,rho,B,T)
 integer,        intent(in) :: ierrlist(:)
 logical,        intent(in) :: fatal_only
 real, optional, intent(in) :: rho,B,T

 return
end subroutine nicil_translate_error
!----------------------------------------------------------------------!
pure subroutine nicil_update_nimhd(icall,eta_ohm,eta_hall,eta_ambi,Bfield,rho,T,nden_save, &
                                   ierrlist,data_out,fdg_in)
 integer,          intent(in)    :: icall
 integer,          intent(inout) :: ierrlist(n_warn)
 real,             intent(out)   :: eta_ohm,eta_hall,eta_ambi
 real,             intent(in)    :: Bfield,rho,T
 real,             intent(inout) :: nden_save(:)
 real,   optional, intent(in)    :: fdg_in(:)
 real,   optional, intent(out)   :: data_out(n_data_out)

 eta_ohm  = 0.
 eta_hall = 0.
 eta_ambi = 0.
 if (present(data_out)) data_out = 0.

 return
end subroutine nicil_update_nimhd
!-----------------------------------------------------------------------
pure subroutine nicil_get_dudt_nimhd(dudtnonideal,eta_ohm,eta_ambi,rho,J,B)
 real, intent(out) :: dudtnonideal
 real, intent(in)  :: J(3),B(3)
 real, intent(in)  :: eta_ohm,eta_ambi,rho

 dudtnonideal = 0.

 return
end subroutine nicil_get_dudt_nimhd
!-----------------------------------------------------------------------
pure subroutine nicil_get_dt_nimhd(dtohm,dthall,dtambi,h,eta_ohm,eta_hall,eta_ambi)
 real, intent(in)  :: h,eta_ohm,eta_hall,eta_ambi
 real, intent(out) :: dtohm,dthall,dtambi

 dtohm  = huge(dtohm)
 dthall = huge(dthall)
 dtambi = huge(dtambi)

 return
end subroutine nicil_get_dt_nimhd
!-----------------------------------------------------------------------
pure subroutine nicil_get_halldrift(eta_hall,Bx,By,Bz,jcurrent,vdrift)
 real,    intent(in)    :: eta_hall,Bx,By,Bz,jcurrent(3)
 real,    intent(out)   :: vdrift(3)

 vdrift = 0.

 return
end subroutine nicil_get_halldrift
!-----------------------------------------------------------------------
pure subroutine nicil_get_ambidrift(eta_ambi,Bx,By,Bz,jcurrent,vdrift)
 real,    intent(in)    :: eta_ambi,Bx,By,Bz,jcurrent(3)
 real,    intent(out)   :: vdrift(3)

 vdrift = 0.

 return
end subroutine nicil_get_ambidrift
!-----------------------------------------------------------------------
pure subroutine nicil_get_vion(eta_hall,eta_ambi,vx,vy,vz,Bx,By,Bz,jcurrent,vion,ierrlist,vdrift_out)
 integer,           intent(inout) :: ierrlist(:)
 real,              intent(in)    :: eta_hall,eta_ambi,vx,vy,vz,Bx,By,Bz,jcurrent(3)
 real,              intent(out)   :: vion(3)
 real,    optional, intent(out)   :: vdrift_out(3)

 vion = 0.
 if (present(vdrift_out)) vdrift_out = 0.

 return
end subroutine nicil_get_vion
!-----------------------------------------------------------------------
pure subroutine nimhd_get_jcbcb(jcbcb,jcb,jcurrent,Bx,By,Bz,B1)
 real, intent(out) :: jcb(3), jcbcb(3)
 real, intent(in)  :: jcurrent(3)
 real, intent(in)  :: Bx, By, Bz, B1

 jcb   = 0.
 jcbcb = 0.

 return
end subroutine nimhd_get_jcbcb
!-----------------------------------------------------------------------
pure subroutine nimhd_get_dBdt(dBnonideal,eta_ohm,eta_hall,eta_ambi,jcurrent,jcb,jcbcb,dxr1,dyr1,dzr1)
 real, intent(out) :: dBnonideal(3)
 real, intent(in)  :: jcurrent(3),jcb(3),jcbcb(3)
 real, intent(in)  :: eta_ohm,eta_hall,eta_ambi,dxr1,dyr1,dzr1

 dBnonideal = 0.

 return
end subroutine nimhd_get_dBdt
!-----------------------------------------------------------------------
pure subroutine nimhd_get_DcrossR(DcrossR,D_in,dx,dy,dz,eta)
 real, intent (in)  :: D_in(3)
 real, intent (out) :: DcrossR(3)
 real, intent (in)  :: dx, dy, dz, eta

 DcrossR = 0.

 return
end subroutine nimhd_get_DcrossR
!----------------------------------------------------------------------!
end module nicil
