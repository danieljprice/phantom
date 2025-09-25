!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module options
!
! Sets default values of input parameters
!  these are overwritten by reading from the input file or
!  by setting them in the setup routine
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: damping, dim, eos, kernel, part, timestep, units,
!   viscosity
!
 use eos,             only:ieos,iopacity_type,use_var_comp ! so this is available via options module
 use damping,         only:idamp ! so this is available via options module
 use dim,             only:curlv ! make available from options module
 use mcfost_utils,    only:use_mcfost,use_mcfost_stellar_parameters
 use radiation_utils, only:implicit_radiation,limit_radiation_flux,implicit_radiation_store_drad
 use shock_capturing, only:alpha,alphamax,alphau,alphaB,beta,disc_viscosity,ireconav
 implicit none
!
! these are parameters which may be changed by the user
! and read from the input file
!
 integer, public :: nfulldump,nmaxdumps,iexternalforce
 real, public :: tolh,rkill
 real(kind=4), public :: twallmax

 integer, public :: ishock_heating,ipdv_heating,icooling,iresistive_heating

 ! div B cleaning
 real, public :: psidecayfac, overcleanfac

! additional .ev data
 logical, public :: calc_erot
! final maximum density
 real,    public :: rhofinal_cgs,rhofinal1

! dust method
 logical, public :: use_dustfrac, use_hybrid, use_porosity

 ! pressure on sinks
 logical, public :: need_pressure_on_sinks

! library use
 logical, public :: write_files

 public :: set_default_options

 ! options from lower-level modules that can also be imported via options module
 public :: ieos,idamp
 public :: iopacity_type
 public :: use_var_comp  ! use variable composition
 public :: curlv
 public :: use_mcfost,use_mcfost_stellar_parameters
 public :: implicit_radiation,limit_radiation_flux,implicit_radiation_store_drad
 public :: alpha,alphamax,alphau,alphaB,ireconav,beta

 private

contains

subroutine set_default_options
 use timestep,        only:set_defaults_timestep
 use part,            only:hfact,Bextx,Bexty,Bextz,ien_type,ien_entropy
 use viscosity,       only:set_defaults_viscosity
 use dim,             only:gr,do_radiation,isothermal
 use kernel,          only:hfact_default
 use eos,             only:polyk2
 use units,           only:set_units
 use mcfost_utils,    only:set_defaults_mcfost
 use radiation_utils, only:set_defaults_radiation
 use shock_capturing, only:set_defaults_shock_capturing
 use dynamic_dtmax,   only:set_defaults_dynamic_dtmax

 ! Default timestepping options
 call set_defaults_timestep

 ! Default dynamic dtmax options
 call set_defaults_dynamic_dtmax

 ! Reset units
 call set_units()

 ! Miscellaneous parameters
 nmaxdumps = -1
 twallmax  = 0.0             ! maximum wall time for run, in seconds
 nfulldump = 10              ! frequency of writing full dumps
 hfact     = hfact_default   ! smoothing length in units of average particle spacing
 Bextx     = 0.              ! external magnetic field
 Bexty     = 0.
 Bextz     = 0.
 tolh      = 1.e-4           ! tolerance on h iterations
 iexternalforce = 0          ! external forces
 if (gr) iexternalforce = 1
 calc_erot = .false.         ! To allow rotational energies to be printed to .ev
 rhofinal_cgs = 0.           ! Final maximum density (0 == ignored)

 ! equation of state
 if (.not.isothermal) then
    ieos = 2
 else
    ieos = 1
 endif
 ishock_heating     = 1
 ipdv_heating       = 1
 iresistive_heating = 1
 icooling           = 0
 ien_type           = 0
 if (gr) ien_type   = ien_entropy
 polyk2             = 0. ! only used for ieos=8
 use_var_comp = .false.  ! variable composition

 ! shock capturing
 call set_defaults_shock_capturing

 ! div B cleaning (MHD only)
 psidecayfac       = 1.0     ! ratio of parabolic to hyperbolic cleaning
 overcleanfac      = 1.0     ! factor by which to increase cleaning speed for div B cleaning

 ! physical viscosity
 call set_defaults_viscosity

 ! mcfost
 call set_defaults_mcfost

 ! damping
 idamp = 0

 ! radius outside which we kill particles
 rkill             = -1.

 ! dust method
 use_dustfrac = .false.

 ! radiation
 call set_defaults_radiation

 ! pressure on sinks
 need_pressure_on_sinks = .false.

 ! enable/disable writing output files
 write_files = .true.

end subroutine set_default_options

end module options
