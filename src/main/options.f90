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
! :Runtime parameters:
!   - calc_erot : *include E_rot in the ev_file*
!   - curlv     : *output curl v in dump files*
!   - track_lum : *write du/dt to dump files (for a "lightcurve")*
!
! :Dependencies: damping, dim, eos, infile_utils, io_control, kernel,
!   mcfost_utils, part, radiation_utils, shock_capturing, timestep, units,
!   viscosity
!
 use eos,             only:ieos,icooling,iopacity_type,use_var_comp ! so this is available via options module
 use damping,         only:idamp ! so this is available via options module
 use dim,             only:curlv,track_lum ! make available from options module
 use part,            only:tolh  ! make available from options module
 use mcfost_utils,    only:use_mcfost,use_mcfost_stellar_parameters
 use radiation_utils, only:implicit_radiation,limit_radiation_flux,implicit_radiation_store_drad
 use shock_capturing, only:alpha,alphamax,alphau,alphaB,beta,disc_viscosity,ireconav
 use io_control,      only:nfulldump
 implicit none
!
! these are parameters which may be changed by the user
! and read from the input file
!
 integer, public :: iexternalforce

! additional .ev data
 logical, public :: calc_erot

! dust method
 logical, public :: use_dustfrac, use_hybrid, use_porosity

 ! pressure on sinks
 logical, public :: need_pressure_on_sinks

! library use
 logical, public :: write_files

 public :: set_default_options,write_options_output,read_options_output

 ! options from lower-level modules that can also be imported via options module
 public :: ieos,icooling,idamp,tolh
 public :: iopacity_type
 public :: use_var_comp  ! use variable composition
 public :: curlv
 public :: use_mcfost,use_mcfost_stellar_parameters
 public :: implicit_radiation,limit_radiation_flux,implicit_radiation_store_drad
 public :: alpha,alphamax,alphau,alphaB,ireconav,beta
 public :: nfulldump

 private

contains

subroutine set_default_options
 use timestep,        only:set_defaults_timestep
 use part,            only:hfact,Bextx,Bexty,Bextz,tolh
 use viscosity,       only:set_defaults_viscosity
 use dim,             only:gr,do_radiation,isothermal
 use kernel,          only:hfact_default
 use eos,             only:set_defaults_eos
 use units,           only:set_units
 use mcfost_utils,    only:set_defaults_mcfost
 use radiation_utils, only:set_defaults_radiation
 use shock_capturing, only:set_defaults_shock_capturing
 use io_control,      only:set_defaults_iocontrol

 ! Default timestepping options
 call set_defaults_timestep

 ! Default io control options
 call set_defaults_iocontrol

 ! Reset units
 call set_units()

 ! Miscellaneous parameters
 hfact     = hfact_default   ! smoothing length in units of average particle spacing
 tolh      = 1.e-4           ! tolerance on h iterations
 Bextx     = 0.              ! external magnetic field
 Bexty     = 0.
 Bextz     = 0.
 iexternalforce = 0          ! external forces
 if (gr) iexternalforce = 1
 calc_erot = .false.         ! To allow rotational energies to be printed to .ev

 ! equation of state
 call set_defaults_eos

 ! shock capturing
 call set_defaults_shock_capturing

 ! physical viscosity
 call set_defaults_viscosity

 ! mcfost
 call set_defaults_mcfost

 ! damping
 idamp = 0

 ! dust method
 use_dustfrac = .false.

 ! radiation
 call set_defaults_radiation

 ! pressure on sinks
 need_pressure_on_sinks = .false.

 ! enable/disable writing output files
 write_files = .true.

end subroutine set_default_options

!-----------------------------------------------------------------------
!+
!  Write options controlling optional output
!+
!-----------------------------------------------------------------------
subroutine write_options_output(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(curlv,'curlv','output curl v in dump files',iunit)
 call write_inopt(track_lum,'track_lum','write du/dt to dump files (for a "lightcurve")',iunit)
 if (calc_erot) call write_inopt(calc_erot,'calc_erot','include E_rot in the ev_file',iunit)

end subroutine write_options_output

!-----------------------------------------------------------------------
!+
!  Read options controlling optional output
!+
!-----------------------------------------------------------------------
subroutine read_options_output(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 ! Non-compulsory: keep current values as defaults
 call read_inopt(curlv,    'curlv',    db,errcount=nerr,default=curlv)
 call read_inopt(track_lum,'track_lum',db,errcount=nerr,default=track_lum)
 call read_inopt(calc_erot,'calc_erot',db,errcount=nerr,default=calc_erot)

end subroutine read_options_output

end module options
