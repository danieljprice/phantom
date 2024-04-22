!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
 use eos,     only:ieos,iopacity_type,use_var_comp ! so this is available via options module
 use damping, only:idamp ! so this is available via options module
 implicit none
!
! these are parameters which may be changed by the user
! and read from the input file
!
 real, public :: avdecayconst
 integer, public :: nfulldump,nmaxdumps,iexternalforce
 real, public :: tolh,damp,rkill
 integer, parameter :: sp = 4 ! single precision
 real(kind=sp), public :: twallmax

! artificial viscosity, thermal conductivity, resistivity

 real, public :: alpha,alphau,beta
 real, public :: alphamax
 real, public :: alphaB, psidecayfac, overcleanfac, hdivbbmax_max
 integer, public :: ishock_heating,ipdv_heating,icooling,iresistive_heating
 integer, public :: ireconav

! additional .ev data
 logical, public :: calc_erot
! final maximum density
 real,    public :: rhofinal_cgs,rhofinal1

! dust method
 logical, public :: use_dustfrac, use_hybrid, use_porosity

! mcfost
 logical, public :: use_mcfost, use_Voronoi_limits_file, use_mcfost_stellar_parameters, mcfost_computes_Lacc
 logical, public :: mcfost_uses_PdV, mcfost_dust_subl
 integer, public :: ISM
 real(kind=sp), public :: mcfost_keep_part
 character(len=80), public :: Voronoi_limits_file

 ! radiation
 logical, public :: exchange_radiation_energy, limit_radiation_flux, implicit_radiation
 logical, public :: implicit_radiation_store_drad

 public :: set_default_options
 public :: ieos,idamp
 public :: iopacity_type
 public :: use_var_comp  ! use variable composition

 private

contains

subroutine set_default_options
 use timestep,  only:set_defaults_timestep
 use part,      only:hfact,Bextx,Bexty,Bextz,mhd,maxalpha,ien_type,ien_entropy
 use viscosity, only:set_defaults_viscosity
 use dim,       only:maxp,maxvxyzu,nalpha,gr,do_radiation
 use kernel,    only:hfact_default
 use eos,       only:polyk2
 use units,     only:set_units

 ! Default timsteps
 call set_defaults_timestep

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
 if (maxvxyzu==4) then
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

 ! artificial viscosity
 if (maxalpha>0 .and. maxalpha==maxp) then
    if (nalpha >= 2) then
       alpha = 0.0 ! Cullen-Dehnen switch
    else
       alpha = 0.1 ! Morris-Monaghan switch
    endif
 else
    alpha = 1.
 endif
 alphamax = 1.0
 call set_defaults_viscosity
 idamp = 0

 ! artificial thermal conductivity
 alphau = 1.
 if (gr) alphau = 0.1
 ireconav = -1

 ! artificial resistivity (MHD only)
 alphaB            = 1.0
 psidecayfac       = 1.0     ! psi decay factor (MHD only)
 overcleanfac      = 1.0     ! factor to increase signal velocity for (only) time steps and psi cleaning
 hdivbbmax_max     = 1.0     ! if > overcleanfac, then use B/(h*|div B|) as a coefficient for dtclean;
 !                           ! this is the max value allowed; test suggest =512 for magnetised colliding flows
 beta              = 2.0     ! beta viscosity term
 avdecayconst      = 0.1     ! decay time constant for viscosity switches

 ! radius outside which we kill particles
 rkill             = -1.

 ! dust method
 use_dustfrac = .false.

 ! mcfost
 use_mcfost = .false.
 use_mcfost_stellar_parameters = .false.
 mcfost_computes_Lacc = .false.
 mcfost_dust_subl = .false.
 mcfost_uses_PdV = .true.
 mcfost_keep_part = 0.999_sp
 ISM = 0

 ! radiation
 if (do_radiation) then
    exchange_radiation_energy = .true.
    limit_radiation_flux = .true.
    iopacity_type = 1
    implicit_radiation = .false.
 else
    exchange_radiation_energy = .false.
    limit_radiation_flux = .false.
    iopacity_type = 0
    implicit_radiation = .false.
 endif
 implicit_radiation_store_drad = .false.

 ! variable composition
 use_var_comp = .false.

end subroutine set_default_options

end module options
