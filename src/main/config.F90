!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dim
!
! Module to determine storage based on compile-time configuration
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none
#include "../../build/phantom-version.h"
 integer, parameter, public :: phantom_version_major = PHANTOM_VERSION_MAJOR
 integer, parameter, public :: phantom_version_minor = PHANTOM_VERSION_MINOR
 integer, parameter, public :: phantom_version_micro = PHANTOM_VERSION_MICRO
 character(len=*), parameter, public :: phantom_version_string = PHANTOM_VERSION_STRING

 public

 character(len=80), parameter :: &
    tagline='Phantom v'//phantom_version_string//' (c) 2007-2025 The Authors'

 ! maximum number of particles
 integer :: maxp = 0 ! memory not allocated initially
#ifdef MAXP
 integer, parameter :: maxp_hard = MAXP
#else
 integer, parameter :: maxp_hard = 5200000
#endif

 ! maximum number of point masses
#ifdef MAXPTMASS
 integer, parameter :: maxptmass = MAXPTMASS
#else
 integer, parameter :: maxptmass = 1000
#endif
 integer, parameter :: nsinkproperties = 22

 logical :: store_sf_ptmass = .false.

 ! storage of thermal energy or not
#ifdef ISOTHERMAL
 integer, parameter :: maxvxyzu = 3
 logical, parameter :: isothermal = .true.
#else
 integer, parameter :: maxvxyzu = 4
 logical, parameter :: isothermal = .false.
#endif

 integer :: maxTdust = 0
 logical :: store_dust_temperature = .false.
#ifdef SINK_RADIATION
 logical, parameter :: sink_radiation = .true.
#else
 logical, parameter :: sink_radiation = .false.
#endif

 ! maximum allowable number of neighbours (safest=maxp)
#ifdef MAXNEIGH
 integer, parameter :: maxneigh = MAXNEIGH
#else
 integer, parameter :: maxneigh = maxp_hard
#endif

! maxmimum storage in linklist
 integer         :: ncellsmax
 integer(kind=8) :: ncellsmaxglobal

!------
! Dust
!------
 integer :: maxp_dustfrac = 0
 integer :: maxp_growth = 0
#ifdef DUST
 logical, parameter :: use_dust = .true.

#ifdef MAXDUSTLARGE
 integer, parameter :: maxdustlarge = MAXDUSTLARGE
#else
 integer, parameter :: maxdustlarge = 11
#endif
#ifdef MAXDUSTSMALL
 integer, parameter :: maxdustsmall = MAXDUSTSMALL
#else
 integer, parameter :: maxdustsmall = 11
#endif

#ifdef DUSTGROWTH
 logical, parameter :: use_dustgrowth = .true.
#else
 logical, parameter :: use_dustgrowth = .false.
#endif
#else
 logical, parameter :: use_dust = .false.
 ! integer, parameter :: ndustfluids = 0
 ! integer, parameter :: ndusttypes = 1 ! to avoid seg faults
 integer, parameter :: maxdustlarge = 1
 integer, parameter :: maxdustsmall = 1
 logical, parameter :: use_dustgrowth = .false.
#endif
 integer, parameter :: maxdusttypes = maxdustsmall + maxdustlarge

 ! kdtree
 integer, parameter :: minpart = 10

 integer :: maxprad = 0

 integer, parameter :: radensumforce      = 1, &
                       radenxpartvecforce = 7, &
                       radensumden        = 3, &
                       radenxpartvetden   = 1
#ifdef RADIATION
 logical, parameter :: do_radiation = .true.
#else
 logical, parameter :: do_radiation = .false.
#endif
 ! rhosum
 integer, parameter :: maxrhosum = 39 + &
                                   maxdustlarge - 1 + &
                                   radensumden

 ! fsum
 integer, parameter :: fsumvars = 25 ! Number of scalars in fsum
 integer, parameter :: fsumarrs = 5  ! Number of arrays  in fsum
 integer, parameter :: maxfsum  = fsumvars + &                  ! Total number of values
                                  fsumarrs*(maxdusttypes-1) + &
                                  radensumforce

! xpartveci
 integer, parameter :: maxxpartvecidens = 14 + radenxpartvetden

 integer, parameter :: maxxpartvecvars = 63 ! Number of scalars in xpartvec
 integer, parameter :: maxxpartvecarrs = 2  ! Number of arrays in xpartvec
 integer, parameter :: maxxpartvecGR   = 33 ! Number of GR values in xpartvec (1 for dens, 16 for gcov, 16 for gcon)
 integer, parameter :: maxxpartveciforce = maxxpartvecvars + &              ! Total number of values
                                           maxxpartvecarrs*(maxdusttypes-1) + &
                                           radenxpartvecforce + &
                                           maxxpartvecGR

#ifdef MPI
 logical, parameter :: mpi = .true.
#else
 logical, parameter :: mpi = .false.
#endif

 ! storage for artificial viscosity switch
 integer :: maxalpha = 0
#ifdef DISC_VISCOSITY
 integer, parameter :: nalpha = 0
#else
#ifdef CONST_AV
 integer, parameter :: nalpha = 0
#else
#ifdef USE_MORRIS_MONAGHAN
 integer, parameter :: nalpha = 1
#else
 integer, parameter :: nalpha = 3
#endif
#endif
#endif

 ! optional storage of curl v
#ifdef CURLV
 integer, parameter :: ndivcurlv = 4
#else
 integer, parameter :: ndivcurlv = 1
#endif
 ! storage of velocity derivatives
 integer :: maxdvdx = 0  ! set to maxp when memory allocated

 ! periodic boundaries
#ifdef PERIODIC
 logical, parameter :: periodic = .true.
#else
 logical, parameter :: periodic = .false.
#endif

 ! Maximum number of particle types
 !
 integer, parameter :: maxtypes = 8 + 2*maxdustlarge - 1

 !
 ! Number of dimensions, where it is needed
 ! (Phantom is hard wired to ndim=3 in a lot of
 !  places; changing this does NOT change the
 !  code dimensionality, it just allows routines
 !  to be written in a way that are agnostic to
 !  the number of dimensions)
 !
 integer, parameter :: ndim = 3


!-----------------
! KROME chemistry
!-----------------
 integer :: maxp_krome = 0
#ifdef KROME
 logical, parameter :: use_krome = .true.
#else
 logical, parameter :: use_krome = .false.
#endif

!-----------------
! Magnetic fields
!-----------------
 integer :: maxmhd = 0
#ifdef MHD
 logical, parameter :: mhd = .true.
#else
 logical, parameter :: mhd = .false.
#endif
 integer, parameter :: maxBevol  = 4  ! size of B-arrays (Bx,By,Bz,psi)
 integer, parameter :: ndivcurlB = 4

! Non-ideal MHD
! if fast_divcurlB=true, then divcurlB is calculated simultaneous with density which leads to a race condition and errors (typically less than a percent)
! divcurlB is only used as diagnostics & divergence cleaning in ideal MHD, so fast_divcurlB=true is reasonable
! divcurlB is used to update the non-ideal terms, so fast_divcurlB=false is required for accuracy (especially if there will be jumps in density)
 integer :: maxmhdni = 0
#ifdef NONIDEALMHD
 logical, parameter :: mhd_nonideal    = .true.
 logical, parameter :: fast_divcurlB   = .false.
 integer, parameter :: n_nden_phantom  = 13      ! number density of chemical species, electrons & n_grains; defined in nicil == 11+2*na
#else
 logical, parameter :: mhd_nonideal    = .false.
 logical, parameter :: fast_divcurlB   = .true.
 integer, parameter :: n_nden_phantom  = 0
#endif
 logical            :: calculate_density  = .true.  ! do not toggle; initialised for efficiency
 logical            :: calculate_divcurlB = .true.  ! do not toggle; initialised for efficiency

!--------------------
! H2 Chemistry
!--------------------
 integer :: maxp_h2 = 0
 integer, parameter :: nabundances = 5
 logical :: h2chemistry = .false.

!--------------------
! Self-gravity
!--------------------
 integer :: maxgrav = 0
#ifdef GRAVITY
 logical, parameter :: gravity = .true.
 integer, parameter :: ngradh = 2
#else
 logical, parameter :: gravity = .false.
 integer, parameter :: ngradh = 1
#endif

!--------------------
! General relativity
!--------------------
 integer :: maxgr = 0
#ifdef GR
 logical, parameter :: gr = .true.
 integer, parameter :: maxptmassgr = maxptmass
#else
 logical, parameter :: gr = .false.
 integer, parameter :: maxptmassgr = 0
#endif

!---------------------
! Numerical relativity
!---------------------
#ifdef NR
 logical, parameter :: nr = .true.
#else
 logical, parameter :: nr = .false.
#endif

!--------------------
! Supertimestepping
!--------------------
 integer :: maxsts = 1

!--------------------
! Dust formation
!--------------------
 logical :: do_nucleation  = .false.
 logical :: update_muGamma = .false.
 integer :: itau_alloc     = 0
 integer :: itauL_alloc    = 0
 integer :: inucleation    = 0
 !number of elements considered in the nucleation chemical network
 integer, parameter :: nElements = 10
#ifdef DUST_NUCLEATION
 logical :: nucleation = .true.
#else
 logical :: nucleation = .false.
#endif
 integer :: maxp_nucleation = 0

!--------------------
! MCFOST library
!--------------------
#ifdef MCFOST
 logical, parameter :: compiled_with_mcfost = .true.
#else
 logical, parameter :: compiled_with_mcfost = .false.
#endif

!--------------------
! Light curve stuff
!--------------------
 integer :: maxlum = 0
#ifdef LIGHTCURVE
 logical, parameter :: lightcurve = .true.
#else
 logical, parameter :: lightcurve = .false.
#endif

!--------------------
! logical for bookkeeping
!--------------------
#ifdef INJECT_PARTICLES
 logical, parameter :: inject_parts = .true.
#else
 logical, parameter :: inject_parts = .false.
#endif

!--------------------
! Adaptive Particle Refinement (APR)
!--------------------
#ifdef APR
 logical, parameter :: use_apr = .true.
 integer, parameter :: apr_maxlevel = 10
#else
 logical, parameter :: use_apr = .false.
 integer, parameter :: apr_maxlevel = 0
#endif
 integer :: maxp_apr = 0

!--------------------
! Sink in tree methods
!--------------------
#ifdef SINKTREE
 logical, parameter :: use_sinktree = .true.
#else
 logical, parameter :: use_sinktree = .false.
#endif
 integer :: maxpsph = 0

!--------------------
! individual timesteps
!--------------------
#ifdef IND_TIMESTEPS
 logical, parameter :: ind_timesteps = .true.
#else
 logical, parameter :: ind_timesteps = .false.
#endif

 !--------------------
 ! Analysis array sizes
 !--------------------
 integer :: maxan = 0
 integer :: maxmhdan = 0
 integer :: maxdustan = 0
 integer :: maxgran = 0
 integer :: maxindan = 0

 !--------------------
 ! Phase and gradh sizes - inconsistent with everything else, but keeping to original logic
 !--------------------
 integer :: maxphase = 0
 integer :: maxgradh = 0

 !--------------------
 ! a place to store the number of the dumpfile; required for restart dumps
 !--------------------
 integer :: idumpfile = 0

contains

subroutine update_max_sizes(n,ntot)
 integer,                   intent(in) :: n
 integer(kind=8), optional, intent(in) :: ntot

 maxp = n
 if (use_apr) then
    maxp = 4*n
    maxp_apr = maxp
 endif

 maxpsph = maxp
 if (use_sinktree) then
    maxp = n+maxptmass
 endif

 if (use_krome) maxp_krome = maxp

 if (h2chemistry) maxp_h2 = maxp

#ifdef SINK_RADIATION
 store_dust_temperature = .true.
#endif

 if (store_dust_temperature) maxTdust = maxp
 if (do_nucleation) maxp_nucleation = maxp

#ifdef NCELLSMAX
 ncellsmax       = NCELLSMAX
 ncellsmaxglobal = NCELLSMAX
#else
 ncellsmax = 2*maxp
 if (present(ntot)) then
    ncellsmaxglobal = 2*ntot
 else
    ncellsmaxglobal = ncellsmax
 endif
#endif

 if (use_dust) then
    maxp_dustfrac = maxp
    if (use_dustgrowth) maxp_growth = maxp
 endif

#ifdef DISC_VISCOSITY
 maxalpha = 0
#else
#ifdef CONST_AV
 maxalpha = 0
#else
#ifdef USE_MORRIS_MONAGHAN
 maxalpha = maxp
#else
 maxalpha = maxp
#endif
#endif
#endif

 if (mhd) then
    maxmhd = maxp
    if (mhd_nonideal) maxmhdni = maxp
 endif

 if (gravity) maxgrav = maxp
 if (gr) maxgr = maxp

#ifdef STS_TIMESTEPS
#ifdef IND_TIMESTEPS
 maxsts = maxp
#endif
#endif

#if LIGHTCURVE
 maxlum = maxp
#endif

#ifndef ANALYSIS
 maxan = maxp
 maxmhdan = maxmhd
 maxdustan = maxp_dustfrac
 maxgran = maxgr
#endif

 if (ind_timesteps) maxindan = maxan

 if (do_radiation) then
    maxprad = maxp
    maxlum = maxp
 endif
! Very convoluted, but follows original logic...
 maxphase = maxan
 maxgradh = maxan
 maxdvdx = maxan

end subroutine update_max_sizes

end module dim
