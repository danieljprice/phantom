!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
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
 character(len=80), parameter :: &  ! module version
    modid="$Id$"

 public

 character(len=80), parameter :: &
    tagline='Phantom v'//phantom_version_string//' (c) 2007-2020 The Authors'

 ! maximum number of particles
 integer :: maxp = 0 ! memory not allocated initially
#ifdef MAXP
 integer, parameter :: maxp_hard = MAXP
#else
 integer, parameter :: maxp_hard = 1200000
#endif

 ! maximum number of point masses
#ifdef MAXPTMASS
 integer, parameter :: maxptmass = MAXPTMASS
#else
 integer, parameter :: maxptmass = 1000
#endif
 integer, parameter :: nsinkproperties = 17

 ! storage of thermal energy or not
#ifdef ISOTHERMAL
 integer, parameter :: maxvxyzu = 3
#else
 integer, parameter :: maxvxyzu = 4
#endif

 ! storage of temperature
 integer :: maxtemp = 0
#ifdef STORE_TEMPERATURE
 logical, parameter :: store_temperature = .true.
#else
 logical, parameter :: store_temperature = .false.
#endif

 integer :: maxTdust = 0
#ifdef SINK_RADIATION
 logical, parameter :: sink_radiation = .true.
 logical, parameter :: store_dust_temperature = .true.
#else
 logical, parameter :: sink_radiation = .false.
 logical, parameter :: store_dust_temperature = .false.
#endif

 ! maximum allowable number of neighbours (safest=maxp)
#ifdef MAXNEIGH
 integer, parameter :: maxneigh = MAXNEIGH
#else
 integer, parameter :: maxneigh = maxp_hard
#endif

! maxmimum storage in linklist
 integer :: ncellsmax

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

 integer, parameter :: &
 radensumforce      = 1,&
 radenxpartvecforce = 7,&
 radensumden        = 3,&
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
 integer, parameter :: fsumvars = 20 ! Number of scalars in fsum
 integer, parameter :: fsumarrs = 5  ! Number of arrays  in fsum
 integer, parameter :: maxfsum  = fsumvars + &                  ! Total number of values
                                  fsumarrs*(maxdusttypes-1) + &
                                  radensumforce

! xpartveci
 integer, parameter :: maxxpartvecidens = 14 + radenxpartvetden

 integer, parameter :: maxxpartvecvars = 57 ! Number of scalars in xpartvec
 integer, parameter :: maxxpartvecarrs = 2  ! Number of arrays in xpartvec
 integer, parameter :: maxxpartvecGR   = 33 ! Number of GR values in xpartvec (1 for dens, 16 for gcov, 16 for gcon)
 integer, parameter :: maxxpartveciforce = maxxpartvecvars + &              ! Total number of values
                                           maxxpartvecarrs*(maxdusttypes-1) + &
                                           radenxpartvecforce + &
                                           maxxpartvecGR

 ! cell storage
 integer, parameter :: maxprocs = 32
#ifdef STACKSIZE
 integer, parameter :: stacksize = STACKSIZE
#else
 integer, parameter :: stacksize = 200000
#endif

 ! storage for artificial viscosity switch
 integer :: maxalpha = 0
#ifdef DISC_VISCOSITY
 integer, parameter :: nalpha = 1
#else
#ifdef CONST_AV
 integer, parameter :: nalpha = 1
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
 integer, parameter :: maxtypes = 7 + 2*maxdustlarge - 1

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
 logical, parameter :: store_gamma = .true.
#else
 logical, parameter :: store_gamma = .false.
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
 logical, parameter :: mhd_racecondition = .false. ! true/false: will calculate divcurlB in/after densityiterate permitting/preventing a race condition
 integer, parameter :: maxBevol  = 4  ! size of B-arrays (Bx,By,Bz,psi)
 integer, parameter :: ndivcurlB = 4

! non-ideal MHD
 integer :: maxmhdni = 0
#ifdef MHD
#ifdef NONIDEALMHD
 logical, parameter :: mhd_nonideal = .true.
 integer, parameter :: n_nden_phantom  = 13  ! number density of chemical species, electrons & n_grains; defined in nicil == 11+2*na
#else
 logical, parameter :: mhd_nonideal = .false.
 integer, parameter :: n_nden_phantom  = 0
#endif
#else
 logical, parameter :: mhd_nonideal = .false.
 integer, parameter :: n_nden_phantom  = 0
#endif

!--------------------
! H2 Chemistry
!--------------------
 integer :: maxp_h2 = 0
#ifdef H2CHEM
 logical, parameter :: h2chemistry = .true.
#else
 logical, parameter :: h2chemistry = .false.
#endif
 integer, parameter :: nabundances = 5

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
#else
 logical, parameter :: gr = .false.
#endif

!--------------------
! Gravitational wave strain
!--------------------
#ifdef GWS
 logical, parameter :: gws = .true.
#else
 logical, parameter :: gws = .false.
#endif

!--------------------
! Supertimestepping
!--------------------
 integer :: maxsts = 1

!--------------------
! Wind cooling
!--------------------
#if defined(WIND) || !defined (ISOTHERMAL)
 logical :: windcooling = .true.
#else
 logical :: windcooling = .false.
#endif

!--------------------
! Dust formation
!--------------------
#ifdef NUCLEATION
 integer :: maxsp = maxp_hard
#else
 integer :: maxsp = 0
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
 logical, parameter :: particles_are_injected = .true.
#else
 logical, parameter :: particles_are_injected = .false.
#endif

 !--------------------
 ! Analysis array sizes
 !--------------------
 integer :: maxan = 0
 integer :: maxmhdan = 0
 integer :: maxdustan = 0
 integer :: maxgran = 0

 !--------------------
 ! Phase and gradh sizes - inconsistent with everything else, but keeping to original logic
 !--------------------
 integer :: maxphase = 0
 integer :: maxgradh = 0

contains

subroutine update_max_sizes(n)
 integer, intent(in) :: n

 maxp = n

#ifdef KROME
 maxp_krome = maxp
#endif

#ifdef SINK_RADIATION
 maxTdust = maxp
#endif

#ifdef STORE_TEMPERATURE
 maxtemp = maxp
#endif

#ifdef NCELLSMAX
 ncellsmax = NCELLSMAX
#else
 ncellsmax = 2*maxp
#endif

#ifdef DUST
 maxp_dustfrac = maxp
#ifdef DUSTGROWTH
 maxp_growth = maxp
#endif
#endif

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

#ifdef MHD
 maxmhd = maxp
#ifdef NONIDEALMHD
 maxmhdni = maxp
#endif
#endif

#ifdef H2CHEM
 maxp_h2 = maxp
#endif

#ifdef GRAVITY
 maxgrav = maxp
#endif

#ifdef GR
 maxgr = maxp
#endif

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

#ifdef RADIATION
 maxprad = maxp
#endif
! Very convoluted, but follows original logic...
 maxphase = maxan
 maxgradh = maxan
 maxdvdx = maxan

end subroutine update_max_sizes

end module dim
