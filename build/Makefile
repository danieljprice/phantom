#--------------------------------------------------------------------------!
# The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
# Copyright (c) 2007-2018 The Authors (see AUTHORS)                        !
# See LICENCE file for usage and distribution conditions                   !
# http://users.monash.edu.au/~dprice/phantom                               !
#--------------------------------------------------------------------------!
#+
#  The Phantom Makefile
#
#  DESCRIPTION:
#   This is the main Makefile for all of the code and utilities
#   Compiler settings are grouped under the SYSTEM variable while
#   compile-time settings for different problems are grouped under
#   the SETUP variable
#
#  OWNER: Daniel Price
#
#  $Id: 2788b71b1c08e560e77dce9849c5cb24a668f4b9 $
#+
#--------------------------------------------------------------------------

.KEEP_STATE:

PHANTOM_VERSION_MAJOR=1
PHANTOM_VERSION_MINOR=3
PHANTOM_VERSION_MICRO=0
VERSION=$(PHANTOM_VERSION_MAJOR).$(PHANTOM_VERSION_MINOR).$(PHANTOM_VERSION_MICRO)

KNOWN_SYSTEM=no
SHELL = /bin/bash
VPATH = ${RUNDIR} ../src/main ../src/utils ../src/setup ../src/tests ../src/lib/NICIL/src
BINDIR= ../bin
UNAME=${shell uname}
#----------------------------------------------------------------
# Sensible defaults for phantom configuration
#----------------------------------------------------------------
CONFIG        = config.F90
SETUPFILE     = setup_unifdis.F90
MODFILE       = moddump_default.f90
ANALYSIS      = analysis_dtheader.f90
MULTIRUNFILE  = multirun.f90
LIVE_ANALYSIS = no
LINKLIST      = dtype_kdtree.F90 kdtree.F90 linklist_kdtree.F90
#
# can comment out the following lines and instead set
# the parameters as environment variables
#
ifndef DOUBLEPRECISION
  DOUBLEPRECISION= yes
endif
ifndef EDITOR
  EDITOR= emacs
endif
ifndef OPENMP
  OPENMP= yes
endif
ifndef SPLASH_DIR
   SPLASH_DIR=${shell if [ -d $$HOME/splash ]; then echo $$HOME/splash; fi}
endif
# MPI= yes
#
# endian can be "BIG", "LITTLE" or anything else which has no effect
#
# ENDIAN= default
#
CC = gcc
CCFLAGS = -O5
LIBCXX = -lstdc++
#FPPFLAGS=
LDFLAGS=
SRCPHOTO=

#----------------------------------------------------------------
# here follows specific configuration options used
# for various types of simulations
#
# preprocessor options are as follows:
#
# -DPERIODIC            ! periodic boundaries
# -DIND_TIMESTEPS       ! individual particle timesteps
# -DSTS_TIMESTEPS       ! super-timestepping
# -DDISC_VISCOSITY      ! use artificial disc viscosity ( nu \propto alpha_sph cs h
#                       ! and calculated for both approaching and receding particles
# -DDRIVING             ! use turbulence driving
# -DMHD                 ! magnetic fields
# -DNONIDEALMHD         ! non-ideal magnetic fields including ionisation; uses NICIL
#
# -DSORT                ! whether or not to sort particles
#
# -DPHOTO               ! turn on the photoevaporation
# -DLIGHTCURVE          ! lightcurve estimation
#----------------------------------------------------------------
#
# Check for obsolete setups and replace with generic version
#
#----------------------------------------------------------------
ifeq ($(SETUP), HLTau) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= HLTau
    override SETUP=disc
endif
ifeq ($(SETUP), dustyHLTau) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= dustyHLTau
    override SETUP=dustydisc
endif
ifeq ($(SETUP), mcfost) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP=mcfost
    override SETUP=disc
    MCFOST=yes
endif
ifeq ($(SETUP), planets) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= planets
    override SETUP=disc
endif
ifeq ($(SETUP), binarydisc) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= binarydisc
    override SETUP=disc
endif
ifeq ($(SETUP), dustybinarydisc) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= dustybinarydisc
    override SETUP=dustydisc
endif
ifeq ($(SETUP), Lense-Thirring) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= Lense-Thirring
    override SETUP=disc
endif
ifeq ($(SETUP), warp) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= warp
    override SETUP=disc
endif
ifeq ($(SETUP), rndisc) # [buildbot skip]
    OBSOLETE_SETUP=yes
    OLDSETUP= rndisc
    override SETUP=lightcurvedisc
endif

#----------------------------------------------------------------
# Current code setup options
#----------------------------------------------------------------
ifeq ($(SETUP), empty)
#   empty setup for external-driver simulation
    SETUPFILE= setup_empty.f90
    IND_TIMESTEPS=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), wddisc)
    ISOTHERMAL=yes
    SETUPFILE=setup_wddisc.f90
    KNOWN_SETUP=yes
    DUST=yes
endif

ifeq ($(SETUP), asteroidwind)
    SETUPFILE=setup_asteroidwind.f90
    SRCINJECT=inject_asteroidwind.f90
    IND_TIMESTEPS=yes
    CONST_AV=yes
    ISOTHERMAL=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), galdisc)
#   galactic disc simulations
    IND_TIMESTEPS=yes
    H2CHEM=yes
    ISOTHERMAL=no
    GRAVITY=no
    MHD=no
    SETUPFILE= setup_galdisc.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), galdiscmhd)
#   galactic disc simulations with magnetic fields
    IND_TIMESTEPS=yes
    H2CHEM=no
    ISOTHERMAL=yes
    GRAVITY=no
    MHD=yes
    SETUPFILE= setup_galdisc.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), turbdrive)
#   driven turbulence
    ifeq ($(IND_TIMESTEPS), yes)
       FPPFLAGS= -DPERIODIC -DCORRECT_BULK_MOTION -DSTIR_FROM_FILE
    else
       FPPFLAGS= -DPERIODIC # -DCORRECT_MEAN_FORCE
    endif
    SETUPFILE= setup_unifdis.F90
    SRCTURB= forcing.F90
    MULTIRUNFILE= multirun_mach.f90
    KNOWN_SETUP=yes
    CURLV=yes
    ISOTHERMAL=yes
endif

ifeq ($(SETUP), turbsubsonic)
#   subsonic, driven turbulence
    FPPFLAGS= -DPERIODIC -DSTIR_FROM_FILE -DCORRECT_MEAN_FORCE
    SETUPFILE= setup_unifdis.F90
    SRCTURB= forcing.F90
    MULTIRUNFILE= multirun_mach.f90
    KNOWN_SETUP=yes
    KERNEL=quintic
    ISOTHERMAL=yes
    IND_TIMESTEPS=no
    CURLV=yes
    USE_STRAIN_TENSOR=yes
endif

ifeq ($(SETUP), taylorgreen)
#   Taylor-Green vortex problem
    FPPFLAGS= -DPERIODIC -DCURLV -DUSE_STRAIN_TENSOR
    SETUPFILE= setup_taylorgreen.f90
    ISOTHERMAL=yes
    KNOWN_SETUP=yes
    KERNEL=quintic
    MODFILE= moddump_taylorgreen.f90
    IND_TIMESTEPS=no
endif

ifeq ($(SETUP), turb)
#   driven supersonic turbulence (hydro, mhd, dusty)
    FPPFLAGS      = -DPERIODIC -DCORRECT_BULK_MOTION -DSTIR_FROM_FILE
    SETUPFILE     = setup_turb.F90
    SRCTURB       = forcing.F90
    IND_TIMESTEPS = yes
    KNOWN_SETUP   = yes
    ISOTHERMAL    = yes
    CURLV         = yes
    MHD           = no
    DUST          = no
    USE_STRAIN_TENSOR = yes
endif

ifeq ($(SETUP), wd)
    SETUPFILE     = setup_star.f90
    IND_TIMESTEPS = no
    MHD           = no
    GRAVITY       = yes
    ISOTHERMAL    = no
    KNOWN_SETUP   = yes
    STORE_TEMP    = yes
    MODFILE       = moddump_binarystar.f90
endif

ifeq ($(SETUP), photoevap)
# Mark Hutchison
    FPPFLAGS= -DDISC_VISCOSITY -DPHOTO
    SETUPFILE= setup_photoevap.f90
    ANALYSIS= analysis_disc.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    SRCPHOTO=photoevap.f90
endif

ifeq ($(SETUP), disc)
#   locally isothermal gas disc
    ifeq ($(DISC_VISCOSITY), no)
        FPPFLAGS=
    else
        FPPFLAGS= -DDISC_VISCOSITY
    endif
    SETUPFILE= setup_disc.f90
    ANALYSIS= analysis_disc.f90
    ISOTHERMAL=yes
    KNOWN_SETUP=yes
    MULTIRUNFILE= multirun.f90
    IND_TIMESTEPS=yes
endif

ifeq ($(SETUP), adiabaticdisc)
#   adiabatic disc
    FPPFLAGS= -DDISC_VISCOSITY
    SETUPFILE= setup_disc.f90
    ANALYSIS= analysis_disc.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    ISOTHERMAL=no
    MULTIRUNFILE= multirun.f90
endif

ifeq ($(SETUP), lightcurvedisc)
#   adiabatic disc with lightcurve
    FPPFLAGS= -DDISC_VISCOSITY -DLIGHTCURVE
    SETUPFILE= setup_disc.f90
    ANALYSIS= analysis_disc.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    ISOTHERMAL=no
    MULTIRUNFILE= multirun.f90
endif

ifeq ($(SETUP), gwdisc)
#   locally isothermal discs
    FPPFLAGS= -DDISC_VISCOSITY
    SETUPFILE= setup_gwdisc.f90
    ANALYSIS= analysis_disc.f90
    MAXP=2000000
    IND_TIMESTEPS=yes
    ISOTHERMAL=yes
    KNOWN_SETUP=yes
    MULTIRUNFILE= multirun.f90
    SRCPOT= ${SRCPOTS:extern_binary.f90=extern_binary_gw.f90}
endif

ifeq ($(SETUP), nshwdisc)
#   disc around a neutron star
    FPPFLAGS= -DDISC_VISCOSITY -DPRDRAG
    SETUPFILE= setup_nsdisc.f90
    ANALYSIS= analysis_disc.f90
    MODFILE= moddump_changemass.f90
    ISOTHERMAL=yes
    IND_TIMESTEPS=yes
    KNOWN_SETUP=yes
    NCELLSMAX=3*maxp
endif

ifeq ($(SETUP), prtest)
#   simple test of prdrag
    FPPFLAGS=
    SETUPFILE= setup_prtest.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), binarydiscMFlow)
#   binarydiscMFlow setup
    FPPFLAGS= -DMFLOW -DUSE_STRAIN_TENSOR #-DVMFLOW  #-DDISC_VISCOSITY# -DBINPOS
    SETUPFILE= setup_disc.f90
    ANALYSIS= analysis_disc_MFlow.f90
#    ANALYSIS= analysis_binarydisc.f90
    MAXP=1000000
    ISOTHERMAL=yes
    CURLV=yes
    KNOWN_SETUP=yes
    LIVE_ANALYSIS=no
    IND_TIMESTEPS=yes
    MODFILE= moddump_removeparticles_cylinder.f90 #moddump_addpartfortest.f90
endif

ifeq ($(SETUP), planetdisc)
#   planet disc interaction with fixed planet orbit
    FPPFLAGS=
    SETUPFILE= setup_planetdisc.f90
    ISOTHERMAL=yes
    IND_TIMESTEPS=yes
    CURLV=yes
    KNOWN_SETUP=yes
    ANALYSIS=analysis_disc.f90
endif

ifeq ($(SETUP), planetatm)
#   disc interaction with fixed planet orbit + atmosphere
    FPPFLAGS=
    SETUPFILE= setup_disc.f90
    ISOTHERMAL=yes
    IND_TIMESTEPS=yes
    CURLV=yes
    KNOWN_SETUP=yes
    ANALYSIS=analysis_disc.f90
endif

ifeq ($(SETUP), torus)
#   MRI torus
    FPPFLAGS=
    SETUPFILE= setup_torus.f90
    ANALYSIS= analysis_torus.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), galcen)
#   galactic centre
    FPPFLAGS=
    SETUPFILE= setup_galcen_stars.f90
    SRCINJECT= inject_galcen_winds.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), quebec)
    SETUPFILE = setup_quebec.f90
    GRAVITY = yes
    KNOWN_SETUP = yes
    MODFILE = moddump_binarystar.f90
endif

ifeq ($(SETUP), tde)
#   tidal disruption simulations
    SETUPFILE= setup_star.f90
    ANALYSIS=analysis_tde.f90
    GRAVITY=yes
    ISOTHERMAL=yes
    MODFILE=moddump_tidal.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), polytrope)
#   single (option 2) or binary (option 3) polytrope test
    SETUPFILE= setup_star.f90
    ANALYSIS= density_profiles.o analysis_polytropes.f90
    GRAVITY=yes
    ISOTHERMAL=yes
    MODFILE=moddump_binarystar.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), neutronstar)
#   neutron star (use option 4)
    SETUPFILE= setup_star.f90
    ISOTHERMAL=yes
    GRAVITY=no     #since external force being used
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), sphereinbox)
#   sphere-in-box setup
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_sphereinbox.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), shock)
#   sod shock tube test
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_shock.F90
    DOUBLEPRECISION=yes
    KERNEL=quintic
    ISOTHERMAL=no
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), dustyshock)
#   sod shock tube test
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_shock.F90
    DUST=yes
    DOUBLEPRECISION=yes
    KERNEL=quintic
    ISOTHERMAL=no
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), mhdshock)
#   Ryu & Brio-Wu shock tube tests
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_shock.F90
    DOUBLEPRECISION=yes
    MHD=yes
    KERNEL=quintic
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), nimhdshock)
#   non-ideal mhd standing and C shock tests
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_shock.F90
    DOUBLEPRECISION=yes
    MHD=yes
    IND_TIMESTEPS=no
    STS_TIMESTEPS=no
    NONIDEALMHD=yes
    KERNEL=WendlandC4
    ISOTHERMAL=yes
    MAXP=6000000
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), dustydisc)
#   locally isothermal dusty discs
    FPPFLAGS= -DDISC_VISCOSITY
    SETUPFILE= setup_disc.f90
    MODFILE= moddump_dustadd.f90
    ISOTHERMAL=yes
    DUST=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    ANALYSIS=analysis_dustydisc.f90
endif

ifeq ($(SETUP), growingdisc)
#   locally isothermal dusty discs with growth and fragmentation
    FPPFLAGS= -DDISC_VISCOSITY
    SETUPFILE= setup_disc.f90
    MODFILE= moddump_dustadd.f90
    ISOTHERMAL=yes
    DUST=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    DUSTGROWTH = yes
    ANALYSIS=analysis_dustydisc.f90
endif

ifeq ($(SETUP), dustybox)
#   dust in a box
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_dustybox.f90
    MODFILE= moddump_dustadd.f90
    ISOTHERMAL=yes
    DUST=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    ANALYSIS= analysis_trackbox.f90
endif

ifeq ($(SETUP), dustysedov)
#   Sedov blast wave test with dust
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_dustysedov.f90
    MODFILE= moddump_dustadd.f90
    DUST=yes
    KNOWN_SETUP=yes
    #IND_TIMESTEPS=no
endif

ifeq ($(SETUP), dustywave)
#   dust in a box
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_wave.f90
    MODFILE= moddump_dustadd.f90
    DUST=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    ANALYSIS= analysis_trackbox.f90
endif

ifeq ($(SETUP), wave)
#   linear wave
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_wave.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    KERNEL=quintic
endif

ifeq ($(SETUP), wavedamp)
#   Wave damping test as per Choi et al (2009)
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_wavedamp.f90
    ISOTHERMAL=yes
    NONIDEALMHD=yes
    MHD=yes
    KNOWN_SETUP=yes
    KERNEL=WendlandC4
    IND_TIMESTEPS=no
    STS_TIMESTEPS=no
    ANALYSIS = analysis_bzrms.f90
    DEBUG=no
endif

ifeq ($(SETUP), sedov)
#   Sedov blast wave test
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_sedov.f90
    IND_TIMESTEPS=yes
    KNOWN_SETUP=yes
    MAXP=2100000
endif

ifeq ($(SETUP), blob)
#   Blob evaporation problem
    FPPFLAGS= -DPERIODIC#-DSORT
    SETUPFILE= setup_blob.f90
    DOUBLEPRECISION= no
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), kh)
#   Kelvin-Helmholtz problem
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_kh.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), mhdrotor)
#   MHD rotor problem
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_mhdrotor.f90
    MHD=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), jadvect)
#   MHD current loop advection problem
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_jadvect.f90
    MHD=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), alfven)
#   MHD circularly polarised Alfven wave problem
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_alfvenwave.f90
    MHD=yes
    KNOWN_SETUP=yes
    KERNEL=quintic
endif

ifeq ($(SETUP), orstang)
#   Orszag-Tang vortex
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_orstang.f90
    MHD=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), balsarakim)
#   Balsara-Kim 2004
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_unifdis.f90
    MHD=yes
    KNOWN_SETUP=yes
    SRCINJECT=inject_sne.f90
    KERNEL=quintic
    IND_TIMESTEPS=yes
    H2CHEM=yes
endif

ifeq ($(SETUP), mhdvortex)
#   Balsara (2004) MHD vortex
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_mhdvortex.f90
    MHD=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), mhdsine)
#   MHD sine wave
    FPPFLAGS= -DPERIODIC
    SETUPFILE= setup_mhdsine.f90
    NONIDEALMHD=no
    MHD=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), mhdblast)
#   MHD blast wave test
    SETUPFILE= setup_mhdblast.f90
    PERIODIC=yes
    MHD=yes
    ISOTHERMAL=no
    IND_TIMESTEPS=no
    KNOWN_SETUP=yes
    MAXP=3000000
endif

ifeq ($(SETUP), mhdwave)
#   propagating isolated MHD wave
    SETUPFILE= setup_mhdwave.f90
    PERIODIC=yes
    MHD=yes
    ISOTHERMAL=no
    IND_TIMESTEPS=no
    KNOWN_SETUP=yes
    MAXP=3000000
endif

ifeq ($(SETUP), cluster)
#   cluster formation (setup)
    FPPFLAGS=
    SETUPFILE= velfield_fromcubes.f90 setup_cluster.f90
    MODFILE= moddump_default.f90
    ANALYSIS= phantom_pdfs.o analysis_MWpdf.f90 #analysis_sinkmass.f90
    ISOTHERMAL=yes
    MHD=no
    GRAVITY=yes
    IND_TIMESTEPS=yes
    KNOWN_SETUP=yes
    MAXPTMASS=1000
    MAXP=3500000
endif

ifeq ($(SETUP), binary)
#   binary setup
    FPPFLAGS= -DCONST_AV
    SRCINJECT= inject_rochelobe.f90
    SETUPFILE= setup_binary.f90
    #SETUPFILE= setup_chinchen.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), common)
#   binary setup
    FPPFLAGS=
    SETUPFILE= setup_common.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), star)
#   import stellar model from 1D stellar evolution code
#   use option 5 of setup_star
    FPPFLAGS=
    SETUPFILE= setup_star.f90
    MODFILE= moddump_binary.f90
    ANALYSIS= ${SRCNIMHD} utils_summary.o utils_omp.o ptmass.o energies.o analysis_common_envelope.F90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    MAXP=10000000
    GRAVITY=yes
    MHD=no
endif

ifeq ($(SETUP), isowind)
#   wind setup (isothermal spherical wind from a sink particle)
    FPPFLAGS=
    ISOTHERMAL=yes
    SETUPFILE= setup_wind.f90
    SRCINJECT= icosahedron.f90 utils_inject.f90 inject_wind.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    STS_TIMESTEPS=no
endif

ifeq ($(SETUP), wind)
#   wind setup (adiabatic supersonic spherical wind from a sink particle)
    SETUPFILE= setup_wind.f90
    SRCINJECT= icosahedron.f90 utils_inject.f90 inject_wind.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    STS_TIMESTEPS=no
endif

ifeq ($(SETUP), bowen)
# bowen wind setup (pulsating dusty wind)
    FPPFLAGS= -DBOWEN
    SETUPFILE= setup_wind.f90
    SRCINJECT= icosahedron.f90 utils_inject.f90 bowen_dust.F90 inject_wind.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    STS_TIMESTEPS=no
endif

ifeq ($(SETUP), BHL)
# Bondi-Hoyle-Lyttleton setup
    SETUPFILE= setup_BHL.f90
    SRCINJECT= inject_BHL.f90
    FPPFLAGS= -DNOWARNRESTRICTEDHJUMP
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
endif

ifeq ($(SETUP), jet)
#   Jet simulation from Price, Tricco & Bate (2012)
    FPPFLAGS= -DPERIODIC -DGRAVITY
    SETUPFILE= setup_sphereinbox.f90
    #ANALYSIS= analysis_jet.f90
    ANALYSIS= ${SRCNIMHD} analysis_protostar_environ.F90
    ISOTHERMAL=yes
    MHD=yes
    IND_TIMESTEPS=yes
    KNOWN_SETUP=yes
    DUST=no
endif

ifeq ($(SETUP), jetnimhd)
#   Simulation from Wurster, Price & Bate (2016,2017) et seq
    SETUPFILE= setup_sphereinbox.f90
    ANALYSIS= ${SRCNIMHD} analysis_protostar_environ.F90
    PERIODIC=yes
    GRAVITY=yes
    ISOTHERMAL=yes
    MHD=yes
    NONIDEALMHD=yes
    IND_TIMESTEPS=yes
    STS_TIMESTEPS=yes
    MODFILE=moddump_CoM.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), sgdisc)
#   self-gravitating disc
    IND_TIMESTEPS=yes
    GRAVITY=yes
    SETUPFILE= setup_disc.f90
#   ANALYSIS = ${LINKLIST} utils_omp.F90 utils_summary.F90 ptmass.F90 analysis_clumpfind.F90
    ANALYSIS = analysis_disc_stresses.f90
#    ANALYSIS = ${LINKLIST} analysis_getneighbours.f90
    DIMFILE= dim_disc.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), dustysgdisc)
#   self-gravitating dustydisc
    SETUPFILE= setup_disc.f90
    GRAVITY=yes
    DUST=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    ANALYSIS=analysis_dustydisc.f90
endif

ifeq ($(SETUP), sgdiscsector)
#   self-gravitating disc sector
    IND_TIMESTEPS=yes
    GRAVITY=yes
    SETUPFILE= setup_discsector.f90
    SRCINJECT= inject_keplerianshear.f90
    DIMFILE= dim_disc.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), dustsettle)
#   dust settling test from PL15
    SETUPFILE= setup_dustsettle.f90
    DUST=yes
    PERIODIC=yes
    ISOTHERMAL=yes
    MODFILE=moddump_dustadd.f90
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
endif

ifeq ($(SETUP), test)
#   default setup for tests
    FPPFLAGS= -DPERIODIC
    KNOWN_SETUP=yes
    CONST_ARTRES=yes
    CURLV=yes
    MHD=yes
    DUST=yes
    KERNEL=cubic
endif

ifeq ($(SETUP), test2)
#   default setup for tests
    FPPFLAGS= -DDISC_VISCOSITY
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
    USE_STRAIN_TENSOR=yes
endif

ifeq ($(SETUP), testcyl)
#   default setup for tests
    FPPFLAGS= -DDISC_VISCOSITY
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    CONST_ARTRES=yes
    CURLV=yes
endif

ifeq ($(SETUP), testkd)
#   default setup for tests
    FPPFLAGS= -DPERIODIC
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    CONST_ARTRES=yes
    CURLV=yes
    MHD=yes
endif

ifeq ($(SETUP), testgrav)
#   self-gravity unit tests
    FPPFLAGS= -DGRAVITY
    KNOWN_SETUP=yes
    CONST_ARTRES=yes
    CURLV=yes
endif

ifeq ($(SETUP), testdust)
#   dust unit tests
    PERIODIC=yes
    DUST=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
endif

ifeq ($(SETUP), testgrowth)
#   dust growth unit tests
    PERIODIC=yes
    DUST=yes
    DUSTGROWTH=yes
    KNOWN_SETUP=yes
    IND_TIMESTEPS=no
endif

ifeq ($(SETUP), testnimhd)
#   non-ideal MHD (+ boundary particle + super-timestepping) unit tests
    PERIODIC=yes
    ISOTHERMAL=yes
    NONIDEALMHD=yes
    MHD=yes
    KERNEL=WendlandC4
    IND_TIMESTEPS=no
    STS_TIMESTEPS=yes
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), testlum)
#   Lense-Thirring setup
    FPPFLAGS= -DLIGHTCURVE -DUSE_STRAIN_TENSOR
    KNOWN_SETUP=yes
    IND_TIMESTEPS=yes
    ISOTHERMAL=no
endif

ifeq ($(SETUP), default)
    KNOWN_SETUP=yes
    SETUPFILE= setup_unifdis.F90
    FPPFLAGS= -DPERIODIC
    DUST=yes
endif

#ifeq ($(SETUP), clumpybox) # [buildbot skip]
#    KNOWN_SETUP=yes
#    ANALYSIS = ${LINKLIST} analysis_velocitydispersion_vs_scale.f90
#    ANALYSIS = adaptivemesh.f90 interpolate3D_amr.f90 analysis_pdfs_dhf.f90
#    SETUPFILE= setup_fromgrid.f90
#endif

ifeq ($(SETUP), galaxies)
#   Read in data created from Wurster&Thacker(2013a,b)
    SETUPFILE= setup_galaxies.f90
    ANALYSIS=analysis_GalMerger.f90
    CONST_AV=no
    ISOTHERMAL=no
    IND_TIMESTEPS=yes
    GRAVITY=yes
    MAXP=2600000
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), nsmerger)
#   Model a neutron star merger; use option 6
    SETUPFILE= setup_star.f90
    ISOTHERMAL=yes
    IND_TIMESTEPS=no
    GRAVITY=yes
    MODFILE=moddump_binarystar.f90
    ANALYSIS=analysis_NSmerger.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), evrard)
#   models the Evrard collapse; use option 7
    SETUPFILE= setup_star.f90
    ISOTHERMAL=no
    GRAVITY=yes
    ANALYSIS=analysis_sphere.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), tokamak)
#   tokamak torus setup
    PERIODIC=no
    ISOTHERMAL=yes
    SETUPFILE= setup_tokamak.f90
    KNOWN_SETUP=yes
endif

ifeq ($(SETUP), sfmcfost) # [buildbot skip]
#   live feedback from mcfost in star formation calculation
    PERIODIC=yes
    SETUPFILE= setup_sphereinbox.f90
    ANALYSIS= analysis_mcfost.f90
    LIVE_ANALYSIS=yes
    ISOTHERMAL=no
    KNOWN_SETUP=yes
    MULTIRUNFILE= multirun.f90
    IND_TIMESTEPS=yes
    GRAVITY=yes
    MCFOST=yes
endif

ifeq ($(SETUP), mcfostcmdline) # [buildbot skip]
#   live feedback from mcfost, superseded by mcfost setup
    FPPFLAGS= -DDISC_VISCOSITY
    SETUPFILE= setup_disc.f90
    ANALYSIS= analysis_mcfostcmdline.f90
    LIVE_ANALYSIS=yes
    ISOTHERMAL=no
    KNOWN_SETUP=yes
    MULTIRUNFILE= multirun.f90
    IND_TIMESTEPS=yes
endif

ifndef SETUPFILE
    SETUPFILE= setup_unifdis.F90
endif

ifndef SRCNIMHD
    SRCNIMHD = nicil.F90 nicil_supplement.F90
endif

ifndef SRCDUST
    SRCDUST = dust.F90 growth.F90
endif

#ifndef SRCGROWTH
#    SRCGROWTH = growth.F90
#endif

ifeq ($(MCFOST), yes)
#   live feedback from mcfost
    ANALYSIS= analysis_mcfost.f90
    LIVE_ANALYSIS=yes
    ISOTHERMAL=no
    MCFOST=yes

    FPPFLAGS+= -DMCFOST
    LDFLAGS+= -L$(MCFOST_DIR)/src -lmcfost $(LIBCXX)
#    LDFLAGS+=  ~/mcfost/src/libtmp/*.o  ~/mcfost/src/*.o  $(LIBCXX)
endif


#----------------------------------------------------------------
ifeq ($(SYSTEM),cray)
    FC=ftn
    FFLAGS=-Oaggress -Ovector3 -Oipa4
    DBLFLAG= -s real64
    CC=cc
    CCFLAGS=-O3
    KNOWN_SYSTEM=yes
ifeq ($(MAP),yes)
    LDFLAGS+=-dynamic -L/${ALLINEA_DIR}/allinea -lmap-sampler -Wl,--eh-frame-hdr
    FFLAGS+= -G2
endif
endif

ifeq ($(SYSTEM),daint)
    FC=ftn
    FFLAGS= -O3 -inline-factor=500 -dynamic -mcmodel=medium -heap-arrays -shared-intel -warn uninitialized -warn unused -warn truncated_source
    DBLFLAG= -r8
    DEBUGFLAG= -check all -WB -traceback -g -debug all
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
    OMPFLAGS = -openmp
    CC = icc
    CCFLAGS = -O3 -mcmodel=medium
    LIBCXX = -cxxlib
    QSYS = slurm
ifeq ($(MPI),daint)
    USEMPI=yes
    FPPFLAGS += -DMPI
endif
endif

ifeq ($(SYSTEM), msg)
    include Makefile_defaults_ifort
    QSYS = sge
    QSHELL = tcsh
    ifeq ($(OPENMP),yes)
       QPE = smp
       NOMP = '$$NSLOTS'
       ifndef NPAR
          NPAR = '4-32'
       endif
    endif
    ifeq ($(MPI),yes)
       QPE = mpi
       ifeq ($(OPENMP),yes)
            QPE = mqu4
            NOMP = 4
       endif
    endif
    #QEXTRA='-l dpod=true -q mqu2'
    #HDF5=yes
#    HDF5ROOT=/opt/sw/hdf5-1.8.0/
endif

ifeq ($(SYSTEM), m2)
#   MASSIVE facility: massive.org.au
    include Makefile_defaults_ifort
    QSYS = pbs
    ifeq ($(OPENMP),yes)
       NOMP='12'
    else
       NOMP='1'
    endif
    QNODES='nodes='$(NMPI)':ppn='$(NOMP)
    WALLTIME='500:00:00'
endif

ifeq ($(SYSTEM), g2)
#   gstar facility
#   Note: gstar has nomp=12; sstar has nomp=16
    include Makefile_defaults_ifort
    QSYS = pbs
    ifeq ($(OPENMP),yes)
       NOMP='16'
    else
       NOMP='1'
    endif
    QNAME='sstar'
    QNODES='nodes='$(NMPI)':ppn='$(NOMP)
    WALLTIME='168:00:00'
    MPIEXEC='mpiexec -npernode 1'
endif

ifeq ($(SYSTEM), ozstar)
#   ozstar facility
    include Makefile_defaults_ifort
    OMPFLAGS=-qopenmp
    NOMP=32
    #QNAME='skylake'
    QSYS = slurm
    WALLTIME='168:00:00'
endif

ifeq ($(SYSTEM), monarch)
    include Makefile_defaults_ifort
    OMPFLAGS=-qopenmp -qopt-report
    QSYS = slurm
    QPROJECT='p01'
    WALLTIME='100:59:59'
    QPARTITION='comp'
endif

ifeq ($(SYSTEM), monarchpsxe)
    include Makefile_defaults_ifort
    QSYS = slurm
    QPROJECT='p01'
endif

ifeq ($(SYSTEM), raijin)
#   raijin (NCI machine)
    include Makefile_defaults_ifort
    #MPI=intel
    FFLAGS= -O3 -mcmodel=medium -shared-intel -ip -axSSE2,SSSE3,SSE4.1,SSE4.2,AVX -inline-factor=500 -warn uninitialized -warn unused -warn truncated_source
    DEBUGFLAG= -check all -WB -traceback -g -fpe0 -fpstkchk
    CCFLAGS= -O3 -ip
    QSYS= pbs
    #PBSRESUBMIT=yes
    NOMP=16
    ifeq ($(MPI),yes)
       NPAR=32
    endif
    QPROJECT='pt4'
    QNAME='normal'
    WALLTIME='48:00:00'
    MPIEXEC='mpiexec -npernode 1'
    QNODES='ncpus='$(NPAR)
    QEXTRA='-l other=hyperthread'
endif

ifeq ($(SYSTEM), gfortran)
    include Makefile_defaults_gfortran
ifneq ($(UNAME), Darwin)
    FFLAGS+= -mcmodel=medium
endif
endif

ifeq ($(SYSTEM), gfortranOSX)  # for use with mac gfortran (5.3.0, 7.3.0 tested)
    include Makefile_defaults_gfortran
endif

ifeq ($(SYSTEM), gfortran44)
    include Makefile_defaults_gfortran
    FC= gfortran -gdwarf-2
    FFLAGS= -O3 -Wall -frecord-marker=4 -finline-functions-called-once -finline-limit=1500 -funroll-loops -ftree-vectorize
    DEBUGFLAG= -g -frange-check -ffpe-trap=invalid,denormal -finit-real=nan -finit-integer=nan -fbacktrace
endif

ifeq ($(SYSTEM), gfortran47)
    include Makefile_defaults_gfortran
    FC= gfortran-mp-4.7 -gdwarf-2
    FFLAGS= -Wall -m64 -O3 -ffast-math -funroll-loops -ftree-loop-linear \
            -finline-functions-called-once \
            -fomit-frame-pointer -finline-limit=3000 --param min-vect-loop-bound=2
    DEBUGFLAG= -Wextra -g -frange-check -fcheck=all -ffpe-trap=denormal -finit-real=nan -finit-integer=nan -fbacktrace
endif

ifeq ($(SYSTEM), complexity)
#   complexity.leicester.dirac.ac.uk
    include Makefile_defaults_ifort
    FFLAGS= -O3 -xhost -ipo -mcmodel=medium -shared-intel -warn uninitialized \
            -warn unused -warn truncated_source
    DEBUGFLAG= -check all -WB -traceback -g -fpe0 -fp-stack-check
    CCFLAGS = -O3 -ipo -mcmodel=medium
    QSYS=pbs
    QNAME=q64
    WALLTIME='48:00:00'
endif

ifeq ($(SYSTEM), isca)
    # local cluster at the University of Exeter
    include Makefile_defaults_ifort
    FFLAGS= -O3 -axAVX -mcmodel=medium \
            -warn uninitialized -warn truncated_source\
            -warn interfaces -nogen-interfaces
    OMPFLAGS= -qopenmp
    DEBUGFLAG= -check all -traceback -g -fpe0 -fp-stack-check -heap-arrays -O0
    QNAME=pq
    WALLTIME='168:00:00'
endif

ifeq ($(SYSTEM), skylake)
# HPCs Skylake cluster at Cambridge
    include Makefile_defaults_ifort
    FFLAGS= -O3 -mcmodel=medium -shared-intel -warn uninitialized -warn unused -warn \
            truncated_source -xCORE-AVX512 -ipo
    OMPFLAGS = -qopenmp
    CCFLAGS = -O3 -mcmodel=medium -xCORE-AVX512 -ipo
    QSYS = slurm
    QPROJECT='DIRAC-DP005-CPU'
    WALLTIME='36:00:00'
endif

ifeq ($(SYSTEM), ifort)
    include Makefile_defaults_ifort
endif

ifeq ($(SYSTEM), ifortmac)
    include Makefile_defaults_ifort
    FFLAGS= -O3 -xhost -shared-intel -warn uninitialized \
            -warn unused -warn truncated_source -Wl,-rpath,/opt/intel/lib
    DEBUGFLAG= -check all -WB -traceback -g -fpe0 -fp-stack-check
endif

ifeq ($(SYSTEM), ifortgcc)
    include Makefile_defaults_ifort
    CC = gcc
    CCFLAGS = -O3
endif

ifeq ($(SYSTEM), hydra)
# this configuration works for the hydra cluster http://www.mpcdf.mpg.de/services/computing/hydra
    include Makefile_defaults_ifort
    FFLAGS= -O3 -xavx -ip -mcmodel=medium -shared-intel -warn uninitialized \
            -warn unused -warn truncated_source
    DEBUGFLAG= -check all -WB -traceback -g -fpe0 -fp-stack-check
    CCFLAGS = -O3 -ipo -mcmodel=medium
endif

ifeq ($(SYSTEM), lyoccf)
# LIO CCF cluster
    include Makefile_defaults_ifort
    FFLAGS= -O3 -ftz -xavx -cpp -sox -fno-alias -fno-fnalias \
            -no-prec-div -no-prec-sqrt -align all -warn uninitialized \
            -warn unused -warn truncated_source
    LIBCXX = -cxxlib
endif

# Set some default files if not defined above
ifdef MAXP
   FPPFLAGS += -DMAXP=${MAXP}
endif
ifdef MAXPTMASS
   FPPFLAGS += -DMAXPTMASS=${MAXPTMASS}
endif
ifdef MAXNEIGH
   FPPFLAGS += -DMAXNEIGH=${MAXNEIGH}
endif
ifdef NCELLSMAX
   FPPFLAGS += -DNCELLSMAX=${NCELLSMAX}
endif
ifdef STACKSIZE
   FPPFLAGS += -DSTACKSIZE=${STACKSIZE}
endif
# Set other optional flags depending on settings

ifeq ($(DEBUG), yes)
    FFLAGS += ${DEBUGFLAG}
    FFLAGS := $(FFLAGS:-O3=-O0)
    FFLAGS := $(FFLAGS:-ipo= )
endif

ifeq ($(ENDIAN), BIG)
    FFLAGS += ${ENDIANFLAGBIG}
endif

ifeq ($(ENDIAN), LITTLE)
    FFLAGS += ${ENDIANFLAGLITTLE}
endif

ifeq ($(OPENMP), yes)
    FFLAGS += ${OMPFLAGS}
endif

ifeq ($(SORT), yes)
    FPPFLAGS += -DSORT
endif

ifeq ($(PERIODIC), yes)
    FPPFLAGS += -DPERIODIC
endif

ifeq ($(GRAVITY), yes)
    FPPFLAGS += -DGRAVITY
endif

ifeq ($(ISOTHERMAL), yes)
    FPPFLAGS += -DISOTHERMAL
endif

ifeq ($(STORE_TEMP), yes)
    FPPFLAGS += -DSTORE_TEMPERATURE
endif

ifeq ($(MHD), yes)
    FPPFLAGS += -DMHD
endif

ifeq ($(DUST), yes)
    FPPFLAGS += -DDUST
    ifndef KERNEL
       KERNEL=quintic
    endif
endif

ifdef MAXDUSTSMALL
   FPPFLAGS += -DMAXDUSTSMALL=${MAXDUSTSMALL}
endif
ifdef MAXDUSTLARGE
   FPPFLAGS += -DMAXDUSTLARGE=${MAXDUSTLARGE}
endif

ifeq ($(DUSTGROWTH), yes)
    FPPFLAGS += -DDUSTGROWTH
endif

ifeq ($(NONIDEALMHD), yes)
    FPPFLAGS += -DNONIDEALMHD
endif

ifeq ($(H2CHEM), yes)
    FPPFLAGS += -DH2CHEM
endif

ifeq ($(DISC_VISCOSITY), yes)
    FPPFLAGS += -DDISC_VISCOSITY
endif

ifeq ($(CONST_AV), yes)
    FPPFLAGS += -DCONST_AV
endif

ifeq ($(MORRIS_MONAGHAN), yes)
    FPPFLAGS += -DUSE_MORRIS_MONAGHAN
endif

ifeq ($(CONST_ARTRES), yes)
    FPPFLAGS += -DCONST_ARTRES
endif

ifeq ($(CURLV), yes)
    FPPFLAGS += -DCURLV
endif

ifeq ($(USE_STRAIN_TENSOR), yes)
    FPPFLAGS += -DUSE_STRAIN_TENSOR
endif

ifeq ($(IND_TIMESTEPS), yes)
    FPPFLAGS += -DIND_TIMESTEPS
endif

ifeq ($(STS_TIMESTEPS), yes)
    FPPFLAGS += -DSTS_TIMESTEPS
endif

ifeq ($(CMACIONIZE), yes)
    FPPFLAGS += -DCMACIONIZE
endif

ifeq ($(DEBUG), yes)
    FFLAGS += -DDEBUG
endif

ifdef SRCTURB
    FPPFLAGS += -DDRIVING
endif

#
# kernel choice
#
ifndef SRCKERNEL
ifdef KERNEL
   SRCKERNEL= kernel_${KERNEL}.f90
else
   SRCKERNEL= kernel_cubic.f90
   KERNEL=cubic
endif
endif

#
# can turn particle injection off
# by setting INJECT_PARTICLES=no
# on command line. Otherwise on
# if injection module selected
#
ifeq ($(INJECT_PARTICLES), no)
   SRCINJECT=
else
ifdef SRCINJECT
    FPPFLAGS += -DINJECT_PARTICLES
endif
endif

ifdef LIGHTCURVE
    FPPFLAGS += -DLIGHTCURVE
endif

# do double precision flag last (append only to FFLAGS)

ZZFFLAGS := ${FFLAGS}
ifeq ($(DOUBLEPRECISION), yes)
    FFLAGS += ${DBLFLAG}
endif

ifeq ($(ANALYSISONLY), yes)
    FPPFLAGS += -DANALYSIS
endif

ifeq ($(LIVE_ANALYSIS), yes)
    FPPFLAGS += -DLIVE_ANALYSIS
    SRCAN = $(ANALYSIS)
else
    SRCAN=
endif

#
# MPI flavour (mostly openmpi these days)
#
ifeq ($(MPI), yes)
    FC= mpif90 `mpif90 --showme:compile`
    CC= mpicc `mpicc --showme:compile`
    LDFLAGS+= `mpif90 --showme:link`
    FPPFLAGS += -DMPI
    USEMPI=yes
endif

ifeq ($(MPI), openmpi)
    FC= openmpif90 `openmpif90 --showme:compile`
    LDFLAGS+= `openmpif90 --showme:link`
    FPPFLAGS += -DMPI
    USEMPI=yes
endif

ifeq ($(MPI), zen)
    FC= mpif90
    LDFLAGS+= -lmpi -lmpiif
    FPPFLAGS += -DMPI
    USEMPI=yes
endif

ifeq ($(MPI), psxe)
    FC= mpiifort
    LDFLAGS+= `mpiifort--showme:link`
    FPPFLAGS += -DMPI
    USEMPI=yes
endif

ifeq ($(MPI), mpifort)
    FC= mpifort
    FPPFLAGS += -DMPI
    USEMPI=yes
endif

ifeq ($(MPI), intel)
    FC= mpif90
    FPPFLAGS += -DMPI
    USEMPI=yes
endif

ifeq ($(USEMPI), yes)
    RUNMPI=$(MPIEXEC)
else
    RUNMPI=
endif

#
# HDF5 libraries (if required)
#
# Requires two directories:
#   - include for Fortran .mod files
#   - lib for the shared library .so files
#
# Often both directories are under one root,
# e.g. HDF5ROOT= /usr/local/opt/hdf5
# In this case just set HDF5ROOT for your machine.
#
# However, sometimes these directories are separate,
# then you must set both HDF5INCLUDE and HDF5LIB.
#
ifeq ($(HDF5), yes)
ifeq (X$(HDF5ROOT), X)
    HDF5ROOT= /usr/local/opt/hdf5
endif
ifeq (X$(HDF5INCLUDE), X)
    HDF5INCLUDE= $(HDF5ROOT)/include
endif
ifeq (X$(HDF5LIB), X)
    HDF5LIB= $(HDF5ROOT)/lib
endif
    FFLAGS+= -I$(HDF5INCLUDE)
    CCFLAGS+= -I$(HDF5INCLUDE)
    LDFLAGS+= -L$(HDF5LIB) -lhdf5 -lhdf5_fortran
    FPPFLAGS+= -DHDF5
endif

IDFLAGS=$(FPPFLAGS)
ifeq ($(DEBUG), yes)
    IDFLAGS += -DDEBUG
endif
#
# select domain decomposition type
#
DOMAIN= mpi_domain.F90
OBJDIR=obj

# define the implicit rule to make a .o file from a .f90 file

.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .c .f

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) ${FPP_PREFIX} $(FPPFLAGS) $< -o $@

%.o : %.c
	$(CC) -c $(CCFLAGS) $< -o $@

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

# these are the sources common to all compilations
ifeq (X$(SRCPOTS), X)
SRCPOTS= extern_corotate.f90 \
         extern_binary.f90 \
         extern_spiral.f90 \
         extern_lensethirring.f90 \
         extern_gnewton.F90 \
         lumin_nsdisc.F90 extern_prdrag.F90 \
         extern_Bfield.f90 \
         extern_neutronstar.f90 \
         extern_staticsine.f90 \
         extern_gwinspiral.f90 \
         externalforces.F90
endif
ifeq (X$(SRCPOT), X)
SRCPOT=${SRCPOTS}
endif

SRCTESTS=${TEST_FASTMATH} test_kernel.f90 test_dust.F90 test_growth.F90 test_nonidealmhd.F90 test_gravity.F90 \
         test_derivs.F90 test_cooling.f90 test_eos.f90 test_externf.f90 test_rwdump.f90 \
         test_step.F90 test_indtstep.F90 test_setdisc.F90 \
         test_link.F90 test_kdtree.F90 test_ptmass.F90 test_luminosity.F90\
         test_gnewton.F90 test_corotate.f90 test_geometry.f90 \
         test_sedov.F90

ifeq (X$(SRCTEST), X)
SRCTEST=${SRCTESTS}
endif

SRCCHEM= coolfunc.f90 fs_data.f90 mol_data.f90 utils_spline.f90 h2cooling.f90 h2chem.f90 cooling.f90

SRCMESA= eos_mesa_microphysics.F90 eos_mesa.f90
SRCEOS = ${SRCMESA} eos_helmholtz.f90 eos.F90

ifeq ($(HDF5), yes)
SRCREADWRITE_DUMPS= utils_hdf5.F90 utils_dumpfiles_hdf5.F90 readwrite_dumps_hdf5.F90
else
SRCREADWRITE_DUMPS= readwrite_dumps.F90
endif

SOURCES= physcon.f90 ${CONFIG} ${SRCKERNEL} io.F90 units.f90 boundary.f90 \
         mpi_utils.F90 dtype_kdtree.F90 utils_omp.F90 utils_cpuinfo.f90 \
         utils_allocate.f90 \
         utils_mathfunc.f90 part.F90 ${DOMAIN} utils_timing.f90 mpi_balance.F90 \
         commons.f90 timestep.f90 utils_dumpfiles.f90 utils_indtimesteps.F90 utils_infiles.f90 \
         utils_sort.f90 utils_supertimestep.F90 utils_tables.f90 \
         utils_sphNG.f90 utils_vectors.f90 utils_datafiles.f90 datafiles.f90 \
         gitinfo.f90 ${SRCFASTMATH} random.f90 checkoptions.F90 ${SRCEOS} \
         set_binary.f90 set_flyby.f90 viscosity.f90 options.f90 centreofmass.f90 ${SRCPOT} damping.f90 \
         set_disc.F90 partinject.F90 utils_filenames.f90 utils_summary.F90 ${SRCCHEM} \
         directsum.f90 prompting.f90 ${SRCDUST} set_dust.F90 set_dust_options.f90 \
         mpi_dens.F90 mpi_force.F90 stack.F90 mpi_derivs.F90 kdtree.F90 linklist_kdtree.F90 ${SRCTURB} \
         ${SRCNIMHD} ${SRCPHOTO} ${SRCINJECT} memory.F90 ${SRCREADWRITE_DUMPS}  \
         quitdump.f90 ptmass.F90 \
         readwrite_infile.F90 dens.F90 force.F90 deriv.F90 energies.F90 sort_particles.F90 \
         evwrite.F90 step_leapfrog.F90 writeheader.F90 ${SRCAN} step_supertimestep.F90 mf_write.f90 evolve.F90 \
         geometry.f90 stretchmap.f90 density_profiles.f90 set_unifdis.f90 set_slab.f90 \
         set_sphere.f90 set_vfield.f90 ${SETUPFILE} checksetup.F90 utils_testsuite.f90 \
         ${SRCTEST} \
         testsuite.F90 initial.F90 leastsquares.f90 solvelinearsystem.f90

OBJECTS1 = $(SOURCES:.f90=.o)
OBJECTS = $(OBJECTS1:.F90=.o)

.PHONY: phantom
phantom: checksystem checkparams $(OBJECTS) phantom.o
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJECTS) phantom.o $(LDFLAGS)
ifeq ($(UNAME), Darwin)
	dsymutil $(BINDIR)/$@
endif

	@sh ../scripts/phantom_version_gen.sh "$(IDFLAGS)"
	@echo ""
	@echo "The Phantom is here (in $(BINDIR)/phantom)"
	@echo ""

#----------------------------------------------------
# generic target for compiling ALL phantom utilities
# this is used in the nightly build checks
#
utils: phantomsetup phantomanalysis \
       multirun phantom_moddump \
       phantom2divv phantom2divb \
       diffdumps showheader ev2mdot phantomevcompare acc2ang \
       phantom2sphNG phantom2gadget testbinary \
       sfutils phantom2pdf phantom2pdf-amr \
       phantom2struct phantom2power libphantom

cleanutils: cleansetup cleananalysis \
            cleanmultirun cleantestbinary cleanmoddump \
            cleanphantom2divv cleanphantom2divb \
            cleandiffdumps cleanev2mdot cleanacc2ang \
            cleanp2s cleanphantom2gadget \
            cleansfutils cleanp2p cleanphantom2pdf-amr \
            cleanphantom2struct cleanphantom2power cleanphantomevcompare cleanlibphantom

#--------------------------------------------------------------
# edit target opens current setup module in the default editor
#
edit: checksetup
	$(EDITOR) ../src/setup/$(SETUPFILE)

#----------------------------------------------------
# these are the sources for anything which uses the readwrite_dumps module
#
SRCDUMP= physcon.f90 ${CONFIG} ${SRCKERNEL} io.F90 units.f90 boundary.f90 mpi_utils.F90 \
         utils_infiles.f90 dtype_kdtree.f90 utils_allocate.f90 part.F90 ${DOMAIN} kdtree.F90 linklist_kdtree.F90 \
         utils_dumpfiles.f90 utils_vectors.f90 utils_mathfunc.f90 \
         utils_datafiles.f90 utils_filenames.f90 datafiles.f90 gitinfo.f90 \
         centreofmass.f90 \
         ${SRCEOS} ${SRCPOT} ${SRCPHOTO} \
         memory.F90 \
         utils_sphNG.f90 \
         commons.f90 timestep.f90 ${SRCFASTMATH} checkoptions.F90 \
         viscosity.f90 options.f90 prompting.f90 ${SRCDUST} \
         ${SRCREADWRITE_DUMPS}
OBJDUMP1= $(SRCDUMP:.f90=.o)
OBJDUMP= $(OBJDUMP1:.F90=.o)

# make first file required for compiling utilities depend on math flags
# to ensure that this is always up to date before compiling anything else.
physcon.o: .make_mathflags

#----------------------------------------------------
# these are the sources for phantom setup utility
#
SRCSETUP= utils_omp.F90 utils_sort.f90 utils_timing.f90 utils_summary.F90 \
          utils_tables.f90 random.f90 mpi_balance.F90 set_dust.F90 set_dust_options.f90 set_binary.f90 set_flyby.f90 \
          utils_indtimesteps.F90 partinject.F90 stack.F90 mpi_dens.F90 mpi_force.F90 mpi_derivs.F90 \
          ${SRCTURB} ${SRCNIMHD} ${SRCCHEM} \
          ptmass.F90 energies.F90 \
          geometry.f90 stretchmap.f90 density_profiles.f90 \
          set_unifdis.f90 set_slab.f90 set_sphere.f90 set_disc.F90 \
          set_vfield.f90 sort_particles.F90 ${SRCINJECT} \
          ${SETUPFILE} checksetup.F90 \
          set_Bfield.f90 damping.f90 readwrite_infile.f90

OBJSETUP1= $(SRCSETUP:.f90=.o)
OBJSETUP= $(OBJDUMP) $(OBJSETUP1:.F90=.o) phantomsetup.o

.PHONY: phantomsetup
phantomsetup: setup

setup: checksystem checkparams $(OBJSETUP)
	$(FC) $(FFLAGS) -o $(BINDIR)/phantomsetup $(OBJSETUP) $(LDFLAGS) $(LIBS)
	@echo ""
	@echo "Phantom setup built"
	@echo ""

cleansetup:
	rm -f $(BINDIR)/phantomsetup

config.o: phantom-version.h

phantom-version.h:
	@echo "creating $@"
	@echo "#define PHANTOM_VERSION_MAJOR $(PHANTOM_VERSION_MAJOR)" > $@
	@echo "#define PHANTOM_VERSION_MINOR $(PHANTOM_VERSION_MINOR)" >> $@
	@echo "#define PHANTOM_VERSION_MICRO $(PHANTOM_VERSION_MICRO)" >> $@
	@echo "#define PHANTOM_VERSION_STRING \"$(VERSION)\"" >> $@

#----------------------------------------------------
# run checkhdf5 for HDF5=yes
#
utils_outputhdf5.o: checkhdf5

#----------------------------------------------------
# these are the sources for the phantom2grid utility
#
ifdef HDF5
OBJP2G= $(OBJDUMP) hdf5utils.o write_grid_hdf5.o interpolate3D.o phantom2grid.o
else
OBJP2G= $(OBJDUMP) interpolate3D.o phantom2grid.o
endif

write_grid_hdf5.o: checkhdf5

.PHONY: phantom2grid
phantom2grid:
	${MAKE} ANALYSISONLY=yes phantom2gridfake

phantom2gridfake: checksystem checkparams $(OBJP2G)
	$(FC) $(FFLAGS) -o $(BINDIR)/phantom2grid $(OBJP2G) $(LDFLAGS)
	@echo ""
	@echo "Phantom2grid: we are here to help you"
	@echo ""

cleanp2g:
	rm -f $(BINDIR)/phantom2grid

#------------------------------------------------------------
# these are the sources for the phantom2pdf utility
# to compute Probability Density Functions from Phantom data
#
OBJP2PDF= $(OBJDUMP) asciiutils.o pdfs.o interpolate3D.o rhomach.o phantom2pdf.o

.PHONY: phantom2pdf
phantom2pdf:
	${MAKE} ANALYSISONLY=yes phantom2pdffake

phantom2pdffake: checksystem checkparams $(OBJP2PDF)
	$(FC) $(FFLAGS) -o $(BINDIR)/phantom2pdf $(OBJP2PDF) $(LDFLAGS)
	@echo ""
	@echo "Phantom2pdf: we are Probably Dramatically Fun"
	@echo ""

cleanp2p:
	rm -f $(BINDIR)/phantom2pdf

pdfs.o: checksplash $(SPLASH_DIR)/src/pdfs.f90
	$(FC) $(FFLAGS) -o $@ -c $(SPLASH_DIR)/src/pdfs.f90

# In case you need the old pdfs.f90 module located in phantom/src/utils/
# rather than the on located in splash. (e.g. analysis_MWpdf.f90 requires the
# phantom version)
phantom_pdfs.o: ../src/utils/pdfs.f90
	$(FC) $(FFLAGS) -o $@ -c $<

asciiutils.o: checksplash $(SPLASH_DIR)/src/asciiutils.f90
	$(FC) $(FFLAGS) -o $@ -c $(SPLASH_DIR)/src/asciiutils.f90

# these are the sources for the grid2pdf utility

ifdef HDF5
OBJG2PDF= io.o utils_filenames.o asciiutils.o \
          hdf5utils.o read_grid_hdf5.o write_grid_hdf5.o io_grid.o pdfs.o rhomach.o grid2pdf.o
else
OBJG2PDF= io.o utils_filenames.o asciiutils.o \
          io_grid.o pdfs.o rhomach.o grid2pdf.o
endif

.PHONY: grid2pdf
grid2pdf: checksys checkparams $(OBJG2PDF)
	@echo "objects are $(OBJG2PDF)"
	$(FC) $(FFLAGS) -o $(BINDIR)/grid2pdf $(OBJG2PDF) $(LDFLAGS)
	@echo ""
	@echo "Grid2pdf: we are Possibly Dangerously Fanatical"
	@echo ""

cleang2p:
	rm -f $(BINDIR)/grid2pdf

#------------------------------------------------------
# Probability Distribution Functions via adaptive mesh
#
.PHONY: phantom2pdf-amr
phantom2pdf-amr:
	${MAKE} phantomanalysis ANALYSIS="adaptivemesh.f90 interpolate3D_amr.F90 asciiutils.f90 pdfs.f90 analysis_pdfs.f90"\
        ANALYSISBIN=$@ ANALYSISONLY=yes

cleanphantom2pdf-amr:
	rm -f $(BINDIR)/phantom2struct

analysis_pdfs.o: interpolate3D_amr.o adaptivemesh.o
interpolate3D_amr.o: adaptivemesh.o

#----------------------------------------------------
# these are the sources for the phantom2power utility
#
SRCP2P= utils_omp.F90 fft3d.f90 fftpack.f \
        power3d.f90 icosahedron.f90 powerspec.f90 \
        interpolate3D.F90 asciiutils.f90 pdfs.f90 phantom2power.f90
OBJP2P1= $(SRCP2P:.f90=.o)
OBJP2P2= $(OBJP2P1:.F90=.o)
OBJP2P= $(OBJDUMP) $(OBJP2P2:.f=.o)

.PHONY: phantom2power
phantom2power:
	${MAKE} ANALYSISONLY=yes phantom2powerfake

phantom2powerfake: checksystem checkparams $(OBJP2P)
	$(FC) $(FFLAGS) -o $(BINDIR)/phantom2power $(OBJP2P)
	@echo ""
	@echo "Phantom2power: we have the power to help you"
	@echo ""

cleanphantom2power:
	rm -f $(BINDIR)/phantom2power

#----------------------------------------------------
# these are the sources for the phantom_moddump utility
#
OBJMOD1 = utils_omp.F90 utils_summary.f90 utils_indtimesteps.F90 \
          utils_tables.f90 checksetup.f90 set_Bfield.f90 \
          partinject.F90 random.f90 set_disc.F90 set_dust.F90 set_binary.f90 ${SRCINJECT} \
          ${SRCTURB} ${SRCNIMHD} ${SRCCHEM} \
          density_profiles.f90 ptmass.F90 damping.f90 readwrite_infile.f90 ${MODFILE:.f90=.o}
OBJMOD2 = ${OBJMOD1:.F90=.o}
OBJMOD = ${OBJMOD2:.f90=.o}
OBJDA= ${OBJDUMP} ${OBJMOD} phantom_moddump.o

phantom_moddump: checksystem checkparams $(OBJDA)
	@echo ""
	@echo "phantom_moddump: we are here to help you"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/phantommoddump $(OBJDA) $(LDFLAGS)

moddump: phantom_moddump

cleanmoddump:
	rm -f $(BINDIR)/phantommoddump

#----------------------------------------------------
# these are the sources for the phantomanalysis utility
#
OBJAN1= ${ANALYSIS:.f90=.o}
OBJAN2= ${OBJAN1:.F90=.o}
OBJAN= ${OBJAN2:.f=.o}
OBJA= utils_sort.o utils_tables.o leastsquares.o solvelinearsystem.o \
      ${OBJDUMP} utils_disc.o set_dust.o set_binary.o ${OBJAN}

ifndef ANALYSISBIN
ANALYSISBIN=phantomanalysis
endif

.PHONY: phantomanalysis
phantomanalysis: checksystem checkparams $(OBJA) phantomanalysis.o
	@echo ""
	@echo "phantomanalysis: we live to serve you"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$(ANALYSISBIN) $(OBJA) phantomanalysis.o $(LDFLAGS)

analysis: phantomanalysis


cleananalysis:
	rm -f $(BINDIR)/phantomanalysis

.PHONY: libphantom
SRCLIB=icosahedron.f90 libphantom-evolve.F90  libphantom-splash.f90  libphantom.F90
OBJLIB1=${SRCLIB:.f90=.o}
OBJLIB=${OBJLIB1:.F90=.o}
libphantom: checksystem checkparams
	${MAKE} phantom ${OBJLIB} SETUP=${SETUP} FFLAGS="${FFLAGS} -fPIC"
	$(FC) -shared -fPIC $(FFLAGS) $(FPPFLAGS) $(DBLFLAG) ${OBJLIB} ${OBJECTS} $(LDFLAGS) -o $(BINDIR)/libphantom.so

cleanlibphantom:
	rm -f $(BINDIR)/libphantom.so

.PHONY: pyanalysis
pyanalysis: libphantom

#------------------------------------------------------
# Various utilities for computing structure functions
# and manipulating the resulting output
#
.PHONY: phantom2struct
phantom2struct:
	${MAKE} phantomanalysis ANALYSIS="utils_timing.f90 io_structurefn.f90 random.f90 struct_part.f90 analysis_structurefn.f90"\
        ANALYSISBIN=$@ ANALYSISONLY=yes

cleanphantom2struct:
	rm -f $(BINDIR)/phantom2struct

# conversion between structure function file formats
.PHONY: struct2struct
STRUCT2STRUCTOBJ= utils_filenames.o io_structurefn.o struct2struct.o
struct2struct: checksys checkparams ${STRUCT2STRUCTOBJ}
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ${STRUCT2STRUCTOBJ}

cleanstruct2struct:
	rm -f $(BINDIR)/struct2struct

# time average of structure function files
.PHONY: time_average_struct time_average_sf
TIMEAVERAGESFOBJ=utils_filenames.o io_structurefn.o time_average_sf.o
time_average_sf: time_average_struct
time_average_struct: checksys checkparams ${TIMEAVERAGESFOBJ}
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ${TIMEAVERAGESFOBJ}

cleantime_average_struct:
	rm -f $(BINDIR)/time_average_struct

# structure function slope calculation
.PHONY: get_struct_slope get_struct_slope
GETSLOPESFOBJ=utils_filenames.o io_structurefn.o leastsquares.o get_struct_slope.o
get_slope_sf: get_struct_slope
get_struct_slope: checksys checkparams ${GETSLOPESFOBJ}
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ${GETSLOPESFOBJ}

cleanget_struct_slope:
	rm -f $(BINDIR)/time_average_struct

sfutils: structutils
structutils: time_average_sf struct2struct get_slope_sf

cleansfutils: cleanstructutils
cleanstructutils: cleantime_average_struct cleanstruct2struct cleanget_struct_slope

#----------------------------------------------------
# utility to calculate divv from a dump file
# compile using all phantom files
#
phantom2divv: checksys checkparams $(OBJECTS) phantom2divv.o
	@echo ""
	@echo "phantom2divv: divergence is beautiful"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJECTS) phantom2divv.o

cleanphantom2divv:
	rm -f $(BINDIR)/phantom2divv

#----------------------------------------------------
# utility to calculate divB & curlB from a dump file
# compile using all phantom files
#
phantom2divb: checksys checkparams $(OBJECTS) phantom2divb.o
	@echo ""
	@echo "phantom2divb: divergence should be eradicated"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJECTS) phantom2divb.o

cleanphantom2divb:
	rm -f $(BINDIR)/phantom2divb

#----------------------------------------------------
# these are the sources for the diffdumps utility
#
diffdumps: checksys checkparams $(OBJDUMP) utils_testsuite.o diffdumps.o
	@echo ""
	@echo "diffdumps: we welcome you"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) utils_testsuite.o diffdumps.o

cleandiffdumps:
	rm -f $(BINDIR)/phantom2divb

#----------------------------------------------------
# these are the sources for the phantom2sphNG utility
#
phantom2sphNG: checksystem checkparams $(OBJDUMP) phantom2sphNG.o
	@echo ""
	@echo "phantom2sphNG: now why would you want to do that?"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) phantom2sphNG.o

p2s: phantom2sphNG

cleanp2s:
	rm -f $(BINDIR)/phantom2sphNG

#----------------------------------------------------
# these are the sources for the phantom2sphNG utility
#
phantom2gadget: checksystem checkparams $(OBJDUMP) phantom2gadget.o
	@echo ""
	@echo "phantom2gadget: now why would you want to do that?"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) phantom2gadget.o

p2g: phantom2gadget

cleanphantom2gadget:
	rm -f $(BINDIR)/phantom2gadget

#----------------------------------------------------
# these are the sources for the phantom2mcfost utility
#
.PHONY: phantom2mcfost
phantom2mcfost: checkmcfost
	${MAKE} phantomanalysis ANALYSIS="analysis_mcfost.f90"\
        ANALYSISBIN=$@ ANALYSISONLY=yes LDFLAGS="-L$(MCFOST_DIR)/src -lmcfost $(LIBCXX)"

analysis_mcfost.o: analysis_mcfost.f90
	$(FC) -c $(FFLAGS) -I$(MCFOST_DIR)/src $< -o $@

analysis_mcfost.o: checkmcfost

cleanphantom2mcfost:
	rm -f $(BINDIR)/phantom2mcfost

#----------------------------------------------------
# these are the sources for the phantom2hdf5 utility
#
ifeq ($(PHANTOM2HDF5), yes)
    FPPFLAGS += -DPHANTOM2HDF5
endif
checkhdf5flags:
   ifneq ($(HDF5), yes)
	@echo "-----------------------------"
	@echo "Need to compile with HDF5=yes"
	@echo "-----------------------------"
	${MAKE} err
   endif
   ifneq ($(PHANTOM2HDF5), yes)
	@echo "-------------------------------------"
	@echo "Need to compile with PHANTOM2HDF5=yes"
	@echo "-------------------------------------"
	${MAKE} err
   endif

.PHONY: phantom2hdf5
phantom2hdf5: checksystem checkparams checkhdf5flags $(OBJDUMP) readwrite_dumps.o phantom2hdf5.o
	@echo ""
	@echo "phantom2hdf5: welcome to the future"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUMP) readwrite_dumps.o phantom2hdf5.o $(LDFLAGS)

p2h: phantom2hdf5

cleanphantom2hdf5:
	rm -f $(BINDIR)/phantom2hdf5

#----------------------------------------------------
# utility to rewrite .ev files using a common header
#
SRCEV=utils_infiles.f90 utils_evfiles.f90 prompting.f90 phantomevcompare.f90
OBJEVC1 = ${SRCEV:.f90=.o}
OBJEVC = ${OBJEVC1:.F90=.o}

.PHONY: phantomevcompare
phantomevcompare: $(OBJEVC)
	@echo ""
	@echo "phantomevcompare: let the graphing begin!"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJEVC)

cleanphantomevcompare:
	rm -f $(BINDIR)/phantomevcompare

#----------------------------------------------------
# these are the sources for the multirun utility
#
SRCMULT = physcon.f90 ${CONFIG} ${SRCKERNEL} io.F90 mpi_utils.F90 ${SRCFASTMATH} \
          units.f90 boundary.f90 utils_allocate.f90 part.F90 timestep.f90 commons.f90 \
          utils_filenames.f90 utils_mathfunc.f90 utils_vectors.f90 utils_omp.F90 utils_datafiles.f90 datafiles.f90 \
          viscosity.f90 options.f90 damping.f90 ${SRCEOS} \
          utils_infiles.f90 utils_dumpfiles.f90 utils_summary.f90 centreofmass.f90 \
          ${SRCCHEM} ${DOMAIN} ${SRCPOT} ptmass.F90 ${LINKLIST} ${SRCTURB} \
          prompting.f90 ${SRCDUST} ${SRCNIMHD} readwrite_infile.f90 ${MULTIRUNFILE}
OBJM1 = ${SRCMULT:.f90=.o}
OBJMULT = ${OBJM1:.F90=.o}

multirun: checksystem checkparams $(OBJMULT)
	@echo ""
	@echo "multirun: your hope is our desire"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJMULT)

cleanmultirun:
	rm -f $(BINDIR)/multirun

#----------------------------------------------------
# utility to plot orbits based on orbital elements (a,e,i,o,w,f)
#
SRCBIN = prompting.f90 utils_datafiles.f90 datafiles.f90 ${CONFIG} physcon.f90 io.F90 \
         mpi_utils.F90 utils_allocate.f90 part.F90 mpi_domain.F90 set_binary.f90 test_binary.f90 testbinary.f90
OBJBIN1 = ${SRCBIN:.f90=.o}
OBJBIN = ${OBJBIN1:.F90=.o}

.PHONY: testbinary

testbin: testbinary

testbinary: checksys checkparams $(OBJBIN)
	@echo ""
	@echo "test_binary: may your orbits orbit"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/testbinary $(OBJBIN)

cleantestbinary:
	rm -f $(BINDIR)/testbinary

#----------------------------------------------------
# check for anything that depends on HDF5
#
checkhdf5:
   ifeq (X${HDF5ROOT}, X)
	@echo; echo "ERROR: HDF5ROOT should be set before compiling with HDF5 utilities"; echo; ${MAKE} err;
   else
	@if [ -d $$HDF5ROOT ]; then echo; echo "HDF5ROOT=$$HDF5ROOT"; echo; else echo; echo "ERROR: Directory given by HDF5ROOT=$$HDF5ROOT does not exist"; echo; ${MAKE} err; fi;
   endif

#----------------------------------------------------
# these are the sources for the plot_kernel utility
#

OBJPLOTK= physcon.o ${SRCKERNEL:.f90=.o} giza-fortran.o plot_kernel.o

plotkernel: checksys checkparams checksplash $(OBJPLOTK)
	@echo ""
	@echo "plot_kernel: may your kernels be normalised"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJPLOTK) $(LDFLAGS) -L$(SPLASH_DIR)/giza/lib -lgiza

plot_kernel.o: ${SRCKERNEL:.f90=.o}
#giza-fortran.o: ${SPLASH_DIR}/giza/src/$@
#	$(FC) $(FFLAGS) -o $@ -c ${SPLASH_DIR}/giza/interface/giza-fortran.F90

cleanplotkernel:
	rm -f $(BINDIR)/plotkernel

#----------------------------------------------------
# these are the sources for the showheader utility
#
SRCSHOWHEADER= utils_dumpfiles.f90 showheader.f90
OBJSHOWHEADER= $(SRCSHOWHEADER:.f90=.o)
showheader: checksys $(OBJSHOWHEADER)
	@echo ""
	@echo "showheader: show me the header!"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJSHOWHEADER)

#----------------------------------------------------
# these are the sources for the evol_dustywaves utility
#
SRCDUSTEVOL= cubicsolve.f90 dustywaves.f90 evol_dustywaves.f90
OBJDUSTEVOL= $(SRCDUSTEVOL:.f90=.o)

evol_dustywaves: checksys $(OBJDUSTEVOL)
	@echo ""
	@echo "dusty wave .ev solutions^TM: All the energy you need."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJDUSTEVOL)

#----------------------------------------------------
# these are the sources for the ev2mdot utility
#
.PHONY: ev2mdot
ev2mdot: checksys utils_filenames.o utils_infiles.o utils_evfiles.o ev2mdot.o
	@echo ""
	@echo "ev2mdot: Accretion rates R us."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ ev2mdot.o utils_filenames.o utils_evfiles.o utils_infiles.o

cleanev2mdot:
	rm -f $(BINDIR)/ev2mdot

#----------------------------------------------------
# these are the sources for the acc2ang utility
#
.PHONY: acc2ang
acc2ang: checksys acc2ang.o
	@echo ""
	@echo "acc2ang: Accreted ang. mom. R us."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ acc2ang.o

cleanacc2ang:
	rm -f $(BINDIR)/acc2ang

#---------------------------
# sources for the mass_flow utility
#
OBJMF1 = ${ANALYSIS:.f90=.o}
OBJMF2 = ${OBJMF1:.F90=.o}
OBJMF = ${OBJMF2:.f=.o}
OBJM= utils_sort.o leastsquares.o solvelinearsystem.o ${OBJDUMP} ${OBJMF} set_binary.o mf_write.o

.PHONY: mflow
mflow: checksys $(OBJM)  mflow.o ev2mdot lombperiod
	@echo ""
	@echo "mflow: mass flow R us."
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@  $(OBJM) mflow.o

.PHONY:lombperiod
lombperiod: powerspectrums.o lombperiod.o
	$(FC) $(FFLAGS) -o $(BINDIR)/$@  lombperiod.o powerspectrums.o

#----------------------------------------------------
# these are the sources for the combinedustdumps utility
#
OBJCDD= ${OBJECTS} combinedustdumps.o

combinedustdumps: checksys checkparams $(OBJCDD)
	@echo ""
	@echo "combinedustdumps: many grains make light work"
	@echo ""
	$(FC) $(FFLAGS) -o $(BINDIR)/$@ $(OBJCDD) $(LDFLAGS)

cleancombinedustdumps:
	rm -f $(BINDIR)/combinedustdumps



#----------------------------------------------------
# target to write appropriate queue submission script
#
ifndef QSYS
   QSYS=pbs
endif
ifndef WALLTIME
   WALLTIME='1000:00:00'
endif
ifndef MAXMEM
   MAXMEM='16G'
endif
ifeq ($(OPENMP),yes)
 ifndef NOMP
   ifdef OMP_NUM_THREADS
      NOMP=$(OMP_NUM_THREADS)
   else
      NOMP=2
   endif
 endif
 ifndef OMP_SCHEDULE
    OMP_SCHEDULE=dynamic
 endif
 ifndef QPE
    QPE=omp
 endif
 ifndef NPAR
    NPAR=$(NOMP)
 endif
endif
ifeq ($(USEMPI),yes)
 ifndef NMPI
    NMPI=8
 endif
 ifndef QPE
    QPE=mpi
 endif
 ifndef NPAR
    NPAR=$(NMPI)
 endif
 ifndef MPIEXEC
    MPIEXEC=mpiexec -np ${NMPI}
 endif
else
 ifndef NMPI
    NMPI=1
 endif
endif
ifndef OUTFILE
  ifeq ($(QSYS),sge)
    OUTFILE=$(INFILE)'.sgeout'
  else
    ifeq ($(QSYS),pbs)
       OUTFILE=$(INFILE)'.pbsout'
    else
       OUTFILE=$(INFILE)'.qout'
    endif
  endif
endif
ifndef MAILTO
   MAILTO=`git config --get user.email`
endif
GETLOG='`grep logfile "$(INFILE)" | sed "s/logfile =//g" | sed "s/\\!.*//g" | sed "s/\s//g"`'

ifndef CMD
CMD='./phantom $(INFILE) >& $$outfile'
endif

.PHONY: qscript

qscript:
    ifneq ($(KNOWN_SYSTEM), yes)
	@echo "Error: qscript needs known SYSTEM variable set"
	@${MAKE} err;
    endif
    ifndef INFILE
	@echo
	@echo "Usage: make qscript INFILE=infile"
	@echo
	@${MAKE} err;
    endif
    # set default values for variables not set
    ifeq ($(QSHELL),tcsh)
	@echo '#!/bin/tcsh'
    else
	@echo '#!/bin/bash'
    endif
    ifeq ($(QSYS),sge)
	@echo '## Sun Grid Engine Script, created by "make qscript" '`date`
        ifeq ($(QSHELL),tcsh)
	    @echo '#$$ -S /bin/tcsh'
        else
	    @echo '#$$ -S /bin/bash'
        endif
	@echo '#$$ -cwd'
	@echo '#$$ -N '`../scripts/randomword.pl`
	@echo '#$$ -o '$(OUTFILE)' -j y'
	@echo '#$$ -l h_rt='$(WALLTIME)
	@echo '#$$ -l h_vmem='$(MAXMEM)
        ifdef MAILTO
	   @echo '#$$ -m ae'
	   @echo '#$$ -M '$(MAILTO)
        endif
        ifdef QPE
	   @echo '#$$ -pe '$(QPE) $(NPAR)
        endif
        ifdef QEXTRA
	   @echo '#$$ '$(QEXTRA)
        endif
	@echo
	@echo 'echo "SGE: HOSTS   = "`cat $$PE_HOSTFILE`'
	@echo 'echo "SGE: NHOSTS  = $$NHOSTS"'
	@echo 'echo "SGE: NSLOTS  = $$NSLOTS"'
	@echo 'echo "SGE: NQUEUES = $$NQUEUES"'
    else ifeq ($(QSYS),pbs)
	@echo '## PBS Job Submission Script, created by "make qscript" '`date`
        ifdef QNODES
	   @echo '#PBS -l '$(QNODES)
        else
           ifeq ($(SYSTEM),zen)
	      @echo '#PBS -l nodes='$(NMPI)':ppn=8:StandardMem'
           else
	      @echo '#PBS -l nodes='$(NMPI)':ppn='$(NOMP)
           endif
        endif
        ifdef JOBNAME
	   @echo '#PBS -N '$(JOBNAME)
        else
	   @echo '#PBS -N '`../scripts/randomword.pl`
        endif
        ifdef QNAME
	   @echo '#PBS -q '$(QNAME)
        endif
        ifdef QPROJECT
	   @echo '#PBS -P '$(QPROJECT)
        endif
	@echo '#PBS -o '$(OUTFILE)
	@echo '#PBS -j oe'
        ifdef MAILTO
	   @echo '#PBS -m e'
	   @echo '#PBS -M '$(MAILTO)
        endif
	@echo '#PBS -l walltime='$(WALLTIME)
	@echo '#PBS -l mem='$(MAXMEM)
        ifdef QEXTRA
	   @echo '#PBS '$(QEXTRA)
        endif
	@echo '## phantom jobs can be restarted:'
	@echo '#PBS -r y'
        ifeq ($(PBSRESUBMIT),yes)
             ifeq ($(QSHELL),tcsh)
	          $(error error: resubmittable scripts require bash, cannot use QSHELL=tcsh);
             endif
	     @echo '#PBS -v NJOBS,NJOB'
	     @echo
	     @echo '#------------------------------------------------------------------------------'
	     @echo '# this is a self-resubmitting PBS script'
	     @echo '# use qsub -v NJOBS=10 <scriptname> to submit'
	     @echo '# with an appropriate value for NJOBS'
	     @echo '#'
	     @echo '# These variables are assumed to be set:'
	     @echo '#   NJOBS is the total number of jobs in a sequence of jobs (defaults to 1)'
	     @echo '#   NJOB is the number of the previous job in the sequence (defaults to 0)'
	     @echo '#------------------------------------------------------------------------------'
	     @echo 'if [ X$$NJOBS == X ]; then'
	     @echo '    echo "NJOBS (total number of jobs in sequence) is not set - defaulting to 1"'
	     @echo '    export NJOBS=1'
	     @echo 'fi'
	     @echo 'if [ X$$NJOB == X ]; then'
	     @echo '    echo "NJOB (previous job number in sequence) is not set - defaulting to 0"'
	     @echo '    export NJOB=0'
	     @echo 'fi'
	     @echo '#'
	     @echo '# Quick termination of job sequence - look for a file called STOP_SEQUENCE'
	     @echo '#'
	     @echo 'if [ -f $$PBS_O_WORKDIR/STOP_SEQUENCE ]; then'
	     @echo '    echo  "Terminating sequence after $$NJOB jobs"'
	     @echo '    exit 0'
	     @echo 'fi'
	     @echo '#'
	     @echo '# Increment the counter to get current job number'
	     @echo '#'
	     @echo 'NJOB=$$(($$NJOB+1))'
	     @echo '#'
	     @echo '# Are we in an incomplete job sequence - more jobs to run ?'
	     @echo '#'
	     @echo 'if [ $$NJOB -lt $$NJOBS ]; then'
	     @echo '    #'
	     @echo '    # Now submit the next job'
	     @echo '    #'
	     @echo '    NEXTJOB=$$(($$NJOB+1))'
	     @echo '    echo "Submitting job number $$NEXTJOB in sequence of $$NJOBS jobs"'
	     @echo '    qsub -z -W depend=afterany:$$PBS_JOBID $$0'
	     @echo 'else'
	     @echo '    echo "Running last job in sequence of $NJOBS jobs"'
	     @echo 'fi'
#	     @echo '#'
#	     @echo '# File manipulation prior to job commencing, eg. clean up previous output files,'
#	     @echo '# check for consistency of checkpoint files, ...'
#	     @echo '#'
#	     @echo 'if [ $$NJOB -gt 1 ]; then'
#	     @echo '   echo " "'
#	     @echo '   # .... USER INSERTION HERE '
#	     @echo 'fi'
	     @echo '#------------------------------------------------------------------------------'
        endif
	@echo
	@echo 'cd $$PBS_O_WORKDIR'
	@echo 'echo "PBS_O_WORKDIR is $$PBS_O_WORKDIR"'
	@echo 'echo "PBS_JOBNAME is $$PBS_JOBNAME"'
	@echo 'env | grep PBS'
	@echo 'cat $$PBS_NODEFILE > nodefile'
    else
        ifdef QNODES
	   @echo '#SBATCH --nodes='$(QNODES)
        else
	   @echo '#SBATCH --nodes='$(NMPI)' --ntasks='$(NOMP)
        endif
	@echo '#SBATCH --cpus-per-task=1'
        ifdef JOBNAME
	   @echo '#SBATCH --job-name='$(JOBNAME)
        else
	   @echo '#SBATCH --job-name='`../scripts/randomword.pl`
        endif
        ifdef QNAME
	   @echo '#SBATCH --queue='$(QNAME)
        endif
        ifdef QPROJECT
	   @echo '#SBATCH --account='$(QPROJECT)
        endif
        ifdef QPARTITION
	   @echo '#SBATCH --partition='$(QPARTITION)
        endif
	@echo '#SBATCH --output='$(OUTFILE)
        ifdef MAILTO
	   @echo '#SBATCH --mail-type=BEGIN'
	   @echo '#SBATCH --mail-type=FAIL'
	   @echo '#SBATCH --mail-type=END'
	   @echo '#SBATCH --mail-user='$(MAILTO)
        endif
	@echo '#SBATCH --time=0-'$(WALLTIME)
	@echo '#SBATCH --mem='$(MAXMEM)
        ifdef QEXTRA
	   @echo '#SBATCH '$(QEXTRA)
        endif
    endif
	@echo 'echo "HOSTNAME = $$HOSTNAME"'
	@echo 'echo "HOSTTYPE = $$HOSTTYPE"'
	@echo 'echo Time is `date`'
	@echo 'echo Directory is `pwd`'
	@echo
    ifeq ($(QSHELL),tcsh)
	@echo 'limit stacksize unlimited'
    else
	@echo 'ulimit -s unlimited'
    endif
    #-- set openMP environment variables
    ifeq ($(OPENMP),yes)
        ifeq ($(QSHELL),tcsh)
	   @echo 'setenv OMP_SCHEDULE "'$(OMP_SCHEDULE)'"'
	   @echo 'setenv OMP_NUM_THREADS '$(NOMP)
	   @echo 'setenv OMP_STACKSIZE 1024m'
        else
	   @echo 'export OMP_SCHEDULE="'$(OMP_SCHEDULE)'"'
	   @echo 'export OMP_NUM_THREADS='$(NOMP)
	   @echo 'export OMP_STACKSIZE=1024m'
        endif
    endif
	@echo
    #-- add lines specific to particular machines
    ifeq ($(SYSTEM),msg)
        ifeq ($(QSHELL),bash)
	   @echo 'source /etc/profile'
	   @echo 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}'
        else
	   @echo 'setenv LD_LIBRARY_PATH '${LD_LIBRARY_PATH}
        endif
	@cat ~/.modules
    endif
	@echo
    #--final line is code execution
	@echo 'echo "starting phantom run..."'
    ifeq ($(QSHELL),tcsh)
	@echo 'setenv outfile '$(GETLOG)
    else
	@echo 'export outfile='$(GETLOG)
    endif
	@echo 'echo "writing output to $$outfile"'
    ifeq ($(USEMPI),yes)
	@echo $(MPIEXEC)' '$(CMD)
    else
	@echo $(CMD)
    endif
    ifeq ($(PBSRESUBMIT),yes)
	@echo
	@echo '#------------------------------------------------------------------------------'
	@echo '# Not expected to reach this point in general but if we do, check that all '
	@echo '# is OK.  If the job command exited with an error, terminate the job'
	@echo '#'
	@echo 'errstat=$$?'
	@echo 'if [ $$errstat -ne 0 ]; then'
	@echo '    # A brief nap so PBS kills us in normal termination. Prefer to '
	@echo '    # be killed by PBS if PBS detected some resource excess'
	@echo '    sleep 5  '
	@echo '    echo "Job number $$NJOB returned an error status $$errstat - stopping job sequence."'
	@echo '    touch $$PBS_O_WORKDIR/STOP_SEQUENCE'
	@echo '    exit $$errstat'
	@echo 'fi'
	@echo '#------------------------------------------------------------------------------'
    endif

#----------------------------------------------------
# unit test for block limits
#
test1: checksystem checkparams $(OBJDUMP) test_blocklimits.o
	$(FC) $(FFLAGS) -o $(BINDIR)/test1 $(OBJDUMP) test_blocklimits.o

#----------------------------------------------------
# run test suite
#
.PHONY: test test2 testcyl testgrav testall
test:
	${MAKE} SETUP=test && $(RUNMPI) $(BINDIR)/phantom test

test2:
	${MAKE} SETUP=test2 && $(RUNMPI) $(BINDIR)/phantom test

testkd:
	${MAKE} SETUP=testkd && $(RUNMPI) $(BINDIR)/phantom test

testcyl:
	${MAKE} SETUP=testcyl && $(RUNMPI) $(BINDIR)/phantom test

testgrav:
	${MAKE} SETUP=testgrav && $(RUNMPI) $(BINDIR)/phantom test gravity

testdust:
	${MAKE} SETUP=testdust && $(RUNMPI) $(BINDIR)/phantom test dust

testgrowth:
	${MAKE} SETUP=testgrowth && $(RUNMPI) $(BINDIR)/phantom test growth

testnimhd:
	${MAKE} SETUP=testnimhd && $(RUNMPI) $(BINDIR)/phantom test nimhd

testall: test test2 testcyl testgrav

#----------------------------------------------------
# this is a utility to test the fast sqrt functions
# to see if they are faster than the native calls
# if so, then the appropriate pre-processor flags
# are added
#
.PHONY: .make_mathflags .make_nofastmath getmathflags checkmath
ifndef FASTSQRT
   FASTSQRT=${shell if [ -e .make_nofastmath ]; then echo no; fi}
endif

ifeq ($(FASTSQRT), no)
   OBJTESTMATH=
   FASTMATH=no
else
   OBJTESTMATH= random.o io.o fastmath.o mpi_utils.o test_fastmath.o getmathflags.o
   FASTMATH=${shell if [ -e .make_mathflags ]; then cat .make_mathflags; fi}
endif

.make_mathflags: checksys $(OBJTESTMATH)
     ifeq ($(FASTSQRT), no)
	@touch .make_mathflags
     else
	@if [ ! -e $@ ]; then \
	    $(FC) $(FFLAGS) -o $(BINDIR)/getmathflags $(OBJTESTMATH) || ${MAKE} fastmathlinkerr; \
	    $(BINDIR)/getmathflags > .make_mathflags; \
	fi
     endif

ifeq ($(FASTMATH), yes)
   SRCFASTMATH=fastmath.o
   TEST_FASTMATH=test_fastmath.F90
   FPPFLAGS+=-DFINVSQRT
else
   SRCFASTMATH=
   TEST_FASTMATH=
endif

fastmath.o: fastmath.f90
	$(FC) $(FFLAGS) -o $@ -c $< || ${MAKE} fastmathlinkerr
test_fastmath.o: test_fastmath.F90
	$(FC) $(FFLAGS) -o $@ -c $< || ${MAKE} fastmathlinkerr
getmathflags.o: getmathflags.f90
	$(FC) $(FFLAGS) -o $@ -c $< || ${MAKE} fastmathlinkerr

fastmathlinkerr:
	@echo "***********************************************************************"
	@echo "*** ERROR linking fastsqrt stuff (requires Fortran->C call)         ***"
	@echo "*** Type make again to ignore this and compile without it           ***"
	@echo "***********************************************************************"
	@touch .make_mathflags
	@touch .make_nofastmath
	${MAKE} err;

#----------------------------------------------------

LASTSYSTEM = ${shell if [ -e .make_lastsystem ]; then cat .make_lastsystem; fi}
LASTSETUP = ${shell if [ -e .make_lastsetup ]; then cat .make_lastsetup; fi}
LASTFPPFLAGS = ${shell if [ -e .make_lastfppflags ]; then cat .make_lastfppflags; fi}
LASTFFLAGS = ${shell if [ -e .make_lastfflags ]; then cat .make_lastfflags; fi}

.PHONY: checksystem checkparams checksplash checksys

checksystem: checksys checksetup

checksys:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling Phantom v$(PHANTOM_VERSION_MAJOR).$(PHANTOM_VERSION_MINOR).$(PHANTOM_VERSION_MICRO) for $(SYSTEM) system..........."
	@echo ""
        ifneq ($(SYSTEM),$(LASTSYSTEM))
	    @echo system changed from ${LASTSYSTEM} to ${SYSTEM}
	    @${MAKE} clean
	    @${MAKE} cleanmathflags
        endif
	@echo $(SYSTEM) > .make_lastsystem
   else
	@echo ""
	@echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=> set the environment variable SYSTEM to one listed "
	@echo "   in build/Makefile and try again"
	@echo ""
	@${MAKE} compilers
	@${MAKE} err;
   endif

checksetup:
   ifeq ($(OBSOLETE_SETUP), yes)
	@echo "make: WARNING: value of SETUP = $(OLDSETUP) is obsolete..."
	@echo "=> setting SETUP = $(SETUP)"
	@echo
   endif
   ifeq ($(KNOWN_SETUP), yes)
	@echo "Using options for "$(SETUP)" setup"
	@echo ""
        ifneq ($(SETUP),$(LASTSETUP))
	    @echo setup changed from ${LASTSETUP} to ${SETUP}
	    @${MAKE} clean
        endif
	@echo $(SETUP) > .make_lastsetup
   else
	@echo "setup '$(SETUP)' not recognised..."
	@echo ""
	@echo "Please set SETUP to one listed in build/Makefile"
	@echo ""
	@echo " e.g.:"
	@echo " make SETUP=sedov"
	@echo " make SETUP=disc"
	@echo " make SETUP=turbdrive"
	@echo ""
	@echo " or:"
	@echo " export SETUP=sedov"
	@echo " make"
	@echo ""
	@echo "You may also wish to consider the following compile-time options:"
	@echo ""
	@echo " DEBUG=yes/no"
	@echo " DOUBLEPRECISION=yes/no"
	@echo " OPENMP=yes/no"
	@echo " ENDIAN=BIG/LITTLE"
	@echo ""
	@${MAKE} err;
   endif

checkparams:
	@echo "Using $(KERNEL) kernel"
   ifeq ($(DEBUG), yes)
	@echo "Debugging flags are ON"
   endif
   ifeq ($(DOUBLEPRECISION), yes)
	@echo "Flags set for DOUBLE PRECISION"
   else
	@echo "Flags set for SINGLE PRECISION"
   endif
   ifeq ($(OPENMP), yes)
	@echo "Compiling in PARALLEL (OpenMP)"
   else
	@echo "Compiling in SERIAL"
   endif
   ifeq ($(ENDIAN), BIG)
	@echo "Flags set for conversion to BIG endian"
   endif
   ifeq ($(ENDIAN), LITTLE)
	@echo "Flags set for conversion to LITTLE endian"
   endif
   ifneq ($(FPPFLAGS),$(LASTFPPFLAGS))
	@echo 'pre-processor flags changed from "'${LASTFPPFLAGS}'" to "'${FPPFLAGS}'"'
	@${MAKE} clean;
	#for x in ../src/*/*.F90; do y=`basename $$x`; rm -f $${y/.F90/.o}; done
   endif
	@echo "Preprocessor flags are "${FPPFLAGS}
	@echo "${FPPFLAGS}" > .make_lastfppflags
   ifneq ($(FFLAGS),$(LASTFFLAGS))
	@echo 'Fortran flags changed from "'${LASTFFLAGS}'" to "'${FFLAGS}'"'
	@${MAKE} clean;
   endif
	@echo "Fortran flags are "${FFLAGS}
	@echo "${FFLAGS}" > .make_lastfflags

checksplash:
   ifneq ("X$(SPLASH_DIR)","X")
	@echo; echo "Compiling SPLASH source files from "$(SPLASH_DIR); echo
   else
	@echo; echo "ERROR: cannot find SPLASH directory needed for some source files - try \"export SPLASH_DIR=${HOME}/splash\""; echo
   endif

checkmcfost:
   ifneq ("X$(MCFOST_DIR)","X")
	@echo; echo "MCFOST directory is "$(MCFOST_DIR); echo;
   else
	@echo; echo "ERROR: cannot find MCFOST directory for linking - set this using MCFOST_DIR"; echo; ${MAKE} err
   endif

giza-fortran.o : $(SPLASH_DIR)/giza/interface/giza-fortran.F90 $(SPLASH_DIR)/giza/lib/libgiza.a
	$(FC) $(FFLAGS) -I$(SPLASH_DIR)/giza/include/ -c $< -o $@

compilers:
	@echo "I suggest one of the following, based on detected Fortran compilers..."; echo;
	@if type -p ifort > /dev/null; then echo "make SYSTEM=ifort"; fi;
	@if type -p pathf90 > /dev/null; then echo "make SYSTEM=pathf90"; fi;
	@if type -p pgf90 > /dev/null; then echo "make SYSTEM=pgf90"; fi;
	@if type -p xlf90_r > /dev/null; then echo "make SYSTEM=ukaff1a [uses xlf90_r]"; fi;
	@if type -p gfortran > /dev/null; then echo "make SYSTEM=gfortran"; fi;
	@if type -p g95 > /dev/null; then echo "make SYSTEM=g95"; fi;
	@echo "(end of possible selections)"; echo;

#----------------------------------------------------
# target to automatically include dependencies in Makefile
# relies on the g95 compiler being present
# (does not have to be used for the main compilation)

depends: clean checksetup
	#@echo '*********************************************************************************'
	#@echo 'First run of Makefile -- creating dependency lines using gfortran, writing to .depends'
	#@echo '*********************************************************************************'
	#@gfortran -M -cpp -c ../src/*/*.*90 > .depends
	#@echo '*************************************************************************'
	#@echo 'If no errors above, then Makefile dependencies were created successfully '
	#@echo ' -- be sure to run "make depends" again if you alter code dependencies'
	#@echo '*************************************************************************'
	#@${MAKE} clean

.depends:
	@if type -p gfortran; then touch .depends; ${MAKE} --quiet SETUP=test depends; else echo "warning: no gfortran so dependencies not calculated"; touch .depends; fi;

include .depends

getdims:
	@echo $(MAXP)

err:
	$(error aborting);

clean:
	rm -f *.o *.mod phantom-version.h

cleanall: clean cleanmathflags
	cd $(BINDIR); rm -f phantom phantomsetup

cleandist: clean cleanall
	rm -f .make_lastsystem .make_lastsetup .make_lastfppflags .depends

cleanmathflags:
	rm -f .make_mathflags bin/getmathflags
