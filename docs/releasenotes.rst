Release notes
=============

v2024.0.0 - 29th Jan 2024
-------------------------

Physics
~~~~~~~
- ability to use numerical relativity backend with phantom (`Magnall et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023PhRvD.108j3534M/abstract>`__; #480)
- further improvements to implicit radiation scheme (thanks to Mike Lau and Ryosuke Hirai; #406,#438,#441,#452,#455,#458,#474)
- further improvements to wind injection and cooling modules (thanks to Lionel Siess, Mats Esseldeurs, Silke Maes and Jolien Malfait; #392,)
- J2 potential due to oblateness implemented for sink particles (#289)
- external potential implemented for geopotential model, to test J2 potential (#289)
- implemented Loren/Bate implicit scheme for drag with dust-as-particles (thanks to Stephane Michoulier, #428,#436)
- dynamic boundary conditions, allowing box with expanding boundaries (thanks to James Wurster; #416)
- bug fix in generalised Farris equation of state (thanks to Nicolas Cuello; #433)

Setup
~~~~~
- major reorganisation of star setup into separate module, can now setup and relax one or more stars in several different setups, allowing one-shot-setup-and-relax for common envelopes, binary stars and tidal disruption events (#405,#407,#413)
- new hierarchical system setup: can now setup an arbitrary number of point masses or stars in hierarchical systems (thanks to Simone Ceppi; #401,#426; see `Ceppi et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.514..906C/abstract>`__)
- relaxation process for stars is restartable, works automatically (#414, #417)
- can setup unbound parabolic and hyperbolic orbits using the standard 6-parameter orbital elements (#443,#448; #302)
- use m1 and m2 in the binary disc setup instead of primary mass and mass ratio (#431)
- new "wind tunnel" setup and injection module (thanks to Mike Lau; #470)
- new "solar system" setup for placing solar system planets and minor bodies by downloading their published orbital elements (#430)
- bugs fixed with asteroid wind setup (#463)
- bug fix with units in GR tidal disruption event setup (#432)
- bug fix with initial velocities in disc setup with self-gravity and dust, properly compute enclosed mass for both gas and dust (thanks to Cristiano Longarini; #427)
- bug fix with turbulent stirring setup (thanks to Terry Tricco; #449)

Analysis/moddump utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~
- cleanup and further enhancements to common envelope analysis routines (thanks to Miguel Gonzalez-Bolivar; #467,#462)
- moddump_sink displays correct value of sink luminosity (#439)
- analysis routine for radio emission from tidal disruption events (thanks to Fitz Hu; #472)
- new analysis routine to compute time of dust formation (`Bermudez-Bustamante et al. 2023 <https://ui.adsabs.harvard.edu/abs/2024arXiv240103644B/abstract>`__)

Other
~~~~~
- github actions workflow now checks that running phantom on the .in file for one timestep succeeds following setup procedure
- github actions workflow checks compilation of phantom+mcfost
- phantom is now enforced to compile without any compiler warnings with gfortran on the master branch
- further work to reduce ugly ifdefs in phantom codebase (#55)
- various bugs with uninitialised variables fixed; all setups now checked with DEBUG=yes


v2023.0.0 - 10th Mar 2023
-------------------------

Physics
~~~~~~~
- Dust nucleation: chemical network allows for self-consistent (carbon-rich) dust formation in simulations once the temperature drops below the condensation temperature, e.g. in AGB star winds (`Siess et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022A%26A...667A..75S/abstract>`__). Just set DUST_NUCLEATION=yes
- Implementation of implicit flux limited diffusion solver for radiation based on `Whitehouse & Bate 2004 <https://ui.adsabs.harvard.edu/abs/2004MNRAS.353.1078W>`__, `Whitehouse, Bate & Monaghan (2005) <https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1367W>`__ and `Bate & Keto (2015) <http://adsabs.harvard.edu/abs/2015MNRAS.449.2643B>`__. The implementation is a direct port of the solver used in Matthew Bate's sphNG code, generously contributed by Matthew and ported to phantom by Mike Lau and Daniel Price.
- Compute and store a mass accretion rate for each sink particle based on the accreted mass over one dtmax interval. This can be used in MCFOST to provide accretion luminosity from each sink particle (`Borchert et al. 2022a <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510L..37B>`__, `2022b <https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.4436B>`__)
- Allow for use of PdV work and shock heating as source terms in MCFOST when performing live-coupled simulations with phantom+MCFOST (thanks to Elli Borchert and Sahl Rowther, #344; see appendix of `Borchert et al. 2022b <https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.4436B>`__)
- New equation of state that implements gas plus radiation pressure plus recombination of Hydrogen and Helium (`Lau et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.4436B>`__; #267)
- New stratified disc locally isothermal equation of state added (thanks to Caitlyn Hardiman; #268, #278, #308)
- GR simulations now work with the gas+radiation equation of state by implementing an option to evolve "s" as the entropy instead of P/rho^gamma (thanks to Fitz Hu; #250, #324)
- New spherical raytracer for sink particle radiation / dust acceleration, particularly in the context of stellar winds (thanks to Mats Esseldeurs, Lionel Siess, Ward Homan)
- option to use reconstruction on velocities to reduce dissipation away from shocks in the GR code (where the standard Cullen & Dehnen shock switch cannot be used)
- added conservation check which stops the MHD code if hdivB/B is too large AND hdivbbmax_max == 1 (thanks to James Wurster)
- added an override to force sink creation if density is too high; also an addition step in the barotropic equation of state at very high densities (thanks to James Wurster)
- option for heating from sink particles, which adds heats all particles within the softening radius of the sink according to the sink particle luminosity (#216)
- mass-weighted interpolations of dust-gas quantities are now default in dust growth (#377; thanks to Stephane Michoulier)
- option for Aeolian erosion in dust growth module (#365, #372; thanks to Stephane Michoulier)
- Added possibility of varying beta cooling with radius as a power law (#361; thanks to Cristiano Longarini)
- option for damping boundary conditions for comparing phantom to simulations of planet-disc interaction with grid codes (#351)


Setup
~~~~~
- can now set up hierarchical quadruple star systems in setup_disc (thanks to Amena Faruqi; #355)
- setup and simple cooling function added for cooling shock problem from `Creasey et al. (2011) <https://ui.adsabs.harvard.edu/abs/2011MNRAS.415.3706C>`__
- can now give --maxp=1e7 flag to phantomsetup instead of compiling with MAXP= which is now obsolete
- major reorganisation of setup_star into neat subroutines with high level functions (#297, also #203)
- dumps written during relax_star now usable as starting file
- various improvements to setting up stars with variable composition (thanks to Mike Lau)
- can now setup triple stars in the stellar wind setup (thanks to Lionel Siess), see e.g. `Maes et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...653A..25M/abstract>`__, `Malfait et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021A&A...652A..51M>`__
- major overhaul of setup_sphereinbox to include turbulence and one-fluid dust; and many more optional input variables (thanks to James Wurster)
- option for particle shuffling in the Sedov blast wave setup (thanks to James Wurster)
- simplified setup options when adding planets in setup_disc
- sensible physical units chosen for special relativistic shock tubes (thanks to Fitz Hu)
- added versatility to setup of power-law size distribution when setting up dust (thanks to Mark Hutchison)
- new SETUP options for isothermal (dusty) self gravitating disc setups (thanks to Cristiano Longarini)
- fixed default cooling in disc_setup (#114; thanks to Benedetta Veronesi)

Bugs
~~~~
- Bugs fixed with using dump files from sphNG to phantom (thanks to Alison Young; #343)
- Bug fix with hsoft=0 when reading MESA file without softening in star setup
- Various bug fixes when setting up Bonnor-Ebert density profiles in the sphere-in-box setup (thanks to James Wurster; #303)
- Bug fix with automated download of data files
- Bug fix in implicit cooling when du/dt goes to zero (#328; thanks to Lionel Siess)
- Various bug fixes with dust nucleation and cooling (Siess)
- Bug fix with artificial conductivity when employing the Minkowski metric in General Relativity
- various issues compiling phantom with MCFOST fixed (thanks to Christophe Pinte; #199)
- seg fault in dump file read utilities fixed if attempting to read an array that has not been allocated
- various bugs with phantom2hdf5 fixed (thanks to Stephen Nielson #368)
- Bug fix reading the default star data files from the data/ directory
- Bug fix in `Farris et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...783..134F>`__ equation of state (#282; thanks to Enrico Ragusa)
- Bug fix for dustfrac and dust-to-gas ratio in dustydisc setup (#273; thanks to Mark Hutchison)
- Bug fix in the initialisation of dustfrac in setup_disc
- Fixed critical bug in calculation of Teff in stellar wind setups (thanks to Lionel Siess)
- Fix the makefiles qscript target to write correct slurm scripts for MPI jobs (#269)
- Bug fix with sink particle creation in MPI (Chan, Liptai via ADACS; #234)
- Bug fix with Koyama & Inutuska cooling; works now for both implicit & explicit (thanks to James Wurster)
- Bug fixes in gravitational wave inspiral with star+sink (thanks to Martina Toscani; #367)
- libtool error fixed, use ar rcs to create libraries

Utils
~~~~~
- New moddump utility for importing sphNG dumps containing sink particles into phantom (thanks to Alison Young)
- New moddump to add a flyby to an evolved simulation (sink particle in parabolic orbit); thanks to Cristiano Longarini
- improved common envelope analysis routines (Lau, Gonzalez, Nielson; #334)
- moddump_sink allows modification of all sink particles in the simulation
- better error messages in moddump_binary if not enough memory is allocated to add a second star (Lau)
- cleanup of moddump_binary to improve code clarity (#257)
- moddump_binary can be used to set up binary of two stars with sink-particle stellar cores (thanks to Mike Lau; #362)
- moddump_rotate to add solid body rotation to a sphere of gas (Lau; #290)
- diffdumps returns a non-zero exit code if files differ (Chan, Liptai via ADACS; #237)


Performance
~~~~~~~~~~~
- major optimisation of MPI communication to avoid bottleneck of openMP code (#310; Chan, Liptai via ADACS)
- optimisation of particle balance between MPI threads (Chan, Liptai via ADACS; #316)
- timing information written in the log file for local and remote parts of density and force (Chan, Liptai via ADACS; #271)
- various MPI and OpenMP memory allocation optimisations and bug fixes (Chan, Liptai via ADACS; #209, #262; #243)
- Optimisations to reduce unnecessary calls when compiling with `MPI=yes` but running with only 1 MPI task (Chan, Liptai via ADACS; #259)

Other
~~~~~
- Switched off the automatic decrease of dtmax if the time between dumps is too large (#342)
- added option to create restart dumps if we go > 24h without a dump (#352; thanks to James Wurster)
- better help for SETUP= flag in Makefile
- configuration added for Flatiron cluster (SYSTEM=rusty and SYSTEM=popeye; thanks to Mike Lau)
- further work to remove unnecessary ifdefs (#55)
- Added MPI unit tests to the testsuite (Chan, Liptai via ADACS; #220, #222, #229, #235, #217, #322)
- major reorganisation of cooling modules; added cooling_solver, cooling_functions and other modules
- bots script can be run as a pre-commit action (Chan, Liptai via ADACS; #223, #317)
- Makefile split into Makefile_setups, Makefile_systems and Makefile_qscripts to avoid clutter (Liptai via ADACS; #261; see #253)
- Timing hierarchy drawn in a nicer way using a tree diagram (Chan via ADACS; #254)
- makefile exit codes are propagated through to calling scripts (#256)
- test suite is now also run using ifort on github runners (Chan, Liptai via ADACS; #228)
- github actions checks on pull requests are now run in parallel (Chan, Liptai via ADACS; #224)
- if dt is too small, exit in step with useful information rather than in get_ibin (thanks to James Wurster)
- option to run bots on staged files only (#213)

Documentation
~~~~~~~~~~~~~
- added list of pre-cooked setups (SETUP=blah) to docs
- added list of all equation of state options (#311)
- additional documentation on the file format specification
- Documentation added regarding Sarracen
- Machine-specific instructions added for Kennedy (St. Andrews) and DiAL
- Documentation for self-gravitating and gravitationally unstable disc setups (thanks to Cristiano Longarini)


v2022.0.0 - 17th Jan 2022
-------------------------

Physics
~~~~~~~
- Option for gravitational wave emission in quadrupole approximation from any simulation (`Toscani et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510..992T/abstract>`__)
- Further improvements to wind injection/line cooling/dust formation (`Siess et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022A%26A...667A..75S/abstract>`__)
- Ideal + radiation + H/He ionisation equation of state (Lau, Hirai)
- Allow for variable composition (X, Z, mu) in stars (Lau, Hirai)
- Radiative feedback implemented via MCFOST based on sink particle Mdot (`Borchert et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510L..37B/abstract>`__)
- Sink particles can now merge (thanks to James Wurster; #172)
- Option for thermal energy floor / minimum temperature (Wurster)
- Fixes/improvements to implicit cooling (Wurster)
- Updated NICIL library for non-ideal MHD coefficients to v2.1 (Wuster; #115)

Setup
~~~~~
- Further improvements to automated relax-star procedure and to setup_star in general (See Appendix C of `Lau et al. 2022 <https://ui.adsabs.harvard.edu/abs/2021arXiv211100923L/abstract>`__)
- Real star profiles allowed in GR tidal disruption event setup and moddump (Hu, Sharma)
- Set up for an hierarchical triple system embedded in a circum-triple disc (`Ceppi et al. 2022 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.514..906C/abstract>`__; `2023 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.520.5817C/abstract>`__; #102, #110)
- Firehose setup added for testing tidal disruption flows

Bugs
~~~~
- Bug fixed where showarrays utility did not work with single precision files (#164)
- Bug fix with particle IDs tracking with MPI (Chan, Liptai via ADACS)
- Bug fix with particle waking with MPI (Chan, Liptai via ADACS)
- Fix missing sink force reduction during initial setup (Chan, Liptai via ADACS)
- Fix reading integer arrays from native phantom files (Chan, Liptai via ADACS)
- Bug fix with seg fault in test suite during sink particle creation (#132)
- Bug fixes with molecular line cooling (Homan)
- Bug fix with timestep during particle injection (Wurster)
- Bug fixes with disc setup (Ragusa)

Utils
~~~~~
- improved common envelope analysis routines (Lau)
- some issues with hdf5 read/write fixed (Chan)
- diffdumps utility now works with MPI
- import/export to Kepler 1D stellar evolution code (Sharma, Heger)
- bug fixes in dustydisc analysis
- fix unit conversion of distance and mass in moddump dustadd.f90 (Longarini)

Other
~~~~~
- entire build and test suite now checked during continuous integration/ github workflows (Chan, Liptai via ADACS)
- fixed warnings regarding temporary array creation when compiling with ifort
- compiler warnings fixed

v2021.0.0 - 25th Jan 2021
-------------------------

Physics
~~~~~~~
- General relativistic hydrodynamics in Kerr, Schwarzschild and Minkowski metrics (`Liptai & Price 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..819L/abstract>`__)
- Major improvements to wind injection/line cooling/dust formation (contributed by Lionel Siess)
- Interface with KROME chemistry library for chemistry+cooling (contributed by Ward Homan)
- Multigrain dust-as-particles (i.e. multiple large grain species) now works (Mentiplay et al. 2020)
- Overdamping problem for small grains fixed when dust is simulated with particles (`Price & Laibe 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3929P/abstract>`__)
- Stepinski-Valageas dust growth algorithm works with both dust-as-mixture and dust-as-particles (`Vericel et al. 2020 <https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.2318V/abstract>`__)
- Preliminary implementation of flux limited diffusion radiation hydro, explicit timestepping only (Biriukov, Borchert)
- Added "ideal + radiation" equation of state (Lau)
- Various improvements to asteroid wind injection modules (Trevascus, Nealon, see `Trevascus et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021MNRAS.505L..21T/abstract>`__)
- gravitational wave inspiral via external force works with sink particles (Toscani)
- gravitational wave emission computed automatically using Quadrupole approximation (Toscani)
- NICIL library for non-ideal MHD diffusion coefficients updated to v1.2.6 (Wurster)

Setup
~~~~~
- Major improvements to setup procedure when mapping MESA stars into phantom, including ability to replace core with softened point mass particle (Lau, Hirai, Gonzalez, de Marco, Reichardt)
- automated relaxation of stellar profiles in phantomsetup using asynchronous shifting (relax-o-matic), similar to `Diehl et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015PASA...32...48D/abstract>`__.
- added random-but-symmetric option to set_sphere, giving arbitrary density profile with centre of mass exactly at origin
- Various setup routines for GR simulations, e.g. setup_grtde for tidal disruption problems (Liptai et al. 2019)
- Dust growth setups (growingdisc,testgrowth)
- Shocktube setup includes special relativistic shock tubes (`Liptai & Price 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..819L/abstract>`__), radiative shocks (Borchert, Biriukov) and dusty shocks with multiple grain sizes (Mentiplay et al. 2020). Also added ability to smooth initial shock front if desired (c.f. Mentiplay et al. 2020)
- Ability to set up initial density profile as Bonnor-Ebert sphere in star formation setups (Wurster)
- Disc setup with dust now shows the percentage of particles not satisfying the terminal velocity approximation (Ragusa)

Bugs
~~~~
- Various bug fixes with radiation hydrodynamics with flux-limited diffusion (Borchert, Biriukov)
- Various issues with live phantom-mcfost simulations fixed (Pinte)
- Various issues with multigrain dust calculations fixed (Mentiplay)
- Various issues with dust growth fixed (Vericel)
- Fixed bug with artificial conductivity being incorrect when non-ideal equations of state were used (Lau, Hirai)
- Bug fix with sink particles not crossing periodic boundaries (Wurster)
- Now check for dead particles present in dump files and remove them
- bug fixes with phantom2pdf_amr for computing volume-weighted probability density functions
- bug fix in analysis_disc regarding where the origin is assumed to be (Nealon)
- bug fix with memory allocation for dvdx, possibly meaning shock viscosity switch was not applied properly

Utils
~~~~~
- splitpart and mergepart utilities added for splitting and merging particles, can be used to continue a simulation at a lower/higher resolution (Nealon, Wurster, Price)
- growth_to_mcfost utility added for radiative transfer post-processing of simulations with dust growth (Vericel)
- major improvements to analysis_common_envelope (Lau, de Marco)
- various issues with phantom2hdf5 utility fixed (Mentiplay, Pinte)
- moddump_sink can be used to modify various sink particle properties by hand (Pinte, Lau)
- analysis_tde for analysing GR tidal disruption calculations (Liptai)
- ev2dot utility for taking derivative of any column in a .ev file (Liptai)
- evcut, evhead, evcat utilities for manipulating/combining .ev files (Liptai)
- combinedustdumps utility for stacking dust-gas simulations performed with single grain sizes (Mentiplay, Price)

Build
~~~~~
- code compiled into more modular and re-usable libraries (libsetup, libphantom)
- phantomtest is now compiled as a separate binary to phantom, where phantomtest depends on phantom but not the other way around
- phantomsetup now compiles using libsetup to keep dependencies clean

Other
~~~~~
- Added rkill option to kill particles outside a certain radius, useful for simulations with particle injection (Veronesi)
- get_derivs_global routine simplifies a lot of code in the test suite
- Remaining static memory allocation removed, phantom itself no longer requires MAXP= flag to increase the particle number beyond 10^6. However, this remains necessary in phantomsetup.
- migrated repositories and continuous integration tests to github
- simplified code due to pressure now being stored on particles, use "conservative to primitive" routine to convert conserved variables to primitive variables
- automated documentation of code modules via sphinx-fortran

Performance
~~~~~~~~~~~
- pressure, temperature and sound speed are now stored on particles, removing the need to call the equation of state routine on neighbours. This improves performance of simulations using tabulated equations of state. Equation of state is now only called once per timestep.


v1.4.0 - 20th January 2020 - 1b48489
------------------------------------

Physics
~~~~~~~

-  Working implementation of dust growth using Stepinski-Valageas 1997
   model (Vericel)
-  updated MCFOST interface for live radiation calculations
-  further improvements to Roche Lobe injection (Worpel)
-  Fixed issue of initial violent response of inner disc - no longer
   correct orbital velocities for surface density turnover in inner disc
-  default alpha_AV is 1.0 instead of 0.1 when using CONST_AV = yes
-  warning added about particles with zero sound speed
-  preliminary work to incorporate Shen (2012) equation of state for
   Neutron stars
-  Support for multi grain dust with multiple large grain species
-  (non-ideal MHD) updated nicil cosmic ray ionisation library to V1.2.6

Setup
~~~~~

-  Bug fixes with disc setup routines
-  Default npart is 10^6 in disc setup
-  better warnings about validity of terminal velocity approximation
-  moved default settings for particle arrays into init_part routine
-  cluster setup reads/writes .setup file

Bugs
~~~~

-  Issue with zero grain sizes upon restart fixed, now checked for in
   checksetup
-  Issue with one fluid setups not working on stable branch fixed
-  Numerous bugs fixed with dust growth implementation (Vericel)
-  git version info prints correctly when running test suite
-  now call update_externalforce before checksetup is run to avoid
   problem with extern_binary
-  Default units changed in galaxies setup to avoid momentum
   conservation warning
-  bug fixes for barotropic ieos=8
-  bug fix with fatal error for particles with energy equal to zero (now a warning)
-  (pyphantom) Added try statements to avoid errors when loading utherm, temperature and bxyz
-  (ptmass) bug fix in bookkeeping of why sink was not created
-  (test_derivs) more precise test of artificial viscosity terms for DISC_VISCOSITY=no,
-  passes test suite when KERNEL=quintic
-  MPI thread-safe downloading of datafiles
-  BUG FIX with memory allocation for dvdx; possibly affecting viscosity switch if DISC_VISCOSITY=no

Performance
~~~~~~~~~~~

-  Improved parallelisation of root node construction in kdtree build

Build
~~~~~

-  Nightly code performance (openMP only) now checked automatically

Utils
~~~~~

-  read_array_from_file in utils_dumpfiles can be used to read real*4
   arrays not read during read_dumpfile (e.g. luminosity)
-  kernels script updated to Python 3
-  several python scripts (evcat,evcut,evhead,ev2dot) added for messing
   around with .ev files (#, Liptai)
-  phantom2hdf5 added to convert dump files to hdf5 format (Mentiplay,
   Liptai)
-  moddump to remove particles inside/outside some radius (Vericel)
-  disc analysis utility now assumes that the disc is around the first
   sink if sinks are present
-  combinedustdumps utility to stack different grain sizes from
   single-grain calculations now works with automatic memory allocation

Other
~~~~~

-  less verbose output during memory allocation
-  update_test_scores routine used to avoid repeated code in test suite
-  optional HDF5 output for easy reading of dump files in Python via
   Plonk (Mentiplay, Liptai)
-  automatic correction of “if(” to “if (” by format-bot


v1.3.0 - 22 Feb 2019 - 4d45cb3
------------------------------

Physics
~~~~~~~

-  Multigrain dust simulations with multiple large grains now possible (Mentiplay). This complements the multigrain method used for small grains, but simulating small and large grain populations simultaneously is not yet fully functional
- Further updates to dust growth algorithms (Vericel)
-  Much improved wind injection routines (Price, Siess)
- Improvements to Roche lobe injection module (Worpel)
- Injection modules can now provide an additional timestep constraint where needed
-  One fluid dust uses method of `Ballabio et al.  (2018) <http://ui.adsabs.harvard.edu/abs/2018MNRAS.477 .2766B>`__ to prevent negative dust fractions
-  can now set a maximum density after which the simulation will end, also dtmax will dynamically decrease/increase if density increases too rapidly (Wurster)
-  removed obsolete and unused etamhd fixed resistivity variable
- reduced timestep from physical viscosity force by factor 0.4: this has been found to lead to much better convergence of disc simulations that use this method (Nixon)

Bugs
~~~~

- bug fix with momentum conservation in two fluid dust-gas drag when ISOTHERMAL=yes
- array bounds error in analysis_tde fixed
- bugfix in read options for externbinary module

Tests
~~~~~

-  test for momentum and energy conservation in two fluid dust-gas drag
-  code performance is now checked nightly against a suite of benchmarks
-  sends error code to system if a fatal error happens (Pinte)
-  added check on the conservation of angular momentum with dust/gas

Setup
~~~~~

-  Binary disc setup uses Farris et al. (2014) locally isothermal equation of state for discs around more than one star
-  Disc setup routine modularised and made more general (Mentiplay)
- gwdisc setup now allows disc inclination (`Pereira et al. 2019 <http://ui.adsabs.harvard.edu/abs/2019MNRAS.4 84...31P>`__)
-  setup_star given fairly major restructure so logic is clearer; more cleanly split interactive from non-interactive parts
-  Flyby setup updated with the following roll angle convention: incl=0 => prograde orbit (disc and perturber anti-clockwise; incl=180 => retrograde orbit (disc anti- and perturber clockwise). See `Cuello et al. 2019 <http://ui.adsabs.harvard.edu/abs/2019MNRAS.483.4114CL>`__
-  minor fixes to dustyshock and dustywave setups (Hutchison)
- binary_w in setup_disc is now 270 degrees by default
- asteroidwind setup added
- added option to setup a settled dusty disc, working with both one and 2 fluid (Dipierro)

Build
~~~~~

-  version number and git sha now written to dump file headers
- memory is now allocated at runtime for main arrays in Phantom (Chan). This avoids the need to recompile with MAXP= when you change the particle number.  Only applies to main phantom binary at present, not to phantomsetup.
- many compiler warnings fixed
- cleanup of evolve module
- obsolete preprocessor flags -DSORT_RADIUS_INI T and -DDUSTFRAC deleted
-  you can now supply JOBNAME= when making job scripts with make qscript, otherwise it continues to choose delightful random words

Analysis
~~~~~~~~

-  Multigrain post-processing works properly with MCFOST
- phantomevcompare will not duplicate data when merging files
-  further integration with MCFOST
- analysis disc planet prints the effective tilt between the inner and outer disc (Nealon)
-  disc analysis now defaults to sorting particles by cylindrical radius - this should fix any discrepancies that may have been occurring.  Deliberately made it very hard not to chose this option (Nealon)
-  disc analysis now returns the total angular momentum components as well (Nealon)
- precession files: these can now be made even if the first file input is not the first file of the simulation (Nealon)
- utils_disc now handles an eccentric disc - bins are defined by semi-major axis, not by radius (Nealon)
- analysis_dustydis c
- Added check Ltot!=0 to prevent NaNs in the output (Ragusa)
- moddump_extenddis c implemented to extend an existing disc simulation in radius (Nealon)
-  disc scale height now calculated from particle positions but works perfectly with a warped disc (Nealon)

Other
~~~~~

-  phantom outputs helpful error message if .setup file is given on command line instead of .in file


v1.2.0 - 20 Jun 2018 - d339b10
------------------------------

This release corresponds to the accepted version of the Phantom paper (v2 on arXiv). Changes compared to v1.1.0:

Physics
~~~~~~~

- Multigrain dust algorithm implemented `(Hutchison, Price & Laibe 2018) <http://ui.adsabs.harvard.edu/abs/2018MNRAS.476.2186H>`__

Build
~~~~~

- SYSTEM=ozstar added


v1.1.0 - 5 Apr 2018
-------------------

Physics
~~~~~~~

-  Helmholtz equation of state implemented (Tricco)
- preliminary work on dust growth (Vericel)

Bugs
~~~~

-  bug fix with magnetic fields on boundary particles
-  bug fix with incorrect fatal error on centre of mass non-conservation
-  angular momentum now conserved during sink particle accretion (#17, Wurster)
- issues with git-lfs fixed
- bug fix with write of B-field to small dump files

Tests
~~~~~

-  setupbot: Nightly checks that phantomsetup does not require unspecified user input

Setup
~~~~~

-  better defaults in several setups so we pass setupbot checks
- set_slab utility routine added for 2D-in-3D setups

Build
~~~~~

- SYSTEM=raijin added


v1.0 - 13 Mar 2018
------------------

Physics
~~~~~~~

-  working MPI implementation (Chan)
-  more robust algorithm for one fluid dust (Ballabio+ 2018)
-  dust algorithm (one fluid/two fluid) chosen at runtime not compile time
-  particle waking with individual timesteps re-implemented (Wurster; 45fae9b)
-  universal disc setup routine (Mentiplay)
-  setup added for flyby simulations (Mentiplay, Cuello)
-  CO cooling implemented (Glover)
-  magnetic field evolves B/rho rather than B (Tricco, Price)
-  stellar wind routine works out-of-the-box (Toupin)
-  improvemements to Galactic Centre winds and cooling (Russell, Price)
-  NICIL updated to v1.2.3 (Wurster)

Bugs
~~~~

-  bug with drag in two fluid dust-gas when hj > hi fixed (Dipierro)
-  updates/bug fixes to MESA Equation of state tabulation
-  bug fix with energy conservation with softened sink particles
-  bug fix with self-gravity + multiple particle types

Tests
~~~~~

-  nightly checks for non-ideal MHD added
-  self gravity checked for all particle types
-  testsuite checked nightly with MPI


v0.9 - 14 Feb 2017
------------------


This is the first public release of Phantom, alongside arXiv paper.

Contains:

-  hydro
-  sink particles
-  self-gravity
-  MHD
-  dust (two fluid and one fluid)
-  ISM chemistry and cooling
-  physical viscosity
-  non-ideal MHD
-  external forces including corotating frame, Lense-Thirring
   precession, P-R drag, fixed binary
