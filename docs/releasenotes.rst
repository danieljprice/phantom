Release notes
=============

v2021.0.0 - 25th Jan 2021
-------------------------

Physics
~~~~~~~
- General relativistic hydrodynamics in Kerr, Schwarzschild and Minkowski metrics (`Liptai & Price 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..819L/abstract>`__)
- Major improvements to wind injection/line cooling/dust formation (contributed by Lionel Siess)
- Interface with KROME chemistry library for chemistry+cooling (contributed by Ward Homan)
- Multigrain dust-as-particles (i.e. multiple large grain species) now works (Mentiplay et al. 2020)
- Overdamping problem for small grains fixed when dust is simulated with particles (`Price & Laibe 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.495.3929P/abstract>`__)
- Stepinski-Valageas dust growth algorithm works with both dust-as-mixture and dust-as-particles (Vericel et al. 2020)
- Preliminary implementation of flux limited diffusion radiation hydro, explicit timestepping only (Biriukov, Borchert)
- Added "ideal + radiation" equation of state (Lau)
- Various improvements to asteroid wind injection modules (Trevascus, Nealon)
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
- code performance is now checked nightly against a `suite of benchmarks <https://bitbucket.org/danielprice/phantom-benchmarks>`__
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
