Setting up stars and tidal disruption events
============================================

Setting up and relaxing a star
------------------------------

First, follow the usual procedure for initiating a new simulation with
phantom. We’ll use the “star” setup, but you can also use the
“polytrope” or “neutronstar” configurations (the first two use self-gravity
for the star, the last one uses an external potential). For tidal disruption
events in general relativity use“grtde”. That is:

make a new directory and write a local Makefile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ mkdir star
   $ cd star
   $ ~/phantom/scripts/writemake.sh star > Makefile

compile phantom and phantomsetup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ make
   $ make setup
   $ ls
   Makefile    phantom*    phantomsetup*

run phantomsetup
~~~~~~~~~~~~~~~~

::

   $./phantomsetup star

   star.setup not found: using interactive setup

   1) Uniform density profile
   2) Polytrope
   3) Density vs r from ascii file
   4) KEPLER star from file
   5) MESA star from file
   6) Piecewise polytrope
   7) Evrard collapse
  Enter which density profile to use ([1:7], default=1): 2
  Setting up Polytrope
  Enter mass unit (e.g. solarm,jupiterm,earthm) (blank="blank",default="solarm"):
  Enter distance unit (e.g. au,pc,kpc,0.1pc) (blank="blank",default="solarr"):
  Enter the approximate number of particles in the sphere ([0:], default=100000):
  Enter the desired EoS for setup (default=2):
  Enter gamma (adiabatic index) ([1.000:7.000], default=1.667):
  Enter the mass of the star (code units) ([0.000:], default=1.000):
  Enter the radius of the star (code units) ([0.000:], default=1.000):
  Relax star automatically during setup? (default=no): y
   Writing star.setup
  STOP please check and edit .setup file and rerun phantomsetup

once you have answered the questions once, they are written to a file
called star.setup, which you can use for subsequent runs of phantomsetup
without having to answer prompts (see below)

relax the star
~~~~~~~~~~~~~~

Open the star.setup file and make sure the parameter “relax_star” to True::

                   relax_star =       T    ! artificial damping of velocities (if on, v=0 initially)

Then run phantomsetup::

   ./phantomsetup star.setup

Which will generate a sequence of snapshots of the relaxation process::

   RELAX-A-STAR-O-MATIC: Etherm:  0.499     Epot: -0.854     R*:   1.00
       WILL stop WHEN: dens error <   1.00% AND Ekin/Epot <   1.000E-07 OR Iter=0

   -------->   TIME =    0.000    : full dump written to file relax_00000   <--------

    Relaxing star: Iter   1/1000, dens error: 20.39%, R*:  0.978     Ekin/Epot:  3.157E-03
    Relaxing star: Iter   2/1000, dens error: 15.91%, R*:  0.978     Ekin/Epot:  4.858E-03
    Relaxing star: Iter   3/1000, dens error: 13.05%, R*:  0.978     Ekin/Epot:  3.303E-03
    Relaxing star: Iter   4/1000, dens error: 11.11%, R*:  0.977     Ekin/Epot:  2.276E-03
    Relaxing star: Iter   5/1000, dens error:  9.69%, R*:  0.975     Ekin/Epot:  1.655E-03

   -------->   TIME =   0.3756E-01: full dump written to file relax_00001   <--------
    Relaxing star: Iter   6/ 500, dens error: 22.17%, R*:   1.03     Ekin/Epot:  1.280E-03

This process will stop when either the density error is less than 1\% and the ratio of kinetic to potential energy
is less than the specified tolerance, OR the number of iterations reaches the specified maximum.

Once complete, you should obtain a relaxed initial conditions snapshot with the correct density and pressure profiles::

   -------->   TIME =    0.000    : full dump written to file star_00000.tmp   <--------


    input file poly.in written successfully.
    To start the calculation, use:

    ./phantom star.in

evolve the star
~~~~~~~~~~~~~~~

If you really want to, you can evolve the star in isolation with the regular code to double-check that the relaxation process worked::

    ./phantom star.in

check the output
~~~~~~~~~~~~~~~~

::

   splash star_0*

(if you have version 2 of splash, the relevant command is "ssplash")

Putting the star on an orbit for a tidal disruption event
---------------------------------------------------------

If you used the “tde” or "grtde" setup then simply compile :doc:`moddump <moddump>`::

   $ make moddump

otherwise you need to specify the tidal moddump file::

   $ make moddump MODFILE=moddump_tidal.f90

Then run moddump on your relaxed star::

   $ ./phantommoddump star_00000 tde 0.0
   ...
   ...
   ...
    writing moddump params file tde.tdeparams
     Edit tde.tdeparams and rerun phantommoddump

When you first run this, a “tde.tdeparams” file will be created. Edit
this to set the star on your desired orbit, and then rerun
phantommoddump::

   # parameters file for a TDE phantommodump
                   beta =       1.000    ! penetration factor
                     mh =   1.000E+06    ! mass of black hole (code units)
                     ms =       1.000    ! mass of star       (code units)
                     rs =       1.000    ! radius of star     (code units)
                  theta =       0.000    ! stellar rotation with respect to x-axis (in degrees)
                    phi =       0.000    ! stellar rotation with respect to y-axis (in degrees)
                     r0 =        490.    ! starting distance

After this you can simply run phantom::

   $ ./phantom tde.in

Adding a magnetic field to the star
-----------------------------------

compile phantommoddump
~~~~~~~~~~~~~~~~~~~~~~

The module used to compile this utility is specified using MODFILE= in
phantom/build/Makefile. The default for the “polytrope” setup is
currently moddump_spheres.f90::

   MODFILE=moddump_spheres.f90

Change this to moddump_default.f90. You can do this temporarily on the
command line by compiling phantommoddump as follows::

   make moddump MODFILE=moddump_default.f90 MHD=yes

run phantommoddump
~~~~~~~~~~~~~~~~~~

::

   $ ./phantommoddump
   PhantomSPH: (c) 2007-2023 The Authors

    Usage: moddump dumpfilein dumpfileout [time] [outformat]

in our case we want::

   ./phantommoddump star_00010 magstar_00000

which will give some errors::

    ERROR! MHD arrays not found in Phantom dump file: got            0

but then prompt you to add magnetic fields::

   add/reset magnetic fields? (default=no): yes

you can follow the prompts to add uniform magnetic fields using this
routine.

now implement something decent in src/setup/set_Bfield.f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

you can either use the pre-cooked magnetic field setups in this routine,
or you can just make a new :doc:`moddump <moddump>` module that sets up the magnetic field in a custom way.
