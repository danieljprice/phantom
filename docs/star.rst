Setting up stars and tidal disruption events
============================================

Setting up and relaxing a star
------------------------------

First, follow the usual procedure for initiating a new simulation with
phantom. We’ll use the “polytrope” setup, but you can also use the
“star” or “neutronstar” configurations (the first two use self-gravity
for the star, the last one uses an external potential). For TDEs use
“tde”. That is:

make a new directory and write a local Makefile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ mkdir star
   $ cd star
   $ ~/phantom/scripts/writemake.sh polytrope > Makefile

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
     nprocs =            1
    ERROR opening star.setup

   Case  1 Uniform density sphere
   Case  2 Polytrope
   Case  3 Binary polytrope
   Case  4 neutron star from file
   Case  5 Red giant (Macquarie)
   Case  6 neutron star with piecewise polytrope EOS
   Case  7 Evrard collapse
   Case  8 KEPLER star from file
   Case  9 Helmholtz free energy eos star
   Enter which setup to use ([1:9], default=1):
   Setting up Uniform density sphere
   Enter mass unit (e.g. solarm,jupiterm,earthm) (blank="blank",default="solarm"):
   Enter distance unit (e.g. au,pc,kpc,0.1pc) (blank="blank",default="solarr"):
   Enter the approximate number of particles in the sphere ([0:666666], default=10000):
   Enter the radius of the star (code units) ([0.000:], default=1.000):
   Enter the mass of the star (code units) ([0.000:], default=1.000):
   Enter the Adiabatic index ([1.000:], default=1.667):
   Enter polyk (sound speed .or. constant in EOS calculation) ([0.000:], default=0.5000):
     set_sphere: Iterating to form sphere with approx        10000  particles
    set_sphere: Iterations complete: added      10659 particles in sphere

once you have answered the questions once, they are written to a file
called star.setup, which you can use for subsequent runs of phantomsetup
without having to answer prompts

relax the star
~~~~~~~~~~~~~~

Open the star.in file and set the parameter “damp” to 0.03:

::

                   damp =       0.030    ! artificial damping of velocities (if on, v=0 initially)

Also change dtmax if you want more frequent output:

::

                  dtmax =       1.000    ! time between dumps

Then run phantom:

::

   ./phantom star.in

In the log file you should see:

::

    Polytropic equation of state: P =   0.424304*rho^  1.666667

check the output
~~~~~~~~~~~~~~~~

::

   ssplash star_0*

Putting the star on an orbit for a tidal disruption event
---------------------------------------------------------

If you used the “tde” setup then simply compile moddump:

::

   $ make moddump

otherwise you need to specify the tidal moddump file

::

   $ make moddump MODFILE=moddump_tidal.f90

Then run moddump on your relaxed star

::

   $ ./phantommoddump star_00010 tde 0.0
   ...
   ...
   ...
    writing moddump params file tde.tdeparams
     Edit tde.tdeparams and rerun phantommoddump

When you first run this, a “tde.tdeparams” file will be created. Edit
this to set the star on your desired orbit, and then rerun
phantommoddump.

::

   # parameters file for a TDE phantommodump
                   beta =       1.000    ! penetration factor
                     mh =   1.000E+06    ! mass of black hole (code units)
                     ms =       1.000    ! mass of star       (code units)
                     rs =       1.000    ! radius of star     (code units)
                  theta =       0.000    ! stellar rotation with respect to x-axis (in degrees)
                    phi =       0.000    ! stellar rotation with respect to y-axis (in degrees)
                     r0 =        490.    ! starting distance

After this you can simply run phantom

::

   $ ./phantom tde.in

Adding a magnetic field to the star
-----------------------------------

compile phantommoddump
~~~~~~~~~~~~~~~~~~~~~~

The module used to compile this utility is specified using MODFILE= in
phantom/build/Makefile. The default for the “polytrope” setup is
currently moddump_spheres.f90

::

   MODFILE=moddump_spheres.f90

Change this to moddump_default.f90. You can do this temporarily on the
command line by compiling phantommoddump as follows:

::

   $ make moddump MODFILE=moddump_default.f90 MHD=yes

run phantommoddump
~~~~~~~~~~~~~~~~~~

::

   $ ./phantommoddump
   PhantomSPH: (c) 2007-2017 The Authors

    Usage: moddump dumpfilein dumpfileout [time] [outformat]

in our case we want:

::

   ./phantommoddump star_00010 magstar_00000

which will give some errors:

::

    ERROR! MHD arrays not found in Phantom dump file: got            0

but then prompt you to add magnetic fields:

::

   add/reset magnetic fields? (default=no): yes

you can follow the prompts to add uniform magnetic fields using this
routine.

now implement something decent in src/setup/set_Bfield.f90
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

you can either implement a general magnetic field setup in this routine,
or you can just make a new moddump module that sets up the magnetic
field in a custom way.
