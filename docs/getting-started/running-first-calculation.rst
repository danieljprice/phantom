Running your first calculation
==============================

First steps
-----------

After obtaining a fresh :doc:`git clone <gitinfo>`, enter
the root-level phantom directory:

::

   cd phantom

You will need to specify some additional environment variables to run
properly:

::

   export OMP_SCHEDULE="dynamic"
   export OMP_STACKSIZE=512M
   ulimit -s unlimited

**Put these commands in your ~/.bashrc file** so they are set every time
you login. The stacksize (ulimit -s) and OMP_STACKSIZE need to be set
because the neighbour caches in phantom are private to each thread, and
hence stored on the per-thread stack – this means the storage can exceed
the default openMP stack size setting, usually causing a seg-fault.

Choosing a compiler
---------------------

Fortran is a compiled language, for which many compilers exist. The free
compiler is called gfortran, but the Intel Fortran Compiler (ifort) is the main
commercial compiler and typically runs phantom ~20% faster than gfortran.

You can choose the compiler and other machine-specific configurations by specifying 
a SYSTEM variable corresponding to one of those listed in `phantom/build/Makefile_systems <https://github.com/danieljprice/phantom/blob/master/build/Makefile_systems>`__. 
You can do this by setting an environment variable (e.g. in bash/sh):

::

   export SYSTEM=ifort

The best way is to put the line above into your .profile/.bashrc so
that it is always set for the machine that you’re using. Alternatively you 
can specify this by including it on the command line every time you compile the code:

::

   make SYSTEM=ifort

A list of current options for the SYSTEM variable is given below, please contribute
if you have configured phantom for your local cluster so that others can benefit:

.. include:: systems-list.rst

Running an example calculation
------------------------------

To run phantom on an example problem, first create a directory for the calculation (best to do this SOMEWHERE ELSE, i.e. NOT in a subdirectory of phantom)::

   mkdir blast; cd blast

Then use the writemake script in the phantom/scripts directory to write a local Makefile::

   ~/phantom/scripts/writemake.sh sedov > Makefile

where “sedov” is the :doc:`name of a SETUP variable in phantom/build/Makefile_setups </user-guide/setups>`
(this argument is optional, but convenient as it means phantom when compiled in this directory will always compile for this setup). Then you should have::

   $ ls
   Makefile

Start by compiling both phantom and the phantomsetup utility::

   make; make setup

which produces::

   $ ls
   Makefile phantom* phantomsetup*

Run phantomsetup with the name you want to give the calculation::

   ./phantomsetup blast

which, after any prompting that the setup routine might perform,
produces both an initial dump file and a runtime options file::

   $ ls
   Makefile        blast.in        blast_00000.tmp     phantom*        phantomsetup*

The dump file (blast_00000.tmp) is a binary file that can be read by
`splash <http://users.monash.edu.au/~dprice/splash>`__. The .tmp appended to the filename is because phantomsetup does not compute the density, so the smoothing lengths and densities in the file are at this stage just guesses.

The input file (blast.in) contains all of the :doc:`runtime configuration
options </user-guide/infile>`. It’s fairly self-explanatory, but probably the main things to note are the end time and the time between dumps::

                   tmax =      0.2000    ! end time
                  dtmax =      0.0100    ! time between dumps

Essentially the code will run to time tmax in code units, writing a dump file at intervals of dtmax. So in the above example we will get 201 dump files produced in total. The other setting to note is the frequency of “full dumps”::

              nfulldump =          10    ! full dump every n dumps

This means that only 1 file in 10 is a “full” or “restart” dump from which the code can be restarted. Dump files in between contain only the particle positions and smoothing lengths, i.e. just enough information to make a movie but without wasting disk space.

The basic physics that is controllable at runtime (any physics that affects memory storage in phantom will require selection of compile-time options also) is contained within the block::

   # options controlling hydrodynamics, artificial dissipation
                   ieos =           2    ! eqn of state (1=isoth; 2=adiab; 3/4=locally iso (sphere/cyl); 5=two phase)
                  alpha =      0.0000    ! MINIMUM shock viscosity parameter (max = 1.0)
                 alphau =      1.0000    ! shock conductivity parameter
                   beta =      2.0000    ! beta viscosity
           avdecayconst =      0.1000    ! decay time constant for viscosity switches
                   damp =      0.0000    ! artificial damping of velocities (if on, v=0 initially)
           ipdv_heating =           1    ! heating from PdV work (0=off, 1=on)
         ishock_heating =           1    ! shock heating (0=off, 1=on)

To be able to use phantom effectively, you need to know enough about SPH
to know what these do. I suggest reading `Price
(2012) <http://ui.adsabs.harvard.edu/abs/2012JCoPh.231..759P>`__ as a first
step.

To run the code, just run phantom with the name of the input file::

   ./phantom blast.in

Note that the first thing that the code does is to compute density, and hence replaces the .tmp file with a “real” dump file::


   -------->   TIME =     0.0000: full dump written to file blast_00000   <--------


    input file blast.in written successfully.

    ---> DELETING temporary dump file blast_00000.tmp <---

Also, note that the input file (blast.in) is automatically updated every time a full dump is written. This means that if you enter the same command again::

   ./phantom blast.in

...then the calculation just picks up from the last full dump file written.

Visualising the output
----------------------

That’s what `splash <http://users.monash.edu.au/~dprice/splash>`__ is for! Use splash to look at the dump files produced by phantom::

   splash blast_0* -r 6 -dev /xw

For the Sedov example shown above, there’s even an exact solution included in splash (use o7 from the splash menu to change to spherical coordinates, then plot density as a function of radius, then use o8 to plot the Sedov exact solution).

.. important::
   If you have v2.x or earlier of splash,
   type ssplash instead of splash to read :doc:`the phantom native binary format</user-guide/dumpfile>`.

The .ev files, which are just ascii files containing global quantities as a function of time::

   $ more blast01.ev
   # [01        time]   [02        ekin]   [03      etherm]   [04        emag]   [05        epot]   [06        etot]   [07      totmom]   [08      angtot]   [09     rho max]   [10     rho ave]   [11          dt]   [12   totentrop]   [13     rmsmach]   [14        vrms]   [15        xcom]   [16        ycom]   [17        zcom]   [18   alpha max]
   0.0000000000E+00   0.0000000000E+00   1.0000000000E+00   0.0000000000E+00   0.0000000000E+00   1.0000000000E+00   0.0000000000E+00   0.0000000000E+00   1.0008253226E+00   1.0008253226E+00   0.0000000000E+00   6.6630010871E-01   0.0000000000E+00   0.0000000000E+00  -1.4823076577E-21  -8.1592566227E-18   1.2378986725E-14   1.0000000000E+00
   ...

The .ev files can be visualised using any standard plotting tool. For example you can use splash with the -ev (or
-e) option::

   splash -e blast*.ev

where column labels should be read automatically from the header of the .ev file.

For more detailed analysis of :doc:`Phantom dump files </user-guide/dumpfile>`, write yourself an analysis module for the :doc:`phantomanalysis </user-guide/analysis>` utility. Analysis modules exist for many common tasks, including interpolating to a 3D grid (both fixed and AMR), computing PDFs, structure functions and power spectra, getting disc surface density profiles, and converting to other formats.
