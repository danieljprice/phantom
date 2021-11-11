Compile-time configuration
==========================

Phantom uses a mix of compile-time and run-time configuration. The
:doc:`run-time configuration <infile>` is specified in the input file
(blah.in). Compile-time options are as follows:

Simple things can be changed on the command line. To change the maximum
number of particles use:

::

   make SETUP=disc MAXP=10000000

Setup block
~~~~~~~~~~~

The compile-time configuration of Phantom is specified using the SETUP
block in phantom/build/Makefile. For example the default disc setup is
“disc”:

::

   ifeq ($(SETUP), disc)
   #   locally isothermal gas disc
       FPPFLAGS= -DDISC_VISCOSITY
       SETUPFILE= setup_disc.f90
       ANALYSIS= analysis_disc.f90
       ISOTHERMAL=yes
       KNOWN_SETUP=yes
       MULTIRUNFILE= multirun.f90
       IND_TIMESTEPS=yes
   endif

The choice of C preprocessor flags (FPPFLAGS=) can be used to specify
particular physics. Otherwise these are specified using variables in
this block as follows:

Code modules
------------

+-----------------+-----------------+-----------------+-----------------+
| *Variable*      | *Setting*       | *Default value* | *Description*   |
+=================+=================+=================+=================+
| SETUPFILE       | .f90 file(s)    | setup_unifdis.F | The setup       |
|                 |                 | 90              | routine and any |
|                 |                 |                 | auxiliary       |
|                 |                 |                 | routines needed |
|                 |                 |                 | by phantomsetup |
+-----------------+-----------------+-----------------+-----------------+
| LINKLIST        | .f90 file(s)    | linklist.F90    | The neighbour   |
|                 |                 |                 | finding         |
|                 |                 |                 | algorithm       |
|                 |                 |                 | (fixed grid,    |
|                 |                 |                 | fixed           |
|                 |                 |                 | cylindrical     |
|                 |                 |                 | grid or kdtree) |
+-----------------+-----------------+-----------------+-----------------+
| ANALYSIS        | .f90 file(s)    | analysis_dthead | (optional) The  |
|                 |                 | er.f90          | analysis        |
|                 |                 |                 | routine and any |
|                 |                 |                 | auxiliary       |
|                 |                 |                 | routines used   |
|                 |                 |                 | by the          |
|                 |                 |                 | phantomanalysis |
|                 |                 |                 | utility         |
+-----------------+-----------------+-----------------+-----------------+
| SRCTURB         | .f90 file(s)    | forcing_ouproce | (optional)      |
|                 |                 | ss.f90          | Turbulence      |
|                 |                 |                 | driving routine |
|                 |                 |                 | (triggers       |
|                 |                 |                 | -DDRIVING)      |
+-----------------+-----------------+-----------------+-----------------+
| SRCINJECT       | .f90 file(s)    | inject_rochelob | (optional)      |
|                 |                 | e.f90           | Module handling |
|                 |                 |                 | particle        |
|                 |                 |                 | injection       |
|                 |                 |                 | (triggers       |
|                 |                 |                 | -DINJECT_PARTIC |
|                 |                 |                 | LES)            |
+-----------------+-----------------+-----------------+-----------------+
| MODFILE         | .f90 file(s)    | moddump.f90     | (optional)      |
|                 |                 |                 | Routine used by |
|                 |                 |                 | moddump utility |
|                 |                 |                 | (to modify an   |
|                 |                 |                 | existing dump   |
|                 |                 |                 | file)           |
+-----------------+-----------------+-----------------+-----------------+

Code performance and accuracy
-----------------------------

+-----------------+-----------------+-----------------+-----------------+
| *Variable*      | *Setting*       | *Default value* | *Description*   |
+=================+=================+=================+=================+
| MPI             | yes/no/openmpi/ | no              | compile with    |
|                 | zen/apac        |                 | MPI             |
|                 |                 |                 | parallelisation |
+-----------------+-----------------+-----------------+-----------------+
| OPENMP          | yes/no          | yes             | compile with    |
|                 |                 |                 | openMP          |
|                 |                 |                 | parallelisation |
+-----------------+-----------------+-----------------+-----------------+
| IND_TIMESTEPS   | yes/no          | no              | use individual  |
|                 |                 |                 | timesteps or    |
|                 |                 |                 | not             |
+-----------------+-----------------+-----------------+-----------------+
| DOUBLEPRECISION | yes/no          | yes             | use 8-byte      |
|                 |                 |                 | reals (no=4     |
|                 |                 |                 | byte)           |
+-----------------+-----------------+-----------------+-----------------+
| DEBUG           | yes/no          | no              | turns on        |
|                 |                 |                 | debugging flags |
|                 |                 |                 | (slow)          |
+-----------------+-----------------+-----------------+-----------------+

Memory usage
------------

+-----------------+-----------------+-----------------+-----------------+
| *Variable*      | *Setting*       | *Default value* | *Description*   |
+=================+=================+=================+=================+
| MAXP            | integer         | 1000000         | maximum number  |
|                 |                 |                 | of particles    |
|                 |                 |                 | (array size)    |
+-----------------+-----------------+-----------------+-----------------+
| MAXPTMASS       | integer         | 2               | maximum number  |
|                 |                 |                 | of point mass   |
|                 |                 |                 | particles       |
|                 |                 |                 | (array size)    |
+-----------------+-----------------+-----------------+-----------------+
| MAXNEIGH        | integer         | same as maxp    | maximum size of |
|                 |                 |                 | neighbour       |
|                 |                 |                 | arrays (default |
|                 |                 |                 | is safest but   |
|                 |                 |                 | can be lower to |
|                 |                 |                 | save memory)    |
+-----------------+-----------------+-----------------+-----------------+
| NCELLSMAX       | integer         | same as maxp    | maximum number  |
|                 |                 |                 | of cells in     |
|                 |                 |                 | fixed-grid      |
|                 |                 |                 | neighbour       |
|                 |                 |                 | finding         |
+-----------------+-----------------+-----------------+-----------------+

Physics
-------

+-----------------+-----------------+-----------------+-----------------+
| *Variable*      | *Setting*       | *Default value* | *Description*   |
+=================+=================+=================+=================+
| PERIODIC        | yes/no          | no              | periodic        |
|                 |                 |                 | boundaries      |
+-----------------+-----------------+-----------------+-----------------+
| ISOTHERMAL      | yes/no          | no              | determines      |
|                 |                 |                 | whether or not  |
|                 |                 |                 | to store        |
|                 |                 |                 | thermal energy  |
+-----------------+-----------------+-----------------+-----------------+
| GRAVITY         | yes/no          | no              | use             |
|                 |                 |                 | self-gravity or |
|                 |                 |                 | not             |
+-----------------+-----------------+-----------------+-----------------+
| MHD             | yes/no          | no              | use             |
|                 |                 |                 | magnetohydrodyn |
|                 |                 |                 | amics           |
|                 |                 |                 | or not          |
+-----------------+-----------------+-----------------+-----------------+
| DUST            | yes/no          | no              | use dust        |
|                 |                 |                 | algorithms or   |
|                 |                 |                 | not             |
+-----------------+-----------------+-----------------+-----------------+
| H2CHEM          | yes/no          | no              | use H2          |
|                 |                 |                 | chemistry or    |
|                 |                 |                 | not             |
+-----------------+-----------------+-----------------+-----------------+
| DISC_VISCOSITY  | yes/no          | no              | apply           |
|                 |                 |                 | artificial      |
|                 |                 |                 | viscosity to    |
|                 |                 |                 | both            |
|                 |                 |                 | approaching and |
|                 |                 |                 | receding        |
|                 |                 |                 | particles and   |
|                 |                 |                 | multiply by     |
|                 |                 |                 | h/rij           |
+-----------------+-----------------+-----------------+-----------------+
| CONST_AV        | yes/no          | no              | use a constant  |
|                 |                 |                 | artificial      |
|                 |                 |                 | viscosity       |
|                 |                 |                 | parameter       |
|                 |                 |                 | instead of the  |
|                 |                 |                 | Morris &        |
|                 |                 |                 | Monaghan switch |
+-----------------+-----------------+-----------------+-----------------+
| CONST_ARTRES    | yes/no          | no              | use a constant  |
|                 |                 |                 | artificial      |
|                 |                 |                 | resistivity     |
|                 |                 |                 | parameter       |
|                 |                 |                 | instead of the  |
|                 |                 |                 | Tricco & Price  |
|                 |                 |                 | switch (MHD     |
|                 |                 |                 | only)           |
+-----------------+-----------------+-----------------+-----------------+
| CURLV           | yes/no          | no              | store curl v    |
|                 |                 |                 | and write it to |
|                 |                 |                 | full dump files |
+-----------------+-----------------+-----------------+-----------------+
| USE_STRAIN_TENS | yes/no          | no              | determines      |
| OR              |                 |                 | whether or not  |
|                 |                 |                 | strain tensor   |
|                 |                 |                 | is stored, and  |
|                 |                 |                 | therefore       |
|                 |                 |                 | whether         |
|                 |                 |                 | physical        |
|                 |                 |                 | viscosity is    |
|                 |                 |                 | done using two  |
|                 |                 |                 | first           |
|                 |                 |                 | derivatives or  |
|                 |                 |                 | two second      |
|                 |                 |                 | derivatives     |
|                 |                 |                 | (see Lodato &   |
|                 |                 |                 | Price 2010)     |
+-----------------+-----------------+-----------------+-----------------+
| DUSTGROWTH      | yes/no          | no              | use dust growth |
|                 |                 |                 | (and/or         |
|                 |                 |                 | fragmentation)  |
|                 |                 |                 | prescription    |
|                 |                 |                 | from Stepinski  |
|                 |                 |                 | & Valageas      |
|                 |                 |                 | (1997) for      |
|                 |                 |                 | two-fluid       |
|                 |                 |                 | algorithm or    |
|                 |                 |                 | not             |
+-----------------+-----------------+-----------------+-----------------+

Examples
--------

For example, to set individual timesteps on the command line:

::

   make IND_TIMESTEPS=yes

or put this in the SETUP block:

::

   ifeq ($(SETUP), isodisc)
       FPPFLAGS= -DDISC_VISCOSITY -DSORT_RADIUS_INIT -DISOTHERMAL
       ...
       IND_TIMESTEPS=yes

or add it as a preprocessor flag:

::

   ifeq ($(SETUP), isodisc)
       FPPFLAGS= -DDISC_VISCOSITY -DSORT_RADIUS_INIT -DISOTHERMAL -DIND_TIMESTEPS
       ...
