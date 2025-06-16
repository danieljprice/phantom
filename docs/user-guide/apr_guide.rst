Splitting/merging particles
===========================

Splitting or merging particles can be used to increase (or decrease) the resolution of your simulation. The
simplest form of this is to split (merge) a whole simulation into a higher (lower) resolution simulation. A good use
case for this is to evolve an accretion disc for a long time at low resolution to settle initial transients, 
then split to high resolution for a smaller number of orbits.

More complex splitting/merging can be achieved by using Adaptive Particle Refinement (APR) (see :ref:`apr`).

Splitting a whole simulation
----------------------------
A global split/merge can be achieved using the ``splitpart`` and ``mergepart`` utilities. These are actually
`moddump <moddump.rst>`_ modules but just compiled with different executable names. You can split a simulation into a higher
resolution simulation as follows::

  make splitpart
  ./splitpart dump_01000 split_01000 --nchild=13

This will split each particle in the simulation into 13, placed in a regular lattice. The ``--nchild`` option
specifies the number of children to split into. If you choose a number other than 13 the children will be
be placed randomly rather than on a regular lattice.

Merging a whole simulation
--------------------------
You can merge particles as follows::

  make mergepart
  ./mergepart dump_01000 merge_01000 --nchild=2

which will produce a snapshot with npart/nchild particles. The merging procedure uses the tree structure
of the particles to merge them into the correct number of particles, as described in Nealon & Price (2025).

.. _apr:

Running phantom with APR
========================

Adaptive Particle Refinement (APR) allows you to arbitrarily set regions of your simulation to have different local resolutions.
A description of the method and our implementation can be found in Nealon & Price (2025) and
you should cite this work if you are using APR.

The rules of this method are:
 1.	A mass factor of 2 is allowed between adjacent refinement levels (e.g. one parent 0> two children)
 2.	Regions are spherical but can be nested together like a layered onion
 3.	You cannot have both refinement and derefinement in the same simulation, e.g. the base refinement level has to be either the maximum or minimum
 4.	Try to ensure that particles have several sound crossing times between subsequent split/merge procedures to reduce noise


Compiling Phantom with APR
--------------------------
Use your usual setup routine but with APR=yes. For example::

     mkdir mycalc
     cd mycalc
     ~/phantom/scripts/writemake.sh disc > Makefile
     export APR=yes
     make setup
     make

Initial conditions with APR
-------------------------
APR does not affect any setups except SETUP=star (but for the relax procedure only). When you check your .in file it will have apr options at the bottom:

::

   apr_max =           3    ! number of additional refinement levels (3 -> 2x resolution)

This option sets how many *extra* levels of refinement you want. Each level corresponds to a factor of 2 in mass.
To increase the spatial resolution by a factor of 2 you will need to use 3 levels here::

  ref_dir =           1    ! increase (1) or decrease (-1) resolution

This chooses whether you are refining or derefining the simulation. If you are refining then the base resolution level will be 0.
If you are derefining the base level will be apr_max::

  apr_type =           2    ! 1: static, 2: moving sink, 3: create clumps

Here you choose what kind of region you want. Current options include:
 1.	A position fixed in space
 2.	Tracking a particular sink particle
 3.	Tracking a gravitationally bound clump (under development).
Depending on what you choose here you will get additional options to describe the properties of the region you selected.
You may need to re-run to get the right options if you alter apr_type. To add your own new region you can edit the apr_region.f90 file.
Note for now that we only allow spherical regions::

  apr_rad =         5.    ! radius of innermost region

This chooses the radius of the core of the refinement zone, where the highest/lowest refinement will be. If you think about
your refinement zone as an onion, this is the core of the onion::

  apr_drad =         10.    ! size of step to next region

This chooses the step width of the nested regions - the width of the shells as you step up your resolution.

Running with APR
--------------------
When APR is implemented and being used Phantom prints out the following statement in the log file::

    Adapative particle refinement is ON

Additionally, because the particle numbers are changing each step you should see the number of
particles being updated during the steps e.g.::

> step 2 / 16 t = 20.92159 dt = 0.072 moved 502 in 0.058 cpu-s < | np = 37902 |
> step 4 / 16 t = 21.06588 dt = 0.072 moved 1792 in 0.070 cpu-s <
> step 6 / 16 t = 21.21017 dt = 0.072 moved 319 in 0.058 cpu-s <
> step 8 / 16 t = 21.35445 dt = 0.072 moved 6175 in 0.097 cpu-s <
> step 10 / 16 t = 21.49874 dt = 0.072 moved 442 in 0.057 cpu-s < | np = 37901 |
> step 12 / 16 t = 21.64303 dt = 0.072 moved 1283 in 0.064 cpu-s <
> step 14 / 16 t = 21.78732 dt = 0.072 moved 476 in 0.058 cpu-s < | np = 37900 |
> step 16 / 16 t = 21.93160 dt = 0.072 moved 37860 in 0.31 cpu-s <

Plotting with APR
--------------------
APR is natively read by both splash and sarracen. The easiest way to check the method is working exactly
as you expect is to scatter plot (not render!) the mass of the particles in your simulation. This will
show you where the refinement zone is and you can confirm the geometry, its evolution as well as the
refinement direction.

Analysis with APR
--------------------
No analysis files that ship with Phantom have been updated to accommodate APR. To do this yourself, any
time you define a particle mass from the massoftype array you will need to edit it to read::

  if (use_apr) then
     pmassi = aprmassoftype(iamtypei,apr_level(i))
  else
     pmassi = massoftype(iamtypei)
  endif

This relies on the apr_level, aprmassoftype and use_apr which can be included with::

  use dim,  only::use_apr
  use part, only::apr_level,aprmassoftype

Note that apr_level is integer(kind=1).
