Moddump recipes
===============

This page lists common tasks that modify an existing Phantom dump file using
built-in modules in ``src/utils/``. For the generic moddump workflow (``.in``
file handling, writing your own module), see :doc:`moddump`.

Compile a specific module from your run directory (with a local Makefile from
``~/phantom/scripts/writemake.sh``)::

   make moddump MODFILE=moddump_addflyby.f90

Run in the **same directory as the input ``.in`` file** (see :doc:`moddump`)::

   ./phantommoddump dump_in dump_out [time]

Add a flyby to an evolved disc
------------------------------

Use this when you already have a relaxed or partially evolved disc and want to
add an **unbound** perturber on a **parabolic orbit** (eccentricity = 1),
without rebuilding the initial conditions from scratch.

This differs from a :doc:`flyby set up at t=0 </examples/disc>` via
``phantomsetup``, which creates both sinks and the disc together.

**Module:** ``moddump_addflyby.f90``

::

   cd ~/runs/mydisc
   make moddump MODFILE=moddump_addflyby.f90
   ./phantommoddump relax_00150 flyby_00000 0.

You will be prompted for:

- mass and accretion radius of the perturber
- distance of minimum approach ``dma``
- initial distance in units of ``dma`` (``n0``; typically :math:`\gg 1`)
- orbit orientation (ascending node angle, inclination)
- suggested time between dumps as a fraction of the flyby time

The module reads the primary sink mass from the dump, recentres the system, and
prints a suggested ``dtmax`` (in code units) for the output ``flyby.in``. Edit
``flyby.in`` as needed, then run::

   ./phantom flyby.in

Orbital element conventions are described in :doc:`/physics/orbits`.

Add dust after gas relaxation
-----------------------------

**Module:** ``moddump_dustadd.f90``

Used in the :doc:`dust settling test </examples/dustsettle>`: relax a gas-only
disc, then add dust particles before the science run.

::

   make moddump MODFILE=moddump_dustadd.f90
   ./phantommoddump relax_00150 disc_00000 0.

Extend a disc in radius
-----------------------

**Module:** ``moddump_extenddisc.f90``

Extends an existing disc outward, assuming the same surface density profile.

::

   make moddump MODFILE=moddump_extenddisc.f90
   ./phantommoddump disc_00100 disc_extended_00000 0.

Add sinks or planets
--------------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Purpose
   * - ``moddump_addsink``
     - Add a sink particle (interactive)
   * - ``moddump_addplanets``
     - Add embedded planets in a disc
   * - ``moddump_sinkbinary``
     - Add a binary sink pair
   * - ``moddump_binary``
     - Convert single star to binary
   * - ``moddump_binarystar``
     - Split relaxed star into binary components

Trim or edit the particle distribution
--------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Module
     - Purpose
   * - ``moddump_removeparticles_radius``
     - Remove particles inside or outside a sphere
   * - ``moddump_removeparticles_cylinder``
     - Remove particles outside a cylinder
   * - ``moddump_recenter`` / ``moddump_CoM``
     - Recentre on centre of mass
   * - ``moddump_faceon``
     - Realign disc face-on, move origin
   * - ``moddump_rescale``
     - Change code units of the dump
   * - ``moddump_changemass``
     - Change particle masses
   * - ``moddump_disc``
     - Disc warp and magnetic field

Change resolution
-----------------

Use the standalone binaries (or the moddump wrappers that call the same logic)::

   make splitpart
   ./splitpart dump_00000 dump_hi_00000

::

   make mergepart
   ./mergepart dump_00000 dump_lo_00000

Import external dumps
---------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Module
     - Purpose
   * - ``moddump_sphNG2phantom``
     - Import sphNG dump
   * - ``moddump_sphNG2phantom_disc``
     - sphNG disc with sinks
   * - ``moddump_sphNG2phantom_addBfield``
     - Import and add B-field

Inject particles from another simulation
----------------------------------------

See :doc:`inject_sim` for the ``inject_sim`` moddump module and the associated
``.in`` options (``start_dump``, ``final_dump``, ``r_inject``).

Other useful moddumps
---------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Module
     - Purpose
   * - ``moddump_sink``
     - Edit sink properties interactively
   * - ``moddump_sink2gas``
     - Replace sinks with resolved gas
   * - ``moddump_growthtomultigrain``
     - Convert dust-growth dump for MCFOST
   * - ``moddump_LTE_to_rad`` / ``moddump_rad_to_LTE``
     - Switch radiation mode
   * - ``moddump_infall``
     - Add infalling cloud
   * - ``moddump_rotate``
     - Solid-body rotation on a sphere

A full list of modules is in ``src/utils/moddump_*.f90``. To compile every
moddump for testing::

   make allmoddumps

Post-processing with phantomanalysis
------------------------------------

To **analyse** dumps without modifying them, use ``phantomanalysis`` with a
built-in ``analysis_*`` module (see :doc:`analysis`). Common disc examples:

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Module
     - Purpose
   * - ``analysis_disc``
     - Azimuthally averaged disc profiles
   * - ``analysis_dustydisc``
     - Disc analysis for ``dustydisc`` setups
   * - ``analysis_disc_planet``
     - Disc with embedded planet
   * - ``analysis_disc_stresses``
     - Stresses and alpha

::

   make analysis ANALYSIS=analysis_disc.f90
   ./phantomanalysis dump_00100

Appendix: built-in module index
-------------------------------

All moddump plugins live in ``src/utils/`` as ``moddump_*.f90``. All analysis
plugins live as ``analysis_*.f90``. Read the module header in each file for
details. To compile a specific analysis module::

   make analysis ANALYSIS=analysis_dustydisc.f90

**moddump modules:** ``moddump_addflyby``, ``moddump_addplanets``, ``moddump_addsink``,
``moddump_binary``, ``moddump_binarystar``, ``moddump_changemass``, ``moddump_CoM``,
``moddump_disc``, ``moddump_dustadd``, ``moddump_dustsinkfrac``, ``moddump_extenddisc``,
``moddump_faceon``, ``moddump_growthtomultigrain``, ``moddump_infall``,
``moddump_LTE_to_rad``, ``moddump_mergepart``, ``moddump_messupSPH``,
``moddump_perturbgas``, ``moddump_polytrope``, ``moddump_radiotde``,
``moddump_rad_to_LTE``, ``moddump_recalcuT``, ``moddump_recenter``,
``moddump_removeparticles_cylinder``, ``moddump_removeparticles_radius``,
``moddump_rescale``, ``moddump_rotate``, ``moddump_sink``, ``moddump_sink2gas``,
``moddump_sinkbinary``, ``moddump_splitpart``, ``moddump_sphNG2phantom``,
``moddump_sphNG2phantom_addBfield``, ``moddump_sphNG2phantom_disc``,
``moddump_taylorgreen``, ``moddump_tdesink``, ``moddump_temp``, ``moddump_tidal``,
``moddump_torus``, plus ``moddump_default`` (template).

**analysis modules (selection):** ``analysis_disc``, ``analysis_dustydisc``,
``analysis_disc_planet``, ``analysis_disc_stresses``, ``analysis_disc_mag``,
``analysis_disc_MFlow``, ``analysis_dustmass``, ``analysis_energies``,
``analysis_angmom``, ``analysis_ptmass``, ``analysis_structurefn``,
``analysis_clumpfind``, ``analysis_mcfost``, ``analysis_tde`` — and others in
``src/utils/analysis_*.f90``.
