Non-ideal magnetohydrodynamics
==============================

Non-ideal MHD (Ohmic resistivity, ambipolar diffusion, Hall effect) is implemented
via the NICIL library (Wurster; see also the
`Phantom code paper <http://ui.adsabs.harvard.edu/abs/2018PASA...35...31P>`__, Section 2.11).
This page describes how to **compile and run** non-ideal calculations; it does not
repeat the governing equations or test-problem definitions from the paper.

Requirements
------------

Non-ideal terms require **both** MHD and NICIL at compile time::

   make MHD=yes NONIDEALMHD=yes
   make setup MHD=yes NONIDEALMHD=yes
   make MHD=yes NONIDEALMHD=yes

Check the build summary line includes ``non-ideal`` (see :doc:`/user-guide/config`).

Useful SETUP presets
--------------------

+------------------+----------------------------------------------------------+
| ``writemake.sh`` | Use case                                                 |
+==================+==========================================================+
| ``nimhdshock``   | Standing / C-shock tubes; constant or NICIL coefficients |
+------------------+----------------------------------------------------------+
| ``jetnimhd``     | Collapsing magnetized core; self-gravity + non-ideal MHD |
+------------------+----------------------------------------------------------+
| ``testnimhd``    | Unit tests and wave-damping tests (developers)           |
+------------------+----------------------------------------------------------+

All of these use ``KERNEL=WendlandC4`` in ``Makefile_setups`` (recommended for
non-ideal shock tests in the paper).

Shock-tube test (``nimhdshock``)
-------------------------------

Good first sanity check; corresponds to paper Section 4.8 / 2.11 tests::

   mkdir -p ~/runs/nimhdshock
   cd ~/runs/nimhdshock
   ~/phantom/scripts/writemake.sh nimhdshock > Makefile
   make setup
   make
   ./phantomsetup shock

Edit ``shock.setup`` to choose the case (interactive setup), which non-ideal terms
are active (``use_ohm``, ``use_hall``, ``use_ambi``), and coefficient values
(``C_OR``, ``C_HE``, ``C_AD`` in the setup file for constant-coefficient runs).

The generated ``shock.in`` contains a ``# options controlling non-ideal MHD`` block
when NICIL is compiled in. Important runtime switches include:

- ``use_ohm``, ``use_hall``, ``use_ambi`` — which terms to include
- ``eta_constant`` — use fixed coefficients vs. NICIL ionisation calculation
- ``eta_const_type`` — how constant coefficients are specified (see module
  ``nicil_supplement`` in ``src/main/nicil_supplement.f90``)
- ``Cdt_diff``, ``Cdt_hall`` — timestep limiter factors for diffusive and Hall terms

Run::

   ./phantom shock.in

Compare with the paper figures or run a lower-resolution copy via the test suite
(``make SETUP=test phantomtest && ./bin/phantomtest nonidealmhd``; see
:doc:`/developer-guide/testing`).

Star formation with non-ideal MHD (``jetnimhd``)
------------------------------------------------

Preset for magnetized, self-gravitating, periodic-box star-formation style setups
(Wurster, Price & Bate 2016, 2017)::

   ~/phantom/scripts/writemake.sh jetnimhd > Makefile
   make setup
   make
   ./phantomsetup sphereinbox

``setup_sphereinbox.f90`` places a dense sphere in a lower-density medium (or a
uniform box if ``density_contrast=1`` for decaying turbulence). Edit ``sphereinbox.setup``
for magnetic field strength (``Bzero`` or mass-to-flux ratio), turbulence
(``rms_mach``), Bonnor–Ebert options, and sinks.

Non-ideal coefficients are set in the ``.in`` file (NICIL block) after setup.
Self-gravity and super-timestepping for non-ideal terms are handled automatically
in this preset.

Ideal MHD benchmarks (comparison)
---------------------------------

The code paper ideal MHD tests (e.g. Orszag–Tang vortex) use **ideal** MHD, not
``NONIDEALMHD``. Reproduce them from the
`phantom-examples <https://github.com/phantomSPH/phantom-examples>`_ repository,
for example::

   git clone https://github.com/phantomSPH/phantom-examples
   cd phantom-examples/mhd/orstang

as described in :doc:`/getting-started/flatiron`.

Runtime output
--------------

Summary output may list limiting timesteps from non-ideal terms
(``dtohmic``, ``dthall``, ``dtambipolar``). Dump headers can store per-particle
``eta_nimhd`` arrays when configured.

Developer checks
----------------

- Full non-ideal unit tests: ``make SETUP=test`` then ``./bin/phantomtest nonidealmhd``
- CI runs these on every pull request (:doc:`/developer-guide/testing`)

Further reading
---------------

- :doc:`/user-guide/infile` — general runtime options
- :doc:`/user-guide/config` — ``MHD`` compile flag
- :doc:`turbulence` — driven turbulence (often combined with MHD in molecular clouds)
- :doc:`/examples/star` and :doc:`/examples/relaxation` — stellar initial conditions
- Code paper Sections 2.11 and 4.8; Wurster et al. NICIL papers (cited in releasenotes)
