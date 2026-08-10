Driven supersonic turbulence
============================

The turbulent driving module implements an Ornstein–Uhlenbeck stirring pattern
(Federrath; see `Price & Federrath 2010 <https://ui.adsabs.harvard.edu/abs/2010MNRAS.405.2135P>`__,
`Tricco, Price & Federrath 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.1179T>`__).
The algorithms and tests are described in the
`Phantom code paper <http://ui.adsabs.harvard.edu/abs/2018PASA...35...31P>`__ (Section 2.4);
this page explains how to **run** turbulence calculations.

Which SETUP to use
------------------

+------------------+----------------------------------------------------------+
| ``writemake.sh`` | Use case                                                 |
+==================+==========================================================+
| ``turb``         | Hydro or MHD turbulence in a periodic box (main recipe)  |
+------------------+----------------------------------------------------------+
| ``dustyturb``    | Same as ``turb`` with dust-as-mixture enabled             |
+------------------+----------------------------------------------------------+
| ``turbdrive``    | Uniform periodic box + driving; used in                   |
|                  | :doc:`density` (zero-step density on external particles) |
+------------------+----------------------------------------------------------+

All three compile with ``PERIODIC=yes`` and ``DRIVING=yes`` (stirring compiled in).

Basic hydro turbulence run
--------------------------

Create a run directory outside the code tree::

   mkdir -p ~/runs/turb
   cd ~/runs/turb
   ~/phantom/scripts/writemake.sh turb > Makefile
   make setup
   make

Run the setup utility (interactive prompts, then edit ``turb.setup`` if needed)::

   ./phantomsetup turb

Default geometry is a unit periodic cube with a cubic particle lattice
(``ilattice=1``; close-packed ``ilattice=2`` is closer to a relaxed distribution).
Default units mimic a molecular cloud (3 pc box, low density, 0.2 km/s sound speed);
see ``setup_turb.f90`` for defaults.

Run the simulation::

   ./phantom turb.in

MHD turbulence
--------------

The ``turb`` preset does **not** enable MHD by default. Rebuild with::

   make MHD=yes
   make setup MHD=yes
   make MHD=yes

Then run ``phantomsetup`` again. Set ``Bz_0`` in ``turb.setup`` for a uniform field in
:math:`z`. The setup prints ``MHD turbulence w/uniform field in z-direction`` when
MHD is active.

Dusty turbulence
----------------

Use the pre-cooked dusty preset::

   ~/phantom/scripts/writemake.sh dustyturb > Makefile
   make setup
   make
   ./phantomsetup turb

This uses ``setup_turb.f90`` with ``DUST=yes`` and the dust-as-mixture method (see
:doc:`/physics/dust`). Configure grain bins in ``turb.setup`` as for other dusty setups.

Stirring parameters (``.in`` file)
----------------------------------

Driving is controlled at runtime when ``DRIVING`` was enabled at compile time.
Key options (see also the code paper Table of parameters, Section 2.4):

+------------------+------------------------------------------+
| Variable         | Role                                     |
+==================+==========================================+
| ``istir``        | Turn stirring on (1) or off (0)           |
+------------------+------------------------------------------+
| ``st_energy``    | Energy input per mode ($E_{\rm m}$)      |
+------------------+------------------------------------------+
| ``st_decay``     | Correlation / decay time ($t_{\rm decay}$)|
+------------------+------------------------------------------+
| ``st_solweight`` | Solenoidal weight ($w$; 1 = solenoidal)  |
+------------------+------------------------------------------+
| ``st_stirmin``   | Minimum driving wavenumber ($k_{\rm min}$) |
+------------------+------------------------------------------+
| ``st_stirmax``   | Maximum driving wavenumber ($k_{\rm max}$) |
+------------------+------------------------------------------+
| ``st_dtfreq``    | Interval between stirring updates        |
+------------------+------------------------------------------+
| ``st_amplfac``   | Overall amplitude multiplier             |
+------------------+------------------------------------------+
| ``st_seed``      | Random number seed for the forcing field |
+------------------+------------------------------------------+

A ``forcing.dat`` file may be written/read to continue an identical forcing pattern
across restarts.

Suggested starting values are written into a fresh ``turb.in`` by ``phantomsetup``
(e.g. ``tmax``, ``dtmax``, ``nfulldump`` for several outputs per crossing time).
Tune ``st_energy`` and ``st_decay`` together: the module sets
``st_OUvar = sqrt(st_energy/st_decay)``.

Resolution and timesteps
------------------------

- Use ``hfact`` and kernel choice as in :doc:`/user-guide/infile` (the paper recommends
  quintic kernels for some dusty tests; default ``turb`` uses the setup default kernel).
- ``C_cour`` and ``C_force`` default to 0.3 and 0.25 (paper test suite values).
- Individual timesteps are enabled in the ``turb`` preset (``IND_TIMESTEPS=yes``).

Analysis and comparison
-----------------------

- Visualise dumps with `splash <https://github.com/danieljprice/splash>`_.
- Measure :math:`\nabla\cdot\mathbf{v}` and :math:`\nabla\times\mathbf{v}` with
  ``phantom2divv`` (see :doc:`/user-guide/utils`).
- For standardised MHD benchmark problems from the code paper (e.g. Orszag–Tang),
  see the external `phantom-examples <https://github.com/phantomSPH/phantom-examples>`__
  repository (also mentioned in :doc:`/getting-started/flatiron`).

Further reading
---------------

- :doc:`/user-guide/config` — compile-time flags (``PERIODIC``, ``DRIVING``, ``MHD``)
- :doc:`/user-guide/setups` — full SETUP list
- :doc:`density` — periodic box without the dedicated ``turb`` setup file
- Code paper Section 2.4 (physics); Section 5.1 (supersonic turbulence application)
