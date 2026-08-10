Dust
=====

Dust can be modelled in phantom in two ways:
- as an evolving dust fraction on a single set of mixture particles (Laibe & Price 2014a,b,c, Price & Laibe 2015, Hutchison et al. 2018)
- as a separate set of particles (Laibe & Price 2012a,b; Price & Laibe 2020; Mentiplay et al. 2020)

Dust-as-mixture is best for well coupled dust species (typically < 1 mm grain sizes) and 
with the Ballabio flux limiter ON is fast in terms of computational cost, but is not accurate for highly decoupled dust species.
Dust-as-particles is accurate for both coupled and decoupled dust species (typically > 1 mm grain sizes), but is slow for small grain sizes.

Dust method (``dust_method``)
-----------------------------

Can be chosen at setup time in ``disc.setup`` (or via interactive ``phantomsetup``).

.. list-table::
   :header-rows: 1
   :widths: 8 22 35 35

   * - Value
     - Name
     - When to use
     - Caveats
   * - 1
     - One fluid (dust fraction on gas particles)
     - Many grain sizes as a mixture; fast planet-disc models
     - ``ilimitdustflux`` in ``.in`` can limit spurious diffusion; poor for very large, strongly decoupled grains
   * - 2
     - Two fluid (separate dust particle sets)
     - Drag, back-reaction, settling, streaming tests
     - ``drag_implicit`` available; highly coupled species may limit dt
   * - 3
     - Hybrid (small = mixture, large = particles)
     - Intended to combine methods
     - **Not widely tested:** no reason why it shouldn't work, but has not been tested thoroughly

Choose method 1 or 2 unless you are actively developing or testing the hybrid capability.

Grain size representation
-------------------------

**Single grain size** (``ndusttypesinp = 1`` or one large-grain species in
two-fluid mode): simplest setup; appropriate for settling tests, drag
experiments, and problems where a single representative size is enough.

**Multigrain mixture** (``dust_method = 1``, multiple bins): log-spaced sizes
between ``smincgs`` and ``smaxcgs`` with power-law index ``sindex`` (e.g. MRN
= 3.5). Mass is distributed across bins automatically. More bins increase
cost and complexity; use enough bins to resolve the physics you care about,
not necessarily the full MRN range at high resolution.

**Dust growth** (``SETUP=growingdisc``): grain sizes evolve during the run.
See :doc:`/examples/dustgrowth`.

**After setup:** add dust to a relaxed gas dump with :doc:`/user-guide/moddump-recipes`
(``moddump_dustadd``), or convert growth output with ``moddump_growthtomultigrain``.

Spatial dust distribution
-------------------------

Controlled in ``disc.setup`` via ``isetdust``:

- **0** — dust density follows the gas profile
- **1** — custom dust surface density and scale height
- **2** — same structure as gas but independent inner/outer radii
  (``R_indust``, ``R_outdust``)

For one-fluid simulations, the runtime option ``ilimitdustflux`` in the ``.in``
file applies the Ballabio et al. (2018) flux limiter to stop short timesteps in
 the outer disc from slowing the calculation down. It is only available for ``dust_method = 1``. Set
it in ``disc.in`` after setup, for example::

   ilimitdustflux = T

How to speed up the dust-as-particles method with implicit drag
---------------------------------------------------------------

To use implicit drag, you must compile with individual timesteps disabled, e.g.::

   make IND_TIMESTEPS=no
   make setup IND_TIMESTEPS=no
   make

Then set in the ``.in`` file (written automatically for two-fluid setups)::

   drag_implicit = T

See the ``dust`` module runtime parameters in ``src/main/dust.f90`` and setup-time
dust options in :doc:`/api/setup/set_dust_options`.

Recommended compile presets
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - ``writemake.sh``
     - Typical use
   * - ``disc``
     - Gas disc, isothermal
   * - ``dustydisc``
     - Disc + dust (default analysis module)
   * - ``growingdisc``
     - Disc + dust growth
   * - ``dustyisosgdisc``
     - Self-gravitating dusty disc
   * - ``isosgdisc``
     - Self-gravitating gas disc

Full list: :doc:`/user-guide/setups` and :doc:`/user-guide/setups-list`.

Further reading
---------------

- :doc:`/examples/disc` — interactive setup walkthrough
- :doc:`/examples/dustsettle` and :doc:`/examples/dustgrowth`
- :doc:`/user-guide/moddump-recipes` — modify dumps (flyby, dustadd, extend disc)
- :doc:`/faq` — common questions on viscosity and dust
