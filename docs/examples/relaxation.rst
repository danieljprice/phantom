
Stellar relaxation
==================


When mapping 1D profile from MESA to a 3D code such as Phantom, a relaxation procedure is required to ensure that the particles are in hydrostatic equilibrium.

In phantom this is achieved by a process called ``relax-o-matic``
(`Lau et al., 2022 <https://doi.org/10.1093/mnras/stac049>`_ - **please see its appendix C for full details of how it works**).


To run ``relax-o-matic`` after setting up a star,
just follow the instructions during setup (from ``phantomsetup``) and type ``yes`` when prompted "Relax star automatically during setup?";

Then, modify the ``.setup`` file and fill in the relaxation options, which are pretty self-explanatory::

	relax_star        =           T    ! relax star automatically during setup
	tol_ekin          =   1.000E-07    ! tolerance on ekin/epot to stop relaxation
	tol_dens          =       1.000    ! % error in density to stop relaxation
	maxits            =        1000    ! maximum number of relaxation iterations
	write_rho_to_file =           F    ! write density profile to file

Run the ``phantomsetup`` again to start the auto relaxation process.
You are now good to go!

If interruptted during the relaxation, run ``phantomsetup`` again and it will pick up from where it left off automatically.

If in doubt, you can always run the star in isolation for a few years to see if relaxation procedure worked.

Rotating stars
--------------

Set ``relax = T`` and specify a signed spin period for each gas star in the
``.setup`` file, for example::

    spin_period1 = 24 hr
    spin_period2 = -48 hr

The period accepts the usual time units; a bare number is in code units.
Positive periods rotate counterclockwise as seen from positive z, negative
periods clockwise. The default, ``0 hr``, leaves the non-rotating workflow
unchanged, including old setup files without this option.
For a binary with prescribed individual spins, leave ``corotate = F`` to avoid
adding a second spin contribution when placing the stars on their orbit.

The existing asynchronous particle shifts include centrifugal acceleration
about z. The entropy profile is retained as a function of mass coordinate;
the relaxed rotating star need not have the original spherical density profile.
Convergence and ``relax*.ev`` use the residual shifting kinetic energy, not
the physical rotational kinetic energy. Returned models and relaxation dumps
contain the physical rigid-spin velocities, ready for evolution without damping.

This option currently requires Newtonian gravity without external forces.
It is not a new runtime damping mode. Excessively rapid rotation may prevent
convergence or cause mass shedding. Remove existing ``relax*`` dumps before
changing the requested spin or stellar profile, then check the released star
in isolation before using it in a collision.

Density-weighted stellar particles
----------------------------------

To generate a star from a radial profile with approximately uniform initial
particle spacing and individually varying masses, select the optional
density-weighted hexagonal close-packed distribution::

    particle_lattice1 = weighted_hcp
    hcp_density_exponent1 = 0.
    np1 = 100000
    relax = T
    spin_period1 = 24 hr

Use ``KERNEL=WendlandC4`` for these stellar models, together with
``MPI=no`` and ``IND_TIMESTEPS=no``. Retain ``hfact=1.5`` in the runtime
``.in`` file: too few neighbours can leave unresolved density errors
where particle masses change sharply. Check convergence with resolution
and neighbour count before interpreting mixing at such interfaces.

Each star can have its own requested particle count. The actual count is the
number of accepted lattice sites, not necessarily exactly ``npN``. Select
``weighted_hcp`` for all gas stars in a multiple-star setup. The default
distribution and existing equal-mass setups are unchanged.

The single-exponent construction follows the parent-star prescription in
StarSmasher (revision ``140a5553f31c37c2d0d879728296c9a9afa82146``).
The exponent ``q`` controls number density and mass: ``n`` is proportional
to ``rho**q`` and particle masses to ``rho**(1-q)``. Thus ``q=0`` gives
uniformly spaced initial sites and ``q=1`` gives equal masses. Intermediate
values trade envelope resolution against core resolution. This generates
new particles from the profile; it does not import another code's dump.

The lattice is an initial condition, not a constraint during evolution.
Particle masses and count remain fixed, but spacing and smoothing lengths
change as the gas moves. Runtime smoothing lengths use a separate normalized
neighbour-count kernel with a flat inner core and a cumulative cubic-spline
taper, following StarSmasher's default ``gflag=1``. Mass density is independently
summed with the selected SPH kernel. Pressure work and softened self-gravity
include the corresponding variable-smoothing-length terms. Artificial
viscosity and conduction use a monotone kernel gradient rather than the
pressure correction.

The initial surface prescription uses ``hcp_nnoptN=23``,
``hcp_surface_factorN=1.3``, ``hcp_hfloorN=0`` and
``hcp_surface_bufferN=-1`` (the automatic two-initial-smoothing-length cut,
limited by the 85-percent enclosed-mass radius). A positive surface buffer
overrides that cut. These are generation controls, not the runtime neighbour
constraint; ``hfact`` remains independent. Inadequately resolved point cores
are rejected rather than silently removed.

MESA profile conversion uses the existing Phantom solar units by default.
For a profile using a different convention, set ``hcp_profile_mass_unitN``
(grams) and ``hcp_profile_radius_unitN`` (centimetres) explicitly. Importing a
MESA profile does not select the MESA equation of state.

Continuous masses currently require Newtonian gas-only hydrodynamics,
``MPI=no``, ``APR=no`` and ``IND_TIMESTEPS=no``; OpenMP is supported.
MHD, radiation transport, dust, sinks, injection and gravitational-wave
diagnostics are not supported in this mode. Radiation pressure in the EOS
is distinct from radiation transport. The existing asynchronous automatic relaxer
supports both static and rigidly rotating models; there is no additional
runtime damping mode. Verify each released stellar model in isolation.

Full dumps store masses, smoothing lengths and mass density in full precision
and support restart. Small dumps retain masses but have reduced precision
geometry and are not restart files. Legacy dumps without individual masses
continue to use the existing mass-per-type representation.
Do not restart variable-mass dumps with older Phantom binaries that lack
individual-mass support.
