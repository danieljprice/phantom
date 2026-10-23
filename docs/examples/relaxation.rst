
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
