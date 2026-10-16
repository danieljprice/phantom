Density-weighted stellar models
==============================

``particle_lattice = weighted_hcp`` generates a star directly from a radial
profile using antipodal hexagonal close-packed sites and a natural cubic
radial map. It follows the single-exponent construction in
`StarSmasher <https://github.com/jalombar/starsmasher/tree/140a5553f31c37c2d0d879728296c9a9afa82146>`_.
No particle snapshot from another code is required.

For ``hcp_density_exponent = q``, the initial number density scales as
``rho**q`` and the particle masses as ``rho**(1-q)``. Thus ``q=0`` gives
uniform spatial sampling with smaller masses in the envelope, whereas
``q=1`` gives equal masses. The accepted lattice count approximates ``np``;
particle masses and count then remain fixed during evolution. Uniform
spacing is an initial condition, not a constraint on the flow.

Build with ``SETUP=star MPI=no APR=no IND_TIMESTEPS=no DOUBLEPRECISION=yes``.
The supported runtime is Newtonian gas hydrodynamics with shared timesteps
and optional OpenMP. APR, MPI, individual timesteps, dust, MHD, radiation,
injection, sinks and gravitational-wave diagnostics are rejected for
continuous masses. This mode changes neither the EOS nor its data tables.

Generation
----------

The public ``data/star_data_files/weighted_hcp_demo.setup`` example uses an
analytic n=1.5 polytrope and an ideal-gas EOS. From the repository root::

   make SYSTEM=gfortran SETUP=star MPI=no APR=no IND_TIMESTEPS=no \
        OPENMP=yes DOUBLEPRECISION=yes phantom setup phantomtest
   mkdir hcp_demo
   cp data/star_data_files/weighted_hcp_demo.setup hcp_demo/demo.setup
   cd hcp_demo
   ../bin/phantomsetup demo

For other profiles, select the usual stellar profile options and add::

   particle_lattice1 = weighted_hcp
   hcp_density_exponent1 = 0.
   np1 = 100000
   relax = F

Use this lattice for every gas star in a multi-star setup. The optional
``hcp_nnopt`` (23), ``hcp_surface_factor`` (1.3) and ``hcp_hfloor`` (0)
set the initial support and smoothing estimates, not the runtime neighbour
constraint. The support edge is the larger of the 85-percent enclosed-mass
radius and ``R - 2*h_c``, where
``h_c = 0.5*R*(hcp_surface_factor*hcp_nnopt/np)**(1/3)``.
Set a positive ``hcp_surface_buffer`` to replace the factor 2; ``-1``
retains the default. An unresolved point core is rejected rather than omitted.

MESA profile conversion normally uses Phantom's solar mass and radius.
``hcp_profile_mass_unit`` and ``hcp_profile_radius_unit`` override the
profile's solar-unit convention in grams and centimetres, independently
of the simulation units. For older profiles these may be ``1.9892e33``
and ``6.9598e10``. These options also take the star-number suffix.

Relaxation and spin
-------------------

The setup relax-o-matic assumes equal masses and is not used here.
Instead generate with ``relax=F`` and evolve the isolated star with::

   idamp = 4
   damp_time_s = 3200.
   damp_stop_s = 32000.
   damp_spin_x_s = 0.
   damp_spin_y_s = 0.
   damp_spin_z_s = 0.

The drag acceleration is
``-(v - V_cm - omega cross (r - R_cm))/tau``. Times are physical seconds;
spin components are radians per second. A zero stop time keeps damping on.
Zero spin gives static relaxation. Nonzero spin drives rigid rotation in
the inertial frame, without adding a separate centrifugal acceleration.
The mass-weighted centre and bulk velocity are updated at each force call.
Damping removes energy without depositing heat and can exchange angular
momentum with the imposed rotation. It is intended for one isolated star.

To assign a z-axis spin to an existing static dump, build ``moddump`` with
``MODFILE=moddump_rotate.f90`` and set ``omega`` in its ``.moddump`` file
or through ``--omega=...``. This utility uses code-time units, retains bulk
velocity, and replaces residual velocities by rigid spin about the mass
centre. The result still needs rotating relaxation. Finally set ``idamp=0``
and check free stability, radii, residual motion and conservation before
using either star in an encounter.

Encounter and output
--------------------

Build ``moddump`` with ``MODFILE=moddump_stellar_encounter.f90`` to combine
two isolated, compatible full dumps on an incoming hyperbolic orbit::

   make SYSTEM=gfortran SETUP=star MPI=no APR=no IND_TIMESTEPS=no \
        OPENMP=yes DOUBLEPRECISION=yes MODFILE=moddump_stellar_encounter.f90 \
        MODDUMPBIN=phantom_stellar_encounter moddump
   phantom_stellar_encounter primary_00000 encounter_00000 0 --maxp=200000 \
        --secondary=secondary_00000 --separation=4 --pericentre=1.25 --ecc=1.6679

On the first invocation the utility writes ``encounter.moddump``; rerun
after checking its options. Distances are in common code units and
``arg_peri`` is in degrees (default 180). The stars must initially be
separate and ``ecc>1``. Units, EOS and smoothing factor must match;
also check composition and EOS parameters in the two run inputs.
Only bulk positions and velocities change: masses, thermal states and
internal spin are retained. Primary IDs are retained and secondary IDs
are offset above the largest primary ID. Damping is disabled on assembly.

Full dumps preserve individual masses, density and smoothing lengths at
the build's precision and can be restarted. Small dumps retain real8 masses
but single-precision geometry and density; they are not restart files.
The Python dump reader also reads the explicit real8 mass block.
An explicit ``--maxp=N`` can reserve extra tree capacity for a centrally
concentrated star without adding particles; repeat it on restart.

The smoothing constraint uses number density, while mass density is the
mass-weighted kernel sum. Pressure work includes the corresponding grad-h
correction. Viscosity and artificial thermal conduction use monotone kernel
gradients instead, conserving pairwise exchange without reversing diffusion.
``C_thermal`` (default 0.05) limits fractional internal-energy losses via
the timestep, without clipping the conservative energy derivative.
Choose ``alphau`` explicitly when comparing thermal treatments: setting
it to zero disables artificial conduction, not viscous shock heating.

``phantomtest hcp damp part dump directsum`` covers construction, mass
transport and restart, rigid-spin damping, kernel density and pairwise
momentum, torque, energy and dissipation checks. These do not establish
resolution convergence or stability for arbitrary stellar profiles or EOS.

References
----------

Gaburov, Lombardi and Portegies Zwart (2010), MNRAS 402, 105--126,
appendix A, describes the non-equal-mass number-density formulation.
The generator is attributed to the pinned StarSmasher construction above;
the isolated spin drag follows its ``nrelax=1`` velocity prescription,
not its separate energy-conserving relaxation procedure.
