.. _set_orbit:

Setting up orbits in phantom
============================

The ``set_orbit`` module provides a generic interface for setting up two-body orbits using various input parameter sets. It is used by many setup routines (e.g. binary, hierarchical systems) to define orbital initial conditions.

The options available in the ``.setup`` file depend on the chosen input type (``itype``).

Input Types
-----------

The typical block in the .setup file looks like this::

   # orbit 
           itype_binary =           2    ! orbital elements (0=aeiOwf,1=flyby,2=Orbit Reconstructor,3=dx,dv,4=x,v)
              binary_dx =     346.000    ! observed dx at t=tmax (code units or e.g. 1 au)
              binary_dy =    -242.000    ! observed dy at t=tmax (code units or e.g. 1 au)
              binary_dz =    -366.667    ! [guessed] dz at t=tmax (code units or e.g. 1 au)
             binary_dvx =       0.200    ! [guessed] dvx at t=tmax (code units or e.g. 1 km/s)
             binary_dvy =      -0.040    ! [guessed] dvy at t=tmax (code units or e.g. 1 km/s)
             binary_dvz =       0.000    ! observed dvz at t=tmax (code units or e.g. 1 km/s)
               binary_d =        1200    ! separation at t=0 if unbound


The control parameter is ``itype`` (or with a prefix, e.g., ``itype_binary``):

*   **0**: Standard Keplerian elements (a, e, i, Omega, w, f)
*   **1**: Flyby parameters (rp, d, e, i, Omega, w)
*   **2**: Orbit Reconstructor (reconstruct orbital parameters to match observed dx,dv at t=tmax)
*   **3**: Relative position and velocity vectors (dx, dv)
*   **4**: Absolute position and velocity vectors (x1, v1, x2, v2)

Detailed Options
----------------

Type 0: Keplerian Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the default mode.

*   ``a``: Semi-major axis. You can specify units (e.g., ``100 au``).
    *   If ``e < 1``, you can alternatively specify the **orbital period** (e.g., ``1000 yrs``).
    *   If ``e > 1`` (hyperbolic), ``a`` should be negative (convention ``a = -GM/2E``).
*   ``e``: Eccentricity (>= 0).
*   ``i``: Inclination in degrees.
*   ``O``: Longitude of ascending node (Omega) in degrees.
*   ``w``: Argument of periapsis (little omega) in degrees.
*   ``f``: Initial true anomaly in degrees (0 = periastron, 180 = apastron).

Type 1: Flyby Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

Useful for setting up unbound encounters.

*   ``rp``: Pericentre distance (closest approach). Can specify units (e.g. ``10 au``).
*   ``d``: Initial separation distance. Can specify units (e.g. ``100 au``).
*   ``e``: Eccentricity.
*   ``i``: Inclination in degrees.
*   ``O``: Longitude of ascending node in degrees.
*   ``w``: Argument of periapsis in degrees.

Type 2: Orbit Reconstructor
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This mode is designed to find an orbit that matches a specific observed separation and velocity at a future time. This is useful for reconstructing orbits of observed systems where only projected quantities are known.

*   ``dx``, ``dy``, ``dz``: Separation vector at the time of observation. ``dz`` is typically a guessed parameter if deprojection is needed.
*   ``dvx``, ``dvy``, ``dvz``: Relative velocity vector at the time of observation. ``dvx``, ``dvy`` are typically guessed parameters.
*   ``d``: Initial separation at t=0 (for unbound orbits).

Type 3: Relative Position and Velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Directly specify the relative Cartesian vectors.

*   ``dx``, ``dy``, ``dz``: Relative position vector components (r2 - r1).
*   ``dvx``, ``dvy``, ``dvz``: Relative velocity vector components (v2 - v1).

Type 4: Absolute Positions and Velocities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify exact coordinates for both bodies. Note that this overrides any centre-of-mass settings typically handled by the setup routine.

*   ``x1``, ``y1``, ``z1``: Position components of body 1.
*   ``vx1``, ``vy1``, ``vz1``: Velocity components of body 1.
*   ``x2``, ``y2``, ``z2``: Position components of body 2.
*   ``vx2``, ``vy2``, ``vz2``: Velocity components of body 2.

Units
-----

Most distance and velocity parameters accept unit strings. For example:

*   ``100 au``
*   ``1 pc``
*   ``5.2 km/s``
*   ``1000 yrs`` (for period input in Type 0)

Prefixes
--------

When used in setup files (like ``binary.setup``), these options will often have a prefix or label to distinguish between multiple orbits (e.g., in a hierarchical triple). Common prefixes include:

*   ``binary_a``, ``binary_e``...
*   ``orbit1_a``, ``orbit2_a``...

Check your specific ``.setup`` file for the exact variable names.

