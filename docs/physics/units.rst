Units in Phantom
================

Phantom solves the equations of hydrodynamics in code units. Depending on the physics, the results of
the simulation can be compared to observations merely by rescaling the units of length, time, and mass
appropriately in post-processing.

Hydrodynamics
-------------

Code units can be understood by writing out the equations of hydrodynamics:

.. math::

    \frac{d\rho}{dt} + \rho \nabla \cdot \mathbf{v} = 0

    \frac{d\mathbf{v}}{dt} = -\frac{\nabla P}{\rho}

    \frac{du}{dt} = -\frac{P}{\rho} \nabla \cdot \mathbf{v}

where :math:`P` is the pressure, :math:`\rho` is the density, 
:math:`\mathbf{v}` is the velocity and :math:`u` is the specific internal energy.

Now, consider rescaling the units of length, time, and mass by arbitrary 
constants :math:`T`, :math:`L`, and :math:`M`.
Using these, we define dimensionless variables for pressure, density, velocity, and 
internal energy according to

.. math::

    P' = P / \left(\frac{M}{L^3} \frac{L^2}{T^2}\right)

    \rho' = \rho / \frac{M}{L^3}

    v' = v / \frac{L}{T}

    u' = u / \frac{L^2}{T^2}

If we then rewrite the equations in terms of these dimensionless variables (e.g. using :math:`dt = T dt'`), we find

.. math::

    \frac{d\rho'}{dt'} + \rho' \nabla' \cdot \mathbf{v'} = 0

    \frac{d\mathbf{v'}}{dt'} = -\frac{\nabla' P'}{\rho'}

    \frac{du'}{dt'} = -\frac{P'}{\rho'} \nabla' \cdot \mathbf{v'}.

We see that the equations are the same, but with the new units. Hence
we find that the equations of hydrodynamics are scale-free.

Hydrodynamics with gravity
--------------------------

The situation changes slightly if we add more physics to the simulation. 
For example, if we add gravity from a central point mass, the equations of motion become

.. math::

    \frac{d\mathbf{v'}}{dt'} = -\frac{\nabla' P'}{\rho'} - \frac{G M}{r'^2} \mathbf{\hat{r}'}

which is no longer scale-free because of the gravitational constant. However, we can
rescale the constant :math:`G` to unity to make the equations scale-free again. 
Typically we rescale the time unit to 

.. math::

    T = \sqrt{\frac{L^3}{G M}}

which means that in code units, :math:`G=1`. This is the default in Phantom, but
the choice of mass and length units remain arbitrary and for many setups
can be set as a string in the .setup file:

.. code-block:: text

    # units
            dist_unit =          au    ! distance unit (e.g. au,pc,kpc,0.1pc)
            mass_unit =      solarm    ! mass unit (e.g. solarm,jupiterm,earthm)

Again, the only reason units are needed is to interpret the results. They are not
used in the calculation itself, except sometimes in the setup of the simulation if the
user wants to enter parameters in physical units.

The unit scalings (``umass``, ``udist``, ``utime``) are written in the header of the
dump files for use by SPLASH and other post-processing tools.

Computing temperature, or using equations of state that depend on temperature
------------------------------------------------------------------------------

The main need for physical units is to compute the temperature of the gas. Again, 
the equations of hydrodynamics only care about the internal energy. For an ideal gas
one can relate temperature to internal energy using:

.. math::

    u = \frac{3}{2} \frac{k_B T}{\mu m_{\mathrm{H}}}

where :math:`k_B` is the Boltzmann constant, :math:`\mu` is the mean molecular weight,
and :math:`m_{\mathrm{H}}` is the mass of the hydrogen atom. Notice that this involves two more
physical constants: the Boltzmann constant :math:`k_B` and the mass of the hydrogen atom
:math:`m_{\mathrm{H}}`. So computing temperature from internal energy breaks the scale-free nature
of the equations. However the temperature is not actually used in the calculation, just the
internal energy (or the sound speed) and is only needed for interpretation.

If we include additional physics, such as cooling or heating, or 
radiation pressure in the equation of state, then the temperature *is* needed in the
calculation. In this case, the equations are no longer scale-free because the unit conversions
are used in the actual simulation to compute the temperature, not just for post-processing.

Geometric units for special and general relativity
--------------------------------------------------

In relativity the equations involve the speed of light :math:`c` as well as
the gravitational constant :math:`G`. Here we choose units such that :math:`c=1`
and :math:`G=1` which leaves only one unit freedom, typically the mass unit.

Hence when using phantom with GR=yes, only the mass unit is changeable:

.. code-block:: text

    # units
            mass_unit =      solarm    ! mass unit (e.g. solarm,jupiterm,earthm)

Typically one sets the mass unit to the black hole mass when using the Schwarzschild
or Kerr metric such that :math:`M=1`, but this is not required.

Geometric units have a nice interpretation in terms of the Schwarzschild radius
and the event horizon radius. With G=c=1:

- velocity is in units of the speed of light,
- length is in units of :math:`GM/c^2` 
- internal energy is in units of :math:`c^2`
- the unit of time is :math:`\sqrt{G M / c^3}`

For example, the Schwarzschild radius of a black hole in geometric units
is :math:`R_S = 2 G M / c^2 = 2M` and the last stable orbit for a test particle is :math:`R_H = 3 R_S = 6M`,
and if :math:`M=1` these are just at :math:`r=2` and :math:`r=6` in code units.

Practical example of how to rescale your simulation to arbitrary units
----------------------------------------------------------------------
Unit rescaling can be better understood with a practical example. Let's imagine we performed
a simulation of a 1 Jupiter-mass circumbinary disc orbiting a pair of binary stars consisting
of a 0.5 solar mass primary star and a 0.5 solar mass secondary, on a circular orbit with a semi-major
axis of 10 au. How should we rescale our results to re-interpret the simulation as a circumbinary disc
around binary supermassive black holes in the centre of a galaxy?

This is where the dimensionless units come in. In the first case we can imagine choosing units of
:math:`M=1` solar mass, :math:`L=1` au and, since the simulation involves gravity, time units
such that :math:`G=1`. You may already notice the similarity to Kepler's third law 
(:math:`P = 2\pi \sqrt{a^3/(GM)}`, where :math:`P` is the period and :math:`a` is the semi-major axis) in the time unit
Indeed, we find that with our choice of code units the time unit is

.. math::

    T = \sqrt{\frac{L^3}{GM}} = \sqrt{\frac{(1.5 \times 10^{13} \mathrm{cm})^3}{6.67 \times 10^{-8} \times 2 \times 10^{33} \mathrm{g}}} = 5 \times 10^6 \mathrm{s} = 58 \mathrm{days} = \frac{1}{2\pi} \mathrm{yrs}

So when we look at the raw simulation outputs (use --code flag in splash to plot in code units) we see
that the mass of our two point masses are :math:`0.5` in code units, while
our binary is separated by 10 code units of length, and completes one orbit at :math:`t=2\pi`.

Reinterpreting this simulation for two SMBHs is thus trivial. We just have to choose different
values of :math:`M` and :math:`L` when interpreting the results. For example, if we choose
:math:`M=10^6` solar masses and :math:`L=1 kpc` then nothing about the simulation changes, we just
interpret :math:`M=0.5` as :math:`0.5 \times 10^6 M_\odot`, a separation of 10 in code units as 10 kpc,
and a time of :math:`t=2\pi` in code units as

.. math::

    t = T \times t_\mathrm{code} = \sqrt{\frac{L^3}{GM}} 2\pi = \sqrt{\frac{(10 \mathrm{kpc})^3}{G 10^6 M_\odot}} = 14,873 \mathrm{yr}

which, as previously, is the orbital period of our binary divided by :math:`2\pi` (this occurs because our mass
unit is equal to the total mass of the binary).

Similarly, we can reinterpret the internal energy for the two different cases, using

.. math::

    u_\mathrm{physical} = L^2/T^2 u_\mathrm{code}

Where :math:`M` and :math:`T` and :math:`L` are written in the code as ``umass``, ``utime`` and ``udist``.
The units module also provides conversion factors for specific variables, like ``unit_ergg`` which is :math:`L^2/T^2`.

This would then give very different values for temperature in the black hole case compared to the stellar binary,
even though the underlying dimensionless internal energy in the two simulations is the same.

Unit conversion routines
------------------------

The units module (src/main/units.f90) provides several utility routines for converting between code units 
and physical units. These are particularly useful when setting up simulations or analyzing results.

The main routines are:

1. ``in_code_units(string,ierr,unit_type)``
   Converts a string like '10.*days' or '10*au' into code units. For example:
   
   .. code-block:: fortran
   
      real :: time
      integer :: ierr
      time = in_code_units('10.*days',ierr,'time')
   
   This will convert 10 days into code units. The routine handles:

   - Time units (s, min, hr, day, yr, Myr)
   - Length units (cm, m, km, au, pc, kpc, Mpc)
   - Mass units (g, solarm, earthm, jupiterm)
   - Density units (g/cm^3, kg/m^3)
   - Velocity units (cm/s, m/s, km/s, c)
   - Mass accretion rate units (g/s, Msun/yr)

2. ``in_units(val,unitstring)``
   Converts a value from code units to physical units. For example:
   
   .. code-block:: fortran
   
      real :: mass_in_msun
      mass_in_msun = in_units(mass,'solarm')
   
   This will convert a mass value from code units to solar masses. The routine supports the same unit types as ``in_code_units``.

3. ``select_unit(string,unit,ierr,unit_type)``
   A lower-level routine that recognizes unit strings and returns the appropriate conversion factor. This is used internally by the other routines.

   These routines are useful when:

   - Setting up initial conditions in physical units
   - Converting simulation results to physical units for comparison with observations
   - Reading parameters from input files that are specified in physical units
   - Writing output in physical units for analysis

The unit conversion routines handle common astronomical units automatically, and the conversion factors are updated whenever the code units are changed (e.g., when using geometric units for relativity). 


