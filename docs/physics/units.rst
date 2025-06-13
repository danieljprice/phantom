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

    P' = P \frac{M}{L^3} \frac{L^2}{T^2} = \frac{P M}{L T^2}

    \rho' = \rho \frac{M}{L^3}

    v' = v \frac{L}{T}

    u' = u \frac{L^2}{T^2}

If we then rewrite the equations in terms of these dimensionless variables, we find

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

    T = \sqrt{G M/L^3}

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

The unit scalings (umass, udist, utime) are written in the header of the dump files for use by SPLASH and other
post-processing tools.

Computing temperature, or using equations of state that depend on temperature
------------------------------------------------------------------------------

The main need for physical units is to compute the temperature of the gas. Again, 
the equations of hydrodynamics only care about the internal energy. For an ideal gas
one can relate temperature to internal energy using:

.. math::

    u = \frac{3}{2} \frac{k_B T}{\mu m_H}

where :math:`k_B` is the Boltzmann constant, :math:`\mu` is the mean molecular weight,
and :math:`m_H` is the mass of the hydrogen atom. Notice that this involves two more
physical constants: the Boltzmann constant :math:`k_B` and the mass of the hydrogen atom
:math:`m_H`. So computing temperature from internal energy breaks the scale-free nature
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
and :math:`G=1` which leaves only one unit freedom, typically the mass unit

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

For example, the Schwarzschild radius of a black hole
is :math:`R_S = 2 G M / c^2` and the last stable orbit for a test particle is :math:`R_H = 3 R_S`.
In geometric units, these are :math:`R_S = 2 M` and :math:`R_H = 6 M`.






