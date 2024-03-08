Composition tracking in phantom
=====================================

This file documents how to track chemical abundances and mixing in phantom.

Composition tracking with fixed abundances
--------------------------------------------
Tracking the chemical composition of the gas in phantom with fixed
abundances is straightforward, since phantom is a Lagrangian code and
the particle identifiers are preserved throughout the simulation.

In the non-MPI code without particle injection, the particles are 
always written to the dump files in the same order, so the particle
id is simply the particle index in the dump file. In the MPI code,
the particle id is stored in the 'iorig' array in the dump file.

Hence composition tracking can done as a post-processing step.

Writing a .comp file
~~~~~~~~~~~~~~~~~~~~~
In practice, the composition of the gas can be tracked by writing a .comp file
either once for the entire simulation, or one per dump file. The .comp file
is a simple ascii file with one row per particle, where the columns are the abundances.
You should also give the first line a header with labels for each element::

    # h1, he3, he4, c12, n14, o16, ne20, mg24, si28, s32, ar36, ca40, ti44, cr48, fe52, ni56
    7.1119142075E-01   9.3180507858E-05   2.7341075111E-01  ...
    7.1119141814E-01   9.3180815248E-05   2.7341075311E-01  ...
    7.0179434509E-01   3.6580185766E-06   2.8264795926E-01  ...
    ...

If the phantom dump is called `foo_0000`, the .comp file should be called `foo_0000.comp`
or simply `foo.comp` if the composition is constant on each particle.

This file is automatically read by splash and used to create extra columns in the visualisation.

See also
--------

- :doc:`Setting up stars and tidal disruption events </examples/star>`
