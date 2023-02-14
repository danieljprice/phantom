Sink particles in phantom
=========================================================

Sink particles (Bate et al. 1995) are used to represent point masses (usually stars or planets) in phantom calculations.
These particles are treated separately to the SPH particles and have the following
properties:

- sink-sink gravity is always computed by direct sum
- sink-gas gravity is mutual and computed by direct sum
- sinks evolve on a separate timestep to the gas
- sinks can accrete gas in a momentum and angular-momentum conserving manner
- can be used to represent stars, or the cores of stars, or planets
- there is no density field computed on sink particles
- no gravitational softening is applied to sink-sink interaction by default
- can emit winds, radiation or other feedback onto the simulation

Sink particles compared to external forces
-------------------------------------------
The main advantage of a sink particle over a fixed external potential
is that they satisfy the conservation laws. For example, a disc of gas
particles orbiting in a 1/r potential will not conserve momentum, but
a disc of gas particles orbiting a sink particle will.

Sink particles compared to fixed binary potential
--------------------------------------------------
A circumbinary disc around a sink particle pair will similarly satisfy the conservation
of linear and angular momentum, unlike a prescribed binary potential. The
backreaction of the gas disc onto the sink means that the binary can naturally
gain or lose angular momentum to the disc, and will gain mass due to accretion.

Sink particle properties
-------------------------
As well as position, velocity, mass and acceleration, sink particle have the following properties:

- hacc: accretion radius, inside of which gas particles are tested for accretion
- hsoft: softening length for gas-sink interaction (zero by default)
- macc: total accreted mass (to avoid round-off error, as this is often a small fraction of total mass)

The full list of extended properties (extracted from `part.F90 <https://github.com/danieljprice/phantom/blob/master/src/main/part.F90>`__) is as follows:

.. include:: sink-properties.rst

These are stored in the following arrays:

- xyzmh_ptmass (positions, mass, accretion radius and extended properties)
- vxyz_ptmass (velocities)
- fxyz_ptmass (accelerations)
- dsdt_ptmass (time derivative of spin, where oblateness is used)

Where to look in the code
--------------------------
Sink particle functionality is implemented in the `ptmass module <https://github.com/danieljprice/phantom/blob/master/src/main/ptmass.F90>`__.
