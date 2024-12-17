Radiation hydrodynamics in phantom
=====================================

There are many ways to skin a cat. And so, there are many ways to include the effects of radiation in phantom. 
These are outlined below, in order of increasing complexity.

Adiabatic vs isothermal equations of state
-------------------------------------------
The most basic way to bound the effects of heating and cooling is to compare
simulations with an adiabatic equation of state (ieos=2), where all energy is assumed to be trapped:

    :math:`P = (\gamma-1)\rho u`

and an isothermal equation of state (ieos=1, with code compiled using ISOTHERMAL=yes):

    :math:`P = c_{\rm s}^2 \rho`

where all energy is radiated away.  See :doc:`Equations of state in phantom <eos>`.

The difference between two simulations using the above methods is a measure of the importance of radiative cooling.
Importantly in neither of the above cases is the temperature of the gas actually computed, though one is free to do so.
One can simply interpret the sound speed or internal energy in terms of temperature after the fact
by choosing code units, but the simulations themselves do not care about the units and are scale-free.

For accretion discs using a strictly isothermal equation of state is a bad idea
so we typically prescribe a fixed temperature profile (e.g. T is a prescribed function of radius)
rather than T=const. There are various options for this in phantom, see :doc:`Equations of state in phantom <eos>` and :doc:`Accretion discs </examples/disc>`.
As above, the temperature is *not actually computed*, but rather the *sound speed* is prescribed.

Radiation pressure
-------------------
The next level of sophistication, assuming radiation is perfectly trapped (optically thick), is to include the effect
of radiation pressure using :doc:`the gas+radiation equation of state (ieos=12) <eos>`.
Here, we assume that the internal energy represents the total specific energy, ie. 
the sum of the gas thermal energy and the radiation energy.
More specifically:

   :math:`u = \frac{3}{2} \frac{k_B T}{\mu m_H} + a T^4`

This equation can be solved backwards for temperature, and the resulting temperature used to
compute the pressure:

   :math:`P = \frac{k_b T}{\mu m_H} + \frac13 a T^4` 

See :doc:`documentation for ieos=12 <eos>`. One important difference
compared to the adiabatic vs. isothermal equation of states is that real physical constants
are needed to compute P from u, and so the units of the simulation are no longer arbitrary.

The main effect of radiation pressure is to reduce the effective adiabatic index of the gas.
One can see this because when the gas pressure term dominates, we have

   :math:`P = \frac23 \rho u`

implying :math:`\gamma=5/3`, whereas when the radiation pressure term dominates, we have

   :math:`P = \frac13 \rho u`

implying :math:`\gamma=4/3`.

Ionization and recombination
------------------------------
One level up from this is to also include ionization and recombination using the
:doc:`the gas+rad+rec equation of state (ieos=15) <eos>`, which utilises a similar
approach to the gas+radiation equation of state above but also solves a Saha equation
for Hydrogen and Helium ionization/recombination. One can also achieve similar physics
using the :doc:`tabulated MESA equation of state (ieos=10) <eos>` but with less control over
which ionizations to switch on/off. One also should be careful that the simulation remains
within the density and temperature ranges of the table.

Radiation transport with flux-limited diffusion
-------------------------------------------------
Next level, we can allow the radiation to diffuse. The implicit flux-limited diffusion
scheme implemented in phantom is from Whitehouse & Bate (2004) and Whitehouse, 
Bate and Monaghan (2005). Flux limited diffusion is correct when the gas is optically
thick but makes an approximation to limit the radiation propagation to the speed of
light in optically thin regions, but still treats propagation of radiation as diffusion
in these regions, which does not capture things like shadowing.

In the code this formulation differs by splitting the gas and radiation pressure
into separate contributions. The gas specific internal energy is stored in the vxyzu array,
while the radiation specific internal energy is stored in the rad array.

Hence the total internal energy per unit mass in this case is::

    ui = vxyzu(4,i) + rad(iradxi,i)

Hence in this case **you need to set both the gas and radiation internal energies during
the setup of the simulation**, or when injecting particles.

As a first step, you can reproduce the gas + radiation equation of state above
by assuming an infinite opacity, by setting the appropriate flag in the input file::

    iopacity_type = 0       ! opacity method (0=inf,1=mesa,2=constant,-1=preserve)

More generally, one can set the opacity as a function of density and temperature from 
the MESA tables (iopacity_type=1) which then allows radiation to diffuse. 
Ionization/recombination is included in these tables, but at the moment not 
molecular hydrogen formation.

For modelling stars other than those powered by contraction, including the 
leakage of radiation in this way requires you to supply a heating source, otherwise 
the star will just steadily cool. For red giants with a sink particle core a simple procedure
is to supply a constant luminosity input from the core (:doc:`sink heating <sinks>`). This is experimental.
Another option would be to include a nuclear burning network (please somebody contribute this).


Polytropic radiative cooling approximation
-------------------------------------------
This is an alternative to computationally expensive radiative transfer in regimes where radiative cooling is 
important. This method estimates the optical depth for each particle and its equilibrium temperature. From 
these the new temperature and internal energy is updated at each timestep. The method implemented here is 
the "modified Lombardi" method of `Young et al. (2024) <https://ui.adsabs.harvard.edu/abs/2024MNRAS.531.1746Y>`__,
 which was based on Stamatellos et al. (2007) and Lombardi et al. (2015). This method is designed for 
 self-gravitating discs around a central star. Stellar heating is included from the most massive sink particle.
 
Use icooling = 9 and ieos = 23 to use the tabulated equation of state which has the opacity tables required
for the cooling calculation. The additional parameters are::

	EOS_file =  myeos.dat   ! File containing tabulated EOS values
	Lstar =       0.440    ! Luminosity of host star for calculating Tmin (Lsun)
	Tfloor =       5.000    ! temperature floor (K); on if > 0
	
	
N.B. This version does not currently include radiation from more than one sink particle. This version does not 
couple with flux-limited diffusion at the moment.


Irradiation from stars with phantom + MCFOST
---------------------------------------------
In regimes where the radiation diffusion time is relatively short, anything not inside stars
or where temperatures are set by external irradiation, a better approach is to use
the :doc:`coupled version of phantom and MCFOST </external-utilities/mcfost>`. 

In this procedure we call MCFOST at discrete intervals (set by the dtmax parameter in the .in file)
which emits and propagates photons until radiative equilibrium is reached. This is a 
good approximation if the time to reach radiative equilibrium is shorter than the time interval
between calls to MCFOST, which is true for example in most protoplanetary discs.

One can also include PdV work and shock heating contributions in the calculation of
radiative equilibrium, so this allows for shock heating from the gas as well as heating
from central stars. The (dust) temperature we receive back from MCFOST is then simply used
as the gas temperature, and is kept fixed on particles between calls to MCFOST.

MCFOST by default assumes the only source of opacity is dust, and if dust is not used in
the simulation will assume that dust is 1% of the gas in order to compute opacities. If the
simulation includes dust species (with :doc:`DUST=yes </user-guide/config>`) or dust formation (with :doc:`NUCLEATION=yes </user-guide/config>`) then
opacities will be computed from the dust information in the simulation. See :doc:`Using Phantom with MCFOST </external-utilities/mcfost>`

Using MCFOST for gas radiative transfer, using the atomic line transfer capabilities outlined
by `Tessore et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021A&A...647A..27T>`__ is also possible but more experimental.

See also
--------

- :doc:`Equations of state available in Phantom <eos>`
- :doc:`MCFOST </external-utilities/mcfost>`
