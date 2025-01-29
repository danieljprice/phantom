Simulation of self-gravitating and gravitationally unstable accretion discs
============================================================================

In PHANTOM, it is possible to perform simulations of accretion discs taking into account the role of the disc self-gravity and, possibly, to trigger gravitational instability. 

Self-gravity in (vertically isothermal) accretion discs
--------------------------------------------------------
Normally in the *disc* and *dustydisc* environments the disc self-gravity is not taken into account. Indeed, usually the disc to star mass ratio is so low that the disc contribution to the gravitational potential is negligible. 
The setups *isosgdisc* and *dustyisosgdisc* allow simulating vertically isothermal discs with self-gravity and disc viscosity ($α_{SS}$), and work as the *disc* and *dustydisc* ones. 

Gravitationally unstable accretion discs
------------------------------------------
The environments *sgdisc* and *dustysgdisc* allow the user to simulate a gravitationally unstable accretion disc. Gravitational instability is triggered by cooling, using an adiabatic equation of state and without disc viscosity. Usually, a cooling law it is prescribed, and the simplest one has been proposed by Gammie (2001): the cooling timescale $t_{cool}$ is assumed to be proportional to the dynamical time of the disc $t_{dyn}$ so that $t_{cool} = β t_{dyn}$. 

In the *setup* file there are the disc parameters, and in the *input* file it is possible to prescribe the cooling law. In particular, the variables to pay attention to are:

- ieos = 2 , to choose an adiabatic equation of state
- icooling = 3, to choose a β cooling prescription with constant β (icooling = 7 prescribes a varying β cooling with the radius)
- beta_cool = #, to choose the value of β cooling

To use the radiative cooling approximation of Young et al. (2024) use:

- ieos = 24
- icooling = 9

See :doc:`Radiation hydrodynamics in phantom </physics/radiation>`.
