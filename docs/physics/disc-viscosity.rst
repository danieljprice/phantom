Disc viscosity in Phantom
==========================

What is the Shakura-Sunyaev viscosity?
--------------------------------------
In viscous accretion disc theory (Pringle 1981), axisymmetric 
accretion discs evolve according to the following equation for the surface density:

.. math::
    
    \frac{\partial \Sigma}{\partial t} = \frac{3}{R} \frac{\partial}{\partial R} \left[ R^{1/2} \frac{\partial}{\partial R} \left( \nu \Sigma R^{1/2} \right) \right]

where :math:`\Sigma` is the surface density (mass per unit area),
:math:`R` is the radial coordinate (distance from the central object),
and :math:`\nu` is an effective viscosity coefficient (area per unit time).

Shakura and Sunyaev (1973) proposed the prescription

.. math::

    \nu = \alpha c_{\rm s} H

where :math:`\alpha` is a dimensionless parameter, :math:`c_s` is the sound speed,
and :math:`H` is the disc scale height. The idea is that the "largest eddy" would
be of order the scale height and that the turbulence motions would be of order
the sound speed.

In this model the mass flow onto the central object is entirely driven by
the effective viscosity, which is assumed to represent some kind of 
underlying microscopic turbulence.

Do I need to add a disc viscosity to my disc?
---------------------------------------------
Phantom is a 3D code so depending on the problem there may already be
non-Keplerian motions present in your disc that will drive the accretion flow
without requiring a viscosity "prescription". Common examples are:

- Discs can be gravitational unstable, driving large-scale gravitational turbulence

- Accretion can be driven by anything that makes large-scale spiral arms,
  driven by planets, binaries, external companions, or infall from the environment  

- The disc may be unstable to the magnetorotational instability (MRI) 
  if you include magnetic fields (MHD=yes)

In these cases you don't need to add an effective viscosity, and should
simply attempt to minimise the numerical viscosity in the code, by using 
the minimum dissipation required for shock capturing, namely::

   # shock capturing
               alpha =       0.000    ! MINIMUM shock viscosity parameter
            alphamax =       1.000    ! MAXIMUM shock viscosity parameter
                beta =       2.000    ! non-linear shock viscosity parameter

Adding Shakura-Sunyaev viscosity by modifying the shock-capturing viscosity
---------------------------------------------------------------------------
If you are simulating an axisymmetric disc and want the results to be
independent of the resolution, you can explicitly add a disc viscosity.

The default way of including disc viscosity in SPH is to modify the
shock-capturing viscosity to mimic a physical viscosity. Details are given in
`Lodato & Price
(2010) <http://ui.adsabs.harvard.edu/abs/2010MNRAS.405.1212L>`__ or in the
Phantom paper. The Shakura-Sunyaev alpha is related to alpha_AV
according to

.. math::
    
    \alpha_{\rm SS} = 1/10 \alpha_{\rm AV} \frac{\langle h \rangle}{H}

where :math:`\langle h \rangle` is the mean smoothing length in an annulus and :math:`H` is the disc scale
height at that annulus.

In the .setup file for disc simulations this translation will be computed
automatically::
    
    # options for gas accretion disc
              alphaSS =       0.005    ! desired alphaSS (0 for minimal needed for shock capturing)

so you enter the "alpha_SS" value here and the code will compute
the minimum "alpha_AV" value required to achieve it, and will also
set the ``disc_viscosity`` flag to true in the .in file::

    # shock capturing
               alpha =  0.15836804    ! MINIMUM shock viscosity parameter
            alphamax =       1.000    ! MAXIMUM shock viscosity parameter
                beta =       2.000    ! non-linear shock viscosity parameter
      disc_viscosity =           T    ! use cs, multiply by h/|rij| and apply to approaching/receding

You are only guaranteed to achieve a constant alpha_SS if the ratio of h/<H> is constant,
which can be achieved by carefully choosing the temperature and density 
profiles of the disc (see Lodato & Price 2010).

Adding Shakura-Sunyaev viscosity as a Navier-Stokes viscosity
--------------------------------------------------------------
The other way to explicitly add viscosity is via the Navier-Stokes equation.
Here you can explicitly set a value for the kinematic viscosity, :math:`\nu`::

   # options controlling physical viscosity
           irealvisc =           0    ! physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)
          shearparam =       0.100    ! magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)
            bulkvisc =       0.000    ! magnitude of bulk viscosity

Option number 2 attempts to prescribe a constant alpha_SS. The problem is that
we do not in general know the local scale height, defined as

.. math::

    H = \frac{c_{\rm s}}{\Omega(R)}

where :math:`\Omega` is the orbital frequency. With ``irealvisc=2``, the code sets
the kinematic viscosity to

.. math::

    \nu = \alpha_{\rm SS} \frac{c_{\rm s}^2}{\Omega(R)}

where :math:`\alpha_{\rm SS}` is the desired Shakura-Sunyaev alpha,
and :math:`\Omega(R) = \sqrt{GM_*/R^3}` is the orbital frequency, a *prescribed*
function of the radial coordinate.

Impoartantly, this implies ONE star at the coordinate origin (the mass is read
from the mass1 parameter from the central potential, or from the first sink particle).

Currently this option does NOT give a sensible answer if you discs around multiple 
stars or sink particles that are not centred on the origin.

This is why in planet-disc interaction comparisons, it is common
to simply prescribe a constant value for :math:`\nu` itself (``irealvisc=1``).

What is the lowest alpha viscosity for a disc that Phantom can simulate?
------------------------------------------------------------------------

1. The lowest meaningful value of alpha_AV to "feel" an imposed viscosity
   is 0.1 (see e.g. Meru & Bate 2010; lower than this and you will get 
   spreading independent of the value of alpha_AV due to particle jitter).

2. For a disc simulated with 10^6 SPH particles we typically resolve the
   scale height by 2-5 resolution lengths, giving h/H ~ 0.2 – 0.5.

3. So therefore the lowest alpha_SS you can simulate at this resolution
   is 2 - 5 x 10^-3. To reach a lower alpha_SS, you have to lower /H,
   which means using more particles.

Further reading
---------------

- :doc:`/examples/disc` — interactive setup walkthrough

