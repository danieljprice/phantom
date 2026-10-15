Ideal MHD
=========

Phantom has mature implementations of both ideal and non-ideal magnetohydrodynamics (MHD)
following the SPMHD formalism of Børve et al. (2001), Price & Monaghan (2004a,b, 2005),
Price (2012), Tricco & Price (2012), Tricco, Price & Bate (2016) and Tsukamoto (2026).
The non-ideal MHD implementation is described in Wurster et al. (2014, 2016).
Full details are given in the phantom paper (Price et al. 2018).

Magnetic field variables
-------------------------

Phantom evolves the magnetic field as :math:`\mathbf{B}/\rho`
according to the Lagrangian form of the induction equation, coupled to a scalar field
:math:`\psi` that is used to control the divergence of the magnetic field

.. math::

   \frac{{\rm d}}{{\rm d}t}\left(\frac{\mathbf{B}}{\rho}\right) = \left(\frac{\mathbf{B}}{\rho}\cdot\nabla\right)\mathbf{v} - \frac{\nabla \psi}{\rho},

where :math:`\psi` evolves according to

.. math::

   \frac{{\rm d}}{{\rm d}t}\left(\frac{\psi}{c_{\rm h}}\right) = -c_{\rm h} \nabla\cdot\mathbf{B} - \frac{\psi}{c_{\rm h}\tau} - \frac12 \frac{\psi}{c_{\rm h}} \nabla\cdot\mathbf{v}.

The code evolves :math:`\psi/c_{\rm h}` rather than :math:`\psi`, so that the
cleaning speed can vary in space and time without breaking energy conservation
(Tricco, Price & Bate 2016). Together these are stored in the ``Bevol`` array, as
:math:`\mathbf{B}/\rho` in ``Bevol(1:3)`` and :math:`\psi/c_{\rm h}` in ``Bevol(4)``.

The cleaning speed is the fast MHD wave speed, multiplied by ``overcleanfac``,

.. math::

   c_{\rm h} = f_{\rm clean}\sqrt{v_{\rm A}^2 + c_{\rm s}^2},

where :math:`v_{\rm A}` is the Alfvén speed, :math:`c_{\rm s}` is the sound speed
and :math:`f_{\rm clean}` is ``overcleanfac`` (default 1). Setting
``overcleanfac > 1`` cleans faster at the cost of shorter timesteps.
The damping timescale is

.. math::

   \tau = \frac{h}{\sigma c_{\rm h}},

where :math:`h` is the smoothing length and :math:`\sigma` is ``psidecayfac``
(default 1), which sets the ratio of parabolic to hyperbolic cleaning.

This is the constrained hyperbolic/parabolic cleaning of Tricco & Price (2012) and
Tricco, Price & Bate (2016). The constraint is that the total energy of the
magnetic and :math:`\psi` fields,

.. math::

   E = \sum_a m_a \left(\frac{B_a^2}{2\mu_0\rho_a} + \frac{\psi_a^2}{2\mu_0\rho_a c_{{\rm h},a}^2}\right),

is exactly conserved by the cleaning in the absence of the parabolic term, and can
only decrease when it is present, since then

.. math::

   \frac{{\rm d}E}{{\rm d}t} = -\sum_a \frac{m_a \psi_a^2}{\mu_0 \rho_a c_{{\rm h},a}^2 \tau_a} \leq 0.

Enforcing the divergence-free condition
----------------------------------------

Since Oct 2026 the code defaults are to use the above hyperbolic/parabolic cleaning but to also apply
a projection-to-zero at each synchronisation (``use_divb_proj = T`` in the runtime ``.in`` file).
This is performed at intervals of ``dtmax`` to ensure that divergence errors cannot build up over time.

Since the projection is tied to ``dtmax``, use ``nout`` to project more often
than dumps are written: for example ``dtmax`` ten times smaller with
``nout = 10`` projects ten times per dump.

The projection iterates until the root mean square of :math:`h|\nabla\cdot\mathbf{B}|/|\mathbf{B}|`
is below :math:`10^{-5}`, and the worst
1% of the divergence accumulated since the previous projection is reduced
by :math:`0.1/N`, where :math:`N` is the number of steps in between.
``divb_tol`` (default ``1e-8``) is a relative residual limit, and
the solve stops if that limit is reached first.

``divb_maxit`` (default 1000) is a safety cap on the
conjugate-gradient iteration. If it is reached first, the correction obtained
so far is still applied and :math:`\psi` is still set to zero.

The projection requires extra memory for the cached neighbour list.
If you don't have enough memory ("out of memory" from ``divb_proj``), set
``use_divb_proj = F`` to rely on cleaning alone.

Currently the div B projection does not work with MPI.

Non-ideal MHD
=============

Phantom implements all three non-ideal MHD terms (Ohmic resistivity, ambipolar
diffusion and the Hall effect) in the strong coupling approximation (Wurster et al. 2014, 2016).
These are controlled by the runtime parameters ``use_ohmic``, ``use_ambipolar`` and ``use_hall``.

Currently all of these are handled with explicit timestepping which means you can get
short timesteps and hence slow simulations if the non-ideal MHD coefficients become large.

Non-ideal MHD coefficients can either be specified as constants or calculated
self-consistently by the in-built NICIL library (Wurster 2016).
NICIL computes the ionisation balance including cosmic ray, thermal ionisation, and charged
dust grain species at a given density and temperature and returns the corresponding values of the
ambipolar, Hall and Ohmic diffusion parameters. See Wurster (2016) for details.