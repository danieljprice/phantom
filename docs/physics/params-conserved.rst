Conservation checks
==============================

Phantom verifies the main conservation laws at every time step. While these quantities
cannot be exactly conserved, Phantom requires that the conservation error remains below a certain threshold.
The limits on conservation errors enforced by the code are listed in the table below, as well as the references 
to justify them and where to find them in code. Since Phantom will stop if any of these thresholds are reached, 
you may wish to turn these checks off (e.g. when testing something). 
To do so, set the I_WILL_NOT_PUBLISH_CRAP environment variable to “yes”:

.. code-block:: bash

    export I_WILL_NOT_PUBLISH_CRAP=yes

You can also change the value of the thresholds where they are defined in the code, although this is not recommended 
unless you know what you are doing. These checks are turned off automatically in some setups that use more advanced 
features in Phantom, more details on this below.

.. include:: params-conserved-list.rst

Some of the conservation checks listed above are automatically turned off when using some specific features of Phantom:

- Total energy conservation is checked only when using an adiabatic EoS with NO cooling processes and ALL the heating terms (PdV work, shock heating, resistive heating when using MHD) 
- Linear momentum conservation is not checked as long as boundary particles are used.
- Angular momentum is checked only if no boundaries are used, and no non-radial forces (iexternalforce>1) are used.
- Energy, momentum and angular momentum conservation are not checked if particle injection, dynamic boundaries, or APR (adaptive particle refinement) are used.