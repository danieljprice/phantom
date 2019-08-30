Running a simulation with dust growth and fragmentation
=======================================================

The dust growth and fragmentation algorithm is described in `Vericel et
al. (in prep) <https://media.giphy.com/media/XIqCQx02E1U9W/giphy.>`__.

.. important::

 Dust growth is not restrained to single star simulations!
 You can use dustgrowth around binaries, you can put planets in your
 disc, you can tilt and warp everything if you want and many more. Be
 creative :).

If you find a bug, please report it either on the phantom slack channel
or by sending me an email at arnaud.vericel@univ-lyon1.fr

The algorithm is tested automatically by the nightly test suite, but you
can check the test manually in your phantom repository:

::

   cd ~/phantom; make testgrowth

Dust growth is a compile time option and you have two ways of setting it
up.

You can import the Makefile using the dustydisc setup and set
``DUSTGROWTH`` to yes before compiling:

::

   ~/phantom/scripts/writemake.sh dustydisc > Makefile
   export DUSTGROWTH=yes
   make; make setup

You can directly import the Makefile using the “growingdisc” setup:

::

   ~/phantom/scripts/writemake.sh growingdisc > Makefile
   make; make setup

When you make the setup file ``./phantomsetup ilovegrowth``, dust growth
adds an option block after the dust section:

::

   # options for growth and fragmentation of dust
                  ifrag =           1    ! fragmentation of dust (0=off,1=on,2=Kobayashi)
                  isnow =           0    ! snow line (0=off,1=position based,2=temperature based)
                  rsnow =        100.    ! snow line position in AU
                  Tsnow =        150.    ! snow line condensation temperature in K
                  vfrag =         15.    ! uniform fragmentation threshold in m/s
                vfragin =       5.000    ! inward fragmentation threshold in m/s
               vfragout =         15.    ! inward fragmentation threshold in m/s
           grainsizemin =       0.005    ! minimum allowed grain size in cm

Here’s a brief description of each of them (remember that they are
described in more details in `Vericel et al. (in
prep) <https://imgflip.com/i/389twd>`__)

::

                  ifrag =           1    ! fragmentation of dust (0=off,1=on,2=Kobayashi)

switch choosing between pure growth (= 0), growth and harsh
fragmentation (= 1) or growth and smoother fragmentation (= 2)

::

                  isnow =           0    ! snow line (0=off,1=position based,2=temperature based)

switch setting up a snow line by distance to the star (= 1) or by its
corresponding temperature (= 2)

::

                  rsnow =        100.    ! snow line position in AU

snow line distance to the star if isnow = 1

::

                  Tsnow =        150.    ! snow line condensation temperature in K

snow line corresponding temperature to the star if isnow = 2

::

                  vfrag =         15.    ! uniform fragmentation threshold in m/s

uniform fragmentation threshold velocity above which fragmentation
occurs (if ifrag > 0 and isnow = 0)

::

                vfragin =       5.000    ! inward fragmentation threshold in m/s
               vfragout =         15.    ! inward fragmentation threshold in m/s

inner and outer fragmentation threshold velocities above which
fragmentation occurs (if ifrag > 0 and isnow > 0)

::

           grainsizemin =       0.005    ! minimum allowed grain size in cm

minimum size you allow dust grains to reach if ifrag > 0

After re-executing the setup file ``./phantomsetup ilovegrowth.setup``,
``ilovegrowth.in`` contains an extra dust growth block only with the
relevant informations, e.g here:

::

   # options controlling growth
                  ifrag =           1    ! dust fragmentation (0=off,1=on,2=Kobayashi)
           grainsizemin =       0.005    ! minimum grain size in cm
                  isnow =           0    ! snow line (0=off,1=position based,2=temperature based)
                  vfrag =         15.    ! uniform fragmentation threshold in m/s

The relative velocity between dust particles in computed using the
viscosity parameter α which is stored in the ``shearparam`` parameter in
the physical viscosity block of the infile.

::

   # options controlling physical viscosity
              irealvisc =           0    ! physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)
             shearparam =       0.010    ! magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)

Independently of the use of physical viscosity or not
(``irealvisc = 0 or 2``), the algorithm will use ``shearparam`` as the
value for α. However, if you set ``irealvisc`` to 1 and set a constant
viscosity, dustgrowth will be unavailable. In the case where you have no
physical viscosity, make sure ``shearparam`` is equal to the ``alphaSS``
that you specified during the setup. This is important for consistency.

After launching the simulation ``./phantom ilovegrowth.in``, the
dumpfiles will contain the grain sizes, grain densities, the ratio
Vrel/Vfrag (if ifrag > 0, else gives Vrel/(1 m/s)) and the Stokes number
St.

Tips to perform great dust growth simulations
---------------------------------------------

Set more gas particles than dust particles. I typically recommend a
ratio of 10 gas particles per dust particle (e.g. 10^6 gas and 10^5
dust).

Let the gas relax before injecting dust. This is pretty important,
especially for small grains. To do that, the Makefile that you imported
from growingdisc
``~/phantom/scripts/writemake.sh growingdisc > Makefile`` has a moddump
utility that you can compile using ``make moddump`` (:doc:`see
here <moddump>`).

I also suggest you use another moddump to remove the few particles that
eventually got ejected at large distances from the star. This will also
accelerate your simulation. You can make that utility with
``make moddump MODFILE=moddump_removeparticles_radius.f90``.

If you consider fragmentation, ``grainsizemin`` should not be too
small or else tiny grains will dictate the timestepping and make your
simulation ridiculously slow. I typically recommend to set that minimum
to 10 to 50 μm.

**Have fun :) and make sure to cite the paper** `Vericel et al. (in
prep) <https://imgflip.com/i/38bw62>`__
