Tidal disruption events
=======================

Setting up and relaxing a star
------------------------------

First, follow the :doc:`usual procedure for initiating a new simulation with
phantom </getting-started/running-first-calculation>`. We’ll use the “:doc:`grtde </user-guide/setups>`” setup.

make a new directory and write a local Makefile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    mkdir tde
    cd tde
    ~/phantom/scripts/writemake.sh grtde > Makefile

If you prefer Newtonian gravity use :doc:`tde </user-guide/setups>` instead of :doc:`grtde </user-guide/setups>`

compile phantom and phantomsetup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ make
   $ make setup
   $ ls
   Makefile    phantom*    phantomsetup*

run phantomsetup
~~~~~~~~~~~~~~~~

::

   ./phantomsetup tde

 writing setup options file tde.setup
  Edit tde.setup and rerun phantomsetup


relax the star
~~~~~~~~~~~~~~

Open the tde.setup file and make sure the parameter “relax_star” to True::

                   relax_star =       T    ! artificial damping of velocities (if on, v=0 initially)

Then run phantomsetup::

   ./phantomsetup tde.setup

Which will generate a sequence of snapshots of the relaxation process::

   RELAX-A-STAR-O-MATIC: Etherm:  0.499     Epot: -0.854     R*:   1.00
       WILL stop WHEN: dens error <   1.00% AND Ekin/Epot <   1.000E-07 OR Iter=0

   -------->   TIME =    0.000    : full dump written to file relax_00000   <--------

    Relaxing star: Iter   1/1000, dens error: 20.39%, R*:  0.978     Ekin/Epot:  3.157E-03
    Relaxing star: Iter   2/1000, dens error: 15.91%, R*:  0.978     Ekin/Epot:  4.858E-03
    Relaxing star: Iter   3/1000, dens error: 13.05%, R*:  0.978     Ekin/Epot:  3.303E-03
    Relaxing star: Iter   4/1000, dens error: 11.11%, R*:  0.977     Ekin/Epot:  2.276E-03
    Relaxing star: Iter   5/1000, dens error:  9.69%, R*:  0.975     Ekin/Epot:  1.655E-03

   -------->   TIME =   0.3756E-01: full dump written to file relax_00001   <--------
    Relaxing star: Iter   6/ 500, dens error: 22.17%, R*:   1.03     Ekin/Epot:  1.280E-03

This process will stop when either the density error is less than 1\% and the ratio of kinetic to potential energy
is less than the specified tolerance, OR the number of iterations reaches the specified maximum.

Once complete, you should obtain a relaxed initial conditions snapshot with the correct density and pressure profiles::

   -------->   TIME =    0.000    : full dump written to file tde_00000.tmp   <--------


    input file tde.in written successfully.
    To start the calculation, use:

    ./phantom tde.in


check the relaxed stellar profile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   splash relax1_* -y density -x r

If you have the latest version of splash, this should automatically plot the analytic profile for comparison based on the contents of the relax1.profile file.

Putting the star on an orbit for a tidal disruption event
---------------------------------------------------------
This should be done automatically using the “tde” or "grtde" setup and it controlled by the following set of options in the .setup file ::

   # options for orbit around black hole
      provide_params =           F    ! manually specify the position and velocity of the star(s)
                beta =       5.000    ! penetration factor
              ecc_bh =       0.800    ! eccentricity (1 for parabolic)
            theta_bh =       0.000    ! inclination of orbit (degrees)
             norbits =       5.000    ! number of orbits
       dumpsperorbit =         100    ! number of dumps per orbit

Running the simulation
----------------------

After this you can simply run phantom::

   ./phantom tde.in
