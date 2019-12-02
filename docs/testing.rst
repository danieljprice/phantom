Making sure your changes are safe to commit
===========================================

There are two automated checking procedures that you can use to make
sure that changes you are about to commit do not break the code. These
are both run nightly anyway (and will email you the results), but here
is how to make sure that your changes do not break the code for others.

Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

You can run the whole test suite using:

::

   make test

This is just a shortcut for the following sequence of commands:

::

   make SETUP=test
   ./bin/phantom test

The test suite is also run nightly to check for regressions. You can run
the complete nightly test yourself using the testbot wrapper script:

::

   cd phantom/scripts
   ./testbot.sh

Running selected parts of the test suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run just part of the test suite by giving an additional argument
as follows:

::

   make SETUP=test && ./bin/phantom test derivs

A non-exhaustive list of possible arguments are as follows:

+-----------------------------------+-----------------------------------+
| argument                          | description                       |
+===================================+===================================+
| derivs                            | unit tests of all derivative      |
|                                   | terms                             |
+-----------------------------------+-----------------------------------+
| derivshydro                       | unit tests of hydrodynamic        |
|                                   | derivatives                       |
+-----------------------------------+-----------------------------------+
| derivsmhd                         | unit tests of the                 |
|                                   | magnetohydrodynamics derivatives  |
+-----------------------------------+-----------------------------------+
| derivsav                          | unit tests of the artificial      |
|                                   | viscosity terms                   |
+-----------------------------------+-----------------------------------+
| derivsvisc                        | unit tests of the physical        |
|                                   | viscosity terms                   |
+-----------------------------------+-----------------------------------+
| derivsind                         | unit tests related to individual  |
|                                   | particle timesteps                |
+-----------------------------------+-----------------------------------+
| derivscontrast                    | unit tests of the hydro           |
|                                   | derivatives in a setup with       |
|                                   | density contrast                  |
+-----------------------------------+-----------------------------------+
| link                              | tests the linklist/neighbour      |
|                                   | finding modules                   |
+-----------------------------------+-----------------------------------+
| step                              | performs checks on the            |
|                                   | timestepping routine              |
+-----------------------------------+-----------------------------------+
| ptmass                            | performs checks of the sink       |
|                                   | particle/point mass module        |
+-----------------------------------+-----------------------------------+
| externf                           | performs tests of the external    |
|                                   | force module                      |
+-----------------------------------+-----------------------------------+
| sedov                             | performs a sedov blast wave test  |
+-----------------------------------+-----------------------------------+
| indtstep                          | performs unit tests of various    |
|                                   | utilities to do with individual   |
|                                   | timestepping                      |
+-----------------------------------+-----------------------------------+
| gravity                           | perform tests of self-gravity     |
|                                   | terms                             |
+-----------------------------------+-----------------------------------+
| dump                              | perform unit tests of write       |
|                                   | to/read from dump files           |
+-----------------------------------+-----------------------------------+
| kernel                            | performs unit tests of the kernel |
|                                   | module                            |
+-----------------------------------+-----------------------------------+

The buildbot
~~~~~~~~~~~~

The buildbot runs nightly and checks that the code compiles in all of
the possible SETUP configurations in the Makefile. You can run this
yourself as follows:

::

   cd phantom/scripts
   ./buildbot.sh

Because some of the setups are compiled with large static arrays which
may exceed the memory limits of your local machine, you can optionally
check only those setups that have idim less than some value, e.g.:

::

   cd phantom/scripts
   ./buildbot.sh 1000000

which checks only SETUPs with idim set to 1 million particles or less in
the dim file.
