Getting your code to pass the github actions
============================================

On every pull request a sequence of continuous integration tests 
are performed to check that code is safe to merge into master.
The scripts in the .github/workflows directory are as follows:

- build: checks that phantom, phantomsetup, phantomanalysis and phantommoddump compile with every possible SETUP= flag
- test: runs the test suite [see below]
- mpi: runs the test suite with MPI [see below]
- mcfost: compiles and links phantom+mcfost

Running the test suite on your own machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run the test suite using:

::

   make test

This is just a shortcut for the following sequence of commands:

::

   make SETUP=test phantomtest
   ./bin/phantomtest

You can run the complete testsuite yourself using the testbot wrapper script:

::

   cd phantom/scripts
   ./testbot.sh

Running selected parts of the test suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run just part of the test suite by giving an additional argument
as follows:

::

   make SETUP=test phantomtest && ./bin/phantomtest derivs

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

which checks only SETUPs with maxp set to 1 million particles or fewer.

FAQ on common reasons for failure
=================================
We enforce the following policies in merging to the master branch:

1. Code must compile and run with and without DEBUG=yes
2. Code must compile with ability to change the precision of reals to real*4 (this is enforced in SETUP=blob)
3. Code must compile with no warnings when compiled with gfortran (enforced with NOWARN=yes which adds the -Werror flag)
4. Testsuite must work with and without MPI, i.e. compile with MPI=yes


