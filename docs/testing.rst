Getting your code to pass the github actions
============================================

On every pull request a sequence of continuous integration tests
are performed to check that code is safe to merge into master.
The scripts in the `.github/workflows <https://github.com/danieljprice/phantom/tree/master/.github/workflows>`_ directory are as follows:

- |build|_: checks that phantom, phantomsetup, phantomanalysis and phantommoddump compile with every possible SETUP= flag
- |test|_: runs the test suite [see below]
- |mpi|_: runs the test suite with MPI [see below]
- |mcfost|_: compiles and links phantom+mcfost

.. |build| image:: https://github.com/danieljprice/phantom/actions/workflows/build.yml/badge.svg
.. _build: https://github.com/danieljprice/phantom/actions/workflows/build.yml

.. |test| image:: https://github.com/danieljprice/phantom/actions/workflows/test.yml/badge.svg
.. _test: https://github.com/danieljprice/phantom/actions/workflows/test.yml

.. |mpi| image:: https://github.com/danieljprice/phantom/actions/workflows/mpi.yml/badge.svg
.. _mpi: https://github.com/danieljprice/phantom/actions/workflows/mpi.yml

.. |mcfost| image:: https://github.com/danieljprice/phantom/actions/workflows/mcfost.yml/badge.svg
.. _mcfost: https://github.com/danieljprice/phantom/actions/workflows/mcfost.yml

Running the test suite on your own machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run the test suite using::

   make test

This is just a shortcut for the following sequence of commands::

   make SETUP=test phantomtest
   ./bin/phantomtest

You can run the complete testsuite yourself using the testbot wrapper script::

   cd phantom/scripts
   ./testbot.sh

Running selected parts of the test suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can run just part of the test suite by giving an additional argument
as follows::

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

The buildbot also runs in `an action <https://github.com/danieljprice/phantom/actions>`_ and checks that the code compiles in :doc:`all of
the possible SETUP configurations in the Makefile <setups>`. You can run this
offline as follows::

   cd phantom/scripts
   ./buildbot.sh

If you want to check only those SETUPS that were failing in the actions,
edit the buildbot.sh script and override the allsetups= line, e.g::

   allsetups='disc star'

Common reasons for failure
~~~~~~~~~~~~~~~~~~~~~~~~~~~
We enforce the following policies in merging to the master branch:

1. Code must compile and run with and without DEBUG=yes
2. Code must compile with ability to change the precision of reals to real*4 (this is enforced in SETUP=blob)
3. Code must compile with no warnings when compiled with gfortran (enforced with NOWARN=yes which adds the -Werror flag)
4. Testsuite must work with and without MPI, i.e. compile with MPI=yes

How to reproduce the github build environment offline
======================================================
Just occasionally it is hard to reproduce a failure in the actions. It *is*
possible to recreate the github actions environment offline, using a Docker container.
I suggest to do this *only* as a last resort. The recommended steps are as follows:

Running the actions locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Install `Docker <https://docs.docker.com/desktop/install/mac-install/>`_
2. Install `act <https://github.com/nektos/act>`_
3. run the pull_request workflow

::

   act pull_request

Checking the phantom build that is failing manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you just want to check things manually but in the same environment
as used in the actions, try the following::

1. Install [Docker](https://docs.docker.com/desktop/install/mac-install/)
2. Install Docker command line tools

::

    brew install docker

3. Install the ubuntu-latest image in Docker, e.g. by typing in a terminal [around 6.5Gb download]

::

   docker pull nektos/act-environments-ubuntu:18.04

4. Run the image, and proceed to run phantom build checks manually

::

   git clone https://github.com/danieljprice/phantom
   mkdir -p runs/mydisc
   cd runs/mydisc
   ~/phantom/scripts/writemake.sh disc > Makefile
   export DEBUG=yes
   export PHANTOM_DIR=~/phantom
   make
   make setup
   make analysis
   make moddump
   ./phantomsetup disc
   ./phantomsetup disc
   ./phantomsetup disc
   ./phantom disc
