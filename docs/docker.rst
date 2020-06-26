Running phantom build/tests in a Docker image (as used by bitbucket pipelines)
==============================================================================

To reproduce/debug build results shown on the bitbucket commits page,
you can run Phantom in the same Docker image on your own machine

Install Docker.
---------------

Visit `docker.com <http://docker.com>`__, click “Get Docker” and install
for your OS.

Pull the Docker image.
----------------------

Open a terminal and type:

::

   docker pull conradchan/phantom

Run an interactive session
--------------------------

Type:

::

   docker run -it conradchan/phantom

Install and run Phantom
-----------------------

Follow the :doc:`usual steps <testing>` to run the test suite. Simplest would
be

::

   git clone https://bitbucket.org/danielprice/phantom

followed by

::

   cd phantom
   export SYSTEM=gfortran
   export OMP_STACKSIZE=512M
   make test

To run exactly what is done in the pipelines, look at the commands in
the bitbucket-pipelines.yml file in the root directory of Phantom.

Running pipelines on your own fork of Phantom
---------------------------------------------

To run the pipelines on your own fork, just click “enable pipelines” in
your bitbucket settings.

You will then need to set environment variables, under “Settings”,
“Pipelines”, “Environment variables”

See
https://bitbucket.org/danielprice/phantom/admin/addon/admin/pipelines/repository-variables
