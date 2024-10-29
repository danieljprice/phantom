PhantomNR
=========

Using PhantomNR to simulate general relativistic hydrodynamics on dynamical spacetimes 
--------------------------------------------------------------------------------------

About phantomNR
~~~~~~~~~~~~~~~

`phantomNR <https://github.com/spencermagnall/phantomNR>`__ is
an extension to the General Relativistic Smoothed Particle Hydrodynamics code Phantom,
that allows for the evolution of relativistic fluids with evolving spacetime metrics. 
This is acomplished via coupling with the numerical relativity framework Einstein Toolkit (ET).
phantomNR's current usage is as a fully relativistic  N-Body code for the simulation of inhomogenous
cosmologies (see `Magnall et al. 2023 <https://ui.adsabs.harvard.edu/abs/2023arXiv230715194M/abstract>`__). 
Einstein Toolkit acts as a "driver" for both the spacetime evolution, and the hydrodynamic evolution. 
As a consquence, simulations are started and mointered entirely within ET, and are setup using a .par 
parameter file which describes the parameters of the simulation. In addition, phantomNR also requires
particle information, which is provided via the standard phantom dump file.   


Compilation and linking
~~~~~~~~~~~~~~~~~~~~~~~ 
You will first need to compile phantom and phantomsetup 
using the flrw setup

::

   scripts/writemake.sh flrw > Makefile 

   make; make setup

which compiles the libphantom.a static library which is 
required for linking and the phantom and phantomsetup binaries. 

You will also need to set the include directory of phantom in Einstein Toolkit
e.g: 

::

   PHANTOM_DIR = /Users/smag0001/phantom/phantomET/bin

Generating a phantom dump file from phantom setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Particles can be setup using phantomsetup in two ways:

1. **Using a regular .setup file** 
   e.g ./phantomsetup flrw.setup will produce a dump file and .in file using an interactive setup routine.
   

3. **Using a .par file** By appending .setup to the end of an Einstein Toolkit parameter file, phantomsetup 
   will automatically read in (most) relevant quantities about the simulation setup and generate an appropriate
   distribution of particles


Troubleshooting
---------------

**Issue**: Large Constraint Violations



**Solution**: Generally, this is indicative of a mismatch between the spacetime setup by Einstein Toolkit
and the particle distribution which is setup by Phantom. A large raw constraint violation, may not always be indicative
of a poorly initialised setup however. It is important to check the relative constraint violations (TODO insert equations)

In many cases, a poor initial constraint is simply a consquence of not setting spacetime and consistently (e.g phi=1e-4 for particles, but phi=1e-6 for spacetime).
We reccomend that the .in and dumpfiles are generated using the .par file of Einstein Toolkit to alleviate this issue.

Constraint violations may also occur due to a low particle and/or grid resolution 




Using phantomNR on Ozstar/NT 
-------------------------------


