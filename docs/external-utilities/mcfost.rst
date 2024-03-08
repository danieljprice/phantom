MCFOST
======

Using MCFOST to generate synthetic images from Phantom data
-----------------------------------------------------------

About MCFOST
~~~~~~~~~~~~

`MCFOST <http://ipag.osug.fr/~pintec/mcfost/docs/html/index.html>`__ is
a Monte-Carlo radiative transfer code written by Christophe Pinte,
François Menard and Gaspar Duchêne. Christophe Pinte is the lead
developer. MCFOST has a dedicated Phantom interface, meaning it can read
phantom data directly. MCFOST is particularly suited to particle data
because the radiative transfer is performed on a Voronoi mesh which can
be built around the particles, rather than requiring interpolation to a
fixed mesh.

Using MCFOST to produce synthetic images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MCFOST can produce synthetic images similar to what would be seen by
telescopes like ALMA and VLT-SPHERE. It uses the distribution of dust
and gas in Phantom, and model spectra for the illumination from sink
particles, to predict the observed flux. MCFOST can operate in three
modes on Phantom data:

1. **Gas-only calculations**. MCFOST assumes a fixed dust-to-gas ratio
   and a grain model to compute opacities
2. **Gas-dust calculations with single (large) grain species**. MCFOST
   uses gas particle information with a fixed dust-to-gas ratio and a
   grain model to compute the emission/opacity caused by small dust
   grains, and interpolates between this and the distribution of large
   grains to predict flux at longer wavelengths
3. **Multigrain gas-dust calculations**. MCFOST uses the entire grain
   size distribution simulated by Phantom to compute observed flux

Producing dust continuum images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For examples, see `Mentiplay et al.
(2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484L.130M/>`__. See
also `Dipierro et al.
(2015) <http://ui.adsabs.harvard.edu/abs/2015MNRAS.453L..73D>`__ and
`Ragusa et al.
(2017) <http://ui.adsabs.harvard.edu/abs/2017MNRAS.464.1449R>`__ for
earlier examples using RADMC3D

::

   mcfost <para_file> -phantom <phantom_dump>

For more information, read the `MCFOST
documentation <http://ipag.osug.fr/~pintec/mcfost/docs/html/mcfost+phantom.html>`__

Producing synthetic line emission
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MCFOST can use the velocity and density information to compute synthetic
emission in various spectral lines, e.g. 12CO, 13CO, C18O, HCO+ and
more. For examples, see `Pinte et al.
(2018) <http://ui.adsabs.harvard.edu/abs/2018ApJ...860L..13P>`__, `Price et
al. (2018) <http://ui.adsabs.harvard.edu/abs/2018MNRAS.477.1270P>`__

Moment maps
~~~~~~~~~~~

Once you have the cubes, you just need to integrate over velocities to
get the various moments:

-  Moment_0 = int F(x,y) dv
-  Moment_1 = int F(x,y) v dv
-  Moment_2 = int F(x,y) v**2 dv

Producing scattered light images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can produce scattered light images using the -img flag along with
the wavelength of interest. For example, to produce a synthetic image
at 1.6 microns::

   mcfost discA.para -phantom disc_00100 -img 1.6

See examples in `Mentiplay et al.
(2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484L.130M/>`__,
`Nealon et al.
(2019) <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.4951N>`__

Troubleshooting
---------------

**Issue**: When trying to process a dump, get the following error:

::

   Reading phantom density file: disc_00012

   *** ERROR - default real size wrong ***

   Error code =           4
   default real size wrong

**Solution**: To save disk space, phantom outputs both ‘small dumps’ and
‘large dumps’ by default. Only the \`large dumps’ can be read by MCFOST.
By default, these are the ones with filenames ending in 0 (the frequency
is set by the “nfulldumps = 10” parameter in the phantom .in file).

Using MCFOST to set temperatures in a live calculation
------------------------------------------------------

Compiling phantom with MCFOST in a live calculation is a bit more tricky
as you need to have MCFOST compiled from source.

You first need to compile libmcfost:

::

   mkdir ~/mcfost-install
   export MCFOST_INSTALL=~/mcfost-install
   export MCFOST_GIT=1
   cd mcfost/src
   make all

then simply set MCFOST=yes when compiling PHANTOM.

Using Phantom+MCFOST on Ozstar
-------------------------------
There is a copy of mcfost and libmcfost.a compiled in /fred/oz015/cpinte/mcfost

To compile phantom with mcfost on ozstar using this pre-compiled version, you will need::

   module load ifort/2018
   ~/phantom/scripts/writemake.sh mcfost > Makefile
   make setup MCFOST_LIBS=/fred/oz015/cpinte/mcfost/src/libtmp
   make MCFOST_LIBS=/fred/oz015/cpinte/mcfost/src/libtmp
   ./phantomsetup disc
   
To run the code with MCFOST you will need::

   export MCFOST_UTILS=/fred/oz015/mcfost/utils
   ./phantom disc

You will also need a disc.para file

Using Phantom+MCFOST on Mac OS with mcfost installed using homebrew
--------------------------------------------------------------------------
A simple way to install mcfost from source on Mac OS is to use the homebrew package::

  brew tap danieljprice/all
  brew install mcfost

This will install mcfost into /usr/local/bin, libmcfost.a into /usr/local/lib/ 
and mcfost2phantom.mod into /usr/local/include. You can then compile phantom 
linked against MCFOST by overriding the linker flags as follows::

   ~/phantom/scripts/writemake.sh disc > Makefile
   make MCFOST=yes PREFIX=/opt/homebrew LIBCXX=-lc++
   make setup MCFOST=yes PREFIX=/opt/homebrew LIBCXX=-lc++
   ./phantomsetup disc

Or, if you have an intel Mac::

   ~/phantom/scripts/writemake.sh disc > Makefile
   make MCFOST=yes PREFIX=/usr/local LIBCXX=-lc++
   make setup MCFOST=yes PREFIX=/usr/local LIBCXX=-lc++
   ./phantomsetup disc


To run the code with MCFOST you will need to create a directory where MCFOST utilities can be installed::

   mkdir -p ~/mcfost-utils/
   export MCFOST_UTILS=~/mcfost-utils
   ./phantom disc

You will also need a disc.para file

