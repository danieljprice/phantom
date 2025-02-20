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

Retrieve a sample parameter file using::

   mcfost -get_para

Copy it to a more sensible filename::

   cp ref4.1.para disc.para 

Then edit the parameter file::

   nano disc.para

In particular you should pay attention to the following parameters::

   #Number of photon packages
   1.28e7                  nbr_photons_eq_th  [set to at least 10 times the number of SPH particles]
   ...

Then proceed, specifying the phantom dump file name you want to make images of::

   mcfost disc.para -phantom disc_00100

The first step computes the temperature, so you then need to make an image at the 
desired wavelength in micron. For example, for a wavelength of 1mm one would use::

   mcfost disc.para -phantom disc_00100 -img 1000

For more information, read the `MCFOST
documentation <http://ipag.osug.fr/~pintec/mcfost/docs/html/mcfost+phantom.html>`__

Producing synthetic line emission
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MCFOST can use the velocity and density information to compute synthetic
emission in various spectral lines, e.g. 12CO, 13CO, C18O, HCO+ and
more. For examples, see `Pinte et al.
(2018) <http://ui.adsabs.harvard.edu/abs/2018ApJ...860L..13P>`__, `Price et
al. (2018) <http://ui.adsabs.harvard.edu/abs/2018MNRAS.477.1270P>`__

The basic procedure is to add the -mol flag::

   mcfost disc.para -phantom disc_00100 -mol

which with default settings should produce a cube of CO line emission.

Moment maps
~~~~~~~~~~~

Once you have the cubes, you just need to integrate over velocities to
get the various moments:

-  Moment_0 = int F(x,y) dv
-  Moment_1 = int F(x,y) v dv
-  Moment_2 = int F(x,y) v**2 dv

This is best done using `pymcfost <https://github.com/cpinte/pymcfost>`__

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


Compiling and running Phantom+MCFOST on Ozstar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Compiling and running Phantom+MCFOST on Mac OS with mcfost installed using homebrew
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


Runtime options for phantom+MCFOST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, when using MCFOST, you should NOT let the temperature evolve between MCFOST calls, hence
the following options should be switched off when MCFOST is activated::

              ipdv_heating =          0    ! heating from PdV work (0=off, 1=on)
           ishock_heating =           0    ! shock heating (0=off, 1=on)

This is because we assume radiative equilibrium at all times, so the temperature is set by the
balance between heating and cooling, and this is computed by MCFOST, not phantom. The temperature
is updated every dtmax.

After compiling phantom+MCFOST as above, you should also find several new options appearing in the .in file::

         use_mcfost =            T    ! use the mcfost library

use this option to switch the call to MCFOST on or off. Beware that the code is compiled with energy
STORED, so running with use_mcfost = F will revert to an ADIABATIC equation of state, but where u=const
on particles if ipdv_heating and ishoc_heating are off (this is not the same as the locally isothermal
equation of state used in normal simulations of discs).

::

    use_mcfost_stars =           F    ! Fix the stellar parameters to mcfost values or update using sink mass

either use the stellar spectra in the MCFOST .para file, or look up spectra based on Siess+2000 isochrones
based on the mass of each sink particle. You should manually set the stellar parameters in the .para file
if you are trying to model a known source (e.g. HD 142527).

::

    mcfost_computes_Lacc =           F    ! Should mcfost compute the accretion luminosity

Accretion luminosity adds an additional radiation source based assuming mass accreted by each sink particle
is converted into radiation on the stellar surface. This is emitted as a blackbody with temperature set by dividing
the accretion luminosity by 4*pi*R^2, where R is the stellar radius (set in the .para file).

::

     mcfost_uses_PdV =           T    ! Should mcfost use the PdV work and shock heating?

The only source of photons in MCFOST by default is from stars (ie. sink particles). If you want to include heating
from shocks and PdV work, you should set this to T. This will add the pdV work and shock heating as source terms
in the Monte Carlo radiative transfer. Recall that when using MCFOST we are assuming radiative equilibrium at
all times, so the temperature is set by the balance between heating and cooling. See Figure A1 in 
`Borchert et al. 2022b <https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.4436B>`__ for an example of the effect 
of PdV work and shock heating on the temperature structure of a disc. Typically it is small.

::

    mcfost_keep_part =       0.999    ! Fraction of particles to keep for MCFOST

MCFOST throws away very distant particles by default when constructing the Voronoi mesh. Set this to 1.0 to keep all particles.

::

                 ISM =           0    ! ISM heating : 0 -> no ISM radiation field, 1 -> ProDiMo, 2 -> Bate & Keto

include additional source of UV from the background interstellar medium, so there is some low temperature even if
no sink particles are present in the simulation

::

    mcfost_dust_subl =           F    ! Should mcfost do dust sublimation (experimental!)

attempts to remove dust in regions where the temperature exceeds the sublimation temperature (1500K). This is experimental.
