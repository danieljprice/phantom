Setting up a softened star
=========================================================

Hydrodynamical simulations of giant stars are computationally demanding, due to
the wide range of dynamical times between the core and the enveleope. Typically,
this problem is alleviated by replacing the stellar core with a point particle
that only interacts gravitationally with the rest of the star via a softened
potential. The stellar profile inside some softening radius is modified to 
preserve hydrostatic equilibrium.

Here, we show how to use Phantom to set up a softened star from a MESA stellar
profile. We name our work directory "star" and make phantomsetup inside:

::

    $ mkdir star
    $ cd star
    $ $PHANTOM_DIR/scripts/writemake.sh star > Makefile
    $ make setup
    $ ls
    Makefile    phantomsetup*

To set up a star called "mystar", run phantomsetup, ``./phantomsetup mystar``,
and answer the prompts. We would like to create a softened star from a MESA
profile, and so we select

::

    Case  5 MESA star from file

and choose the desired EoS when prompted. Currently, stars can only be softened
with (i) Adiabatic/polytropic EoS (ieos = 1), (ii) Ideal gas plus radiation EoS
(ieos = 12), and (iii) MESA EoS (ieos = 10). For (i), you will also be prompted
the adiabatic index, ``gamma``.

You will also be prompted the path to your MESA profile:

::

    Enter file name containing density profile (blank="blank",default="P12_Phantom_Profile.data"):

By default, Phantom looks for the MESA headers "mass_grams", "rho", "cell_specific_IE",
"radius_cm", "pressure", "temperature", "x_mass_fraction_H", and "y_mass_fraction_He",
in any order and reads the data, assumming they are unlogged and in cgs units. If 
the headers of your MESA profile are not in this form, you will have to change them.

After providing the path to your MESA profile, answer "yes" to soften the profile,
and you will be presented with three options:

::

    Soften the core density profile and add a sink particle core? (default=no): yes
    Options for core softening:

    1. Specify radius of density softening
    2. Specify mass of sink particle core (not recommended)
    3. Specify both radius of density softening length and mass
    of sink particle core (if you do not know what you are
    doing, you will obtain a poorly softened profile)
    Select option above : ([1:3], default=1): 

Phantom's core softening procedure requires two parameters: the softening radius, h,
and the core mass, mcore. The softening radius is the radius below which the stellar
profile is replaced with the softened profile. Option 1 is the recommended option,
which allows the user to specify the softening radius, and automatically solves for
a value for the core mass that would give a "nice" softened profile.

The core softening procedure will write a data file containing the softened profile.
Specify the name of the file in the prompt shown below:

::

    Enter output file name of cored stellar profile: (blank="blank",default="mysoftenedstar.dat"): 


If you would like to relax the star automatically, answer yes to the following prompt:

::

    Relax star automatically during setup? (default=no):

At the end of the procedure, a temporary dump file, ``mystar_00000.tmp`` is written.
In this dump, gas particles are mapped to the softened profile, and a sink particle
core with mass mcore and softening length hsoft = h/2 is also added. A setup file,
``mystar.setup`` is also written, which saves the setup options and allows the entire
setup procedure to be performed by running ``./phantomsetup mystar.setup`` without
going through interactive prompting.


Setting up a softened star using the .setup file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The relevant options in the setup file for core softening are

::

    # core softening and sink stellar core options
            isoftcore =           T    ! Soften the core of an input MESA profile
                 ieos =           2    ! 1=isothermal,2=adiabatic,10=MESA,12=idealplusrad
                gamma =       1.667    ! adiabatic index
        isofteningopt =           3    ! 1=supply hsoft, 2=supply mcore, 3=supply both
   unsoftened_profile =  profile51.data   ! Path to MESA profile for softening
       outputfilename =  mysoftenedstar.dat   ! Output path for softened MESA profile
                hdens =         15.    ! Radius of core softening
                mcore =       3.150    ! Mass of sink particle stellar core
                hsoft =       7.500    ! Softening length of sink particle stellar core
            isinkcore =           T    ! Add a sink particle stellar core

In this example,

* The adiabatic EoS has been chosen (ieos = 2), and so the adiabatic index has also
  been specified (gamma = 1.667).

* Both the softening radius and core mass have been specified (isofteningopt = 3).
  The core mass is 3.15 Msun, while the softening radius is 15 Rsun.

* phantomsetup will read the original mesa profile, "profile51.data", and write
  a softened profile, "mysoftenedstar.dat".