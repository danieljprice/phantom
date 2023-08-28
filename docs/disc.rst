Protoplanetary discs
============================================

We consider the following examples below:

1. Circumbinary disc
2. Flyby interaction of star with protoplanetary disc
3. Protoplanetary disc with dust, gas and planets
4. Self-gravitating disc

Circumbinary disc
------------------

make a new directory and write a local Makefile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure your directory is NOT a subdirectory of the code::

   $ mkdir -p ~/runs/mydisc
   $ cd ~/runs/mydisc
   $ ~/phantom/scripts/writemake.sh disc > Makefile

compile phantom and phantomsetup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   $ make
   $ make setup
   $ ls
   Makefile    phantom*    phantomsetup*

run phantomsetup
~~~~~~~~~~~~~~~~

The setup procedure asks some basic questions to get the structure of
the .setup file correct. Answer 2 for the number of stars::

   ./phantomsetup disc

    -----------------------------------------------------------------

         Welcome to the New Disc Setup

    -----------------------------------------------------------------
     disc.setup not found: using interactive setup

    ===========================
    +++  CENTRAL OBJECT(S)  +++
    ===========================
    Do you want to use sink particles or an external potential?
     0=potential
     1=sinks
     ([0:1], default=1):
    How many sinks? ([1:], default=1): 2
    Do you want the binary orbit to be bound (elliptic) or unbound (parabolic/    hyperbolic) [flyby]?
     0=bound
     1=unbound
     ([0:1], default=0):

    =================
    +++  DISC(S)  +++
    =================
    Do you want a circumbinary disc? (default=yes):
    Do you want a circumprimary disc? (default=no):
    Do you want a circumsecondary disc? (default=no):
    How do you want to set the gas disc mass?
     0=total disc mass
     1=mass within annulus
     2=surface density normalisation
     3=surface density at reference radius
     4=minimum Toomre Q
     ([0:4], default=0):
    Do you want to exponentially taper the outer gas disc profile? (default=no):
    Do you want to warp the disc? (default=no):

    =================
    +++  PLANETS  +++
    =================
    How many planets? ([0:9], default=0):

    ================
    +++  OUTPUT  +++
    ================
    Enter time between dumps as fraction of binary period ([0.000:], default=0.1000):
    Enter number of orbits to simulate ([0:], default=100):

     writing setup options file disc.setup

     >>> please edit disc.setup to set parameters for your problem then rerun  phantomsetup <<<

After answering the questions, this will create a file called sim.setup which contains setup options. Open this file in your favourite text editor to proceed...

edit the .setup file and rerun phantomsetup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After editing the .setup file, run phantomsetup again::

   ./phantomsetup disc

You should see output along the lines of::

   reading setup options from disc.setup
   opening database from disc.setup with 38 entries

   writing setup options file disc.setup

   Central objects represented by two sinks
   Primary mass:         1.00    solarm
   Binary mass ratio:   0.200
   Accretion Radius 1:   1.00    au
   Accretion Radius 2:  0.500    au

    ---------- binary parameters -----------
    primary mass     :   1.00000
    secondary mass   :  0.200000
    mass ratio m2/m1 :  0.200000
    reduced mass     :  0.166667
   ...
   setting ieos=3 for locally isothermal disc around origin
   dust

   ...
   # gas disc parameters - this file is NOT read by setup
                   R_in =         25.    ! inner disc boundary
                  R_ref =         25.    ! reference radius
               R_out =        125.    ! outer disc boundary

   ...
   -------->   TIME =    0.000    : full dump written to file disc_00000.tmp   <--------

    input file disc.in written successfully.
    To start the calculation, use:

    ./phantom disc.in

The above procedure prints a .discparams file (in the above example would be
called disc.discparams) that contains some of the parameters used to
initialise the disc setup.

For a circumbinary disc the equation of state is set to a vertically isothermal equation of state (ieos=3) where the radius is taken with respect to *the coordinate origin*. See :doc:`Equations of state available in Phantom <eos>`

check the .in file and proceed to run phantom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    ./phantom disc.in


Flyby interaction of star with protoplanetary disc
--------------------------------------------------
The procedure for a flyby interaction is as above but answering the questions differently::

   ./phantomsetup flyby

   How many sinks? ([1:], default=1): 2
   Do you want the binary orbit to be bound (elliptic) or unbound (parabolic/hyperbolic) [flyby]?
    0=bound
    1=unbound
    ([0:1], default=0): 1

   =================
   +++  DISC(S)  +++
   =================
   Do you want a circumprimary disc? (default=yes): yes
   Do you want a circumsecondary disc? (default=no): no

which produces::

 writing setup options file flyby.setup

 Central object represented by a sink at the system origin with a perturber sink
   Primary mass:         1.00    solarm
   Perturber mass:       1.00    solarm
   Accretion Radius 1:   1.00    au
   Accretion Radius 2:   1.00    au

  ---------- flyby parameters -----------
  primary mass            :    1.00
  secondary mass          :    1.00
  mass ratio              :    1.00

For a circumprimary disc the equation of state is set to ieos=6, such that the radius is taken with respect to the first sink particle in the simulation. See :doc:`Equations of state available in Phantom <eos>`

The Farris et al. (2014) :doc:`equation of state <eos>` (ieos=14 for a binary or ieos=13 if there are more than two stars) is also useful for a flyby simulation if one does not want to have excessively cold material around the secondary


Protoplanetary disc with embedded planets
-----------------------------------------------
To add planets to a protoplanetary disc simulation, simply amend the line specifying the number of planets you want in the disc::

   # set planets
              nplanets =           3    ! number of planets

and re-run phantomsetup, which will add the missing parameters to the .setup file::

   ./phantomsetup disc

after editing the .setup file, proceed to run phantomsetup again::

   ./phantomsetup disc

and finally proceed to run phantom::

   ./phantom disc.in

Protoplanetary disc with dust, gas and planets
-----------------------------------------------

compile phantom and phantomsetup with SETUP=dustydisc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To add dust you just need to compile with DUST=yes. Rather than having to remember to type 'make DUST=yes' and 'make setup DUST=yes' it's easier to use
the pre-cooked setup configuration *dustydisc* for this::

   ~/phantom/scripts/writemake.sh dustydisc > Makefile
   make setup
   make

run phantomsetup and decide which dust method to use
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The default is to add a single species of dust as a separate set of SPH particles::

   $ ./phantomsetup disc

   ...
   ==============
   +++  DUST  +++
   ==============
   Which dust method do you want? (1=one fluid,2=two fluid,3=Hybrid) ([1:3], default=2):
   Enter total dust to gas ratio ([0.000:], default=0.1000E-01):
   How many large grain sizes do you want? ([1:11], default=1):
   How do you want to set the dust density profile?
    0=equal to the gas
    1=custom
    2=equal to the gas, but with unique cutoffs
    ([0:2], default=0):

setup the desired grain size distribution for multigrain simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simplest method for simulating a range of grain sizes is to use
the dust-as-mixture method, where up to 11 grain sizes are allowed by default, simply edit the .setup file as follows::

   # options for dust
            dust_method =           1    ! dust method (1=one fluid,2=two    fluid,3=Hybrid)
            dust_to_gas =       0.010    ! dust to gas ratio
          ndusttypesinp =          11    ! number of grain sizes
      ilimitdustfluxinp =           T    ! limit dust diffusion using Ballabio et    al. (2018)
             igrainsize =           0    ! grain size distribution (0=log-   space,1=manually)
          igrainsizelog =           0    ! select parameters to fix    (0=smin,smax|1=s1,sN|2=s1,logds|3=sN,logds|4=s1,sN,logds)
                smincgs =   1.000E-04    ! min grain size (in cm)
                smaxcgs =       1.000    ! max grain size (in cm)
                 sindex =       3.500    ! grain size power-law index (e.g. MRN =    3.5)
             igraindens =           0    ! grain density input (0=equal,1=manually)
           graindensinp =       3.000    ! intrinsic grain density (in g/cm^3)
               isetdust =           0    ! how to set dust density profile (0=equal    to gas,1=custom,2=equal to gas with cutoffs)

then run phantomsetup again::

    $ ./phantomsetup disc

The 'limit dust diffusion' makes the simulation inaccurate for the very largest grains but ensures that the simulations do not become prohibitively slow by ensuring that decoupled dust species do not control the simulation timestep. If you want to simulate such species accurately and cheaply you should add these species using separate sets of dust particles.

check the .in file and proceed to run phantom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Assuming setup has completed correctly, you can run phantom as previously::

    $ ./phantom disc.in
