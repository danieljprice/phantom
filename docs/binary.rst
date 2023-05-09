Binary stars and common envelope evolution
============================================

Using SETUP=binary
------------------
The one-stop-shop to setup a binary star simulation is to use SETUP=binary::

   ~/phantom/scripts/writemake.sh binary > Makefile
   make setup
   ./phantomsetup sim

giving::

  -----------------------------------------------------------------
   Welcome to the Ultimate Binary Setup
  -----------------------------------------------------------------

   writing setup options file binary.setup
    Edit sim.setup and rerun phantomsetup


This will create a file called sim.setup which contains setup options. Open this file in
your favourite text editor to proceed...


Two sink particles in orbit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simplest is to setup two sink particles in a binary (iprofile=0), by amending the iprofile flags::

   # options for star 1
              iprofile1 =           0    ! 0=Sink,1=Unif,2=Poly,3=Dens,4=KEPL,5=MESA,6=Pie

   # options for star 2
              iprofile2 =           0    ! 0=Sink,1=Unif,2=Poly,3=Dens,4=KEPL,5=MESA,6=Pie

Then run phantomsetup again to rewrite the required options::

 $ ./phantomsetup sim

which will give::

  ERROR: hacc1 not found
  ERROR: hacc2 not found
   2 error(s) during read of setup file: re-writing...
  writing setup options file sim.setup
   Edit sim.setup and rerun phantomsetup

Run this again to finish the setup::

   ./phantomsetup sim

This creates a sim.in file for the main code. You should edit the tmax and dtmax
to give the desired finishing time (tmax) and time between snapshots (dtmax)::

   cat sim.in

               dumpfile =  sim_00000.tmp  ! dump file to start from
   ...
                   tmax =     200000.    ! end time
                  dtmax =         50.    ! time between dumps
   ...
                  alpha =       0.000    ! MINIMUM art. viscosity parameter (max = 1.0)
                 alphau =       1.000    ! art. conductivity parameter
   ...
          icreate_sinks =           0    ! allow automatic sink particle creation
                  f_acc =       0.500    ! particles < f_acc*h_acc accreted without checks


After editing the .in file, proceed to run the simulation::

   make
   ./phantom sim

and have a look at the outputs with splash::

   splash sim_0*


Two polytropes in orbit
~~~~~~~~~~~~~~~~~~~~~~~~~
The default is to setup-and-relax two polytropes and place them in orbit. For this edit
your sim.setup file to give::

  # options for star 1
             iprofile1 =           2    ! 0=Sink,1=Unif,2=Poly,3=Dens,4=KEPL,5=MESA,6=Pie
                Mstar1 =       1.000    ! mass of star1
                Rstar1 =       1.000    ! radius of star1
                   np1 =        1000    ! number of particles

  # options for star 2
             iprofile2 =           2    ! 0=Sink,1=Unif,2=Poly,3=Dens,4=KEPL,5=MESA,6=Pie
                Mstar2 =       1.000    ! mass of star2
                Rstar2 =       1.000    ! radius of star2

Then run phantomsetup again to rewrite the required options::

 $ ./phantomsetup sim

This time you should see the automated relax-a-star procedure kick in::

    RELAX-A-STAR-O-MATIC: Etherm:  0.463     Epot: -0.822     R*:   1.00
       WILL stop WHEN: dens error <   1.00% AND Ekin/Epot <   1.000E-07 OR Iter=1000
    Relaxing star: Iter   1/1000, dens error: 13.72%, R*:  0.924     Ekin/Epot:  3.398E-03
    Relaxing star: Iter   2/1000, dens error: 10.97%, R*:  0.915     Ekin/Epot:  4.673E-03
    ...

As previously, you can then just proceed to run the simulation after editing the sim.in file
::

   ./phantom sim.in

Two stars from MESA profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To use stellar profiles from the MESA code, select iprofile=5 in your .setup file
and enter the name of the ascii data file containing the input profile
(most files produced by MESA just work...)::

  # options for star 1
            iprofile1 =           5    ! 0=Sink,1=Unif,2=Poly,3=Dens,4=KEPL,5=MESA,6=Pie
       input_profile1 =  P12_Phantom_Profile.data   ! Path to input profile
           isoftcore1 =           0    ! 0=no core softening, 1=cubic, 2=const. entropy
           isinkcore1 =           F    ! Add a sink particle stellar core
                  np1 =     1000000    ! number of particles

  # options for star 2
            iprofile2 =           5    ! 0=Sink,1=Unif,2=Poly,3=Dens,4=KEPL,5=MESA,6=Pie
       input_profile2 =  P12_Phantom_Profile.data   ! Path to input profile
           isoftcore2 =           0    ! 0=no core softening, 1=cubic, 2=const. entropy
           isinkcore2 =           F    ! Add a sink particle stellar core

Notice that you do not get to set the particle resolution for the second star,
since the mass of the particles is fixed by the mass and particle number in star 1.

Replacing dense stellar cores with sink particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the options above you have the option to remove the dense core of the star
which causes small timesteps in the code, and replace it with a softened point
mass. The default option for this is isoftcore1=2 and isinkcore1=1.

For more details, see :doc:`Setting up a softened star <softstar>`


Using SETUP=star and moddump_binary
------------------------------------
See :doc:`Setting up stars and tidal disruption events <star>` for the older two-step procedure. The options available are
identical, but with a bit more flexibility and without
having to re-run the relaxation procedure over and over again.
