Running the dust settling test
==============================

This test is described in `Price & Laibe
2015 <http://ui.adsabs.harvard.edu/abs/2015MNRAS.451.5332P>`__

Set up the problem as usual with SETUP=dustsettle:

::

   mkdir runs/dustsettle
   cd runs/dustsettle
   ~/phantom/scripts/writemake.sh dustsettle > Makefile
   make setup
   make
   make moddump

To run the test you first need to relax the disc atmosphere:

::

   ./phantomsetup relax

I recommend you do NOT add dust at this stage (i.e. set dtg=0). Just run
the relaxation with gas only:

::

   ./phantom relax.in

PL15 recommended you relax the disc for ~15 orbits. By default each dump
file is output every 0.1 orbits, so 15 orbits corresponds to dump file
number 150.

Then use moddump to add dust:

::

   ./phantommoddump relax_00150 disc_00000 0.

and run the dust settling test:

::

   ./phantom disc.in

If you keep the same dtmax as in relax.in then every 10 orbits
corresponds to disc_00100, disc_00200 etc. which are the times shown in
Figures 7-9 of the PL15 paper.
