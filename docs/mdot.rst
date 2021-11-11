Computing accretion rates from Phantom data
===========================================

There are several ways to do this.

Use the .ev file output
-----------------------

For a single or binary system the cumulative accreted mass onto the
first two sink particles, or of the two stars if using an external
binary potential, is printed to the .ev file. The advantage of using the
.ev file is that the data is written frequently – every few timesteps or
so, so you get a very fine-grained output compared to just using the
dump files.

To get accretion rates from the .ev data, use the ev2mdot utility. To
compile it, use:

::

   make ev2mdot

which compiles the following program into the phantom/bin directory:

::

   $ cd ~/phantom/bin/
   $ ls ev*
   ev2mdot*

then run this with the .ev file(s) as arguments:

::

   $ ev2mdot disc*.ev
   disc01.ev --> disc.mdot
   disc02.ev --> disc.mdot
   disc03.ev --> disc.mdot
   disc04.ev --> disc.mdot
   disc05.ev --> disc.mdot
   disc06.ev --> disc.mdot
   disc07.ev --> disc.mdot

the .mdot files are regular ascii file with columns as follows:

::

   time mdot macc mdot1 macc1 mdot2 macc2

where mdot is the instantaneous mass accretion rate and macc the
cumulative mass accreted, and the 1 and 2 referring to the values on the
two stars (if present).

ev2mdot assumes the shortest possible time interval by default, but to
get smoother accretion rates you can average over a longer time by
specifying an optional argument:

::

   $ ev2mdot 0.1 disc*.ev
   disc01.ev --> disc01.mdot
   ...

Write an analysis module
------------------------

More generally, to have full access to all of the sink particle
information in the dump files, :doc:`write yourself a module for the
phantomanalysis utility <analysis>`. Then you can just import the sink
particle arrays directly and perform whatever analysis you desire.

You can access the mass accreted by a sink particle by first importing
the main sink arrays:

::

   subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
    use part,     only:xyzmh_ptmass, nptmass, imacc

then get the accreted mass using the imacc tag. For example, the
accreted mass on sink particle 1 would be

::

      mass_accreted_by_sink1 = xyzmh_ptmass(imacc,1)

converting to physical units
----------------------------

To get an accretion rate in Msun/year, import the mass and time units as
follows:

::

    use units, only:umass,utime

multiplying by the relevant unit converts to cgs units. For example, to
get the accreted mass in g we would use

::

      mass_accreted_by_sink1 = xyzmh_ptmass(imacc,1)*umass

to get this in solar masses, just divide by the mass of the Sun

::

      use physcon, only:solarm
      ...
      mass_accreted_by_sink1 = xyzmh_ptmass(imacc,1)*umass/solarm

the same procedure applies to the time scaling, i.e.

::

      use physcon, only:years
      ...
      time_in_years = time_in_code*utime/years
