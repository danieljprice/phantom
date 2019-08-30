imorbel
=======

About
-----

`imorbel <https://github.com/drgmk/imorbel>`__ is a tool for fitting
orbits from observations of short orbital arcs in the plane of the sky.
This can be used to sample trial orbits for accretion discs around
binary stars with phantom.

An example of this is the paper by `Price et al.
(2018) <http://adsabs.harvard.edu/abs/2018MNRAS.477.1270P>`__ where we
used imorbel to generate possible orbits to model the HD142527
circumbinary disc

Using the code
--------------

Download imorbel `from github <https://github.com/drgmk/imorbel>`__. For
example:

::

   $ git clone https://github.com/drgmk/imorbel.git

To run imorbel you need the following python3 packages installed via
your python package manager

::

   corner
   emcee

otherwise you will get errors like

::

   Traceback (most recent call last):
     File "./imorbel.py", line 4, in <module>
       import corner
   ModuleNotFoundError: No module named 'corner'

I installed these with anaconda using “pip install corner” and “pip
install emcee”, but the method depends on how you installed python.

You can generate orbits by entering dates, positions on the sky (in RA
and Dec), masses of the two stars and the distance. For HD142527, with
data from `Lacour et
al. 2016 <http://adsabs.harvard.edu/abs/2016A%26A...590A..90L>`__ the
command we used was:

::

   python3 ./imorbel.py --date 2012-3-11 2013-3-17 2013-4-11 2013-7-14 2014-4-11 2014-5-12 --sep .0897 .082 .0863 .0825 .0797 .0772 --e_sep 0.0026 0.0021 0.0019 0.0011 0.0056 0.0006 --pa 133.1 126.3 126.6 123.8 119.5 116.6 --e_pa 1.9 1.6 1.4 1.2 8.7 0.5 --mass 2.2 --e_mass 0.2 --distance 156 --e_distance 7 --Qmax 90 --nzvz 500 --nelem 10000 --pickle-samples 

To sample orbits interactively simply add the –interactive flag:

::

   python3 ./imorbel.py --interactive --date 2012-3-11 2013-3-17 2013-4-11 2013-7-14 2014-4-11 2014-5-12 --sep .0897 .082 .0863 .0825 .0797 .0772 --e_sep 0.0026 0.0021 0.0019 0.0011 0.0056 0.0006 --pa 133.1 126.3 126.6 123.8 119.5 116.6 --e_pa 1.9 1.6 1.4 1.2 8.7 0.5 --mass 2.2 --e_mass 0.2 --distance 156 --e_distance 7 --Qmax 90 --nzvz 500 --nelem 10000 --pickle-samples 

or with the Claudi+2018 data:

::

   python3 ./imorbel.py --interactive --date 2015-05-13 2015-07-03 2016-03-26 2016-06-13 2017-05-16 2018-04-14 --sep 0.069 0.0655 0.060 0.061 0.0485 0.044 --e_sep  0.002 0.0004 0.002 0.002 0.0005 0.001 --pa 110.2 106.2 97.1 96.3 77.8 55.4 --e_pa 0.5 0.4 0.5 0.5 0.2 0.4 --mass 2.1 --e_mass 0.1 --distance 156 --e_distance 7 --Qmax 90 --nzvz 500 --nelem 10000 --pickle-samples

You could also do it on `Grant Kennedy’s
website <http://camd21.ast.cam.ac.uk/~grant/imorbel/err.php?date=2012-3-11+2013-3-17+2013-4-11+2013-7-14+2014-4-11+2014-5-12&sep=.0897+.082+.0863+.0825+.0797+.0772&e_sep=0.0026+0.0021+0.0019+0.0011+0.0056+0.0006&pa=133.1+126.3+126.6+123.8+119.5+116.6&e_pa=1.9+1.6+1.4+1.2+8.7+0.5&mass=2.2&e_mass=0.2&distance=156&e_distance=8&qmax=90&nzvz=500&nelem=100000>`__
(slightly older version of code), though this won’t give the samples

Using trial orbits in phantom
-----------------------------

For interactively selected orbits, a list of orbital elements will be
printed to the screen:

::

   aeiOwf: 33.6000324392704 0.4853543141484768 119.29326149651105 339.30649589409865 200.12769110596912 35.01567603911597
   aeiOwf: 11.821472499230216 0.37672804809788407 141.12984550892793 306.18736185593275 306.426539881919 236.10163083841906
   aeiOwf: 3280.9876328902697 0.9883468413889207 106.52285476254282 181.57356080235368 67.97524118829074 11.699046192710135
   aeiOwf: 35.33855513483261 0.3114925057433448 116.03817713220647 176.54161957756938 60.19326045361112 11.080508657264627

To use these in phantom simply insert these as rows in a file called
\`orbits.txt’. For example:

::

   $ cat orbits.txt
   31.3582167882 0.739160120129 131.339014199 44.9459443034 27.8888589917 249.265639093
   34.2737927615 0.496467238691 119.354619111 159.201554668 19.9803890039 35.0395196638
   28.8429090404 0.398147851957 120.376279213 340.336584475 201.528494675 33.781760767
   26.4715486983 0.241370019165 119.93067425 349.778360603 218.016179643 25.9373105301
   33.6611161237 0.499347757899 119.709251375 157.599190237 19.2024362133 33.8823996072
   38.9483985735 0.611146033329 120.286389327 19.254039125 354.024378053 268.315878271

then compile phantom with SETUP=disc and proceed to run phantomsetup as
usual

::

   $ mkdir mydisc
   $ cd mydisc
   $ ~/phantom/scripts/writemake.sh disc > Makefile
   $ make setup

giving something like:

::

   $ ./phantomsetup disc
     nprocs =            1

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
   How many sinks? ([1:2], default=1): 2
   Do you want the binary orbit to be bound (elliptic) or unbound (parabolic/hyperbolic) [flyby]?
    0=bound
    1=unbound
    ([0:1], default=0): 0

   =================
   +++  DISC(S)  +++
   =================
   Do you want a circumbinary disc? (default=yes): 

   << press enter for remaining questions >>

Then if orbits.txt exists, simply run phantomsetup again using:

::

   $ cp ../orbits.txt .
   $ ls
   disc.setup  Makefile    orbits.dat  phantomsetup*
   $ ./phantomsetup disc
    writing setup options file discA.setup

    writing setup options file discB.setup

    writing setup options file discC.setup

    writing setup options file discD.setup

    writing setup options file discE.setup

    writing setup options file discF.setup

    writing setup options file discG.setup

    writing setup options file discH.setup

this will automatically read each orbit and generate a separate .setup
file for each orbit. For example:

::

   $ diff discA.setup discB.setup 
   18,23c18,23
   <             binary_a =  31.3582168    ! binary semi-major axis
   <             binary_e =  0.73916012    ! binary eccentricity
   <             binary_i =  131.339014    ! i, inclination (deg)
   <             binary_O =  44.9459443    ! Omega, PA of ascending node (deg)
   <             binary_w =   27.888859    ! w, argument of periapsis (deg)
   <             binary_f =  249.265639    ! f, initial true anomaly (deg,180=apastron)
   ---
   >             binary_a =  38.9483986    ! binary semi-major axis
   >             binary_e =  0.611146033  ! binary eccentricity
   >             binary_i =  120.286389    ! i, inclination (deg)
   >             binary_O =  19.2540391    ! Omega, PA of ascending node (deg)
   >             binary_w =  354.024378    ! w, argument of periapsis (deg)
   >             binary_f =  268.315878    ! f, initial true anomaly (deg,180=apastron)

Initiating a series of calculations for each orbit
--------------------------------------------------

You can then use the run-setups script in phantom to generate
directories for each of the .setup files:

::

   $ ~/phantom/scripts/run-setups.sh disc?.setup
   creating directory discA
   moving setup file discA.setup into discA/
   entering directory discA
   writing discA/Makefile
   writing discA/run.qscript
   ...

giving a series of directories that are ready to go:

::

   $ ls
   Makefile    discA/      discC/      discE/      discG/      orbits.dat
   disc.setup  discB/      discD/      discF/      discH/      phantomsetup*

If you are on a machine where “make qscript” works, each directory
should contain a job submission script, so you can submit all the jobs
using, for example:

::

   $ for x in disc?; do cd $x; sbatch run.qscript; cd ..; done
