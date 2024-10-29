Analysis of Phantom output
==========================

Visualisation of Phantom output
-------------------------------

Dump files
~~~~~~~~~~

That's what `splash <https://github.com/danieljprice/splash>`_ is for! Use splash to interactively inspect the dump
files produced by phantom, e.g.:

.. code-block:: bash

   splash dump_0*

or for a column density rendering:

.. code-block:: bash

   splash -r 6 dump_0*

and press Enter for "Hollywood mode". You can also use splash to convert to ascii:

.. code-block:: bash

   splash to ascii dump_0*

or to interpolate to a 3D grid:

.. code-block:: bash

   splash to grid dump_0*

To make a movie, just give "/png" as the output:

.. code-block:: bash

   splash -r 6 dump_0* -dev /png

Then use the ffmpeg script in the splash/scripts directory to convert the png files into an mp4 movie:

.. code-block:: bash

   ~/splash/scripts/movie.sh mymovie.mp4 splash

More details in the `splash documentation <https://splash-viz.readthedocs.io>`_

For more detailed analysis of Phantom dump files, write yourself an
analysis module for the phantomanalysis utility as described below.
Analysis modules exist for many common tasks, including interpolating to
a 3D grid (both fixed and AMR), computing PDFs, structure functions and
power spectra, getting disc surface density profiles, and converting to
other formats, and it is simple to write your own.

Global quantities
~~~~~~~~~~~~~~~~~

The .ev files, which are just ascii files containing
global quantities as a function of time, can be visualised using asplash
with the -ev (or -e) option:

.. code-block:: bash

   splash -ev dump*.ev

Global quantities not in the .ev file can also be obtained using the
splash calc utility, e.g.:

.. code-block:: bash

   splash calc max dump_0*

which produces a file containing the maximum of each quantity in the
dump files as a function of time.

Customised analysis
-------------------

Reading phantom dump files
~~~~~~~~~~~~~~~~~~~~~~~~~~

The short answer is, do not under any circumstances attempt to do this
yourself! (If you need convincing, just have a quick look at how long
the read_data_sphNG.f90 file in splash is). The best way to read/analyse
phantom dumps, aside from using splash to visualise them, is to use the
built-in phantomanalysis utility (described below), or the
:doc:`sarracen </external-utilities/sarracen>` python package. A full description of
the data format and how to read it can be found :doc:`here <dumpfile>`.

Sarracen
~~~~~~~~

- See :doc:`How to analyse and visualise phantom data with sarracen </external-utilities/sarracen>`

Phantomanalysis
~~~~~~~~~~~~~~~

Compile the phantomanalysis utility using:

.. code-block:: bash

   make analysis

which compiles the phantomanalysis binary using the analysis module you
specified in `build/Makefile_setups <https://github.com/danieljprice/phantom/blob/master/build/Makefile_setups>`__:

.. code-block:: make

   ifeq ($(SETUP), isodisc)
       ...
       SETUPFILE= setup_disc.f90
       ANALYSIS= analysis_disc.f90
       ...

giving

.. code-block:: bash

   $ ls
   phantomanalysis*

which you can then run on a series of snapshots

.. code-block:: bash

   $ ./phantomanalysis dump_0*

Phantomanalysis is a simple wrapper that reads all of the dump files on the command line in sequence and calls the analysis routine specified in the ANALYSIS variable, in this case analysis_disc.f90. For a list of pre-built analysis tools, see the :doc:`list of Phantom
utilities <utils>`.

Compiling your own phantomanalysis module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also select the module on the command line using, for example

.. code-block:: bash

   make analysis ANALYSIS=analysis_blah.f90

You can then write an analysis_blah.f90 to do whatever it is you want,
even if you what you want is:

- something completely trivial (see for example `analysis_dtheader.f90 <https://github.com/danieljprice/phantom/blob/master/src/utils/analysis_dtheader.f90>`__ which just compares the time from each dump file with the time in the previous dump file); or
- conversion to another format; or
- actually performing some analysis (e.g. `analysis_disc.f90 <https://github.com/danieljprice/phantom/blob/master/src/utils/analysis_disc.f90>`__ which bins particles into rings for comparison with 1D alpha-disc evolution calculations).

The call to analysis passes the most useful information on the particles
(positions, velocities, thermal energy, particle masses and numbers of
particles). **Any remaining information can also be accessed via the
usual phantom modules**. For example, you can access sink particle
arrays using:

.. code-block:: fortran

   use part, only:xyzmh_ptmass,vxyz_ptmass


Converting to another format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apart from writing a short analysis module, you can also use the convert
utility in splash. For example, to convert all files to ascii format
(not recommended, they’ll be huge):

.. code-block:: bash

   splash to ascii blast_0*

You can also convert to other code formats, e.g.:

.. code-block:: bash

   splash to gadget blast_0*

Analysis with pyanalysis
~~~~~~~~~~~~~~~~~~~~~~~~

An alternative method for analysis in python is to compile the phantom
pyanalysis utility using:

.. code-block:: bash

   make pyanalysis

which compiles the libphantom library giving

.. code-block:: bash

   $ ls
   libphantom.so* libanalysis.py

Now you can import the PhantomAnalysis class from libanalysis.py. This
can be done interactively in iPython or in a Python script

.. code-block:: python

   In [1]: from libanalysis import PhantomAnalysis as pa

and create an instance of this class with a phantom dumpfile

.. code-block:: python

   In [2]: dumpfile = 'blast_00000'

   In [3]: dump = pa(dumpfile)

This loads the dumpfile and places particle quantities into numpy
arrays. These quantities are accessible as attributes of the
PhantomAnalysis class. For example

.. code-block:: python

   In [4]: print dump.npart
   125000

   In [5]: print dump.xyzh
   [[-0.49  -0.47  -0.45  ...,  0.45   0.47   0.49 ]
    [-0.49  -0.49  -0.49  ...,  0.49   0.49   0.49 ]
    [-0.49  -0.49  -0.49  ...,  0.49   0.49   0.49 ]
    [ 0.024  0.024  0.024 ...,  0.024  0.024  0.024]]

   In [6]: print dump.vxyz
   [[ 0.  0.  0. ...,  0.  0.  0.]
    [ 0.  0.  0. ...,  0.  0.  0.]
    [ 0.  0.  0. ...,  0.  0.  0.]]

   In [7]: print dump.utherm
   [ 0.  0.  0. ...,  0.  0.  0.]

List of variables

-  time
-  hfact
-  massofgas
-  units (dictionary) {‘udist’, ‘umass’, ‘utime’, ‘udens’, ‘umagfd’}
-  npart
-  xyzh
-  vxyz
-  utherm
-  nptmass
-  ptmass_xyzmh
-  ptmass_vxyz
-  ptmass_spinxyz

Loading phantom HDF5 dumps into python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To get yourself HDF5 dumpfiles, have a look at :doc:`Running phantom with HDF5 output <hdf5>`.

Import h5py and load the dumpfile

.. code-block:: python

   In [1]: import h5py

   In [2]: f = h5py.File('disc_00000.h5','r')

List the main containers in the file

.. code-block:: python

   In [3]: list(f.keys())
   Out[3]: ['header', 'particles', 'sinks']

List the particle arrays that are available

.. code-block:: python

   In [4]: list(f['particles'].keys())
   Out[4]: ['divv', 'dt', 'h', 'itype', 'pressure', 'vxyz', 'xyz']

Extract the ``xyz`` array from the file

.. code-block:: python

   In [5]: f['particles']['xyz'].value
   Out[5]:
   array([[ -6.05266606,  -6.66164664,  -0.34922808],
          [  2.55540523,  17.91264485,   0.52264339],
          [ 15.26729989,  -6.75512839,  -0.70489168],
          ...,
          [ -9.45331138,   1.34188609,   0.69513828],
          [ 12.67824199,   3.35761305,  -0.39397658],
          [-11.34601204,   0.75837632,   0.6858956 ]])

See `h5py docs <http://docs.h5py.org/en/stable/quick.html>`__ for more information
