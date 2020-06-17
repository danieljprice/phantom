Analysis of Phantom output
==========================

Visualisation of Phantom output
-------------------------------

Dump files
~~~~~~~~~~

That's what splash is for! Use ssplash to look at the dump
files produced by phantom, e.g.:

.. code-block:: bash

   ssplash dump_0*

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

   asplash -e dump*.ev

where to label the columns properly, set the following environment
variable:

.. code-block:: bash

   export ASPLASH_COLUMNSFILE=~/phantom/scripts/columns

Global quantities not in the .ev file can also be obtained using the
splash calc utility, e.g.:

.. code-block:: bash

   ssplash calc max dump_0*

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
phantomanalysis utility, see below.

Phantomanalysis
~~~~~~~~~~~~~~~

Compile the phantomanalysis utility using:

.. code-block:: bash

   make analysis

which compiles the phantomanalysis binary using the analysis module you
specified in phantom/build/Makefile:

.. code-block:: make

   ifeq ($(SETUP), isodisc)
       ...
       LINKLIST= linklist_cyl.F90
       SETUPFILE= setup_disc.f90
       ANALYSIS= analysis_disc.f90
       ...

giving

.. code-block:: bash

   $ ls
   phantomanalysis*

Phantomanalysis is a simple wrapper that reads all of the dump files on
the command line in sequence and calls the analysis routine specified in
the analysis_blah.f90 module.

You can then write an analysis_blah.f90 to do whatever it is you want,
even if you what you want is: \* something completely trivial (see for
example analysis_dtheader.f90 which just compares the time from each
dump file with the time in the previous dump file); or \* conversion to
another format; or \* actually performing some analysis
(e.g. analysis_disc.f90 which bins particles into rings for comparison
with 1D alpha-disc evolution calculations).

The call to analysis passes the most useful information on the particles
(positions, velocities, thermal energy, particle masses and numbers of
particles). \**Any remaining information can also be accessed via the
usual phantom modules**. For example, you can access sink particle
arrays using:

.. code-block:: fortran

   use part, only:xyzmh_ptmass,vxyz_ptmass

For a list of pre-built analysis tools, see the list of Phantom
utilities :doc:`utils`.

Converting to another format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apart from writing a short analysis module, you can also use the convert
utility in splash. For example, to convert all files to ascii format
(not recommended, they’ll be huge):

.. code-block:: bash

   ssplash to ascii blast_0*

To avoid precision loss, you will need to ensure that splash is compiled
in double precision (use make DOUBLEPRECISION=yes when compiling splash)

Analysis with Python
~~~~~~~~~~~~~~~~~~~~

Compile the phantom pyanalysis utility using:

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
