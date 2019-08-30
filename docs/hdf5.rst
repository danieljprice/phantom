Reading and writing Phantom output files in HDF5 format
=======================================================

From
`Wikipedia <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__:

   Hierarchical Data Format (HDF) is a set of file formats (HDF4, HDF5)
   designed to store and organize large amounts of data. Originally
   developed at the National Center for Supercomputing Applications, it
   is supported by The HDF Group, a non-profit corporation whose mission
   is to ensure continued development of HDF5 technologies and the
   continued accessibility of data stored in HDF.

HDF5 has the following nice features.

-  It is widely available.
-  It has bindings in C and Fortran.
-  It has command line tools for reading data.
-  It has Python packages to read data into NumPy arrays.
-  It has compression built-in.

Compiling and installing HDF5 libraries
---------------------------------------

.. note::

 On most supercomputing clusters HDF5 is available, e.g. Ozstar.
 The following instructions are for a machine you control where the HDF5
 libraries are not yet available.

Using a package manager
~~~~~~~~~~~~~~~~~~~~~~~

On macOS you can install HDF5 with Homebrew.

::

   brew install hdf5

The shared object library and include files are at
``/usr/local/opt/hdf5``. Use this directory as ``HDF5ROOT`` (see below).

On Ubuntu 18.04, for example, you can install HDF5 with apt.

::

   sudo apt install libhdf5-serial-dev

The location of the library is then
``/usr/lib/x86_64-linux-gnu/hdf5/serial``. Use this directory as
``HDF5ROOT`` (see below).

Compiling
~~~~~~~~~

We use the Fortran bindings in Phantom. This means that you must use the
same compiler to compile Phantom as was used for the HDF5 library on
your machine. Typically this is GCC (gfortran) *not* Intel Fortran
(ifort). If you want to compile Phantom with ifort you must compile HDF5
yourself.

First download the source from
https://www.hdfgroup.org/downloads/hdf5/source-code/. Then run the
following code. It may take a while.

::

   # Set compilers to Intel.
   export CC=icc
   export F9X=ifort
   export CXX=icpc

   # Extract tar ball. I assume you downloaded the tar.gz file.
   tar -zxvf hdf5-1.10.5.tar.gz
   cd hdf5-1.10.5

   # Configure and Make
   ./configure --prefix=/usr/local/hdf5 --enable-fortran --enable-cxx
   make

   # Run tests
   make check

   # Install to /usr/local/hdf5
   sudo mkdir /usr/local/hdf5
   # You need ownership over the directory to which you will install
   # Replace <user> appropriately
   sudo chown <user> /usr/local/hdf5
   make install

Compiling Phantom
-----------------

Writing HDF5 output is a compile time option and requires access to the
Fortran HDF5 library. To compile for HDF5 output set ``HDF5_DIR``, for
example if HDF5 was installed with Homebrew on macOS

::

   HDF5_DIR=/usr/local/opt/hdf5

or if it was installed with APT on Ubuntu

::

   HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial

Then compile with

::

   make HDF5=yes HDF5ROOT=$HDF5_DIR

The variable ``HDF5ROOT`` specifies the location of the HDF5 library.

.. note::

 You **may** need to add to ``LD_LIBRARY_PATH`` (in Linux) or
 ``DYLD_LIBRARY_PATH`` (in macOS) to point to the HDF5 library location.
 For example, on Linux with HDF5 compiled with ifort and installed to
 ``/usr/local/hdf5``

 ``export LD_LIBRARY_PATH=/usr/local/hdf5/lib:$LD_LIBRARY_PATH``

Ozstar
~~~~~~

On Ozstar you need to make sure that the OpenMPI and HDF5 modules are
loaded. The variable ``HDF5_DIR`` gives the location of the HDF5 library
once the HDF5 module is loaded.

::

   module load iccifort/2018.1.163-gcc-6.4.0
   module load openmpi/3.0.0
   module load hdf5/1.10.1

Then when you compile Phantom use ``HDF5_DIR`` for ``HDF5ROOT``:

::

   make SYSTEM=ozstar HDF5=yes HDF5ROOT=$HDF5_DIR phantom setup

Note that you must have the HDF5 module loaded when running phantom,
phantomsetup, etc. So make sure to put ``module load hdf5/1.10.1`` in
your Slurm job file.

Converting standard output files to HDF5 format with phantom2hdf5
-----------------------------------------------------------------

``phantom2hdf5`` is a utility that can convert standard Phantom dump
files to HDF5 format.

To compile it, type

::

   make phantom2hdf5 HDF5=yes PHANTOM2HDF5=yes

Recall that you may need to set ``HDF5ROOT`` too.

Now pass a file (or a list of files) to the converter

::

   ./phantom2hdf5 dump_00*

Which returns an HDF5 version of each dumpfile

::

   $ ls
   dump_00000     dump_00001     dump_00002     dump_00003
   dump_00000.h5  dump_00001.h5  dump_00002.h5  dump_00003.h5
   ...

Reading Phantom HDF5 dump files in Python
-----------------------------------------

You can now read the data from the dump file with the command line tools
available with HDF5 or with the Python package h5py.

Command line
~~~~~~~~~~~~

To see all the available datasets:

::

   h5ls -r dump_00000.h5

This produces output like

::

   /                        Group
   /header                  Group
   /header/Bextx            Dataset {SCALAR}
   /header/Bexty            Dataset {SCALAR}
   /header/Bextz            Dataset {SCALAR}
   /header/C_cour           Dataset {SCALAR}
   /header/C_force          Dataset {SCALAR}
   /header/RK2              Dataset {SCALAR}
   /header/alpha            Dataset {SCALAR}
   /header/alphaB           Dataset {SCALAR}
   /header/alphau           Dataset {SCALAR}
   /header/angtot_in        Dataset {SCALAR}
   /header/dtmax            Dataset {SCALAR}
   /header/dum              Dataset {SCALAR}
   /header/etot_in          Dataset {SCALAR}
   /header/fileident        Dataset {SCALAR}
   /header/gamma            Dataset {SCALAR}
   /header/get_conserv      Dataset {SCALAR}
   /header/graindens        Dataset {2}
   /header/grainsize        Dataset {2}
   /header/hfact            Dataset {SCALAR}
   /header/idust            Dataset {SCALAR}
   /header/ieos             Dataset {SCALAR}
   /header/iexternalforce   Dataset {SCALAR}
   /header/isink            Dataset {SCALAR}
   /header/majorv           Dataset {SCALAR}
   /header/massoftype       Dataset {7}
   /header/mdust_in         Dataset {2}
   /header/microv           Dataset {SCALAR}
   /header/minorv           Dataset {SCALAR}
   /header/nblocks          Dataset {SCALAR}
   /header/ndustlarge       Dataset {SCALAR}
   /header/ndustsmall       Dataset {SCALAR}
   /header/npartoftype      Dataset {7}
   /header/nparttot         Dataset {SCALAR}
   /header/nptmass          Dataset {SCALAR}
   /header/ntypes           Dataset {SCALAR}
   /header/polyk2           Dataset {SCALAR}
   /header/qfacdisc         Dataset {SCALAR}
   /header/rhozero          Dataset {SCALAR}
   /header/time             Dataset {SCALAR}
   /header/tolh             Dataset {SCALAR}
   /header/totmom_in        Dataset {SCALAR}
   /header/udist            Dataset {SCALAR}
   /header/umagfd           Dataset {SCALAR}
   /header/umass            Dataset {SCALAR}
   /header/utime            Dataset {SCALAR}
   /header/xmax             Dataset {SCALAR}
   /header/xmin             Dataset {SCALAR}
   /header/ymax             Dataset {SCALAR}
   /header/ymin             Dataset {SCALAR}
   /header/zmax             Dataset {SCALAR}
   /header/zmin             Dataset {SCALAR}
   /particles               Group
   /particles/divv          Dataset {10250000}
   /particles/dt            Dataset {10250000}
   /particles/h             Dataset {10250000}
   /particles/itype         Dataset {10250000}
   /particles/pressure      Dataset {10250000}
   /particles/vxyz          Dataset {10250000, 3}
   /particles/xyz           Dataset {10250000, 3}
   /sinks                   Group
   /sinks/h                 Dataset {4}
   /sinks/hsoft             Dataset {4}
   /sinks/m                 Dataset {4}
   /sinks/maccreted         Dataset {4}
   /sinks/spinxyz           Dataset {4, 3}
   /sinks/tlast             Dataset {4}
   /sinks/vxyz              Dataset {4, 3}
   /sinks/xyz               Dataset {4, 3}

You can access a particular value like

::

   h5dump -d "/header/npartoftype" dump_00000.h5

This produces output like

::

   HDF5 "dump_00000.h5" {
   DATASET "/header/npartoftype" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 7 ) / ( 7 ) }
      DATA {
      (0): 10000000, 250000, 0, 0, 0, 0, 0
      }
   }
   }

Python with h5py
~~~~~~~~~~~~~~~~

The Python package h5py comes with Anaconda. Alternatively you can
install it with pip or Conda.

::

   conda install h5py

To read a dump file

::

   >>> import h5py
   >>> f = h5py.File('dump_00000.h5')

Then you can access datasets like

::

   >>> f['particles/xyz'][:]

Plonk
~~~~~

Plonk is a Python package for analysis and visualisation of Phantom
data. It is open source and available at
https://github.com/dmentipl/plonk.

It uses compiled Fortran from Splash to do interpolation from the data
to a pixel grid.
