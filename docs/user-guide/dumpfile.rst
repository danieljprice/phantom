The phantom native binary format
================================

The phantom native file format is a self-describing binary file format that can be written/read in native Fortran code with no dependent libraries. This makes phantom easy to compile. The file format is common to phantom and also the sphNG code by Matthew Bate.

Why you should not try to read this format directly
----------------------------------------------------
It is NOT recommended to delve into gritty details of the file format yourself. Instead, phantom provides several plug-and-play utilities that allow you to perform :doc:`post-simulation analysis <analysis>` or :doc:`modify dump files <moddump>`. These utilities will read the file for you and hand you the necessary information needed.

splash
~~~~~~
For SPH simulations the raw data is not so useful since to visualise
fields in a smooth manner one needs to use the SPH kernel. This is the
purpose of `splash <https://github.com/danieljprice/splash>`_, to enable you to produce smooth plots and visualisations
from the code. It reads the raw data files and gives you plots and visualisations::

  splash file_00000 -r 6 --movie

sarracen
~~~~~~~~
:doc:`Sarracen </external-utilities/sarracen>` is a package with similar functionality to splash but done in Python::

   import sarracen
   sdf = sarracen.read_phantom('file_00000')
   sdf.render('rho')

the sdf object is a pandas dataframe that you can use to access the raw data in the snapshot.

phantomanalysis
~~~~~~~~~~~~~~~~
phantomanalysis is a phantom utility which allows you
to write a plug-in module (e.g. `src/utils/analysis_CoM.f90 <https://github.com/danieljprice/phantom/blob/master/src/utils/analysis_CoM.f90>`__)
to perform the analysis you desire. In your run directory you can compile this as follows::

   make analysis
  ./phantomanalysis file_00000

More details can be found :doc:`here <analysis>`. Other example analysis modules can be found in the `src/utils <https://github.com/danieljprice/phantom/tree/master/src/utils>`__ directory.

phantommoddump
~~~~~~~~~~~~~~~~
phantommoddump is a phantom utility which allows you
to write a plug-in module (e.g. `src/utils/moddump_sink.f90 <https://github.com/danieljprice/phantom/blob/master/src/utils/moddump_sink.f90>`__)
to perform modifications on the data files. In your run directory you can compile this as follows::

   make moddump
  ./phantommoddump file_00000 fileout_00000 0.

More details can be found :doc:`here <analysis>`. Other example moddump plugins can be found in the  `src/utils <https://github.com/danieljprice/phantom/tree/master/src/utils>`__ directory.

showheader
~~~~~~~~~~
There are two simple utilities provided as part of phantom that can show you the basic contents of the output file. These are showheader and showarrays
Compile `showheader <https://github.com/danieljprice/phantom/blob/master/src/utils/showheader.f90>`__ in your run directory as follows::

  $ make showheader
   $ ./showheader file_00000

   FT:Phantom:2021.0.0:63a3980 (hydro+1dust): 11/11/2021 15:52:02.3
    nparttot            10000000
    ntypes                    28
    npartoftype         10000000
    npartoftype                0
    npartoftype                0
    npartoftype                0
    ...

showarrays
~~~~~~~~~~
Compile `showarrays <https://github.com/danieljprice/phantom/blob/master/src/utils/showarrays.f90>`__  in your run directory as follows::

  $ make showarrays
   $ ./showarrays file_00000

   FT:Phantom:2021.0.0:63a3980 (hydro+1dust): 11/11/2021 15:52:02.3
   :: nblocks =   1 array lengths per block =  2
   Block   1 array block   1: size   10000000
    iorig           int*4
    x               real    [  -353.47121609956514        367.94104018200619       -76.374084753256568        19.191506229777517        316.09769022614961       ...
    y               real    [  -84.985590881333579        97.586921713523182       -198.72303085298060        210.31220931816267       -69.913275828622105       ...
    z               real    [   68.441574469373009       -59.149183743236037        2.6592414939047884       -1.4559665587079171       -13.662210030961006       ...
    ...
    Block   1 array block   2: size          2
    x               real    [ -0.25337264432514489        67.435846805220919      ]
    y               real    [  0.18615987303660178       -86.722973955883475      ]


Reading datafiles into Python (recommended way)
------------------------------------------------
The simplest way to read the raw datafiles into Python is to use the
:doc:`sarracen </external-utilities/sarracen>` package:

   - https://github.com/ttricco/sarracen

Reading datafiles into Python (other ways)
------------------------------------------------
Another possibility to read the raw datafiles into Python you can use the
`pre-cooked python script <https://github.com/danieljprice/phantom/blob/master/scripts/readPhantomDump.py>`__ for this::

  phantom/scripts/readPhantomDump.py

where basic usage is:

.. code-block:: python

   dump = read_dump('file_00000')
   print (dump)

You can also use splash to convert to other formats which can be easily read into Python::

  splash to ascii file_00000

And finally you can output rendered pixel maps from splash for final plotting in Python::

  splash -o ascii file_00000 -r 6 -dev /png

which will produce a .pix file containing the raw image (that matches what is shown in splash.png)::

  > writing pixel map to file file_00000_logcolumndensitygcm2_proj.pix ...OK

which can be plotted with::

  python ~/splash/scripts/plot_pix.py file_00000_logcolumndensitygcm2_proj.pix

An alternative to splash with similar functionality in Python is Plonk. The main
limitation is that currently Plonk cannot read the native file format
and requires conversion to hdf5 format first.

A portable Fortran module for reading the datafiles
----------------------------------------------------
If you are STILL not satisfied with the above pre-cooked utilities, we provide
a portable, dependency-free Fortran module (`src/main/utils_dumpfiles.f90 <https://github.com/danieljprice/phantom/blob/master/src/main/utils_dumpfiles.f90>`__)
that can be incorporated into the source code of other codes. This provides
high level functionality that can be used to read the file format.
This library is what is used within Phantom itself to read/write the data files
but is not specific to Phantom in any way.

reading the file header
~~~~~~~~~~~~~~~~~~~~~~~~
The basic operations to read the file header would be:

.. code-block:: fortran

  use dump_utils
  integer :: iu,ierr,nblocks,narraylengths
  character(len=lenid)  :: fileid
  character(len=lentag) :: tag
  type(dump_h)          :: hdr
  integer(kind=8)       :: ilen(4)
  integer               :: nums(ndatatypes,4)
  logical               :: got_x,got_y,got_z,match

  iu = 12    ! the Fortran unit code you want to open the file on

  call open_dumpfile_r(iu,'file_00000',fileid,ierr)

  call read_header(iu,hdr,.true.,ierr)

which returns the file header into a derived type (struct) called hdr. Header
variables (of any type) can then be extracted from the header by calling the
subroutine `extract`:

.. code-block:: fortran

  call extract('nptmass',nptmass,hdr,ierr)

reading arrays
~~~~~~~~~~~~~~
To read arrays efficiently (i.e. in order they were written without skipping
around in the file) you need a certain parsing structure of nested loops but
can then just call the read_arrays routine to extract array information:

.. code-block:: fortran

   call extract('nblocks',nblocks,hdr,ierr)
   call free_header(hdr,ierr)

   read (iu, *) number
   narraylengths = number/nblocks
   narraylengths = nblocks/narraysperblock

   do iblock=1,nblocks
      call read_block_header(narraylengths,ilen,nums,idisk1,ierr)
      do iarr=1,narraylengths
         do k=1,ndatatypes
            do i=1,nums(k,iarr)
               read(iu,*) tag
               call read_array(x,'x',got_x,k,1,ilen,0,iu,tag,match,ierr)
               call read_array(y,'y',got_y,k,1,ilen,0,iu,tag,match,ierr)
               call read_array(z,'z',got_y,k,1,ilen,0,iu,tag,match,ierr)
               if (.not.match) read(iu,*)  ! skip array
            enddo
         enddo
      enddo
   enddo

   close(iu)


reading a single array
~~~~~~~~~~~~~~~~~~~~~~
If efficiency is not your concern, a high level function called `read_array_from_file`
is provided, which can be used to extract a single column from a file:

.. code-block:: fortran

   use dumputils, only:read_array_from_file
   real, allocatable :: x(n)
   integer :: iu

   iu = 12
   call read_array_from_file(iu,'dump_00000','x',x,ierr)


If you still really want to know how the file format works
-----------------------------------------------------------
The basic format is a global header followed by a series of self-describing
blocks. Each block contains a set of arrays of the same length with one
of eight possible data types.

The opening gambit
~~~~~~~~~~~~~~~~~~~
The file is a Fortran binary file. Each 'write' statement in Fortran
writes a 4-byte tag at the beginning and end. In other languages you
will need to read these tags and can use them to decide the length of the line.

The first line consists of magic numbers designed for sanity checking::

  <4 bytes>i1,r1,i2,iversion,i3<4 bytes>

There are three integers: i1, i2 and i3 which should be equal to 060769,
060878 and 690706, respectively. i1 is written in the *default integer* kind,
while r1 is a real number that should be equal to i2 that is written in
the *default real kind*. Successful read of these numbers should be used to
decide whether:

- the file is corrupt or written in the wrong endian(if i1 is not read correctly)
- the default real kind is 4-bytes or 8-bytes (if i2 is read correctly, or not)
- the default integer kind is 4-bytes or 8-bytes (if i3 is read correctly, or not)
- the version of the file format (incremented only if the file format becomes backwards incompatible).

Typically iversion=1 although in (very) old files this
number did not exist and so reading one of the magic numbers here corresponds to a file with iversion=0.

The second line contains a 100-character file identifier::

  <4 bytes>fileid<4 bytes>

Typically in dumps written by phantom this contains code version and date information::

  FT:Phantom:2021.0.0:63a3980 (hydro+1dust): 11/11/2021 15:52:02.3

The first letter of the file id indicates if the file is a 'full dump' (F) or 'small dump' (S).
The second letter (T) indicates the file is written in the 'tagged' format, where printed labels
are written prior to each array being written to the file.

The predefined data types
~~~~~~~~~~~~~~~~~~~~~~~~~
The eight pre-defined data types are, in order:

  i) default integer
  ii) 1-byte integer (integer*1)
  iii) 2-byte integer (integer*2)
  iv) 4-byte integer (integer*4)
  v) 8-byte integer (integer*8)
  vi) default real
  vii) 4-byte real (real*4)
  viii) 8-byte real (real*8)

The 'default integer' and 'default real' are floating precision types
which can be 4-byte or 8-byte. Their type is determined by reading the
magic numbers on the first line of the file. The least used type is
the 2-byte integer, and it is possible that the meaning of this one
may be changed at some stage.

The global header
~~~~~~~~~~~~~~~~~
The global header is a simple loop over the 8 predefined data types, where for each type we write::

 loop i=1,8
    <4 bytes>nvars<4 bytes>
    <4 bytes>tags(1:nvars)<4 bytes>
    <4 bytes>vals(1:nvals)<4 bytes>
 end loop

where:

- 'nvars' is a 4-byte integer
- tags is an array of strings, each tag is 16 characters in length
- vals is an array of variables of the specified type

Typically the maximum number of variables in the header is small,
e.g. no more than 256, although this is not required.

The array blocks
~~~~~~~~~~~~~~~~
Following the global header is a series of array blocks. Each block
contains a series of arrays of the same length that is defined in the
block header. Typically a block would contain a set of arrays
for a subset of the particles (e.g. from an MPI domain decomposition).
For a non-MPI phantom simulation there usually only two
blocks, containing the arrays for "normal" particles and the arrays
for sink particles. For an MPI simulation with N threads these blocks
would be repeated N times.

Importantly, the headers for each block are written prior to the
blocks themselves, to enable efficient memory allocation. So following
the file header the next lines are::

  <4 bytes>nblocks<4 bytes>
  loop i=1,nblocks
       <4 bytes>n(i),nums(1:8,i)<4 bytes>
  end loop

where n is an 8-byte integer containing the array length for each block
and nums(1:8) is a 4-byte integer specifying the number of arrays of
each of the 8 predefined data types that are written in that block.

The above lines are immediately followed by the arrays themselves, written as::

   loop i=1,nblocks
      <4 bytes>tag<4 bytes>
      <4 bytes>integer_array(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>int1_array(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>int4_array(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>real_array(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>real_array2(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>real_array3(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>real4_array(1:n(i))<4 bytes>
      <4 bytes>tag<4 bytes>
      <4 bytes>real8_array(1:n(i))<4 bytes>
   end loop

where tag is a 16-character string containing the label for the array,
n(i) is the array length read from the block header and the type of the
array variables are written in order of their data type. In the above
example the block header would have corresponded to::

   <4 bytes>n(i),1,1,0,1,0,3,1,1<4 bytes>

since we wrote 1 integer array, 1 integer*1, 0 integer*2, 1 integer*4,
0 integer*8, 3 real arrays, 1 real*4 array and 1 real*8 array.

...and that's the data format. But just use the routines :)
