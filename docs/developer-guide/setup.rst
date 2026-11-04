Initial conditions / writing your own setup routine
---------------------------------------------------

Steps are:

1. check that a :doc:`setup routine for your problem </user-guide/setups>` does not already exist.
2. copy one of the existing setups (e.g. setup_unifdis.f90) to a new
   file (e.g. setup_myproblem.f90).
3. edit this file.
4. set this as the SETUPFILE= within your SETUP block in
   build/phantom/Makefile_setups.
5. compile phantomsetup (“make setup” in the run directory).

Setting up additional particle arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The important arrays that require filling are explicitly passed to the
setpart routine:

::

   subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time)

Other arrays can be setup by importing them via the usual code modules.
For example, to set up sink particles, import the arrays with:

::

   use part, only:xyzmh_ptmass, vxyz_ptmass

and set these to the appropriate values.

Setting default values for runtime options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since phantomsetup is also responsible for (re)writing the initial
runtime options file, it is also possible to provide default options for
runtime parameters in the setup routine. Just import the relevant
parameters from the relevant module:

::

   use options, only:alpha,alphau

These variables can be edited as usual in the .in file at runtime, but
is a good way of giving sensible defaults (note however, that setting
runtime options in your setup routine will cause them to be overwritten
in the .in file every time you run phantomsetup; to prevent this, you
can declare the default values only if the .in file does not exist).

The .setup files
~~~~~~~~~~~~~~~~

Several of the existing setup_*.f90 subroutines write .setup files whose
values are only required only to initialise the problem. For these setup
files, when phantomsetup is run for the first time with a given file
model prefix, then it will interactively ask for the required values; on
all subsequent runs, phantomsetup will read these values from the .setup
file. If options are to be added to these setup files, then they must be
added in the interactive question section, the write_setupfile
subroutine, and the read_setupfile subroutine (all three sections will
be present in the subroutine).

Utility routines
~~~~~~~~~~~~~~~~

There are several routines in the :doc:`libsetup </api/setup/index>` library to assist with setting up the
particle positions:

+---------------------------------------------+-----------------------------------------------+
| module/subroutine                           | purpose                                       |
+=============================================+===============================================+
| :doc:`set_unifdis </api/setup/set_unifdis>` | sets up a uniform particle distribution with  |
|                                             | particles set on various lattice options      |
+---------------------------------------------+-----------------------------------------------+
| :doc:`set_sphere </api/setup/set_sphere>`   | sets up a uniform density sphere              |
|                                             | (interface to set_unifdis)                    |
+---------------------------------------------+-----------------------------------------------+
| :doc:`set_disc </api/setup/set_disc>`       | sets up a single accretion disc with given    |
|                                             | surface density and temperature profiles      |
+---------------------------------------------+-----------------------------------------------+
| :doc:`set_binary </api/setup/set_binary>`   | sets up a binary consisting of two point mass |
|                                             | particles with all 6 orbital elements         |
+---------------------------------------------+-----------------------------------------------+
| :doc:`set_slab </api/setup/set_slab>`       | sets up a thin slab for 2D tests done in 3D   |
+---------------------------------------------+-----------------------------------------------+

Try to use these wherever possible. This makes the setup routines much
simpler and means less cut-and-pasting between otherwise similar setups.
For example, it is perfectly valid to have a setup routine that sets up
a single disc with different parameterisations of Rin, Rout, but you
should achieve this by a simple setup_blah.f90 file that just calls
set_disc.

Internal labels
~~~~~~~~~~~~~~~

Try to refer to particle types and array components by their
parameterised labels in the part module. That is, use

::

   npartoftype(igas) = 1000000

rather than

::

   npartoftype(1) = 1000000

Similarly, you should NEVER set or use iphase directly, instead you
should use the interface routines in the part module (isetphase,
iactive, iamgas, iamdust, iamtype) to set/extract the relevant
information
