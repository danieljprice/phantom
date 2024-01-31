Modifying phantom dump files
============================

The right way to edit the particle information in a phantom dump file is
to use the moddump utility. Compile this from your runtime directory
(with a local Makefile written using ~/phantom/scripts/writemake.sh) as
follows:

::

   make moddump

This builds a binary called phantommoddump, for which you can write a
moddump module to do what you want. The default moddump version is just
a “do nothing” operation that you can edit to do what you want. Make
your own copy:

::

   cd phantom/src/utils/
   cp moddump_default.f90 moddump_mine.f90

Then specify the name of the moddump module in your SETUP block as
follows:

::

   ifeq ($(SETUP),blah)
       FPPFLAGS=...
       ...
       MODFILE=moddump_mine.f90

Now remake phantommoddump as above, and run it:

::

   $  ./phantommoddump
    PhantomSPH: (c) 2007-2018 The Authors

    Usage: moddump dumpfilein dumpfileout [time] [outformat]

So for example to edit the dump file dump_00100 and reset the time in
the modified file to zero, you would use:

::

   ./phantommoddump dump_00100 newdump_00000 0.0

Handling of input file options
------------------------------

For convenience, phantommoddump also creates a .in file. The rules for
this are:

1. First read the .in file with the prefix of the INPUT file if it
   exists (e.g. dump.in in the above example)
2. Then read the .in file with the prefix of the OUTPUT file in case it
   already exists (e.g. newdump.in in the above example)
3. Write a *new* .in file with the prefix of the OUTPUT file
   (e.g. newdump.in)

So the correct usage for phantommoddump is to perform the moddump
operation **in the same directory as the .in file corresponding to the
prefix of the INPUT file**. This way, the output .in file will be
identical, modified as desired by changes in your moddump module.
Otherwise, you will get a .in file with the code defaults which may be
different to the settings used in your previous calculation.

For other useful tools, see the full :doc:`list of Phantom
utilities <utils>`.
