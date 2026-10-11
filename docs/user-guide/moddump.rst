Modifying phantom dump files
============================

The right way to edit the particle information in a phantom dump file is
to use the moddump utility. Compile this from your runtime directory
(with a local Makefile written using ~/phantom/scripts/writemake.sh) as
follows::

   make moddump

This builds a binary called phantommoddump, for which you can write a
moddump module to do what you want.

Using built-in moddump modules
------------------------------

Phantom ships with many ready-made modules in ``src/utils/`` (e.g.
``moddump_addflyby``, ``moddump_dustadd``). To use one, specify it when
building::

   make moddump MODFILE=moddump_addflyby.f90

Then run::

   ./phantommoddump dump_in dump_out [time]

Task-oriented recipes (flyby on evolved disc, add dust, extend disc, etc.)
are in :doc:`moddump-recipes`.

Moddump parameter files
-----------------------

Modules with options store them in a ``.moddump`` file named after the
prefix of the OUTPUT dump. For example::

   ./phantommoddump relax_00150 flyby_00000 0.

uses ``flyby.moddump``:

- if the file does not exist, you are prompted for the options, your answers
  are saved to ``flyby.moddump`` and phantommoddump stops, so you can check
  the file before rerunning (as for phantomsetup and .setup files)
- if the file exists, the options are read from it without prompting and the
  dump is modified, so you can edit the file and rerun to repeat or adjust the
  operation
- if options are missing from the file, they are added with default values and
  phantommoddump stops so you can check them before rerunning

The header of the ``.moddump`` file records which dump file was modified and
the time given to the output dump.

Any option in the ``.moddump`` file can also be set on the command line as
``--option=value``, for example::

   ./splitpart dump_01000 split_01000 --nchild=2

Options given this way override the values in the file, and the file is
rewritten with the values that were used. If there is no file yet, the options
are used with the defaults instead of prompting, the file is written and
phantommoddump stops, so run the same command again to modify the dump.

Writing your own moddump module
-------------------------------

The default moddump version is just a “do nothing” operation that you can
edit to do what you want. Make your own copy::

   cd phantom/src/utils/
   cp moddump_default.f90 moddump_mine.f90

If your module has options, read and write them using the ``moddump_utils``
module, in the same way as the ``read_options_*`` and ``write_options_*``
routines in the main code. Call ``get_moddump_options`` at the start of
``modify_dump``, and supply routines to prompt for, write and read the options::

   subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
    use moddump_utils, only:get_moddump_options
    ...
    call get_moddump_options(read_options_moddump,write_options_moddump,moddump_interactive)
    ...
   end subroutine modify_dump

   subroutine moddump_interactive()
    use prompting, only:prompt

    call prompt('Enter angular frequency omega',omega)

   end subroutine moddump_interactive

   subroutine write_options_moddump(iunit)
    use moddump_utils, only:write_inopt
    integer, intent(in) :: iunit

    call write_inopt(omega,'omega','angular frequency of solid body rotation',iunit)

   end subroutine write_options_moddump

   subroutine read_options_moddump(db,nerr)
    use moddump_utils, only:inopts,read_inopt
    type(inopts), intent(inout) :: db(:)
    integer,      intent(inout) :: nerr

    call read_inopt(omega,'omega',db,errcount=nerr)

   end subroutine read_options_moddump

The ``.moddump`` file itself (opening, closing and the header) is handled by
``get_moddump_options``. See ``src/utils/moddump_rotate.f90`` for a complete
example. Modules without options need none of this.

Then specify the name of the moddump module in your SETUP block as
follows::

   ifeq ($(SETUP),blah)
       FPPFLAGS=...
       ...
       MODFILE=moddump_mine.f90

Now remake phantommoddump as above, and run it::

   $  ./phantommoddump
    PhantomSPH: (c) 2007-2026 The Authors

    Usage: moddump dumpfilein dumpfileout [time] [outformat]

So for example to edit the dump file dump_00100 and reset the time in
the modified file to zero, you would use::

   ./phantommoddump dump_00100 newdump_00000 0.0

Handling of input file options
------------------------------

For convenience, phantommoddump also creates a .in file. The rules for
this are:

1. First read the .in file with the prefix of the INPUT file if it
   exists (e.g. dump.in in the above example)
2. Then read the .in file with the prefix of the OUTPUT file in case it
   already exists (e.g. newdump.in in the above example)
3. Write a *new* .in file with the prefix of the OUTPUT file
   (e.g. newdump.in)

So the correct usage for phantommoddump is to perform the moddump
operation **in the same directory as the .in file corresponding to the
prefix of the INPUT file**. This way, the output .in file will be
identical, modified as desired by changes in your moddump module.
Otherwise, you will get a .in file with the code defaults which may be
different to the settings used in your previous calculation.

For other useful tools, see the full :doc:`list of Phantom
utilities <utils>` and :doc:`moddump-recipes`.
