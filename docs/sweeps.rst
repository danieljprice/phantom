Performing parameter sweeps the easy way
========================================

See `phantom-config <https://github.com/dmentipl/phantom-config>`__ for
a Python package to parse, convert, modify, and generate Phantom config
files, i.e. .in and .setup files.

Loops in the .in file
---------------------

Phantom input files have some nice features, namely a “loop syntax” that
can be used to perform parameter sweeps on ANY real or integer variable
defined in the input file. For example, to perform a set of calculations
using different values of the artificial viscosity parameter, edit the
relevant line in the .in file (in this example called blast.in) from

::

                  alpha =      1.0000    ! art. viscosity parameter

to

::

                  alpha = 0.1 to 1.0 step 0.1    ! art. viscosity parameter

then run either phantom or phantomsetup on the modified input file:

::

   ./phantom blast.in

and you get the following output:

::

    __
   /  \
   |  LOOP DETECTED IN VARIABLE alpha in input file on line 29
   \--->
   WRITING   10 INPUT FILES VARYING alpha
   ---------------blast_alpha_0.10.in---------------
                  alpha =      0.1000    ! art. viscosity parameter
   ---------------blast_alpha_0.20.in---------------
                  alpha =      0.2000    ! art. viscosity parameter
   ---------------blast_alpha_0.30.in---------------
                  alpha =      0.3000    ! art. viscosity parameter
   ---------------blast_alpha_0.40.in---------------
                  alpha =      0.4000    ! art. viscosity parameter
   ---------------blast_alpha_0.50.in---------------
                  alpha =      0.5000    ! art. viscosity parameter
   ---------------blast_alpha_0.60.in---------------
                  alpha =      0.6000    ! art. viscosity parameter
   ---------------blast_alpha_0.70.in---------------
                  alpha =      0.7000    ! art. viscosity parameter
   ---------------blast_alpha_0.80.in---------------
                  alpha =      0.8000    ! art. viscosity parameter
   ---------------blast_alpha_0.90.in---------------
                  alpha =      0.9000    ! art. viscosity parameter
   ---------------blast_alpha_1.00.in---------------
                  alpha =      1.0000    ! art. viscosity parameter

and you should have a sequence of input files in the current directory,
including the original “master” file:

::

   $ ls *.in
   blast.in        blast_alpha_0.30.in blast_alpha_0.60.in blast_alpha_0.90.in
   blast_alpha_0.10.in blast_alpha_0.40.in blast_alpha_0.70.in blast_alpha_1.00.in
   blast_alpha_0.20.in blast_alpha_0.50.in blast_alpha_0.80.in

To set up a sequence of calculations, use the stagerun.sh script in the
phantom/scripts directory:

::

   $ ~/phantom/scripts/stageruns.sh blast_alpha*.in

This creates a directory corresponding to each input file, writes a
local Makefile in each directory using the writemake.sh script, compiles
phantom in each directory, and even writes a pbs/sge job submission
script using the customisable “[[qscript|make qscript]]” functionality
in the Phantom Makefile (customise the output of this script via the
SYSTEM block in build/phantom/Makefile).

Loop syntax reference
---------------------

The example given above shows the basic loop syntax with linear steps.
You can also use logarithmic steps:

::

                  alpha = 0.1 to 1.0 logstep 0.1    ! art. viscosity parameter

which increments alpha in logarithmic intervals:

::


   /  \
   |  LOOP DETECTED IN VARIABLE alpha in input file on line 29
   \--->
   WRITING   11 INPUT FILES VARYING alpha
   ---------------blast_alpha_0.10.in---------------
                  alpha =      0.1000    !  art. viscosity parameter
   ---------------blast_alpha_0.13.in---------------
                  alpha =      0.1259    !  art. viscosity parameter
   ---------------blast_alpha_0.16.in---------------
                  alpha =      0.1585    !  art. viscosity parameter
   ---------------blast_alpha_0.20.in---------------
                  alpha =      0.1995    !  art. viscosity parameter
   ---------------blast_alpha_0.25.in---------------
                  alpha =      0.2512    !  art. viscosity parameter
   ---------------blast_alpha_0.32.in---------------
                  alpha =      0.3162    !  art. viscosity parameter
   ---------------blast_alpha_0.40.in---------------
                  alpha =      0.3981    !  art. viscosity parameter
   ---------------blast_alpha_0.50.in---------------
                  alpha =      0.5012    !  art. viscosity parameter
   ---------------blast_alpha_0.63.in---------------
                  alpha =      0.6310    !  art. viscosity parameter
   ---------------blast_alpha_0.79.in---------------
                  alpha =      0.7943    !  art. viscosity parameter
   ---------------blast_alpha_1.00.in---------------
                  alpha =      1.0000    !  art. viscosity parameter

You can also specify the number of steps rather than the step interval:

::

                  alpha = 0.1 to 1.0 nstep 3    ! art. viscosity parameter

producing:

::

   /  \
   |  LOOP DETECTED IN VARIABLE alpha in input file on line 29
   \--->
   WRITING    3 INPUT FILES VARYING alpha
   ---------------blast_alpha_0.10.in---------------
                  alpha =      0.1000    !  art. viscosity parameter
   ---------------blast_alpha_0.55.in---------------
                  alpha =      0.5500    !  art. viscosity parameter
   ---------------blast_alpha_1.00.in---------------
                  alpha =      1.0000    !  art. viscosity parameter

And similarly the number of logarithmic steps:

::

                  alpha = 0.1 to 1.0 nlogstep 3    ! art. viscosity parameter

giving:

::

    __
   /  \
   |  LOOP DETECTED IN VARIABLE alpha in input file on line 29
   \--->
   WRITING    3 INPUT FILES VARYING alpha
   ---------------blast_alpha_0.10.in---------------
                  alpha =      0.1000    !  art. viscosity parameter
   ---------------blast_alpha_0.32.in---------------
                  alpha =      0.3162    !  art. viscosity parameter
   ---------------blast_alpha_1.00.in---------------
                  alpha =      1.0000    !  art. viscosity parameter

Integer variables
-----------------

Loop syntax also works with integer quantities in the input file. For
example:

::

         ishock_heating =  0 to 1 step 1  ! shock heating (0=off, 1=on)

produces:

::

    __
   /  \
   |  LOOP DETECTED IN VARIABLE ishock_heating in input file on line 35
   \--->
   ---------------blast_ishock_heating_0.in---------------
         ishock_heating =           0    !  shock heating (0=off, 1=on)
   ---------------blast_ishock_heating_1.in---------------
         ishock_heating =           1    !  shock heating (0=off, 1=on)
