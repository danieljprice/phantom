Configuring/writing SLURM/PBS/Sun Grid Engine scripts from the Makefile
=======================================================================

To write a SLURM, PBS or Sun Grid Engine script in the run directory,
use

::

   make qscript INFILE=blah.in > run.job

where blah.in is the name of the Phantom input file. This writes a batch
submission script that can be used to submit a Phantom calculation to a
job submission queue:

::

   sbatch run.job

or, using PBS:

::

   qsub run.job

The output of the command above can be configured in the SYSTEM block in
build/Makefile, meaning that you can configure the script according to
the queueing system of the machine that you are running Phantom on.

Here is an example from the Monarch machine, which uses the SLURM
queueing engine:

::

   ifeq ($(SYSTEM), monarch)
       FC= ifort
       ...
       QSYS = slurm
       QPROJECT='p01'
       WALLTIME='100:59:59'
   endif

where we have specified our project ID and the maximum wall time limit
allowed by the local admins.

Below is a second example, showing the configuration for the vayu
machine, which uses the PBS queueing system:

::

   if ($(SYSTEM), vayu)
       FC=ifort
      ...
       QSYS= pbs
       MAIL=daniel.price@monash.edu
       PBSRESUBMIT=yes
       ifneq ($(MPI),yes)
          NOMP=8
          NPAR='8:8'
       endif
   endif

Finally, here is another example from the Monash Sun Grid, which uses
the Sun Grid Engine:

::

   ifeq ($(SYSTEM), msg)
       ...
       QSYS = sge
       QSHELL = tcsh
       MAIL = daniel.price@monash.edu
       ifeq ($(OPENMP),yes)
          QPE = smp
          NOMP = '$$NSLOTS'
          ifndef NPAR
             NPAR = '4-32'
          endif
       endif
       ifeq ($(MPI),yes)
          QPE = mpi
          ifeq ($(OPENMP),yes)
               QPE = mqu4
               NOMP = 4
          endif
       endif
   endif

In the above, $NPAR is the number of requested slots, whereas $NSLOTS is
a variable set by the Sun Grid Engine to be the actual number of cpus
allocated.

A full description of the options are given below:

+-----------------------+-----------------------+-----------------------+
| option                | settings              | description           |
+=======================+=======================+=======================+
| QSYS=pbs              | slurm, pbs or sge     | specifies whether a   |
|                       |                       | SLURM, PBS or Sun     |
|                       |                       | Grid Engine script    |
|                       |                       | should be written.    |
+-----------------------+-----------------------+-----------------------+
| QSHELL=tcsh           | tcsh or bash          | specifies whether or  |
|                       |                       | not the submission    |
|                       |                       | script and runtime    |
|                       |                       | shell is tcsh or bash |
+-----------------------+-----------------------+-----------------------+
| WALLTIME=‘500:00:00’  | dd:hh:mm:ss           | default hard walltime |
|                       |                       | limit for job scripts |
+-----------------------+-----------------------+-----------------------+
| PBSRESUBMIT=yes       | yes or no             | specifies whether or  |
|                       |                       | not to write a        |
|                       |                       | self-resubmitting PBS |
|                       |                       | script (e.g. if the   |
|                       |                       | queue has runtime     |
|                       |                       | limits)               |
+-----------------------+-----------------------+-----------------------+
| NPAR=4-32             | number or range       | specifies requested   |
|                       |                       | number of cpus (fixed |
|                       |                       | for pbs scripts, can  |
|                       |                       | be a range for sge    |
|                       |                       | scripts)              |
+-----------------------+-----------------------+-----------------------+
| NOMP=16               | a number              | specifies number of   |
|                       |                       | openMP threads in job |
|                       |                       | scripts               |
+-----------------------+-----------------------+-----------------------+
| QPE=smp               | string dependent on   | specifies parallel    |
|                       | queue config          | environment in sge    |
|                       |                       | scripts: use qconf    |
|                       |                       | -spl to see the       |
|                       |                       | parallel environments |
|                       |                       | setup for your        |
|                       |                       | queueing system.      |
+-----------------------+-----------------------+-----------------------+
| QNODES=‘nodes=1:ppn=’ | a string              | specifies any node    |
| $(NOMP)               |                       | configuration         |
|                       |                       | specific to the       |
|                       |                       | machine               |
+-----------------------+-----------------------+-----------------------+
