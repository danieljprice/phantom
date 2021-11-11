Tips and tricks for running on Pawsey machines
==============================================

Magnus
------

Use:

::

   export SYSTEM=cray

for the Cray Fortran compiler

When compiling on the head node, sometimes you get

::

   /bin/bash: line 1: 16605 Illegal instruction     ../bin/phantom test

which is caused by different instruction sets (Haswell vs. Sandy Bridge)

::

   module swap craype-haswell craype-sandybridge

Apparently magnus-1 and magnus-2 are Sandybridge while magnus-3 and
magnus-4 are Haswell.

Allinea MAP profiling tools
---------------------------

::

   module load forge

you need to first create the libraries needed to run the map tool:

::

   cd $HOME
   mkdir allinea
   cd allinea
   make-profiler-libraries

then phantom needs to be compiled with some extra flags in
build/Makefile:

::

   FFLAGS+=-G2
   LDFLAGS+= -L$HOME/allinea -lmap ... < and some other stuff >

then to run Allinea on openMP code in the run script need:

::

   export ALLINEA_NO_MPI_AUTODETECT=1
   map --profile aprun -n 1 -d 24 ./phantom disc.in

CrayPAT profiling tools
-----------------------

::

   module load perftools/6.1.0

First compile Phantom using SYSTEM=cray (i.e. using CRAY Fortran
compiler):

::

   make SYSTEM=cray

Then use:

::

   pat_build -o phantom-pat

then submit to queue using the execute command

::

   <rest of job script>
   aprun -n 1 -d 24 ./phantom+pat disc.in

then should output phantom+pat-.xf. Then run

::

   pat_report <name of .xf file> > sampling.txt

This produces sampling.txt, also .apa and .ap2 files. Latter is for
viewing using Apprentice2:

::

   app2 phantom*.ap2

Typical job script for scaling and profiling tests
--------------------------------------------------

::

   #!/bin/bash --login
   #SBATCH --account=director2006
   #SBATCH --nodes=1
   #SBATCH --job-name=scaling
   #SBATCH --output=disc.in.pbsout
   #SBATCH --mail-type=FAIL
   #SBATCH --mail-type=END
   #SBATCH --mail-user=daniel.price@monash.edu
   #SBATCH --time=03:00:00
   #SBATCH --mem=10000M

   #module load intel
   #module load mpt
   #export SYSTEM=ifort

   ulimit -s unlimited
   export OMP_SCHEDULE="dynamic"
   ##export OMP_NUM_THREADS=24
   export OMP_STACKSIZE=1024M
   export LD_LIBRARY_PATH=/scratch/director2006/rnealon/runs/replicate_magnus/allinea:$LD_LIBRARY_PATH
   export ALLINEA_NO_MPI_AUTODETECT=1

   echo "starting phantom run..."
   export outfile=`grep logfile "disc.in" | sed "s/logfile =//g" | sed "s/\\!.*//g" | sed "s/\s//g"`
   echo "writing output to $outfile"

   for x in 1 2 4 12 24; do
      export OMP_NUM_THREADS=$x
      map --profile aprun -n 1 -d $OMP_NUM_THREADS ./phantom blast.in >& $outfile
      cp blast01.log $x\_thread.log
      cp blast.in.s blast.in
   done

Using intel instead of cray Fortran compiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   module swap PrgEnv-cray PrgEnv-intel

Also PrgEnv-gnu for Gnu compilers. This makes ftn->ifort and cc->icc.
Also links various libraries appropriate to the compiler. MPI code also
compiles with just ftn and cc instead of mpif90 and mpicc.

Getting an interactive node
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   salloc -p debugq
