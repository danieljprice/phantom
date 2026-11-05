Running phantom with MPI
========================

Phantom uses hybrid MPI and openMP parallelisation: the particles are
divided between MPI tasks, and each task runs several openMP threads.
On a single desktop machine openMP alone is fine. On a cluster, and
always if you want more than one node, compile with MPI=yes and run a
few MPI tasks per node, each with several threads.

Installing MPI
--------------
To compile phantom with MPI, you will need an MPI distribution installed. On your local machine you can install openMPI:

Ubuntu
~~~~~~

::

    sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

Mac OS
~~~~~~

::

    brew install openmpi

On a cluster, load the MPI module that goes with your compiler (e.g. gompi for gfortran with openMPI).

Compiling Phantom with MPI
--------------------------
Follow the usual procedure but with MPI=yes, i.e.:

::

     mkdir mycalc
     cd mycalc
     ~/phantom/scripts/writemake.sh disc > Makefile
     export MPI=yes
     make setup
     make

Running the code with MPI
-------------------------
This is the same as usual but pre-pended by "mpiexec -np 2", where 2 is the number of MPI tasks.
The number of openMP threads per task is set as usual, so the following runs 2 tasks with 4 threads each:

::

   export OMP_NUM_THREADS=4
   mpiexec -np 2 --map-by slot:PE=4 --bind-to core ./phantom disc.in

With openMPI the ``--map-by slot:PE=4`` gives each task four cores; otherwise each task is bound
to a single core and its threads share it.

How many tasks and threads?
---------------------------
- Use 8-16 threads per task, and enough tasks to fill the nodes, e.g. 4 tasks of 16 threads on a 64-core node.
  On OzSTAR's 64-core milan nodes, 4 tasks of 16 threads is faster than openMP alone on the same node,
  and 2-4 nodes are faster again;
- any number of tasks works (it does not need to be a power of 2);
- bind the threads of each task to its cores with::

     export OMP_PROC_BIND=close OMP_PLACES=cores

- each task allocates room for up to four times its share of the particles,
  so with four or more tasks a run needs about four times the memory of the particle arrays in openMP alone;
- with openMP alone on a whole node, spread the memory over all of the node's NUMA domains with::

     numactl --interleave=all ./phantom disc.in

  (``make qscript`` adds this to bash scripts that do not use MPI or a custom ``CMD``). On OzSTAR's milan nodes this is 10-35% faster
  with 64 threads. Do not use it with MPI, where each task's memory is best kept close to its cores.

The particles are divided between the tasks automatically, balancing the work. Forces and conservation are the same as with openMP alone, to round-off.
At each full dump the log file reports the time each task spends in each phase and the data it exchanges with other tasks.

Writing dumps with MPI
----------------------
On a Lustre file system, all MPI tasks write a full dump at the same time where MPI-IO can use direct I/O and phantom was compiled with gfortran or Intel Fortran (ifort or ifx). ``make qscript MPI=yes`` puts the settings for this in the job script; otherwise add them to yours::

    export I_MPI_EXTRA_FILESYSTEM=1              # Intel MPI
    export OMPI_MCA_io_romio321_priority=100     # Open MPI 4.x
    export OMPI_MCA_io_romio341_priority=100     # Open MPI 5.x

It helps most when the run directory is striped over many storage targets, e.g. ``lfs setstripe -c 16 -S 1M .`` before the dumps are written. Otherwise the tasks take turns, as they also do with ``export PHANTOM_DUMP_PARALLEL=no``. The file is the same either way.

So far this has been tested with Open MPI and Intel MPI only. Other MPI libraries (e.g. MPICH, MVAPICH2, Cray MPICH) are untested: they write at the same time only if their MPI-IO reports direct I/O on Lustre, and otherwise take turns.

Running on a cluster
--------------------
On a supercomputing cluster, write your queue script with MPI=yes, giving the number of MPI tasks (NMPI) and
openMP threads per task (NOMP). For example, with slurm, 8 tasks of 16 threads (two 64-core nodes):

::

    make qscript MPI=yes NMPI=8 NOMP=16 INFILE=disc.in > run.q
    sbatch run.q

which runs the code with ``srun --cpus-per-task=16 ./phantom disc.in``. Set MAXMEM to the memory per node you need.
See the guides for :doc:`OzSTAR </getting-started/ozstar>` and :doc:`other clusters </getting-started/clusters>`.
