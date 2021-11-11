Running phantom with MPI
========================

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
This is the same as usual but pre-pended by "mpiexec -np 4" where 4 is the number of MPI threads

::

   ./mpiexec -np 4 ./phantom disc.in
   
Running on a cluster
--------------------
On a supercomputing cluster, you can amend your queue script as necessary by writing it with MPI=yes. 
For example, on the gadi machine with 2 MPI threads and 48 openMP threads (96 cpus in total) you would use:

::

    make qscript MPI=yes NOMP=48 NMPI=2 > run.q
    qsub run.q
