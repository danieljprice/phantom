Getting started on the NCI supercomputer (Australian National Supercomputing Facility)
======================================================================================

Apply for an account at http://nci.org.au

If you are in Daniel Price’s research group, request to join project “fu7”

Log in 
-------

::

   ssh -Y USER@gadi.nci.org.au

Configure your environment
------------------

First edit your .bashrc file in your favourite text editor::

   vi ~/.bashrc

Mine has::

   export SYSTEM=nci
   export OMP_STACKSIZE=512M
   export PATH=$PATH:/scratch/fu7/splash/bin
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/fu7/splash/giza/lib
   export MAXP=2000000
   ulimit -s unlimited
   source ~/.modules
  
If you are using phantom+mcfost, you will need the following lines::
   
   export MCFOST_DIR=/scratch/fu7/mcfost-src/
   export MCFOST_AUTO_UPDATE=0
   export MCFOST_INSTALL=/scratch/fu7/mcfost/
   export MCFOST_GIT=1
   export MCFOST_NO_XGBOOST=1
   export MCFOST_LIBS=/scratch/fu7/mcfost
   export MCFOST_UTILS=/scratch/fu7/mcfost/utils

Then relevant modules in your .modules file::

   vi ~/.modules

Mine contains::

   module load intel-compiler
   module load intel-mpi
   module load git

where the last line is needed for git's large file storage (LFS) to work.

Finally, make a shortcut to the /scratch filesystem::

   cd /scratch/fu7
   mkdir $USER
   cd
   ln -s /scratch/fu7/$USER runs
   cd runs
   pwd -P

Get phantom
-----------

Clone a copy of phantom into your home directory and install git-lfs::

   $ cd $HOME
   $ git clone https://github.com/danieljprice/phantom.git
   $ cd phantom
   $ git lfs install

Run a calculation
------------------
   
then make a subdirectory for the name of the calculation you want to run
(e.g. tde)::

   $ cd; cd runs
   $ mkdir tde
   $ cd tde
   $ ~/phantom/scripts/writemake.sh grtde > Makefile
   $ make setup
   $ make

then run phantomsetup to create your initial conditions::

   $ ./phantomsetup tde
   (just press enter to all the questions to get the default)

To run the code, you need to write a pbs script. You can get an
example by typing “make qscript”::

   $ make qscript INFILE=tde.in JOBNAME=myrun > run.qscript

should produce something like::

  $ cat run.qscript
  #!/bin/bash
  ## PBS Job Submission Script, created by "make qscript" Tue Mar 31 12:32:08 AEDT 2020
  #PBS -l ncpus=48
  #PBS -N myrun
  #PBS -q normal
  #PBS -P fu7
  #PBS -o tde.in.pbsout
  #PBS -j oe
  #PBS -m e
  #PBS -M daniel.price@monash.edu
  #PBS -l walltime=48:00:00
  #PBS -l mem=16G
  #PBS -l other=hyperthread
  ## phantom jobs can be restarted:
  #PBS -r y

  cd $PBS_O_WORKDIR
  echo "PBS_O_WORKDIR is $PBS_O_WORKDIR"
  echo "PBS_JOBNAME is $PBS_JOBNAME"
  env | grep PBS
  cat $PBS_NODEFILE > nodefile
  echo "HOSTNAME = $HOSTNAME"
  echo "HOSTTYPE = $HOSTTYPE"
  echo Time is `date`
  echo Directory is `pwd`

  ulimit -s unlimited
  export OMP_SCHEDULE="dynamic"
  export OMP_NUM_THREADS=48
  export OMP_STACKSIZE=1024m

  echo "starting phantom run..."
  export outfile=`grep logfile "tde.in" | sed "s/logfile =//g" | sed "s/\\!.*//g" | sed "s/\s//g"`
  echo "writing output to $outfile"
  ./phantom tde.in >& $outfile

You can then proceed to submit the job to the queue using::

  qsub run.qscript

Check the status using::

  qstat -u $USER


more info
---------

For more information on the actual machine `read the
userguide <https://opus.nci.org.au/display/Help/Preparing+for+Gadi>`__
