Getting started on the NCI supercomputer (Australian National Supercomputing Facility)
=======================================================================================

Apply for an account at http://nci.org.au

If you are in Daniel Price’s research group, request to join project “wk74”

Log in
-------

Please read the :doc:`general instructions for how to log in/out and copy files to/from a remote cluster <clusters>`
::

   ssh -Y USER@gadi.nci.org.au

Configure your environment
---------------------------

First edit your .bashrc file in your favourite text editor::

   vi ~/.bashrc

Mine has::

   export SYSTEM=nci
   export PROJECT=wk74
   export OMP_STACKSIZE=512M
   export OMP_NUM_THREADS=32
   export SPLASH_DIR=/g/data/$PROJECT/splash
   export PATH=$PATH:$SPLASH_DIR/bin
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SPLASH_DIR/giza/lib
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/hdf5/1.10.5/lib
   export MAXP=2000000
   ulimit -s unlimited
   source ~/.modules

If you are using phantom+mcfost, you will need the following lines::

   export MCFOST_DIR=/g/data/$PROJECT/mcfost-src/
   export MCFOST_AUTO_UPDATE=0
   export MCFOST_INSTALL=/g/data/$PROJECT/mcfost/
   export MCFOST_GIT=1
   export MCFOST_NO_XGBOOST=1
   export MCFOST_LIBS=/g/data/$PROJECT/mcfost
   export MCFOST_UTILS=/g/data/$PROJECT/mcfost/utils
   export HDF5ROOT=/apps/hdf5/1.10.5/lib

Then relevant modules in your .modules file::

   vi ~/.modules

Mine contains::

   module load intel-compiler
   module load intel-mpi
   module load git
   module load hdf5

where the last line is needed for git's large file storage (LFS) to work.

Finally, make a shortcut to the /g/data filesystem::

   cd /g/data/$PROJECT
   mkdir $USER
   cd
   ln -s /g/data/$PROJECT/$USER runs
   cd runs
   pwd -P

Get phantom
------------

Clone a copy of phantom into your home directory::

   $ cd $HOME
   $ git clone https://github.com/danieljprice/phantom.git
   $ cd phantom
   
and tell git who you are::

   $ git config --global user.name "Joe Bloggs"
   $ git config --global user.email "joe.bloggs@monash.edu"


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

   $ make qscript INFILE=tde.in JOBNAME=myrun > run.q

should produce something like::

  $ cat run.q
  #!/bin/bash
  ## PBS Job Submission Script, created by "make qscript" Tue Mar 31 12:32:08 AEDT 2020
  #PBS -l ncpus=48
  #PBS -N myrun
  #PBS -q normal
  #PBS -P wk74
  #PBS -o tde.in.pbsout
  #PBS -j oe
  #PBS -m e
  #PBS -M daniel.price@monash.edu
  #PBS -l walltime=48:00:00
  #PBS -l mem=16G
  #PBS -l storage=gdata/wk74
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

  qsub run.q

Check the status using::

  qstat -u $USER

How to keep your job running for more than 48 hours
----------------------------------------------------

Often you will want to keep your calculation going for longer than the 48-hour maximum queue limit.
To achieve this you can just submit another job with the same script
that depends on completion of the previous job

First find out the job id of the job you have already submitted::

  $ qstat
  Job id                 Name             User              Time Use S Queue
  ---------------------  ---------------- ----------------  -------- - -----
  18780261.gadi-pbs      croc             abc123            402:48:0 R normal-exec

Then submit another job that depends on this one::

   qsub -W depend=afterany:18780261 run.q

The job will remain in the queue until the previous job completes. Then
when the new job starts phantom will just carry on where it left off.

A more sophisticated version of the above can be achieved by generating your
PBS script with PBSRESUBMIT=yes::

  make qscript INFILE=tde.in PBSRESUBMIT=yes > run.q

you can check the details of this using::

  cat run.q

and submit your script using::

  qsub -v NJOBS=10 run.q

which will automagically submit 10 jobs to the queue, each depending on completion of the previous job.

how to not annoy everybody else
---------------------------------
Do not fill the disk quota! Use a mix of small and full dumps where possible and set dtmax to a reasonable value to avoid generating large numbers of unnecessary large files.

For how to move the results of your calculations off gadi see :doc:`here </user-guide/data-curation>`

how to use splash to make movies without your job getting killed
-----------------------------------------------------------------
If you try to make a sequence of images using splash on the login node
(e.g. by typing /png or file.png at the device prompt), your job will get killed
due to the runtime limits::

  Graphics device/type (? to see list, default /xw):/png
  ...
  Killed
  
A simple workaround for this is to launch N instances of splash using a bash loop::

  $ for x in dump_0*; do echo $x; done
  
Then replace "echo $x" with the relevant splash command::

  $ for x in dump_0*; do splash -r 6 -dev $x.png $x; done

If you still get prompts that need answers you can follow the procedure `here <https://splash-viz.readthedocs.io/en/latest/other.html#reading-processing-data-into-images-without-having-to-answer-prompts>`, or simply list the answers to the prompts in a file (here called answers.txt) and use::

  $ for x in dump_0*; do splash -r 6 -dev $x.png $x < answers.txt; done

this way each process is short and your movie-making can proceed without getting killed.

more info
----------
See :doc:`general instructions for how to log in/out and copy files to/from a remote cluster <clusters>`

For more information on the actual machine `read the
userguide <https://opus.nci.org.au/display/Help/Preparing+for+Gadi>`__
