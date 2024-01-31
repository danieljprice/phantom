Running phantom on a remote cluster (general instructions)
======================================================================================

The following are some general tips and tricks for working on a remote computing cluster. See the userguide for guidelines specific to particular machines.

Logging in to a remote machine
------------------------------

Once you have an account on the remote machine, with a username, proceed as follows.

In the following I will assume your username is "me101". Login using the "secure shell" command::

   $ ssh -Y me101@gadi.nci.org.au

The '-Y' is to enable forwarding of X-Windows connection (necessary for splash remote visualisation).

Usually I make an alias for ssh so the -Y gets added automatically.
That is, add the following line to your .bashrc file::

   $ alias ssh='ssh -Y'

Finally, it is helpful to not have to type your password every time you log in.
For instructions on how to do this, google for "ssh key exchange" and follow
some of the excellent online tutorials.

Copying files to/from a remote machine
--------------------------------------

Copy files using the scp command. For example, to copy the file "myfile.in" TO the remote machine::

   $ scp myfile.in me101@gadi.nci.org.au:

Pay special attention to the colon, otherwise you will just get a local file called me101@gadi.nci.org.au! If you want to specify the path on the remote host, you
can give the path relative to your home directory after the colon, e.g.::

   $ scp myfile.in me101@gadi.nci.org.au:runs/phantom/mysim/

To copy FROM the remote machine, just issue the same command in reverse::

   $ scp me101@gadi.nci.org.au:runs/phantom/mysim/myfile.in .

where the '.' indicates 'copy to the current directory'.

To ease the process of copying files to/from remote machines, I usually
set up some bash helper functions. To do this add the following lines
to your ~/.bashrc::

 gadiget ()
  {
    scp me101@gadi.nci.org.au:$* .;
  }

  gadiput ()
  {
    scp $* me101@gadi.nci.org.au:;
  }

then you can simply type the following to retrieve files from the remote machine::

   $ gadiget runs/phantom/mysim/myfile.in

or, to copy files onto your home directory on the remote machine::

   $ gadiput myfile.in


Environment variables
---------------------
To compile and run phantom it is necessary to set several
"environment variables". These are variables which are set using
the "export" command in bash, e.g.::

   $ export MYNAME=Daniel

then you can query them using the "echo" command, e.g.::

   $ echo $NAME
    Daniel

You can also list *all* your current environment variable settings
using the "env" command, e.g.::

   $ env
   ...
   NAME=Daniel

Modules
-------
Most supercomputing facilities allow you to load certain compilers or programs that are not installed by default using the "module" command.

You can search for relevant modules (i.e. missing programs) using "module avail", then load them using "module load"::

  $ module avail
   $ module load intel-compiler

Then you can list the modules you have loaded using "module list"::

  $ module list
   Currently Loaded Modulefiles:
    1) pbs   2) intel-compiler/2020.1.217   3) ffmpeg/4.1.3

Setting environment variables and modules every time you log in
----------------------------------------------------------------

To set environment variables and aliases every time you login, simply place the relevant commands in the .bashrc file in your home directory. For example, edit this file using your favourite text editor::

   $ vi ~/.bashrc

Then add the necessary environment variables needed for phantom (and splash) to run correctly. For example, mine contains::

   export SYSTEM=nci
   export OMP_STACKSIZE=512M
   export PATH=$PATH:/scratch/fu7/splash/bin
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/fu7/splash/giza/lib
   source ~/.modules

The .modules file is a file I used to list my "module load" commands::

   vi ~/.modules

where mine currently contains::

   module load intel-compiler
   module load ffmpeg

If you use the "make qscript" functionality in phantom, then your
job submission script will automatically load the relevant modules.

Configuring git on your cluster
-------------------------------
The first time you use git on a particular machine, you need to specify your name and email address::

   $ git config user.name "Daniel Price"
   $ git config user.email "me101@gmail.com"

These are important because they are what is written on your commit log messages when you make commits to the git repo.

Get phantom
-----------

Clone a copy of phantom into your home directory::

   $ git clone https://github.com/danieljprice/phantom

Distinguishing between /home and /scratch filesystems
-------------------------------------------------------
Most supercomputing facilities only give you a small quota
of filespace that is backed up, i.e. your "home" space. Calculations
that output large files are expected to use more temporary disk
space (i.e. "scratch" space). So this means you should check out
and compile code in your home space, but perform the actual
calculations in scratch filespace. I usually configure this by
making some simple shortcuts (using the "ln -s" command).

For example, I make a subdirectory called "runs" from my home directory which is really a shortcut to the /scratch filesystem::

   cd /scratch/fu7
   mkdir $USER
   cd
   ln -s /scratch/fu7/$USER runs
   cd runs
   pwd -P

You can see that it is a shortcut by typing "ls -l" in your home space::

   $ ls -l
   lrwxrwxrwx 1 me101 pt4        19 Dec 13 14:01 runs -> /scratch/fu7/me101

Creating your initial conditions
--------------------------------

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

Submitting a job to the cluster
--------------------------------
To run the code, you need to write a pbs or slurm script. You can get an example by typing “make qscript”::

   $ make qscript INFILE=tde.in JOBNAME=myrun > run.q

should produce something like::

  $ cat run.q
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

  qsub run.q

Check the status using::

  qstat -u $USER

Delete the job using::

  qdel <jobid>

where the *jobid* is the id for the job listed in qstat

Slurm vs pbs, and configuring the phantom Makefile for your cluster
-------------------------------------------------------------------
Slurm is the main alternative to PBS for managing job submission. The functionality is similar but the commands are different. You can configure the "make qscript" command in phantom to use either by making your own "SYSTEM" block in `phantom/build/Makefile_systems <https://github.com/danieljprice/phantom/blob/master/build/Makefile_systems>`__ and specifying QSYS=slurm or QSYS=pbs. For example, add the following lines to `phantom/build/Makefile_systems <https://github.com/danieljprice/phantom/blob/master/build/Makefile_systems>`__::

  ifeq ($(SYSTEM), mycluster)
      include Makefile_defaults_ifort
      QSYS = slurm
      QPROJECT='p01'
  endif

Then set your SYSTEM environment variable appropriately and you should get a slurm script instead of a PBS script::

  $ export SYSTEM=mycluster
  $ make qscript INFILE=tde.in

For slurm the relevant commands are then::

  sbatch run.q

and check status using::

  squeue

Monitoring your job output
--------------------------
The pbs script written by "make qscript" ensures that the output
of your calculation, that would normally go to the screen, instead
is output to a file called "tde01.log", which is automatically updated to "tde02.log" if you restart the calculation. You can monitor the output if the calculation while it is running in the queue using the "tail -f" command::

   $ tail -f tde01.log

Type ctrl-c to quit the "tail -f". Obviously you can also look
at the dump files as they arrive using splash::

   $ splash tde_0*

A common problem is to have forgotten to type "ssh -Y", which will give you the following error::

  Graphics device/type (? to see list, default /xw):
  %giza - ERROR - _giza_open_device_xw: Connection to the X server could not be made
   ERROR opening plotting device

To fix this you need to log out and log back in again using "ssh -Y". If you are using the Windows Linux Subsystem you will also need to have "Xming" running in the background.

Further useful tips
-------------------
It is always useful to improve your knowledge of shell scripting.
For example, I would recommend learning how to use loops, e.g.::

  $ cd phantom; for x in src/main/*.*90; do echo $x; done

and string replacement::

  $ export NAME=Daniel; echo ${NAME/Dan/Span}

You should also learn how to use the "grep" command to search for particular lines of code, e.g.::

  $ grep "call deriv" ~/phantom/src/*/*.*90
