Getting started on DiAL (Data Intensive@Leicester, DiRAC cluster at Leicester)
==============================================================================

Apply for an account
--------------------

Follow the instructions for how to apply for an account:
https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/dirac/request-access

First time you log in
---------------------

Log in with your DiRAC username:

::

   $ ssh -Y USERNAME@dial.dirac.ac.uk

Show available software

::

   $ module avail

Load intel compiler and SPLASH (currently Phantom works best with this
compiler)

::

   $ module load intel/compilers/18.0.3
   $ module load splash

Add the above two commands and the following three commands to your
~/.bashrc:

::

   export SYSTEM=ifort
   ulimit -s unlimited
   export OMP_STACKSIZE=512M

Now, when you login again, all these should be set automatically.

Get Phantom
~~~~~~~~~~~

Clone a copy of Phantom into your home directory:

::

   $ cd ~
   $ git clone https://USERNAME@bitbucket.org/danielprice/phantom.git

Performing a calculation
------------------------

Make a directory (within your user folder) for runs:

::

   $ cd /scratch/PROJECT/USERNAME
   $ mkdir runs
   $ cd runs

then make a subdirectory for the name of the calculation you want to run
(e.g. shock)

::

   $ mkdir shock
   $ cd shock
   $ ~/phantom/scripts/writemake.sh shock > Makefile
   $ make setup
   $ make

then run phantomsetup to create your initial conditions

::

   $ ./phantomsetup shock
   (just press enter to all the questions to get the default)

To run the code, you need to write a submission script. You can get an
example by typing “make qscript”:

::

   $ make qscript INFILE=shock.in > run.qscript

should produce something like

::

   $ cat run.qscript
   #!/bin/bash
   ## PBS Job Submission Script, created by "make qscript" Tue Aug 21 23:17:04 BST 2018
   #PBS -l nodes=1:ppn=36
   #PBS -N [enter a job name here]
   #PBS -A [enter an account number here]
   #PBS -o disc.in.pbsout
   #PBS -j oe
   #PBS -m ea
   #PBS -M [enter your email here]
   #PBS -l walltime=48:00:00
   #PBS -l mem=16G
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

   module load intel/compilers/18
   module load gcc/7.3
   ulimit -s unlimited
   export OMP_SCHEDULE="dynamic"
   export OMP_NUM_THREADS=36
   export OMP_STACKSIZE=1024m


   echo "starting phantom run..."
   export outfile=`grep logfile "disc.in" | sed "s/logfile =//g" | sed "s/\\!.*//g" | sed "s/\s//g"`
   echo "writing output to $outfile"
   ./phantom disc.in >& $outfile

You will need to enter the email destination, job name and account
number as required (without the brackets). You can then submit this to
the queue using

::

   $ qsub run.qscript
   22054.master.cm.cluster

and check status using the qstat command and your username, e.g.

::

   $ qstat -u [your username]
   Job ID                  Username  Queue    Jobname  SessID  NDS   TSK   Memory      Time  S  Time
   22054.master.cm.cluste  [......]  dirac25x dp005    6678     1     36     16gb  01:00:00  Q       

When it has started, you can follow what the calculation is doing by
looking at the .log file:

::

   $ tail -f shock01.log

(press ctrl-c to quit the tail -f command). You should obtain a series
of dump files:

::

   $ ls
   shock_00000
   shock_00001
   shock_00002

which you can view with splash (use the “ssplash” binary to view this
format):

::

   $ ssplash shock_0*

You can also check conserved quantities by plotting things in the .ev
file. The first line of the file shows you what each column is:

::

   $ head shock01.ev

and you can plot these columns using “asplash -ev” or any other program
for plotting ascii files, like gnu plot:

::

   $ asplash -ev *.ev

Acknowledgements
----------------

If you use DiAL in a publication it should be acknowledged with the
following text (from their website):

This work was performed using the DiRAC Data Intensive service at
Leicester, operated by the University of Leicester IT Services, which
forms part of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The
equipment was funded by BEIS capital funding via STFC capital grants
ST/K000373/1 and ST/R002363/1 and STFC DiRAC Operations grant
ST/R001014/1. DiRAC is part of the National e-Infrastructure.

More info
---------

More info is available on the DiAL website:
https://www2.le.ac.uk/offices/itservices/ithelp/services/hpc/dirac
