Getting started on Monarch (Monash campus cluster)
==================================================

Apply for an account
--------------------

Follow the instructions for how to apply for an account:
https://confluence.apps.monash.edu/display/monarch/How+to+Access+MonARCH

If you are in Daniel Price’s research group, click on “Join existing
project” and use the “pMona0001” project

First time you log in
---------------------

Log in:

::

   $ ssh -Y USERNAME@monarch.erc.monash.edu

show available software

::

   $ module avail

load intel compilers and recent git

::

   $ module load intel
   $ module load git/2.19.0

initialise git-lfs

::

   $ git lfs install

get phantom
~~~~~~~~~~~

Clone a copy of phantom into your home directory

::

   $ git clone https://USERNAME@bitbucket.org/danielprice/phantom.git

get splash
~~~~~~~~~~

Finally, install splash in your home directory by following the
instructions on the `splash home
page <http://users.monash.edu.au/~dprice/splash/>`__

I put the “module load” commands in a file called ~/.modules which
contains the modules I want loaded every time I log in. For example:

::

   $ cat .modules
   module load intel/2018u3
   module load gcc/5.4.0
   module load git/2.19.0

Then, add the following lines to your ~/.bashrc

::

   source ~/.modules
   export SYSTEM=monarch
   ulimit -s unlimited
   export OMP_STACKSIZE=512M
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/splash/giza/lib
   export PATH=${PATH}:${HOME}/splash/bin

Now, when you login again, all these should be set automatically.

Performing a calculation
------------------------

I usually make a soft link / shortcut called “runs” pointing to the
directory where I want to run my calculations:

::

   $ ln -s pMona0001/joebloggs/runs runs
   $ cd runs
   $ pwd -P
   /mnt/lustre/projects/pMona0001/dprice/runs

then make a subdirectory for the name of the calculation you want to run
(e.g. cluster)

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

To run the code, you need to write a slurm script. You can get an
example by typing “make qscript”:

::

   $ make qscript INFILE=shock.in > run.qscript

should produce something like

::

   $ cat run.qscript
   $ make qscript INFILE=shock.in
   #!/bin/bash
   #SBATCH --nodes=1 --ntasks=16
   #SBATCH --cpus-per-task=1
   #SBATCH --job-name=buck-jump
   #SBATCH --account=p01
   #SBATCH --output=shock.in.qout
   #SBATCH --mail-type=BEGIN
   #SBATCH --mail-type=FAIL
   #SBATCH --mail-type=END
   #SBATCH --mail-user=daniel.price@monash.edu
   #SBATCH --time=0-100:59:59
   #SBATCH --mem=16G
   echo "HOSTNAME = $HOSTNAME"
   echo "HOSTTYPE = $HOSTTYPE"
   echo Time is `date`
   echo Directory is `pwd`

   ulimit -s unlimited
   export OMP_SCHEDULE="dynamic"
   export OMP_NUM_THREADS=16
   export OMP_STACKSIZE=1024m

   echo "starting phantom run..."
   export outfile=`grep logfile "shock.in" | sed "s/logfile =//g" | sed "s/\\!.*//g" | sed "s/\s//g"`
   echo "writing output to $outfile"
   ./phantom shock.in >& $outfile

You can then submit this to the queue using

::

   $ sbatch run.qscript
   Submitted batch job 2162704

and check status using

::

   $ squeue
            2162702     medium extended price  R 2-12:13:46      1 hs9
            2162703     medium extended price  R 2-12:13:46      1 hs9
            2162704     medium extended price  R 2-12:13:46      1 hs9

You can follow what the calculation is doing by looking at the .log
file:

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
