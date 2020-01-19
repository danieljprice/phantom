How to run Phantom on OzStar
============================

See also older instructions for :doc:`gstar <g2>`.

Apply for an account
--------------------

https://supercomputing.swin.edu.au/account-management/new_account_request

If you are in Daniel Price’s research group, from your account
management page, request “join project” and select “oz015 - Price/Pinte
research group”

First time you log in
---------------------

::

   $ ssh -Y USERNAME@ozstar.swin.edu.au

show available software

::

   $ module avail

load intel compilers, git, git-lfs and splash

::

   $ module load ifort/2016.2.181-gcc-6.4.0
   $ module load git/2.16.0
   $ module load git-lfs/2.4.0
   $ module load splash/2.7.0

Get phantom
~~~~~~~~~~~

Clone a copy of phantom into your home directory

::

   $ git clone https://USERNAME@bitbucket.org/danielprice/phantom.git

Set your username and email address
-----------------------------------

Ensure that your name and email address are set, as follows:

::

   cd phantom
   git config --global user.name "Joe Bloggs"
   git config --global user.email "joe.bloggs@monash.edu"

Please use your full name in the format above, as this is what appears
in the commit logs (and in the AUTHORS file).

Initialise git-lfs
------------------

Ensure you have done “module load git-lfs” as above, then type:

::

   git lfs install

edit your .bashrc file
----------------------

I put the “module load” commands in a file called ~/.modules which
contains the modules I want loaded every time I log in. For example:

::

   $ cat .modules
   module load ifort/2016.2.181-gcc-6.4.0
   module load git/2.16.0
   module load git-lfs/2.4.0
   module load splash/2.7.0

Then, add the following lines to your ~/.bashrc

::

   source ~/.modules
   export SYSTEM=ozstar
   ulimit -s unlimited
   export OMP_STACKSIZE=512M

Now, when you login again, all these should be set automatically.

Performing a calculation
------------------------

You should *not* perform calculations in your home space - this is for
code and small files. Calculations should be run in the “project” area
in /fred/PROJECT_NAME/USERNAME

I usually make a soft link / shortcut called “runs” pointing to the
directory where I want to run my calculations:

::

   $ cd /fred/oz015
   $ mkdir USERNAME
   $ cd
   $ ln -s /fred/oz015/USERNAME runs
   $ cd runs
   $ pwd -P
   /fred/oz015/USERNAME

then make a subdirectory for the name of the calculation you want to run
(e.g. shock)

::

   $ mkdir shock
   $ cd shock
   $ ~/phantom/scripts/writemake.sh shock > Makefile
   $ make setup
   $ make
   $ ./phantomsetup myshock

To run the code, you need to write a slurm script. You can get an
example by typing “make qscript”:

::

   $ make qscript INFILE=myshock.in > run.q

should produce something like

::

   $ cat run.q
   #!/bin/bash
   #SBATCH --nodes=1 --ntasks=32
   #SBATCH --cpus-per-task=1
   #SBATCH --job-name=anesthesimeter
   #SBATCH --output=myshock.in.qout
   #SBATCH --mail-type=BEGIN
   #SBATCH --mail-type=FAIL
   #SBATCH --mail-type=END
   #SBATCH --mail-user=daniel.price@monash.edu
   #SBATCH --time=0-168:00:00
   #SBATCH --mem=16G
   echo "HOSTNAME = $HOSTNAME"
   echo "HOSTTYPE = $HOSTTYPE"
   echo Time is `date`
   echo Directory is `pwd`

   ulimit -s unlimited
   export OMP_SCHEDULE="dynamic"
   export OMP_NUM_THREADS=32
   export OMP_STACKSIZE=1024m


   echo "starting phantom run..."
   export outfile=`grep logfile "myshock.in" | sed "s/logfile =//g" | sed "s/\\!.*//g" | sed "s/\s//g"`
   echo "writing output to $outfile"
   ./phantom myshock.in >& $outfile

You can then submit this to the queue using

::

   $ sbatch run.q
   Submitted batch job 245936

and check status using

::

   $ squeue -u dprice
                JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
               245936   skylake phonator   dprice PD       0:00      1 (Resources)

splash on OzStar
~~~~~~~~~~~~~~~~

There is a version of splash you can get by loading the relevant module
(module load splash). If you want a more recent version there is a
version that gets regularly updated in the shared project folder
(/fred/oz015/splash):

::

   /fred/oz015/splash/bin/ssplash

You can add this directory in your path by putting the following lines
in your ~/.bashrc file:

::

   export PATH=/fred/oz015/splash/bin:${PATH}

more info
~~~~~~~~~

For more information on the actual machine `read the
userguide <https://supercomputing.swin.edu.au>`__

getting your job to run quickly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

check the online job monitor!
