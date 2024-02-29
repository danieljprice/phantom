How to run Phantom on OzStar
============================

See also general instructions for :doc:`running phantom on a remote cluster <clusters>`.

Apply for an account
--------------------

https://supercomputing.swin.edu.au/account-management/new_account_request

If you are in Daniel Price’s research group, from your account
management page, request “join project”:

https://supercomputing.swin.edu.au/account-management/project_join_request

and select “oz015 - Price/Pinte research group”

First time you log in
---------------------

Replace USERNAME below with your username::

   $ ssh -Y USERNAME@ozstar.swin.edu.au

show available software::

   $ module avail

load intel compilers, git and splash. These might have different names
to the below, but should be similar::

   $ module load intel-compilers/2023.0.0
   $ module load ffmpeg/5.1.2
   $ module load gompi/2023a
   $ module load hdf5/1.14.0

Get phantom
~~~~~~~~~~~

Clone a copy of phantom into your home directory::

   $ git clone https://github.com/danieljprice/phantom.git

Set your username and email address
-----------------------------------

Ensure that your name and email address are set, as follows::

   cd phantom
   git config --global user.name "Joe Bloggs"
   git config --global user.email "joe.bloggs@monash.edu"

Please use your full name in the format above, as this is what appears
in the commit logs (and in the AUTHORS file).

edit your .bashrc file
----------------------

I put the “module load” commands in a file called ~/.modules which
contains the modules I want loaded every time I log in. For example::

   $ cat .modules
   module load intel-compilers/2023.0.0
   module load ffmpeg/5.1.2
   module load gompi/2023a
   module load hdf5/1.14.0

Then, add the following lines to your ~/.bashrc::

   source ~/.modules
   export SYSTEM=ozstar
   ulimit -s unlimited
   export OMP_STACKSIZE=512M

Now, when you login again, all these should be set automatically.

Performing a calculation
------------------------

You should *not* perform calculations in your home space - this is for
code and small files. Calculations should be run in the “project” area
in /fred/PROJECT_NAME/$USER

I usually make a soft link / shortcut called “runs” pointing to the
directory where I want to run my calculations::

   $ cd /fred/oz015
   $ mkdir $USER
   $ cd
   $ ln -s /fred/oz015/$USER runs
   $ cd runs
   $ pwd -P
   /fred/oz015/USERNAME

then make a subdirectory for the name of the calculation you want to run
(e.g. shock)::

   $ mkdir shock
   $ cd shock
   $ ~/phantom/scripts/writemake.sh shock > Makefile
   $ make setup
   $ make
   $ ./phantomsetup myshock

To run the code, you need to write a slurm script. You can get an
example by typing “make qscript”::

   $ make qscript INFILE=myshock.in > run.q

should produce something like::

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

You can then submit this to the queue using::

   $ sbatch run.q
   Submitted batch job 245936

and check status using::

   $ squeue -u dprice
                JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
               245936   skylake phonator   dprice PD       0:00      1 (Resources)

splash on OzStar
~~~~~~~~~~~~~~~~

There is a version of splash you can get by loading the relevant module
(module load splash). If you want a more recent version there is a
version that gets regularly updated in the shared project folder
(/fred/oz015/splash)::

   /fred/oz015/splash/bin/splash

You can add this directory in your path by putting the following lines
in your ~/.bashrc file::

   export PATH=/fred/oz015/splash/bin:${PATH}
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/fred/oz015/splash/giza/lib

getting your job to run quickly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

first check the `online job monitor <https://supercomputing.swin.edu.au/monitor/>`__, click on "Future" and check for available nodes with the largest number of cpus available (typically either 16 or 32). The "skylake" queue is the default. If you notice spare nodes on other queues, e.g. sstar or gstar you can request this queue via your job submission script, e.g.::

     #SBATCH --nodes=1 --ntasks=16
     ...
     #SBATCH --partition=sstar
     ...
     export OMP_NUM_THREADS=16

where as above you also need to adjust the number of cpus you are requesting to fit the node size. In the sstar queue, the default nodes have only 16 cpus: as the job can only run on one node, you need to either request 16 cpus in your job submission script as above, or request the single 32 core node in sstar using ::

     #SBATCH --nodes=1 --ntasks=32
     ...
     #SBATCH --partition=sstar
     #SBATCH -C largemem
     ...
     export OMP_NUM_THREADS=32     

getting your job to restart automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ozstar has fairly generous queue limits (168 hrs) but if you want your job to automatically continue longer than this (use with caution), you can submit a second job to the queue that depends on the first one completing, e.g.::

    $ sbatch run.q
    Submitted batch job 21300377
    
    $ sbatch --dependency=afterany:21300377 run.q
    Submitted batch job 21300378

You should then be able to see two jobs in the queue, with one waiting on the other to finish::

   $ squeue -u $USER
             21300378   skylake  enemata   dprice PD       0:00      1 (Dependency) 
             21300377   skylake  enemata   dprice R        0:10      1 (john110) 

You can use either "afterok" to start the next job only if the first job completed successfully, or "afterany" to restart when the previous job terminates for any reason.


more info
~~~~~~~~~

For more information on the actual machine `read the
userguide <https://supercomputing.swin.edu.au>`__

See also general instructions for :doc:`running phantom on a remote cluster <clusters>`.


