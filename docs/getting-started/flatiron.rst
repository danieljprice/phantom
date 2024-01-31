Getting started on Rusty (Flatiron cluster)
===========================================

See also :doc:`general instructions for running phantom on a remote cluster <clusters>`.

We assume you already have a Flatiron username and account

First time you log in
---------------------

Make sure you log in with the -Y flag to enable X-Windows forwarding::

   ssh -Y -p 61022 <USER>@gateway.flatironinstitute.org
   ssh -Y rusty

show available software::

   module avail

load intel compilers and splash
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   module load intel-oneapi-compilers/2023.0.0
   module use ~yjiang/modulefiles
   module load splash/3.8.3

get phantom
~~~~~~~~~~~~
Clone a copy of phantom into your home directory::

   git clone https://github.com/danieljprice/phantom.git

set your username and email address
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ensure that your name and email address are set, as follows::

   cd phantom
   git config --global user.name "Joe Bloggs"
   git config --global user.email "joe.bloggs@monash.edu"

Please use your full name in the format above, as this is what appears
in the commit logs (and in the AUTHORS file).

edit your .bashrc file
~~~~~~~~~~~~~~~~~~~~~~

I put the “module load” commands in a file called ~/.modules which
contains the modules I want every time I log in. For example::

   $ cat .modules
   module load intel-oneapi-compilers/2023.0.0

Then, add the following lines to your ~/.bashrc::

   source ~/.modules
   export SYSTEM=rusty
   ulimit -s unlimited
   export OMP_STACKSIZE=512M
   export OMP_SCHEDULE=dynamic

Now, when you login again, all these should be set automatically.

Performing a calculation
------------------------

You should *not* perform calculations in your home space - this is for
code and small files. Calculations should be run in the “ceph” area
in /mnt/ceph/users/$USER/

On other machines I usually make a soft link / shortcut called “runs” pointing to the directory where I want to run my calculations::

   $ cd
   $ ln -s /mnt/ceph/users/$USER runs
   $ cd runs
   $ pwd -P
   /mnt/ceph/users/USERNAME

However on the Flatiron machines there is already a shortcut called "ceph" in your homespace.

Then make a subdirectory for the name of the calculation you want to run
(e.g. shock)::

   mkdir shock
   cd shock
   ~/phantom/scripts/writemake.sh shock > Makefile
   make setup
   make
   ./phantomsetup shock

To run the code, you need to write a slurm script. You can get an
example by typing “make qscript”::

   make qscript NOMP=10 INFILE=shock.in > run.q

should produce something like::

   $ cat run.q
   #!/bin/bash
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=10
   #SBATCH --job-name=audiencia
   #SBATCH --partition=gen
   #SBATCH --output=shock.in.qout
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
   export OMP_NUM_THREADS=10
   export OMP_STACKSIZE=1024m

   echo "starting phantom run..."
   export outfile=`grep logfile "shock.in" | sed "s/logfile =//g" | sed "s/   \\!.*//g" | sed "s/\s//g"`
   echo "writing output to $outfile"
   ./phantom shock.in >& $outfile

You can then submit this to the "temp" queue using::

   $ sbatch -p temp run.q
   Submitted batch job 2547013

check status using::

   $ squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           2547013      temp audienci   dprice  R       0:02      1 worker3109

Once the job is running, follow the output log using the tail -f command::

   tail -f shock01.log

splash on rusty
----------------

There is a version of splash you can get by loading the relevant module::

   module use ~yjiang/modulefiles
   module load splash/3.8.3

Alternatively, you can also install a local copy of splash in your home space::

   cd
   git clone https://github.com/danieljprice/splash
   cd splash; git clone https://github.com/danieljprice/giza
   make withgiza SYSTEM=ifort

You can add this directory in your path by putting the following lines
in your ~/.bashrc file::

   export SPLASH_DIR=/mnt/home/dprice/splash
   export PATH=${SPLASH_DIR}/bin:${PATH}
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SPLASH_DIR}/giza/lib

visualising the shock problem
------------------------------
Simple plots with splash can be launched from the command line, e.g.::

  splash shock_0* -y density -x x

then press space to step through the snapshots. You can plot the exact solution
for the shock tube by switching it on from the o) menu::

  Please enter your selection now (y axis or option):o
  ------------- particle plot options -------------------
  8) plot exact solution                 (  3 )
  ...
  enter option ([0:9], default=0): 8
  ...
  3) shock tube
  ...
  enter exact solution to plot ([0:18], default=3):
  plotting exact solution number  3
  no file sod.shk
   >> read sod.setup
   rhoL, rho_R =    1.0000000000000000       0.12500000000000000
   pr_L, pr_R  =    1.0000000000000000       0.10000000000000001
   v_L,  v_R   =    0.0000000000000000        0.0000000000000000
  set adiabatic gamma manually? (no=read from dumps) (default=no):

then press 's' from the main menu to save, and remake your plot::

  Please enter your selection now (y axis or option):6
  (x axis) (default=1):
  (render) (0=none) ([0:16], default=0):
  Graphics device/type (? to see list, default /xw):

To make an offline version of the same plot from the command line, use::

  splash shock_0* -y density -x x -dev shock.pdf

read the `splash userguide <https://splash-viz.readthedocs.org>`__ for more

more interesting examples
-------------------------
To proceed to a more interesting calculation, just change the name of the :doc:`SETUP
parameter <setups>` when you created the Makefile in the run directory, as per
the :doc:`examples <examples>`::

  cd ~/ceph
  mkdir disc
  cd disc
  ~/phantom/scripts/writemake.sh disc > Makefile
  make setup; make

followed by::

  ./phantomsetup disc

after editing the .setup file as desired, then finish the setup and run the calculation::

  ./phantomsetup disc
  make qscript NOMP=10 INFILE=disc.in > run.q
  sbatch -p temp run.q

...and you're off and rolling.

more examples of test problems, e.g. Orszag-Tang vortex
----------------------------------------------------------
To run your own version of standardised test problems shown in the `code paper <http://adsabs.harvard.edu/abs/2018PASA...35...31P>`_,
you can clone the `phantom-examples <https://github.com/phantomSPH/phantom-examples>`_ repository::

  cd ~/ceph
  git clone https://github.com/phantomSPH/phantom-examples
  cd phantom-examples/mhd/orstang
  make setup; make

followed by::

  ./phantomsetup orstang256.setup
  make qscript NOMP=10 INFILE=orstang256.in > run.q
  sbatch -p temp run.q

more info
---------

For more information on the actual machine `read the
userguide <https://wiki.flatironinstitute.org/SCC/Hardware/Rusty>`__
