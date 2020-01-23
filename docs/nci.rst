Getting started on the NCI supercomputer (Australian National Supercomputing Facility)
======================================================================================

Apply for an account at http://nci.org.au

If you are in Daniel Price’s research group, request to join project “fu7”

Log in 
-------

::

   ssh -Y USER@gadi.nci.org.au

Make a shortcut to the /short filesystem
----------------------------------------

::

   cd /scratch/fu7
   mkdir $USER
   cd
   ln -s /scratch/fu7/$USER runs
   cd runs
   pwd -P

Edit your ~/.bashrc file
------------------

Mine has:

::

   export SYSTEM=nci
   export OMP_STACKSIZE=512M
   export PATH=$PATH:$HOME/splash/bin
   export MAXP=2000000
   source ~/.modules

Put relevant modules in your .modules file
------------------------------------------

Mine contains:

::

   module load intel-compiler
   module load intel-mpi

more info
---------

For more information on the actual machine `read the
userguide <https://opus.nci.org.au/display/Help/Preparing+for+Gadi>`__
