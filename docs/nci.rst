Getting started on the raijin supercomputer
===========================================

Apply for an account at http://nci.org.au

If you are in Daniel Price’s research group, request to join projects
“fu7” and “pt4”

Log in to raijin
----------------

::

   ssh -Y USERNAME@raijin.nci.org.au

Make a shortcut to the /short filesystem
----------------------------------------

::

   cd /short/pt4
   mkdir USERNAME
   cd
   ln -s /short/pt4/USERNAME runs
   cd runs
   pwd -P

Edit your .profile
------------------

Mine has:

::

   export SYSTEM=raijin
   export OMP_STACKSIZE=512M
   export PATH=$PATH:$HOME/splash/bin
   export MAXP=2000000
   source ~/.modules

Put relevant modules in your .modules file
------------------------------------------

Mine contains:

::

   module load intel-fc
   module load intel-mpi

more info
---------

For more information on the actual machine `read the
userguide <https://opus.nci.org.au/display/Help/Raijin+User+Guide>`__
