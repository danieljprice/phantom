How to set up and run a common envelope binary simulation
=========================================================

Use SETUP=star:

::

   ~/phantom/scripts/writemake.sh star > Makefile

::

   make setup
   ./phantomsetup star

Will create star.in. Set the following options:

::

                   tmax =       4000.    ! end time
                  dtmax =         50.    ! time between dumps
   ...
                  alpha =       1.000    ! MINIMUM art. viscosity parameter (max = 1.0)
                 alphau =       1.000    ! art. conductivity parameter
   ...
                   damp =       0.030    ! artificial damping of velocities (if on, v=0 initially)
   ...
         ishock_heating =           0    ! shock heating (0=off, 1=on)
   ...
          icreate_sinks =           1    ! allow automatic sink particle creation
           rho_crit_cgs =         10.    ! density above which sink particles are created (g/cm^3)
                 r_crit =   1.000E+04    ! critical radius for point mass creation (no new sinks < r_crit from existing sink)
                  h_acc =       0.030    ! accretion radius for new sink particles
                  hsoft =       1.000    ! softening length for sink particles (Plummer)
                  f_acc =       1.000    ! particles < f_acc*h_acc accreted without checks

Run this to completion:

::

   make
   ./phantom star.in

Then compile phantommoddump:

::

   make moddump

and run it:

::

   ./phantommoddump star_00080 binary_00000.tmp 0.0

Now you’re ready to run the binary calculation. Use

::

   ./phantomsetup binary.in

to create a binary.in file, and then change the file as follows

::

               dumpfile =  binary_00000.tmp  ! dump file to start from
   ...
                   tmax =     200000.    ! end time
                  dtmax =         50.    ! time between dumps     **These are fairly arbitrary. Depends on desired simulation**
   ...
                  alpha =       0.100    ! MINIMUM art. viscosity parameter (max = 1.0)
                 alphau =       1.000    ! art. conductivity parameter
   ...
          icreate_sinks =           0    ! allow automatic sink particle creation
                  f_acc =       0.500    ! particles < f_acc*h_acc accreted without checks

Once changed, run

::

   ./phantom binary.in

to start simulation. If anybody reads this, just note that the above
changes to binary.in were made without any REAL knowledge of how this
might affect the end result. Do not try this at home.
