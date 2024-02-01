A brief guide to Phantom runtime options
========================================

Probably the main thing you need to understand when running simulations with Phantom is what the input options mean. Here is an attempt to explain the most important ones. Let's start with a typical Phantom input file:

::

   # Runtime options file for Phantom, written 13/12/2012 08:56:08.2
   # Options not present assume their default values
   # This file is updated automatically after a full dump

   # job name
                logfile =  myrun01.log   ! file to which output is directed
               dumpfile =  myrun_00000.tmp   ! dump file to start from

   # options controlling run time and input/output
                   tmax =      2.0000    ! end time
                  dtmax =      0.0100    ! time between dumps
                   nmax =     1000000    ! maximum number of timesteps (0=just get derivs and stop)
                   nout =          -1    ! number of steps between dumps (-ve=ignore)
              nmaxdumps =          -1    ! stop after n full dumps (-ve=ignore)
               twallmax =       00:00    ! maximum wall time (hh:mm, 00:00=ignore)
              nfulldump =          10    ! full dump every n dumps
               iverbose =           0    ! verboseness of output to log file (0-5)

   # options controlling accuracy
                 C_cour =      0.3000    ! Courant number
                C_force =      0.2500    ! dt_force number
                   tolv =   1.000E-02    ! tolerance on v iterations in timestepping
                  hfact =      1.2000    ! h in units of particle spacing [h = hfact(m/rho)^(1/3)]
                   tolh =   1.000E-04    ! tolerance on h-rho iterations
      restartonshortest =           F    ! restart with all particles on shortest timestep

   # options controlling hydrodynamics, artificial dissipation
                   ieos =           2    ! eqn of state (1=isoth; 2=adiab; 3/4=locally iso (sphere/cyl); 5=two phase)
                  alpha =      0.1000    ! MINIMUM art. viscosity parameter (max = 1.0)
                 alphau =      1.0000    ! art. conductivity parameter
                   beta =      2.0000    ! beta viscosity
           avdecayconst =      0.1000    ! decay time constant for viscosity switches
                   damp =      0.0000    ! artificial damping of velocities (if on, v=0 initially)
           ipdv_heating =           1    ! heating from PdV work (0=off, 1=on)
         ishock_heating =           1    ! shock heating (0=off, 1=on)

   # options relating to external forces
         iexternalforce =           0    ! external force to apply (1=1/r^2,2=1/rcyl^2,3=binary,4=prdrag,8=spiral,9=lt,10=ns)

   # options controlling physical viscosity
              irealvisc =           0    ! physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)
             shearparam =      0.1000    ! magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)
               bulkvisc =      0.0000    ! magnitude of bulk viscosity

Now, to deconstruct:

options controlling run time and input/output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

logfile
-------

::

                logfile =  myrun01.log   ! file to which output is directed

Specifies the name of the output log. Actually this is ignored as by default the output is directed to stdout.

dumpfile
--------

::

               dumpfile =  myrun_00000.tmp   ! dump file to start from

Specifies the name of the file containing particle positions etc. from which to start the calculation. The .tmp indicates that the file is a temporary file created by phantomsetup which has *not* had the density calculated. This means that the smoothing lengths, and hence densities, are just a guess.

tmax
----

::

                   tmax =      2.0000    ! end time

Maximum time in code units for the simulation, i.e. the code will stop when execution when t reaches this value.

dtmax
-----

::

                  dtmax =      0.0100    ! time between dumps

Interval (in code units of time) between snapshot files (“dumps”). Dump files can either be full dumps - containing everything needed to restart a calculation, or “small” dumps, containing just the basic information needed to perform a visualisation.

nmax
----

::

                   nmax =     1000000    ! maximum number of timesteps (0=just get derivs and stop)

The maximum number of timesteps that the code will take. The main use of this option is to set this to zero, meaning that the code will just compute the density and/or forces, replace and delete the .tmp file, and stop.

nout
----

::

                   nout =          -1    ! number of steps between dumps (-ve=ignore)

Can be used to specify that snapshots should be written after a set number of timesteps rather than at intervals of dtmax. Mainly useful for debugging, or if the timescales in the simulation are changing rapidly.  Best to set dtmax=tmax if this option is used, otherwise undefined behaviour may result.

nmaxdumps
---------

::

              nmaxdumps =          -1    ! stop after n full dumps (-ve=ignore)

Stop the code after completing n full dumps (i.e. nmaxdump intervals of dtmax*nfulldump). Can be used instead of specifying tmax, e.g. that we want 1000 dump files in the simulation.

twallmax
--------

::

               twallmax =       00:00    ! maximum wall time (hh:mm, 00:00=ignore)

Maximum wall time for the simulation. Set this if you are running in a queue with strict walltime limits. The code predicts whether or not a full dump will be written within the remaining walltime, and if not, stops execution, meaning that cpu time is not wasted evolving data that will never be written to disk.

nfulldump
---------

::

              nfulldump =          10    ! full dump every n dumps

The above would specify that every 10th dump would be a full (“restart”) dump. Use nfulldump=1 to make every dump file a full dump, meaning that the code can be restarted from every output file. Small dumps are mainly useful for saving disk space but retaining the ability to make a movie with reasonable time-resolution.

iverbose
--------

::

               iverbose =           0    ! verboseness of output to log file (0-5)

Specifies the verboseness level of the output to the log file.  iverbose=0 is the default, iverbose=1 gives some extra information, while iverbose=2 gives timing and neighbour statistics for every timestep.

options controlling accuracy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C_cour
------

::

                 C_cour =      0.3000    ! Courant number

Timestep constraint from Courant condition is :math:`dt < C_{cour} h / v_{sig}`. 0.3 is roughly the maximum for stability.

C_force
-------

::

                C_force =      0.2500    ! dt_force number

Timestep constraint based on force condition :math:`dt < C_{force} \sqrt(h/|f|)`. Does not usually constrain timestep.

tolv
----

::

                   tolv =   1.000E-02    ! tolerance on v iterations in timestepping

This is related to the treatment of velocity-dependent forces in the leapfrog integrator. The basic leapfrog scheme is as follows:

::

   x^1 = x^0 + dt*v^0 + dt/2 * f^0
   f^1 = f(x^1)
   v^1 = v_0 + dt/2 * (f^0 + f^1)

…which becomes implicit if the force depends on velocity, i.e. f^1 = f(x^1, v^1). In SPH the only velocity-dependent part of the force is the viscosity terms (external forces that depend on v are handled separately). To deal with this we make a prediction, i.e.:

::

   x^1 = x^0 + dt*v^0 + dt/2 * f^0
   v^* = v^0 + dt*f^0
   f^* = f(x^1, v^*)
   v^1 = v_0 + dt/2 * (f^0 + f^*)

then if v^1 and the prediction (v^*) differ by more than tolv, we iterate the force evaluation. Since this is a waste (it is as expensive as doing another timestep), the integrator also reduces the timestep such that this condition should not occur (this is when the timestep is controlled by the “error” condition). A value of tolv > 100 means that tolv will be ignored, i.e. no iterations will ever be taken.
hfact
-----

::

                  hfact =      1.2000    ! h in units of particle spacing [h = hfact(m/rho)^(1/3)]

Specifies the smoothing length in units of the particle spacing, and hence both the resolution length and the number of neighbours. The latter is proportional to hfact^3 in 3D. The formula relating hfact to the number of neighbours is given in Price (2012). Note that hfact is referred to as eta in that paper.

tolh
----

::

                   tolh =   1.000E-04    ! tolerance on h-rho iterations

Smoothing length and density are iterated self-consistently since they are mutually dependent until the conditio (h - h_prev)/h_0 < tolh is satisfied. Essentially, this is the tolerance to which the smoothing length and density are interchangeable, since only the smoothing length is stored in Phantom. Thus think of it as the error you are prepared to allow in the density (with higher errors in low density regions). In other words, don’t increase this, but decreasing it is OK and does not substantially increase the cost since density iterations are cheap in Phantom.

restartonshortest
-----------------

::

      restartonshortest =           F    ! restart with all particles on shortest timestep

just specifies that all particles should be put on the shortest timestep when the simulation restarts. Use this if the calculation has crashed for some reason and you want to try to bravely carry on.

options controlling hydrodynamics, artificial dissipation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ieos
----

::

                   ieos =           2    ! eqn of state (1=isoth; 2=adiab; 3/4=locally iso (sphere/cyl); 5=two phase)

specifies the equation of state. See eos.f90 for details.

alpha
-----

::

                  alpha =      0.1000    ! MINIMUM art. viscosity parameter (max = 1.0)

the main artificial viscosity parameter – either the fixed constant value if ialphaind=0 (specified in the dim file), or the minimum value if ialphaind=idim. Note that which of these has been compiled is indicated by the comment.

alphau
------

::

                 alphau =      1.0000    ! art. conductivity parameter

artificial conductivity parameter. Generally safe to leave this as 1.0 since we use the signal velocity for conductivity that depends on the velocity difference. Read Price (2008) for why artificial conductivity is important. Note that Phantom employs the vsig=|vab.r\| signal velocity (appropriate for simulations involving gravity), as noted in the footnote of that paper, in Price (2012), and in Wadsley et al.  (2008).

beta
----

::

                   beta =      2.0000    ! beta viscosity

beta viscosity parameter. Mainly important for preventing particle interpenetration. Important to set this higher (to 4.0) if the Mach number exceeds ~5. Note that beta should generally be ON even if no alpha viscosity is applied, otherwise particles can penetrate each other (e.g. going through the midplane of a disc).

avdecayconst
------------

::

           avdecayconst =      0.1000    ! decay time constant for viscosity switches

Constant in the decay timescale for the viscosity, conductivity and resistivity switches. Roughly gives inverse of how many smoothing lengths away from a shock that the viscosity will take to be damped away. No reason to change this particularly.

damp
----

::

                   damp =      0.0000    ! artificial damping of velocities (if on, v=0 initially)

specifies a damping of the form

::

       dv/dt = f - damp*v

i.e., should be a number between 0.0 and 1.0 specifying the fraction of the kinetic energy to be removed at every timestep. Use this when trying to relax the initial conditions of a simulation by damping the particles into a relaxed configuration. A good value when doing this is in the range 0.02-0.05.

ipdv_heating
------------

::

           ipdv_heating =           1    ! heating from PdV work (0=off, 1=on)

specifies whether or not the PdV work term du/dt = -P/rho*(div v) is included in the evolution of thermal energy. Not relevant if an isothermal equation of state is used.

ishock_heating
--------------

::

         ishock_heating =           1    ! shock heating (0=off, 1=on)

specifies whether or not the energy dissipation from the artificial viscosity term is included in the evolution of thermal energy (du/dt).  Not relevant if an isothermal equation of state is used.
