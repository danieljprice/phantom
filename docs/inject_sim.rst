
Injecting particles from existing simulations to new simulations
=========================================================

Initial setup
-------------

To ensure the particle mass and units are consistent in both existing & new simulations, 
   it is recommended to use 'phantommoddump' with existing simulations to setup new simulations

::

   make SRCINJECT=inject_sim.f90; make moddump SRCINJECT=inject_sim.f90
   ./phantommoddump YOUR_EXISTING_SIMULATION YOUR_NEW_SIMULATION TIME

'phantommodump' might produce a parameter file depending on the setup, 
   in that case one would need to run 
::
   
   ./phantommoddump YOUR_EXISTING_SIMULATION YOUR_NEW_SIMULATION TIME' 

one more time after setting up the parameters
 
At the end of these instructions, an initial dump of the new simulaton and a .in file are created.

::

Content of the .in file
--------------------------

Options controlling particle injection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   # options controlling particle injection
       start_dump =  'dump_00000'   ! dumpfile to start for injection
         r_inject =     5.000E+14   ! radius to inject tde outflow (in cm)
       final_dump =  'dump_02000'   ! stop injection after this dump

Here’s a brief description of each of them

::

       start_dump =  'dump_00000'   ! dumpfile to start for injection

set the dump start to inject. The code will check the start_dump time and start injection when the time is reached in new simulations
Once a dump is used by injection, the dump number will automatically increased by 1. The new dump is written to .in file once a full dump is saved

If the dumps are in a different directory, 
 
::

       start_dump =  'PATH/TO/YOUR/OTHER/DIR/dump_00000'   ! dumpfile to start for injection

can read dumps from other directory. The path needs to be the RELATIVE path to the working directory
!!!--------------------------------------!!!
NOTE: qotation marks are NECESSARY with path
!!!--------------------------------------!!!

::

         r_inject =     5.000E+14   ! radius to inject tde outflow (in cm)

set the radius for inject. For TDE outflow specifically, once a particle pass this radius from inside to outside in the existing simulations, it is injected to the new simulations

::

       final_dump =  'dump_02000'   ! stop injection after this dump

set the dump to stop injection. The injection dump number keep increasing by 1 after each injection and will stop once reaching this set final_dump.
If there is a PATH in start_dump, it is NECESSARY in final_dump as well.

